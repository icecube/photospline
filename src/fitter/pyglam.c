#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>

#include "glam.h"

/* Sparse GLAM Python interface */

static PyObject *pybox(PyObject *self, PyObject *args);
static PyObject *pyrho(PyObject *self, PyObject *args);
static PyObject *pyfit(PyObject *self, PyObject *args, PyObject *kw);
static PyObject *pygrideval(PyObject *self, PyObject *args);
static PyObject *pynnls(PyObject *self, PyObject *args, PyObject *kw);

static cholmod_sparse *construct_penalty(struct splinetable* out,
    PyObject* py_penalty, double smooth, int monodim, cholmod_common* c);

static PyObject *splinetable_mod;

PyDoc_STRVAR(pynnls__doc__, "nnls(AtA, Atb, algorithm=\"block\") => x");

static PyMethodDef methods[] = {
  { "box", pybox, METH_VARARGS },
  { "rho", pyrho, METH_VARARGS },
  { "fit", (PyObject *(*)(PyObject *, PyObject *))pyfit,
      METH_VARARGS | METH_KEYWORDS },
  { "grideval", pygrideval, METH_VARARGS },
  { "nnls", (PyObject *(*)(PyObject *, PyObject *))pynnls,
      METH_VARARGS | METH_KEYWORDS, pynnls__doc__ },
  { NULL, NULL }
};

void initspglam(void)
{
  import_array();

  /* Look up the splinetable class */
  splinetable_mod =
      PyImport_ImportModule("photospline.splinetable");

  if (splinetable_mod == NULL) {
    PyErr_SetString(PyExc_ImportError,
        "Could not import splinetable module");
    return;
  }
  Py_InitModule("spglam", methods);
}

static cholmod_sparse *
numpy2d_to_sparse(PyArrayObject *a, cholmod_common *c)
{
  cholmod_dense ad;
  cholmod_sparse *sp, *spt;

  ad.nrow = a->dimensions[1];
  ad.ncol = a->dimensions[0];
  ad.nzmax = ad.nrow * ad.ncol;
  ad.d = ad.nrow;
  ad.x = a->data;
  ad.z = NULL;
  ad.xtype = CHOLMOD_REAL;
  ad.dtype = CHOLMOD_DOUBLE;

  sp = cholmod_l_dense_to_sparse(&ad, 1, c);
  
  /* Correct for row-major/column-major ordering issues */
  spt = cholmod_l_transpose(sp, 1, c);
  cholmod_l_free_sparse(&sp, c);

  return spt;
}

static PyArrayObject *
numpy_sparse_to_2d(cholmod_sparse *a, cholmod_common *c)
{
  npy_intp dimensions[2];
  PyArrayObject *out;
  cholmod_dense *ad;
  cholmod_sparse *at;

  dimensions[0] = a->nrow;
  dimensions[1] = a->ncol;
  out = (PyArrayObject *)PyArray_SimpleNew(2, dimensions,
      PyArray_DOUBLE);

  at = cholmod_l_transpose(a, 1, c); /* Fix row-major/column-major */
  ad = cholmod_l_sparse_to_dense(at, c);
  cholmod_l_free_sparse(&at, c);
  memcpy(out->data, ad->x, sizeof(double) * ad->nrow * ad->ncol);
  cholmod_l_free_dense(&ad, c);

  return out;
}

static int
numpynd_to_ndsparse(PyObject *in, struct ndsparse *out)
{
  PyArrayObject *inar;
  int i, j, elements, coord, currow;
  unsigned int *moduli;

  inar = (PyArrayObject *)PyArray_ContiguousFromObject(in,
      PyArray_DOUBLE, 1, INT_MAX);
  if (inar == NULL)
    return -1;

  out->rows = 0;
  out->ndim = inar->nd;
  out->ranges = malloc(sizeof(unsigned)*out->ndim);
  out->i = malloc(sizeof(int *)*out->ndim);
  moduli = malloc(sizeof(unsigned)*out->ndim);

  elements = 1;
  for (i = 0; i < inar->nd; i++) {
    out->ranges[i] = inar->dimensions[i];
    elements *= inar->dimensions[i];
  }

  moduli[out->ndim-1] = 1;
  for (i = out->ndim-2; i >= 0; i--)
    moduli[i] = moduli[i+1]*out->ranges[i+1];
  
  /* Find how many non-zeros we have */
  for (i = 0; i < elements; i++) {
    if (((double *)(inar->data))[i] != 0.0)
      out->rows++;
  }
  assert(out->rows>0);

  /* Init coordinates */
  for (i = 0; i < inar->nd; i++)
    out->i[i] = malloc(sizeof(int)*out->rows);
  out->x = malloc(sizeof(double)*out->rows);

  currow = 0;
  for (i = 0; i < elements; i++)  {
    if (((double *)(inar->data))[i] == 0)
      continue;
    out->x[currow] = ((double *)(inar->data))[i];

    coord = i;
    for (j = 0; j < inar->nd; j++) {
      out->i[j][currow] = coord / moduli[j];
      coord = coord % moduli[j];
    }

    currow++;
  }
  Py_DECREF(inar);
  free(moduli);

  return 0;
}

static PyArrayObject *
numpy_ndsparse_to_ndarray(struct ndsparse *a)
{
  double *x;
  PyArrayObject *out;
  unsigned int moduli[a->ndim];
  npy_intp dimensions[a->ndim];
  int i, j, k, elements;

  /* Change the type of a->ranges to pacify numpy */
  for (i = 0; i < a->ndim; i++)
    dimensions[i] = a->ranges[i];

  out = (PyArrayObject *)PyArray_SimpleNew(a->ndim, dimensions,
      PyArray_DOUBLE);

  moduli[a->ndim-1] = 1;
  for (i = a->ndim-2; i >= 0; i--)
    moduli[i] = moduli[i+1]*a->ranges[i+1];

  /* Find and initialize the data array to zero */
  x = (double *)out->data;
  elements = 1;
  for (i = 0; i < a->ndim; i++) 
    elements *= a->ranges[i];
  memset(x, 0, sizeof(double)*elements);

  for (i = 0; i < a->rows; i++) {
    k = 0;
    for (j = 0; j < a->ndim; j++)
      k += a->i[j][i]*moduli[j];
    
    x[k] = a->x[i];
  }

  return out;
}

static PyObject *pybox(PyObject *self, PyObject *args)
{
  PyObject *xa, *xb;
  PyArrayObject *a, *b, *result_array;
  cholmod_common c;
  cholmod_sparse *am, *bm, *result;

  if (!PyArg_ParseTuple(args, "OO", &xa, &xb))
    return NULL;

  a = (PyArrayObject *)PyArray_ContiguousFromObject(xa,
      PyArray_DOUBLE, 2, 2);
  b = (PyArrayObject *)PyArray_ContiguousFromObject(xb,
      PyArray_DOUBLE, 2, 2);

  if (a == NULL || b == NULL) {
    PyErr_SetString(PyExc_ValueError,
        "arrays must be two-dimensional and of type float");
    return NULL;
  }

  cholmod_l_start(&c);

  am = numpy2d_to_sparse(a, &c);
  bm = numpy2d_to_sparse(b, &c);

  Py_DECREF(a);
  Py_DECREF(b);

  result = box(am, bm, &c);
  cholmod_l_free_sparse(&am, &c);
  cholmod_l_free_sparse(&bm, &c);

  if (result == NULL) {
    PyErr_SetString(PyExc_ValueError,
        "arrays must have compatible dimensions");
    cholmod_l_finish(&c);
    return NULL;
  }

  result_array = numpy_sparse_to_2d(result, &c);

  cholmod_l_free_sparse(&result, &c);
  cholmod_l_finish(&c);

  return PyArray_Return(result_array);
}

void print_ndsparse_py(struct ndsparse *a) {
  PyObject *py = (PyObject *)numpy_ndsparse_to_ndarray(a);
  PyObject_Print(py, stdout, 0);
  Py_DECREF(py);
}

void printndsparse(struct ndsparse *a) {
  double *x;
  int moduli[a->ndim];
  int i, j, k, elements;

  elements = 1;
  for (i = 0; i < a->ndim; i++) 
    elements *= a->ranges[i];
  moduli[a->ndim-1] = 1;
  for (i = a->ndim-2; i >= 0; i--)
    moduli[i] = moduli[i+1]*a->ranges[i+1];

  for (i = 0; i < a->ndim; i++)
    printf("Dimension %d: %d\n",i,a->ranges[i]);

  x = calloc(elements,sizeof(double));

  for (i = 0; i < a->rows; i++) {
    k = 0;
    for (j = 0; j < a->ndim; j++)
      k += a->i[j][i]*moduli[j];
    
    x[k] = a->x[i];
  }

  for (i = 0; i < elements; i++) {
    if (a->ndim > 0 && i % moduli[0] == 0 && i != 0) printf("\n");
    if (a->ndim > 1 && i % moduli[1] == 0 && i != 0) printf("\n");
    printf("%lf\t",x[i]);
  }
  printf("\n");

  free(x);
}

static PyObject *pyrho(PyObject *self, PyObject *args)
{
  /* This takes the reverse calling order in Python:
   * (matrix, ndarray, dim) */

  PyObject *xa, *xb;
  PyArrayObject *a, *result_array;
  struct ndsparse nd;
  cholmod_sparse *am;
  cholmod_common c;
  int dim, i, err;

  if (!PyArg_ParseTuple(args, "OOi", &xa, &xb, &dim))
    return NULL;

  a = (PyArrayObject *)PyArray_ContiguousFromObject(xa,
      PyArray_DOUBLE, 2, 2);
  if (a == NULL || numpynd_to_ndsparse(xb, &nd) != 0) {
    if (a != NULL) {
      Py_DECREF(a);
    }
  
    PyErr_SetString(PyExc_ValueError,
        "could not decode arrays");
    return NULL;
  }

  cholmod_l_start(&c);

  am = numpy2d_to_sparse(a, &c);
  err = slicemultiply(&nd, am, dim, &c);

  cholmod_l_free_sparse(&am, &c);
  cholmod_l_finish(&c);

  if (err == 0) {
    result_array = numpy_ndsparse_to_ndarray(&nd);
  } else {
    PyErr_SetString(PyExc_ValueError,
        "Dimensions do not match");
    result_array = NULL;
  }

  for (i = 0; i < nd.ndim; i++)
    free(nd.i[i]);
  free(nd.x); free(nd.ranges);

  return (PyObject *)result_array;
}

static PyObject *pyfit(PyObject *self, PyObject *args, PyObject *kw)
{
  PyObject *z, *w, *coords, *knots, *periods, *order,
      *penorder_py, *result;
  PyArrayObject *result_arr, *z_arr;
  struct ndsparse data;
  struct splinetable out;
  double *data_arr, *weights;
  double **c_coords;
  cholmod_common c;
  int *penorder = NULL;
  int *moduli;
  npy_intp* naxes;
  double smooth;
  int i, j, k, elements, err;
  int monodim = -1;
  char *keywordargs[] = {"z", "w", "coords", "knots", "order", "smooth",
      "periods", "penalties", "monodim", NULL};

  /* Initialize a few things to NULL */
  moduli = NULL;
  weights = NULL;
  z_arr = NULL;
  c_coords = NULL;
  result = NULL;
  periods = NULL;
  penorder_py = NULL;
  smooth = 1;
  memset(&out, 0, sizeof(out));

  /* Parse our arguments from Python land */
  if (!PyArg_ParseTupleAndKeywords(args, kw, "OOOOO|dOOi", keywordargs,
      &z, &w, &coords, &knots, &order, &smooth, &periods, &penorder_py,
      &monodim))
    return NULL;

  /* Parse weights first to avoid storing data with 0 weight */
  err = numpynd_to_ndsparse(w, &data);
  if (err == 0)
    z_arr = (PyArrayObject *)PyArray_ContiguousFromObject(z,
        PyArray_DOUBLE, data.ndim, data.ndim);

  if (err != 0 || z_arr == NULL) {
    PyErr_SetString(PyExc_ValueError,
        "could not decode arrays");
    return NULL;
  }

  /* Now validate the input */

  result = NULL;
  for (i = 0; i < data.ndim; i++) {
    if (z_arr->dimensions[i] != data.ranges[i]) {
      Py_DECREF(z_arr);
      PyErr_SetString(PyExc_ValueError,
          "weight and data array dimensions do not match");
      goto exit;
    }
  }
  if (data.ndim == 0) {
    Py_DECREF(z_arr);
	PyErr_SetString(PyExc_ValueError,
	  "data must be at least one-dimensional");
	goto exit;
  }

  /* Set up the data array */
  moduli = malloc(sizeof(int) * data.ndim);
  moduli[data.ndim-1] = 1;
  for (i = data.ndim-2; i >= 0; i--)
    moduli[i] = moduli[i+1]*data.ranges[i+1];

  data_arr = calloc(data.rows,sizeof(double));
  for (i = 0; i < data.rows; i++) {
    k = 0;
    for (j = 0; j < data.ndim; j++)
      k += moduli[j]*data.i[j][i];
    data_arr[i] = ((double *)(z_arr->data))[k];
  }
  free(moduli);

  /* Swap the data into the ndsparse structure to satisfy glamfit's
   * calling convention */
  weights = data.x;
  data.x = data_arr;

  /* We don't need the Python structure anymore */
  Py_DECREF(z_arr);

  /* Check knot and coords for consistency */
  if (!PySequence_Check(knots) || data.ndim != PySequence_Length(knots)) {
    PyErr_SetString(PyExc_TypeError,
        "knots must be a sequence with one row for each dimension");
    goto exit;
  }
  if (!PySequence_Check(coords) || data.ndim != PySequence_Length(coords)) {
    PyErr_SetString(PyExc_TypeError,
        "coord must be a sequence with one row for each dimension");
    goto exit;
  }

  /* Start setting up the spline table */
  out.ndim = data.ndim;

  out.knots = calloc(out.ndim,sizeof(double *));
  out.nknots = calloc(out.ndim,sizeof(long));
  for (i = 0; i < PySequence_Length(knots); i++) {
    PyArrayObject *knot_vec;
    knot_vec = (PyArrayObject *)PyArray_ContiguousFromObject(
        PySequence_GetItem(knots, i),
        PyArray_DOUBLE, 1, 1);

    if (knot_vec == NULL) {
      PyErr_SetString(PyExc_TypeError,
          "knots cannot be read as arrays");
      goto exit;
    }

    out.nknots[i] = knot_vec->dimensions[0];
    out.knots[i] = calloc(out.nknots[i], sizeof(double));
    memcpy(out.knots[i], knot_vec->data,
        out.nknots[i] * sizeof(double));

    Py_DECREF(knot_vec);
  }

  out.order = calloc(out.ndim,sizeof(int));
  penorder = calloc(out.ndim,sizeof(int));
  penorder[0] = 2;
  {

    if (PySequence_Check(order)) {
      for (i = 0; i < out.ndim; i++)
        out.order[i] = PyLong_AsLong(PySequence_GetItem(
            order,i));
    } else {
      out.order[0] = PyLong_AsLong(order);
      for (i = 1; i < out.ndim; i++)
        out.order[i] = out.order[0];
    }
  }
  
  c_coords = malloc(sizeof(double *)*out.ndim);
  for (i = 0; i < out.ndim; i++)
    c_coords[i] = NULL;
  for (i = 0; i < PySequence_Length(coords); i++) {
    PyArrayObject *coord_vec;
    coord_vec = (PyArrayObject *)PyArray_ContiguousFromObject(
        PySequence_GetItem(coords, i),
        PyArray_DOUBLE, 1, 1);

    if (coord_vec == NULL) {
      PyErr_SetString(PyExc_TypeError,
          "coords cannot be read as arrays");
      goto exit;
    }

    if (coord_vec->dimensions[0] != data.ranges[i]) {
      PyErr_SetString(PyExc_ValueError,
          "wrong number of coords");
      Py_DECREF(coord_vec);
      goto exit;
    }

    c_coords[i] = calloc(data.ranges[i], sizeof(double));
    memcpy(c_coords[i], coord_vec->data,
        data.ranges[i] * sizeof(double));

    Py_DECREF(coord_vec);
  }

  /* Do the fit */
  cholmod_l_start(&c);
  
  if (penorder_py == NULL) {
    /* do a fit with default penalties, weighted with smoothness */
    glamfit(&data, weights, c_coords, &out, smooth, out.order,
        penorder, monodim, 1, &c);
  } else {
    /* do a fit with arbitrary linear combinations of penalties */
    cholmod_sparse* penalty = construct_penalty(&out, penorder_py,
        smooth, monodim, &c);
    glamfit_complex(&data, weights, c_coords, &out, out.order,
        penalty, monodim, 1, &c);
          cholmod_l_free_sparse(&penalty, &c);
  }
  cholmod_l_finish(&c);

  /* Now process the splinetable into a numpy array */
  elements = 1;
  naxes = malloc(sizeof(npy_intp) * out.ndim);
  for (i = 0; i < out.ndim; i++) {
    elements *= out.naxes[i];
    naxes[i] = out.naxes[i];
  }
  result_arr = (PyArrayObject *)PyArray_SimpleNew(out.ndim, naxes,
      PyArray_FLOAT);

  if (result_arr != NULL) {
    PyObject *splinetable_cls;

    memcpy(result_arr->data, out.coefficients,
        elements*sizeof(float));

    /* Look up the splinetable class */
    splinetable_cls = PyObject_GetAttrString(splinetable_mod,
        "SplineTable");
    if (splinetable_cls == NULL) {
      PyErr_SetString(PyExc_ImportError,
          "Could not find spline table class");
      goto exit;
    }
    result = PyInstance_New(splinetable_cls, NULL, NULL);
    if (result == NULL) {
      PyErr_SetString(PyExc_ImportError,
          "Could not instantiate spline table class");
      goto exit;
    }

    PyObject_SetAttrString(result, "order", order);
    PyObject_SetAttrString(result, "knots", knots);
    if (periods == NULL) {
      npy_intp dim;

      dim = out.ndim;
      periods = PyArray_ZEROS(1, &dim, PyArray_INT, 0);
      PyObject_SetAttrString(result, "periods", periods);
      Py_DECREF(periods);
    } else {
      PyObject_SetAttrString(result, "periods", periods);
    }
    PyObject_SetAttrString(result, "coefficients",
        (PyObject *)result_arr);

    Py_DECREF(result_arr);
    Py_DECREF(splinetable_cls);
  }

exit:
  for (i = 0; i < data.ndim; i++)
    free(data.i[i]);
  free(data.x); free(data.ranges);
  if (weights) free(weights);
  if (c_coords != NULL) {
    for (i = 0; i < data.ndim; i++)
      if (c_coords[i])
	    free(c_coords[i]);
    free(c_coords);
  }
  if (penorder)
    free(penorder);
  if (out.periods)
    free(out.periods);
  if (out.knots) {
    for (i = 0; i < out.ndim; i++)
      free(out.knots[i]);
    free(out.knots);
  }
  if (out.nknots)
    free(out.nknots); 
  if (out.order)
    free(out.order);
  if (out.coefficients) {
    free(out.naxes);
    free(out.coefficients);
  }

  return (PyObject *)result;
}

/* Construct a penalty term from a dictionary of penalty orders */
static cholmod_sparse*
construct_penalty(struct splinetable* out, PyObject* py_penalty,
    double smooth, int monodim, cholmod_common* c)
{
  long *nsplines;
  cholmod_sparse *penalty;
  size_t sidelen;
  int i,j;
  long order;
  PyObject *key,*value,*py_scale;
  Py_ssize_t ppos = 0;
  double scale;
  
  nsplines = calloc(out->ndim,sizeof(long));
  
  sidelen = 1;
  for (i = 0; i < out->ndim; i++) {
    nsplines[i] = out->nknots[i] - out->order[i] - 1;
    sidelen *= nsplines[i];
  }
  
  penalty = cholmod_l_spzeros(sidelen, sidelen, 1, CHOLMOD_REAL,
      c);

  while (PyDict_Next(py_penalty, &ppos, &key, &value)) {
    order = PyInt_AsLong(key);
    if (PySequence_Check(value)) {
      for (j = 0; j < out->ndim; j++) {
        py_scale = PySequence_GetItem(value,j);
        scale = PyFloat_AsDouble(py_scale);
        Py_DECREF(py_scale);
        penalty = add_penalty_term(nsplines,
            out->knots[j], out->ndim, j,
            out->order[j], order, scale,
            j == monodim, penalty, c);
      }
    } else {
      scale = (py_penalty == Py_None) ? smooth :
          PyFloat_AsDouble(value);
      for (j = 0; j < out->ndim; j++) {
        penalty = add_penalty_term(nsplines,
            out->knots[j], out->ndim, j,
            out->order[j], order, scale,
            j == monodim, penalty, c);
      }  
    }
  }
  
  free(nsplines);
  return penalty;
}

static PyObject *
pygrideval(PyObject *self, PyObject *args)
{
  PyObject *table, *coords, *knots, *order_obj, *coeff;
  PyArrayObject *result;
  struct ndsparse nd;
  cholmod_common c;
  long order;
  int i;
  
  if (!PyArg_ParseTuple(args, "OO", &table, &coords))
    return NULL;

  knots = PyObject_GetAttrString(table, "knots");
  if (!PySequence_Check(knots) || !PySequence_Check(coords)) {
    PyErr_SetString(PyExc_TypeError,
      "Knots or coords not a sequence");
    return NULL;
  }

  order_obj = PyObject_GetAttrString(table, "order");

  coeff = PyObject_GetAttrString(table, "coefficients");
  numpynd_to_ndsparse(coeff, &nd);
  Py_DECREF(coeff);

  cholmod_l_start(&c);

  for (i = 0; i < PySequence_Length(knots); i++) {
    PyArrayObject *coord_vec, *knots_vec;
    cholmod_sparse *basis, *basist;

    coord_vec = (PyArrayObject *)PyArray_ContiguousFromObject(
        PySequence_GetItem(coords, i),
        PyArray_DOUBLE, 1, 1);
    knots_vec = (PyArrayObject *)PyArray_ContiguousFromObject(
        PySequence_GetItem(knots, i),
        PyArray_DOUBLE, 1, 1);
    if (PySequence_Check(order_obj))
      order = PyLong_AsLong(PySequence_GetItem(order_obj,i));
    else
      order = PyLong_AsLong(order_obj);
    
    basis = bsplinebasis((double *)knots_vec->data,
        knots_vec->dimensions[0], (double *)coord_vec->data,
        coord_vec->dimensions[0], order, &c);
    basist = cholmod_l_transpose(basis, 1, &c);
    cholmod_l_free_sparse(&basis, &c);

    Py_DECREF(coord_vec);
    Py_DECREF(knots_vec);

    slicemultiply(&nd, basist, i, &c);
  
    cholmod_l_free_sparse(&basist, &c);
  }
  Py_DECREF(order_obj);

  cholmod_l_finish(&c);

  result = numpy_ndsparse_to_ndarray(&nd);
  for (i = 0; i < nd.ndim; i++)
    free(nd.i[i]);
  free(nd.i);
  free(nd.x);
  free(nd.ranges);

  return ((PyObject *)result);
}

static PyObject *pynnls(PyObject *self, PyObject *args, PyObject *kw)
{
  PyObject *xa, *xb;
  char *algorithm = "block";
  enum nnls_algorithm {BLOCK, BLOCK_UPDOWN, BLOCK3} alg;
  PyArrayObject *a, *b, *result_array;
  npy_intp dimensions[2];
  cholmod_common c;
  cholmod_sparse *am;
  cholmod_dense *bm, *result;

  char *keywordargs[] = {"A", "y", "nnls", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kw, "OO|s", keywordargs, &xa, &xb, &algorithm))
    return NULL;

  if (!strcmp(algorithm, "block3"))
    alg = BLOCK3;
  else if (!strcmp(algorithm, "block"))
    alg = BLOCK;
  else if (!strcmp(algorithm, "block_updown"))
    alg = BLOCK_UPDOWN;
  else {
    PyErr_SetString(PyExc_ValueError,
        "algorithm must be one of 'block', 'block_updown', 'block3'");
    return (NULL);
  }

  a = (PyArrayObject *)PyArray_ContiguousFromObject(xa,
      PyArray_DOUBLE, 2, 2);
  b = (PyArrayObject *)PyArray_ContiguousFromObject(xb,
      PyArray_DOUBLE, 1, 1);

  if (a == NULL || b == NULL) {
    PyErr_SetString(PyExc_ValueError,
        "arrays must be two-dimensional and of type float");
    return NULL;
  }

  cholmod_l_start(&c);

  am = numpy2d_to_sparse(a, &c);
  bm = cholmod_l_allocate_dense(b->dimensions[0], 1, b->dimensions[0],
      CHOLMOD_REAL, &c);
  memcpy(bm->x, b->data, b->dimensions[0]*sizeof(double));

  Py_DECREF(a);
  Py_DECREF(b);

  switch (alg) {
    case BLOCK:
      result = nnls_normal_block(am, bm, 1, &c);
      break;
    case BLOCK_UPDOWN:
      result = nnls_normal_block_updown(am, bm, 1, &c);
      break;
    case BLOCK3:
      result = nnls_normal_block3(am, bm, 1, &c);
      break;
  }
  cholmod_l_free_dense(&bm, &c);

  if (result == NULL) {
    PyErr_SetString(PyExc_ValueError,
        "arrays must have compatible dimensions");
    cholmod_l_finish(&c);
    return NULL;
  }

  dimensions[0] = am->nrow;
  cholmod_l_free_sparse(&am, &c);

  result_array = (PyArrayObject *)PyArray_SimpleNew(1, dimensions,
      PyArray_DOUBLE);
  memcpy(result_array->data, result->x, dimensions[0]*sizeof(double));

  cholmod_l_free_dense(&result, &c);
  cholmod_l_finish(&c);

  return PyArray_Return(result_array);
}

