#include <Python.h>
#ifdef HAVE_NUMPY
	#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
	#include <numpy/ndarrayobject.h>
#endif
#include <stdio.h>
#include <limits.h>

#include <photospline/cinter/splinetable.h>
//#include <photospline/glam.h>

typedef struct{
	PyObject_HEAD
	struct splinetable table;
} pysplinetable;

static void
pysplinetable_dealloc(pysplinetable* self){
    splinetable_free(&self->table);
}

static PyObject*
pysplinetable_new(PyTypeObject* type, PyObject* args, PyObject* kwds){
	pysplinetable* self;

	self = (pysplinetable*)type->tp_alloc(type, 0);
	if(self)
		splinetable_init(&self->table);

	return (PyObject*)self;
}

static int
pysplinetable_init(pysplinetable* self, PyObject* args, PyObject* kwds){
	//base initialization
	splinetable_init(&self->table);
	
    static char* kwlist[] = {"path", NULL};
	char* path=NULL;
	
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s", kwlist, &path))
        return -1;
	
	if(path){ //try to read from the specified file
		int status=readsplinefitstable(path, &self->table);
		if(status!=0){
			PyErr_SetString(PyExc_IOError, "Unable to read spline from input file");
			return(-2);
		}
	}
	
    return 0;
}

static int
pysplinetable_print(pysplinetable* self, FILE* fp, int flags){
	uint32_t ndim=splinetable_ndim(&self->table);
	fprintf(fp,"Splinetable with %u dimension",ndim);
	if(ndim!=1)
		fprintf(fp,"s");
	return(0);
}

static PyObject*
pysplinetable_write(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"path", NULL};
	
	char* path=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &path))
		return(NULL);
	
	int status=writesplinefitstable(path,&self->table);
	if(status!=0){
		PyErr_SetString(PyExc_IOError, "Unable to write spline to output file");
		return(NULL);
	}
	
	return(Py_None);
}

static PyObject*
pysplinetable_get_aux_value(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"key", NULL};
	
	char* key=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &key))
		return(NULL);
	
	const char* value=splinetable_get_key(&self->table, key);
	if(value==NULL){
		PyErr_SetString(PyExc_KeyError, "Key not found");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("s",value);
	return(result);
}

//OMIT: read_key (users can do their own casting/parsing in python)
//TODO: write_key

#ifdef HAVE_NUMPY
static PyObject*
pysplinetable_getknots(pysplinetable *self, void *closure)
{
	PyObject *list = PyTuple_New(splinetable_ndim(&self->table));
	for (int dim = 0; dim < splinetable_ndim(&self->table); dim++) {
		npy_intp nknots = splinetable_nknots(&self->table, dim);
		PyObject *knots = PyArray_SimpleNewFromData(1, &nknots, NPY_DOUBLE, (void*)splinetable_knots(&self->table, dim));
		PyArray_CLEARFLAGS((PyArrayObject*)knots, NPY_ARRAY_WRITEABLE);
		PyTuple_SetItem(list, dim, knots);
	}
	
	return list;
}
#endif

static PyObject*
pysplinetable_getorder(pysplinetable* self, void *closure){
	PyObject *list = PyTuple_New(splinetable_ndim(&self->table));
	for (int dim = 0; dim < splinetable_ndim(&self->table); dim++) {
		PyTuple_SetItem(list, dim, PyInt_FromLong(splinetable_order(&self->table, dim)));
	}
	
	return list;
}

static PyObject*
pysplinetable_getextents(pysplinetable* self, void *closure){
	PyObject *list = PyTuple_New(splinetable_ndim(&self->table));
	for (int dim = 0; dim < splinetable_ndim(&self->table); dim++) {
		PyTuple_SetItem(list, dim, PyTuple_Pack(2,
		    PyFloat_FromDouble(splinetable_lower_extent(&self->table,dim)),
		    PyFloat_FromDouble(splinetable_upper_extent(&self->table,dim))));
	}
	
	return list;
}

#ifdef HAVE_NUMPY
static PyObject*
pysplinetable_getcoeffcients(pysplinetable* self, void *closure){
	unsigned int ndim = splinetable_ndim(&self->table);
	assert (ndim > 0);
	npy_intp dims[ndim];
	npy_intp strides[ndim];
	for (unsigned int i = 0; i < ndim; i++) {
		dims[i] = splinetable_ncoeffs(&self->table, i);
		strides[i] = sizeof(float)*splinetable_stride(&self->table, i);
	}
	return PyArray_New(&PyArray_Type, ndim, dims, NPY_FLOAT, strides,
	    (void*)splinetable_coefficients(&self->table), sizeof(float),
	    NPY_ARRAY_CARRAY_RO, NULL);
}
#endif

static PyObject*
pysplinetable_getndim(pysplinetable* self, void *closure){
	return PyInt_FromLong(splinetable_ndim(&self->table));
}

static PyObject*
pysplinetable_searchcenters(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"x", NULL};
	
	PyObject* pyx=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &pyx))
		return(NULL);
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	uint32_t ndim=splinetable_ndim(&self->table);
	if(xlen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of x must match the table dimension");
		return(NULL);
	}
	
	//unpack x
	//assume these are arbitrary sequences, not numpy arrays
	//TODO: should be possible to make a more efficient case for numpy arrays
	{
		//a small amount of evil
		double x[ndim];
		int centers[ndim];
		
		for(unsigned int i=0; i!=ndim; i++){
			PyObject* xi=PySequence_GetItem(pyx,i);
			x[i]=PyFloat_AsDouble(xi);
			Py_DECREF(xi); //done with this
			//printf("x[%u]=%lf\n",i,x[i]);
		}
		
		int status=tablesearchcenters(&self->table,x,centers);
		if(!status){
			PyErr_SetString(PyExc_ValueError, "tablesearchcenters failed");
			return(NULL);
		}
		
		PyObject* result=PyTuple_New(ndim);
		for(unsigned int i=0; i!=ndim; i++)
			PyTuple_SetItem(result,i,Py_BuildValue("i",centers[i]));
		return(result);
	}
}

static PyObject*
pysplinetable_evaluate(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"x", "centers", "derivatives", NULL};
	
	PyObject* pyx=NULL;
	PyObject* pycenters=NULL;
	unsigned long derivatives=0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|k", kwlist, &pyx, &pycenters, &derivatives))
		return(NULL);
	
	if(!PySequence_Check(pyx)){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pycenters)){
		PyErr_SetString(PyExc_ValueError, "centers must be a sequence");
		return(NULL);
	}
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	Py_ssize_t centerslen=PySequence_Length(pycenters);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	if(centerslen==-1){
		PyErr_SetString(PyExc_ValueError, "centers must be a sequence");
		return(NULL);
	}
	uint32_t ndim=splinetable_ndim(&self->table);
	if(xlen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of x must match the table dimension");
		return(NULL);
	}
	if(centerslen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of centers must match the table dimension");
		return(NULL);
	}
	
	unsigned int invalid_derivative_mask=~((1<<ndim)-1);
	if((derivatives&invalid_derivative_mask)!=0){
		PyErr_SetString(PyExc_ValueError, "Bits beyond the table dimension must not be set in derivatives");
		return(NULL);
	}
	
	//unpack x and centers
	//assume these are arbitrary sequences, not numpy arrays
	//TODO: should be possible to make a more efficient case for numpy arrays
	{
		//a small amount of evil
		double x[ndim];
		int centers[ndim];
		
		for(unsigned int i=0; i!=ndim; i++){
			PyObject* xi=PySequence_GetItem(pyx,i);
			x[i]=PyFloat_AsDouble(xi);
			Py_DECREF(xi); //done with this
			PyObject* centeri=PySequence_GetItem(pycenters,i);
			centers[i]=PyInt_AsLong(centeri);
			Py_DECREF(centeri); //done with this
		}
		
		double result=ndsplineeval(&self->table,x,centers,derivatives);
		return(Py_BuildValue("d",result));
	}
}

//attempts to do search centers and ndsplineeval in one step
static PyObject*
pysplinetable_evaluate_simple(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"x", "derivatives", NULL};
	
	PyObject* pyx=NULL;
	unsigned long derivatives=0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|k", kwlist, &pyx, &derivatives))
		return(NULL);
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	uint32_t ndim=splinetable_ndim(&self->table);
	if(xlen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of x must match the table dimension");
		return(NULL);
	}
	
	unsigned int invalid_derivative_mask=~((1<<ndim)-1);
	if((derivatives&invalid_derivative_mask)!=0){
		PyErr_SetString(PyExc_ValueError, "Bits beyond the table dimension must not be set in derivatives");
		return(NULL);
	}
	
	//unpack x
	//assume these are arbitrary sequences, not numpy arrays
	//TODO: should be possible to make a more efficient case for numpy arrays
	{
		//a small amount of evil
		double x[ndim];
		int centers[ndim];
		
		for(unsigned int i=0; i!=ndim; i++){
			PyObject* xi=PySequence_GetItem(pyx,i);
			x[i]=PyFloat_AsDouble(xi);
			Py_DECREF(xi); //done with this
			//printf("x[%u]=%lf\n",i,x[i]);
		}
		
		int status=tablesearchcenters(&self->table,x,centers);
		if(!status){
			PyErr_SetString(PyExc_ValueError, "tablesearchcenters failed");
			return(NULL);
		}
		
		double result=ndsplineeval(&self->table,x,centers,derivatives);
		return(Py_BuildValue("d",result));
	}
}

static PyObject*
pysplinetable_evaluate_gradient(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"x", "centers", NULL};
	
	PyObject* pyx=NULL;
	PyObject* pycenters=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, &pyx, &pycenters))
		return(NULL);
	
	if(!PySequence_Check(pyx)){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pycenters)){
		PyErr_SetString(PyExc_ValueError, "centers must be a sequence");
		return(NULL);
	}
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	Py_ssize_t centerslen=PySequence_Length(pycenters);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	if(centerslen==-1){
		PyErr_SetString(PyExc_ValueError, "centers must be a sequence");
		return(NULL);
	}
	uint32_t ndim=splinetable_ndim(&self->table);
	if(xlen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of x must match the table dimension");
		return(NULL);
	}
	if(centerslen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of centers must match the table dimension");
		return(NULL);
	}
	
	//unpack x and centers
	//assume these are arbitrary sequences, not numpy arrays
	//TODO: should be possible to make a more efficient case for numpy arrays
	{
		//a small amount of evil
		double x[ndim];
		int centers[ndim];
		double evaluates[ndim+1];
		
		for(unsigned int i=0; i!=ndim; i++){
			PyObject* xi=PySequence_GetItem(pyx,i);
			x[i]=PyFloat_AsDouble(xi);
			Py_DECREF(xi); //done with this
			PyObject* centeri=PySequence_GetItem(pycenters,i);
			centers[i]=PyInt_AsLong(centeri);
			Py_DECREF(centeri); //done with this
		}
		
		ndsplineeval_gradient(&self->table,x,centers,evaluates);
		
		PyObject* result=PyTuple_New(ndim+1);
		for(unsigned int i=0; i!=ndim+1; i++)
			PyTuple_SetItem(result,i,Py_BuildValue("d",evaluates[i]));
		return(result);
	}
}

static PyObject*
pysplinetable_deriv2(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"x", "centers", "derivatives", NULL};
	
	PyObject* pyx=NULL;
	PyObject* pycenters=NULL;
	unsigned long derivatives=0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOk", kwlist, &pyx, &pycenters, &derivatives))
		return(NULL);
	
	if(!PySequence_Check(pyx)){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pycenters)){
		PyErr_SetString(PyExc_ValueError, "centers must be a sequence");
		return(NULL);
	}
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	Py_ssize_t centerslen=PySequence_Length(pycenters);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	if(centerslen==-1){
		PyErr_SetString(PyExc_ValueError, "centers must be a sequence");
		return(NULL);
	}
	uint32_t ndim=splinetable_ndim(&self->table);
	if(xlen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of x must match the table dimension");
		return(NULL);
	}
	if(centerslen!=ndim){
		PyErr_SetString(PyExc_ValueError, "Length of centers must match the table dimension");
		return(NULL);
	}
	
	unsigned int invalid_derivative_mask=~((1<<ndim)-1);
	if((derivatives&invalid_derivative_mask)!=0){
		PyErr_SetString(PyExc_ValueError, "Bits beyond the table dimension must not be set in derivatives");
		return(NULL);
	}
	
	//unpack x and centers
	//assume these are arbitrary sequences, not numpy arrays
	//TODO: should be possible to make a more efficient case for numpy arrays
	{
		//a small amount of evil
		double x[ndim];
		int centers[ndim];
		
		for(unsigned int i=0; i!=ndim; i++){
			PyObject* xi=PySequence_GetItem(pyx,i);
			x[i]=PyFloat_AsDouble(xi);
			Py_DECREF(xi); //done with this
			PyObject* centeri=PySequence_GetItem(pycenters,i);
			centers[i]=PyInt_AsLong(centeri);
			Py_DECREF(centeri); //done with this
		}
		
		double result=ndsplineeval_deriv2(&self->table,x,centers,derivatives);
		return(Py_BuildValue("d",result));
	}
}

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
//TODO: fitting!
#endif //PHOTOSPLINE_INCLUDES_SPGLAM

//TODO: sampling? //requires C impl

//TODO: permutation

static PyGetSetDef pysplinetable_properties[] = {
	{"order", (getter)pysplinetable_getorder, NULL, "Order of spline in each dimension", NULL},
#ifdef HAVE_NUMPY
	{"knots", (getter)pysplinetable_getknots, NULL, "Knot vectors for each dimension", NULL},
	{"coefficients", (getter)pysplinetable_getcoeffcients, NULL, "Spline coefficients", NULL},
#endif
	{"extents", (getter)pysplinetable_getextents, NULL, "Range of support in each dimension", NULL},
	{"ndim", (getter)pysplinetable_getndim, NULL, "Number of dimensions", NULL},
	{NULL}
};

static PyMethodDef pysplinetable_methods[] = {
	{"write", (PyCFunction)pysplinetable_write, METH_KEYWORDS,
	 "Write the spline to a FITS file at the given path"},
	{"aux_value", (PyCFunction)pysplinetable_get_aux_value, METH_KEYWORDS,
	 "Get the value associated with an auxilliary key"},
	{"search_centers", (PyCFunction)pysplinetable_searchcenters, METH_KEYWORDS,
	 "Look up the basis function indices corresponding to a set of coordinates"},
	{"evaluate", (PyCFunction)pysplinetable_evaluate, METH_KEYWORDS,
	 "Evaluate the spline at a set of coordinates or its derivatives in the given dimensions"},
	{"evaluate_simple", (PyCFunction)pysplinetable_evaluate_simple, METH_KEYWORDS,
	 "Evaluate the spline at a set of coordinates or its derivatives in the given dimensions"},
	{"evaluate_gradient", (PyCFunction)pysplinetable_evaluate_gradient, METH_KEYWORDS,
	 "Evaluate the spline and all of its derivatives at a set of coordinates"},
	{"deriv2", (PyCFunction)pysplinetable_deriv2, METH_KEYWORDS,
	 "Evaluate the second derivative of the spline in the given dimensions"},
	
    {NULL}  /* Sentinel */
};

static PyTypeObject pysplinetableType = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"pyphotospline.Splinetable", /*tp_name*/
	sizeof(pysplinetable),     /*tp_basicsize*/
	0,                         /*tp_itemsize*/
	(destructor)pysplinetable_dealloc, /*tp_dealloc*/
	(printfunc)pysplinetable_print, /*tp_print*/
	0,                         /*tp_getattr*/
	0,                         /*tp_setattr*/
	0,                         /*tp_compare*/
	0,                         /*tp_repr*/
	0,                         /*tp_as_number*/
	0,                         /*tp_as_sequence*/
	0,                         /*tp_as_mapping*/
	0,                         /*tp_hash */
	(ternaryfunc)pysplinetable_evaluate_simple, /*tp_call*/
	0,                         /*tp_str*/
	0,                         /*tp_getattro*/
	0,                         /*tp_setattro*/
	0,                         /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Tensor product B-spline", /* tp_doc */
	0,		                   /* tp_traverse */
    0,		                   /* tp_clear */
    0,		                   /* tp_richcompare */
    0,		                   /* tp_weaklistoffset */
    0,		                   /* tp_iter */
    0,		                   /* tp_iternext */
    pysplinetable_methods,     /* tp_methods */
    0,                         /* tp_members */
    pysplinetable_properties,  /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)pysplinetable_init, /* tp_init */
    0,                         /* tp_alloc */
    pysplinetable_new,         /* tp_new */
};

static PyMethodDef photospline_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initpyphotospline(void){
	PyObject* module;
	
	pysplinetableType.tp_new = PyType_GenericNew;
	if (PyType_Ready(&pysplinetableType) < 0)
		return;
	
	module = Py_InitModule3("pyphotospline", photospline_methods,
	                   "A package for fitting gridded data to tensor-product "
	                   "B-spline surfaces and evaluating those surfaces");
	
	Py_INCREF(&pysplinetableType);
	PyModule_AddObject(module, "SplineTable", (PyObject *)&pysplinetableType);
	
#ifdef HAVE_NUMPY
	import_array();
#endif
}
