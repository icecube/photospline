#include <Python.h>
#ifdef HAVE_NUMPY
	#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
	#include <numpy/ndarrayobject.h>
#endif
#include <stdio.h>
#include <limits.h>

#include <photospline/splinetable-mod.h>

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
typedef struct{
	PyObject_HEAD
	photospline::ndsparse* data;
} pyndsparse;

static void
pyndsparse_dealloc(pyndsparse* self){
	delete(self->data);
}

static PyObject*
pyndsparse_new(PyTypeObject* type, PyObject* args, PyObject* kwds){
	pyndsparse* self;
	
	self = (pyndsparse*)type->tp_alloc(type, 0);
	if(self){
		try{
			self->data=new photospline::ndsparse();
		}catch(std::exception& ex){
			PyErr_SetString(PyExc_Exception,
			                (std::string("Unable to allocate ndsparse: ")+ex.what()).c_str());
			return(NULL);
		}catch(...){
			PyErr_SetString(PyExc_Exception, "Unable to allocate ndsparse");
			return(NULL);
		}
	}
	
	return (PyObject*)self;
}

static int
pyndsparse_init(pyndsparse* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"rows", "ndim", NULL};
	unsigned long long rows, ndim;
	
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "KK", (char**)kwlist, &rows, &ndim))
        return -1;
	
	try{
		self->data=new photospline::ndsparse(rows,ndim);
	}catch(std::exception& ex){
		PyErr_SetString(PyExc_Exception,
						(std::string("Unable to allocate ndsparse: ")+ex.what()).c_str());
		return(-2);
	}catch(...){
		PyErr_SetString(PyExc_Exception, "Unable to allocate ndsparse");
		return(-2);
	}
	
    return 0;
}

static int
pyndsparse_print(pyndsparse* self, FILE* fp, int flags){
	fprintf(fp,"ndsparse with %zu dimension",self->data->ndim);
	if(self->data->ndim!=1)
		fprintf(fp,"s");
	fprintf(fp," and space for %zu entr",self->data->rows);
	if(self->data->ndim!=1)
		fprintf(fp,"ies");
	else
		fprintf(fp,"y");
	fprintf(fp," (%zu entr",self->data->entriesInserted);
	if(self->data->entriesInserted!=1)
		fprintf(fp,"ies filled)");
	else
		fprintf(fp,"y filled)");
	return(0);
}

static PyObject*
pyndsparse_insert(pyndsparse* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"value", "indices", NULL};
	
	double value;
	PyObject* pyindices=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "dO", (char**)kwlist,
									 &value, &pyindices))
		return(NULL);
	
	if(!PySequence_Check(pyindices)){
		PyErr_SetString(PyExc_ValueError, "indices must be a sequence");
		return(NULL);
	}
	Py_ssize_t len=PySequence_Length(pyindices);
	if(len!=self->data->ndim){
		PyErr_SetString(PyExc_ValueError, "Length of indices must match the ndsparse dimension");
		return(NULL);
	}
	
	unsigned int indices[self->data->ndim];
	for(unsigned int i=0; i!=self->data->ndim; i++){
		PyObject* idx=PySequence_GetItem(pyindices,i);
		indices[i]=PyInt_AsLong(idx);
		Py_DECREF(idx); //done with this
	}
	self->data->insertEntry(value,indices);
	
	Py_INCREF(Py_None);
	return Py_None;
}

static PyMethodDef pyndsparse_methods[] = {
	{"insert", (PyCFunction)pyndsparse_insert, METH_KEYWORDS,
		"Insert a data point"},
	{NULL}
};

static PyTypeObject pyndsparseType = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"pyphotospline.ndsparse", /*tp_name*/
	sizeof(pyndsparse),     /*tp_basicsize*/
	0,                         /*tp_itemsize*/
	(destructor)pyndsparse_dealloc, /*tp_dealloc*/
	(printfunc)pyndsparse_print, /*tp_print*/
	0,                         /*tp_getattr*/
	0,                         /*tp_setattr*/
	0,                         /*tp_compare*/
	0,                         /*tp_repr*/
	0,                         /*tp_as_number*/
	0,                         /*tp_as_sequence*/
	0,                         /*tp_as_mapping*/
	0,                         /*tp_hash */
	0,                         /*tp_call*/
	0,                         /*tp_str*/
	0,                         /*tp_getattro*/
	0,                         /*tp_setattro*/
	0,                         /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "", /* tp_doc */
	0,		                   /* tp_traverse */
    0,		                   /* tp_clear */
    0,		                   /* tp_richcompare */
    0,		                   /* tp_weaklistoffset */
    0,		                   /* tp_iter */
    0,		                   /* tp_iternext */
    pyndsparse_methods,        /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)pyndsparse_init, /* tp_init */
    0,                         /* tp_alloc */
    pyndsparse_new,         /* tp_new */
};
#endif //#ifdef PHOTOSPLINE_INCLUDES_SPGLAM

typedef struct{
	PyObject_HEAD
	photospline::splinetable<>* table;
} pysplinetable;

static void
pysplinetable_dealloc(pysplinetable* self){
	delete(self->table);
}

static PyObject*
pysplinetable_new(PyTypeObject* type, PyObject* args, PyObject* kwds){
	pysplinetable* self;

	self = (pysplinetable*)type->tp_alloc(type, 0);
	if(self){
		try{
			self->table=new photospline::splinetable<>();
		}catch(std::exception& ex){
			PyErr_SetString(PyExc_Exception,
			                (std::string("Unable to allocate spline table: ")+ex.what()).c_str());
			return(NULL);
		}catch(...){
			PyErr_SetString(PyExc_Exception, "Unable to allocate spline table");
			return(NULL);
		}
	}

	return (PyObject*)self;
}

static int
pysplinetable_init(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"path", NULL};
	char* path=NULL;
	
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", (char**)kwlist, &path))
        return -1;
	
	try{
		self->table=new photospline::splinetable<>(path);
	}catch(std::exception& ex){
		PyErr_SetString(PyExc_Exception,
						(std::string("Unable to allocate spline table: ")+ex.what()).c_str());
		return(-2);
	}catch(...){
		PyErr_SetString(PyExc_Exception, "Unable to allocate spline table");
		return(-2);
	}
	
    return 0;
}

static int
pysplinetable_print(pysplinetable* self, FILE* fp, int flags){
	uint32_t ndim=self->table->get_ndim();
	fprintf(fp,"Splinetable with %u dimension",ndim);
	if(ndim!=1)
		fprintf(fp,"s");
	return(0);
}

static PyObject*
pysplinetable_write(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"path", NULL};
	
	char* path=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", (char**)kwlist, &path))
		return(NULL);
	
	try{
		self->table->write_fits(path);
	}catch(std::exception& ex){
		PyErr_SetString(PyExc_Exception,ex.what());
		return(NULL);
	}
	
	return(Py_None);
}

static PyObject*
pysplinetable_get_aux_value(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"key", NULL};
	
	char* key=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", (char**)kwlist, &key))
		return(NULL);
	
	const char* value=self->table->get_aux_value(key);
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
	PyObject* list = PyTuple_New(self->table->get_ndim());
	for (uint32_t dim = 0; dim < self->table->get_ndim(); dim++) {
		npy_intp nknots = self->table->get_nknots(dim);
#ifndef NPY_1_7_API_VERSION //old numpy
		PyObject* knots = PyArray_New(&PyArray_Type, 1, &nknots, NPY_DOUBLE,
		            NULL, (void*)self->table->get_knots(dim), sizeof(double),
		            NPY_CARRAY_RO,NULL);
#else //newer numpy
		PyObject* knots = PyArray_SimpleNewFromData(1, &nknots, NPY_DOUBLE, (void*)self->table->get_knots(dim));
		PyArray_CLEARFLAGS((PyArrayObject*)knots, NPY_ARRAY_WRITEABLE);
#endif
		((PyArrayObject*)knots)->base=(PyObject*)self;
		Py_INCREF(self);
		PyTuple_SetItem(list, dim, knots);
	}
	
	return list;
}
#endif

static PyObject*
pysplinetable_getorder(pysplinetable* self, void *closure){
	PyObject *list = PyTuple_New(self->table->get_ndim());
	for (int dim = 0; dim < self->table->get_ndim(); dim++) {
		PyTuple_SetItem(list, dim, PyInt_FromLong(self->table->get_order(dim)));
	}
	
	return list;
}

static PyObject*
pysplinetable_getextents(pysplinetable* self, void *closure){
	PyObject *list = PyTuple_New(self->table->get_ndim());
	for (int dim = 0; dim < self->table->get_ndim(); dim++) {
		PyTuple_SetItem(list, dim, PyTuple_Pack(2,
		    PyFloat_FromDouble(self->table->lower_extent(dim)),
		    PyFloat_FromDouble(self->table->upper_extent(dim))));
	}
	
	return list;
}

#ifdef HAVE_NUMPY
static PyObject*
pysplinetable_getcoeffcients(pysplinetable* self, void *closure){
	unsigned int ndim = self->table->get_ndim();
	assert (ndim > 0);
	npy_intp dims[ndim];
	npy_intp strides[ndim];
	for (unsigned int i = 0; i < ndim; i++) {
		dims[i] = self->table->get_ncoeffs(i);
		strides[i] = sizeof(float)*self->table->get_stride(i);
	}
	PyObject* arr=PyArray_New(&PyArray_Type, ndim, dims, NPY_FLOAT, strides,
	    (void*)self->table->get_coefficients(), sizeof(float),
#ifndef NPY_1_7_API_VERSION //old numpy
	    NPY_CARRAY_RO,
#else //newer numpy
	    NPY_ARRAY_CARRAY_RO,
#endif
	    NULL);
	((PyArrayObject*)arr)->base=(PyObject*)self;
	Py_INCREF(self);
	return arr;
}
#endif

static PyObject*
pysplinetable_getndim(pysplinetable* self, void *closure){
	return PyInt_FromLong(self->table->get_ndim());
}

static PyObject*
pysplinetable_searchcenters(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"x", NULL};
	
	PyObject* pyx=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", (char**)kwlist, &pyx))
		return(NULL);
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	uint32_t ndim=self->table->get_ndim();
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
		
		if(!self->table->searchcenters(x,centers)){
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
	static const char* kwlist[] = {"x", "centers", "derivatives", NULL};
	
	PyObject* pyx=NULL;
	PyObject* pycenters=NULL;
	unsigned long derivatives=0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|k", (char**)kwlist, &pyx, &pycenters, &derivatives))
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
	uint32_t ndim=self->table->get_ndim();
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
		
		double result=self->table->ndsplineeval(x,centers,derivatives);
		return(Py_BuildValue("d",result));
	}
}

//attempts to do search centers and ndsplineeval in one step
static PyObject*
pysplinetable_evaluate_simple(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"x", "derivatives", NULL};
	
	PyObject* pyx=NULL;
	unsigned long derivatives=0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|k", (char**)kwlist, &pyx, &derivatives))
		return(NULL);
	
	Py_ssize_t xlen=PySequence_Length(pyx);
	
	if(xlen==-1){
		PyErr_SetString(PyExc_ValueError, "x must be a sequence");
		return(NULL);
	}
	uint32_t ndim=self->table->get_ndim();
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
		
		if(!self->table->searchcenters(x,centers)){
			PyErr_SetString(PyExc_ValueError, "tablesearchcenters failed");
			return(NULL);
		}
		
		double result=self->table->ndsplineeval(x,centers,derivatives);
		return(Py_BuildValue("d",result));
	}
}

static PyObject*
pysplinetable_evaluate_gradient(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"x", "centers", NULL};
	
	PyObject* pyx=NULL;
	PyObject* pycenters=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", (char**)kwlist, &pyx, &pycenters))
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
	uint32_t ndim=self->table->get_ndim();
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
		
		self->table->ndsplineeval_gradient(x,centers,evaluates);
		
		PyObject* result=PyTuple_New(ndim+1);
		for(unsigned int i=0; i!=ndim+1; i++)
			PyTuple_SetItem(result,i,Py_BuildValue("d",evaluates[i]));
		return(result);
	}
}

static PyObject*
pysplinetable_deriv2(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"x", "centers", "derivatives", NULL};
	
	PyObject* pyx=NULL;
	PyObject* pycenters=NULL;
	unsigned long derivatives=0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOk", (char**)kwlist, &pyx, &pycenters, &derivatives))
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
	uint32_t ndim=self->table->get_ndim();
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
		
		double result=self->table->ndsplineeval_deriv2(x,centers,derivatives);
		return(Py_BuildValue("d",result));
	}
}

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
static PyObject*
pyphotospline_glam_fit(PyObject* self, PyObject* args, PyObject* kwds);
#endif

//TODO: sampling?

static PyObject*
pysplinetable_permute(pysplinetable* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"permutation", NULL};
	
	PyObject* pypermutation=NULL;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", (char**)kwlist, &pypermutation))
		return(NULL);
	
	if(!PySequence_Check(pypermutation)){
		PyErr_SetString(PyExc_ValueError, "permutation must be a sequence");
		return(NULL);
	}
	
	try{
		std::vector<size_t> permutation;
		for(unsigned int i=0; i!=(unsigned)PySequence_Length(pypermutation); i++){
			PyObject* idx=PySequence_GetItem(pypermutation,i);
			permutation.push_back(PyInt_AsSsize_t(idx));
			Py_DECREF(idx); //done with this
		}
		self->table->permuteDimensions(permutation);
	}catch(std::exception& ex){
		PyErr_SetString(PyExc_ValueError, ex.what());
		return(NULL);
	}
	
	Py_INCREF(Py_None);
	return(Py_None);
}

static PyGetSetDef pysplinetable_properties[] = {
	{(char*)"order", (getter)pysplinetable_getorder, NULL, (char*)"Order of spline in each dimension", NULL},
#ifdef HAVE_NUMPY
	{(char*)"knots", (getter)pysplinetable_getknots, NULL, (char*)"Knot vectors for each dimension", NULL},
	{(char*)"coefficients", (getter)pysplinetable_getcoeffcients, NULL, (char*)"Spline coefficients", NULL},
#endif
	{(char*)"extents", (getter)pysplinetable_getextents, NULL, (char*)"Range of support in each dimension", NULL},
	{(char*)"ndim", (getter)pysplinetable_getndim, NULL, (char*)"Number of dimensions", NULL},
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
	{"permute_dimensions", (PyCFunction)pysplinetable_permute, METH_KEYWORDS,
	 "Permute the dimensions of an existing spline table"},
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
#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
	{"glam_fit", (PyCFunction)pyphotospline_glam_fit, METH_KEYWORDS,
	 "Fit a spline table to data"},
#endif
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
	
#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
	pyndsparseType.tp_new = PyType_GenericNew;
	if (PyType_Ready(&pyndsparseType) < 0)
		return;
#endif
	
	module = Py_InitModule3("pyphotospline", photospline_methods,
	                   "A package for fitting gridded data to tensor-product "
	                   "B-spline surfaces and evaluating those surfaces");
	
	Py_INCREF(&pysplinetableType);
	PyModule_AddObject(module, "SplineTable", (PyObject*)&pysplinetableType);
#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
	Py_INCREF(&pyndsparseType);
	PyModule_AddObject(module, "ndsparse", (PyObject*)&pyndsparseType);
#endif
	
#ifdef HAVE_NUMPY
	import_array();
#endif
}

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
static PyObject*
pyphotospline_glam_fit(PyObject* self, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"data", "weights", "coordinates", "order",
		"knots", "smoothing", "penaltyOrder", "monodim", "verbose", NULL};
	
	PyObject* pydata=NULL;
	PyObject* pyweights=NULL;
	PyObject* pycoordinates=NULL;
	PyObject* pyorder=NULL;
	PyObject* pyknots=NULL;
	PyObject* pysmoothing=NULL;
	PyObject* pyporder=NULL;
	unsigned long monodim=photospline::splinetable<>::no_monodim;
	unsigned char verbose=1;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOOO|kb", (char**)kwlist,
									 &pydata, &pyweights, &pycoordinates,
									 &pyorder, &pyknots, &pysmoothing, &pyporder,
									 &monodim, &verbose))
		return(NULL);
	
	//Check for valid inputs
	
	//data should be an ndsparse
	if(!PyObject_IsInstance(pydata,(PyObject*)&pyndsparseType)){
		PyErr_SetString(PyExc_ValueError, "data must be an ndsparse");
		return(NULL);
	}
	const ndsparse& data=*(((pyndsparse*)pydata)->data);
	
	//various things should be sequences
	if(!PySequence_Check(pyweights)){
		PyErr_SetString(PyExc_ValueError, "weights must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pycoordinates)){
		PyErr_SetString(PyExc_ValueError, "coordinates must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pyorder)){
		PyErr_SetString(PyExc_ValueError, "order must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pyknots)){
		PyErr_SetString(PyExc_ValueError, "knots must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pysmoothing)){
		PyErr_SetString(PyExc_ValueError, "smoothing must be a sequence");
		return(NULL);
	}
	if(!PySequence_Check(pyporder)){
		PyErr_SetString(PyExc_ValueError, "penaltyOrder must be a sequence");
		return(NULL);
	}
	
	//length of weights should match number of rows in data
	if(PySequence_Length(pyweights)!=data.rows){
		PyErr_SetString(PyExc_ValueError, "weights must be have the same length as data");
		return(NULL);
	}
	
	//length of coordinates should match dimension of data
	if(PySequence_Length(pycoordinates)!=data.ndim){
		PyErr_SetString(PyExc_ValueError, "coordinates must be have the same length as the dimension of data");
		return(NULL);
	}
	
	//length of order should match dimension of data
	if(PySequence_Length(pyorder)!=data.ndim){
		PyErr_SetString(PyExc_ValueError, "order must be have the same length as the dimension of data");
		return(NULL);
	}
	
	//length of knots should match dimension of data
	if(PySequence_Length(pyknots)!=data.ndim){
		PyErr_SetString(PyExc_ValueError, "knots must be have the same length as the dimension of data");
		return(NULL);
	}
	
	//coordinates should be a sequence of sequences
	for(unsigned int i=0; i!=(unsigned)PySequence_Length(pycoordinates); i++){
		PyObject* pycoordinates_i=PySequence_GetItem(pycoordinates,i);
		if(!PySequence_Check(pycoordinates_i)){
			PyErr_SetString(PyExc_ValueError, "all entries in coordinates must be sequences");
			return(NULL);
		}
		Py_DECREF(pycoordinates_i);
	}
	
	//knots should be a sequence of sequences
	for(unsigned int i=0; i!=(unsigned)PySequence_Length(pyknots); i++){
		PyObject* pyknots_i=PySequence_GetItem(pyknots,i);
		if(!PySequence_Check(pyknots_i)){
			PyErr_SetString(PyExc_ValueError, "all entries in knots must be sequences");
			return(NULL);
		}
		Py_DECREF(pyknots);
	}
	
	//If possible we do not want to copy the input data, but it may be necessary.
	//To avoid an explosion of type combinations, though, we will always pass
	//array_views for most of the arguments to fit, and in case we do need to
	//copy these vectors will be used for their backing stores.
	std::vector<double> weights_store(0);
	std::vector<std::vector<double>> coordinates_store(0);
	std::vector<uint32_t> order_store(0);
	std::vector<std::vector<double>> knots_store(0);
	std::vector<double> smoothing_store(0);
	std::vector<uint32_t> porder_store(0);
	//These views will always be the actual arguments
	using photospline::detail::array_view;
	array_view<double> weights;
	std::vector<array_view<double>> coordinates;
	array_view<uint32_t> order;
	std::vector<array_view<double>> knots;
	array_view<double> smoothing;
	array_view<uint32_t> porder;
	
	//We can avoid copying input data if it is in a numpy array which has the
	//correct type and does not have funny layout.
#define compatible_numpy_array(obj,element_type) (\
	PyArray_Check(obj) \
	&& PyArray_EquivTypenums(PyArray_DESCR(obj)->type_num,element_type) \
	&& (PyArray_FLAGS(obj)&NPY_C_CONTIGUOUS) \
	&& (PyArray_FLAGS(obj)&NPY_ALIGNED) \
	)
	
	//Extract inputs
	
	// extract weights
#ifdef HAVE_NUMPY
	if(compatible_numpy_array(pyweights,NPY_DOUBLE))
		weights.reset((double*)PyArray_DATA(pyweights),data.rows);
	else //note sneaky line break across #endif
#endif
	{ //have to copy
		weights_store.resize(data.rows);
		for(unsigned int i=0; i!=data.rows; i++){
			PyObject* weight_i=PySequence_GetItem(pyweights,i);
			weights_store[i]=PyFloat_AsDouble(weight_i);
			Py_DECREF(weight_i); //done with this
		}
		weights.reset(weights_store.data(),data.rows);
	}
	
	// extract coordinates
	coordinates.resize(data.ndim);
	for(unsigned int j=0; j!=data.ndim; j++){
		PyObject* pycoordinates_i=PySequence_GetItem(pycoordinates,j);
		Py_ssize_t ncoordinates=PySequence_Length(pycoordinates_i);
		assert(ncoordinates>0);
		
#ifdef HAVE_NUMPY
		if(compatible_numpy_array(pycoordinates_i,NPY_DOUBLE))
			coordinates[j].reset((double*)PyArray_DATA(pycoordinates_i),ncoordinates);
		else //note sneaky line break across #endif
#endif
		{ //have to copy
			coordinates_store.resize(j+1);
			coordinates_store[j].resize(ncoordinates);
			for(unsigned int i=0; i!=(size_t)ncoordinates; i++){
				PyObject* coordinate_i=PySequence_GetItem(pycoordinates_i,i);
				coordinates_store[j][i]=PyFloat_AsDouble(coordinate_i);
				Py_DECREF(coordinate_i); //done with this
			}
			coordinates[j].reset(coordinates_store[j].data(),ncoordinates);
		}
		
		Py_DECREF(pycoordinates_i); //done with this
	}
	
	// extract order
#ifdef HAVE_NUMPY
	//TODO: Numpy's UInt may not be the same as uint32_t
	if(compatible_numpy_array(pyorder,NPY_UINT))
		order.reset((uint32_t*)PyArray_DATA(pyorder),data.ndim);
	else //note sneaky line break across #endif
#endif
	{ //have to copy
		order_store.resize(data.ndim);
		for(unsigned int i=0; i!=data.ndim; i++){
			PyObject* pyorder_i=PySequence_GetItem(pyorder,i);
			order_store[i]=PyInt_AsUnsignedLongMask(pyorder_i);
			Py_DECREF(pyorder_i); //done with this
		}
		order.reset(order_store.data(),data.ndim);
	}
	
	// extract knots
	knots.resize(data.ndim);
	for(unsigned int j=0; j!=data.ndim; j++){
		PyObject* pyknots_i=PySequence_GetItem(pyknots,j);
		Py_ssize_t nknots=PySequence_Length(pyknots_i);
		assert(nknots>0);
		
#ifdef HAVE_NUMPY
		if(compatible_numpy_array(pyknots_i,NPY_DOUBLE))
			knots[j].reset((double*)PyArray_DATA(pyknots_i),nknots);
		else //note sneaky line break across #endif
#endif
		{ //have to copy
			knots_store.resize(j+1);
			knots_store[j].resize(nknots);
			for(unsigned int i=0; i!=(size_t)nknots; i++){
				PyObject* knot_i=PySequence_GetItem(pyknots_i,i);
				knots_store[j][i]=PyFloat_AsDouble(knot_i);
				Py_DECREF(knot_i); //done with this
			}
			knots[j].reset(knots_store[j].data(),nknots);
		}
		
		Py_DECREF(pyknots_i); //done with this
	}
	
	// extract smoothing
#ifdef HAVE_NUMPY
	if(compatible_numpy_array(pysmoothing,NPY_DOUBLE))
		smoothing.reset((double*)PyArray_DATA(pysmoothing),data.ndim);
	else //note sneaky line break across #endif
#endif
	{ //have to copy
		smoothing_store.resize(data.rows);
		for(unsigned int i=0; i!=data.ndim; i++){
			PyObject* smoothing_i=PySequence_GetItem(pysmoothing,i);
			smoothing_store[i]=PyFloat_AsDouble(smoothing_i);
			Py_DECREF(smoothing_i); //done with this
		}
		smoothing.reset(smoothing_store.data(),data.ndim);
	}
	
	// extract porder
#ifdef HAVE_NUMPY
	//TODO: Numpy's UInt may not be the same as uint32_t
	if(compatible_numpy_array(pyporder,NPY_UINT))
		porder.reset((uint32_t*)PyArray_DATA(pyporder),data.ndim);
	else //note sneaky line break across #endif
#endif
	{ //have to copy
		porder_store.resize(data.ndim);
		for(unsigned int i=0; i!=data.ndim; i++){
			PyObject* pyporder_i=PySequence_GetItem(pyporder,i);
			porder_store[i]=PyInt_AsUnsignedLongMask(pyporder_i);
			Py_DECREF(pyporder_i); //done with this
		}
		porder.reset(porder_store.data(),data.ndim);
	}
	
#undef compatible_numpy_array
	
	//make a spline
	std::unique_ptr<pysplinetable,void(*)(pysplinetable*)>
	spline((pysplinetable*)pysplinetable_new(&pysplinetableType,Py_None,Py_None),
	                         &pysplinetable_dealloc);
	
	//Finally do the fit
	try{
		spline->table->fit(data, weights, coordinates, order, knots, smoothing, porder, monodim, verbose);
	}catch(std::logic_error& ex){
		PyErr_SetString(PyExc_ValueError, ex.what());
		return(NULL);
	}catch(std::exception& ex){
		PyErr_SetString(PyExc_Exception, ex.what());
		return(NULL);
	}
	
	return((PyObject*)spline.release());
}
#endif //PHOTOSPLINE_INCLUDES_SPGLAM
