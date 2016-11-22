#include <Python.h>
//#include <numpy/ndarrayobject.h>
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
	printf("pysplinetable_new\n");
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
pysplinetable_ndim(pysplinetable* self){
	PyObject* result=Py_BuildValue("I",splinetable_ndim(&self->table));
	return(result);
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

//TODO: get_key
//OMIT: read_key (uses can do their own casting/parsing in python)
//TODO: write_key

static PyObject*
pysplinetable_order(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("I",splinetable_order(&self->table,dim));
	return(result);
}

static PyObject*
pysplinetable_nknots(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("I",splinetable_nknots(&self->table,dim));
	return(result);
}

//TODO: knots

static PyObject*
pysplinetable_knot(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", "idx", NULL};
	
	unsigned int dim,idx;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "II", kwlist, &dim, &idx))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	uint32_t nknots=splinetable_nknots(&self->table,dim);
	if(idx>=nknots){
		PyErr_SetString(PyExc_ValueError, "Knot index out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("d",splinetable_knot(&self->table,dim,idx));
	return(result);
}

static PyObject*
pysplinetable_lower_extent(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("d",splinetable_lower_extent(&self->table,dim));
	return(result);
}

static PyObject*
pysplinetable_upper_extent(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("d",splinetable_upper_extent(&self->table,dim));
	return(result);
}

static PyObject*
pysplinetable_period(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("d",splinetable_period(&self->table,dim));
	return(result);
}

static PyObject*
pysplinetable_ncoeffs(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("L",splinetable_ncoeffs(&self->table,dim));
	return(result);
}

static PyObject*
pysplinetable_total_ncoeffs(pysplinetable* self){
	PyObject* result=Py_BuildValue("L",splinetable_total_ncoeffs(&self->table));
	return(result);
}

static PyObject*
pysplinetable_stride(pysplinetable* self, PyObject* args, PyObject* kwds){
	static char* kwlist[] = {"dim", NULL};
	
	unsigned int dim;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "I", kwlist, &dim))
		return(NULL);
	
	uint32_t ndim=splinetable_ndim(&self->table);
	if(dim>=ndim){
		PyErr_SetString(PyExc_ValueError, "Dimension out of range");
		return(NULL);
	}
	
	PyObject* result=Py_BuildValue("L",splinetable_stride(&self->table,dim));
	return(result);
}

//TODO: coefficients

//TODO: evaluation!

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
//TODO: fitting!
#endif //PHOTOSPLINE_INCLUDES_SPGLAM

//TODO: sampling? //requires C impl

//TODO: permutation

static PyMethodDef pysplinetable_methods[] = {
	{"write", (PyCFunction)pysplinetable_write, METH_KEYWORDS,
	 "Write the spline to a FITS file at the given path"},
	{"ndim", (PyCFunction)pysplinetable_ndim, METH_NOARGS,
	 "Return the number of dimensions the table has"},
	{"order", (PyCFunction)pysplinetable_order, METH_KEYWORDS,
	 "Return the order of the spline in the given dimension"},
	{"nknots", (PyCFunction)pysplinetable_nknots, METH_KEYWORDS,
	 "Return the number of knots the spline has in the given dimension"},
	{"knot", (PyCFunction)pysplinetable_nknots, METH_KEYWORDS,
	 "Return the number of position of the knot with the given index in the given dimension"},
	{"lower_extent", (PyCFunction)pysplinetable_lower_extent, METH_KEYWORDS,
	 "Return the minimum extent of the spline in the given dimension"},
	{"upper_extent", (PyCFunction)pysplinetable_upper_extent, METH_KEYWORDS,
	 "Return the maximum extent of the spline in the given dimension"},
	{"period", (PyCFunction)pysplinetable_period, METH_KEYWORDS,
	 "Return the period of the spline in the given dimension"},
	{"ncoeffs", (PyCFunction)pysplinetable_ncoeffs, METH_KEYWORDS,
	 "Return the number of coefficients the spline has in the given dimension"},
	{"total_ncoeffs", (PyCFunction)pysplinetable_total_ncoeffs, METH_NOARGS,
	 "Return the total number of coefficients the table has"},
	{"stride", (PyCFunction)pysplinetable_stride, METH_KEYWORDS,
	 "Return the stride in the spline coefficients in the given dimension"},
	
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
	0,                         /*tp_call*/
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
    0,                         /* tp_getset */
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
	PyModule_AddObject(module, "Splinetable", (PyObject *)&pysplinetableType);
}

