#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>

#include <photospline/glam.h>

/* Sparse GLAM Python interface */

static PyObject *pybox(PyObject *self, PyObject *args);
static PyObject *pyrho(PyObject *self, PyObject *args);
static PyObject *pyfit(PyObject *self, PyObject *args, PyObject *kw);
static PyObject *pygrideval(PyObject *self, PyObject *args);
static PyObject *pynnls(PyObject *self, PyObject *args, PyObject *kw);

