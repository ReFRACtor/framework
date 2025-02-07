// There is a bunch of cookie cutter stuff needed to initialize a
// module for python, particularly if you want this to work with both
// python 2 and 3.
//
// This header sets this stuff up. Before including this, you should
// #define the PYTHON_MODULE_NAME, which is the name of the library
// python will load.
//
// There is often module specific code. This initialization will call
// the function module_init, passing the module as an argument. You
// should define this function before including.
//
// So for example:
//
// #define PYTHON_MODULE_NAME _swig_wrap
// #include "python_lib_init.h"
//
// static void module_init(PyObject* module)
// {
//    INIT_MODULE(package, "_swig_std", INIT_FUNC(swig_std));
// }
//
// Initializing cython code is almost the same, so we include it in
// this header. Just define DO_CYTHON for the difference. Rather than
// putting all the modules into a parent package, we put them into a
// python list. 

#include <Python.h>
#include <iostream>

// The double macros such as XINIT_NAME_PYTHON3 and INIT_NAME_PYTHON3
// are needed to get PYTHON_MODULE_NAME to expand before the macro
// call. See
// https://gcc.gnu.org/onlinedocs/cpp/Argument-Prescan.html#Argument-Prescan
// for details on these somewhat arcane rules for the c preprocessor
#define CONCATE(s,t) s ## t
#define XINIT_NAME_PYTHON3(x) CONCATE(PyInit_,x)
#define INIT_NAME_PYTHON3 XINIT_NAME_PYTHON3(PYTHON_MODULE_NAME)
#define XINIT_NAME_PYTHON2(x) CONCATE(init,x)
#define INIT_NAME_PYTHON2 XINIT_NAME_PYTHON2(PYTHON_MODULE_NAME)
#define str(s) #s
#define XPYTHON_MODULE_NAME_STRING(s) str(s)
#define PYTHON_MODULE_NAME_STRING XPYTHON_MODULE_NAME_STRING(PYTHON_MODULE_NAME)

static void module_init(PyObject* module);

// Python 2 and 3 do strings differently, so we have a simple macro to 
// keep from having lots of ifdefs spread around. See 
// https://wiki.python.org/moin/PortingExtensionModulesToPy3k for details on
// this.
#if PY_MAJOR_VERSION > 2
#define Text_FromUTF8(str) PyUnicode_FromString(str)
#else
#define Text_FromUTF8(str) PyString_FromString(str)
#endif

// Python 2 and 3 have different name for their swig init functions
#if PY_MAJOR_VERSION > 2
#define INIT_FUNC(S) PyInit_ ## S
#define INIT_TYPE PyObject *
#define INIT_MODULE init_extension_module3
#else
#define INIT_FUNC(S) init ## S
#define INIT_TYPE void
#define INIT_MODULE init_extension_module2
#endif

extern "C" {
#if PY_MAJOR_VERSION > 2
  PyObject * INIT_NAME_PYTHON3(void);
#else
  void INIT_NAME_PYTHON2(void);
#endif
}

// Used throughout SWIG wrapper, define here because it is convenient.
std::string parse_python_exception() {
  PyObject *type = NULL, *value = NULL, *tb = NULL;
  std::string ret = "Python error that I can't parse";
  PyErr_Fetch(&type, &value, &tb);
  PyObject * temp_bytes = PyUnicode_AsEncodedString(value, "ASCII", 
						    "ignore");
  if(temp_bytes) {
    ret = PyBytes_AS_STRING(temp_bytes); // Borrowed pointer
    Py_DECREF(temp_bytes);
  }
  // Try to get a traceback if we can
  PyErr_NormalizeException(&type, &value, &tb);
  PyObject* mod = PyImport_ImportModule("traceback");
  PyObject* err_str_list = NULL;
  if(tb) {
    err_str_list = PyObject_CallMethodObjArgs(mod,
	      Text_FromUTF8("format_exception"), type, value, tb, NULL);
  }
  if(err_str_list) {
    PyObject* err_str = 
      PyObject_CallMethodObjArgs(Text_FromUTF8(""),
				 Text_FromUTF8("join"), 
				 err_str_list, NULL);
    if(err_str) {
        PyObject * temp_bytes = PyUnicode_AsEncodedString(err_str, "ASCII", 
	"strict");
	if(temp_bytes) {
	  ret = PyBytes_AS_STRING(temp_bytes); // Borrowed pointer
	  Py_DECREF(temp_bytes);
	}
    }
    Py_XDECREF(err_str);
  }
  Py_XDECREF(mod);
  Py_XDECREF(err_str_list);
  Py_XDECREF(type);
  Py_XDECREF(value);
  Py_XDECREF(tb);
  return ret;
}

#if PY_MAJOR_VERSION > 2
// Version for python 3
static void init_extension_module3(PyObject* package, const char *modulename,
				  PyObject * (*initfunction)(void)) {
  PyObject *module = initfunction();
  PyObject *module_dic = PyImport_GetModuleDict();
  PyDict_SetItem(module_dic, Text_FromUTF8(modulename), module);
#ifdef DO_CYTHON
  // For cython, we create list rather than add to a package
  PyList_Append(package, Text_FromUTF8(modulename));
#else  
  if(PyModule_AddObject(package, (char *)modulename, module)) {
    std::cerr << "Initialisation in PyImport_AddObject failed for module "
	      << modulename << "\n";
    return;
  }
  Py_INCREF(module);
#endif
}
#else 
// Version for python 2
static void init_extension_module2(PyObject* package, const char *modulename,
				  void (*initfunction)(void)) {
  PyObject *module = PyImport_AddModule((char *)modulename);
  if(!module) {
    std::cerr << "Initialisation in PyImport_AddModule failed for module "
	      << modulename << "\n";
    return;
  }
#ifdef DO_CYTHON
  // For cython, we create list rather than add to a package
  PyList_Append(package, Text_FromUTF8(modulename));
#else  
  if(PyModule_AddObject(package, (char *)modulename, module)) {
    std::cerr << "Initialisation in PyImport_AddObject failed for module "
	      << modulename << "\n";
    return;
  }
  Py_INCREF(module);
#endif
  initfunction();
}
#endif


// This next blob of code comes from 
// https://wiki.python.org/moin/PortingExtensionModulesToPy3k

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m, PyObject *args) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef wrap_methods[] = {
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {NULL, NULL, METH_NOARGS, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int wrap_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int wrap_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        PYTHON_MODULE_NAME_STRING,
        NULL,
        sizeof(struct module_state),
        wrap_methods,
        NULL,
        wrap_traverse,
        wrap_clear,
        NULL
};

#define INITERROR return NULL

PyObject *
INIT_NAME_PYTHON3(void)

#else
#define INITERROR return

void
  INIT_NAME_PYTHON2(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule(PYTHON_MODULE_NAME_STRING, wrap_methods);
#endif

    if (module == NULL) {
        std::cerr << "Initialization failed\n";
        INITERROR;
    }
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException(PYTHON_MODULE_NAME_STRING ".Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    module_init(module);
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}



