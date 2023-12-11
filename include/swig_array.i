// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
#ifndef SWIG_MODULE_ALREADY_DONE
#if SWIG_VERSION < 0x040000  
%module(directors="1", allprotected="1") foo
#else
%module(moduleimport="from ._swig_wrap import $module", directors="1", allprotected="1") foo
#endif
#endif

%{
#define DO_IMPORT_ARRAY
%}
%include "swig_array_inc.i"
%include <std_vector.i>

//--------------------------------------------------------------
// Before using numpy, we need to call the numpy supplied 
// function 'import_array' which imports the module and does set up.
//--------------------------------------------------------------

%init {
    import_array();
}

%{

//--------------------------------------------------------------
// Return numpy module
//--------------------------------------------------------------

PyObject* numpy_module()
{
  static PyObject* mod = 0;
  if(!mod)
    mod = PyImport_ImportModule("numpy");
  return mod;
}

//--------------------------------------------------------------
// Return numpy.float64
//--------------------------------------------------------------

PyObject* numpy_dot_float64()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "float64");
  return res;
}

PyObject* numpy_dot_float32()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "float32");
  return res;
}
 
PyObject* numpy_dot_int32()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "int32");
  return res;
}

PyObject* numpy_dot_uint32()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "uint32");
  return res;
}

PyObject* numpy_dot_int16()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "int16");
  return res;
}

PyObject* numpy_dot_uint16()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "uint16");
  return res;
}

PyObject* numpy_dot_int8()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "int8");
  return res;
}

PyObject* numpy_dot_uint8()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "uint8");
  return res;
}
 
PyObject* numpy_dot_bool()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "bool");
  return res;
}

PyObject* numpy_dot_object()
{
  static PyObject* res = 0;
  if(!res)
    res = PyObject_GetAttrString(numpy_module(), "object");
  return res;
}

%}

// Allow conversion to a binary String in the target language
%include "cdata.i"

namespace blitz {
// Define blitz::Array for use in Python
template<class T, int D> class Array  {
public:
  // These functions aren't normally used, because typemaps
  // automatically map from blitz::Array to numpy or narry. But leave
  // in place for helping with other languages.
  Array(int e1);
  Array(int e1, int e2, blitz::GeneralArrayStorage<D> storage = 
	blitz::FortranArray<D>());
  Array(int e1, int e2, int e3, blitz::GeneralArrayStorage<D> storage = 
	blitz::FortranArray<D>());
  Array(int e1, int e2, int e3, int e4, blitz::GeneralArrayStorage<D> storage = 
	blitz::FortranArray<D>());
  T* data();
  int size() const;
  %extend {
     T read(int i1) {return (*$self)(i1);}
     T read(int i1,int i2) {return (*$self)(i1,i2);}
     T read(int i1, int i2, int i3) {return (*$self)(i1,i2,i3);}
     T read(int i1, int i2, int i3, int i4) {return (*$self)(i1,i2,i3,i4);}
     void write(int i1, T val) {(*$self)(i1) = val;}
     void write(int i1,int i2, T val) {(*$self)(i1,i2) = val;}
     void write(int i1, int i2, int i3, T val) {(*$self)(i1,i2,i3) = val;}
     void write(int i1, int i2, int i3, int i4, T val) 
     {(*$self)(i1,i2,i3,i4) = val;}
     void* datav() { return (void*) $self->data(); }
     int shape0() { return $self->shape()[0]; }
     int shape1() { return $self->shape()[1]; }
     int shape2() { return $self->shape()[2]; }
     int shape3() { return $self->shape()[3]; }
  }
};

// Define blitz::Range for use in Python
enum { fromStart = Blitz::INT_MIN, toEnd = Blitz::INT_MIN };
class Range {
public:
  Range();
  explicit Range(int slicePosition);
  Range(int first, int last, int stride=1);
  int first(int lowRange = 0) const;
  int last(int highRange = 0) const;
  unsigned length(int =0) const;
  int stride() const;
  bool isAscendingContiguous() const;
  void setRange(int first, int last, int stride=1);
  static Range all();
  bool isUnitStride() const;
};

}

//--------------------------------------------------------------
// Convert to numpy. Note that there is a complication in the 
// lifetime of the pointed to array. numpy can't take ownership
// of the memory in the blitz::Array, since it wasn't allocated
// by python. Instead, numpy just points to the memory. To ensure
// that the blitz::Array memory isn't freeded, we also stash a
// python object wrapping around the blitz::Array that holds onto
// the object. This gets placed in a special area set up by numpy
// exactly for this purpose called "BASE". When the numpy array 
// get deleted, it also deletes the numpy. If this is the only
// reference to the blitz::Array memory, then the memory gets
// cleaned up then.
//--------------------------------------------------------------

%define %blitz_to_numpy(TYPE, DIM, from_obj, to_obj)
    // Copy out dimensions and stride from blitz array
  npy_intp dims[DIM], stride[DIM];
  for(int i = 0; i < DIM; ++i) {
    dims[i] = from_obj->extent(i);
    // Note numpy stride is in terms of bytes, while blitz in in terms
    // of type T.
    stride[i] = from_obj->stride(i) * sizeof(TYPE);
  }

  // Create new numpy object using Numpy C API
  to_obj = PyArray_New(&PyArray_Type, DIM, dims, type_to_npy<TYPE >(), 
			stride, from_obj->data(), 0, 0, 0);
  blitz::Array<TYPE, DIM>* t = new blitz::Array<TYPE, DIM>(*from_obj);
  // Stash pointer to original blitz array as detailed above
  PyArray_SetBaseObject((PyArrayObject*) to_obj, 
			SWIG_NewPointerObj(SWIG_as_voidptr(t), 
				   $descriptor(blitz::Array<TYPE, DIM>*), 					   SWIG_POINTER_NEW | SWIG_POINTER_OWN ));
%enddef

//************************************************************
// Type map to use python type numpy as input and output
//************************************************************

#ifdef SWIGPYTHON
//--------------------------------------------------------------
// Swig doesn't have typemap templates, so we define a macro to
// do this for each type and dimension, and then call the macro
// below to set this up for a range of types and sizes.
// The PRECEDENCE is the order that we check for a match with
// overloaded functions. The actual value doesn't matter too much,
// just make sure it is different from any of the type check
// in swig, and that it is different for different array types
//--------------------------------------------------------------

%define %array_template(NMTYPE, TYPE,DIM, PRECEDENCE)

//--------------------------------------------------------------
// Convert to numpy. See description above for lifetime issues.
//--------------------------------------------------------------

%typemap(out) blitz::Array<TYPE, DIM> {
  // Treat as pointer for the purposes of the macro
  %blitz_to_numpy(TYPE, DIM, (&$1), $result);
}

//--------------------------------------------------------------
// Convert to numpy. See description above for lifetime issues.
//--------------------------------------------------------------

%typemap(out) const blitz::Array<TYPE, DIM>& {
  %blitz_to_numpy(TYPE, DIM, $1, $result);
}

%typemap(out) blitz::Array<TYPE, DIM>& {
  npy_intp dims[DIM], stride[DIM];
  for(int i = 0; i < DIM; ++i) {
    dims[i] = $1->extent(i);
    // Note numpy stride is in terms of bytes, while blitz in in terms
    // of type T.
    stride[i] = $1->stride(i) * sizeof(TYPE);
  }
  $result = PyArray_New(&PyArray_Type, DIM, dims, type_to_npy<TYPE>(), 
			stride, $1->data(), 0, 0, 0);
  PyArray_UpdateFlags((PyArrayObject*)$result, NPY_ARRAY_WRITEABLE);
  blitz::Array<TYPE, DIM>* t = new blitz::Array<TYPE, DIM>(*$1);
  PyArray_SetBaseObject((PyArrayObject*)$result, 
			SWIG_NewPointerObj(SWIG_as_voidptr(t), 
				   $descriptor(blitz::Array<TYPE, DIM>*), 					   SWIG_POINTER_NEW | SWIG_POINTER_OWN ));
}

%typemap(out) std::vector<blitz::Array<TYPE, DIM> > {
    
    $result = PyList_New($1.size());
    for(int idx = 0; idx < (int) $1.size(); idx++) {
        PyObject* element_np_arr;
        %blitz_to_numpy(TYPE, DIM, (&$1.at(idx)), element_np_arr);
        PyList_SetItem($result, idx, element_np_arr);
    }
}

//--------------------------------------------------------------
/// Handle multiple array returns
//--------------------------------------------------------------

%typemap(in, numinputs=0) blitz::Array<TYPE, DIM>& OUTPUT (blitz::Array<TYPE, DIM> temp) {
   $1 = &temp;
}

%typemap(in, numinputs=0) blitz::Array<TYPE, DIM>& OUTPUT1 (blitz::Array<TYPE, DIM> temp) {
   $1 = &temp;
}

%typemap(in, numinputs=0) blitz::Array<TYPE, DIM>& OUTPUT2 (blitz::Array<TYPE, DIM> temp) {
   $1 = &temp;
}


%typemap(argout) blitz::Array<TYPE, DIM>& OUTPUT {
   PyObject *res;
   %blitz_to_numpy(TYPE, DIM, $1, res);
   $result = SWIG_AppendOutput($result, res);
}

%typemap(argout) blitz::Array<TYPE, DIM>& OUTPUT1 {
  PyObject *res;
  %blitz_to_numpy(TYPE, DIM, $1, res);
  $result = SWIG_AppendOutput($result, res);
}

%typemap(argout) blitz::Array<TYPE, DIM>& OUTPUT2 {
  PyObject *res;
  %blitz_to_numpy(TYPE, DIM, $1, res);
  $result = SWIG_AppendOutput($result, res);
}

//--------------------------------------------------------------
// Convert any type first to a numpy array (doesn't copy if 
// already a numpy array), and then set blitz array to point to
// this.
//--------------------------------------------------------------

%typemap(in) const blitz::Array<TYPE, DIM>& (blitz::Array<TYPE, DIM> a, PythonObject numpy) 
{
  int res = SWIG_ConvertPtr($input, (void**)(&$1), $descriptor(blitz::Array<TYPE, DIM>*), 
			    %convertptr_flags);
  if(!SWIG_IsOK(res)) {
    numpy.obj = to_numpy<TYPE >($input);
    if(!numpy.obj) {
      SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)numpy.obj) !=DIM) {
      SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
      return NULL;
    }
    a.reference(to_blitz_array<TYPE, DIM>(numpy));
    $1 = &a;
  }
}

//--------------------------------------------------------------
// Version that forces a copy of data
//--------------------------------------------------------------

%typemap(in) const blitz::Array<TYPE, DIM>& FORCE_COPY (blitz::Array<TYPE, DIM> a, PythonObject numpy) 
{
  int res = SWIG_ConvertPtr($input, (void**)(&$1), $descriptor(blitz::Array<TYPE, DIM>*), 
			    %convertptr_flags);
  if(!SWIG_IsOK(res)) {
    numpy.obj = to_numpy<TYPE >($input);
    if(!numpy.obj) {
      SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
    }
    if(PyArray_NDIM((PyArrayObject*)numpy.obj) !=DIM) {
      SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
    }
    a.reference(to_blitz_array<TYPE, DIM>(numpy).copy());
    $1 = &a;
  }
}

//--------------------------------------------------------------
// Convert any type first to a numpy array (doesn't copy if 
// already a numpy array), and then set blitz array to point to
// this.
//--------------------------------------------------------------

%typemap(in) blitz::Array<TYPE, DIM> (PythonObject numpy) 
{
  numpy.obj = to_numpy<TYPE >($input);
  if(!numpy.obj) {
    SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
  }
  if(PyArray_NDIM((PyArrayObject*)numpy.obj) !=DIM) {
    SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
  }
  $1 = to_blitz_array<TYPE, DIM>(numpy);
}

//--------------------------------------------------------------
// Handle conversion in directors
//--------------------------------------------------------------

%typemap(directorout) blitz::Array<TYPE, DIM> (PythonObject numpy) 
{
  PythonObject t(to_numpy<TYPE >($input));
  if(!t.obj) {
    SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
  }
  if(PyArray_NDIM((PyArrayObject*)t.obj) !=DIM) {
    SWIG_Error(SWIG_TypeError, "in method '$symname', expecting type  Array<TYPE,DIM>");
  }
  $result.reference(to_blitz_array<TYPE, DIM>(t).copy());
}

%typemap(directorin) const blitz::Array<TYPE, DIM>& 
{
  npy_intp dims[DIM], stride[DIM];
  for(int i = 0; i < DIM; ++i) {
    dims[i] = $1.extent(i);
    // Note numpy stride is in terms of bytes, while blitz in in terms
    // of type T.
    stride[i] = $1.stride(i) * sizeof(TYPE);
  }
  PyObject* res = PyArray_New(&PyArray_Type, DIM, dims, type_to_npy<TYPE >(), 
			      stride, const_cast<TYPE*>($1.data()), 0, 0, 0);
  blitz::Array<TYPE, DIM>* t = new blitz::Array<TYPE, DIM>($1);
  PyArray_SetBaseObject((PyArrayObject*)res, 
			SWIG_NewPointerObj(SWIG_as_voidptr(t), 
				   $descriptor(blitz::Array<TYPE, DIM>*), 					   SWIG_POINTER_NEW | SWIG_POINTER_OWN ));
  $input = res;
}

//--------------------------------------------------------------
// Check if object can be converted to a blitz::Array.
//--------------------------------------------------------------

%typecheck(PRECEDENCE) blitz::Array<TYPE, DIM>, const blitz::Array<TYPE, DIM>& {
  PythonObject t(to_numpy<TYPE >($input));
  $1 = (t.obj && PyArray_NDIM((PyArrayObject*)t.obj) ==DIM ? 1 : 0);
}

//--------------------------------------------------------------
// Convert a list of numpy arrays into a vector of blitz arrays
//--------------------------------------------------------------

%typecheck(PRECEDENCE) std::vector<blitz::Array<TYPE, DIM> >, const std::vector<blitz::Array<TYPE, DIM> > {
  // Here we are just checking that the item is an iterable, we are not checking that
  // each element can be converted to a blitz array. Leave that as a failure during the typemap conversion
  PyObject *iterator = PyObject_GetIter($input);
  $1 = iterator != NULL;
}

%typecheck(PRECEDENCE) std::vector<blitz::Array<TYPE, DIM> >&, const std::vector<blitz::Array<TYPE, DIM> >& {
  // Here we are just checking that the item is an iterable, we are not checking that
  // each element can be converted to a blitz array. Leave that as a failure during the typemap conversion
  PyObject *iterator = PyObject_GetIter($input);
  $1 = iterator != NULL;
}

%typemap(in) std::vector<blitz::Array<TYPE, DIM> > (std::vector<blitz::Array<TYPE, DIM> > arr_vec), const std::vector<blitz::Array<TYPE, DIM> > (std::vector<blitz::Array<TYPE, DIM> > arr_vec)
{ 
  try {
    iter_to_vector_of_arrays<TYPE, DIM>($input, arr_vec);
  } catch (ArrayConversionException& exc) {
    SWIG_exception_fail(SWIG_TypeError, exc.what());
  }

  $1 = arr_vec;
}

%typemap(in) std::vector<blitz::Array<TYPE, DIM> >& (std::vector<blitz::Array<TYPE, DIM> > arr_vec), const std::vector<blitz::Array<TYPE, DIM> >& (std::vector<blitz::Array<TYPE, DIM> > arr_vec)
{ 
  try {
    iter_to_vector_of_arrays<TYPE, DIM>($input, arr_vec);
  } catch (ArrayConversionException& exc) {
    SWIG_exception_fail(SWIG_TypeError, exc.what());
  }

  $1 = &arr_vec;
}

//--------------------------------------------------------------
// Create a name usable by Python for each evaluated blitz array type
//--------------------------------------------------------------

%template(BlitzArray_ ## NMTYPE ## _ ## DIM) blitz::Array<TYPE, DIM>;

%enddef

//--------------------------------------------------------------
// Evaluate the array template for multiple dimension sizes
//--------------------------------------------------------------

%define %array_all_template(NMTYPE, TYPE, PRECEDENCE)

%array_template(NMTYPE, TYPE, 1, PRECEDENCE ## 0);
%array_template(NMTYPE, TYPE, 2, PRECEDENCE ## 1);
%array_template(NMTYPE, TYPE, 3, PRECEDENCE ## 2);
%array_template(NMTYPE, TYPE, 4, PRECEDENCE ## 3);
%array_template(NMTYPE, TYPE, 5, PRECEDENCE ## 4);
%array_template(NMTYPE, TYPE, 6, PRECEDENCE ## 5);
%array_template(NMTYPE, TYPE, 7, PRECEDENCE ## 6);
%array_template(NMTYPE, TYPE, 8, PRECEDENCE ## 7);

%enddef

//--------------------------------------------------------------
// Evaluate array templates for various types, for multiple dimension sizes
//--------------------------------------------------------------

%array_all_template(double, double, 115);
%array_all_template(bool, bool, 116);
%array_all_template(int, int, 117);
%array_all_template(uint, unsigned int, 118);
%array_all_template(float, float, 119);
%array_all_template(short, short int, 120);
%array_all_template(ushort, unsigned short int, 121);
%array_all_template(char, char, 122);
%array_all_template(uchar, unsigned char, 123);
#endif  // end SWIGPYTHON

//--------------------------------------------------------------
// List of things "import *" will include
//--------------------------------------------------------------

%pythoncode %{
__all__ = []  
for t in ("double", "bool", "int", "uint", "float", "short", "ushort", "char", "uchar"):
    for d in range(1,8+1):
        __all__.append("BlitzArray_%s_%d" % (t, d))
%}
