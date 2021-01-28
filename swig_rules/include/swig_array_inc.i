// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// Stuff to include everywhere for using blitz::array. We also have 
// swig_array.i which includes stuff compiled in one spot.

%{
// Don't want to use threads with ruby
//#undef _REENTRANT

// Forward declaration
template<typename P_type> class PythonMemoryBlockReference;
// This is a very evil kludge. We use one of the function names in
// MemoryBlockReference
// to sneak in a friend declaration so we can access the internal block_
// variable. This is because we don't want to edit the actual blitz header.
#define blockLength() _fake() {return 0;}         \
  friend class PythonMemoryBlockReference<T_type>; \
  sizeType blockLength()
#include <blitz/memblock.h>
#undef blockLength
  
#include <blitz/array.h>
#include <blitz/range.h>
#define PY_ARRAY_UNIQUE_SYMBOL geocal_ARRAY_API
#ifndef DO_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif
// See https://github.com/numpy/numpy/issues/3008 for explanation of
// this.
// We'll have to update this as the numpy API increases
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <exception>
#include <vector>

//--------------------------------------------------------------
// Helper class for python that holds an object and when deleted
// decrements the reference to it.
//--------------------------------------------------------------

class PythonObject {
public:
  PythonObject(PyObject* Obj = 0) : obj(Obj) {}
  ~PythonObject() { Py_XDECREF(obj); }
  PyObject* obj;
  operator PyObject*() {return obj;}
};

//--------------------------------------------------------------
/// Special memory block for blitz. This doesn't actually own
/// the data, rather it owns a reference to a PythonObject that
/// keeps a reference to the numpy PythonObject so it doesn't
/// go away on the python side while we are using it on the
/// C++ side. 
//--------------------------------------------------------------

template<typename P_type>
class PythonMemoryBlock : public blitz::MemoryBlock<P_type> {
public:
  typedef P_type T_type;
  PythonMemoryBlock(PyObject* numpy_obj)
    : blitz::MemoryBlock<P_type>(PyArray_NBYTES((PyArrayObject*) numpy_obj),
                 (T_type *) PyArray_DATA((PyArrayObject*) numpy_obj)),
      python_obj(numpy_obj)
  {
    Py_XINCREF(python_obj); 
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
    std::cout << "PythonMemoryBlock: have reference to numpy object data at " << blitz::MemoryBlock<P_type>::data() << "\n"
              << "   numpy python object " << python_obj << "\n"
              << "   numpy python object ref count (after incrementing) " << Py_REFCNT(python_obj) <<"\n";
#endif    
  }
  virtual ~PythonMemoryBlock()
  {
    // Don't actually want to free real data. Instead, we just
    // have the python_obj go out of scope to clean up our reference
    // to it.
    blitz::MemoryBlock<P_type>::dataBlockAddress() = 0;
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
    std::cout << "PythonMemoryBlock: removing reference to numpy object data at " << blitz::MemoryBlock<P_type>::data() << "\n"
              << "   numpy python object " << python_obj << "\n"
              << "   numpy python object ref count (before decrementing) " << Py_REFCNT(python_obj) <<"\n";
#endif    
    Py_XDECREF(python_obj);
  }
private:
  PyObject* python_obj;
};

//--------------------------------------------------------------
/// Don't actually need a special MemoryBlockReference, except
/// that this is the only way to change the block of an existing
/// MemoryBlockReference (e.g., a blitz::Array). So this briefly
/// exists, uses itself to change the block, and then disappears.
//--------------------------------------------------------------

template<typename P_type>
class PythonMemoryBlockReference : public blitz::MemoryBlockReference<P_type> {
public:
  typedef P_type T_type;
  template<int N_rank>
  PythonMemoryBlockReference(blitz::Array<T_type, N_rank>& a,
                           PyObject* numpy_obj)
    : blitz::MemoryBlockReference<T_type>(0, a.data(),
                                          blitz::neverDeleteData)
  {
    blitz::MemoryBlockReference<T_type>::block_ =
      new PythonMemoryBlock<T_type>(numpy_obj);
    a.changeBlock(*this);
  }
};
 
PyObject* numpy_module();
PyObject* numpy_dot_float64();
PyObject* numpy_dot_float32();
PyObject* numpy_dot_int32();
PyObject* numpy_dot_uint32();
PyObject* numpy_dot_int16();
PyObject* numpy_dot_uint16();
PyObject* numpy_dot_int8();
PyObject* numpy_dot_uint8();
PyObject* numpy_dot_bool();

//--------------------------------------------------------------
// Helper routines to map a template type to the code numpy uses
// for that type.
//--------------------------------------------------------------

template<class T> int type_to_npy();
template<> inline int type_to_npy<double>() {return NPY_DOUBLE;}
template<> inline int type_to_npy<float>() {return NPY_FLOAT;}
template<> inline int type_to_npy<int>() {return NPY_INT;}
template<> inline int type_to_npy<unsigned int>() {return NPY_UINT;}
template<> inline int type_to_npy<short int>() {return NPY_SHORT;}
template<> inline int type_to_npy<unsigned short int>() {return NPY_USHORT;}
template<> inline int type_to_npy<char>() {return NPY_BYTE;}
template<> inline int type_to_npy<unsigned char>() {return NPY_UBYTE;}
template<> inline int type_to_npy<bool>() {return NPY_BOOL;}

inline char* npy_type_str(int type_enum)
{
    switch (type_enum) {
        case NPY_DOUBLE:
            return "double";
        case NPY_FLOAT:
            return "float";
        case NPY_INT:
            return "int";
        case NPY_UINT:
            return "unsigned int";
        case NPY_SHORT:
            return "short int";
        case NPY_USHORT:
            return "unsigned short int";
        case NPY_BYTE:
            return "char";
        case NPY_UBYTE:
            return "unsigned char";
        case NPY_BOOL:
            return "bool";
        default:
            return "Unknown";
    }
}

//--------------------------------------------------------------
// Custom exception for throwing errors we encounter when
// converting values that can be caught by typemaps and used
// to output the actual error message
//--------------------------------------------------------------

struct ArrayConversionException : public std::exception
{
  ArrayConversionException(const std::string& msg) : message(msg) {}
  virtual ~ArrayConversionException() throw() {}
  const char * what () const throw ()
  {
    return message.c_str();
  }

  std::string message;
};

//--------------------------------------------------------------
// Use the numpy command "asarray" to convert various python 
// objects to a numpy object. This may return null, if the 
// "asarray" fails. 
//--------------------------------------------------------------

template<class T> PyObject* to_numpy(PyObject* obj);

template<> inline PyObject* to_numpy<double>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                             PyString_FromString("asarray"), 
                                             obj, numpy_dot_float64(), NULL);
  // Don't worry about errors , since we just return a null
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<float>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                             PyString_FromString("asarray"), 
                                             obj, numpy_dot_float32(), NULL);
  // Don't worry about errors , since we just return a null
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<bool>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_bool(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<int>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_int32(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<unsigned int>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_uint32(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<short int>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_int16(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<unsigned short int>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_uint16(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<char>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_int8(), NULL);
  PyErr_Clear();
  return res;
}

template<> inline PyObject* to_numpy<unsigned char>(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(numpy_module(), 
                                    PyString_FromString("asarray"), 
                                    obj, numpy_dot_uint8(), NULL);
  PyErr_Clear();
  return res;
}
 
//--------------------------------------------------------------
// Convert a numpy array to a blitz::Array. The numpy should 
// already be the right data type before calling these (you can
// call to_numpy, if that is convenient). The underlying data is 
// still owned by the numpy object, so you need to make sure that
// the numpy object doesn't get deleted until you are done with
// the blitz::Array.
//
// If this fails, we throw an exception.
//--------------------------------------------------------------

template<class T, int D> 
inline blitz::Array<T, D> to_blitz_array(PyObject* numpy_obj)
{
  PyArrayObject* numpy = (PyArrayObject*) numpy_obj;
  if(PyArray_NDIM(numpy) != D) {
    std::stringstream err;
    err << "Error converting PyObject to Blitz array. "
        << "Dimension of array passed: " << PyArray_NDIM(numpy) 
        << " is not the expected size: " << D
        << " for blitz::Array<" << npy_type_str(type_to_npy<T>()) << ", " << D << ">";
    throw std::runtime_error(err.str().c_str());
  }
  if(PyArray_TYPE(numpy) != type_to_npy<T>()) {
    std::stringstream err;
    err << "Error converting PyObject to Blitz array. "
        << "Type of array passed: " << npy_type_str(PyArray_TYPE(numpy))
        << " is not the expected type: " << npy_type_str(type_to_npy<T>())
        << " for blitz::Array<" << npy_type_str(type_to_npy<T>()) << ", " << D << ">";
    throw std::runtime_error(err.str().c_str());
  }
  blitz::TinyVector<int, D> shape, stride;
  for(int i = 0; i < D; ++i) {
    shape(i) = PyArray_DIM(numpy, i);
    // Note numpy stride is in terms of bytes, while blitz in in terms
    // of type T.
    stride(i) = PyArray_STRIDE(numpy, i) / sizeof(T);
    if((int) (stride(i) * sizeof(T)) != (int) PyArray_STRIDE(numpy, i)) {
      throw std::runtime_error("Error converting PyObject to Blitz array. blitz::Array can't handle strides that aren't an even multiple of sizeof(T)");
    }
  }
  blitz::Array<T, D> a((T*)PyArray_DATA(numpy), shape, stride, 
                       blitz::neverDeleteData);
  // Stash a reference to numpy_obj in array, so it doesn't disappear
  // while the blitz::Array still exists
  PythonMemoryBlockReference<T> br(a, numpy_obj);
  return a;
}

template<class TYPE, int DIM> 
inline void iter_to_vector_of_arrays(PyObject *sequence, std::vector<blitz::Array<TYPE, DIM> >& arr_vec)
{

  if (PyUnicode_Check(sequence)) {
    throw ArrayConversionException("iter_to_vector_of_arrays: unicode objects are not supported");
  }

  PyObject *iterator = PyObject_GetIter(sequence);
  PyObject *item;

  if (iterator == NULL) {
    throw ArrayConversionException("iter_to_vector_of_arrays: passed object is not iterable");
  }

  int iter_index = 0;
  while ((item = PyIter_Next(iterator))) {
    PythonObject numpy;
    numpy.obj = to_numpy<TYPE>(item);
    if(!numpy.obj) {
      std::stringstream err_msg;
      err_msg << "iter_to_vector_of_arrays: object is not a numpy object at index: " << iter_index;
      throw ArrayConversionException(err_msg.str());
    }
    if(PyArray_NDIM((PyArrayObject*) numpy.obj) != DIM) {
      std::stringstream err_msg;
      err_msg << "iter_to_vector_of_arrays: incorrect dimension of numpy object at index: " << iter_index << ", expected dimension: " << DIM;
      throw ArrayConversionException(err_msg.str());
    }

    blitz::Array<TYPE, DIM> item_arr(to_blitz_array<TYPE, DIM>(numpy));
    arr_vec.push_back(item_arr);

    Py_DECREF(item);
    iter_index++;
  }

  Py_DECREF(iterator);

  if (PyErr_Occurred()) {
     throw ArrayConversionException("iter_to_vector_of_arrays: an error occured iterating over input object");
  }
}

%}

