// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "swig_boost_array_inc.i"
%include "swig_array_inc.i"

namespace boost {
template<class T, int D> class array {
public:
#ifdef SWIGPYTHON
  %rename(_size) size();
  int size();
  %pythoncode %{
    @property
    def size(self):
       return self._size()

    def set(self, v):
       '''Set the value to the given value. The value should be a container
       that we can enumerate over to fill in the value of the array'''
       for i, t in enumerate(v):
          if(i >= self.size):
              raise RuntimeError("Value must be exactly %d in size" % self.size)
          self[i] = t
       if(i != self.size - 1):
           raise RuntimeError("Value must be exactly %d in size" % self.size)

    def __iter__(self):
       for i in range(self.size):
           yield self[i]

    def __reduce__(self):
       return _new_from_set, (self.__class__, list(self))
  %}
  %extend {
    T __getitem__(int i) {return (*$self)[i];}
    void __setitem__(int i, const T& V) {(*$self)[i] = V;}
    std::string print_to_string() const
      { 
	std::ostringstream os;
	BOOST_FOREACH(T i, *$self)
	  os << i << " ";
	return os.str();
      }
   }
#else
   int size();
   %extend {
     T read(int i) {return (*$self)[i];}
     void write(int i, const T& V) {(*$self)[i] = V;}
   }
#endif
};

#ifdef SWIGPYTHON

//--------------------------------------------------------------
// Swig doesn't have typemap templates, so we define a macro to
// do this for each type and length, and then call the macro
// below to set this up for a range of types and sizes.
// The PRECEDENCE is the order that we check for a match with
// overloaded functions. The actual value doesn't matter too much,
// just make sure it is different from any of the type check
// in swig, and that it is different for different array types
//--------------------------------------------------------------
  
%define %boost_array_template(NAME,TYPE, LEN, PRECEDENCE)
%template(NAME) boost::array<TYPE, LEN>;

%typemap(in) const boost::array<TYPE, LEN>& (boost::array<TYPE, LEN> a, PythonObject numpy) 
{
  numpy.obj = to_numpy<TYPE >($input);
  if(!numpy.obj)
    return NULL;
  a = to_boost_array<TYPE, LEN>(numpy);
  $1 = &a;
}

%typemap(in) boost::array<TYPE, LEN> (PythonObject numpy) 
{
  numpy.obj = to_numpy<TYPE >($input);
  if(!numpy.obj)
    return NULL;
  $1 = to_boost_array<TYPE, LEN>(numpy);
}

%typecheck(PRECEDENCE) boost::array<TYPE, LEN>, const boost::array<TYPE, LEN>& {
  PythonObject t(to_numpy<TYPE >($input));
  if(!t.obj || PyArray_NDIM((PyArrayObject*) t.obj) != 1)
    return 0;
  blitz::Array<T, 1> b = to_blitz_array<TYPE, 1>(t.obj);
  return (b.rows() == LEN ? 1 : 0);
}
  
%enddef

%define %python_attribute_boost_array(NAME, TYPE, LEN)
   %extend {
    blitz::Array<TYPE, 1> _ ## NAME() const {
      blitz::Array<TYPE, 1> res(LEN);
      for(int i = 0; i < LEN; ++i)
        res(i) = $self->NAME[i];
      return res;
    }
    void _ ## NAME(const blitz::Array<TYPE, 1>& V) {
      if(V.rows() != LEN)
	throw std::runtime_error("Array not expeced size");
      for(int i = 0; i < LEN; ++i)
        $self->NAME[i] = V(i);
    }
    }
%pythoncode {
@property
def NAME(self):
    return self._ ## NAME()

@NAME.setter
def NAME(self, value):
  self._ ## NAME(value)
}
%enddef

#endif
}
%boost_array_template(Array_double_20, double, 20, 1241)
%boost_array_template(Array_double_12, double, 12, 1242)
%boost_array_template(Array_double_14, double, 14, 1243)
%boost_array_template(Array_double_3, double, 3, 1244)
%boost_array_template(Array_bool_20, bool, 20, 1245)


