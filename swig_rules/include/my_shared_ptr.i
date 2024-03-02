// There appears to be a bug in the shared_ptr handler of SWIG, as of
// version 3.0.12. This is "boost_shared_ptr.i" from
// share/swig/3.0.12/python/, with some fixes.
//
// Also merge in change in the git repository
// (https://github.com/swig/swig) that haven't made it into the SWIG
// released software yet. This is from
// d9cac176f6a03b9a57527ef5c1b5659b8857242a.
//
// The fixes are related to proper handling of shared_ptr returned by
// directors.
%include <shared_ptr.i>

%{
// See DirectorNotes.md for discussion of this.
#include "python_ref_ptr_cleanup.h"
%}
// Set SHARED_PTR_DISOWN to $disown if required, for example
// #define SHARED_PTR_DISOWN $disown
#if !defined(SHARED_PTR_DISOWN)
#define SHARED_PTR_DISOWN 0
#endif

%fragment("SWIG_null_deleter_python", "header", fragment="SWIG_null_deleter") {
%#define SWIG_NO_NULL_DELETER_SWIG_BUILTIN_INIT
}

// Language specific macro implementing all the customisations for handling the smart pointer
%define SWIG_SHARED_PTR_TYPEMAPS(CONST, TYPE...)

// %naturalvar is as documented for member variables
%naturalvar TYPE;
%naturalvar SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >;

// destructor wrapper customisation
%feature("unref") TYPE
//"if (debug_shared) { cout << \"deleting use_count: \" << (*smartarg1).use_count() << \" [\" << (boost::get_deleter<SWIG_null_deleter>(*smartarg1) ? std::string(\"CANNOT BE DETERMINED SAFELY\") : ( (*smartarg1).get() ? (*smartarg1)->getValue() : std::string(\"NULL PTR\") )) << \"]\" << endl << flush; }\n"
                               "(void)arg1; delete smartarg1;"

// Typemap customisations...

// plain value
%typemap(in) CONST TYPE (void *argp, int res = 0) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (!argp) {
    %argument_nullref("$type", $symname, $argnum);
  } else {
    $1 = *(%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *)->get());
    if (newmem & SWIG_CAST_NEW_MEMORY) delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
  }
}
%typemap(out) CONST TYPE {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(new $1_ltype(($1_ltype &)$1));
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

%typemap(varin) CONST TYPE {
  void *argp = 0;
  int newmem = 0;
  int res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %variable_fail(res, "$type", "$name");
  }
  if (!argp) {
    %variable_nullref("$type", "$name");
  } else {
    $1 = *(%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *)->get());
    if (newmem & SWIG_CAST_NEW_MEMORY) delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
  }
}
%typemap(varout) CONST TYPE {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(new $1_ltype(($1_ltype &)$1));
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

// plain pointer
// Note: $disown not implemented by default as it will lead to a memory leak of the shared_ptr instance
%typemap(in) CONST TYPE * (void  *argp = 0, int res = 0, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > tempshared, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartarg = 0) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SHARED_PTR_DISOWN | %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    tempshared = *%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    $1 = %const_cast(tempshared.get(), $1_ltype);
  } else {
    smartarg = %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    $1 = %const_cast((smartarg ? smartarg->get() : 0), $1_ltype);
  }
}

%typemap(out, fragment="SWIG_null_deleter_python") CONST TYPE * {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = $1 ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >($1 SWIG_NO_NULL_DELETER_$owner) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), $owner | SWIG_POINTER_OWN));
}

%typemap(varin) CONST TYPE * {
  void *argp = 0;
  int newmem = 0;
  int res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %variable_fail(res, "$type", "$name");
  }
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > tempshared;
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartarg = 0;
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    tempshared = *%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    $1 = %const_cast(tempshared.get(), $1_ltype);
  } else {
    smartarg = %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    $1 = %const_cast((smartarg ? smartarg->get() : 0), $1_ltype);
  }
}
%typemap(varout, fragment="SWIG_null_deleter_python") CONST TYPE * {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = $1 ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >($1 SWIG_NO_NULL_DELETER_0) : 0;
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

// plain reference
%typemap(in) CONST TYPE & (void  *argp = 0, int res = 0, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > tempshared) {
  int newmem = 0;
  // Added mms
  // First check to see if all ready pointer type
  TYPE *ptr;
  res = SWIG_ConvertPtrAndOwn($input, (void**)(&ptr), $descriptor(TYPE *), %convertptr_flags, &newmem);
  if (SWIG_IsOK(res)) {
    $1 = ptr;
  } else {
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (!argp) { %argument_nullref("$type", $symname, $argnum); }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    tempshared = *%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    $1 = %const_cast(tempshared.get(), $1_ltype);
  } else {
    $1 = %const_cast(%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *)->get(), $1_ltype);
  }
  }
}
%typemap(out, fragment="SWIG_null_deleter_python") CONST TYPE & {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >($1 SWIG_NO_NULL_DELETER_$owner);
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

%typemap(varin) CONST TYPE & {
  void *argp = 0;
  int newmem = 0;
  int res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %variable_fail(res, "$type", "$name");
  }
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > tempshared;
  if (!argp) {
    %variable_nullref("$type", "$name");
  }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    tempshared = *%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    $1 = *%const_cast(tempshared.get(), $1_ltype);
  } else {
    $1 = *%const_cast(%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *)->get(), $1_ltype);
  }
}
%typemap(varout, fragment="SWIG_null_deleter_python") CONST TYPE & {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(&$1 SWIG_NO_NULL_DELETER_0);
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

// plain pointer by reference
// Note: $disown not implemented by default as it will lead to a memory leak of the shared_ptr instance
%typemap(in) TYPE *CONST& (void  *argp = 0, int res = 0, $*1_ltype temp = 0, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > tempshared) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SHARED_PTR_DISOWN | %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    tempshared = *%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    delete %reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
    temp = %const_cast(tempshared.get(), $*1_ltype);
  } else {
    temp = %const_cast(%reinterpret_cast(argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *)->get(), $*1_ltype);
  }
  $1 = &temp;
}
%typemap(out, fragment="SWIG_null_deleter_python") TYPE *CONST& {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = *$1 ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(*$1 SWIG_NO_NULL_DELETER_$owner) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

%typemap(varin) TYPE *CONST& %{
#error "varin typemap not implemented"
%}
%typemap(varout) TYPE *CONST& %{
#error "varout typemap not implemented"
%}

// shared_ptr by value
%typemap(in) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > (void *argp, int res = 0) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (argp) $1 = *(%reinterpret_cast(argp, $&ltype));
  if (newmem & SWIG_CAST_NEW_MEMORY) delete %reinterpret_cast(argp, $&ltype);
  // Added mms
  // Special handling if this is a director class.
  // See DirectorNotes.md for discussion of this.
  Swig::Director* dp = dynamic_cast<Swig::Director*>($1.get());
  if(dp) {
    $1.reset($1.get(), PythonRefPtrCleanup(dp->swig_get_self()));
  }
}
%typemap(out) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = $1 ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >($1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

// Added mms
%typemap(directorout,noblock=1,implicitconv=1) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > (void * swig_argp, int swig_res = 0) {
  int newmem = 0;
  swig_res = SWIG_ConvertPtrAndOwn($input,&swig_argp,$&descriptor, %convertptr_flags | %implicitconv_flag, &newmem);
  if (!SWIG_IsOK(swig_res)) {
    %dirout_fail(swig_res,"$type");
  }
  $result = *(%reinterpret_cast(swig_argp, $&ltype));
  // Special handling if this is a director class. In that case, we
  // don't own the underlying python object. Instead,
  // we tell python we have a reference to the underlying object, and
  // when this gets destroyed we decrement the reference to the python
  // object. 
  Swig::Director* dp = dynamic_cast<Swig::Director*>($result.get());
  if(dp) {
    $result.reset($result.get(), PythonRefPtrCleanup(dp->swig_get_self()));
    }
  if (newmem & SWIG_CAST_NEW_MEMORY) %delete(%reinterpret_cast(swig_argp, $&ltype));
}

%typemap(directorout,noblock=1,implicitconv=1) CONST TYPE (void * swig_argp, int swig_res = 0) {
  int newmem = 0;
  swig_res = SWIG_ConvertPtrAndOwn($input,&swig_argp,$descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags | %implicitconv_flag, &newmem);
  if (!SWIG_IsOK(swig_res)) {
    %dirout_fail(swig_res,"$type");
  }
  $result = *(%reinterpret_cast(swig_argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *)->get());
  if (newmem & SWIG_CAST_NEW_MEMORY) delete %reinterpret_cast(swig_argp, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *);
}

%typemap(directorin) CONST TYPE& {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(($1_ltype)&$1_name, SWIG_null_deleter());
  $input = SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN);
}

%typemap(directorin) CONST TYPE {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(&$1_name, SWIG_null_deleter());
  $input = SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN);
}

%typemap(varin) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > {
  int newmem = 0;
  void *argp = 0;
  int res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %variable_fail(res, "$type", "$name");
  }
  $1 = argp ? *(%reinterpret_cast(argp, $&ltype)) : SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE >();
  if (newmem & SWIG_CAST_NEW_MEMORY) delete %reinterpret_cast(argp, $&ltype);
}
%typemap(varout) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = $1 ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >($1) : 0;
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}


// shared_ptr by reference
%typemap(in) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > & (void *argp, int res = 0, $*1_ltype tempshared, $*1_ltype temp2shared) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    if (argp) tempshared = *%reinterpret_cast(argp, $ltype);
    delete %reinterpret_cast(argp, $ltype);
    $1 = &tempshared;
  } else {
    $1 = (argp) ? %reinterpret_cast(argp, $ltype) : &tempshared;
  }
  // Added mms
  // Special handling if this is a director class.
  // See DirectorNotes.md for discussion of this.
  Swig::Director* dp = dynamic_cast<Swig::Director*>($1->get());
  if(dp) {
    temp2shared.reset($1->get(), PythonRefPtrCleanup(dp->swig_get_self()));
    $1 = &temp2shared;
  }
}
%typemap(out) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > & {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = *$1 ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(*$1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

// shared_ptr by reference, where we pass the underlying object
// pointed to in python for a director object. This is the same for a
// shared pointer that is not a director, but in the case of the
// director we pass back the shared_ptr contained in the underlying
// python object. Note that this should *not* be used for anything
// that we keep a copy for (e.g., use to set a shared_ptr in an
// object). The primary used of this is weak_ptr - since we normally
// return a *different* shared_ptr (although pointing to the same
// object) this breaks the weak_ptr semantics. See
// DirectorNotes.md for a bit more of a discussion of this,
// and the DirectorExampleWeakPtr unit test in swig-rules-skeleton
// repository pythnon unit tests for an example illustrating the
// issue.
//
// Note that it is important that SHARED_PTR_NO_OWN not actually be copied to a persistent 
// shared_ptr. If you do, then you are likely to get a seg fault - this is one of those
// "giving you enough rope to hang yourself" sort of thing. This
// handles a corner case, but by *breaking* in a controlled way the
// lifetime management that we have in place.


%typemap(in) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > & SHARED_PTR_NO_OWN (void *argp, int res = 0, $*1_ltype tempshared) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    if (argp) tempshared = *%reinterpret_cast(argp, $ltype);
    delete %reinterpret_cast(argp, $ltype);
    $1 = &tempshared;
  } else {
    $1 = (argp) ? %reinterpret_cast(argp, $ltype) : &tempshared;
  }
  // Added mms
  // Special handling if this is a director class.
  // See DirectorNotes.md for discussion of this.
  //
  // Note that it is important that SHARED_PTR_NO_OWN not actually be copied to a persistent 
  // shared_ptr. If you do, then you are likely to get a seg fault - this is one of those
  // "giving you enough rope to hang yourself" sort of thing. This
  // handles a corner case, but by *breaking* in a controlled way the
  // lifetime management that we have in place.
  
  Swig::Director* dp = dynamic_cast<Swig::Director*>($1->get());
  if(dp) {
    PyObject* this_pobj = PyObject_GetAttr($input, PyString_FromString("this"));
    res = SWIG_ConvertPtrAndOwn(this_pobj, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
    Py_DECREF(this_pobj);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (newmem & SWIG_CAST_NEW_MEMORY) {
      if (argp) tempshared = *%reinterpret_cast(argp, $ltype);
      delete %reinterpret_cast(argp, $ltype);
      $1 = &tempshared;
    } else {
      $1 = (argp) ? %reinterpret_cast(argp, $ltype) : &tempshared;
    }
  }    
}

%typemap(varin) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > & %{
#error "varin typemap not implemented"
%}
%typemap(varout) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > & %{
#error "varout typemap not implemented"
%}

// shared_ptr by pointer
%typemap(in) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > * (void *argp, int res = 0, $*1_ltype tempshared) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (newmem & SWIG_CAST_NEW_MEMORY) {
    if (argp) tempshared = *%reinterpret_cast(argp, $ltype);
    delete %reinterpret_cast(argp, $ltype);
    $1 = &tempshared;
  } else {
    $1 = (argp) ? %reinterpret_cast(argp, $ltype) : &tempshared;
  }
}
%typemap(out) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > * {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = ($1 && *$1) ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(*$1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
  if ($owner) delete $1;
}

%typemap(varin) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > * %{
#error "varin typemap not implemented"
%}
%typemap(varout) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > * %{
#error "varout typemap not implemented"
%}

// shared_ptr by pointer reference
%typemap(in) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *& (void *argp, int res = 0, SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > tempshared, $*1_ltype temp = 0) {
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  if (argp) tempshared = *%reinterpret_cast(argp, $*ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete %reinterpret_cast(argp, $*ltype);
  temp = &tempshared;
  $1 = &temp;
}
%typemap(out) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *& {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = (*$1 && **$1) ? new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(**$1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

%typemap(varin) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *& %{
#error "varin typemap not implemented"
%}
%typemap(varout) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *& %{
#error "varout typemap not implemented"
%}

// Typecheck typemaps
// Note: SWIG_ConvertPtr with void ** parameter set to 0 instead of using SWIG_ConvertPtrAndOwn, so that the casting
// function is not called thereby avoiding a possible smart pointer copy constructor call when casting up the inheritance chain.
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER, equivalent="TYPE *", noblock=1)
                      TYPE CONST,
                      TYPE CONST &,
                      TYPE CONST *,
                      TYPE *CONST&,
                      SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >,
                      SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > &,
                      SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *,
                      SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *& {
  int res = SWIG_ConvertPtr($input, 0, $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), 0);
  $1 = SWIG_CheckState(res);
}

// Multiple Output typemap

%typemap(in,numinputs=0) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > &OUTPUT (SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > temp) {
  $1 = &temp;
}

%typemap(argout) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > &OUTPUT {
  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *smartresult = new SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >(*$1);
  %append_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< TYPE > *), SWIG_POINTER_OWN));
}

// various missing typemaps - If ever used (unlikely) ensure compilation error rather than runtime bug
%typemap(in) CONST TYPE[], CONST TYPE[ANY], CONST TYPE (CLASS::*) %{
#error "typemaps for $1_type not available"
%}
%typemap(out) CONST TYPE[], CONST TYPE[ANY], CONST TYPE (CLASS::*) %{
#error "typemaps for $1_type not available"
%}

%typemap(doctype) SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >,
                  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > &,
                  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *,
                  SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE > *&
  %{TYPE%}


%template() SWIG_SHARED_PTR_QNAMESPACE::shared_ptr< CONST TYPE >;


%enddef
