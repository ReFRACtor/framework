namespace FullPhysics {
// This contains documentation for use by doxygen.
/****************************************************************//**
 \page python_doxygen Python specific documentation
 
  ReFRACtor uses [SWIG](http://www.swig.org/) to build its Python
  bindings from interface (.i) files.

  By default in SWIG, type information is global across all SWIG
  modules loaded, and this can cause type conflicts between modules
  that were not designed to work together. To solve this in ReFRACtor
  , we define SWIG_TYPE_TABLE and set it to "refractor".

  There is a caveat that only other modules compiled with
  SWIG_TYPE_TABLE set to "refractor" will share ReFRACtor type
  information.

  Therefore, if you build a child/derived project that intends to use
  or extend the ReFRACtor Python modules, you MUST set
  SWIG_TYPE_TABLE to "refractor" and ensure your project is with
  compiled with -DSWIG_TYPE_TABLE=refractor.

  For additional information see Section 16.4 [The SWIG runtime code]
  (http://www.swig.org/Doc3.0/Modules.html#Modules_nn2) from the
  [SWIG documentation](http://www.swig.org/Doc3.0/Contents.html)

*******************************************************************/
}
