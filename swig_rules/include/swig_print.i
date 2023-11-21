//--------------------------------------------------------------
// We use the generic name "print_to_string" in our classes, but
// this corresponds to the __str__ object in python. So go ahead
// and rename everything if we are in python
//--------------------------------------------------------------

%rename(__str__) print_to_string;

// Ignore the print function. Python can't really override this
// effectively, instead it should override desc. The print function
// interacts incorrectly with our ostream swig rules in swig_iostream
//%rename(print_desc) print;
%rename ("$ignore") print;

// Trying to override print in python through a director doesn't
// really work very well. So we instroduce a desc function that can be
// overriden.
%define %director_python_print()
%extend {
  virtual std::string desc() const { return "SwigObject"; }
  virtual void print(std::ostream& Os) const { Os << desc(); }
}
%enddef

