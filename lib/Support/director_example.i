// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "director_example.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::DirectorExample);
%fp_shared_ptr(FullPhysics::DirectorExampleUser);
%fp_shared_ptr(FullPhysics::DirectorExampleWeakPtr);

// Note we need the template before its first use, so we have all
// the typemaps in place.
%template(Vector_DirectorExample) std::vector<boost::shared_ptr<FullPhysics::DirectorExample> >;
%template(Vector_Vector_DirectorExample) std::vector<std::vector<boost::shared_ptr<FullPhysics::DirectorExample> > >;
%template(Vector_Vector_int) std::vector<std::vector<int> >;

namespace FullPhysics {
%feature("director") DirectorExample;

class DirectorExample : public GenericObject {
public:
  DirectorExample(int V, const std::string& Name);
  std::string print_to_string() const;
  virtual std::string desc() const;
  virtual int func(int v) const;
  int value_to_add() const;
  void value_to_add(int v);
  const std::string& name() const;
  %pickle_serialization();
protected:  
  // Only meant for swig serialization
  DirectorExample();
};

class DirectorExampleWeakPtr : public GenericObject {
public:
  // See DirectorNotes.md on the used of SHARED_PTR_NO_OWN
  DirectorExampleWeakPtr(const boost::shared_ptr<DirectorExample>& SHARED_PTR_NO_OWN);
  std::string print_to_string() const;
  boost::shared_ptr<DirectorExample> a() const;
  void a(const boost::shared_ptr<DirectorExample>& SHARED_PTR_NO_OWN);
  bool expired() const;
  long use_count() const;
  %pickle_serialization();
};
  
class DirectorExampleUser : public GenericObject {
public:
  DirectorExampleUser();
  std::string print_to_string() const;
  const std::vector<boost::shared_ptr<DirectorExample> >&
    value_vec() const;
  void set_vec(const std::vector<boost::shared_ptr<DirectorExample> >& V);
  const std::vector<std::vector<boost::shared_ptr<DirectorExample> > >&
    value_vec_vec() const;
  void set_vec_vec(const std::vector<std::vector<boost::shared_ptr<DirectorExample> > >& V);
  const boost::shared_ptr<DirectorExample>&
    value() const;
  const std::vector<int>& value_vec_int() const;
  void set_vec_int(const std::vector<int>& V);
  const std::vector<std::vector<int> >& value_vec_vec_int() const;
  void set_vec_vec_int(const std::vector<std::vector<int> >& V);
  void set_value(const boost::shared_ptr<DirectorExample>& V);
  int apply_func(int v) const;
  int apply_func(int v, int ind) const;
  int apply_func(int v, int ind, int ind2) const;
  %pickle_serialization();
};

}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(director_example, DirectorExample)
  
// List of things "import *" will include
%python_export("DirectorExample", "DirectorExampleUser",
	       "DirectorExampleWeakPtr",
	       "Vector_DirectorExample", "Vector_Vector_DirectorExample",
	       "Vector_Vector_int");

  
