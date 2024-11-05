#ifndef DIRECTOR_EXAMPLE_H
#define DIRECTOR_EXAMPLE_H
#include "printable.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include <iostream>

namespace FullPhysics {
/****************************************************************//**
  This is a simple class that we can use to test having a derived
  class in python to override the director functions.
*******************************************************************/

class DirectorExample : public Printable<DirectorExample> {
public:
//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  DirectorExample(int value_to_add, const std::string& name)
    : value_to_add_(value_to_add),
      name_(name)
  {
    std::cerr << "Creating C++ " << name_ << "\n";
  }
  virtual ~DirectorExample()
  {
    std::cerr << "Destroying C++ " << name_ << "\n";
  }
  virtual int func(int v) const { return v + value_to_add_; }
  virtual void print(std::ostream& Os) const;
  int value_to_add() const { return value_to_add_;}
  void value_to_add(int v) { value_to_add_ = v;}
  const std::string& name() const { return name_;}
protected:
  // Directors with swig seems to need to have this available
  DirectorExample() {}
private:
  int value_to_add_;
  std::string name_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  A weak pointer to a DirectorExample, used to test proper handling of
  weak pointers and also pointer tracking in the serialization code.
*******************************************************************/
class DirectorExampleWeakPtr : public Printable<DirectorExampleWeakPtr> {
public:
  DirectorExampleWeakPtr(const boost::shared_ptr<DirectorExample>& V)
    : a_(V) {}
  virtual ~DirectorExampleWeakPtr() {}
  void print(std::ostream& Os) const
  {
    Os << "----------------------\n";
    Os << "DirectorExampleWeakPtr\n";
    if(a())
      Os << *a();
    else
      Os << "None";
    Os << "\n----------------------\n";
  }
  boost::shared_ptr<DirectorExample> a() const { return a_.lock();}
  void a(const boost::shared_ptr<DirectorExample>& V) { a_ = V;}
  bool expired() const {return a_.expired();}
  long use_count() const {return a_.use_count();}
private:
  boost::weak_ptr<DirectorExample> a_;
  DirectorExampleWeakPtr() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  Sample class to use DirectorExample, can use to make sure python
  can be passed correctly.
*******************************************************************/
class DirectorExampleUser : public Printable<DirectorExampleUser> {
public:
  DirectorExampleUser()
  {
    std::cerr << "Creating DirectorExampleUser\n";
  }
  virtual ~DirectorExampleUser()
  {
    std::cerr << "Destroying DirectorExampleUser\n";
  }
  virtual void print(std::ostream& Os) const { Os << "DirectorExampleUser"; }
  const std::vector<boost::shared_ptr<DirectorExample> >&
  value_vec() const { return value_vec_; }
  void set_vec(const std::vector<boost::shared_ptr<DirectorExample> >& V)
  {value_vec_ = V; }
  const std::vector<std::vector<boost::shared_ptr<DirectorExample> > >&
  value_vec_vec() const { return value_vec2_; }
  void set_vec_vec(const std::vector<std::vector<boost::shared_ptr<DirectorExample> > >& V)
  {value_vec2_ = V; }
  const std::vector<int>&
  value_vec_int() const { return t_; }
  void set_vec_int(const std::vector<int>& V)
  {t_ = V; }
  const std::vector<std::vector<int> >&
  value_vec_vec_int() const { return t2_; }
  void set_vec_vec_int(const std::vector<std::vector<int> >& V)
  {t2_ = V; }
  const boost::shared_ptr<DirectorExample>&
  value() const { return value_; }
  void set_value(const boost::shared_ptr<DirectorExample>& V)
  {value_ = V; }
  int apply_func(int v) const { return value_->func(v);}
  int apply_func(int v, int ind) const { return value_vec_[ind]->func(v);}
  int apply_func(int v, int ind, int ind2) const { return value_vec2_[ind][ind2]->func(v);}
private:
  std::vector<boost::shared_ptr<DirectorExample> >  value_vec_;
  std::vector<std::vector<boost::shared_ptr<DirectorExample> > > value_vec2_;
  boost::shared_ptr<DirectorExample> value_;
  std::vector<int> t_;
  std::vector<std::vector<int> > t2_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
  
}
FP_EXPORT_KEY(DirectorExample);
FP_EXPORT_KEY(DirectorExampleWeakPtr);
FP_EXPORT_KEY(DirectorExampleUser);
#endif

