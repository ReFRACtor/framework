#ifndef LIDORT_INTERFACE_COMMON_H
#define LIDORT_INTERFACE_COMMON_H

#include <iostream>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "fp_exception.h"


/* This file was auto-generated */

#define FORTRAN_TRUE_INT 1

#define BYTE_SIZE_ERROR_CHECK(var_name, c_size, f_size) \
  if(c_size != f_size) { \
    std::stringstream err_msg; \
    err_msg << "Size of C variable: " << c_size \
            << " for " << var_name \
            << " does not match size of Fortran variable: " << f_size; \
    throw Exception(err_msg.str()); \
  }

namespace FullPhysics {

class Lidort_Structure : public Printable<Lidort_Structure> {
public:
  Lidort_Structure() : fortran_type_c(0), owns_pointer(true) {}
  Lidort_Structure(void* allocated_f_type_c) : fortran_type_c(allocated_f_type_c), owns_pointer(false) {}
  void* fortran_type_ptr() { return fortran_type_c; }

  virtual void print(std::ostream &output_stream) const = 0;

protected:
  void *fortran_type_c;
  bool owns_pointer;
};

}
#endif