// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "pca_eigensolver.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::PCAEigenSolver);
%fp_shared_ptr(FullPhysics::PCAEigenSolverGeneric);
%fp_shared_ptr(FullPhysics::PCAEigenSolverFortran);

namespace FullPhysics {

class PCAEigenSolver : public GenericObject {
public:
  virtual std::vector<blitz::Array<double, 1> > data_mean() const = 0;
  virtual std::vector<blitz::Array<double, 2> > eof_properties() const;
  virtual blitz::Array<double, 2> principal_components() const;
  virtual blitz::Array<double, 2> correction_2m(
        const blitz::Array<double, 1>& lidort_mean, const blitz::Array<double, 1>& twostream_mean,
        const blitz::Array<double, 2>& lidort_plus, const blitz::Array<double, 2>& twostream_plus,
        const blitz::Array<double, 2>& lidort_minus, const blitz::Array<double, 2>& twostream_minus);
  virtual blitz::Array<double, 2> correction_3m(
        const blitz::Array<double, 1>& lidort_mean,
        const blitz::Array<double, 1>& twostream_mean,
        const blitz::Array<double, 1>& first_order_man,
        const blitz::Array<double, 2>& lidort_plus,
        const blitz::Array<double, 2>& twostream_plus,
        const blitz::Array<double, 2>& first_order_plus,
        const blitz::Array<double, 2>& lidort_minus,
        const blitz::Array<double, 2>& twostream_minus,
        const blitz::Array<double, 2>& first_order_minus);
};

class PCAEigenSolverGeneric : public PCAEigenSolver {
public:
  PCAEigenSolverGeneric(const std::vector<blitz::Array<double, 2> >&
                        gridded_data, int num_eofs);
  virtual std::vector<blitz::Array<double, 1> > data_mean() const;
};

class PCAEigenSolverFortran : public PCAEigenSolver {
public:
  PCAEigenSolverFortran(const blitz::Array<double, 3>& gridded_data,
                        int num_eofs);
  virtual std::vector<blitz::Array<double, 1> > data_mean() const;
};

}

