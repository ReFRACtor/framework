#ifndef PCA_EIGENSOLVER_H
#define PCA_EIGENSOLVER_H

#include <vector>
#include <blitz/array.h>

namespace FullPhysics {

class PCAEigenSolver {
public:

    ~PCAEigenSolver() = default;

    virtual std::vector<blitz::Array<double, 1> > data_mean() const = 0;
    virtual blitz::Array<double, 2> eofs() const = 0;
    virtual blitz::Array<double, 2> principal_components() const = 0;
};

class PCAEigenSolverGeneric : public PCAEigenSolver {
public:
    
    PCAEigenSolverGeneric(const std::vector<blitz::Array<double, 2> >& gridded_data, int num_eofs);

    virtual std::vector<blitz::Array<double, 1> > data_mean() const { return atmos_mean_; }
    virtual blitz::Array<double, 2> eofs() const { return eofs_; }
    virtual blitz::Array<double, 2> principal_components() const { return prin_comps_; };

private:

    void solve(const std::vector<blitz::Array<double, 2> >& gridded_data, int num_eofs);

    // Result variables
    std::vector<blitz::Array<double, 1> > atmos_mean_;
    blitz::Array<double, 2> eofs_;
    blitz::Array<double, 2> prin_comps_;

};

class PCAEigenSolverFortran : public PCAEigenSolver {
public:
 
    PCAEigenSolverFortran(const blitz::Array<double, 3>& gridded_data, int num_eofs);

    virtual std::vector<blitz::Array<double, 1> > data_mean() const;
    virtual blitz::Array<double, 2> eofs() const { return eofs_; }
    virtual blitz::Array<double, 2> principal_components() const { return prin_comps_; };

private:

    void solve(const blitz::Array<double, 3>& gridded_data, int num_eofs);

    // Result variables
    blitz::Array<double, 2> atmos_mean_;
    blitz::Array<double, 2> eofs_;
    blitz::Array<double, 2> prin_comps_;

};

}

#endif
