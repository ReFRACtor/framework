#ifndef PCA_EIGENSOLVER_H
#define PCA_EIGENSOLVER_H

#include <vector>
#include <blitz/array.h>

namespace FullPhysics {

class PCAEigenSolver {
public:

    ~PCAEigenSolver() = default;

    /// Returns the lograithmic mean of the input data averaged over their second dimension
    /// Each array in the vector will have the size of the first dimension of the data
    /// given to the class originally in the corresponding data index.
    virtual std::vector<blitz::Array<double, 1> > data_mean() const = 0;

    /// Returns the EOF perturbations of the input data's mean values for the +/- EOFs properties
    /// Each array in the vector is sized NUM_DATA * PLUS_MINUS * NUM_EOFS*2 
    /// Where NUM_DATA is the size of the original data for the corresponding data index
    /// Where PLUS_MINUS is always 2. Index 0 is the plus (+) EOFs, index 1 is the minus (1) EOFs 
    virtual std::vector<blitz::Array<double, 3> > data_perturbations() const;

    /// Returns the EOF properties for the packed data
    /// Dimensions are NUM_EOF * NUM_PACKED
    virtual std::vector<blitz::Array<double, 2> > eof_properties() const = 0;

    /// Returns the prinicpal components for each EOF 
    /// Dimensions are NUM_EOF * NUM_POINTS
    virtual blitz::Array<double, 2> principal_components() const = 0;

    /// Compute the radiance correction factor given outsputs from the various RT and using the prinicpal compontents computed
    /// by this class
    virtual blitz::Array<double, 2> correction(
        const blitz::Array<double, 1>& lidort_mean, const blitz::Array<double, 1>& twostream_mean, const blitz::Array<double, 1>& first_order_man,
        const blitz::Array<double, 2>& lidort_plus, const blitz::Array<double, 2>& twostream_plus, const blitz::Array<double, 2>& first_order_plus,
        const blitz::Array<double, 2>& lidort_minus, const blitz::Array<double, 2>& twostream_minus, const blitz::Array<double, 2>& first_order_minus);

};

class PCAEigenSolverGeneric : public PCAEigenSolver {
public:
    
    PCAEigenSolverGeneric(const std::vector<blitz::Array<double, 2> >& gridded_data, int num_eofs);

    virtual std::vector<blitz::Array<double, 1> > data_mean() const { return atmos_mean_; }
    virtual std::vector<blitz::Array<double, 2> > eof_properties() const { return eofs_; }
    virtual blitz::Array<double, 2> principal_components() const { return prin_comps_; };

private:

    void solve(const std::vector<blitz::Array<double, 2> >& gridded_data, int num_eofs);

    // Result variables
    std::vector<blitz::Array<double, 1> > atmos_mean_;
    std::vector<blitz::Array<double, 2> > eofs_;
    blitz::Array<double, 2> prin_comps_;

};

class PCAEigenSolverFortran : public PCAEigenSolver {
public:
 
    PCAEigenSolverFortran(const blitz::Array<double, 3>& gridded_data, int num_eofs);

    virtual std::vector<blitz::Array<double, 1> > data_mean() const;
    virtual std::vector<blitz::Array<double, 2> > eof_properties() const;
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
