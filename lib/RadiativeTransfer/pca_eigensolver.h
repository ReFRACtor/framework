#ifndef PCA_EIGENSOLVER_H
#define PCA_EIGENSOLVER_H

#include <vector>
#include <blitz/array.h>

namespace FullPhysics {

class PCAEigenSolver {
public:
    
    PCAEigenSolver(std::vector<blitz::Array<double, 2> >& gridded_data, int num_eofs);

    std::vector<blitz::Array<double, 1> > data_mean() const { return atmos_mean_; }
    blitz::Array<double, 2> eofs() const { return eofs_; }
    blitz::Array<double, 2> principal_components() const { return prin_comps_; };


private:

    void solve(std::vector<blitz::Array<double, 2> >& gridded_data, int num_eofs);

    // Result variables
    std::vector<blitz::Array<double, 1> > atmos_mean_;
    blitz::Array<double, 2> eofs_;
    blitz::Array<double, 2> prin_comps_;

};

}

#endif
