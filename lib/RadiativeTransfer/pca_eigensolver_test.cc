#include "unit_test_support.h"
#include "global_fixture.h"

#include "pca_optical_properties.h"
#include "pca_eigensolver.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pca_eigensolver, GlobalFixture)

BOOST_AUTO_TEST_CASE(compare_with_fortran)
{
    std::string opt_prop_filename = test_data_dir() + "in/pca/pca_optical_properties.h5";
    auto opt_props = PCAOpticalPropertiesFile(opt_prop_filename);

    int num_eofs = 4;

    // Pack data for our generic adaptation class
    std::vector<Array<double, 2> > gridded_data_g;
    gridded_data_g.push_back(opt_props.total_optical_depth());
    gridded_data_g.push_back(opt_props.single_scattering_albedo());

    auto eigen_g = PCAEigenSolverGeneric(gridded_data_g, num_eofs);

    // Pack data for the fortran adaptation class
    // Get sizes from gridded_data_g to avoid rereading arrays from disk
    int num_layers = gridded_data_g[0].rows();
    int num_points = gridded_data_g[1].cols();
    Array<double, 3> gridded_data_f(2, num_layers, num_points);

    auto all = Range::all();
    gridded_data_f(0, all, all) = opt_props.total_optical_depth();
    gridded_data_f(1, all, all) = opt_props.single_scattering_albedo();

    auto eigen_f = PCAEigenSolverFortran(gridded_data_f, num_eofs);

    int num_var = gridded_data_g.size();
    BOOST_CHECK_MATRIX_CLOSE_TOL(eigen_f.principal_components(), eigen_g.principal_components(), 1e-10);

    for (int ivar = 0; ivar < num_var; ivar++) {
        BOOST_CHECK_MATRIX_CLOSE_TOL(eigen_f.eof_properties()[ivar], eigen_g.eof_properties()[ivar], 1e-10);
        BOOST_CHECK_MATRIX_CLOSE_TOL(eigen_f.data_mean()[ivar], eigen_g.data_mean()[ivar], 1e-10);

        // These should always match since the routine is defined in the base class
        BOOST_CHECK_MATRIX_CLOSE_TOL(eigen_f.data_perturbations()[ivar], eigen_g.data_perturbations()[ivar], 1e-10);
    }
}

BOOST_AUTO_TEST_SUITE_END()
