#include "unit_test_support.h"
#include "pca_fixture.h"

#include "pca_eigensolver.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void pca_2m_correction(int *max_eofs_2p1, int *Max_geoms, int *n_eofs, int *npoints, int *ngeoms, double *PrinComps, double *intensity_LD_bin, double *intensity_2S_bin, double *Intensity_Corrfacs);
    void pca_3m_correction(int *max_eofs_2p1, int *Max_geoms, int *n_eofs, int *npoints, int *ngeoms, double *PrinComps, double *intensity_LD_bin, double *intensity_2S_bin, double *intensity_FO_bin, double *Intensity_Corrfacs);
}

BOOST_FIXTURE_TEST_SUITE(pca_eigensolver, PcaFixture)

BOOST_AUTO_TEST_CASE(compare_with_fortran)
{
    Range ra = Range::all();

    int num_eofs = 4;

    std::vector<Array<double, 2> > gridded_data_g = gridded_data_generic();

    auto eigen_g = PCAEigenSolverGeneric(gridded_data_g, num_eofs);

    Array<double, 3> gridded_data_f = gridded_data_fortran();

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

BOOST_AUTO_TEST_CASE(one_point)
{
    Range ra = Range::all();

    int num_eofs = 4;
    int num_points_used = 1; // Check for errors when using only 1 data point

    std::vector<Array<double, 2> > gridded_data_g = gridded_data_generic(num_points_used);

    auto eigen_g = PCAEigenSolverGeneric(gridded_data_g, num_eofs);

    Array<double, 3> gridded_data_f = gridded_data_fortran(num_points_used);

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

BOOST_AUTO_TEST_CASE(correction)
{

    Range ra = Range::all();

    int num_eofs = 4;

    std::vector<Array<double, 2> > gridded_data_g = gridded_data_generic();

    auto eigen_g = PCAEigenSolverGeneric(gridded_data_g, num_eofs);

    // Use some real values extracted from the offline code, the numbers are
    // not that important given we will be comparing against the fortran code, but
    // its a good source for something with real variations instead of just
    // using some random numbers
    Array<double, 2> intensity_lidort(num_eofs*2+1, 1);
    intensity_lidort = 2.4064995192916996E-005, 2.3588829644955006E-005, 2.4566762449535553E-005, 2.6114228387539818E-005, 2.2182206896868496E-005, 2.5399849528192935E-005, 2.2847515948604265E-005, 2.4658610808848054E-005, 2.3495405676822065E-005;

    Array<double, 2> intensity_twostream(num_eofs*2+1, 1);
    intensity_twostream =  1.1095953790896775E-005, 1.0865986217343071E-005, 1.1340814520655838E-005, 1.2034765021963400E-005, 1.0233404237759563E-005, 1.1735165461935445E-005, 1.0514491219450307E-005, 1.1362391060291389E-005, 1.0841184599708410E-005;

    Array<double, 2> intensity_first_order(num_eofs*2+1, 1);
    intensity_first_order =  6.9918313186754838E-004, 7.0016319085633699E-004, 6.9817588710241810E-004, 7.0011172556585317E-004, 6.9811788370752247E-004, 7.0822910409202525E-004, 6.9018242119303706E-004, 6.9168810758034576E-004, 7.0673637868372610E-004;

    // Extract values from the packed eof arrays used by the fortran code into the seperated values we use in C++
    int nstokes = 1;

    Array<double, 1> lidort_mean(nstokes);
    Array<double, 1> twostream_mean(nstokes);
    Array<double, 1> first_order_mean(nstokes);
 
    Array<double, 2> lidort_plus(num_eofs, nstokes);
    Array<double, 2> twostream_plus(num_eofs, nstokes);
    Array<double, 2> first_order_plus(num_eofs, nstokes);
 
    Array<double, 2> lidort_minus(num_eofs, nstokes);
    Array<double, 2> twostream_minus(num_eofs, nstokes);
    Array<double, 2> first_order_minus(num_eofs, nstokes);

    lidort_mean(0) = intensity_lidort(0);
    twostream_mean(0) = intensity_twostream(0);
    first_order_mean(0) = intensity_first_order(0);

    for (int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
        int plus_idx = 2 * eof_idx + 1;
        int minus_idx = plus_idx + 1;
 
        lidort_plus(eof_idx, 0) = intensity_lidort(plus_idx);
        lidort_minus(eof_idx, 0) = intensity_lidort(minus_idx);
 
        twostream_plus(eof_idx, 0) = intensity_twostream(plus_idx);
        twostream_minus(eof_idx, 0) = intensity_twostream(minus_idx);
 
        first_order_plus(eof_idx, 0) = intensity_first_order(plus_idx);
        first_order_minus(eof_idx, 0) = intensity_first_order(minus_idx);
    }

    // Run our correction implementation and compare against fortran code
    int ngeoms = 1;

    Array<double, 2> prin_comps(num_eofs, num_data_points, ColumnMajorArray<2>());
    prin_comps = eigen_g.principal_components();

    blitz::Array<double, 2> calc_correction_3m(num_data_points, nstokes);
    calc_correction_3m = eigen_g.correction_3m(lidort_mean, twostream_mean, first_order_mean, lidort_plus, twostream_plus, first_order_plus, lidort_minus, twostream_minus, first_order_minus);

    blitz::Array<double, 2> calc_correction_2m(num_data_points, nstokes);
    calc_correction_2m = eigen_g.correction_2m(lidort_mean, twostream_mean, lidort_plus, twostream_plus, lidort_minus, twostream_minus);

    blitz::Array<double, 2> expt_correction_3m(num_data_points, ngeoms, ColumnMajorArray<2>());
    blitz::Array<double, 2> expt_correction_2m(num_data_points, ngeoms, ColumnMajorArray<2>());

    int num_eofs_2p1 = num_eofs * 2 + 1;

    pca_3m_correction(&num_eofs_2p1, &ngeoms, &num_eofs, &num_data_points, &ngeoms, prin_comps.dataFirst(), intensity_lidort.dataFirst(), intensity_twostream.dataFirst(), intensity_first_order.dataFirst(), expt_correction_3m.dataFirst());
    pca_2m_correction(&num_eofs_2p1, &ngeoms, &num_eofs, &num_data_points, &ngeoms, prin_comps.dataFirst(), intensity_lidort.dataFirst(), intensity_twostream.dataFirst(), expt_correction_2m.dataFirst());

    BOOST_CHECK_MATRIX_CLOSE_TOL(expt_correction_3m(Range::all(), 0), calc_correction_3m(Range::all(), 0), 1e-10);
    BOOST_CHECK_MATRIX_CLOSE_TOL(expt_correction_2m(Range::all(), 0), calc_correction_2m(Range::all(), 0), 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
