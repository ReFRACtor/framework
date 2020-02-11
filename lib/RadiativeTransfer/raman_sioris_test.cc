#include "unit.h"
#include "unit_test_support.h"

#include "hdf_file.h"

extern "C" {
    void get_raman(int *nz, int *nw, double *sza, double *vza, double *sca, double *albedo, bool *do_upwelling, double *ts, double *rhos, double *wave, double *sol, double *taus, double *rspec, bool *problems);
}

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(raman_sioris, GlobalFixture)

BOOST_AUTO_TEST_CASE(offline_data)
{
    // Unit test data was converted from test_raman_input.dat file from MUSES OMI code in idl-retrieve
    // that was originally provided by Xiong Liu
    HdfFile offline_data(test_data_dir() + "/in/raman_sioris/offline_test_data.h5");

    int nw = offline_data.read_field<int>("/nw");
    int nz = offline_data.read_field<int>("/nz");
    double sza = offline_data.read_field<double>("/sza");
    double vza = offline_data.read_field<double>("/vza");
    double sca = offline_data.read_field<double>("/sca");

    Array<double, 1> ts = offline_data.read_field<double, 1>("/ts");
    Array<double, 1> rhos = offline_data.read_field<double, 1>("/rhos");

    Array<double, 1> wave = offline_data.read_field<double, 1>("/wave");
    Array<double, 1> sol = offline_data.read_field<double, 1>("/sol");

    Array<double, 2> taus_c = offline_data.read_field<double, 2>("/taus");

    // Convert to column major ordering
    Array<double, 2> taus_f(to_fortran(taus_c));

    bool do_upwelling = true;
    double albedo = 0.5;

    Array<double, 1> rspec_calc(nw, blitz::ColumnMajorArray<1>());
    bool problems = false;

    get_raman(&nz, &nw, &sza, &vza, &sca, &albedo, &do_upwelling, ts.dataFirst(), rhos.dataFirst(), wave.dataFirst(), sol.dataFirst(), taus_f.dataFirst(), rspec_calc.dataFirst(), &problems);

    Array<double, 1> rspec_expt = offline_data.read_field<double, 1>("rspec");

    BOOST_CHECK_MATRIX_CLOSE_TOL(rspec_expt, rspec_expt, 1e-16);

}

BOOST_AUTO_TEST_SUITE_END()
