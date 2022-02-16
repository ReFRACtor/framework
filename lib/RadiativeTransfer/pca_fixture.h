#include "global_fixture.h"
#include "hdf_file.h"
#include "unit_test_support.h"

#include "optical_properties_wrt_rt.h"

using namespace FullPhysics;
using namespace blitz;

class PcaFixture : public GlobalFixture {
public:

    PcaFixture() {
        Range ra = Range::all();

        boost::shared_ptr<HdfFile> prop_data(new HdfFile(test_data_dir() + "in/pca/optical_properties.h5"));

        Array<double, 2> ray_sca = prop_data->read_field<double, 2>("rayleigh_optical_depth");
        Array<double, 3> gas_abs = prop_data->read_field<double, 3>("gas_optical_depth");
        Array<double, 3> aer_ext = prop_data->read_field<double, 3>("aerosol_extinction");
        Array<double, 3> aer_sca = prop_data->read_field<double, 3>("aerosol_scattering");

        Array<double, 2> total_od = prop_data->read_field<double, 2>("total_optical_depth");
        Array<double, 2> total_ssa = prop_data->read_field<double, 2>("single_scattering_albedo");

        for(int dat_idx = 0; dat_idx < ray_sca.rows(); dat_idx++) {
            boost::shared_ptr<OpticalPropertiesWrtRt> dat_props(new OpticalPropertiesWrtRt());
            ArrayAd<double, 1> ray_od_ad(ray_sca(dat_idx, ra));
            ArrayAd<double, 2> gas_od_ad(gas_abs(dat_idx, ra, ra));
            ArrayAd<double, 2> aer_ext_ad(aer_ext(dat_idx, ra, ra));
            ArrayAd<double, 2> aer_sca_ad(aer_sca(dat_idx, ra, ra));

            // Empty phase functions value
            std::vector<ArrayAd<double, 3> > pf;
            for(int aer_idx = 0; aer_idx < aer_ext.depth(); aer_idx++) {
                ArrayAd<double, 3> empty(1, ray_od_ad.rows(), 6, 0);
                pf.push_back(empty);
            }

            dat_props->initialize(ray_od_ad, gas_od_ad, aer_ext_ad, aer_sca_ad, pf);

            // Double check that the values from the input file matches what we compute
            BOOST_CHECK_MATRIX_CLOSE(total_od(dat_idx, ra), dat_props->total_optical_depth().value());
            BOOST_CHECK_MATRIX_CLOSE_TOL(total_ssa(dat_idx, ra), dat_props->total_single_scattering_albedo().value(), 1e-6);

            opt_props.push_back(dat_props);
        }

        num_layers = opt_props[0]->number_layers();
        num_data_points = opt_props.size();

        grid_total_od.resize(num_layers, num_data_points);
        grid_total_ssa.resize(num_layers, num_data_points);

        for(int dat_idx = 0 ; dat_idx < num_data_points;  dat_idx++) {
            grid_total_od(ra, dat_idx) = opt_props[dat_idx]->total_optical_depth().value();
            grid_total_ssa(ra, dat_idx) = opt_props[dat_idx]->total_single_scattering_albedo().value();
        }
    }

    std::vector<Array<double, 2> > gridded_data_generic(int num_points_ret = 0)
    {
        if (num_points_ret < 1) {
            num_points_ret = num_data_points;
        }

        Range ra = Range::all();
        Range r_points(0, num_points_ret-1);

        // Pack data for our generic adaptation class
        std::vector<Array<double, 2> > gridded_data_g;

        gridded_data_g.push_back(grid_total_od(ra, r_points));
        gridded_data_g.push_back(grid_total_ssa(ra, r_points));

        return gridded_data_g;
    }

    Array<double, 3> gridded_data_fortran(int num_points_ret = 0)
    {
        if (num_points_ret < 1) {
            num_points_ret = num_data_points;
        }

        // Pack data for the fortran adaptation class
        // Get sizes from gridded_data_g to avoid rereading arrays from disk
        Array<double, 3> gridded_data_f(2, num_layers, num_points_ret);

        Range ra = Range::all();
        Range r_points(0, num_points_ret-1);

        gridded_data_f(0, ra, ra) = grid_total_od(ra, r_points);
        gridded_data_f(1, ra, ra) = grid_total_ssa(ra, r_points);

        return gridded_data_f;
    }

    std::vector<boost::shared_ptr<OpticalPropertiesWrtRt> > opt_props;

    Array<double, 2> grid_total_od;
    Array<double, 2> grid_total_ssa;

    int num_layers;
    int num_data_points;

};
