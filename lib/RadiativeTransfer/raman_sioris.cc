#include "raman_sioris.h"

#include "unit.h"

// for to_fortran
#include "linear_algebra.h"

using namespace blitz;
using namespace FullPhysics;

extern "C" {
    void get_raman(int *nz, int *nw, double *sza, double *vza, double *sca, double *albedo, bool *do_upwelling, const double *ts, const double *rhos, const double *wave, const double *sol, const double *taus, double *rspec, bool *problems);
}

Array<double, 1> FullPhysics::compute_raman_sioris(double solar_zenith, double viewing_zenith, double scattering_angle, double albedo, bool do_upwelling, const Array<double, 1> &temperature_layers, const Array<double, 1>& air_number_density, const SpectralDomain &grid, const Array<double, 1> &solar_irradiance, const Array<double, 2> &total_optical_depth)
{

    int num_layers = temperature_layers.rows(); // nz
    int num_points = grid.data().rows(); // nw

    if (air_number_density.rows() != num_layers) {
        Exception err;
        err << "Number of air_number_density layers: " << air_number_density.rows() 
            << " does not match the number of layers specified by the temperature layers: " << num_layers;
        throw err;
    }
    
    if (total_optical_depth.cols() != num_layers) {
        Exception err;
        err << "Number of total_optical_depth layers: " << total_optical_depth.rows() 
            << " does not match the number of layers specified by the temperature layers: " << num_layers;
        throw err;
    }

    if(solar_irradiance.rows() != num_points) {
        Exception err;
        err << "Number of solar_irradiance grid points: " << solar_irradiance.rows() 
            << " does not match the number of grid points: " << num_points;
        throw err;
    }

    if(total_optical_depth.rows() != num_points) {
        Exception err;
        err << "Number of total_optical_depth grid points: " << total_optical_depth.rows() 
            << " does not match the number of grid points: " << num_points;
        throw err;
    }

    Array<double, 1> nm_grid = grid.convert_wave(units::nm);

    Array<double, 2> total_optical_depth_f = to_fortran(total_optical_depth);

    Array<double, 1> raman_spec(num_points);

    bool problems = false;

    get_raman(&num_layers, &num_points, &solar_zenith, &viewing_zenith, &scattering_angle, &albedo, &do_upwelling, temperature_layers.dataFirst(), air_number_density.dataFirst(), nm_grid.dataFirst(), solar_irradiance.dataFirst(), total_optical_depth_f.dataFirst(), raman_spec.dataFirst(), &problems);

    if (problems) {
        throw Exception("raman_sioris fortran code encountered an internal issue.");
    }

    return raman_spec;

}
