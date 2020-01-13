#include "ostream_pad.h"

#include "pca_rt.h"

#include "optical_properties_wrt_rt.h"
#include "optical_properties_pca.h"

using namespace FullPhysics;
using namespace blitz;

PCARt::PCARt(const boost::shared_ptr<AtmosphereStandard>& Atm,
             const std::string Primary_absorber,
             const PCABinning::Method Bin_method, const int Num_bins,
             const int Num_eofs,
             const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
             const blitz::Array<double, 1>& Sza, 
             const blitz::Array<double, 1>& Zen, 
             const blitz::Array<double, 1>& Azm,
             int Number_streams, 
             int Number_moments, 
             bool do_solar_sources, 
             bool do_thermal_emission) 
: RadiativeTransferFixedStokesCoefficient(Stokes_coef), atm(Atm), primary_absorber(Primary_absorber), 
  bin_method(Bin_method), num_bins(Num_bins), num_eofs(Num_eofs)
{
    lidort_rt.reset(new LidortRt(Atm, Stokes_coef, Sza, Zen, Azm, false, Number_streams, Number_moments, true, do_solar_sources, do_thermal_emission, do_thermal_emission));

    twostream_rt.reset(new TwostreamRt(Atm, Stokes_coef, Sza, Zen, Azm, true, do_solar_sources, do_thermal_emission));

    first_order_rt.reset(new FirstOrderRt(Atm, Stokes_coef, Sza, Zen, Azm, Number_streams, Number_moments, do_solar_sources, do_thermal_emission));

}

void PCARt::compute_bins(const SpectralDomain& Spec_domain, int Spec_index) const
{
    // Compute optical properties for the whole band
    pca_opt.clear();

    for (int dom_idx = 0; dom_idx < Spec_domain.data().rows(); dom_idx++) {
        boost::shared_ptr<OpticalPropertiesWrtRt> point_opt_prop(new OpticalPropertiesWrtRt());
        point_opt_prop->initialize(DoubleWithUnit(Spec_domain.data()(dom_idx), Spec_domain.units()), Spec_index,
                                   atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());
        pca_opt.push_back(point_opt_prop);
    }

    int primary_abs_index = atm->absorber_ptr()->gas_index(primary_absorber);

    // Compute bins
    pca_bin.reset(new PCABinning(pca_opt, bin_method, num_bins, primary_abs_index));
}

blitz::Array<double, 2> PCARt::stokes(const SpectralDomain& Spec_domain, int Spec_index) const
{

    // Our final output from this routine
    Array<double, 2> stokes_output(Spec_domain.data().rows(), number_stokes());

    // Convenience since other RTs take wavenumbers as input
    Array<double, 1> wavenumbers(Spec_domain.wavenumber());

    // Need to convert aerosol from atmosphere to AerosolOptical for PCA optical properties packing
    boost::shared_ptr<AerosolOptical> aerosol_optical = boost::dynamic_pointer_cast<AerosolOptical>(atm->aerosol_ptr());

    if(!aerosol_optical) {
        throw Exception("Failed to convert aerosol class to AerosolOptical");
    }

    compute_bins(Spec_domain, Spec_index);

    auto bins = pca_bin->bin_indexes();
    Array<int, 1> num_bin_points = pca_bin->num_bin_points();

    // Extract optical properties intermediate variable values for each bin
    int num_gas = atm->absorber_ptr()->number_species();
    int num_aerosol = atm->aerosol_ptr()->number_particle();
    int num_layers = pca_opt[0]->number_layers();

    // Compute values for each bin
    pca_bin_solvers.clear();
    for (int bin_idx = 0; bin_idx < bins.size(); bin_idx++) {

        // Skip empty bins
        if (num_bin_points(bin_idx) == 0) {
            continue;
        }

        // Solve PCA solution for intermediate variable optical properties for the bin

        // Copy out the intermediate variable values for the spectral domain values in the current bin
        // This requires us to copy into a new array since the bin locations won't be contiguous
        std::vector<blitz::Array<double, 2> > data_list;

        ArrayAd<double, 2> packed_0 = OpticalPropertiesPca::pack(pca_opt[bins[bin_idx](0)]);
        int num_pack_var = packed_0.cols();

        for (int var_idx = 0; var_idx < num_pack_var; var_idx++) {
            blitz::Array<double, 2> var_packed(num_layers, num_bin_points(bin_idx));

            var_packed(Range::all(), 0) = packed_0.value()(Range::all(), var_idx);

            data_list.push_back(var_packed);
        }

        for(int dom_idx = 1; dom_idx < num_bin_points(bin_idx); dom_idx++) {
            ArrayAd<double, 2> point_packed = OpticalPropertiesPca::pack(pca_opt[bins[bin_idx](dom_idx)]);

            for (int var_idx = 0; var_idx < num_pack_var; var_idx++) {
                data_list[var_idx](Range::all(), dom_idx) = point_packed.value()(Range::all(), var_idx);
            }
        }

        // Create bin averaged spectral point
        Array<double, 1> bin_wns(num_bin_points(bin_idx));
        for(int dom_idx = 0; dom_idx < num_bin_points(bin_idx); dom_idx++) {
            bin_wns(dom_idx) = wavenumbers(bins[bin_idx](dom_idx));
        }
        double avg_wn = mean(bin_wns);

        // Create PCA solution from binned intermediate variables
        boost::shared_ptr<PCAEigenSolver> pca_solver;
        pca_solver.reset(new PCAEigenSolverGeneric(data_list, num_eofs));
        pca_bin_solvers.push_back(pca_solver);
        
        // Create EOF intermediate values for mean and plus and minus each eof
        // We need to translate back from the eigen solvers vectors to arrays for each intermediate variable index
        
        // Compute mean IV values
        Array<double, 2> bin_mean_iv(num_layers, num_pack_var);
        for (int var_idx = 0; var_idx < num_pack_var; var_idx++) {
            bin_mean_iv(Range::all(), var_idx) = exp(pca_solver->data_mean()[var_idx]);
        }
        
        // Compute bin mean intensity from each RT 
        ArrayAd<double, 2> bin_mean_iv_ad(bin_mean_iv);
        boost::shared_ptr<OpticalPropertiesPca> bin_mean_opt_props(new OpticalPropertiesPca(bin_mean_iv_ad, avg_wn, aerosol_optical, num_gas, num_aerosol));

        Array<double, 1> lidort_mean(lidort_rt->stokes_single_wn(avg_wn, Spec_index, bin_mean_opt_props));
        Array<double, 1> twostream_mean(twostream_rt->stokes_single_wn(avg_wn, Spec_index, bin_mean_opt_props));
        Array<double, 1> first_order_mean(first_order_rt->stokes_single_wn(avg_wn, Spec_index, bin_mean_opt_props));

        // Compute bin mean plus eof and mean minus eof intensity from each RT
        Array<double, 2> lidort_plus(num_eofs, number_stokes());
        Array<double, 2> twostream_plus(num_eofs, number_stokes());
        Array<double, 2> first_order_plus(num_eofs, number_stokes());

        Array<double, 2> lidort_minus(num_eofs, number_stokes());
        Array<double, 2> twostream_minus(num_eofs, number_stokes());
        Array<double, 2> first_order_minus(num_eofs, number_stokes());

        for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
            Array<double, 2> bin_plus_iv(num_layers, num_pack_var);
            Array<double, 2> bin_minus_iv(num_layers, num_pack_var);

            for (int var_idx = 0; var_idx < num_pack_var; var_idx++) {
                bin_plus_iv(Range::all(), var_idx) = 
                    exp(pca_solver->data_mean()[var_idx] + pca_solver->eof_properties()[var_idx](Range::all(), eof_idx));

                bin_minus_iv(Range::all(), var_idx) = 
                    exp(pca_solver->data_mean()[var_idx] - pca_solver->eof_properties()[var_idx](Range::all(), eof_idx));
            }

            ArrayAd<double, 2> bin_plus_iv_ad(bin_plus_iv);
            ArrayAd<double, 2> bin_minus_iv_ad(bin_minus_iv);

            boost::shared_ptr<OpticalPropertiesPca> bin_plus_opt_props(new OpticalPropertiesPca(bin_plus_iv_ad, avg_wn, aerosol_optical, num_gas, num_aerosol));
            boost::shared_ptr<OpticalPropertiesPca> bin_minus_opt_props(new OpticalPropertiesPca(bin_minus_iv_ad, avg_wn, aerosol_optical, num_gas, num_aerosol));

            lidort_plus(eof_idx, Range::all()) = lidort_rt->stokes_single_wn(avg_wn, Spec_index, bin_plus_opt_props);
            lidort_minus(eof_idx, Range::all()) = lidort_rt->stokes_single_wn(avg_wn, Spec_index, bin_minus_opt_props);

            twostream_plus(eof_idx, Range::all()) = twostream_rt->stokes_single_wn(avg_wn, Spec_index, bin_plus_opt_props);
            twostream_minus(eof_idx, Range::all()) = twostream_rt->stokes_single_wn(avg_wn, Spec_index, bin_minus_opt_props);

            first_order_plus(eof_idx, Range::all()) = first_order_rt->stokes_single_wn(avg_wn, Spec_index, bin_plus_opt_props);
            first_order_minus(eof_idx, Range::all()) = first_order_rt->stokes_single_wn(avg_wn, Spec_index, bin_minus_opt_props);

        }

        // Compute correction factors
        Array<double, 2> bin_corrections(num_bin_points(bin_idx), number_stokes());
        bin_corrections = pca_solver->correction(
                lidort_mean, twostream_mean, first_order_mean,
                lidort_plus, twostream_plus, first_order_plus,
                lidort_minus, twostream_minus, first_order_minus);

        // Compute 2stream and first order for all points in the bin and use correction to get final values
        for(int dom_idx = 0; dom_idx < num_bin_points(bin_idx); dom_idx++) {
            int grid_idx = bins[bin_idx](dom_idx);
            int dom_wn = wavenumbers(grid_idx);

            Array<double, 1> twostream_full(twostream_rt->stokes_single_wn(dom_wn, Spec_index, pca_opt[grid_idx]));
            Array<double, 1> first_order_full(first_order_rt->stokes_single_wn(dom_wn, Spec_index, pca_opt[grid_idx]));

            stokes_output(grid_idx, Range::all()) = 
                bin_corrections(dom_idx, Range::all()) * (twostream_full(Range::all()) + first_order_full(Range::all()));
        }
    }

    return stokes_output;
}

ArrayAd<double, 2> PCARt::stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const
{
    return ArrayAd<double, 2>(stokes(Spec_domain, Spec_index));
} 

void PCARt::print(std::ostream& Os, bool Short_form) const 
{
    OstreamPad opad(Os, "    ");
    Os << "PCARt\n";
    Os << "  ";
    RadiativeTransferFixedStokesCoefficient::print(opad, Short_form);
    opad.strict_sync();
    Os << "\nPrimary absorber: " << primary_absorber << "\n"
       << "Number EOFs: " << num_eofs << "\n";

    Os << "\n  LIDORT RT:\n";
    lidort_rt->print(opad, Short_form);
    opad.strict_sync();

    Os << "\n  2stream RT:\n";
    twostream_rt->print(opad, Short_form);
    opad.strict_sync();

    Os << "\n  First Order RT:\n";
    first_order_rt->print(opad, Short_form);
    opad.strict_sync();
}
