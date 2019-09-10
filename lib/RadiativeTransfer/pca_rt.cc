#include "pca_rt.h"

using namespace FullPhysics;

PcaRt::PcaRt(const boost::shared_ptr<RtAtmosphere>& Atm,
             const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
             const blitz::Array<double, 1>& Sza, 
             const blitz::Array<double, 1>& Zen, 
             const blitz::Array<double, 1>& Azm,
             int Number_streams, 
             int Number_moments, 
             bool do_solar_sources = true, 
             bool do_thermal_emission = false,
             const std::string Primary_absorber,
             const int Num_eofs) : atm(Atm), primary_absorber(Primary_absorber), num_eofs(Num_eofs)
{
    lidort_rt.reset(LidortRt(Atm, Stokes_coef, Sza, Zen, Azm, false, Number_streams, Number_moments, true, do_solar_sources, do_termal_emission, do_thermal_emission));

    twostream_rt.reset(TwostreamRt(Atm, Stokes_coeff, Sza, Zen, Azm, true, do_solar_sources, do_thermal_emission));

    first_order_rt.reset(FirstOrderRt(Atm, Stokes_coeff, Sza, Zen, Azm, Number_streams, Number_moments, do_solar_sources, do_thermal_emission));

}

blitz::Array<double, 2> PcaRt::stokes(const SpectralDomain& Spec_domain, int Spec_index) const
{

    // Our final output from this routine
    Array<double, 2> stokes_output(Spec_domain.rows(), number_stokes());

    // Convenience since other RTs take wavenumbers as input
    Array<double, 1> wavenumbers(Spec_domain.wavenumber());

    // Compute optical properties for the whole band
    PCAOpticalPropertiesAtmosphere pca_opt(atm, Spec_domain, Spec_index, primary_absorber);

    // Compute bins
    PcaBinning pca_bin(pca_opt, 9); // Only implemented binning method requires 9 bins always
    auto bins = pca_bin.bin_indexes();

    // Compute values for each bin
    for (int bin_idx = 0; bin_idx < bins.size(); bin_idx++) {
        int num_bin_points = bins[var_idx].rows();

        // Solve PCA solution for intermediate variable optical properties for the bin
        Array<double, 3> bin_inp_iv(pca_opt.intermediate_variable());

        // Copy out the intermediate variable values for the spectral domain values in the current bin
        // This requires us to copy into a new array since the bin locations won't be contiguous
        std::vector<blitz::Array<double, 2> > data_list;
        for (int var_idx = 0; var_idx < bin_inp_iv.cols(); var_idx++) {
            Array<double, 2> var_bin_vals(bin_inp_iv.rows(), num_bin_points);
            for(int dom_idx = 0; dom_idx < num_bin_points; dom_idx++) {
                var_bin_bals(Range::all(), dom_idx) = bin_inp_iv(Range::all(), var_idx, bins[var_idx](dom_idx));
            }
            data_list.push_back(data_list);
        }

        // Create bin averaged spectral point
        Array<double, 1> bin_wns(num_bin_points);
        for(int dom_idx = 0; dom_idx < num_bin_points; dom_idx++) {
            bin_wns(dom_idx) = wavenumbers(bins[bin_idx](dom_idx));
        }
        double avg_wn = bin_wns / num_bin_points;

        // Create PCA solution from binned intermediate variables
        PcaEigenSolverGeneric pca_solver(data_list, num_eofs); // 4 = # EOFs, TBD make this configrable
        
        // Create EOF intermediate values for mean and plus and minus each eof
        // We need to translate back from the eigen solvers vectors to arrays for each intermediate variable index
        
        // Compute mean IV values
        Array<double, 2> bin_mean_iv(bin_inp_iv.rows(), bin_inp_iv.cols());
        for (int var_idx = 0; var_idx < bin_inp_iv.cols(); var_idx++) {
            bin_mean_iv(Range::all(), var_idx) = pca_solver.data_mean()[var_idx];
        }
        
        // Compute bin mean intensity from each RT 
        ArrayAd<double, 2> bin_mean_iv_ad(bin_mean_iv)

        Array<double, 1> lidort_mean(lidort_rt->stokes(avg_wn, Spec_index, bin_mean_iv_ad));
        Array<double, 1> twostream_mean(twostream_rt->stokes(avg_wn, Spec_index, bin_mean_iv_ad));
        Array<double, 1> first_order_mean(first_order_rt->stokes(avg_wn, Spec_index, bin_mean_iv_ad));

        // Compute bin mean plus eof and mean minus eof intensity from each RT
        Array<double, 2> lidort_plus(num_eofs, number_stokes());
        Array<double, 2> twostream_plus(num_eofs, number_stokes());
        Array<double, 2> first_order_plus(num_eofs, number_stokes());

        Array<double, 2> lidort_minus(num_eofs, number_stokes());
        Array<double, 2> twostream_minus(num_eofs, number_stokes());
        Array<double, 2> first_order_minus(num_eofs, number_stokes());

        for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
            Array<double, 2> bin_plus_iv(bin_inp_iv.rows(), bin_inp_iv.cols());
            Array<double, 2> bin_minus_iv(bin_inp_iv.rows(), bin_inp_iv.cols());

            for (int var_idx = 0; var_idx < bin_inp_iv.cols(); var_idx++) {
                bin_plus_iv(eof_idx, Range::all(), var_idx) = 
                    pca_solver.data_mean()[var_idx] + pca_solver.eof_properties()[var_idx](Range::all(), eof_idx);

                bin_plus_iv(eof_idx, Range::all(), var_idx) = 
                    pca_solver.data_mean()[var_idx] - pca_solver.eof_properties()[var_idx](Range::all(), eof_idx);
            }

            ArrayAd<double, 2> bin_plus_iv_ad(bin_plus_iv);
            ArrayAd<double, 2> bin_minus_iv_ad(bin_minus_iv);

            lidort_plus(eof_idx, Range::all()) = lidort_rt->stokes(avg_wn, Spec_index, bin_plus_iv_ad);
            lidort_minus(eof_idx, Range::all()) = lidort_rt->stokes(avg_wn, Spec_index, bin_minus_iv_ad);

            twostream_plus(eof_idx, Range::all()) = twostream_rt->stokes(avg_wn, Spec_index, bin_plus_iv_ad);
            twostream_minus(eof_idx, Range::all()) = twostream_rt->stokes(avg_wn, Spec_index, bin_minus_iv_ad);

            first_order_plus(eof_idx, Range::all()) = first_order_rt->stokes(avg_wn, Spec_index, bin_plus_iv_ad);
            first_order_minus(eof_idx, Range::all()) = first_order_rt->stokes(avg_wn, Spec_index, bin_minus_iv_ad);

        }

        // Compute correction factors
        Array<double, 2> bin_corrections(num_bin_points, number_stokes());
        bin_corrections = pca_solver->correction(
                lidort_mean, twostream_mean, first_order_man,
                lidort_plus, twostream_plus, first_order_plus,
                lidort_minus, twostream_minus, first_order_minus);

        // Compute 2stream and first order for all points in the bin and use correction to get final values
        for(int dom_idx = 0; dom_idx < num_bin_points; dom_idx++) {
            int grid_idx = bins[bin_idx](dom_idx);
            int dom_wn = wavenumbers(grid_idx);

            Array<double, 1> twostream_full(twostream_rt->stokes(dom_wn, Spec_index, bin_inp_iv));
            Array<double, 1> first_order_full(first_order_rt->stokes(dom_wn, Spec_index, bin_inp_iv));

            stokes_output(grid_idx, Range::all()) = 
                bin_corrections(dom_idx, Range::all()) * twostream_full(Range::all()) + first_order_full(Range::all());
        }
    }

}

ArrayAd<double, 2> PcaRt::stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const
{
} 
