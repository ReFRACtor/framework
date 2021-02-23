#include <fstream>
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

#include "pca_rt.h"

#include "optical_properties_wrt_rt.h"

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
             bool do_thermal_emission,
             bool do_3M_correction) 
: RadiativeTransferFixedStokesCoefficient(Stokes_coef), atm(Atm), primary_absorber(Primary_absorber), 
  bin_method(Bin_method), num_bins(Num_bins), num_eofs(Num_eofs), do_3m_correction(do_3M_correction)
{

    // Need to convert aerosol from atmosphere to AerosolOptical for PCA optical properties packing
    aerosol_optical = boost::dynamic_pointer_cast<AerosolOptical>(atm->aerosol_ptr());

    // Only throw an error if the original pointer is not already null
    if(atm->aerosol_ptr() && !aerosol_optical) {
        throw Exception("Failed to convert aerosol class to AerosolOptical");
    }

    lidort_rt.reset(new LidortRt(Atm, Stokes_coef, Sza, Zen, Azm, false, Number_streams, Number_moments, true, do_solar_sources, do_thermal_emission, do_thermal_emission));

    twostream_rt.reset(new TwostreamRt(Atm, Stokes_coef, Sza, Zen, Azm, true, do_solar_sources, do_thermal_emission));

    first_order_rt.reset(new FirstOrderRt(Atm, Stokes_coef, Sza, Zen, Azm, Number_streams, Number_moments, do_solar_sources, do_thermal_emission));

    num_gas = atm->absorber_ptr()->number_species();
    num_layer = atm->pressure_ptr()->number_layer();

    if (aerosol_optical) {
        num_aerosol = aerosol_optical->number_particle();
    } else {
        num_aerosol = 0;
    }

    num_packed_var = 2 + num_aerosol;
}

void PCARt::clear_pca_objects() const
{
    pca_opt.clear();
    pca_bin.reset();
    pca_solvers.clear();
}

void PCARt::compute_bins(const SpectralDomain& Spec_domain, int Spec_index) const
{
    
    boost::optional<std::string> progress_message("Optical Properties: Channel " + boost::lexical_cast<std::string>(Spec_index + 1));
    boost::shared_ptr<boost::progress_display> progress_bar = progress_display(Spec_domain.data(), progress_message);

    // Compute optical properties for the whole band
    for (int dom_idx = 0; dom_idx < Spec_domain.data().rows(); dom_idx++) {
        boost::shared_ptr<OpticalPropertiesWrtRt> point_opt_prop(new OpticalPropertiesWrtRt());
        point_opt_prop->initialize(DoubleWithUnit(Spec_domain.data()(dom_idx), Spec_domain.units()), Spec_index,
                                   atm->absorber_ptr(), atm->rayleigh_ptr(), atm->aerosol_ptr());
        pca_opt.push_back(point_opt_prop);

        if(progress_bar) *progress_bar += 1;
    }

    int primary_abs_index = atm->absorber_ptr()->gas_index(primary_absorber);

    // Compute bins
    pca_bin.reset(new PCABinning(pca_opt, bin_method, num_bins, primary_abs_index));
}

const boost::shared_ptr<PCAEigenSolver> PCARt::compute_bin_solution(const Array<int, 1>& data_indexes) const
{
    // Solve PCA solution for intermediate variable optical properties for the bin

    std::vector<blitz::Array<double, 2> > data_list;

    int num_point = data_indexes.rows();
    
    // Allocate an array for each packed variable
    for (int var_idx = 0; var_idx < num_packed_var; var_idx++) {
        blitz::Array<double, 2> var_packed(num_layer, num_point);
        data_list.push_back(var_packed);
    }

    // Copy out the intermediate variable values for the spectral domain values in the current bin
    // This requires us to copy into a new array since the bin locations won't be contiguous
    for(int dom_idx = 0; dom_idx < num_point; dom_idx++) {
        ArrayAd<double, 2> point_packed = OpticalPropertiesPca::pack(pca_opt[data_indexes(dom_idx)]);

        for (int var_idx = 0; var_idx < num_packed_var; var_idx++) {
            data_list[var_idx](Range::all(), dom_idx) = point_packed.value()(Range::all(), var_idx);
        }
    }

    // Create PCA solution from binned intermediate variables
    boost::shared_ptr<PCAEigenSolver> bin_solver;
    bin_solver.reset(new PCAEigenSolverGeneric(data_list, num_eofs));
    
    // Save PCA solver to class attribute to be accessible for debugging
    pca_solvers.push_back(bin_solver);

    return bin_solver;
}

const boost::shared_ptr<PCABinOpticalProperties> PCARt::compute_bin_optical_props(boost::shared_ptr<PCAEigenSolver> pca_solver, double bin_wn) const
{
    // Create EOF intermediate values for mean and plus and minus each eof
    //
    // We need to translate back from the eigen solvers vectors to arrays for each intermediate variable index
 
    boost::shared_ptr<PCABinOpticalProperties> bin_opt_props(new PCABinOpticalProperties());

    // Compute mean IV values
    Array<double, 2> bin_mean_iv(num_layer, num_packed_var);
    for (int var_idx = 0; var_idx < num_packed_var; var_idx++) {
        bin_mean_iv(Range::all(), var_idx) = exp(pca_solver->data_mean()[var_idx]);
    }
    ArrayAd<double, 2> bin_mean_iv_ad(bin_mean_iv);
    
    // Compute bin mean intensity from each RT 
    bin_opt_props->mean.reset(new OpticalPropertiesPca(bin_mean_iv_ad, bin_wn, aerosol_optical, num_gas, num_aerosol));

    for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
        Array<double, 2> bin_plus_iv(num_layer, num_packed_var);
        Array<double, 2> bin_minus_iv(num_layer, num_packed_var);

        for (int var_idx = 0; var_idx < num_packed_var; var_idx++) {
            bin_plus_iv(Range::all(), var_idx) = 
                exp(pca_solver->data_mean()[var_idx] + pca_solver->eof_properties()[var_idx](Range::all(), eof_idx));

            bin_minus_iv(Range::all(), var_idx) = 
                exp(pca_solver->data_mean()[var_idx] - pca_solver->eof_properties()[var_idx](Range::all(), eof_idx));
        }

        ArrayAd<double, 2> bin_plus_iv_ad(bin_plus_iv);
        ArrayAd<double, 2> bin_minus_iv_ad(bin_minus_iv);

        boost::shared_ptr<OpticalPropertiesPca> plus_opt_props(new OpticalPropertiesPca(bin_plus_iv_ad, bin_wn, aerosol_optical, num_gas, num_aerosol));
        boost::shared_ptr<OpticalPropertiesPca> minus_opt_props(new OpticalPropertiesPca(bin_minus_iv_ad, bin_wn, aerosol_optical, num_gas, num_aerosol));
    
        bin_opt_props->eof_plus.push_back(plus_opt_props);
        bin_opt_props->eof_minus.push_back(minus_opt_props);
    }

    return bin_opt_props;
}

const double PCARt::bin_effective_wavenumber(Array<double, 1> &win_wavenumbers, int bin_index) const
{
    // Create bin averaged spectral point
    int num_bin_points = pca_bin->num_bin_points()(bin_index);
    Array<int, 1> bin_indexes(pca_bin->bin_indexes()[bin_index]);

    double avg_wn = 0;
    for(int dom_idx = 0; dom_idx < num_bin_points; dom_idx++) {
        avg_wn += win_wavenumbers(bin_indexes(dom_idx));
    }
    avg_wn /= num_bin_points;

    return avg_wn;
}

Array<double, 2> PCARt::compute_bin_correction_factors(boost::shared_ptr<PCAEigenSolver>& pca_solver, boost::shared_ptr<PCABinOpticalProperties>& bin_opt_props, double bin_wn, int channel_index) const
{
    Array<double, 1> lidort_mean(lidort_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->mean));
    Array<double, 1> twostream_mean(twostream_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->mean));

    // Compute bin mean plus eof and mean minus eof intensity from each RT
    Array<double, 2> lidort_plus(num_eofs, number_stokes());
    Array<double, 2> twostream_plus(num_eofs, number_stokes());

    Array<double, 2> lidort_minus(num_eofs, number_stokes());
    Array<double, 2> twostream_minus(num_eofs, number_stokes());

    Array<double, 1> first_order_mean;
    Array<double, 2> first_order_plus;
    Array<double, 2> first_order_minus;

    if (do_3m_correction) { 
        first_order_mean.reference(first_order_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->mean));
        first_order_plus.resize(num_eofs, number_stokes());
        first_order_minus.resize(num_eofs, number_stokes());
    }

    bool fo_captured = false;
    for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
        lidort_plus(eof_idx, Range::all()) = lidort_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->eof_plus[eof_idx]);
        lidort_minus(eof_idx, Range::all()) = lidort_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->eof_minus[eof_idx]);

        twostream_plus(eof_idx, Range::all()) = twostream_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->eof_plus[eof_idx]);
        twostream_minus(eof_idx, Range::all()) = twostream_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->eof_minus[eof_idx]);

        first_order_plus(eof_idx, Range::all()) = first_order_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->eof_plus[eof_idx]);
        first_order_minus(eof_idx, Range::all()) = first_order_rt->stokes_single_wn(bin_wn, channel_index, bin_opt_props->eof_minus[eof_idx]);
    }

    // Compute correction factors
    Array<double, 2> bin_corrections;
    if (do_3m_correction) {
        bin_corrections.reference(
            pca_solver->correction_3m(
                lidort_mean, twostream_mean, first_order_mean,
                lidort_plus, twostream_plus, first_order_plus,
                lidort_minus, twostream_minus, first_order_minus)
            );
    } else {
        bin_corrections.reference(
            pca_solver->correction_2m(
                lidort_mean, twostream_mean,
                lidort_plus, twostream_plus,
                lidort_minus, twostream_minus)
            );
     }

    return bin_corrections;
}

blitz::Array<double, 2> PCARt::stokes(const SpectralDomain& Spec_domain, int Spec_index) const
{
    FunctionTimer ft(timer.function_timer(true));

    // Clear out debugging and per call objects
    clear_pca_objects();

    // Our final output from this routine
    Array<double, 2> stokes_output(Spec_domain.data().rows(), number_stokes());

    // Convenience since other RTs take wavenumbers as input
    Array<double, 1> wavenumbers(Spec_domain.wavenumber());

    // Create optical property bins from full range optical properties
    compute_bins(Spec_domain, Spec_index);

    std::vector<blitz::Array<int, 1> > bins = pca_bin->bin_indexes();

    boost::optional<std::string> progress_message("PCA RT: Channel " + boost::lexical_cast<std::string>(Spec_index + 1));
    boost::shared_ptr<boost::progress_display> progress_bar = progress_display(Spec_domain.data(), progress_message);

    // Compute values for each bin
    for (int bin_idx = 0; bin_idx < bins.size(); bin_idx++) {
        int num_bin_points = pca_bin->num_bin_points()(bin_idx);

        // Skip empty bins
        if (num_bin_points == 0) {
            continue;
        }

        // We need an effective wavenumber for aerosol property computations
        double bin_wn = bin_effective_wavenumber(wavenumbers, bin_idx);

        // Solve PCA solution from optical properties for the bin
        boost::shared_ptr<PCAEigenSolver> pca_solver = compute_bin_solution(bins[bin_idx]);
        
        // Compute the binned optical properties
        boost::shared_ptr<PCABinOpticalProperties> bin_opt_props = compute_bin_optical_props(pca_solver, bin_wn);

        // Compute per wavenumber correction factors from RT runs of binned optical properties
        Array<double, 2> bin_corrections = compute_bin_correction_factors(pca_solver, bin_opt_props, bin_wn, Spec_index);
        
        // Compute 2stream and first order for all points in the bin and use correction to get final values
        for(int dom_idx = 0; dom_idx < num_bin_points; dom_idx++) {
            int grid_idx = bins[bin_idx](dom_idx);
            int dom_wn = wavenumbers(grid_idx);

            Array<double, 1> twostream_full(twostream_rt->stokes_single_wn(dom_wn, Spec_index, pca_opt[grid_idx]));
            Array<double, 1> first_order_full(first_order_rt->stokes_single_wn(dom_wn, Spec_index, pca_opt[grid_idx]));

            // The difference here is the 2M correction only applies to the 2stream values
            if (do_3m_correction) {
                stokes_output(grid_idx, Range::all()) = 
                    bin_corrections(dom_idx, Range::all()) * (twostream_full(Range::all()) + first_order_full(Range::all()));
            } else {
                stokes_output(grid_idx, Range::all()) = 
                    bin_corrections(dom_idx, Range::all()) * twostream_full(Range::all()) + first_order_full(Range::all());
            }

            if (progress_bar) *progress_bar += 1;
        }
    }

    Logger::info() << atm->timer_info();

    return stokes_output;
}

ArrayAd<double, 2> PCARt::stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const
{
    FunctionTimer ft(timer.function_timer(true));

    // Clear out debugging and per call objects
    clear_pca_objects();

    // Our final output from this routine
    ArrayAd<double, 2> stokes_output(Spec_domain.data().rows(), number_stokes(), 0);

    // Convenience since other RTs take wavenumbers as input
    Array<double, 1> wavenumbers(Spec_domain.wavenumber());

    // Create optical property bins from full range optical properties
    compute_bins(Spec_domain, Spec_index);

    std::vector<blitz::Array<int, 1> > bins = pca_bin->bin_indexes();

    boost::optional<std::string> progress_message("PCA RT + Jacobian: Channel " + boost::lexical_cast<std::string>(Spec_index + 1));
    boost::shared_ptr<boost::progress_display> progress_bar = progress_display(Spec_domain.data(), progress_message);

    // Compute values for each bin
    for (int bin_idx = 0; bin_idx < bins.size(); bin_idx++) {
        int num_bin_points = pca_bin->num_bin_points()(bin_idx);

        // Skip empty bins
        if (num_bin_points == 0) {
            continue;
        }

        // We need an effective wavenumber for aerosol property computations
        double bin_wn = bin_effective_wavenumber(wavenumbers, bin_idx);

        // Solve PCA solution from optical properties for the bin
        boost::shared_ptr<PCAEigenSolver> pca_solver = compute_bin_solution(bins[bin_idx]);
        
        // Compute the binned optical properties
        boost::shared_ptr<PCABinOpticalProperties> bin_opt_props = compute_bin_optical_props(pca_solver, bin_wn);

        // Compute per wavenumber correction factors from RT runs of binned optical properties
        Array<double, 2> bin_corrections = compute_bin_correction_factors(pca_solver, bin_opt_props, bin_wn, Spec_index);
        
        // Compute 2stream and first order for all points in the bin and use correction to get final values
        for(int dom_idx = 0; dom_idx < num_bin_points; dom_idx++) {
            int grid_idx = bins[bin_idx](dom_idx);
            int dom_wn = wavenumbers(grid_idx);

            ArrayAd<double, 1> twostream_full(twostream_rt->stokes_and_jacobian_single_wn(dom_wn, Spec_index, pca_opt[grid_idx]));
            ArrayAd<double, 1> first_order_full(first_order_rt->stokes_and_jacobian_single_wn(dom_wn, Spec_index, pca_opt[grid_idx]));

            // The difference here is the 2M correction only applies to the 2stream values
            if (do_3m_correction) {
                stokes_output.value()(grid_idx, Range::all()) = 
                    bin_corrections(dom_idx, Range::all()) * (twostream_full.value()(Range::all()) + first_order_full.value()(Range::all()));
            } else {
                stokes_output.value()(grid_idx, Range::all()) = 
                    bin_corrections(dom_idx, Range::all()) * twostream_full.value()(Range::all()) + first_order_full.value()(Range::all());
            }
            
            if (stokes_output.number_variable() == 0) {
                stokes_output.resize_number_variable(twostream_full.number_variable());
            }

            // Corrections can be applied for each jacobian independently according to Vijay
            for (int jac_idx = 0; jac_idx < twostream_full.number_variable(); jac_idx++) {
                if (do_3m_correction) {
                    stokes_output.jacobian()(grid_idx, Range::all(), jac_idx) = 
                        bin_corrections(dom_idx, Range::all()) * (twostream_full.jacobian()(Range::all(), jac_idx) + first_order_full.jacobian()(Range::all(), jac_idx));
                } else {
                    stokes_output.jacobian()(grid_idx, Range::all(), jac_idx) = 
                        bin_corrections(dom_idx, Range::all()) * twostream_full.jacobian()(Range::all(), jac_idx) + first_order_full.jacobian()(Range::all(), jac_idx);
                }
            }

            if (progress_bar) *progress_bar += 1;
        }
    }

    Logger::info() << atm->timer_info();

    return stokes_output;
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
