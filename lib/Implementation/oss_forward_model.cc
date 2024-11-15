#include "oss_forward_model.h"
#include <boost/fusion/include/any.hpp>


using namespace FullPhysics;
using namespace blitz;

OssForwardModel::OssForwardModel(std::vector<boost::shared_ptr<AbsorberVmr>>& Vmr,
        const boost::shared_ptr<Pressure>& Pressure_,
        const boost::shared_ptr<Temperature>& Temperature_,
        const boost::shared_ptr<SurfaceTemperature>& Skin_temperature,
        const boost::shared_ptr<GroundPiecewise>& Ground_,
        DoubleWithUnit Obs_zen_ang, DoubleWithUnit Sol_zen_ang,
        DoubleWithUnit Lat, DoubleWithUnit Surf_alt, bool Lambertian,
        const std::string& Sel_file, const std::string& Od_file, const std::string& Sol_file,
        const std::string& Fix_file, const std::string& Ch_sel_file,
        std::vector<boost::shared_ptr<SpectralDomain>> Channel_domains,
        int Max_chans) :
            vmr(Vmr), pressure(Pressure_), temperature(Temperature_),
            skin_temperature(Skin_temperature), ground(Ground_), obs_zen_ang(Obs_zen_ang),
            sol_zen_ang(Sol_zen_ang), lat(Lat), surf_alt(Surf_alt), lambertian(Lambertian),
            sel_file(Sel_file), sel_file_sz(Sel_file.length()), od_file(Od_file),
            od_file_sz(Od_file.length()), sol_file(Sol_file), sol_file_sz(Sol_file.length()),
            fix_file(Fix_file), fix_file_sz(Fix_file.length()),ch_sel_file(Ch_sel_file),
            ch_sel_file_sz(Ch_sel_file.length()), channel_domains(Channel_domains),
            max_chans(Max_chans), is_setup(false) {
}

void OssForwardModel::setup_grid() {

    std::vector<std::string> gas_names;
    std::vector<std::string> gas_jacobian_names;

    for (int gas_index = 0; gas_index < vmr.size(); gas_index++) {
        gas_names.push_back(vmr[gas_index]->gas_name());
        if (retrieval_flags && retrieval_flags->gas_levels()[gas_index].rows()) {
            gas_jacobian_names.push_back(vmr[gas_index]->gas_name());
        }
    }

    int num_vert_lev = pressure->number_level();
    int num_surf_points = ground->spectral_points().value.rows();

    /* TODO: Clouds disabled */
    float min_extinct_cld = 999.0;

    fixed_inputs = boost::make_shared<OssFixedInputs>(gas_names, gas_jacobian_names, sel_file, od_file, sol_file, fix_file,
            ch_sel_file, num_vert_lev, num_surf_points, min_extinct_cld, channel_domains, max_chans);
    oss_master = boost::make_shared<OssMasters>(fixed_inputs);
    oss_master->init();
    center_spectral_point.units = oss_master->fixed_outputs->full_spectral_point.units;
    center_spectral_point.value.resize(oss_master->fixed_outputs->full_spectral_point.value.rows());
    center_spectral_point.value = cast<double>(oss_master->fixed_outputs->full_spectral_point.value);
    is_setup = true;

}

Spectrum OssForwardModel::radiance(int sensor_index, bool skip_jacobian) const {
    if (!is_setup) {
        throw Exception("Call setup_grid() to initialize before asking for radiances.");
    }
    ArrayAdWithUnit<double, 1> pressure_grid = pressure->pressure_grid().convert(units::mbar);
    Array<float, 1> oss_pressure(pressure_grid.value.rows());
    for (int i = 0; i < oss_pressure.rows(); i++) {
        int level_index = (pressure_grid.rows() - 1) - i;
        oss_pressure(i) = static_cast<float>(pressure_grid.value.value()(level_index));
    }

    ArrayAdWithUnit<double, 1> temperature_grid =
            temperature->temperature_grid(*pressure.get()).convert(units::K);
    Array<float, 1> oss_temperature(temperature_grid.value.rows());
    for (int i = 0; i < oss_temperature.rows(); i++) {
        int level_index = (temperature_grid.rows() - 1) - i;
        oss_temperature(i) = static_cast<float>(temperature_grid.value.value()(level_index));
    }

    float oss_skin_temp = skin_temperature->surface_temperature(sensor_index).convert(units::K).value.value();

    Array<float, 2> vmr_gas(vmr.size(), pressure->number_level());
    for (int gas_index = 0; gas_index < vmr.size(); gas_index++) {
        vmr_gas(gas_index, Range::all()) = cast<float>(vmr[gas_index]->vmr_grid(*pressure.get()).value());
    }

    ArrayWithUnit<double, 1> surface_grid = ground->spectral_points().convert(units::inv_cm);
    Array<float, 1> oss_surface_grid(cast<float>(surface_grid.value(Range::all())));
    Array<float, 1> oss_emiss(surface_grid.value.rows());
    for (int point_index = 0; point_index < surface_grid.value.rows(); point_index++) {
        oss_emiss(point_index) = static_cast<float>(ground->surface_parameter(
                surface_grid.value(point_index), sensor_index).value()(0));
    }
    Array<float, 1> oss_refl(1.0 - oss_emiss);

    /* TODO: Clouds disabled */
    float scale_cld = 0.0;
    float pressure_cld = 0.0;
    int num_cld = 2; // Note in example: "dummy space, following the setup in main_ir.f90"
    Array<float, 1> ext_cld(num_cld);
    ext_cld = 0;
    Array<float, 1> cld_grid(num_cld);
    cld_grid = 0;

    float oss_obs_zen_ang = static_cast<float>(obs_zen_ang.convert(units::deg).value);
    float oss_sol_zen_ang = static_cast<float>(sol_zen_ang.convert(units::deg).value);
    float oss_lat = static_cast<float>(lat.convert(units::deg).value);
    float oss_surf_alt = static_cast<float>(surf_alt.convert(units::m).value);
    boost::shared_ptr<OssModifiedInputs> modified_inputs = boost::make_shared<OssModifiedInputs>(
            oss_pressure, oss_temperature, oss_skin_temp, vmr_gas, oss_emiss, oss_refl,
            scale_cld, pressure_cld, ext_cld, oss_surface_grid,cld_grid,
            oss_obs_zen_ang, oss_sol_zen_ang, oss_lat, oss_surf_alt, lambertian);

    boost::shared_ptr<OssModifiedOutputs> modified_outputs(oss_master->run_fwd_model(sensor_index, modified_inputs));
    cached_outputs = modified_outputs;
    Array<double, 1> rad(cast<double>(modified_outputs->y.value(Range::all())));

    int num_state_variables = 0;
    ArrayAd<double, 1> res(rad.rows(), num_state_variables);
    if (retrieval_flags) {
        num_state_variables += retrieval_flags->num_total_flags();

        res.resize(rad.rows(), num_state_variables);
        /* Fill in jacobians according to retrieval flags */
        int sv_idx = 0;
        for (const int& temp_level : retrieval_flags->temp_levels()) {
            // std::cout << "Retrieving temp level " << temp_level << "\n";
            res.jacobian()(Range::all(), sv_idx) = modified_outputs->xk_temp.value(Range::all(), temp_level);
            sv_idx++;
        }

        Array<bool, 1> skin_temp_sensors(retrieval_flags->skin_temp_sensors());
        for(int skin_temp_sensor_idx = 0; skin_temp_sensor_idx < skin_temp_sensors.rows(); skin_temp_sensor_idx++) {
            if(skin_temp_sensors(skin_temp_sensor_idx)) {
                // If we are not currently modeling the radiance for a given sensor channel
                // then set the jacobians to zero since it still is part of the retrieval
                // vector structure
                if(sensor_index == skin_temp_sensor_idx) {
                    res.jacobian()(Range::all(), sv_idx) = modified_outputs->xk_tskin.value;
                } else {
                    res.jacobian()(Range::all(), sv_idx) = 0.0;
                }
                sv_idx++;
            }
        }

        int gas_jacob_index = 0;
        for (const Array<int, 1>& gas_levels_for_gas : retrieval_flags->gas_levels()) {
            if(gas_levels_for_gas.rows()) {
                for (const int& gas_level : gas_levels_for_gas) {
                    // std::cout << "Retrieving level " << gas_level << " of gas jacobian " << gas_jacob_index << "\n";
                    res.jacobian()(Range::all(), sv_idx) = modified_outputs->xk_out_gas.value(gas_jacob_index,
                            Range::all(), gas_level);
                    sv_idx++;
                }
                gas_jacob_index++;
            }
        }

        for (const int& surf_point : retrieval_flags->emissivity_flags()) {
            // std::cout << "Retrieving emissivity surface point " << surf_point << "\n";
            res.jacobian()(Range::all(), sv_idx) = modified_outputs->xk_em.value(Range::all(), surf_point);
            sv_idx++;
        }
        for (const int& surf_point : retrieval_flags->reflectivity_flags()) {
            // std::cout << "Retrieving reflectivity surface point " << surf_point << "\n";
            res.jacobian()(Range::all(), sv_idx) = modified_outputs->xk_rf.value(Range::all(), surf_point);
            sv_idx++;
        }

        if (sv_idx != num_state_variables) {
            throw Exception("Number of counted SV elements and actually assigned elements did not match");
        }

    }
    res.value() = rad;
    Spectrum convolved_spec(spectral_domain(sensor_index), SpectralRange(res, Unit("W / (m^2 sr cm^{-1})")));
    notify_spectrum_update(convolved_spec, "convolved", sensor_index);
    return convolved_spec;
}

void OssForwardModel::setup_retrieval(const boost::shared_ptr<OssRetrievalFlags>& Retrieval_flags) {
    if (is_setup) {
        throw Exception("Too late to setup retrieval. setup_retrieval() first then setup_grid()");
    }
    retrieval_flags = Retrieval_flags;

}


void OssForwardModel::notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int sensor_index) const
{
  if (olist.size() > 0)
    notify_update_do(boost::make_shared<NamedSpectrum>(updated_spec, spec_name, sensor_index));
}
