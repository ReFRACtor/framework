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
            max_chans(Max_chans) {
    std::vector<std::string> gas_names;
    std::vector<std::string> gas_jacobian_names;

    for (int gas_index = 0; gas_index < vmr.size(); gas_index++) {
        gas_names.push_back(vmr[gas_index]->gas_name());
        gas_jacobian_names.push_back(vmr[gas_index]->gas_name());

        // Disabled because can not call ->state_used before items are attached to state vector
        // which wil be after this consturctor has been called
        //if (any(vmr[gas_index]->state_used())) {
        //    gas_jacobian_names.push_back(vmr[gas_index]->gas_name());
        //}
    }

    int num_vert_lev = pressure->number_level();
    int num_surf_points = ground->spectral_points().value.rows();

    /* TODO: Clouds disabled */
    float min_extinct_cld = 999.0;

    fixed_inputs = boost::make_shared<OssFixedInputs>(gas_names, gas_jacobian_names, sel_file, od_file, sol_file, fix_file,
            ch_sel_file, num_vert_lev, num_surf_points, min_extinct_cld, channel_domains, max_chans);
    oss_master = boost::make_shared<OssMasters>(fixed_inputs);
}

void OssForwardModel::setup_grid() {
    oss_master->init();
    center_spectral_point.units = oss_master->fixed_outputs->center_spectral_point.units;
    center_spectral_point.value.resize(oss_master->fixed_outputs->center_spectral_point.value.rows());
    center_spectral_point.value = cast<double>(oss_master->fixed_outputs->center_spectral_point.value);
}

Spectrum OssForwardModel::radiance(int channel_index, bool skip_jacobian) const {
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

    float oss_skin_temp = skin_temperature->surface_temperature(channel_index).convert(units::K).value.value();

    Array<float, 2> vmr_gas(vmr.size(), pressure->number_level());
    for (int gas_index = 0; gas_index < vmr.size(); gas_index++) {
        vmr_gas(gas_index, Range::all()) = cast<float>(vmr[gas_index]->vmr_grid(*pressure.get()).value());
    }

    ArrayWithUnit<double, 1> surface_grid = ground->spectral_points().convert(units::inv_cm);
    Array<float, 1> oss_surface_grid(cast<float>(surface_grid.value(Range::all())));
    Array<float, 1> oss_emiss(surface_grid.value.rows());
    for (int point_index = 0; point_index < surface_grid.value.rows(); point_index++) {
        oss_emiss(point_index) = static_cast<float>(ground->surface_parameter(
                surface_grid.value(point_index), channel_index).value()(0));
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

    boost::shared_ptr<OssModifiedOutputs> modified_outputs(oss_master->run_fwd_model(channel_index, modified_inputs));
    cached_outputs = modified_outputs;
    Array<double, 1> rad(cast<double>(modified_outputs->y.value(Range::all())));

    /* TODO: Add jacobian to SpectralRange */
    return Spectrum(spectral_domain(channel_index), SpectralRange(rad, Unit("W / cm^2 / sr / cm^-1")));
}

void OssForwardModel::setup_retrieval(OssRetrievalFlags Retrieval_flags) const {
    retrieval_flags = boost::make_shared<OssRetrievalFlags>(Retrieval_flags);
}
