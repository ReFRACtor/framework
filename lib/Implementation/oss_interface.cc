#include <numeric>

#include "oss_interface.h"
#include "string_vector_to_char.h"

using namespace FullPhysics;
using namespace blitz;

void OssMasters::init(double spec_domain_tolerance) {
    int num_gases = fixed_inputs->gas_names.size();
    boost::shared_ptr<StringVectorToChar> oss_gases = fixed_inputs->oss_gas_names;
    int len_gas_substr = oss_gases->substrlen;
    std::string oss_gas_str = oss_gases->c_str;

    int num_jacob = fixed_inputs->gas_jacobian_names.size();
    boost::shared_ptr<StringVectorToChar> oss_gas_jacob = fixed_inputs->oss_gas_jacobian_names;
    int len_jacob_substr = oss_gas_jacob->substrlen;
    std::string oss_gas_jacob_str = oss_gas_jacob->c_str;

    int sel_file_sz = fixed_inputs->sel_file.length();
    int od_file_sz  = fixed_inputs->od_file.length();
    int sol_file_sz = fixed_inputs->sol_file.length();
    int fix_file_sz = fixed_inputs->fix_file.length();
    int ch_sel_file_sz = fixed_inputs->ch_sel_file.length();
    int num_vert_lev = fixed_inputs->num_vert_lev;
    int num_surf_points = fixed_inputs->num_surf_points;
    float min_extinct_cld = fixed_inputs->min_extinct_cld.convert(Unit("km^-1")).value;

    ArrayWithUnit<float, 1> full_spectral_point =
            ArrayWithUnit<float, 1>(blitz::Array<float, 1>(fixed_inputs->max_chans),
            Unit("Wavenumbers"));
    int num_chan;

    // Note std::string guaranteed to be contiguous since C++11 per 21.4.1.5 in ISO C++11 std.
    // This allows &(std::string())[0] to give you a non-const char * reference in contrast to c_str()
    cppinitwrapper(num_gases, len_gas_substr, &oss_gas_str[0],
            num_jacob, len_jacob_substr, &oss_gas_jacob_str[0],
            fixed_inputs->sel_file.c_str(), sel_file_sz,
            fixed_inputs->od_file.c_str(), od_file_sz,
            fixed_inputs->sol_file.c_str(), sol_file_sz,
            fixed_inputs->fix_file.c_str(), fix_file_sz,
            fixed_inputs->ch_sel_file.c_str(), ch_sel_file_sz,
            num_vert_lev, num_surf_points,
            min_extinct_cld, num_chan,
            full_spectral_point.value.data(), fixed_inputs->max_chans);

    full_spectral_point.value.resizeAndPreserve(num_chan);

    std::vector<std::vector<int>> spec_pt_idxs;
    if (fixed_inputs->channel_domains.size()) {
        for (int sensor_num = 0; sensor_num < fixed_inputs->channel_domains.size(); sensor_num++) {
            std::vector<int> sensor_spec_pt_idx;
            for (double wavenumber : fixed_inputs->channel_domains[sensor_num]->wavenumber()) {
                // Search for wavenumber in full_spectral_point and add to spec_pt_idxs for sensor
                bool found = false;
                for (int full_idx = 0; full_idx < full_spectral_point.rows(); full_idx++) {

                    if (std::abs(full_spectral_point.value(full_idx) - wavenumber) < spec_domain_tolerance) {
                        found = true;
                        // +1 because fortran interface expects channel idx to begin at 1
                        sensor_spec_pt_idx.push_back(full_idx + 1);
                    }
                }
                if (!found) {
                    std::stringstream select_err;
                    select_err << "Requested wavenumber " << wavenumber << " not found within full list.";
                    throw Exception(select_err.str());
                }
            }
            spec_pt_idxs.push_back(sensor_spec_pt_idx);
            int oss_chan_idx;
            int num_sensor_spec_pt = sensor_spec_pt_idx.size();
            cpploadchanselect(num_sensor_spec_pt, sensor_spec_pt_idx.data(), oss_chan_idx);
            // std::cout << "Channel set loaded for sensor " << sensor_num << " at OSS channel idx " << oss_chan_idx << "\n";
        }
    } else { // Add all spectral points to a single sensor
        std::vector<int> sensor_spec_pt_idx(num_chan);
        // fortran interface says channel idx starts at 1
        std::iota(sensor_spec_pt_idx.begin(), sensor_spec_pt_idx.end(), 1);
        int oss_chan_idx;
        cpploadchanselect(num_chan, sensor_spec_pt_idx.data(), oss_chan_idx);
        // std::cout << "All channels loaded at OSS channel idx " << oss_chan_idx << "\n";
    }

    fixed_outputs = boost::make_shared<OssFixedOutputs>(num_chan, full_spectral_point);

}

boost::shared_ptr<OssModifiedOutputs> OssMasters::run_fwd_model(int Channel_index,
        boost::shared_ptr<OssModifiedInputs> Modified_inputs) {
    int num_vert_lev = fixed_inputs->num_vert_lev;
    int num_surf_points = fixed_inputs->num_surf_points;
    int num_gas = fixed_inputs->gas_names.size();
    float skin_temp = Modified_inputs->skin_temp.convert(units::K).value;
    int num_surf_grid = Modified_inputs->surf_grid.rows();
    float scale_cld = Modified_inputs->scale_cld.convert(units::dimensionless).value;
    float pressure_cld = Modified_inputs->pressure_cld.convert(units::mbar).value;
    int num_cld = Modified_inputs->cld_grid.rows();
    float obs_zen_ang = Modified_inputs->obs_zen_ang.convert(units::deg).value;
    float sol_zen_ang = Modified_inputs->sol_zen_ang.convert(units::deg).value;
    float lat = Modified_inputs->lat.convert(units::deg).value;
    float surf_alt = Modified_inputs->surf_alt.convert(units::m).value;
    int is_lambertian = Modified_inputs->lambertian;
    int num_gas_jacob = fixed_inputs->gas_jacobian_names.size();

    Channel_index++; // OSS channel indexes start at 1
    cppsetchanselect(Channel_index);
    int num_chan;
    cppgetcountchannel(Channel_index, num_chan);

    /* Outputs */
    blitz::Array<float, 1> y = blitz::Array<float, 1>(num_chan);
    blitz::Array<float, 2> xk_temp = blitz::Array<float, 2>(num_chan, num_vert_lev);
    blitz::Array<float, 1> xk_tskin = blitz::Array<float, 1>(num_chan);
    blitz::Array<float, 3> xk_out_gas = blitz::Array<float, 3>(num_gas_jacob, num_chan, num_vert_lev);
    blitz::Array<float, 2> xk_em = blitz::Array<float, 2>(num_chan, num_surf_points);
    blitz::Array<float, 2> xk_rf = blitz::Array<float, 2>(num_chan, num_surf_points);
    blitz::Array<float, 1> xk_cldln_pres = blitz::Array<float, 1>(num_chan);
    blitz::Array<float, 2> xk_cldln_ext = blitz::Array<float, 2>(num_chan, num_cld);

    /* Converted inputs */
    blitz::Array<float, 1> pressure = Modified_inputs->pressure.convert(units::mbar).value;
    blitz::Array<float, 1> temp = Modified_inputs->temp.convert(units::K).value;
    blitz::Array<float, 2> vmr_gas = Modified_inputs->vmr_gas.convert(units::dimensionless).value;
    blitz::Array<float, 1> emis = Modified_inputs->emis.convert(units::dimensionless).value;
    blitz::Array<float, 1> refl = Modified_inputs->refl.convert(units::dimensionless).value;
    blitz::Array<float, 1> ext_cld = Modified_inputs->ext_cld.convert(Unit("km^-1")).value;
    blitz::Array<float, 1> surf_grid = Modified_inputs->surf_grid.convert(units::inv_cm).value;
    blitz::Array<float, 1> cld_grid = Modified_inputs->cld_grid.convert(units::inv_cm).value;

    cppfwdwrapper(num_vert_lev, num_gas, pressure.data(),
            temp.data(), skin_temp, vmr_gas.data(), num_surf_grid, emis.data(), refl.data(),
            scale_cld, pressure_cld, num_cld, ext_cld.data(), surf_grid.data(), cld_grid.data(),
            obs_zen_ang, sol_zen_ang, lat, surf_alt, is_lambertian, num_gas_jacob, num_chan, y.data(),
            xk_temp.data(), xk_tskin.data(), xk_out_gas.data(), xk_em.data(), xk_rf.data(),
            xk_cldln_pres.data(), xk_cldln_ext.data());

    return boost::make_shared<OssModifiedOutputs>(y, xk_temp, xk_tskin, xk_out_gas, xk_em, xk_rf, xk_cldln_pres, xk_cldln_ext);
}

