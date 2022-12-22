/* This file was auto-generated */

%include "fp_common.i"

%{
#include "spurr_interface_types.h"
%}

%fp_shared_ptr(FullPhysics::Spurr_Type_Structure);
%fp_shared_ptr(FullPhysics::Spurr_Pars_Base);

%fp_shared_ptr(FullPhysics::Spurr_Brdf_Lin_Sup_Inputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Lin_Sup_Outputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Sup_Inputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Sup_Outputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Input_Exception_Handling_Base);
%fp_shared_ptr(FullPhysics::Spurr_Brdf_Output_Exception_Handling_Base);
%fp_shared_ptr(FullPhysics::Sleave_Sup_Inputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Lincontrol_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Linoptical_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Lininputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Lincontrol_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Lininputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Linatmos_Base);
%fp_shared_ptr(FullPhysics::Spurr_Linsurf_Base);
%fp_shared_ptr(FullPhysics::Spurr_Linoutputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lin_Sup_Brdf_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lin_Sup_Sleave_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lin_Sup_Ss_Surf_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lin_Sup_Ss_Base);
%fp_shared_ptr(FullPhysics::Spurr_Lin_Sup_Inout_Base);
%fp_shared_ptr(FullPhysics::Spurr_Main_Outputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Wladjusted_Outputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Exception_Handling_Base);
%fp_shared_ptr(FullPhysics::Spurr_Input_Exception_Handling_Base);
%fp_shared_ptr(FullPhysics::Spurr_Outputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Sup_Brdf_Base);
%fp_shared_ptr(FullPhysics::Spurr_Sup_Sleave_Base);
%fp_shared_ptr(FullPhysics::Spurr_Sup_Ss_Base);
%fp_shared_ptr(FullPhysics::Spurr_Sup_Inout_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Boolean_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Control_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Sunrays_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Uservalues_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Chapman_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Optical_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Write_Base);
%fp_shared_ptr(FullPhysics::Spurr_Fixed_Inputs_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Boolean_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Control_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Sunrays_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Uservalues_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Chapman_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Optical_Base);
%fp_shared_ptr(FullPhysics::Spurr_Modified_Inputs_Base);
 

namespace FullPhysics {

%nodefaultctor Spurr_Type_Structure;

class Spurr_Type_Structure {
public:
    std::string print_to_string() const;
};

%nodefaultctor Spurr_Pars_Base;

class Spurr_Pars_Base {
public:
  virtual const double& bigexp() const = 0;
  virtual const int& bpdfndvi_idx() const = 0;
  virtual const int& bpdfsoil_idx() const = 0;
  virtual const int& bpdfvegn_idx() const = 0;
  virtual const int& coxmunk_idx() const = 0;
  virtual const double& deg_to_rad() const = 0;
  virtual const int& dnidx() const = 0;
  virtual const double& eps3() const = 0;
  virtual const double& eps4() const = 0;
  virtual const double& eps5() const = 0;
  virtual const double& four() const = 0;
  virtual const double& half() const = 0;
  virtual const int& hapke_idx() const = 0;
  virtual const double& hopital_tolerance() const = 0;
  virtual const int& lambertian_idx() const = 0;
  virtual const int& lidense_idx() const = 0;
  virtual const int& lisparse_idx() const = 0;
  virtual const int& max_allstrms() const = 0;
  virtual const int& max_allstrms_p1() const = 0;
  virtual const int& max_atmoswfs() const = 0;
  virtual const int& max_brdf_kernels() const = 0;
  virtual const int& max_brdf_parameters() const = 0;
  virtual const int& max_directions() const = 0;
  virtual const int& max_geometries() const = 0;
  virtual const int& max_messages() const = 0;
  virtual const int& max_msrs_muquad() const = 0;
  virtual const int& max_msrs_phiquad() const = 0;
  virtual const int& max_partlayers() const = 0;
  virtual const int& max_sleavewfs() const = 0;
  virtual const int& max_surfacewfs() const = 0;
  virtual const double& max_tau_qpath() const = 0;
  virtual const double& max_tau_spath() const = 0;
  virtual const double& max_tau_upath() const = 0;
  virtual const int& max_taylor_terms() const = 0;
  virtual const int& max_thermal_coeffs() const = 0;
  virtual const int& max_user_levels() const = 0;
  virtual const int& max_user_obsgeoms() const = 0;
  virtual const int& max_user_relazms() const = 0;
  virtual const int& max_user_streams() const = 0;
  virtual const int& maxbandtotal() const = 0;
  virtual const int& maxbeams() const = 0;
  virtual const int& maxbrdf_idx() const = 0;
  virtual const int& maxfinelayers() const = 0;
  virtual const int& maxfourier() const = 0;
  virtual const int& maxlayers() const = 0;
  virtual const int& maxmoments() const = 0;
  virtual const int& maxmoments_input() const = 0;
  virtual const int& maxsthalf_brdf() const = 0;
  virtual const int& maxstreams() const = 0;
  virtual const int& maxstreams_2() const = 0;
  virtual const int& maxstreams_brdf() const = 0;
  virtual const int& maxstreams_p1() const = 0;
  virtual const int& maxstreams_scaling() const = 0;
  virtual const int& maxtotal() const = 0;
  virtual const double& minus_one() const = 0;
  virtual const double& minus_two() const = 0;
  virtual const int& modfresnel_idx() const = 0;
  virtual const int& newcmglint_idx() const = 0;
  virtual const double& omega_smallnum() const = 0;
  virtual const double& one() const = 0;
  virtual const double& onep5() const = 0;
  virtual const double& pi2() const = 0;
  virtual const double& pi4() const = 0;
  virtual const double& pie() const = 0;
  virtual const double& pio2() const = 0;
  virtual const double& pio4() const = 0;
  virtual const double& quarter() const = 0;
  virtual const int& rahman_idx() const = 0;
  virtual const int& rossthick_idx() const = 0;
  virtual const int& rossthin_idx() const = 0;
  virtual const int& roujean_idx() const = 0;
  virtual const int& rtkhotspot_idx() const = 0;
  virtual const double& smallnum() const = 0;
  virtual const int& snowbrdf_idx() const = 0;
  virtual const double& taylor_small() const = 0;
  virtual const double& three() const = 0;
  virtual const double& two() const = 0;
  virtual const int& upidx() const = 0;
  virtual const double& zero() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;

};


%nodefaultctor Spurr_Brdf_Lin_Sup_Inputs_Base;

class Spurr_Brdf_Lin_Sup_Inputs_Base : public Spurr_Type_Structure {
public:

  virtual const bool bs_do_bsavalue_wf() const = 0;
  virtual void bs_do_bsavalue_wf(const bool& bs_do_bsavalue_wf_in) = 0;
  
  virtual const blitz::Array<bool, 1> bs_do_kernel_factor_wfs() const = 0;
  virtual void bs_do_kernel_factor_wfs(const blitz::Array<bool, 1>& bs_do_kernel_factor_wfs_in) = 0;
  
  virtual const blitz::Array<bool, 2> bs_do_kernel_params_wfs() const = 0;
  virtual void bs_do_kernel_params_wfs(const blitz::Array<bool, 2>& bs_do_kernel_params_wfs_in) = 0;
  
  virtual const blitz::Array<bool, 1> bs_do_kparams_derivs() const = 0;
  virtual void bs_do_kparams_derivs(const blitz::Array<bool, 1>& bs_do_kparams_derivs_in) = 0;
  
  virtual const bool bs_do_windspeed_wf() const = 0;
  virtual void bs_do_windspeed_wf(const bool& bs_do_windspeed_wf_in) = 0;
  
  virtual const bool bs_do_wsavalue_wf() const = 0;
  virtual void bs_do_wsavalue_wf(const bool& bs_do_wsavalue_wf_in) = 0;
  
  virtual const int& bs_n_kernel_factor_wfs() const = 0;
  virtual void bs_n_kernel_factor_wfs(const int& bs_n_kernel_factor_wfs_in) = 0;
  
  virtual const int& bs_n_kernel_params_wfs() const = 0;
  virtual void bs_n_kernel_params_wfs(const int& bs_n_kernel_params_wfs_in) = 0;
  
  virtual const int& bs_n_surface_wfs() const = 0;
  virtual void bs_n_surface_wfs(const int& bs_n_surface_wfs_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Brdf_Lin_Sup_Outputs_Base;

class Spurr_Brdf_Lin_Sup_Outputs_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * bs_ls_brdf_f
   * bs_ls_brdf_f_0
   * bs_ls_dbounce_brdfunc
   * bs_ls_emissivity
   * bs_ls_user_brdf_f
   * bs_ls_user_brdf_f_0
   * bs_ls_user_emissivity
   
  */
};

%nodefaultctor Spurr_Brdf_Sup_Inputs_Base;

class Spurr_Brdf_Sup_Inputs_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<double, 1>& bs_beam_szas() const = 0;
  virtual void bs_beam_szas(const blitz::Array<double, 1>& bs_beam_szas_in) = 0;
  
  virtual const blitz::Array<double, 1>& bs_brdf_factors() const = 0;
  virtual void bs_brdf_factors(const blitz::Array<double, 1>& bs_brdf_factors_in) = 0;
  
  virtual const std::vector< std::string > bs_brdf_names() const = 0;
  
  virtual const blitz::Array<double, 2>& bs_brdf_parameters() const = 0;
  virtual void bs_brdf_parameters(const blitz::Array<double, 2>& bs_brdf_parameters_in) = 0;
  
  virtual const double& bs_bsa_value() const = 0;
  virtual void bs_bsa_value(const double& bs_bsa_value_in) = 0;
  
  virtual const bool bs_do_brdf_surface() const = 0;
  virtual void bs_do_brdf_surface(const bool& bs_do_brdf_surface_in) = 0;
  
  virtual const bool bs_do_bsa_scaling() const = 0;
  virtual void bs_do_bsa_scaling(const bool& bs_do_bsa_scaling_in) = 0;
  
  virtual const bool bs_do_directbounce_only() const = 0;
  virtual void bs_do_directbounce_only(const bool& bs_do_directbounce_only_in) = 0;
  
  virtual const bool bs_do_doublet_geometry() const = 0;
  virtual void bs_do_doublet_geometry(const bool& bs_do_doublet_geometry_in) = 0;
  
  virtual const bool bs_do_facetisotropy() const = 0;
  virtual void bs_do_facetisotropy(const bool& bs_do_facetisotropy_in) = 0;
  
  virtual const bool bs_do_foamoption() const = 0;
  virtual void bs_do_foamoption(const bool& bs_do_foamoption_in) = 0;
  
  virtual const bool bs_do_glintshadow() const = 0;
  virtual void bs_do_glintshadow(const bool& bs_do_glintshadow_in) = 0;
  
  virtual const bool bs_do_glitter_msrcorr() const = 0;
  virtual void bs_do_glitter_msrcorr(const bool& bs_do_glitter_msrcorr_in) = 0;
  
  virtual const bool bs_do_glitter_msrcorr_dbonly() const = 0;
  virtual void bs_do_glitter_msrcorr_dbonly(const bool& bs_do_glitter_msrcorr_dbonly_in) = 0;
  
  virtual const bool bs_do_newcmglint() const = 0;
  virtual void bs_do_newcmglint(const bool& bs_do_newcmglint_in) = 0;
  
  virtual const bool bs_do_shadow_effect() const = 0;
  virtual void bs_do_shadow_effect(const bool& bs_do_shadow_effect_in) = 0;
  
  virtual const bool bs_do_solar_sources() const = 0;
  virtual void bs_do_solar_sources(const bool& bs_do_solar_sources_in) = 0;
  
  virtual const bool bs_do_surface_emission() const = 0;
  virtual void bs_do_surface_emission(const bool& bs_do_surface_emission_in) = 0;
  
  virtual const bool bs_do_user_obsgeoms() const = 0;
  virtual void bs_do_user_obsgeoms(const bool& bs_do_user_obsgeoms_in) = 0;
  
  virtual const bool bs_do_user_streams() const = 0;
  virtual void bs_do_user_streams(const bool& bs_do_user_streams_in) = 0;
  
  virtual const bool bs_do_wsa_scaling() const = 0;
  virtual void bs_do_wsa_scaling(const bool& bs_do_wsa_scaling_in) = 0;
  
  virtual const bool bs_do_wsabsa_output() const = 0;
  virtual void bs_do_wsabsa_output(const bool& bs_do_wsabsa_output_in) = 0;
  
  virtual const int& bs_glitter_msrcorr_nmuquad() const = 0;
  virtual void bs_glitter_msrcorr_nmuquad(const int& bs_glitter_msrcorr_nmuquad_in) = 0;
  
  virtual const int& bs_glitter_msrcorr_nphiquad() const = 0;
  virtual void bs_glitter_msrcorr_nphiquad(const int& bs_glitter_msrcorr_nphiquad_in) = 0;
  
  virtual const int& bs_glitter_msrcorr_order() const = 0;
  virtual void bs_glitter_msrcorr_order(const int& bs_glitter_msrcorr_order_in) = 0;
  
  virtual const blitz::Array<bool, 1> bs_lambertian_kernel_flag() const = 0;
  virtual void bs_lambertian_kernel_flag(const blitz::Array<bool, 1>& bs_lambertian_kernel_flag_in) = 0;
  
  virtual const int& bs_n_brdf_kernels() const = 0;
  virtual void bs_n_brdf_kernels(const int& bs_n_brdf_kernels_in) = 0;
  
  virtual const blitz::Array<int, 1>& bs_n_brdf_parameters() const = 0;
  virtual void bs_n_brdf_parameters(const blitz::Array<int, 1>& bs_n_brdf_parameters_in) = 0;
  
  virtual const int& bs_n_user_doublets() const = 0;
  virtual void bs_n_user_doublets(const int& bs_n_user_doublets_in) = 0;
  
  virtual const int& bs_n_user_obsgeoms() const = 0;
  virtual void bs_n_user_obsgeoms(const int& bs_n_user_obsgeoms_in) = 0;
  
  virtual const int& bs_n_user_relazms() const = 0;
  virtual void bs_n_user_relazms(const int& bs_n_user_relazms_in) = 0;
  
  virtual const int& bs_n_user_streams() const = 0;
  virtual void bs_n_user_streams(const int& bs_n_user_streams_in) = 0;
  
  virtual const int& bs_nbeams() const = 0;
  virtual void bs_nbeams(const int& bs_nbeams_in) = 0;
  
  virtual const int& bs_nstreams() const = 0;
  virtual void bs_nstreams(const int& bs_nstreams_in) = 0;
  
  virtual const int& bs_nstreams_brdf() const = 0;
  virtual void bs_nstreams_brdf(const int& bs_nstreams_brdf_in) = 0;
  
  virtual const double& bs_salinity() const = 0;
  virtual void bs_salinity(const double& bs_salinity_in) = 0;
  
  virtual const blitz::Array<double, 1>& bs_user_angles_input() const = 0;
  virtual void bs_user_angles_input(const blitz::Array<double, 1>& bs_user_angles_input_in) = 0;
  
  virtual const blitz::Array<double, 2>& bs_user_doublets() const = 0;
  virtual void bs_user_doublets(const blitz::Array<double, 2>& bs_user_doublets_in) = 0;
  
  virtual const blitz::Array<double, 2>& bs_user_obsgeoms() const = 0;
  virtual void bs_user_obsgeoms(const blitz::Array<double, 2>& bs_user_obsgeoms_in) = 0;
  
  virtual const blitz::Array<double, 1>& bs_user_relazms() const = 0;
  virtual void bs_user_relazms(const blitz::Array<double, 1>& bs_user_relazms_in) = 0;
  
  virtual const double& bs_wavelength() const = 0;
  virtual void bs_wavelength(const double& bs_wavelength_in) = 0;
  
  virtual const blitz::Array<int, 1>& bs_which_brdf() const = 0;
  virtual void bs_which_brdf(const blitz::Array<int, 1>& bs_which_brdf_in) = 0;
  
  virtual const blitz::Array<double, 1>& bs_winddir() const = 0;
  virtual void bs_winddir(const blitz::Array<double, 1>& bs_winddir_in) = 0;
  
  virtual const double& bs_windspeed() const = 0;
  virtual void bs_windspeed(const double& bs_windspeed_in) = 0;
  
  virtual const double& bs_wsa_value() const = 0;
  virtual void bs_wsa_value(const double& bs_wsa_value_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Brdf_Sup_Outputs_Base;

class Spurr_Brdf_Sup_Outputs_Base : public Spurr_Type_Structure {
public:

  virtual const double& bs_bsa_calculated() const = 0;
  virtual void bs_bsa_calculated(const double& bs_bsa_calculated_in) = 0;
  
  virtual const blitz::Array<double, 1>& bs_bsa_kernels() const = 0;
  virtual void bs_bsa_kernels(const blitz::Array<double, 1>& bs_bsa_kernels_in) = 0;
  
  virtual const double& bs_wsa_calculated() const = 0;
  virtual void bs_wsa_calculated(const double& bs_wsa_calculated_in) = 0;
  
  virtual const blitz::Array<double, 1>& bs_wsa_kernels() const = 0;
  virtual void bs_wsa_kernels(const blitz::Array<double, 1>& bs_wsa_kernels_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * bs_brdf_f
   * bs_brdf_f_0
   * bs_dbounce_brdfunc
   * bs_emissivity
   * bs_user_brdf_f
   * bs_user_brdf_f_0
   * bs_user_emissivity
   
  */
};

%nodefaultctor Spurr_Brdf_Input_Exception_Handling_Base;

class Spurr_Brdf_Input_Exception_Handling_Base : public Spurr_Type_Structure {
public:

  virtual const std::vector< std::string > bs_inputactions() const = 0;
  
  virtual const std::vector< std::string > bs_inputmessages() const = 0;
  
  virtual const int& bs_ninputmessages() const = 0;
  virtual void bs_ninputmessages(const int& bs_ninputmessages_in) = 0;
  
  virtual const int& bs_status_inputread() const = 0;
  virtual void bs_status_inputread(const int& bs_status_inputread_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Brdf_Output_Exception_Handling_Base;

class Spurr_Brdf_Output_Exception_Handling_Base : public Spurr_Type_Structure {
public:

  virtual const int& bs_noutputmessages() const = 0;
  virtual void bs_noutputmessages(const int& bs_noutputmessages_in) = 0;
  
  virtual const std::vector< std::string > bs_outputmessages() const = 0;
  
  virtual const int& bs_status_output() const = 0;
  virtual void bs_status_output(const int& bs_status_output_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Sleave_Sup_Inputs_Base;

class Sleave_Sup_Inputs_Base : public Spurr_Type_Structure {
public:

  virtual const bool sl_azimuthdep() const = 0;
  virtual void sl_azimuthdep(const bool& sl_azimuthdep_in) = 0;
  
  virtual const blitz::Array<double, 1>& sl_beam_szas() const = 0;
  virtual void sl_beam_szas(const blitz::Array<double, 1>& sl_beam_szas_in) = 0;
  
  virtual const double& sl_chlorconc() const = 0;
  virtual void sl_chlorconc(const double& sl_chlorconc_in) = 0;
  
  virtual const bool sl_do_doublet_geometry() const = 0;
  virtual void sl_do_doublet_geometry(const bool& sl_do_doublet_geometry_in) = 0;
  
  virtual const bool sl_do_exact() const = 0;
  virtual void sl_do_exact(const bool& sl_do_exact_in) = 0;
  
  virtual const bool sl_do_exactonly() const = 0;
  virtual void sl_do_exactonly(const bool& sl_do_exactonly_in) = 0;
  
  virtual const bool sl_do_facetisotropy() const = 0;
  virtual void sl_do_facetisotropy(const bool& sl_do_facetisotropy_in) = 0;
  
  virtual const bool sl_do_fluorescence() const = 0;
  virtual void sl_do_fluorescence(const bool& sl_do_fluorescence_in) = 0;
  
  virtual const bool sl_do_foamoption() const = 0;
  virtual void sl_do_foamoption(const bool& sl_do_foamoption_in) = 0;
  
  virtual const bool sl_do_fourier_output() const = 0;
  virtual void sl_do_fourier_output(const bool& sl_do_fourier_output_in) = 0;
  
  virtual const bool sl_do_glintshadow() const = 0;
  virtual void sl_do_glintshadow(const bool& sl_do_glintshadow_in) = 0;
  
  virtual const bool sl_do_isotropic() const = 0;
  virtual void sl_do_isotropic(const bool& sl_do_isotropic_in) = 0;
  
  virtual const bool sl_do_roughsurface() const = 0;
  virtual void sl_do_roughsurface(const bool& sl_do_roughsurface_in) = 0;
  
  virtual const bool sl_do_sleaving() const = 0;
  virtual void sl_do_sleaving(const bool& sl_do_sleaving_in) = 0;
  
  virtual const bool sl_do_solar_sources() const = 0;
  virtual void sl_do_solar_sources(const bool& sl_do_solar_sources_in) = 0;
  
  virtual const bool sl_do_user_obsgeoms() const = 0;
  virtual void sl_do_user_obsgeoms(const bool& sl_do_user_obsgeoms_in) = 0;
  
  virtual const bool sl_do_user_streams() const = 0;
  virtual void sl_do_user_streams(const bool& sl_do_user_streams_in) = 0;
  
  virtual const double& sl_fl_amplitude755() const = 0;
  virtual void sl_fl_amplitude755(const double& sl_fl_amplitude755_in) = 0;
  
  virtual const bool sl_fl_do_datagaussian() const = 0;
  virtual void sl_fl_do_datagaussian(const bool& sl_fl_do_datagaussian_in) = 0;
  
  virtual const blitz::Array<int, 1>& sl_fl_epoch() const = 0;
  virtual void sl_fl_epoch(const blitz::Array<int, 1>& sl_fl_epoch_in) = 0;
  
  virtual const blitz::Array<double, 2>& sl_fl_inputgaussians() const = 0;
  virtual void sl_fl_inputgaussians(const blitz::Array<double, 2>& sl_fl_inputgaussians_in) = 0;
  
  virtual const double& sl_fl_latitude() const = 0;
  virtual void sl_fl_latitude(const double& sl_fl_latitude_in) = 0;
  
  virtual const double& sl_fl_longitude() const = 0;
  virtual void sl_fl_longitude(const double& sl_fl_longitude_in) = 0;
  
  virtual const double& sl_fl_wavelength() const = 0;
  virtual void sl_fl_wavelength(const double& sl_fl_wavelength_in) = 0;
  
  virtual const int& sl_n_user_doublets() const = 0;
  virtual void sl_n_user_doublets(const int& sl_n_user_doublets_in) = 0;
  
  virtual const int& sl_n_user_obsgeoms() const = 0;
  virtual void sl_n_user_obsgeoms(const int& sl_n_user_obsgeoms_in) = 0;
  
  virtual const int& sl_n_user_relazms() const = 0;
  virtual void sl_n_user_relazms(const int& sl_n_user_relazms_in) = 0;
  
  virtual const int& sl_n_user_streams() const = 0;
  virtual void sl_n_user_streams(const int& sl_n_user_streams_in) = 0;
  
  virtual const int& sl_nbeams() const = 0;
  virtual void sl_nbeams(const int& sl_nbeams_in) = 0;
  
  virtual const int& sl_nstreams() const = 0;
  virtual void sl_nstreams(const int& sl_nstreams_in) = 0;
  
  virtual const double& sl_salinity() const = 0;
  virtual void sl_salinity(const double& sl_salinity_in) = 0;
  
  virtual const blitz::Array<double, 1>& sl_user_angles_input() const = 0;
  virtual void sl_user_angles_input(const blitz::Array<double, 1>& sl_user_angles_input_in) = 0;
  
  virtual const blitz::Array<double, 2>& sl_user_doublets() const = 0;
  virtual void sl_user_doublets(const blitz::Array<double, 2>& sl_user_doublets_in) = 0;
  
  virtual const blitz::Array<double, 2>& sl_user_obsgeoms() const = 0;
  virtual void sl_user_obsgeoms(const blitz::Array<double, 2>& sl_user_obsgeoms_in) = 0;
  
  virtual const blitz::Array<double, 1>& sl_user_relazms() const = 0;
  virtual void sl_user_relazms(const blitz::Array<double, 1>& sl_user_relazms_in) = 0;
  
  virtual const double& sl_wavelength() const = 0;
  virtual void sl_wavelength(const double& sl_wavelength_in) = 0;
  
  virtual const blitz::Array<double, 1>& sl_winddir() const = 0;
  virtual void sl_winddir(const blitz::Array<double, 1>& sl_winddir_in) = 0;
  
  virtual const double& sl_windspeed() const = 0;
  virtual void sl_windspeed(const double& sl_windspeed_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Lincontrol_Base;

class Spurr_Fixed_Lincontrol_Base : public Spurr_Type_Structure {
public:

  virtual const std::vector< std::string > ts_columnwf_names() const = 0;
  
  virtual const blitz::Array<bool, 1> ts_layer_vary_flag() const = 0;
  virtual void ts_layer_vary_flag(const blitz::Array<bool, 1>& ts_layer_vary_flag_in) = 0;
  
  virtual const blitz::Array<int, 1>& ts_layer_vary_number() const = 0;
  virtual void ts_layer_vary_number(const blitz::Array<int, 1>& ts_layer_vary_number_in) = 0;
  
  virtual const int& ts_n_sleave_wfs() const = 0;
  virtual void ts_n_sleave_wfs(const int& ts_n_sleave_wfs_in) = 0;
  
  virtual const int& ts_n_surface_wfs() const = 0;
  virtual void ts_n_surface_wfs(const int& ts_n_surface_wfs_in) = 0;
  
  virtual const int& ts_n_totalcolumn_wfs() const = 0;
  virtual void ts_n_totalcolumn_wfs(const int& ts_n_totalcolumn_wfs_in) = 0;
  
  virtual const std::vector< std::string > ts_profilewf_names() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Linoptical_Base;

class Spurr_Fixed_Linoptical_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<double, 2>& ts_l_deltau_vert_input() const = 0;
  virtual void ts_l_deltau_vert_input(const blitz::Array<double, 2>& ts_l_deltau_vert_input_in) = 0;
  
  virtual const blitz::Array<double, 2>& ts_l_omega_total_input() const = 0;
  virtual void ts_l_omega_total_input(const blitz::Array<double, 2>& ts_l_omega_total_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Lininputs_Base;

class Spurr_Fixed_Lininputs_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Fixed_Lincontrol_Base& fixed_lincontrol_base() = 0;
  virtual const Spurr_Fixed_Lincontrol_Base& fixed_lincontrol_base() const = 0;
  
  virtual Spurr_Fixed_Linoptical_Base& fixed_linoptical_base() = 0;
  virtual const Spurr_Fixed_Linoptical_Base& fixed_linoptical_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Lincontrol_Base;

class Spurr_Modified_Lincontrol_Base : public Spurr_Type_Structure {
public:

  virtual const bool ts_do_atmos_lbbf() const = 0;
  virtual void ts_do_atmos_lbbf(const bool& ts_do_atmos_lbbf_in) = 0;
  
  virtual const bool ts_do_atmos_linearization() const = 0;
  virtual void ts_do_atmos_linearization(const bool& ts_do_atmos_linearization_in) = 0;
  
  virtual const bool ts_do_column_linearization() const = 0;
  virtual void ts_do_column_linearization(const bool& ts_do_column_linearization_in) = 0;
  
  virtual const bool ts_do_linearization() const = 0;
  virtual void ts_do_linearization(const bool& ts_do_linearization_in) = 0;
  
  virtual const bool ts_do_profile_linearization() const = 0;
  virtual void ts_do_profile_linearization(const bool& ts_do_profile_linearization_in) = 0;
  
  virtual const bool ts_do_simulation_only() const = 0;
  virtual void ts_do_simulation_only(const bool& ts_do_simulation_only_in) = 0;
  
  virtual const bool ts_do_sleave_wfs() const = 0;
  virtual void ts_do_sleave_wfs(const bool& ts_do_sleave_wfs_in) = 0;
  
  virtual const bool ts_do_surface_lbbf() const = 0;
  virtual void ts_do_surface_lbbf(const bool& ts_do_surface_lbbf_in) = 0;
  
  virtual const bool ts_do_surface_linearization() const = 0;
  virtual void ts_do_surface_linearization(const bool& ts_do_surface_linearization_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Lininputs_Base;

class Spurr_Modified_Lininputs_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Modified_Lincontrol_Base& modified_lincontrol_base() = 0;
  virtual const Spurr_Modified_Lincontrol_Base& modified_lincontrol_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Linatmos_Base;

class Spurr_Linatmos_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_abbwfs_fluxes
   * ts_abbwfs_jacobians
   
  */
};

%nodefaultctor Spurr_Linsurf_Base;

class Spurr_Linsurf_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_flux_diffuse_surfwf
   * ts_ls_layer_mssts
   * ts_ls_surf_mssts
   * ts_sbbwfs_fluxes
   * ts_sbbwfs_jacobians
   * ts_surfacewf
   
  */
};

%nodefaultctor Spurr_Linoutputs_Base;

class Spurr_Linoutputs_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Linatmos_Base& linatmos_base() = 0;
  virtual const Spurr_Linatmos_Base& linatmos_base() const = 0;
  
  virtual Spurr_Linsurf_Base& linsurf_base() = 0;
  virtual const Spurr_Linsurf_Base& linsurf_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Lin_Sup_Brdf_Base;

class Spurr_Lin_Sup_Brdf_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_ls_brdf_f
   * ts_ls_brdf_f_0
   * ts_ls_emissivity
   * ts_ls_exactdb_brdfunc
   * ts_ls_user_brdf_f
   * ts_ls_user_brdf_f_0
   * ts_ls_user_emissivity
   
  */
};

%nodefaultctor Spurr_Lin_Sup_Sleave_Base;

class Spurr_Lin_Sup_Sleave_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_lssl_slterm_f_0
   * ts_lssl_slterm_isotropic
   * ts_lssl_slterm_userangles
   * ts_lssl_user_slterm_f_0
   
  */
};

%nodefaultctor Spurr_Lin_Sup_Ss_Surf_Base;

class Spurr_Lin_Sup_Ss_Surf_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_surfacewf_db
   
  */
};

%nodefaultctor Spurr_Lin_Sup_Ss_Base;

class Spurr_Lin_Sup_Ss_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Lin_Sup_Ss_Surf_Base& linsup_ss_surf_base() = 0;
  virtual const Spurr_Lin_Sup_Ss_Surf_Base& linsup_ss_surf_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Lin_Sup_Inout_Base;

class Spurr_Lin_Sup_Inout_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Lin_Sup_Brdf_Base& linsup_brdf_base() = 0;
  virtual const Spurr_Lin_Sup_Brdf_Base& linsup_brdf_base() const = 0;
  
  virtual Spurr_Lin_Sup_Sleave_Base& linsup_sleave_base() = 0;
  virtual const Spurr_Lin_Sup_Sleave_Base& linsup_sleave_base() const = 0;
  
  virtual Spurr_Lin_Sup_Ss_Base& linsup_ss_base() = 0;
  virtual const Spurr_Lin_Sup_Ss_Base& linsup_ss_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Main_Outputs_Base;

class Spurr_Main_Outputs_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<int, 1>& ts_fourier_saved() const = 0;
  virtual void ts_fourier_saved(const blitz::Array<int, 1>& ts_fourier_saved_in) = 0;
  
  virtual const blitz::Array<double, 2>& ts_lostrans() const = 0;
  virtual void ts_lostrans(const blitz::Array<double, 2>& ts_lostrans_in) = 0;
  
  virtual const int& ts_n_geometries() const = 0;
  virtual void ts_n_geometries(const int& ts_n_geometries_in) = 0;
  
  virtual const blitz::Array<double, 2>& ts_pathgeoms() const = 0;
  virtual void ts_pathgeoms(const blitz::Array<double, 2>& ts_pathgeoms_in) = 0;
  
  virtual const double& ts_planetary_sbterm() const = 0;
  virtual void ts_planetary_sbterm(const double& ts_planetary_sbterm_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_solarbeam_boatrans() const = 0;
  virtual void ts_solarbeam_boatrans(const blitz::Array<double, 1>& ts_solarbeam_boatrans_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_albmed_fluxes
   * ts_albmed_user
   * ts_contribs
   * ts_dnflux_direct
   * ts_flux_diffuse
   * ts_layer_mssts
   * ts_planetary_transterm
   * ts_surf_mssts
   * ts_trnmed_fluxes
   * ts_trnmed_user
   
  */
};

%nodefaultctor Spurr_Wladjusted_Outputs_Base;

class Spurr_Wladjusted_Outputs_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_wladjusted_direct
   * ts_wladjusted_f_ords_0
   * ts_wladjusted_f_user_0
   * ts_wladjusted_isotropic
   
  */
};

%nodefaultctor Spurr_Exception_Handling_Base;

class Spurr_Exception_Handling_Base : public Spurr_Type_Structure {
public:

  virtual const std::vector< std::string > ts_actions() const = 0;
  
  virtual const std::vector< std::string > ts_checkmessages() const = 0;
  
  virtual const std::string ts_message() const = 0;
  
  virtual const int& ts_ncheckmessages() const = 0;
  virtual void ts_ncheckmessages(const int& ts_ncheckmessages_in) = 0;
  
  virtual const int& ts_status_calculation() const = 0;
  virtual void ts_status_calculation(const int& ts_status_calculation_in) = 0;
  
  virtual const int& ts_status_inputcheck() const = 0;
  virtual void ts_status_inputcheck(const int& ts_status_inputcheck_in) = 0;
  
  virtual const std::string ts_trace_1() const = 0;
  
  virtual const std::string ts_trace_2() const = 0;
  
  virtual const std::string ts_trace_3() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Input_Exception_Handling_Base;

class Spurr_Input_Exception_Handling_Base : public Spurr_Type_Structure {
public:

  virtual const std::vector< std::string > ts_inputactions() const = 0;
  
  virtual const std::vector< std::string > ts_inputmessages() const = 0;
  
  virtual const int& ts_ninputmessages() const = 0;
  virtual void ts_ninputmessages(const int& ts_ninputmessages_in) = 0;
  
  virtual const int& ts_status_inputread() const = 0;
  virtual void ts_status_inputread(const int& ts_status_inputread_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Outputs_Base;

class Spurr_Outputs_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Main_Outputs_Base& main_outputs_base() = 0;
  virtual const Spurr_Main_Outputs_Base& main_outputs_base() const = 0;
  
  virtual Spurr_Exception_Handling_Base& exception_handling_base() = 0;
  virtual const Spurr_Exception_Handling_Base& exception_handling_base() const = 0;
  
  virtual Spurr_Wladjusted_Outputs_Base& wladjusted_outputs_base() = 0;
  virtual const Spurr_Wladjusted_Outputs_Base& wladjusted_outputs_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Sup_Brdf_Base;

class Spurr_Sup_Brdf_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_brdf_f
   * ts_brdf_f_0
   * ts_emissivity
   * ts_exactdb_brdfunc
   * ts_user_brdf_f
   * ts_user_brdf_f_0
   * ts_user_emissivity
   
  */
};

%nodefaultctor Spurr_Sup_Sleave_Base;

class Spurr_Sup_Sleave_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_slterm_f_0
   * ts_slterm_isotropic
   * ts_slterm_userangles
   * ts_user_slterm_f_0
   
  */
};

%nodefaultctor Spurr_Sup_Ss_Base;

class Spurr_Sup_Ss_Base : public Spurr_Type_Structure {
public:

  virtual void print(std::ostream &output_stream) const = 0;
  
  /* 

   Common methods in derived classes with mismatching signatures: 
   * ts_contribs_ss
   
  */
};

%nodefaultctor Spurr_Sup_Inout_Base;

class Spurr_Sup_Inout_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Sup_Brdf_Base& sup_brdf_base() = 0;
  virtual const Spurr_Sup_Brdf_Base& sup_brdf_base() const = 0;
  
  virtual Spurr_Sup_Sleave_Base& sup_sleave_base() = 0;
  virtual const Spurr_Sup_Sleave_Base& sup_sleave_base() const = 0;
  
  virtual Spurr_Sup_Ss_Base& sup_ss_base() = 0;
  virtual const Spurr_Sup_Ss_Base& sup_ss_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Boolean_Base;

class Spurr_Fixed_Boolean_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<bool, 1> ts_do_albtrn_media() const = 0;
  virtual void ts_do_albtrn_media(const blitz::Array<bool, 1>& ts_do_albtrn_media_in) = 0;
  
  virtual const bool ts_do_boa_illumination() const = 0;
  virtual void ts_do_boa_illumination(const bool& ts_do_boa_illumination_in) = 0;
  
  virtual const bool ts_do_dnwelling() const = 0;
  virtual void ts_do_dnwelling(const bool& ts_do_dnwelling_in) = 0;
  
  virtual const bool ts_do_fluorescence() const = 0;
  virtual void ts_do_fluorescence(const bool& ts_do_fluorescence_in) = 0;
  
  virtual const bool ts_do_fullrad_mode() const = 0;
  virtual void ts_do_fullrad_mode(const bool& ts_do_fullrad_mode_in) = 0;
  
  virtual const bool ts_do_mssts() const = 0;
  virtual void ts_do_mssts(const bool& ts_do_mssts_in) = 0;
  
  virtual const bool ts_do_plane_parallel() const = 0;
  virtual void ts_do_plane_parallel(const bool& ts_do_plane_parallel_in) = 0;
  
  virtual const bool ts_do_planetary_problem() const = 0;
  virtual void ts_do_planetary_problem(const bool& ts_do_planetary_problem_in) = 0;
  
  virtual const bool ts_do_sl_isotropic() const = 0;
  virtual void ts_do_sl_isotropic(const bool& ts_do_sl_isotropic_in) = 0;
  
  virtual const bool ts_do_surface_emission() const = 0;
  virtual void ts_do_surface_emission(const bool& ts_do_surface_emission_in) = 0;
  
  virtual const bool ts_do_surface_leaving() const = 0;
  virtual void ts_do_surface_leaving(const bool& ts_do_surface_leaving_in) = 0;
  
  virtual const bool ts_do_tf_iteration() const = 0;
  virtual void ts_do_tf_iteration(const bool& ts_do_tf_iteration_in) = 0;
  
  virtual const bool ts_do_thermal_emission() const = 0;
  virtual void ts_do_thermal_emission(const bool& ts_do_thermal_emission_in) = 0;
  
  virtual const bool ts_do_toa_contribs() const = 0;
  virtual void ts_do_toa_contribs(const bool& ts_do_toa_contribs_in) = 0;
  
  virtual const bool ts_do_toa_illumination() const = 0;
  virtual void ts_do_toa_illumination(const bool& ts_do_toa_illumination_in) = 0;
  
  virtual const bool ts_do_upwelling() const = 0;
  virtual void ts_do_upwelling(const bool& ts_do_upwelling_in) = 0;
  
  virtual const bool ts_do_water_leaving() const = 0;
  virtual void ts_do_water_leaving(const bool& ts_do_water_leaving_in) = 0;
  
  virtual const bool ts_do_wladjusted_output() const = 0;
  virtual void ts_do_wladjusted_output(const bool& ts_do_wladjusted_output_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Control_Base;

class Spurr_Fixed_Control_Base : public Spurr_Type_Structure {
public:

  virtual const double& ts_asymtx_tolerance() const = 0;
  virtual void ts_asymtx_tolerance(const double& ts_asymtx_tolerance_in) = 0;
  
  virtual const double& ts_boa_illumination() const = 0;
  virtual void ts_boa_illumination(const double& ts_boa_illumination_in) = 0;
  
  virtual const double& ts_fourier_accuracy() const = 0;
  virtual void ts_fourier_accuracy(const double& ts_lidort_accuracy_in) = 0;
  
  virtual const int& ts_n_thermal_coeffs() const = 0;
  virtual void ts_n_thermal_coeffs(const int& ts_n_thermal_coeffs_in) = 0;
  
  virtual const int& ts_nfinelayers() const = 0;
  virtual void ts_nfinelayers(const int& ts_nfinelayers_in) = 0;
  
  virtual const int& ts_nlayers() const = 0;
  virtual void ts_nlayers(const int& ts_nlayers_in) = 0;
  
  virtual const int& ts_nstreams() const = 0;
  virtual void ts_nstreams(const int& ts_nstreams_in) = 0;
  
  virtual const int& ts_taylor_order() const = 0;
  virtual void ts_taylor_order(const int& ts_taylor_order_in) = 0;
  
  virtual const double& ts_tf_criterion() const = 0;
  virtual void ts_tf_criterion(const double& ts_tf_criterion_in) = 0;
  
  virtual const int& ts_tf_maxiter() const = 0;
  virtual void ts_tf_maxiter(const int& ts_tf_maxiter_in) = 0;
  
  virtual const double& ts_toa_illumination() const = 0;
  virtual void ts_toa_illumination(const double& ts_toa_illumination_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Sunrays_Base;

class Spurr_Fixed_Sunrays_Base : public Spurr_Type_Structure {
public:

  virtual const double& ts_flux_factor() const = 0;
  virtual void ts_flux_factor(const double& ts_flux_factor_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Uservalues_Base;

class Spurr_Fixed_Uservalues_Base : public Spurr_Type_Structure {
public:

  virtual const int& ts_n_user_levels() const = 0;
  virtual void ts_n_user_levels(const int& ts_n_user_levels_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Chapman_Base;

class Spurr_Fixed_Chapman_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<int, 1>& ts_finegrid() const = 0;
  virtual void ts_finegrid(const blitz::Array<int, 1>& ts_finegrid_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_height_grid() const = 0;
  virtual void ts_height_grid(const blitz::Array<double, 1>& ts_height_grid_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_pressure_grid() const = 0;
  virtual void ts_pressure_grid(const blitz::Array<double, 1>& ts_pressure_grid_in) = 0;
  
  virtual const double& ts_rfindex_parameter() const = 0;
  virtual void ts_rfindex_parameter(const double& ts_rfindex_parameter_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_temperature_grid() const = 0;
  virtual void ts_temperature_grid(const blitz::Array<double, 1>& ts_temperature_grid_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Optical_Base;

class Spurr_Fixed_Optical_Base : public Spurr_Type_Structure {
public:

  virtual const double& ts_atmos_wavelength() const = 0;
  virtual void ts_atmos_wavelength(const double& ts_atmos_wavelength_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_deltau_vert_input() const = 0;
  virtual void ts_deltau_vert_input(const blitz::Array<double, 1>& ts_deltau_vert_input_in) = 0;
  
  virtual const double& ts_lambertian_albedo() const = 0;
  virtual void ts_lambertian_albedo(const double& ts_lambertian_albedo_in) = 0;
  
  virtual const double& ts_surface_bb_input() const = 0;
  virtual void ts_surface_bb_input(const double& ts_surface_bb_input_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_thermal_bb_input() const = 0;
  virtual void ts_thermal_bb_input(const blitz::Array<double, 1>& ts_thermal_bb_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Write_Base;

class Spurr_Fixed_Write_Base : public Spurr_Type_Structure {
public:

  virtual const bool ts_do_debug_write() const = 0;
  virtual void ts_do_debug_write(const bool& ts_do_debug_write_in) = 0;
  
  virtual const bool ts_do_write_fourier() const = 0;
  virtual void ts_do_write_fourier(const bool& ts_do_write_fourier_in) = 0;
  
  virtual const bool ts_do_write_input() const = 0;
  virtual void ts_do_write_input(const bool& ts_do_write_input_in) = 0;
  
  virtual const bool ts_do_write_results() const = 0;
  virtual void ts_do_write_results(const bool& ts_do_write_results_in) = 0;
  
  virtual const bool ts_do_write_scenario() const = 0;
  virtual void ts_do_write_scenario(const bool& ts_do_write_scenario_in) = 0;
  
  virtual const std::string ts_fourier_write_filename() const = 0;
  
  virtual const std::string ts_input_write_filename() const = 0;
  
  virtual const std::string ts_results_write_filename() const = 0;
  
  virtual const std::string ts_scenario_write_filename() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Fixed_Inputs_Base;

class Spurr_Fixed_Inputs_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Fixed_Boolean_Base& fixed_boolean_base() = 0;
  virtual const Spurr_Fixed_Boolean_Base& fixed_boolean_base() const = 0;
  
  virtual Spurr_Fixed_Chapman_Base& fixed_chapman_base() = 0;
  virtual const Spurr_Fixed_Chapman_Base& fixed_chapman_base() const = 0;
  
  virtual Spurr_Fixed_Control_Base& fixed_control_base() = 0;
  virtual const Spurr_Fixed_Control_Base& fixed_control_base() const = 0;
  
  virtual Spurr_Fixed_Optical_Base& fixed_optical_base() = 0;
  virtual const Spurr_Fixed_Optical_Base& fixed_optical_base() const = 0;
  
  virtual Spurr_Fixed_Sunrays_Base& fixed_sunrays_base() = 0;
  virtual const Spurr_Fixed_Sunrays_Base& fixed_sunrays_base() const = 0;
  
  virtual Spurr_Fixed_Uservalues_Base& fixed_uservalues_base() = 0;
  virtual const Spurr_Fixed_Uservalues_Base& fixed_uservalues_base() const = 0;
  
  virtual Spurr_Fixed_Write_Base& fixed_write_base() = 0;
  virtual const Spurr_Fixed_Write_Base& fixed_write_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Boolean_Base;

class Spurr_Modified_Boolean_Base : public Spurr_Type_Structure {
public:

  virtual const bool ts_do_additional_mvout() const = 0;
  virtual void ts_do_additional_mvout(const bool& ts_do_additional_mvout_in) = 0;
  
  virtual const bool ts_do_bvp_telescoping() const = 0;
  virtual void ts_do_bvp_telescoping(const bool& ts_do_bvp_telescoping_in) = 0;
  
  virtual const bool ts_do_chapman_function() const = 0;
  virtual void ts_do_chapman_function(const bool& ts_do_chapman_function_in) = 0;
  
  virtual const bool ts_do_deltam_scaling() const = 0;
  virtual void ts_do_deltam_scaling(const bool& ts_do_deltam_scaling_in) = 0;
  
  virtual const bool ts_do_double_convtest() const = 0;
  virtual void ts_do_double_convtest(const bool& ts_do_double_convtest_in) = 0;
  
  virtual const bool ts_do_doublet_geometry() const = 0;
  virtual void ts_do_doublet_geometry(const bool& ts_do_doublet_geometry_in) = 0;
  
  virtual const bool ts_do_external_wleave() const = 0;
  virtual void ts_do_external_wleave(const bool& ts_do_external_wleave_in) = 0;
  
  virtual const bool ts_do_focorr() const = 0;
  virtual void ts_do_focorr(const bool& ts_do_focorr_in) = 0;
  
  virtual const bool ts_do_focorr_external() const = 0;
  virtual void ts_do_focorr_external(const bool& ts_do_focorr_external_in) = 0;
  
  virtual const bool ts_do_focorr_nadir() const = 0;
  virtual void ts_do_focorr_nadir(const bool& ts_do_focorr_nadir_in) = 0;
  
  virtual const bool ts_do_focorr_outgoing() const = 0;
  virtual void ts_do_focorr_outgoing(const bool& ts_do_focorr_outgoing_in) = 0;
  
  virtual const bool ts_do_mvout_only() const = 0;
  virtual void ts_do_mvout_only(const bool& ts_do_mvout_only_in) = 0;
  
  virtual const bool ts_do_observation_geometry() const = 0;
  virtual void ts_do_observation_geometry(const bool& ts_do_observation_geometry_in) = 0;
  
  virtual const bool ts_do_rayleigh_only() const = 0;
  virtual void ts_do_rayleigh_only(const bool& ts_do_rayleigh_only_in) = 0;
  
  virtual const bool ts_do_refractive_geometry() const = 0;
  virtual void ts_do_refractive_geometry(const bool& ts_do_refractive_geometry_in) = 0;
  
  virtual const bool ts_do_solar_sources() const = 0;
  virtual void ts_do_solar_sources(const bool& ts_do_solar_sources_in) = 0;
  
  virtual const bool ts_do_solution_saving() const = 0;
  virtual void ts_do_solution_saving(const bool& ts_do_solution_saving_in) = 0;
  
  virtual const bool ts_do_sscorr_truncation() const = 0;
  virtual void ts_do_sscorr_truncation(const bool& ts_do_sscorr_truncation_in) = 0;
  
  virtual const bool ts_do_thermal_transonly() const = 0;
  virtual void ts_do_thermal_transonly(const bool& ts_do_thermal_transonly_in) = 0;
  
  virtual const bool ts_do_user_streams() const = 0;
  virtual void ts_do_user_streams(const bool& ts_do_user_streams_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Control_Base;

class Spurr_Modified_Control_Base : public Spurr_Type_Structure {
public:

  virtual const int& ts_nmoments_input() const = 0;
  virtual void ts_nmoments_input(const int& ts_nmoments_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Sunrays_Base;

class Spurr_Modified_Sunrays_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<double, 1>& ts_beam_szas() const = 0;
  virtual void ts_beam_szas(const blitz::Array<double, 1>& ts_beam_szas_in) = 0;
  
  virtual const int& ts_nbeams() const = 0;
  virtual void ts_nbeams(const int& ts_nbeams_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Uservalues_Base;

class Spurr_Modified_Uservalues_Base : public Spurr_Type_Structure {
public:

  virtual const double& ts_geometry_specheight() const = 0;
  virtual void ts_geometry_specheight(const double& ts_geometry_specheight_in) = 0;
  
  virtual const int& ts_n_user_doublets() const = 0;
  virtual void ts_n_user_doublets(const int& ts_n_user_doublets_in) = 0;
  
  virtual const int& ts_n_user_obsgeoms() const = 0;
  virtual void ts_n_user_obsgeoms(const int& ts_n_user_obsgeoms_in) = 0;
  
  virtual const int& ts_n_user_relazms() const = 0;
  virtual void ts_n_user_relazms(const int& ts_n_user_relazms_in) = 0;
  
  virtual const int& ts_n_user_streams() const = 0;
  virtual void ts_n_user_streams(const int& ts_n_user_streams_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_user_angles_input() const = 0;
  virtual void ts_user_angles_input(const blitz::Array<double, 1>& ts_user_angles_input_in) = 0;
  
  virtual const blitz::Array<double, 2>& ts_user_doublets() const = 0;
  virtual void ts_user_doublets(const blitz::Array<double, 2>& ts_user_doublets_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_user_levels() const = 0;
  virtual void ts_user_levels(const blitz::Array<double, 1>& ts_user_levels_in) = 0;
  
  virtual const blitz::Array<double, 2>& ts_user_obsgeoms_input() const = 0;
  virtual void ts_user_obsgeoms_input(const blitz::Array<double, 2>& ts_user_obsgeoms_input_in) = 0;
  
  virtual const blitz::Array<double, 1>& ts_user_relazms() const = 0;
  virtual void ts_user_relazms(const blitz::Array<double, 1>& ts_user_relazms_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Chapman_Base;

class Spurr_Modified_Chapman_Base : public Spurr_Type_Structure {
public:

  virtual const double& ts_earth_radius() const = 0;
  virtual void ts_earth_radius(const double& ts_earth_radius_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Optical_Base;

class Spurr_Modified_Optical_Base : public Spurr_Type_Structure {
public:

  virtual const blitz::Array<double, 1>& ts_omega_total_input() const = 0;
  virtual void ts_omega_total_input(const blitz::Array<double, 1>& ts_omega_total_input_in) = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};

%nodefaultctor Spurr_Modified_Inputs_Base;

class Spurr_Modified_Inputs_Base : public Spurr_Type_Structure {
public:

  virtual Spurr_Modified_Boolean_Base& modified_boolean_base() = 0;
  virtual const Spurr_Modified_Boolean_Base& modified_boolean_base() const = 0;
  
  virtual Spurr_Modified_Chapman_Base& modified_chapman_base() = 0;
  virtual const Spurr_Modified_Chapman_Base& modified_chapman_base() const = 0;
  
  virtual Spurr_Modified_Control_Base& modified_control_base() = 0;
  virtual const Spurr_Modified_Control_Base& modified_control_base() const = 0;
  
  virtual Spurr_Modified_Optical_Base& modified_optical_base() = 0;
  virtual const Spurr_Modified_Optical_Base& modified_optical_base() const = 0;
  
  virtual Spurr_Modified_Sunrays_Base& modified_sunrays_base() = 0;
  virtual const Spurr_Modified_Sunrays_Base& modified_sunrays_base() const = 0;
  
  virtual Spurr_Modified_Uservalues_Base& modified_uservalues_base() = 0;
  virtual const Spurr_Modified_Uservalues_Base& modified_uservalues_base() const = 0;
  
  virtual void print(std::ostream &output_stream) const = 0;
  
};



}