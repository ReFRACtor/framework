// This file was auto-generated

%include "fp_common.i"

%include "spurr_interface_types.i"

%{
#include "vlidort_interface_types.h"
%}

%fp_shared_ptr(FullPhysics::VLidort_Pars);

%fp_shared_ptr(FullPhysics::VBrdf_Linsup_Inputs);
%fp_shared_ptr(FullPhysics::VBrdf_Linsup_Outputs);
%fp_shared_ptr(FullPhysics::VBrdf_Sup_Inputs);
%fp_shared_ptr(FullPhysics::VBrdf_Sup_Outputs);
%fp_shared_ptr(FullPhysics::VBrdf_Input_Exception_Handling);
%fp_shared_ptr(FullPhysics::VBrdf_Output_Exception_Handling);
%fp_shared_ptr(FullPhysics::VSleave_Sup_Inputs);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Lincontrol);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Linoptical);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Lininputs);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Lincontrol);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Lininputs);
%fp_shared_ptr(FullPhysics::VLidort_Lincol);
%fp_shared_ptr(FullPhysics::VLidort_Linprof);
%fp_shared_ptr(FullPhysics::VLidort_Linatmos);
%fp_shared_ptr(FullPhysics::VLidort_Linsurf);
%fp_shared_ptr(FullPhysics::VLidort_Linoutputs);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Brdf);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Sleave);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Ss_Col);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Ss_Prof);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Ss_Surf);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Ss);
%fp_shared_ptr(FullPhysics::VLidort_Linsup_Inout);
%fp_shared_ptr(FullPhysics::VLidort_Main_Outputs);
%fp_shared_ptr(FullPhysics::VLidort_Wladjusted_Outputs);
%fp_shared_ptr(FullPhysics::VLidort_Exception_Handling);
%fp_shared_ptr(FullPhysics::VLidort_Input_Exception_Handling);
%fp_shared_ptr(FullPhysics::VLidort_Outputs);
%fp_shared_ptr(FullPhysics::VLidort_Sup_Brdf);
%fp_shared_ptr(FullPhysics::VLidort_Sup_Sleave);
%fp_shared_ptr(FullPhysics::VLidort_Sup_Ss);
%fp_shared_ptr(FullPhysics::VLidort_Sup_Inout);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Boolean);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Control);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Sunrays);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Uservalues);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Chapman);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Optical);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Write);
%fp_shared_ptr(FullPhysics::VLidort_Fixed_Inputs);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Boolean);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Control);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Sunrays);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Uservalues);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Chapman);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Optical);
%fp_shared_ptr(FullPhysics::VLidort_Modified_Inputs);
 

namespace FullPhysics {

%nodefaultctor VLidort_Pars;

class VLidort_Pars : public Spurr_Pars_Base {
public:

  const double& bigexp() const;
  const int& bpdfndvi_idx() const;
  const int& bpdfsoil_idx() const;
  const int& bpdfvegn_idx() const;
  const int& coxmunk_idx() const;
  const double& deg_to_rad() const;
  const int& dnidx() const;
  const double& eps3() const;
  const double& eps4() const;
  const double& eps5() const;
  const double& four() const;
  const int& gisscoxmunk_cri_idx() const;
  const int& gisscoxmunk_idx() const;
  const double& half() const;
  const int& hapke_idx() const;
  const double& hopital_tolerance() const;
  const int& lambertian_idx() const;
  const int& lidense_idx() const;
  const int& lisparse_idx() const;
  const int& max_allstrms() const;
  const int& max_allstrms_p1() const;
  const int& max_atmoswfs() const;
  const int& max_brdf_kernels() const;
  const int& max_brdf_parameters() const;
  const int& max_directions() const;
  const int& max_geometries() const;
  const int& max_messages() const;
  const int& max_msrs_muquad() const;
  const int& max_msrs_phiquad() const;
  const int& max_partlayers() const;
  const int& max_psols() const;
  const int& max_scatpsols() const;
  const int& max_sleavewfs() const;
  const int& max_surfacewfs() const;
  const int& max_szangles() const;
  const double& max_tau_qpath() const;
  const double& max_tau_spath() const;
  const double& max_tau_upath() const;
  const int& max_taylor_terms() const;
  const int& max_thermal_coeffs() const;
  const int& max_user_levels() const;
  const int& max_user_obsgeoms() const;
  const int& max_user_relazms() const;
  const int& max_user_streams() const;
  const int& max_user_vzangles() const;
  const int& max_ustrmstks() const;
  const int& maxbandtotal() const;
  const int& maxbeams() const;
  const int& maxbrdf_idx() const;
  const int& maxevalues() const;
  const int& maxfinelayers() const;
  const int& maxfourier() const;
  const int& maxlayers() const;
  const int& maxmoments() const;
  const int& maxmoments_input() const;
  const int& maxsthalf_brdf() const;
  const int& maxstokes() const;
  const int& maxstokes_sq() const;
  const int& maxstreams() const;
  const int& maxstreams_2() const;
  const int& maxstreams_21() const;
  const int& maxstreams_brdf() const;
  const int& maxstreams_p1() const;
  const int& maxstreams_p2() const;
  const int& maxstreams_scaling() const;
  const int& maxstrmstks() const;
  const int& maxstrmstks_2() const;
  const int& maxstrmstks_21() const;
  const int& maxstrmstks_p1() const;
  const int& maxstrmstks_p2() const;
  const int& maxstrmstks_p4() const;
  const int& maxtotal() const;
  const double& minus_one() const;
  const double& minus_two() const;
  const int& modfresnel_idx() const;
  const int& newcmglint_idx() const;
  const int& newgcmglint_idx() const;
  const double& omega_smallnum() const;
  const double& one() const;
  const double& onep5() const;
  const double& pi2() const;
  const double& pi4() const;
  const double& pie() const;
  const double& pio2() const;
  const double& pio4() const;
  const double& quarter() const;
  const int& rahman_idx() const;
  const int& rossthick_idx() const;
  const int& rossthin_idx() const;
  const int& roujean_idx() const;
  const int& rtkhotspot_idx() const;
  const double& smallnum() const;
  const int& snowbrdf_idx() const;
  const double& taylor_large() const;
  const double& taylor_small() const;
  const double& three() const;
  const double& two() const;
  const int& upidx() const;
  const int& vlidort_debug() const;
  const int& vlidort_errunit() const;
  const int& vlidort_funit() const;
  const int& vlidort_info() const;
  const int& vlidort_inunit() const;
  const int& vlidort_resunit() const;
  const int& vlidort_scenunit() const;
  const int& vlidort_serious() const;
  const int& vlidort_success() const;
  const std::string vlidort_version_number() const;
  const int& vlidort_warning() const;
  const double& zero() const;
  static VLidort_Pars& instance();

  virtual void print(std::ostream &output_stream) const;

  VLidort_Pars();

};


class VBrdf_Linsup_Inputs : public Spurr_Brdf_Lin_Sup_Inputs_Base {
public:
  VBrdf_Linsup_Inputs();
  VBrdf_Linsup_Inputs(const VBrdf_Linsup_Inputs& src);
  ~VBrdf_Linsup_Inputs();

  const bool bs_do_bsavalue_wf() const;
  void bs_do_bsavalue_wf(const bool& bs_do_bsavalue_wf_in);
  
  const blitz::Array<bool, 1> bs_do_kernel_factor_wfs() const;
  void bs_do_kernel_factor_wfs(const blitz::Array<bool, 1>& bs_do_kernel_factor_wfs_in);
  
  const blitz::Array<bool, 2> bs_do_kernel_params_wfs() const;
  void bs_do_kernel_params_wfs(const blitz::Array<bool, 2>& bs_do_kernel_params_wfs_in);
  
  const blitz::Array<bool, 1> bs_do_kparams_derivs() const;
  void bs_do_kparams_derivs(const blitz::Array<bool, 1>& bs_do_kparams_derivs_in);
  
  const bool bs_do_windspeed_wf() const;
  void bs_do_windspeed_wf(const bool& bs_do_windspeed_wf_in);
  
  const bool bs_do_wsavalue_wf() const;
  void bs_do_wsavalue_wf(const bool& bs_do_wsavalue_wf_in);
  
  const int& bs_n_kernel_factor_wfs() const;
  void bs_n_kernel_factor_wfs(const int& bs_n_kernel_factor_wfs_in);
  
  const int& bs_n_kernel_params_wfs() const;
  void bs_n_kernel_params_wfs(const int& bs_n_kernel_params_wfs_in);
  
  const int& bs_n_surface_wfs() const;
  void bs_n_surface_wfs(const int& bs_n_surface_wfs_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VBrdf_Linsup_Outputs : public Spurr_Brdf_Lin_Sup_Outputs_Base {
public:
  VBrdf_Linsup_Outputs();
  VBrdf_Linsup_Outputs(const VBrdf_Linsup_Outputs& src);
  ~VBrdf_Linsup_Outputs();

  const blitz::Array<double, 5>& bs_ls_brdf_f() const;
  void bs_ls_brdf_f(const blitz::Array<double, 5>& bs_ls_brdf_f_in);
  
  const blitz::Array<double, 5>& bs_ls_brdf_f_0() const;
  void bs_ls_brdf_f_0(const blitz::Array<double, 5>& bs_ls_brdf_f_0_in);
  
  const blitz::Array<double, 5>& bs_ls_dbounce_brdfunc() const;
  void bs_ls_dbounce_brdfunc(const blitz::Array<double, 5>& bs_ls_dbounce_brdfunc_in);
  
  const blitz::Array<double, 3>& bs_ls_emissivity() const;
  void bs_ls_emissivity(const blitz::Array<double, 3>& bs_ls_emissivity_in);
  
  const blitz::Array<double, 5>& bs_ls_user_brdf_f() const;
  void bs_ls_user_brdf_f(const blitz::Array<double, 5>& bs_ls_user_brdf_f_in);
  
  const blitz::Array<double, 5>& bs_ls_user_brdf_f_0() const;
  void bs_ls_user_brdf_f_0(const blitz::Array<double, 5>& bs_ls_user_brdf_f_0_in);
  
  const blitz::Array<double, 3>& bs_ls_user_emissivity() const;
  void bs_ls_user_emissivity(const blitz::Array<double, 3>& bs_ls_user_emissivity_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VBrdf_Sup_Inputs : public Spurr_Brdf_Sup_Inputs_Base {
public:
  VBrdf_Sup_Inputs();
  VBrdf_Sup_Inputs(const VBrdf_Sup_Inputs& src);
  ~VBrdf_Sup_Inputs();

  const blitz::Array<double, 1>& bs_beam_szas() const;
  void bs_beam_szas(const blitz::Array<double, 1>& bs_beam_szas_in);
  
  const blitz::Array<double, 1>& bs_brdf_factors() const;
  void bs_brdf_factors(const blitz::Array<double, 1>& bs_brdf_factors_in);
  
  const std::vector< std::string > bs_brdf_names() const;
  
  const blitz::Array<double, 2>& bs_brdf_parameters() const;
  void bs_brdf_parameters(const blitz::Array<double, 2>& bs_brdf_parameters_in);
  
  const double& bs_bsa_value() const;
  void bs_bsa_value(const double& bs_bsa_value_in);
  
  const bool bs_do_brdf_surface() const;
  void bs_do_brdf_surface(const bool& bs_do_brdf_surface_in);
  
  const bool bs_do_bsa_scaling() const;
  void bs_do_bsa_scaling(const bool& bs_do_bsa_scaling_in);
  
  const bool bs_do_directbounce_only() const;
  void bs_do_directbounce_only(const bool& bs_do_directbounce_only_in);
  
  const bool bs_do_doublet_geometry() const;
  void bs_do_doublet_geometry(const bool& bs_do_doublet_geometry_in);
  
  const bool bs_do_facetisotropy() const;
  void bs_do_facetisotropy(const bool& bs_do_facetisotropy_in);
  
  const bool bs_do_foamoption() const;
  void bs_do_foamoption(const bool& bs_do_foamoption_in);
  
  const bool bs_do_glintshadow() const;
  void bs_do_glintshadow(const bool& bs_do_glintshadow_in);
  
  const bool bs_do_glitter_msrcorr() const;
  void bs_do_glitter_msrcorr(const bool& bs_do_glitter_msrcorr_in);
  
  const bool bs_do_glitter_msrcorr_dbonly() const;
  void bs_do_glitter_msrcorr_dbonly(const bool& bs_do_glitter_msrcorr_dbonly_in);
  
  const bool bs_do_newcmglint() const;
  void bs_do_newcmglint(const bool& bs_do_newcmglint_in);
  
  const bool bs_do_newgcmglint() const;
  void bs_do_newgcmglint(const bool& bs_do_newgcmglint_in);
  
  const bool bs_do_shadow_effect() const;
  void bs_do_shadow_effect(const bool& bs_do_shadow_effect_in);
  
  const bool bs_do_solar_sources() const;
  void bs_do_solar_sources(const bool& bs_do_solar_sources_in);
  
  const bool bs_do_surface_emission() const;
  void bs_do_surface_emission(const bool& bs_do_surface_emission_in);
  
  const bool bs_do_user_obsgeoms() const;
  void bs_do_user_obsgeoms(const bool& bs_do_user_obsgeoms_in);
  
  const bool bs_do_user_streams() const;
  void bs_do_user_streams(const bool& bs_do_user_streams_in);
  
  const bool bs_do_wsa_scaling() const;
  void bs_do_wsa_scaling(const bool& bs_do_wsa_scaling_in);
  
  const bool bs_do_wsabsa_output() const;
  void bs_do_wsabsa_output(const bool& bs_do_wsabsa_output_in);
  
  const int& bs_glitter_msrcorr_nmuquad() const;
  void bs_glitter_msrcorr_nmuquad(const int& bs_glitter_msrcorr_nmuquad_in);
  
  const int& bs_glitter_msrcorr_nphiquad() const;
  void bs_glitter_msrcorr_nphiquad(const int& bs_glitter_msrcorr_nphiquad_in);
  
  const int& bs_glitter_msrcorr_order() const;
  void bs_glitter_msrcorr_order(const int& bs_glitter_msrcorr_order_in);
  
  const blitz::Array<bool, 1> bs_lambertian_kernel_flag() const;
  void bs_lambertian_kernel_flag(const blitz::Array<bool, 1>& bs_lambertian_kernel_flag_in);
  
  const int& bs_n_brdf_kernels() const;
  void bs_n_brdf_kernels(const int& bs_n_brdf_kernels_in);
  
  const blitz::Array<int, 1>& bs_n_brdf_parameters() const;
  void bs_n_brdf_parameters(const blitz::Array<int, 1>& bs_n_brdf_parameters_in);
  
  const int& bs_n_user_doublets() const;
  void bs_n_user_doublets(const int& bs_n_user_doublets_in);
  
  const int& bs_n_user_obsgeoms() const;
  void bs_n_user_obsgeoms(const int& bs_n_user_obsgeoms_in);
  
  const int& bs_n_user_relazms() const;
  void bs_n_user_relazms(const int& bs_n_user_relazms_in);
  
  const int& bs_n_user_streams() const;
  void bs_n_user_streams(const int& bs_n_user_streams_in);
  
  const int& bs_nbeams() const;
  void bs_nbeams(const int& bs_nbeams_in);
  
  const int& bs_nstokes() const;
  void bs_nstokes(const int& bs_nstokes_in);
  
  const int& bs_nstreams() const;
  void bs_nstreams(const int& bs_nstreams_in);
  
  const int& bs_nstreams_brdf() const;
  void bs_nstreams_brdf(const int& bs_nstreams_brdf_in);
  
  const double& bs_salinity() const;
  void bs_salinity(const double& bs_salinity_in);
  
  const blitz::Array<double, 1>& bs_user_angles_input() const;
  void bs_user_angles_input(const blitz::Array<double, 1>& bs_user_angles_input_in);
  
  const blitz::Array<double, 2>& bs_user_doublets() const;
  void bs_user_doublets(const blitz::Array<double, 2>& bs_user_doublets_in);
  
  const blitz::Array<double, 2>& bs_user_obsgeoms() const;
  void bs_user_obsgeoms(const blitz::Array<double, 2>& bs_user_obsgeoms_in);
  
  const blitz::Array<double, 1>& bs_user_relazms() const;
  void bs_user_relazms(const blitz::Array<double, 1>& bs_user_relazms_in);
  
  const double& bs_wavelength() const;
  void bs_wavelength(const double& bs_wavelength_in);
  
  const blitz::Array<int, 1>& bs_which_brdf() const;
  void bs_which_brdf(const blitz::Array<int, 1>& bs_which_brdf_in);
  
  const blitz::Array<double, 1>& bs_winddir() const;
  void bs_winddir(const blitz::Array<double, 1>& bs_winddir_in);
  
  const double& bs_windspeed() const;
  void bs_windspeed(const double& bs_windspeed_in);
  
  const double& bs_wsa_value() const;
  void bs_wsa_value(const double& bs_wsa_value_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VBrdf_Sup_Outputs : public Spurr_Brdf_Sup_Outputs_Base {
public:
  VBrdf_Sup_Outputs();
  VBrdf_Sup_Outputs(const VBrdf_Sup_Outputs& src);
  ~VBrdf_Sup_Outputs();

  const blitz::Array<double, 4>& bs_brdf_f() const;
  void bs_brdf_f(const blitz::Array<double, 4>& bs_brdf_f_in);
  
  const blitz::Array<double, 4>& bs_brdf_f_0() const;
  void bs_brdf_f_0(const blitz::Array<double, 4>& bs_brdf_f_0_in);
  
  const double& bs_bsa_calculated() const;
  void bs_bsa_calculated(const double& bs_bsa_calculated_in);
  
  const blitz::Array<double, 1>& bs_bsa_kernels() const;
  void bs_bsa_kernels(const blitz::Array<double, 1>& bs_bsa_kernels_in);
  
  const blitz::Array<double, 4>& bs_dbounce_brdfunc() const;
  void bs_dbounce_brdfunc(const blitz::Array<double, 4>& bs_dbounce_brdfunc_in);
  
  const blitz::Array<double, 2>& bs_emissivity() const;
  void bs_emissivity(const blitz::Array<double, 2>& bs_emissivity_in);
  
  const blitz::Array<double, 4>& bs_user_brdf_f() const;
  void bs_user_brdf_f(const blitz::Array<double, 4>& bs_user_brdf_f_in);
  
  const blitz::Array<double, 4>& bs_user_brdf_f_0() const;
  void bs_user_brdf_f_0(const blitz::Array<double, 4>& bs_user_brdf_f_0_in);
  
  const blitz::Array<double, 2>& bs_user_emissivity() const;
  void bs_user_emissivity(const blitz::Array<double, 2>& bs_user_emissivity_in);
  
  const double& bs_wsa_calculated() const;
  void bs_wsa_calculated(const double& bs_wsa_calculated_in);
  
  const blitz::Array<double, 1>& bs_wsa_kernels() const;
  void bs_wsa_kernels(const blitz::Array<double, 1>& bs_wsa_kernels_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VBrdf_Input_Exception_Handling : public Spurr_Brdf_Input_Exception_Handling_Base {
public:
  VBrdf_Input_Exception_Handling();
  VBrdf_Input_Exception_Handling(const VBrdf_Input_Exception_Handling& src);
  ~VBrdf_Input_Exception_Handling();

  const std::vector< std::string > bs_inputactions() const;
  
  const std::vector< std::string > bs_inputmessages() const;
  
  const int& bs_ninputmessages() const;
  void bs_ninputmessages(const int& bs_ninputmessages_in);
  
  const int& bs_status_inputread() const;
  void bs_status_inputread(const int& bs_status_inputread_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VBrdf_Output_Exception_Handling : public Spurr_Brdf_Output_Exception_Handling_Base {
public:
  VBrdf_Output_Exception_Handling();
  VBrdf_Output_Exception_Handling(const VBrdf_Output_Exception_Handling& src);
  ~VBrdf_Output_Exception_Handling();

  const int& bs_noutputmessages() const;
  void bs_noutputmessages(const int& bs_noutputmessages_in);
  
  const std::vector< std::string > bs_outputmessages() const;
  
  const int& bs_status_output() const;
  void bs_status_output(const int& bs_status_output_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VSleave_Sup_Inputs : public Sleave_Sup_Inputs_Base {
public:
  VSleave_Sup_Inputs();
  VSleave_Sup_Inputs(const VSleave_Sup_Inputs& src);
  ~VSleave_Sup_Inputs();

  const bool sl_azimuthdep() const;
  void sl_azimuthdep(const bool& sl_azimuthdep_in);
  
  const blitz::Array<double, 1>& sl_beam_szas() const;
  void sl_beam_szas(const blitz::Array<double, 1>& sl_beam_szas_in);
  
  const double& sl_chlorconc() const;
  void sl_chlorconc(const double& sl_chlorconc_in);
  
  const bool sl_do_doublet_geometry() const;
  void sl_do_doublet_geometry(const bool& sl_do_doublet_geometry_in);
  
  const bool sl_do_exact() const;
  void sl_do_exact(const bool& sl_do_exact_in);
  
  const bool sl_do_exactonly() const;
  void sl_do_exactonly(const bool& sl_do_exactonly_in);
  
  const bool sl_do_facetisotropy() const;
  void sl_do_facetisotropy(const bool& sl_do_facetisotropy_in);
  
  const bool sl_do_fluorescence() const;
  void sl_do_fluorescence(const bool& sl_do_fluorescence_in);
  
  const bool sl_do_foamoption() const;
  void sl_do_foamoption(const bool& sl_do_foamoption_in);
  
  const bool sl_do_fourier_output() const;
  void sl_do_fourier_output(const bool& sl_do_fourier_output_in);
  
  const bool sl_do_glintshadow() const;
  void sl_do_glintshadow(const bool& sl_do_glintshadow_in);
  
  const bool sl_do_isotropic() const;
  void sl_do_isotropic(const bool& sl_do_isotropic_in);
  
  const bool sl_do_roughsurface() const;
  void sl_do_roughsurface(const bool& sl_do_roughsurface_in);
  
  const bool sl_do_sleaving() const;
  void sl_do_sleaving(const bool& sl_do_sleaving_in);
  
  const bool sl_do_solar_sources() const;
  void sl_do_solar_sources(const bool& sl_do_solar_sources_in);
  
  const bool sl_do_user_obsgeoms() const;
  void sl_do_user_obsgeoms(const bool& sl_do_user_obsgeoms_in);
  
  const bool sl_do_user_streams() const;
  void sl_do_user_streams(const bool& sl_do_user_streams_in);
  
  const double& sl_fl_amplitude755() const;
  void sl_fl_amplitude755(const double& sl_fl_amplitude755_in);
  
  const bool sl_fl_do_datagaussian() const;
  void sl_fl_do_datagaussian(const bool& sl_fl_do_datagaussian_in);
  
  const blitz::Array<int, 1>& sl_fl_epoch() const;
  void sl_fl_epoch(const blitz::Array<int, 1>& sl_fl_epoch_in);
  
  const blitz::Array<double, 2>& sl_fl_inputgaussians() const;
  void sl_fl_inputgaussians(const blitz::Array<double, 2>& sl_fl_inputgaussians_in);
  
  const double& sl_fl_latitude() const;
  void sl_fl_latitude(const double& sl_fl_latitude_in);
  
  const double& sl_fl_longitude() const;
  void sl_fl_longitude(const double& sl_fl_longitude_in);
  
  const double& sl_fl_wavelength() const;
  void sl_fl_wavelength(const double& sl_fl_wavelength_in);
  
  const int& sl_n_user_doublets() const;
  void sl_n_user_doublets(const int& sl_n_user_doublets_in);
  
  const int& sl_n_user_obsgeoms() const;
  void sl_n_user_obsgeoms(const int& sl_n_user_obsgeoms_in);
  
  const int& sl_n_user_relazms() const;
  void sl_n_user_relazms(const int& sl_n_user_relazms_in);
  
  const int& sl_n_user_streams() const;
  void sl_n_user_streams(const int& sl_n_user_streams_in);
  
  const int& sl_nbeams() const;
  void sl_nbeams(const int& sl_nbeams_in);
  
  const int& sl_nstokes() const;
  void sl_nstokes(const int& sl_nstokes_in);
  
  const int& sl_nstreams() const;
  void sl_nstreams(const int& sl_nstreams_in);
  
  const double& sl_salinity() const;
  void sl_salinity(const double& sl_salinity_in);
  
  const blitz::Array<double, 1>& sl_user_angles_input() const;
  void sl_user_angles_input(const blitz::Array<double, 1>& sl_user_angles_input_in);
  
  const blitz::Array<double, 2>& sl_user_doublets() const;
  void sl_user_doublets(const blitz::Array<double, 2>& sl_user_doublets_in);
  
  const blitz::Array<double, 2>& sl_user_obsgeoms() const;
  void sl_user_obsgeoms(const blitz::Array<double, 2>& sl_user_obsgeoms_in);
  
  const blitz::Array<double, 1>& sl_user_relazms() const;
  void sl_user_relazms(const blitz::Array<double, 1>& sl_user_relazms_in);
  
  const std::string sl_vsleave_datapath() const;
  
  const double& sl_wavelength() const;
  void sl_wavelength(const double& sl_wavelength_in);
  
  const blitz::Array<double, 1>& sl_winddir() const;
  void sl_winddir(const blitz::Array<double, 1>& sl_winddir_in);
  
  const double& sl_windspeed() const;
  void sl_windspeed(const double& sl_windspeed_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Lincontrol : public Spurr_Fixed_Lincontrol_Base {
public:
  VLidort_Fixed_Lincontrol();
  VLidort_Fixed_Lincontrol(const VLidort_Fixed_Lincontrol& src);
  ~VLidort_Fixed_Lincontrol();

  const std::vector< std::string > ts_columnwf_names() const;
  
  const blitz::Array<bool, 1> ts_layer_vary_flag() const;
  void ts_layer_vary_flag(const blitz::Array<bool, 1>& ts_layer_vary_flag_in);
  
  const blitz::Array<int, 1>& ts_layer_vary_number() const;
  void ts_layer_vary_number(const blitz::Array<int, 1>& ts_layer_vary_number_in);
  
  const int& ts_n_sleave_wfs() const;
  void ts_n_sleave_wfs(const int& ts_n_sleave_wfs_in);
  
  const int& ts_n_surface_wfs() const;
  void ts_n_surface_wfs(const int& ts_n_surface_wfs_in);
  
  const int& ts_n_totalcolumn_wfs() const;
  void ts_n_totalcolumn_wfs(const int& ts_n_totalcolumn_wfs_in);
  
  const int& ts_n_totalprofile_wfs() const;
  void ts_n_totalprofile_wfs(const int& ts_n_totalprofile_wfs_in);
  
  const std::vector< std::string > ts_profilewf_names() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Linoptical : public Spurr_Fixed_Linoptical_Base {
public:
  VLidort_Fixed_Linoptical();
  VLidort_Fixed_Linoptical(const VLidort_Fixed_Linoptical& src);
  ~VLidort_Fixed_Linoptical();

  const blitz::Array<double, 2>& ts_l_deltau_vert_input() const;
  void ts_l_deltau_vert_input(const blitz::Array<double, 2>& ts_l_deltau_vert_input_in);
  
  const blitz::Array<double, 4>& ts_l_fmatrix_dn() const;
  void ts_l_fmatrix_dn(const blitz::Array<double, 4>& ts_l_fmatrix_dn_in);
  
  const blitz::Array<double, 4>& ts_l_fmatrix_up() const;
  void ts_l_fmatrix_up(const blitz::Array<double, 4>& ts_l_fmatrix_up_in);
  
  const blitz::Array<double, 4>& ts_l_greekmat_total_input() const;
  void ts_l_greekmat_total_input(const blitz::Array<double, 4>& ts_l_greekmat_total_input_in);
  
  const blitz::Array<double, 2>& ts_l_omega_total_input() const;
  void ts_l_omega_total_input(const blitz::Array<double, 2>& ts_l_omega_total_input_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Lininputs : public Spurr_Fixed_Lininputs_Base {
public:
  VLidort_Fixed_Lininputs();
  VLidort_Fixed_Lininputs(const VLidort_Fixed_Lininputs& src);
  ~VLidort_Fixed_Lininputs();

  const VLidort_Fixed_Lincontrol& cont() const;
  const Spurr_Fixed_Lincontrol_Base& fixed_lincontrol_base() const;
  
  const VLidort_Fixed_Linoptical& optical() const;
  const Spurr_Fixed_Linoptical_Base& fixed_linoptical_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Lincontrol : public Spurr_Modified_Lincontrol_Base {
public:
  VLidort_Modified_Lincontrol();
  VLidort_Modified_Lincontrol(const VLidort_Modified_Lincontrol& src);
  ~VLidort_Modified_Lincontrol();

  const bool ts_do_atmos_lbbf() const;
  void ts_do_atmos_lbbf(const bool& ts_do_atmos_lbbf_in);
  
  const bool ts_do_atmos_linearization() const;
  void ts_do_atmos_linearization(const bool& ts_do_atmos_linearization_in);
  
  const bool ts_do_column_linearization() const;
  void ts_do_column_linearization(const bool& ts_do_column_linearization_in);
  
  const bool ts_do_linearization() const;
  void ts_do_linearization(const bool& ts_do_linearization_in);
  
  const bool ts_do_profile_linearization() const;
  void ts_do_profile_linearization(const bool& ts_do_profile_linearization_in);
  
  const bool ts_do_simulation_only() const;
  void ts_do_simulation_only(const bool& ts_do_simulation_only_in);
  
  const bool ts_do_sleave_wfs() const;
  void ts_do_sleave_wfs(const bool& ts_do_sleave_wfs_in);
  
  const bool ts_do_surface_lbbf() const;
  void ts_do_surface_lbbf(const bool& ts_do_surface_lbbf_in);
  
  const bool ts_do_surface_linearization() const;
  void ts_do_surface_linearization(const bool& ts_do_surface_linearization_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Lininputs : public Spurr_Modified_Lininputs_Base {
public:
  VLidort_Modified_Lininputs();
  VLidort_Modified_Lininputs(const VLidort_Modified_Lininputs& src);
  ~VLidort_Modified_Lininputs();

  const VLidort_Modified_Lincontrol& mcont() const;
  const Spurr_Modified_Lincontrol_Base& modified_lincontrol_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Lincol : public Spurr_Type_Structure {
public:
  VLidort_Lincol();
  VLidort_Lincol(const VLidort_Lincol& src);
  ~VLidort_Lincol();

  const blitz::Array<double, 3>& ts_albmed_fluxes_colwf() const;
  void ts_albmed_fluxes_colwf(const blitz::Array<double, 3>& ts_albmed_fluxes_colwf_in);
  
  const blitz::Array<double, 3>& ts_albmed_user_colwf() const;
  void ts_albmed_user_colwf(const blitz::Array<double, 3>& ts_albmed_user_colwf_in);
  
  const blitz::Array<double, 5>& ts_columnwf() const;
  void ts_columnwf(const blitz::Array<double, 5>& ts_columnwf_in);
  
  const blitz::Array<double, 4>& ts_dnflux_direct_colwf() const;
  void ts_dnflux_direct_colwf(const blitz::Array<double, 4>& ts_dnflux_direct_colwf_in);
  
  const blitz::Array<double, 4>& ts_dnmeanst_direct_colwf() const;
  void ts_dnmeanst_direct_colwf(const blitz::Array<double, 4>& ts_dnmeanst_direct_colwf_in);
  
  const blitz::Array<double, 5>& ts_flux_diffuse_colwf() const;
  void ts_flux_diffuse_colwf(const blitz::Array<double, 5>& ts_flux_diffuse_colwf_in);
  
  const blitz::Array<double, 4>& ts_lc_layer_mssts() const;
  void ts_lc_layer_mssts(const blitz::Array<double, 4>& ts_lc_layer_mssts_in);
  
  const blitz::Array<double, 3>& ts_lc_lostrans() const;
  void ts_lc_lostrans(const blitz::Array<double, 3>& ts_lc_lostrans_in);
  
  const blitz::Array<double, 3>& ts_lc_surf_mssts() const;
  void ts_lc_surf_mssts(const blitz::Array<double, 3>& ts_lc_surf_mssts_in);
  
  const blitz::Array<double, 5>& ts_meanst_diffuse_colwf() const;
  void ts_meanst_diffuse_colwf(const blitz::Array<double, 5>& ts_meanst_diffuse_colwf_in);
  
  const blitz::Array<double, 1>& ts_planetary_sbterm_colwf() const;
  void ts_planetary_sbterm_colwf(const blitz::Array<double, 1>& ts_planetary_sbterm_colwf_in);
  
  const blitz::Array<double, 3>& ts_planetary_transterm_colwf() const;
  void ts_planetary_transterm_colwf(const blitz::Array<double, 3>& ts_planetary_transterm_colwf_in);
  
  const blitz::Array<double, 3>& ts_transbeam_colwf() const;
  void ts_transbeam_colwf(const blitz::Array<double, 3>& ts_transbeam_colwf_in);
  
  const blitz::Array<double, 3>& ts_trnmed_fluxes_colwf() const;
  void ts_trnmed_fluxes_colwf(const blitz::Array<double, 3>& ts_trnmed_fluxes_colwf_in);
  
  const blitz::Array<double, 3>& ts_trnmed_user_colwf() const;
  void ts_trnmed_user_colwf(const blitz::Array<double, 3>& ts_trnmed_user_colwf_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linprof : public Spurr_Type_Structure {
public:
  VLidort_Linprof();
  VLidort_Linprof(const VLidort_Linprof& src);
  ~VLidort_Linprof();

  const blitz::Array<double, 4>& ts_albmed_fluxes_profwf() const;
  void ts_albmed_fluxes_profwf(const blitz::Array<double, 4>& ts_albmed_fluxes_profwf_in);
  
  const blitz::Array<double, 4>& ts_albmed_user_profwf() const;
  void ts_albmed_user_profwf(const blitz::Array<double, 4>& ts_albmed_user_profwf_in);
  
  const blitz::Array<double, 5>& ts_dnflux_direct_profwf() const;
  void ts_dnflux_direct_profwf(const blitz::Array<double, 5>& ts_dnflux_direct_profwf_in);
  
  const blitz::Array<double, 5>& ts_dnmeanst_direct_profwf() const;
  void ts_dnmeanst_direct_profwf(const blitz::Array<double, 5>& ts_dnmeanst_direct_profwf_in);
  
  const blitz::Array<double, 6>& ts_flux_diffuse_profwf() const;
  void ts_flux_diffuse_profwf(const blitz::Array<double, 6>& ts_flux_diffuse_profwf_in);
  
  const blitz::Array<double, 5>& ts_lp_layer_mssts() const;
  void ts_lp_layer_mssts(const blitz::Array<double, 5>& ts_lp_layer_mssts_in);
  
  const blitz::Array<double, 3>& ts_lp_lostrans() const;
  void ts_lp_lostrans(const blitz::Array<double, 3>& ts_lp_lostrans_in);
  
  const blitz::Array<double, 4>& ts_lp_surf_mssts() const;
  void ts_lp_surf_mssts(const blitz::Array<double, 4>& ts_lp_surf_mssts_in);
  
  const blitz::Array<double, 6>& ts_meanst_diffuse_profwf() const;
  void ts_meanst_diffuse_profwf(const blitz::Array<double, 6>& ts_meanst_diffuse_profwf_in);
  
  const blitz::Array<double, 2>& ts_planetary_sbterm_profwf() const;
  void ts_planetary_sbterm_profwf(const blitz::Array<double, 2>& ts_planetary_sbterm_profwf_in);
  
  const blitz::Array<double, 4>& ts_planetary_transterm_profwf() const;
  void ts_planetary_transterm_profwf(const blitz::Array<double, 4>& ts_planetary_transterm_profwf_in);
  
  const blitz::Array<double, 6>& ts_profilewf() const;
  void ts_profilewf(const blitz::Array<double, 6>& ts_profilewf_in);
  
  const blitz::Array<double, 4>& ts_transbeam_profwf() const;
  void ts_transbeam_profwf(const blitz::Array<double, 4>& ts_transbeam_profwf_in);
  
  const blitz::Array<double, 4>& ts_trnmed_fluxes_profwf() const;
  void ts_trnmed_fluxes_profwf(const blitz::Array<double, 4>& ts_trnmed_fluxes_profwf_in);
  
  const blitz::Array<double, 4>& ts_trnmed_user_profwf() const;
  void ts_trnmed_user_profwf(const blitz::Array<double, 4>& ts_trnmed_user_profwf_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linatmos : public Spurr_Linatmos_Base {
public:
  VLidort_Linatmos();
  VLidort_Linatmos(const VLidort_Linatmos& src);
  ~VLidort_Linatmos();

  const blitz::Array<double, 5>& ts_abbwfs_fluxes() const;
  void ts_abbwfs_fluxes(const blitz::Array<double, 5>& ts_abbwfs_fluxes_in);
  
  const blitz::Array<double, 5>& ts_abbwfs_jacobians() const;
  void ts_abbwfs_jacobians(const blitz::Array<double, 5>& ts_abbwfs_jacobians_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsurf : public Spurr_Linsurf_Base {
public:
  VLidort_Linsurf();
  VLidort_Linsurf(const VLidort_Linsurf& src);
  ~VLidort_Linsurf();

  const blitz::Array<double, 5>& ts_flux_diffuse_surfwf() const;
  void ts_flux_diffuse_surfwf(const blitz::Array<double, 5>& ts_flux_diffuse_surfwf_in);
  
  const blitz::Array<double, 4>& ts_ls_layer_mssts() const;
  void ts_ls_layer_mssts(const blitz::Array<double, 4>& ts_ls_layer_mssts_in);
  
  const blitz::Array<double, 3>& ts_ls_surf_mssts() const;
  void ts_ls_surf_mssts(const blitz::Array<double, 3>& ts_ls_surf_mssts_in);
  
  const blitz::Array<double, 5>& ts_meanst_diffuse_surfwf() const;
  void ts_meanst_diffuse_surfwf(const blitz::Array<double, 5>& ts_meanst_diffuse_surfwf_in);
  
  const blitz::Array<double, 4>& ts_sbbwfs_fluxes() const;
  void ts_sbbwfs_fluxes(const blitz::Array<double, 4>& ts_sbbwfs_fluxes_in);
  
  const blitz::Array<double, 4>& ts_sbbwfs_jacobians() const;
  void ts_sbbwfs_jacobians(const blitz::Array<double, 4>& ts_sbbwfs_jacobians_in);
  
  const blitz::Array<double, 5>& ts_surfacewf() const;
  void ts_surfacewf(const blitz::Array<double, 5>& ts_surfacewf_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linoutputs : public Spurr_Linoutputs_Base {
public:
  VLidort_Linoutputs();
  VLidort_Linoutputs(const VLidort_Linoutputs& src);
  ~VLidort_Linoutputs();

  const VLidort_Linatmos& atmos() const;
  const Spurr_Linatmos_Base& linatmos_base() const;
  
  const VLidort_Lincol& col() const;
  
  
  const VLidort_Linprof& prof() const;
  
  
  const VLidort_Linsurf& surf() const;
  const Spurr_Linsurf_Base& linsurf_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Brdf : public Spurr_Lin_Sup_Brdf_Base {
public:
  VLidort_Linsup_Brdf();
  VLidort_Linsup_Brdf(const VLidort_Linsup_Brdf& src);
  ~VLidort_Linsup_Brdf();

  const blitz::Array<double, 5>& ts_ls_brdf_f() const;
  void ts_ls_brdf_f(const blitz::Array<double, 5>& ts_ls_brdf_f_in);
  
  const blitz::Array<double, 5>& ts_ls_brdf_f_0() const;
  void ts_ls_brdf_f_0(const blitz::Array<double, 5>& ts_ls_brdf_f_0_in);
  
  const blitz::Array<double, 3>& ts_ls_emissivity() const;
  void ts_ls_emissivity(const blitz::Array<double, 3>& ts_ls_emissivity_in);
  
  const blitz::Array<double, 5>& ts_ls_exactdb_brdfunc() const;
  void ts_ls_exactdb_brdfunc(const blitz::Array<double, 5>& ts_ls_exactdb_brdfunc_in);
  
  const blitz::Array<double, 5>& ts_ls_user_brdf_f() const;
  void ts_ls_user_brdf_f(const blitz::Array<double, 5>& ts_ls_user_brdf_f_in);
  
  const blitz::Array<double, 5>& ts_ls_user_brdf_f_0() const;
  void ts_ls_user_brdf_f_0(const blitz::Array<double, 5>& ts_ls_user_brdf_f_0_in);
  
  const blitz::Array<double, 3>& ts_ls_user_emissivity() const;
  void ts_ls_user_emissivity(const blitz::Array<double, 3>& ts_ls_user_emissivity_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Sleave : public Spurr_Lin_Sup_Sleave_Base {
public:
  VLidort_Linsup_Sleave();
  VLidort_Linsup_Sleave(const VLidort_Linsup_Sleave& src);
  ~VLidort_Linsup_Sleave();

  const blitz::Array<double, 5>& ts_lssl_slterm_f_0() const;
  void ts_lssl_slterm_f_0(const blitz::Array<double, 5>& ts_lssl_slterm_f_0_in);
  
  const blitz::Array<double, 3>& ts_lssl_slterm_isotropic() const;
  void ts_lssl_slterm_isotropic(const blitz::Array<double, 3>& ts_lssl_slterm_isotropic_in);
  
  const blitz::Array<double, 5>& ts_lssl_slterm_userangles() const;
  void ts_lssl_slterm_userangles(const blitz::Array<double, 5>& ts_lssl_slterm_userangles_in);
  
  const blitz::Array<double, 5>& ts_lssl_user_slterm_f_0() const;
  void ts_lssl_user_slterm_f_0(const blitz::Array<double, 5>& ts_lssl_user_slterm_f_0_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Ss_Col : public Spurr_Type_Structure {
public:
  VLidort_Linsup_Ss_Col();
  VLidort_Linsup_Ss_Col(const VLidort_Linsup_Ss_Col& src);
  ~VLidort_Linsup_Ss_Col();

  const blitz::Array<double, 4>& ts_columnwf_db() const;
  void ts_columnwf_db(const blitz::Array<double, 4>& ts_columnwf_db_in);
  
  const blitz::Array<double, 5>& ts_columnwf_ss() const;
  void ts_columnwf_ss(const blitz::Array<double, 5>& ts_columnwf_ss_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Ss_Prof : public Spurr_Type_Structure {
public:
  VLidort_Linsup_Ss_Prof();
  VLidort_Linsup_Ss_Prof(const VLidort_Linsup_Ss_Prof& src);
  ~VLidort_Linsup_Ss_Prof();

  const blitz::Array<double, 5>& ts_profilewf_db() const;
  void ts_profilewf_db(const blitz::Array<double, 5>& ts_profilewf_db_in);
  
  const blitz::Array<double, 6>& ts_profilewf_ss() const;
  void ts_profilewf_ss(const blitz::Array<double, 6>& ts_profilewf_ss_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Ss_Surf : public Spurr_Lin_Sup_Ss_Surf_Base {
public:
  VLidort_Linsup_Ss_Surf();
  VLidort_Linsup_Ss_Surf(const VLidort_Linsup_Ss_Surf& src);
  ~VLidort_Linsup_Ss_Surf();

  const blitz::Array<double, 4>& ts_surfacewf_db() const;
  void ts_surfacewf_db(const blitz::Array<double, 4>& ts_surfacewf_db_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Ss : public Spurr_Lin_Sup_Ss_Base {
public:
  VLidort_Linsup_Ss();
  VLidort_Linsup_Ss(const VLidort_Linsup_Ss& src);
  ~VLidort_Linsup_Ss();

  const VLidort_Linsup_Ss_Col& col() const;
  
  
  const VLidort_Linsup_Ss_Prof& prof() const;
  
  
  const VLidort_Linsup_Ss_Surf& surf() const;
  const Spurr_Lin_Sup_Ss_Surf_Base& linsup_ss_surf_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Linsup_Inout : public Spurr_Lin_Sup_Inout_Base {
public:
  VLidort_Linsup_Inout();
  VLidort_Linsup_Inout(const VLidort_Linsup_Inout& src);
  ~VLidort_Linsup_Inout();

  const VLidort_Linsup_Brdf& brdf() const;
  const Spurr_Lin_Sup_Brdf_Base& linsup_brdf_base() const;
  
  const VLidort_Linsup_Sleave& sleave() const;
  const Spurr_Lin_Sup_Sleave_Base& linsup_sleave_base() const;
  
  const VLidort_Linsup_Ss& ss() const;
  const Spurr_Lin_Sup_Ss_Base& linsup_ss_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Main_Outputs : public Spurr_Main_Outputs_Base {
public:
  VLidort_Main_Outputs();
  VLidort_Main_Outputs(const VLidort_Main_Outputs& src);
  ~VLidort_Main_Outputs();

  const blitz::Array<double, 2>& ts_albmed_fluxes() const;
  void ts_albmed_fluxes(const blitz::Array<double, 2>& ts_albmed_fluxes_in);
  
  const blitz::Array<double, 2>& ts_albmed_user() const;
  void ts_albmed_user(const blitz::Array<double, 2>& ts_albmed_user_in);
  
  const blitz::Array<double, 3>& ts_contribs() const;
  void ts_contribs(const blitz::Array<double, 3>& ts_contribs_in);
  
  const blitz::Array<double, 3>& ts_dnflux_direct() const;
  void ts_dnflux_direct(const blitz::Array<double, 3>& ts_dnflux_direct_in);
  
  const blitz::Array<double, 3>& ts_dnmeanst_direct() const;
  void ts_dnmeanst_direct(const blitz::Array<double, 3>& ts_dnmeanst_direct_in);
  
  const blitz::Array<double, 4>& ts_flux_diffuse() const;
  void ts_flux_diffuse(const blitz::Array<double, 4>& ts_flux_diffuse_in);
  
  const blitz::Array<int, 1>& ts_fourier_saved() const;
  void ts_fourier_saved(const blitz::Array<int, 1>& ts_fourier_saved_in);
  
  const blitz::Array<double, 3>& ts_layer_mssts() const;
  void ts_layer_mssts(const blitz::Array<double, 3>& ts_layer_mssts_in);
  
  const blitz::Array<double, 2>& ts_lostrans() const;
  void ts_lostrans(const blitz::Array<double, 2>& ts_lostrans_in);
  
  const blitz::Array<double, 4>& ts_meanst_diffuse() const;
  void ts_meanst_diffuse(const blitz::Array<double, 4>& ts_meanst_diffuse_in);
  
  const int& ts_n_geometries() const;
  void ts_n_geometries(const int& ts_n_geometries_in);
  
  const blitz::Array<double, 2>& ts_pathgeoms() const;
  void ts_pathgeoms(const blitz::Array<double, 2>& ts_pathgeoms_in);
  
  const double& ts_planetary_sbterm() const;
  void ts_planetary_sbterm(const double& ts_planetary_sbterm_in);
  
  const blitz::Array<double, 2>& ts_planetary_transterm() const;
  void ts_planetary_transterm(const blitz::Array<double, 2>& ts_planetary_transterm_in);
  
  const blitz::Array<double, 1>& ts_solarbeam_boatrans() const;
  void ts_solarbeam_boatrans(const blitz::Array<double, 1>& ts_solarbeam_boatrans_in);
  
  const blitz::Array<double, 4>& ts_stokes() const;
  void ts_stokes(const blitz::Array<double, 4>& ts_stokes_in);
  
  const blitz::Array<double, 2>& ts_surf_mssts() const;
  void ts_surf_mssts(const blitz::Array<double, 2>& ts_surf_mssts_in);
  
  const blitz::Array<int, 1>& ts_sza_offsets() const;
  void ts_sza_offsets(const blitz::Array<int, 1>& ts_sza_offsets_in);
  
  const blitz::Array<int, 1>& ts_szd_offsets() const;
  void ts_szd_offsets(const blitz::Array<int, 1>& ts_szd_offsets_in);
  
  const blitz::Array<double, 2>& ts_trnmed_fluxes() const;
  void ts_trnmed_fluxes(const blitz::Array<double, 2>& ts_trnmed_fluxes_in);
  
  const blitz::Array<double, 2>& ts_trnmed_user() const;
  void ts_trnmed_user(const blitz::Array<double, 2>& ts_trnmed_user_in);
  
  const blitz::Array<int, 2>& ts_vza_offsets() const;
  void ts_vza_offsets(const blitz::Array<int, 2>& ts_vza_offsets_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Wladjusted_Outputs : public Spurr_Wladjusted_Outputs_Base {
public:
  VLidort_Wladjusted_Outputs();
  VLidort_Wladjusted_Outputs(const VLidort_Wladjusted_Outputs& src);
  ~VLidort_Wladjusted_Outputs();

  const blitz::Array<double, 4>& ts_wladjusted_direct() const;
  void ts_wladjusted_direct(const blitz::Array<double, 4>& ts_wladjusted_direct_in);
  
  const blitz::Array<double, 4>& ts_wladjusted_f_ords_0() const;
  void ts_wladjusted_f_ords_0(const blitz::Array<double, 4>& ts_wladjusted_f_ords_0_in);
  
  const blitz::Array<double, 4>& ts_wladjusted_f_user_0() const;
  void ts_wladjusted_f_user_0(const blitz::Array<double, 4>& ts_wladjusted_f_user_0_in);
  
  const blitz::Array<double, 2>& ts_wladjusted_isotropic() const;
  void ts_wladjusted_isotropic(const blitz::Array<double, 2>& ts_wladjusted_isotropic_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Exception_Handling : public Spurr_Exception_Handling_Base {
public:
  VLidort_Exception_Handling();
  VLidort_Exception_Handling(const VLidort_Exception_Handling& src);
  ~VLidort_Exception_Handling();

  const std::vector< std::string > ts_actions() const;
  
  const std::vector< std::string > ts_checkmessages() const;
  
  const std::string ts_message() const;
  
  const int& ts_ncheckmessages() const;
  void ts_ncheckmessages(const int& ts_ncheckmessages_in);
  
  const int& ts_status_calculation() const;
  void ts_status_calculation(const int& ts_status_calculation_in);
  
  const int& ts_status_inputcheck() const;
  void ts_status_inputcheck(const int& ts_status_inputcheck_in);
  
  const std::string ts_trace_1() const;
  
  const std::string ts_trace_2() const;
  
  const std::string ts_trace_3() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Input_Exception_Handling : public Spurr_Input_Exception_Handling_Base {
public:
  VLidort_Input_Exception_Handling();
  VLidort_Input_Exception_Handling(const VLidort_Input_Exception_Handling& src);
  ~VLidort_Input_Exception_Handling();

  const std::vector< std::string > ts_inputactions() const;
  
  const std::vector< std::string > ts_inputmessages() const;
  
  const int& ts_ninputmessages() const;
  void ts_ninputmessages(const int& ts_ninputmessages_in);
  
  const int& ts_status_inputread() const;
  void ts_status_inputread(const int& ts_status_inputread_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Outputs : public Spurr_Outputs_Base {
public:
  VLidort_Outputs();
  VLidort_Outputs(const VLidort_Outputs& src);
  ~VLidort_Outputs();

  const VLidort_Main_Outputs& main() const;
  const Spurr_Main_Outputs_Base& main_outputs_base() const;
  
  const VLidort_Exception_Handling& status() const;
  const Spurr_Exception_Handling_Base& exception_handling_base() const;
  
  const VLidort_Wladjusted_Outputs& wlout() const;
  const Spurr_Wladjusted_Outputs_Base& wladjusted_outputs_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Sup_Brdf : public Spurr_Sup_Brdf_Base {
public:
  VLidort_Sup_Brdf();
  VLidort_Sup_Brdf(const VLidort_Sup_Brdf& src);
  ~VLidort_Sup_Brdf();

  const blitz::Array<double, 4>& ts_brdf_f() const;
  void ts_brdf_f(const blitz::Array<double, 4>& ts_brdf_f_in);
  
  const blitz::Array<double, 4>& ts_brdf_f_0() const;
  void ts_brdf_f_0(const blitz::Array<double, 4>& ts_brdf_f_0_in);
  
  const blitz::Array<double, 2>& ts_emissivity() const;
  void ts_emissivity(const blitz::Array<double, 2>& ts_emissivity_in);
  
  const blitz::Array<double, 4>& ts_exactdb_brdfunc() const;
  void ts_exactdb_brdfunc(const blitz::Array<double, 4>& ts_exactdb_brdfunc_in);
  
  const blitz::Array<double, 4>& ts_user_brdf_f() const;
  void ts_user_brdf_f(const blitz::Array<double, 4>& ts_user_brdf_f_in);
  
  const blitz::Array<double, 4>& ts_user_brdf_f_0() const;
  void ts_user_brdf_f_0(const blitz::Array<double, 4>& ts_user_brdf_f_0_in);
  
  const blitz::Array<double, 2>& ts_user_emissivity() const;
  void ts_user_emissivity(const blitz::Array<double, 2>& ts_user_emissivity_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Sup_Sleave : public Spurr_Sup_Sleave_Base {
public:
  VLidort_Sup_Sleave();
  VLidort_Sup_Sleave(const VLidort_Sup_Sleave& src);
  ~VLidort_Sup_Sleave();

  const blitz::Array<double, 4>& ts_slterm_f_0() const;
  void ts_slterm_f_0(const blitz::Array<double, 4>& ts_slterm_f_0_in);
  
  const blitz::Array<double, 2>& ts_slterm_isotropic() const;
  void ts_slterm_isotropic(const blitz::Array<double, 2>& ts_slterm_isotropic_in);
  
  const blitz::Array<double, 4>& ts_slterm_userangles() const;
  void ts_slterm_userangles(const blitz::Array<double, 4>& ts_slterm_userangles_in);
  
  const blitz::Array<double, 4>& ts_user_slterm_f_0() const;
  void ts_user_slterm_f_0(const blitz::Array<double, 4>& ts_user_slterm_f_0_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Sup_Ss : public Spurr_Sup_Ss_Base {
public:
  VLidort_Sup_Ss();
  VLidort_Sup_Ss(const VLidort_Sup_Ss& src);
  ~VLidort_Sup_Ss();

  const blitz::Array<double, 3>& ts_contribs_ss() const;
  void ts_contribs_ss(const blitz::Array<double, 3>& ts_contribs_ss_in);
  
  const blitz::Array<double, 3>& ts_stokes_db() const;
  void ts_stokes_db(const blitz::Array<double, 3>& ts_stokes_db_in);
  
  const blitz::Array<double, 4>& ts_stokes_ss() const;
  void ts_stokes_ss(const blitz::Array<double, 4>& ts_stokes_ss_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Sup_Inout : public Spurr_Sup_Inout_Base {
public:
  VLidort_Sup_Inout();
  VLidort_Sup_Inout(const VLidort_Sup_Inout& src);
  ~VLidort_Sup_Inout();

  const VLidort_Sup_Brdf& brdf() const;
  const Spurr_Sup_Brdf_Base& sup_brdf_base() const;
  
  const VLidort_Sup_Sleave& sleave() const;
  const Spurr_Sup_Sleave_Base& sup_sleave_base() const;
  
  const VLidort_Sup_Ss& ss() const;
  const Spurr_Sup_Ss_Base& sup_ss_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Boolean : public Spurr_Fixed_Boolean_Base {
public:
  VLidort_Fixed_Boolean();
  VLidort_Fixed_Boolean(const VLidort_Fixed_Boolean& src);
  ~VLidort_Fixed_Boolean();

  const blitz::Array<bool, 1> ts_do_albtrn_media() const;
  void ts_do_albtrn_media(const blitz::Array<bool, 1>& ts_do_albtrn_media_in);
  
  const bool ts_do_boa_illumination() const;
  void ts_do_boa_illumination(const bool& ts_do_boa_illumination_in);
  
  const bool ts_do_dnwelling() const;
  void ts_do_dnwelling(const bool& ts_do_dnwelling_in);
  
  const bool ts_do_fluorescence() const;
  void ts_do_fluorescence(const bool& ts_do_fluorescence_in);
  
  const bool ts_do_fourier0_nstokes2() const;
  void ts_do_fourier0_nstokes2(const bool& ts_do_fourier0_nstokes2_in);
  
  const bool ts_do_fullrad_mode() const;
  void ts_do_fullrad_mode(const bool& ts_do_fullrad_mode_in);
  
  const bool ts_do_lambertian_surface() const;
  void ts_do_lambertian_surface(const bool& ts_do_lambertian_surface_in);
  
  const bool ts_do_mssts() const;
  void ts_do_mssts(const bool& ts_do_mssts_in);
  
  const bool ts_do_plane_parallel() const;
  void ts_do_plane_parallel(const bool& ts_do_plane_parallel_in);
  
  const bool ts_do_planetary_problem() const;
  void ts_do_planetary_problem(const bool& ts_do_planetary_problem_in);
  
  const bool ts_do_sl_isotropic() const;
  void ts_do_sl_isotropic(const bool& ts_do_sl_isotropic_in);
  
  const bool ts_do_specialist_option_1() const;
  void ts_do_specialist_option_1(const bool& ts_do_specialist_option_1_in);
  
  const bool ts_do_specialist_option_2() const;
  void ts_do_specialist_option_2(const bool& ts_do_specialist_option_2_in);
  
  const bool ts_do_specialist_option_3() const;
  void ts_do_specialist_option_3(const bool& ts_do_specialist_option_3_in);
  
  const bool ts_do_surface_emission() const;
  void ts_do_surface_emission(const bool& ts_do_surface_emission_in);
  
  const bool ts_do_surface_leaving() const;
  void ts_do_surface_leaving(const bool& ts_do_surface_leaving_in);
  
  const bool ts_do_tf_iteration() const;
  void ts_do_tf_iteration(const bool& ts_do_tf_iteration_in);
  
  const bool ts_do_thermal_emission() const;
  void ts_do_thermal_emission(const bool& ts_do_thermal_emission_in);
  
  const bool ts_do_toa_contribs() const;
  void ts_do_toa_contribs(const bool& ts_do_toa_contribs_in);
  
  const bool ts_do_toa_illumination() const;
  void ts_do_toa_illumination(const bool& ts_do_toa_illumination_in);
  
  const bool ts_do_upwelling() const;
  void ts_do_upwelling(const bool& ts_do_upwelling_in);
  
  const bool ts_do_water_leaving() const;
  void ts_do_water_leaving(const bool& ts_do_water_leaving_in);
  
  const bool ts_do_wladjusted_output() const;
  void ts_do_wladjusted_output(const bool& ts_do_wladjusted_output_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Control : public Spurr_Fixed_Control_Base {
public:
  VLidort_Fixed_Control();
  VLidort_Fixed_Control(const VLidort_Fixed_Control& src);
  ~VLidort_Fixed_Control();

  const double& ts_asymtx_tolerance() const;
  void ts_asymtx_tolerance(const double& ts_asymtx_tolerance_in);
  
  const double& ts_boa_illumination() const;
  void ts_boa_illumination(const double& ts_boa_illumination_in);
  
  const int& ts_n_thermal_coeffs() const;
  void ts_n_thermal_coeffs(const int& ts_n_thermal_coeffs_in);
  
  const int& ts_nfinelayers() const;
  void ts_nfinelayers(const int& ts_nfinelayers_in);
  
  const int& ts_nlayers() const;
  void ts_nlayers(const int& ts_nlayers_in);
  
  const int& ts_nlayers_cutoff() const;
  void ts_nlayers_cutoff(const int& ts_nlayers_cutoff_in);
  
  const int& ts_nlayers_noms() const;
  void ts_nlayers_noms(const int& ts_nlayers_noms_in);
  
  const int& ts_nstokes() const;
  void ts_nstokes(const int& ts_nstokes_in);
  
  const int& ts_nstreams() const;
  void ts_nstreams(const int& ts_nstreams_in);
  
  const int& ts_taylor_order() const;
  void ts_taylor_order(const int& ts_taylor_order_in);
  
  const double& ts_tf_criterion() const;
  void ts_tf_criterion(const double& ts_tf_criterion_in);
  
  const int& ts_tf_maxiter() const;
  void ts_tf_maxiter(const int& ts_tf_maxiter_in);
  
  const double& ts_toa_illumination() const;
  void ts_toa_illumination(const double& ts_toa_illumination_in);
  
  const double& ts_fourier_accuracy() const;
  void ts_fourier_accuracy(const double& ts_vlidort_accuracy_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Sunrays : public Spurr_Fixed_Sunrays_Base {
public:
  VLidort_Fixed_Sunrays();
  VLidort_Fixed_Sunrays(const VLidort_Fixed_Sunrays& src);
  ~VLidort_Fixed_Sunrays();

  const double& ts_flux_factor() const;
  void ts_flux_factor(const double& ts_flux_factor_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Uservalues : public Spurr_Fixed_Uservalues_Base {
public:
  VLidort_Fixed_Uservalues();
  VLidort_Fixed_Uservalues(const VLidort_Fixed_Uservalues& src);
  ~VLidort_Fixed_Uservalues();

  const int& ts_n_user_levels() const;
  void ts_n_user_levels(const int& ts_n_user_levels_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Chapman : public Spurr_Fixed_Chapman_Base {
public:
  VLidort_Fixed_Chapman();
  VLidort_Fixed_Chapman(const VLidort_Fixed_Chapman& src);
  ~VLidort_Fixed_Chapman();

  const blitz::Array<int, 1>& ts_finegrid() const;
  void ts_finegrid(const blitz::Array<int, 1>& ts_finegrid_in);
  
  const blitz::Array<double, 1>& ts_height_grid() const;
  void ts_height_grid(const blitz::Array<double, 1>& ts_height_grid_in);
  
  const blitz::Array<double, 1>& ts_pressure_grid() const;
  void ts_pressure_grid(const blitz::Array<double, 1>& ts_pressure_grid_in);
  
  const double& ts_rfindex_parameter() const;
  void ts_rfindex_parameter(const double& ts_rfindex_parameter_in);
  
  const blitz::Array<double, 1>& ts_temperature_grid() const;
  void ts_temperature_grid(const blitz::Array<double, 1>& ts_temperature_grid_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Optical : public Spurr_Fixed_Optical_Base {
public:
  VLidort_Fixed_Optical();
  VLidort_Fixed_Optical(const VLidort_Fixed_Optical& src);
  ~VLidort_Fixed_Optical();

  const double& ts_atmos_wavelength() const;
  void ts_atmos_wavelength(const double& ts_atmos_wavelength_in);
  
  const blitz::Array<double, 1>& ts_deltau_vert_input() const;
  void ts_deltau_vert_input(const blitz::Array<double, 1>& ts_deltau_vert_input_in);
  
  const blitz::Array<double, 3>& ts_fmatrix_dn() const;
  void ts_fmatrix_dn(const blitz::Array<double, 3>& ts_fmatrix_dn_in);
  
  const blitz::Array<double, 3>& ts_fmatrix_up() const;
  void ts_fmatrix_up(const blitz::Array<double, 3>& ts_fmatrix_up_in);
  
  const blitz::Array<double, 3>& ts_greekmat_total_input() const;
  void ts_greekmat_total_input(const blitz::Array<double, 3>& ts_greekmat_total_input_in);
  
  const double& ts_lambertian_albedo() const;
  void ts_lambertian_albedo(const double& ts_lambertian_albedo_in);
  
  const double& ts_surface_bb_input() const;
  void ts_surface_bb_input(const double& ts_surface_bb_input_in);
  
  const blitz::Array<double, 1>& ts_thermal_bb_input() const;
  void ts_thermal_bb_input(const blitz::Array<double, 1>& ts_thermal_bb_input_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Write : public Spurr_Fixed_Write_Base {
public:
  VLidort_Fixed_Write();
  VLidort_Fixed_Write(const VLidort_Fixed_Write& src);
  ~VLidort_Fixed_Write();

  const bool ts_do_debug_write() const;
  void ts_do_debug_write(const bool& ts_do_debug_write_in);
  
  const bool ts_do_write_fourier() const;
  void ts_do_write_fourier(const bool& ts_do_write_fourier_in);
  
  const bool ts_do_write_input() const;
  void ts_do_write_input(const bool& ts_do_write_input_in);
  
  const bool ts_do_write_results() const;
  void ts_do_write_results(const bool& ts_do_write_results_in);
  
  const bool ts_do_write_scenario() const;
  void ts_do_write_scenario(const bool& ts_do_write_scenario_in);
  
  const std::string ts_fourier_write_filename() const;
  
  const std::string ts_input_write_filename() const;
  
  const std::string ts_results_write_filename() const;
  
  const std::string ts_scenario_write_filename() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Fixed_Inputs : public Spurr_Fixed_Inputs_Base {
public:
  VLidort_Fixed_Inputs();
  VLidort_Fixed_Inputs(const VLidort_Fixed_Inputs& src);
  ~VLidort_Fixed_Inputs();

  const VLidort_Fixed_Boolean& f_bool() const;
  const Spurr_Fixed_Boolean_Base& fixed_boolean_base() const;
  
  const VLidort_Fixed_Chapman& chapman() const;
  const Spurr_Fixed_Chapman_Base& fixed_chapman_base() const;
  
  const VLidort_Fixed_Control& cont() const;
  const Spurr_Fixed_Control_Base& fixed_control_base() const;
  
  const VLidort_Fixed_Optical& optical() const;
  const Spurr_Fixed_Optical_Base& fixed_optical_base() const;
  
  const VLidort_Fixed_Sunrays& sunrays() const;
  const Spurr_Fixed_Sunrays_Base& fixed_sunrays_base() const;
  
  const VLidort_Fixed_Uservalues& userval() const;
  const Spurr_Fixed_Uservalues_Base& fixed_uservalues_base() const;
  
  const VLidort_Fixed_Write& write() const;
  const Spurr_Fixed_Write_Base& fixed_write_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Boolean : public Spurr_Modified_Boolean_Base {
public:
  VLidort_Modified_Boolean();
  VLidort_Modified_Boolean(const VLidort_Modified_Boolean& src);
  ~VLidort_Modified_Boolean();

  const bool ts_do_additional_mvout() const;
  void ts_do_additional_mvout(const bool& ts_do_additional_mvout_in);
  
  const bool ts_do_bvp_telescoping() const;
  void ts_do_bvp_telescoping(const bool& ts_do_bvp_telescoping_in);
  
  const bool ts_do_chapman_function() const;
  void ts_do_chapman_function(const bool& ts_do_chapman_function_in);
  
  const bool ts_do_classical_solution() const;
  void ts_do_classical_solution(const bool& ts_do_classical_solution_in);
  
  const bool ts_do_deltam_scaling() const;
  void ts_do_deltam_scaling(const bool& ts_do_deltam_scaling_in);
  
  const bool ts_do_double_convtest() const;
  void ts_do_double_convtest(const bool& ts_do_double_convtest_in);
  
  const bool ts_do_doublet_geometry() const;
  void ts_do_doublet_geometry(const bool& ts_do_doublet_geometry_in);
  
  const bool ts_do_external_wleave() const;
  void ts_do_external_wleave(const bool& ts_do_external_wleave_in);
  
  const bool ts_do_focorr() const;
  void ts_do_focorr(const bool& ts_do_focorr_in);
  
  const bool ts_do_focorr_external() const;
  void ts_do_focorr_external(const bool& ts_do_focorr_external_in);
  
  const bool ts_do_focorr_nadir() const;
  void ts_do_focorr_nadir(const bool& ts_do_focorr_nadir_in);
  
  const bool ts_do_focorr_outgoing() const;
  void ts_do_focorr_outgoing(const bool& ts_do_focorr_outgoing_in);
  
  const bool ts_do_mvout_only() const;
  void ts_do_mvout_only(const bool& ts_do_mvout_only_in);
  
  const bool ts_do_observation_geometry() const;
  void ts_do_observation_geometry(const bool& ts_do_observation_geometry_in);
  
  const bool ts_do_rayleigh_only() const;
  void ts_do_rayleigh_only(const bool& ts_do_rayleigh_only_in);
  
  const bool ts_do_refractive_geometry() const;
  void ts_do_refractive_geometry(const bool& ts_do_refractive_geometry_in);
  
  const bool ts_do_solar_sources() const;
  void ts_do_solar_sources(const bool& ts_do_solar_sources_in);
  
  const bool ts_do_solution_saving() const;
  void ts_do_solution_saving(const bool& ts_do_solution_saving_in);
  
  const bool ts_do_sscorr_truncation() const;
  void ts_do_sscorr_truncation(const bool& ts_do_sscorr_truncation_in);
  
  const bool ts_do_sscorr_usefmat() const;
  void ts_do_sscorr_usefmat(const bool& ts_do_sscorr_usefmat_in);
  
  const bool ts_do_thermal_transonly() const;
  void ts_do_thermal_transonly(const bool& ts_do_thermal_transonly_in);
  
  const bool ts_do_user_streams() const;
  void ts_do_user_streams(const bool& ts_do_user_vzangles_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Control : public Spurr_Modified_Control_Base {
public:
  VLidort_Modified_Control();
  VLidort_Modified_Control(const VLidort_Modified_Control& src);
  ~VLidort_Modified_Control();

  const int& ts_nmoments_input() const;
  void ts_nmoments_input(const int& ts_ngreek_moments_input_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Sunrays : public Spurr_Modified_Sunrays_Base {
public:
  VLidort_Modified_Sunrays();
  VLidort_Modified_Sunrays(const VLidort_Modified_Sunrays& src);
  ~VLidort_Modified_Sunrays();

  const int& ts_nbeams() const;
  void ts_nbeams(const int& ts_n_szangles_in);
  
  const blitz::Array<double, 1>& ts_beam_szas() const;
  void ts_beam_szas(const blitz::Array<double, 1>& ts_szangles_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Uservalues : public Spurr_Modified_Uservalues_Base {
public:
  VLidort_Modified_Uservalues();
  VLidort_Modified_Uservalues(const VLidort_Modified_Uservalues& src);
  ~VLidort_Modified_Uservalues();

  const double& ts_geometry_specheight() const;
  void ts_geometry_specheight(const double& ts_geometry_specheight_in);
  
  const int& ts_n_user_doublets() const;
  void ts_n_user_doublets(const int& ts_n_user_doublets_in);
  
  const int& ts_n_user_obsgeoms() const;
  void ts_n_user_obsgeoms(const int& ts_n_user_obsgeoms_in);
  
  const int& ts_n_user_relazms() const;
  void ts_n_user_relazms(const int& ts_n_user_relazms_in);
  
  const int& ts_n_user_streams() const;
  void ts_n_user_streams(const int& ts_n_user_vzangles_in);
  
  const blitz::Array<double, 2>& ts_user_doublets() const;
  void ts_user_doublets(const blitz::Array<double, 2>& ts_user_doublets_in);
  
  const blitz::Array<double, 1>& ts_user_levels() const;
  void ts_user_levels(const blitz::Array<double, 1>& ts_user_levels_in);
  
  const blitz::Array<double, 2>& ts_user_obsgeoms_input() const;
  void ts_user_obsgeoms_input(const blitz::Array<double, 2>& ts_user_obsgeoms_input_in);
  
  const blitz::Array<double, 1>& ts_user_relazms() const;
  void ts_user_relazms(const blitz::Array<double, 1>& ts_user_relazms_in);
  
  const blitz::Array<double, 1>& ts_user_angles_input() const;
  void ts_user_angles_input(const blitz::Array<double, 1>& ts_user_vzangles_input_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Chapman : public Spurr_Modified_Chapman_Base {
public:
  VLidort_Modified_Chapman();
  VLidort_Modified_Chapman(const VLidort_Modified_Chapman& src);
  ~VLidort_Modified_Chapman();

  const double& ts_earth_radius() const;
  void ts_earth_radius(const double& ts_earth_radius_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Optical : public Spurr_Modified_Optical_Base {
public:
  VLidort_Modified_Optical();
  VLidort_Modified_Optical(const VLidort_Modified_Optical& src);
  ~VLidort_Modified_Optical();

  const blitz::Array<double, 1>& ts_omega_total_input() const;
  void ts_omega_total_input(const blitz::Array<double, 1>& ts_omega_total_input_in);
  
  virtual void print(std::ostream &output_stream) const;
};

class VLidort_Modified_Inputs : public Spurr_Modified_Inputs_Base {
public:
  VLidort_Modified_Inputs();
  VLidort_Modified_Inputs(const VLidort_Modified_Inputs& src);
  ~VLidort_Modified_Inputs();

  const VLidort_Modified_Boolean& mbool() const;
  const Spurr_Modified_Boolean_Base& modified_boolean_base() const;
  
  const VLidort_Modified_Chapman& mchapman() const;
  const Spurr_Modified_Chapman_Base& modified_chapman_base() const;
  
  const VLidort_Modified_Control& mcont() const;
  const Spurr_Modified_Control_Base& modified_control_base() const;
  
  const VLidort_Modified_Optical& moptical() const;
  const Spurr_Modified_Optical_Base& modified_optical_base() const;
  
  const VLidort_Modified_Sunrays& msunrays() const;
  const Spurr_Modified_Sunrays_Base& modified_sunrays_base() const;
  
  const VLidort_Modified_Uservalues& muserval() const;
  const Spurr_Modified_Uservalues_Base& modified_uservalues_base() const;
  
  virtual void print(std::ostream &output_stream) const;
};



}