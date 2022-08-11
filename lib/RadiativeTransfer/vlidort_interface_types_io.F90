module vlidort_interface_types_io

use iso_c_binding
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

! Links to type: "vbrdf_linsup_inputs" from module: "vbrdf_linsup_inputs_def_m" in file: "vbrdf_lin_sup_inputs_def.f90"
! Allocs and initializes type
subroutine vbrdf_linsup_inputs_c_write(lun, fortran_type_c) bind(C)
  use vbrdf_linsup_inputs_def_m, only : vbrdf_linsup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vbrdf_linsup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_linsup_inputs_f_write(lun, fortran_type_f)

end subroutine vbrdf_linsup_inputs_c_write

subroutine vbrdf_linsup_inputs_f_write(lun, fortran_type_f) 
  use vbrdf_linsup_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_linsup_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_do_bsavalue_wf
  write(UNIT=lun) fortran_type_f%bs_do_kernel_factor_wfs
  write(UNIT=lun) fortran_type_f%bs_do_kernel_params_wfs
  write(UNIT=lun) fortran_type_f%bs_do_kparams_derivs
  write(UNIT=lun) fortran_type_f%bs_do_windspeed_wf
  write(UNIT=lun) fortran_type_f%bs_do_wsavalue_wf
  write(UNIT=lun) fortran_type_f%bs_n_kernel_factor_wfs
  write(UNIT=lun) fortran_type_f%bs_n_kernel_params_wfs
  write(UNIT=lun) fortran_type_f%bs_n_surface_wfs
  
end subroutine vbrdf_linsup_inputs_f_write

subroutine vbrdf_linsup_inputs_c_read(lun, fortran_type_c) bind(C)
  use vbrdf_linsup_inputs_def_m, only : vbrdf_linsup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vbrdf_linsup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_linsup_inputs_f_read(lun, fortran_type_f)

end subroutine vbrdf_linsup_inputs_c_read

subroutine vbrdf_linsup_inputs_f_read(lun, fortran_type_f) 
  use vbrdf_linsup_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_linsup_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_do_bsavalue_wf
  read(UNIT=lun) fortran_type_f%bs_do_kernel_factor_wfs
  read(UNIT=lun) fortran_type_f%bs_do_kernel_params_wfs
  read(UNIT=lun) fortran_type_f%bs_do_kparams_derivs
  read(UNIT=lun) fortran_type_f%bs_do_windspeed_wf
  read(UNIT=lun) fortran_type_f%bs_do_wsavalue_wf
  read(UNIT=lun) fortran_type_f%bs_n_kernel_factor_wfs
  read(UNIT=lun) fortran_type_f%bs_n_kernel_params_wfs
  read(UNIT=lun) fortran_type_f%bs_n_surface_wfs
  
end subroutine vbrdf_linsup_inputs_f_read

! Links to type: "vbrdf_linsup_outputs" from module: "vbrdf_linsup_outputs_def_m" in file: "vbrdf_lin_sup_outputs_def.f90"
! Allocs and initializes type
subroutine vbrdf_linsup_outputs_c_write(lun, fortran_type_c) bind(C)
  use vbrdf_linsup_outputs_def_m, only : vbrdf_linsup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vbrdf_linsup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_linsup_outputs_f_write(lun, fortran_type_f)

end subroutine vbrdf_linsup_outputs_c_write

subroutine vbrdf_linsup_outputs_f_write(lun, fortran_type_f) 
  use vbrdf_linsup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_linsup_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_ls_brdf_f
  write(UNIT=lun) fortran_type_f%bs_ls_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_ls_dbounce_brdfunc
  write(UNIT=lun) fortran_type_f%bs_ls_emissivity
  write(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f
  write(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_ls_user_emissivity
  
end subroutine vbrdf_linsup_outputs_f_write

subroutine vbrdf_linsup_outputs_c_read(lun, fortran_type_c) bind(C)
  use vbrdf_linsup_outputs_def_m, only : vbrdf_linsup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vbrdf_linsup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_linsup_outputs_f_read(lun, fortran_type_f)

end subroutine vbrdf_linsup_outputs_c_read

subroutine vbrdf_linsup_outputs_f_read(lun, fortran_type_f) 
  use vbrdf_linsup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_linsup_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_ls_brdf_f
  read(UNIT=lun) fortran_type_f%bs_ls_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_ls_dbounce_brdfunc
  read(UNIT=lun) fortran_type_f%bs_ls_emissivity
  read(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f
  read(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_ls_user_emissivity
  
end subroutine vbrdf_linsup_outputs_f_read

! Links to type: "vbrdf_sup_inputs" from module: "vbrdf_sup_inputs_def_m" in file: "vbrdf_sup_inputs_def.f90"
! Allocs and initializes type
subroutine vbrdf_sup_inputs_c_write(lun, fortran_type_c) bind(C)
  use vbrdf_sup_inputs_def_m, only : vbrdf_sup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vbrdf_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_sup_inputs_f_write(lun, fortran_type_f)

end subroutine vbrdf_sup_inputs_c_write

subroutine vbrdf_sup_inputs_f_write(lun, fortran_type_f) 
  use vbrdf_sup_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_sup_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_beam_szas
  write(UNIT=lun) fortran_type_f%bs_brdf_factors
  write(UNIT=lun) fortran_type_f%bs_brdf_names
  write(UNIT=lun) fortran_type_f%bs_brdf_parameters
  write(UNIT=lun) fortran_type_f%bs_bsa_value
  write(UNIT=lun) fortran_type_f%bs_do_brdf_surface
  write(UNIT=lun) fortran_type_f%bs_do_bsa_scaling
  write(UNIT=lun) fortran_type_f%bs_do_directbounce_only
  write(UNIT=lun) fortran_type_f%bs_do_doublet_geometry
  write(UNIT=lun) fortran_type_f%bs_do_facetisotropy
  write(UNIT=lun) fortran_type_f%bs_do_foamoption
  write(UNIT=lun) fortran_type_f%bs_do_glintshadow
  write(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr
  write(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr_dbonly
  write(UNIT=lun) fortran_type_f%bs_do_newcmglint
  write(UNIT=lun) fortran_type_f%bs_do_newgcmglint
  write(UNIT=lun) fortran_type_f%bs_do_shadow_effect
  write(UNIT=lun) fortran_type_f%bs_do_solar_sources
  write(UNIT=lun) fortran_type_f%bs_do_surface_emission
  write(UNIT=lun) fortran_type_f%bs_do_user_obsgeoms
  write(UNIT=lun) fortran_type_f%bs_do_user_streams
  write(UNIT=lun) fortran_type_f%bs_do_wsa_scaling
  write(UNIT=lun) fortran_type_f%bs_do_wsabsa_output
  write(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nmuquad
  write(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nphiquad
  write(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_order
  write(UNIT=lun) fortran_type_f%bs_lambertian_kernel_flag
  write(UNIT=lun) fortran_type_f%bs_n_brdf_kernels
  write(UNIT=lun) fortran_type_f%bs_n_brdf_parameters
  write(UNIT=lun) fortran_type_f%bs_n_user_doublets
  write(UNIT=lun) fortran_type_f%bs_n_user_obsgeoms
  write(UNIT=lun) fortran_type_f%bs_n_user_relazms
  write(UNIT=lun) fortran_type_f%bs_n_user_streams
  write(UNIT=lun) fortran_type_f%bs_nbeams
  write(UNIT=lun) fortran_type_f%bs_nstokes
  write(UNIT=lun) fortran_type_f%bs_nstreams
  write(UNIT=lun) fortran_type_f%bs_nstreams_brdf
  write(UNIT=lun) fortran_type_f%bs_salinity
  write(UNIT=lun) fortran_type_f%bs_user_angles_input
  write(UNIT=lun) fortran_type_f%bs_user_doublets
  write(UNIT=lun) fortran_type_f%bs_user_obsgeoms
  write(UNIT=lun) fortran_type_f%bs_user_relazms
  write(UNIT=lun) fortran_type_f%bs_wavelength
  write(UNIT=lun) fortran_type_f%bs_which_brdf
  write(UNIT=lun) fortran_type_f%bs_winddir
  write(UNIT=lun) fortran_type_f%bs_windspeed
  write(UNIT=lun) fortran_type_f%bs_wsa_value
  
end subroutine vbrdf_sup_inputs_f_write

subroutine vbrdf_sup_inputs_c_read(lun, fortran_type_c) bind(C)
  use vbrdf_sup_inputs_def_m, only : vbrdf_sup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vbrdf_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_sup_inputs_f_read(lun, fortran_type_f)

end subroutine vbrdf_sup_inputs_c_read

subroutine vbrdf_sup_inputs_f_read(lun, fortran_type_f) 
  use vbrdf_sup_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_sup_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_beam_szas
  read(UNIT=lun) fortran_type_f%bs_brdf_factors
  read(UNIT=lun) fortran_type_f%bs_brdf_names
  read(UNIT=lun) fortran_type_f%bs_brdf_parameters
  read(UNIT=lun) fortran_type_f%bs_bsa_value
  read(UNIT=lun) fortran_type_f%bs_do_brdf_surface
  read(UNIT=lun) fortran_type_f%bs_do_bsa_scaling
  read(UNIT=lun) fortran_type_f%bs_do_directbounce_only
  read(UNIT=lun) fortran_type_f%bs_do_doublet_geometry
  read(UNIT=lun) fortran_type_f%bs_do_facetisotropy
  read(UNIT=lun) fortran_type_f%bs_do_foamoption
  read(UNIT=lun) fortran_type_f%bs_do_glintshadow
  read(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr
  read(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr_dbonly
  read(UNIT=lun) fortran_type_f%bs_do_newcmglint
  read(UNIT=lun) fortran_type_f%bs_do_newgcmglint
  read(UNIT=lun) fortran_type_f%bs_do_shadow_effect
  read(UNIT=lun) fortran_type_f%bs_do_solar_sources
  read(UNIT=lun) fortran_type_f%bs_do_surface_emission
  read(UNIT=lun) fortran_type_f%bs_do_user_obsgeoms
  read(UNIT=lun) fortran_type_f%bs_do_user_streams
  read(UNIT=lun) fortran_type_f%bs_do_wsa_scaling
  read(UNIT=lun) fortran_type_f%bs_do_wsabsa_output
  read(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nmuquad
  read(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nphiquad
  read(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_order
  read(UNIT=lun) fortran_type_f%bs_lambertian_kernel_flag
  read(UNIT=lun) fortran_type_f%bs_n_brdf_kernels
  read(UNIT=lun) fortran_type_f%bs_n_brdf_parameters
  read(UNIT=lun) fortran_type_f%bs_n_user_doublets
  read(UNIT=lun) fortran_type_f%bs_n_user_obsgeoms
  read(UNIT=lun) fortran_type_f%bs_n_user_relazms
  read(UNIT=lun) fortran_type_f%bs_n_user_streams
  read(UNIT=lun) fortran_type_f%bs_nbeams
  read(UNIT=lun) fortran_type_f%bs_nstokes
  read(UNIT=lun) fortran_type_f%bs_nstreams
  read(UNIT=lun) fortran_type_f%bs_nstreams_brdf
  read(UNIT=lun) fortran_type_f%bs_salinity
  read(UNIT=lun) fortran_type_f%bs_user_angles_input
  read(UNIT=lun) fortran_type_f%bs_user_doublets
  read(UNIT=lun) fortran_type_f%bs_user_obsgeoms
  read(UNIT=lun) fortran_type_f%bs_user_relazms
  read(UNIT=lun) fortran_type_f%bs_wavelength
  read(UNIT=lun) fortran_type_f%bs_which_brdf
  read(UNIT=lun) fortran_type_f%bs_winddir
  read(UNIT=lun) fortran_type_f%bs_windspeed
  read(UNIT=lun) fortran_type_f%bs_wsa_value
  
end subroutine vbrdf_sup_inputs_f_read

! Links to type: "vbrdf_sup_outputs" from module: "vbrdf_sup_outputs_def_m" in file: "vbrdf_sup_outputs_def.f90"
! Allocs and initializes type
subroutine vbrdf_sup_outputs_c_write(lun, fortran_type_c) bind(C)
  use vbrdf_sup_outputs_def_m, only : vbrdf_sup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vbrdf_sup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_sup_outputs_f_write(lun, fortran_type_f)

end subroutine vbrdf_sup_outputs_c_write

subroutine vbrdf_sup_outputs_f_write(lun, fortran_type_f) 
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_sup_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_brdf_f
  write(UNIT=lun) fortran_type_f%bs_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_bsa_calculated
  write(UNIT=lun) fortran_type_f%bs_bsa_kernels
  write(UNIT=lun) fortran_type_f%bs_dbounce_brdfunc
  write(UNIT=lun) fortran_type_f%bs_emissivity
  write(UNIT=lun) fortran_type_f%bs_user_brdf_f
  write(UNIT=lun) fortran_type_f%bs_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_user_emissivity
  write(UNIT=lun) fortran_type_f%bs_wsa_calculated
  write(UNIT=lun) fortran_type_f%bs_wsa_kernels
  
end subroutine vbrdf_sup_outputs_f_write

subroutine vbrdf_sup_outputs_c_read(lun, fortran_type_c) bind(C)
  use vbrdf_sup_outputs_def_m, only : vbrdf_sup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vbrdf_sup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_sup_outputs_f_read(lun, fortran_type_f)

end subroutine vbrdf_sup_outputs_c_read

subroutine vbrdf_sup_outputs_f_read(lun, fortran_type_f) 
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_sup_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_brdf_f
  read(UNIT=lun) fortran_type_f%bs_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_bsa_calculated
  read(UNIT=lun) fortran_type_f%bs_bsa_kernels
  read(UNIT=lun) fortran_type_f%bs_dbounce_brdfunc
  read(UNIT=lun) fortran_type_f%bs_emissivity
  read(UNIT=lun) fortran_type_f%bs_user_brdf_f
  read(UNIT=lun) fortran_type_f%bs_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_user_emissivity
  read(UNIT=lun) fortran_type_f%bs_wsa_calculated
  read(UNIT=lun) fortran_type_f%bs_wsa_kernels
  
end subroutine vbrdf_sup_outputs_f_read

! Links to type: "vbrdf_input_exception_handling" from module: "vbrdf_sup_outputs_def_m" in file: "vbrdf_sup_outputs_def.f90"
! Allocs and initializes type
subroutine vbrdf_input_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use vbrdf_sup_outputs_def_m, only : vbrdf_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vbrdf_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_input_exception_handling_f_write(lun, fortran_type_f)

end subroutine vbrdf_input_exception_handling_c_write

subroutine vbrdf_input_exception_handling_f_write(lun, fortran_type_f) 
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_input_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_inputactions
  write(UNIT=lun) fortran_type_f%bs_inputmessages
  write(UNIT=lun) fortran_type_f%bs_ninputmessages
  write(UNIT=lun) fortran_type_f%bs_status_inputread
  
end subroutine vbrdf_input_exception_handling_f_write

subroutine vbrdf_input_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use vbrdf_sup_outputs_def_m, only : vbrdf_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vbrdf_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_input_exception_handling_f_read(lun, fortran_type_f)

end subroutine vbrdf_input_exception_handling_c_read

subroutine vbrdf_input_exception_handling_f_read(lun, fortran_type_f) 
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_input_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_inputactions
  read(UNIT=lun) fortran_type_f%bs_inputmessages
  read(UNIT=lun) fortran_type_f%bs_ninputmessages
  read(UNIT=lun) fortran_type_f%bs_status_inputread
  
end subroutine vbrdf_input_exception_handling_f_read

! Links to type: "vbrdf_output_exception_handling" from module: "vbrdf_sup_outputs_def_m" in file: "vbrdf_sup_outputs_def.f90"
! Allocs and initializes type
subroutine vbrdf_output_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use vbrdf_sup_outputs_def_m, only : vbrdf_output_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vbrdf_output_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_output_exception_handling_f_write(lun, fortran_type_f)

end subroutine vbrdf_output_exception_handling_c_write

subroutine vbrdf_output_exception_handling_f_write(lun, fortran_type_f) 
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_output_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_noutputmessages
  write(UNIT=lun) fortran_type_f%bs_outputmessages
  write(UNIT=lun) fortran_type_f%bs_status_output
  
end subroutine vbrdf_output_exception_handling_f_write

subroutine vbrdf_output_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use vbrdf_sup_outputs_def_m, only : vbrdf_output_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vbrdf_output_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vbrdf_output_exception_handling_f_read(lun, fortran_type_f)

end subroutine vbrdf_output_exception_handling_c_read

subroutine vbrdf_output_exception_handling_f_read(lun, fortran_type_f) 
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vbrdf_output_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_noutputmessages
  read(UNIT=lun) fortran_type_f%bs_outputmessages
  read(UNIT=lun) fortran_type_f%bs_status_output
  
end subroutine vbrdf_output_exception_handling_f_read

! Links to type: "vsleave_sup_inputs" from module: "vsleave_sup_inputs_def_m" in file: "vsleave_sup_inputs_def.f90"
! Allocs and initializes type
subroutine vsleave_sup_inputs_c_write(lun, fortran_type_c) bind(C)
  use vsleave_sup_inputs_def_m, only : vsleave_sup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vsleave_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vsleave_sup_inputs_f_write(lun, fortran_type_f)

end subroutine vsleave_sup_inputs_c_write

subroutine vsleave_sup_inputs_f_write(lun, fortran_type_f) 
  use vsleave_sup_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vsleave_sup_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%sl_azimuthdep
  write(UNIT=lun) fortran_type_f%sl_beam_szas
  write(UNIT=lun) fortran_type_f%sl_chlorconc
  write(UNIT=lun) fortran_type_f%sl_do_doublet_geometry
  write(UNIT=lun) fortran_type_f%sl_do_exact
  write(UNIT=lun) fortran_type_f%sl_do_exactonly
  write(UNIT=lun) fortran_type_f%sl_do_facetisotropy
  write(UNIT=lun) fortran_type_f%sl_do_fluorescence
  write(UNIT=lun) fortran_type_f%sl_do_foamoption
  write(UNIT=lun) fortran_type_f%sl_do_fourier_output
  write(UNIT=lun) fortran_type_f%sl_do_glintshadow
  write(UNIT=lun) fortran_type_f%sl_do_isotropic
  write(UNIT=lun) fortran_type_f%sl_do_roughsurface
  write(UNIT=lun) fortran_type_f%sl_do_sleaving
  write(UNIT=lun) fortran_type_f%sl_do_solar_sources
  write(UNIT=lun) fortran_type_f%sl_do_user_obsgeoms
  write(UNIT=lun) fortran_type_f%sl_do_user_streams
  write(UNIT=lun) fortran_type_f%sl_fl_amplitude755
  write(UNIT=lun) fortran_type_f%sl_fl_do_datagaussian
  write(UNIT=lun) fortran_type_f%sl_fl_epoch
  write(UNIT=lun) fortran_type_f%sl_fl_inputgaussians
  write(UNIT=lun) fortran_type_f%sl_fl_latitude
  write(UNIT=lun) fortran_type_f%sl_fl_longitude
  write(UNIT=lun) fortran_type_f%sl_fl_wavelength
  write(UNIT=lun) fortran_type_f%sl_n_user_doublets
  write(UNIT=lun) fortran_type_f%sl_n_user_obsgeoms
  write(UNIT=lun) fortran_type_f%sl_n_user_relazms
  write(UNIT=lun) fortran_type_f%sl_n_user_streams
  write(UNIT=lun) fortran_type_f%sl_nbeams
  write(UNIT=lun) fortran_type_f%sl_nstokes
  write(UNIT=lun) fortran_type_f%sl_nstreams
  write(UNIT=lun) fortran_type_f%sl_salinity
  write(UNIT=lun) fortran_type_f%sl_user_angles_input
  write(UNIT=lun) fortran_type_f%sl_user_doublets
  write(UNIT=lun) fortran_type_f%sl_user_obsgeoms
  write(UNIT=lun) fortran_type_f%sl_user_relazms
  write(UNIT=lun) fortran_type_f%sl_vsleave_datapath
  write(UNIT=lun) fortran_type_f%sl_wavelength
  write(UNIT=lun) fortran_type_f%sl_winddir
  write(UNIT=lun) fortran_type_f%sl_windspeed
  
end subroutine vsleave_sup_inputs_f_write

subroutine vsleave_sup_inputs_c_read(lun, fortran_type_c) bind(C)
  use vsleave_sup_inputs_def_m, only : vsleave_sup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vsleave_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vsleave_sup_inputs_f_read(lun, fortran_type_f)

end subroutine vsleave_sup_inputs_c_read

subroutine vsleave_sup_inputs_f_read(lun, fortran_type_f) 
  use vsleave_sup_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vsleave_sup_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%sl_azimuthdep
  read(UNIT=lun) fortran_type_f%sl_beam_szas
  read(UNIT=lun) fortran_type_f%sl_chlorconc
  read(UNIT=lun) fortran_type_f%sl_do_doublet_geometry
  read(UNIT=lun) fortran_type_f%sl_do_exact
  read(UNIT=lun) fortran_type_f%sl_do_exactonly
  read(UNIT=lun) fortran_type_f%sl_do_facetisotropy
  read(UNIT=lun) fortran_type_f%sl_do_fluorescence
  read(UNIT=lun) fortran_type_f%sl_do_foamoption
  read(UNIT=lun) fortran_type_f%sl_do_fourier_output
  read(UNIT=lun) fortran_type_f%sl_do_glintshadow
  read(UNIT=lun) fortran_type_f%sl_do_isotropic
  read(UNIT=lun) fortran_type_f%sl_do_roughsurface
  read(UNIT=lun) fortran_type_f%sl_do_sleaving
  read(UNIT=lun) fortran_type_f%sl_do_solar_sources
  read(UNIT=lun) fortran_type_f%sl_do_user_obsgeoms
  read(UNIT=lun) fortran_type_f%sl_do_user_streams
  read(UNIT=lun) fortran_type_f%sl_fl_amplitude755
  read(UNIT=lun) fortran_type_f%sl_fl_do_datagaussian
  read(UNIT=lun) fortran_type_f%sl_fl_epoch
  read(UNIT=lun) fortran_type_f%sl_fl_inputgaussians
  read(UNIT=lun) fortran_type_f%sl_fl_latitude
  read(UNIT=lun) fortran_type_f%sl_fl_longitude
  read(UNIT=lun) fortran_type_f%sl_fl_wavelength
  read(UNIT=lun) fortran_type_f%sl_n_user_doublets
  read(UNIT=lun) fortran_type_f%sl_n_user_obsgeoms
  read(UNIT=lun) fortran_type_f%sl_n_user_relazms
  read(UNIT=lun) fortran_type_f%sl_n_user_streams
  read(UNIT=lun) fortran_type_f%sl_nbeams
  read(UNIT=lun) fortran_type_f%sl_nstokes
  read(UNIT=lun) fortran_type_f%sl_nstreams
  read(UNIT=lun) fortran_type_f%sl_salinity
  read(UNIT=lun) fortran_type_f%sl_user_angles_input
  read(UNIT=lun) fortran_type_f%sl_user_doublets
  read(UNIT=lun) fortran_type_f%sl_user_obsgeoms
  read(UNIT=lun) fortran_type_f%sl_user_relazms
  read(UNIT=lun) fortran_type_f%sl_vsleave_datapath
  read(UNIT=lun) fortran_type_f%sl_wavelength
  read(UNIT=lun) fortran_type_f%sl_winddir
  read(UNIT=lun) fortran_type_f%sl_windspeed
  
end subroutine vsleave_sup_inputs_f_read

! Links to type: "vlidort_fixed_lincontrol" from module: "vlidort_lininputs_def_m" in file: "vlidort_lin_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_lincontrol_c_write(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_fixed_lincontrol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_lincontrol_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_lincontrol_c_write

subroutine vlidort_fixed_lincontrol_f_write(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_lincontrol), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_columnwf_names
  write(UNIT=lun) fortran_type_f%ts_layer_vary_flag
  write(UNIT=lun) fortran_type_f%ts_layer_vary_number
  write(UNIT=lun) fortran_type_f%ts_n_sleave_wfs
  write(UNIT=lun) fortran_type_f%ts_n_surface_wfs
  write(UNIT=lun) fortran_type_f%ts_n_totalcolumn_wfs
  write(UNIT=lun) fortran_type_f%ts_n_totalprofile_wfs
  write(UNIT=lun) fortran_type_f%ts_profilewf_names
  
end subroutine vlidort_fixed_lincontrol_f_write

subroutine vlidort_fixed_lincontrol_c_read(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_fixed_lincontrol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_lincontrol_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_lincontrol_c_read

subroutine vlidort_fixed_lincontrol_f_read(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_lincontrol), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_columnwf_names
  read(UNIT=lun) fortran_type_f%ts_layer_vary_flag
  read(UNIT=lun) fortran_type_f%ts_layer_vary_number
  read(UNIT=lun) fortran_type_f%ts_n_sleave_wfs
  read(UNIT=lun) fortran_type_f%ts_n_surface_wfs
  read(UNIT=lun) fortran_type_f%ts_n_totalcolumn_wfs
  read(UNIT=lun) fortran_type_f%ts_n_totalprofile_wfs
  read(UNIT=lun) fortran_type_f%ts_profilewf_names
  
end subroutine vlidort_fixed_lincontrol_f_read

! Links to type: "vlidort_fixed_linoptical" from module: "vlidort_lininputs_def_m" in file: "vlidort_lin_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_linoptical_c_write(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_fixed_linoptical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_linoptical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_linoptical_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_linoptical_c_write

subroutine vlidort_fixed_linoptical_f_write(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_linoptical), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_l_deltau_vert_input
  write(UNIT=lun) fortran_type_f%ts_l_fmatrix_dn
  write(UNIT=lun) fortran_type_f%ts_l_fmatrix_up
  write(UNIT=lun) fortran_type_f%ts_l_greekmat_total_input
  write(UNIT=lun) fortran_type_f%ts_l_omega_total_input
  
end subroutine vlidort_fixed_linoptical_f_write

subroutine vlidort_fixed_linoptical_c_read(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_fixed_linoptical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_linoptical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_linoptical_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_linoptical_c_read

subroutine vlidort_fixed_linoptical_f_read(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_linoptical), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_l_deltau_vert_input
  read(UNIT=lun) fortran_type_f%ts_l_fmatrix_dn
  read(UNIT=lun) fortran_type_f%ts_l_fmatrix_up
  read(UNIT=lun) fortran_type_f%ts_l_greekmat_total_input
  read(UNIT=lun) fortran_type_f%ts_l_omega_total_input
  
end subroutine vlidort_fixed_linoptical_f_read

! Links to type: "vlidort_fixed_lininputs" from module: "vlidort_lininputs_def_m" in file: "vlidort_lin_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_lininputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_fixed_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_lininputs_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_lininputs_c_write

subroutine vlidort_fixed_lininputs_f_write(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_lininputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_fixed_lincontrol), pointer :: cont_lcl  
  type(vlidort_fixed_linoptical), pointer :: optical_lcl  
  
  ! Get pointer to types
  cont_lcl => fortran_type_f%cont
  optical_lcl => fortran_type_f%optical
  
  call vlidort_fixed_lincontrol_f_write(lun, cont_lcl)
  call vlidort_fixed_linoptical_f_write(lun, optical_lcl)
  
end subroutine vlidort_fixed_lininputs_f_write

subroutine vlidort_fixed_lininputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_fixed_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_lininputs_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_lininputs_c_read

subroutine vlidort_fixed_lininputs_f_read(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_lininputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_fixed_lincontrol), pointer :: cont_lcl  
  type(vlidort_fixed_linoptical), pointer :: optical_lcl  
  
  ! Get pointer to types
  cont_lcl => fortran_type_f%cont
  optical_lcl => fortran_type_f%optical
  
  call vlidort_fixed_lincontrol_f_read(lun, cont_lcl)
  call vlidort_fixed_linoptical_f_read(lun, optical_lcl)
  
end subroutine vlidort_fixed_lininputs_f_read

! Links to type: "vlidort_modified_lincontrol" from module: "vlidort_lininputs_def_m" in file: "vlidort_lin_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_lincontrol_c_write(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_modified_lincontrol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_lincontrol_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_lincontrol_c_write

subroutine vlidort_modified_lincontrol_f_write(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_lincontrol), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_atmos_lbbf
  write(UNIT=lun) fortran_type_f%ts_do_atmos_linearization
  write(UNIT=lun) fortran_type_f%ts_do_column_linearization
  write(UNIT=lun) fortran_type_f%ts_do_linearization
  write(UNIT=lun) fortran_type_f%ts_do_profile_linearization
  write(UNIT=lun) fortran_type_f%ts_do_simulation_only
  write(UNIT=lun) fortran_type_f%ts_do_sleave_wfs
  write(UNIT=lun) fortran_type_f%ts_do_surface_lbbf
  write(UNIT=lun) fortran_type_f%ts_do_surface_linearization
  
end subroutine vlidort_modified_lincontrol_f_write

subroutine vlidort_modified_lincontrol_c_read(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_modified_lincontrol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_lincontrol_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_lincontrol_c_read

subroutine vlidort_modified_lincontrol_f_read(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_lincontrol), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_atmos_lbbf
  read(UNIT=lun) fortran_type_f%ts_do_atmos_linearization
  read(UNIT=lun) fortran_type_f%ts_do_column_linearization
  read(UNIT=lun) fortran_type_f%ts_do_linearization
  read(UNIT=lun) fortran_type_f%ts_do_profile_linearization
  read(UNIT=lun) fortran_type_f%ts_do_simulation_only
  read(UNIT=lun) fortran_type_f%ts_do_sleave_wfs
  read(UNIT=lun) fortran_type_f%ts_do_surface_lbbf
  read(UNIT=lun) fortran_type_f%ts_do_surface_linearization
  
end subroutine vlidort_modified_lincontrol_f_read

! Links to type: "vlidort_modified_lininputs" from module: "vlidort_lininputs_def_m" in file: "vlidort_lin_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_lininputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_modified_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_lininputs_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_lininputs_c_write

subroutine vlidort_modified_lininputs_f_write(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_lininputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_modified_lincontrol), pointer :: mcont_lcl  
  
  ! Get pointer to types
  mcont_lcl => fortran_type_f%mcont
  
  call vlidort_modified_lincontrol_f_write(lun, mcont_lcl)
  
end subroutine vlidort_modified_lininputs_f_write

subroutine vlidort_modified_lininputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_lininputs_def_m, only : vlidort_modified_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_lininputs_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_lininputs_c_read

subroutine vlidort_modified_lininputs_f_read(lun, fortran_type_f) 
  use vlidort_lininputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_lininputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_modified_lincontrol), pointer :: mcont_lcl  
  
  ! Get pointer to types
  mcont_lcl => fortran_type_f%mcont
  
  call vlidort_modified_lincontrol_f_read(lun, mcont_lcl)
  
end subroutine vlidort_modified_lininputs_f_read

! Links to type: "vlidort_lincol" from module: "vlidort_linoutputs_def_m" in file: "vlidort_lin_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_lincol_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_lincol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_lincol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_lincol_f_write(lun, fortran_type_f)

end subroutine vlidort_lincol_c_write

subroutine vlidort_lincol_f_write(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_lincol), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_albmed_fluxes_colwf
  write(UNIT=lun) fortran_type_f%ts_albmed_user_colwf
  write(UNIT=lun) fortran_type_f%ts_columnwf
  write(UNIT=lun) fortran_type_f%ts_dnflux_direct_colwf
  write(UNIT=lun) fortran_type_f%ts_dnmeanst_direct_colwf
  write(UNIT=lun) fortran_type_f%ts_flux_diffuse_colwf
  write(UNIT=lun) fortran_type_f%ts_lc_layer_mssts
  write(UNIT=lun) fortran_type_f%ts_lc_lostrans
  write(UNIT=lun) fortran_type_f%ts_lc_surf_mssts
  write(UNIT=lun) fortran_type_f%ts_meanst_diffuse_colwf
  write(UNIT=lun) fortran_type_f%ts_planetary_sbterm_colwf
  write(UNIT=lun) fortran_type_f%ts_planetary_transterm_colwf
  write(UNIT=lun) fortran_type_f%ts_transbeam_colwf
  write(UNIT=lun) fortran_type_f%ts_trnmed_fluxes_colwf
  write(UNIT=lun) fortran_type_f%ts_trnmed_user_colwf
  
end subroutine vlidort_lincol_f_write

subroutine vlidort_lincol_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_lincol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_lincol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_lincol_f_read(lun, fortran_type_f)

end subroutine vlidort_lincol_c_read

subroutine vlidort_lincol_f_read(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_lincol), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_albmed_fluxes_colwf
  read(UNIT=lun) fortran_type_f%ts_albmed_user_colwf
  read(UNIT=lun) fortran_type_f%ts_columnwf
  read(UNIT=lun) fortran_type_f%ts_dnflux_direct_colwf
  read(UNIT=lun) fortran_type_f%ts_dnmeanst_direct_colwf
  read(UNIT=lun) fortran_type_f%ts_flux_diffuse_colwf
  read(UNIT=lun) fortran_type_f%ts_lc_layer_mssts
  read(UNIT=lun) fortran_type_f%ts_lc_lostrans
  read(UNIT=lun) fortran_type_f%ts_lc_surf_mssts
  read(UNIT=lun) fortran_type_f%ts_meanst_diffuse_colwf
  read(UNIT=lun) fortran_type_f%ts_planetary_sbterm_colwf
  read(UNIT=lun) fortran_type_f%ts_planetary_transterm_colwf
  read(UNIT=lun) fortran_type_f%ts_transbeam_colwf
  read(UNIT=lun) fortran_type_f%ts_trnmed_fluxes_colwf
  read(UNIT=lun) fortran_type_f%ts_trnmed_user_colwf
  
end subroutine vlidort_lincol_f_read

! Links to type: "vlidort_linprof" from module: "vlidort_linoutputs_def_m" in file: "vlidort_lin_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_linprof_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linprof

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linprof), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linprof_f_write(lun, fortran_type_f)

end subroutine vlidort_linprof_c_write

subroutine vlidort_linprof_f_write(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linprof), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_albmed_fluxes_profwf
  write(UNIT=lun) fortran_type_f%ts_albmed_user_profwf
  write(UNIT=lun) fortran_type_f%ts_dnflux_direct_profwf
  write(UNIT=lun) fortran_type_f%ts_dnmeanst_direct_profwf
  write(UNIT=lun) fortran_type_f%ts_flux_diffuse_profwf
  write(UNIT=lun) fortran_type_f%ts_lp_layer_mssts
  write(UNIT=lun) fortran_type_f%ts_lp_lostrans
  write(UNIT=lun) fortran_type_f%ts_lp_surf_mssts
  write(UNIT=lun) fortran_type_f%ts_meanst_diffuse_profwf
  write(UNIT=lun) fortran_type_f%ts_planetary_sbterm_profwf
  write(UNIT=lun) fortran_type_f%ts_planetary_transterm_profwf
  write(UNIT=lun) fortran_type_f%ts_profilewf
  write(UNIT=lun) fortran_type_f%ts_transbeam_profwf
  write(UNIT=lun) fortran_type_f%ts_trnmed_fluxes_profwf
  write(UNIT=lun) fortran_type_f%ts_trnmed_user_profwf
  
end subroutine vlidort_linprof_f_write

subroutine vlidort_linprof_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linprof

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linprof), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linprof_f_read(lun, fortran_type_f)

end subroutine vlidort_linprof_c_read

subroutine vlidort_linprof_f_read(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linprof), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_albmed_fluxes_profwf
  read(UNIT=lun) fortran_type_f%ts_albmed_user_profwf
  read(UNIT=lun) fortran_type_f%ts_dnflux_direct_profwf
  read(UNIT=lun) fortran_type_f%ts_dnmeanst_direct_profwf
  read(UNIT=lun) fortran_type_f%ts_flux_diffuse_profwf
  read(UNIT=lun) fortran_type_f%ts_lp_layer_mssts
  read(UNIT=lun) fortran_type_f%ts_lp_lostrans
  read(UNIT=lun) fortran_type_f%ts_lp_surf_mssts
  read(UNIT=lun) fortran_type_f%ts_meanst_diffuse_profwf
  read(UNIT=lun) fortran_type_f%ts_planetary_sbterm_profwf
  read(UNIT=lun) fortran_type_f%ts_planetary_transterm_profwf
  read(UNIT=lun) fortran_type_f%ts_profilewf
  read(UNIT=lun) fortran_type_f%ts_transbeam_profwf
  read(UNIT=lun) fortran_type_f%ts_trnmed_fluxes_profwf
  read(UNIT=lun) fortran_type_f%ts_trnmed_user_profwf
  
end subroutine vlidort_linprof_f_read

! Links to type: "vlidort_linatmos" from module: "vlidort_linoutputs_def_m" in file: "vlidort_lin_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_linatmos_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linatmos

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linatmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linatmos_f_write(lun, fortran_type_f)

end subroutine vlidort_linatmos_c_write

subroutine vlidort_linatmos_f_write(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linatmos), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_abbwfs_fluxes
  write(UNIT=lun) fortran_type_f%ts_abbwfs_jacobians
  
end subroutine vlidort_linatmos_f_write

subroutine vlidort_linatmos_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linatmos

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linatmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linatmos_f_read(lun, fortran_type_f)

end subroutine vlidort_linatmos_c_read

subroutine vlidort_linatmos_f_read(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linatmos), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_abbwfs_fluxes
  read(UNIT=lun) fortran_type_f%ts_abbwfs_jacobians
  
end subroutine vlidort_linatmos_f_read

! Links to type: "vlidort_linsurf" from module: "vlidort_linoutputs_def_m" in file: "vlidort_lin_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_linsurf_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linsurf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsurf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsurf_f_write(lun, fortran_type_f)

end subroutine vlidort_linsurf_c_write

subroutine vlidort_linsurf_f_write(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsurf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_flux_diffuse_surfwf
  write(UNIT=lun) fortran_type_f%ts_ls_layer_mssts
  write(UNIT=lun) fortran_type_f%ts_ls_surf_mssts
  write(UNIT=lun) fortran_type_f%ts_meanst_diffuse_surfwf
  write(UNIT=lun) fortran_type_f%ts_sbbwfs_fluxes
  write(UNIT=lun) fortran_type_f%ts_sbbwfs_jacobians
  write(UNIT=lun) fortran_type_f%ts_surfacewf
  
end subroutine vlidort_linsurf_f_write

subroutine vlidort_linsurf_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linsurf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsurf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsurf_f_read(lun, fortran_type_f)

end subroutine vlidort_linsurf_c_read

subroutine vlidort_linsurf_f_read(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsurf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_flux_diffuse_surfwf
  read(UNIT=lun) fortran_type_f%ts_ls_layer_mssts
  read(UNIT=lun) fortran_type_f%ts_ls_surf_mssts
  read(UNIT=lun) fortran_type_f%ts_meanst_diffuse_surfwf
  read(UNIT=lun) fortran_type_f%ts_sbbwfs_fluxes
  read(UNIT=lun) fortran_type_f%ts_sbbwfs_jacobians
  read(UNIT=lun) fortran_type_f%ts_surfacewf
  
end subroutine vlidort_linsurf_f_read

! Links to type: "vlidort_linoutputs" from module: "vlidort_linoutputs_def_m" in file: "vlidort_lin_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_linoutputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linoutputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linoutputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linoutputs_f_write(lun, fortran_type_f)

end subroutine vlidort_linoutputs_c_write

subroutine vlidort_linoutputs_f_write(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linoutputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_linatmos), pointer :: atmos_lcl  
  type(vlidort_lincol), pointer :: col_lcl  
  type(vlidort_linprof), pointer :: prof_lcl  
  type(vlidort_linsurf), pointer :: surf_lcl  
  
  ! Get pointer to types
  atmos_lcl => fortran_type_f%atmos
  col_lcl => fortran_type_f%col
  prof_lcl => fortran_type_f%prof
  surf_lcl => fortran_type_f%surf
  
  call vlidort_linatmos_f_write(lun, atmos_lcl)
  call vlidort_lincol_f_write(lun, col_lcl)
  call vlidort_linprof_f_write(lun, prof_lcl)
  call vlidort_linsurf_f_write(lun, surf_lcl)
  
end subroutine vlidort_linoutputs_f_write

subroutine vlidort_linoutputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linoutputs_def_m, only : vlidort_linoutputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linoutputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linoutputs_f_read(lun, fortran_type_f)

end subroutine vlidort_linoutputs_c_read

subroutine vlidort_linoutputs_f_read(lun, fortran_type_f) 
  use vlidort_linoutputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linoutputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_linatmos), pointer :: atmos_lcl  
  type(vlidort_lincol), pointer :: col_lcl  
  type(vlidort_linprof), pointer :: prof_lcl  
  type(vlidort_linsurf), pointer :: surf_lcl  
  
  ! Get pointer to types
  atmos_lcl => fortran_type_f%atmos
  col_lcl => fortran_type_f%col
  prof_lcl => fortran_type_f%prof
  surf_lcl => fortran_type_f%surf
  
  call vlidort_linatmos_f_read(lun, atmos_lcl)
  call vlidort_lincol_f_read(lun, col_lcl)
  call vlidort_linprof_f_read(lun, prof_lcl)
  call vlidort_linsurf_f_read(lun, surf_lcl)
  
end subroutine vlidort_linoutputs_f_read

! Links to type: "vlidort_linsup_brdf" from module: "vlidort_linsup_brdf_def_m" in file: "vlidort_lin_sup_brdf_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_brdf_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_brdf_def_m, only : vlidort_linsup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_brdf_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_brdf_c_write

subroutine vlidort_linsup_brdf_f_write(lun, fortran_type_f) 
  use vlidort_linsup_brdf_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_brdf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_ls_brdf_f
  write(UNIT=lun) fortran_type_f%ts_ls_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_ls_emissivity
  write(UNIT=lun) fortran_type_f%ts_ls_exactdb_brdfunc
  write(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f
  write(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_ls_user_emissivity
  
end subroutine vlidort_linsup_brdf_f_write

subroutine vlidort_linsup_brdf_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_brdf_def_m, only : vlidort_linsup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_brdf_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_brdf_c_read

subroutine vlidort_linsup_brdf_f_read(lun, fortran_type_f) 
  use vlidort_linsup_brdf_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_brdf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_ls_brdf_f
  read(UNIT=lun) fortran_type_f%ts_ls_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_ls_emissivity
  read(UNIT=lun) fortran_type_f%ts_ls_exactdb_brdfunc
  read(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f
  read(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_ls_user_emissivity
  
end subroutine vlidort_linsup_brdf_f_read

! Links to type: "vlidort_linsup_sleave" from module: "vlidort_linsup_sleave_def_m" in file: "vlidort_lin_sup_sleave_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_sleave_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_sleave_def_m, only : vlidort_linsup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_sleave_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_sleave_c_write

subroutine vlidort_linsup_sleave_f_write(lun, fortran_type_f) 
  use vlidort_linsup_sleave_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_sleave), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_lssl_slterm_f_0
  write(UNIT=lun) fortran_type_f%ts_lssl_slterm_isotropic
  write(UNIT=lun) fortran_type_f%ts_lssl_slterm_userangles
  write(UNIT=lun) fortran_type_f%ts_lssl_user_slterm_f_0
  
end subroutine vlidort_linsup_sleave_f_write

subroutine vlidort_linsup_sleave_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_sleave_def_m, only : vlidort_linsup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_sleave_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_sleave_c_read

subroutine vlidort_linsup_sleave_f_read(lun, fortran_type_f) 
  use vlidort_linsup_sleave_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_sleave), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_lssl_slterm_f_0
  read(UNIT=lun) fortran_type_f%ts_lssl_slterm_isotropic
  read(UNIT=lun) fortran_type_f%ts_lssl_slterm_userangles
  read(UNIT=lun) fortran_type_f%ts_lssl_user_slterm_f_0
  
end subroutine vlidort_linsup_sleave_f_read

! Links to type: "vlidort_linsup_ss_col" from module: "vlidort_linsup_ss_def_m" in file: "vlidort_lin_sup_ss_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_ss_col_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss_col

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_ss_col), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_col_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_col_c_write

subroutine vlidort_linsup_ss_col_f_write(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss_col), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_columnwf_db
  write(UNIT=lun) fortran_type_f%ts_columnwf_ss
  
end subroutine vlidort_linsup_ss_col_f_write

subroutine vlidort_linsup_ss_col_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss_col

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_ss_col), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_col_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_col_c_read

subroutine vlidort_linsup_ss_col_f_read(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss_col), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_columnwf_db
  read(UNIT=lun) fortran_type_f%ts_columnwf_ss
  
end subroutine vlidort_linsup_ss_col_f_read

! Links to type: "vlidort_linsup_ss_prof" from module: "vlidort_linsup_ss_def_m" in file: "vlidort_lin_sup_ss_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_ss_prof_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss_prof

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_ss_prof), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_prof_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_prof_c_write

subroutine vlidort_linsup_ss_prof_f_write(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss_prof), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_profilewf_db
  write(UNIT=lun) fortran_type_f%ts_profilewf_ss
  
end subroutine vlidort_linsup_ss_prof_f_write

subroutine vlidort_linsup_ss_prof_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss_prof

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_ss_prof), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_prof_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_prof_c_read

subroutine vlidort_linsup_ss_prof_f_read(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss_prof), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_profilewf_db
  read(UNIT=lun) fortran_type_f%ts_profilewf_ss
  
end subroutine vlidort_linsup_ss_prof_f_read

! Links to type: "vlidort_linsup_ss_surf" from module: "vlidort_linsup_ss_def_m" in file: "vlidort_lin_sup_ss_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_ss_surf_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss_surf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_ss_surf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_surf_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_surf_c_write

subroutine vlidort_linsup_ss_surf_f_write(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss_surf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_surfacewf_db
  
end subroutine vlidort_linsup_ss_surf_f_write

subroutine vlidort_linsup_ss_surf_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss_surf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_ss_surf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_surf_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_surf_c_read

subroutine vlidort_linsup_ss_surf_f_read(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss_surf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_surfacewf_db
  
end subroutine vlidort_linsup_ss_surf_f_read

! Links to type: "vlidort_linsup_ss" from module: "vlidort_linsup_ss_def_m" in file: "vlidort_lin_sup_ss_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_ss_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_c_write

subroutine vlidort_linsup_ss_f_write(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_linsup_ss_col), pointer :: col_lcl  
  type(vlidort_linsup_ss_prof), pointer :: prof_lcl  
  type(vlidort_linsup_ss_surf), pointer :: surf_lcl  
  
  ! Get pointer to types
  col_lcl => fortran_type_f%col
  prof_lcl => fortran_type_f%prof
  surf_lcl => fortran_type_f%surf
  
  call vlidort_linsup_ss_col_f_write(lun, col_lcl)
  call vlidort_linsup_ss_prof_f_write(lun, prof_lcl)
  call vlidort_linsup_ss_surf_f_write(lun, surf_lcl)
  
end subroutine vlidort_linsup_ss_f_write

subroutine vlidort_linsup_ss_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_ss_def_m, only : vlidort_linsup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_ss_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_ss_c_read

subroutine vlidort_linsup_ss_f_read(lun, fortran_type_f) 
  use vlidort_linsup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_ss), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_linsup_ss_col), pointer :: col_lcl  
  type(vlidort_linsup_ss_prof), pointer :: prof_lcl  
  type(vlidort_linsup_ss_surf), pointer :: surf_lcl  
  
  ! Get pointer to types
  col_lcl => fortran_type_f%col
  prof_lcl => fortran_type_f%prof
  surf_lcl => fortran_type_f%surf
  
  call vlidort_linsup_ss_col_f_read(lun, col_lcl)
  call vlidort_linsup_ss_prof_f_read(lun, prof_lcl)
  call vlidort_linsup_ss_surf_f_read(lun, surf_lcl)
  
end subroutine vlidort_linsup_ss_f_read

! Links to type: "vlidort_linsup_inout" from module: "vlidort_linsup_inout_def_m" in file: "vlidort_lin_sup_def.f90"
! Allocs and initializes type
subroutine vlidort_linsup_inout_c_write(lun, fortran_type_c) bind(C)
  use vlidort_linsup_inout_def_m, only : vlidort_linsup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_linsup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_inout_f_write(lun, fortran_type_f)

end subroutine vlidort_linsup_inout_c_write

subroutine vlidort_linsup_inout_f_write(lun, fortran_type_f) 
  use vlidort_linsup_inout_def_m
  use vlidort_linsup_brdf_def_m
  use vlidort_linsup_ss_def_m
  use vlidort_linsup_sleave_def_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_inout), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_linsup_brdf), pointer :: brdf_lcl  
  type(vlidort_linsup_sleave), pointer :: sleave_lcl  
  type(vlidort_linsup_ss), pointer :: ss_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  sleave_lcl => fortran_type_f%sleave
  ss_lcl => fortran_type_f%ss
  
  call vlidort_linsup_brdf_f_write(lun, brdf_lcl)
  call vlidort_linsup_sleave_f_write(lun, sleave_lcl)
  call vlidort_linsup_ss_f_write(lun, ss_lcl)
  
end subroutine vlidort_linsup_inout_f_write

subroutine vlidort_linsup_inout_c_read(lun, fortran_type_c) bind(C)
  use vlidort_linsup_inout_def_m, only : vlidort_linsup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_linsup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_linsup_inout_f_read(lun, fortran_type_f)

end subroutine vlidort_linsup_inout_c_read

subroutine vlidort_linsup_inout_f_read(lun, fortran_type_f) 
  use vlidort_linsup_inout_def_m
  use vlidort_linsup_brdf_def_m
  use vlidort_linsup_ss_def_m
  use vlidort_linsup_sleave_def_m
  
  integer, intent(in) :: lun
  type(vlidort_linsup_inout), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_linsup_brdf), pointer :: brdf_lcl  
  type(vlidort_linsup_sleave), pointer :: sleave_lcl  
  type(vlidort_linsup_ss), pointer :: ss_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  sleave_lcl => fortran_type_f%sleave
  ss_lcl => fortran_type_f%ss
  
  call vlidort_linsup_brdf_f_read(lun, brdf_lcl)
  call vlidort_linsup_sleave_f_read(lun, sleave_lcl)
  call vlidort_linsup_ss_f_read(lun, ss_lcl)
  
end subroutine vlidort_linsup_inout_f_read

! Links to type: "vlidort_main_outputs" from module: "vlidort_outputs_def_m" in file: "vlidort_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_main_outputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_main_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_main_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_main_outputs_f_write(lun, fortran_type_f)

end subroutine vlidort_main_outputs_c_write

subroutine vlidort_main_outputs_f_write(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_main_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_albmed_fluxes
  write(UNIT=lun) fortran_type_f%ts_albmed_user
  write(UNIT=lun) fortran_type_f%ts_contribs
  write(UNIT=lun) fortran_type_f%ts_dnflux_direct
  write(UNIT=lun) fortran_type_f%ts_dnmeanst_direct
  write(UNIT=lun) fortran_type_f%ts_flux_diffuse
  write(UNIT=lun) fortran_type_f%ts_fourier_saved
  write(UNIT=lun) fortran_type_f%ts_layer_mssts
  write(UNIT=lun) fortran_type_f%ts_lostrans
  write(UNIT=lun) fortran_type_f%ts_meanst_diffuse
  write(UNIT=lun) fortran_type_f%ts_n_geometries
  write(UNIT=lun) fortran_type_f%ts_pathgeoms
  write(UNIT=lun) fortran_type_f%ts_planetary_sbterm
  write(UNIT=lun) fortran_type_f%ts_planetary_transterm
  write(UNIT=lun) fortran_type_f%ts_solarbeam_boatrans
  write(UNIT=lun) fortran_type_f%ts_stokes
  write(UNIT=lun) fortran_type_f%ts_surf_mssts
  write(UNIT=lun) fortran_type_f%ts_sza_offsets
  write(UNIT=lun) fortran_type_f%ts_szd_offsets
  write(UNIT=lun) fortran_type_f%ts_trnmed_fluxes
  write(UNIT=lun) fortran_type_f%ts_trnmed_user
  write(UNIT=lun) fortran_type_f%ts_vza_offsets
  
end subroutine vlidort_main_outputs_f_write

subroutine vlidort_main_outputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_main_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_main_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_main_outputs_f_read(lun, fortran_type_f)

end subroutine vlidort_main_outputs_c_read

subroutine vlidort_main_outputs_f_read(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_main_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_albmed_fluxes
  read(UNIT=lun) fortran_type_f%ts_albmed_user
  read(UNIT=lun) fortran_type_f%ts_contribs
  read(UNIT=lun) fortran_type_f%ts_dnflux_direct
  read(UNIT=lun) fortran_type_f%ts_dnmeanst_direct
  read(UNIT=lun) fortran_type_f%ts_flux_diffuse
  read(UNIT=lun) fortran_type_f%ts_fourier_saved
  read(UNIT=lun) fortran_type_f%ts_layer_mssts
  read(UNIT=lun) fortran_type_f%ts_lostrans
  read(UNIT=lun) fortran_type_f%ts_meanst_diffuse
  read(UNIT=lun) fortran_type_f%ts_n_geometries
  read(UNIT=lun) fortran_type_f%ts_pathgeoms
  read(UNIT=lun) fortran_type_f%ts_planetary_sbterm
  read(UNIT=lun) fortran_type_f%ts_planetary_transterm
  read(UNIT=lun) fortran_type_f%ts_solarbeam_boatrans
  read(UNIT=lun) fortran_type_f%ts_stokes
  read(UNIT=lun) fortran_type_f%ts_surf_mssts
  read(UNIT=lun) fortran_type_f%ts_sza_offsets
  read(UNIT=lun) fortran_type_f%ts_szd_offsets
  read(UNIT=lun) fortran_type_f%ts_trnmed_fluxes
  read(UNIT=lun) fortran_type_f%ts_trnmed_user
  read(UNIT=lun) fortran_type_f%ts_vza_offsets
  
end subroutine vlidort_main_outputs_f_read

! Links to type: "vlidort_wladjusted_outputs" from module: "vlidort_outputs_def_m" in file: "vlidort_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_wladjusted_outputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_wladjusted_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_wladjusted_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_wladjusted_outputs_f_write(lun, fortran_type_f)

end subroutine vlidort_wladjusted_outputs_c_write

subroutine vlidort_wladjusted_outputs_f_write(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_wladjusted_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_wladjusted_direct
  write(UNIT=lun) fortran_type_f%ts_wladjusted_f_ords_0
  write(UNIT=lun) fortran_type_f%ts_wladjusted_f_user_0
  write(UNIT=lun) fortran_type_f%ts_wladjusted_isotropic
  
end subroutine vlidort_wladjusted_outputs_f_write

subroutine vlidort_wladjusted_outputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_wladjusted_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_wladjusted_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_wladjusted_outputs_f_read(lun, fortran_type_f)

end subroutine vlidort_wladjusted_outputs_c_read

subroutine vlidort_wladjusted_outputs_f_read(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_wladjusted_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_wladjusted_direct
  read(UNIT=lun) fortran_type_f%ts_wladjusted_f_ords_0
  read(UNIT=lun) fortran_type_f%ts_wladjusted_f_user_0
  read(UNIT=lun) fortran_type_f%ts_wladjusted_isotropic
  
end subroutine vlidort_wladjusted_outputs_f_read

! Links to type: "vlidort_exception_handling" from module: "vlidort_outputs_def_m" in file: "vlidort_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_exception_handling_f_write(lun, fortran_type_f)

end subroutine vlidort_exception_handling_c_write

subroutine vlidort_exception_handling_f_write(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_actions
  write(UNIT=lun) fortran_type_f%ts_checkmessages
  write(UNIT=lun) fortran_type_f%ts_message
  write(UNIT=lun) fortran_type_f%ts_ncheckmessages
  write(UNIT=lun) fortran_type_f%ts_status_calculation
  write(UNIT=lun) fortran_type_f%ts_status_inputcheck
  write(UNIT=lun) fortran_type_f%ts_trace_1
  write(UNIT=lun) fortran_type_f%ts_trace_2
  write(UNIT=lun) fortran_type_f%ts_trace_3
  
end subroutine vlidort_exception_handling_f_write

subroutine vlidort_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_exception_handling_f_read(lun, fortran_type_f)

end subroutine vlidort_exception_handling_c_read

subroutine vlidort_exception_handling_f_read(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_actions
  read(UNIT=lun) fortran_type_f%ts_checkmessages
  read(UNIT=lun) fortran_type_f%ts_message
  read(UNIT=lun) fortran_type_f%ts_ncheckmessages
  read(UNIT=lun) fortran_type_f%ts_status_calculation
  read(UNIT=lun) fortran_type_f%ts_status_inputcheck
  read(UNIT=lun) fortran_type_f%ts_trace_1
  read(UNIT=lun) fortran_type_f%ts_trace_2
  read(UNIT=lun) fortran_type_f%ts_trace_3
  
end subroutine vlidort_exception_handling_f_read

! Links to type: "vlidort_input_exception_handling" from module: "vlidort_outputs_def_m" in file: "vlidort_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_input_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_input_exception_handling_f_write(lun, fortran_type_f)

end subroutine vlidort_input_exception_handling_c_write

subroutine vlidort_input_exception_handling_f_write(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_input_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_inputactions
  write(UNIT=lun) fortran_type_f%ts_inputmessages
  write(UNIT=lun) fortran_type_f%ts_ninputmessages
  write(UNIT=lun) fortran_type_f%ts_status_inputread
  
end subroutine vlidort_input_exception_handling_f_write

subroutine vlidort_input_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_input_exception_handling_f_read(lun, fortran_type_f)

end subroutine vlidort_input_exception_handling_c_read

subroutine vlidort_input_exception_handling_f_read(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_input_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_inputactions
  read(UNIT=lun) fortran_type_f%ts_inputmessages
  read(UNIT=lun) fortran_type_f%ts_ninputmessages
  read(UNIT=lun) fortran_type_f%ts_status_inputread
  
end subroutine vlidort_input_exception_handling_f_read

! Links to type: "vlidort_outputs" from module: "vlidort_outputs_def_m" in file: "vlidort_outputs_def.f90"
! Allocs and initializes type
subroutine vlidort_outputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_outputs_f_write(lun, fortran_type_f)

end subroutine vlidort_outputs_c_write

subroutine vlidort_outputs_f_write(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_main_outputs), pointer :: main_lcl  
  type(vlidort_exception_handling), pointer :: status_lcl  
  type(vlidort_wladjusted_outputs), pointer :: wlout_lcl  
  
  ! Get pointer to types
  main_lcl => fortran_type_f%main
  status_lcl => fortran_type_f%status
  wlout_lcl => fortran_type_f%wlout
  
  call vlidort_main_outputs_f_write(lun, main_lcl)
  call vlidort_exception_handling_f_write(lun, status_lcl)
  call vlidort_wladjusted_outputs_f_write(lun, wlout_lcl)
  
end subroutine vlidort_outputs_f_write

subroutine vlidort_outputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_outputs_def_m, only : vlidort_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_outputs_f_read(lun, fortran_type_f)

end subroutine vlidort_outputs_c_read

subroutine vlidort_outputs_f_read(lun, fortran_type_f) 
  use vlidort_outputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_main_outputs), pointer :: main_lcl  
  type(vlidort_exception_handling), pointer :: status_lcl  
  type(vlidort_wladjusted_outputs), pointer :: wlout_lcl  
  
  ! Get pointer to types
  main_lcl => fortran_type_f%main
  status_lcl => fortran_type_f%status
  wlout_lcl => fortran_type_f%wlout
  
  call vlidort_main_outputs_f_read(lun, main_lcl)
  call vlidort_exception_handling_f_read(lun, status_lcl)
  call vlidort_wladjusted_outputs_f_read(lun, wlout_lcl)
  
end subroutine vlidort_outputs_f_read

! Links to type: "vlidort_sup_brdf" from module: "vlidort_sup_brdf_def_m" in file: "vlidort_sup_brdf_def.f90"
! Allocs and initializes type
subroutine vlidort_sup_brdf_c_write(lun, fortran_type_c) bind(C)
  use vlidort_sup_brdf_def_m, only : vlidort_sup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_sup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_brdf_f_write(lun, fortran_type_f)

end subroutine vlidort_sup_brdf_c_write

subroutine vlidort_sup_brdf_f_write(lun, fortran_type_f) 
  use vlidort_sup_brdf_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_brdf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_brdf_f
  write(UNIT=lun) fortran_type_f%ts_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_emissivity
  write(UNIT=lun) fortran_type_f%ts_exactdb_brdfunc
  write(UNIT=lun) fortran_type_f%ts_user_brdf_f
  write(UNIT=lun) fortran_type_f%ts_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_user_emissivity
  
end subroutine vlidort_sup_brdf_f_write

subroutine vlidort_sup_brdf_c_read(lun, fortran_type_c) bind(C)
  use vlidort_sup_brdf_def_m, only : vlidort_sup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_sup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_brdf_f_read(lun, fortran_type_f)

end subroutine vlidort_sup_brdf_c_read

subroutine vlidort_sup_brdf_f_read(lun, fortran_type_f) 
  use vlidort_sup_brdf_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_brdf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_brdf_f
  read(UNIT=lun) fortran_type_f%ts_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_emissivity
  read(UNIT=lun) fortran_type_f%ts_exactdb_brdfunc
  read(UNIT=lun) fortran_type_f%ts_user_brdf_f
  read(UNIT=lun) fortran_type_f%ts_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_user_emissivity
  
end subroutine vlidort_sup_brdf_f_read

! Links to type: "vlidort_sup_sleave" from module: "vlidort_sup_sleave_def_m" in file: "vlidort_sup_sleave_def.f90"
! Allocs and initializes type
subroutine vlidort_sup_sleave_c_write(lun, fortran_type_c) bind(C)
  use vlidort_sup_sleave_def_m, only : vlidort_sup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_sup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_sleave_f_write(lun, fortran_type_f)

end subroutine vlidort_sup_sleave_c_write

subroutine vlidort_sup_sleave_f_write(lun, fortran_type_f) 
  use vlidort_sup_sleave_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_sleave), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_slterm_f_0
  write(UNIT=lun) fortran_type_f%ts_slterm_isotropic
  write(UNIT=lun) fortran_type_f%ts_slterm_userangles
  write(UNIT=lun) fortran_type_f%ts_user_slterm_f_0
  
end subroutine vlidort_sup_sleave_f_write

subroutine vlidort_sup_sleave_c_read(lun, fortran_type_c) bind(C)
  use vlidort_sup_sleave_def_m, only : vlidort_sup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_sup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_sleave_f_read(lun, fortran_type_f)

end subroutine vlidort_sup_sleave_c_read

subroutine vlidort_sup_sleave_f_read(lun, fortran_type_f) 
  use vlidort_sup_sleave_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_sleave), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_slterm_f_0
  read(UNIT=lun) fortran_type_f%ts_slterm_isotropic
  read(UNIT=lun) fortran_type_f%ts_slterm_userangles
  read(UNIT=lun) fortran_type_f%ts_user_slterm_f_0
  
end subroutine vlidort_sup_sleave_f_read

! Links to type: "vlidort_sup_ss" from module: "vlidort_sup_ss_def_m" in file: "vlidort_sup_ss_def.f90"
! Allocs and initializes type
subroutine vlidort_sup_ss_c_write(lun, fortran_type_c) bind(C)
  use vlidort_sup_ss_def_m, only : vlidort_sup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_sup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_ss_f_write(lun, fortran_type_f)

end subroutine vlidort_sup_ss_c_write

subroutine vlidort_sup_ss_f_write(lun, fortran_type_f) 
  use vlidort_sup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_ss), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_contribs_ss
  write(UNIT=lun) fortran_type_f%ts_stokes_db
  write(UNIT=lun) fortran_type_f%ts_stokes_ss
  
end subroutine vlidort_sup_ss_f_write

subroutine vlidort_sup_ss_c_read(lun, fortran_type_c) bind(C)
  use vlidort_sup_ss_def_m, only : vlidort_sup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_sup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_ss_f_read(lun, fortran_type_f)

end subroutine vlidort_sup_ss_c_read

subroutine vlidort_sup_ss_f_read(lun, fortran_type_f) 
  use vlidort_sup_ss_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_ss), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_contribs_ss
  read(UNIT=lun) fortran_type_f%ts_stokes_db
  read(UNIT=lun) fortran_type_f%ts_stokes_ss
  
end subroutine vlidort_sup_ss_f_read

! Links to type: "vlidort_sup_inout" from module: "vlidort_sup_inout_def_m" in file: "vlidort_sup_def.f90"
! Allocs and initializes type
subroutine vlidort_sup_inout_c_write(lun, fortran_type_c) bind(C)
  use vlidort_sup_inout_def_m, only : vlidort_sup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_sup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_inout_f_write(lun, fortran_type_f)

end subroutine vlidort_sup_inout_c_write

subroutine vlidort_sup_inout_f_write(lun, fortran_type_f) 
  use vlidort_sup_inout_def_m
  use vlidort_sup_brdf_def_m
  use vlidort_sup_ss_def_m
  use vlidort_sup_sleave_def_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_inout), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_sup_brdf), pointer :: brdf_lcl  
  type(vlidort_sup_sleave), pointer :: sleave_lcl  
  type(vlidort_sup_ss), pointer :: ss_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  sleave_lcl => fortran_type_f%sleave
  ss_lcl => fortran_type_f%ss
  
  call vlidort_sup_brdf_f_write(lun, brdf_lcl)
  call vlidort_sup_sleave_f_write(lun, sleave_lcl)
  call vlidort_sup_ss_f_write(lun, ss_lcl)
  
end subroutine vlidort_sup_inout_f_write

subroutine vlidort_sup_inout_c_read(lun, fortran_type_c) bind(C)
  use vlidort_sup_inout_def_m, only : vlidort_sup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_sup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_sup_inout_f_read(lun, fortran_type_f)

end subroutine vlidort_sup_inout_c_read

subroutine vlidort_sup_inout_f_read(lun, fortran_type_f) 
  use vlidort_sup_inout_def_m
  use vlidort_sup_brdf_def_m
  use vlidort_sup_ss_def_m
  use vlidort_sup_sleave_def_m
  
  integer, intent(in) :: lun
  type(vlidort_sup_inout), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_sup_brdf), pointer :: brdf_lcl  
  type(vlidort_sup_sleave), pointer :: sleave_lcl  
  type(vlidort_sup_ss), pointer :: ss_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  sleave_lcl => fortran_type_f%sleave
  ss_lcl => fortran_type_f%ss
  
  call vlidort_sup_brdf_f_read(lun, brdf_lcl)
  call vlidort_sup_sleave_f_read(lun, sleave_lcl)
  call vlidort_sup_ss_f_read(lun, ss_lcl)
  
end subroutine vlidort_sup_inout_f_read

! Links to type: "vlidort_fixed_boolean" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_boolean_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_boolean_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_boolean_c_write

subroutine vlidort_fixed_boolean_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_boolean), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_albtrn_media
  write(UNIT=lun) fortran_type_f%ts_do_boa_illumination
  write(UNIT=lun) fortran_type_f%ts_do_dnwelling
  write(UNIT=lun) fortran_type_f%ts_do_fluorescence
  write(UNIT=lun) fortran_type_f%ts_do_fourier0_nstokes2
  write(UNIT=lun) fortran_type_f%ts_do_fullrad_mode
  write(UNIT=lun) fortran_type_f%ts_do_lambertian_surface
  write(UNIT=lun) fortran_type_f%ts_do_mssts
  write(UNIT=lun) fortran_type_f%ts_do_plane_parallel
  write(UNIT=lun) fortran_type_f%ts_do_planetary_problem
  write(UNIT=lun) fortran_type_f%ts_do_sl_isotropic
  write(UNIT=lun) fortran_type_f%ts_do_specialist_option_1
  write(UNIT=lun) fortran_type_f%ts_do_specialist_option_2
  write(UNIT=lun) fortran_type_f%ts_do_specialist_option_3
  write(UNIT=lun) fortran_type_f%ts_do_surface_emission
  write(UNIT=lun) fortran_type_f%ts_do_surface_leaving
  write(UNIT=lun) fortran_type_f%ts_do_tf_iteration
  write(UNIT=lun) fortran_type_f%ts_do_thermal_emission
  write(UNIT=lun) fortran_type_f%ts_do_toa_contribs
  write(UNIT=lun) fortran_type_f%ts_do_toa_illumination
  write(UNIT=lun) fortran_type_f%ts_do_upwelling
  write(UNIT=lun) fortran_type_f%ts_do_water_leaving
  write(UNIT=lun) fortran_type_f%ts_do_wladjusted_output
  
end subroutine vlidort_fixed_boolean_f_write

subroutine vlidort_fixed_boolean_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_boolean_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_boolean_c_read

subroutine vlidort_fixed_boolean_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_boolean), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_albtrn_media
  read(UNIT=lun) fortran_type_f%ts_do_boa_illumination
  read(UNIT=lun) fortran_type_f%ts_do_dnwelling
  read(UNIT=lun) fortran_type_f%ts_do_fluorescence
  read(UNIT=lun) fortran_type_f%ts_do_fourier0_nstokes2
  read(UNIT=lun) fortran_type_f%ts_do_fullrad_mode
  read(UNIT=lun) fortran_type_f%ts_do_lambertian_surface
  read(UNIT=lun) fortran_type_f%ts_do_mssts
  read(UNIT=lun) fortran_type_f%ts_do_plane_parallel
  read(UNIT=lun) fortran_type_f%ts_do_planetary_problem
  read(UNIT=lun) fortran_type_f%ts_do_sl_isotropic
  read(UNIT=lun) fortran_type_f%ts_do_specialist_option_1
  read(UNIT=lun) fortran_type_f%ts_do_specialist_option_2
  read(UNIT=lun) fortran_type_f%ts_do_specialist_option_3
  read(UNIT=lun) fortran_type_f%ts_do_surface_emission
  read(UNIT=lun) fortran_type_f%ts_do_surface_leaving
  read(UNIT=lun) fortran_type_f%ts_do_tf_iteration
  read(UNIT=lun) fortran_type_f%ts_do_thermal_emission
  read(UNIT=lun) fortran_type_f%ts_do_toa_contribs
  read(UNIT=lun) fortran_type_f%ts_do_toa_illumination
  read(UNIT=lun) fortran_type_f%ts_do_upwelling
  read(UNIT=lun) fortran_type_f%ts_do_water_leaving
  read(UNIT=lun) fortran_type_f%ts_do_wladjusted_output
  
end subroutine vlidort_fixed_boolean_f_read

! Links to type: "vlidort_fixed_control" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_control_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_control_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_control_c_write

subroutine vlidort_fixed_control_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_control), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_asymtx_tolerance
  write(UNIT=lun) fortran_type_f%ts_boa_illumination
  write(UNIT=lun) fortran_type_f%ts_n_thermal_coeffs
  write(UNIT=lun) fortran_type_f%ts_nfinelayers
  write(UNIT=lun) fortran_type_f%ts_nlayers
  write(UNIT=lun) fortran_type_f%ts_nlayers_cutoff
  write(UNIT=lun) fortran_type_f%ts_nlayers_noms
  write(UNIT=lun) fortran_type_f%ts_nstokes
  write(UNIT=lun) fortran_type_f%ts_nstreams
  write(UNIT=lun) fortran_type_f%ts_taylor_order
  write(UNIT=lun) fortran_type_f%ts_tf_criterion
  write(UNIT=lun) fortran_type_f%ts_tf_maxiter
  write(UNIT=lun) fortran_type_f%ts_toa_illumination
  write(UNIT=lun) fortran_type_f%ts_vlidort_accuracy
  
end subroutine vlidort_fixed_control_f_write

subroutine vlidort_fixed_control_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_control_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_control_c_read

subroutine vlidort_fixed_control_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_control), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_asymtx_tolerance
  read(UNIT=lun) fortran_type_f%ts_boa_illumination
  read(UNIT=lun) fortran_type_f%ts_n_thermal_coeffs
  read(UNIT=lun) fortran_type_f%ts_nfinelayers
  read(UNIT=lun) fortran_type_f%ts_nlayers
  read(UNIT=lun) fortran_type_f%ts_nlayers_cutoff
  read(UNIT=lun) fortran_type_f%ts_nlayers_noms
  read(UNIT=lun) fortran_type_f%ts_nstokes
  read(UNIT=lun) fortran_type_f%ts_nstreams
  read(UNIT=lun) fortran_type_f%ts_taylor_order
  read(UNIT=lun) fortran_type_f%ts_tf_criterion
  read(UNIT=lun) fortran_type_f%ts_tf_maxiter
  read(UNIT=lun) fortran_type_f%ts_toa_illumination
  read(UNIT=lun) fortran_type_f%ts_vlidort_accuracy
  
end subroutine vlidort_fixed_control_f_read

! Links to type: "vlidort_fixed_sunrays" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_sunrays_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_sunrays_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_sunrays_c_write

subroutine vlidort_fixed_sunrays_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_sunrays), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_flux_factor
  
end subroutine vlidort_fixed_sunrays_f_write

subroutine vlidort_fixed_sunrays_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_sunrays_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_sunrays_c_read

subroutine vlidort_fixed_sunrays_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_sunrays), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_flux_factor
  
end subroutine vlidort_fixed_sunrays_f_read

! Links to type: "vlidort_fixed_uservalues" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_uservalues_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_uservalues_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_uservalues_c_write

subroutine vlidort_fixed_uservalues_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_uservalues), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_n_user_levels
  
end subroutine vlidort_fixed_uservalues_f_write

subroutine vlidort_fixed_uservalues_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_uservalues_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_uservalues_c_read

subroutine vlidort_fixed_uservalues_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_uservalues), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_n_user_levels
  
end subroutine vlidort_fixed_uservalues_f_read

! Links to type: "vlidort_fixed_chapman" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_chapman_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_chapman_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_chapman_c_write

subroutine vlidort_fixed_chapman_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_chapman), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_finegrid
  write(UNIT=lun) fortran_type_f%ts_height_grid
  write(UNIT=lun) fortran_type_f%ts_pressure_grid
  write(UNIT=lun) fortran_type_f%ts_rfindex_parameter
  write(UNIT=lun) fortran_type_f%ts_temperature_grid
  
end subroutine vlidort_fixed_chapman_f_write

subroutine vlidort_fixed_chapman_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_chapman_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_chapman_c_read

subroutine vlidort_fixed_chapman_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_chapman), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_finegrid
  read(UNIT=lun) fortran_type_f%ts_height_grid
  read(UNIT=lun) fortran_type_f%ts_pressure_grid
  read(UNIT=lun) fortran_type_f%ts_rfindex_parameter
  read(UNIT=lun) fortran_type_f%ts_temperature_grid
  
end subroutine vlidort_fixed_chapman_f_read

! Links to type: "vlidort_fixed_optical" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_optical_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_optical_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_optical_c_write

subroutine vlidort_fixed_optical_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_optical), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_atmos_wavelength
  write(UNIT=lun) fortran_type_f%ts_deltau_vert_input
  write(UNIT=lun) fortran_type_f%ts_fmatrix_dn
  write(UNIT=lun) fortran_type_f%ts_fmatrix_up
  write(UNIT=lun) fortran_type_f%ts_greekmat_total_input
  write(UNIT=lun) fortran_type_f%ts_lambertian_albedo
  write(UNIT=lun) fortran_type_f%ts_surface_bb_input
  write(UNIT=lun) fortran_type_f%ts_thermal_bb_input
  
end subroutine vlidort_fixed_optical_f_write

subroutine vlidort_fixed_optical_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_optical_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_optical_c_read

subroutine vlidort_fixed_optical_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_optical), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_atmos_wavelength
  read(UNIT=lun) fortran_type_f%ts_deltau_vert_input
  read(UNIT=lun) fortran_type_f%ts_fmatrix_dn
  read(UNIT=lun) fortran_type_f%ts_fmatrix_up
  read(UNIT=lun) fortran_type_f%ts_greekmat_total_input
  read(UNIT=lun) fortran_type_f%ts_lambertian_albedo
  read(UNIT=lun) fortran_type_f%ts_surface_bb_input
  read(UNIT=lun) fortran_type_f%ts_thermal_bb_input
  
end subroutine vlidort_fixed_optical_f_read

! Links to type: "vlidort_fixed_write" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_write_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_write

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_write), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_write_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_write_c_write

subroutine vlidort_fixed_write_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_write), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_debug_write
  write(UNIT=lun) fortran_type_f%ts_do_write_fourier
  write(UNIT=lun) fortran_type_f%ts_do_write_input
  write(UNIT=lun) fortran_type_f%ts_do_write_results
  write(UNIT=lun) fortran_type_f%ts_do_write_scenario
  write(UNIT=lun) fortran_type_f%ts_fourier_write_filename
  write(UNIT=lun) fortran_type_f%ts_input_write_filename
  write(UNIT=lun) fortran_type_f%ts_results_write_filename
  write(UNIT=lun) fortran_type_f%ts_scenario_write_filename
  
end subroutine vlidort_fixed_write_f_write

subroutine vlidort_fixed_write_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_write

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_write), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_write_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_write_c_read

subroutine vlidort_fixed_write_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_write), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_debug_write
  read(UNIT=lun) fortran_type_f%ts_do_write_fourier
  read(UNIT=lun) fortran_type_f%ts_do_write_input
  read(UNIT=lun) fortran_type_f%ts_do_write_results
  read(UNIT=lun) fortran_type_f%ts_do_write_scenario
  read(UNIT=lun) fortran_type_f%ts_fourier_write_filename
  read(UNIT=lun) fortran_type_f%ts_input_write_filename
  read(UNIT=lun) fortran_type_f%ts_results_write_filename
  read(UNIT=lun) fortran_type_f%ts_scenario_write_filename
  
end subroutine vlidort_fixed_write_f_read

! Links to type: "vlidort_fixed_inputs" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_fixed_inputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_fixed_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_inputs_f_write(lun, fortran_type_f)

end subroutine vlidort_fixed_inputs_c_write

subroutine vlidort_fixed_inputs_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_fixed_boolean), pointer :: bool_lcl  
  type(vlidort_fixed_chapman), pointer :: chapman_lcl  
  type(vlidort_fixed_control), pointer :: cont_lcl  
  type(vlidort_fixed_optical), pointer :: optical_lcl  
  type(vlidort_fixed_sunrays), pointer :: sunrays_lcl  
  type(vlidort_fixed_uservalues), pointer :: userval_lcl  
  type(vlidort_fixed_write), pointer :: write_lcl  
  
  ! Get pointer to types
  bool_lcl => fortran_type_f%bool
  chapman_lcl => fortran_type_f%chapman
  cont_lcl => fortran_type_f%cont
  optical_lcl => fortran_type_f%optical
  sunrays_lcl => fortran_type_f%sunrays
  userval_lcl => fortran_type_f%userval
  write_lcl => fortran_type_f%write
  
  call vlidort_fixed_boolean_f_write(lun, bool_lcl)
  call vlidort_fixed_chapman_f_write(lun, chapman_lcl)
  call vlidort_fixed_control_f_write(lun, cont_lcl)
  call vlidort_fixed_optical_f_write(lun, optical_lcl)
  call vlidort_fixed_sunrays_f_write(lun, sunrays_lcl)
  call vlidort_fixed_uservalues_f_write(lun, userval_lcl)
  call vlidort_fixed_write_f_write(lun, write_lcl)
  
end subroutine vlidort_fixed_inputs_f_write

subroutine vlidort_fixed_inputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_fixed_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_fixed_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_fixed_inputs_f_read(lun, fortran_type_f)

end subroutine vlidort_fixed_inputs_c_read

subroutine vlidort_fixed_inputs_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_fixed_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_fixed_boolean), pointer :: bool_lcl  
  type(vlidort_fixed_chapman), pointer :: chapman_lcl  
  type(vlidort_fixed_control), pointer :: cont_lcl  
  type(vlidort_fixed_optical), pointer :: optical_lcl  
  type(vlidort_fixed_sunrays), pointer :: sunrays_lcl  
  type(vlidort_fixed_uservalues), pointer :: userval_lcl  
  type(vlidort_fixed_write), pointer :: write_lcl  
  
  ! Get pointer to types
  bool_lcl => fortran_type_f%bool
  chapman_lcl => fortran_type_f%chapman
  cont_lcl => fortran_type_f%cont
  optical_lcl => fortran_type_f%optical
  sunrays_lcl => fortran_type_f%sunrays
  userval_lcl => fortran_type_f%userval
  write_lcl => fortran_type_f%write
  
  call vlidort_fixed_boolean_f_read(lun, bool_lcl)
  call vlidort_fixed_chapman_f_read(lun, chapman_lcl)
  call vlidort_fixed_control_f_read(lun, cont_lcl)
  call vlidort_fixed_optical_f_read(lun, optical_lcl)
  call vlidort_fixed_sunrays_f_read(lun, sunrays_lcl)
  call vlidort_fixed_uservalues_f_read(lun, userval_lcl)
  call vlidort_fixed_write_f_read(lun, write_lcl)
  
end subroutine vlidort_fixed_inputs_f_read

! Links to type: "vlidort_modified_boolean" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_boolean_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_boolean_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_boolean_c_write

subroutine vlidort_modified_boolean_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_boolean), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_additional_mvout
  write(UNIT=lun) fortran_type_f%ts_do_bvp_telescoping
  write(UNIT=lun) fortran_type_f%ts_do_chapman_function
  write(UNIT=lun) fortran_type_f%ts_do_classical_solution
  write(UNIT=lun) fortran_type_f%ts_do_deltam_scaling
  write(UNIT=lun) fortran_type_f%ts_do_double_convtest
  write(UNIT=lun) fortran_type_f%ts_do_doublet_geometry
  write(UNIT=lun) fortran_type_f%ts_do_external_wleave
  write(UNIT=lun) fortran_type_f%ts_do_focorr
  write(UNIT=lun) fortran_type_f%ts_do_focorr_external
  write(UNIT=lun) fortran_type_f%ts_do_focorr_nadir
  write(UNIT=lun) fortran_type_f%ts_do_focorr_outgoing
  write(UNIT=lun) fortran_type_f%ts_do_mvout_only
  write(UNIT=lun) fortran_type_f%ts_do_observation_geometry
  write(UNIT=lun) fortran_type_f%ts_do_rayleigh_only
  write(UNIT=lun) fortran_type_f%ts_do_refractive_geometry
  write(UNIT=lun) fortran_type_f%ts_do_solar_sources
  write(UNIT=lun) fortran_type_f%ts_do_solution_saving
  write(UNIT=lun) fortran_type_f%ts_do_sscorr_truncation
  write(UNIT=lun) fortran_type_f%ts_do_sscorr_usefmat
  write(UNIT=lun) fortran_type_f%ts_do_thermal_transonly
  write(UNIT=lun) fortran_type_f%ts_do_user_vzangles
  
end subroutine vlidort_modified_boolean_f_write

subroutine vlidort_modified_boolean_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_boolean_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_boolean_c_read

subroutine vlidort_modified_boolean_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_boolean), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_additional_mvout
  read(UNIT=lun) fortran_type_f%ts_do_bvp_telescoping
  read(UNIT=lun) fortran_type_f%ts_do_chapman_function
  read(UNIT=lun) fortran_type_f%ts_do_classical_solution
  read(UNIT=lun) fortran_type_f%ts_do_deltam_scaling
  read(UNIT=lun) fortran_type_f%ts_do_double_convtest
  read(UNIT=lun) fortran_type_f%ts_do_doublet_geometry
  read(UNIT=lun) fortran_type_f%ts_do_external_wleave
  read(UNIT=lun) fortran_type_f%ts_do_focorr
  read(UNIT=lun) fortran_type_f%ts_do_focorr_external
  read(UNIT=lun) fortran_type_f%ts_do_focorr_nadir
  read(UNIT=lun) fortran_type_f%ts_do_focorr_outgoing
  read(UNIT=lun) fortran_type_f%ts_do_mvout_only
  read(UNIT=lun) fortran_type_f%ts_do_observation_geometry
  read(UNIT=lun) fortran_type_f%ts_do_rayleigh_only
  read(UNIT=lun) fortran_type_f%ts_do_refractive_geometry
  read(UNIT=lun) fortran_type_f%ts_do_solar_sources
  read(UNIT=lun) fortran_type_f%ts_do_solution_saving
  read(UNIT=lun) fortran_type_f%ts_do_sscorr_truncation
  read(UNIT=lun) fortran_type_f%ts_do_sscorr_usefmat
  read(UNIT=lun) fortran_type_f%ts_do_thermal_transonly
  read(UNIT=lun) fortran_type_f%ts_do_user_vzangles
  
end subroutine vlidort_modified_boolean_f_read

! Links to type: "vlidort_modified_control" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_control_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_control_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_control_c_write

subroutine vlidort_modified_control_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_control), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_ngreek_moments_input
  
end subroutine vlidort_modified_control_f_write

subroutine vlidort_modified_control_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_control_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_control_c_read

subroutine vlidort_modified_control_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_control), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_ngreek_moments_input
  
end subroutine vlidort_modified_control_f_read

! Links to type: "vlidort_modified_sunrays" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_sunrays_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_sunrays_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_sunrays_c_write

subroutine vlidort_modified_sunrays_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_sunrays), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_n_szangles
  write(UNIT=lun) fortran_type_f%ts_szangles
  
end subroutine vlidort_modified_sunrays_f_write

subroutine vlidort_modified_sunrays_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_sunrays_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_sunrays_c_read

subroutine vlidort_modified_sunrays_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_sunrays), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_n_szangles
  read(UNIT=lun) fortran_type_f%ts_szangles
  
end subroutine vlidort_modified_sunrays_f_read

! Links to type: "vlidort_modified_uservalues" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_uservalues_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_uservalues_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_uservalues_c_write

subroutine vlidort_modified_uservalues_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_uservalues), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_geometry_specheight
  write(UNIT=lun) fortran_type_f%ts_n_user_doublets
  write(UNIT=lun) fortran_type_f%ts_n_user_obsgeoms
  write(UNIT=lun) fortran_type_f%ts_n_user_relazms
  write(UNIT=lun) fortran_type_f%ts_n_user_vzangles
  write(UNIT=lun) fortran_type_f%ts_user_doublets
  write(UNIT=lun) fortran_type_f%ts_user_levels
  write(UNIT=lun) fortran_type_f%ts_user_obsgeoms_input
  write(UNIT=lun) fortran_type_f%ts_user_relazms
  write(UNIT=lun) fortran_type_f%ts_user_vzangles_input
  
end subroutine vlidort_modified_uservalues_f_write

subroutine vlidort_modified_uservalues_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_uservalues_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_uservalues_c_read

subroutine vlidort_modified_uservalues_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_uservalues), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_geometry_specheight
  read(UNIT=lun) fortran_type_f%ts_n_user_doublets
  read(UNIT=lun) fortran_type_f%ts_n_user_obsgeoms
  read(UNIT=lun) fortran_type_f%ts_n_user_relazms
  read(UNIT=lun) fortran_type_f%ts_n_user_vzangles
  read(UNIT=lun) fortran_type_f%ts_user_doublets
  read(UNIT=lun) fortran_type_f%ts_user_levels
  read(UNIT=lun) fortran_type_f%ts_user_obsgeoms_input
  read(UNIT=lun) fortran_type_f%ts_user_relazms
  read(UNIT=lun) fortran_type_f%ts_user_vzangles_input
  
end subroutine vlidort_modified_uservalues_f_read

! Links to type: "vlidort_modified_chapman" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_chapman_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_chapman_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_chapman_c_write

subroutine vlidort_modified_chapman_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_chapman), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_earth_radius
  
end subroutine vlidort_modified_chapman_f_write

subroutine vlidort_modified_chapman_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_chapman_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_chapman_c_read

subroutine vlidort_modified_chapman_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_chapman), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_earth_radius
  
end subroutine vlidort_modified_chapman_f_read

! Links to type: "vlidort_modified_optical" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_optical_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_optical_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_optical_c_write

subroutine vlidort_modified_optical_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_optical), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_omega_total_input
  
end subroutine vlidort_modified_optical_f_write

subroutine vlidort_modified_optical_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_optical_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_optical_c_read

subroutine vlidort_modified_optical_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_optical), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_omega_total_input
  
end subroutine vlidort_modified_optical_f_read

! Links to type: "vlidort_modified_inputs" from module: "vlidort_inputs_def_m" in file: "vlidort_inputs_def.f90"
! Allocs and initializes type
subroutine vlidort_modified_inputs_c_write(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(vlidort_modified_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_inputs_f_write(lun, fortran_type_f)

end subroutine vlidort_modified_inputs_c_write

subroutine vlidort_modified_inputs_f_write(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_modified_boolean), pointer :: mbool_lcl  
  type(vlidort_modified_chapman), pointer :: mchapman_lcl  
  type(vlidort_modified_control), pointer :: mcont_lcl  
  type(vlidort_modified_optical), pointer :: moptical_lcl  
  type(vlidort_modified_sunrays), pointer :: msunrays_lcl  
  type(vlidort_modified_uservalues), pointer :: muserval_lcl  
  
  ! Get pointer to types
  mbool_lcl => fortran_type_f%mbool
  mchapman_lcl => fortran_type_f%mchapman
  mcont_lcl => fortran_type_f%mcont
  moptical_lcl => fortran_type_f%moptical
  msunrays_lcl => fortran_type_f%msunrays
  muserval_lcl => fortran_type_f%muserval
  
  call vlidort_modified_boolean_f_write(lun, mbool_lcl)
  call vlidort_modified_chapman_f_write(lun, mchapman_lcl)
  call vlidort_modified_control_f_write(lun, mcont_lcl)
  call vlidort_modified_optical_f_write(lun, moptical_lcl)
  call vlidort_modified_sunrays_f_write(lun, msunrays_lcl)
  call vlidort_modified_uservalues_f_write(lun, muserval_lcl)
  
end subroutine vlidort_modified_inputs_f_write

subroutine vlidort_modified_inputs_c_read(lun, fortran_type_c) bind(C)
  use vlidort_inputs_def_m, only : vlidort_modified_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(vlidort_modified_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call vlidort_modified_inputs_f_read(lun, fortran_type_f)

end subroutine vlidort_modified_inputs_c_read

subroutine vlidort_modified_inputs_f_read(lun, fortran_type_f) 
  use vlidort_inputs_def_m
  use vlidort_pars_m
  
  integer, intent(in) :: lun
  type(vlidort_modified_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(vlidort_modified_boolean), pointer :: mbool_lcl  
  type(vlidort_modified_chapman), pointer :: mchapman_lcl  
  type(vlidort_modified_control), pointer :: mcont_lcl  
  type(vlidort_modified_optical), pointer :: moptical_lcl  
  type(vlidort_modified_sunrays), pointer :: msunrays_lcl  
  type(vlidort_modified_uservalues), pointer :: muserval_lcl  
  
  ! Get pointer to types
  mbool_lcl => fortran_type_f%mbool
  mchapman_lcl => fortran_type_f%mchapman
  mcont_lcl => fortran_type_f%mcont
  moptical_lcl => fortran_type_f%moptical
  msunrays_lcl => fortran_type_f%msunrays
  muserval_lcl => fortran_type_f%muserval
  
  call vlidort_modified_boolean_f_read(lun, mbool_lcl)
  call vlidort_modified_chapman_f_read(lun, mchapman_lcl)
  call vlidort_modified_control_f_read(lun, mcont_lcl)
  call vlidort_modified_optical_f_read(lun, moptical_lcl)
  call vlidort_modified_sunrays_f_read(lun, msunrays_lcl)
  call vlidort_modified_uservalues_f_read(lun, muserval_lcl)
  
end subroutine vlidort_modified_inputs_f_read


end module vlidort_interface_types_io