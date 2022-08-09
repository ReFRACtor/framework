
module VBRDF_LINSUP_MASTERS_M_WRAP

use iso_c_binding
use vbrdf_linsup_masters_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine vbrdf_linsup_masters_m_vbrdf_lin_inputmaster_wrap (filnam_in_len, &
                                                              filnam_in, &
                                                              vbrdf_sup_in_in, &
                                                              vbrdf_linsup_in_in, &
                                                              vbrdf_sup_inputstatus_in) bind(C)
  use vlidort_pars_m
  use vbrdf_findpar_m
  use vbrdf_sup_inputs_def_m
  use vbrdf_sup_outputs_def_m
  use vbrdf_linsup_inputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_linsup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_linsup_inputs), pointer :: vbrdf_linsup_in_lcl
  type(vbrdf_input_exception_handling), pointer :: vbrdf_sup_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_linsup_in_in, vbrdf_linsup_in_lcl)
  call c_f_pointer(vbrdf_sup_inputstatus_in, vbrdf_sup_inputstatus_lcl)

  call vbrdf_lin_inputmaster(filnam_lcl, &
                             vbrdf_sup_in_lcl, &
                             vbrdf_linsup_in_lcl, &
                             vbrdf_sup_inputstatus_lcl)

end subroutine vbrdf_linsup_masters_m_vbrdf_lin_inputmaster_wrap

subroutine vbrdf_linsup_masters_m_vbrdf_lin_mainmaster_wrap (do_debug_restoration_in, &
                                                             nmoments_input_in, &
                                                             vbrdf_sup_in_in, &
                                                             vbrdf_linsup_in_in, &
                                                             vbrdf_sup_out_in, &
                                                             vbrdf_linsup_out_in, &
                                                             vbrdf_sup_outputstatus_in) bind(C)
  use vlidort_pars_m
  use vbrdf_sup_inputs_def_m
  use vbrdf_sup_outputs_def_m
  use vbrdf_linsup_inputs_def_m
  use vbrdf_linsup_outputs_def_m
  use vbrdf_sup_aux_m
  use vbrdf_sup_kernels_m
  use vbrdf_linsup_kernels_m
  use vbrdf_sup_routines_m
  use vbrdf_linsup_routines_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_restoration_in
  integer(c_int), intent(in) :: nmoments_input_in
  type(c_ptr), intent(in) :: vbrdf_sup_in_in
  type(c_ptr), intent(in) :: vbrdf_linsup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_out_in
  type(c_ptr), intent(out) :: vbrdf_linsup_out_in
  type(c_ptr), intent(out) :: vbrdf_sup_outputstatus_in

  ! Local variables
  logical(kind=4) :: do_debug_restoration_lcl
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_linsup_inputs), pointer :: vbrdf_linsup_in_lcl
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vbrdf_linsup_outputs), pointer :: vbrdf_linsup_out_lcl
  type(vbrdf_output_exception_handling), pointer :: vbrdf_sup_outputstatus_lcl

  ! Convert input arguments
  do_debug_restoration_lcl = do_debug_restoration_in
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_linsup_in_in, vbrdf_linsup_in_lcl)
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vbrdf_linsup_out_in, vbrdf_linsup_out_lcl)
  call c_f_pointer(vbrdf_sup_outputstatus_in, vbrdf_sup_outputstatus_lcl)

  call vbrdf_lin_mainmaster(do_debug_restoration_lcl, &
                            nmoments_input_in, &
                            vbrdf_sup_in_lcl, &
                            vbrdf_linsup_in_lcl, &
                            vbrdf_sup_out_lcl, &
                            vbrdf_linsup_out_lcl, &
                            vbrdf_sup_outputstatus_lcl)

end subroutine vbrdf_linsup_masters_m_vbrdf_lin_mainmaster_wrap

end module VBRDF_LINSUP_MASTERS_M_WRAP

module VBRDF_SUP_MASTERS_M_WRAP

use iso_c_binding
use vbrdf_sup_masters_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine vbrdf_sup_masters_m_vbrdf_inputmaster_wrap (filnam_in_len, &
                                                       filnam_in, &
                                                       vbrdf_sup_in_in, &
                                                       vbrdf_sup_inputstatus_in) bind(C)
  use vlidort_pars_m
  use vbrdf_findpar_m
  use vbrdf_sup_inputs_def_m
  use vbrdf_sup_outputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_input_exception_handling), pointer :: vbrdf_sup_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_sup_inputstatus_in, vbrdf_sup_inputstatus_lcl)

  call vbrdf_inputmaster(filnam_lcl, &
                         vbrdf_sup_in_lcl, &
                         vbrdf_sup_inputstatus_lcl)

end subroutine vbrdf_sup_masters_m_vbrdf_inputmaster_wrap

subroutine vbrdf_sup_masters_m_vbrdf_mainmaster_wrap (do_debug_restoration_in, &
                                                      nmoments_input_in, &
                                                      vbrdf_sup_in_in, &
                                                      vbrdf_sup_out_in, &
                                                      vbrdf_sup_outputstatus_in) bind(C)
  use vlidort_pars_m
  use vbrdf_sup_inputs_def_m
  use vbrdf_sup_outputs_def_m
  use vbrdf_sup_aux_m
  use vbrdf_sup_kernels_m
  use vbrdf_sup_routines_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_restoration_in
  integer(c_int), intent(in) :: nmoments_input_in
  type(c_ptr), intent(in) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_out_in
  type(c_ptr), intent(out) :: vbrdf_sup_outputstatus_in

  ! Local variables
  logical(kind=4) :: do_debug_restoration_lcl
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vbrdf_output_exception_handling), pointer :: vbrdf_sup_outputstatus_lcl

  ! Convert input arguments
  do_debug_restoration_lcl = do_debug_restoration_in
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vbrdf_sup_outputstatus_in, vbrdf_sup_outputstatus_lcl)

  call vbrdf_mainmaster(do_debug_restoration_lcl, &
                        nmoments_input_in, &
                        vbrdf_sup_in_lcl, &
                        vbrdf_sup_out_lcl, &
                        vbrdf_sup_outputstatus_lcl)

end subroutine vbrdf_sup_masters_m_vbrdf_mainmaster_wrap

end module VBRDF_SUP_MASTERS_M_WRAP

module VLIDORT_INPUTS_M_WRAP

use iso_c_binding
use vlidort_inputs_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine v_inputs_m_v_brdf_sup_init_wrap (vlidort_sup_in) bind(C)
  use vlidort_pars_m
  use vlidort_sup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_sup_in

  ! Local variables
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)

  call vlidort_brdf_sup_init(vlidort_sup_lcl)

end subroutine v_inputs_m_v_brdf_sup_init_wrap

subroutine v_inputs_m_v_input_master_wrap (filnam_in_len, &
                                           filnam_in, &
                                           vlidort_fixin_in, &
                                           vlidort_modin_in, &
                                           vlidort_inputstatus_in) bind(C)
  use vlidort_pars_m
  use vlidort_inputs_def_m
  use vlidort_outputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: vlidort_fixin_in
  type(c_ptr), intent(out) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_input_exception_handling), pointer :: vlidort_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_inputstatus_in, vlidort_inputstatus_lcl)

  call vlidort_input_master(filnam_lcl, &
                            vlidort_fixin_lcl, &
                            vlidort_modin_lcl, &
                            vlidort_inputstatus_lcl)

end subroutine v_inputs_m_v_input_master_wrap

subroutine v_inputs_m_v_sleave_sup_init_wrap (vlidort_sup_in) bind(C)
  use vlidort_pars_m
  use vlidort_sup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_sup_in

  ! Local variables
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)

  call vlidort_sleave_sup_init(vlidort_sup_lcl)

end subroutine v_inputs_m_v_sleave_sup_init_wrap

subroutine v_inputs_m_v_ss_sup_init_wrap (vlidort_sup_in) bind(C)
  use vlidort_pars_m
  use vlidort_sup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_sup_in

  ! Local variables
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)

  call vlidort_ss_sup_init(vlidort_sup_lcl)

end subroutine v_inputs_m_v_ss_sup_init_wrap

subroutine v_inputs_m_v_sup_init_wrap (vlidort_sup_in) bind(C)
  use vlidort_sup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_sup_in

  ! Local variables
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)

  call vlidort_sup_init(vlidort_sup_lcl)

end subroutine v_inputs_m_v_sup_init_wrap

end module VLIDORT_INPUTS_M_WRAP

module VLIDORT_MASTERS_M_WRAP

use iso_c_binding
use vlidort_masters_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine v_masters_m_v_master_wrap (do_debug_input_in, &
                                      vlidort_fixin_in, &
                                      vlidort_modin_in, &
                                      vlidort_sup_in, &
                                      vlidort_out_in) bind(C)
  use vlidort_pars_m
  use vlidort_inputs_def_m
  use vlidort_sup_inout_def_m
  use vlidort_outputs_def_m
  use vlidort_work_def_m
  use vlidort_inputs_m
  use vlidort_geometry_m
  use vlidort_miscsetups_m
  use vlidort_thermalsup_m
  use vlidort_converge_m
  use vlidort_pack_m
  use vlidort_writemodules_m
  use vlidort_vfo_interface_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_input_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in

  ! Local variables
  logical(kind=4) :: do_debug_input_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl

  ! Convert input arguments
  do_debug_input_lcl = do_debug_input_in
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)

  call vlidort_master(do_debug_input_lcl, &
                      vlidort_fixin_lcl, &
                      vlidort_modin_lcl, &
                      vlidort_sup_lcl, &
                      vlidort_out_lcl)

end subroutine v_masters_m_v_master_wrap

end module VLIDORT_MASTERS_M_WRAP

module VLIDORT_L_INPUTS_M_WRAP

use iso_c_binding
use vlidort_l_inputs_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine v_l_inputs_m_v_brdf_linsup_init_wrap (vlidort_linsup_in) bind(C)
  use vlidort_pars_m
  use vlidort_linsup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_linsup_in

  ! Local variables
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)

  call vlidort_brdf_linsup_init(vlidort_linsup_lcl)

end subroutine v_l_inputs_m_v_brdf_linsup_init_wrap

subroutine v_l_inputs_m_v_l_input_master_wrap (filnam_in_len, &
                                               filnam_in, &
                                               vlidort_fixin_in, &
                                               vlidort_modin_in, &
                                               vlidort_linfixin_in, &
                                               vlidort_linmodin_in, &
                                               vlidort_inputstatus_in) bind(C)
  use vlidort_pars_m
  use vlidort_inputs_def_m
  use vlidort_lininputs_def_m
  use vlidort_outputs_def_m
  use vlidort_inputs_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: vlidort_fixin_in
  type(c_ptr), intent(out) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_linfixin_in
  type(c_ptr), intent(out) :: vlidort_linmodin_in
  type(c_ptr), intent(out) :: vlidort_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_input_exception_handling), pointer :: vlidort_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_inputstatus_in, vlidort_inputstatus_lcl)

  call vlidort_l_input_master(filnam_lcl, &
                              vlidort_fixin_lcl, &
                              vlidort_modin_lcl, &
                              vlidort_linfixin_lcl, &
                              vlidort_linmodin_lcl, &
                              vlidort_inputstatus_lcl)

end subroutine v_l_inputs_m_v_l_input_master_wrap

subroutine v_l_inputs_m_v_linsup_init_wrap (vlidort_linsup_in) bind(C)
  use vlidort_linsup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_linsup_in

  ! Local variables
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)

  call vlidort_linsup_init(vlidort_linsup_lcl)

end subroutine v_l_inputs_m_v_linsup_init_wrap

subroutine v_l_inputs_m_v_sleave_linsup_init_wrap (vlidort_linsup_in) bind(C)
  use vlidort_pars_m
  use vlidort_linsup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_linsup_in

  ! Local variables
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)

  call vlidort_sleave_linsup_init(vlidort_linsup_lcl)

end subroutine v_l_inputs_m_v_sleave_linsup_init_wrap

subroutine v_l_inputs_m_v_ss_linsup_init_wrap (vlidort_linsup_in) bind(C)
  use vlidort_pars_m
  use vlidort_linsup_inout_def_m

  ! Arguments
  type(c_ptr), intent(inout) :: vlidort_linsup_in

  ! Local variables
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)

  call vlidort_ss_linsup_init(vlidort_linsup_lcl)

end subroutine v_l_inputs_m_v_ss_linsup_init_wrap

end module VLIDORT_L_INPUTS_M_WRAP

module VLIDORT_LCS_MASTERS_M_WRAP

use iso_c_binding
use vlidort_lcs_masters_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine v_lcs_masters_m_v_lcs_master_wrap (do_debug_input_in, &
                                              vlidort_fixin_in, &
                                              vlidort_modin_in, &
                                              vlidort_sup_in, &
                                              vlidort_out_in, &
                                              vlidort_linfixin_in, &
                                              vlidort_linmodin_in, &
                                              vlidort_linsup_in, &
                                              vlidort_linout_in) bind(C)
  use vlidort_pars_m
  use vlidort_inputs_def_m
  use vlidort_sup_inout_def_m
  use vlidort_outputs_def_m
  use vlidort_lininputs_def_m
  use vlidort_linsup_inout_def_m
  use vlidort_linoutputs_def_m
  use vlidort_work_def_m
  use vlidort_linwork_def_m
  use vlidort_inputs_m
  use vlidort_geometry_m
  use vlidort_miscsetups_m
  use vlidort_thermalsup_m
  use vlidort_converge_m
  use vlidort_pack_m
  use vlidort_writemodules_m
  use vlidort_vfo_lcs_interface_m
  use vlidort_l_inputs_m
  use vlidort_l_pack_m
  use vlidort_l_writemodules_m
  use vlidort_l_thermalsup_m
  use vlidort_lc_pack_m
  use vlidort_lc_miscsetups_m
  use vlidort_lcs_converge_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_input_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in
  type(c_ptr), intent(in) :: vlidort_linfixin_in
  type(c_ptr), intent(inout) :: vlidort_linmodin_in
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_linout_in

  ! Local variables
  logical(kind=4) :: do_debug_input_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_linoutputs), pointer :: vlidort_linout_lcl

  ! Convert input arguments
  do_debug_input_lcl = do_debug_input_in
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_linout_in, vlidort_linout_lcl)

  call vlidort_lcs_master(do_debug_input_lcl, &
                          vlidort_fixin_lcl, &
                          vlidort_modin_lcl, &
                          vlidort_sup_lcl, &
                          vlidort_out_lcl, &
                          vlidort_linfixin_lcl, &
                          vlidort_linmodin_lcl, &
                          vlidort_linsup_lcl, &
                          vlidort_linout_lcl)

end subroutine v_lcs_masters_m_v_lcs_master_wrap

end module VLIDORT_LCS_MASTERS_M_WRAP

module VLIDORT_LPS_MASTERS_M_WRAP

use iso_c_binding
use vlidort_lps_masters_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine v_lps_masters_m_v_lps_master_wrap (do_debug_input_in, &
                                              vlidort_fixin_in, &
                                              vlidort_modin_in, &
                                              vlidort_sup_in, &
                                              vlidort_out_in, &
                                              vlidort_linfixin_in, &
                                              vlidort_linmodin_in, &
                                              vlidort_linsup_in, &
                                              vlidort_linout_in) bind(C)
  use vlidort_pars_m
  use vlidort_inputs_def_m
  use vlidort_sup_inout_def_m
  use vlidort_outputs_def_m
  use vlidort_lininputs_def_m
  use vlidort_linsup_inout_def_m
  use vlidort_linoutputs_def_m
  use vlidort_work_def_m
  use vlidort_linwork_def_m
  use vlidort_inputs_m
  use vlidort_geometry_m
  use vlidort_miscsetups_m
  use vlidort_thermalsup_m
  use vlidort_converge_m
  use vlidort_pack_m
  use vlidort_writemodules_m
  use vlidort_vfo_lps_interface_m
  use vlidort_l_inputs_m
  use vlidort_l_pack_m
  use vlidort_l_writemodules_m
  use vlidort_l_thermalsup_m
  use vlidort_lp_pack_m
  use vlidort_lp_miscsetups_m
  use vlidort_lps_converge_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_input_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in
  type(c_ptr), intent(in) :: vlidort_linfixin_in
  type(c_ptr), intent(inout) :: vlidort_linmodin_in
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_linout_in

  ! Local variables
  logical(kind=4) :: do_debug_input_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_linoutputs), pointer :: vlidort_linout_lcl

  ! Convert input arguments
  do_debug_input_lcl = do_debug_input_in
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_linout_in, vlidort_linout_lcl)

  call vlidort_lps_master(do_debug_input_lcl, &
                          vlidort_fixin_lcl, &
                          vlidort_modin_lcl, &
                          vlidort_sup_lcl, &
                          vlidort_out_lcl, &
                          vlidort_linfixin_lcl, &
                          vlidort_linmodin_lcl, &
                          vlidort_linsup_lcl, &
                          vlidort_linout_lcl)

end subroutine v_lps_masters_m_v_lps_master_wrap

end module VLIDORT_LPS_MASTERS_M_WRAP

module VLIDORT_VBRDF_SUP_ACCESSORIES_M_WRAP

use iso_c_binding
use vlidort_vbrdf_sup_accessories_m
use vlidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine v_vbrdf_sup_accessories_m_set_v_vbrdf_inputs_wrap (vbrdf_sup_out_in, &
                                                              vlidort_fixin_in, &
                                                              vlidort_modin_in, &
                                                              vlidort_sup_in) bind(C)
  use vbrdf_sup_outputs_def_m
  use vlidort_pars_m
  use vlidort_io_defs_m
  use vlidort_sup_inout_def_m

  ! Arguments
  type(c_ptr), intent(in) :: vbrdf_sup_out_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(in) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in

  ! Local variables
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)

  call set_vlidort_vbrdf_inputs(vbrdf_sup_out_lcl, &
                                vlidort_fixin_lcl, &
                                vlidort_modin_lcl, &
                                vlidort_sup_lcl)

end subroutine v_vbrdf_sup_accessories_m_set_v_vbrdf_inputs_wrap

subroutine v_vbrdf_sup_accessories_m_v_vbrdf_input_check_wrap (vbrdf_sup_in_in, &
                                                               vlidort_fixin_in, &
                                                               vlidort_modin_in, &
                                                               vlidort_vbrdfcheck_status_in) bind(C)
  use vlidort_pars_m
  use vbrdf_sup_inputs_def_m
  use vlidort_inputs_def_m
  use vlidort_outputs_def_m

  ! Arguments
  type(c_ptr), intent(in) :: vbrdf_sup_in_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(in) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_vbrdfcheck_status_in

  ! Local variables
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_exception_handling), pointer :: vlidort_vbrdfcheck_status_lcl

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_vbrdfcheck_status_in, vlidort_vbrdfcheck_status_lcl)

  call vlidort_vbrdf_input_check(vbrdf_sup_in_lcl, &
                                 vlidort_fixin_lcl, &
                                 vlidort_modin_lcl, &
                                 vlidort_vbrdfcheck_status_lcl)

end subroutine v_vbrdf_sup_accessories_m_v_vbrdf_input_check_wrap

subroutine v_vbrdf_sup_accessories_m_v_vbrdf_input_check_error_wrap (errorfile_in_len, &
                                                                     errorfile_in, &
                                                                     vlidort_vbrdfcheck_status_in) bind(C)
  use vbrdf_sup_aux_m
  use vlidort_outputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: errorfile_in_len
  character(kind=c_char) , intent(inout) :: errorfile_in(errorfile_in_len+1)
  type(c_ptr), intent(in) :: vlidort_vbrdfcheck_status_in

  ! Local variables
  character(kind=c_char, len=errorfile_in_len) :: errorfile_lcl
  integer :: len_idx
  type(vlidort_exception_handling), pointer :: vlidort_vbrdfcheck_status_lcl

  ! Convert input arguments
  do len_idx = 1, errorfile_in_len
    errorfile_lcl(len_idx:len_idx) = &
      errorfile_in(len_idx)
  end do
  call c_f_pointer(vlidort_vbrdfcheck_status_in, vlidort_vbrdfcheck_status_lcl)

  call vlidort_vbrdf_input_check_error(errorfile_lcl, &
                                       vlidort_vbrdfcheck_status_lcl)

end subroutine v_vbrdf_sup_accessories_m_v_vbrdf_input_check_error_wrap

end module VLIDORT_VBRDF_SUP_ACCESSORIES_M_WRAP
