module VBRDF_LINSUP_MASTERS_M_IO

use iso_c_binding
use vlidort_interface_types_io
use vlidort_pars_m
use vbrdf_findpar_m
use vbrdf_sup_inputs_def_m
use vbrdf_sup_outputs_def_m
use vbrdf_linsup_inputs_def_m
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

! This module was auto-generated 

implicit none

contains

subroutine vbrdf_linsup_masters_m_read_wrap (filename_in, filename_in_len, vbrdf_sup_in_in, &
                                             vbrdf_linsup_in_in, &
                                             vbrdf_sup_inputstatus_in, &
                                             vbrdf_sup_out_in, &
                                             vbrdf_linsup_out_in, &
                                             vbrdf_sup_outputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_linsup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: vbrdf_sup_out_in
  type(c_ptr), intent(out) :: vbrdf_linsup_out_in
  type(c_ptr), intent(out) :: vbrdf_sup_outputstatus_in

  ! Local variables
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_linsup_inputs), pointer :: vbrdf_linsup_in_lcl
  type(vbrdf_input_exception_handling), pointer :: vbrdf_sup_inputstatus_lcl
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vbrdf_linsup_outputs), pointer :: vbrdf_linsup_out_lcl
  type(vbrdf_output_exception_handling), pointer :: vbrdf_sup_outputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_linsup_in_in, vbrdf_linsup_in_lcl)
  call c_f_pointer(vbrdf_sup_inputstatus_in, vbrdf_sup_inputstatus_lcl)
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vbrdf_linsup_out_in, vbrdf_linsup_out_lcl)
  call c_f_pointer(vbrdf_sup_outputstatus_in, vbrdf_sup_outputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vbrdf_linsup_masters_m_read(filename_lcl, vbrdf_sup_in_lcl, &
                                   vbrdf_linsup_in_lcl, &
                                   vbrdf_sup_inputstatus_lcl, &
                                   vbrdf_sup_out_lcl, &
                                   vbrdf_linsup_out_lcl, &
                                   vbrdf_sup_outputstatus_lcl)

end subroutine vbrdf_linsup_masters_m_read_wrap

subroutine vbrdf_linsup_masters_m_read (filename, vbrdf_sup_in_in, &
                                        vbrdf_linsup_in_in, &
                                        vbrdf_sup_inputstatus_in, &
                                        vbrdf_sup_out_in, &
                                        vbrdf_linsup_out_in, &
                                        vbrdf_sup_outputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vbrdf_sup_inputs), intent(inout), pointer :: vbrdf_sup_in_in
  type(vbrdf_linsup_inputs), intent(inout), pointer :: vbrdf_linsup_in_in
  type(vbrdf_input_exception_handling), intent(inout), pointer :: vbrdf_sup_inputstatus_in
  type(vbrdf_sup_outputs), intent(inout), pointer :: vbrdf_sup_out_in
  type(vbrdf_linsup_outputs), intent(inout), pointer :: vbrdf_linsup_out_in
  type(vbrdf_output_exception_handling), intent(inout), pointer :: vbrdf_sup_outputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vbrdf_sup_inputs_f_read(500, vbrdf_sup_in_in)
  call vbrdf_linsup_inputs_f_read(500, vbrdf_linsup_in_in)
  call vbrdf_input_exception_handling_f_read(500, vbrdf_sup_inputstatus_in)
  call vbrdf_sup_outputs_f_read(500, vbrdf_sup_out_in)
  call vbrdf_linsup_outputs_f_read(500, vbrdf_linsup_out_in)
  call vbrdf_output_exception_handling_f_read(500, vbrdf_sup_outputstatus_in)
  close(500)

end subroutine vbrdf_linsup_masters_m_read

subroutine vbrdf_linsup_masters_m_write_wrap (filename_in, filename_in_len, vbrdf_sup_in_in, &
                                              vbrdf_linsup_in_in, &
                                              vbrdf_sup_inputstatus_in, &
                                              vbrdf_sup_out_in, &
                                              vbrdf_linsup_out_in, &
                                              vbrdf_sup_outputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_linsup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: vbrdf_sup_out_in
  type(c_ptr), intent(out) :: vbrdf_linsup_out_in
  type(c_ptr), intent(out) :: vbrdf_sup_outputstatus_in

  ! Local variables
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_linsup_inputs), pointer :: vbrdf_linsup_in_lcl
  type(vbrdf_input_exception_handling), pointer :: vbrdf_sup_inputstatus_lcl
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vbrdf_linsup_outputs), pointer :: vbrdf_linsup_out_lcl
  type(vbrdf_output_exception_handling), pointer :: vbrdf_sup_outputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_linsup_in_in, vbrdf_linsup_in_lcl)
  call c_f_pointer(vbrdf_sup_inputstatus_in, vbrdf_sup_inputstatus_lcl)
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vbrdf_linsup_out_in, vbrdf_linsup_out_lcl)
  call c_f_pointer(vbrdf_sup_outputstatus_in, vbrdf_sup_outputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vbrdf_linsup_masters_m_write(filename_lcl, vbrdf_sup_in_lcl, &
                                    vbrdf_linsup_in_lcl, &
                                    vbrdf_sup_inputstatus_lcl, &
                                    vbrdf_sup_out_lcl, &
                                    vbrdf_linsup_out_lcl, &
                                    vbrdf_sup_outputstatus_lcl)

end subroutine vbrdf_linsup_masters_m_write_wrap

subroutine vbrdf_linsup_masters_m_write (filename, vbrdf_sup_in_in, &
                                         vbrdf_linsup_in_in, &
                                         vbrdf_sup_inputstatus_in, &
                                         vbrdf_sup_out_in, &
                                         vbrdf_linsup_out_in, &
                                         vbrdf_sup_outputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vbrdf_sup_inputs), intent(inout), pointer :: vbrdf_sup_in_in
  type(vbrdf_linsup_inputs), intent(inout), pointer :: vbrdf_linsup_in_in
  type(vbrdf_input_exception_handling), intent(inout), pointer :: vbrdf_sup_inputstatus_in
  type(vbrdf_sup_outputs), intent(inout), pointer :: vbrdf_sup_out_in
  type(vbrdf_linsup_outputs), intent(inout), pointer :: vbrdf_linsup_out_in
  type(vbrdf_output_exception_handling), intent(inout), pointer :: vbrdf_sup_outputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vbrdf_sup_inputs_f_write(500, vbrdf_sup_in_in)
  call vbrdf_linsup_inputs_f_write(500, vbrdf_linsup_in_in)
  call vbrdf_input_exception_handling_f_write(500, vbrdf_sup_inputstatus_in)
  call vbrdf_sup_outputs_f_write(500, vbrdf_sup_out_in)
  call vbrdf_linsup_outputs_f_write(500, vbrdf_linsup_out_in)
  call vbrdf_output_exception_handling_f_write(500, vbrdf_sup_outputstatus_in)
  close(500)

end subroutine vbrdf_linsup_masters_m_write


 
end module VBRDF_LINSUP_MASTERS_M_IO

module VBRDF_SUP_MASTERS_M_IO

use iso_c_binding
use vlidort_interface_types_io
use vlidort_pars_m
use vbrdf_findpar_m
use vbrdf_sup_inputs_def_m
use vbrdf_sup_outputs_def_m
use vlidort_pars_m
use vbrdf_sup_inputs_def_m
use vbrdf_sup_outputs_def_m
use vbrdf_sup_aux_m
use vbrdf_sup_kernels_m
use vbrdf_sup_routines_m

! This module was auto-generated 

implicit none

contains

subroutine vbrdf_sup_masters_m_read_wrap (filename_in, filename_in_len, vbrdf_sup_in_in, &
                                          vbrdf_sup_inputstatus_in, &
                                          vbrdf_sup_out_in, &
                                          vbrdf_sup_outputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: vbrdf_sup_out_in
  type(c_ptr), intent(out) :: vbrdf_sup_outputstatus_in

  ! Local variables
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_input_exception_handling), pointer :: vbrdf_sup_inputstatus_lcl
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vbrdf_output_exception_handling), pointer :: vbrdf_sup_outputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_sup_inputstatus_in, vbrdf_sup_inputstatus_lcl)
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vbrdf_sup_outputstatus_in, vbrdf_sup_outputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vbrdf_sup_masters_m_read(filename_lcl, vbrdf_sup_in_lcl, &
                                vbrdf_sup_inputstatus_lcl, &
                                vbrdf_sup_out_lcl, &
                                vbrdf_sup_outputstatus_lcl)

end subroutine vbrdf_sup_masters_m_read_wrap

subroutine vbrdf_sup_masters_m_read (filename, vbrdf_sup_in_in, &
                                     vbrdf_sup_inputstatus_in, &
                                     vbrdf_sup_out_in, &
                                     vbrdf_sup_outputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vbrdf_sup_inputs), intent(inout), pointer :: vbrdf_sup_in_in
  type(vbrdf_input_exception_handling), intent(inout), pointer :: vbrdf_sup_inputstatus_in
  type(vbrdf_sup_outputs), intent(inout), pointer :: vbrdf_sup_out_in
  type(vbrdf_output_exception_handling), intent(inout), pointer :: vbrdf_sup_outputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vbrdf_sup_inputs_f_read(500, vbrdf_sup_in_in)
  call vbrdf_input_exception_handling_f_read(500, vbrdf_sup_inputstatus_in)
  call vbrdf_sup_outputs_f_read(500, vbrdf_sup_out_in)
  call vbrdf_output_exception_handling_f_read(500, vbrdf_sup_outputstatus_in)
  close(500)

end subroutine vbrdf_sup_masters_m_read

subroutine vbrdf_sup_masters_m_write_wrap (filename_in, filename_in_len, vbrdf_sup_in_in, &
                                           vbrdf_sup_inputstatus_in, &
                                           vbrdf_sup_out_in, &
                                           vbrdf_sup_outputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vbrdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: vbrdf_sup_out_in
  type(c_ptr), intent(out) :: vbrdf_sup_outputstatus_in

  ! Local variables
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vbrdf_input_exception_handling), pointer :: vbrdf_sup_inputstatus_lcl
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vbrdf_output_exception_handling), pointer :: vbrdf_sup_outputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vbrdf_sup_inputstatus_in, vbrdf_sup_inputstatus_lcl)
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vbrdf_sup_outputstatus_in, vbrdf_sup_outputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vbrdf_sup_masters_m_write(filename_lcl, vbrdf_sup_in_lcl, &
                                 vbrdf_sup_inputstatus_lcl, &
                                 vbrdf_sup_out_lcl, &
                                 vbrdf_sup_outputstatus_lcl)

end subroutine vbrdf_sup_masters_m_write_wrap

subroutine vbrdf_sup_masters_m_write (filename, vbrdf_sup_in_in, &
                                      vbrdf_sup_inputstatus_in, &
                                      vbrdf_sup_out_in, &
                                      vbrdf_sup_outputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vbrdf_sup_inputs), intent(inout), pointer :: vbrdf_sup_in_in
  type(vbrdf_input_exception_handling), intent(inout), pointer :: vbrdf_sup_inputstatus_in
  type(vbrdf_sup_outputs), intent(inout), pointer :: vbrdf_sup_out_in
  type(vbrdf_output_exception_handling), intent(inout), pointer :: vbrdf_sup_outputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vbrdf_sup_inputs_f_write(500, vbrdf_sup_in_in)
  call vbrdf_input_exception_handling_f_write(500, vbrdf_sup_inputstatus_in)
  call vbrdf_sup_outputs_f_write(500, vbrdf_sup_out_in)
  call vbrdf_output_exception_handling_f_write(500, vbrdf_sup_outputstatus_in)
  close(500)

end subroutine vbrdf_sup_masters_m_write


 
end module VBRDF_SUP_MASTERS_M_IO

module VLIDORT_INPUTS_M_IO

use iso_c_binding
use vlidort_interface_types_io
use vlidort_pars_m
use vlidort_sup_inout_def_m
use vlidort_pars_m
use vlidort_inputs_def_m
use vlidort_outputs_def_m
use vlidort_pars_m
use vlidort_sup_inout_def_m
use vlidort_pars_m
use vlidort_sup_inout_def_m
use vlidort_sup_inout_def_m

! This module was auto-generated 

implicit none

contains

subroutine v_inputs_m_read_wrap (filename_in, filename_in_len, vlidort_sup_in, &
                                 vlidort_fixin_in, &
                                 vlidort_modin_in, &
                                 vlidort_inputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_fixin_in
  type(c_ptr), intent(out) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_inputstatus_in

  ! Local variables
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_input_exception_handling), pointer :: vlidort_inputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_inputstatus_in, vlidort_inputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_inputs_m_read(filename_lcl, vlidort_sup_lcl, &
                             vlidort_fixin_lcl, &
                             vlidort_modin_lcl, &
                             vlidort_inputstatus_lcl)

end subroutine v_inputs_m_read_wrap

subroutine vlidort_inputs_m_read (filename, vlidort_sup_in, &
                                  vlidort_fixin_in, &
                                  vlidort_modin_in, &
                                  vlidort_inputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_input_exception_handling), intent(inout), pointer :: vlidort_inputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_sup_inout_f_read(500, vlidort_sup_in)
  call vlidort_fixed_inputs_f_read(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_read(500, vlidort_modin_in)
  call vlidort_input_exception_handling_f_read(500, vlidort_inputstatus_in)
  close(500)

end subroutine vlidort_inputs_m_read

subroutine v_inputs_m_write_wrap (filename_in, filename_in_len, vlidort_sup_in, &
                                  vlidort_fixin_in, &
                                  vlidort_modin_in, &
                                  vlidort_inputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_fixin_in
  type(c_ptr), intent(out) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_inputstatus_in

  ! Local variables
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_input_exception_handling), pointer :: vlidort_inputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_inputstatus_in, vlidort_inputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_inputs_m_write(filename_lcl, vlidort_sup_lcl, &
                              vlidort_fixin_lcl, &
                              vlidort_modin_lcl, &
                              vlidort_inputstatus_lcl)

end subroutine v_inputs_m_write_wrap

subroutine vlidort_inputs_m_write (filename, vlidort_sup_in, &
                                   vlidort_fixin_in, &
                                   vlidort_modin_in, &
                                   vlidort_inputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_input_exception_handling), intent(inout), pointer :: vlidort_inputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_sup_inout_f_write(500, vlidort_sup_in)
  call vlidort_fixed_inputs_f_write(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_write(500, vlidort_modin_in)
  call vlidort_input_exception_handling_f_write(500, vlidort_inputstatus_in)
  close(500)

end subroutine vlidort_inputs_m_write


 
end module VLIDORT_INPUTS_M_IO

module VLIDORT_MASTERS_M_IO

use iso_c_binding
use vlidort_interface_types_io
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

! This module was auto-generated 

implicit none

contains

subroutine v_masters_m_read_wrap (filename_in, filename_in_len, vlidort_fixin_in, &
                                  vlidort_modin_in, &
                                  vlidort_sup_in, &
                                  vlidort_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in

  ! Local variables
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_masters_m_read(filename_lcl, vlidort_fixin_lcl, &
                              vlidort_modin_lcl, &
                              vlidort_sup_lcl, &
                              vlidort_out_lcl)

end subroutine v_masters_m_read_wrap

subroutine vlidort_masters_m_read (filename, vlidort_fixin_in, &
                                   vlidort_modin_in, &
                                   vlidort_sup_in, &
                                   vlidort_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_outputs), intent(inout), pointer :: vlidort_out_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_fixed_inputs_f_read(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_read(500, vlidort_modin_in)
  call vlidort_sup_inout_f_read(500, vlidort_sup_in)
  call vlidort_outputs_f_read(500, vlidort_out_in)
  close(500)

end subroutine vlidort_masters_m_read

subroutine v_masters_m_write_wrap (filename_in, filename_in_len, vlidort_fixin_in, &
                                   vlidort_modin_in, &
                                   vlidort_sup_in, &
                                   vlidort_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in

  ! Local variables
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_masters_m_write(filename_lcl, vlidort_fixin_lcl, &
                               vlidort_modin_lcl, &
                               vlidort_sup_lcl, &
                               vlidort_out_lcl)

end subroutine v_masters_m_write_wrap

subroutine vlidort_masters_m_write (filename, vlidort_fixin_in, &
                                    vlidort_modin_in, &
                                    vlidort_sup_in, &
                                    vlidort_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_outputs), intent(inout), pointer :: vlidort_out_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_fixed_inputs_f_write(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_write(500, vlidort_modin_in)
  call vlidort_sup_inout_f_write(500, vlidort_sup_in)
  call vlidort_outputs_f_write(500, vlidort_out_in)
  close(500)

end subroutine vlidort_masters_m_write


 
end module VLIDORT_MASTERS_M_IO

module VLIDORT_L_INPUTS_M_IO

use iso_c_binding
use vlidort_interface_types_io
use vlidort_pars_m
use vlidort_linsup_inout_def_m
use vlidort_pars_m
use vlidort_inputs_def_m
use vlidort_lininputs_def_m
use vlidort_outputs_def_m
use vlidort_inputs_m
use vlidort_linsup_inout_def_m
use vlidort_pars_m
use vlidort_linsup_inout_def_m
use vlidort_pars_m
use vlidort_linsup_inout_def_m

! This module was auto-generated 

implicit none

contains

subroutine v_l_inputs_m_read_wrap (filename_in, filename_in_len, vlidort_linsup_in, &
                                   vlidort_fixin_in, &
                                   vlidort_modin_in, &
                                   vlidort_linfixin_in, &
                                   vlidort_linmodin_in, &
                                   vlidort_inputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_fixin_in
  type(c_ptr), intent(out) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_linfixin_in
  type(c_ptr), intent(out) :: vlidort_linmodin_in
  type(c_ptr), intent(out) :: vlidort_inputstatus_in

  ! Local variables
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_input_exception_handling), pointer :: vlidort_inputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_inputstatus_in, vlidort_inputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_l_inputs_m_read(filename_lcl, vlidort_linsup_lcl, &
                               vlidort_fixin_lcl, &
                               vlidort_modin_lcl, &
                               vlidort_linfixin_lcl, &
                               vlidort_linmodin_lcl, &
                               vlidort_inputstatus_lcl)

end subroutine v_l_inputs_m_read_wrap

subroutine vlidort_l_inputs_m_read (filename, vlidort_linsup_in, &
                                    vlidort_fixin_in, &
                                    vlidort_modin_in, &
                                    vlidort_linfixin_in, &
                                    vlidort_linmodin_in, &
                                    vlidort_inputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_linsup_inout), intent(inout), pointer :: vlidort_linsup_in
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_fixed_lininputs), intent(inout), pointer :: vlidort_linfixin_in
  type(vlidort_modified_lininputs), intent(inout), pointer :: vlidort_linmodin_in
  type(vlidort_input_exception_handling), intent(inout), pointer :: vlidort_inputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_linsup_inout_f_read(500, vlidort_linsup_in)
  call vlidort_fixed_inputs_f_read(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_read(500, vlidort_modin_in)
  call vlidort_fixed_lininputs_f_read(500, vlidort_linfixin_in)
  call vlidort_modified_lininputs_f_read(500, vlidort_linmodin_in)
  call vlidort_input_exception_handling_f_read(500, vlidort_inputstatus_in)
  close(500)

end subroutine vlidort_l_inputs_m_read

subroutine v_l_inputs_m_write_wrap (filename_in, filename_in_len, vlidort_linsup_in, &
                                    vlidort_fixin_in, &
                                    vlidort_modin_in, &
                                    vlidort_linfixin_in, &
                                    vlidort_linmodin_in, &
                                    vlidort_inputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_fixin_in
  type(c_ptr), intent(out) :: vlidort_modin_in
  type(c_ptr), intent(out) :: vlidort_linfixin_in
  type(c_ptr), intent(out) :: vlidort_linmodin_in
  type(c_ptr), intent(out) :: vlidort_inputstatus_in

  ! Local variables
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_input_exception_handling), pointer :: vlidort_inputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_inputstatus_in, vlidort_inputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_l_inputs_m_write(filename_lcl, vlidort_linsup_lcl, &
                                vlidort_fixin_lcl, &
                                vlidort_modin_lcl, &
                                vlidort_linfixin_lcl, &
                                vlidort_linmodin_lcl, &
                                vlidort_inputstatus_lcl)

end subroutine v_l_inputs_m_write_wrap

subroutine vlidort_l_inputs_m_write (filename, vlidort_linsup_in, &
                                     vlidort_fixin_in, &
                                     vlidort_modin_in, &
                                     vlidort_linfixin_in, &
                                     vlidort_linmodin_in, &
                                     vlidort_inputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_linsup_inout), intent(inout), pointer :: vlidort_linsup_in
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_fixed_lininputs), intent(inout), pointer :: vlidort_linfixin_in
  type(vlidort_modified_lininputs), intent(inout), pointer :: vlidort_linmodin_in
  type(vlidort_input_exception_handling), intent(inout), pointer :: vlidort_inputstatus_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_linsup_inout_f_write(500, vlidort_linsup_in)
  call vlidort_fixed_inputs_f_write(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_write(500, vlidort_modin_in)
  call vlidort_fixed_lininputs_f_write(500, vlidort_linfixin_in)
  call vlidort_modified_lininputs_f_write(500, vlidort_linmodin_in)
  call vlidort_input_exception_handling_f_write(500, vlidort_inputstatus_in)
  close(500)

end subroutine vlidort_l_inputs_m_write


 
end module VLIDORT_L_INPUTS_M_IO

module VLIDORT_LCS_MASTERS_M_IO

use iso_c_binding
use vlidort_interface_types_io
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

! This module was auto-generated 

implicit none

contains

subroutine v_lcs_masters_m_read_wrap (filename_in, filename_in_len, vlidort_fixin_in, &
                                      vlidort_modin_in, &
                                      vlidort_sup_in, &
                                      vlidort_out_in, &
                                      vlidort_linfixin_in, &
                                      vlidort_linmodin_in, &
                                      vlidort_linsup_in, &
                                      vlidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in
  type(c_ptr), intent(in) :: vlidort_linfixin_in
  type(c_ptr), intent(inout) :: vlidort_linmodin_in
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_linout_in

  ! Local variables
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_linoutputs), pointer :: vlidort_linout_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_linout_in, vlidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_lcs_masters_m_read(filename_lcl, vlidort_fixin_lcl, &
                                  vlidort_modin_lcl, &
                                  vlidort_sup_lcl, &
                                  vlidort_out_lcl, &
                                  vlidort_linfixin_lcl, &
                                  vlidort_linmodin_lcl, &
                                  vlidort_linsup_lcl, &
                                  vlidort_linout_lcl)

end subroutine v_lcs_masters_m_read_wrap

subroutine vlidort_lcs_masters_m_read (filename, vlidort_fixin_in, &
                                       vlidort_modin_in, &
                                       vlidort_sup_in, &
                                       vlidort_out_in, &
                                       vlidort_linfixin_in, &
                                       vlidort_linmodin_in, &
                                       vlidort_linsup_in, &
                                       vlidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_outputs), intent(inout), pointer :: vlidort_out_in
  type(vlidort_fixed_lininputs), intent(inout), pointer :: vlidort_linfixin_in
  type(vlidort_modified_lininputs), intent(inout), pointer :: vlidort_linmodin_in
  type(vlidort_linsup_inout), intent(inout), pointer :: vlidort_linsup_in
  type(vlidort_linoutputs), intent(inout), pointer :: vlidort_linout_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_fixed_inputs_f_read(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_read(500, vlidort_modin_in)
  call vlidort_sup_inout_f_read(500, vlidort_sup_in)
  call vlidort_outputs_f_read(500, vlidort_out_in)
  call vlidort_fixed_lininputs_f_read(500, vlidort_linfixin_in)
  call vlidort_modified_lininputs_f_read(500, vlidort_linmodin_in)
  call vlidort_linsup_inout_f_read(500, vlidort_linsup_in)
  call vlidort_linoutputs_f_read(500, vlidort_linout_in)
  close(500)

end subroutine vlidort_lcs_masters_m_read

subroutine v_lcs_masters_m_write_wrap (filename_in, filename_in_len, vlidort_fixin_in, &
                                       vlidort_modin_in, &
                                       vlidort_sup_in, &
                                       vlidort_out_in, &
                                       vlidort_linfixin_in, &
                                       vlidort_linmodin_in, &
                                       vlidort_linsup_in, &
                                       vlidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in
  type(c_ptr), intent(in) :: vlidort_linfixin_in
  type(c_ptr), intent(inout) :: vlidort_linmodin_in
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_linout_in

  ! Local variables
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_linoutputs), pointer :: vlidort_linout_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_linout_in, vlidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_lcs_masters_m_write(filename_lcl, vlidort_fixin_lcl, &
                                   vlidort_modin_lcl, &
                                   vlidort_sup_lcl, &
                                   vlidort_out_lcl, &
                                   vlidort_linfixin_lcl, &
                                   vlidort_linmodin_lcl, &
                                   vlidort_linsup_lcl, &
                                   vlidort_linout_lcl)

end subroutine v_lcs_masters_m_write_wrap

subroutine vlidort_lcs_masters_m_write (filename, vlidort_fixin_in, &
                                        vlidort_modin_in, &
                                        vlidort_sup_in, &
                                        vlidort_out_in, &
                                        vlidort_linfixin_in, &
                                        vlidort_linmodin_in, &
                                        vlidort_linsup_in, &
                                        vlidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_outputs), intent(inout), pointer :: vlidort_out_in
  type(vlidort_fixed_lininputs), intent(inout), pointer :: vlidort_linfixin_in
  type(vlidort_modified_lininputs), intent(inout), pointer :: vlidort_linmodin_in
  type(vlidort_linsup_inout), intent(inout), pointer :: vlidort_linsup_in
  type(vlidort_linoutputs), intent(inout), pointer :: vlidort_linout_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_fixed_inputs_f_write(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_write(500, vlidort_modin_in)
  call vlidort_sup_inout_f_write(500, vlidort_sup_in)
  call vlidort_outputs_f_write(500, vlidort_out_in)
  call vlidort_fixed_lininputs_f_write(500, vlidort_linfixin_in)
  call vlidort_modified_lininputs_f_write(500, vlidort_linmodin_in)
  call vlidort_linsup_inout_f_write(500, vlidort_linsup_in)
  call vlidort_linoutputs_f_write(500, vlidort_linout_in)
  close(500)

end subroutine vlidort_lcs_masters_m_write


 
end module VLIDORT_LCS_MASTERS_M_IO

module VLIDORT_LPS_MASTERS_M_IO

use iso_c_binding
use vlidort_interface_types_io
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

! This module was auto-generated 

implicit none

contains

subroutine v_lps_masters_m_read_wrap (filename_in, filename_in_len, vlidort_fixin_in, &
                                      vlidort_modin_in, &
                                      vlidort_sup_in, &
                                      vlidort_out_in, &
                                      vlidort_linfixin_in, &
                                      vlidort_linmodin_in, &
                                      vlidort_linsup_in, &
                                      vlidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in
  type(c_ptr), intent(in) :: vlidort_linfixin_in
  type(c_ptr), intent(inout) :: vlidort_linmodin_in
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_linout_in

  ! Local variables
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_linoutputs), pointer :: vlidort_linout_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_linout_in, vlidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_lps_masters_m_read(filename_lcl, vlidort_fixin_lcl, &
                                  vlidort_modin_lcl, &
                                  vlidort_sup_lcl, &
                                  vlidort_out_lcl, &
                                  vlidort_linfixin_lcl, &
                                  vlidort_linmodin_lcl, &
                                  vlidort_linsup_lcl, &
                                  vlidort_linout_lcl)

end subroutine v_lps_masters_m_read_wrap

subroutine vlidort_lps_masters_m_read (filename, vlidort_fixin_in, &
                                       vlidort_modin_in, &
                                       vlidort_sup_in, &
                                       vlidort_out_in, &
                                       vlidort_linfixin_in, &
                                       vlidort_linmodin_in, &
                                       vlidort_linsup_in, &
                                       vlidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_outputs), intent(inout), pointer :: vlidort_out_in
  type(vlidort_fixed_lininputs), intent(inout), pointer :: vlidort_linfixin_in
  type(vlidort_modified_lininputs), intent(inout), pointer :: vlidort_linmodin_in
  type(vlidort_linsup_inout), intent(inout), pointer :: vlidort_linsup_in
  type(vlidort_linoutputs), intent(inout), pointer :: vlidort_linout_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_fixed_inputs_f_read(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_read(500, vlidort_modin_in)
  call vlidort_sup_inout_f_read(500, vlidort_sup_in)
  call vlidort_outputs_f_read(500, vlidort_out_in)
  call vlidort_fixed_lininputs_f_read(500, vlidort_linfixin_in)
  call vlidort_modified_lininputs_f_read(500, vlidort_linmodin_in)
  call vlidort_linsup_inout_f_read(500, vlidort_linsup_in)
  call vlidort_linoutputs_f_read(500, vlidort_linout_in)
  close(500)

end subroutine vlidort_lps_masters_m_read

subroutine v_lps_masters_m_write_wrap (filename_in, filename_in_len, vlidort_fixin_in, &
                                       vlidort_modin_in, &
                                       vlidort_sup_in, &
                                       vlidort_out_in, &
                                       vlidort_linfixin_in, &
                                       vlidort_linmodin_in, &
                                       vlidort_linsup_in, &
                                       vlidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(inout) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(out) :: vlidort_out_in
  type(c_ptr), intent(in) :: vlidort_linfixin_in
  type(c_ptr), intent(inout) :: vlidort_linmodin_in
  type(c_ptr), intent(inout) :: vlidort_linsup_in
  type(c_ptr), intent(out) :: vlidort_linout_in

  ! Local variables
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vlidort_outputs), pointer :: vlidort_out_lcl
  type(vlidort_fixed_lininputs), pointer :: vlidort_linfixin_lcl
  type(vlidort_modified_lininputs), pointer :: vlidort_linmodin_lcl
  type(vlidort_linsup_inout), pointer :: vlidort_linsup_lcl
  type(vlidort_linoutputs), pointer :: vlidort_linout_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vlidort_out_in, vlidort_out_lcl)
  call c_f_pointer(vlidort_linfixin_in, vlidort_linfixin_lcl)
  call c_f_pointer(vlidort_linmodin_in, vlidort_linmodin_lcl)
  call c_f_pointer(vlidort_linsup_in, vlidort_linsup_lcl)
  call c_f_pointer(vlidort_linout_in, vlidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_lps_masters_m_write(filename_lcl, vlidort_fixin_lcl, &
                                   vlidort_modin_lcl, &
                                   vlidort_sup_lcl, &
                                   vlidort_out_lcl, &
                                   vlidort_linfixin_lcl, &
                                   vlidort_linmodin_lcl, &
                                   vlidort_linsup_lcl, &
                                   vlidort_linout_lcl)

end subroutine v_lps_masters_m_write_wrap

subroutine vlidort_lps_masters_m_write (filename, vlidort_fixin_in, &
                                        vlidort_modin_in, &
                                        vlidort_sup_in, &
                                        vlidort_out_in, &
                                        vlidort_linfixin_in, &
                                        vlidort_linmodin_in, &
                                        vlidort_linsup_in, &
                                        vlidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vlidort_outputs), intent(inout), pointer :: vlidort_out_in
  type(vlidort_fixed_lininputs), intent(inout), pointer :: vlidort_linfixin_in
  type(vlidort_modified_lininputs), intent(inout), pointer :: vlidort_linmodin_in
  type(vlidort_linsup_inout), intent(inout), pointer :: vlidort_linsup_in
  type(vlidort_linoutputs), intent(inout), pointer :: vlidort_linout_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vlidort_fixed_inputs_f_write(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_write(500, vlidort_modin_in)
  call vlidort_sup_inout_f_write(500, vlidort_sup_in)
  call vlidort_outputs_f_write(500, vlidort_out_in)
  call vlidort_fixed_lininputs_f_write(500, vlidort_linfixin_in)
  call vlidort_modified_lininputs_f_write(500, vlidort_linmodin_in)
  call vlidort_linsup_inout_f_write(500, vlidort_linsup_in)
  call vlidort_linoutputs_f_write(500, vlidort_linout_in)
  close(500)

end subroutine vlidort_lps_masters_m_write


 
end module VLIDORT_LPS_MASTERS_M_IO

module VLIDORT_VBRDF_SUP_ACCESSORIES_M_IO

use iso_c_binding
use vlidort_interface_types_io
use vbrdf_sup_outputs_def_m
use vlidort_pars_m
use vlidort_io_defs_m
use vlidort_sup_inout_def_m
use vlidort_pars_m
use vbrdf_sup_inputs_def_m
use vlidort_inputs_def_m
use vlidort_outputs_def_m
use vbrdf_sup_aux_m
use vlidort_outputs_def_m

! This module was auto-generated 

implicit none

contains

subroutine v_vbrdf_sup_accessories_m_read_wrap (filename_in, filename_in_len, vbrdf_sup_out_in, &
                                                vlidort_fixin_in, &
                                                vlidort_modin_in, &
                                                vlidort_sup_in, &
                                                vbrdf_sup_in_in, &
                                                vlidort_vbrdfcheck_status_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vbrdf_sup_out_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(in) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(in) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vlidort_vbrdfcheck_status_in

  ! Local variables
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vlidort_exception_handling), pointer :: vlidort_vbrdfcheck_status_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vlidort_vbrdfcheck_status_in, vlidort_vbrdfcheck_status_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_vbrdf_sup_accessories_m_read(filename_lcl, vbrdf_sup_out_lcl, &
                                            vlidort_fixin_lcl, &
                                            vlidort_modin_lcl, &
                                            vlidort_sup_lcl, &
                                            vbrdf_sup_in_lcl, &
                                            vlidort_vbrdfcheck_status_lcl)

end subroutine v_vbrdf_sup_accessories_m_read_wrap

subroutine vlidort_vbrdf_sup_accessories_m_read (filename, vbrdf_sup_out_in, &
                                                 vlidort_fixin_in, &
                                                 vlidort_modin_in, &
                                                 vlidort_sup_in, &
                                                 vbrdf_sup_in_in, &
                                                 vlidort_vbrdfcheck_status_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vbrdf_sup_outputs), intent(inout), pointer :: vbrdf_sup_out_in
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vbrdf_sup_inputs), intent(inout), pointer :: vbrdf_sup_in_in
  type(vlidort_exception_handling), intent(inout), pointer :: vlidort_vbrdfcheck_status_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vbrdf_sup_outputs_f_read(500, vbrdf_sup_out_in)
  call vlidort_fixed_inputs_f_read(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_read(500, vlidort_modin_in)
  call vlidort_sup_inout_f_read(500, vlidort_sup_in)
  call vbrdf_sup_inputs_f_read(500, vbrdf_sup_in_in)
  call vlidort_exception_handling_f_read(500, vlidort_vbrdfcheck_status_in)
  close(500)

end subroutine vlidort_vbrdf_sup_accessories_m_read

subroutine v_vbrdf_sup_accessories_m_write_wrap (filename_in, filename_in_len, vbrdf_sup_out_in, &
                                                 vlidort_fixin_in, &
                                                 vlidort_modin_in, &
                                                 vlidort_sup_in, &
                                                 vbrdf_sup_in_in, &
                                                 vlidort_vbrdfcheck_status_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: vbrdf_sup_out_in
  type(c_ptr), intent(in) :: vlidort_fixin_in
  type(c_ptr), intent(in) :: vlidort_modin_in
  type(c_ptr), intent(inout) :: vlidort_sup_in
  type(c_ptr), intent(in) :: vbrdf_sup_in_in
  type(c_ptr), intent(out) :: vlidort_vbrdfcheck_status_in

  ! Local variables
  type(vbrdf_sup_outputs), pointer :: vbrdf_sup_out_lcl
  type(vlidort_fixed_inputs), pointer :: vlidort_fixin_lcl
  type(vlidort_modified_inputs), pointer :: vlidort_modin_lcl
  type(vlidort_sup_inout), pointer :: vlidort_sup_lcl
  type(vbrdf_sup_inputs), pointer :: vbrdf_sup_in_lcl
  type(vlidort_exception_handling), pointer :: vlidort_vbrdfcheck_status_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(vbrdf_sup_out_in, vbrdf_sup_out_lcl)
  call c_f_pointer(vlidort_fixin_in, vlidort_fixin_lcl)
  call c_f_pointer(vlidort_modin_in, vlidort_modin_lcl)
  call c_f_pointer(vlidort_sup_in, vlidort_sup_lcl)
  call c_f_pointer(vbrdf_sup_in_in, vbrdf_sup_in_lcl)
  call c_f_pointer(vlidort_vbrdfcheck_status_in, vlidort_vbrdfcheck_status_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call vlidort_vbrdf_sup_accessories_m_write(filename_lcl, vbrdf_sup_out_lcl, &
                                             vlidort_fixin_lcl, &
                                             vlidort_modin_lcl, &
                                             vlidort_sup_lcl, &
                                             vbrdf_sup_in_lcl, &
                                             vlidort_vbrdfcheck_status_lcl)

end subroutine v_vbrdf_sup_accessories_m_write_wrap

subroutine vlidort_vbrdf_sup_accessories_m_write (filename, vbrdf_sup_out_in, &
                                                  vlidort_fixin_in, &
                                                  vlidort_modin_in, &
                                                  vlidort_sup_in, &
                                                  vbrdf_sup_in_in, &
                                                  vlidort_vbrdfcheck_status_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(vbrdf_sup_outputs), intent(inout), pointer :: vbrdf_sup_out_in
  type(vlidort_fixed_inputs), intent(inout), pointer :: vlidort_fixin_in
  type(vlidort_modified_inputs), intent(inout), pointer :: vlidort_modin_in
  type(vlidort_sup_inout), intent(inout), pointer :: vlidort_sup_in
  type(vbrdf_sup_inputs), intent(inout), pointer :: vbrdf_sup_in_in
  type(vlidort_exception_handling), intent(inout), pointer :: vlidort_vbrdfcheck_status_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call vbrdf_sup_outputs_f_write(500, vbrdf_sup_out_in)
  call vlidort_fixed_inputs_f_write(500, vlidort_fixin_in)
  call vlidort_modified_inputs_f_write(500, vlidort_modin_in)
  call vlidort_sup_inout_f_write(500, vlidort_sup_in)
  call vbrdf_sup_inputs_f_write(500, vbrdf_sup_in_in)
  call vlidort_exception_handling_f_write(500, vlidort_vbrdfcheck_status_in)
  close(500)

end subroutine vlidort_vbrdf_sup_accessories_m_write


 
end module VLIDORT_VBRDF_SUP_ACCESSORIES_M_IO

