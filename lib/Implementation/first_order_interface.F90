
module FO_DTGEOMETRY_MASTER_M_WRAP

use iso_c_binding
use fo_dtgeometry_master_m
use fo_geometry_routines_m

! This module was auto-generated 

implicit none

contains

subroutine fo_dtgeometry_master_m_fo_dtgeometry_master_wrap (maxgeoms, &
                                                             maxlayers, &
                                                             maxfine, &
                                                             do_planpar, &
                                                             do_enhanced_ps, &
                                                             ngeoms, &
                                                             nlayers, &
                                                             nfine, &
                                                             dtr, &
                                                             eradius, &
                                                             heights, &
                                                             alpha_boa, &
                                                             donadir, &
                                                             docrit, &
                                                             acrit, &
                                                             extinc, &
                                                             raycon, &
                                                             radii, &
                                                             alpha, &
                                                             cota, &
                                                             nfinedivs, &
                                                             xfine, &
                                                             wfine, &
                                                             csqfine, &
                                                             cotfine, &
                                                             alphafine, &
                                                             radiifine, &
                                                             ncrit, &
                                                             radcrit, &
                                                             cotcrit, &
                                                             mu1, &
                                                             fail, &
                                                             message_len, &
                                                             message, &
                                                             trace_len, &
                                                             trace) bind(C)
  use fo_geometry_routines_m

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxfine
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), intent(in) :: nfine
  real(c_double), intent(in) :: dtr
  real(c_double), intent(in) :: eradius
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: heights
  real(c_double), dimension(MAXGEOMS), intent(inout) :: alpha_boa
  logical(c_bool), dimension(MAXGEOMS), intent(inout) :: donadir
  logical(c_bool), intent(inout) :: docrit
  real(c_double), intent(in) :: acrit
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinc
  real(c_double), dimension(MAXGEOMS), intent(inout) :: raycon
  real(c_double), dimension(0:MAXLAYERS), intent(inout) :: radii
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: alpha
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: cota
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(inout) :: nfinedivs
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: csqfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: cotfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: alphafine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: radiifine
  integer(c_int), dimension(MAXGEOMS), intent(out) :: ncrit
  real(c_double), dimension(MAXGEOMS), intent(out) :: radcrit
  real(c_double), dimension(MAXGEOMS), intent(out) :: cotcrit
  real(c_double), dimension(MAXGEOMS), intent(out) :: mu1
  logical(c_bool), intent(out) :: fail
  integer(c_int), intent(in) :: message_len
  character(kind=c_char) , intent(inout) :: message(message_len+1)
  integer(c_int), intent(in) :: trace_len
  character(kind=c_char) , intent(inout) :: trace(trace_len+1)

  ! Local variables
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXGEOMS) :: donadir_lcl
  logical(kind=4) :: docrit_lcl
  logical(kind=4) :: fail_lcl
  character(kind=c_char, len=message_len) :: message_lcl
  integer :: len_idx
  character(kind=c_char, len=trace_len) :: trace_lcl

  ! Convert input arguments
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  donadir_lcl = donadir
  docrit_lcl = docrit

  call fo_dtgeometry_master(maxgeoms, &
                            maxlayers, &
                            maxfine, &
                            do_planpar_lcl, &
                            do_enhanced_ps_lcl, &
                            ngeoms, &
                            nlayers, &
                            nfine, &
                            dtr, &
                            eradius, &
                            heights, &
                            alpha_boa, &
                            donadir_lcl, &
                            docrit_lcl, &
                            acrit, &
                            extinc, &
                            raycon, &
                            radii, &
                            alpha, &
                            cota, &
                            nfinedivs, &
                            xfine, &
                            wfine, &
                            csqfine, &
                            cotfine, &
                            alphafine, &
                            radiifine, &
                            ncrit, &
                            radcrit, &
                            cotcrit, &
                            mu1, &
                            fail_lcl, &
                            message_lcl, &
                            trace_lcl)

  ! Convert output arguments
  donadir = donadir_lcl
  docrit = docrit_lcl
  fail = fail_lcl
  do len_idx = 1, message_len
    message(len_idx) = &
      message_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(message_lcl(:))+1
  message(len_idx) = c_null_char
  do len_idx = 1, trace_len
    trace(len_idx) = &
      trace_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(trace_lcl(:))+1
  trace(len_idx) = c_null_char

end subroutine fo_dtgeometry_master_m_fo_dtgeometry_master_wrap

end module FO_DTGEOMETRY_MASTER_M_WRAP

module FO_SSGEOMETRY_MASTER_M_WRAP

use iso_c_binding
use fo_ssgeometry_master_m
use fo_geometry_generic_m
use fo_geometry_routines_m

! This module was auto-generated 

implicit none

contains

subroutine fo_ssgeometry_master_m_fo_ssgeometry_master_wrap (maxgeoms, &
                                                             maxszas, &
                                                             maxvzas, &
                                                             maxazms, &
                                                             maxlayers, &
                                                             maxfine, &
                                                             do_obsgeom, &
                                                             do_chapman, &
                                                             do_planpar, &
                                                             do_enhanced_ps, &
                                                             ngeoms, &
                                                             nszas, &
                                                             nvzas, &
                                                             nazms, &
                                                             nlayers, &
                                                             nfine, &
                                                             dtr, &
                                                             pie, &
                                                             vsign, &
                                                             eradius, &
                                                             heights, &
                                                             obsgeom_boa, &
                                                             alpha_boa, &
                                                             theta_boa, &
                                                             phi_boa, &
                                                             donadir, &
                                                             docrit, &
                                                             acrit, &
                                                             extinc, &
                                                             raycon, &
                                                             radii, &
                                                             alpha, &
                                                             cota, &
                                                             nfinedivs, &
                                                             xfine, &
                                                             wfine, &
                                                             csqfine, &
                                                             cotfine, &
                                                             alphafine, &
                                                             radiifine, &
                                                             ncrit, &
                                                             radcrit, &
                                                             cotcrit, &
                                                             mu0, &
                                                             mu1, &
                                                             cosscat, &
                                                             chapfacs, &
                                                             sunpaths, &
                                                             ntraverse, &
                                                             sunpathsfine, &
                                                             ntraversefine, &
                                                             fail, &
                                                             message_len, &
                                                             message, &
                                                             trace_len, &
                                                             trace) bind(C)
  use fo_geometry_generic_m
  use fo_geometry_routines_m

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxszas
  integer(c_int), intent(in) :: maxvzas
  integer(c_int), intent(in) :: maxazms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxfine
  logical(c_bool), intent(in) :: do_obsgeom
  logical(c_bool), intent(in) :: do_chapman
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nszas
  integer(c_int), intent(in) :: nvzas
  integer(c_int), intent(in) :: nazms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), intent(in) :: nfine
  real(c_double), intent(in) :: dtr
  real(c_double), intent(in) :: pie
  real(c_double), intent(in) :: vsign
  real(c_double), intent(in) :: eradius
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: heights
  real(c_double), dimension(MAXGEOMS, 3), intent(inout) :: obsgeom_boa
  real(c_double), dimension(MAXVZAS), intent(inout) :: alpha_boa
  real(c_double), dimension(MAXSZAS), intent(inout) :: theta_boa
  real(c_double), dimension(MAXAZMS), intent(inout) :: phi_boa
  logical(c_bool), dimension(MAXGEOMS), intent(inout) :: donadir
  logical(c_bool), intent(inout) :: docrit
  real(c_double), intent(in) :: acrit
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinc
  real(c_double), dimension(MAXGEOMS), intent(inout) :: raycon
  real(c_double), dimension(0:MAXLAYERS), intent(inout) :: radii
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: alpha
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: cota
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(inout) :: nfinedivs
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: csqfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: cotfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: alphafine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: radiifine
  integer(c_int), dimension(MAXGEOMS), intent(out) :: ncrit
  real(c_double), dimension(MAXGEOMS), intent(out) :: radcrit
  real(c_double), dimension(MAXGEOMS), intent(out) :: cotcrit
  real(c_double), dimension(MAXGEOMS), intent(out) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(out) :: mu1
  real(c_double), dimension(MAXGEOMS), intent(out) :: cosscat
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(out) :: chapfacs
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(out) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(out) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(out) :: ntraversefine
  logical(c_bool), intent(out) :: fail
  integer(c_int), intent(in) :: message_len
  character(kind=c_char) , intent(inout) :: message(message_len+1)
  integer(c_int), intent(in) :: trace_len
  character(kind=c_char) , intent(inout) :: trace(trace_len+1)

  ! Local variables
  logical(kind=4) :: do_obsgeom_lcl
  logical(kind=4) :: do_chapman_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXGEOMS) :: donadir_lcl
  logical(kind=4) :: docrit_lcl
  logical(kind=4) :: fail_lcl
  character(kind=c_char, len=message_len) :: message_lcl
  integer :: len_idx
  character(kind=c_char, len=trace_len) :: trace_lcl

  ! Convert input arguments
  do_obsgeom_lcl = do_obsgeom
  do_chapman_lcl = do_chapman
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  donadir_lcl = donadir
  docrit_lcl = docrit

  call fo_ssgeometry_master(maxgeoms, &
                            maxszas, &
                            maxvzas, &
                            maxazms, &
                            maxlayers, &
                            maxfine, &
                            do_obsgeom_lcl, &
                            do_chapman_lcl, &
                            do_planpar_lcl, &
                            do_enhanced_ps_lcl, &
                            ngeoms, &
                            nszas, &
                            nvzas, &
                            nazms, &
                            nlayers, &
                            nfine, &
                            dtr, &
                            pie, &
                            vsign, &
                            eradius, &
                            heights, &
                            obsgeom_boa, &
                            alpha_boa, &
                            theta_boa, &
                            phi_boa, &
                            donadir_lcl, &
                            docrit_lcl, &
                            acrit, &
                            extinc, &
                            raycon, &
                            radii, &
                            alpha, &
                            cota, &
                            nfinedivs, &
                            xfine, &
                            wfine, &
                            csqfine, &
                            cotfine, &
                            alphafine, &
                            radiifine, &
                            ncrit, &
                            radcrit, &
                            cotcrit, &
                            mu0, &
                            mu1, &
                            cosscat, &
                            chapfacs, &
                            sunpaths, &
                            ntraverse, &
                            sunpathsfine, &
                            ntraversefine, &
                            fail_lcl, &
                            message_lcl, &
                            trace_lcl)

  ! Convert output arguments
  donadir = donadir_lcl
  docrit = docrit_lcl
  fail = fail_lcl
  do len_idx = 1, message_len
    message(len_idx) = &
      message_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(message_lcl(:))+1
  message(len_idx) = c_null_char
  do len_idx = 1, trace_len
    trace(len_idx) = &
      trace_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(trace_lcl(:))+1
  trace(len_idx) = c_null_char

end subroutine fo_ssgeometry_master_m_fo_ssgeometry_master_wrap

end module FO_SSGEOMETRY_MASTER_M_WRAP
