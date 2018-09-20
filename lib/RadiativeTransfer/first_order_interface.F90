
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

module FO_SCALARSS_RTCALCS_I_M_WRAP

use iso_c_binding
use fo_scalarss_rtcalcs_i_m

! This module was auto-generated 

implicit none

contains

subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_dn_wrap (maxgeoms, &
                                                          maxlayers, &
                                                          maxfine, &
                                                          max_user_levels, &
                                                          do_deltam_scaling, &
                                                          do_planpar, &
                                                          do_regular_ps, &
                                                          do_enhanced_ps, &
                                                          donadir, &
                                                          ngeoms, &
                                                          nlayers, &
                                                          nfinedivs, &
                                                          n_user_levels, &
                                                          user_levels, &
                                                          extinction, &
                                                          deltaus, &
                                                          exactscat_dn, &
                                                          flux, &
                                                          mu1, &
                                                          ncrit, &
                                                          radcrit, &
                                                          cotcrit, &
                                                          xfine, &
                                                          wfine, &
                                                          csqfine, &
                                                          cotfine, &
                                                          raycon, &
                                                          radii, &
                                                          cota, &
                                                          sunpaths, &
                                                          ntraverse, &
                                                          sunpathsfine, &
                                                          ntraversefine, &
                                                          intensity_dn, &
                                                          cumsource_dn) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_regular_ps
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXGEOMS), intent(in) :: donadir
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: exactscat_dn
  real(c_double), intent(in) :: flux
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  integer(c_int), dimension(MAXGEOMS), intent(in) :: ncrit
  real(c_double), dimension(MAXGEOMS), intent(in) :: radcrit
  real(c_double), dimension(MAXGEOMS), intent(in) :: cotcrit
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: csqfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: cotfine
  real(c_double), dimension(MAXGEOMS), intent(in) :: raycon
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: radii
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: cota
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dn
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_dn

  ! Local variables
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_regular_ps_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXGEOMS) :: donadir_lcl

  ! Convert input arguments
  do_deltam_scaling_lcl = do_deltam_scaling
  do_planpar_lcl = do_planpar
  do_regular_ps_lcl = do_regular_ps
  do_enhanced_ps_lcl = do_enhanced_ps
  donadir_lcl = donadir

  call ss_integral_i_dn(maxgeoms, &
                        maxlayers, &
                        maxfine, &
                        max_user_levels, &
                        do_deltam_scaling_lcl, &
                        do_planpar_lcl, &
                        do_regular_ps_lcl, &
                        do_enhanced_ps_lcl, &
                        donadir_lcl, &
                        ngeoms, &
                        nlayers, &
                        nfinedivs, &
                        n_user_levels, &
                        user_levels, &
                        extinction, &
                        deltaus, &
                        exactscat_dn, &
                        flux, &
                        mu1, &
                        ncrit, &
                        radcrit, &
                        cotcrit, &
                        xfine, &
                        wfine, &
                        csqfine, &
                        cotfine, &
                        raycon, &
                        radii, &
                        cota, &
                        sunpaths, &
                        ntraverse, &
                        sunpathsfine, &
                        ntraversefine, &
                        intensity_dn, &
                        cumsource_dn)

end subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_dn_wrap

subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap (maxgeoms, &
                                                          maxlayers, &
                                                          maxfine, &
                                                          max_user_levels, &
                                                          do_deltam_scaling, &
                                                          do_planpar, &
                                                          do_regular_ps, &
                                                          do_enhanced_ps, &
                                                          donadir, &
                                                          ngeoms, &
                                                          nlayers, &
                                                          nfinedivs, &
                                                          n_user_levels, &
                                                          user_levels, &
                                                          reflec, &
                                                          extinction, &
                                                          deltaus, &
                                                          exactscat_up, &
                                                          flux, &
                                                          mu0, &
                                                          mu1, &
                                                          ncrit, &
                                                          xfine, &
                                                          wfine, &
                                                          csqfine, &
                                                          cotfine, &
                                                          raycon, &
                                                          cota, &
                                                          sunpaths, &
                                                          ntraverse, &
                                                          sunpathsfine, &
                                                          ntraversefine, &
                                                          intensity_up, &
                                                          intensity_db, &
                                                          cumsource_up) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_regular_ps
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXGEOMS), intent(in) :: donadir
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  real(c_double), dimension(MAXGEOMS), intent(in) :: reflec
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: exactscat_up
  real(c_double), intent(in) :: flux
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  integer(c_int), dimension(MAXGEOMS), intent(in) :: ncrit
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: csqfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: cotfine
  real(c_double), dimension(MAXGEOMS), intent(in) :: raycon
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: cota
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_db
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_up

  ! Local variables
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_regular_ps_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXGEOMS) :: donadir_lcl

  ! Convert input arguments
  do_deltam_scaling_lcl = do_deltam_scaling
  do_planpar_lcl = do_planpar
  do_regular_ps_lcl = do_regular_ps
  do_enhanced_ps_lcl = do_enhanced_ps
  donadir_lcl = donadir

  call ss_integral_i_up(maxgeoms, &
                        maxlayers, &
                        maxfine, &
                        max_user_levels, &
                        do_deltam_scaling_lcl, &
                        do_planpar_lcl, &
                        do_regular_ps_lcl, &
                        do_enhanced_ps_lcl, &
                        donadir_lcl, &
                        ngeoms, &
                        nlayers, &
                        nfinedivs, &
                        n_user_levels, &
                        user_levels, &
                        reflec, &
                        extinction, &
                        deltaus, &
                        exactscat_up, &
                        flux, &
                        mu0, &
                        mu1, &
                        ncrit, &
                        xfine, &
                        wfine, &
                        csqfine, &
                        cotfine, &
                        raycon, &
                        cota, &
                        sunpaths, &
                        ntraverse, &
                        sunpathsfine, &
                        ntraversefine, &
                        intensity_up, &
                        intensity_db, &
                        cumsource_up)

end subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap

subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_updn_wrap (maxgeoms, &
                                                            maxlayers, &
                                                            maxfine, &
                                                            max_user_levels, &
                                                            do_upwelling, &
                                                            do_dnwelling, &
                                                            do_deltam_scaling, &
                                                            do_planpar, &
                                                            do_regular_ps, &
                                                            do_enhanced_ps, &
                                                            donadir, &
                                                            ngeoms, &
                                                            nlayers, &
                                                            nfinedivs, &
                                                            n_user_levels, &
                                                            user_levels, &
                                                            reflec, &
                                                            extinction, &
                                                            deltaus, &
                                                            exactscat_up, &
                                                            exactscat_dn, &
                                                            flux, &
                                                            mu0, &
                                                            mu1, &
                                                            ncrit, &
                                                            radcrit, &
                                                            cotcrit, &
                                                            xfine, &
                                                            wfine, &
                                                            csqfine, &
                                                            cotfine, &
                                                            raycon, &
                                                            radii, &
                                                            cota, &
                                                            sunpaths_up, &
                                                            ntraverse_up, &
                                                            sunpathsfine_up, &
                                                            ntraversefine_up, &
                                                            sunpaths_dn, &
                                                            ntraverse_dn, &
                                                            sunpathsfine_dn, &
                                                            ntraversefine_dn, &
                                                            intensity_up, &
                                                            intensity_db, &
                                                            cumsource_up, &
                                                            intensity_dn, &
                                                            cumsource_dn) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_regular_ps
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXGEOMS), intent(in) :: donadir
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  real(c_double), dimension(MAXGEOMS), intent(in) :: reflec
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: exactscat_up
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: exactscat_dn
  real(c_double), intent(in) :: flux
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  integer(c_int), dimension(MAXGEOMS), intent(in) :: ncrit
  real(c_double), dimension(MAXGEOMS), intent(in) :: radcrit
  real(c_double), dimension(MAXGEOMS), intent(in) :: cotcrit
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: csqfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: cotfine
  real(c_double), dimension(MAXGEOMS), intent(in) :: raycon
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: radii
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: cota
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_up
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse_up
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_up
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_up
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_dn
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse_dn
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_dn
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_dn
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_db
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dn
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_dn

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_dnwelling_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_regular_ps_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXGEOMS) :: donadir_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_deltam_scaling_lcl = do_deltam_scaling
  do_planpar_lcl = do_planpar
  do_regular_ps_lcl = do_regular_ps
  do_enhanced_ps_lcl = do_enhanced_ps
  donadir_lcl = donadir

  call ss_integral_i_updn(maxgeoms, &
                          maxlayers, &
                          maxfine, &
                          max_user_levels, &
                          do_upwelling_lcl, &
                          do_dnwelling_lcl, &
                          do_deltam_scaling_lcl, &
                          do_planpar_lcl, &
                          do_regular_ps_lcl, &
                          do_enhanced_ps_lcl, &
                          donadir_lcl, &
                          ngeoms, &
                          nlayers, &
                          nfinedivs, &
                          n_user_levels, &
                          user_levels, &
                          reflec, &
                          extinction, &
                          deltaus, &
                          exactscat_up, &
                          exactscat_dn, &
                          flux, &
                          mu0, &
                          mu1, &
                          ncrit, &
                          radcrit, &
                          cotcrit, &
                          xfine, &
                          wfine, &
                          csqfine, &
                          cotfine, &
                          raycon, &
                          radii, &
                          cota, &
                          sunpaths_up, &
                          ntraverse_up, &
                          sunpathsfine_up, &
                          ntraversefine_up, &
                          sunpaths_dn, &
                          ntraverse_dn, &
                          sunpathsfine_dn, &
                          ntraversefine_dn, &
                          intensity_up, &
                          intensity_db, &
                          cumsource_up, &
                          intensity_dn, &
                          cumsource_dn)

end subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_updn_wrap

end module FO_SCALARSS_RTCALCS_I_M_WRAP
