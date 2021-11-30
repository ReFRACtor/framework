
module FO_DTWPGEOMETRY_MASTER_M_WRAP

use iso_c_binding
use fo_dtwpgeometry_master_m
use fo_wpgeometry_routines_m

! This module was auto-generated 

implicit none

contains

subroutine fo_dtwpgeometry_master_m_fo_dtwpgeometry_master_wrap (maxgeoms, &
                                                                 maxlayers, &
                                                                 maxpartials, &
                                                                 maxfine, &
                                                                 dtr, &
                                                                 eradius, &
                                                                 do_upwelling, &
                                                                 do_planpar, &
                                                                 do_enhanced_ps, &
                                                                 do_partials, &
                                                                 ngeoms, &
                                                                 nlayers, &
                                                                 npartials, &
                                                                 nfine, &
                                                                 partial_layeridx, &
                                                                 heights, &
                                                                 alpha_boa, &
                                                                 partial_heights, &
                                                                 mu1, &
                                                                 radii, &
                                                                 losw_paths, &
                                                                 alpha, &
                                                                 sina, &
                                                                 cosa, &
                                                                 radii_p, &
                                                                 losp_paths, &
                                                                 alpha_p, &
                                                                 sina_p, &
                                                                 cosa_p, &
                                                                 nfinedivs, &
                                                                 xfine, &
                                                                 wfine, &
                                                                 radiifine, &
                                                                 alphafine, &
                                                                 sinfine, &
                                                                 cosfine, &
                                                                 nfinedivs_p, &
                                                                 xfine_p, &
                                                                 wfine_p, &
                                                                 radiifine_p, &
                                                                 alphafine_p, &
                                                                 sinfine_p, &
                                                                 cosfine_p, &
                                                                 fail, &
                                                                 message_len, &
                                                                 message, &
                                                                 trace_len, &
                                                                 trace) bind(C)
  use fo_wpgeometry_routines_m

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  real(c_double), intent(in) :: dtr
  real(c_double), intent(in) :: eradius
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), intent(in) :: do_partials
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), intent(in) :: npartials
  integer(c_int), intent(in) :: nfine
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: heights
  real(c_double), dimension(MAXGEOMS), intent(inout) :: alpha_boa
  real(c_double), dimension(MAXPARTIALS), intent(in) :: partial_heights
  real(c_double), dimension(MAXGEOMS), intent(out) :: mu1
  real(c_double), dimension(0:MAXLAYERS), intent(inout) :: radii
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(out) :: losw_paths
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: alpha
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: sina
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: cosa
  real(c_double), dimension(MAXPARTIALS), intent(inout) :: radii_p
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: losp_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: alpha_p
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: sina_p
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: cosa_p
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(inout) :: nfinedivs
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: radiifine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: alphafine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: sinfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: cosfine
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: nfinedivs_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: radiifine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: alphafine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: sinfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: cosfine_p
  logical(c_bool), intent(out) :: fail
  integer(c_int), intent(in) :: message_len
  character(kind=c_char) , intent(inout) :: message(message_len+1)
  integer(c_int), intent(in) :: trace_len
  character(kind=c_char) , intent(inout) :: trace(trace_len+1)

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: fail_lcl
  character(kind=c_char, len=message_len) :: message_lcl
  integer :: len_idx
  character(kind=c_char, len=trace_len) :: trace_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_partials_lcl = do_partials

  call fo_dtwpgeometry_master(maxgeoms, &
                              maxlayers, &
                              maxpartials, &
                              maxfine, &
                              dtr, &
                              eradius, &
                              do_upwelling_lcl, &
                              do_planpar_lcl, &
                              do_enhanced_ps_lcl, &
                              do_partials_lcl, &
                              ngeoms, &
                              nlayers, &
                              npartials, &
                              nfine, &
                              partial_layeridx, &
                              heights, &
                              alpha_boa, &
                              partial_heights, &
                              mu1, &
                              radii, &
                              losw_paths, &
                              alpha, &
                              sina, &
                              cosa, &
                              radii_p, &
                              losp_paths, &
                              alpha_p, &
                              sina_p, &
                              cosa_p, &
                              nfinedivs, &
                              xfine, &
                              wfine, &
                              radiifine, &
                              alphafine, &
                              sinfine, &
                              cosfine, &
                              nfinedivs_p, &
                              xfine_p, &
                              wfine_p, &
                              radiifine_p, &
                              alphafine_p, &
                              sinfine_p, &
                              cosfine_p, &
                              fail_lcl, &
                              message_lcl, &
                              trace_lcl)

  ! Convert output arguments
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

end subroutine fo_dtwpgeometry_master_m_fo_dtwpgeometry_master_wrap

end module FO_DTWPGEOMETRY_MASTER_M_WRAP

module FO_SSWPGEOMETRY_MASTER_M_WRAP

use iso_c_binding
use fo_sswpgeometry_master_m
use fo_wpgeometry_routines_m

! This module was auto-generated 

implicit none

contains

subroutine fo_sswpgeometry_master_m_fo_sswpgeometry_master_wrap (maxgeoms, &
                                                                 maxszas, &
                                                                 maxvzas, &
                                                                 maxazms, &
                                                                 maxlayers, &
                                                                 maxpartials, &
                                                                 maxfine, &
                                                                 do_obsgeom, &
                                                                 do_doublet, &
                                                                 do_chapman, &
                                                                 do_planpar, &
                                                                 do_enhanced_ps, &
                                                                 do_partials, &
                                                                 ngeoms, &
                                                                 nszas, &
                                                                 nvzas, &
                                                                 nazms, &
                                                                 nlayers, &
                                                                 nfine, &
                                                                 npartials, &
                                                                 partial_layeridx, &
                                                                 dtr, &
                                                                 pie, &
                                                                 vsign, &
                                                                 eradius, &
                                                                 nv_offset, &
                                                                 na_offset, &
                                                                 nd_offset, &
                                                                 heights, &
                                                                 partial_heights, &
                                                                 obsgeom_boa, &
                                                                 alpha_boa, &
                                                                 theta_boa, &
                                                                 phi_boa, &
                                                                 donadir, &
                                                                 raycon, &
                                                                 mu0, &
                                                                 mu1, &
                                                                 cosscat, &
                                                                 radii, &
                                                                 losw_paths, &
                                                                 alpha, &
                                                                 sina, &
                                                                 cosa, &
                                                                 sunpaths, &
                                                                 ntraverse, &
                                                                 chapfacs, &
                                                                 theta_all, &
                                                                 radii_p, &
                                                                 losp_paths, &
                                                                 alpha_p, &
                                                                 sina_p, &
                                                                 cosa_p, &
                                                                 sunpaths_p, &
                                                                 ntraverse_p, &
                                                                 chapfacs_p, &
                                                                 nfinedivs, &
                                                                 xfine, &
                                                                 wfine, &
                                                                 radiifine, &
                                                                 alphafine, &
                                                                 sinfine, &
                                                                 cosfine, &
                                                                 sunpathsfine, &
                                                                 ntraversefine, &
                                                                 nfinedivs_p, &
                                                                 xfine_p, &
                                                                 wfine_p, &
                                                                 radiifine_p, &
                                                                 alphafine_p, &
                                                                 sinfine_p, &
                                                                 cosfine_p, &
                                                                 sunpathsfine_p, &
                                                                 ntraversefine_p, &
                                                                 fail, &
                                                                 message_len, &
                                                                 message, &
                                                                 trace_len, &
                                                                 trace) bind(C)
  use fo_wpgeometry_routines_m

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxszas
  integer(c_int), intent(in) :: maxvzas
  integer(c_int), intent(in) :: maxazms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  logical(c_bool), intent(in) :: do_obsgeom
  logical(c_bool), intent(in) :: do_doublet
  logical(c_bool), intent(in) :: do_chapman
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), intent(in) :: do_partials
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nszas
  integer(c_int), intent(in) :: nvzas
  integer(c_int), intent(in) :: nazms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), intent(in) :: nfine
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), intent(in) :: dtr
  real(c_double), intent(in) :: pie
  real(c_double), intent(in) :: vsign
  real(c_double), intent(in) :: eradius
  integer(c_int), dimension(MAXSZAS), intent(in) :: nv_offset
  integer(c_int), dimension(MAXSZAS, MAXVZAS), intent(in) :: na_offset
  integer(c_int), dimension(MAXSZAS), intent(in) :: nd_offset
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: heights
  real(c_double), dimension(MAXPARTIALS), intent(in) :: partial_heights
  real(c_double), dimension(MAXGEOMS, 3), intent(inout) :: obsgeom_boa
  real(c_double), dimension(MAXVZAS), intent(inout) :: alpha_boa
  real(c_double), dimension(MAXSZAS), intent(inout) :: theta_boa
  real(c_double), dimension(MAXAZMS), intent(inout) :: phi_boa
  logical(c_bool), dimension(MAXGEOMS), intent(inout) :: donadir
  real(c_double), dimension(MAXGEOMS), intent(inout) :: raycon
  real(c_double), dimension(MAXGEOMS), intent(out) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(out) :: mu1
  real(c_double), dimension(MAXGEOMS), intent(out) :: cosscat
  real(c_double), dimension(0:MAXLAYERS), intent(inout) :: radii
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(out) :: losw_paths
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: alpha
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: sina
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(inout) :: cosa
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(out) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(out) :: chapfacs
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: theta_all
  real(c_double), dimension(MAXPARTIALS), intent(inout) :: radii_p
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: losp_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: alpha_p
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: sina_p
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: cosa_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(out) :: sunpaths_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: ntraverse_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(out) :: chapfacs_p
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(inout) :: nfinedivs
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: radiifine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: alphafine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: sinfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(inout) :: cosfine
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(out) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(out) :: ntraversefine
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(inout) :: nfinedivs_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: radiifine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: alphafine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: sinfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(inout) :: cosfine_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(out) :: sunpathsfine_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(out) :: ntraversefine_p
  logical(c_bool), intent(out) :: fail
  integer(c_int), intent(in) :: message_len
  character(kind=c_char) , intent(inout) :: message(message_len+1)
  integer(c_int), intent(in) :: trace_len
  character(kind=c_char) , intent(inout) :: trace(trace_len+1)

  ! Local variables
  logical(kind=4) :: do_obsgeom_lcl
  logical(kind=4) :: do_doublet_lcl
  logical(kind=4) :: do_chapman_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4), dimension(MAXGEOMS) :: donadir_lcl
  logical(kind=4) :: fail_lcl
  character(kind=c_char, len=message_len) :: message_lcl
  integer :: len_idx
  character(kind=c_char, len=trace_len) :: trace_lcl

  ! Convert input arguments
  do_obsgeom_lcl = do_obsgeom
  do_doublet_lcl = do_doublet
  do_chapman_lcl = do_chapman
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_partials_lcl = do_partials
  donadir_lcl = donadir

  call fo_sswpgeometry_master(maxgeoms, &
                              maxszas, &
                              maxvzas, &
                              maxazms, &
                              maxlayers, &
                              maxpartials, &
                              maxfine, &
                              do_obsgeom_lcl, &
                              do_doublet_lcl, &
                              do_chapman_lcl, &
                              do_planpar_lcl, &
                              do_enhanced_ps_lcl, &
                              do_partials_lcl, &
                              ngeoms, &
                              nszas, &
                              nvzas, &
                              nazms, &
                              nlayers, &
                              nfine, &
                              npartials, &
                              partial_layeridx, &
                              dtr, &
                              pie, &
                              vsign, &
                              eradius, &
                              nv_offset, &
                              na_offset, &
                              nd_offset, &
                              heights, &
                              partial_heights, &
                              obsgeom_boa, &
                              alpha_boa, &
                              theta_boa, &
                              phi_boa, &
                              donadir_lcl, &
                              raycon, &
                              mu0, &
                              mu1, &
                              cosscat, &
                              radii, &
                              losw_paths, &
                              alpha, &
                              sina, &
                              cosa, &
                              sunpaths, &
                              ntraverse, &
                              chapfacs, &
                              theta_all, &
                              radii_p, &
                              losp_paths, &
                              alpha_p, &
                              sina_p, &
                              cosa_p, &
                              sunpaths_p, &
                              ntraverse_p, &
                              chapfacs_p, &
                              nfinedivs, &
                              xfine, &
                              wfine, &
                              radiifine, &
                              alphafine, &
                              sinfine, &
                              cosfine, &
                              sunpathsfine, &
                              ntraversefine, &
                              nfinedivs_p, &
                              xfine_p, &
                              wfine_p, &
                              radiifine_p, &
                              alphafine_p, &
                              sinfine_p, &
                              cosfine_p, &
                              sunpathsfine_p, &
                              ntraversefine_p, &
                              fail_lcl, &
                              message_lcl, &
                              trace_lcl)

  ! Convert output arguments
  donadir = donadir_lcl
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

end subroutine fo_sswpgeometry_master_m_fo_sswpgeometry_master_wrap

end module FO_SSWPGEOMETRY_MASTER_M_WRAP

module FO_SCALARSS_RTCALCS_I_M_WRAP

use iso_c_binding
use fo_scalarss_rtcalcs_i_m

! This module was auto-generated 

implicit none

contains

subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_dn_wrap (maxgeoms, &
                                                          maxlayers, &
                                                          maxpartials, &
                                                          maxfine, &
                                                          maxmoments_input, &
                                                          max_user_levels, &
                                                          do_deltam_scaling, &
                                                          do_phasfunc, &
                                                          do_partials, &
                                                          do_planpar, &
                                                          do_enhanced_ps, &
                                                          flux, &
                                                          do_sources_dn, &
                                                          do_sources_dn_p, &
                                                          ngeoms, &
                                                          nlayers, &
                                                          nfinedivs, &
                                                          nmoments_input, &
                                                          n_user_levels, &
                                                          user_levels, &
                                                          npartials, &
                                                          nfinedivs_p, &
                                                          partial_outindex, &
                                                          partial_outflag, &
                                                          partial_layeridx, &
                                                          extinction, &
                                                          deltaus, &
                                                          omega, &
                                                          truncfac, &
                                                          phasmoms, &
                                                          phasfunc_dn, &
                                                          mu1, &
                                                          legpoly_dn, &
                                                          losw_paths, &
                                                          losp_paths, &
                                                          xfine, &
                                                          wfine, &
                                                          sunpaths, &
                                                          ntraverse, &
                                                          sunpathsfine, &
                                                          ntraversefine, &
                                                          xfine_p, &
                                                          wfine_p, &
                                                          sunpaths_p, &
                                                          ntraverse_p, &
                                                          sunpathsfine_p, &
                                                          ntraversefine_p, &
                                                          intensity_dn, &
                                                          cumsource_dn, &
                                                          lostrans_dn) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: maxmoments_input
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_phasfunc
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  real(c_double), intent(in) :: flux
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: nmoments_input
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT), intent(in) :: phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_dn
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_dn
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dn
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_dn
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_dn

  ! Local variables
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_phasfunc_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_deltam_scaling_lcl = do_deltam_scaling
  do_phasfunc_lcl = do_phasfunc
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  partial_outflag_lcl = partial_outflag

  call ss_integral_i_dn(maxgeoms, &
                        maxlayers, &
                        maxpartials, &
                        maxfine, &
                        maxmoments_input, &
                        max_user_levels, &
                        do_deltam_scaling_lcl, &
                        do_phasfunc_lcl, &
                        do_partials_lcl, &
                        do_planpar_lcl, &
                        do_enhanced_ps_lcl, &
                        flux, &
                        do_sources_dn_lcl, &
                        do_sources_dn_p_lcl, &
                        ngeoms, &
                        nlayers, &
                        nfinedivs, &
                        nmoments_input, &
                        n_user_levels, &
                        user_levels, &
                        npartials, &
                        nfinedivs_p, &
                        partial_outindex, &
                        partial_outflag_lcl, &
                        partial_layeridx, &
                        extinction, &
                        deltaus, &
                        omega, &
                        truncfac, &
                        phasmoms, &
                        phasfunc_dn, &
                        mu1, &
                        legpoly_dn, &
                        losw_paths, &
                        losp_paths, &
                        xfine, &
                        wfine, &
                        sunpaths, &
                        ntraverse, &
                        sunpathsfine, &
                        ntraversefine, &
                        xfine_p, &
                        wfine_p, &
                        sunpaths_p, &
                        ntraverse_p, &
                        sunpathsfine_p, &
                        ntraversefine_p, &
                        intensity_dn, &
                        cumsource_dn, &
                        lostrans_dn)

end subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_dn_wrap

subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap (maxgeoms, &
                                                          maxlayers, &
                                                          maxpartials, &
                                                          maxfine, &
                                                          maxmoments_input, &
                                                          max_user_levels, &
                                                          do_deltam_scaling, &
                                                          do_phasfunc, &
                                                          do_surface_leaving, &
                                                          do_water_leaving, &
                                                          do_partials, &
                                                          do_planpar, &
                                                          do_enhanced_ps, &
                                                          flux, &
                                                          do_sources_up, &
                                                          do_sources_up_p, &
                                                          ngeoms, &
                                                          nlayers, &
                                                          nfinedivs, &
                                                          nmoments_input, &
                                                          n_user_levels, &
                                                          user_levels, &
                                                          npartials, &
                                                          nfinedivs_p, &
                                                          partial_outindex, &
                                                          partial_outflag, &
                                                          partial_layeridx, &
                                                          extinction, &
                                                          deltaus, &
                                                          omega, &
                                                          truncfac, &
                                                          phasmoms, &
                                                          phasfunc_up, &
                                                          reflec, &
                                                          slterm, &
                                                          mu0, &
                                                          mu1, &
                                                          legpoly_up, &
                                                          losw_paths, &
                                                          losp_paths, &
                                                          xfine, &
                                                          wfine, &
                                                          sunpaths, &
                                                          ntraverse, &
                                                          sunpathsfine, &
                                                          ntraversefine, &
                                                          xfine_p, &
                                                          wfine_p, &
                                                          sunpaths_p, &
                                                          ntraverse_p, &
                                                          sunpathsfine_p, &
                                                          ntraversefine_p, &
                                                          intensity_up, &
                                                          intensity_db, &
                                                          cumsource_up, &
                                                          cumtrans, &
                                                          lostrans_up) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: maxmoments_input
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_phasfunc
  logical(c_bool), intent(in) :: do_surface_leaving
  logical(c_bool), intent(in) :: do_water_leaving
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  real(c_double), intent(in) :: flux
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: nmoments_input
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT), intent(in) :: phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_up
  real(c_double), dimension(MAXGEOMS), intent(in) :: reflec
  real(c_double), dimension(MAXGEOMS), intent(in) :: slterm
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_up
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_db
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: cumtrans
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_up

  ! Local variables
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_phasfunc_lcl
  logical(kind=4) :: do_surface_leaving_lcl
  logical(kind=4) :: do_water_leaving_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_deltam_scaling_lcl = do_deltam_scaling
  do_phasfunc_lcl = do_phasfunc
  do_surface_leaving_lcl = do_surface_leaving
  do_water_leaving_lcl = do_water_leaving
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  partial_outflag_lcl = partial_outflag

  call ss_integral_i_up(maxgeoms, &
                        maxlayers, &
                        maxpartials, &
                        maxfine, &
                        maxmoments_input, &
                        max_user_levels, &
                        do_deltam_scaling_lcl, &
                        do_phasfunc_lcl, &
                        do_surface_leaving_lcl, &
                        do_water_leaving_lcl, &
                        do_partials_lcl, &
                        do_planpar_lcl, &
                        do_enhanced_ps_lcl, &
                        flux, &
                        do_sources_up_lcl, &
                        do_sources_up_p_lcl, &
                        ngeoms, &
                        nlayers, &
                        nfinedivs, &
                        nmoments_input, &
                        n_user_levels, &
                        user_levels, &
                        npartials, &
                        nfinedivs_p, &
                        partial_outindex, &
                        partial_outflag_lcl, &
                        partial_layeridx, &
                        extinction, &
                        deltaus, &
                        omega, &
                        truncfac, &
                        phasmoms, &
                        phasfunc_up, &
                        reflec, &
                        slterm, &
                        mu0, &
                        mu1, &
                        legpoly_up, &
                        losw_paths, &
                        losp_paths, &
                        xfine, &
                        wfine, &
                        sunpaths, &
                        ntraverse, &
                        sunpathsfine, &
                        ntraversefine, &
                        xfine_p, &
                        wfine_p, &
                        sunpaths_p, &
                        ntraverse_p, &
                        sunpathsfine_p, &
                        ntraversefine_p, &
                        intensity_up, &
                        intensity_db, &
                        cumsource_up, &
                        cumtrans, &
                        lostrans_up)

end subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap

subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_updn_wrap (maxgeoms, &
                                                            maxlayers, &
                                                            maxpartials, &
                                                            maxfine, &
                                                            maxmoments_input, &
                                                            max_user_levels, &
                                                            do_upwelling, &
                                                            do_dnwelling, &
                                                            do_deltam_scaling, &
                                                            do_phasfunc, &
                                                            do_surface_leaving, &
                                                            do_water_leaving, &
                                                            do_partials, &
                                                            do_planpar, &
                                                            do_enhanced_ps, &
                                                            do_sources_up, &
                                                            do_sources_up_p, &
                                                            do_sources_dn, &
                                                            do_sources_dn_p, &
                                                            ngeoms, &
                                                            nlayers, &
                                                            nfinedivs, &
                                                            nmoments_input, &
                                                            n_user_levels, &
                                                            user_levels, &
                                                            npartials, &
                                                            nfinedivs_p, &
                                                            partial_outindex, &
                                                            partial_outflag, &
                                                            partial_layeridx, &
                                                            flux, &
                                                            extinction, &
                                                            deltaus, &
                                                            omega, &
                                                            truncfac, &
                                                            phasmoms, &
                                                            phasfunc_up, &
                                                            phasfunc_dn, &
                                                            reflec, &
                                                            slterm, &
                                                            mu0, &
                                                            mu1, &
                                                            legpoly_up, &
                                                            legpoly_dn, &
                                                            losw_paths, &
                                                            losp_paths, &
                                                            xfine_up, &
                                                            wfine_up, &
                                                            sunpaths_up, &
                                                            ntraverse_up, &
                                                            sunpathsfine_up, &
                                                            ntraversefine_up, &
                                                            xfine_dn, &
                                                            wfine_dn, &
                                                            sunpaths_dn, &
                                                            ntraverse_dn, &
                                                            sunpathsfine_dn, &
                                                            ntraversefine_dn, &
                                                            xfine_up_p, &
                                                            wfine_up_p, &
                                                            sunpaths_up_p, &
                                                            ntraverse_up_p, &
                                                            sunpathsfine_up_p, &
                                                            ntraversefine_up_p, &
                                                            xfine_dn_p, &
                                                            wfine_dn_p, &
                                                            sunpaths_dn_p, &
                                                            ntraverse_dn_p, &
                                                            sunpathsfine_dn_p, &
                                                            ntraversefine_dn_p, &
                                                            intensity_up, &
                                                            intensity_db, &
                                                            cumsource_up, &
                                                            cumtrans, &
                                                            lostrans_up, &
                                                            intensity_dn, &
                                                            cumsource_dn, &
                                                            lostrans_dn) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: maxmoments_input
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_phasfunc
  logical(c_bool), intent(in) :: do_surface_leaving
  logical(c_bool), intent(in) :: do_water_leaving
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: nmoments_input
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), intent(in) :: flux
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT), intent(in) :: phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_up
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_dn
  real(c_double), dimension(MAXGEOMS), intent(in) :: reflec
  real(c_double), dimension(MAXGEOMS), intent(in) :: slterm
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_up
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_dn
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_up
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse_up
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_up
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_dn
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse_dn
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_dn
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_dn
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_up_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_up_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_up_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_dn_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_dn_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_dn_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_db
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: cumtrans
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dn
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_dn
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_dn

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_dnwelling_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_phasfunc_lcl
  logical(kind=4) :: do_surface_leaving_lcl
  logical(kind=4) :: do_water_leaving_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_deltam_scaling_lcl = do_deltam_scaling
  do_phasfunc_lcl = do_phasfunc
  do_surface_leaving_lcl = do_surface_leaving
  do_water_leaving_lcl = do_water_leaving
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  partial_outflag_lcl = partial_outflag

  call ss_integral_i_updn(maxgeoms, &
                          maxlayers, &
                          maxpartials, &
                          maxfine, &
                          maxmoments_input, &
                          max_user_levels, &
                          do_upwelling_lcl, &
                          do_dnwelling_lcl, &
                          do_deltam_scaling_lcl, &
                          do_phasfunc_lcl, &
                          do_surface_leaving_lcl, &
                          do_water_leaving_lcl, &
                          do_partials_lcl, &
                          do_planpar_lcl, &
                          do_enhanced_ps_lcl, &
                          do_sources_up_lcl, &
                          do_sources_up_p_lcl, &
                          do_sources_dn_lcl, &
                          do_sources_dn_p_lcl, &
                          ngeoms, &
                          nlayers, &
                          nfinedivs, &
                          nmoments_input, &
                          n_user_levels, &
                          user_levels, &
                          npartials, &
                          nfinedivs_p, &
                          partial_outindex, &
                          partial_outflag_lcl, &
                          partial_layeridx, &
                          flux, &
                          extinction, &
                          deltaus, &
                          omega, &
                          truncfac, &
                          phasmoms, &
                          phasfunc_up, &
                          phasfunc_dn, &
                          reflec, &
                          slterm, &
                          mu0, &
                          mu1, &
                          legpoly_up, &
                          legpoly_dn, &
                          losw_paths, &
                          losp_paths, &
                          xfine_up, &
                          wfine_up, &
                          sunpaths_up, &
                          ntraverse_up, &
                          sunpathsfine_up, &
                          ntraversefine_up, &
                          xfine_dn, &
                          wfine_dn, &
                          sunpaths_dn, &
                          ntraverse_dn, &
                          sunpathsfine_dn, &
                          ntraversefine_dn, &
                          xfine_up_p, &
                          wfine_up_p, &
                          sunpaths_up_p, &
                          ntraverse_up_p, &
                          sunpathsfine_up_p, &
                          ntraversefine_up_p, &
                          xfine_dn_p, &
                          wfine_dn_p, &
                          sunpaths_dn_p, &
                          ntraverse_dn_p, &
                          sunpathsfine_dn_p, &
                          ntraversefine_dn_p, &
                          intensity_up, &
                          intensity_db, &
                          cumsource_up, &
                          cumtrans, &
                          lostrans_up, &
                          intensity_dn, &
                          cumsource_dn, &
                          lostrans_dn)

end subroutine fo_scalarss_rtcalcs_i_m_ss_integral_i_updn_wrap

end module FO_SCALARSS_RTCALCS_I_M_WRAP

module FO_SCALARSS_RTCALCS_ILPS_M_WRAP

use iso_c_binding
use fo_scalarss_rtcalcs_ilps_m

! This module was auto-generated 

implicit none

contains

subroutine fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_dn_wrap (maxgeoms, &
                                                                maxlayers, &
                                                                maxpartials, &
                                                                maxfine, &
                                                                maxmoments_input, &
                                                                max_user_levels, &
                                                                max_atmoswfs, &
                                                                do_deltam_scaling, &
                                                                do_phasfunc, &
                                                                do_partials, &
                                                                do_planpar, &
                                                                do_enhanced_ps, &
                                                                flux, &
                                                                do_sources_dn, &
                                                                do_sources_dn_p, &
                                                                do_profilewfs, &
                                                                lvaryflags, &
                                                                lvarynums, &
                                                                lvarymoms, &
                                                                ngeoms, &
                                                                nlayers, &
                                                                nfinedivs, &
                                                                nmoments_input, &
                                                                n_user_levels, &
                                                                user_levels, &
                                                                npartials, &
                                                                nfinedivs_p, &
                                                                partial_outindex, &
                                                                partial_outflag, &
                                                                partial_layeridx, &
                                                                extinction, &
                                                                deltaus, &
                                                                omega, &
                                                                truncfac, &
                                                                phasmoms, &
                                                                phasfunc_dn, &
                                                                l_extinction, &
                                                                l_deltaus, &
                                                                l_omega, &
                                                                l_truncfac, &
                                                                l_phasmoms, &
                                                                l_phasfunc_dn, &
                                                                mu1, &
                                                                legpoly_dn, &
                                                                losw_paths, &
                                                                losp_paths, &
                                                                xfine, &
                                                                wfine, &
                                                                sunpaths, &
                                                                ntraverse, &
                                                                sunpathsfine, &
                                                                ntraversefine, &
                                                                xfine_p, &
                                                                wfine_p, &
                                                                sunpaths_p, &
                                                                ntraverse_p, &
                                                                sunpathsfine_p, &
                                                                ntraversefine_p, &
                                                                intensity_dn, &
                                                                lp_jacobians_dn, &
                                                                lostrans_dn, &
                                                                lp_lostrans_dn) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: maxmoments_input
  integer(c_int), intent(in) :: max_user_levels
  integer(c_int), intent(in) :: max_atmoswfs
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_phasfunc
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  real(c_double), intent(in) :: flux
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  logical(c_bool), intent(in) :: do_profilewfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: lvaryflags
  integer(c_int), dimension(MAXLAYERS), intent(in) :: lvarynums
  logical(c_bool), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: lvarymoms
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: nmoments_input
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT), intent(in) :: phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_dn
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_extinction
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltaus
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, MAX_ATMOSWFS), intent(in) :: l_phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS, MAX_ATMOSWFS), intent(in) :: l_phasfunc_dn
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_dn
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dn
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dn
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_dn
  real(c_double), dimension(MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_lostrans_dn

  ! Local variables
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_phasfunc_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4) :: do_profilewfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: lvaryflags_lcl
  logical(kind=4), dimension(MAXLAYERS, MAX_ATMOSWFS) :: lvarymoms_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_deltam_scaling_lcl = do_deltam_scaling
  do_phasfunc_lcl = do_phasfunc
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  do_profilewfs_lcl = do_profilewfs
  lvaryflags_lcl = lvaryflags
  lvarymoms_lcl = lvarymoms
  partial_outflag_lcl = partial_outflag

  call ss_integral_ilps_dn(maxgeoms, &
                           maxlayers, &
                           maxpartials, &
                           maxfine, &
                           maxmoments_input, &
                           max_user_levels, &
                           max_atmoswfs, &
                           do_deltam_scaling_lcl, &
                           do_phasfunc_lcl, &
                           do_partials_lcl, &
                           do_planpar_lcl, &
                           do_enhanced_ps_lcl, &
                           flux, &
                           do_sources_dn_lcl, &
                           do_sources_dn_p_lcl, &
                           do_profilewfs_lcl, &
                           lvaryflags_lcl, &
                           lvarynums, &
                           lvarymoms_lcl, &
                           ngeoms, &
                           nlayers, &
                           nfinedivs, &
                           nmoments_input, &
                           n_user_levels, &
                           user_levels, &
                           npartials, &
                           nfinedivs_p, &
                           partial_outindex, &
                           partial_outflag_lcl, &
                           partial_layeridx, &
                           extinction, &
                           deltaus, &
                           omega, &
                           truncfac, &
                           phasmoms, &
                           phasfunc_dn, &
                           l_extinction, &
                           l_deltaus, &
                           l_omega, &
                           l_truncfac, &
                           l_phasmoms, &
                           l_phasfunc_dn, &
                           mu1, &
                           legpoly_dn, &
                           losw_paths, &
                           losp_paths, &
                           xfine, &
                           wfine, &
                           sunpaths, &
                           ntraverse, &
                           sunpathsfine, &
                           ntraversefine, &
                           xfine_p, &
                           wfine_p, &
                           sunpaths_p, &
                           ntraverse_p, &
                           sunpathsfine_p, &
                           ntraversefine_p, &
                           intensity_dn, &
                           lp_jacobians_dn, &
                           lostrans_dn, &
                           lp_lostrans_dn)

end subroutine fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_dn_wrap

subroutine fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_up_wrap (maxgeoms, &
                                                                maxlayers, &
                                                                maxpartials, &
                                                                maxfine, &
                                                                maxmoments_input, &
                                                                max_user_levels, &
                                                                max_atmoswfs, &
                                                                max_surfacewfs, &
                                                                max_sleavewfs, &
                                                                do_deltam_scaling, &
                                                                do_phasfunc, &
                                                                do_surface_leaving, &
                                                                do_water_leaving, &
                                                                do_partials, &
                                                                do_planpar, &
                                                                do_enhanced_ps, &
                                                                flux, &
                                                                do_sources_up, &
                                                                do_sources_up_p, &
                                                                do_profilewfs, &
                                                                do_surfacewfs, &
                                                                do_sleavewfs, &
                                                                n_reflecwfs, &
                                                                n_sleavewfs, &
                                                                n_surfacewfs, &
                                                                lvaryflags, &
                                                                lvarynums, &
                                                                lvarymoms, &
                                                                ngeoms, &
                                                                nlayers, &
                                                                nfinedivs, &
                                                                nmoments_input, &
                                                                n_user_levels, &
                                                                user_levels, &
                                                                npartials, &
                                                                nfinedivs_p, &
                                                                partial_outindex, &
                                                                partial_outflag, &
                                                                partial_layeridx, &
                                                                extinction, &
                                                                deltaus, &
                                                                omega, &
                                                                truncfac, &
                                                                phasmoms, &
                                                                phasfunc_up, &
                                                                reflec, &
                                                                slterm, &
                                                                l_extinction, &
                                                                l_deltaus, &
                                                                l_omega, &
                                                                l_truncfac, &
                                                                l_phasmoms, &
                                                                l_phasfunc_up, &
                                                                ls_reflec, &
                                                                lssl_slterm, &
                                                                mu0, &
                                                                mu1, &
                                                                legpoly_up, &
                                                                losw_paths, &
                                                                losp_paths, &
                                                                xfine, &
                                                                wfine, &
                                                                sunpaths, &
                                                                ntraverse, &
                                                                sunpathsfine, &
                                                                ntraversefine, &
                                                                xfine_p, &
                                                                wfine_p, &
                                                                sunpaths_p, &
                                                                ntraverse_p, &
                                                                sunpathsfine_p, &
                                                                ntraversefine_p, &
                                                                intensity_up, &
                                                                intensity_db, &
                                                                lp_jacobians_up, &
                                                                lp_jacobians_db, &
                                                                ls_jacobians_db, &
                                                                cumtrans, &
                                                                lostrans_up, &
                                                                lp_cumtrans, &
                                                                lp_lostrans_up) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: maxmoments_input
  integer(c_int), intent(in) :: max_user_levels
  integer(c_int), intent(in) :: max_atmoswfs
  integer(c_int), intent(in) :: max_surfacewfs
  integer(c_int), intent(in) :: max_sleavewfs
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_phasfunc
  logical(c_bool), intent(in) :: do_surface_leaving
  logical(c_bool), intent(in) :: do_water_leaving
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  real(c_double), intent(in) :: flux
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  logical(c_bool), intent(in) :: do_profilewfs
  logical(c_bool), intent(in) :: do_surfacewfs
  logical(c_bool), intent(in) :: do_sleavewfs
  integer(c_int), intent(in) :: n_reflecwfs
  integer(c_int), intent(in) :: n_sleavewfs
  integer(c_int), intent(in) :: n_surfacewfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: lvaryflags
  integer(c_int), dimension(MAXLAYERS), intent(in) :: lvarynums
  logical(c_bool), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: lvarymoms
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: nmoments_input
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT), intent(in) :: phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_up
  real(c_double), dimension(MAXGEOMS), intent(in) :: reflec
  real(c_double), dimension(MAXGEOMS), intent(in) :: slterm
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_extinction
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltaus
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, MAX_ATMOSWFS), intent(in) :: l_phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS, MAX_ATMOSWFS), intent(in) :: l_phasfunc_up
  real(c_double), dimension(MAXGEOMS, MAX_SURFACEWFS), intent(in) :: ls_reflec
  real(c_double), dimension(MAXGEOMS, MAX_SLEAVEWFS), intent(in) :: lssl_slterm
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_up
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_db
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_db
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAX_SURFACEWFS), intent(out) :: ls_jacobians_db
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: cumtrans
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_cumtrans
  real(c_double), dimension(MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_lostrans_up

  ! Local variables
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_phasfunc_lcl
  logical(kind=4) :: do_surface_leaving_lcl
  logical(kind=4) :: do_water_leaving_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4) :: do_profilewfs_lcl
  logical(kind=4) :: do_surfacewfs_lcl
  logical(kind=4) :: do_sleavewfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: lvaryflags_lcl
  logical(kind=4), dimension(MAXLAYERS, MAX_ATMOSWFS) :: lvarymoms_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_deltam_scaling_lcl = do_deltam_scaling
  do_phasfunc_lcl = do_phasfunc
  do_surface_leaving_lcl = do_surface_leaving
  do_water_leaving_lcl = do_water_leaving
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  do_profilewfs_lcl = do_profilewfs
  do_surfacewfs_lcl = do_surfacewfs
  do_sleavewfs_lcl = do_sleavewfs
  lvaryflags_lcl = lvaryflags
  lvarymoms_lcl = lvarymoms
  partial_outflag_lcl = partial_outflag

  call ss_integral_ilps_up(maxgeoms, &
                           maxlayers, &
                           maxpartials, &
                           maxfine, &
                           maxmoments_input, &
                           max_user_levels, &
                           max_atmoswfs, &
                           max_surfacewfs, &
                           max_sleavewfs, &
                           do_deltam_scaling_lcl, &
                           do_phasfunc_lcl, &
                           do_surface_leaving_lcl, &
                           do_water_leaving_lcl, &
                           do_partials_lcl, &
                           do_planpar_lcl, &
                           do_enhanced_ps_lcl, &
                           flux, &
                           do_sources_up_lcl, &
                           do_sources_up_p_lcl, &
                           do_profilewfs_lcl, &
                           do_surfacewfs_lcl, &
                           do_sleavewfs_lcl, &
                           n_reflecwfs, &
                           n_sleavewfs, &
                           n_surfacewfs, &
                           lvaryflags_lcl, &
                           lvarynums, &
                           lvarymoms_lcl, &
                           ngeoms, &
                           nlayers, &
                           nfinedivs, &
                           nmoments_input, &
                           n_user_levels, &
                           user_levels, &
                           npartials, &
                           nfinedivs_p, &
                           partial_outindex, &
                           partial_outflag_lcl, &
                           partial_layeridx, &
                           extinction, &
                           deltaus, &
                           omega, &
                           truncfac, &
                           phasmoms, &
                           phasfunc_up, &
                           reflec, &
                           slterm, &
                           l_extinction, &
                           l_deltaus, &
                           l_omega, &
                           l_truncfac, &
                           l_phasmoms, &
                           l_phasfunc_up, &
                           ls_reflec, &
                           lssl_slterm, &
                           mu0, &
                           mu1, &
                           legpoly_up, &
                           losw_paths, &
                           losp_paths, &
                           xfine, &
                           wfine, &
                           sunpaths, &
                           ntraverse, &
                           sunpathsfine, &
                           ntraversefine, &
                           xfine_p, &
                           wfine_p, &
                           sunpaths_p, &
                           ntraverse_p, &
                           sunpathsfine_p, &
                           ntraversefine_p, &
                           intensity_up, &
                           intensity_db, &
                           lp_jacobians_up, &
                           lp_jacobians_db, &
                           ls_jacobians_db, &
                           cumtrans, &
                           lostrans_up, &
                           lp_cumtrans, &
                           lp_lostrans_up)

end subroutine fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_up_wrap

subroutine fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_updn_wrap (maxgeoms, &
                                                                  maxlayers, &
                                                                  maxpartials, &
                                                                  maxfine, &
                                                                  maxmoments_input, &
                                                                  max_user_levels, &
                                                                  max_atmoswfs, &
                                                                  max_surfacewfs, &
                                                                  max_sleavewfs, &
                                                                  do_upwelling, &
                                                                  do_dnwelling, &
                                                                  do_deltam_scaling, &
                                                                  do_phasfunc, &
                                                                  do_surface_leaving, &
                                                                  do_water_leaving, &
                                                                  do_partials, &
                                                                  do_planpar, &
                                                                  do_enhanced_ps, &
                                                                  do_sources_up, &
                                                                  do_sources_up_p, &
                                                                  do_sources_dn, &
                                                                  do_sources_dn_p, &
                                                                  do_profilewfs, &
                                                                  do_surfacewfs, &
                                                                  do_sleavewfs, &
                                                                  n_reflecwfs, &
                                                                  n_sleavewfs, &
                                                                  n_surfacewfs, &
                                                                  lvaryflags, &
                                                                  lvarynums, &
                                                                  lvarymoms, &
                                                                  ngeoms, &
                                                                  nlayers, &
                                                                  nfinedivs, &
                                                                  nmoments_input, &
                                                                  n_user_levels, &
                                                                  user_levels, &
                                                                  npartials, &
                                                                  nfinedivs_p, &
                                                                  partial_outindex, &
                                                                  partial_outflag, &
                                                                  partial_layeridx, &
                                                                  flux, &
                                                                  extinction, &
                                                                  deltaus, &
                                                                  omega, &
                                                                  truncfac, &
                                                                  phasmoms, &
                                                                  phasfunc_up, &
                                                                  phasfunc_dn, &
                                                                  reflec, &
                                                                  slterm, &
                                                                  ls_reflec, &
                                                                  lssl_slterm, &
                                                                  l_extinction, &
                                                                  l_deltaus, &
                                                                  l_omega, &
                                                                  l_truncfac, &
                                                                  l_phasmoms, &
                                                                  l_phasfunc_up, &
                                                                  l_phasfunc_dn, &
                                                                  mu0, &
                                                                  mu1, &
                                                                  legpoly_up, &
                                                                  legpoly_dn, &
                                                                  losw_paths, &
                                                                  losp_paths, &
                                                                  xfine_up, &
                                                                  wfine_up, &
                                                                  sunpaths_up, &
                                                                  ntraverse_up, &
                                                                  sunpathsfine_up, &
                                                                  ntraversefine_up, &
                                                                  xfine_dn, &
                                                                  wfine_dn, &
                                                                  sunpaths_dn, &
                                                                  ntraverse_dn, &
                                                                  sunpathsfine_dn, &
                                                                  ntraversefine_dn, &
                                                                  xfine_up_p, &
                                                                  wfine_up_p, &
                                                                  sunpaths_up_p, &
                                                                  ntraverse_up_p, &
                                                                  sunpathsfine_up_p, &
                                                                  ntraversefine_up_p, &
                                                                  xfine_dn_p, &
                                                                  wfine_dn_p, &
                                                                  sunpaths_dn_p, &
                                                                  ntraverse_dn_p, &
                                                                  sunpathsfine_dn_p, &
                                                                  ntraversefine_dn_p, &
                                                                  intensity_up, &
                                                                  intensity_db, &
                                                                  lp_jacobians_up, &
                                                                  lp_jacobians_db, &
                                                                  ls_jacobians_db, &
                                                                  lostrans_up, &
                                                                  lp_lostrans_up, &
                                                                  intensity_dn, &
                                                                  lp_jacobians_dn, &
                                                                  cumtrans, &
                                                                  lp_cumtrans, &
                                                                  lostrans_dn, &
                                                                  lp_lostrans_dn) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: maxmoments_input
  integer(c_int), intent(in) :: max_user_levels
  integer(c_int), intent(in) :: max_atmoswfs
  integer(c_int), intent(in) :: max_surfacewfs
  integer(c_int), intent(in) :: max_sleavewfs
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_phasfunc
  logical(c_bool), intent(in) :: do_surface_leaving
  logical(c_bool), intent(in) :: do_water_leaving
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  logical(c_bool), intent(in) :: do_profilewfs
  logical(c_bool), intent(in) :: do_surfacewfs
  logical(c_bool), intent(in) :: do_sleavewfs
  integer(c_int), intent(in) :: n_reflecwfs
  integer(c_int), intent(in) :: n_sleavewfs
  integer(c_int), intent(in) :: n_surfacewfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: lvaryflags
  integer(c_int), dimension(MAXLAYERS), intent(in) :: lvarynums
  logical(c_bool), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: lvarymoms
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: nmoments_input
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), intent(in) :: flux
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT), intent(in) :: phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_up
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: phasfunc_dn
  real(c_double), dimension(MAXGEOMS), intent(in) :: reflec
  real(c_double), dimension(MAXGEOMS), intent(in) :: slterm
  real(c_double), dimension(MAXGEOMS, MAX_SURFACEWFS), intent(in) :: ls_reflec
  real(c_double), dimension(MAXGEOMS, MAX_SLEAVEWFS), intent(in) :: lssl_slterm
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_extinction
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltaus
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_truncfac
  real(c_double), dimension(MAXLAYERS, 0:MAXMOMENTS_INPUT, MAX_ATMOSWFS), intent(in) :: l_phasmoms
  real(c_double), dimension(MAXLAYERS, MAXGEOMS, MAX_ATMOSWFS), intent(in) :: l_phasfunc_up
  real(c_double), dimension(MAXLAYERS, MAXGEOMS, MAX_ATMOSWFS), intent(in) :: l_phasfunc_dn
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu0
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_up
  real(c_double), dimension(0:MAXMOMENTS_INPUT, MAXGEOMS), intent(in) :: legpoly_dn
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_up
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse_up
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_up
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn
  real(c_double), dimension(0:MAXLAYERS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_dn
  integer(c_int), dimension(0:MAXLAYERS, MAXGEOMS), intent(in) :: ntraverse_dn
  real(c_double), dimension(MAXLAYERS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_dn
  integer(c_int), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_dn
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_up_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_up_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_up_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXGEOMS), intent(in) :: sunpaths_dn_p
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: ntraverse_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: sunpathsfine_dn_p
  integer(c_int), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: ntraversefine_dn_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_db
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_db
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAX_SURFACEWFS), intent(out) :: ls_jacobians_db
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_up
  real(c_double), dimension(MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_lostrans_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dn
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dn
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: cumtrans
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_cumtrans
  real(c_double), dimension(MAXGEOMS, MAXLAYERS), intent(out) :: lostrans_dn
  real(c_double), dimension(MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_lostrans_dn

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_dnwelling_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_phasfunc_lcl
  logical(kind=4) :: do_surface_leaving_lcl
  logical(kind=4) :: do_water_leaving_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4) :: do_profilewfs_lcl
  logical(kind=4) :: do_surfacewfs_lcl
  logical(kind=4) :: do_sleavewfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: lvaryflags_lcl
  logical(kind=4), dimension(MAXLAYERS, MAX_ATMOSWFS) :: lvarymoms_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_deltam_scaling_lcl = do_deltam_scaling
  do_phasfunc_lcl = do_phasfunc
  do_surface_leaving_lcl = do_surface_leaving
  do_water_leaving_lcl = do_water_leaving
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  do_profilewfs_lcl = do_profilewfs
  do_surfacewfs_lcl = do_surfacewfs
  do_sleavewfs_lcl = do_sleavewfs
  lvaryflags_lcl = lvaryflags
  lvarymoms_lcl = lvarymoms
  partial_outflag_lcl = partial_outflag

  call ss_integral_ilps_updn(maxgeoms, &
                             maxlayers, &
                             maxpartials, &
                             maxfine, &
                             maxmoments_input, &
                             max_user_levels, &
                             max_atmoswfs, &
                             max_surfacewfs, &
                             max_sleavewfs, &
                             do_upwelling_lcl, &
                             do_dnwelling_lcl, &
                             do_deltam_scaling_lcl, &
                             do_phasfunc_lcl, &
                             do_surface_leaving_lcl, &
                             do_water_leaving_lcl, &
                             do_partials_lcl, &
                             do_planpar_lcl, &
                             do_enhanced_ps_lcl, &
                             do_sources_up_lcl, &
                             do_sources_up_p_lcl, &
                             do_sources_dn_lcl, &
                             do_sources_dn_p_lcl, &
                             do_profilewfs_lcl, &
                             do_surfacewfs_lcl, &
                             do_sleavewfs_lcl, &
                             n_reflecwfs, &
                             n_sleavewfs, &
                             n_surfacewfs, &
                             lvaryflags_lcl, &
                             lvarynums, &
                             lvarymoms_lcl, &
                             ngeoms, &
                             nlayers, &
                             nfinedivs, &
                             nmoments_input, &
                             n_user_levels, &
                             user_levels, &
                             npartials, &
                             nfinedivs_p, &
                             partial_outindex, &
                             partial_outflag_lcl, &
                             partial_layeridx, &
                             flux, &
                             extinction, &
                             deltaus, &
                             omega, &
                             truncfac, &
                             phasmoms, &
                             phasfunc_up, &
                             phasfunc_dn, &
                             reflec, &
                             slterm, &
                             ls_reflec, &
                             lssl_slterm, &
                             l_extinction, &
                             l_deltaus, &
                             l_omega, &
                             l_truncfac, &
                             l_phasmoms, &
                             l_phasfunc_up, &
                             l_phasfunc_dn, &
                             mu0, &
                             mu1, &
                             legpoly_up, &
                             legpoly_dn, &
                             losw_paths, &
                             losp_paths, &
                             xfine_up, &
                             wfine_up, &
                             sunpaths_up, &
                             ntraverse_up, &
                             sunpathsfine_up, &
                             ntraversefine_up, &
                             xfine_dn, &
                             wfine_dn, &
                             sunpaths_dn, &
                             ntraverse_dn, &
                             sunpathsfine_dn, &
                             ntraversefine_dn, &
                             xfine_up_p, &
                             wfine_up_p, &
                             sunpaths_up_p, &
                             ntraverse_up_p, &
                             sunpathsfine_up_p, &
                             ntraversefine_up_p, &
                             xfine_dn_p, &
                             wfine_dn_p, &
                             sunpaths_dn_p, &
                             ntraverse_dn_p, &
                             sunpathsfine_dn_p, &
                             ntraversefine_dn_p, &
                             intensity_up, &
                             intensity_db, &
                             lp_jacobians_up, &
                             lp_jacobians_db, &
                             ls_jacobians_db, &
                             lostrans_up, &
                             lp_lostrans_up, &
                             intensity_dn, &
                             lp_jacobians_dn, &
                             cumtrans, &
                             lp_cumtrans, &
                             lostrans_dn, &
                             lp_lostrans_dn)

end subroutine fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_updn_wrap

end module FO_SCALARSS_RTCALCS_ILPS_M_WRAP

module FO_SCALARSS_SPHERFUNCS_M_WRAP

use iso_c_binding
use fo_scalarss_spherfuncs_m

! This module was auto-generated 

implicit none

contains

subroutine fo_scalarss_spherfuncs_m_fo_scalarss_spherfuncs_wrap (maxmoments, &
                                                                 maxgeoms, &
                                                                 nmoments, &
                                                                 ngeoms, &
                                                                 starter, &
                                                                 do_spherfunc, &
                                                                 cosscat, &
                                                                 df1, &
                                                                 df2, &
                                                                 ss_pleg) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxmoments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: nmoments
  integer(c_int), intent(in) :: ngeoms
  logical(c_bool), intent(inout) :: starter
  logical(c_bool), intent(in) :: do_spherfunc
  real(c_double), dimension(MAXGEOMS), intent(in) :: cosscat
  real(c_double), dimension(MAXMOMENTS), intent(inout) :: df1
  real(c_double), dimension(MAXMOMENTS), intent(inout) :: df2
  real(c_double), dimension(0:MAXMOMENTS, MAXGEOMS), intent(out) :: ss_pleg

  ! Local variables
  logical(kind=4) :: starter_lcl
  logical(kind=4) :: do_spherfunc_lcl

  ! Convert input arguments
  starter_lcl = starter
  do_spherfunc_lcl = do_spherfunc

  call fo_scalarss_spherfuncs(maxmoments, &
                              maxgeoms, &
                              nmoments, &
                              ngeoms, &
                              starter_lcl, &
                              do_spherfunc_lcl, &
                              cosscat, &
                              df1, &
                              df2, &
                              ss_pleg)

  ! Convert output arguments
  starter = starter_lcl

end subroutine fo_scalarss_spherfuncs_m_fo_scalarss_spherfuncs_wrap

end module FO_SCALARSS_SPHERFUNCS_M_WRAP

module FO_THERMAL_RTCALCS_I_M_WRAP

use iso_c_binding
use fo_thermal_rtcalcs_i_m

! This module was auto-generated 

implicit none

contains

subroutine fo_thermal_rtcalcs_i_m_dte_integral_i_dn_wrap (maxgeoms, &
                                                          maxlayers, &
                                                          maxpartials, &
                                                          maxfine, &
                                                          max_user_levels, &
                                                          do_thermset, &
                                                          do_deltam_scaling, &
                                                          do_partials, &
                                                          do_planpar, &
                                                          do_enhanced_ps, &
                                                          do_sources_dn, &
                                                          do_sources_dn_p, &
                                                          ngeoms, &
                                                          nlayers, &
                                                          nfinedivs, &
                                                          n_user_levels, &
                                                          user_levels, &
                                                          npartials, &
                                                          nfinedivs_p, &
                                                          partial_outindex, &
                                                          partial_outflag, &
                                                          partial_layeridx, &
                                                          bb_input, &
                                                          extinction, &
                                                          deltaus, &
                                                          omega, &
                                                          truncfac, &
                                                          mu1, &
                                                          losw_paths, &
                                                          losp_paths, &
                                                          xfine, &
                                                          wfine, &
                                                          hfine, &
                                                          xfine_p, &
                                                          wfine_p, &
                                                          hfine_p, &
                                                          intensity_dta_dn, &
                                                          cumsource_dn, &
                                                          tcom1) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(inout) :: do_thermset
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: bb_input
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_dn
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_dn
  real(c_double), dimension(MAXLAYERS, 2), intent(inout) :: tcom1

  ! Local variables
  logical(kind=4) :: do_thermset_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_thermset_lcl = do_thermset
  do_deltam_scaling_lcl = do_deltam_scaling
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  partial_outflag_lcl = partial_outflag

  call dte_integral_i_dn(maxgeoms, &
                         maxlayers, &
                         maxpartials, &
                         maxfine, &
                         max_user_levels, &
                         do_thermset_lcl, &
                         do_deltam_scaling_lcl, &
                         do_partials_lcl, &
                         do_planpar_lcl, &
                         do_enhanced_ps_lcl, &
                         do_sources_dn_lcl, &
                         do_sources_dn_p_lcl, &
                         ngeoms, &
                         nlayers, &
                         nfinedivs, &
                         n_user_levels, &
                         user_levels, &
                         npartials, &
                         nfinedivs_p, &
                         partial_outindex, &
                         partial_outflag_lcl, &
                         partial_layeridx, &
                         bb_input, &
                         extinction, &
                         deltaus, &
                         omega, &
                         truncfac, &
                         mu1, &
                         losw_paths, &
                         losp_paths, &
                         xfine, &
                         wfine, &
                         hfine, &
                         xfine_p, &
                         wfine_p, &
                         hfine_p, &
                         intensity_dta_dn, &
                         cumsource_dn, &
                         tcom1)

  ! Convert output arguments
  do_thermset = do_thermset_lcl

end subroutine fo_thermal_rtcalcs_i_m_dte_integral_i_dn_wrap

subroutine fo_thermal_rtcalcs_i_m_dte_integral_i_up_wrap (maxgeoms, &
                                                          maxlayers, &
                                                          maxpartials, &
                                                          maxfine, &
                                                          max_user_levels, &
                                                          do_thermset, &
                                                          do_deltam_scaling, &
                                                          do_partials, &
                                                          do_planpar, &
                                                          do_enhanced_ps, &
                                                          do_sources_up, &
                                                          do_sources_up_p, &
                                                          ngeoms, &
                                                          nlayers, &
                                                          nfinedivs, &
                                                          n_user_levels, &
                                                          user_levels, &
                                                          npartials, &
                                                          nfinedivs_p, &
                                                          partial_outindex, &
                                                          partial_outflag, &
                                                          partial_layeridx, &
                                                          bb_input, &
                                                          surfbb, &
                                                          user_emissivity, &
                                                          extinction, &
                                                          deltaus, &
                                                          omega, &
                                                          truncfac, &
                                                          mu1, &
                                                          losw_paths, &
                                                          losp_paths, &
                                                          xfine, &
                                                          wfine, &
                                                          hfine, &
                                                          xfine_p, &
                                                          wfine_p, &
                                                          hfine_p, &
                                                          intensity_dta_up, &
                                                          intensity_dts, &
                                                          cumsource_up, &
                                                          tcom1, &
                                                          lostrans_up, &
                                                          lostrans_up_p) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(inout) :: do_thermset
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: bb_input
  real(c_double), intent(in) :: surfbb
  real(c_double), dimension(MAXGEOMS), intent(in) :: user_emissivity
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dts
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_up
  real(c_double), dimension(MAXLAYERS, 2), intent(inout) :: tcom1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(out) :: lostrans_up
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: lostrans_up_p

  ! Local variables
  logical(kind=4) :: do_thermset_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_thermset_lcl = do_thermset
  do_deltam_scaling_lcl = do_deltam_scaling
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  partial_outflag_lcl = partial_outflag

  call dte_integral_i_up(maxgeoms, &
                         maxlayers, &
                         maxpartials, &
                         maxfine, &
                         max_user_levels, &
                         do_thermset_lcl, &
                         do_deltam_scaling_lcl, &
                         do_partials_lcl, &
                         do_planpar_lcl, &
                         do_enhanced_ps_lcl, &
                         do_sources_up_lcl, &
                         do_sources_up_p_lcl, &
                         ngeoms, &
                         nlayers, &
                         nfinedivs, &
                         n_user_levels, &
                         user_levels, &
                         npartials, &
                         nfinedivs_p, &
                         partial_outindex, &
                         partial_outflag_lcl, &
                         partial_layeridx, &
                         bb_input, &
                         surfbb, &
                         user_emissivity, &
                         extinction, &
                         deltaus, &
                         omega, &
                         truncfac, &
                         mu1, &
                         losw_paths, &
                         losp_paths, &
                         xfine, &
                         wfine, &
                         hfine, &
                         xfine_p, &
                         wfine_p, &
                         hfine_p, &
                         intensity_dta_up, &
                         intensity_dts, &
                         cumsource_up, &
                         tcom1, &
                         lostrans_up, &
                         lostrans_up_p)

  ! Convert output arguments
  do_thermset = do_thermset_lcl

end subroutine fo_thermal_rtcalcs_i_m_dte_integral_i_up_wrap

subroutine fo_thermal_rtcalcs_i_m_dte_integral_i_updn_wrap (maxgeoms, &
                                                            maxlayers, &
                                                            maxpartials, &
                                                            maxfine, &
                                                            max_user_levels, &
                                                            do_upwelling, &
                                                            do_dnwelling, &
                                                            do_thermset, &
                                                            do_deltam_scaling, &
                                                            do_partials, &
                                                            do_planpar, &
                                                            do_enhanced_ps, &
                                                            do_sources_up, &
                                                            do_sources_up_p, &
                                                            do_sources_dn, &
                                                            do_sources_dn_p, &
                                                            ngeoms, &
                                                            nlayers, &
                                                            nfinedivs, &
                                                            n_user_levels, &
                                                            user_levels, &
                                                            npartials, &
                                                            nfinedivs_p, &
                                                            partial_outindex, &
                                                            partial_outflag, &
                                                            partial_layeridx, &
                                                            bb_input, &
                                                            surfbb, &
                                                            user_emissivity, &
                                                            extinction, &
                                                            deltaus, &
                                                            omega, &
                                                            truncfac, &
                                                            mu1, &
                                                            losw_paths, &
                                                            losp_paths, &
                                                            xfine_up, &
                                                            wfine_up, &
                                                            hfine_up, &
                                                            xfine_dn, &
                                                            wfine_dn, &
                                                            hfine_dn, &
                                                            xfine_up_p, &
                                                            wfine_up_p, &
                                                            hfine_up_p, &
                                                            xfine_dn_p, &
                                                            wfine_dn_p, &
                                                            hfine_dn_p, &
                                                            intensity_dta_up, &
                                                            intensity_dts, &
                                                            intensity_dta_dn, &
                                                            cumsource_up, &
                                                            cumsource_dn, &
                                                            tcom1, &
                                                            lostrans_up, &
                                                            lostrans_up_p) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(inout) :: do_thermset
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: bb_input
  real(c_double), intent(in) :: surfbb
  real(c_double), dimension(MAXGEOMS), intent(in) :: user_emissivity
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine_dn
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_dn_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dts
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_dn
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_up
  real(c_double), dimension(0:MAXLAYERS, MAXGEOMS), intent(out) :: cumsource_dn
  real(c_double), dimension(MAXLAYERS, 2), intent(inout) :: tcom1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(out) :: lostrans_up
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: lostrans_up_p

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_dnwelling_lcl
  logical(kind=4) :: do_thermset_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_thermset_lcl = do_thermset
  do_deltam_scaling_lcl = do_deltam_scaling
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  partial_outflag_lcl = partial_outflag

  call dte_integral_i_updn(maxgeoms, &
                           maxlayers, &
                           maxpartials, &
                           maxfine, &
                           max_user_levels, &
                           do_upwelling_lcl, &
                           do_dnwelling_lcl, &
                           do_thermset_lcl, &
                           do_deltam_scaling_lcl, &
                           do_partials_lcl, &
                           do_planpar_lcl, &
                           do_enhanced_ps_lcl, &
                           do_sources_up_lcl, &
                           do_sources_up_p_lcl, &
                           do_sources_dn_lcl, &
                           do_sources_dn_p_lcl, &
                           ngeoms, &
                           nlayers, &
                           nfinedivs, &
                           n_user_levels, &
                           user_levels, &
                           npartials, &
                           nfinedivs_p, &
                           partial_outindex, &
                           partial_outflag_lcl, &
                           partial_layeridx, &
                           bb_input, &
                           surfbb, &
                           user_emissivity, &
                           extinction, &
                           deltaus, &
                           omega, &
                           truncfac, &
                           mu1, &
                           losw_paths, &
                           losp_paths, &
                           xfine_up, &
                           wfine_up, &
                           hfine_up, &
                           xfine_dn, &
                           wfine_dn, &
                           hfine_dn, &
                           xfine_up_p, &
                           wfine_up_p, &
                           hfine_up_p, &
                           xfine_dn_p, &
                           wfine_dn_p, &
                           hfine_dn_p, &
                           intensity_dta_up, &
                           intensity_dts, &
                           intensity_dta_dn, &
                           cumsource_up, &
                           cumsource_dn, &
                           tcom1, &
                           lostrans_up, &
                           lostrans_up_p)

  ! Convert output arguments
  do_thermset = do_thermset_lcl

end subroutine fo_thermal_rtcalcs_i_m_dte_integral_i_updn_wrap

end module FO_THERMAL_RTCALCS_I_M_WRAP

module FO_THERMAL_RTCALCS_ILPS_M_WRAP

use iso_c_binding
use fo_thermal_rtcalcs_ilps_m

! This module was auto-generated 

implicit none

contains

subroutine fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_dn_wrap (maxgeoms, &
                                                                maxlayers, &
                                                                maxpartials, &
                                                                maxfine, &
                                                                max_user_levels, &
                                                                max_atmoswfs, &
                                                                do_thermset, &
                                                                do_deltam_scaling, &
                                                                do_partials, &
                                                                do_planpar, &
                                                                do_enhanced_ps, &
                                                                do_sources_dn, &
                                                                do_sources_dn_p, &
                                                                do_profilewfs, &
                                                                lvaryflags, &
                                                                lvarynums, &
                                                                ngeoms, &
                                                                nlayers, &
                                                                nfinedivs, &
                                                                n_user_levels, &
                                                                user_levels, &
                                                                npartials, &
                                                                nfinedivs_p, &
                                                                partial_outindex, &
                                                                partial_outflag, &
                                                                partial_layeridx, &
                                                                bb_input, &
                                                                extinction, &
                                                                deltaus, &
                                                                omega, &
                                                                truncfac, &
                                                                l_extinction, &
                                                                l_deltaus, &
                                                                l_omega, &
                                                                l_truncfac, &
                                                                mu1, &
                                                                losw_paths, &
                                                                losp_paths, &
                                                                xfine, &
                                                                wfine, &
                                                                hfine, &
                                                                xfine_p, &
                                                                wfine_p, &
                                                                hfine_p, &
                                                                intensity_dta_dn, &
                                                                lp_jacobians_dta_dn, &
                                                                tcom1, &
                                                                l_tcom1) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  integer(c_int), intent(in) :: max_atmoswfs
  logical(c_bool), intent(inout) :: do_thermset
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  logical(c_bool), intent(in) :: do_profilewfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: lvaryflags
  integer(c_int), dimension(MAXLAYERS), intent(in) :: lvarynums
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: bb_input
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_extinction
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltaus
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_truncfac
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_dn
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dta_dn
  real(c_double), dimension(MAXLAYERS, 2), intent(inout) :: tcom1
  real(c_double), dimension(MAXLAYERS, 2, MAX_ATMOSWFS), intent(inout) :: l_tcom1

  ! Local variables
  logical(kind=4) :: do_thermset_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4) :: do_profilewfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: lvaryflags_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_thermset_lcl = do_thermset
  do_deltam_scaling_lcl = do_deltam_scaling
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  do_profilewfs_lcl = do_profilewfs
  lvaryflags_lcl = lvaryflags
  partial_outflag_lcl = partial_outflag

  call dte_integral_ilps_dn(maxgeoms, &
                            maxlayers, &
                            maxpartials, &
                            maxfine, &
                            max_user_levels, &
                            max_atmoswfs, &
                            do_thermset_lcl, &
                            do_deltam_scaling_lcl, &
                            do_partials_lcl, &
                            do_planpar_lcl, &
                            do_enhanced_ps_lcl, &
                            do_sources_dn_lcl, &
                            do_sources_dn_p_lcl, &
                            do_profilewfs_lcl, &
                            lvaryflags_lcl, &
                            lvarynums, &
                            ngeoms, &
                            nlayers, &
                            nfinedivs, &
                            n_user_levels, &
                            user_levels, &
                            npartials, &
                            nfinedivs_p, &
                            partial_outindex, &
                            partial_outflag_lcl, &
                            partial_layeridx, &
                            bb_input, &
                            extinction, &
                            deltaus, &
                            omega, &
                            truncfac, &
                            l_extinction, &
                            l_deltaus, &
                            l_omega, &
                            l_truncfac, &
                            mu1, &
                            losw_paths, &
                            losp_paths, &
                            xfine, &
                            wfine, &
                            hfine, &
                            xfine_p, &
                            wfine_p, &
                            hfine_p, &
                            intensity_dta_dn, &
                            lp_jacobians_dta_dn, &
                            tcom1, &
                            l_tcom1)

  ! Convert output arguments
  do_thermset = do_thermset_lcl

end subroutine fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_dn_wrap

subroutine fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_up_wrap (maxgeoms, &
                                                                maxlayers, &
                                                                maxpartials, &
                                                                maxfine, &
                                                                max_user_levels, &
                                                                max_atmoswfs, &
                                                                max_surfacewfs, &
                                                                do_thermset, &
                                                                do_deltam_scaling, &
                                                                do_partials, &
                                                                do_planpar, &
                                                                do_enhanced_ps, &
                                                                do_sources_up, &
                                                                do_sources_up_p, &
                                                                do_profilewfs, &
                                                                do_surfacewfs, &
                                                                lvaryflags, &
                                                                lvarynums, &
                                                                n_surfacewfs, &
                                                                ngeoms, &
                                                                nlayers, &
                                                                nfinedivs, &
                                                                n_user_levels, &
                                                                user_levels, &
                                                                npartials, &
                                                                nfinedivs_p, &
                                                                partial_outindex, &
                                                                partial_outflag, &
                                                                partial_layeridx, &
                                                                bb_input, &
                                                                surfbb, &
                                                                user_emissivity, &
                                                                ls_user_emissivity, &
                                                                extinction, &
                                                                deltaus, &
                                                                omega, &
                                                                truncfac, &
                                                                l_extinction, &
                                                                l_deltaus, &
                                                                l_omega, &
                                                                l_truncfac, &
                                                                mu1, &
                                                                losw_paths, &
                                                                losp_paths, &
                                                                xfine, &
                                                                wfine, &
                                                                hfine, &
                                                                xfine_p, &
                                                                wfine_p, &
                                                                hfine_p, &
                                                                intensity_dta_up, &
                                                                intensity_dts, &
                                                                lp_jacobians_dta_up, &
                                                                lp_jacobians_dts_up, &
                                                                ls_jacobians_dts, &
                                                                tcom1, &
                                                                l_tcom1, &
                                                                lostrans_up, &
                                                                lostrans_up_p, &
                                                                l_lostrans_up, &
                                                                l_lostrans_up_p) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  integer(c_int), intent(in) :: max_atmoswfs
  integer(c_int), intent(in) :: max_surfacewfs
  logical(c_bool), intent(inout) :: do_thermset
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  logical(c_bool), intent(in) :: do_profilewfs
  logical(c_bool), intent(in) :: do_surfacewfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: lvaryflags
  integer(c_int), dimension(MAXLAYERS), intent(in) :: lvarynums
  integer(c_int), intent(in) :: n_surfacewfs
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: bb_input
  real(c_double), intent(in) :: surfbb
  real(c_double), dimension(MAXGEOMS), intent(in) :: user_emissivity
  real(c_double), dimension(MAXGEOMS, MAX_SURFACEWFS), intent(in) :: ls_user_emissivity
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_extinction
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltaus
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_truncfac
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dts
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dta_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dts_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAX_SURFACEWFS), intent(out) :: ls_jacobians_dts
  real(c_double), dimension(MAXLAYERS, 2), intent(inout) :: tcom1
  real(c_double), dimension(MAXLAYERS, 2, MAX_ATMOSWFS), intent(inout) :: l_tcom1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(out) :: lostrans_up
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: lostrans_up_p
  real(c_double), dimension(MAXLAYERS, MAXGEOMS, MAX_ATMOSWFS), intent(out) :: l_lostrans_up
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS, MAX_ATMOSWFS), intent(out) :: l_lostrans_up_p

  ! Local variables
  logical(kind=4) :: do_thermset_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4) :: do_profilewfs_lcl
  logical(kind=4) :: do_surfacewfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: lvaryflags_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_thermset_lcl = do_thermset
  do_deltam_scaling_lcl = do_deltam_scaling
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  do_profilewfs_lcl = do_profilewfs
  do_surfacewfs_lcl = do_surfacewfs
  lvaryflags_lcl = lvaryflags
  partial_outflag_lcl = partial_outflag

  call dte_integral_ilps_up(maxgeoms, &
                            maxlayers, &
                            maxpartials, &
                            maxfine, &
                            max_user_levels, &
                            max_atmoswfs, &
                            max_surfacewfs, &
                            do_thermset_lcl, &
                            do_deltam_scaling_lcl, &
                            do_partials_lcl, &
                            do_planpar_lcl, &
                            do_enhanced_ps_lcl, &
                            do_sources_up_lcl, &
                            do_sources_up_p_lcl, &
                            do_profilewfs_lcl, &
                            do_surfacewfs_lcl, &
                            lvaryflags_lcl, &
                            lvarynums, &
                            n_surfacewfs, &
                            ngeoms, &
                            nlayers, &
                            nfinedivs, &
                            n_user_levels, &
                            user_levels, &
                            npartials, &
                            nfinedivs_p, &
                            partial_outindex, &
                            partial_outflag_lcl, &
                            partial_layeridx, &
                            bb_input, &
                            surfbb, &
                            user_emissivity, &
                            ls_user_emissivity, &
                            extinction, &
                            deltaus, &
                            omega, &
                            truncfac, &
                            l_extinction, &
                            l_deltaus, &
                            l_omega, &
                            l_truncfac, &
                            mu1, &
                            losw_paths, &
                            losp_paths, &
                            xfine, &
                            wfine, &
                            hfine, &
                            xfine_p, &
                            wfine_p, &
                            hfine_p, &
                            intensity_dta_up, &
                            intensity_dts, &
                            lp_jacobians_dta_up, &
                            lp_jacobians_dts_up, &
                            ls_jacobians_dts, &
                            tcom1, &
                            l_tcom1, &
                            lostrans_up, &
                            lostrans_up_p, &
                            l_lostrans_up, &
                            l_lostrans_up_p)

  ! Convert output arguments
  do_thermset = do_thermset_lcl

end subroutine fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_up_wrap

subroutine fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_updn_wrap (maxgeoms, &
                                                                  maxlayers, &
                                                                  maxpartials, &
                                                                  maxfine, &
                                                                  max_user_levels, &
                                                                  max_atmoswfs, &
                                                                  max_surfacewfs, &
                                                                  do_upwelling, &
                                                                  do_dnwelling, &
                                                                  do_thermset, &
                                                                  do_deltam_scaling, &
                                                                  do_partials, &
                                                                  do_planpar, &
                                                                  do_enhanced_ps, &
                                                                  do_sources_up, &
                                                                  do_sources_up_p, &
                                                                  do_sources_dn, &
                                                                  do_sources_dn_p, &
                                                                  do_profilewfs, &
                                                                  do_surfacewfs, &
                                                                  lvaryflags, &
                                                                  lvarynums, &
                                                                  n_surfacewfs, &
                                                                  ngeoms, &
                                                                  nlayers, &
                                                                  nfinedivs, &
                                                                  n_user_levels, &
                                                                  user_levels, &
                                                                  npartials, &
                                                                  nfinedivs_p, &
                                                                  partial_outindex, &
                                                                  partial_outflag, &
                                                                  partial_layeridx, &
                                                                  bb_input, &
                                                                  surfbb, &
                                                                  user_emissivity, &
                                                                  ls_user_emissivity, &
                                                                  extinction, &
                                                                  deltaus, &
                                                                  omega, &
                                                                  truncfac, &
                                                                  l_extinction, &
                                                                  l_deltaus, &
                                                                  l_omega, &
                                                                  l_truncfac, &
                                                                  mu1, &
                                                                  losw_paths, &
                                                                  losp_paths, &
                                                                  xfine_up, &
                                                                  wfine_up, &
                                                                  hfine_up, &
                                                                  xfine_dn, &
                                                                  wfine_dn, &
                                                                  hfine_dn, &
                                                                  xfine_up_p, &
                                                                  wfine_up_p, &
                                                                  hfine_up_p, &
                                                                  xfine_dn_p, &
                                                                  wfine_dn_p, &
                                                                  hfine_dn_p, &
                                                                  intensity_dta_up, &
                                                                  intensity_dts, &
                                                                  intensity_dta_dn, &
                                                                  lp_jacobians_dta_up, &
                                                                  lp_jacobians_dts_up, &
                                                                  ls_jacobians_dts, &
                                                                  lp_jacobians_dta_dn, &
                                                                  tcom1, &
                                                                  l_tcom1, &
                                                                  lostrans_up, &
                                                                  lostrans_up_p, &
                                                                  l_lostrans_up, &
                                                                  l_lostrans_up_p) bind(C)

  ! Arguments
  integer(c_int), intent(in) :: maxgeoms
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxpartials
  integer(c_int), intent(in) :: maxfine
  integer(c_int), intent(in) :: max_user_levels
  integer(c_int), intent(in) :: max_atmoswfs
  integer(c_int), intent(in) :: max_surfacewfs
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(inout) :: do_thermset
  logical(c_bool), intent(in) :: do_deltam_scaling
  logical(c_bool), intent(in) :: do_partials
  logical(c_bool), intent(in) :: do_planpar
  logical(c_bool), intent(in) :: do_enhanced_ps
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_up
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_up_p
  logical(c_bool), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: do_sources_dn
  logical(c_bool), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: do_sources_dn_p
  logical(c_bool), intent(in) :: do_profilewfs
  logical(c_bool), intent(in) :: do_surfacewfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: lvaryflags
  integer(c_int), dimension(MAXLAYERS), intent(in) :: lvarynums
  integer(c_int), intent(in) :: n_surfacewfs
  integer(c_int), intent(in) :: ngeoms
  integer(c_int), intent(in) :: nlayers
  integer(c_int), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: nfinedivs
  integer(c_int), intent(in) :: n_user_levels
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: user_levels
  integer(c_int), intent(in) :: npartials
  integer(c_int), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: nfinedivs_p
  integer(c_int), dimension(MAX_USER_LEVELS), intent(in) :: partial_outindex
  logical(c_bool), dimension(MAX_USER_LEVELS), intent(in) :: partial_outflag
  integer(c_int), dimension(MAXPARTIALS), intent(in) :: partial_layeridx
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: bb_input
  real(c_double), intent(in) :: surfbb
  real(c_double), dimension(MAXGEOMS), intent(in) :: user_emissivity
  real(c_double), dimension(MAXGEOMS, MAX_SURFACEWFS), intent(in) :: ls_user_emissivity
  real(c_double), dimension(MAXLAYERS), intent(in) :: extinction
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltaus
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega
  real(c_double), dimension(MAXLAYERS), intent(in) :: truncfac
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_extinction
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltaus
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_truncfac
  real(c_double), dimension(MAXGEOMS), intent(in) :: mu1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(in) :: losw_paths
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(in) :: losp_paths
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine_up
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn
  real(c_double), dimension(MAXLAYERS, MAXFINE, MAXGEOMS), intent(in) :: hfine_dn
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_up_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: xfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: wfine_dn_p
  real(c_double), dimension(MAXPARTIALS, MAXFINE, MAXGEOMS), intent(in) :: hfine_dn_p
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dts
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS), intent(out) :: intensity_dta_dn
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dta_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAXLAYERS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dts_up
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAX_SURFACEWFS), intent(out) :: ls_jacobians_dts
  real(c_double), dimension(MAX_USER_LEVELS, MAXGEOMS, MAX_ATMOSWFS), intent(out) :: lp_jacobians_dta_dn
  real(c_double), dimension(MAXLAYERS, 2), intent(inout) :: tcom1
  real(c_double), dimension(MAXLAYERS, 2, MAX_ATMOSWFS), intent(inout) :: l_tcom1
  real(c_double), dimension(MAXLAYERS, MAXGEOMS), intent(out) :: lostrans_up
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS), intent(out) :: lostrans_up_p
  real(c_double), dimension(MAXLAYERS, MAXGEOMS, MAX_ATMOSWFS), intent(out) :: l_lostrans_up
  real(c_double), dimension(MAXPARTIALS, MAXGEOMS, MAX_ATMOSWFS), intent(out) :: l_lostrans_up_p

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_dnwelling_lcl
  logical(kind=4) :: do_thermset_lcl
  logical(kind=4) :: do_deltam_scaling_lcl
  logical(kind=4) :: do_partials_lcl
  logical(kind=4) :: do_planpar_lcl
  logical(kind=4) :: do_enhanced_ps_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_up_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_up_p_lcl
  logical(kind=4), dimension(MAXLAYERS, MAXGEOMS) :: do_sources_dn_lcl
  logical(kind=4), dimension(MAXPARTIALS, MAXGEOMS) :: do_sources_dn_p_lcl
  logical(kind=4) :: do_profilewfs_lcl
  logical(kind=4) :: do_surfacewfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: lvaryflags_lcl
  logical(kind=4), dimension(MAX_USER_LEVELS) :: partial_outflag_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_thermset_lcl = do_thermset
  do_deltam_scaling_lcl = do_deltam_scaling
  do_partials_lcl = do_partials
  do_planpar_lcl = do_planpar
  do_enhanced_ps_lcl = do_enhanced_ps
  do_sources_up_lcl = do_sources_up
  do_sources_up_p_lcl = do_sources_up_p
  do_sources_dn_lcl = do_sources_dn
  do_sources_dn_p_lcl = do_sources_dn_p
  do_profilewfs_lcl = do_profilewfs
  do_surfacewfs_lcl = do_surfacewfs
  lvaryflags_lcl = lvaryflags
  partial_outflag_lcl = partial_outflag

  call dte_integral_ilps_updn(maxgeoms, &
                              maxlayers, &
                              maxpartials, &
                              maxfine, &
                              max_user_levels, &
                              max_atmoswfs, &
                              max_surfacewfs, &
                              do_upwelling_lcl, &
                              do_dnwelling_lcl, &
                              do_thermset_lcl, &
                              do_deltam_scaling_lcl, &
                              do_partials_lcl, &
                              do_planpar_lcl, &
                              do_enhanced_ps_lcl, &
                              do_sources_up_lcl, &
                              do_sources_up_p_lcl, &
                              do_sources_dn_lcl, &
                              do_sources_dn_p_lcl, &
                              do_profilewfs_lcl, &
                              do_surfacewfs_lcl, &
                              lvaryflags_lcl, &
                              lvarynums, &
                              n_surfacewfs, &
                              ngeoms, &
                              nlayers, &
                              nfinedivs, &
                              n_user_levels, &
                              user_levels, &
                              npartials, &
                              nfinedivs_p, &
                              partial_outindex, &
                              partial_outflag_lcl, &
                              partial_layeridx, &
                              bb_input, &
                              surfbb, &
                              user_emissivity, &
                              ls_user_emissivity, &
                              extinction, &
                              deltaus, &
                              omega, &
                              truncfac, &
                              l_extinction, &
                              l_deltaus, &
                              l_omega, &
                              l_truncfac, &
                              mu1, &
                              losw_paths, &
                              losp_paths, &
                              xfine_up, &
                              wfine_up, &
                              hfine_up, &
                              xfine_dn, &
                              wfine_dn, &
                              hfine_dn, &
                              xfine_up_p, &
                              wfine_up_p, &
                              hfine_up_p, &
                              xfine_dn_p, &
                              wfine_dn_p, &
                              hfine_dn_p, &
                              intensity_dta_up, &
                              intensity_dts, &
                              intensity_dta_dn, &
                              lp_jacobians_dta_up, &
                              lp_jacobians_dts_up, &
                              ls_jacobians_dts, &
                              lp_jacobians_dta_dn, &
                              tcom1, &
                              l_tcom1, &
                              lostrans_up, &
                              lostrans_up_p, &
                              l_lostrans_up, &
                              l_lostrans_up_p)

  ! Convert output arguments
  do_thermset = do_thermset_lcl

end subroutine fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_updn_wrap

end module FO_THERMAL_RTCALCS_ILPS_M_WRAP
