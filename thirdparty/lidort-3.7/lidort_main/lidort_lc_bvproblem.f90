! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.7 F90                               #
! #  Release Date :   June 2014                             #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! #       Surface-leaving, BRDF Albedo-scaling     (3.7)    # 
! #       Taylor series, BBF Jacobians, ThreadSafe (3.7)    #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #  For linearizations involving atmospheric parameters        #
! #            LC_BVP_SOLUTION_MASTER                           #
! #            LC_BVP_COLUMN_SETUP                              #
! #            LC_BVP_SURFACE_SETUP                             #
! #            LC_BEAMSOLUTION_NEQK                             #
! #                                                             #
! #  For linearizations (atmospheric) using telescoped BVP      #
! #            LC_BVPTEL_SOLUTION_MASTER                        #
! #            LC_BVPTEL_COLUMN_SETUP                           #
! #                                                             #
! ###############################################################

module lidort_lc_bvproblem

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  LIDORT Use dependencies

   USE lidort_aux,      only : DGBTRS, DGETRS
   USE lidort_Taylor_m, only : TAYLOR_SERIES_L_1

!private
public

contains

SUBROUTINE LC_BVP_SOLUTION_MASTER &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                      & ! Input, FLags
             DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER_COMPONENT, IBEAM,  & ! Input, Flags and order
             NLAYERS,NTOTAL,N_SUB,N_SUP,NSTREAMS,NSTREAMS_2,N_WEIGHTFUNCS, & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
             SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
             BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,         & ! Input, BVP Bandmat
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous solution 
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + Thermal
             NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,         & ! Output - Linearized Constants + PI
             STATUS, MESSAGE, TRACE )                                        ! Output - Exception handling

!  Linearization of the Boundary Problem Solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)   :: TAYLOR_ORDER

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER_COMPONENT, IBEAM

!  Number of layers, BVP control

      INTEGER  , intent(in)  ::  NLAYERS, NTOTAL, N_SUB, N_SUP

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  Factor, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Direct beam solutions

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER  , intent(in)  ::  LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(in)  ::  MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Homogeneous and thermal solution variables
!  ------------------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized Beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local help variables 

      INTEGER    ::  I, I1, Q, N, AA

!  initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Linearization of the regular BVP case
!  =====================================

!  Set up the column vectors for Column linearizations
!  ---------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'

      CALL LC_BVP_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                      & ! Input, FLags
             DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER_COMPONENT, IBEAM,  & ! Input, Flags and order
             NLAYERS, NTOTAL, NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS,         & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
             SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
             L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS  ( 'n', NTOTAL, N_SUB, N_SUP, N_WEIGHTFUNCS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT,  COL2_WF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_WEIGHTFUNCS for DGBTRS call in LC_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(I,N,Q) = COL2_WF(LCONMASK(I,N),Q)
              PCON(I,N,Q) = COL2_WF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_WEIGHTFUNCS for 1-layer: DGETRS call in LC_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_WEIGHTFUNCS
            NCON(I,N,Q) = SCOL2_WF(I,Q)
            PCON(I,N,Q) = SCOL2_WF(I1,Q)
          ENDDO
        ENDDO

      ENDIF

!  linearized BVP results
!  ======================

!  Associated quantities

      DO N = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          DO AA = 1, NSTREAMS
             DO Q = 1, N_WEIGHTFUNCS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  finish

      RETURN
END SUBROUTINE LC_BVP_SOLUTION_MASTER

!

SUBROUTINE LC_BVP_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                      & ! Input, FLags
             DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IPARTIC,          & ! Input, Flags and order
             NLAYERS, NTOTAL, NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS,         & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
             SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,    & ! Input, Surface Stuff
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
             L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  Linearized column vector setup (bulk property linearization)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS, NTOTAL

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER, IPARTIC

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  Factor, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Direct beam solutions

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Homogeneous
!  -----------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized Beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BBV Problem

      REAL(fpk), intent(out) :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  help variables

      INTEGER    :: Q, AA, N, N1, I, I1, CM, C0, K, NV, M, IB
      REAL(fpk)  :: CPOS, CNEG, L_HOM, L_PARTIC
      REAL(fpk)  :: L_BEAM, L_HOMD, L_HOMU, FAC, FAC3
      LOGICAL    :: MODIFIED_BOUNDARY

!  Local linearized reflectance arrays

      REAL(fpk)  :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  initialise
!  ----------

      M  = FOURIER
      IB = IPARTIC

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, MAX_ATMOSWFS
          COL2_WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

      DO NV = 1, NLAYERS
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, NSTREAMS_2
            DO Q = 1, N_WEIGHTFUNCS
              L_WUPPER(I,NV,Q) = L_T_WUPPER(I,NV,Q)
              L_WLOWER(I,NV,Q) = L_T_WLOWER(I,NV,Q)
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS_2
            DO Q = 1, N_WEIGHTFUNCS
              L_WUPPER(I,NV,Q) = ZERO
              L_WLOWER(I,NV,Q) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Top of first layer (TOA), UPPER boundary condition
!  --------------------------------------------------

      N = 1

!  Get the linearized beam solution for the first layer
!   Adds to the thermal solution if already present

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LC_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS, M, IB, N,        & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
      ENDIF

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO I = 1, NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          L_PARTIC = - L_WUPPER(I,N,Q)
          L_HOM    = ZERO
          DO AA = 1, NSTREAMS
            CPOS = L_XPOS(I,AA,N,Q)
            CNEG =   T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
            L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
          ENDDO
          COL2_WF(I,Q) = L_PARTIC - L_HOM
        ENDDO
      ENDDO
 
!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer
!   Adds to the thermal solution if already present

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LC_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS, M, IB, N1,       & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
        ENDIF

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTREAMS_2
          CM = C0 + I
          DO Q = 1, N_WEIGHTFUNCS
            L_HOMD = ZERO
            L_HOMU = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(I,AA,N1,Q)
              CNEG =   T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + &
                     L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
              L_HOMU = L_HOMU + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
              CNEG = L_XNEG(I,AA,N,Q)
              CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + &
                     L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I,AA,N)
              L_HOMD = L_HOMD + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            L_PARTIC      = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            L_HOM         = L_HOMU - L_HOMD
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
          ENDDO
        ENDDO

!  End layer

      ENDDO

!  LOWER layer 
!  -----------

      N = NLAYERS
      MODIFIED_BOUNDARY = .TRUE.

!  get the linearized downward-reflected term

      CALL LC_BVP_SURFACE_SETUP                            &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BOUNDARY, FOURIER,           & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_WEIGHTFUNCS, & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output 

!  Compute the solution

      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      DO I = 1, NSTREAMS
        CM = C0 + I
        I1 = I + NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)
          L_HOM    = ZERO
          DO AA = 1, NSTREAMS
            CPOS =   T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
                   L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
            CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
            CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
            L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
          ENDDO
          COL2_WF(CM,Q) = - L_PARTIC - L_HOM
        ENDDO
      ENDDO

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          FAC = - DIRECT_BEAM(I,IPARTIC)
          DO Q = 1, N_WEIGHTFUNCS
            L_BEAM = ZERO
            DO K = 1, NLAYERS
              FAC3 = FAC * DELTAU_SLANT(N,K,IPARTIC)
              L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K) * FAC3
            ENDDO
            COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
          ENDDO
        ENDDO
      ENDIF

!  Single layer

      IF ( NLAYERS .eq. 1 ) THEN
        DO I = 1, NTOTAL
          DO Q = 1, N_WEIGHTFUNCS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  finish

      RETURN
END SUBROUTINE LC_BVP_COLUMN_SETUP

!

SUBROUTINE LC_BVP_SURFACE_SETUP                            &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, FOURIER_COMPONENT,     & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_WEIGHTFUNCS, & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output 

!  Linearized surface reflectance terms

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  input arguments
!  ---------------

!  BRDF flag

      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE

!  Number of streams and layers

      INTEGER  , intent(in)  ::  NSTREAMS, NLAYERS

!  Number of weighting functions

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier component

      INTEGER  , intent(in)  ::  FOURIER_COMPONENT

!  Flag for type of boundary condition

      LOGICAL  , intent(in)  ::  MODIFIED_BCL4

!  overall surface flag and surface factor, albedo

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

      REAL(fpk), intent(out) :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      REAL(fpk)  :: PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HSP_U, HSM_U, REFL_P, REFL_B, REFL_M
      REAL(fpk)  :: FACTOR
      INTEGER    :: AA, J, Q, N, I, M

!  Initial section
!  ---------------

!  Always zero the result to start

      DO I = 1, NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          R2_L_PARTIC(I,Q)    = ZERO
          DO AA = 1, NSTREAMS
            R2_L_HOMP(I,AA,Q) = ZERO
            R2_L_HOMM(I,AA,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  FOurier component

      M = FOURIER_COMPONENT

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

      N = NLAYERS

!  Particular integral parts

      DO J = 1, NSTREAMS
        DO Q = 1, N_WEIGHTFUNCS
          PV_W(J,Q) = L_WLOWER(J,N,Q) * QUAD_STRMWTS(J)
        ENDDO
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BCL4 ) THEN
        DO J = 1, NSTREAMS
          DO AA = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              HSP_U = L_XPOS(J,AA,N,Q) *   T_DELT_EIGEN(AA,N) + &
                        XPOS(J,AA,N)   * L_T_DELT_EIGEN(AA,N,Q)
              HSM_U = L_XNEG(J,AA,N,Q)
              HV_P(J,AA,Q) = QUAD_STRMWTS(J)*HSP_U
              HV_M(J,AA,Q) = QUAD_STRMWTS(J)*HSM_U
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Integrated Downward reflection (Calculation, Lambertian case)
!  -------------------------------------------------------------

     IF ( .not. DO_BRDF_SURFACE ) THEN

!  amplitude

        FACTOR = SURFACE_FACTOR * ALBEDO

!  Only a solution for the isotropic part.

        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO Q = 1, N_WEIGHTFUNCS

!  Particular solution

            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + PV_W(J,Q)
            ENDDO
            REFL_B = REFL_B * FACTOR
            DO I = 1, NSTREAMS
              R2_L_PARTIC(I,Q) = R2_L_PARTIC(I,Q) + REFL_B
            ENDDO

!  Homogeneous solutions (only for modified BC)

            IF ( MODIFIED_BCL4 ) THEN
              DO AA = 1, NSTREAMS
                REFL_P = ZERO
                REFL_M = ZERO
                DO J = 1, NSTREAMS
                  REFL_P = REFL_P + HV_P(J,AA,Q)
                  REFL_M = REFL_M + HV_M(J,AA,Q)
                ENDDO
                REFL_P = REFL_P * FACTOR
                REFL_M = REFL_M * FACTOR
                DO I = 1, NSTREAMS
                  R2_L_HOMP(I,AA,Q) = R2_L_HOMP(I,AA,Q) + REFL_P
                  R2_L_HOMM(I,AA,Q) = R2_L_HOMM(I,AA,Q) + REFL_M
                ENDDO
              ENDDO
            ENDIF

!  end parameter loop

          ENDDO
        ENDIF

!  Integrated Downward reflection (Calculation, Bidirectional case)
!  ----------------------------------------------------------------

      ELSE

        DO I = 1, NSTREAMS
          DO Q = 1, N_WEIGHTFUNCS

!  particular solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

            REFL_B = ZERO
            DO J = 1, NSTREAMS
!              REFL_B = REFL_B + PV_W(J,Q) * BRDF_F(M,J,I)
              REFL_B = REFL_B + PV_W(J,Q) * BRDF_F(M,I,J)
            ENDDO
            REFL_B = REFL_B * SURFACE_FACTOR
            R2_L_PARTIC(I,Q) = R2_L_PARTIC(I,Q) + REFL_B

!  homogeneous solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

            IF ( MODIFIED_BCL4 ) THEN
              DO AA = 1, NSTREAMS
                REFL_P = ZERO
                REFL_M = ZERO
                DO J = 1, NSTREAMS
!                  REFL_P = REFL_P + HV_P(J,AA,Q) * BRDF_F(M,J,I)
!                  REFL_M = REFL_M + HV_M(J,AA,Q) * BRDF_F(M,J,I)
                  REFL_P = REFL_P + HV_P(J,AA,Q) * BRDF_F(M,I,J)
                  REFL_M = REFL_M + HV_M(J,AA,Q) * BRDF_F(M,I,J)
                ENDDO
                REFL_P = REFL_P * SURFACE_FACTOR
                REFL_M = REFL_M * SURFACE_FACTOR
                R2_L_HOMP(I,AA,Q) = R2_L_HOMP(I,AA,Q) + REFL_P
                R2_L_HOMM(I,AA,Q) = R2_L_HOMM(I,AA,Q) + REFL_M
              ENDDO
            ENDIF

!  end parameter and stream loops

          ENDDO
        ENDDO

!  End BRDF clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_BVP_SURFACE_SETUP

!

SUBROUTINE LC_BEAMSOLUTION_NEQK &
       ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flag and order
         NSTREAMS, NSTREAMS_2, N_WEIGHTFUNCS, M, IB, N,        & ! Input, Numbers
         DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
         INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
         T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
         GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
         ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
         L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N, due to column.
!   This is the bulk property linearization

!  Completely rewritten, 11 October 2013
!    * Final version with Taylor series expansions
!    * Removal of LC_GAMMA

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier number, beam index, layer index

      INTEGER  , intent(in)  ::  M, IB, N

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvalues, Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ===============

!  Local linearized Green's function multipliers

      REAL(fpk)  :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: AA, I, I1, Q
      REAL(fpk)  :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(fpk)  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT
      REAL(fpk)  :: L_WDEL, L_ZDEL, L_LAM, L_KEG, L_DELTA, L_AST, L_BST, L_MULT(MAX_ATMOSWFS)
      REAL(fpk)  :: LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)

!  Green's function solution
!  =========================

!  No linearized particular solution beyond the cutoff layer. ALSO -
!  Nothing if layer is inactive (Does not depend on solution saving)
!    Exit Immediately (Solutions have already been zeroed).

!      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N)) .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) RETURN
      IF (.NOT.DO_LAYER_SCATTERING(M,N) .OR. (N .GT.LAYER_PIS_CUTOFF(IB))) RETURN

!  Linearizations of optical depth integrations (Linearized Green function multipliers)
!  ------------------------------------------------------------------------------------

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

      DO AA = 1, NSTREAMS

         ZDEL  = T_DELT_EIGEN(AA,N)
         ZWDEL = ZDEL * WDEL
         AST   = CONST * ATERM_SAVE(AA,N) 
         BST   = CONST * BTERM_SAVE(AA,N) 

!  Downwelling, Make allowances for Taylor series

        IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
           EPS   = GAMMA_M(AA,N)
           DELTA = DELTAU_VERT(N)
           LAM   = AVERAGE_SECANT(N,IB)
           DO Q = 1, N_WEIGHTFUNCS
              L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
              L_KEG  = L_KEIGEN(AA,N,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KEG, L_LAM, WDEL, LAM, L_MULT(Q) )
           ENDDO
        ELSE
           MULT = ( ZDEL - WDEL ) / GAMMA_M(AA,N)
           DO Q = 1, N_WEIGHTFUNCS
              L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
              L_WDEL =  LC_T_DELT_MUBAR (N,IB,Q) ; L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
              L_MULT(Q) = ( ( L_ZDEL - L_WDEL ) - MULT * (L_LAM - L_KEG) ) / GAMMA_M(AA,N)
           ENDDO
        ENDIF
        DO Q = 1, N_WEIGHTFUNCS
           L_AST =  LC_INITIAL_TRANS(N,IB,Q)  + L_ATERM_SAVE(AA,N,Q)
           L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * L_AST + L_MULT(Q) * AST
        ENDDO

!  Upwelling

        MULT = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
        DO Q = 1, N_WEIGHTFUNCS
           L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
           L_WDEL =  LC_T_DELT_MUBAR (N,IB,Q) ; L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
           L_MULT(Q) = - ( L_ZDEL*WDEL + L_WDEL*ZDEL + MULT*(L_LAM+L_KEG) ) / GAMMA_P(AA,N)
        ENDDO
        DO Q = 1, N_WEIGHTFUNCS
           L_BST =  LC_INITIAL_TRANS(N,IB,Q)  + L_BTERM_SAVE(AA,N,Q)
           L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * L_BST + L_MULT(Q) * BST
        ENDDO

!  End discrete ordinate loop

      ENDDO

!  Set linearized form of particular integral at boundaries
!  --------------------------------------------------------

      DO Q = 1, N_WEIGHTFUNCS
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           S_P_U = ZERO ; S_P_L = ZERO
           S_M_U = ZERO ; S_M_L = ZERO
           DO AA = 1, NSTREAMS
              S_P_U = S_P_U + L_GFUNC_UP(AA,Q) * XPOS(I1,AA,N) + GFUNC_UP(AA,N) * L_XPOS(I1,AA,N,Q)
              S_M_U = S_M_U + L_GFUNC_UP(AA,Q) * XPOS(I,AA,N)  + GFUNC_UP(AA,N) * L_XPOS(I,AA,N,Q)
              S_P_L = S_P_L + L_GFUNC_DN(AA,Q) * XPOS(I,AA,N)  + GFUNC_DN(AA,N) * L_XPOS(I,AA,N,Q)
              S_M_L = S_M_L + L_GFUNC_DN(AA,Q) * XPOS(I1,AA,N) + GFUNC_DN(AA,N) * L_XPOS(I1,AA,N,Q)
           ENDDO
           LBSOL(I,Q,1)  = S_P_U ; LBSOL(I1,Q,1) = S_M_U
           LBSOL(I1,Q,2) = S_M_L ; LBSOL(I,Q,2)  = S_P_L
        ENDDO
      ENDDO

!  Add to existing solution

      DO Q = 1, N_WEIGHTFUNCS
         DO I = 1, NSTREAMS_2
            L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
            L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
         ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_BEAMSOLUTION_NEQK

!

SUBROUTINE LC_BVPTEL_SOLUTION_MASTER &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER,                    & ! Input, Flags, Order
             IBEAM, NLAYERS,  NSTREAMS, NSTREAMS_2, N_LAYER_WFS,            & ! Input, Numbers
             N_SUB, N_SUP, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,  & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,                  & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam Quantities
             BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,                        & ! Input, BVProblem
             T_DELT_DISORDS, L_T_DELT_DISORDS,  LCON_XVEC, MCON_XVEC,       & ! Input, Extinction, non active layers 
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                          & ! Input, Homogeneous + constants
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,  & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,          & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                      & ! Input, Linearized Homogeneous
             L_ATERM_SAVE, L_BTERM_SAVE,                                    & ! Input, Linearized Greens
             L_WUPPER, L_WLOWER, NCON, PCON, NCON_XVEC, PCON_XVEC,          & ! Output
             STATUS, MESSAGE, TRACE )                                         ! Output

!  Linearization of the Telescoped Boundary Problem Solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                              ZERO, LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_LAYER_WFS


!  BVProblem Band matrix control

      INTEGER  , intent(in)  ::  N_SUB, N_SUP

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2       (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDTELMAT2 (MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT    (MAXSTREAMS_2)

!  discrete ordinate factors (BVP telescoping, solutions saving)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local help variables 

      INTEGER    :: I, I1, K, K1, Q, M
      INTEGER    :: NS, N, N1, NAF, NAL, AA, C0, CM, CP
      REAL(fpk)  :: SHOM, L_HOM1, L_HOM2

!  initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Wf count

      M = FOURIER
      Q = 0

!  Set up linearized BVP, column vector  B' where AX = B'
!  ======================================================

!  Bulk: Compute the main column B' where AX = B'

      CALL LC_BVPTEL_COLUMN_SETUP &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                            & ! Input, Flags and order
             NLAYERS, NSTREAMS, NSTREAMS_2, N_LAYER_WFS, M, IBEAM,         & ! Input, Numbers
             ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,               & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,                 & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous solution
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous 
             L_ATERM_SAVE, L_BTERM_SAVE,                                   & ! Input, Linearized Greens 
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                      ! Output

!       if ( fourier_component .eq. 3 .and. ibeam.eq.4 ) then
!          do n = 1, N_BVTELMATRIX_SIZE
!             write(*,*)n,coltel2_wf(n,1),coltel2_wf(n,2)
!          enddo
!       endif

!  Solve linearized BVP, several active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS  ( 'n', N_BVTELMATRIX_SIZE, N_SUB, N_SUP, N_LAYER_WFS, &
                       BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, COLTEL2_WF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = 'DGBTRS call in L_BVPTEL_BEAMSOLUTION_MASTER'
         STATUS  = LIDORT_SERIOUS
         RETURN
        ENDIF

!  Set linearized integration constants, active layers

        C0 = - NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         C0 = C0 + NSTREAMS_2
         DO I = 1, NSTREAMS
          CM = C0 + I
          CP = CM + NSTREAMS
          DO Q = 1, N_LAYER_WFS
           NCON(I,N,Q) = COLTEL2_WF(CM,Q)
           PCON(I,N,Q) = COLTEL2_WF(CP,Q)
          ENDDO
         ENDDO
        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NSTREAMS_2, N_LAYER_WFS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in LC_BVPTEL_BEAMSOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set linearized integration constants for active layer

        N = ACTIVE_LAYERS(1)
        DO K = 1, NSTREAMS
         K1 = K + NSTREAMS 
         DO Q = 1, N_LAYER_WFS
           NCON(K,N,Q) = SCOL2_WF(K,Q)
           PCON(K,N,Q) = SCOL2_WF(K1,Q)
         ENDDO
        ENDDO

!  end clause for backsubstitution

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            DO Q = 1, N_LAYER_WFS
              NCON_XVEC(I,K,N,Q) = NCON(K,N,Q) * XPOS(I,K,N)
              PCON_XVEC(I,K,N,Q) = PCON(K,N,Q) * XNEG(I,K,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required, first layer is always active

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              L_HOM1 = NCON_XVEC(I1,AA,NAF,Q) + LCON(AA,NAF)*L_XPOS(I1,AA,NAF,Q)
              L_HOM2 = T_DELT_EIGEN(AA,NAF) * ( PCON_XVEC(I1,AA,NAF,Q)   + &
                         MCON(AA,NAF)*L_XNEG(I1,AA,NAF,Q) )              + &
                         L_T_DELT_EIGEN(AA,NAF,Q) * MCON_XVEC(I1,AA,NAF)
              SHOM = SHOM + L_HOM1 + L_HOM2
            ENDDO
            PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
            NCON(I,N1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards.
!   Additional linearizations if you are passing through the varying layer.

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            NCON(I,N,Q) = ZERO
            PCON(I,N,Q) = T_DELT_DISORDS(I,N1) * PCON(I,N1,Q) &
                      + L_T_DELT_DISORDS(I,N1,Q) * MCON(I,N1)
          ENDDO
        ENDDO
      ENDDO     

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!   -- PCON values are zero (no upwelling radiation)
!   -- NCON values propagated downwards from bottom of last active layer

!  layer immediately below Last active layer
!    .... Require linearized solutions at bottom of last active layer
!    .... Additional linearizations in the last active layer, always active

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN
        N1 = NAL + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              L_HOM2 = PCON_XVEC(I,AA,NAL,Q) + MCON(AA,NAL)*L_XNEG(I,AA,NAL,Q)
              L_HOM1 = T_DELT_EIGEN(AA,NAL) * &
                ( NCON_XVEC(I,AA,NAL,Q) + LCON(AA,NAL)*L_XPOS(I,AA,NAL,Q) ) + &
                  L_T_DELT_EIGEN(AA,NAL,Q) * LCON_XVEC(I,AA,NAL)
              SHOM = SHOM + L_HOM1 + L_HOM2
            ENDDO
            NCON(I,N1,Q) = L_WLOWER(I,NAL,Q) + SHOM
            PCON(I,N1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  other layers to bottom of medium: propagate downwards.
!   Additional variation if you are passing through the varying layer.

      DO N = NAL + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_LAYER_WFS
            PCON(I,N,Q) = ZERO
            NCON(I,N,Q) = T_DELT_DISORDS(I,N1) * NCON(I,N1,Q) &
                      + L_T_DELT_DISORDS(I,N1,Q) * LCON(I,N1)
          ENDDO
        ENDDO
      ENDDO

!  Associated quantities for inactive layers
!  -----------------------------------------

!  atmosphere layers with no scattering

      DO N = 1, NLAYERS
        IF ( N .LT. NAF .OR. N.GT.NAL ) THEN
          DO I = 1, NSTREAMS_2
            DO AA = 1, NSTREAMS
             DO Q = 1, N_LAYER_WFS
              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
             ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
 
!  finish

      RETURN
END SUBROUTINE LC_BVPTEL_SOLUTION_MASTER

!

SUBROUTINE LC_BVPTEL_COLUMN_SETUP &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NLAYERS, NSTREAMS, NSTREAMS_2, N_LAYER_WFS, M, IBEAM, & ! Input, Numbers
             ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,       & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                  & ! Input, Beam Quantities
             T_DELT_EIGEN, XPOS, XNEG, LCON, MCON,                         & ! Input, Homogeneous solution
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                     & ! Input, Linearized Homogeneous 
             L_ATERM_SAVE, L_BTERM_SAVE,                                   & ! Input, Linearized Greens 
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                      ! Output

!  Column setup for the linearized telescoped BVP

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                              MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_LAYER_WFS

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  M, IBEAM

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(fpk), intent(out) :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER     :: Q,AA,N,N1,NS,I,I1,CM,C0,NAF
      REAL(fpk)   :: CPOS, CNEG, L_HOM, L_BEAM, L_HOMD, L_HOMU

!  Try this safety-first zeroing

      DO I = 1, NSTREAMS_2
        DO Q = 1, N_LAYER_WFS
          DO N = 1, NLAYERS
            L_WUPPER(I,N,Q) = ZERO
            L_WLOWER(I,N,Q) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Get the linearized solutions for all active layers
!    Always need this, regardless of number of active layers

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        CALL LC_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                      & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N_LAYER_WFS, M, IBEAM, N,         & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,           & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,            & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,   & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS,   & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                   & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,     & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                      ! Output
      ENDDO

!  Go to special case for only 1 active layer

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  zero column vector

       DO I = 1, N_BVTELMATRIX_SIZE
         DO Q = 1, MAX_ATMOSWFS
           COLTEL2_WF(I,Q) = ZERO
         ENDDO
       ENDDO

!  top of first active layer, first boundary condition
!  ---------------------------------------------------

       NS = 1
       N = ACTIVE_LAYERS(NS)
       NAF = N

!  require homogeneous and beam solution linearizations

       DO I = 1, NSTREAMS
         DO Q = 1, N_LAYER_WFS
           L_BEAM = - L_WUPPER(I,N,Q)
           L_HOM  = ZERO
            DO AA = 1, NSTREAMS
             CPOS = L_XPOS(I,AA,N,Q)
             CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) +  &
                  L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
             L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
           ENDDO
           COLTEL2_WF(I,Q) = L_BEAM - L_HOM
         ENDDO
       ENDDO

!  Intermediate boundaries between active layers
!  ---------------------------------------------

       DO NS = 1, NLAYERS_TEL - 1

!  offsets

         N  = ACTIVE_LAYERS(NS)
         N1 = N + 1
         C0 = NS*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer N1

         DO I = 1, NSTREAMS_2
           CM = C0 + I
           DO Q = 1, N_LAYER_WFS
             L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
             L_HOMU = ZERO
             L_HOMD = ZERO
             DO AA = 1, NSTREAMS
               CPOS = L_XPOS(I,AA,N1,Q)
               CNEG =   T_DELT_EIGEN(AA,N1)   * L_XNEG(I,AA,N1,Q) + &
                      L_T_DELT_EIGEN(AA,N1,Q) *   XNEG(I,AA,N1)
               L_HOMU = L_HOMU + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
               CNEG = L_XNEG(I,AA,N,Q)
               CPOS =   T_DELT_EIGEN(AA,N)    * L_XPOS(I,AA,N,Q)  +  &
                      L_T_DELT_EIGEN(AA,N,Q)  *   XPOS(I,AA,N)
               L_HOMD = L_HOMD + LCON(AA,N)  * CPOS + MCON(AA,N)  * CNEG
             ENDDO
             L_HOM            = L_HOMU - L_HOMD
             COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
            ENDDO
         ENDDO

!  End loop over intermediate active layer boundaries

       ENDDO

!  Final boundary, bottom of lowest active layer
!  ---------------------------------------------

       NS = NLAYERS_TEL
       N  = ACTIVE_LAYERS(NS)      
       C0 = (NS-1)*NSTREAMS_2 + NSTREAMS

!  last active layer is varying

       DO I = 1, NSTREAMS
         CM = C0 + I
         I1 = I + NSTREAMS
         DO Q = 1, N_LAYER_WFS
           L_BEAM = L_WLOWER(I1,N,Q)
           L_HOM  = ZERO
           DO AA = 1, NSTREAMS
             CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + &
                  L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
             CNEG = L_XNEG(I1,AA,N,Q)
             L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
           ENDDO
           COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
         ENDDO
       ENDDO

!  Continuation point for the single-active-layer case

      ELSE

!  zero column vector
!   Probably not needed

       DO I = 1, NSTREAMS_2
         DO Q = 1, MAX_ATMOSWFS
           SCOL2_WF(I,Q) = ZERO
         ENDDO
       ENDDO

!  top of active layer
!  -------------------

       NS = 1
       N = ACTIVE_LAYERS(NS)

!  layer that is varying,
!    require homogeneous and beam solution linearizations

       DO I = 1, NSTREAMS
         DO Q = 1, N_LAYER_WFS
           L_BEAM = - L_WUPPER(I,N,Q)
           L_HOM  = ZERO
           DO AA = 1, NSTREAMS
            CPOS = L_XPOS(I,AA,N,Q)
             CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + &
                  L_T_DELT_EIGEN(AA,N,Q) *   XNEG(I,AA,N)
             L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
           ENDDO
           SCOL2_WF(I,Q) = L_BEAM - L_HOM
         ENDDO
       ENDDO

!  Bottom of active layer
!  ----------------------

!  active layer is varying layer

       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         DO Q = 1, N_LAYER_WFS
           L_BEAM = L_WLOWER(I1,N,Q)
           L_HOM  = ZERO
           DO AA = 1, NSTREAMS
             CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) +  &
                  L_T_DELT_EIGEN(AA,N,Q) *   XPOS(I1,AA,N)
             CNEG = L_XNEG(I1,AA,N,Q)
             L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
           ENDDO
           SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
         ENDDO
       ENDDO

!  End layer clause

      ENDIF

!  finish

      RETURN
END SUBROUTINE LC_BVPTEL_COLUMN_SETUP

!  End

end module lidort_lc_bvproblem
