
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ##########################################################
! #                                                        #
! #      PUBLIC Subroutines in this Module                 #
! #                                                        #
! #         SURFACEWF_BVP_SOLUTION,      LS BVP Master     #
! #         SURFACEWF_BVPTEL_SOLUTION,   LS BVPTel Master  #
! #         SURFACEWF_POSTPROCESS_MASTER LS PP Master      #
! #                                                        #
! #      PRIVATE                                           #
! #                                                        #
! #         BOA_SURFACEWF                                  #
! #                                                        #
! #         UPUSER_SURFACEWF, calling...                   #
! #            LS_WHOLELAYER_STERM_UP                      #
! #            LS_PARTLAYER_STERM_UP                       #
! #                                                        #
! #         DNUSER_SURFACEWF, calling...                   #
! #            LS_WHOLELAYER_STERM_DN                      #
! #            LS_PARTLAYER_STERM_DN                       #
! #                                                        #
! #         VLIDORT_LS_INTEGRATED_OUTPUT, Calling...       #
! #            QUADSURFACEWF_LEVEL_UP                      #
! #            QUADSURFACEWF_LEVEL_DN                      #
! #            QUADSURFACEWF_OFFGRID_UP                    #
! #            QUADSURFACEWF_OFFGRID_DN                    #
! #                                                        #
! ##########################################################

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO constructions with IF blocks
!     * Generalized surface treatment with telescoped BVP
!     * Rearrangement of subroutines

!  Version 2.8.1. April-May 2019.
!   4/9/19. Introduce Water-leaving control
!   4/9/19. Introduce variability for adjusted water-leaving transmittance (handle only)

!  1/31/21. Version 2.8.3. Some changes.
!    -- BRDF arrays are defined locally for each Fourier
!    -- Introduce DO_MSSTS flag for linearized MSST output
!    -- Drop LOCAL_UM_START; all UM loops start with 1
!    -- Use the post processing mask system, don't need separate OBSGEOM calculation

      MODULE vlidort_ls_wfsurface_m

      PRIVATE :: BOA_SURFACEWF, UPUSER_SURFACEWF, DNUSER_SURFACEWF, VLIDORT_LS_INTEGRATED_OUTPUT,              &
                 LS_WHOLELAYER_STERM_UP, LS_PARTLAYER_STERM_UP, LS_WHOLELAYER_STERM_DN, LS_PARTLAYER_STERM_DN, &
                 QUADSURFACEWF_LEVEL_UP, QUADSURFACEWF_OFFGRID_UP, QUADSURFACEWF_LEVEL_DN, QUADSURFACEWF_OFFGRID_DN

      PUBLIC  :: SURFACEWF_BVP_SOLUTION,    &
                 SURFACEWF_BVPTEL_SOLUTION, &
                 SURFACEWF_POSTPROCESS_MASTER

      CONTAINS

      SUBROUTINE SURFACEWF_BVP_SOLUTION ( FOURIER,&
              DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, DO_LAMBERTIAN_SURFACE, & ! Input
              DO_WATER_LEAVING, IBEAM, NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS, & ! Input
              NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,         & ! Input
              MUELLER_INDEX, K_REAL, K_COMPLEX, QUAD_STRMWTS,                     & ! Input
              SURFACE_FACTOR, ATMOS_ATTN, LS_TRANS_ATMOS_FINAL, SL_QUADTERM,      & ! Input
              SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0,                      & ! Input
              T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WLOWER,                         & ! Input
              BANDMAT2, IPIVOT, SMAT2, SIPIVOT, LCON, MCON,                       & ! Input
              NCON_SWF, PCON_SWF, STATUS, MESSAGE, TRACE )                          ! Output

!  Regular-BVP solution for the surface weighting functions
!   4/9/19. Introduce variability for adjusted water-leaving transmittance (handle only)

!  1/31/21. Version 2.8.3. BRDF arrays LS_BRDF_F, LS_BRDF_F_0 are defined locally for each Fourier

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAX_SURFACEWFS,         &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ, &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO

      USE lapack_tools_m, only : DGBTRS, DGETRS

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::          DO_WATER_LEAVING   ! New for 2.8.1, 4/9/19

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  other derived numbers

      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

!  Bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Surface stuff
!    -- 1/31/21. Version 2.8.3.  LS_BRDF_F, LS_BRDF_F_0 defined locally, dropo the MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: SURFBB
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F     ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F_0   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )

!  4/9/19. Adjusted water-leaving inputs

      DOUBLE PRECISION, INTENT (IN) :: SL_QUADTERM(MAXSTREAMS,MAXBEAMS,MAXSTOKES)
      DOUBLE PRECISION, INTENT (IN) :: LS_TRANS_ATMOS_FINAL(MAXBEAMS,MAX_SURFACEWFS)

!  RTE solution variables

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  BVP Matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  Attenuation to the surface

      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )

!  outputs
!  -------

!  Linearized Solution constants of integration

      DOUBLE PRECISION, intent(out) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, intent(out) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  Exception handling

      INTEGER      , intent(inout) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

!  local variables
!  ---------------

!  Column vectors for solving BCs

      DOUBLE PRECISION  :: COL2_SWF    (MAXTOTAL,MAX_SURFACEWFS)
      DOUBLE PRECISION  :: SCOL2_SWF   (MAXSTREAMS_2,MAX_SURFACEWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  other help variables

      INTEGER ::          N, IB, Q, M, NSTOKES_ACTUAL
      INTEGER ::          K, KO1, K0, K1, K2
      INTEGER ::          I, J, O1, O2, OM, IR, IROW, IROW1, IROW_S, IROW1_S, C0, CM
      DOUBLE PRECISION :: L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I, REFL_B
      DOUBLE PRECISION :: FACTOR_BRDF, REFL_ATTN, AWF_DIRECT, AWF_EMISS
      DOUBLE PRECISION :: H_1,  H_2,  H_1_S,  H_2_S,  H1,  H2,  DBR
      DOUBLE PRECISION :: H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI, EMISS_VAR
      DOUBLE PRECISION :: H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

!  Help arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  Ground level boundary condition
!  -------------------------------

!  Initialize exception handling

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Initialize

      IB = IBEAM
      M  = FOURIER
      FACTOR_BRDF = SURFACE_FACTOR

!  initialise. Vitally necessary. R. Spurr and V. Natraj, 20 january 2006

      COL2_SWF(1:NTOTAL,1:N_SURFACE_WFS) = ZERO

!  If this is the surface

      N  = NLAYERS
      C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

!  Save some quantities
!  --------------------

!  Only want 1 component for Lambertian

      NSTOKES_ACTUAL = NSTOKES
      IF ( DO_LAMBERTIAN_SURFACE ) NSTOKES_ACTUAL = 1

!  (This is a repetition of earlier code and could be stored)

!  start loops

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES_ACTUAL

!  Beam

          PV_W(J,O1) = WLOWER(J,O1,N) * QUAD_STRMWTS(J)

!  real homogeneous solution contributions

          DO K = 1, K_REAL(N)
            H1 = SOLA_XPOS(J,O1,K,N)
            H2 = SOLB_XNEG(J,O1,K,N)
            HV_P(J,O1,K) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K) = QUAD_STRMWTS(J)*H2
          ENDDO

!  Complex homogeneous solution contributions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            HV_P(J,O1,K1) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K1,N)
            HV_P(J,O1,K2) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K2,N)
            HV_M(J,O1,K1) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K1,N)
            HV_M(J,O1,K2) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K2,N)
          ENDDO

!  End loops

        ENDDO
      ENDDO

!  ===============================
!  1. Set up the BVP Column Vector
!  ===============================

!  Lambertian case
!  ===============

!  Skip if not flagged

! replace GOTO 598 by IF Clause
!      IF ( .not. DO_LAMBERTIAN_SURFACE ) goto 998

      IF ( DO_LAMBERTIAN_SURFACE ) THEN

!  Only 1 weighting function. Q = 1.

        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  OLD CODE, 10 January 2011
!  Beam contribution for this Kernel
!            L_BEAM = R2_BEAM(I,O1)
!  Real homogeneous solutions for this Kernel
!            L_HOM_R = ZERO
!            DO K = 1, K_REAL(N)
!              L_HOM_R = L_HOM_R &
!                + LCON(K,N) * R2_HOMP(I,O1,K) * T_DELT_EIGEN(K,N) + MCON(K,N) * R2_HOMM(I,O1,K)
!            ENDDO
!  Complex homogeneous solutions for this Kernel
!            L_HOM_CR = ZERO
!            KO1 = K_REAL(N) + 1
!            DO K = 1, K_COMPLEX(N)
!              K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
!              H1R =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K1,N) - R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K2,N)
!              H1I =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K2,N) + R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K1,N)
!              L_HOM_CR =  L_HOM_CR + LCON(K1,N) * H1R  - LCON(K2,N) * H1I &
!                                   + MCON(K1,N) * R2_HOMM(I,O1,K1) - MCON(K2,N) * R2_HOMM(I,O1,K2)
!            ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!  Beam contribution for this Kernel, only for NStokes = 1
!  @@@ Rob fix 17jan11, this L_BEAM statement was outside the IF loop
!            L_BEAM = REFL_B * SURFACE_FACTOR

              L_BEAM = ZERO
              IF ( O1 .eq. 1 ) THEN
                REFL_B = SUM(PV_W(1:NSTREAMS,O1))
                L_BEAM = REFL_B * SURFACE_FACTOR
              ENDIF

!  Real homogeneous solutions contribution to this Kernel

              L_HOM_R = ZERO
              IF ( O1 .eq. 1 ) THEN
                DO K = 1, K_REAL(N)
                  H_1 = ZERO ; H_2 = ZERO
                  DO J = 1, NSTREAMS
                    H_1 = H_1 + HV_P(J,O1,K) * SURFACE_FACTOR
                    H_2 = H_2 + HV_M(J,O1,K) * SURFACE_FACTOR
                  ENDDO
                  L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) + MCON(K,N) * H_2
                ENDDO
              ENDIF

!  Complex homogeneous solutions for this Kernel

              L_HOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              IF ( O1 .eq. 1 ) THEN
                DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1 + 1
                  H_1_CR = ZERO ; H_2_CR = ZERO
                  H_1_CI = ZERO ; H_2_CI = ZERO
                  DO J = 1, NSTREAMS
                    H_1_CR = H_1_CR + HV_P(J,O1,K1) * SURFACE_FACTOR
                    H_2_CR = H_2_CR + HV_M(J,O1,K1) * SURFACE_FACTOR
                    H_1_CI = H_1_CI + HV_P(J,O1,K2) * SURFACE_FACTOR
                    H_2_CI = H_2_CI + HV_M(J,O1,K2) * SURFACE_FACTOR
                  ENDDO
                  H1R =   H_1_CR * T_DELT_EIGEN(K1,N) - H_1_CI * T_DELT_EIGEN(K2,N)
                  H1I =   H_1_CR * T_DELT_EIGEN(K2,N) + H_1_CI * T_DELT_EIGEN(K1,N)
                  L_HOM_CR = L_HOM_CR + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I  &
                                      + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI
                ENDDO
              ENDIF

!  Final contribution

              COL2_SWF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  End Streams and Stokes loops

            ENDDO
          ENDDO

!  Add direct beam variation of albedo

          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            O1 = 1
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              IROW = IR + O1
              CM   = C0 + IROW
              COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + ATMOS_ATTN(IB)
            ENDDO
          ENDIF

!  4/9/19. New section for adjusted waterleaving contribution

          if ( FOURIER.eq.0 .and. DO_WATER_LEAVING ) then
             O1 = 1
             DO I = 1, NSTREAMS
                IR = NSTOKES*(I-1)
                IROW = IR + O1
                CM   = C0 + IROW
                COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + LS_TRANS_ATMOS_FINAL(IB,Q) * SL_QUADTERM(I,IB,O1)
             ENDDO
          endif
        
!  @@ Rob Fix: Use ATMOS_ATTN instead of Direct Beam (Unnormalized WF!!)
!              COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + DIRECT_BEAM(I,IB,O1)
!              OOnly need O1 = 1 component (Lambertian)
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          DO I = 1, NSTREAMS
!            IR = NSTOKES*(I-1)
!            DO O1 = 1, NSTOKES
!              IROW = IR + O1 ; CM   = C0 + IROW
!              COL2_SWF(CM,Q) =  COL2_SWF(CM,Q) + DIRECT_BEAM(I,IB,O1)
!            ENDDO
!          ENDDO
!        ENDIF

!  If surface emission, include emissivity variation. This code added for Version 2.4RT

          IF ( DO_INCLUDE_SURFEMISS ) THEN
            O1   = 1
            EMISS_VAR = SURFBB
            DO I = 1, NSTREAMS
              IR   = NSTOKES*(I-1)
              IROW = IR + O1
              CM   = C0 + IROW
              COL2_SWF(CM,Q) = COL2_SWF(CM,Q) - EMISS_VAR
            ENDDO
          ENDIF

!  @@@ Rob Fix FORMER CODE in F77 Version
!             Changed after input from V. Natraj 11/14/2010.
!  @@@ Rob Fix, exclude this line
!          IF ( DO_SOLAR_SOURCES ) EMISS_VAR = SURFBB * PI4
!  @@@ Rob Fix, do not multiply EMISS_VAR by albedo
!          EMISS_VAR = LAMBERTIAN_ALBEDO * EMISS_VAR

!  Copy for the single layer case

          IF ( NLAYERS .EQ. 1 ) THEN
            DO N = 1, NTOTAL
              SCOL2_SWF(N,Q) = COL2_SWF(N,Q)
            ENDDO
          ENDIF

!  End parameter loop

        ENDDO

!  BRDF boundary conditions
!  ========================

      ELSE

!  Continuation point. Removed, Version 2.8
!998   continue

!  Diffuse scatter contributions
!  -----------------------------

!  Start weighting function loop

        DO Q = 1, N_SURFACE_WFS

!  start loops

          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW

!  Beam contribution for this Kernel
!  @@@@@@@@@@@ Rob Fix, 2/9/11, (I,J) instead of (J,I)
!    -- 1/31/21. Version 2.8.3.  LS_BRDF_F, defined locally, drop the "M" Fourier index

              REFL_B = ZERO
              DO J = 1, NSTREAMS
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
!                  REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,OM,J,I)
                  REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,OM,I,J)
                ENDDO
              ENDDO
              L_BEAM = REFL_B * FACTOR_BRDF

!  Real homogeneous solutions contribution to this Kernel
!  @@@@@@@@@@@ Rob Fix, 2/9/11, (I,J) instead of (J,I)
!    -- 1/31/21. Version 2.8.3.  LS_BRDF_F, defined locally, drop the "M" Fourier index

              L_HOM_R = ZERO
              DO K = 1, K_REAL(N)
                H_1 = ZERO ; H_2 = ZERO
                DO J = 1, NSTREAMS
                  H_1_S = ZERO ; H_2_S = ZERO
                  DO O2 = 1, NSTOKES
                    OM  = MUELLER_INDEX(O1,O2)
!                    H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,OM,J,I)
!                    H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,OM,J,I)
                    H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,OM,I,J)
                    H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,OM,I,J)
                  ENDDO
                  H_1 = H_1 + H_1_S
                  H_2 = H_2 + H_2_S
                ENDDO
                H_1 = FACTOR_BRDF * H_1
                H_2 = FACTOR_BRDF * H_2
                L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) + MCON(K,N) * H_2
              ENDDO
    
!  homogeneous complex solutions
!    -- 1/31/21. Version 2.8.3.  LS_BRDF_F, defined locally, drop the "M" Fourier index

              L_HOM_CR = ZERO
              KO1 = K_REAL(N) + 1
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1 + 1
                H_1_CR = ZERO ; H_2_CR = ZERO
                H_1_CI = ZERO ; H_2_CI = ZERO
                DO J = 1, NSTREAMS
                  H_1_S_CR = ZERO ; H_2_S_CR = ZERO
                  H_1_S_CI = ZERO ; H_2_S_CI = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    DBR = LS_BRDF_F(Q,OM,J,I)
                    H_1_S_CR = H_1_S_CR + HV_P(J,O2,K1) * DBR
                    H_2_S_CR = H_2_S_CR + HV_M(J,O2,K1) * DBR
                    H_1_S_CI = H_1_S_CI + HV_P(J,O2,K2) * DBR
                    H_2_S_CI = H_2_S_CI + HV_M(J,O2,K2) * DBR
                  ENDDO
                  H_1_CR = H_1_CR + H_1_S_CR
                  H_2_CR = H_2_CR + H_2_S_CR
                  H_1_CI = H_1_CI + H_1_S_CI
                  H_2_CI = H_2_CI + H_2_S_CI
                ENDDO
                H_1_CR = FACTOR_BRDF  * H_1_CR
                H_1_CI = FACTOR_BRDF  * H_1_CI
                H_2_CR = FACTOR_BRDF  * H_2_CR
                H_2_CI = FACTOR_BRDF  * H_2_CI
                H1R =   H_1_CR * T_DELT_EIGEN(K1,N) - H_1_CI * T_DELT_EIGEN(K2,N)
                H1I =   H_1_CR * T_DELT_EIGEN(K2,N) + H_1_CI * T_DELT_EIGEN(K1,N)
                L_HOM_CR = L_HOM_CR + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I &
                                    + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI
              ENDDO
    
!  Final contribution

              COL2_SWF(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

!  End loops

            ENDDO
          ENDDO

!  Direct beam reflection
!    -- 1/31/21. Version 2.8.3.  LS_BRDF_F_0, defined locally, drop the "M" Fourier index

          IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            REFL_ATTN = ATMOS_ATTN(IB)
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM   = C0 + IROW
                OM = MUELLER_INDEX(O1,1)
                AWF_DIRECT = REFL_ATTN * LS_BRDF_F_0(Q,OM,I,IB)
                COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + AWF_DIRECT
              ENDDO
            ENDDO
          ENDIF

!  4/9/19. New section for adjusted waterleaving contribution

          if ( FOURIER.eq.0 .and. DO_WATER_LEAVING ) then
             O1 = 1
             DO I = 1, NSTREAMS
                IR = NSTOKES*(I-1)
                IROW = IR + O1
                CM   = C0 + IROW
                COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + LS_TRANS_ATMOS_FINAL(IB,Q) * SL_QUADTERM(I,IB,O1)
             ENDDO
          endif
        
!  If surface emission, include emissivity variation, BRDF surface

          IF ( DO_INCLUDE_SURFEMISS ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM   = C0 + IROW
                AWF_EMISS = SURFBB * LS_EMISSIVITY(Q,O1,I)
                COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + AWF_EMISS
              ENDDO
            ENDDO
          ENDIF

!  Copy for the single layer case

          IF ( NLAYERS .EQ. 1 ) THEN
             DO N = 1, NTOTAL
              SCOL2_SWF(N,Q) = COL2_SWF(N,Q)
            ENDDO
          ENDIF

!  End parameter loop

        ENDDO

!  End BRDF vs. Lambertian clause

      ENDIF

!  ================
!  2. SOLVE THE BVP
!  ================

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_SWF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in SURFACEWF_BVP_SOLUTION'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, all layer

        DO Q = 1, N_SURFACE_WFS
          DO N = 1, NLAYERS
            C0 = (N-1)*NSTKS_NSTRMS_2
            DO K = 1, K_REAL(N)
              IROW = K
              IROW1 = IROW + NSTKS_NSTRMS
              NCON_SWF(Q,K,N) = COL2_SWF(C0+IROW,Q)
              PCON_SWF(Q,K,N) = COL2_SWF(C0+IROW1,Q)
            ENDDO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              IROW    = K + K_REAL(N)
              IROW1   = IROW + NSTKS_NSTRMS
              IROW_S  = IROW + K_COMPLEX(N)
              IROW1_S = IROW_S + NSTKS_NSTRMS
              NCON_SWF(Q,K1,N) = COL2_SWF(C0+IROW,   Q)
              NCON_SWF(Q,K2,N) = COL2_SWF(C0+IROW_S, Q)
              PCON_SWF(Q,K1,N) = COL2_SWF(C0+IROW1,  Q)
              PCON_SWF(Q,K2,N) = COL2_SWF(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_SWF

        CALL DGETRS &
           ( 'N', NTOTAL, N_SURFACE_WFS, SMAT2, MAXSTRMSTKS_2, &
              SIPIVOT, SCOL2_SWF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (Reg. 1 layer) in SURFACEWF_BVP_SOLUTION'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, 1 layer

        DO Q = 1, N_SURFACE_WFS
          N = 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_SWF(Q,K,N) = SCOL2_SWF(IROW,Q)
            PCON_SWF(Q,K,N) = SCOL2_SWF(IROW1,Q)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_SWF(Q,K1,N) = SCOL2_SWF(IROW,   Q)
            NCON_SWF(Q,K2,N) = SCOL2_SWF(IROW_S, Q)
            PCON_SWF(Q,K1,N) = SCOL2_SWF(IROW1,  Q)
            PCON_SWF(Q,K2,N) = SCOL2_SWF(IROW1_S,Q)
          ENDDO
        ENDDO

!  end clause

      ENDIF

!  End subroutine

      RETURN
      END SUBROUTINE SURFACEWF_BVP_SOLUTION

!

      SUBROUTINE SURFACEWF_BVPTEL_SOLUTION &
            ( DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFACE,  FOURIER, IBEAM,     & ! Input Flags/indices
              NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,                      & ! Input Numbers
              NSTKS_NSTRMS, NSTKS_NSTRMS_2,                                   & ! Input
              N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,                 & ! BVPTel Control
              N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,                   & ! BVPTel Control
              MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, QUAD_STRMWTS, & ! Input
              LS_BRDF_F, LS_BRDF_F_0, T_DELT_DISORDS, ATMOS_ATTN,             & ! Inputs
              T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WLOWER, LCON, MCON,         & ! Input
              BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,                         & ! Input BVPTel
              NCON_SWF, PCON_SWF, STATUS, MESSAGE, TRACE )                      ! Output

!  Regular-BVP solution for the surface weighting functions

!  1/31/21. Version 2.8.3. BRDF arrays LS_BRDF_F, LS_BRDF_F_0 are defined locally for each Fourier

      USE VLIDORT_pars_m, only : MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAX_SURFACEWFS, &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, MAXSTOKES_SQ,     &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE

      USE lapack_tools_m, only : DGBTRS, DGETRS

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFACE

!  Basic control numbers

      INTEGER, INTENT (IN) ::          FOURIER, IBEAM
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  other derived numbers

      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Telescoped BVP Control

      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::          NLAYERS_TEL
      INTEGER, INTENT (IN) ::          ACTIVE_LAYERS ( MAXLAYERS )

!  Bookkeeping

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  discrete ordinate stuff

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  Surface stuff
!    -- 1/31/21. Version 2.8.3. LS_BRDF_F, LS_BRDF_F_0 defined locally, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )

!  RTE solution variables

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  BVP Matrices

      DOUBLE PRECISION, INTENT (IN) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOTTEL   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2   ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTRMSTKS_2 )

!  Attenuation to the surface

      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )

!  outputs
!  -------

!  Linearized Solution constants of integration

      DOUBLE PRECISION, intent(out) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, intent(out) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  Exception handling

      INTEGER      , intent(inout) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

!  local variables
!  ---------------

!  Column vectors for solving BCs

      DOUBLE PRECISION  :: COLTEL2_SWF    (MAXTOTAL,MAX_SURFACEWFS)
      DOUBLE PRECISION  :: SCOL2_SWF      (MAXSTREAMS_2,MAX_SURFACEWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  other help variables

      INTEGER ::          N, N1, IB, Q, M
      INTEGER ::          NAF, NAL, NAL1, NS
      INTEGER ::          K, KO1, K0, K1, K2
      INTEGER ::          I, I1, J, O1, O2, OM, IR, IROW, IROW1, IROW_S, IROW1_S, C0, CM, ICOW, IC
      DOUBLE PRECISION :: L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I, REFL_B, L_HOM1, L_HOM2, SHOM, SHOM_R
      DOUBLE PRECISION :: FACTOR_BRDF, REFL_ATTN, AWF_DIRECT
      DOUBLE PRECISION :: H_1,  H_2,  H_1_S,  H_2_S,  H1,  H2,  DBR
      DOUBLE PRECISION :: H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI
      DOUBLE PRECISION :: H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI
      DOUBLE PRECISION :: SHOM_CR, L_HOM1CR, L_HOM2CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, PXR1, NXR2, PXR2

!  Help arrays

      DOUBLE PRECISION :: PV_W ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: HV_P ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: HV_M ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: CUMTRANS(MAXSTREAMS)

!  Ground level boundary condition
!  -------------------------------

!  Initialize exception handling

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Initialize

      IB = IBEAM
      M  = FOURIER
      FACTOR_BRDF = SURFACE_FACTOR

!  Initialise. Vitally necessary. R. Spurr and V. Natraj, 20 january 2006

      COLTEL2_SWF(1:N_BVTELMATRIX_SIZE,1:MAX_SURFACEWFS) = ZERO

!  Initialise Ground level boundary condition

      M  = FOURIER
      NS = NLAYERS_TEL
      N = ACTIVE_LAYERS(NS)
      C0 = (NLAYERS_TEL-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      IB = IBEAM

!  Cumulative transmittance again

      CUMTRANS(1:NSTREAMS) = ONE
      DO N1 = NLAYERS, N+1, -1
         CUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
      ENDDO

!  Bottom boundary condition
!  =========================

!  (This is a repetition of earlier code and could be stored)

!  Start loops

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES

!  Beam

          PV_W(J,O1) = WLOWER(J,O1,N) * QUAD_STRMWTS(J)

!  Real homogeneous solution contributions

          DO K = 1, K_REAL(N)
            H1 = SOLA_XPOS(J,O1,K,N)
            H2 = SOLB_XNEG(J,O1,K,N)
            HV_P(J,O1,K) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K) = QUAD_STRMWTS(J)*H2
          ENDDO

!  Complex homogeneous solution contributions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            HV_P(J,O1,K1) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K1,N)
            HV_P(J,O1,K2) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K2,N)
            HV_M(J,O1,K1) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K1,N)
            HV_M(J,O1,K2) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K2,N)
          ENDDO

!  End loops

        ENDDO
      ENDDO

!  BRDF boundary conditions
!  ========================

!  Diffuse scatter contributions
!  -----------------------------

!  Start weighting function loop

      DO Q = 1, N_SURFACE_WFS

!  Start loops

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW

!  Beam contribution for this Kernel
!  @@@@@@@@@@@ Rob Fix, 2/9/11, (I,J) instead of (J,I)
!    -- 1/31/21. Version 2.8.3. LS_BRDF_F defined locally , drop Fourier index "M"

            REFL_B = ZERO
            DO J = 1, NSTREAMS
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
!                REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,OM,J,I)
                REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,OM,I,J)
              ENDDO
            ENDDO
            L_BEAM = REFL_B * FACTOR_BRDF

!  Real homogeneous solutions contribution to this Kernel
!  @@@@@@@@@@@ Rob Fix, 2/9/11, (I,J) instead of (J,I)
!    -- 1/31/21. Version 2.8.3. LS_BRDF_F defined locally , drop Fourier index "M"

            L_HOM_R = ZERO
            DO K = 1, K_REAL(N)
              H_1 = ZERO ; H_2 = ZERO
              DO J = 1, NSTREAMS
                H_1_S = ZERO ; H_2_S = ZERO
                DO O2 = 1, NSTOKES
                  OM  = MUELLER_INDEX(O1,O2)
!                  H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,OM,J,I)
!                  H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,OM,J,I)
                  H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,OM,I,J)
                  H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,OM,I,J)
                ENDDO
                H_1 = H_1 + H_1_S
                H_2 = H_2 + H_2_S
              ENDDO
              H_1 = FACTOR_BRDF * H_1
              H_2 = FACTOR_BRDF * H_2
              L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) + MCON(K,N) * H_2
            ENDDO

!  Homogeneous complex solutions
!    -- 1/31/21. Version 2.8.3. LS_BRDF_F defined locally , drop Fourier index "M"

            L_HOM_CR = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1 + 1
              H_1_CR = ZERO ; H_2_CR = ZERO
              H_1_CI = ZERO ; H_2_CI = ZERO
              DO J = 1, NSTREAMS
                H_1_S_CR = ZERO ; H_2_S_CR = ZERO
                H_1_S_CI = ZERO ; H_2_S_CI = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  DBR = LS_BRDF_F(Q,OM,J,I)
                  H_1_S_CR = H_1_S_CR + HV_P(J,O2,K1) * DBR
                  H_2_S_CR = H_2_S_CR + HV_M(J,O2,K1) * DBR
                  H_1_S_CI = H_1_S_CI + HV_P(J,O2,K2) * DBR
                  H_2_S_CI = H_2_S_CI + HV_M(J,O2,K2) * DBR
                ENDDO
                H_1_CR = H_1_CR + H_1_S_CR
                H_2_CR = H_2_CR + H_2_S_CR
                H_1_CI = H_1_CI + H_1_S_CI
                H_2_CI = H_2_CI + H_2_S_CI
              ENDDO
              H_1_CR = FACTOR_BRDF  * H_1_CR
              H_1_CI = FACTOR_BRDF  * H_1_CI
              H_2_CR = FACTOR_BRDF  * H_2_CR
              H_2_CI = FACTOR_BRDF  * H_2_CI
              H1R =   H_1_CR * T_DELT_EIGEN(K1,N) - H_1_CI * T_DELT_EIGEN(K2,N)
              H1I =   H_1_CR * T_DELT_EIGEN(K2,N) + H_1_CI * T_DELT_EIGEN(K1,N)
              L_HOM_CR = L_HOM_CR + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I &
                                  + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI
            ENDDO

!  Final contribution

            COLTEL2_SWF(CM,Q) = ( L_BEAM + L_HOM_R + L_HOM_CR ) * CUMTRANS(I)

!  End loops

          ENDDO
        ENDDO

!  Direct beam reflection
!    -- 1/31/21. Version 2.8.3. LS_BRDF_F_0 defined locally , drop Fourier index "M"

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          REFL_ATTN = ATMOS_ATTN(IB)
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              OM = MUELLER_INDEX(O1,1)
              AWF_DIRECT = CUMTRANS(I) * REFL_ATTN * LS_BRDF_F_0(Q,OM,I,IB)
              COLTEL2_SWF(CM,Q) = COLTEL2_SWF(CM,Q) + AWF_DIRECT
            ENDDO
          ENDDO
        ENDIF

!  Copy for the single layer case

        IF ( NLAYERS_TEL .EQ. 1 ) THEN
           DO N = 1, NSTKS_NSTRMS_2
            SCOL2_SWF(N,Q) = COLTEL2_SWF(N,Q)
          ENDDO
        ENDIF

!  End parameter loop

      ENDDO

!  ================
!  2. SOLVE THE BVP
!  ================

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COLTEL2_SWF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUBDIAG, &
              N_BVTELMATRIX_SUPDIAG, N_SURFACE_WFS, &
              BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
              COLTEL2_SWF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in SURFACEWF_BVPTEL_SOLUTION'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, all ACTIVE LAYERS

        DO Q = 1, N_SURFACE_WFS
          C0 = - NSTKS_NSTRMS_2
          DO NS = 1, NLAYERS_TEL
            N = ACTIVE_LAYERS(NS)
            C0 = C0 + NSTKS_NSTRMS_2
            DO K = 1, K_REAL(N)
              IROW = K
              IROW1 = IROW + NSTKS_NSTRMS
              NCON_SWF(Q,K,N) = COLTEL2_SWF(C0+IROW,Q)
              PCON_SWF(Q,K,N) = COLTEL2_SWF(C0+IROW1,Q)
            ENDDO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              IROW    = K + K_REAL(N)
              IROW1   = IROW + NSTKS_NSTRMS
              IROW_S  = IROW + K_COMPLEX(N)
              IROW1_S = IROW_S + NSTKS_NSTRMS
              NCON_SWF(Q,K1,N) = COLTEL2_SWF(C0+IROW,   Q)
              NCON_SWF(Q,K2,N) = COLTEL2_SWF(C0+IROW_S, Q)
              PCON_SWF(Q,K1,N) = COLTEL2_SWF(C0+IROW1,  Q)
              PCON_SWF(Q,K2,N) = COLTEL2_SWF(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

        NAL = ACTIVE_LAYERS(1)

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_SWF

        CALL DGETRS &
           ( 'N', NSTKS_NSTRMS_2, N_SURFACE_WFS, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              SCOL2_SWF, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (Reg. 1 layer) in SURFACEWF_BVPTEL_SOLUTION'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, 1 layer

        DO Q = 1, N_SURFACE_WFS
          DO K = 1, K_REAL(NAL)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_SWF(Q,K,NAL) = SCOL2_SWF(IROW,Q)
            PCON_SWF(Q,K,NAL) = SCOL2_SWF(IROW1,Q)
          ENDDO
          KO1 = K_REAL(NAL) + 1
          DO K = 1, K_COMPLEX(NAL)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(NAL)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(NAL)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_SWF(Q,K1,NAL) = SCOL2_SWF(IROW,   Q)
            NCON_SWF(Q,K2,NAL) = SCOL2_SWF(IROW_S, Q)
            PCON_SWF(Q,K1,NAL) = SCOL2_SWF(IROW1,  Q)
            PCON_SWF(Q,K2,NAL) = SCOL2_SWF(IROW1_S,Q)
          ENDDO
        ENDDO

!  End clause

      ENDIF

!  ==============================================================
!  3. Set Linearized integration constants for non-active layers
!  ==============================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  ------------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required, first layer is always active

!mick fix 9/19/2017 - rotated indices for NCON_SWF & PCON_SWF in several locations in this section

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1

!  Start stream, stokes and parameter loops

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            ICOW = IC + O1

            DO Q = 1, N_SURFACE_WFS

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(NAF)
                !L_HOM1 = NCON_SWF(K,NAF,Q) * SOLA_XPOS(I1,O1,K,NAF)
                !L_HOM2 = T_DELT_EIGEN(K,NAF) * PCON_SWF(K,NAF,Q) * SOLB_XNEG(I1,O1,K,NAF)
                L_HOM1 = NCON_SWF(Q,K,NAF) * SOLA_XPOS(I1,O1,K,NAF)
                L_HOM2 = T_DELT_EIGEN(K,NAF) * PCON_SWF(Q,K,NAF) * SOLB_XNEG(I1,O1,K,NAF)
                SHOM_R = SHOM_R + L_HOM1 + L_HOM2
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(NAF) + 1
              DO K = 1, K_COMPLEX(NAF)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                !NXR1  = NCON_SWF(K1,NAF,Q) * SOLA_XPOS(I1,O1,K1,NAF) - NCON_SWF(K2,NAF,Q) * SOLA_XPOS(I1,O1,K2,NAF)
                !PXR1  = PCON_SWF(K1,NAF,Q) * SOLB_XNEG(I1,O1,K1,NAF) - PCON_SWF(K2,NAF,Q) * SOLB_XNEG(I1,O1,K2,NAF)
                !PXR2  = PCON_SWF(K1,NAF,Q) * SOLB_XNEG(I1,O1,K2,NAF) + PCON_SWF(K2,NAF,Q) * SOLB_XNEG(I1,O1,K1,NAF)
                NXR1  = NCON_SWF(Q,K1,NAF) * SOLA_XPOS(I1,O1,K1,NAF) - NCON_SWF(Q,K2,NAF) * SOLA_XPOS(I1,O1,K2,NAF)
                PXR1  = PCON_SWF(Q,K1,NAF) * SOLB_XNEG(I1,O1,K1,NAF) - PCON_SWF(Q,K2,NAF) * SOLB_XNEG(I1,O1,K2,NAF)
                PXR2  = PCON_SWF(Q,K1,NAF) * SOLB_XNEG(I1,O1,K2,NAF) + PCON_SWF(Q,K2,NAF) * SOLB_XNEG(I1,O1,K1,NAF)
                L_HOM1CR = NXR1 
                L_HOM2CR =  T_DELT_EIGEN(K1,NAF) * PXR1 - T_DELT_EIGEN(K2,NAF)   * PXR2
                SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
              ENDDO

!  Total

              SHOM = SHOM_R + SHOM_CR
              !PCON_SWF(ICOW,N1,Q) = SHOM
              !NCON_SWF(ICOW,N1,Q) = ZERO
              PCON_SWF(Q,ICOW,N1) = SHOM
              NCON_SWF(Q,ICOW,N1) = ZERO

!  End loops

            ENDDO
          ENDDO
        ENDDO

!  End active layer varying

      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards
!    Additional linearizations if you are passing through the varying layer

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            !NCON_SWF(IROW,N,1:N_SURFACE_WFS) = ZERO
            !PCON_SWF(IROW,N,1:N_SURFACE_WFS) = T_DELT_DISORDS(I,N1) * PCON_SWF(IROW,N,N_SURFACE_WFS)
            NCON_SWF(1:N_SURFACE_WFS,IROW,N) = ZERO

!mick fix 9/19/2017 - changed PCON_SWF indices on RHS from (N_SURFACE_WFS,IROW,N) to (1:N_SURFACE_WFS,IROW,N1)
            PCON_SWF(1:N_SURFACE_WFS,IROW,N) = T_DELT_DISORDS(I,N1) * PCON_SWF(1:N_SURFACE_WFS,IROW,N1)
          ENDDO
        ENDDO
      ENDDO

!  Transmittance layers BELOW active layer(s)
!  ------------------------------------------

!  This section substantially revised for Version 2.8
!   - General surface treatment for the Telescpied BVP

!       ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- NCON values are  propagated downwards from bottom of last active layer
!   -- PCON values also propagated downwards, BUT only present if surface condition
!  1.   Require linearized solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

!mick fix 9/19/2017 - rotated indices for NCON_SWF & PCON_SWF in several locations in this section

      NAL = ACTIVE_LAYERS(NLAYERS_TEL) ; NAL1 = NAL + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  N-constants, always required

!  Start stream, stokes, parameter loops

        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            DO Q = 1, N_SURFACE_WFS

!  Real homogeneous solutions
!   Needs checking------------------------- 22 March 2007

              SHOM_R = ZERO
              DO K = 1, K_REAL(NAL)
                !NXR  = NCON_SWF(K,NAL,Q) * SOLA_XPOS(I,O1,K,NAL)
                !PXR  = PCON_SWF(K,NAL,Q) * SOLB_XNEG(I,O1,K,NAL)
                NXR  = NCON_SWF(Q,K,NAL) * SOLA_XPOS(I,O1,K,NAL)
                PXR  = PCON_SWF(Q,K,NAL) * SOLB_XNEG(I,O1,K,NAL)
                SHOM_R = SHOM_R + PXR + T_DELT_EIGEN(K,NAL) * NXR
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              KO1 = K_REAL(NAL) + 1
              DO K = 1, K_COMPLEX(NAL)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                !NXR1  =   NCON_SWF(K1,NAL,Q) * SOLA_XPOS(I,O1,K1,NAL) - NCON_SWF(K2,NAL,Q) * SOLA_XPOS(I,O1,K2,NAL)
                !NXR2  =   NCON_SWF(K1,NAL,Q) * SOLA_XPOS(I,O1,K2,NAL) + NCON_SWF(K2,NAL,Q) * SOLA_XPOS(I,O1,K1,NAL)
                !PXR1  =   PCON_SWF(K1,NAL,Q) * SOLB_XNEG(I,O1,K1,NAL) - PCON_SWF(K2,NAL,Q) * SOLB_XNEG(I,O1,K2,NAL)
                NXR1  =   NCON_SWF(Q,K1,NAL) * SOLA_XPOS(I,O1,K1,NAL) - NCON_SWF(Q,K2,NAL) * SOLA_XPOS(I,O1,K2,NAL)
                NXR2  =   NCON_SWF(Q,K1,NAL) * SOLA_XPOS(I,O1,K2,NAL) + NCON_SWF(Q,K2,NAL) * SOLA_XPOS(I,O1,K1,NAL)
                PXR1  =   PCON_SWF(Q,K1,NAL) * SOLB_XNEG(I,O1,K1,NAL) - PCON_SWF(Q,K2,NAL) * SOLB_XNEG(I,O1,K2,NAL)
                L_HOM2CR = PXR1
                L_HOM1CR = T_DELT_EIGEN(K1,NAL) * NXR1 - T_DELT_EIGEN(K2,NAL) * NXR2
                SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
              ENDDO

!  Total

              SHOM = SHOM_R + SHOM_CR
              !NCON_SWF(IROW,NAL1,Q) = SHOM
              NCON_SWF(Q,IROW,NAL1) = SHOM

!  End loops

            ENDDO 
          ENDDO
        ENDDO

!  Other layers to bottom of medium: propagate downwards.
!   Additional variation, since you are passing through varying layers.

        DO N = NAL + 2, NLAYERS
          N1 = N - 1
          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              !NCON_SWF(IROW,N,1:N_SURFACE_WFS) = T_DELT_DISORDS(I,N1) * NCON_SWF(IROW,N,1:N_SURFACE_WFS)
!mick fix 9/19/2017 - changed NCON_SWF indices on RHS from (1:N_SURFACE_WFS,IROW,N) to (1:N_SURFACE_WFS,IROW,N1)
              NCON_SWF(1:N_SURFACE_WFS,IROW,N) = T_DELT_DISORDS(I,N1) * NCON_SWF(1:N_SURFACE_WFS,IROW,N1)
            ENDDO
          ENDDO
        ENDDO

!  P-Constants need to be determined if there is a surface condition. Otherwise zero.

        IF ( DO_INCLUDE_SURFACE ) THEN

!  Start stream, stokes and parameter loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              ICOW = IC + O1
              DO Q = 1, N_SURFACE_WFS

!  Real homogeneous solutions

                SHOM_R = ZERO
                DO K = 1, K_REAL(NAL)
                  !NXR  = NCON_SWF(K,NAL,Q) * SOLA_XPOS(I1,O1,K,NAL)
                  !PXR  = PCON_SWF(K,NAL,Q) * SOLB_XNEG(I1,O1,K,NAL)
                  NXR  = NCON_SWF(Q,K,NAL) * SOLA_XPOS(I1,O1,K,NAL)
                  PXR  = PCON_SWF(Q,K,NAL) * SOLB_XNEG(I1,O1,K,NAL)
                  SHOM_R = SHOM_R + NXR + T_DELT_EIGEN(K,NAL) * PXR
                ENDDO

!  Complex homogeneous solutions

                SHOM_CR = ZERO
                KO1 = K_REAL(NAL) + 1
                DO K = 1, K_COMPLEX(NAL)
                  K0 = 2*K - 2
                  K1 = KO1 + K0
                  K2 = K1  + 1
                  !NXR1  = NCON_SWF(K1,NAL,Q) * SOLA_XPOS(I1,O1,K1,NAL) - NCON_SWF(K2,NAL,Q) * SOLA_XPOS(I1,O1,K2,NAL)
                  !PXR1  = PCON_SWF(K1,NAL,Q) * SOLB_XNEG(I1,O1,K1,NAL) - PCON_SWF(K2,NAL,Q) * SOLB_XNEG(I1,O1,K2,NAL)
                  !PXR2  = PCON_SWF(K1,NAL,Q) * SOLB_XNEG(I1,O1,K2,NAL) + PCON_SWF(K2,NAL,Q) * SOLB_XNEG(I1,O1,K1,NAL)
                  NXR1  = NCON_SWF(Q,K1,NAL) * SOLA_XPOS(I1,O1,K1,NAL) - NCON_SWF(Q,K2,NAL) * SOLA_XPOS(I1,O1,K2,NAL)
                  PXR1  = PCON_SWF(Q,K1,NAL) * SOLB_XNEG(I1,O1,K1,NAL) - PCON_SWF(Q,K2,NAL) * SOLB_XNEG(I1,O1,K2,NAL)
                  PXR2  = PCON_SWF(Q,K1,NAL) * SOLB_XNEG(I1,O1,K2,NAL) + PCON_SWF(Q,K2,NAL) * SOLB_XNEG(I1,O1,K1,NAL)
                  L_HOM1CR = NXR1
                  L_HOM2CR = T_DELT_EIGEN(K1,NAL) * PXR1 - T_DELT_EIGEN(K2,NAL) * PXR2
                  SHOM_CR = SHOM_CR + L_HOM1CR + L_HOM2CR
                ENDDO

!  Total

                SHOM = SHOM_R + SHOM_CR
                !PCON_SWF(ICOW,NAL1,Q) = SHOM / T_DELT_DISORDS(I,NAL1)
                PCON_SWF(Q,ICOW,NAL1) = SHOM / T_DELT_DISORDS(I,NAL1)

!  End loops

              ENDDO
            ENDDO
          ENDDO

!  Other constants propagated

          DO N = NAL + 2, NLAYERS
            N1 = N - 1
            DO I = 1, NSTREAMS
              I1  = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                DO Q = 1, N_SURFACE_WFS
                  !PCON_SWF(ICOW,N,Q) = PCON_SWF(ICOW,N1,Q)  / T_DELT_DISORDS(I,N)
                  PCON_SWF(Q,ICOW,N) = PCON_SWF(Q,ICOW,N1)  / T_DELT_DISORDS(I,N)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  Otherwise all P-constants are zero
        
        ELSE

          DO N = NAL1, NLAYERS
            DO I = 1, NSTREAMS
              I1  = I + NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                DO Q = 1, N_SURFACE_WFS
                  !PCON_SWF(ICOW,N,Q) = ZERO
                  PCON_SWF(Q,ICOW,N) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  End general surface treatment

        ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  End subroutine

      RETURN
      END SUBROUTINE SURFACEWF_BVPTEL_SOLUTION

!

      SUBROUTINE SURFACEWF_POSTPROCESS_MASTER ( &
          DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_MSSTS,                       & ! Input flags (general)
          DO_DBCORRECTION, DO_INCLUDE_MVOUTPUT, DO_OBSERVATION_GEOMETRY,               & ! Input flags (general)
          DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,               & ! Input flags (thermal)
          DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,             & ! Input flags (surface)
          FOURIER, IBEAM, ECHT_NSTOKES, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS, & ! input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          & ! Input Level output control
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                & ! input partlayer control
          N_DIRECTIONS, WHICH_DIRECTIONS, MUELLER_INDEX,                               & ! Input bookkeeping
          FLUX_MULTIPLIER, FLUXVEC, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,        & ! input Flux + Quads
          SL_USERTERM, ALBEDO, USER_BRDF_F, SURFBB, LS_USER_EMISSIVITY, LS_EMISSIVITY, & ! input Surface stuff
          LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL,           & ! input Surface stuff
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,    & ! Input transmittances (User/Disords)
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                       & ! Input Homog. solutions
          T_UTUP_EIGEN, T_UTDN_EIGEN, UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,      & ! Input Homog. Solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,  UT_HMULT_DU, UT_HMULT_DD,       & ! Input Homog.Multipliers
          ATMOS_ATTN, STOKES_DOWNSURF, NCON_SWF, PCON_SWF,                             & ! BOA downwelling, Linearized BVP
          SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F, MINT_SURFACEWF, FLUX_SURFACEWF)  ! Output

!  1/31/21. Version 2.8.3. Several minor changes
!    -- Drop LOCAL_UM_START; all UM loops start with 1
!    -- Use the post processing mask system, don't need separate OBSGEOM calculation
!    -- Introduce DO_MSSTS flag, calculate linearized MSSTS output LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F
!    -- BRDF arrays USER_BRDF_F, LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0 all defined locally each Fourier

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, MAX_SZANGLES,             &
                                  MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_DIRECTIONS, MAX_SURFACEWFS, &
                                  MAXSTOKES_SQ, MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, UPIDX, DNIDX, ZERO 


      IMPLICIT NONE

!  Flags

!      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT  !   removed version 2.8

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      
!  4/9/19. Direct surface inclusion flags. replaces DO_INCLUDE_DIRECTBEAM

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTRF
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTSL
!      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!  Indices

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Post-processing mask,replace LOCAL_UM_START and N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Whole layer bookkeeping

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  Partial layer bookkeeping

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Bookkeeping

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Flux and Quadrature

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Surface reflection
!    incident quadrature streams, reflected user streams
!    -- BRDF array USER_BRDF_F  defined locally each Fourier, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: SURFBB

!  4/9/19. Adjusted water-leaving inputs

      DOUBLE PRECISION, INTENT (IN) :: SL_USERTERM(MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES)
      DOUBLE PRECISION, INTENT (IN) :: LS_TRANS_ATMOS_FINAL(MAXBEAMS,MAX_SURFACEWFS)

!  Linearized surface BRDF and emissivity
!    -- BRDF arrays defined locally each Fourier, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )

!  Homogeneous solution

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  User-angle homogeneous solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  homogeneous solution multiplier

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

! BOA source contributions

      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )

!  linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  Linearized surface output
!  -------------------------

      DOUBLE PRECISION, INTENT (INOUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  Local variables
!  ---------------

!  Linearized BOA terms

      DOUBLE PRECISION :: LS_BOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_BOA_THTONLY_SOURCE ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

!  Other local variables

      INTEGER ::          LUM, UM

!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling weighting functions
!  -----------------------------

      IF ( DO_UPWELLING ) THEN

!  Get the surface term (LS_BOA_SOURCE).
!  4/19/19. Version 2.8.1. Introduce INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL (two directbeam flags)
!  4/19/19. Version 2.8.1. Introduce SL_USERTERM and LS_TRANS_ATMOS_FINAL, to handle waterleaving linearization

!  1/31/21. Version 2.8.3.
!     -- Use post-processing mask N_PPSTREAMS, PPSTREAM_MASK, replace LOCAL_UM_START, drop DO_OBSERVATION_GEOMETRY
!     -- Define all BRDF inputs locally for each Fourier
!mick fix 1/5/2021 - added ECHT_NSTOKES to inputs

        CALL BOA_SURFACEWF ( &
          DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_DBCORRECTION, DO_LAMBERTIAN_SURFACE, & ! Input flags (general/surface) 
          DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,                & ! Input flags (thermal)
          DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, FOURIER, IBEAM,                     & ! Input flags, indices
          ECHT_NSTOKES, NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,                      & ! input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK,                                                   & ! input Control numbers
          MUELLER_INDEX, FLUXVEC, SURFACE_FACTOR, QUAD_STRMWTS,                         & ! input Bookkeeping + Quads
          ALBEDO, USER_BRDF_F, SL_USERTERM, SURFBB, LS_BRDF_F, LS_TRANS_ATMOS_FINAL,    & ! input Surface stuff
          LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_USER_EMISSIVITY, LS_EMISSIVITY,          & ! input Surface stuff
          T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                        & ! RTE Homog. solutions
          ATMOS_ATTN, STOKES_DOWNSURF, NCON_SWF, PCON_SWF,                              & ! BOA downwelling, Linearized BVP
          LS_BOA_SOURCE, LS_BOA_THTONLY_SOURCE )                                          ! Output

!  Upwelling surface weighting function field - post processing

!  1/31/21. Version 2.8.3. 
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F), now output
!mick fix 1/5/2021 - added ECHT_NSTOKES to inputs

        CALL UPUSER_SURFACEWF ( &
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_THERMAL_TRANSONLY, DO_MSSTS, IBEAM, & ! Input Flags
          ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,                    & ! Input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK,                                                      & ! Input Control numbers
          FLUX_MULTIPLIER, K_REAL, K_COMPLEX, UTAU_LEVEL_MASK_UP,                          & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                    & ! input partlayer control
          T_DELT_USERM, T_UTUP_USERM, UHOM_UPDN, UHOM_UPUP, LS_BOA_SOURCE,                 & ! Input User RT Solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, NCON_SWF, PCON_SWF,                  & ! Input Multipliers/Linearized-BVP
          SURFACEWF_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F )                                   ! Output

!mick temp fix 9/7/2012 - ELSE added, Zeroes the output if Upwelling flag not set
!  -- 1/31/21. Version 2.8.3. Use post processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present)
!                  - replaced NSTOKES with ECHT_NSTOKES here 
!                  - initialized LS_LAYER_MSSTS_F & LS_SURF_MSSTS_F 

      ELSE
         DO LUM = 1, N_PPSTREAMS
            !UM = PPSTREAM_MASK(LUM,IBEAM)
            SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,LUM,IBEAM,1:ECHT_NSTOKES,UPIDX) = ZERO
         ENDDO
         LS_LAYER_MSSTS_F(IBEAM,1:ECHT_NSTOKES,1:NLAYERS,1:N_SURFACE_WFS) = ZERO
         LS_SURF_MSSTS_F (IBEAM,1:ECHT_NSTOKES,1:N_SURFACE_WFS ) = ZERO
      ENDIF

!  Downwelling surface weighting function field - post processing

      IF ( DO_DNWELLING ) THEN

!  1/31/21. Version 2.8.3. 
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F), now output
!mick fix 1/5/2021 - added ECHT_NSTOKES & NLAYERS to inputs

        CALL DNUSER_SURFACEWF ( &
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_THERMAL_TRANSONLY, DO_MSSTS, & ! Input Flags
          IBEAM, ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,      & ! Input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK,                                               & ! Input Control numbers
          FLUX_MULTIPLIER, K_REAL, K_COMPLEX, UTAU_LEVEL_MASK_DN,                   & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,             & ! input partlayer control
          T_DELT_USERM, T_UTDN_USERM, UHOM_DNDN, UHOM_DNUP,                         & ! Input RT Solutions
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD, NCON_SWF, PCON_SWF,           & ! Input Multipliers/IntegConsts
          SURFACEWF_F, LS_LAYER_MSSTS_F )                                             ! Output

!mick temp fix 9/7/2012 - ELSE added, Zeroes the output if Downwelling flag not set
!  -- 1/31/21. Version 2.8.3. Use post processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present)
!                  - replaced NSTOKES with ECHT_NSTOKES here 
!                  - initialized LS_LAYER_MSSTS_F 

      ELSE
         DO LUM = 1, N_PPSTREAMS
            !UM = PPSTREAM_MASK(LUM,IBEAM)
            SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,LUM,IBEAM,1:ECHT_NSTOKES,DNIDX) = ZERO
         ENDDO
         LS_LAYER_MSSTS_F(IBEAM,1:ECHT_NSTOKES,1:NLAYERS,1:N_SURFACE_WFS) = ZERO
      ENDIF

!  Mean value output
!  -----------------

!  Version 2.8. Drop DO_QUAD_OUTPUT flag
!  1/31/21. Version 2.8.3. No Changes here

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL VLIDORT_LS_INTEGRATED_OUTPUT ( &
          DO_THERMAL_TRANSONLY, IBEAM,                                    & ! Input flag and index
          NSTOKES, NSTREAMS, NLAYERS, N_DIRECTIONS, N_SURFACE_WFS,        & ! Input numbers
          N_USER_LEVELS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          & ! Input level control
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input partials control
          WHICH_DIRECTIONS, K_REAL, K_COMPLEX, FLUX_MULTIPLIER,           & ! Input Bookkeeping and Flux
          QUAD_WEIGHTS, QUAD_STRMWTS, T_DELT_DISORDS, T_DISORDS_UTUP,     & ! Input Quads/Disord-Transm
          SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, & ! input Homog. Soln/Transm.
          LS_BOA_THTONLY_SOURCE, NCON_SWF, PCON_SWF,                      & ! Input BOA/Linearized constants
          MINT_SURFACEWF, FLUX_SURFACEWF )                                  ! Output
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACEWF_POSTPROCESS_MASTER

!

      SUBROUTINE BOA_SURFACEWF ( &
          DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_DBCORRECTION, DO_LAMBERTIAN_SURFACE, & ! Input flags (general/surface) 
          DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,                & ! Input flags (thermal)
          DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, FOURIER, IBEAM,                     & ! Input flags, indices
          ECHT_NSTOKES, NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,                      & ! input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK,                                                   & ! input Control numbers
          MUELLER_INDEX, FLUXVEC, SURFACE_FACTOR, QUAD_STRMWTS,                         & ! input Bookkeeping + Quads
          ALBEDO, USER_BRDF_F, SL_USERTERM, SURFBB, LS_BRDF_F, LS_TRANS_ATMOS_FINAL,    & ! input Surface stuff
          LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_USER_EMISSIVITY, LS_EMISSIVITY,          & ! input Surface stuff
          T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                        & ! RTE Homog. solutions
          ATMOS_ATTN, STOKES_DOWNSURF, NCON_SWF, PCON_SWF,                              & ! BOA downwelling, Linearized BVP
          LS_BOA_SOURCE, LS_BOA_THTONLY_SOURCE )                                          ! Output

!  4/19/19. Version 2.8.1. Introduce INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL (two directbeam flags)
!  4/19/19. Version 2.8.1. Introduce SL_USERTERM and LS_TRANS_ATMOS_FINAL, to handle waterleaving linearization

!  1/31/21. Version 2.8.3. Two Changes
!     -- Use post-processing mask N_PPSTREAMS, PPSTREAM_MASK, replace LOCAL_UM_START, drop DO_OBSERVATION_GEOMETRY
!     -- Define all BRDF inputs locally for each Fourier

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS, &
                                  MAXSTOKES_SQ, MAXSTREAMS_2, MAXSTRMSTKS, MAXEVALUES, ZERO 

      IMPLICIT NONE
      

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
      
!  4/9/19. Direct surface inclusion flags. Replaces DO_INCLUDE_DIRECTBEAM

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTRF
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTSL

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL

!  Indices

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Post-processing mask

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Bookkeeping, Quads

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR

!  Surface inputs
!   -- 1/31/21. Version 2.8.3. BRDF array is defined locally, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: SURFBB

!  4/9/19, Version 2.8.1. Adjusted water-leaving inputs

      DOUBLE PRECISION, INTENT (IN)  :: SL_USERTERM(MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES)
      DOUBLE PRECISION, INTENT (IN)  :: LS_TRANS_ATMOS_FINAL(MAXBEAMS,MAX_SURFACEWFS)

!  Linearized surface BRDF and emissivity
!   -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )

!  Homogeneous solution

      INTEGER, INTENT (IN) ::          K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

! BOA source contributions

      DOUBLE PRECISION, INTENT (IN) :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )

!  linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  OUTPUT
!  ======

      DOUBLE PRECISION, INTENT (OUT) :: LS_BOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (OUT) :: LS_BOA_THTONLY_SOURCE ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

!  Local variables
!  ---------------

      LOGICAL ::          DO_QTHTONLY
      INTEGER ::          UM, I, J, N, IB, O1, O2, M, OM, Q, LUM
      INTEGER ::          K, KO1, K0, K1, K2
      DOUBLE PRECISION :: INTEGRAND ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      DOUBLE PRECISION :: H1, H2, NXR, PXR, NXR1, NXR2, PXR1

!  Initialise
!  ----------

!  Shorthand

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      M   = FOURIER

!  Special flag. removed QUAD_OUTPUT, Version 2.8
!      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND. ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND. DO_INCLUDE_MVOUTPUT )

!  initialise derivative of BOA source function
!   -- 1/31/21. Version 2.8.3. .Use post-processing mask
!mick fix 1/5/2021 - replaced NSTOKES with ECHT_NSTOKES here 

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SURFACE_WFS
            LS_BOA_SOURCE(Q,UM,1:ECHT_NSTOKES) = ZERO
         ENDDO
      ENDDO

!  Thermal tranmsittance only, special term
!mick fix 1/5/2021 - replaced NSTOKES with ECHT_NSTOKES here

      IF ( DO_QTHTONLY ) THEN
         DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
               LS_BOA_THTONLY_SOURCE(Q,I,1:ECHT_NSTOKES) = ZERO
            ENDDO
         ENDDO
      ENDIF

!  Skip diffuse-field variation for thermal transmittance-only
!   GOTO statement removed, Version 2.8 8/4/16
!      IF ( DO_THERMAL_TRANSONLY .OR. .NOT. DO_USER_STREAMS ) GO TO 599

      IF ( .not.DO_THERMAL_TRANSONLY .AND. DO_USER_STREAMS ) THEN

!  Contribution due to derivatives of BV constants
!  -----------------------------------------------

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

!  start loops

        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

              SUM_R = ZERO
              DO K = 1, K_REAL(N)
                 NXR = NCON_SWF(Q,K,N) * SOLA_XPOS(I,O1,K,N)
                 PXR = PCON_SWF(Q,K,N) * SOLB_XNEG(I,O1,K,N)
                 SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
              ENDDO

!  Complex solutions

              SUM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                       - NCON_SWF(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
                NXR2 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                       + NCON_SWF(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
                PXR1 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                       - PCON_SWF(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
                H1 =  NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
                H2 =  PXR1
                SUM_CR = SUM_CR + H1 + H2
              ENDDO

!  Final result

              INTEGRAND(Q,I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

!  end loops

            ENDDO
          ENDDO
        ENDDO

!  integrated reflectance term
!  ---------------------------

!  integrate Lambertian case, same for all user-streams
!   -- 1/31/21. Version 2.8.3. Use post-processing mask

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          IF ( FOURIER.EQ.0 ) THEN
            O1 = 1 ;  Q = 1
            REFLEC = SURFACE_FACTOR * ALBEDO * SUM(INTEGRAND(Q,1:NSTREAMS,O1))
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               LS_BOA_SOURCE(Q,UM,O1) = REFLEC
            ENDDO
          ENDIF
        ENDIF

!  BRDF case
!   -- 1/31/21. Version 2.8.3. Use post-processing mask
!   -- 1/31/21. Version 2.8.3. BRDF array is defined locally, drop Fourier index "M"

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
           DO Q = 1, N_SURFACE_WFS
              DO LUM = 1, N_PPSTREAMS
                 UM = PPSTREAM_MASK(LUM,IB)
                 DO O1 = 1, NSTOKES
                    REFLEC = ZERO
                    DO J = 1, NSTREAMS
                       S_REFLEC = ZERO
                       DO O2 = 1, NSTOKES
                          OM = MUELLER_INDEX(O1,O2)
                          S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * USER_BRDF_F(OM,UM,J)
                       ENDDO
                       REFLEC = REFLEC + S_REFLEC
                    ENDDO
                    LS_BOA_SOURCE(Q,UM,O1) = REFLEC * SURFACE_FACTOR
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

!  Contributions due to direct variation of kernel parameter
!  ---------------------------------------------------------

!  Lambertian (this is the albedo)
!   -- 1/31/21. Version 2.8.3. Use post-processing mask

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          IF ( FOURIER .EQ. 0 ) THEN
            O1 = 1 ; Q  = 1
!  @@@ Rob Fix. Leave out Albedo ( Un-normalzied WFs )
!           REFLEC = SUM(STOKES_DOWNSURF(1:NSTREAMS,O1)) * ALBEDO * SURFACE_FACTOR
            REFLEC = SUM(STOKES_DOWNSURF(1:NSTREAMS,O1)) * SURFACE_FACTOR
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
            ENDDO
          ENDIF
        ENDIF

!   Non-Lambertian (need derivative of BRDF Fourier term)
!   -- 1/31/21. Version 2.8.3. Use post-processing mask
!   -- 1/31/21. Version 2.8.3. BRDF array is defined locally, drop Fourier index "M"

        IF ( .not.DO_LAMBERTIAN_SURFACE ) THEN
           DO Q = 1, N_SURFACE_WFS
              DO LUM = 1, N_PPSTREAMS
                 UM = PPSTREAM_MASK(LUM,IB)
                 DO O1 = 1, NSTOKES
                    REFLEC = ZERO
                    DO J = 1, NSTREAMS
                       S_REFLEC = ZERO
                       DO O2 = 1, NSTOKES
                          OM = MUELLER_INDEX(O1,O2)
                          S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * LS_USER_BRDF_F(Q,OM,UM,J)
                       ENDDO
                       REFLEC = REFLEC + S_REFLEC
                    ENDDO
                    REFLEC = REFLEC * SURFACE_FACTOR
                    LS_BOA_SOURCE(Q,UM,O1)=LS_BOA_SOURCE(Q,UM,O1) + REFLEC
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

!  End clause for scattering solutions (replaces GOTO 599)

      ENDIF

!  Continuation point for avoiding diffuse field computation
! 599  continue   !  Removed, Version 2.8

!  Thermal tranmsittance and Integrated output (quadrature terms)
!     This is the reflectance term
!   -- 1/31/21. Version 2.8.3. BRDF array is defined locally, drop Fourier index "M"

      IF ( DO_QTHTONLY ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1 ; Q = 1
          REFLEC = SUM(STOKES_DOWNSURF(1:NSTREAMS,O1))
          REFLEC = SURFACE_FACTOR * REFLEC
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(Q,I,O1) = LS_BOA_THTONLY_SOURCE(Q,I,O1)  + REFLEC
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) * LS_BRDF_F(Q,OM,I,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = REFLEC * SURFACE_FACTOR
                LS_BOA_THTONLY_SOURCE(Q,I,O1) = LS_BOA_THTONLY_SOURCE(Q,I,O1)  + REFLEC
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Add linearization term for variation of direct beam reflectance
!  ---------------------------------------------------------------

!   -- 1/31/21. Version 2.8.3. Use post-processing mask
!   -- 1/31/21. Version 2.8.3. BRDF array is defined locally, drop Fourier index "M"

      IF ( DO_USER_STREAMS ) THEN
         IF ( DO_INCLUDE_DIRECTRF.and..not.DO_DBCORRECTION ) THEN
            REFL_ATTN = ATMOS_ATTN(IB)
            IF ( DO_LAMBERTIAN_SURFACE.AND.M.EQ.0) THEN
               O1 = 1 ; Q  = 1
               REFLEC = REFL_ATTN   ! REFLEC = REFL_ATTN * ALBEDO  !   Formerly, was normalized
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
               ENDDO
            ELSE IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
               DO Q = 1, N_SURFACE_WFS
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     DO O1 = 1, NSTOKES
                        REFLEC = ZERO
                        DO O2 = 1, NSTOKES
                           OM = MUELLER_INDEX(O1,O2)
                           REFLEC = REFLEC + LS_USER_BRDF_F_0(Q,OM,UM,IB) * FLUXVEC(O2)
                        ENDDO
                        REFLEC = REFLEC * REFL_ATTN
                        LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF    ! @@@ If (DO_USER_STREAM) clause, Added January 2011

!  Add emissivity variation at user defined angles.
!  ------------------------------------------------

!    Apparenly only present for Fourier zero
!  (expression for emissivity variation follows from Kirchhoff's law)
!  @@@ Rob Fix, add DO_USER_STREAMS to If Block

!   -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_USER_STREAMS ) THEN
         IF ( DO_INCLUDE_SURFEMISS .and. M.eq.0 .and..not.DO_MSMODE_THERMAL ) THEN
            IF ( DO_LAMBERTIAN_SURFACE ) THEN
               O1 = 1 ; Q  = 1
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) - SURFBB
               ENDDO
            ELSE
               DO Q = 1, N_SURFACE_WFS
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     DO O1 = 1, NSTOKES
                        LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) - SURFBB * LS_USER_EMISSIVITY(Q,O1,UM)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF    ! @@@ If (DO_USER_STREAM) clause, Added January 2011

!  Add surface-leaving term linearization if flagged
!   4/9/19. This is new for the water-leaving linearization. Only for Fourier zero
!   -- 1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_USER_STREAMS ) THEN
         IF ( DO_INCLUDE_DIRECTSL .and. M.eq.0 ) THEN
            DO Q = 1, N_SURFACE_WFS
               REFLEC = LS_TRANS_ATMOS_FINAL(IB,Q)
               DO O1 = 1, NSTOKES
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC * SL_USERTERM(UM,IB,O1)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

!  Thermal Transmittance only (For Integrated product)
!  ---------------------------------------------------

!    Surface emissivity term

      IF ( DO_INCLUDE_SURFEMISS .and. DO_QTHTONLY ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1 ; Q = 1
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(Q,I,O1) = LS_BOA_THTONLY_SOURCE(Q,I,O1) -  SURFBB
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                LS_BOA_THTONLY_SOURCE(Q,I,O1) = LS_BOA_THTONLY_SOURCE(Q,I,O1) - SURFBB * LS_EMISSIVITY(Q,O1,I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BOA_SURFACEWF

!

      SUBROUTINE UPUSER_SURFACEWF ( &
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_THERMAL_TRANSONLY, DO_MSSTS, IBEAM, & ! Input Flags
          ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,                    & ! Input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK,                                                      & ! Input Control numbers
          FLUX_MULTIPLIER, K_REAL, K_COMPLEX, UTAU_LEVEL_MASK_UP,                          & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                    & ! input partlayer control
          T_DELT_USERM, T_UTUP_USERM, UHOM_UPDN, UHOM_UPUP, LS_BOA_SOURCE,                 & ! Input User RT Solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, NCON_SWF, PCON_SWF,                  & ! Input Multipliers/Linearized-BVP
          SURFACEWF_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F )                                   ! Output

!  1/31/21. Version 2.8.3. 
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- LInearizations of MSST functions (LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F), now output
!     -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS, MAX_USER_LEVELS, &
                                  MAX_SURFACEWFS, MAXSTRMSTKS, MAXEVALUES, UPIDX, ZERO 

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!  Indices

      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Post-processing mask,replace LOCAL_UM_START and N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Bookkeeping and Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Level/partial output control

      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  user stream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  User homog solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Linearized surface term

      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, UT, IB, LUM

      DOUBLE PRECISION :: LS_CUMUL_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_LAYER_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_FINAL_SOURCE

!  Local indices

      IB  = IBEAM

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops
!                  - replaced NSTOKES with ECHT_NSTOKES here 
!                  - initialized LS_LAYER_MSSTS_F & LS_SURF_MSSTS_F

      !IF ( DO_USER_STREAMS ) THEN
      !   DO UTA = 1, N_USER_LEVELS
      !      DO LUM = 1, N_PPSTREAMS
      !         UM = PPSTREAM_MASK(LUM,IB)
      !         DO Q = 1, N_SURFACE_WFS
      !            SURFACEWF_F(Q,UTA,UM,IB,1:NSTOKES,UPIDX) = ZERO
      !        ENDDO
      !      ENDDO
      !   ENDDO
      !ENDIF

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            DO UTA = 1, N_USER_LEVELS
               DO Q = 1, N_SURFACE_WFS
                  SURFACEWF_F(Q,UTA,LUM,IB,1:ECHT_NSTOKES,UPIDX) = ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      LS_LAYER_MSSTS_F(IB,1:ECHT_NSTOKES,1:NLAYERS,1:N_SURFACE_WFS) = ZERO
      LS_SURF_MSSTS_F (IB,1:ECHT_NSTOKES,1:N_SURFACE_WFS ) = ZERO

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to the BOA sum
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_SURFACE_WFS
               LS_CUMUL_SOURCE(Q,UM,1:NSTOKES) = LS_BOA_SOURCE(Q,UM,1:NSTOKES)
            ENDDO
         ENDDO
      ENDIF

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST surface source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

      IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) then
         DO O1 = 1, NSTOKES
            DO Q = 1, N_SURFACE_WFS
               LS_SURF_MSSTS_F(IB,O1,Q) = FLUX_MULTIPLIER * LS_BOA_SOURCE(Q,UM,O1)
            ENDDO
         ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LS_WHOLELAYER_STERM_UP ( &
              DO_THERMAL_TRANSONLY, IB, N, NSTOKES,      & ! input flags/indices
              N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK, & ! input numbers
              K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,   & ! Input User-homog solutions
              HMULT_1, HMULT_2, NCON_SWF, PCON_SWF,      & ! Input multipliers and lin-constants
              LS_LAYER_SOURCE )                            ! Output

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
               DO O1 = 1, NSTOKES
                  DO Q = 1, N_SURFACE_WFS
                     LS_LAYER_MSSTS_F(IB,O1,N,Q) = FLUX_MULTIPLIER * LS_LAYER_SOURCE(Q,IB,O1)
                  ENDDO
               ENDDO
            ENDIF

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                     LS_CUMUL_SOURCE(Q,UM,O1) = LS_LAYER_SOURCE(Q,UM,O1) + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
            ENDDO

!  End layer loop

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LS_PARTLAYER_STERM_UP ( &
              DO_THERMAL_TRANSONLY, IB, UT, N, NSTOKES,     & ! input flags/indices
              N_SURFACE_WFS,N_PPSTREAMS, PPSTREAM_MASK,     & ! input numbers
              K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,      & ! Input User Homog solutions
              UT_HMULT_UU, UT_HMULT_UD, NCON_SWF, PCON_SWF, & ! Input multipliers and Lin.constants
              LS_LAYER_SOURCE )                               ! Output

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                     LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,UM,O1) + T_UTUP_USERM(UT,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                     SURFACEWF_F(Q,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * LS_FINAL_SOURCE
                  ENDDO
               ENDDO
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!  1/31/21. Version 2.8.3. Use post-processing mask

           IF ( DO_USER_STREAMS ) THEN
              DO LUM = 1, N_PPSTREAMS
                 UM = PPSTREAM_MASK(LUM,IB)
                 DO Q = 1, N_SURFACE_WFS
                    DO O1 = 1, NSTOKES
                       SURFACEWF_F(Q,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,UM,O1)
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF

!  end level output assignation

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE UPUSER_SURFACEWF

!

      SUBROUTINE DNUSER_SURFACEWF ( &
          DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_THERMAL_TRANSONLY, DO_MSSTS, & ! Input Flags
          IBEAM, ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,      & ! Input Control numbers
          N_PPSTREAMS, PPSTREAM_MASK,                                               & ! Input Control numbers
          FLUX_MULTIPLIER, K_REAL, K_COMPLEX, UTAU_LEVEL_MASK_DN,                   & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,             & ! input partlayer control
          T_DELT_USERM, T_UTDN_USERM, UHOM_DNDN, UHOM_DNUP,                         & ! Input RT Solutions
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD, NCON_SWF, PCON_SWF,           & ! Input Multipliers/IntegConsts
          SURFACEWF_F, LS_LAYER_MSSTS_F )                                             ! Output

!  1/31/21. Version 2.8.3. 
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- LInearizations of MSST functions (LS_LAYER_MSSTS_F), now output
!     -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_USER_VZANGLES, MAX_SZANGLES, MAX_DIRECTIONS, MAX_USER_LEVELS, &
                                  MAX_SURFACEWFS, MAXSTRMSTKS, MAXEVALUES, DNIDX, ZERO 

      IMPLICIT NONE

!  Inputs
!  ======

!  flags

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!  Indices

      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Post-processing mask

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Level/partial output control

      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  bookkeeping and Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  user stream transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  User homog solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, UT, IB, LUM

      DOUBLE PRECISION :: LS_CUMUL_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_LAYER_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_TOA_SOURCE   ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LS_FINAL_SOURCE

!  Initialise

      IB  = IBEAM

!  Zero all Fourier component output
!     mick fix 1/8/2013 - changed limit of do loop  !DO UM = 1, LOCAL_UM_START
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops
!                  - replaced NSTOKES with ECHT_NSTOKES here 
!                  - initialized LS_LAYER_MSSTS_F

      !IF ( DO_USER_STREAMS ) THEN
      !   DO UTA = 1, N_USER_LEVELS
      !      DO LUM = 1, N_PPSTREAMS
      !         UM = PPSTREAM_MASK(LUM,IB)
      !         DO Q = 1, N_SURFACE_WFS
      !            SURFACEWF_F(Q,UTA,UM,IB,1:NSTOKES,DNIDX) = ZERO
      !         ENDDO
      !      ENDDO
      !   ENDDO
      !ENDIF

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            DO UTA = 1, N_USER_LEVELS
               DO Q = 1, N_SURFACE_WFS
                  SURFACEWF_F(Q,UTA,LUM,IB,1:ECHT_NSTOKES,DNIDX) = ZERO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      LS_LAYER_MSSTS_F(IB,1:ECHT_NSTOKES,1:NLAYERS,1:N_SURFACE_WFS) = ZERO

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms
!  1/31/21. Version 2.8.3. Use post-processing mask

      IF ( DO_USER_STREAMS ) THEN
         DO Q = 1, N_SURFACE_WFS
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO O1 = 1, NSTOKES
                  LS_TOA_SOURCE(Q,UM,O1)   = ZERO
                  LS_CUMUL_SOURCE(Q,UM,O1) = LS_TOA_SOURCE(Q,UM,O1)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LS_WHOLELAYER_STERM_DN ( &
              DO_THERMAL_TRANSONLY, IB, N, NSTOKES,      & ! input flags/indices
              N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK, & ! input numbers
              K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,   & ! Input User-homog solutions
              HMULT_1, HMULT_2, NCON_SWF, PCON_SWF,      & ! Input multipliers and lin-constants
              LS_LAYER_SOURCE )                            ! Output

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
               DO O1 = 1, NSTOKES
                  DO Q = 1, N_SURFACE_WFS
                     LS_LAYER_MSSTS_F(IB,O1,N,Q) = FLUX_MULTIPLIER * LS_LAYER_SOURCE(Q,IB,O1)
                  ENDDO
               ENDDO
            ENDIF

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                     LS_CUMUL_SOURCE(Q,UM,O1) = LS_LAYER_SOURCE(Q,UM,O1) + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
            ENDDO

!  End layer loop

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LS_PARTLAYER_STERM_DN ( &
              DO_THERMAL_TRANSONLY, IB, UT, N, NSTOKES,     & ! input flags/indices
              N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK,    & ! input numbers
              K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,      & ! Input User Homog solutions
              UT_HMULT_DU, UT_HMULT_DD, NCON_SWF, PCON_SWF, & ! Input multipliers and Lin.constants
              LS_LAYER_SOURCE )                               ! Output

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SURFACE_WFS
                  DO O1 = 1, NSTOKES
                     LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,UM,O1) + T_UTDN_USERM(UT,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                     SURFACEWF_F(Q,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * LS_FINAL_SOURCE
                  ENDDO
               ENDDO
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term
!  1/31/21. Version 2.8.3. Use post-processing mask

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SURFACE_WFS
                   DO O1 = 1, NSTOKES
                     SURFACEWF_F(Q,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,UM,O1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE DNUSER_SURFACEWF

!

      SUBROUTINE LS_WHOLELAYER_STERM_UP ( &
        DO_THERMAL_TRANSONLY, IBEAM, N, NSTOKES,   & ! Input indices
        N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK, & ! input numbers
        K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,   & ! Input User-homog solutions
        HMULT_1, HMULT_2, NCON_SWF, PCON_SWF,      & ! Input multipliers and lin-constants
        LS_LAYERSOURCE )                             ! Output

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_SURFACEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO 

      IMPLICIT NONE

!  Inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Indices

      INTEGER, INTENT (IN) ::          N, IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  User homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, IB, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      KO1 = K_REAL(N) + 1
      IB  = IBEAM

!  Thermal transmittance only
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                 LS_LAYERSOURCE(Q,UM,O1) = ZERO
              ENDDO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and Stokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SWF(Q,K,N)*UHOM_UPDN(UM,O1,K,N)
                  PUXR = PCON_SWF(Q,K,N)*UHOM_UPUP(UM,O1,K,N)
                  H1 =  NUXR * HMULT_2(K,UM,N)
                  H2 =  PUXR * HMULT_1(K,UM,N)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SWF(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                          - NCON_SWF(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SWF(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                          + NCON_SWF(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SWF(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                          - PCON_SWF(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SWF(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                          + PCON_SWF(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
                  H1 = NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N)
                  H2 = PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish
 
      RETURN
      END SUBROUTINE LS_WHOLELAYER_STERM_UP

!

      SUBROUTINE LS_WHOLELAYER_STERM_DN ( &
        DO_THERMAL_TRANSONLY, IBEAM, N, NSTOKES,   & ! Input indices
        N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK, & ! input numbers
        K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,   & ! Input User-homog solutions
        HMULT_1, HMULT_2, NCON_SWF, PCON_SWF,      & ! Input multipliers and lin-constants
        LS_LAYERSOURCE )                             ! Output

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_SURFACEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO 

      IMPLICIT NONE

!  Inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Indices

      INTEGER, INTENT (IN) ::          N, IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  User homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, IB, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Local indices

      KO1 = K_REAL(N) + 1
      IB  = IBEAM

!  Thermal transmittance only
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                 LS_LAYERSOURCE(Q,UM,O1) = ZERO
              ENDDO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and Stokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SWF(Q,K,N)*UHOM_DNDN(UM,O1,K,N)
                  PUXR = PCON_SWF(Q,K,N)*UHOM_DNUP(UM,O1,K,N)
                  H1 =  NUXR * HMULT_1(K,UM,N)
                  H2 =  PUXR * HMULT_2(K,UM,N)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SWF(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                          - NCON_SWF(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SWF(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                          + NCON_SWF(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SWF(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                          - PCON_SWF(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SWF(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                          + PCON_SWF(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
                  H1 = NUXR1 * HMULT_1(K1,UM,N) - NUXR2 * HMULT_1(K2,UM,N)
                  H2 = PUXR1 * HMULT_2(K1,UM,N) - PUXR2 * HMULT_2(K2,UM,N)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_WHOLELAYER_STERM_DN

!

      SUBROUTINE LS_PARTLAYER_STERM_UP ( &
        DO_THERMAL_TRANSONLY, IBEAM, UT, N, NSTOKES,  & ! Input indices
        N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK,    & ! input numbers
        K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,      & ! Input User Homog solutions
        UT_HMULT_UU, UT_HMULT_UD, NCON_SWF, PCON_SWF, & ! Input multipliers and Lin.constants
        LS_LAYERSOURCE )                                ! Output

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_SURFACEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO 

      IMPLICIT NONE

!  Inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Indices

      INTEGER, INTENT (IN) ::          N, IBEAM, UT

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing mask

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  User homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, IB, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Local indices

      KO1 = K_REAL(N) + 1
      IB  = IBEAM

!  Thermal transmittance only
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                 LS_LAYERSOURCE(Q,UM,O1) = ZERO
              ENDDO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loops over user angles, weighting functions and Stokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SWF(Q,K,N) * UHOM_UPDN(UM,O1,K,N)
                  PUXR = PCON_SWF(Q,K,N) * UHOM_UPUP(UM,O1,K,N)
                  H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
                  H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SWF(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                          - NCON_SWF(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SWF(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                          + NCON_SWF(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SWF(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                          - PCON_SWF(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SWF(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                          + PCON_SWF(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
                  H1 =   NUXR1 * UT_HMULT_UD(K1,UM,UT) - NUXR2 * UT_HMULT_UD(K2,UM,UT)
                  H2 =   PUXR1 * UT_HMULT_UU(K1,UM,UT) - PUXR2 * UT_HMULT_UU(K2,UM,UT)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_PARTLAYER_STERM_UP

!

      SUBROUTINE LS_PARTLAYER_STERM_DN ( &
        DO_THERMAL_TRANSONLY, IBEAM, UT, N, NSTOKES,  & ! Input indices
        N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK,    & ! input numbers
        K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,      & ! Input User Homog solutions
        UT_HMULT_DU, UT_HMULT_DD, NCON_SWF, PCON_SWF, & ! Input multipliers and Lin.constants
        LS_LAYERSOURCE )                                ! Output

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_SURFACEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO 

      IMPLICIT NONE

!  Inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  Indices

      INTEGER, INTENT (IN) ::          N, IBEAM, UT

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing mask

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  User homog solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYERSOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, IB, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Local indices

      KO1 = K_REAL(N) + 1
      IB  = IBEAM

!  Thermal transmittance only
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                 LS_LAYERSOURCE(Q,UM,O1) = ZERO
              ENDDO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loops over user angles, weighting functions and Stokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SWF(Q,K,N) * UHOM_DNDN(UM,O1,K,N)
                  PUXR = PCON_SWF(Q,K,N) * UHOM_DNUP(UM,O1,K,N)
                  H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
                  H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SWF(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                          - NCON_SWF(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SWF(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                          + NCON_SWF(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SWF(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                          - PCON_SWF(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SWF(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                          + PCON_SWF(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
                  H1 =   NUXR1 * UT_HMULT_DD(K1,UM,UT) - NUXR2 * UT_HMULT_DD(K2,UM,UT)
                  H2 =   PUXR1 * UT_HMULT_DU(K1,UM,UT) - PUXR2 * UT_HMULT_DU(K2,UM,UT)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_PARTLAYER_STERM_DN

!

      SUBROUTINE VLIDORT_LS_INTEGRATED_OUTPUT ( &
        DO_THERMAL_TRANSONLY, IBEAM,                                    & ! Input flag and index
        NSTOKES, NSTREAMS, NLAYERS, N_DIRECTIONS, N_SURFACE_WFS,        & ! Input numbers
        N_USER_LEVELS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          & ! Input level control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input partials control
        WHICH_DIRECTIONS, K_REAL, K_COMPLEX, FLUX_MULTIPLIER,           & ! Input Bookkeeping and Flux
        QUAD_WEIGHTS, QUAD_STRMWTS, T_DELT_DISORDS, T_DISORDS_UTUP,     & ! Input Quads/Disord-Transm
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, & ! input Homog. Soln/Transm.
        LS_BOA_THTONLY_SOURCE, NCON_SWF, PCON_SWF,                      & ! Input BOA/Linearized constants
        MINT_SURFACEWF, FLUX_SURFACEWF )                                  ! Output

!  Flux output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!  revision for Version 2.8, 8/3/16
!    Dropped the DO_INCLUDE_MVOUTPUT flag, Cleaned up code.

!  1/31/21. Version 2.8.3. No Changes here

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS,              &
                                  MAX_SZANGLES, MAX_DIRECTIONS, MAX_USER_LEVELS, MAX_SURFACEWFS, &
                                  MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, UPIDX, DNIDX, ZERO, HALF, PI2

      IMPLICIT NONE

!  inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  indices

      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  Level-boundary and partial-layer control

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  bookkeeping

      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!  Homogeneous solution variables

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  homog solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Thermal BOA source

      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_THTONLY_SOURCE ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output 
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  Local quadrature output (for debug)

      DOUBLE PRECISION :: QSURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, Q
      DOUBLE PRECISION :: SMI, SFX

!  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADSURFACEWF_OFFGRID_UP ( &
                DO_THERMAL_TRANSONLY, UTA, UT, N,                  & ! Inputs Flag and Indices
                NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,         & ! Inputs Numbers
                FLUX_MULTIPLIER, T_DELT_DISORDS, T_DISORDS_UTUP,   & ! Input transmittances, flux
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! inputs Homog Solution
                T_UTUP_EIGEN, T_UTDN_EIGEN, LS_BOA_THTONLY_SOURCE, & ! inputs Homog Solution
                NCON_SWF, PCON_SWF,                                & ! inputs Homog Solution
                QSURFACEWF_F )                                       ! Output
            ELSE
              CALL QUADSURFACEWF_LEVEL_UP ( &
                DO_THERMAL_TRANSONLY, UTA, NLEVEL,                 & ! Input flag and indices
                NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,         & ! Inputs Numbers
                FLUX_MULTIPLIER, T_DELT_DISORDS,                   & ! Input transmittances, flux
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! inputs Homog Solution
                T_DELT_EIGEN, LS_BOA_THTONLY_SOURCE,               & ! inputs Homog Solution
                NCON_SWF, PCON_SWF,                                & ! inputs Homog Solution
                QSURFACEWF_F )                                       ! Output
            ENDIF
          ENDDO
        ENDIF

!  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL QUADSURFACEWF_OFFGRID_DN ( &
                DO_THERMAL_TRANSONLY, UTA, UT, N,                  & ! Inputs Flag and Indices
                NSTOKES, NSTREAMS, N_SURFACE_WFS, FLUX_MULTIPLIER, & ! Inputs Numbers
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! inputs Homog Solution
                T_UTUP_EIGEN, T_UTDN_EIGEN, NCON_SWF, PCON_SWF,    & ! inputs Homog Solution
                QSURFACEWF_F )                                       ! Output
            ELSE
              CALL QUADSURFACEWF_LEVEL_DN ( &
                DO_THERMAL_TRANSONLY, UTA, NLEVEL,                 & ! Input flag and indices
                NSTOKES, NSTREAMS, N_SURFACE_WFS, FLUX_MULTIPLIER, & ! Input numbes and flux
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! input Homog Solution
                T_DELT_EIGEN, NCON_SWF, PCON_SWF,                  & ! input Homog Solution
                QSURFACEWF_F )                                       ! output
            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux output (diffuse term only)
!   Dropped the DO_INCLUDE_MVOUTPUTflag from Version 2.8, 8/3/16

        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QSURFACEWF_F(Q,UTA,I,O1)
                SFX = SFX + QUAD_STRMWTS(I) * QSURFACEWF_F(Q,UTA,I,O1)
              ENDDO
              MINT_SURFACEWF(Q,UTA,IBEAM,O1,WDIR) = SMI * HALF
              FLUX_SURFACEWF(Q,UTA,IBEAM,O1,WDIR) = SFX * PI2
            ENDDO
          ENDDO
        ENDDO

!  End directions loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LS_INTEGRATED_OUTPUT

!

      SUBROUTINE QUADSURFACEWF_LEVEL_UP ( &
        DO_THERMAL_TRANSONLY, UTA, NL,                     & ! Input flag and indices
        NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,         & ! Inputs Numbers
        FLUX_MULTIPLIER, T_DELT_DISORDS,                   & ! Input transmittances, flux
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! inputs Homog Solution
        T_DELT_EIGEN, LS_BOA_THTONLY_SOURCE,               & ! inputs Homog Solution
        NCON_SWF, PCON_SWF,                                & ! inputs Homog Solution
        QSURFACEWF_F )                                       ! Output

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_LEVELS, &
                                  MAX_SURFACEWFS, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, ZERO

!  inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  indices

      INTEGER, INTENT (IN) ::          UTA, NL

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  Homogeneous solution variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  homog solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  Thermal BOA source

      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_THTONLY_SOURCE ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, O1, K, KO1, K0, K1, K2, LAY, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, THELP
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1

!  For the lowest level
!  ====================

!  Thermal transmittance-only solution

      IF ( NL.EQ.NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            THELP = LS_BOA_THTONLY_SOURCE(Q,I,O1)
            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  --------------------

      IF ( NL .EQ. NLAYERS ) THEN

!  Offset

        KO1 = K_REAL(NL) + 1

!  parameter loop

        DO Q = 1, N_SURFACE_WFS

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(NL)
                NXR = NCON_SWF(Q,K,NL) * SOLA_XPOS(I1,O1,K,NL)
                PXR = PCON_SWF(Q,K,NL) * SOLB_XNEG(I1,O1,K,NL)
                SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(NL)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1 =   NCON_SWF(Q,K1,NL) * SOLA_XPOS(I1,O1,K1,NL) &
                       - NCON_SWF(Q,K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
                NXR2 =   NCON_SWF(Q,K1,NL) * SOLA_XPOS(I1,O1,K2,NL) &
                       + NCON_SWF(Q,K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
                PXR1 =   PCON_SWF(Q,K1,NL) * SOLB_XNEG(I1,O1,K1,NL) &
                       - PCON_SWF(Q,K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
                H1 =   NXR1 * T_DELT_EIGEN(K1,NL) - NXR2 * T_DELT_EIGEN(K2,NL)
                H2 =  PXR1
                SHOM_CR = SHOM_CR + H1 + H2
              ENDDO

!  collect solution

              QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  Thermal transmittance-only solution

      IF ( NL.NE.NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            THELP = LS_BOA_THTONLY_SOURCE(Q,I,O1)
            DO LAY = NLAYERS, N, -1
              THELP = THELP*T_DELT_DISORDS(I,LAY)
            ENDDO
            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  --------------------

      IF ( NL .NE. NLAYERS ) THEN

!  Offset

        KO1 = K_REAL(N) + 1

!  parameter loop

        DO Q = 1, N_SURFACE_WFS

!  Stokes and streams loops

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR = NCON_SWF(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
                PXR = PCON_SWF(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
                SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                       - NCON_SWF(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
                PXR1 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                       - PCON_SWF(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
                PXR2 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                       + PCON_SWF(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
                H1 =   PXR1 * T_DELT_EIGEN(K1,N) - PXR2 * T_DELT_EIGEN(K2,N)
                H2 =  NXR1
                SHOM_CR = SHOM_CR + H1 + H2
              ENDDO

!  collect solution

              QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_LEVEL_UP

!

      SUBROUTINE QUADSURFACEWF_LEVEL_DN ( &
        DO_THERMAL_TRANSONLY, UTA, NL,                     & ! Input flag and indices
        NSTOKES, NSTREAMS, N_SURFACE_WFS, FLUX_MULTIPLIER, & ! Input numbes and flux
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! input Homog Solution
        T_DELT_EIGEN, NCON_SWF, PCON_SWF,                  & ! input Homog Solution
        QSURFACEWF_F )                                       ! output

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_LEVELS, &
                                  MAX_SURFACEWFS, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  indices

      INTEGER, INTENT (IN) ::          UTA, NL

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Homogeneous solution variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  homog solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1

      N = NL

!  Zero level = TOA no solutions
!  =============================

!  Downwelling weighting function at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              QSURFACEWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For other levels in the atmosphere
!  ==================================

      IF ( NL .NE. 0 .and. DO_THERMAL_TRANSONLY) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              QSURFACEWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      endif

!  Scattering solutions
!  --------------------

      IF ( NL .NE. 0 ) THEN

!  Offset

       KO1 = K_REAL(N) + 1

!  parameter loop

       DO Q = 1, N_SURFACE_WFS

!  Stokes and streams loops

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SWF(Q,K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_SWF(Q,K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
             K0 = 2 * K - 2
             K1 = KO1 + K0
             K2 = K1  + 1
             NXR1 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                    - NCON_SWF(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
             NXR2 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                    + NCON_SWF(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
             PXR1 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                    - PCON_SWF(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
             H1 =   NXR1 * T_DELT_EIGEN(K1,N) &
                  - NXR2 * T_DELT_EIGEN(K2,N)
             H2 =  PXR1
             SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
       ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_LEVEL_DN

!

      SUBROUTINE QUADSURFACEWF_OFFGRID_UP ( &
        DO_THERMAL_TRANSONLY, UTA, UT, N,                  & ! Inputs Flag and Indices
        NSTOKES, NSTREAMS, NLAYERS, N_SURFACE_WFS,         & ! Inputs Numbers
        FLUX_MULTIPLIER, T_DELT_DISORDS, T_DISORDS_UTUP,   & ! Input transmittances, flux
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! inputs Homog Solution
        T_UTUP_EIGEN, T_UTDN_EIGEN, LS_BOA_THTONLY_SOURCE, & ! inputs Homog Solution
        NCON_SWF, PCON_SWF,                                & ! inputs Homog Solution
        QSURFACEWF_F )                                       ! Output

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS,  MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_LEVELS, &
                                  MAX_SURFACEWFS, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  indices

      INTEGER, INTENT (IN) ::          UTA, UT, N

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Discrete ordinate transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!  Homogeneous solution variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  homog solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Thermal BOA source

      DOUBLE PRECISION, INTENT (IN) :: LS_BOA_THTONLY_SOURCE ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2, LAY, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, THELP
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For thermal transmittance-only
!  ------------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            THELP = LS_BOA_THTONLY_SOURCE(Q,I,O1)
            DO LAY = NLAYERS, N+1, -1
              THELP = THELP*T_DELT_DISORDS(I,LAY)
            ENDDO
            THELP = THELP*T_DISORDS_UTUP(I,UT)
            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER * THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop

      DO Q = 1, N_SURFACE_WFS

!  stream directions

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SWF(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
              PXR = PCON_SWF(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
              SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) + PXR * T_UTUP_EIGEN(K,UT)
            ENDDO

!  complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in T_UTDN_EIGEN, T_UTUP_EIGEN

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                     - NCON_SWF(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
              NXR2 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I1,O1,K2,N) &
                     + NCON_SWF(Q,K2,N) * SOLA_XPOS(I1,O1,K1,N)
              PXR1 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                     - PCON_SWF(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
              PXR2 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                     + PCON_SWF(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
              H1 =   NXR1 * T_UTDN_EIGEN(K1,UT) - NXR2 * T_UTDN_EIGEN(K2,UT)
              H2 =   PXR1 * T_UTUP_EIGEN(K1,UT) - PXR2 * T_UTUP_EIGEN(K2,UT)
              SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  set result

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  end O1 and I loops, end Q loop

          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_OFFGRID_UP

!

      SUBROUTINE QUADSURFACEWF_OFFGRID_DN ( &
        DO_THERMAL_TRANSONLY, UTA, UT, N,                  & ! Inputs Flag and Indices
        NSTOKES, NSTREAMS, N_SURFACE_WFS, FLUX_MULTIPLIER, & ! Inputs Numbers
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! inputs Homog Solution
        T_UTUP_EIGEN, T_UTDN_EIGEN, NCON_SWF, PCON_SWF,    & ! inputs Homog Solution
        QSURFACEWF_F )                                       ! Output

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS,  MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_LEVELS, &
                                  MAX_SURFACEWFS, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, ZERO

      IMPLICIT NONE

!  inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY

!  indices

      INTEGER, INTENT (IN) ::          UTA, UT, N

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS

!  Flux

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Homogeneous solution variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  homog solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  Linearized integration constants

      DOUBLE PRECISION, INTENT (IN) :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ------

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For thermal transmittance-only, no contribution
!  -----------------------------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        O1 = 1
        DO Q = 1, N_SURFACE_WFS
          QSURFACEWF_F(Q,UTA,1:NSTREAMS,O1) = ZERO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop

      DO Q = 1, N_SURFACE_WFS

!  stream directions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SWF(Q,K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_SWF(Q,K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) + PXR * T_UTUP_EIGEN(K,UT)
            ENDDO

!  complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in T_UTDN_EIGEN, T_UTUP_EIGEN

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                     - NCON_SWF(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
              NXR2 =   NCON_SWF(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                     + NCON_SWF(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
              PXR1 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                     - PCON_SWF(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
              PXR2 =   PCON_SWF(Q,K1,N) * SOLB_XNEG(I,O1,K2,N) &
                     + PCON_SWF(Q,K2,N) * SOLB_XNEG(I,O1,K1,N)
              H1 =   NXR1 * T_UTDN_EIGEN(K1,UT) - NXR2 * T_UTDN_EIGEN(K2,UT)
              H2 =   PXR1 * T_UTUP_EIGEN(K1,UT) - PXR2 * T_UTUP_EIGEN(K2,UT)
              SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  set result

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  end O1 and I loops, end Q loop

          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE QUADSURFACEWF_OFFGRID_DN

!  End module

      END MODULE vlidort_ls_wfsurface_m

