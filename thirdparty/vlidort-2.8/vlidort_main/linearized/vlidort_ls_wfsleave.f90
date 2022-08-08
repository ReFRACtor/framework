
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
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level PUBLIC routines--------------            #
! #            VLIDORT_LSSL_WFS   (master)                 #
! #                                                        #
! #     High level Jacobian routines  ---------            #
! #            LSSL_UPUSER_SURFACEWF                       #
! #            LSSL_DNUSER_SURFACEWF                       #
! #                                                        #
! #     Post-processing at user angles --------            #
! #            LSSL_WHOLELAYER_STERM_UP                    #
! #            LSSL_WHOLELAYER_STERM_DN                    #
! #            LSSL_PARTLAYER_STERM_UP                     #
! #            LSSL_PARTLAYER_STERM_DN                     #
! #                                                        #
! ##########################################################
! #                                                        #
! #     High level Jacobian routine   ---------            #
! #            LSSL_INTEGRATED_OUTPUT (master)             #
! #                                                        #
! #     Post-processing at quad angles --------            #
! #            LSSL_QUADWF_LEVEL_UP                        #
! #            LSSL_QUADWF_LEVEL_DN                        #
! #            LSSL_QUADWF_OFFGRID_UP                      #
! #            LSSL_QUADWF_OFFGRID_DN                      #
! #                                                        #
! ##########################################################

      MODULE vlidort_ls_wfsleave_m

!  Only 1 routine available publicly to the rest of VLIDORT...
!     --- Version 2.8, the VLIDORT_LSSL_DBCORRECTION "Public" subroutine has been removed
!                      Reason: This correction now done as part of the FO code, version 1.5

!  4/9/19. Code streamlined in several ways
!    - removal of the LSSL_Setups routine (trivial anyway). Now done here      
!    - Use of adjusted water-leaving trans-atmos dependency.
      
!  1/31/21. Version 2.8.3. Some changes.
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier
!    -- Introduce DO_MSSTS flag (input) and linearized MSST output (LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F)
!    -- Drop LOCAL_UM_START; all UM loops start with 1
!    -- Use the post processing mask system, don't need separate OBSGEOM calculation
!    -- Extension to all Fourier components possible with water-leaving.

      PUBLIC  :: VLIDORT_LSSL_WFS

      PRIVATE :: LSSL_UPUSER_SURFACEWF,    &
                 LSSL_DNUSER_SURFACEWF,    &
                 LSSL_INTEGRATED_OUTPUT,   &
                 LSSL_WHOLELAYER_STERM_UP, LSSL_WHOLELAYER_STERM_DN, &
                 LSSL_PARTLAYER_STERM_UP,  LSSL_PARTLAYER_STERM_DN,  & 
                 LSSL_QUADWF_LEVEL_UP,     LSSL_QUADWF_LEVEL_DN,     &
                 LSSL_QUADWF_OFFGRID_UP,   LSSL_QUADWF_OFFGRID_DN

      CONTAINS

      SUBROUTINE VLIDORT_LSSL_WFS ( &
          DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, DO_MSSTS,           & ! input flags (general)
          DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_LAMBERTIAN_SURFACE,             & ! input flags (general)
          DO_WATER_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTSL, ECHT_NSTOKES,    & ! input flags (sleave)
          NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_PPSTREAMS, PPSTREAM_MASK,   & ! Input Numbers
          N_SLEAVE_WFS, N_REFLEC_WFS, N_DIRECTIONS, WHICH_DIRECTIONS,              & ! Input Numbers
          NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,              & ! Input Numbers
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, FLUX_MULTIPLIER,                 & ! Input level out
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Input Partial control
          FOURIER, IBEAM, FLUX_FACTOR, DELTA_FACTOR, SURFACE_FACTOR,               & ! Inputs bookkeeping
          ALBEDO, USER_BRDF_F, TRANS_ATMOS_FINAL, SL_QUADTERM, SL_USERTERM,        & ! Inputs surface
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,            & ! Inputs Lin Sleave
          QUAD_WEIGHTS, QUAD_STRMWTS, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Inputs Quads/BVP
          LSSL_TRANS_ATMOS_FINAL, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Inputs Transmittances
          K_REAL, K_COMPLEX, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,             & ! Input RT solutions
          SOLA_XPOS, SOLB_XNEG, UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,        & ! Input  RT solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,    & ! Input Multipliers
          SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F,                          & ! Output (Main)
          MINT_SURFACEWF, FLUX_SURFACEWF, STATUS, MESSAGE, TRACE )                   ! Output (Flux) and exceptions

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, MAX_USER_LEVELS,    &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_DIRECTIONS,  &
                                 MAXMOMENTS, MAX_SURFACEWFS, MAX_SLEAVEWFS, MAXTOTAL, MAXBANDTOTAL,  &
                                 MAXSTREAMS, MAXSTRMSTKS_2, MAXSTRMSTKS, MAXSTREAMS_2, MAXEVALUES,   &
                                 MAXSTOKES_SQ, ZERO, ONE, PI4, UPIDX, DNIDX, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE LAPACK_TOOLS_m, Only : DGBTRS, DGETRS

      IMPLICIT NONE

!  INPUTS
!  ======

!  Flags

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!  Other flags

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT
      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE         ! Added 11  Sep 2012

      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::          DO_WATER_LEAVING   ! 4/9/19. New
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTSL

!  Indices

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  1/31/21. Version 2.8.3. Post-processing mask,replace LOCAL_UM_START and N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Other numbers

      INTEGER, INTENT (IN) ::          NTOTAL, N_SUBDIAG, N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS, NSTKS_NSTRMS_2

!  Output level/partial control

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Surface inputs
!    -- 1/31/21. Version 2.8.3. BRDF array defined locally for each Fourier, ,drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) :: SURFACE_FACTOR, FLUX_FACTOR, DELTA_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

!  Flux and quadrature

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  4/9/19 Surface-leaving contributions, added

      DOUBLE PRECISION, INTENT (IN) :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  4/9/19. Surface-leaving linearized inputs
!    -- 1/31/21. Version 2.8.3. SLEAVE arrays defined locally for each Fourier, ,drop MAXMOMENTS dimension      

      DOUBLE PRECISION, INTENT (IN) :: LSSL_SLTERM_ISOTROPIC( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LSSL_SLTERM_F_0      ( MAX_SLEAVEWFS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LSSL_USER_SLTERM_F_0 ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  The LSSL Jacobian is currently not enabled, set to zero.

      DOUBLE PRECISION, INTENT (IN) :: TRANS_ATMOS_FINAL     ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LSSL_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_SLEAVEWFS )

!  RT Solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  RT Solution transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  BVProblem inputs

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT  ( MAXSTRMSTKS_2 )

!  user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  User solutions

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Solution multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized surface output

      DOUBLE PRECISION, INTENT (INOUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT)  :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Linearized BOA terms

      DOUBLE PRECISION :: LSSL_BOA_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  4/9/19. Local linearized terms defined internally in this subroutine
!mick fix 9/7/2012 - moved MAX_SLEAVEWFS to 1st dimension

      DOUBLE PRECISION :: LSSL_QUADTERM ( MAX_SLEAVEWFS, MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: LSSL_USERTERM ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  Linearized BVP solution

      DOUBLE PRECISION :: COL2_WFSLEAVE  ( MAXTOTAL, MAX_SLEAVEWFS )
      DOUBLE PRECISION :: SCOL2_WFSLEAVE ( MAXSTRMSTKS_2, MAX_SLEAVEWFS )

      DOUBLE PRECISION :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  Other local variables

      INTEGER ::           INFO, NSTOKES_ACTUAL, IB, M, N, Q, O1
      INTEGER ::           K, K0, K1, K2, KO1, C0, CM, I, UM, LUM, IOFF
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION  :: KS
      CHARACTER (LEN=3) :: CI

!  @ @@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
      INTEGER          :: O2, J, OM
      DOUBLE PRECISION :: INTEGRAND ( MAX_SLEAVEWFS, MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: SUM_R, SUM_CR, REFLEC, S_REFLEC, HELP
      DOUBLE PRECISION :: H1, H2, NXR, PXR, NXR1, NXR2, PXR1, SL
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Initialise status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Proxies

      IB  = IBEAM
      M   = FOURIER
      LUM = 1

!  Actual number of Stokes calculations required

      NSTOKES_ACTUAL = NSTOKES
      IF ( DO_SL_ISOTROPIC ) NSTOKES_ACTUAL = 1

!  Nothing to do if M > 0 and Isotropic

!mick fix - turned off
      !IF ( .not. DO_INCLUDE_DIRECTBEAM  ) RETURN
      IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) RETURN

!  4/9/19. Prepare terms. This used to be in the setups subroutine.
!  ----------------------------------------------------------------
      
!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases
!    -- 1/31/21. Version 2.8.3. SLEAVE arrays LSSL_SLTERM_F_0/LSSL_USER_SLTERM_F_0 defined locally, drop FOURIER "M" index      
!    -- 1/31/21. Version 2.8.3. Use postprocessing mask    

       HELP = FLUX_FACTOR / DELTA_FACTOR
       IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
          DO Q = 1, N_SLEAVE_WFS
             DO O1 = 1, NSTOKES
                SL = LSSL_SLTERM_ISOTROPIC(Q,O1,IB) * HELP
                LSSL_QUADTERM(Q,O1,1:NSTREAMS) = SL
                IF ( DO_USER_STREAMS ) THEN
                   DO LUM = 1, N_PPSTREAMS
                      UM = PPSTREAM_MASK(LUM,IB)
                      LSSL_USERTERM(Q,O1,UM) = SL
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO Q = 1, N_SLEAVE_WFS
             DO O1 = 1, NSTOKES
                DO I = 1, NSTREAMS
                   SL = LSSL_SLTERM_F_0(Q,O1,I,IB) * HELP
                   LSSL_QUADTERM(Q,O1,I) = SL
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                   DO LUM = 1, N_PPSTREAMS
                      UM = PPSTREAM_MASK(LUM,IB)
                      LSSL_USERTERM(Q,O1,UM) = LSSL_USER_SLTERM_F_0(Q,O1,UM,IB) * HELP
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDIF
      
!  BVP solution for perturbed integration constants
!  ------------------------------------------------

!  Compute the main column B' where AX = B'
!     Regular BVP Solution --->  NO TELESCOPING HERE

!  initialize. Vitally necessary

      COL2_WFSLEAVE(1:NTOTAL,1:N_SLEAVE_WFS) = ZERO

!  Last layer : Add direct beam variation due to surface leaving 
!    -Waterleaving case uses linearized adjustment.      

      N  = NLAYERS
      C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
      if ( DO_WATER_LEAVING ) THEN
         DO I = 1, NSTREAMS
            IOFF = C0 + NSTOKES*(I-1)
            DO O1 = 1, NSTOKES_ACTUAL
               CM = IOFF + O1 
               DO Q = 1, N_SLEAVE_WFS
                  COL2_WFSLEAVE(CM,Q) = LSSL_QUADTERM(Q,O1,I)  * TRANS_ATMOS_FINAL(IB) &
                                        + SL_QUADTERM(I,IB,O1) * LSSL_TRANS_ATMOS_FINAL(IB,Q)  
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO I = 1, NSTREAMS
            IOFF = C0 + NSTOKES*(I-1)
            DO O1 = 1, NSTOKES_ACTUAL
               CM = IOFF + O1 
               DO Q = 1, N_SLEAVE_WFS
                COL2_WFSLEAVE(CM,Q) = LSSL_QUADTERM(Q,O1,I) 
               ENDDO
            ENDDO
         ENDDO
      ENDIF

         
!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
         DO Q = 1, N_SLEAVE_WFS
            DO N = 1, NTOTAL
               SCOL2_WFSLEAVE(N,Q) = COL2_WFSLEAVE(N,Q)
            ENDDO
         ENDDO
      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SLEAVE_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_WFSLEAVE, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
           WRITE(CI, '(I3)' ) INFO
           MESSAGE = 'argument i illegal value, for i = '//CI
           TRACE   = 'DGBTRS call (multilayer) in VLIDORT_SLEAVE WFS'
           STATUS  = VLIDORT_SERIOUS ; RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, all layer

        DO N = 1, NLAYERS
           C0 = (N-1)*NSTKS_NSTRMS_2
           DO K = 1, K_REAL(N)
              IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
              DO Q = 1, N_SLEAVE_WFS
                 NCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW,Q)
                 PCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW1,Q)
              ENDDO
           ENDDO
           KO1 = K_REAL(N) + 1
           DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
              IROW    = K + K_REAL(N) ; IROW1   = IROW + NSTKS_NSTRMS
              IROW_S  = IROW + K_COMPLEX(N) ; IROW1_S = IROW_S + NSTKS_NSTRMS
              DO Q = 1, N_SLEAVE_WFS
                 NCON_SLEAVE(Q,K1,N) = COL2_WFSLEAVE(C0+IROW,   Q)
                 NCON_SLEAVE(Q,K2,N) = COL2_WFSLEAVE(C0+IROW_S, Q)
                 PCON_SLEAVE(Q,K1,N) = COL2_WFSLEAVE(C0+IROW1,  Q)
                 PCON_SLEAVE(Q,K2,N) = COL2_WFSLEAVE(C0+IROW1_S,Q)
              ENDDO
           ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFSLEAVE

        CALL DGETRS &
           ( 'N', NTOTAL, N_SLEAVE_WFS, SMAT2, MAXSTRMSTKS_2, &
              SIPIVOT, SCOL2_WFSLEAVE, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
           WRITE(CI, '(I3)' ) INFO
           MESSAGE = 'argument i illegal value, for i = '//CI
           TRACE   = 'DGETRS call (Reg. 1 layer) in VLIDORT_SLEAVE_WFS'
           STATUS  = VLIDORT_SERIOUS ; RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, 1 layer

        N = 1
        DO K = 1, K_REAL(N)
           IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
           DO Q = 1, N_SLEAVE_WFS
              NCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW,Q)
              PCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW1,Q)
           ENDDO
        ENDDO
        KO1 = K_REAL(N) + 1
        DO K = 1, K_COMPLEX(N)
           K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
           IROW    = K + K_REAL(N) ; IROW1   = IROW + NSTKS_NSTRMS
           IROW_S  = IROW + K_COMPLEX(N) ; IROW1_S = IROW_S + NSTKS_NSTRMS
           DO Q = 1, N_SLEAVE_WFS
              NCON_SLEAVE(Q,K1,N) = SCOL2_WFSLEAVE(IROW,   Q)
              NCON_SLEAVE(Q,K2,N) = SCOL2_WFSLEAVE(IROW_S, Q)
              PCON_SLEAVE(Q,K1,N) = SCOL2_WFSLEAVE(IROW1,  Q)
              PCON_SLEAVE(Q,K2,N) = SCOL2_WFSLEAVE(IROW1_S,Q)
           ENDDO
        ENDDO

!  end clause

      ENDIF

!  debug------------------------------------------
!      if ( do_debug_write.and.fourier.eq.0 ) then
!        DO N = 1, NLAYERS
!          DO K = 1, K_REAL(N)
!            write(86,'(3i2,1p6e13.5)')FOURIER,N,K,LCON(K,N), MCON(K,N),
!                                      NCON_SLEAVE(1,K,N),PCON_SLEAVE(1,K,N)
!          ENDDO
!        ENDDO
!      ENDIF

!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling weighting functions
!  -----------------------------

      IF ( DO_UPWELLING .and. DO_USER_STREAMS ) THEN

!  Derivative of BOA source function
!        KS = SURFACE_FACTOR. Bug rob fix 9/9/14. Should be 1.0 always
!   4/9/19. Water-leaving must include derivatives of adjusting transmittance

!  1/31/21. Version 2.8.3. Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS     
 
         KS = one
         IF ( DO_INCLUDE_DIRECTSL .AND. M.EQ.0) THEN
            IF ( DO_WATER_LEAVING ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO O1 = 1, NSTOKES_ACTUAL
                     DO Q = 1, N_SLEAVE_WFS
                        LSSL_BOA_SOURCE(Q,UM,O1) = KS * ( LSSL_USERTERM(Q,O1,UM) * TRANS_ATMOS_FINAL(IB) &
                                                       + SL_USERTERM(UM,IB,O1)   * LSSL_TRANS_ATMOS_FINAL(IB,Q) ) 
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  DO O1 = 1, NSTOKES_ACTUAL
                     LSSL_BOA_SOURCE(1:N_SLEAVE_WFS,UM,O1) = KS * LSSL_USERTERM(1:N_SLEAVE_WFS,O1,UM)
                  ENDDO
               ENDDO
            ENDIF
         ELSE
            LSSL_BOA_SOURCE(1:N_SLEAVE_WFS,:,:) = ZERO
         ENDIF

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@

!  Diffuse Term: Contribution due to derivatives of BVP constants
!  -------------------------------------------------------------

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).dI_DOWN(-j)/dS

!  start loops

        N = NLAYERS
        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES_ACTUAL

!  Real homogeneous solutions

              SUM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,O1,K,N)
                PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,O1,K,N)
                SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
              ENDDO

!  Complex solutions

              SUM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                       - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
                NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                       + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
                PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                       - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
                H1 =  NXR1 * T_DELT_EIGEN(K1,N) &
                     -NXR2 * T_DELT_EIGEN(K2,N)
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

!  Lambertian case, same for all user-streams
!  1/31/21. Version 2.8.3. Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
    
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          IF ( FOURIER.EQ.0 ) THEN
            O1 = 1
            DO Q = 1, N_SLEAVE_WFS
               REFLEC = SURFACE_FACTOR * ALBEDO * SUM(INTEGRAND(Q,1:NSTREAMS,O1))
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  LSSL_BOA_SOURCE(Q,UM,O1) = LSSL_BOA_SOURCE(Q,UM,O1) + REFLEC
               ENDDO
              REFLEC = SURFACE_FACTOR * REFLEC * ALBEDO

            ENDDO
          ENDIF
        ENDIF

!  BRDF case
!    -- 1/31/21. Version 2.8.3. Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS     
!    -- 1/31/21. Version 2.8.3. USER_BRDF_F defined locally for each Fourier, drop "M" index

        IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
           DO Q = 1, N_SLEAVE_WFS
              DO LUM = 1, N_PPSTREAMS
                 UM = PPSTREAM_MASK(LUM,IB)
                 DO O1 = 1, NSTOKES_ACTUAL
                    REFLEC = ZERO
                    DO J = 1, NSTREAMS
                      S_REFLEC = ZERO
                      DO O2 = 1, NSTOKES
                         OM = MAXSTOKES*(O1-1) + O2
                         S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * USER_BRDF_F(OM,UM,J)
                      ENDDO
                      REFLEC = REFLEC + S_REFLEC
                   ENDDO
                   LSSL_BOA_SOURCE(Q,UM,O1) = LSSL_BOA_SOURCE(Q,UM,O1) + REFLEC * SURFACE_FACTOR
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!  Upwelling Surface WF contribution for SLEAVE
!  --------------------------------------------

!  1/31/21. Version 2.8.3. 
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F), now output
!mick fix 1/5/2021 - added ECHT_NSTOKES to inputs

        CALL LSSL_UPUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY, DO_MSSTS, IBEAM, FLUX_MULTIPLIER,    & ! Input flag/beam
          ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SLEAVE_WFS, N_REFLEC_WFS,  & ! Input Numbers
          N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_UP,               & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! input partlayer control
          LSSL_BOA_SOURCE, T_DELT_USERM, T_UTUP_USERM,     & ! Input surface/user transmittances
          K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,         & ! Input user solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,      & ! Input Multipliers
          NCON_SLEAVE, PCON_SLEAVE,                        & ! Linearized BVP
          SURFACEWF_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F )   ! Output

!mick temp fix 9/7/2012 - ELSE added, Zeroes the output if Upwelling flag not set
!  -- 1/31/21. Version 2.8.3. Use post processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present)
!                  - replaced NSTOKES with ECHT_NSTOKES here
!                  - initialized LS_LAYER_MSSTS_F & LS_SURF_MSSTS_F 

      ELSE
         DO LUM = 1, N_PPSTREAMS
            !UM = PPSTREAM_MASK(LUM,IBEAM)
            SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS ,&
                      1:N_USER_LEVELS,LUM,IBEAM,1:ECHT_NSTOKES,UPIDX) = ZERO
         ENDDO
         LS_LAYER_MSSTS_F(IBEAM,1:ECHT_NSTOKES,1:NLAYERS,N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS) = ZERO
         LS_SURF_MSSTS_F (IBEAM,1:ECHT_NSTOKES,N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS) = ZERO
      ENDIF

!  Downwelling Surface WF contribution for SLEAVE
!  ----------------------------------------------

      IF ( DO_DNWELLING ) THEN

!  1/31/21. Version 2.8.3. 
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F), now output
!mick fix 1/5/2021 - added ECHT_NSTOKES & NLAYERS to inputs

        CALL LSSL_DNUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY, DO_MSSTS, IBEAM, FLUX_MULTIPLIER,    & ! Input flag/beam
          ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SLEAVE_WFS, N_REFLEC_WFS, & ! Input Numbers
          N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_DN,               & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! input partlayer control
          T_DELT_USERM, T_UTDN_USERM,                      & ! Input user transmittances
          K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,         & ! Input user solutions
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,      & ! Input Multipliers
          NCON_SLEAVE, PCON_SLEAVE,                        & ! Linearized BVP
          SURFACEWF_F, LS_LAYER_MSSTS_F )                    ! OUTPUT

!mick temp fix 9/7/2012 - ELSE added, Zeroes the output if Upwelling flag not set
!  -- 1/31/21. Version 2.8.3. Use post processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present)
!                  - replaced NSTOKES with ECHT_NSTOKES here 
!                  - initialized LS_LAYER_MSSTS_F

      ELSE
         DO LUM = 1, N_PPSTREAMS
            !UM = PPSTREAM_MASK(LUM,IBEAM)
            SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS ,&
                      1:N_USER_LEVELS,LUM,IBEAM,1:ECHT_NSTOKES,DNIDX) = ZERO
         ENDDO
         LS_LAYER_MSSTS_F(IBEAM,1:ECHT_NSTOKES,1:NLAYERS,N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS) = ZERO
      ENDIF

!  mean value output
!  -----------------

!  Verson 2.8, removed Quad output flag
!      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
!  1/31/21. Version 2.8.3. No Changes here

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL LSSL_INTEGRATED_OUTPUT ( &
          IBEAM, FLUX_MULTIPLIER,  NSTOKES, NSTREAMS,      & ! Input numbers
          NLAYERS, N_USER_LEVELS, N_DIRECTIONS,            & ! Input numbers
          WHICH_DIRECTIONS, N_SLEAVE_WFS, N_REFLEC_WFS,    & ! Input numbers
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          & ! Input Control
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,         & ! Input Control
          PARTLAYERS_LAYERIDX, QUAD_WEIGHTS, QUAD_STRMWTS, & ! Input Control/Quads
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,        & ! Input RT Trans,
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,         & ! Input RT Solutions
          NCON_SLEAVE, PCON_SLEAVE,                        & ! Input BVProblem
          MINT_SURFACEWF, FLUX_SURFACEWF )                   ! OUTPUT
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_LSSL_WFS

!

      SUBROUTINE LSSL_UPUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY, DO_MSSTS, IBEAM, FLUX_MULTIPLIER,    & ! Input flag/beam
          ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SLEAVE_WFS, N_REFLEC_WFS,  & ! Input Numbers
          N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_UP,               & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! input partlayer control
          LSSL_BOA_SOURCE, T_DELT_USERM, T_UTUP_USERM,     & ! Input surface/user transmittances
          K_REAL, K_COMPLEX, UHOM_UPDN, UHOM_UPUP,         & ! Input user solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,      & ! Input Multipliers
          NCON_SLEAVE, PCON_SLEAVE,                        & ! Linearized BVP
          SURFACEWF_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F )   ! Output

!  1/31/21. Version 2.8.3.
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F), now output

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, MAX_USER_LEVELS,   &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_DIRECTIONS, &
                                 MAX_SURFACEWFS, MAX_SLEAVEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO, UPIDX

      IMPLICIT NONE

!  Inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!  Beam and flux

      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Numbers

      INTEGER, INTENT (IN) ::          N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS

!  1/31/21. Version 2.8.3. Post-processing mask,replace LOCAL_UM_START and N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  output control

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Surface term

      DOUBLE PRECISION, INTENT (IN) :: LSSL_BOA_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  RT Solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  BVProblem results

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ======

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  local variables
!  ---------------

!  help

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1, UT, IB, LUM
      DOUBLE PRECISION :: LSSL_FINAL_SOURCE

!  Source arrays

      DOUBLE PRECISION :: LSSL_CUMUL_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_LAYER_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  Proxy

      IB = IBEAM

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops
!                  - replaced NSTOKES with ECHT_NSTOKES here
!                  - initialized LS_LAYER_MSSTS_F & LS_SURF_MSSTS_F

      !DO UTA = 1, N_USER_LEVELS
      !   DO LUM = 1, N_PPSTREAMS
      !      UM = PPSTREAM_MASK(LUM,IB)
      !      DO Q = 1, N_SLEAVE_WFS
      !         Q1 = Q + N_REFLEC_WFS
      !         SURFACEWF_F(Q1,UTA,UM,IB,1:NSTOKES,UPIDX) = ZERO
      !      ENDDO
      !   ENDDO
      !ENDDO

      DO LUM = 1, N_PPSTREAMS
         DO UTA = 1, N_USER_LEVELS
            DO Q = 1, N_SLEAVE_WFS
               Q1 = Q + N_REFLEC_WFS
               SURFACEWF_F(Q1,UTA,LUM,IB,1:ECHT_NSTOKES,UPIDX) = ZERO
            ENDDO
         ENDDO
      ENDDO

      LS_LAYER_MSSTS_F(IB,1:ECHT_NSTOKES,1:NLAYERS,N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS) = ZERO
      LS_SURF_MSSTS_F (IB,1:ECHT_NSTOKES,N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS) = ZERO

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to the BOA sum
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SLEAVE_WFS
            LSSL_CUMUL_SOURCE(Q,UM,1:NSTOKES) = LSSL_BOA_SOURCE(Q,UM,1:NSTOKES)
         ENDDO
      ENDDO

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST surface source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

      IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) then
         DO O1 = 1, NSTOKES
            DO Q = 1, N_SLEAVE_WFS
               Q1 = Q + N_REFLEC_WFS
               LS_SURF_MSSTS_F(IB,O1,Q1) = FLUX_MULTIPLIER * LSSL_BOA_SOURCE(Q,UM,O1)
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

         NUT = NLEVEL + 1
         DO N = NSTART, NUT, -1

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LSSL_WHOLELAYER_STERM_UP ( &
              NSTOKES, IB, N, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,              &
              UHOM_UPUP, UHOM_UPDN, HMULT_1, HMULT_2,                   &
              LSSL_LAYER_SOURCE )

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
               DO O1 = 1, NSTOKES
                  DO Q = 1, N_SLEAVE_WFS
                     Q1 = Q + N_REFLEC_WFS
                     LS_LAYER_MSSTS_F(IB,O1,N,Q1) = FLUX_MULTIPLIER * LSSL_LAYER_SOURCE(Q,IB,O1)
                  ENDDO
               ENDDO
            ENDIF

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  DO O1 = 1, NSTOKES
                     LSSL_CUMUL_SOURCE(Q,UM,O1) = LSSL_LAYER_SOURCE(Q,UM,O1) + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO

!  Offgrid output
!  --------------

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LSSL_PARTLAYER_STERM_UP ( &
              NSTOKES, IB, N, UT, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,                  &
              UHOM_UPUP, UHOM_UPDN, UT_HMULT_UU, UT_HMULT_UD,               &
              LSSL_LAYER_SOURCE )

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  DO O1 = 1, NSTOKES
                     LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM,O1) &
                         + T_UTUP_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                     SURFACEWF_F(Q1,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
                  ENDDO
               ENDDO
            ENDDO

!  Ongrid output
!  -------------

!  1/31/21. Version 2.8.3. Use post-processing mask

         ELSE

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  DO O1 = 1, NSTOKES
                     SURFACEWF_F(Q1,UTA,LUM,IB,O1,UPIDX) = FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
            ENDDO

         ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_UPUSER_SURFACEWF

!

      SUBROUTINE LSSL_DNUSER_SURFACEWF ( &
          DO_OBSERVATION_GEOMETRY, DO_MSSTS, IBEAM, FLUX_MULTIPLIER,    & ! Input flag/beam
          ECHT_NSTOKES, NSTOKES, NLAYERS, N_USER_LEVELS, N_SLEAVE_WFS, N_REFLEC_WFS, & ! Input Numbers
          N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_DN,               & ! Input bookkeeping
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! input partlayer control
          T_DELT_USERM, T_UTDN_USERM,                      & ! Input user transmittances
          K_REAL, K_COMPLEX, UHOM_DNDN, UHOM_DNUP,         & ! Input user solutions
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,      & ! Input Multipliers
          NCON_SLEAVE, PCON_SLEAVE,                        & ! Linearized BVP
          SURFACEWF_F, LS_LAYER_MSSTS_F )                    ! Output

!  1/31/21. Version 2.8.3.
!     -- Use post-processing mask, replace LOCAL_UM_START and N_USER_STREAMS
!     -- or the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F), now output

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, MAX_USER_LEVELS,   &
                                 MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_DIRECTIONS, &
                                 MAX_SURFACEWFS, MAX_SLEAVEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO, DNIDX

      IMPLICIT NONE

!  Inputs
!  ======

!  flag

      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY

!  1/31/21. Version 2.8.3. Add flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!  Beam and flux

      INTEGER, INTENT (IN) ::          IBEAM
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

!  Numbers

      INTEGER, INTENT (IN) ::          ECHT_NSTOKES, NSTOKES
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS 

!  1/31/21. Version 2.8.3. Post-processing mask,replace LOCAL_UM_START and N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  output control

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  user transmittances

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  RT Solutions

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Multipliers

      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  BVProblem results

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ======

      DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 1/31/21. Version 2.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  local variables
!  ---------------

!  help

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1, UT, IB, LUM
      DOUBLE PRECISION :: LSSL_FINAL_SOURCE

!  Source arrays

      DOUBLE PRECISION :: LSSL_TOA_SOURCE  ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_LAYER_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION :: LSSL_CUMUL_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  Proxy

      IB  = IBEAM

!  Zero all Fourier component output
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask
!mick fix 1/5/2021 - removed application of PPSTREAM_MASK (not applied to an array with both
!                    UM & IBEAM present); rearranged order of loops
!                  - replaced NSTOKES with ECHT_NSTOKES here 
!                  - initialized LS_LAYER_MSSTS_F

      !DO UTA = 1, N_USER_LEVELS
      !   DO LUM = 1, N_PPSTREAMS
      !      UM = PPSTREAM_MASK(LUM,IB)
      !      DO Q = 1, N_SLEAVE_WFS
      !         Q1 = Q + N_REFLEC_WFS
      !         SURFACEWF_F(Q1,UTA,UM,IB,1:NSTOKES,DNIDX) = ZERO
      !      ENDDO
      !   ENDDO
      !ENDDO

      DO LUM = 1, N_PPSTREAMS
         DO UTA = 1, N_USER_LEVELS
            DO Q = 1, N_SLEAVE_WFS
               Q1 = Q + N_REFLEC_WFS
               SURFACEWF_F(Q1,UTA,LUM,IB,1:ECHT_NSTOKES,DNIDX) = ZERO
            ENDDO
         ENDDO
      ENDDO

      LS_LAYER_MSSTS_F(IB,1:ECHT_NSTOKES,1:NLAYERS,N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS) = ZERO

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SLEAVE_WFS
            LSSL_TOA_SOURCE(Q,UM,1:NSTOKES)   = ZERO
            LSSL_CUMUL_SOURCE(Q,UM,1:NSTOKES) = LSSL_TOA_SOURCE(Q,UM,1:NSTOKES)
         ENDDO
      ENDDO

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

         NUT = NLEVEL
         DO N = NSTART, NUT

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LSSL_WHOLELAYER_STERM_DN ( &
              NSTOKES, IB, N, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,              &
              UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2,                   &
              LSSL_LAYER_SOURCE )

!  1/31/21. Version 2.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
               DO O1 = 1, NSTOKES
                  DO Q = 1, N_SLEAVE_WFS
                     Q1 = Q + N_REFLEC_WFS
                     LS_LAYER_MSSTS_F(IB,O1,N,Q1) = FLUX_MULTIPLIER * LSSL_LAYER_SOURCE(Q,IB,O1)
                  ENDDO
               ENDDO
            ENDIF

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  DO O1 = 1, NSTOKES
                     LSSL_CUMUL_SOURCE(Q,UM,O1) = LSSL_LAYER_SOURCE(Q,UM,O1) + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
            ENDDO

!  end layer loop

         ENDDO

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

!  1/31/21. Version 2.8.3. Use post-processing mask, drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START

            CALL LSSL_PARTLAYER_STERM_DN ( &
              NSTOKES, IB, N, UT, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
              K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,                  &
              UHOM_DNDN, UHOM_DNUP, UT_HMULT_DD, UT_HMULT_DU,               &
              LSSL_LAYER_SOURCE )

!  1/31/21. Version 2.8.3. Use post-processing mask

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  DO O1 = 1, NSTOKES
                     LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM,O1) &
                         + T_UTDN_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM,O1)
                     SURFACEWF_F(Q1,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
                  ENDDO
               ENDDO
            ENDDO

!  Ongrid output
!  -------------

!  User-defined stream output, just set to the cumulative source term
!  1/31/21. Version 2.8.3. Use post-processing mask

         ELSE

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  DO O1 = 1, NSTOKES
                     SURFACEWF_F(Q1,UTA,LUM,IB,O1,DNIDX) = FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM,O1)
                  ENDDO
               ENDDO
            ENDDO

         ENDIF

!  Check for updating the recursion

         IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_DNUSER_SURFACEWF

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_UP ( &
        NSTOKES, IB, N, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,              &
        UHOM_UPUP, UHOM_UPDN, HMULT_1, HMULT_2,                   &
        LSSL_LAYERSOURCE )

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                 MAX_SLEAVEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          IB, N
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Solution inputs

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  outputs

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Offset

      KO1 = K_REAL(N) + 1

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and Stokes

!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SLEAVE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SLEAVE(Q,K,N)*UHOM_UPDN(UM,O1,K,N)
                  PUXR = PCON_SLEAVE(Q,K,N)*UHOM_UPUP(UM,O1,K,N)
                  H1 =  NUXR * HMULT_2(K,UM,N)
                  H2 =  PUXR * HMULT_1(K,UM,N)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                          - NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                          + NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                          - PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                          + PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
                  H1 =   NUXR1 * HMULT_2(K1,UM,N) - NUXR2 * HMULT_2(K2,UM,N)
                  H2 =   PUXR1 * HMULT_1(K1,UM,N) - PUXR2 * HMULT_1(K2,UM,N)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_UP

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_DN ( &
        NSTOKES, IB, N, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,              &
        UHOM_DNDN, UHOM_DNUP, HMULT_1, HMULT_2,                   &
        LSSL_LAYERSOURCE )

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only :  MAXSTOKES, MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
                                  MAX_SLEAVEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          IB, N
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Solution inputs

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  outputs

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  Offset

      KO1 = K_REAL(N) + 1

!  Homogeneous solutions
!  =====================

!  Loop over user angles, weighting functions and STokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SLEAVE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SLEAVE(Q,K,N)*UHOM_DNDN(UM,O1,K,N)
                  PUXR = PCON_SLEAVE(Q,K,N)*UHOM_DNUP(UM,O1,K,N)
                  H1 =  NUXR * HMULT_1(K,UM,N)
                  H2 =  PUXR * HMULT_2(K,UM,N)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                          - NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                          + NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                          - PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                          + PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
                  H1 =   NUXR1 * HMULT_1(K1,UM,N) - NUXR2 * HMULT_1(K2,UM,N)
                  H2 =   PUXR1 * HMULT_2(K1,UM,N) - PUXR2 * HMULT_2(K2,UM,N)    
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_DN

!

      SUBROUTINE LSSL_PARTLAYER_STERM_UP ( &
        NSTOKES, IB, N, UT, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,                  &
        UHOM_UPUP, UHOM_UPDN, UT_HMULT_UU, UT_HMULT_UD,               &
        LSSL_LAYERSOURCE )

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, &
                                 MAXBEAMS, MAX_SLEAVEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          IB, N, UT
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Solution input variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPDN   ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_UPUP   ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2, NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local range

      KO1 = K_REAL(N) + 1

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SLEAVE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SLEAVE(Q,K,N) * UHOM_UPDN(UM,O1,K,N)
                  PUXR = PCON_SLEAVE(Q,K,N) * UHOM_UPUP(UM,O1,K,N)
                  H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
                  H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N) &
                          - NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N) &
                          + NCON_SLEAVE(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N) &
                          - PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N) &
                          + PCON_SLEAVE(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
                  H1 =   NUXR1 * UT_HMULT_UD(K1,UM,UT) - NUXR2 * UT_HMULT_UD(K2,UM,UT)
                  H2 =   PUXR1 * UT_HMULT_UU(K1,UM,UT) - PUXR2 * UT_HMULT_UU(K2,UM,UT)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_PARTLAYER_STERM_UP

!

      SUBROUTINE LSSL_PARTLAYER_STERM_DN ( &
        NSTOKES, IB, N, UT, N_SLEAVE_WFS, N_PPSTREAMS, PPSTREAM_MASK, &
        K_REAL, K_COMPLEX, NCON_SLEAVE, PCON_SLEAVE,                  &
        UHOM_DNDN, UHOM_DNUP, UT_HMULT_DD, UT_HMULT_DU, &
        LSSL_LAYERSOURCE )

!  1/31/21. Version 2.8.3. Use post-processing mask
!    -- drop DO_OBSERVATION_GEOMETRY, N_USER_STREAMS, LOCAL_UM_START
!    -- Add MAXBEAMS to parameter list

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS, &
                                 MAXBEAMS, MAX_SLEAVEWFS, MAXSTRMSTKS, MAXEVALUES, ZERO

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          IB, N, UT
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS

!  1/31/21. Version 2.8.3. Introduce Post-processing masks

      INTEGER, INTENT (IN) ::          N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Solution input variables

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNDN   ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UHOM_DNUP   ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: LSSL_LAYERSOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          UM, O1, K, KO1, K0, K1, K2, Q, LUM
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local range

      KO1 = K_REAL(N) + 1

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes
!    -- 1/31/21. Version 2.8.3. Use Post-processing mask

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_SLEAVE_WFS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  NUXR = NCON_SLEAVE(Q,K,N) * UHOM_DNDN(UM,O1,K,N)
                  PUXR = PCON_SLEAVE(Q,K,N) * UHOM_DNUP(UM,O1,K,N)
                  H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
                  H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
                  SHOM_R = SHOM_R + H1 + H2
               ENDDO

!  Complex homogeneous solutions
!     Rob, 12/21/10, Changed N -> UT in UT_HMULT_DD, UT_HMULT_DU

               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                  NUXR1 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N) &
                          - NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
                  NUXR2 =   NCON_SLEAVE(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N) &
                          + NCON_SLEAVE(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
                  PUXR1 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N) &
                          - PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
                  PUXR2 =   PCON_SLEAVE(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N) &
                          + PCON_SLEAVE(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
                  H1 =   NUXR1 * UT_HMULT_DD(K1,UM,UT) - NUXR2 * UT_HMULT_DD(K2,UM,UT)
                  H2 =   PUXR1 * UT_HMULT_DU(K1,UM,UT) - PUXR2 * UT_HMULT_DU(K2,UM,UT)
                  SHOM_CR = SHOM_CR + H1 + H2
               ENDDO

!  homogeneous contribution

               LSSL_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

!  End loops over Q, O1 and UM

            ENDDO
         ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_PARTLAYER_STERM_DN

!

      SUBROUTINE LSSL_INTEGRATED_OUTPUT ( &
          IBEAM, FLUX_MULTIPLIER,  NSTOKES, NSTREAMS,      & ! Input numbers
          NLAYERS, N_USER_LEVELS, N_DIRECTIONS,            & ! Input numbers
          WHICH_DIRECTIONS, N_SLEAVE_WFS, N_REFLEC_WFS,    & ! Input numbers
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          & ! Input Control
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,         & ! Input Control
          PARTLAYERS_LAYERIDX, QUAD_WEIGHTS, QUAD_STRMWTS, & ! Input Control/Quads
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,        & ! Input RT Trans,
          K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,         & ! Input RT Solutions
          NCON_SLEAVE, PCON_SLEAVE,                        & ! Input BVProblem
          MINT_SURFACEWF, FLUX_SURFACEWF )                   ! OUTPUT

!  1/31/21. Version 2.8.3. No Changes here

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_LEVELS, MAX_SZANGLES,         &
                                 MAXSTREAMS, MAX_USER_STREAMS, MAX_DIRECTIONS, MAX_SURFACEWFS, MAX_SLEAVEWFS, &
                                 MAXSTRMSTKS, MAXSTREAMS_2, MAXEVALUES, UPIDX, DNIDX, ZERO, HALF, PI2

      IMPLICIT NONE

!  INPUT
!  =====

!  Basic munbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Directional and Level output controls

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  stream and multiplier input

      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Solution variables

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output
!  ======

      DOUBLE PRECISION, INTENT (INOUT) :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Local variables
!  ===============

!  Local quadrature output

      DOUBLE PRECISION :: QSLEAVEWF_F ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, Q, Q1
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
              CALL LSSL_QUADWF_OFFGRID_UP ( &
                   UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS, &
                   FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN, & 
                   K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,     &
                   NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
            ELSE
              CALL LSSL_QUADWF_LEVEL_UP  (  &
                 UTA, NLEVEL, N_SLEAVE_WFS,                &
                 NSTOKES, NSTREAMS, NLAYERS,               &
                 FLUX_MULTIPLIER, T_DELT_EIGEN,            &
                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,  &
                 NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
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
              CALL LSSL_QUADWF_OFFGRID_DN ( &
                   UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS, &
                   FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN, & 
                   K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,     &
                   NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
            ELSE
              CALL LSSL_QUADWF_LEVEL_DN  (  &
                 UTA, NLEVEL, N_SLEAVE_WFS,                &
                 NSTOKES, NSTREAMS,                        &
                 FLUX_MULTIPLIER, T_DELT_EIGEN,            &
                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,  &
                 NCON_SLEAVE, PCON_SLEAVE, QSLEAVEWF_F )
            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux output (diffuse term only). If flagged
!  Version 2.8, Drop the MVOUT flag
!        IF ( DO_INCLUDE_MVOUT ) THEN
!          ...............
!        ENDIF

        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO O1 = 1, NSTOKES
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QSLEAVEWF_F(Q,UTA,I,O1)
                SFX = SFX + QUAD_STRMWTS(I) * QSLEAVEWF_F(Q,UTA,I,O1)
              ENDDO
              MINT_SURFACEWF(Q1,UTA,IBEAM,O1,WDIR) = SMI * HALF
              FLUX_SURFACEWF(Q1,UTA,IBEAM,O1,WDIR) = SFX * PI2
            ENDDO
          ENDDO
        ENDDO

!  End directions loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_INTEGRATED_OUTPUT

!

      SUBROUTINE LSSL_QUADWF_LEVEL_UP (                                     &
        UTA, NL, N_SLEAVE_WFS, NSTOKES, NSTREAMS, NLAYERS, FLUX_MULTIPLIER, &
        T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,              &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

!  1/31/21. Version 2.8.3. No Changes here

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_LEVELS,    &
                                 MAX_SLEAVEWFS, MAXSTRMSTKS, MAXSTREAMS_2, MAXEVALUES, ZERO

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS, NLAYERS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1

!  For the lowest level
!  ====================

      IF ( NL .EQ. NLAYERS ) THEN

!  Offset

        KO1 = K_REAL(NL) + 1

!  parameter loop, Stokes and streams loops

        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(NL)
                NXR = NCON_SLEAVE(Q,K,NL) * SOLA_XPOS(I1,O1,K,NL)
                PXR = PCON_SLEAVE(Q,K,NL) * SOLB_XNEG(I1,O1,K,NL)
                SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(NL)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1 =   NCON_SLEAVE(Q,K1,NL) * SOLA_XPOS(I1,O1,K1,NL) &
                       - NCON_SLEAVE(Q,K2,NL) * SOLA_XPOS(I1,O1,K2,NL)
                NXR2 =   NCON_SLEAVE(Q,K1,NL) * SOLA_XPOS(I1,O1,K2,NL) &
                       + NCON_SLEAVE(Q,K2,NL) * SOLA_XPOS(I1,O1,K1,NL)
                PXR1 =   PCON_SLEAVE(Q,K1,NL) * SOLB_XNEG(I1,O1,K1,NL) &
                       - PCON_SLEAVE(Q,K2,NL) * SOLB_XNEG(I1,O1,K2,NL)
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

      IF ( NL .NE. NLAYERS ) THEN

!  Offset

        KO1 = K_REAL(N) + 1

!  parameter loop, Stokes and streams loops

        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO O1 = 1, NSTOKES

!  Real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(N)
                NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
                PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
                SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
              ENDDO

!  Complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(N)
                K0 = 2 * K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                       - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
                PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                       - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
                PXR2 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                       + PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
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
      END SUBROUTINE LSSL_QUADWF_LEVEL_UP

!

      SUBROUTINE LSSL_QUADWF_LEVEL_DN ( &
        UTA, NL, N_SLEAVE_WFS, NSTOKES, NSTREAMS, FLUX_MULTIPLIER, &
        T_DELT_EIGEN, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,     &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

!  1/31/21. Version 2.8.3. No Changes here

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_USER_LEVELS,    &
                                 MAX_SLEAVEWFS, MAXSTRMSTKS, MAXSTREAMS_2, MAXEVALUES, ZERO

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1

!  Downwelling weighting functions at TOA ( or N = 0 ) are zero
!  ============================================================

      IF ( NL .EQ. 0 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              QSURFACEWF_F(Q,UTA,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Downwelling weighting functions at other levels
!  ===============================================

!  Offset

      N = NL
      KO1 = K_REAL(N) + 1

!  parameter loop,  Stokes and streams loops

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

!  Complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                     - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
              NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                     + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
              PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                     - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
              H1 =   NXR1 * T_DELT_EIGEN(K1,N) - NXR2 * T_DELT_EIGEN(K2,N)
              H2 =  PXR1
              SHOM_CR = SHOM_CR + H1 + H2
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I,O1) = FLUX_MULTIPLIER*(SHOM_R+SHOM_CR)

!  Finish loops

          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_DN

!

      SUBROUTINE LSSL_QUADWF_OFFGRID_UP ( &
        UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS,  &
        FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,  &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,      &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

!  1/31/21. Version 2.8.3. No Changes here

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_LEVELS, &
                                 MAX_SLEAVEWFS, MAXSTRMSTKS, MAXSTREAMS_2, MAXEVALUES, ZERO

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, UT, N, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop, stream and Stokes loops

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,O1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,O1,K,N)
              SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) + PXR * T_UTUP_EIGEN(K,UT)
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I1,O1,K1,N) &
                     - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I1,O1,K2,N)
              NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I1,O1,K2,N) &
                     + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I1,O1,K1,N)
              PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K1,N) &
                     - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K2,N)
              PXR2 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I1,O1,K2,N) &
                     + PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I1,O1,K1,N)
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
      END SUBROUTINE LSSL_QUADWF_OFFGRID_UP

!

      SUBROUTINE LSSL_QUADWF_OFFGRID_DN ( &
        UTA, UT, N, N_SLEAVE_WFS, NSTOKES, NSTREAMS,  &
        FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,  &
        K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,      &
        NCON_SLEAVE, PCON_SLEAVE, QSURFACEWF_F )

!  1/31/21. Version 2.8.3. No Changes here

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_LEVELS, &
                                 MAX_SLEAVEWFS, MAXSTRMSTKS, MAXSTREAMS_2, MAXEVALUES, ZERO

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, UT, N, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTOKES, NSTREAMS
      DOUBLE PRECISION, INTENT (IN) :: FLUX_MULTIPLIER

      DOUBLE PRECISION, INTENT (IN) :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTRMSTKS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: QSURFACEWF_F ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      INTEGER ::          I, O1, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION :: NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Offset

      KO1 = K_REAL(N) + 1

!  parameter loop, stream and Stokes loops

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,O1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,O1,K,N)
              SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              NXR1 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K1,N) &
                     - NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
              NXR2 =   NCON_SLEAVE(Q,K1,N) * SOLA_XPOS(I,O1,K2,N) &
                     + NCON_SLEAVE(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
              PXR1 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K1,N) &
                     - PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
              PXR2 =   PCON_SLEAVE(Q,K1,N) * SOLB_XNEG(I,O1,K2,N) &
                     + PCON_SLEAVE(Q,K2,N) * SOLB_XNEG(I,O1,K1,N)
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
      END SUBROUTINE LSSL_QUADWF_OFFGRID_DN

!  Finish module

      END MODULE vlidort_ls_wfsleave_m

