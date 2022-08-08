
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #           VLIDORT_CONVERGE                                  #
! #           VLIDORT_CONVERGE_OBSGEO                           #
! #           VLIDORT_CONVERGE_DOUBLET                          #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. Separate module for the converge routines
!          ==> Uses Input  type structure VLIDORT_SS directly
!          ==> Uses output type structure VLIDORT_Out, filled directly as needed
!          ==> Addition of new VLIDORT_CONVERGE_DOUBLET subroutine
!          ==> Argument lists for the Converge routines streamlined

module vlidort_converge_m

!  Parameter types

   USE VLIDORT_PARS_m, only : fpk, ONE

!  Dependencies

   USE VLIDORT_Outputs_def_m
   USE VLIDORT_Sup_SS_def_m

!  Everything public

public

contains

      SUBROUTINE VLIDORT_CONVERGE &
       ( DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
         DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS,    & ! Input flags
         NSTOKES, NSTREAMS, NLAYERS, N_OUT_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! Input control numbers        
         IBEAM, FOURIER, N_CONVTESTS, VLIDORT_ACCURACY, VZA_OFFSETS, AZMFAC,       & ! Input numbers, convergence
         N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,                          & ! Input bookkeeping
         STOKES_F, MS_CONTRIBS_F, VLIDORT_SS, VLIDORT_Out,                         & ! Input and output fields
         FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                                  ! Output diagnostics

!  convergence test on the Stokes vector intensity
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  1/31/21. Version 2.8.3. Separate module for the converge routines
!          ==> Takes Input  type structure VLIDORT_SS directly, Replaced STOKES_SS/DB (SS inputs)
!          ==> Takes output type structure VLIDORT_Out, replaces dummy array STOKES (final output) 
!          ==> Use statement for parameter file has been streamlined
!          ==> DO_FOCORR_EXTERNAL has been removed
!          ==> Drop LOCAL_UM_START argument. UM loops now start with 1

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAX_USER_LEVELS, MAX_GEOMETRIES,     &
                                 MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_SZANGLES, &
                                 MAX_DIRECTIONS, ZERO, UPIDX, DNIDX, MAXFOURIER

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE

!  Other flags

      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS

!  Convergence control

      LOGICAL, INTENT (IN)          ::  DO_DOUBLE_CONVTEST
      INTEGER, INTENT (IN)          ::  N_CONVTESTS
      DOUBLE PRECISION, INTENT (IN) ::  VLIDORT_ACCURACY

!  Fourier component and beam.

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Control integers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_RELAZMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  Directional control

      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Bookkeeping: Offsets for geometry indexing

      INTEGER, INTENT (IN) ::           N_OUT_STREAMS
      INTEGER, INTENT (IN) ::           VZA_OFFSETS ( MAX_SZANGLES, MAX_USER_VZANGLES )

!  Local number of azimuths and azimuth factors

      INTEGER, INTENT (IN) ::           LOCAL_N_USERAZM
      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

!  Fourier component inputs

      DOUBLE PRECISION, INTENT (IN) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  Single scatter solutions and Direct Beam results

      TYPE(VLIDORT_Sup_SS), INTENT (IN) :: VLIDORT_SS

!  Modified/output variables
!  -------------------------

!  Converged fields

      TYPE(VLIDORT_Main_Outputs), INTENT (INOUT) :: VLIDORT_Out

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  Local variables
!  ---------------

      INTEGER ::          COUNT, COUNT_A, V, N
      INTEGER ::          I, IDIR, UT, UA, W, O
      DOUBLE PRECISION :: TNEW ( MAXSTOKES ), ACCUR, TOLD ( MAXSTOKES ), TAZM ( MAXSTOKES )

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

        IF ( .not. DO_FOCORR_ALONE ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  DO O = 1, NSTOKES
                    VLIDORT_Out%TS_STOKES(UT,V,O,W) = STOKES_F(UT,I,IBEAM,O,W)
                    !CALL TP7G (FOURIER,UT,I,IBEAM,V,O,W,VLIDORT_Out%TS_STOKES,STOKES_F)
                  ENDDO
                ENDDO
               ENDDO
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IBEAM) = FOURIER
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  DO O = 1, NSTOKES
                    VLIDORT_Out%TS_STOKES(UT,V,O,W) = ZERO
                  ENDDO
                ENDDO
               ENDDO
            ENDDO
          ENDDO
        ENDIF

!  TOA contribution functions (only if flagged)

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    VLIDORT_Out%TS_CONTRIBS(V,O,N) = MS_CONTRIBS_F(I,IBEAM,O,N)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    VLIDORT_Out%TS_CONTRIBS(V,O,N) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 2.8    Much simpler condition

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!        IF ( DO_SSFULL .OR. ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

        IF ( DO_FOCORR ) THEN
         DO IDIR = 1, N_DIRECTIONS
          W = WHICH_DIRECTIONS(IDIR)
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IBEAM,I) + UA
                DO O = 1, NSTOKES
                  !CALL TP7E (FOURIER,UT,V,O,W,VLIDORT_Out%TS_STOKES,VLIDORT_SS%TS_STOKES_SS)
                  VLIDORT_Out%TS_STOKES(UT,V,O,W) = &
                          VLIDORT_Out%TS_STOKES(UT,V,O,W) + VLIDORT_SS%TS_STOKES_SS(UT,V,O,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
         ENDDO
        ENDIF

!    CONTRIBS. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag
!     Version 2.8    Changed Logic for SS terms

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_FOCORR ) THEN
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IBEAM,I) + UA
                DO O = 1, NSTOKES
                  DO N = 1, NLAYERS
                    VLIDORT_Out%TS_CONTRIBS(V,O,N) = &
                            VLIDORT_Out%TS_CONTRIBS(V,O,N) + VLIDORT_SS%TS_CONTRIBS_SS(V,O,N)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the direct beam component if flagged (upwelling only)
!       Convergence on STOKES = STOKES + DBEXACT. Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
!        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
        IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS(IBEAM,I) + UA
                DO O = 1, NSTOKES
                  !CALL TP7F (FOURIER,UT,V,O,VLIDORT_Out%TS_STOKES,VLIDORT_SS%TS_STOKES_DB)
                  VLIDORT_Out%TS_STOKES(UT,V,O,UPIDX) = &
                          VLIDORT_Out%TS_STOKES(UT,V,O,UPIDX) + VLIDORT_SS%TS_STOKES_DB(UT,V,O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_OUT_STREAMS
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  DO O = 1, NSTOKES
                    !CALL TP7G (FOURIER,UT,I,IBEAM,V,O,W,VLIDORT_Out%TS_STOKES,STOKES_F)
                    TOLD(O) = VLIDORT_Out%TS_STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IBEAM,UA,O)*STOKES_F(UT,I,IBEAM,O,W)
                    VLIDORT_Out%TS_STOKES(UT,V,O,W) = TOLD(O) + TAZM(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Examine convergence on intensity only
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AN
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO UA = 1, N_USER_RELAZMS
            COUNT_A = 0
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_OUT_STREAMS
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  DO O = 1, NSTOKES
                    !CALL TP7G (FOURIER,UT,I,IBEAM,V,O,W,VLIDORT_Out%TS_STOKES,STOKES_F)
                    TOLD(O) = VLIDORT_Out%TS_STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IBEAM,UA,O)*STOKES_F(UT,I,IBEAM,O,W)
                    TNEW(O) = TOLD(O) + TAZM(O)
                  ENDDO
                  IF ( TAZM(1) .NE. ZERO ) THEN
                    ACCUR     = ABS(TAZM(1)/TNEW(1))
                    IF ( ACCUR .LT. VLIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  DO O = 1, NSTOKES
                    VLIDORT_Out%TS_STOKES(UT,V,O,W) = TNEW(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component

        IF ( DO_TOA_CONTRIBS ) THEN
          DO UA = 1, LOCAL_N_USERAZM
            DO I = 1, N_OUT_STREAMS
              V = VZA_OFFSETS(IBEAM,I) + UA
              DO N = 1, NLAYERS
                DO O = 1, NSTOKES
                  TOLD(O) = VLIDORT_Out%TS_CONTRIBS(V,O,N)
                  TAZM(O) = AZMFAC(I,IBEAM,UA,O)*MS_CONTRIBS_F(I,IBEAM,O,N)
                  VLIDORT_Out%TS_CONTRIBS(V,O,N) = TOLD(O) + TAZM(O)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER
          ENDIF
        ENDIF

!  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CONVERGE

!

      SUBROUTINE VLIDORT_CONVERGE_DOUBLET ( &
         DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING,           & ! Input flags
         DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, & ! Input flags
         NSTOKES, NSTREAMS, N_OUT_STREAMS, N_USER_LEVELS, IBEAM, FOURIER,     & ! Input numbers
         N_CONVTESTS, VLIDORT_ACCURACY, SZA_DOUBLET_OFFSETS,                  & ! Input Bookkeep, Conv.
         AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS,                              & ! Input Bookkeep, Conv.
         STOKES_F, VLIDORT_SS, VLIDORT_Out,                                   & ! Input/Output fields
         FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                             ! Output Convergence

!  convergence test on the Stokes vector intensity

!  1/31/21. Version 2.8.3. BRAND NEW module for the DO_DOUBLET Convergence option
!          ==> Takes Input/Output structures directly, as with the other routines ion this module 
!          ==> Uses SZA_DOUBLET_OFFSET, Tie azimuth output to user-vza indexing (LUA = 1)
!          ==> Drop the Contribution stuff, drop 1, in this subroutine

      USE VLIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS, MAXLAYERS, &
                                 MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS, MAX_USER_LEVELS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION

      LOGICAL, INTENT (IN) ::           DO_DOUBLE_CONVTEST
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

!  Numbers, basic

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES, NSTREAMS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  Numbers, derived

      INTEGER, INTENT (IN) ::           N_CONVTESTS
      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           N_OUT_STREAMS

!  Bookkeeping, Accuracy and azimuth factors

      INTEGER         , INTENT (IN) ::  WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER         , INTENT (IN) ::  SZA_DOUBLET_OFFSETS ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) ::  VLIDORT_ACCURACY
      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

!  Fourier component input

      DOUBLE PRECISION, INTENT (IN) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Single scatter solutions and Direct Beam results

      TYPE(VLIDORT_Sup_SS), intent(in)     :: VLIDORT_SS

!  modified/output variables
!  -------------------------

!  Converged fields

      TYPE(VLIDORT_Main_Outputs), INTENT (INOUT) :: VLIDORT_Out

!  Convergence testing

      INTEGER, INTENT (INOUT) ::          FOURIER_SAVED ( MAX_SZANGLES )
      INTEGER, INTENT (INOUT) ::          TESTCONV
      LOGICAL, INTENT (INOUT) ::          LOCAL_ITERATION

!  Local variables

      INTEGER ::          LUA = 1
      INTEGER ::          COUNT, V, I, IDIR, UT, W, O
      DOUBLE PRECISION :: TNEW ( MAXSTOKES ), ACCUR, TOLD ( MAXSTOKES ), TAZM ( MAXSTOKES )

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

        IF ( .not. DO_FOCORR_ALONE ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                DO O = 1, NSTOKES
                  VLIDORT_Out%TS_STOKES(UT,V,O,W) = STOKES_F(UT,I,IBEAM,O,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IBEAM) = FOURIER
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                DO O = 1, NSTOKES
                  VLIDORT_Out%TS_STOKES(UT,V,O,W) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  1/31/21. Version 2.8.3. Skip TOA Contribs here

!    STOKES. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 2.8    Much simpler condition

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!        IF ( DO_SSFULL .OR. ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

        IF ( DO_FOCORR ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                DO O = 1, NSTOKES
                  VLIDORT_Out%TS_STOKES(UT,V,O,W) = VLIDORT_Out%TS_STOKES(UT,V,O,W) + VLIDORT_SS%TS_STOKES_SS(UT,V,O,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!    STOKES. Add the direct beam component if flagged (upwelling only)
!       Convergence on STOKES = STOKES + DBEXACT. Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms
        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
!        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

        IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_OUT_STREAMS
              V = SZA_DOUBLET_OFFSETS(IBEAM) + I
              DO O = 1, NSTOKES
                VLIDORT_Out%TS_STOKES(UT,V,O,UPIDX) = VLIDORT_Out%TS_STOKES(UT,V,O,UPIDX) + VLIDORT_SS%TS_STOKES_DB(UT,V,O)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_OUT_STREAMS
                  V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                  DO O = 1, NSTOKES
                    TOLD(O) = VLIDORT_Out%TS_STOKES(UT,V,O,W)
                    TAZM(O) = AZMFAC(I,IBEAM,LUA,O)*STOKES_F(UT,I,IBEAM,O,W)
                    VLIDORT_Out%TS_STOKES(UT,V,O,W) = TOLD(O) + TAZM(O)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  Examine convergence on intensity only
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith)
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_OUT_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                DO O = 1, NSTOKES
                  TOLD(O) = VLIDORT_Out%TS_STOKES(UT,V,O,W)
                  TAZM(O) = AZMFAC(I,IBEAM,LUA,O)*STOKES_F(UT,I,IBEAM,O,W)
                  TNEW(O) = TOLD(O) + TAZM(O)
                ENDDO
                IF ( TAZM(1) .NE. ZERO ) THEN
                  ACCUR     = ABS(TAZM(1)/TNEW(1))
                  IF ( ACCUR .LT. VLIDORT_ACCURACY ) THEN
                    COUNT   = COUNT + 1
                  ENDIF
                ELSE
                  COUNT   = COUNT + 1
                ENDIF
                VLIDORT_Out%TS_STOKES(UT,V,1:NSTOKES,W) = TNEW(1:NSTOKES)
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER
          ENDIF
        ENDIF

!  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CONVERGE_DOUBLET

!

      SUBROUTINE VLIDORT_CONVERGE_OBSGEO &
       ( DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
         DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS,    & ! Input flags
         DO_MSSTS, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, IBEAM, FOURIER,      & ! Input numbers
         N_CONVTESTS, VLIDORT_ACCURACY, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS,    & ! Input Bookkeep, Conv.
         STOKES_F, MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F, VLIDORT_SS,         & ! Input/Output fields
         VLIDORT_Out, FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                      ! Output Convergence

!  convergence test on the Stokes vector intensity
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  1/31/21. Version 2.8.3. Separate module for the converge routines
!          ==> This is the OBSGEO Version, which calculates MSST convergence if flagged (DO_MSSTS set)
!          ==> Takes Input  type structure VLIDORT_SS directly, Replaced STOKES_SS/DB (SS inputs)
!          ==> Takes output type structure VLIDORT_Out, replaces dummy array STOKES (final output) 
!          ==> Use statement for parameter file has been streamlined
!          ==> DO_FOCORR_EXTERNAL has been removed

!  module, dimensions and numbers

      USE VLIDORT_pars_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAX_USER_LEVELS, MAX_GEOMETRIES,     &
                                 MAX_USER_VZANGLES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_SZANGLES, &
                                 MAX_DIRECTIONS, ZERO, UPIDX, DNIDX, MAXFOURIER

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags, oldere versions.
!      LOGICAL, INTENT (IN) ::           DO_SSCORR_NADIR
!      LOGICAL, INTENT (IN) ::           DO_SSCORR_OUTGOING
!      LOGICAL, INTENT (IN) ::           DO_SS_EXTERNAL    ! New 15 March 2012
!      LOGICAL, INTENT (IN) ::           DO_SSFULL

!  1/31/21. Version 2.8.3, removed this flag
!      LOGICAL, INTENT (IN) ::           DO_FOCORR_EXTERNAL

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE

!  Other flags

      LOGICAL, INTENT (IN) ::           DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS

!  1/31/21. Version 2.8.3. Add the MSST flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  Convergence control

      LOGICAL, INTENT (IN)          ::  DO_DOUBLE_CONVTEST
      INTEGER, INTENT (IN)          ::  N_CONVTESTS
      DOUBLE PRECISION, INTENT (IN) ::  VLIDORT_ACCURACY

!  Fourier component and beam.

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Control integers

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  Directional control

      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

! azimuth factors

      DOUBLE PRECISION, INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )

!  Fourier component inputs

      DOUBLE PRECISION, INTENT (IN) ::  STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (IN) ::  MS_CONTRIBS_F &
          ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  1/31/21. Version 2.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      DOUBLE PRECISION, INTENT (IN) :: LAYER_MSSTS_F  ( MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION, INTENT (IN) :: SURF_MSSTS_F   ( MAX_SZANGLES, MAXSTOKES  )

!  Single scatter solutions and Direct Beam results

      TYPE(VLIDORT_Sup_SS), intent(in)     :: VLIDORT_SS

!  modified/output variables
!  -------------------------

!  Converged fields

      TYPE(VLIDORT_Main_Outputs), INTENT (INOUT) :: VLIDORT_Out

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

      INTEGER ::          COUNT, N
      INTEGER ::          IDIR, UT, W, O, LUM, LUA
      DOUBLE PRECISION :: TNEW ( MAXSTOKES ), ACCUR, TOLD ( MAXSTOKES ), TAZM ( MAXSTOKES )

!  Local user indices

      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on STOKES = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on STOKES = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED, do not get calculated)

!  Single scatter calculation alone is initialized to zero here (Version

        IF ( .not. DO_FOCORR_ALONE ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W) = STOKES_F(UT,LUM,IBEAM,O,W)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IBEAM) = FOURIER
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  1/31/21. Version 2.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now Downwelling OR Upwelling.

        VLIDORT_Out%TS_SURF_MSSTS (IBEAM,1:NSTOKES) = ZERO 
        IF ( DO_MSSTS ) THEN
          IF ( DO_UPWELLING ) THEN
            VLIDORT_Out%TS_SURF_MSSTS(IBEAM,1:NSTOKES) = SURF_MSSTS_F(IBEAM,1:NSTOKES)
          ENDIF
          DO N = 1, NLAYERS
            VLIDORT_Out%TS_LAYER_MSSTS(IBEAM,1:NSTOKES,N) = LAYER_MSSTS_F(IBEAM,1:NSTOKES,N)
          ENDDO
        ELSE
          VLIDORT_Out%TS_LAYER_MSSTS(IBEAM,1:NSTOKES,1:NLAYERS) = ZERO
        ENDIF

!  TOA contribution functions (only if flagged)

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO O = 1, NSTOKES
              DO N = 1, NLAYERS
                VLIDORT_Out%TS_CONTRIBS(IBEAM,O,N) = MS_CONTRIBS_F(LUM,IBEAM,O,N)
              ENDDO
            ENDDO
          ELSE
            DO O = 1, NSTOKES
              DO N = 1, NLAYERS
                VLIDORT_Out%TS_CONTRIBS(IBEAM,O,N) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag
!     Version 2.8    some renaming of variables

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!        IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

        IF ( DO_FOCORR ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                !CALL TP7E2 (FOURIER,UT,IBEAM,O,W,VLIDORT_Out%TS_STOKES,VLIDORT_SS%TS_STOKES_SS)
                VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W) = &
                        VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W) + VLIDORT_SS%TS_STOKES_SS(UT,IBEAM,O,W)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!    CONTRIBS. Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag
!     Version 2.8    Changed Logic for SS terms

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_FOCORR ) THEN
            DO O = 1, NSTOKES
              DO N = 1, NLAYERS
                VLIDORT_Out%TS_CONTRIBS(IBEAM,O,N) = &
                        VLIDORT_Out%TS_CONTRIBS(IBEAM,O,N) + VLIDORT_SS%TS_CONTRIBS_SS(IBEAM,O,N)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!    STOKES. Add the direct beam component if flagged (upwelling only)
!       Convergence on STOKES = STOKES + DBEXACT
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

        !IF ( (DO_SSFULL.OR.DO_DBCORRECTION).AND.DO_UPWELLING ) THEN
!        IF ( ((DO_SSFULL.OR.DO_DBCORRECTION).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
        IF ( (DO_FOCORR.OR.DO_DBCORRECTION) .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO O = 1, NSTOKES
              !CALL TP7F2 (FOURIER,UT,IBEAM,O,VLIDORT_Out%TS_STOKES,VLIDORT_SS%TS_STOKES_DB)
              VLIDORT_Out%TS_STOKES(UT,IBEAM,O,UPIDX) = &
                      VLIDORT_Out%TS_STOKES(UT,IBEAM,O,UPIDX) + VLIDORT_SS%TS_STOKES_DB(UT,IBEAM,O)
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each geometry, add Fourier component
!     - for direction, user optical depth

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO O = 1, NSTOKES
                  !CALL TP7G2 (FOURIER,UT,UM,IBEAM,O,W,VLIDORT_Out%TS_STOKES,STOKES_F)
                  TOLD(O) = VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W)
                  TAZM(O) = AZMFAC(LUM,IBEAM,LUA,O)*STOKES_F(UT,LUM,IBEAM,O,W)
                  VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W) = TOLD(O) + TAZM(O)
                ENDDO
              ENDDO
            ENDDO

!  Examine convergence on intensity only
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AN
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO O = 1, NSTOKES
                !CALL TP7G2 (FOURIER,UT,UM,IBEAM,O,W,VLIDORT_Out%TS_STOKES,STOKES_F)
                TOLD(O) = VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W)
                TAZM(O) = AZMFAC(LUM,IBEAM,LUA,O)*STOKES_F(UT,LUM,IBEAM,O,W)
                TNEW(O) = TOLD(O) + TAZM(O)
              ENDDO
              IF ( TAZM(1) .NE. ZERO ) THEN
                ACCUR = DABS(TAZM(1)/TNEW(1))
                IF ( ACCUR .LT. VLIDORT_ACCURACY ) THEN
                  COUNT = COUNT + 1
                ENDIF
              ELSE
                COUNT = COUNT + 1
              ENDIF
              DO O = 1, NSTOKES
                VLIDORT_Out%TS_STOKES(UT,IBEAM,O,W) = TNEW(O)
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
              LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF


!  1/31/21. Version 2.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

        IF ( DO_MSSTS ) THEN
          IF ( DO_UPWELLING ) THEN
            DO O = 1, NSTOKES
              TOLD(O) = VLIDORT_Out%TS_SURF_MSSTS(IBEAM,O)
              TAZM(O) = AZMFAC(LUM,IBEAM,LUA,O)*SURF_MSSTS_F(IBEAM,O)
              VLIDORT_Out%TS_SURF_MSSTS(IBEAM,O) = TOLD(O) + TAZM(O)
            ENDDO
          ENDIF
          DO N = 1, NLAYERS
            DO O = 1, NSTOKES
              TOLD(O) = VLIDORT_Out%TS_LAYER_MSSTS(IBEAM,O,N)
              TAZM(O) = AZMFAC(LUM,IBEAM,LUA,O)*LAYER_MSSTS_F(IBEAM,O,N)
              VLIDORT_Out%TS_LAYER_MSSTS(IBEAM,O,N) = TOLD(O) + TAZM(O)
            ENDDO
          ENDDO
        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component

        IF ( DO_TOA_CONTRIBS ) THEN
          DO N = 1, NLAYERS
            DO O = 1, NSTOKES
              TOLD(O) = VLIDORT_Out%TS_CONTRIBS(IBEAM,O,N)
              TAZM(O) = AZMFAC(LUM,IBEAM,LUA,O)*MS_CONTRIBS_F(LUM,IBEAM,O,N)
              VLIDORT_Out%TS_CONTRIBS(IBEAM,O,N) = TOLD(O) + TAZM(O)
            ENDDO
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER
          ENDIF
        ENDIF

!  For all Fourier, keep saveing the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_CONVERGE_OBSGEO

!  End module

END MODULE vlidort_converge_m

