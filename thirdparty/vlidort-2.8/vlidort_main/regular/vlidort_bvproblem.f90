
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
! #            BVP_MATRIXSETUP_MASTER      (master)             #
! #             BVP_MATRIX_INIT                                 #
! #             BVP_SURFACE_SETUP_HOM                           #
! #             BVP_MATRIX_SETUP                                #
! #             BVP_MATRIX_LUD                                  #
! #                                                             #
! #            BVP_SOLUTION_MASTER      (master)                #
! #             BVP_SURFACE_SETUP_BEAM                          #
! #             BVP_COLUMN_SETUP                                #
! #             BVP_BACKSUB                                     #
! #             BVP_BACKSUB_1  (for use with Media-properties)  #
! #                                                             #
! # Telescoped BVP: Subroutines in this Module                  #
! #                                                             #
! #            BVPTEL_MATRIXSETUP_MASTER      (master)          #
! #             BVPTEL_MATRIX_INIT                              #
! #             BVPTEL_MATRIX_INIT_OMP                          #
! #             BVPTEL_SURFACE_SETUP_HOM                        #
! #             BVPTEL_MATRIX_SETUP                             #
! #             BVPTEL_MATRIX_LUD                               #
! #                                                             #
! #            BVPTEL_SOLUTION_MASTER      (master)             #
! #             BVPTEL_SURFACE_SETUP_BEAM                       #
! #             BVPTEL_COLUMN_SETUP                             #
! #             BVPTEL_BACKSUB                                  #
! #                                                             #
! #            BMAT_ROWMASK    (integer function)               #
! #            BTELMAT_ROWMASK (integer function)               #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. Several Changes
!   -- Green's function treatment for Nstokes = 3
!   -- Doublet post-processing, MSST output for sphericity correction
!   -- Bookkeeping changes (BRDF/SLEAVE arrays defined locally, each Fourier

      MODULE vlidort_bvproblem_m

      PRIVATE
      PUBLIC :: BVP_MATRIXSETUP_MASTER,    &
                BVP_SOLUTION_MASTER,       &
                BVPTEL_MATRIXSETUP_MASTER, &
                BVPTEL_SOLUTION_MASTER,    &
                BVP_SURFACE_SETUP_BEAM,    &
                BVP_COLUMN_SETUP, BVP_BACKSUB_1
      
      CONTAINS

      SUBROUTINE BVP_MATRIXSETUP_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,       & ! Input Flags
        FOURIER, NSTOKES, NSTREAMS, NLAYERS,             & ! Input Numbers
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & ! Input Numbers
        NTOTAL, N_SUBDIAG, N_SUPDIAG, SURFACE_FACTOR,    & ! Input Numbers
        QUAD_STRMWTS, ALBEDO, BRDF_F,                    & ! Input Surface stuff
        MUELLER_INDEX, K_REAL, K_COMPLEX,                & ! Input Bookkeeping
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,              & ! Input RTE stuff
        R2_HOMP, R2_HOMM, AXBID_F,                       & ! Output Surface reflection
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                & ! Output BVP Matrices
        STATUS, MESSAGE, TRACE )                           ! Exception handling

!  1/31/21. Version 2.8.3. BRDF input array defined locally, each Fourier.

!  module, dimensions and numbers
!    -- 1/31/21. Version 2.8.3. Remove MAXMOMENTS dimension.

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, MAXSTRMSTKS_2, &
                                 MAXEVALUES, MAXSTOKES_SQ, MAXBANDTOTAL, MAXTOTAL, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS

!  Version 2.8. Re-arrange argument list

      IMPLICIT NONE

!  Flags

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT (IN) ::            FOURIER
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS

!  Derived numbers

      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NTOTAL
      INTEGER, INTENT (IN) ::            N_SUBDIAG
      INTEGER, INTENT (IN) ::            N_SUPDIAG
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2

!  Surface inputs
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally, each Fourier.

      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  solution bookkeeping

      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::            K_REAL   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            K_COMPLEX ( MAXLAYERS )

!  Solutions to homogeneous RTE

      DOUBLE PRECISION, INTENT (IN) ::   T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Surface reflected solution outputs  

      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  BVP Matrices
  
      DOUBLE PRECISION, INTENT (OUT) ::  BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (OUT) ::           IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) ::  SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::           SIPIVOT ( MAXSTRMSTKS_2 )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER :: STATUS_SUB
      INTEGER :: NMIN ( MAXTOTAL )
      INTEGER :: NMAX ( MAXTOTAL )
      INTEGER :: KALL

!  Initialize exception handling

      STATUS  = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Additional setups for the albedo layer
!   Version 2.8, Rearrange the arguments
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally, each Fourier.

      CALL BVP_SURFACE_SETUP_HOM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,              &
        FOURIER, NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS,      &
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,           &
        MUELLER_INDEX, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
        R2_HOMP, R2_HOMM, AXBID_F )

!  initialize compression matrix (Do this for every Fourier component)
!   Version 2.8, Rearrange the arguments

      CALL BVP_MATRIX_INIT ( &
        NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,    &
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        NMIN, NMAX, KALL, BANDMAT2, SMAT2 )

!  set up boundary values matrix in compressed form (the "A" as in AX=B)
!   Version 2.8, Rearrange the arguments

      CALL BVP_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, FOURIER, NSTOKES, NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & 
        NMIN, NMAX, KALL, K_REAL, K_COMPLEX,                     &
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM,    &
        BANDMAT2, SMAT2 )

!  SVD decomposition of compressed boundary values matrix
!   Version 2.8, Rearrange the arguments

      CALL BVP_MATRIX_LUD ( &
        NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT,      &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVP_MATRIXSETUP_MASTER

!

      SUBROUTINE BVP_MATRIX_INIT ( &
        NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,    &
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        NMIN, NMAX, KALL, BANDMAT2, SMAT2 )

!  Initialise the compressed matrix

      USE VLIDORT_PARS_m, Only: MAXBANDTOTAL, MAXTOTAL, MAXSTRMSTKS_2, ZERO

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  BVP matrices

      INTEGER, INTENT (OUT) ::          NMIN ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          NMAX ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          KALL
      DOUBLE PRECISION, INTENT (OUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  Local variables
!  ---------------

!  compression function

!      INTEGER  BMAT_ROWMASK
!      EXTERNAL BMAT_ROWMASK

!  help variables

      INTEGER :: I, J, N3, JS, JF, IS, LF, L, I1, NMJ, NXJ

!  special case

      IF ( NLAYERS .EQ. 1 ) THEN
!  Enhancement # 7, S. Quesada 7/15/15. installed here 2/8/16
        SMAT2(1:NTOTAL, 1:NTOTAL)  = ZERO
        RETURN
      ENDIF

!  compression row indices, NMIN, NMAX and KALL

 !  Enhancement # 8, S. Quesada 7/15/15. installed here 2/8/16
      NMIN(1:N_SUPDIAG+1) = 1

      DO J = N_SUPDIAG + 2, NTOTAL
        NMIN(J) = J - N_SUPDIAG
      ENDDO
      DO J = 1, NTOTAL - N_SUBDIAG
        NMAX(J) = J + N_SUBDIAG
      ENDDO

 !  Enhancement # 9, S. Quesada 7/15/15. installed here 2/8/16
      NMAX(NTOTAL-N_SUBDIAG+1:NTOTAL) = NTOTAL

      KALL = N_SUBDIAG + N_SUPDIAG + 1

!  Former code to assign Array - now replaced by external function
!    Jukaa Kujanpaa, FMI, August 2005
!      DO I = 1, NTOTAL
!        DO J = 1, NTOTAL
!          IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
!            BMAT_ROWMASK(I,J) = KALL + I - J
!          ENDIF
!        ENDDO
!      ENDDO

!  compression matrix zeroing

      N3 = NSTKS_NSTRMS_2 + NSTKS_NSTRMS
      LF = NLAYERS - 2

!  upper band top

      JS = NSTKS_NSTRMS_2 + 1
      JF = N3 - 1
      DO I = 1, NSTKS_NSTRMS
        DO J = JS, JF + I
          NMJ = NMIN(J); NXJ = NMAX(J)
          BANDMAT2(BMAT_ROWMASK(I,J,NMJ,NXJ,KALL),J) = ZERO
        ENDDO
      ENDDO

!  upper band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 - NSTKS_NSTRMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTKS_NSTRMS_2-1
          I1 = I + IS
          DO J = JS, JF + I
            NMJ = NMIN(J); NXJ = NMAX(J)
            BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band

      DO L = 1, LF
        IS = L*NSTKS_NSTRMS_2 + NSTKS_NSTRMS
        JS = IS - N3 + 1
        JF = IS - NSTKS_NSTRMS
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS + I, JF
            NMJ = NMIN(J); NXJ = NMAX(J)
            BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band bottom

      JS = LF * NSTKS_NSTRMS_2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTKS_NSTRMS
      DO I = 1, NSTKS_NSTRMS
        I1 = I + IS
        DO J = JS + I, JF
          NMJ = NMIN(J); NXJ = NMAX(J)
          BANDMAT2(BMAT_ROWMASK(I1,J,NMJ,NXJ,KALL),J) = ZERO
        ENDDO
      ENDDO

!  finish

      RETURN
      END SUBROUTINE BVP_MATRIX_INIT

!

      SUBROUTINE BVP_SURFACE_SETUP_HOM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,              &
        FOURIER, NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS,      &
        SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,           &
        MUELLER_INDEX, K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, &
        R2_HOMP, R2_HOMM, AXBID_F )

!  This is the Lambertian or BRDF surface routine
!     --> Reflected homogeneous solutions

!  Original construction 2004,  with BRDFs July 26, 2010

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO construction with IF block

!  1/31/21. Version 2.8.3. BRDF input array defined locally, each Fourier.
!    -- remove MAXMOMENTS from parameter list

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, MAXSTRMSTKS_2, &
                                  MAXBANDTOTAL, MAXTOTAL, MAXEVALUES, MAXSTOKES_SQ, zero

      IMPLICIT NONE

!  Flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  Numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS

!  Surface inputs
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally, each Fourier. Remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  RTE solution inputs and bookkeeping

      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Surface reflected solution outputs

      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  local variables
!  ---------------

      INTEGER ::         I,J,O1,K,KO1,K0,K1,K2,M,NL
      DOUBLE PRECISION :: KMULT

!  Initialization
!  ==============

!  Zero total reflected contributions
!    Enhancement # 1, S. Quesada 7/15/15. installed here 2/8/16
      R2_HOMP(1:NSTREAMS, 1:NSTOKES, 1:NSTKS_NSTRMS) = ZERO
      R2_HOMM(1:NSTREAMS, 1:NSTOKES, 1:NSTKS_NSTRMS) = ZERO
      AXBID_F(1:NSTREAMS, 1:NSTREAMS, :)             = ZERO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Fourier component

      M = FOURIER

!  Last layer

      NL = NLAYERS

!  Offset

      KO1 = K_REAL(NL) + 1

!  Return if Fourier component not zero for Lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .GT. 0 ) RETURN

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

      IF ( DO_LAMBERTIAN_SURFACE ) THEN

!  Integrate Downward streams of particular solutions

        KMULT = SURFACE_FACTOR * ALBEDO

!  Homogeneous real solutions
!    Enhancement # 2, S. Quesada 7/15/15. installed here 2/8/16
        DO K = 1, K_REAL(NL)
          DO I = 1, NSTREAMS
            R2_HOMP(I,1,K) = sum( QUAD_STRMWTS(1:NSTREAMS) * SOLA_XPOS(1:NSTREAMS,1,K,NL) ) * KMULT
            R2_HOMM(I,1,K) = sum( QUAD_STRMWTS(1:NSTREAMS) * SOLB_XNEG(1:NSTREAMS,1,K,NL) ) * KMULT
          ENDDO
        ENDDO

!  Homogeneous complex solutions
!    Enhancement # 3, S. Quesada 7/15/15. installed here 2/8/16
        DO K = 1, K_COMPLEX(NL)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          DO I = 1, NSTREAMS
            R2_HOMP(I,1,K1) = sum( QUAD_STRMWTS(1:NSTREAMS) * SOLA_XPOS(1:NSTREAMS,1,K1,NL) ) * KMULT
            R2_HOMM(I,1,K1) = sum( QUAD_STRMWTS(1:NSTREAMS) * SOLB_XNEG(1:NSTREAMS,1,K1,NL) ) * KMULT
            R2_HOMP(I,1,K2) = sum( QUAD_STRMWTS(1:NSTREAMS) * SOLA_XPOS(1:NSTREAMS,1,K2,NL) ) * KMULT
            R2_HOMM(I,1,K2) = sum( QUAD_STRMWTS(1:NSTREAMS) * SOLB_XNEG(1:NSTREAMS,1,K2,NL) ) * KMULT
          ENDDO
        ENDDO

!  Total BRDF case
!  ===============

      ELSE

!  help variable
!    -- Enhancement # 4, S. Quesada 7/15/15. installed here 2/8/16
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally, remove M=Fourier index

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            DO J = 1, NSTREAMS
              AXBID_F(I, J, MUELLER_INDEX(O1,1:NSTOKES)) = QUAD_STRMWTS(J) * BRDF_F(MUELLER_INDEX(O1,1:NSTOKES), I, J)
            ENDDO
          ENDDO
        ENDDO

!  homogeneous real solutions
!    Enhancement # 4, S. Quesada 7/15/15. installed here 2/8/16
        DO K = 1, K_REAL(NL)
          DO O1 = 1, NSTOKES
            DO I = 1, NSTREAMS
              R2_HOMP(I,O1,K) = SURFACE_FACTOR * &
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K, NL) )
              R2_HOMM(I,O1,K) = SURFACE_FACTOR * &
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K, NL) )
            ENDDO
          ENDDO
        ENDDO

!  homogeneous complex solutions
!    Enhancement # 5, S. Quesada 7/15/15. installed here 2/8/16
        DO K = 1, K_COMPLEX(NL)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              R2_HOMP(I,O1,K1) = SURFACE_FACTOR * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K1, NL) )
              R2_HOMM(I,O1,K1) = SURFACE_FACTOR * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K1, NL) )
              R2_HOMP(I,O1,K2) = SURFACE_FACTOR * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K2, NL) )
              R2_HOMM(I,O1,K2) = SURFACE_FACTOR * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K2, NL) )
            ENDDO
          ENDDO
        ENDDO

!  End clause BRDF vs. Lambertian

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_SURFACE_SETUP_HOM

!

      SUBROUTINE BVP_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, FOURIER, NSTOKES, NSTREAMS, NLAYERS, &
        NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,        & 
        NMIN, NMAX, KALL, K_REAL, K_COMPLEX,                     &
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM,    &
        BANDMAT2, SMAT2 )

!  Fills up the compressed matrix directly

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO construction with IF block

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, &
                                 MAXSTRMSTKS_2, MAXEVALUES, MAXBANDTOTAL, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Flag

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE

!  Numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NTOTAL
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2

!  Bookkeeping numbers

      INTEGER, INTENT (IN) ::           NMIN ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           NMAX ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           KALL
      INTEGER, INTENT (IN) ::           K_REAL    ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  RTE solution variables

      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  BVP matrices

      DOUBLE PRECISION, INTENT (INOUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  compression function
!  --------------------

!      INTEGER  BMAT_ROWMASK
!      EXTERNAL BMAT_ROWMASK

!  Local variables
!  ---------------

      INTEGER ::          I, I1, N, N1, C0, CE_OFFSET, O1, M
      INTEGER ::          EP, EM, EPC, EPS, EMS
      INTEGER ::          CP, CM, CEP, CEM, CEPS, CEMS, CEM1, CEP1
      INTEGER ::          ROWEP, ROWEM, ROWEPS, ROWEMS, ROWCEP, ROWCEM
      INTEGER ::          ROWCEPS, ROWCEMS, ROWCEP1, ROWCEM1
      INTEGER ::          IR, IROW, KO11, KO1, K0, K1, K2
      DOUBLE PRECISION :: XPNET, XMNET
      DOUBLE PRECISION :: X1R, X2R, X1I, X2I, X1R_H, X1I_H

!  Initialization

      EP = 0

!  Fourier component

      M = FOURIER

!  General case, NLAYERS > 1
!  =========================

      IF ( NLAYERS .GT. 1 ) THEN

!mick fix 2/17/11 - turn back on initialise (fix under test!)
!  Initialise - makes no difference
!    Enhancement # 10, S. Quesada 7/15/15. installed here 2/8/16

        BANDMAT2(1:MAXBANDTOTAL, 1:NTOTAL) = ZERO

!  top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

        N = 1
        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1

!  real solutions

            DO EP = 1, K_REAL(N)
              EM = EP + NSTKS_NSTRMS
              ROWEP = BMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
              ROWEM = BMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
              BANDMAT2(ROWEP,EP) = SOLA_XPOS(I,O1,EP,N)
              BANDMAT2(ROWEM,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
            ENDDO

!  complex solutions (REWORKED)

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              EMS = EPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

              ROWEP  = BMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
              ROWEM  = BMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
              ROWEPS = BMAT_ROWMASK(IROW,EPS,NMIN(EPS),NMAX(EPS),KALL)
              ROWEMS = BMAT_ROWMASK(IROW,EMS,NMIN(EMS),NMAX(EMS),KALL)

              BANDMAT2(ROWEP,EP)   =   X1R
              BANDMAT2(ROWEPS,EPS) = - X1I
              BANDMAT2(ROWEM,EM)   =   X2R
              BANDMAT2(ROWEMS,EMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  intermediate layer boundaries
!  -----------------------------

        C0 = - NSTKS_NSTRMS
        DO N = 2, NLAYERS
          N1 = N - 1
          KO1  = K_REAL(N) + 1
          KO11 = K_REAL(N1) + 1
          C0        = C0 + NSTKS_NSTRMS_2
          CE_OFFSET = C0 - NSTKS_NSTRMS
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW

!  real solutions for layer above the boundary

              DO EP = 1, K_REAL(N1)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS
                ROWCEP = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
                BANDMAT2(ROWCEP,CEP)   = T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
                BANDMAT2(ROWCEM,CEM)   = SOLB_XNEG(I,O1,EP,N1)
              ENDDO

!  real solutions for layer below the boundary
!   ( Note the change of sign !!! )

              DO EP = 1, K_REAL(N)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS
                CEP1 = CEP + NSTKS_NSTRMS_2
                CEM1 = CEM + NSTKS_NSTRMS_2

                ROWCEP1 = BMAT_ROWMASK(CM,CEP1,NMIN(CEP1),NMAX(CEP1),KALL)
                ROWCEM1 = BMAT_ROWMASK(CM,CEM1,NMIN(CEM1),NMAX(CEM1),KALL)

                BANDMAT2(ROWCEP1,CEP1) = - SOLA_XPOS(I,O1,EP,N)
                BANDMAT2(ROWCEM1,CEM1) = - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
              ENDDO

!  complex solutions for layer above boundary

              DO EPC = 1, K_COMPLEX(N1)
                K0 = 2*EPC - 2
                K1 = KO11 + K0
                K2 = K1  + 1

                EP   = K_REAL(N1) + EPC
                EPS  = K_REAL(N1) + K_COMPLEX(N1) + EPC
                CEP  = CE_OFFSET + EP
                CEM  = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS
                CEMS = CEPS + NSTKS_NSTRMS

                X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1) - SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
                X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1) + SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
                X2R = SOLB_XNEG(I,O1,K1,N1)
                X2I = SOLB_XNEG(I,O1,K2,N1)

                ROWCEP  = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDMAT2(ROWCEP ,CEP)  =   X1R
                BANDMAT2(ROWCEPS,CEPS) = - X1I
                BANDMAT2(ROWCEM ,CEM)  =   X2R
                BANDMAT2(ROWCEMS,CEMS) = - X2I
              ENDDO

!  complex solutions for layer below boundary
!   ( Note the change of sign !!! )

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
                CEM  = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
                CEMS = CEPS + NSTKS_NSTRMS

                X1R = SOLA_XPOS(I,O1,K1,N)
                X1I = SOLA_XPOS(I,O1,K2,N)
                X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
                X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

                ROWCEP  = BMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDMAT2(ROWCEP,CEP)   = - X1R
                BANDMAT2(ROWCEPS,CEPS) =   X1I
                BANDMAT2(ROWCEM,CEM)   = - X2R
                BANDMAT2(ROWCEMS,CEMS) =   X2I
              ENDDO

            ENDDO
          ENDDO
        ENDDO

!  bottom BC (with Surface contributions if flagged)
!  -------------------------------------------------

        N   = NLAYERS
        KO1 = K_REAL(N) + 1
        C0  = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS

!  With surface terms
!  ------------------

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = C0 + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS

                XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
                XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)

                ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

                BANDMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                BANDMAT2(ROWCEM,CEM) = XMNET
              ENDDO

!  Complex solutions

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS
                CEMS = CEPS + NSTKS_NSTRMS

                X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
                X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
                X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
                X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)

                X1R = X1R_H * T_DELT_EIGEN(K1,N) - X1I_H * T_DELT_EIGEN(K2,N)
                X1I = X1R_H * T_DELT_EIGEN(K2,N) + X1I_H * T_DELT_EIGEN(K1,N)

                ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDMAT2(ROWCEP,CEP)   =   X1R
                BANDMAT2(ROWCEPS,CEPS) = - X1I
                BANDMAT2(ROWCEM,CEM)   =   X2R
                BANDMAT2(ROWCEMS,CEMS) = - X2I

              ENDDO

            ENDDO
          ENDDO

!  No surface terms
!  ----------------

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = C0 + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS

                XPNET = SOLA_XPOS(I1,O1,EP,N)
                XMNET = SOLB_XNEG(I1,O1,EP,N)

                ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

                BANDMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                BANDMAT2(ROWCEM,CEM) = XMNET
              ENDDO

!  Complex solutions
!    ----- Deep bug found for X1I, 17 December 2005.

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                CEP  = CE_OFFSET + EP
                CEM  = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS
                CEMS = CEPS + NSTKS_NSTRMS

                X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N) - SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
                X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N) + SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
                X2R = SOLB_XNEG(I1,O1,K1,N)
                X2I = SOLB_XNEG(I1,O1,K2,N)

                ROWCEP  = BMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDMAT2(ROWCEP,CEP)   =   X1R
                BANDMAT2(ROWCEPS,CEPS) = - X1I
                BANDMAT2(ROWCEM,CEM)   =   X2R
                BANDMAT2(ROWCEMS,CEMS) = - X2I
             ENDDO

            ENDDO
          ENDDO

!  End surface clause

        ENDIF

!  special case for 1 layer only
!  =============================

      ELSE

!  top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

        N = 1
        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1

!  real solutions

            DO EP = 1, K_REAL(N)
              EM = EP + NSTKS_NSTRMS
              SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
              SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
            ENDDO

!  REWORKED complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              EMS = EPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

              SMAT2(IROW,EP)  =   X1R
              SMAT2(IROW,EPS) = - X1I
              SMAT2(IROW,EM)  =   X2R
              SMAT2(IROW,EMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  bottom BC (with Surface contributions if flagged)
!  -------------------------------------------------

!   With surface terms

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = NSTKS_NSTRMS + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = EP
                CEM = CEP + NSTKS_NSTRMS
                XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
                XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
                SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                SMAT2(CP,CEM) = XMNET
              ENDDO

!  Complex solutions

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                EM  = EP  + NSTKS_NSTRMS
                CEP = EP
                CEM = EM
                CEPS = EPS
                CEMS = CEPS + NSTKS_NSTRMS
                X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
                X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
                X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
                X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
                X1R = X1R_H * T_DELT_EIGEN(K1,N) - X1I_H * T_DELT_EIGEN(K2,N)
                X1I = X1R_H * T_DELT_EIGEN(K2,N) + X1I_H * T_DELT_EIGEN(K1,N)
                SMAT2(CP,CEP)  =   X1R
                SMAT2(CP,CEPS) = - X1I
                SMAT2(CP,CEM)  =   X2R
                SMAT2(CP,CEMS) = - X2I
              ENDDO

            ENDDO
          ENDDO

!  No Surface terms
!  ----------------

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = NSTKS_NSTRMS + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = EP
                CEM = CEP + NSTKS_NSTRMS
                XPNET = SOLA_XPOS(I1,O1,EP,N)
                XMNET = SOLB_XNEG(I1,O1,EP,N)
                SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                SMAT2(CP,CEM) = XMNET
              ENDDO

!  REWORKED Complex solutions

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                EM  = EP  + NSTKS_NSTRMS
                CEP = EP
                CEM = EM
                CEPS = EPS
                CEMS = CEPS + NSTKS_NSTRMS
                X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
                X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
                X2R = SOLB_XNEG(I1,O1,K1,N)
                X2I = SOLB_XNEG(I1,O1,K2,N)
                SMAT2(CP,CEP)  =   X1R
                SMAT2(CP,CEPS) = - X1I
                SMAT2(CP,CEM)  =   X2R
                SMAT2(CP,CEMS) = - X2I
              ENDDO

            ENDDO
          ENDDO

!  End clause for surface contributions

        ENDIF

!  End NLAYERS > 1 vs NLAYLERS = 1 clause

      ENDIF

!  normal return and finish

      RETURN
      END SUBROUTINE BVP_MATRIX_SETUP

!

      SUBROUTINE BVP_MATRIX_LUD ( &
        NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT,      &
        STATUS, MESSAGE, TRACE )

!  Solves the boundary value problem.

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.

      USE VLIDORT_PARS_m, only : MAXBANDTOTAL, MAXTOTAL, MAXSTRMSTKS_2, VLIDORT_SUCCESS, VLIDORT_SERIOUS
      USE LAPACK_TOOLS_m, only : DGBTRF, DGETRF

      IMPLICIT NONE

!  Numbers

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG

!  BVP Matrices

      DOUBLE PRECISION, INTENT (INOUT) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::            IPIVOT ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::            SIPIVOT ( MAXSTRMSTKS_2 )

!  Exception handling

      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI

!  Intialize Exception handling

      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  LUD the BVP matrix: With compression (multilayers)
!  --------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF &
           ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
             BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_LUD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_LUD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  LUD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF &
           ( NTOTAL, NTOTAL, SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call (nlayers=1) in BVP_MATRIX_LUD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_MATRIX_LUD

!

      SUBROUTINE BVP_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS,       & ! Input Surface Flags
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTBEAM,            & ! Input Surface Flags
        DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IBEAM,                & ! Input illumination flags, indices
        DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                 & ! Input Water-leaving flags
        NSTOKES, NSTREAMS, NLAYERS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,  & ! Numbers
        NTOTAL, N_SUBDIAG, N_SUPDIAG, K_REAL, K_COMPLEX, MUELLER_INDEX,        & ! Numbers, bookkeeping
        TF_MAXITER, TF_CRITERION, FLUX_FACTOR, FLUXVEC, SURFACE_FACTOR,        & ! Input factors/WL control
        QUAD_STRMWTS, ALBEDO, AXBID_F, SURFBB, EMISSIVITY,                     & ! Input surface terms
        SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX,                        & ! Input Sleave/illumination
        BEAM_BOATRANS, LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,   & ! Input Direct-flux
        RF_DIRECT_BEAM, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WUPPER, WLOWER,    & ! Input RTE solutions
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                      & ! Input BVP matrices/pivots
        COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM, R2_BEAM, LCON, MCON,      & ! Modified input/output
        STATUS, MESSAGE, TRACE_1, TRACE_2 )                                      ! Exception handling

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.

!  Major overhaul for 2.8.1, 4/9/19.
!  ---------------------------------

!  BVP_ADJUST_BACKSUB subroutine.
!    - This is an automatic adjustment for self-consistency with Water-leaving situations.
!    - Only done for Fourier = 0, and Water-leaving flag is on, with TF adjustment control
!    - Add Surface-leaving control and terms SLTERM_ISOTROPIC, SLTERM_F_0
!    - Add Direct Flux inputs LOCAL_CSZA, BEAM_CUTOFF, TRANS_MUBAR, INITIAL_TRANS
!    - Additional Output TRANS_ATMOS_FINAL, Modified Direct-radiance outputs

!  Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
!      --- Some rearrangement of subroutine arguments        

!  1/31/21. Version 2.8.3. BRDF/SLEAVE input arrays defined locally, each Fourier
!    -- Remove MAXMOMENTS from parameter list

      USE VLIDORT_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAX_SZANGLES,           &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE

      IMPLICIT NONE

!  Surface flags

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFEMISS
      LOGICAL, intent(in)  ::            DO_SURFACE_LEAVING

!  Illumination flags. New 3/23/19

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_TOAFLUX
      LOGICAL, INTENT (IN) ::            DO_INCLUDE_BOAFLUX

!  Fourier component and beam number

      INTEGER, INTENT (IN) ::            FOURIER
      INTEGER, INTENT (IN) ::            IBEAM

!  4/9/19. Additional Water-leaving control

      LOGICAL  , intent(in)  :: DO_WATER_LEAVING
      LOGICAL  , intent(in)  :: DO_EXTERNAL_WLEAVE
      LOGICAL  , intent(in)  :: DO_TF_ITERATION
      INTEGER  , INTENT (IN) :: TF_MAXITER
      DOUBLE PRECISION, INTENT (IN) :: TF_CRITERION

!  Basic control numbers

      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS

      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2

      INTEGER, INTENT (IN) ::            NTOTAL
      INTEGER, INTENT (IN) ::            N_SUBDIAG
      INTEGER, INTENT (IN) ::            N_SUPDIAG

!  Bookkeeping

      INTEGER, INTENT (IN) ::            K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            K_COMPLEX ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  Version 2p7 input, 2/19/14 (bookkeeping)

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR, FLUXVEC(MAXSTOKES)

!  Surface and TOA/BOA Fluxes (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) ::  BOAFLUX
      DOUBLE PRECISION, INTENT (IN) ::  TOAFLUX
     
!  Surface inputs

      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::   SURFBB

!  4/9/19. Surface-leaving quantities ==> original supplement results (not flux-adjusted)
!  1/31/21. Version 2.8.3. SLEAVE input array defined locally, each Fourier, remove MAXMOMENTS

      LOGICAL         , INTENT (IN) :: DO_SL_ISOTROPIC
      DOUBLE PRECISION, INTENT (IN) :: SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SLTERM_F_0  ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )

!  4/9/19. Reflected Direct beam, Always input here

      DOUBLE PRECISION, intent(in)  :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  Direct-flux inputs (Optical properties, Average-secant parameterization)

      DOUBLE PRECISION, INTENT (IN) :: BEAM_BOATRANS ( MAXBEAMS )
      INTEGER         , INTENT (IN) :: BEAM_CUTOFF   ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR     ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA       ( 0:MAXLAYERS, MAXBEAMS )
      
!  RTE solution inputs

      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  RTE solution inputs

      DOUBLE PRECISION, INTENT (IN) ::   WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  BVP matrix inputs

      DOUBLE PRECISION, INTENT (IN) ::   BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::            IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::   SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::            SIPIVOT ( MAXSTRMSTKS_2 )

!  BVP Column In/Out
!mick fix 7/29/2014 - placed COL2 & SCOL2 in call to make VLIDORT threadsafe
      DOUBLE PRECISION, INTENT (INOUT) ::  COL2  ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) ::  SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

      !  outputs

!  output
!  ------

!  4/9/19 Transmittance output. Only required for water-leaving adjustment.

      DOUBLE PRECISION, intent (inout) :: TRANS_ATMOS_FINAL  ( MAX_SZANGLES )

!  4/9/19. Surface leaving term. This has already been set (non-Waterleaving scenarios)
!             -- Will be iteratively adjusted for Waterleaving case.
      
      DOUBLE PRECISION, intent (inout) :: SL_QUADTERM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  Help array for the reflectance solution

      DOUBLE PRECISION, INTENT (INOUT) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  Solution constants of integration

      DOUBLE PRECISION, INTENT (INOUT) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE, TRACE_1, TRACE_2

!  Local variables
!  ---------------

      INTEGER          :: STATUS_SUB, O1
      LOGICAL          :: DO_ADJUSTED_BVP
      DOUBLE PRECISION :: TFACTOR

!  Include sleaving flag, introduced 4/9/19.

      LOGICAL :: DO_INCLUDE_SLEAVING
      
!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  4/9/19. include sleaving flag (for the Standard BVP)
!      DO_INCLUDE_SLEAVING = ( DO_SURFACE_LEAVING .and. FOURIER.eq.0 ).and. &
!           ( .not.DO_WATER_LEAVING .or. ( DO_WATER_LEAVING .and. DO_EXTERNAL_WLEAVE ) ) 
      
!  5/24/21. Version 2.8.3. Revised Condition for including SLTERMS in the BVProblem
!     -- DO_INCLUDE_SLEAVING now possible with All Fouriers (Water-leaving)
!     -- This should mirror the setting of DO_CALC_SLTERMS in DIRECT_RADIANCE.

!      DO_INCLUDE_SLEAVING = DO_SURFACE_LEAVING
!      if ( DO_SURFACE_LEAVING ) THEN
!         IF ( DO_WATER_LEAVING ) then
!           if ( .not.DO_EXTERNAL_WLEAVE .and. FOURIER .eq. 0 ) THEN
!             DO_INCLUDE_SLEAVING = .false.  !  ADJUSTED BVP solution, Only for M = 0
!           endif
!         ELSE
!           if ( FOURIER.gt.0 ) THEN
!             DO_INCLUDE_SLEAVING = .false.  ! SIF absent for M > 0
!           endif
!         ENDIF
!      ENDIF

      DO_INCLUDE_SLEAVING = .false.
      IF ( DO_SURFACE_LEAVING ) THEN
         IF ( DO_WATER_LEAVING ) THEN
            IF ( DO_EXTERNAL_WLEAVE ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  IF ( FOURIER.eq. 0 ) DO_INCLUDE_SLEAVING = .true.
               ELSE
                  DO_INCLUDE_SLEAVING = .true.
               ENDIF
            ELSE
               IF ( FOURIER.gt.0 ) DO_INCLUDE_SLEAVING = .true.
            ENDIF
         ELSE
            IF ( FOURIER.eq.0 )  DO_INCLUDE_SLEAVING = .true.
         ENDIF
      ENDIF

!  Regular BVP using compressed-band matrices, etc..
!  ==================================================

!  --Additional setups for the albedo layer
!   Version 2.8, Rearrange the argument list

      CALL BVP_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,       &
        FOURIER, NSTOKES, NSTREAMS, NLAYERS,             &
        SURFACE_FACTOR, ALBEDO, QUAD_STRMWTS,            &
        MUELLER_INDEX, WLOWER, AXBID_F,                  &
        R2_BEAM )
      !CALL TP2A (NSTREAMS,NSTOKES,R2_BEAM)

!  --set up Column for solution vector (the "B" as in AX=B)
!   Version 2.8, Rearrange the argument list
!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
!    4/22/19. Control for the surface-leaving contribution (added and refined)
!    1/31/21. Version 2.8.3. SL_QUADTERM is now non-zero for all FOurier Components (water-leaving)

      CALL BVP_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,      & ! Input surface flags
             DO_INCLUDE_SLEAVING, DO_INCLUDE_TOAFLUX,   DO_INCLUDE_BOAFLUX,        & ! Input SL/illum flags
             FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS,           & ! Input numbers
             NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,     & ! Input numbers
             TOAFLUX, BOAFLUX, SURFBB, EMISSIVITY,                 & ! Input TOA/BOA sources
             WUPPER, WLOWER, R2_BEAM, RF_DIRECT_BEAM, SL_QUADTERM, & ! Input solutions
             COL2, SCOL2 )                                           ! Output
      !CALL TP2B (NSTKS_NSTRMS_2,IBEAM,COL2)

!  4/9/19. Water-leaving Adjustment conditions for Fourier zero
!  ------------------------------------------------------------
      
!      - TRUE  if Water-leaving and DO_TF_Iteration, otherwise false
!      - FALSE if Water-leaving and .not. DO_TF_Iteration
!                 ===> SET TRANS_ATMOS_FINAL = 1.0          if  External_Wleave is set
!                 ===> SET TRANS_ATMOS_FINAL = Gordon value if  External_Wleave is NOT. set
!     Gordon's result = Square root of Solar Transmittance

!     1/31/21. Version 2.8.3. Zeroing of TRANS_ATMOS_FINAL only done for Fourier zero

!  5/24/21. Version 2.8.3. No longer using Gordon's result ==> CONSEQUENCE, WATER_LEAVING is:

!     (1) Coupled    , with DO_TF_ITERATION =.true., and TRANS_ATMOS_FINAL iterated (DO_ADJUSTED_BVP) for Fourier 0
!  or (2) Not Coupled, with DO_TF_ITERATION =.false., no adjustment and DO_EXTERNAL_WLEAVE = .true.

      DO_ADJUSTED_BVP =  ( DO_WATER_LEAVING .and. DO_TF_ITERATION .and. FOURIER.eq.0 )
      IF ( FOURIER.eq.0 ) TRANS_ATMOS_FINAL(IBEAM) = ZERO
      IF ( ( DO_WATER_LEAVING .and..not. DO_TF_ITERATION ) .and. FOURIER.eq.0 ) THEN
         IF ( DO_EXTERNAL_WLEAVE ) THEN
            TRANS_ATMOS_FINAL(IBEAM) = ONE
!         ELSE
!            TRANS_ATMOS_FINAL(IBEAM)    = SQRT(BEAM_BOATRANS(IBEAM))
         ENDIF
      ENDIF

!  1/31/21. Version 2.8.3. For Additional Fourier Components > 0, Water-leaving adjustment is necessary
!   -- This adjustment is done in DIRECTRADIANCE now.
!      IF ( DO_WATER_LEAVING .and. DO_TF_ITERATION .and. FOURIER.GT.0 ) then
!         TFACTOR = TRANS_ATMOS_FINAL(IBEAM) ; O1 = 1
!         SL_QUADTERM(1:NSTREAMS,IBEAM,O1) = TFACTOR * SL_QUADTERM(1:NSTREAMS,IBEAM,O1)
!      ENDIF

!   4/9/19 EITHER : Solve the ADJUSTED boundary value problem (Water-leaving flux-adjustment)
!  ==========================================================================================

!  THIS IS ONLY DONE FOR FOURIER = 0

      IF ( DO_ADJUSTED_BVP ) then    
         CALL BVP_ADJUSTED_BACKSUB &
           ( TF_MAXITER, TF_CRITERION, FOURIER, IBEAM,                           & ! TF iteration control, indices
             NSTOKES, NSTREAMS, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,           & ! Input Numbers basic
             NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,        & ! Input Numbers bookkeeping
             FLUX_FACTOR, FLUXVEC, QUAD_STRMWTS, BEAM_BOATRANS,                  & ! Input Quad/Flux
             DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,                      & ! Input surface leaving
             LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,               & ! Input Direct-flux
             SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, WLOWER,                         & ! Input RTE solutions
             IPIVOT, BANDMAT2, SIPIVOT, SMAT2,                                   & ! Input BVP matrices
             TRANS_ATMOS_FINAL, SL_QUADTERM, COL2, SCOL2, LCON, MCON,            & ! Modified Input/Output
             STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Output
         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            STATUS = VLIDORT_SERIOUS ; RETURN
         ENDIF
      ENDIF

!  4/9/19 OR    : Solve the STANDARD boundary problem (back substitution)
!  ======================================================================

      IF ( .not. DO_ADJUSTED_BVP ) then
         CALL BVP_BACKSUB &
            ( FOURIER, IBEAM, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
              NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,       &
              BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,         &
              LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )
         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Regular BVP Solution: back-substitution in BVP_SOLUTION_MASTER'
            STATUS = VLIDORT_SERIOUS ; RETURN
         ENDIF
      ENDIF
      !CALL TP2C (NSTKS_NSTRMS,NLAYERS,LCON,MCON)

!  finish

      RETURN
      END SUBROUTINE BVP_SOLUTION_MASTER

!

      SUBROUTINE BVP_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,       & ! Input flags
        FOURIER, NSTOKES, NSTREAMS, NLAYERS,             & ! Input numbers
        SURFACE_FACTOR, ALBEDO, QUAD_STRMWTS, & ! Input surface
        MUELLER_INDEX, WLOWER, AXBID_F,                  & ! Input solutions
        R2_BEAM )

!  This is the Lambertian or BRDF surface routine
!     --> Reflected beam solutions

!   Original code with BRDFs, July 26, 2010.  RT SOLUTIONS Inc,

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO construction with IF block

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAXSTOKES_SQ, ZERO

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

!  surface

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )

!  RTE stuff

      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER  ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      DOUBLE PRECISION :: REFL_B, REFL_B_S
      INTEGER          :: I, J, O1, O2, OM, NL

!  Initialization
!  ==============

!  Zero total reflected contributions

      R2_BEAM(1:NSTREAMS,1:NSTOKES)    = ZERO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if not the Fourier m = 0 component (Lambertian surface)

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .GT. 0 ) RETURN

!  Lambertian surface
!  ==================

      if ( DO_LAMBERTIAN_SURFACE ) THEN

!  Integrate Downward streams of particular solutions, set reflectance

        REFL_B = DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),WLOWER(1:NSTREAMS,1,NLAYERS))
        REFL_B = REFL_B * SURFACE_FACTOR * ALBEDO
        R2_BEAM(1:NSTREAMS,1) = REFL_B

!  BRDF surface
!  ============

      ELSE

!  Last layer

        NL = NLAYERS

!  Integrate Downward streams of particular solutions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B_S = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                REFL_B_S = REFL_B_S + AXBID_F(I,J,OM) * WLOWER(J,O2,NL)
              ENDDO
              REFL_B = REFL_B + REFL_B_S
            ENDDO
            R2_BEAM(I,O1) = SURFACE_FACTOR * REFL_B
          ENDDO
        ENDDO

!  End clause BRDF vs. Lambertian

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_SURFACE_SETUP_BEAM

!

      SUBROUTINE BVP_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,      & ! Input surface flags
             DO_INCLUDE_SLEAVING, DO_INCLUDE_TOAFLUX,   DO_INCLUDE_BOAFLUX,        & ! Input SL/illum flags
             FOURIER, IBEAM, NSTOKES, NSTREAMS, NLAYERS,           & ! Input numbers
             NSTREAMS_2, NTOTAL, NSTKS_NSTRMS, NSTKS_NSTRMS_2,     & ! Input numbers
             TOAFLUX, BOAFLUX, SURFBB, EMISSIVITY,                 & ! Input TOA/BOA sources
             WUPPER, WLOWER, R2_BEAM, RF_DIRECT_BEAM, SL_QUADTERM, & ! Input solutions
             COL2, SCOL2 )

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO construction with IF block

!   Version 2.8.1, Control for TOA/BOA  illumination added, 3/23/19
!    4/22/19. Control for the surface-leaving contribution (added and refined)

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                                 MAXSTREAMS_2, MAXSTRMSTKS_2, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Surface flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SLEAVING

!  Illumination flags. New 3/23/19

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_TOAFLUX
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_BOAFLUX

!  numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NTOTAL
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2

!  Surface Blackbody and Emissivity, and TOA/BOA Fluxes (new 3/23/19)

      DOUBLE PRECISION, INTENT (IN) ::  BOAFLUX
      DOUBLE PRECISION, INTENT (IN) ::  TOAFLUX
      DOUBLE PRECISION, INTENT (IN) ::  EMISSIVITY ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  SURFBB

!  particular solutions

      DOUBLE PRECISION, INTENT (IN) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Reflected Direct beam

      DOUBLE PRECISION, INTENT (IN) ::  RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  Surface-leaving (actually an input here)

      DOUBLE PRECISION, intent(inout)  :: SL_QUADTERM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER          :: M, I,I1,LAY,LAY1,C0,CM,O1,IR,IROW

!  Proxy (for debug only)

      M = FOURIER

!  General Case, NLAYERS > 1
!  =========================

      IF ( NLAYERS .GT. 1 ) THEN

!  zero column vector

        COL2(1:NTOTAL,IBEAM) = ZERO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

        LAY = 1
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            COL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
          ENDDO
        ENDDO

!  Version 2.8.1, Add illumination (assumed Isotropic) at TOA. 3/23/19

        IF ( DO_INCLUDE_TOAFLUX ) then
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1) ; O1 = 1 ; IROW = IR + O1
            COL2(IROW,IBEAM)   = COL2(IROW,IBEAM) + TOAFLUX
          ENDDO
        ENDIF
       
!  intermediate layer boundaries (will not be done if NLAYERS = 1 )
!  -----------------------------

        DO LAY = 2, NLAYERS
          LAY1 = LAY - 1
          C0 = LAY1*NSTKS_NSTRMS_2 - NSTKS_NSTRMS
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COL2(CM,IBEAM) = WUPPER(I,O1,LAY) - WLOWER(I,O1,LAY1)
            ENDDO
          ENDDO
        ENDDO

!  lowest (surface) boundary (diffuse radiation terms only)
!  --------------------------------------------------------

        LAY = NLAYERS
        C0 = (LAY-1)*NSTKS_NSTRMS_2 + NSTKS_NSTRMS

!  with surface contribueions, include integrated downward reflectances R2_BEAM
!    no surface, similar code excluding integrated reflectance

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
            ENDDO
          ENDDO
        ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          IF ( DO_INCLUDE_SURFACE ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                COL2(CM,IBEAM) = COL2(CM,IBEAM) + RF_DIRECT_BEAM(I,IBEAM,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Add surface-leaving contribution  (only to final level)

        IF ( DO_INCLUDE_SLEAVING ) THEN
           DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                COL2(CM,IBEAM) = COL2(CM,IBEAM) + SL_QUADTERM(I,IBEAM,O1)
              ENDDO
           ENDDO
        ENDIF
         
!  Add thermal emission of ground surface (only to final level)
!  ------------------------------------------------------------

!mick fix 2/14/2012 - changed treatment of emissivity in this NLAYERS > 1
!                     section to be consistent with LIDORT
!mick fix 9/19/2017 - include possibility if polarized emissivity: added inner Stokes loop

!      IF ( DO_INCLUDE_SURFEMISS ) THEN
!        O1 = 1
!        LOCAL_EMISS = ONE - ALBEDO
!        FP_SBB = SURFBB
!        IF ( DO_LAMBERTIAN_SURFACE ) THEN
!          DO I = 1, NSTREAMS
!            CM = C0 + NSTOKES*(I-1) + O1
!            COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * local_emiss
!          ENDDO
!        ELSE
!          DO I = 1, NSTREAMS
!            CM = C0 + NSTOKES*(I-1) + O1
!!  @@@ Rob Fix. Ordering of EMISSIVITY indices wrong
!!           COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(I,O1)
!            COL2(CM,IBEAM) = COL2(CM,IBEAM) + FP_SBB * EMISSIVITY(O1,I)
!          ENDDO
!        ENDIF
!      ENDIF

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          !O1 = 1
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              CM = C0 + NSTOKES*(I-1) + O1
              COL2(CM,IBEAM) = COL2(CM,IBEAM) + SURFBB * EMISSIVITY(O1,I)
            ENDDO
          ENDDO
        ENDIF

!  Version 2.8.1, Add illumination (assumed Isotropic) at BOA. 3/23/19

        IF ( DO_INCLUDE_BOAFLUX ) then
          DO I = 1, NSTREAMS
            O1 = 1 ; CM = C0 + NSTOKES*(I-1) + O1
            COL2(CM,IBEAM)   = COL2(CM,IBEAM) + BOAFLUX
          ENDDO
        ENDIF
       
!  debug

!      do i = 1, ntotal
!        write(500,*)i,COL2(i,IBEAM)
!      enddo
!      pause'hello'

!  special case, NLAYERS = 1
!  =========================

      ELSE

!  zero column vector

        SCOL2(1:NTOTAL,IBEAM) = ZERO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

        LAY = 1
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,LAY)
          ENDDO
        ENDDO

!  Version 2.8.1, Add illumination (assumed isotropic) at TOA. 3/23/19

        IF ( DO_INCLUDE_TOAFLUX ) then
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1) ; O1 = 1 ; IROW = IR + O1
            SCOL2(IROW,IBEAM)   = SCOL2(IROW,IBEAM) + TOAFLUX
          ENDDO
       ENDIF
       
!  lowest (surface) boundary (diffuse radiation terms only)
!  --------------------------------------------------------

!  with non-zero surface, include integrated downward reflectance R2_BEAM
!    no surface, similar code excluding integrated reflectance

        C0 = NSTKS_NSTRMS
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = - WLOWER(I1,O1,LAY)
            ENDDO
          ENDDO
        ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          IF ( DO_INCLUDE_SURFACE ) THEN
            DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + RF_DIRECT_BEAM(I,IBEAM,O1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Add surface-leaving contribution  (only to final level)

        IF ( DO_INCLUDE_SLEAVING ) THEN
           DO I = 1, NSTREAMS
              IR = NSTOKES*(I-1)
              DO O1 = 1, NSTOKES
                IROW = IR + O1
                CM = C0 + IROW
                SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SL_QUADTERM(I,IBEAM,O1)
              ENDDO
           ENDDO
        ENDIF
        
!  Add thermal emission of ground surface (only to final level)
!  ------------------------------------------------------------

!mick fix 2/14/2012 - simplified loop structure here to make it like the
!                     NLAYERS > 1 thermal section above
!mick fix 9/19/2017 - include possibility if polarized emissivity: added inner Stokes loop

!      IF ( DO_INCLUDE_SURFEMISS ) THEN
!        NSTOKES_FIRST = 1
!        DO I = 1, NSTREAMS
!          IR = NSTOKES*(I-1)
!          DO O1 = 1, NSTOKES_FIRST
!            IROW = IR + O1
!            CM = C0 + IROW
!!  @@@ Rob Fix. Ordering of EMISSIVITY indices wrong
!!           SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB*EMISSIVITY(I,O1)
!            SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB*EMISSIVITY(O1,I)
!          ENDDO
!        ENDDO
!      ENDIF

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          !O1 = 1
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              CM = C0 + NSTOKES*(I-1) + O1
              SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB * EMISSIVITY(O1,I)
            ENDDO
          ENDDO
        ENDIF

!  Version 2.8.1, Add illumination (assumed Isotropic) at BOA. 3/23/19

        IF ( DO_INCLUDE_BOAFLUX ) then
          DO I = 1, NSTREAMS
            O1 = 1 ; CM = C0 + NSTOKES*(I-1) + O1
            SCOL2(CM,IBEAM)   = SCOL2(CM,IBEAM) + BOAFLUX
          ENDDO
        ENDIF

!  End clause NLAYERS > 1 vs. NLAYERS = 1

      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVP_COLUMN_SETUP

!

      SUBROUTINE BVP_BACKSUB ( &
        FOURIER, IBEAM, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,       &
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,         &
        LCON, MCON, STATUS, MESSAGE, TRACE )

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.

      USE VLIDORT_PARS_m, only : MAXLAYERS, MAXBEAMS, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS
      USE LAPACK_TOOLS_m, only : DGBTRS, DGETRS

      IMPLICIT NONE

!  indices

      INTEGER, INTENT (IN) ::             FOURIER
      INTEGER, INTENT (IN) ::             IBEAM

!  Numbers

      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NTOTAL
      INTEGER, INTENT (IN) ::             N_SUBDIAG
      INTEGER, INTENT (IN) ::             N_SUPDIAG
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::             K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_COMPLEX ( MAXLAYERS )

!  BVP Matrices

      DOUBLE PRECISION, INTENT (IN) ::    BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::             IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::    SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::             SIPIVOT ( MAXSTRMSTKS_2 )

!  BVP Columns, and results

      DOUBLE PRECISION, INTENT (INOUT) :: COL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Error handling

      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           M, C0, LAY, K, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S
      INTEGER ::           INFO
      CHARACTER (LEN=3) :: CI
      CHARACTER (LEN=2) :: CB

!  mick fix 7/23/2014. Local arrays defined
      DOUBLE PRECISION  :: LOC_COL2  ( MAXTOTAL, 1 )
      DOUBLE PRECISION  :: LOC_SCOL2 ( MAXSTRMSTKS_2, 1 )

!  Proxy (for debug only)

      M = FOURIER

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  BVP back-substitution: With compression (multilayers)
!  =====================================================

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("COL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended
        LOC_COL2(1:NTOTAL,1) = COL2(1:NTOTAL,IBEAM)
        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              LOC_COL2, MAXTOTAL, INFO )
        COL2(1:NTOTAL,IBEAM) = LOC_COL2(1:NTOTAL,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVP_BACKSUB, Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,LAY) = COL2(C0+IROW,IBEAM)
            MCON(K,LAY) = COL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,LAY) = COL2(C0+IROW,   IBEAM)
            LCON(K2,LAY) = COL2(C0+IROW_S, IBEAM)
            MCON(K1,LAY) = COL2(C0+IROW1,  IBEAM)
            MCON(K2,LAY) = COL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO

!  debug---------------------- Removed, Old F77 code

!  Solve the boundary problem: No compression, Single Layer only
!  =============================================================

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NSTKS_NSTRMS_2,1) = SCOL2(1:NSTKS_NSTRMS_2,IBEAM)
        CALL DGETRS ( 'N', NSTKS_NSTRMS_2, 1, &
                      SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                      LOC_SCOL2, MAXSTRMSTKS_2, INFO )
        SCOL2(1:NSTKS_NSTRMS_2,IBEAM) = LOC_SCOL2(1:NSTKS_NSTRMS_2,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (nlayers=1) in BVP_BACKSUB, Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,LAY) = SCOL2(IROW,IBEAM)
          MCON(K,LAY) = SCOL2(IROW1,IBEAM)
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW = K + K_REAL(LAY)
          IROW1 = IROW + NSTKS_NSTRMS
          IROW_S = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,LAY) = SCOL2(IROW,    IBEAM)
          LCON(K2,LAY) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,LAY) = SCOL2(IROW1,   IBEAM)
          MCON(K2,LAY) = SCOL2(IROW1_S, IBEAM)
        ENDDO

!  End clause NLAYERS > 1 vs. NLAYERS = 1

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_BACKSUB

!

      SUBROUTINE BVP_BACKSUB_1 ( &
        NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,           &
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, &
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,   &
        LCON, MCON, STATUS, MESSAGE, TRACE )

!  Solves the boundary value problem. Dedicated routine for the Mediaprops Code

      USE VLIDORT_PARS_m, only : MAXLAYERS, MAXSTRMSTKS, MAXSTRMSTKS_2, &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      USE LAPACK_TOOLS_m, only : DGBTRS, DGETRS
      
      IMPLICIT NONE

!  Control integers

      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NTOTAL
      INTEGER, INTENT (IN) ::             N_SUBDIAG
      INTEGER, INTENT (IN) ::             N_SUPDIAG
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::             K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_COMPLEX ( MAXLAYERS )

!  BVP Matrices and pivots

      DOUBLE PRECISION, INTENT (IN) ::    BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::    SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

      INTEGER, INTENT (IN) ::             IPIVOT  ( MAXTOTAL )
      INTEGER, INTENT (IN) ::             SIPIVOT ( MAXSTRMSTKS_2 )

!  BVP Columns, and results

      DOUBLE PRECISION, INTENT (INOUT) :: COL2  ( MAXTOTAL, 1 )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, 1 )

!  output
!  ------

!  Solution constants of integration

      DOUBLE PRECISION, INTENT (INOUT) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Error handling

      INTEGER          , INTENT (OUT)  :: STATUS
      CHARACTER (LEN=*), INTENT (OUT)  :: MESSAGE, TRACE

!  local variables
!  ---------------

      INTEGER ::           C0, LAY, K, KO1, K0, K1, K2
      INTEGER ::           IROW, IROW1, IROW_S, IROW1_S, INFO
      CHARACTER (LEN=3) :: CI

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  BVP back-substitution: With compression (multilayers)
!  =====================================================

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT,      &
              COL2, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVP_BACKSUB_1 '
          STATUS  = VLIDORT_SERIOUS ; RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTKS_NSTRMS_2
          KO1 = K_REAL(LAY) + 1
          DO K = 1, K_REAL(LAY)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,LAY) = COL2(C0+IROW,1)
            MCON(K,LAY) = COL2(C0+IROW1,1)
          ENDDO
          DO K = 1, K_COMPLEX(LAY)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(LAY)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(LAY)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,LAY) = COL2(C0+IROW,   1)
            LCON(K2,LAY) = COL2(C0+IROW_S, 1)
            MCON(K1,LAY) = COL2(C0+IROW1,  1)
            MCON(K2,LAY) = COL2(C0+IROW1_S,1)
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  =============================================================

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS ( 'N', NSTKS_NSTRMS_2, 1, &
                      SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
                      SCOL2, MAXSTRMSTKS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (nlayers=1) in BVP_BACKSUB_1'
          STATUS  = VLIDORT_SERIOUS ; RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

        LAY = 1
        KO1 = K_REAL(LAY) + 1
        DO K = 1, K_REAL(LAY)
          IROW = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,LAY) = SCOL2(IROW,1)
          MCON(K,LAY) = SCOL2(IROW1,1)
        ENDDO
        DO K = 1, K_COMPLEX(LAY)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW = K + K_REAL(LAY)
          IROW1 = IROW + NSTKS_NSTRMS
          IROW_S = K + K_REAL(LAY) + K_COMPLEX(LAY)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,LAY) = SCOL2(IROW,    1)
          LCON(K2,LAY) = SCOL2(IROW_S,  1)
          MCON(K1,LAY) = SCOL2(IROW1,   1)
          MCON(K2,LAY) = SCOL2(IROW1_S, 1)
        ENDDO

!  End clause NLAYERS > 1 vs. NLAYERS = 1

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVP_BACKSUB_1

!

SUBROUTINE BVP_ADJUSTED_BACKSUB &
        ( TF_MAXITER, TF_CRITERION, FOURIER, IBEAM,                    & ! TF iteration control, indices
          NSTOKES, NSTREAMS, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,    & ! Input Numbers basic
          NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, & ! Input Numbers bookkeeping
          FLUX_FACTOR, FLUXVEC, QUAD_STRMWTS, SOLARBEAM_BOATRANS,      & ! Input Quad/Flux
          DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,               & ! Input surface leaving
          LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,        & ! Input Direct-flux
          SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, WLOWER,                  & ! Input RTE solutions
          IPIVOT, BANDMAT2, SIPIVOT, SMAT2,                            & ! Input BVP matrices
          TRANS_ATMOS_FINAL, SL_QUADTERM, COL2, SCOL2, LCON, MCON,     & ! Modified Input/Output
          STATUS, MESSAGE, TRACE_1, TRACE_2 )                            ! Output

!  Special Fourier = 0 component calculation for the adjusted atmospheric downwelling Flux at BOA

!  1/31/21. Version 2.8.3. SLEAVE input array defined locally, each Fourier, remove MAXMOMENTS from parameter list

      USE VLIDORT_PARS_m, Only : MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_SZANGLES, MAXBANDTOTAL, &
                                  MAXTOTAL, MAXEVALUES, MAXSTREAMS_2, MAXSTRMSTKS, MAXSTRMSTKS_2,        &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO, ONE, TWO, PI2

      USE LAPACK_TOOLS_m, only : DGBTRS, DGETRS

      implicit none

!  Input
!  -----

!  beam number, Fourier index

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  New for version 3.8. 7/6/16. Iteration control

      INTEGER         , INTENT (IN) :: TF_MAXITER
      DOUBLE PRECISION, INTENT (IN) :: TF_CRITERION

!  Input numbers

      INTEGER, INTENT (IN) :: NSTOKES, NSTREAMS, NLAYERS
      INTEGER, INTENT (IN) :: NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Version 2p7 input, 2/19/14 (bookkeeping)

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR, FLUXVEC(MAXSTOKES)

!  Quadrature

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Direct beam input

      DOUBLE PRECISION, INTENT (IN) :: SOLARBEAM_BOATRANS  ( MAXBEAMS )

!  Surface-leaving quantities, these are the original supplement results (not adjusted by tranmsittance)
!    -- 1/31/21. Version 2.8.3. SLEAVE input array defined locally, remove MAXMOMENTS dimension

      LOGICAL         , INTENT (IN) :: DO_SL_ISOTROPIC
      DOUBLE PRECISION, INTENT (IN) :: SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: SLTERM_F_0       ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )

!  Direct-flux inputs (Optical properties, Average-secant parameterization)

      INTEGER         , INTENT (IN) :: BEAM_CUTOFF      ( MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_MUBAR     ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LOCAL_CSZA  ( 0:MAXLAYERS, MAXBEAMS )

!  Homogeneous solution inputs

      INTEGER         , INTENT (IN) :: K_REAL       ( MAXLAYERS )
      INTEGER         , INTENT (IN) :: K_COMPLEX    ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Particular integral inputs

      DOUBLE PRECISION, INTENT (IN) :: WLOWER  ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  BVP Band matrices and pivots

      DOUBLE PRECISION, INTENT (IN) :: BANDMAT2  ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER         , INTENT (IN) :: IPIVOT    ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) :: SMAT2     ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER         , INTENT (IN) :: SIPIVOT   ( MAXSTRMSTKS_2 )

!  Output
!  ======

!  Transmittance output

      DOUBLE PRECISION, INTENT (INOUT) :: TRANS_ATMOS_FINAL  ( MAX_SZANGLES )

!  Adjusted Surface-leaving contribution (actually an output from here)

      DOUBLE PRECISION, INTENT (INOUT) :: SL_QUADTERM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES  )

!  Adjusted BVP column vectors. MODIFIED INPUTS

      DOUBLE PRECISION, intent(inout) :: COL2      ( MAXTOTAL,      MAXBEAMS )
      DOUBLE PRECISION, intent(inout) :: SCOL2     ( MAXSTRMSTKS_2, MAXBEAMS )

!  Solution constants of integration

      DOUBLE PRECISION, intent(inout) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, intent(inout) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_1, TRACE_2

!  LOCAL VARIABLES
!  ===============

!  Saved particular integral at lower boundary

      DOUBLE PRECISION :: WLOWER_SAVED       ( MAXSTRMSTKS_2, MAXSTOKES )

!  BVProblem, local arrays

      DOUBLE PRECISION :: COL2_LOCAL  ( MAXTOTAL, 1 )
      DOUBLE PRECISION :: SCOL2_LOCAL ( MAXSTRMSTKS_2, 1 )
      DOUBLE PRECISION :: LCON_LOCAL  ( MAXSTRMSTKS )
      DOUBLE PRECISION :: MCON_LOCAL  ( MAXSTRMSTKS )

!  Adjusted surface-leaving Direct Beam
!      DOUBLE PRECISION :: SL_DIRECT_BEAM ( MAXSTREAMS, MAXSTOKES )

!  Discrete ordinate downwelling field at the lowest level

      DOUBLE PRECISION :: QSTOKES_F    ( MAXSTREAMS, MAXSTOKES )

!  Transmittances and Fluxes

      DOUBLE PRECISION :: TRANS_ATMOS_IN, TRANS_ATMOS
      DOUBLE PRECISION :: DIRECT_TRANS, DIRECT_FLUX, FLUX_DIRECT(MAXSTOKES), FLUX(MAXSTOKES)

!  Flux multiplier and surface factors

      DOUBLE PRECISION :: FLUX_MULTIPLIER
      DOUBLE PRECISION :: DELTA_FACTOR
      DOUBLE PRECISION :: SURFACE_FACTOR  

!  Helper variables

      character*2      :: CB
      character*3      :: CI
      LOGICAL          :: DO_ITERATION
      INTEGER          :: IR, IROW, IROW1, IROW_S, IROW1_S
      INTEGER          :: K, KO1, K0, K1, K2
      INTEGER          :: JITER, I, O1, C0, C1, CM, N, INFO, STATUS_SUB
      DOUBLE PRECISION :: TFACTOR, SL, CONV, SHOM_R, SHOM_CR
      DOUBLE PRECISION :: LXR, MXR, LXR_CR, LXR_CI, MXR_CR

!  START OF CODE
!  =============

!  Initialize output and exception handling 

      TRANS_ATMOS_FINAL(IBEAM) = one
      STATUS = VLIDORT_SUCCESS
      MESSAGE = ' ' ; TRACE_1 = ' ' ; TRACE_2 = ' '
      
!  SETUP CALCULATION. Only need the saved column using basic direct beam.
!  =================

!  surface reflectance factors, flux multiplier

      SURFACE_FACTOR = TWO
      DELTA_FACTOR   = ONE
      FLUX_MULTIPLIER = DELTA_FACTOR

!  Saved the lowest-layer particular integral

      WLOWER_SAVED(1:NSTREAMS_2,1:NSTOKES) =  WLOWER(1:NSTREAMS_2,1:NSTOKES,NLAYERS) 

!  PART I. ITERATION LOOP to find adjusted transmittance-flux
!  ==========================================================

!  Compute Direct Flux. Should be initialized if direct beam does not reach ground.

      DIRECT_FLUX = zero
      IF ( NLAYERS .LE. BEAM_CUTOFF(IBEAM) ) THEN
         DIRECT_TRANS = INITIAL_TRANS(NLAYERS,IBEAM) * T_DELT_MUBAR(NLAYERS,IBEAM)
         DIRECT_FLUX  = FLUX_FACTOR * DIRECT_TRANS * LOCAL_CSZA(NLAYERS,IBEAM)
         FLUX_DIRECT(1:NSTOKES) = DIRECT_FLUX * FLUXVEC(1:NSTOKES)
      ENDIF

!  initialize loop

      DO_ITERATION = .true. ; JITER = 0 

!  Initial guess - use Square root of Solar Transmittance. (Gordon's result).
!                - If not available set to COS(SZA) as default

      IF ( NLAYERS .LE. BEAM_CUTOFF(IBEAM) ) THEN
         TRANS_ATMOS_IN = SQRT(SOLARBEAM_BOATRANS(IBEAM))
      ELSE
         TRANS_ATMOS_IN = LOCAL_CSZA(NLAYERS,IBEAM)
      ENDIF
!      write(*,*)'start',TRANS_ATMOS_IN, 1.5*DIRECT_TRANS

!  Offset

      C0 = NTOTAL - NSTKS_NSTRMS
      C1 = C0     - NSTKS_NSTRMS

!  Start iteration

      DO WHILE ( DO_ITERATION .and.JITER.lt.TF_MAXITER )

!  Copy the saved column vector. Only the lowest entries will be adjusted
!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("COL2") to a local array with
!         a 2nd dim of size one to use the LAPACK routine DGBTRS again as intended

         IF ( NLAYERS.gt.1 ) then
            COL2_LOCAL(1:NTOTAL,1) = COL2(1:NTOTAL,IBEAM)
         ELSE
            SCOL2_LOCAL(1:NTOTAL,1) = SCOL2(1:NTOTAL,IBEAM)
         ENDIF

!  increase step

         JITER = JITER + 1

!  Find adjusted surface-leaving contribution (unpolarized)
!    -- Zero the remaining polarized entries
!    -- 1/31/21. Version 2.8.3. SLEAVE input array defined locally, remove Fourier index

         O1 = 1 ; TFACTOR = TRANS_ATMOS_IN * FLUX_FACTOR / DELTA_FACTOR
         IF ( DO_SL_ISOTROPIC ) THEN
            SL_QUADTERM(1:NSTREAMS,IBEAM,O1) =  SLTERM_ISOTROPIC(O1,IBEAM) * TFACTOR
         ELSE
            DO I = 1, NSTREAMS
               SL_QUADTERM(I,IBEAM,O1) = SLTERM_F_0(O1,I,IBEAM) * TFACTOR
            ENDDO
         ENDIF
         SL_QUADTERM(1:NSTREAMS,IBEAM,2:NSTOKES) = zero
         
!  Add adjusted surface-leaving contribution to lowest entries in column vector 

         IF ( NLAYERS .gt. 1 ) then
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1)
               DO O1 = 1, NSTOKES
                  IROW = IR + O1 ; CM = C0 + IROW
                  COL2_LOCAL(CM,1) = COL2_LOCAL(CM,1) + SL_QUADTERM(I,IBEAM,O1)
               ENDDO
            ENDDO
         ELSE
            DO I = 1, NSTREAMS
               IR = NSTOKES*(I-1)
               DO O1 = 1, NSTOKES
                  IROW = IR + O1 ; CM = C0 + IROW
                  SCOL2_LOCAL(CM,1) = SCOL2_LOCAL(CM,1) + SL_QUADTERM(I,IBEAM,O1)
               ENDDO
            ENDDO
         ENDIF
         
!  Perform LAPACK substitution (DGBTRS) using RHS column vector COL2

         IF ( NLAYERS .gt. 1 ) then
           CALL DGBTRS &
             ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
                BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_LOCAL, MAXTOTAL, INFO )
         ELSE
            CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTRMSTKS_2, SIPIVOT, SCOL2_LOCAL, MAXSTRMSTKS_2, INFO )
         ENDIF

!  (error tracing)

         IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO ; WRITE(CB, '(I2)' ) IBEAM
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE_1 = 'DGBTRS/DGETRS call in iteration loop, Beam # '//CB
            TRACE_2 = 'Iterated Backsubstitution Error in BVP_ADJUSTED_BACKSUB, called by BVP_SOLUTION_MASTER'
            STATUS  = VLIDORT_SERIOUS ; RETURN
         ENDIF

!  Only interested in lowest layer

          N = NLAYERS
          
!  Set integration constants LCON and MCON for -/+ eigensolutions, Lowest layer only

         IF ( N.gt.1 ) then
            KO1 = K_REAL(N) + 1
            DO K = 1, K_REAL(N)
               IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
               LCON_LOCAL(K) = COL2_LOCAL(C1+IROW,1) ; MCON_LOCAL(K) = COL2_LOCAL(C1+IROW1,1)
            ENDDO
            DO K = 1, K_COMPLEX(N)
               K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
               IROW    = K + K_REAL(N)        ; IROW1   = IROW + NSTKS_NSTRMS
               IROW_S  = IROW + K_COMPLEX(N)  ; IROW1_S = IROW_S + NSTKS_NSTRMS
               LCON_LOCAL(K1) = COL2_LOCAL(C1+IROW, 1) ; LCON_LOCAL(K2) = COL2_LOCAL(C1+IROW_S, 1)
               MCON_LOCAL(K1) = COL2_LOCAL(C1+IROW1,1) ; MCON_LOCAL(K2) = COL2_LOCAL(C1+IROW1_S,1)
            ENDDO
         ELSE
            KO1 = K_REAL(N) + 1
            DO K = 1, K_REAL(N)
               IROW = K ; IROW1 = IROW + NSTKS_NSTRMS
               LCON_LOCAL(K) = SCOL2_LOCAL(IROW,1) ; MCON_LOCAL(K) = SCOL2_LOCAL(IROW1,1)
            ENDDO
            DO K = 1, K_COMPLEX(N)
               K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
               IROW    = K + K_REAL(N)        ; IROW1   = IROW + NSTKS_NSTRMS
               IROW_S  = IROW + K_COMPLEX(N)  ; IROW1_S = IROW_S + NSTKS_NSTRMS
               LCON_LOCAL(K1) = SCOL2_LOCAL(IROW, 1) ; LCON_LOCAL(K2) = SCOL2_LOCAL(IROW_S, 1)
               MCON_LOCAL(K1) = SCOL2_LOCAL(IROW1,1) ; MCON_LOCAL(K2) = SCOL2_LOCAL(IROW1_S,1)
            ENDDO
         ENDIF

!  Compute downwelling field at the surface

         DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
               SHOM_R = ZERO
               DO K = 1, K_REAL(N)
                  LXR  = LCON_LOCAL(K)*SOLA_XPOS(I,O1,K,N)
                  MXR  = MCON_LOCAL(K)*SOLB_XNEG(I,O1,K,N)
                  SHOM_R = SHOM_R + LXR * T_DELT_EIGEN(K,N) + MXR
               ENDDO
               SHOM_CR = ZERO
               DO K = 1, K_COMPLEX(N)
                  K0 = 2*K-2 ; K1 = KO1 + K0 ;  K2 = K1  + 1
                  LXR_CR =  LCON_LOCAL(K1) * SOLA_XPOS(I,O1,K1,N) - LCON_LOCAL(K2) * SOLA_XPOS(I,O1,K2,N)
                  LXR_CI =  LCON_LOCAL(K1) * SOLA_XPOS(I,O1,K2,N) + LCON_LOCAL(K2) * SOLA_XPOS(I,O1,K1,N)
                  MXR_CR =  MCON_LOCAL(K1) * SOLB_XNEG(I,O1,K1,N) - MCON_LOCAL(K2) * SOLB_XNEG(I,O1,K2,N)
                  SHOM_CR = SHOM_CR + LXR_CR*T_DELT_EIGEN(K1,N) - LXR_CI*T_DELT_EIGEN(K2,N) + MXR_CR
               ENDDO
               QSTOKES_F(I,O1) = FLUX_MULTIPLIER * ( WLOWER_SAVED(I,O1) + SHOM_R + SHOM_CR )
            ENDDO
         ENDDO

!  Compute Diffuse Flux, Add direct flux

         DO O1 = 1, NSTOKES
           FLUX(O1) = PI2 * DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QSTOKES_F(1:NSTREAMS,O1))
           FLUX(O1) = FLUX(O1) + FLUX_DIRECT(O1)
        ENDDO
        
!  Compute flux transmittance

         TRANS_ATMOS = FLUX(1) / LOCAL_CSZA(NLAYERS,IBEAM) / FLUX_FACTOR

!  Key debug.....
!         write(*,*)IBEAM,JITER,DIRECT_TRANS,&
!              (FLUX(1)-FLUX_DIRECT(1))/ LOCAL_CSZA(NLAYERS,IBEAM) / FLUX_FACTOR, TRANS_ATMOS_IN, TRANS_ATMOS
        
!  Examine convergence. If not set, update the transmittance, beam-by-beam

         CONV = ABS(ONE-(TRANS_ATMOS_IN/TRANS_ATMOS))
         IF ( CONV .lt. TF_CRITERION ) THEN
            DO_ITERATION = .false.
            TRANS_ATMOS_FINAL(IBEAM) = TRANS_ATMOS
         ELSE
            TRANS_ATMOS_IN = TRANS_ATMOS
         ENDIF
         
!  debug
         
!         write(*,*)ibeam,jiter,SLTERM_ISOTROPIC(IBEAM),SL_DIRECT_BEAM(1),TRANS_ATMOS,CONV,DO_ITERATION

!  Done iteration

      ENDDO

!  error due to lack of convergence

      IF ( DO_ITERATION ) THEN
         message = 'Adjustment iteration failed; TRANS_ATMOS not converging'
         TRACE_1 = 'After iteration loop  in BVP_ADJUSTED_BACKSUB'
         TRACE_2 = 'Error in BVP_ADJUSTED_BACKSUB, called by BVP_SOLUTION_MASTER'
         STATUS = VLIDORT_SERIOUS ; RETURN
      ENDIF

!  PART II - Post-iteration: Adjust SL, update Direct beam, Column vectors, final solution
!  =======================================================================================

!   - Reserve scaling on USERANGLES surface leaving until after FOURIER 0 call. Looks like this
!                      SCALED_SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM,ISCENE) = &
!                             SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) * TFACTOR

!  First, get final value of the adjusted Direct Beam surface-leaving contribution
!    -- 1/31/21. Version 2.8.3. SLEAVE input array defined locally, remove Fourier index

      O1 = 1 ; TFACTOR = TRANS_ATMOS_FINAL(IBEAM) * FLUX_FACTOR / DELTA_FACTOR
      IF ( DO_SL_ISOTROPIC ) THEN
         SL = SLTERM_ISOTROPIC(O1,IBEAM) * TFACTOR
         SL_QUADTERM(1:NSTREAMS,IBEAM,O1) =  SL
      ELSE
         DO I = 1, NSTREAMS
            SL_QUADTERM(I,IBEAM,O1) = SLTERM_F_0(O1,I,IBEAM) * TFACTOR
         ENDDO
      ENDIF

!  Second, update saved  COL2, SCOL2 by adding the surface leaving final contribution

      IF ( NLAYERS .gt. 1 ) then
         DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
               IROW = IR + O1 ; CM = C0 + IROW
               COL2(CM,IBEAM) = COL2(CM,IBEAM) + SL_QUADTERM(I,IBEAM,O1)
            ENDDO
         ENDDO
      ELSE
         DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
               IROW = IR + O1 ; CM = C0 + IROW
               SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SL_QUADTERM(I,IBEAM,O1)
            ENDDO
         ENDDO
      ENDIF

!  Final backsubstitution with the final COL2, get LCON and MCON

      CALL BVP_BACKSUB ( &
           FOURIER, IBEAM, NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
           NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX,       &
           BANDMAT2, IPIVOT, SMAT2, SIPIVOT, COL2, SCOL2,         &
           LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )       ! output

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
         TRACE_2 = 'Final Backsubstitution Error in BVP_ADJUSTED_BACKSUB, called by BVP_SOLUTION_MASTER'
         STATUS = VLIDORT_SERIOUS ; RETURN
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_ADJUSTED_BACKSUB  

!

SUBROUTINE BVPTEL_MATRIXSETUP_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,          & ! Input Flags
        DO_LAYER_SCATTERING, FOURIER, NSTOKES, NSTREAMS,    & ! Input Numbers
        NLAYERS, NSTREAMS_2, NSTKS_NSTRMS_2, NSTKS_NSTRMS,  & ! Input Numbers
        DO_BVTEL_INITIAL, BVTEL_FOURIER,                    & ! BVP Tel Control
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,     & ! BVP Tel Control
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,       & ! BVP Tel Control
        SURFACE_FACTOR, MUELLER_INDEX, K_REAL, K_COMPLEX,   & ! Input Bookkeeping
        QUAD_STRMWTS, ALBEDO, BRDF_F,            & ! Input Surface inputs
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_DELT_DISORDS, & ! Input RTE stuff
        R2_HOMP, R2_HOMM, AXBID_F, CUMTRANS, CUMQUAD,       & ! Output Surface reflection
        BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,             & ! Output BVP Matrices
        STATUS, MESSAGE, TRACE )                              ! Exception handling

!  Sets up the telescoped boundary value problem.
!    Standard case: Fourier > 0. With surface reflection term

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO constructions with IF blocks
!     * Introduce general treatment of surfaces

!   1/31/21. Version 2.8.3. BRDF/SLEAVE input arrays defined locally, remove Fourier index

      USE VLIDORT_PARS_m, Only : MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXSTREAMS_2, MAXSTRMSTKS_2, &
                                 MAXEVALUES, MAXSTOKES_SQ, MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::            DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::            DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
!      LOGICAL, INTENT (IN) ::            DO_SPECIALIST_OPTION_2 ! removed, Version 2.8

!  Numbers

      INTEGER, INTENT (IN) ::            FOURIER
      INTEGER, INTENT (IN) ::            NSTOKES
      INTEGER, INTENT (IN) ::            NSTREAMS
      INTEGER, INTENT (IN) ::            NLAYERS

      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::            NSTREAMS_2
      INTEGER, INTENT (IN) ::            NSTKS_NSTRMS_2

!  BVP Tel Control

      LOGICAL, INTENT (INOUT) ::         DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::         BVTEL_FOURIER
      INTEGER, INTENT (OUT) ::           N_BVTELMATRIX_SIZE
      INTEGER, INTENT (OUT) ::           N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (OUT) ::           N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (OUT) ::           NLAYERS_TEL
      INTEGER, INTENT (OUT) ::           ACTIVE_LAYERS ( MAXLAYERS )

!  Bookkeeping

      DOUBLE PRECISION, INTENT (IN) ::   SURFACE_FACTOR
      INTEGER, INTENT (IN) ::            MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::            K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::            K_COMPLEX ( MAXLAYERS )

!  surface inputs
!     -- 1/31/21. Version 2.8.3. BRDF input array defined locally, remove MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) ::   QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::   ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::   BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::   SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   T_DELT_EIGEN   ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::   T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  BVP Matrices

      DOUBLE PRECISION, INTENT (OUT) ::  BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (OUT) ::           IPIVOTTEL ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT) ::  SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::           SIPIVOT ( MAXSTRMSTKS_2 )

!  surface-reflected solutions

      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) ::  AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  Local output from Surface subroutine

      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION, INTENT (OUT) :: CUMQUAD (MAXSTREAMS)

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  local variables
!  ---------------

      INTEGER :: STATUS_SUB, NLAST

!  local output from initialization

      INTEGER :: NMINTEL ( MAXTOTAL )
      INTEGER :: NMAXTEL ( MAXTOTAL )
      INTEGER :: KALLTEL

!  start code
!  ----------

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  initialize compression matrix (Do this for every Fourier component)
!   Version 2.8,  Rearrange Call Arguments.

      CALL BVPTEL_MATRIX_INIT_OMP ( &
        FOURIER, NLAYERS, NSTKS_NSTRMS, DO_LAYER_SCATTERING,       &
        DO_BVTEL_INITIAL, BVTEL_FOURIER, N_BVTELMATRIX_SIZE,       &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
        ACTIVE_LAYERS, NMINTEL, NMAXTEL, KALLTEL,                  &
        BANDTELMAT2 )

!  Do the surface setup if required
!  --------------------------------

!  Verson 2.7 and earlier ---------
!  - Specialist Option 2 has Telescoping even for Fourier 0, so need
!       to have reflected solutions in lower boundary. Lambertian only

!  Version 2.8 --------------------
!  - Telescoping for higher Fourier components now with a general
!    BRDF surface condition. Formerly, this feature was absent.
!    Specialist option flag no longer needed

!  New Call for Version 2.8
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally

      NLAST = ACTIVE_LAYERS(NLAYERS_TEL)
      CALL BVPTEL_SURFACE_SETUP_HOM  ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, FOURIER,             & ! Flags
            NSTOKES, NSTREAMS, NSTKS_NSTRMS, NLAYERS, NLAST, MUELLER_INDEX, & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,        & ! Input surface
            K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_DISORDS,        & ! RTE solutions
            R2_HOMP, R2_HOMM, AXBID_F, CUMTRANS, CUMQUAD )                    ! Output

!  Old Call
!      IF ( DO_SPECIALIST_OPTION_2 ) THEN
!        IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
!          CALL BVP_SURFACE_SETUP_HOM  ( &
!            DO_INCLUDE_SURFACE, FOURIER, &
!            SURFACE_FACTOR, &
!            NSTOKES, NSTREAMS, &
!            NLAYERS, DO_LAMBERTIAN_SURFACE, &
!            ALBEDO, BRDF_F, &
!            QUAD_STRMWTS, NSTKS_NSTRMS, &
!            MUELLER_INDEX, &
!            K_REAL, K_COMPLEX, &
!            SOLA_XPOS, SOLB_XNEG, &
!            R2_HOMP, R2_HOMM, &
!            AXBID_F )
!        ENDIF
!      ENDIF

!  set up boundary values matrix in compressed form (the "A" as in AX=B)
!   Version 2.8, Remove Argument SPECIALIST OPTION. Rearrange Call Arguments.

!  New  Call

      CALL BVPTEL_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, NSTOKES, NSTREAMS, NLAYERS,              &
        NLAYERS_TEL, ACTIVE_LAYERS, NMINTEL, NMAXTEL, KALLTEL,       &
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM,        &
        BANDTELMAT2, SMAT2 )

!  Old Call
!      CALL BVPTEL_MATRIX_SETUP ( &
!        DO_INCLUDE_SURFACE, &
!        NSTOKES, NSTREAMS, &
!        NLAYERS, DO_SPECIALIST_OPTION_2, &
!        NSTREAMS_2, NSTKS_NSTRMS, &
!        NSTKS_NSTRMS_2, T_DELT_EIGEN, &
!        K_REAL, K_COMPLEX, &
!        SOLA_XPOS, SOLB_XNEG, &
!        NLAYERS_TEL, ACTIVE_LAYERS, &
!        R2_HOMP, R2_HOMM, &
!        NMINTEL, NMAXTEL, KALLTEL, &
!        BANDTELMAT2, SMAT2 )

!  SVD decomposition of compressed boundary values matrix
!   Version 2.8, Rearrange Call Arguments.

      CALL BVPTEL_MATRIX_LUD ( &
        NSTKS_NSTRMS_2, NLAYERS_TEL, N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,    &
        BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,          &
        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  return

      RETURN
      END SUBROUTINE BVPTEL_MATRIXSETUP_MASTER

!

      SUBROUTINE BVPTEL_MATRIX_INIT_OMP ( &
        FOURIER, NLAYERS, NSTKS_NSTRMS, DO_LAYER_SCATTERING,       &
        DO_BVTEL_INITIAL, BVTEL_FOURIER, N_BVTELMATRIX_SIZE,       &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, NLAYERS_TEL, &
        ACTIVE_LAYERS, NMINTEL, NMAXTEL, KALLTEL,                  &
        BANDTELMAT2 )

!  Initialise the compressed matrix
!    Version 2.8. re-arrange arguments. This is the OMP version

      USE VLIDORT_PARS_m, only : MAXMOMENTS, MAXLAYERS, MAXBANDTOTAL, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Line 1 input

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Lines 2-4 output

      LOGICAL, INTENT (INOUT) ::        DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::        BVTEL_FOURIER
      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SIZE

      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (OUT) ::          N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (OUT) ::          NLAYERS_TEL

      INTEGER, INTENT (OUT) ::          ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, INTENT (OUT) ::          NMINTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          NMAXTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::          KALLTEL

!  Line 5, the initialized matrix

      DOUBLE PRECISION, INTENT (OUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )

!  Local variables
!  ---------------

      INTEGER :: I, J, NS, N, N3, LOCAL_FOURIER

!  compression function

!      INTEGER  BTELMAT_ROWMASK
!      EXTERNAL BTELMAT_ROWMASK

!mick fix 7/29/2014 - Changed subroutine to define output vars each time to make VLIDORT 
!                     threadsafe.  Removed the saved vars and the original if block with
!                     its else section in this version of the subroutine.  Also added
!                     if block to save the original fourier component used to perform
!                     these setup calculations.

!  Set up
!  ------

!  Save fourier component used the first time subroutine is called

      IF ( DO_BVTEL_INITIAL ) THEN
        LOCAL_FOURIER = FOURIER
        BVTEL_FOURIER = FOURIER
        DO_BVTEL_INITIAL = .FALSE.
      ELSE
        LOCAL_FOURIER = BVTEL_FOURIER
      ENDIF

!  Determine active layers in atmosphere

      NS = 0
      ACTIVE_LAYERS = 0
      DO N = 1, NLAYERS
       IF ( DO_LAYER_SCATTERING(LOCAL_FOURIER,N) ) THEN
        NS = NS + 1
        ACTIVE_LAYERS(NS) = N
       ENDIF
      ENDDO
      NLAYERS_TEL = NS

!  Size of Reduced BVTEL matrix

      N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS * NLAYERS_TEL

!  Exit if only one active layer

      IF ( NLAYERS_TEL .EQ. 1 ) GO TO 345

!  Number of sub and super diagonals

      N_BVTELMATRIX_SUPDIAG = 3 * NSTKS_NSTRMS - 1
      N_BVTELMATRIX_SUBDIAG = 3 * NSTKS_NSTRMS - 1

!  Compression row indices

      DO J = 1, N_BVTELMATRIX_SUPDIAG + 1
        NMINTEL(J) = 1
      ENDDO
      DO J = N_BVTELMATRIX_SUPDIAG + 2, N_BVTELMATRIX_SIZE
        NMINTEL(J) = J - N_BVTELMATRIX_SUPDIAG
      ENDDO

      DO J = 1, N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG
        NMAXTEL(J) = J + N_BVTELMATRIX_SUBDIAG
      ENDDO
      N3 = N_BVTELMATRIX_SIZE - N_BVTELMATRIX_SUBDIAG + 1
      DO J = N3, N_BVTELMATRIX_SIZE
        NMAXTEL(J) = N_BVTELMATRIX_SIZE
      ENDDO

      KALLTEL = N_BVTELMATRIX_SUBDIAG + N_BVTELMATRIX_SUPDIAG + 1

!  Avoid fancy zeroing - adopt kludge
!   Potential Danger point

      DO I = 1, MAXBANDTOTAL
        DO J = 1, N_BVTELMATRIX_SIZE
          BANDTELMAT2(I,J) = ZERO
        ENDDO
      ENDDO

!  Control

      GO TO 346

!  Control point

 345  CONTINUE

!  Single layer setting

        N_BVTELMATRIX_SIZE = 2 * NSTKS_NSTRMS

 346  CONTINUE

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_INIT_OMP

!

      SUBROUTINE BVPTEL_SURFACE_SETUP_HOM ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, FOURIER,             & ! Flags
            NSTOKES, NSTREAMS, NSTKS_NSTRMS, NLAYERS, NLAST, MUELLER_INDEX, & ! Numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, BRDF_F,        & ! Input surface
            K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_DISORDS,        & ! RTE solutions
            R2_HOMP, R2_HOMM, AXBID_F, CUMTRANS, CUMQUAD )                    ! Output

!  This is the Lambertian or BRDF surface routine for BVP telescoping
!      --> Reflected homogeneous solutions

!  Original version of code  July 26, 2010. Using new BRDF material. RT SOLUTIONS Inc,

!  Version 2.8. August 2016. Programmed 6/27/16, by R. Spurr, RT Solutions Inc.
!     * Use performance-enhanced do-loops, and rearrange argument lists
!     * New Routine: General surface condition for any set of telescoped layers
!     * Formerly: Specialist-option 2 with Lambertian or BRDF surface. Active at ground.
!     * Replaced GOTO construction with IF block

!  1/31/21. Version 2.8.3. BRDF input array defined locally, drop MAXMOMENTS from parameter list

      USE VLIDORT_PARS_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXSTOKES, &
                                 MAXSTOKES_SQ, MAXEVALUES, MAXLAYERS, zero, one

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS, NLAST
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS

!  surface
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally, drop MAXMOMENTS dimension

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  BRDF_F ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )

!  RTE solutions

      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (OUT) :: AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (OUT) :: CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION, INTENT (OUT) :: CUMQUAD (MAXSTREAMS)

!  local variables
!  ---------------

      INTEGER ::          I,J,O1,K,KO1,K0,K1,K2,M,NL,N1
      DOUBLE PRECISION :: KMULT

!  Initialization
!  ==============

!  Zero total reflected contributions

      R2_HOMP(1:NSTREAMS, 1:NSTOKES, 1:NSTKS_NSTRMS) = ZERO
      R2_HOMM(1:NSTREAMS, 1:NSTOKES, 1:NSTKS_NSTRMS) = ZERO
      AXBID_F(1:NSTREAMS, 1:NSTREAMS, :) = ZERO
      CUMTRANS (1:NSTREAMS)              = ZERO
      CUMQUAD  (1:NSTREAMS)              = ZERO

!  Return with Zeroed values if surface flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Proxies (Fourier component,Last active layer

      M = FOURIER
      NL = NLAST

!  Real-valued separation constants Offset

      KO1 = K_REAL(NL) + 1

!  Return if Fourier component not zero for Lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .GT. 0 ) RETURN

!  Cumulative transmittance and help array
      
      CUMTRANS (1:NSTREAMS) = ONE
      do N1 = NL + 1, NLAYERS
        CUMTRANS(1:NSTREAMS) =  CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
      enddo
      CUMQUAD(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * QUAD_STRMWTS(1:NSTREAMS)

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .Eq. 0) THEN

!  Surface albedo factor

        KMULT = SURFACE_FACTOR * ALBEDO

!  Homogeneous real solutions

        DO K = 1, K_REAL(NL)
          DO I = 1, NSTREAMS
            R2_HOMP(I,1,K) = KMULT * CUMTRANS(I) * sum( CUMQUAD(1:NSTREAMS) * SOLA_XPOS(1:NSTREAMS,1,K,NL) )
            R2_HOMM(I,1,K) = KMULT * CUMTRANS(I) * sum( CUMQUAD(1:NSTREAMS) * SOLB_XNEG(1:NSTREAMS,1,K,NL) )
          ENDDO
        ENDDO

!  Homogeneous complex solutions

        KO1 = K_REAL(NL) + 1
        DO K = 1, K_COMPLEX(NL)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          DO I = 1, NSTREAMS
            R2_HOMP(I,1,K1) = KMULT * CUMTRANS(I) * sum( CUMQUAD(1:NSTREAMS) * SOLA_XPOS(1:NSTREAMS,1,K1,NL) )
            R2_HOMM(I,1,K1) = KMULT * CUMTRANS(I) * sum( CUMQUAD(1:NSTREAMS) * SOLB_XNEG(1:NSTREAMS,1,K1,NL) )
            R2_HOMP(I,1,K2) = KMULT * CUMTRANS(I) * sum( CUMQUAD(1:NSTREAMS) * SOLA_XPOS(1:NSTREAMS,1,K2,NL) )
            R2_HOMM(I,1,K2) = KMULT * CUMTRANS(I) * sum( CUMQUAD(1:NSTREAMS) * SOLB_XNEG(1:NSTREAMS,1,K2,NL) )
          ENDDO
        ENDDO

!  Total BRDF case
!  ===============

      ELSE

!  help variable
!    -- 1/31/21. Version 2.8.3. BRDF input array defined locally, drop M=FOURIER index

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            DO J = 1, NSTREAMS
              AXBID_F(I, J, MUELLER_INDEX(O1,1:NSTOKES)) = CUMQUAD(J) * BRDF_F(MUELLER_INDEX(O1,1:NSTOKES), I, J)
            ENDDO
          ENDDO
        ENDDO

!  homogeneous real solutions

        DO K = 1, K_REAL(NL)
          DO O1 = 1, NSTOKES
            DO I = 1, NSTREAMS
              R2_HOMP(I,O1,K) = SURFACE_FACTOR * CUMTRANS(I) * &
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K, NL) )
              R2_HOMM(I,O1,K) = SURFACE_FACTOR * CUMTRANS(I) * &
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K, NL) )
            ENDDO
          ENDDO
        ENDDO

!  homogeneous complex solutions

        DO K = 1, K_COMPLEX(NL)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              R2_HOMP(I,O1,K1) = SURFACE_FACTOR * CUMTRANS(I) * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K1, NL) )
              R2_HOMM(I,O1,K1) = SURFACE_FACTOR * CUMTRANS(I) * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K1, NL) )
              R2_HOMP(I,O1,K2) = SURFACE_FACTOR * CUMTRANS(I) * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLA_XPOS(1:NSTREAMS, 1:NSTOKES, K2, NL) )
              R2_HOMM(I,O1,K2) = SURFACE_FACTOR * CUMTRANS(I) * & 
                sum( AXBID_F(I, 1:NSTREAMS, MUELLER_INDEX(O1,1:NSTOKES)) * SOLB_XNEG(1:NSTREAMS, 1:NSTOKES, K2, NL) )
            ENDDO
          ENDDO
        ENDDO

!  End BRDF clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_SURFACE_SETUP_HOM

!

      SUBROUTINE BVPTEL_MATRIX_SETUP ( &
        DO_INCLUDE_SURFACE, NSTOKES, NSTREAMS, NLAYERS,              &
        NLAYERS_TEL, ACTIVE_LAYERS, NMIN, NMAX, KALL,                &
        NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, K_REAL, K_COMPLEX, &
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM,        &
        BANDTELMAT2, SMAT2 )

!  Version 2.8. August 2016. Programmed 6/27/16, by R. Spurr, RT Solutions Inc.
!     * Use performance-enhanced do-loops, and rearrange argument lists
!     * general surface treatment now implemented
!     * Replaced GOTO construction with IF block

!  Fills up the matrix directly (compressed or 1-layer)

      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXSTREAMS, MAXLAYERS, MAXEVALUES, &
                                 MAXSTREAMS_2, MAXSTRMSTKS_2, MAXBANDTOTAL, MAXTOTAL, ZERO

      IMPLICIT NONE

!  flag and control

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

!  Numbers

      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION, INTENT (IN) ::  R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  telescoping control

      INTEGER, INTENT (IN) ::           NLAYERS_TEL
      INTEGER, INTENT (IN) ::           ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           NMIN ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           NMAX ( MAXTOTAL )
      INTEGER, INTENT (IN) ::           KALL

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (OUT)   :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

!  Local variables
!  ---------------

      INTEGER ::          I, J, I1, N, N1, O1, IR, IROW, NS
      INTEGER ::          EP, EM, EPC, EPS, EMS
      INTEGER ::          CP, CM, CEP, CEM, CEPS, CEMS, CEP1, CEM1
      INTEGER ::          ROWEP, ROWEM, ROWEPS, ROWEMS, ROWCEP, ROWCEM
      INTEGER ::          ROWCEPS, ROWCEMS, ROWCEP1, ROWCEM1
      DOUBLE PRECISION :: XPNET, XMNET, X1R, X2R, X1I, X2I, X1R_H, X1I_H
      INTEGER ::          KO1, K0, K1, K2, KO11
      INTEGER ::          CE_OFFSET, C0

!  compression function

!      INTEGER  BTELMAT_ROWMASK
!      EXTERNAL BTELMAT_ROWMASK

!  Initialization

      EP = 0

!  General Case NLAYERS_TEL > 1
!  ============================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  top BC for first active layer 1: no downward diffuse radiation
!  --------------------------------------------------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)
        KO1 = K_REAL(N) + 1
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1

!  real solutions

            DO EP = 1, K_REAL(N)
              EM = EP + NSTKS_NSTRMS
              ROWEP = BTELMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
              ROWEM = BTELMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
              BANDTELMAT2(ROWEP,EP)  = SOLA_XPOS(I,O1,EP,N)
              BANDTELMAT2(ROWEM,EM)  = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
            ENDDO

!  complex solutions (REWORKED)

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1

              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              EMS = EPS + NSTKS_NSTRMS

              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

              ROWEP  = BTELMAT_ROWMASK(IROW,EP,NMIN(EP),NMAX(EP),KALL)
              ROWEM  = BTELMAT_ROWMASK(IROW,EM,NMIN(EM),NMAX(EM),KALL)
              ROWEPS = BTELMAT_ROWMASK(IROW,EPS,NMIN(EPS),NMAX(EPS),KALL)
              ROWEMS = BTELMAT_ROWMASK(IROW,EMS,NMIN(EMS),NMAX(EMS),KALL)

              BANDTELMAT2(ROWEP,EP)   =   X1R
              BANDTELMAT2(ROWEPS,EPS) = - X1I
              BANDTELMAT2(ROWEM,EM)   =   X2R
              BANDTELMAT2(ROWEMS,EMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  intermediate layer boundaries
!  -----------------------------

        C0 = - NSTKS_NSTRMS
        DO NS = 2, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          N1 = N - 1
          KO1  = K_REAL(N) + 1
          KO11 = K_REAL(N1) + 1
          C0   = C0 + NSTKS_NSTRMS_2
          CE_OFFSET = C0 - NSTKS_NSTRMS
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW

!  real solutions for layer above the boundary

              DO EP = 1, K_REAL(N1)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS

                ROWCEP = BTELMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM = BTELMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)

                BANDTELMAT2(ROWCEP,CEP)   = T_DELT_EIGEN(EP,N1)*SOLA_XPOS(I,O1,EP,N1)
                BANDTELMAT2(ROWCEM,CEM)   = SOLB_XNEG(I,O1,EP,N1)
              ENDDO

!  real solutions for layer below the boundary
!   ( Note the change of sign !!! )

              DO EP = 1, K_REAL(N)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS
                CEP1 = CEP + NSTKS_NSTRMS_2
                CEM1 = CEM + NSTKS_NSTRMS_2

                ROWCEP1 = BTELMAT_ROWMASK(CM,CEP1,NMIN(CEP1),NMAX(CEP1),KALL)
                ROWCEM1 = BTELMAT_ROWMASK(CM,CEM1,NMIN(CEM1),NMAX(CEM1),KALL)

                BANDTELMAT2(ROWCEP1,CEP1) = - SOLA_XPOS(I,O1,EP,N)
                BANDTELMAT2(ROWCEM1,CEM1) = - T_DELT_EIGEN(EP,N)*SOLB_XNEG(I,O1,EP,N)
              ENDDO

!  complex solutions for layer above boundary

              DO EPC = 1, K_COMPLEX(N1)
                K0 = 2*EPC - 2
                K1 = KO11 + K0
                K2 = K1  + 1

                EP  = K_REAL(N1) + EPC
                EPS = K_REAL(N1) + K_COMPLEX(N1) + EPC
                CEP  = CE_OFFSET + EP
                CEM  = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS
                CEMS = CEPS + NSTKS_NSTRMS

                X1R = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K1,N1) - SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K2,N1)
                X1I = SOLA_XPOS(I,O1,K1,N1)*T_DELT_EIGEN(K2,N1) + SOLA_XPOS(I,O1,K2,N1)*T_DELT_EIGEN(K1,N1)
                X2R = SOLB_XNEG(I,O1,K1,N1)
                X2I = SOLB_XNEG(I,O1,K2,N1)

                ROWCEP  = BTELMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BTELMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BTELMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BTELMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDTELMAT2(ROWCEP,CEP)   =   X1R
                BANDTELMAT2(ROWCEPS,CEPS) = - X1I
                BANDTELMAT2(ROWCEM,CEM)   =   X2R
                BANDTELMAT2(ROWCEMS,CEMS) = - X2I
              ENDDO

!  complex solutions for layer below boundary
!   ( Note the change of sign !!! )

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                CEP  = CE_OFFSET + EP + NSTKS_NSTRMS_2
                CEM  = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS + NSTKS_NSTRMS_2
                CEMS = CEPS + NSTKS_NSTRMS

                X1R = SOLA_XPOS(I,O1,K1,N)
                X1I = SOLA_XPOS(I,O1,K2,N)
                X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
                X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)

                ROWCEP  = BTELMAT_ROWMASK(CM,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BTELMAT_ROWMASK(CM,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BTELMAT_ROWMASK(CM,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BTELMAT_ROWMASK(CM,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDTELMAT2(ROWCEP,CEP)   = - X1R
                BANDTELMAT2(ROWCEPS,CEPS) =   X1I
                BANDTELMAT2(ROWCEM,CEM)   = - X2R
                BANDTELMAT2(ROWCEMS,CEMS) =   X2I
              ENDDO

            ENDDO
          ENDDO
        ENDDO

!  bottom BC for Lowest active layer
!  ---------------------------------

!   Old code : Normally no surface additions, except Specialist # 2
!   2.8: General surface condition for telescoped-BVP, 6/27/16 

        N = ACTIVE_LAYERS(NLAYERS_TEL)
        KO1 = K_REAL(N) + 1
        C0  = C0 + NSTKS_NSTRMS_2
        CE_OFFSET = C0 - NSTKS_NSTRMS

!  If there is a surface contribution

!      IF ( DO_INCLUDE_SURFACE.AND.DO_SPECIALIST_OPTION_2 ) THEN -----> OLD CODE
        IF ( DO_INCLUDE_SURFACE ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = C0 + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS

                XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
                XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)

                ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

                BANDTELMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                BANDTELMAT2(ROWCEM,CEM) = XMNET
              ENDDO

!  Complex solutions

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS
                CEMS = CEPS + NSTKS_NSTRMS

                X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
                X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
                X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
                X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)

                X1R = X1R_H * T_DELT_EIGEN(K1,N) - X1I_H * T_DELT_EIGEN(K2,N)
                X1I = X1R_H * T_DELT_EIGEN(K2,N) + X1I_H * T_DELT_EIGEN(K1,N)

                ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BTELMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BTELMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDTELMAT2(ROWCEP,CEP)   =   X1R
                BANDTELMAT2(ROWCEPS,CEPS) = - X1I
                BANDTELMAT2(ROWCEM,CEM)   =   X2R
                BANDTELMAT2(ROWCEMS,CEMS) = - X2I
              ENDDO

!  End stokes and streams loops

            ENDDO
          ENDDO

!  No surface reflections

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = C0 + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = CE_OFFSET + EP
                CEM = CEP + NSTKS_NSTRMS

                XPNET = SOLA_XPOS(I1,O1,EP,N)
                XMNET = SOLB_XNEG(I1,O1,EP,N)

                ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)

                BANDTELMAT2(ROWCEP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                BANDTELMAT2(ROWCEM,CEM) = XMNET
              ENDDO

!  Complex solutions

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1

                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                CEP  = CE_OFFSET + EP
                CEM  = CEP + NSTKS_NSTRMS
                CEPS = CE_OFFSET + EPS
                CEMS = CEPS + NSTKS_NSTRMS

                X1R = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K1,N) - SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K2,N)
                X1I = SOLA_XPOS(I1,O1,K1,N) * T_DELT_EIGEN(K2,N) + SOLA_XPOS(I1,O1,K2,N) * T_DELT_EIGEN(K1,N)
                X2R = SOLB_XNEG(I1,O1,K1,N)
                X2I = SOLB_XNEG(I1,O1,K2,N)

                ROWCEP  = BTELMAT_ROWMASK(CP,CEP,NMIN(CEP),NMAX(CEP),KALL)
                ROWCEM  = BTELMAT_ROWMASK(CP,CEM,NMIN(CEM),NMAX(CEM),KALL)
                ROWCEPS = BTELMAT_ROWMASK(CP,CEPS,NMIN(CEPS),NMAX(CEPS),KALL)
                ROWCEMS = BTELMAT_ROWMASK(CP,CEMS,NMIN(CEMS),NMAX(CEMS),KALL)

                BANDTELMAT2(ROWCEP,CEP)   =   X1R
                BANDTELMAT2(ROWCEPS,CEPS) = - X1I
                BANDTELMAT2(ROWCEM,CEM)   =   X2R
                BANDTELMAT2(ROWCEMS,CEMS) = - X2I
              ENDDO

            ENDDO
          ENDDO

!  End surface Yes/No clause

        ENDIF

!  special case. Only 1 active layer
!  =================================

      ELSE


        N = ACTIVE_LAYERS(1)
        KO1 = K_REAL(N) + 1

!  initialize using the SMAT2 matrix

        DO I = 1, NSTKS_NSTRMS_2
          DO J = 1, NSTKS_NSTRMS_2
            SMAT2(I,J) = ZERO
          ENDDO
        ENDDO

!  top of the layer (downwelling only)
!  -----------------------------------

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1

!  real solutions

            DO EP = 1, K_REAL(N)
              EM = EP + NSTKS_NSTRMS
              SMAT2(IROW,EP) = SOLA_XPOS(I,O1,EP,N)
              SMAT2(IROW,EM) = SOLB_XNEG(I,O1,EP,N)*T_DELT_EIGEN(EP,N)
            ENDDO

!  REWORKED complex solutions

            DO EPC = 1, K_COMPLEX(N)
              K0 = 2*EPC - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              EP  = K_REAL(N) + EPC
              EPS = K_REAL(N) + K_COMPLEX(N) + EPC
              EM  = EP  + NSTKS_NSTRMS
              EMS = EPS + NSTKS_NSTRMS
              X1R = SOLA_XPOS(I,O1,K1,N)
              X1I = SOLA_XPOS(I,O1,K2,N)
              X2R = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K2,N)
              X2I = SOLB_XNEG(I,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLB_XNEG(I,O1,K2,N)*T_DELT_EIGEN(K1,N)
              SMAT2(IROW,EP)  =   X1R
              SMAT2(IROW,EPS) = - X1I
              SMAT2(IROW,EM)  =   X2R
              SMAT2(IROW,EMS) = - X2I
            ENDDO

          ENDDO
        ENDDO

!  bottom BC, with Surface
!  -----------------------

!   Old Code: Only albedo additions for specialist option 2
!   2.8: General surface condition for telescoped-BVP, 6/27/16 

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = NSTKS_NSTRMS + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = EP
                CEM = CEP + NSTKS_NSTRMS
                XPNET = SOLA_XPOS(I1,O1,EP,N) - R2_HOMP(I,O1,EP)
                XMNET = SOLB_XNEG(I1,O1,EP,N) - R2_HOMM(I,O1,EP)
                SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                SMAT2(CP,CEM) = XMNET
              ENDDO

!  Complex solutions

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                EM  = EP  + NSTKS_NSTRMS
                CEP = EP
                CEM = EM
                CEPS = EPS
                CEMS = CEPS + NSTKS_NSTRMS
                X1R_H = SOLA_XPOS(I1,O1,K1,N) - R2_HOMP(I,O1,K1)
                X1I_H = SOLA_XPOS(I1,O1,K2,N) - R2_HOMP(I,O1,K2)
                X2R = SOLB_XNEG(I1,O1,K1,N) - R2_HOMM(I,O1,K1)
                X2I = SOLB_XNEG(I1,O1,K2,N) - R2_HOMM(I,O1,K2)
                X1R = X1R_H * T_DELT_EIGEN(K1,N) - X1I_H * T_DELT_EIGEN(K2,N)
                X1I = X1R_H * T_DELT_EIGEN(K2,N) + X1I_H * T_DELT_EIGEN(K1,N)
                SMAT2(CP,CEP)  =   X1R
                SMAT2(CP,CEPS) = - X1I
                SMAT2(CP,CEM)  =   X2R
                SMAT2(CP,CEMS) = - X2I
              ENDDO

            ENDDO
          ENDDO

!  bottom BC (with no surface)
!  ---------------------------

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CP = NSTKS_NSTRMS + IROW

!  real solutions

              DO EP = 1, K_REAL(N)
                CEP = EP
                CEM = CEP + NSTKS_NSTRMS
                XPNET = SOLA_XPOS(I1,O1,EP,N)
                XMNET = SOLB_XNEG(I1,O1,EP,N)
                SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
                SMAT2(CP,CEM) = XMNET
              ENDDO

!  REWORKED Complex solutions
!    Note to self. The second set of solutions seems correct.

              DO EPC = 1, K_COMPLEX(N)
                K0 = 2*EPC - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                EP  = K_REAL(N) + EPC
                EPS = K_REAL(N) + K_COMPLEX(N) + EPC
                EM  = EP  + NSTKS_NSTRMS
                CEP = EP
                CEM = EM
                CEPS = EPS
                CEMS = CEPS + NSTKS_NSTRMS
                X1R = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K1,N) - SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K2,N)
                X1I = SOLA_XPOS(I1,O1,K1,N)*T_DELT_EIGEN(K2,N) + SOLA_XPOS(I1,O1,K2,N)*T_DELT_EIGEN(K1,N)
                X2R = SOLB_XNEG(I1,O1,K1,N)
                X2I = SOLB_XNEG(I1,O1,K2,N)
                IF ( N.EQ.NLAYERS ) THEN
                  SMAT2(CP,CEP)  =   X1R   ! apparently wrong
                  SMAT2(CP,CEPS) = - X1I   ! apparently wrong
                  SMAT2(CP,CEM)  =   X2R   ! apparently wrong
                  SMAT2(CP,CEMS) = - X2I   ! apparently wrong
                ELSE
                  SMAT2(CP,CEP)  =   X1R  ! tested 29 December 2005
                  SMAT2(CP,CEPS) = - X1I  ! tested 29 December 2005
                  SMAT2(CP,CEM)  =   X2R  ! tested 29 December 2005
                  SMAT2(CP,CEMS) = - X2I  ! tested 29 December 2005
                ENDIF
              ENDDO

            ENDDO
          ENDDO

!  End surface vs. No surface clause

        ENDIF

!  End NLAYERS_TEL > 1 vs. NLAYERS_TEL = 1

      ENDIF

!  normal return and finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_SETUP

!

      SUBROUTINE BVPTEL_MATRIX_LUD ( &
        NSTKS_NSTRMS_2, NLAYERS_TEL, N_BVTELMATRIX_SIZE, &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,    &
        BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,          &
        STATUS, MESSAGE, TRACE )

! TELESCOPED boundary value problem SVD decomposition.
!   Version 2.8, Rearrange Call Arguments.

      USE VLIDORT_PARS_m, only : MAXSTRMSTKS_2, MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS
      USE LAPACK_TOOLS_m, only : DGBTRF, DGETRF

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::             NLAYERS_TEL

      DOUBLE PRECISION, INTENT (INOUT) :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      DOUBLE PRECISION, INTENT (INOUT) :: SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )

      INTEGER, INTENT (OUT) ::            IPIVOTTEL ( MAXTOTAL )
      INTEGER, INTENT (OUT) ::            SIPIVOT ( MAXSTRMSTKS_2 )
      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER           :: INFO
      CHARACTER (LEN=3) :: CI

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SVD the BVPTEL matrix: With compression (multilayers)
!  ----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF &
           ( N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SIZE, &
             N_BVTELMATRIX_SUBDIAG, N_BVTELMATRIX_SUPDIAG, &
             BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_LUD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_LUD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK LU-decomposition for single layer matrix

        CALL DGETRF &
           (  NSTKS_NSTRMS_2, NSTKS_NSTRMS_2, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF (nlayers_tel=1)call in BVPTEL_MATRIX_LUD'
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_MATRIX_LUD

!

      SUBROUTINE BVPTEL_SOLUTION_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,               & ! Flags
        DO_INCLUDE_DIRECTBEAM, FOURIER, IBEAM, NSTOKES, NSTREAMS,& ! Numbers
        NLAYERS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,       & ! Numbers
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,          & ! BVP Tel Control
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,            & ! BVP Tel Control
        MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR,        & ! Input Bookkeeping
        ALBEDO, AXBID_F, CUMTRANS, CUMQUAD,           & ! Input Surface inputs
        WLOWER, WUPPER, DIRECT_BEAM,                             & ! Input RTE stuff
        SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_DELT_DISORDS,      & ! Input RTE stuff
        BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT, COLTEL2, SCOL2,  & ! Input BVProblem
        R2_BEAM, LCON, MCON, STATUS, MESSAGE, TRACE )              ! Output and Status

!  Version 2.8. August 2016.
!     * Use performance-enhanced do-loops
!     * Rearrange argument lists.
!     * Replaced GOTO constructions with IF blocks
!     * Introduce general treatment of surfaces

!  1/31/21. Version 2.8.3. SLEAVE input array defined locally, each Fourier
!    -- NOTE TO SELF, INTRODUCE Higher-order SLEAVE Fourier components

      USE VLIDORT_PARS_m, only : MAXSTREAMS, MAXSTOKES, MAXBEAMS, MAXLAYERS, MAXEVALUES,   &
                                 MAXSTOKES_SQ, MAXSTREAMS_2, MAXSTRMSTKS, MAXSTRMSTKS_2,   &
                                 MAXBANDTOTAL, MAXTOTAL, VLIDORT_SUCCESS, VLIDORT_SERIOUS

      IMPLICIT NONE

!  Flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM
!      LOGICAL, INTENT (IN) ::           DO_SPECIALIST_OPTION_2  ! OLD CODE

!  Numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS

      INTEGER, INTENT (IN) ::           NSTREAMS_2
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2

!  BVP control

      INTEGER, INTENT (IN) ::           NLAYERS_TEL
      INTEGER, INTENT (IN) ::           ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SUBDIAG

!  Surface

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )
      DOUBLE PRECISION, INTENT (IN) ::  CUMQUAD (MAXSTREAMS)
      DOUBLE PRECISION, INTENT (IN) ::  CUMTRANS(MAXSTREAMS)

!  Bookkeeping

      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           K_COMPLEX ( MAXLAYERS )

!  RTE inputs

      DOUBLE PRECISION, INTENT (IN) ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN   ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  BVP Matrix inputs

!      INTEGER, INTENT (IN) ::           IPIVOT ( MAXTOTAL )
      DOUBLE PRECISION, INTENT (IN) ::  SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::           SIPIVOT ( MAXSTRMSTKS_2 )
      DOUBLE PRECISION, INTENT (IN) ::  BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::           IPIVOTTEL ( MAXTOTAL )

!  Column output

!mick fix 7/29/2014 - placed COLTEL2 & SCOL2 in call to make VLIDORT threadsafe
      DOUBLE PRECISION, INTENT (INOUT) :: COLTEL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

      INTEGER ::          STATUS_SUB, NLAST

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Do the surface setup if required
!  --------------------------------

!  - Specialist Option 2 has Telescoping even for Fourier 0, so need
!       to have reflected solutions in the lower boundary, Foremrly
!       Lambertian only , but now BRDF is included
!  - Telescoping for higher Fourier components now with a general
!    BRDF surface condition. Formerly, this feature was absent.

!  New  Call for Version 2.8

      NLAST = ACTIVE_LAYERS(NLAYERS_TEL)
      CALL BVPTEL_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   &  
        FOURIER, NSTOKES, NSTREAMS, NLAST,           &
        SURFACE_FACTOR, ALBEDO, AXBID_F,  &
        MUELLER_INDEX, CUMTRANS, CUMQUAD, WLOWER,    &
        R2_BEAM )

!  Old Code
!      IF ( DO_SPECIALIST_OPTION_2 ) THEN
!        IF ( ACTIVE_LAYERS(NLAYERS_TEL).EQ.NLAYERS ) THEN
!      CALL BVPTEL_SURFACE_SETUP_BEAM ( &
!            DO_INCLUDE_SURFACE, FOURIER, &
!            SURFACE_FACTOR, &
!            NSTOKES, NSTREAMS, &
!            NLAYERS, DO_LAMBERTIAN_SURFACE, &
!            ALBEDO, &
!            QUAD_STRMWTS, MUELLER_INDEX, &
!            WLOWER, AXBID_F, &
!            R2_BEAM )
!        ENDIF
!      ENDIF

!  --set up Column for solution vector (the "B" as in AX=B)
!   Version 2.8, Remove Argument SPECIALIST OPTION. Rearrange Call Arguments.

!  New Version 2.8  Call

      CALL BVPTEL_COLUMN_SETUP                                 &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,         & ! input
            IBEAM, NSTOKES, NSTREAMS, NLAYERS, NSTREAMS_2,     & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
            NSTKS_NSTRMS, NSTKS_NSTRMS_2, T_DELT_DISORDS,      & ! input
            WLOWER, WUPPER, R2_BEAM, DIRECT_BEAM,              & ! Input
            COLTEL2, SCOL2 )                                     ! Output

!  Old Code
!      CALL BVPTEL_COLUMN_SETUP ( &
!        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
!        IBEAM, &
!        NSTOKES, NSTREAMS, &
!        NLAYERS, DO_SPECIALIST_OPTION_2, &
!        NSTREAMS_2, NSTKS_NSTRMS, &
!        NSTKS_NSTRMS_2, &
!        R2_BEAM, DIRECT_BEAM, &
!        WUPPER, WLOWER, &
!        N_BVTELMATRIX_SIZE, NLAYERS_TEL, &
!        ACTIVE_LAYERS, &
!        SCOL2, COLTEL2 )

!  --Solve the boundary problem for this Fourier component (back substit
!   Version 2.8, Rearrange Call Arguments.

!  New Call

      CALL BVPTEL_BACKSUB ( &
        DO_INCLUDE_SURFACE, FOURIER, IBEAM,                       &
        NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,           &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,             &
        K_REAL, K_COMPLEX, T_DELT_DISORDS,                        &
        WUPPER, WLOWER, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,       &
        BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, COLTEL2, SCOL2,   &
        LCON, MCON, STATUS_SUB, MESSAGE, TRACE )

!  Old Code
!      CALL BVPTEL_BACKSUB ( &
!        IBEAM, NSTOKES, &
!        NSTREAMS, NLAYERS, &
!        NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
!        T_DELT_DISORDS, T_DELT_EIGEN, &
!        K_REAL, K_COMPLEX, &
!        SOLA_XPOS, SOLB_XNEG, &
!        WUPPER, WLOWER, &
!        IPIVOT, &
!        SMAT2, SIPIVOT, &
!        N_BVTELMATRIX_SIZE, &
!        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG, &
!        NLAYERS_TEL, ACTIVE_LAYERS, &
!        BANDTELMAT2, IPIVOTTEL, &
!        COLTEL2, &
!        SCOL2, LCON, MCON, &
!        STATUS_SUB, MESSAGE, TRACE )

!  error tracing

      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        RETURN
      ENDIF

!  return

      RETURN
      END SUBROUTINE BVPTEL_SOLUTION_MASTER

!

      SUBROUTINE BVPTEL_SURFACE_SETUP_BEAM ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   &  
        FOURIER, NSTOKES, NSTREAMS, NLAST,           &
        SURFACE_FACTOR, ALBEDO, AXBID_F,  &
        MUELLER_INDEX, CUMTRANS, CUMQUAD, WLOWER,    &
        R2_BEAM )

!  This is the Lambertian or BRDF surface routine for BVP telescoping
!      --> Reflected beam solutions

!  Original version of code  July 26, 2010. Using new BRDF material. RT SOLUTIONS Inc,

!  Version 2.8. Programmed 6/27/16, by R. Spurr, RT Solutions Inc.
!    - General surface condition for any set of telescoped layers. BVP
!    - Formerly: Specialist-option 2 with Lambertian or BRDF surface. Active at ground.
!    - Replaced GOTO constructions with IF blocks

      USE VLIDORT_PARS_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXSTOKES, &
                                 MAXSTOKES_SQ, MAXLAYERS, zero, one


      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_LAMBERTIAN_SURFACE

!  numbers

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAST

!  surface

      DOUBLE PRECISION, INTENT (IN) ::  SURFACE_FACTOR
      DOUBLE PRECISION, INTENT (IN) ::  ALBEDO
      DOUBLE PRECISION, INTENT (IN) ::  AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  RTE solutions

      INTEGER, INTENT (IN) ::           MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  CUMQUAD (MAXSTREAMS)
      DOUBLE PRECISION, INTENT (IN) ::  CUMTRANS(MAXSTREAMS)
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER  ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  local variables
!  ---------------

      DOUBLE PRECISION :: REFL_B, REFL_B_S
      INTEGER          :: I, J, O1, O2, OM

!  Initialization
!  ==============

!  Zero total reflected contributions

      R2_BEAM(1:NSTREAMS,1:NSTOKES)    = ZERO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Return if not the Fourier m = 0 component (Lambertian surface)

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .GT. 0 ) RETURN

!  Lambertian surface
!  ==================

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .Eq. 0) THEN

!  Integrate Downward streams of particular solutions, set reflectance

        REFL_B = DOT_PRODUCT(CUMQUAD(1:NSTREAMS),WLOWER(1:NSTREAMS,1,NLAST))
        REFL_B = REFL_B * SURFACE_FACTOR * ALBEDO
        R2_BEAM(1:NSTREAMS,1) = REFL_B * CUMTRANS(1:NSTREAMS)

!  BRDF surface
!  ------------

      ELSE

!  Integrate Downward streams of particular solutions

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B_S = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                REFL_B_S = REFL_B_S + AXBID_F(I,J,OM) * WLOWER(J,O2,NLAST)
              ENDDO
              REFL_B = REFL_B + REFL_B_S
            ENDDO
           R2_BEAM(I,O1) = SURFACE_FACTOR * CUMTRANS(I) * REFL_B
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_SURFACE_SETUP_BEAM

!

      SUBROUTINE BVPTEL_COLUMN_SETUP &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,         & ! input
            IBEAM, NSTOKES, NSTREAMS, NLAYERS, NSTREAMS_2,     & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
            NSTKS_NSTRMS, NSTKS_NSTRMS_2, T_DELT_DISORDS,      & ! input
            WLOWER, WUPPER, R2_BEAM, DIRECT_BEAM,              & ! Input
            COLTEL2, SCOL2 )                                     ! Output

!  Sets up the telescoped boundary value problem, RHS vector

!  Version 2.8. Programmed 6/27/16, by R. Spurr, RT Solutions Inc.
!    - General surface condition for any set of telescoped layers. BVP
!    - Formerly: Specialist-option 2 with Lambertian or BRDF surface. Active at ground.

      USE VLIDORT_PARS_m, only : MAXSTREAMS, MAXBEAMS, MAXSTOKES, MAXLAYERS, &
                                 MAXSTREAMS_2, MAXSTRMSTKS_2, MAXTOTAL, ZERO, ONE

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (IN) ::           DO_INCLUDE_DIRECTBEAM

!  Numbers

      INTEGER, INTENT (IN) ::           IBEAM
      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           NSTREAMS_2

!  BVP control

      INTEGER, INTENT (IN) ::           N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::           NLAYERS_TEL
      INTEGER, INTENT (IN) ::           ACTIVE_LAYERS ( MAXLAYERS )

!  Misc

      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::           NSTKS_NSTRMS_2
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )

!  RTE stuff

      DOUBLE PRECISION, INTENT (IN) ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  output

      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: COLTEL2 ( MAXTOTAL, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER ::          I, I1, N, N1, NS, C0, CM
      INTEGER ::          IR, O1, IROW
      DOUBLE PRECISION :: CUMTRANS ( MAXSTREAMS )

!  General Case, Number of layers > 1
!  ==================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  zero column vector

         COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM) = ZERO

!  Upper boundary for first active layer: no downward diffuse radiation

        NS = 1
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            COLTEL2(IROW,IBEAM)   = - WUPPER(I,O1,N)
          ENDDO
        ENDDO

!  intermediate layer boundaries

        C0 = - NSTKS_NSTRMS
        DO NS = 2, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          N1 = N - 1
          C0 = C0 + NSTKS_NSTRMS_2
          DO I = 1, NSTREAMS_2
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COLTEL2(CM,IBEAM) = WUPPER(I,O1,N) - WLOWER(I,O1,N1)
            ENDDO
          ENDDO
        ENDDO

!  Lower boundary for last active layer NLAYERS_TEL:
!   Older code only had Lambertian surface for Specialist Option 2
!   Version 2.8. Now general surface contribution for Telescoped BVP

        NS = NLAYERS_TEL
        N  = ACTIVE_LAYERS(NS)
        C0 = C0 + NSTKS_NSTRMS_2

!  This is the basic term, excludes reflected diffuse and direct-beam terms

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS ; IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1 ; CM = C0 + IROW
            COLTEL2(CM,IBEAM) = - WLOWER(I1,O1,N)
          ENDDO
        ENDDO

!  reflected contribution with surface

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS ; IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1 ; CM = C0 + IROW
              COLTEL2(CM,IBEAM) = COLTEL2(CM,IBEAM)  + R2_BEAM(I,O1)
            ENDDO
          ENDDO
        ENDIF

!  Add direct beam solution. CUMTRANS calculation is new, R. Spurr 06/27/16 (from LIDORT)
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
          CUMTRANS(1:NSTREAMS) = ONE
          DO N1 = NLAYERS, N+1, -1
            CUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
          ENDDO
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              COLTEL2(CM,IBEAM) = COLTEL2(CM,IBEAM) + CUMTRANS(I) * DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF

!  Special case, for single layer setup
!  ====================================

      ELSE

!  active layer for telescoped BVP

        N = ACTIVE_LAYERS(1)

!  initialize using the SCOL2 matrix

        SCOL2(1:NSTKS_NSTRMS_2,IBEAM) = ZERO

!  Upper boundary for layer (downwelling only)

        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            SCOL2(IROW,IBEAM)   = - WUPPER(I,O1,N)
          ENDDO
        ENDDO

!  Lower boundary  (upwelling only) for single active layer:
!   Older code only had Lambertian surface for Specialist Option 2
!   Version 2.8. Now general surface contribution for Telescoped BVP

        C0 = NSTKS_NSTRMS

!    This is the basic term, excludes reflected diffuse and direct-beam terms

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM = C0 + IROW
            SCOL2(CM,IBEAM) = - WLOWER(I1,O1,N) + R2_BEAM(I,O1)
          ENDDO
        ENDDO

!    Add diffuse reflected contribution

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = - WLOWER(I1,O1,N) + R2_BEAM(I,O1)
            ENDDO
          ENDDO
        ENDIF

!  Add direct beam solution. CUMTRANS calculation is new, R. Spurr 06/27/16 (from LIDORT)
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
          CUMTRANS(1:NSTREAMS) = ONE
          DO N1 = NLAYERS, N+1, -1
            CUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
          ENDDO
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM = C0 + IROW
              SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + DIRECT_BEAM(I,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF

!  done

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_COLUMN_SETUP

!

      SUBROUTINE BVPTEL_BACKSUB ( &
        DO_INCLUDE_SURFACE, FOURIER, IBEAM,                       &
        NSTOKES, NSTREAMS, NLAYERS, NSTKS_NSTRMS, NSTKS_NSTRMS_2, &
        N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,           &
        N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,             &
        K_REAL, K_COMPLEX, T_DELT_DISORDS,                        &
        WUPPER, WLOWER, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,       &
        BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, COLTEL2, SCOL2,   &
        LCON, MCON, STATUS, MESSAGE, TRACE )

!  Solves the telescoped boundary value problem.
!   Older code only had Lambertian surface for Specialist Option 2

!  Version 2.8. Programmed 6/27/16, by R. Spurr, RT Solutions Inc.
!    - General surface condition for any set of telescoped layers. BVP
!    - Formerly: Specialist-option 2 with Lambertian or BRDF surface. Active at ground.


      USE VLIDORT_PARS_m, only : MAXSTOKES, MAXLAYERS, MAXBEAMS, MAXSTREAMS, MAXSTRMSTKS,         &
                                 MAXSTREAMS_2, MAXEVALUES, MAXSTRMSTKS_2, MAXBANDTOTAL, MAXTOTAL, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, ZERO
      USE LAPACK_TOOLS_m, only : DGBTRS, DGETRS

      IMPLICIT NONE

!  Flag and indices

      LOGICAL, intent(in)  ::             DO_INCLUDE_SURFACE
      INTEGER, INTENT (IN) ::             IBEAM
      INTEGER, INTENT (IN) ::             FOURIER

!  Numbers

      INTEGER, INTENT (IN) ::             NSTOKES
      INTEGER, INTENT (IN) ::             NSTREAMS
      INTEGER, INTENT (IN) ::             NLAYERS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::             NSTKS_NSTRMS_2

!  BVP telescoping control

      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SIZE
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUPDIAG
      INTEGER, INTENT (IN) ::             N_BVTELMATRIX_SUBDIAG
      INTEGER, INTENT (IN) ::             NLAYERS_TEL
      INTEGER, INTENT (IN) ::             ACTIVE_LAYERS ( MAXLAYERS )

!  Misc

      DOUBLE PRECISION, INTENT (IN) ::    T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_REAL ( MAXLAYERS )
      INTEGER, INTENT (IN) ::             K_COMPLEX ( MAXLAYERS )

!  RTE solutions

      DOUBLE PRECISION, INTENT (IN) ::    T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::    WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  BVP matrices and columns
!      INTEGER, INTENT (IN) ::             IPIVOT ( MAXTOTAL )

      DOUBLE PRECISION, INTENT (IN) ::    SMAT2 ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER, INTENT (IN) ::             SIPIVOT ( MAXSTRMSTKS_2 )
      DOUBLE PRECISION, INTENT (IN) ::    BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::             IPIVOTTEL ( MAXTOTAL )

!  Results

      DOUBLE PRECISION, INTENT (INOUT) :: COLTEL2 ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )
      DOUBLE PRECISION, INTENT (INOUT) :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Error handling

      INTEGER, INTENT (OUT) ::            STATUS
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) ::  TRACE

!  local variables
!  ---------------

      INTEGER ::           K, KO1, K0, K1, K2
      INTEGER ::           M, I, I1, N, NS, N1, NAL, O1, IC, ICOW, INFO
      INTEGER ::           C0, IR, IROW, IROW1, IROW_S, IROW1_S
      DOUBLE PRECISION ::  SPAR, SHOM, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION ::  SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION ::  LXR, MXR, LXR_CR, LXR_CI, MXR_CR, MXR_CI
      CHARACTER (LEN=3) :: CI
      CHARACTER (LEN=2) :: CB

      DOUBLE PRECISION  :: LOC_COLTEL2 ( MAXTOTAL, 1 )
      DOUBLE PRECISION  :: LOC_SCOL2   ( MAXSTRMSTKS_2, 1 )

!  start code
!  ----------

!  debug proxy only

      M = FOURIER

!  Intialize Exception handling

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Back-substitution for multi-layer BVP TEL
!  =========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK substitution using RHS column vector COLTEL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("COLTEL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended
        LOC_COLTEL2(1:N_BVTELMATRIX_SIZE,1) = COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM)
        CALL DGBTRS &
           ( 'n', N_BVTELMATRIX_SIZE, &
              N_BVTELMATRIX_SUBDIAG, N_BVTELMATRIX_SUPDIAG, 1, &
              BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
              LOC_COLTEL2, MAXTOTAL, INFO )
        COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM) = LOC_COLTEL2(1:N_BVTELMATRIX_SIZE,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set integration constants for active layers

        C0 = - NSTKS_NSTRMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTKS_NSTRMS_2
          KO1 = K_REAL(N) + 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            LCON(K,N) = COLTEL2(C0+IROW,IBEAM)
            MCON(K,N) = COLTEL2(C0+IROW1,IBEAM)
          ENDDO
          DO K = 1, K_COMPLEX(N)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            LCON(K1,N) = COLTEL2(C0+IROW,   IBEAM)
            LCON(K2,N) = COLTEL2(C0+IROW_S, IBEAM)
            MCON(K1,N) = COLTEL2(C0+IROW1,  IBEAM)
            MCON(K2,N) = COLTEL2(C0+IROW1_S,IBEAM)
          ENDDO
        ENDDO

!  Solve the boundary problem: Single Layer only
!  =============================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  Active layer

        N = ACTIVE_LAYERS(1)
        KO1 = K_REAL(N) + 1

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NSTKS_NSTRMS_2,1) = SCOL2(1:NSTKS_NSTRMS_2,IBEAM)
        CALL DGETRS &
           ( 'N', NSTKS_NSTRMS_2, 1, &
              SMAT2, MAXSTRMSTKS_2, SIPIVOT, &
              LOC_SCOL2, MAXSTRMSTKS_2, INFO )
        SCOL2(1:NSTKS_NSTRMS_2,IBEAM) = LOC_SCOL2(1:NSTKS_NSTRMS_2,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = VLIDORT_SERIOUS
          RETURN
        ENDIF

!  set real constants from the solution vector

        DO K = 1, K_REAL(N)
          IROW  = K
          IROW1 = IROW + NSTKS_NSTRMS
          LCON(K,N) = SCOL2(IROW,IBEAM)
          MCON(K,N) = SCOL2(IROW1,IBEAM)
        ENDDO

!  set complex constants from the solution vector

        DO K = 1, K_COMPLEX(N)
          K0 = 2*K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          IROW    = K + K_REAL(N)
          IROW1   = IROW + NSTKS_NSTRMS
          IROW_S  = K + K_REAL(N) + K_COMPLEX(N)
          IROW1_S = IROW_S + NSTKS_NSTRMS
          LCON(K1,N) = SCOL2(IROW,    IBEAM)
          LCON(K2,N) = SCOL2(IROW_S,  IBEAM)
          MCON(K1,N) = SCOL2(IROW1,   IBEAM)
          MCON(K2,N) = SCOL2(IROW1_S, IBEAM)
        ENDDO

!  end clause for backsubstitution

      ENDIF

!  Set integration constants for non-active layers
!  ===============================================

!  Transmittance layers above active layer
!  ---------------------------------------

!   -- LCON values are zero (no downwelling radiation)
!   -- MCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!    .... Require solutions at top of active layer

      NAL = ACTIVE_LAYERS(1)
      KO1 = K_REAL(NAL) + 1

      IF ( NAL .GT. 1 ) THEN

        N1 = NAL - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          IR = ( I1 - 1 ) * NSTOKES
          IC = ( I - 1  ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            ICOW = IC + O1

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAL)
             LXR = LCON(K,NAL) * SOLA_XPOS(I1,O1,K,NAL)
             MXR = MCON(K,NAL) * SOLB_XNEG(I1,O1,K,NAL)
             HOM1 = LXR
             HOM2 = MXR * T_DELT_EIGEN(K,NAL)
             SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NAL)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR = LCON(K1,NAL) * SOLA_XPOS(I1,O1,K1,NAL) - &
                       LCON(K2,NAL) * SOLA_XPOS(I1,O1,K2,NAL)
              MXR_CR = MCON(K1,NAL) * SOLB_XNEG(I1,O1,K1,NAL) - &
                       MCON(K2,NAL) * SOLB_XNEG(I1,O1,K2,NAL)
              MXR_CI = MCON(K1,NAL) * SOLB_XNEG(I1,O1,K2,NAL) + &
                       MCON(K2,NAL) * SOLB_XNEG(I1,O1,K1,NAL)
              HOM1CR = LXR_CR
              HOM2CR = MXR_CR*T_DELT_EIGEN(K1,NAL) &
                     - MXR_CI*T_DELT_EIGEN(K2,NAL)
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

            SHOM = SHOM_R + SHOM_CR
            SPAR = WUPPER(I1,O1,NAL)
            MCON(ICOW,N1) = SPAR + SHOM
            LCON(ICOW,N1) = ZERO

          ENDDO
        ENDDO

      ENDIF

!  other layers to top, just propagate

      DO N = NAL - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            LCON(IROW,N) = ZERO
            MCON(IROW,N) = T_DELT_DISORDS(I,N1) * MCON(IROW,N1)
          ENDDO
        ENDDO
      ENDDO

!  Transmittance layers below active layer
!  =======================================

!  Version 2.8. general surface treatment for the problem

!     ** Only do this if active scattering is above (not adjacent to) the surface layer
!   -- LCON values are propagated downwards from bottom of last active layer
!   -- MCON values also propagated downwards, BUT only present if surface condition
!  1.   Require solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      KO1 = K_REAL(NAL) + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  L-constants, always required
!  ----------------------------

        N1 = NAL + 1
        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1

!  real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, K_REAL(NAL)
              LXR = LCON(K,NAL) * SOLA_XPOS(I,O1,K,NAL)
              MXR = MCON(K,NAL) * SOLB_XNEG(I,O1,K,NAL)
              HOM1 = LXR * T_DELT_EIGEN(K,NAL)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO

!  complex homogeneous solutions

            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(NAL)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR = LCON(K1,NAL) * SOLA_XPOS(I,O1,K1,NAL) - &
                       LCON(K2,NAL) * SOLA_XPOS(I,O1,K2,NAL)
              LXR_CI = LCON(K1,NAL) * SOLA_XPOS(I,O1,K2,NAL) + &
                       LCON(K2,NAL) * SOLA_XPOS(I,O1,K1,NAL)
              MXR_CR = MCON(K1,NAL) * SOLB_XNEG(I,O1,K1,NAL) - &
                       MCON(K2,NAL) * SOLB_XNEG(I,O1,K2,NAL)
              HOM1CR = LXR_CR*T_DELT_EIGEN(K1,NAL) &
                      -LXR_CI*T_DELT_EIGEN(K2,NAL)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO

!  real part and add particular solution
!    Note to self. The  MINUS sign appears to be correct !  WHY?

!            SHOM = SHOM_R + SHOM_CR   ! apparently wrong
            SHOM = SHOM_R - SHOM_CR
            SPAR = WLOWER(I,O1,NAL)
            LCON(IROW,N1) = SPAR + SHOM

          ENDDO
        ENDDO

!  Remaining streams by propagation

        DO N = NAL + 2, NLAYERS
          N1 = N - 1
          DO I = 1, NSTREAMS
            IR = ( I - 1 ) * NSTOKES
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              LCON(IROW,N)   = LCON(IROW,N1) * T_DELT_DISORDS(I,N1)
            ENDDO
          ENDDO
        ENDDO

!  M-Constants need to be determined if there is a surface condition. Otherwise zero.
!  ----------------------------------------------------------------------------------

        IF ( DO_INCLUDE_SURFACE ) THEN

!  Next layer

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            IC = ( I - 1  ) * NSTOKES
            DO O1 = 1, NSTOKES
              ICOW = IC + O1

!  real homogeneous solutions

              SHOM_R = ZERO
              DO K = 1, K_REAL(NAL)
                LXR = LCON(K,NAL) * SOLA_XPOS(I1,O1,K,NAL)
                MXR = MCON(K,NAL) * SOLB_XNEG(I1,O1,K,NAL)
                HOM1 = LXR
                HOM2 = MXR * T_DELT_EIGEN(K,NAL)
                SHOM_R = SHOM_R + HOM1 + HOM2
              ENDDO

!  complex homogeneous solutions

              SHOM_CR = ZERO
              DO K = 1, K_COMPLEX(NAL)
                K0 = 2*K - 2
                K1 = KO1 + K0
                K2 = K1  + 1
                LXR_CR = LCON(K1,NAL) * SOLA_XPOS(I1,O1,K1,NAL) - &
                         LCON(K2,NAL) * SOLA_XPOS(I1,O1,K2,NAL)
                MXR_CR = MCON(K1,NAL) * SOLB_XNEG(I1,O1,K1,NAL) - &
                         MCON(K2,NAL) * SOLB_XNEG(I1,O1,K2,NAL)
                MXR_CI = MCON(K1,NAL) * SOLB_XNEG(I1,O1,K2,NAL) + &
                         MCON(K2,NAL) * SOLB_XNEG(I1,O1,K1,NAL)
                HOM1CR = LXR_CR
                HOM2CR = MXR_CR*T_DELT_EIGEN(K1,NAL) &
                       - MXR_CI*T_DELT_EIGEN(K2,NAL)
                SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
              ENDDO

!  real part and add particular solution
!    ---Sets Real integration constants (no complex ones)

              SHOM = SHOM_R + SHOM_CR
              SPAR = WUPPER(I1,O1,NAL)
              MCON(ICOW,NAL+1) = ( SPAR + SHOM ) / T_DELT_DISORDS(I,NAL+1)

!  End

            ENDDO
          ENDDO

!  Remaining layers to ground

          DO N = NAL + 2, NLAYERS
            N1 = N - 1
            DO I = 1, NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                MCON(ICOW,N) = MCON(ICOW,N1) / T_DELT_DISORDS(I,N)
              ENDDO
            ENDDO
          ENDDO

!  No surface - Set M constants to zero

        ELSE

          DO N = NAL + 1, NLAYERS
            N1 = N - 1
            DO I = 1, NSTREAMS
              IC = ( I - 1  ) * NSTOKES
              DO O1 = 1, NSTOKES
                ICOW = IC + O1
                MCON(ICOW,N) = ZERO
              ENDDO
            ENDDO
          ENDDO

!  End surface clause

        ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPTEL_BACKSUB

!

      INTEGER FUNCTION BMAT_ROWMASK ( &
        I, J, NMINJ, NMAXJ, KALL )

!  Integer function for Band-matrix compression
!   Replaces Array BMAT_ROWMASK(I,J) which was too large.
!   Replacement by J. Kujanpaa, FMI, August 2005.

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: I
      INTEGER, INTENT (IN) :: J
      INTEGER, INTENT (IN) :: NMINJ
      INTEGER, INTENT (IN) :: NMAXJ
      INTEGER, INTENT (IN) :: KALL

!  Assign the function (Keep the hard-wired Stop for now)
!    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMINJ) .AND. (I.LE.NMAXJ) ) THEN
         BMAT_ROWMASK = KALL + I - J
      ELSE
         BMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BMAT_ROWMASK:',I,J
         STOP ' -- Fatal error, consult with author'
      ENDIF

!  Finish

      RETURN
      END FUNCTION BMAT_ROWMASK

!

      INTEGER FUNCTION BTELMAT_ROWMASK ( &
        I, J, NMINJ, NMAXJ, KALL )

!  Integer function for Band-matrix compression
!   Replaces Array BTELMAT_ROWMASK(I,J) which was too large.
!   Replacement by J. Kujanpaa, FMI, August 2005.

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: I
      INTEGER, INTENT (IN) :: J
      INTEGER, INTENT (IN) :: NMINJ
      INTEGER, INTENT (IN) :: NMAXJ
      INTEGER, INTENT (IN) :: KALL

!  Assign the function (Keep the hard-wired Stop for now)
!    KALL and NMIN,NMAX are assigned in BVP_MATRIX_INIT

      IF ( (I.GE.NMINJ) .AND. (I.LE.NMAXJ) ) THEN
         BTELMAT_ROWMASK = KALL + I - J
      ELSE
         BTELMAT_ROWMASK = 0
         WRITE(*,*) 'ERROR=out of bounds, BTELMAT_ROWMASK:',I,J
         STOP ' -- Fatal error, consult with author'
      ENDIF

!  Finish

      RETURN
      END FUNCTION BTELMAT_ROWMASK

      END MODULE vlidort_bvproblem_m

