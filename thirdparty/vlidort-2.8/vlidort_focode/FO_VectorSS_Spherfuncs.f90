
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
! ###########################################################

! ###############################################################
! #                                                             #
! #              FIRST-ORDER SCALAR/VECTOR MODEL                #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)            #
! #                                                             #
! #  This Version :   1.5.3                                     #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #   Version 1.1,   13 February  2012, First Code              #
! #   Version 1.2,   01 June      2012, Modularization          #
! #   Version 1.3a,  29 October   2012, Obsgeom Multi-geom.     #
! #   Version 1.3b,  24 January   2013, BRDF/SL Supplements     #
! #   Version 1.4,   31 July      2013, Lattice Multi-geom.     #
! #   Version 1.5,   7  July      2016. Use Fmatrix/Phasfunc    #
! #   Version 1.5,   22 August    2016. Partial-layer output.   #
! #   Version 1.5,   30 April     2017. Shakedown completed.    #
! #   Version 1.5.1, 30 September 2019. Revision Thermal DT.    #
! #   Version 1.5.3, 31 March     2021. Doublet geometry.       #
! #                                                             #
! #   FO Version 1.5   coincides (V)LIDORT Version (2.8)3.8     #
! #   FO Version 1.5.1 coincides (V)LIDORT Version (2.8.1)3.8.1 #
! #   FO Version 1.5.3 coincides (V)LIDORT Version (2.8.3)3.8.3 #
! #                                                             #
! ###############################################################

!    ###########################################################
!    #                                                         #
!    # This is Version 1.5.3 of the FO software library.       #
!    # This library comes with the GNU General Public License, #
!    # Version 3.0. Please read this license carefully.        #
!    #                                                         #
!    #      Copyright (c) 2010-2021.                           #
!    #          Robert Spurr, RT Solutions Inc.                #
!    #                                                         #
!    # This file is part of FO CODE Version 1.5.3.             #
!    #                                                         #
!    # FO CODE is free software: you can redistribute it       #
!    # and/or modify it under the terms of the GNU General     #
!    # Public License as published by the Free Software        #
!    # Foundation, either version 3 of the License, or any     #
!    # later version.                                          #
!    #                                                         #
!    # FO CODE is distributed in the hope that it will be      #
!    # useful, but WITHOUT ANY WARRANTY; without even the      #
!    # implied warranty of MERCHANTABILITY or FITNESS FOR A    #
!    # PARTICULAR PURPOSE.  See the GNU General Public License #
!    # for more details.                                       #
!    #                                                         #
!    # You should have received a copy of the GNU General      #
!    # Public License along with FO CODE Version 1.5.3         #
!    # If not, see <http://www.gnu.org/licenses/>.             #
!    #                                                         #
!    ###########################################################

module FO_VectorSS_spherfuncs_m

!  Generalized spherical functions and rotation-angle sines/cosines

!   Extension to multiple geometries, 29 October 2012       (Version 1.3)
!   Extension to Multiple LATTICE geometries, 19 March 2015 (Version 1.4)

!  Version 1.5, 7/6/16. Optional calculation of spherical functions

!  1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Include option for Doublet-geometry.

public

contains

SUBROUTINE FO_VectorSS_spherfuncs &
        ( MAXMOMENTS, MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, DTR,       & ! Inputs
          NMOMENTS, NGEOMS, NSZAS, NVZAS, NAZMS, NSTOKES, VSIGN,      & ! Inputs
          STARTER, DO_OBSGEOM, DO_DOUBLET, DO_SPHERFUNC, DO_SUNLIGHT, & ! Inputs
          NA_OFFSETS, ND_OFFSETS, THETA, ALPHA, PHI, COSSCAT,         & ! Inputs
          ROTATIONS, GSHELP, GENSPHER )                                ! Outputs

!  1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Include option for Doublet-geometry (Add flag and offsets)
!  1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Rearranged the argument list

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  max dimensions

      INTEGER  , intent(in)    :: MAXMOMENTS, MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS

!  numbers and offsets (1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Add ND_OFFSETS)

      INTEGER  , intent(in)    :: NMOMENTS, NGEOMS, NSZAS, NVZAS, NAZMS, NSTOKES
      INTEGER  , intent(in)    :: NA_OFFSETS ( MAXSZAS, MAXVZAS )
      INTEGER  , intent(in)    :: ND_OFFSETS ( MAXSZAS )

!  Geometry control (1/31/21. Version 2.8.3. (Upgrade from 4/15/20). DO_DOUBLET)

      LOGICAL  , intent(in)    :: DO_OBSGEOM
      LOGICAL  , intent(in)    :: DO_DOUBLET

!  control. [Starter flag may be re-set]. 
!     Spherical function flag, Version 1.5, 7/6/16.

      LOGICAL  , intent(inout) :: STARTER
      LOGICAL  , intent(in)    :: DO_SPHERFUNC
      LOGICAL  , intent(in)    :: DO_SUNLIGHT

!  Geometry

      REAL(fpk), intent(in)    :: DTR, THETA(MAXSZAS), ALPHA(MAXVZAS)
      REAL(fpk), intent(in)    :: PHI(MAXAZMS), COSSCAT(MAXGEOMS), VSIGN

!  Outputs
!  =======

      REAL(fpk), intent(inout) :: ROTATIONS(4,MAXGEOMS)
      REAL(fpk), intent(inout) :: GENSPHER(0:MAXMOMENTS,4,MAXGEOMS)
      REAL(fpk), intent(inout) :: GSHELP(7,MAXMOMENTS)

!  Local
!  =====

      integer   :: L, LNEW, LOLD, ITMP, V, NS, NV, NA
      REAL(fpk) :: Dl, Dl1, Dl21, DL41, DLL1
      REAL(fpk) :: SQL4, SQL41, QROOT6, PERLL4, WFACT
      REAL(fpk) :: CTHETA, STHETA, CALPHA, SALPHA, CPHI
      REAL(fpk) :: THETA_R, ALPHA_R, PHI_R, C1, S1, C2, S2
      REAL(fpk) :: SINSCAT, HELP_SINSCAT, MUP1, MUM1, MUS
      REAL(fpk) :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2
      REAL(fpk) :: P00(2), P02(2), P2P2(2), P2M2(2)

!  Parameter numbers

      REAL(fpk), parameter :: ZERO  = 0.0_fpk
      REAL(fpk), parameter :: ONE   = 1.0_fpk
      REAL(fpk), parameter :: TWO   = 2.0_fpk
      REAL(fpk), parameter :: FOUR  = 4.0_fpk
      REAL(fpk), parameter :: SMALL = 1.0e-12_fpk

!  Section 1: Generate help array, only do this once
!  =================================================

!mick fix 4/3/2015 - moved GENSPHER & ROTATIONS initializations outside if block
   GENSPHER = ZERO ; ROTATIONS = ZERO

!Rob fix 7/6/16. Starter calculation only done if spherical functions are needed.
!   IF ( STARTER ) THEN

   IF ( STARTER .and. DO_SPHERFUNC ) THEN
      !GSHELP = ZERO ; GENSPHER = ZERO ; ROTATIONS = ZERO
      GSHELP = ZERO
      DO L = 2, NMOMENTS
         DL = DBLE(L) ; DL1 = DBLE(L-1) ; DL21 = DBLE(2*L-1)
         GSHELP(1,L) = DL21 / DL
         GSHELP(2,L) = DL1  / DL
         if ( nstokes .gt. 1 ) then
            if ( L .EQ. 2 ) THEN
               SQL41 = ZERO
            ELSE
               SQL4  = SQL41 ; DL41  = DL*DL - FOUR
               SQL41 = DSQRT(DL41)
               GSHELP(3,L)  = DL21 / SQL41
               GSHELP(4,L)  = SQL4 / SQL41
               if ( .not. do_sunlight ) then
                  DLL1 = DL * DL1  ;  PERLL4 = ONE / DL1 / DL41
                  GSHELP(5,L) = DL21 * DLL1 * PERLL4
                  GSHELP(6,L) = DL21 * FOUR * PERLL4
                  GSHELP(7,L) = DL  * ( DL1*DL1 - FOUR) * PERLL4
               endif
            ENDIF
         endif
      ENDDO
      STARTER = .false.
   ENDIF

   QROOT6 = - 0.25_fpk * SQRT(6.0_fpk)

!  Section 2. Observational geometry loop
!  ======================================

   if ( do_ObsGeom ) then

!  Start loop

     do V = 1, ngeoms

!  Section 2a. Rotation angles
!  ---------------------------

       mus = COSSCAT(v)
       THETA_R = DTR * THETA(V) ; CTHETA = COS(THETA_R) ; STHETA = SIN(THETA_R)
       ALPHA_R = DTR * ALPHA(V) ; CALPHA = COS(ALPHA_R) ; SALPHA = SIN(ALPHA_R)
       PHI_R   = DTR * PHI(V)   ; CPHI   = COS(PHI_R)

!  THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING
!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

       HELP_SINSCAT = ( ONE - MUS * MUS )
       IF ( HELP_SINSCAT.LE.ZERO ) THEN
         SINSCAT = SMALL
       ELSE
         SINSCAT = SQRT ( HELP_SINSCAT )
       ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

       IF ( DABS(SINSCAT) .LE. SMALL ) THEN
         CSIG1 = ZERO
         CSIG2 = ZERO
       ELSE
         IF ( STHETA .EQ. ZERO ) THEN
           CSIG1 = -  CPHI
         ELSE
           CSIG1 = ( + VSIGN * CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
!           CSIG1 = ( + CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
         ENDIF
         IF ( SALPHA .EQ. ZERO ) THEN
           CSIG2 = -  CPHI
         ELSE
           CSIG2 = ( - CTHETA - VSIGN * CALPHA * MUS ) / SINSCAT / SALPHA
!           CSIG2 = ( - CTHETA - CALPHA * MUS ) / SINSCAT / SALPHA
         ENDIF
       ENDIF
!      write(*,*)'hello',CSIG1,CSIG2

!  DEBUG, 06 OCTOBER 2010. FOUND DOWNWELLING BUG, NEED SEPARATE GEOMETRY
!      WRITE(7,'(I4,2F6.1,1P6E15.7)')N,VSIGN,PHI/DEG_TO_RAD,
!     &      SALPHA, CALPHA, STHETA, CTHETA, CSIG1,CSIG2

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

       IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
       IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE
       IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
       IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)

       CSIG1_2 = TWO * CSIG1
       CSIG2_2 = TWO * CSIG2
       IF ( ABS(CSIG1-ONE) .LT. SMALL ) THEN
         SSIG1 = ZERO
       ELSE
         SSIG1 = SQRT ( ONE - CSIG1 * CSIG1 )
       ENDIF
       IF ( ABS(CSIG2-ONE).LT. SMALL ) THEN
         SSIG2 = ZERO
       ELSE
         SSIG2 = SQRT ( ONE - CSIG2 * CSIG2 )
       ENDIF

!  RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1/S2. H/VDM, EQS. 94-95.

       C1 = CSIG1_2 * CSIG1 - ONE
       C2 = CSIG2_2 * CSIG2 - ONE

!  PHI IS INPUT IN DEGREES

       IF ( PHI(V) .LE. 180.0_fpk ) THEN
         S1 = CSIG1_2 * SSIG1
         S2 = CSIG2_2 * SSIG2
       ELSE
         S1 = -CSIG1_2 * SSIG1
         S2 = -CSIG2_2 * SSIG2
       ENDIF

!  OUTPUT

       ROTATIONS(1,V) = C1
       ROTATIONS(2,V) = S1
       ROTATIONS(3,V) = C2
       ROTATIONS(4,V) = S2

!  Sections 2b/2c are now optional

       IF ( DO_SPHERFUNC ) THEN

!  Section 2b. Legendre only: Nstokes = 1 (scalar case)
!  ----------------------------------------------------

         IF ( NSTOKES .eq. 1 ) THEN
           GENSPHER(0,1,V) = ONE
           GENSPHER(1,1,V) = MUS
           DO L = 2, NMOMENTS
             GENSPHER(L,1,V) = GSHELP(1,L)*GENSPHER(L-1,1,V)*MUS - GSHELP(2,L)*GENSPHER(L-2,1,V)
           ENDDO
         ENDIF

!  Section 2c. : Nstokes > 1 (vector case). General case
!  -----------------------------------------------------

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

         IF ( NSTOKES .gt. 1 ) THEN

           LNEW = 1
           LOLD = 2

           DO L = 0, NMOMENTS

!  FIRST MOMENT

             IF ( L .EQ. 0 ) THEN
               P00(LOLD)  = ONE  ; P00(LNEW)  = ZERO
               P02(LOLD)  = ZERO ; P02(LNEW)  = ZERO
               P2P2(LOLD) = ZERO ; P2P2(LNEW) = ZERO
               P2M2(LOLD) = ZERO ; P2M2(LNEW) = ZERO
             ELSE IF ( L.eq.1 ) then
               P00(LOLD) = MUS
             ELSE
               P00(LOLD) = GSHELP(1,L) * MUS * P00(LNEW) - GSHELP(2,L)*P00(LOLD)
             END IF
             GENSPHER(L,1,V) = P00(LOLD)

!  Other moments

             IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

               MUM1 = ONE - MUS
               MUP1 = ONE + MUS
               P02(LOLD) = QROOT6 * MUM1 * MUP1
               P02(LNEW) = ZERO
               SQL41     = ZERO
               GENSPHER(L,2,V) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

               if ( .not. do_sunlight ) then
                  P2P2(LOLD)= MUP1 * MUP1 / FOUR
                  P2M2(LOLD)= MUM1 * MUM1 / FOUR
                  GENSPHER(L,3,V) = P2P2(LOLD)
                  GENSPHER(L,4,V) = P2M2(LOLD)
               endif

             ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

               P02(LOLD) = GSHELP(3,L) * MUS * P02(LNEW) - GSHELP(4,L) * P02(LOLD)
               GENSPHER(L,2,V) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

               IF ( .not. do_sunlight ) then
                  WFACT = GSHELP(5,L) * MUS - GSHELP(6,L)
                  P2P2(LOLD) = WFACT * P2P2(LNEW) - GSHELP(7,L) * P2P2(LOLD)
                  WFACT = GSHELP(5,L) * MUS + GSHELP(6,L)
                  P2M2(LOLD) = WFACT * P2M2(LNEW) - GSHELP(7,L) * P2M2(LOLD)
                  GENSPHER(L,3,V) = P2P2(LOLD)
                  GENSPHER(L,4,V) = P2M2(LOLD)
                ENDIF

             END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L; PREVENTS SWAPPING OF ENTIRE ARRAYS.

             ITMP = LNEW
             LNEW = LOLD
             LOLD = ITMP

!  END MOMENTS

           END DO

!  ENd Nstokes > 1 clause

         ENDIF

!  End spherfunc option

       ENDIF

!  END GEOMETRIES

     ENDDO

!  End of Observational geometry section

   ENDIF

!  Section 3. Doublet geometry loop
!  ================================

!  1/31/21. Version 2.8.3. (Upgrade from 4/15/20). This section is competely new for Doublet geometry option

   if ( do_Doublet ) then

    do ns = 1, nszas
     THETA_R = DTR * THETA(ns) ; CTHETA = COS(THETA_R) ; STHETA = SIN(THETA_R)
     do nv = 1, nvzas
       ALPHA_R = DTR * ALPHA(nv) ; CALPHA = COS(ALPHA_R) ; SALPHA = SIN(ALPHA_R)
       PHI_R   = DTR * PHI(nv)   ; CPHI   = COS(PHI_R)
       v = nd_offsets(ns) + nv
       mus = COSSCAT(v)

!  Section 3a. Rotation angles
!  ---------------------------

!  THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING
!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

       HELP_SINSCAT = ( ONE - MUS * MUS )
       IF ( HELP_SINSCAT.LE.ZERO ) THEN
         SINSCAT = SMALL
       ELSE
         SINSCAT = SQRT ( HELP_SINSCAT )
       ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

       IF ( DABS(SINSCAT) .LE. SMALL ) THEN
         CSIG1 = ZERO
         CSIG2 = ZERO
       ELSE
         IF ( STHETA .EQ. ZERO ) THEN
           CSIG1 = -  CPHI
         ELSE
           CSIG1 = ( + VSIGN * CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
!           CSIG1 = ( + CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
         ENDIF
         IF ( SALPHA .EQ. ZERO ) THEN
           CSIG2 = -  CPHI
         ELSE
           CSIG2 = ( - CTHETA - VSIGN * CALPHA * MUS ) / SINSCAT / SALPHA
!           CSIG2 = ( - CTHETA - CALPHA * MUS ) / SINSCAT / SALPHA
         ENDIF
       ENDIF
!       write(*,*)'hello',CSIG1,CSIG2

!  DEBUG, 06 OCTOBER 2010. FOUND DOWNWELLING BUG, NEED SEPARATE GEOMETRY
!      WRITE(7,'(I4,2F6.1,1P6E15.7)')N,VSIGN,PHI/DEG_TO_RAD,
!     &      SALPHA, CALPHA, STHETA, CTHETA, CSIG1,CSIG2

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

       IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
       IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE
       IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
       IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)

       CSIG1_2 = TWO * CSIG1
       CSIG2_2 = TWO * CSIG2
       IF ( ABS(CSIG1-ONE) .LT. SMALL ) THEN
         SSIG1 = ZERO
       ELSE
         SSIG1 = SQRT ( ONE - CSIG1 * CSIG1 )
       ENDIF
       IF ( ABS(CSIG2-ONE).LT. SMALL ) THEN
         SSIG2 = ZERO
       ELSE
         SSIG2 = SQRT ( ONE - CSIG2 * CSIG2 )
       ENDIF

!  RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1/S2. H/VDM, EQS. 94-95.

       C1 = CSIG1_2 * CSIG1 - ONE
       C2 = CSIG2_2 * CSIG2 - ONE

!  PHI IS INPUT IN DEGREES

       IF ( PHI(nv) .LE. 180.0_fpk ) THEN
         S1 = CSIG1_2 * SSIG1
         S2 = CSIG2_2 * SSIG2
       ELSE
         S1 = -CSIG1_2 * SSIG1
         S2 = -CSIG2_2 * SSIG2
       ENDIF

!  OUTPUT

       ROTATIONS(1,V) = C1
       ROTATIONS(2,V) = S1
       ROTATIONS(3,V) = C2
       ROTATIONS(4,V) = S2

!  Sections 3b/3c are now optional

       IF ( DO_SPHERFUNC ) THEN

!  Section 3b. Legendre only: Nstokes = 1 (scalar case)
!  ----------------------------------------------------

         IF ( NSTOKES .eq. 1 ) THEN
           GENSPHER(0,1,V) = ONE
           GENSPHER(1,1,V) = MUS
           DO L = 2, NMOMENTS
             GENSPHER(L,1,V) = GSHELP(1,L)*GENSPHER(L-1,1,V)*MUS - GSHELP(2,L)*GENSPHER(L-2,1,V)
           ENDDO
         ENDIF

!  Section 3c. : Nstokes > 1 (vector case). General case
!  -----------------------------------------------------

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

         IF ( NSTOKES .GT. 1 ) THEN

           LNEW = 1
           LOLD = 2
   
           DO L = 0, NMOMENTS

!  FIRST MOMENT

             IF ( L .EQ. 0 ) THEN
               P00(LOLD)  = ONE  ; P00(LNEW)  = ZERO
               P02(LOLD)  = ZERO ; P02(LNEW)  = ZERO
               P2P2(LOLD) = ZERO ; P2P2(LNEW) = ZERO
               P2M2(LOLD) = ZERO ; P2M2(LNEW) = ZERO
             ELSE IF ( L.eq.1 ) then
               P00(LOLD) = MUS
             ELSE
               P00(LOLD) = GSHELP(1,L) * MUS * P00(LNEW) - GSHELP(2,L)*P00(LOLD)
             END IF
             GENSPHER(L,1,V) = P00(LOLD)

!  Other moments

             IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

               MUM1 = ONE - MUS
               MUP1 = ONE + MUS
               P02(LOLD) = QROOT6 * MUM1 * MUP1
               P02(LNEW) = ZERO
               SQL41     = ZERO
               GENSPHER(L,2,V) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

               if ( .not. do_sunlight ) then
                  P2P2(LOLD)= MUP1 * MUP1 / FOUR
                  P2M2(LOLD)= MUM1 * MUM1 / FOUR
                  GENSPHER(L,3,V) = P2P2(LOLD)     
                  GENSPHER(L,4,V) = P2M2(LOLD)
               endif

             ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

               P02(LOLD) = GSHELP(3,L) * MUS * P02(LNEW) - GSHELP(4,L) * P02(LOLD)
               GENSPHER(L,2,V) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

               IF ( .not. do_sunlight ) then
                  WFACT = GSHELP(5,L) * MUS - GSHELP(6,L)
                  P2P2(LOLD) = WFACT * P2P2(LNEW) - GSHELP(7,L) * P2P2(LOLD)
                  WFACT = GSHELP(5,L) * MUS + GSHELP(6,L)
                  P2M2(LOLD) = WFACT * P2M2(LNEW) - GSHELP(7,L) * P2M2(LOLD)
                  GENSPHER(L,3,V) = P2P2(LOLD)
                  GENSPHER(L,4,V) = P2M2(LOLD)
               ENDIF

             END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L; PREVENTS SWAPPING OF ENTIRE ARRAYS.

             ITMP = LNEW
             LNEW = LOLD
             LOLD = ITMP

!  END MOMENTS

           END DO

!  ENd Nstokes > 1 clause

         ENDIF

!  End spherfunc option

       ENDIF

!  END GEOMETRIES

     ENDDO
    ENDDO

!  End DOUBLET Geometry option

   ENDIF

!  Section 4. Lattice geometry loop
!  ================================

!  1/31/21. Version 2.8.3. (Upgrade from 4/15/20). Add do_doublet to if-condition for lattce.

   if ( .not. do_obsgeom .and. .not. do_Doublet ) then

    do ns = 1, nszas
     THETA_R = DTR * THETA(ns) ; CTHETA = COS(THETA_R) ; STHETA = SIN(THETA_R)
     do nv = 1, nvzas
      ALPHA_R = DTR * ALPHA(nv) ; CALPHA = COS(ALPHA_R) ; SALPHA = SIN(ALPHA_R)
      do na = 1, nazms
       PHI_R   = DTR * PHI(na)   ; CPHI   = COS(PHI_R)
       v = na_offsets(ns,nv) + na
       mus = COSSCAT(v)

!  Section 4a. Rotation angles
!  ---------------------------

!  THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING
!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

       HELP_SINSCAT = ( ONE - MUS * MUS )
       IF ( HELP_SINSCAT.LE.ZERO ) THEN
         SINSCAT = SMALL
       ELSE
         SINSCAT = SQRT ( HELP_SINSCAT )
       ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

       IF ( DABS(SINSCAT) .LE. SMALL ) THEN
         CSIG1 = ZERO
         CSIG2 = ZERO
       ELSE
         IF ( STHETA .EQ. ZERO ) THEN
           CSIG1 = -  CPHI
         ELSE
           CSIG1 = ( + VSIGN * CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
!           CSIG1 = ( + CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
         ENDIF
         IF ( SALPHA .EQ. ZERO ) THEN
           CSIG2 = -  CPHI
         ELSE
           CSIG2 = ( - CTHETA - VSIGN * CALPHA * MUS ) / SINSCAT / SALPHA
!           CSIG2 = ( - CTHETA - CALPHA * MUS ) / SINSCAT / SALPHA
         ENDIF
       ENDIF
!       write(*,*)'hello',CSIG1,CSIG2

!  DEBUG, 06 OCTOBER 2010. FOUND DOWNWELLING BUG, NEED SEPARATE GEOMETRY
!      WRITE(7,'(I4,2F6.1,1P6E15.7)')N,VSIGN,PHI/DEG_TO_RAD,
!     &      SALPHA, CALPHA, STHETA, CTHETA, CSIG1,CSIG2

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

       IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
       IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE
       IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
       IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)

       CSIG1_2 = TWO * CSIG1
       CSIG2_2 = TWO * CSIG2
       IF ( ABS(CSIG1-ONE) .LT. SMALL ) THEN
         SSIG1 = ZERO
       ELSE
         SSIG1 = SQRT ( ONE - CSIG1 * CSIG1 )
       ENDIF
       IF ( ABS(CSIG2-ONE).LT. SMALL ) THEN
         SSIG2 = ZERO
       ELSE
         SSIG2 = SQRT ( ONE - CSIG2 * CSIG2 )
       ENDIF

!  RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1/S2. H/VDM, EQS. 94-95.

       C1 = CSIG1_2 * CSIG1 - ONE
       C2 = CSIG2_2 * CSIG2 - ONE

!  PHI IS INPUT IN DEGREES

       IF ( PHI(na) .LE. 180.0_fpk ) THEN
         S1 = CSIG1_2 * SSIG1
         S2 = CSIG2_2 * SSIG2
       ELSE
         S1 = -CSIG1_2 * SSIG1
         S2 = -CSIG2_2 * SSIG2
       ENDIF

!  OUTPUT

       ROTATIONS(1,V) = C1
       ROTATIONS(2,V) = S1
       ROTATIONS(3,V) = C2
       ROTATIONS(4,V) = S2

!  Sections 3b/3c are now optional

       IF ( DO_SPHERFUNC ) THEN

!  Section 4b. Legendre only: Nstokes = 1 (scalar case)
!  ----------------------------------------------------

         IF ( NSTOKES .eq. 1 ) THEN
           GENSPHER(0,1,V) = ONE
           GENSPHER(1,1,V) = MUS
           DO L = 2, NMOMENTS
             GENSPHER(L,1,V) = GSHELP(1,L)*GENSPHER(L-1,1,V)*MUS - GSHELP(2,L)*GENSPHER(L-2,1,V)
           ENDDO
         ENDIF

!  Section 4c. : Nstokes > 1 (vector case). General case
!  -----------------------------------------------------

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

         IF ( NSTOKES .GT. 1 ) THEN

           LNEW = 1
           LOLD = 2
   
           DO L = 0, NMOMENTS

!  FIRST MOMENT

             IF ( L .EQ. 0 ) THEN
               P00(LOLD)  = ONE  ; P00(LNEW)  = ZERO
               P02(LOLD)  = ZERO ; P02(LNEW)  = ZERO
               P2P2(LOLD) = ZERO ; P2P2(LNEW) = ZERO
               P2M2(LOLD) = ZERO ; P2M2(LNEW) = ZERO
             ELSE IF ( L.eq.1 ) then
               P00(LOLD) = MUS
             ELSE
               P00(LOLD) = GSHELP(1,L) * MUS * P00(LNEW) - GSHELP(2,L)*P00(LOLD)
             END IF
             GENSPHER(L,1,V) = P00(LOLD)

!  Other moments

             IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

               MUM1 = ONE - MUS
               MUP1 = ONE + MUS
               P02(LOLD) = QROOT6 * MUM1 * MUP1
               P02(LNEW) = ZERO
               SQL41     = ZERO
               GENSPHER(L,2,V) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

               if ( .not. do_sunlight ) then
                  P2P2(LOLD)= MUP1 * MUP1 / FOUR
                  P2M2(LOLD)= MUM1 * MUM1 / FOUR
                  GENSPHER(L,3,V) = P2P2(LOLD)     
                  GENSPHER(L,4,V) = P2M2(LOLD)
               endif

             ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

               P02(LOLD) = GSHELP(3,L) * MUS * P02(LNEW) - GSHELP(4,L) * P02(LOLD)
               GENSPHER(L,2,V) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

               IF ( .not. do_sunlight ) then
                  WFACT = GSHELP(5,L) * MUS - GSHELP(6,L)
                  P2P2(LOLD) = WFACT * P2P2(LNEW) - GSHELP(7,L) * P2P2(LOLD)
                  WFACT = GSHELP(5,L) * MUS + GSHELP(6,L)
                  P2M2(LOLD) = WFACT * P2M2(LNEW) - GSHELP(7,L) * P2M2(LOLD)
                  GENSPHER(L,3,V) = P2P2(LOLD)
                  GENSPHER(L,4,V) = P2M2(LOLD)
               ENDIF

             END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L; PREVENTS SWAPPING OF ENTIRE ARRAYS.

             ITMP = LNEW
             LNEW = LOLD
             LOLD = ITMP

!  END MOMENTS

           END DO

!  ENd Nstokes > 1 clause

         ENDIF

!  End spherfunc option

       ENDIF

!  END GEOMETRIES

      ENDDO
     ENDDO
    ENDDO

!  End LATTICE Geometry option

   ENDIF

!  Finish

   RETURN
END SUBROUTINE FO_VectorSS_spherfuncs

!  Finish

end module FO_VectorSS_spherfuncs_m


