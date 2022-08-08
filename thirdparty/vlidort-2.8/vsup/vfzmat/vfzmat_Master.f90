
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

module vfzmat_Master_m

!mick mod 9/19/2017 - use new numerical subroutines

   Use vfzmat_Numerical_m, only : GETQUAD2, BSpline, Seval
   Use vfzmat_Rotation_m
   Use vfzmat_PhasMat_m,   only : vfzmat_PhasMat2
   Use vfzmat_DevelopCoeffs_m
   Use vfzmat_ExpandCoeffs_m

public

contains

subroutine vfzmat_Master &
   ( max_moments, max_geoms, max_szas, max_vzas, max_azms, maxlayers,    & ! input  Dimensions (VLIDORT)
     max_InAngles, N_InAngles, InAngles, InFmatrices, Exist_InFmatrices, & ! input  Fmatrices
     do_upwelling, do_dnwelling, do_ObsGeoms, Sunlight,                  & ! input  Flags
     ncoeffs, nlayers, nstokes, n_geoms, n_szas, n_vzas, n_azms,         & ! input  Numbers
     offsets, dtr, szas, vzas, azms, obsgeoms,                           & ! Input  Geometries
     OutFmatrices_up, OutFmatrices_dn, Zmatrices_up, Zmatrices_dn, FMatCoeffs )

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  Purpose
!  -------

!  Stand-alone routine to develop Z-matrices and F-Matrix expansion coefficients,
!    given only a set of F-Matrix inputs on a regular scattering-angle grid

!  Allows for various VLIDORT-based geometry options, upwelling and/or downwelling.

!  Stage 1. (a) Interpolate F-matrices to values implied by geometrical input.
!           (b) Transformation from scattering plane (Fmat) to meridional plane 
!               (Zmat) follows the rotations given in the VLIDORT code.

!  Stage 2. (a) Interpolate F-matrix input to quadrature grid for integration
!           (b) develop Coefficients from integrations using spherical-functions

!  Fmatrix input convention (J is the angle grid on input, N the layer index)
!  ------------------------

!       InFmatrices(N,J,1) = F11(J)
!       InFmatrices(N,J,2) = F22(J)
!       InFmatrices(N,J,3) = F33(J)
!       InFmatrices(N,J,4) = F44(J)
!       InFmatrices(N,J,5) = F12(J)
!       InFmatrices(N,J,6) = F34(J)

!   For Mie scattering, F11 = F22, F33 = F44

!  Convention for using FMatCoeffs output to get VLIDORT Greekmat input
!  --------------------------------------------------------------------

!       GREEKMAT(L,N,1)  --> use + FMatCoeffs(L,N,1)
!       GREEKMAT(L,N,2)  --> use - FMatCoeffs(L,N,5)
!       GREEKMAT(L,N,5)  --> use - FMatCoeffs(L,N,5)
!       GREEKMAT(L,N,6)  --> use + FMatCoeffs(L,N,2)
!       GREEKMAT(L,N,11) --> use + FMatCoeffs(L,N,3)
!       GREEKMAT(L,N,12) --> use - FMatCoeffs(L,N,6)
!       GREEKMAT(L,N,15) --> use + FMatCoeffs(L,N,6)
!       GREEKMAT(L,N,16) --> use + FMatCoeffs(L,N,4)

!    all other entries zero.

!  PATCH Upgrade, 08 November 2019
!  ===============================

!   BSPLINE/Seval routines must be done together inside m-loop
!   formerly, BSPLINE output was only for the m = 6 value.

   implicit none

!  precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER, INTENT(IN)   :: Max_InAngles, max_moments, maxlayers
   INTEGER, INTENT(IN)   :: max_geoms, max_szas, max_vzas, max_azms

!  Directional flags

   LOGICAL, INTENT(IN)   :: do_upwelling, do_dnwelling

!  Flags for use of observational geometry, Sunlight
 
   LOGICAL, INTENT(IN)   :: do_obsgeoms, Sunlight

!  Input F-matrix stuff ( angles and 6 scattering matrix entries )
!  Existence flags introduced 9/19/16

   INTEGER  , INTENT(IN) :: N_InAngles
   REAL(dpk), INTENT(IN) :: InAngles    ( Max_InAngles )
   REAL(dpk), INTENT(IN) :: InFmatrices ( Maxlayers, Max_InAngles, 6 )
   LOGICAL  , INTENT(IN) :: Exist_InFmatrices ( Maxlayers )

!  numbers (general)

   INTEGER, INTENT(IN)   :: nlayers, nstokes, ncoeffs

!  Geometry numbers

   INTEGER, INTENT(IN)   :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER, INTENT(IN)   :: offsets(max_szas,max_vzas)

!  Angles. Convention as for  VLIDORT

   REAL(dpk), INTENT(IN) :: dtr
   REAL(dpk), INTENT(IN) :: szas(max_szas)
   REAL(dpk), INTENT(IN) :: vzas(max_vzas)
   REAL(dpk), INTENT(IN) :: azms(max_azms)
   REAL(dpk), INTENT(IN) :: obsgeoms(max_geoms,3)

!  Output
!  ------

!  Output Fmatrices (Interpolated)

   REAL(dpk), INTENT(OUT) :: OutFmatrices_up   ( Max_Geoms, MaxLayers, 6 )
   REAL(dpk), INTENT(OUT) :: OutFmatrices_dn   ( Max_Geoms, MaxLayers, 6 )

!  Zmatrices

   REAL(dpk), INTENT(OUT) :: Zmatrices_up(max_geoms,maxlayers,4,4)
   REAL(dpk), INTENT(OUT) :: Zmatrices_dn(max_geoms,maxlayers,4,4)

!  Fmatrix coefficients

   REAL(dpk), INTENT(OUT) :: FMatCoeffs(maxlayers,0:max_moments,6)

!  local variables
!  ---------------

!  Parameters

   REAL(dpk), PARAMETER :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   REAL(dpk), PARAMETER :: d_half  = 0.5_dpk, d_two  = 2.0_dpk
   INTEGER  , PARAMETER :: Max_OutAngles = 5000                   ! quadrature

!  Quadrature

   INTEGER   :: N_OutAngles
   REAL(dpk) :: OutAngles  ( Max_OutAngles )
   REAL(dpk) :: OutCosines ( Max_OutAngles )
   REAL(dpk) :: OutWeights ( Max_OutAngles )

!  Local Output Quad Fmatrices (Interpolated)

   REAL(dpk) :: OutFmatrices_Quad ( Max_OutAngles, 6 )

!  Coefficients. Initialized locally inside "Develop" subroutine

   REAL(dpk) :: expcoeffs(0:max_OutAngles,6)

!  rotational angles

   REAL(dpk) :: C1(max_geoms), S1(max_geoms)
   REAL(dpk) :: C2(max_geoms), S2(max_geoms)

!  Scattering angle cosines

   REAL(dpk) :: COSSCAT_up(max_geoms)
   REAL(dpk) :: COSSCAT_dn(max_geoms)

!  Directional sign

   REAL(dpk) :: vsign

!  Help variables

   INTEGER   :: i, ia, ib, k, k1, k2, l, m, n, um, v
   REAL(dpk) :: ctheta, stheta, calpha, salpha, cphi
   REAL(dpk) :: InCosines(Max_InAngles)
   REAL(dpk) :: Local_InFmatrices(Max_InAngles,6)
   REAL(dpk) :: Check_InFmatrices(Max_InAngles,6)
   !REAL(dpk) :: x1,x2,w1,w2,x,y,y2(Max_InAngles,6),yp1,ypn
   REAL(dpk) :: x1,x2,w1,w2,x,y
   REAL(dpk) :: bbas( Max_InAngles ),cbas( Max_InAngles ),dbas( Max_InAngles )

!  Debug

   LOGICAL, PARAMETER :: do_debug_input  = .false.
   LOGICAL, PARAMETER :: check_expansion = .false.

!  VFZMAT input debug

      IF (do_debug_input) CALL vfzmat_debug_input()

!  Initialize
!  ----------

!  Zero output

   OutFmatrices_up = d_zero
   OutFmatrices_dn = d_zero
   Zmatrices_up    = d_zero
   Zmatrices_dn    = d_zero
   FMatCoeffs      = d_zero

!  Quadrature - set number of angles

   N_OutAngles = Max_OutAngles

!  STAGE 1. Develop Z-matrices
!  ===========================

!  1. Get the scattering angle cosines
!  -----------------------------------

   IF ( .not. Do_Obsgeoms ) THEN
     DO IB = 1, n_szas
       ctheta = cos ( szas(ib) * dtr )
       stheta = sin ( szas(ib) * dtr )
       DO UM = 1, n_vzas
         calpha = cos ( vzas(um) * dtr )
         salpha = sin ( vzas(um) * dtr )
         DO IA = 1, n_azms
           cphi = cos ( azms(ia) * dtr )
           V = OFFSETS(IB,UM) + IA
           COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
           COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
         ENDDO
       ENDDO
     ENDDO
   ELSE
     DO V = 1, n_geoms
       CTHETA = cos(obsgeoms(v,1)*dtr)
       STHETA = sin(obsgeoms(v,1)*dtr)
       CALPHA = cos(obsgeoms(v,2)*dtr)
       SALPHA = sin(obsgeoms(v,2)*dtr)
       CPHI   = cos(obsgeoms(v,3)*dtr)
       COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
       COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
     ENDDO
   ENDIF

!  Cosines of input angles.
!    -- Reverse directions for the Splining (Monotonically increasing)

   do k = 1, N_InAngles
     k1 = N_InAngles + 1 - k 
     InCosines(k1) = cos ( InAngles(k) * dtr )
   enddo

!  2. Spline-Interpolate F-matrices to VLIDORT Geometrical grid
!  ------------------------------------------------------------

!  Start layer loop - Only perform if the F-matrices exist in a given layer (9/19/16)

   DO N = 1, nlayers
     if ( Exist_InFmatrices(n) ) then

!  Local F-matrix, reverse order

       do k = 1, N_InAngles
         k1 = N_InAngles + 1 - k 
         Local_InFmatrices(k1,1:6) = InFmatrices(n,k,1:6)
       enddo

!  Start coefficient loop

       do m = 1, 6

!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

         !yp1 = d_zero
         !ypn = ( Local_InFmatrices(N_InAngles,m) - Local_InFmatrices(N_InAngles-1,m) ) / &
         !        ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
         !Call vfzmat_SPLINE(Max_InAngles,InCosines,Local_InFmatrices(:,m),N_InAngles,yp1,ypn,y2(:,m))
         Call BSpline (Max_InAngles,N_InAngles,InCosines,Local_InFmatrices(:,m),bbas,cbas,dbas)

!  Interpolate upwelling

         if ( do_upwelling ) then
           do v = 1, n_geoms
             x = cosscat_up(v)
             !Call vfzmat_SPLINT(Max_InAngles,InCosines,Local_InFmatrices(:,m),y2(:,m),N_InAngles,x,y)
             Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InFmatrices(:,m),bbas,cbas,dbas,y)
             OutFmatrices_up(v,n,m) = y
           enddo
         endif

!  Interpolate downwelling

         if ( do_dnwelling ) then
           do v = 1, n_geoms
             x = cosscat_dn(v)
             !Call vfzmat_SPLINT(Max_InAngles,InCosines,Local_InFmatrices(:,m),y2(:,m),N_InAngles,x,y)
             Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InFmatrices(:,m),bbas,cbas,dbas,y)
             OutFmatrices_dn(v,n,m) = y
           enddo
         endif

!  End coefficient type loop

       enddo

!  End layer loop

     endif
   enddo

!  3. Get the Z-matrices
!  ---------------------

!  upwelling. [ C1/S1/C2/S2 are local, will be overwritten ]

   if ( do_upwelling ) then
      vsign = - d_one
      Call vfzmat_Rotation &
       ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
         do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
         offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
         C1, S1, C2, S2 )
      Call vfzmat_PhasMat2 &
       ( max_geoms, maxlayers, nlayers, nstokes, n_geoms, Sunlight, & ! inputs
         C1, S1, C2, S2, Exist_InFmatrices, OutFmatrices_up,        & ! Inputs
         Zmatrices_up )
   endif

!  Downwelling

   if ( do_dnwelling ) then
      vsign = + d_one
      Call vfzmat_Rotation &
       ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
         do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
         offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
         C1, S1, C2, S2 )
      Call vfzmat_PhasMat2 &
       ( max_geoms, maxlayers, nlayers, nstokes, n_geoms, Sunlight, & ! inputs
         C1, S1, C2, S2, Exist_InFmatrices, OutFmatrices_dn,        & ! Inputs
         Zmatrices_dn )
   endif

!  STAGE 2. Develop Coefficients
!  =============================

!  1. Quadrature
!  -------------

   Call GETQUAD2 ( -d_one, d_one, N_OutAngles, OutCosines, OutWeights )
   do k1 = 1, N_OutAngles / 2
      k2 = N_OutAngles + 1 - k1 
      x1 = OutCosines(k1) ; x2 = OutCosines(k2)
      w1 = OutWeights(k1) ; w2 = OutWeights(k2)
      OutCosines(k1) = x2 ; OutCosines(k2) = x1
      OutWeights(k1) = w2 ; OutWeights(k2) = w1
   enddo
   do k = 1, N_OutAngles
      OutAngles(k) = acos ( OutCosines(k) ) / dtr
   enddo

!  2. Spline-Interpolate F-matrices to Quadrature grid for coefficients
!  --------------------------------------------------------------------

!  Start layer loop - Only perform if the F-matrices exist in a given layer (9/19/16)

   DO N = 1, nlayers
     if ( Exist_InFmatrices(n) ) then

!  Local F-matrix, reverse order

       do k = 1, N_InAngles
         k1 = N_InAngles + 1 - k 
         Local_InFmatrices(k1,1:6) = InFmatrices(n,k,1:6)
       enddo

!  Start coefficient loop

       do m = 1, 6

!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

         !yp1 = d_zero
         !ypn = ( Local_InFmatrices(N_InAngles,m) - Local_InFmatrices(N_InAngles-1,m) ) / &
         !        ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
         !Call vfzmat_SPLINE(Max_InAngles,InCosines,Local_InFmatrices(:,m),N_InAngles,yp1,ypn,y2(:,m))
         Call BSpline (Max_InAngles,N_InAngles,InCosines,Local_InFmatrices(:,m),bbas,cbas,dbas)

!  interpolate

         do k = 1, N_OutAngles
            k1 = N_OutAngles + 1 - k 
            x = OutCosines(k1)
            !Call vfzmat_SPLINT(Max_InAngles,InCosines,Local_InFmatrices(:,m),y2(:,m),N_InAngles,x,y)
            Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InFmatrices(:,m),bbas,cbas,dbas,y)
            OutFmatrices_Quad(k1,m) = y
         enddo

!  End coefficient loop

       enddo

!  Get the coefficients from the DevelopCoeffs routine.
!     Output = Expcoeffs,  which is local array

       Call vfzmat_DevelopCoeffs &
         ( max_OutAngles, ncoeffs, N_OutAngles, &
          OutCosines, OutWeights, OutFmatrices_Quad, Expcoeffs )

!  Set Output
 
       do L = 0, ncoeffs
         FMatCoeffs(N,L,1:6) = Expcoeffs(L,1:6)
       enddo

!  Check workings by calling Expand routine, just for one layer
!    - this should give results very close to the original Fmatrix input.  

       if ( check_expansion ) then
         Call vfzmat_ExpandCoeffs &
           ( max_InAngles, max_OutAngles, n_InAngles, ncoeffs, nstokes, InCosines, expcoeffs, Check_InFmatrices )
         if ( n.eq.1 ) then
           do k = 1, N_InAngles
             k1 = N_InAngles + 1 - k 
             write(771,55)k,k1,InAngles(k),InCosines(k1), InFmatrices(n,k,1:6)
             write(772,55)k,k1,InAngles(k),InCosines(k1), Check_InFmatrices(k,1:6)
           enddo
55         format(2i5,2f10.4,1p6e15.7)
         endif
       endif

!  End layer loop

     endif
   enddo

!  Done

   return

   contains

!  VFZMAT subroutine for debugging input

subroutine vfzmat_debug_input()

   write(*,*)
   write(*,*) 'vfzmat_Master inputs:'

   write(*,*)
   write(*,*) 'max_moments  = ',max_moments
   write(*,*) 'max_geoms    = ',max_geoms
   write(*,*) 'max_szas     = ',max_szas
   write(*,*) 'max_vzas     = ',max_vzas
   write(*,*) 'max_azms     = ',max_azms
   write(*,*) 'maxlayers    = ',maxlayers
   write(*,*) 'max_InAngles = ',max_InAngles
   write(*,*) 'N_InAngles   = ',N_InAngles

   write(*,*)
   do i=1,N_InAngles
     write(*,*) 'i = ',i,' InAngles(i) = ',InAngles(i)
   enddo

   write(*,*)
   do n=1,nlayers
     do i=1,N_InAngles
       do k=1,6
         write(*,*) 'n = ',n,' i = ',i,' k = ',k,' InFmatrices(n,i,k) = ',InFmatrices(n,i,k)
       enddo
     enddo
   enddo

   write(*,*)
   do n=1,nlayers
     write(*,*) 'n = ',n,'Exist_InFmatrices(n) = ',Exist_InFmatrices(n)
   enddo

   write(*,*)
   write(*,*) 'do_upwelling = ',do_upwelling
   write(*,*) 'do_dnwelling = ',do_dnwelling
   write(*,*) 'do_ObsGeoms  = ',do_ObsGeoms
   write(*,*) 'Sunlight     = ',Sunlight

   write(*,*)
   write(*,*) 'ncoeffs      = ',ncoeffs
   write(*,*) 'nlayers      = ',nlayers
   write(*,*) 'nstokes      = ',nstokes
   write(*,*) 'n_szas       = ',n_szas
   write(*,*) 'n_vzas       = ',n_vzas
   write(*,*) 'n_azms       = ',n_azms
   if ( do_ObsGeoms ) write(*,*) 'n_geoms      = ',n_geoms

   write(*,*)
   do ib=1,n_szas
     do i=1,n_vzas
       write(*,*) 'ib = ',ib,' um = ',um,'offsets(ib,um) = ',offsets(ib,um)
     enddo
   enddo

   write(*,*)
   write(*,*) 'dtr = ',dtr

   write(*,*)
   do ib=1,n_szas
     write(*,*) 'ib = ',ib,' szas(ib) = ',szas(ib)
   enddo

   write(*,*)
   do um=1,n_vzas
     write(*,*) 'um = ',um,' vzas(um) = ',vzas(um)
   enddo

   write(*,*)
   do ia=1,n_azms
     write(*,*) 'ia = ',ia,' azms(ia) = ',azms(ia)
   enddo

   write(*,*)
   if ( do_ObsGeoms ) then
     do v=1,n_geoms
       write(*,*) 'v = ',v,' obsgeoms(v,1:3) = ',obsgeoms(v,1:3)
     enddo
   endif

   !stop

end subroutine vfzmat_debug_input

end subroutine vfzmat_Master

!  done module

end module vfzmat_Master_m

