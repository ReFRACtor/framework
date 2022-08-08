
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
! #            vfzmat_PhasMat1                                  #
! #            vfzmat_PhasMat2                                  #
! #                                                             #
! ###############################################################

module vfzmat_PhasMat_m

public

contains

subroutine vfzmat_PhasMat1 &
   ( max_geoms, nstokes, n_geoms, Sunlight, & ! inputs
     C1, S1, C2, S2, Fmatrices,             & ! Inputs
     Zmatrices )

!  Stand-alone routine to develop Z-matrices from F-Matrix inputs for
!  various goemetry combinations, upwelling or downwelling. SINGLE ENTRY.

!  Programmed 19 September 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  Transformation from scattering plan (Fmat) to Meriodional Plane (Zmat)
!   follows the rotations given in the VLIDORT code.

!  Remarks
!  -------

!  1. Angles and numbers should be set up beforehand as for VLIDORT.

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input
!  -----

!  dimensions

   INTEGER, INTENT (IN) :: max_geoms

!  numbers

   INTEGER, INTENT (IN) :: nstokes, n_geoms

!  Flag for Sunlight
 
   LOGICAL, intent(in)   :: Sunlight

!  Rotational variables

   REAL(KIND=dpk), INTENT (IN)   :: C1(max_geoms), S1(max_geoms)
   REAL(KIND=dpk), INTENT (IN)   :: C2(max_geoms), S2(max_geoms)

!  Fmatrices

   REAL(KIND=dpk), INTENT (IN) :: Fmatrices(max_geoms,6)

!  output
!  ------

!  Zmatrices output, initialized globally beforehand

   REAL(KIND=dpk), INTENT (InOut) :: Zmatrices(max_geoms,4,4)

!  local variables
!  ---------------

   INTEGER            :: v, index_11, index_12, index_22, index_33, index_34, index_44 
   REAL    (KIND=dpk) :: HELP2C1, HELP2S1, HELP3C1, HELP3S1

!  New indexing consistent with Tmatrix output

  index_11 = 1
  index_22 = 2
  index_33 = 3
  index_44 = 4
  index_12 = 5
  index_34 = 6

!  Nstokes = 1, just copy F11 to Z(1,1) element and return (no rotation)

   if ( nstokes .eq. 1 ) then
      do v = 1, n_geoms
         Zmatrices(v,1,1) = Fmatrices(v,index_11) 
      enddo
      return
   endif

!  Z_matrices, general case
!  ------------------------
 
!  For sunlight, only need the first column of Z matrix

   IF ( SUNLIGHT ) THEN
      DO V = 1, N_GEOMS
         Zmatrices(V,1,1) =   Fmatrices(V,Index_11)
         Zmatrices(V,2,1) = - Fmatrices(V,Index_12) * C2(V)
         Zmatrices(V,3,1) =   Fmatrices(V,Index_12) * S2(V)
      ENDDO
   ENDIF

!  For general case. Not tested as of 15 April 2005.

   IF ( .NOT. SUNLIGHT ) THEN
      DO V = 1, N_GEOMS
         Zmatrices(V,1,1) =   Fmatrices(V,Index_11)
         Zmatrices(V,2,1) = - Fmatrices(V,Index_12) * C2(V)
         Zmatrices(V,3,1) =   Fmatrices(V,Index_12) * S2(V)
!  This code untested
         HELP2C1 = Fmatrices(V,Index_22) * C1(V)
         HELP2S1 = Fmatrices(V,Index_22) * S1(V)
         HELP3C1 = Fmatrices(V,Index_33) * C1(V)
         HELP3S1 = Fmatrices(V,Index_33) * S1(V)
         Zmatrices(V,1,2) =   Fmatrices(V,Index_12) * C1(V)
         Zmatrices(V,1,3) = - Fmatrices(V,Index_12) * S1(V)
         Zmatrices(V,2,2) =   C2(V) * HELP2C1 - S2(V) * HELP3S1
         Zmatrices(V,2,3) = - C2(V) * HELP2S1 - S2(V) * HELP3C1
         Zmatrices(V,2,4) = - Fmatrices(V,Index_34) * S2(V)
         Zmatrices(V,3,2) =   S2(V) * HELP2C1 + C2(V) * HELP3S1
         Zmatrices(V,3,3) = - S2(V) * HELP2S1 + C2(V) * HELP3C1
         Zmatrices(V,3,4) =   Fmatrices(V,Index_34) * C2(V)
         Zmatrices(V,4,2) = - Fmatrices(V,Index_34) * S1(V)
         Zmatrices(V,4,3) = - Fmatrices(V,Index_34) * C1(V)
         Zmatrices(V,4,4) =   Fmatrices(V,Index_44)
      ENDDO
   ENDIF

!  Done

   RETURN
END SUBROUTINE vfzmat_PhasMat1

!

subroutine vfzmat_PhasMat2 &
   ( max_geoms, maxlayers, nlayers, nstokes, n_geoms,      & ! inputs
     Sunlight, C1, S1, C2, S2, Exist_Fmatrices, Fmatrices, & ! Inputs
     Zmatrices )

!  Stand-alone routine to develop Z-matrices from F-Matrix inputs for
!  various goemetry combinations, upwelling or downwelling.

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  Transformation from scattering plan (Fmat) to Meriodional Plane (Zmat)
!   follows the rotations given in the VLIDORT code.

!  Remarks
!  -------

!  1. Angles and numbers should be set up beforehand as for VLIDORT.
!  2. Existence flag introduced as input, 9/19/16

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input
!  -----

!  dimensions

   INTEGER, INTENT (IN) :: maxlayers, max_geoms

!  numbers

   INTEGER, INTENT (IN) :: nlayers, nstokes, n_geoms

!  Flag for Sunlight
 
   LOGICAL, intent(in)   :: Sunlight

!  Rotational variables

   REAL(KIND=dpk), INTENT (IN)   :: C1(max_geoms), S1(max_geoms)
   REAL(KIND=dpk), INTENT (IN)   :: C2(max_geoms), S2(max_geoms)

!  Fmatrices

   REAL(KIND=dpk), INTENT (IN) :: Fmatrices(max_geoms,maxlayers,6)
   Logical       , INTENT (IN) :: Exist_Fmatrices(maxlayers)

!  output
!  ------

!  Zmatrices output, initialized globally beforehand

   REAL(KIND=dpk), INTENT (InOut) :: Zmatrices(max_geoms,maxlayers,4,4)

!  local variables
!  ---------------

   INTEGER            :: n, v, index_11, index_12, index_22, index_33, index_34, index_44 
   REAL    (KIND=dpk) :: HELP2C1, HELP2S1, HELP3C1, HELP3S1

!  New indexing consistent with Tmatrix output

  index_11 = 1
  index_22 = 2
  index_33 = 3
  index_44 = 4
  index_12 = 5
  index_34 = 6

!  Nstokes = 1, just copy F11 to Z(1,1) element and return (no rotation)

   if ( nstokes .eq. 1 ) then
     do n  = 1, nlayers
       if ( Exist_Fmatrices(n) ) then
         do v = 1, n_geoms
           Zmatrices(v,n,1,1) = Fmatrices(v,n,index_11) 
         enddo
       endif
     enddo
     return
   endif

!  Z_matrices, general case
!  ------------------------
 
!  For sunlight, only need the first column of Z matrix

   IF ( SUNLIGHT ) THEN
     DO N = 1, NLAYERS
       IF ( Exist_Fmatrices(n) ) then
         DO V = 1, N_GEOMS
           Zmatrices(V,N,1,1) =   Fmatrices(V,N,Index_11)
           Zmatrices(V,N,2,1) = - Fmatrices(V,N,Index_12) * C2(V)
           Zmatrices(V,N,3,1) =   Fmatrices(V,N,Index_12) * S2(V)
         ENDDO
       ENDIF
     ENDDO
   ENDIF

!  For general case. Not tested as of 15 April 2005.

   IF ( .NOT. SUNLIGHT ) THEN
     DO N = 1, NLAYERS
       IF ( Exist_Fmatrices(n) ) then
         DO V = 1, N_GEOMS
           Zmatrices(V,N,1,1) =   Fmatrices(V,N,Index_11)
           Zmatrices(V,N,2,1) = - Fmatrices(V,N,Index_12) * C2(V)
           Zmatrices(V,N,3,1) =   Fmatrices(V,N,Index_12) * S2(V)
!  This code untested
           HELP2C1 = Fmatrices(V,N,Index_22) * C1(V)
           HELP2S1 = Fmatrices(V,N,Index_22) * S1(V)
           HELP3C1 = Fmatrices(V,N,Index_33) * C1(V)
           HELP3S1 = Fmatrices(V,N,Index_33) * S1(V)
           Zmatrices(V,N,1,2) =   Fmatrices(V,N,Index_12) * C1(V)
           Zmatrices(V,N,1,3) = - Fmatrices(V,N,Index_12) * S1(V)
           Zmatrices(V,N,2,2) =   C2(V) * HELP2C1 - S2(V) * HELP3S1
           Zmatrices(V,N,2,3) = - C2(V) * HELP2S1 - S2(V) * HELP3C1
           Zmatrices(V,N,2,4) = - Fmatrices(V,N,Index_34) * S2(V)
           Zmatrices(V,N,3,2) =   S2(V) * HELP2C1 + C2(V) * HELP3S1
           Zmatrices(V,N,3,3) = - S2(V) * HELP2S1 + C2(V) * HELP3C1
           Zmatrices(V,N,3,4) =   Fmatrices(V,N,Index_34) * C2(V)
           Zmatrices(V,N,4,2) = - Fmatrices(V,N,Index_34) * S1(V)
           Zmatrices(V,N,4,3) = - Fmatrices(V,N,Index_34) * C1(V)
           Zmatrices(V,N,4,4) =   Fmatrices(V,N,Index_44)
         ENDDO
       ENDIF
     ENDDO
   ENDIF

!  Done

   RETURN
END SUBROUTINE vfzmat_PhasMat2

!  End module

End module vfzmat_PhasMat_m

