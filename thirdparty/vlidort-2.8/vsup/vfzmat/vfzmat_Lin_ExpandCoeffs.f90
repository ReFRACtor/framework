
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

module vfzmat_Lin_ExpandCoeffs_m

public

contains

SUBROUTINE vfzmat_Lin_ExpandCoeffs2 &
   ( maxcoeffs, maxlayers, maxgeoms, maxwfs, nstokes, ncoeffs, nlayers, ngeoms, & ! Inputs
     do_sunlight, Qvary, Qnums, Lvarymoms, L_greekmatvec, genspher,             & ! Inputs
     L_fmatvec )                                                                  ! Outputs

!  Stand-alone routine to define the nonzero entries of the linearized F-Matrix using
!  the generalized spherical function expansion coefficients of the nonzero
!  entries of the Greek matrix which have also been linearized

!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!  Programmed 3/2/20 by M. Christi
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

   IMPLICIT NONE

!  Parameters
!  ----------

!  Precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Input
!  -----

!  Control

   INTEGER,   INTENT(IN) :: maxcoeffs, maxlayers, maxgeoms, maxwfs
   INTEGER,   INTENT(IN) :: nstokes, ncoeffs, nlayers, ngeoms

   LOGICAL,   INTENT(IN) :: do_sunlight

   LOGICAL,   INTENT(IN) :: Qvary ( maxlayers ), Lvarymoms ( maxlayers, maxwfs )
   INTEGER,   INTENT(IN) :: Qnums ( maxlayers )

!  Coefficients of the nonzero entries of the 4x4 Greek Matrix (vector version)

   REAL(dpk), INTENT(IN) :: L_greekmatvec ( maxwfs, 0:maxcoeffs, maxlayers, 16 )

!  Generalized spherical functions.
!    Rotations(1-4) = C1, S1, C2, S2

   REAL(dpk), INTENT(IN) :: genspher ( 0:maxcoeffs, 4, maxgeoms )

!  Output
!  ------

!  Nonzero entries of the F-matrix (vector version)

   REAL(dpk), INTENT(OUT) :: L_fmatvec ( maxwfs, maxlayers, maxgeoms, 6 )

!  Local variables
!  ---------------

   REAL(dpk), PARAMETER :: d_zero  = 0.0_dpk, d_half  = 0.5_dpk

   INTEGER   :: g1, g2, g3, g4, g5, g6, gmask(6), n, q, v
   REAL(dpk) :: L_fmat3, L_fmat4, L_gmsum(0:maxcoeffs), L_gmdif(0:maxcoeffs)

!  Initialization

   L_fmatvec = d_zero

!  Define some local variables

   !Entries in this version of gmask correspond to the following positions
   !in the 4x4 Greek Matrix:
   !          11, 12, 22, 33, 34, 44         
   gmask = (/  1,  2,  6, 11, 12, 16 /)

!  Get F-matrix

   if ( nstokes .eq. 1 ) then

!  Scalar case

     g1 = gmask(1)
     do v = 1, ngeoms
       do n = 1, nlayers
         if ( Qvary(n) ) then
           do q = 1, Qnums(n)
             if ( Lvarymoms(n,q) ) then
               L_fmatvec(q,n,v,1) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
             endif
           enddo
         endif
       enddo
     enddo

   else

!  Vector cases (nstokes > 1)

     if ( do_sunlight ) then

!    Vector - with sunlight

       g1 = gmask(1); g2 = gmask(2)
       do v = 1, ngeoms
         do n = 1, nlayers
           if ( Qvary(n) ) then
             do q = 1, Qnums(n)
               if ( Lvarymoms(n,q) ) then
                 L_fmatvec(q,n,v,1) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
                 L_fmatvec(q,n,v,2) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g2),genspher(0:ncoeffs,2,v))
               endif
             enddo
           endif
         enddo
       enddo

     else

!    Vector - general case with possible linear and circular polarization
!      prepare to use full 4x4 matrix; code introduced but not tested, 05 October 2010

       g1 = gmask(1); g2 = gmask(2); g3 = gmask(3); g4 = gmask(4)
       do v = 1, ngeoms
         do n = 1, nlayers
           if ( Qvary(n) ) then
             do q = 1, Qnums(n)
               if ( Lvarymoms(n,q) ) then
                 L_fmatvec(q,n,v,1) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
                 L_fmatvec(q,n,v,2) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g2),genspher(0:ncoeffs,2,v))

                 L_gmsum(0:ncoeffs) = L_greekmatvec(q,0:ncoeffs,n,g3) + L_greekmatvec(q,0:ncoeffs,n,g4)
                 L_fmat3            = DOT_PRODUCT(L_gmsum(0:ncoeffs),genspher(0:ncoeffs,3,v))

                 L_gmdif(0:ncoeffs) = L_greekmatvec(q,0:ncoeffs,n,g3) - L_greekmatvec(q,0:ncoeffs,n,g4)
                 L_fmat4            = DOT_PRODUCT(L_gmdif(0:ncoeffs),genspher(0:ncoeffs,4,v))

                 L_fmatvec(q,n,v,3) = ( L_fmat3 + L_fmat4 ) * d_half
                 L_fmatvec(q,n,v,4) = ( L_fmat3 - L_fmat4 )

                 if ( nstokes.eq.4 ) then
                   g5 = gmask(5); g6 = gmask(6)
                   L_fmatvec(q,n,v,5) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g5),genspher(0:ncoeffs,2,v))
                   L_fmatvec(q,n,v,6) = DOT_PRODUCT(L_greekmatvec(q,0:ncoeffs,n,g6),genspher(0:ncoeffs,1,v))
                 endif
               endif
             enddo
           endif 
         enddo
       enddo

     endif  !End sunlight IF block
   endif  !End stokes IF block

END SUBROUTINE vfzmat_Lin_ExpandCoeffs2

!  End module

End Module vfzmat_Lin_ExpandCoeffs_m

