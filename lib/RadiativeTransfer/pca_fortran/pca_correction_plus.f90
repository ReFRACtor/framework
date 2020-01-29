module pca_correction_plus_m

!  This Module :

!   Version 1. Original code. 19 November 2011
!   Version 2. Through-differentiation of PCA solution. 1-3 February 2012
!   Version 3. Robustness and Recoding, 2M/3M.  08 November 2012
!   Version 4. Multigeometry.       05 April    2013

!  History of PCA developments:

!   27 January  2012. Alternative EOF Scheme, using Eigenmethods.
!   03 February 2012. Column  Linearization of the EOF correction
!   06 February 2012. Profile Linearization of the EOF correction
!   24 May      2012. Albedo Closure included (new routine) 

!  Subroutines:

!    pca_3M_correction_pplus  --> Correction with profile Jacobian,  LD/2S/FO (3 Intensities)
!    pca_3M_correction_cplus  --> Correction with Column  Jacobian,  LD/2S/FO (3 Intensities)
!    pca_3M_correction_csplus --> Correction with Column and Albedo, LD/2S/FO (3 Intensities)

!    pca_2M_correction_pplus  --> Correction with profile Jacobian,  LD/FO (2 Intensities)
!    pca_2M_correction_cplus  --> Correction with Column  Jacobian,  LD/FO (2 Intensities)
!    pca_2M_correction_csplus --> Correction with Column and Albedo, LD/FO (2 Intensities)

!  Remark:
!   Max_eofs_2p1 = 2 * Max_eofs + 1
!   Max_eofs_4p1 = 4 * Max_eofs + 1      Should be the input for the 3OD routine

!  NOTE (12 April 2013)
!  ====================

!   THERE IS CURRENTLY NO LINEARIZED ROUTINES USING THE 3OD (third-order differencing)

!  precision parameters

!      INTEGER, PARAMETER :: sp     = selected_real_kind(6)
      INTEGER, PARAMETER :: dp     = selected_real_kind(15)

public

contains

subroutine pca_3M_correction_pplus &
          ( max_eofs_2p1, Max_geoms, maxlayers, max_atmoswfs, do_Jacobian, & ! Input dimensioning 
            n_eofs, npoints, ngeoms, nlayers, n_atmoswfs, Kvary, Kpars,    & ! Input control
            PrinComps, LP_PrinComps,                                       & ! Input Principal Components 
            intensity_LD_bin,    intensity_2S_bin,    intensity_FO_bin,    & ! Input (Bin-Intensities)
            LP_Jacobians_LD_bin, LP_Jacobians_2S_bin, LP_Jacobians_FO_bin, & ! Input (Linearized)
            intensity_corrfacs,  LP_Jacobians_corrfacs )                     ! Output (Corrections)

!  Corrections with Profile (layer) Jacobians. 06 February 2012.

      IMPLICIT NONE

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms, maxlayers, max_atmoswfs

!  Control

      LOGICAL      , intent(in) :: do_Jacobian
      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms, nlayers, n_atmoswfs

!  Profile lienarization control

      LOGICAL      , intent(in) :: Kvary(maxlayers)
      INTEGER      , intent(in) :: Kpars(maxlayers)

!  Principal components, and Jacobians thereunto. Pre-allocated.

      REAL(kind=dp), intent(in) :: PrinComps    (n_eofs,npoints)
      REAL(kind=dp), intent(in) :: LP_PrinComps (n_eofs,npoints,nlayers,n_atmoswfs) 

!  Bin-Intensities 
!     (2S_BIN = 2-stream, LD_BIN = LIDORT N-stream, FO_BIN = single-scattering)

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_2S_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,Max_geoms)

!  Linearized Bin_intensities

      REAL(kind=dp), intent(in) :: LP_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,maxlayers,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LP_JACOBIANS_2S_BIN(max_eofs_2p1,Max_geoms,maxlayers,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LP_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,maxlayers,max_atmoswfs)

!  outputs; Intensity corrections and their linearizations. Pre-allocated.

      REAL(kind=dp), intent(out) :: INTENSITY_CORRFACS   (npoints,ngeoms)
      REAL(kind=dp), intent(out) :: LP_JACOBIANS_CORRFACS(npoints,ngeoms,nlayers,n_atmoswfs)

!  local variables
!  ---------------

      INTEGER      :: I, M, MM, MP, Q, K, V

      REAL(kind=dp) :: ID0, IEOF, TERM_1(N_EOFS), TERM_2(N_EOFS), X, XSQ
      REAL(kind=dp) :: IDPLUS, IDMINUS, HI_SS, LO_SS, HI_SS_MP, LO_SS_MP, HI_SS_MM, LO_SS_MM

      REAL(kind=dp) :: LP_ID0(NLAYERS,N_ATMOSWFS), LP_IEOF(NLAYERS,N_ATMOSWFS)
      REAL(kind=dp) :: LP_TERM_1(N_EOFS,NLAYERS,N_ATMOSWFS), LP_TERM_2(N_EOFS,NLAYERS,N_ATMOSWFS)
      REAL(kind=dp) :: LP_IDPLUS, LP_IDMINUS, LP_HI_SS, LP_LO_SS, LP_X, LP_XSQ

      REAL(kind=dp), parameter ::  HALF = 0.5_dp, two = 2.0_dp

!  Zero output
!  -----------

!    PLACEHOLDER - is it really necessary?

!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  mean value
!  ----------

!  Mean-value calculation

      HI_SS = INTENSITY_LD_BIN(1,V)+INTENSITY_FO_BIN(1,V)
      LO_SS = INTENSITY_2S_BIN(1,V)+INTENSITY_FO_BIN(1,V)
      ID0   = LOG ( HI_SS / LO_SS )

!  Linearized mean value calculation

      if ( do_Jacobian ) then
         LP_ID0 = 0.0_dp    ! Safety-first zeroing
         do k = 1, nlayers
            if ( kvary(k) ) then
               do Q = 1, kpars(k)
                  LP_HI_SS    = LP_JACOBIANS_LD_BIN(1,V,K,Q) + LP_JACOBIANS_FO_BIN(1,V,K,Q)
                  LP_LO_SS    = LP_JACOBIANS_2S_BIN(1,V,K,Q) + LP_JACOBIANS_FO_BIN(1,V,K,Q)
                  LP_ID0(K,Q) = ( LP_HI_SS / HI_SS ) - ( LP_LO_SS / LO_SS )
               enddo
            endif
         enddo
      endif

!  Plus/minus values
!  -----------------

      DO M = 1, N_EOFS

!  Indices

         MP = 2*M ; MM = MP + 1

!  Plus/Minus values. Finite Differences TERM_1 and TERM_2

         HI_SS_MP   = INTENSITY_LD_BIN(MP,V)+INTENSITY_FO_BIN(MP,V)
         LO_SS_MP   = INTENSITY_2S_BIN(MP,V)+INTENSITY_FO_BIN(MP,V)
         IDPLUS     = LOG ( HI_SS_MP / LO_SS_MP )
         HI_SS_MM   = INTENSITY_LD_BIN(MM,V)+INTENSITY_FO_BIN(MM,V)
         LO_SS_MM   = INTENSITY_2S_BIN(MM,V)+INTENSITY_FO_BIN(MM,V)
         IDMINUS    = LOG ( HI_SS_MM / LO_SS_MM )
         TERM_1(M)  = ( IDPLUS - IDMINUS )            * HALF
         TERM_2(M)  = ( IDPLUS + IDMINUS - TWO*ID0 )  * HALF

!  Linearizations

         if ( do_Jacobian ) then
            do k = 1, nlayers
               if ( kvary(k) ) then
                  do Q = 1, kpars(k)
                     LP_HI_SS   = LP_JACOBIANS_LD_BIN(MP,V,k,Q) + LP_JACOBIANS_FO_BIN(MP,V,k,Q)
                     LP_LO_SS   = LP_JACOBIANS_2S_BIN(MP,V,k,Q) + LP_JACOBIANS_FO_BIN(MP,V,k,Q)
                     LP_IDPLUS  = ( LP_HI_SS / HI_SS_MP ) - ( LP_LO_SS / LO_SS_MP )
                     LP_HI_SS   = LP_JACOBIANS_LD_BIN(MM,V,k,Q) + LP_JACOBIANS_FO_BIN(MM,V,k,Q)
                     LP_LO_SS   = LP_JACOBIANS_2S_BIN(MM,V,k,Q) + LP_JACOBIANS_FO_BIN(MM,V,k,Q)
                     LP_IDMINUS = ( LP_HI_SS / HI_SS_MM ) - ( LP_LO_SS / LO_SS_MM )
                     LP_TERM_1(M,k,Q)  = ( LP_IDPLUS - LP_IDMINUS )                   * HALF
                     LP_TERM_2(M,k,Q)  = ( LP_IDPLUS + LP_IDMINUS - TWO*LP_ID0(k,Q) ) * HALF
                  enddo
               endif
            enddo
         endif

!  End EOF loop

      ENDDO

!  set correction for all wavelengths
!  ----------------------------------

!  Including Jacobians. Return when done.

      if ( do_Jacobian ) then
         DO I = 1, npoints
            IEOF = ID0 ; LP_IEOF = LP_ID0
            DO M = 1, N_EOFS
               X   = PrinComps(M,I) ; XSQ   = X * X
               IEOF   = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
               do k = 1, nlayers
                  if ( kvary(k) ) then
                     do Q = 1, kpars(k)
                        LP_X = LP_PrinComps(M,I,K,Q) ; LP_XSQ = TWO * X * LP_X
                        LP_IEOF(k,Q) = LP_IEOF(k,Q) + LP_TERM_1(M,k,Q) * X   + TERM_1(M) * LP_X   &
                                                    + LP_TERM_2(M,k,Q) * XSQ + TERM_2(M) * LP_XSQ
                     enddo
                  endif
               enddo
            ENDDO
            IEOF = EXP(IEOF)
            INTENSITY_CORRFACS(I,V) = IEOF
            do k = 1, nlayers
               LP_JACOBIANS_CORRFACS(I,V,k,1:n_atmoswfs) = IEOF * LP_IEOF(k,1:n_atmoswfs)
            enddo
         ENDDO
         RETURN          ! RETURN
      endif

!  No Jacobians, just intensity

      DO I = 1, npoints
         IEOF = ID0 
         DO M = 1, N_EOFS
            X    = PrinComps(M,I) ; XSQ   = X * X
            IEOF = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
         ENDDO
         IEOF = EXP(IEOF)
         INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_3M_correction_pplus


subroutine pca_3M_correction_cplus &
          ( max_eofs_2p1, Max_geoms, max_atmoswfs,                          & ! Input dimensioning 
            do_Jacobian, n_eofs, npoints, ngeoms, n_atmoswfs,               & ! Input control
            PrinComps, LC_PrinComps,                                        & ! Input Principal Components 
            intensity_LD_bin,    intensity_2S_bin,    intensity_FO_bin,     & ! Input (Bin-Intensities)
            LC_Jacobians_LD_bin, LC_Jacobians_2S_bin, LC_Jacobians_FO_bin,  & ! Input (Linearized
            intensity_corrfacs, LC_Jacobians_corrfacs )                       ! Output (Corrections)

!  Column (bulk) Jacobians. 03 February 2012.

      IMPLICIT NONE

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms, max_atmoswfs

!  Control

      LOGICAL      , intent(in) :: do_Jacobian
      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms, n_atmoswfs

!  Principal components, and Jacobians thereunto. Pre-allocated.

      REAL(kind=dp), intent(in) :: PrinComps    (n_eofs,npoints)
      REAL(kind=dp), intent(in) :: LC_PrinComps (n_eofs,npoints,n_atmoswfs) 

!  Bin-Intensities 
!     (2S_BIN = 2-stream, LD_BIN = LIDORT N-stream, FO_BIN = single-scattering)

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_2S_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,Max_geoms)

!  Linearized Bin_intensities

      REAL(kind=dp), intent(in) :: LC_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LC_JACOBIANS_2S_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LC_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)

!  outputs; Intensity corrections and their linearizations. Pre-allocated.

      REAL(kind=dp), intent(out) :: INTENSITY_CORRFACS   (npoints,ngeoms)
      REAL(kind=dp), intent(out) :: LC_JACOBIANS_CORRFACS(npoints,ngeoms,n_atmoswfs)

!  local variables
!  ---------------

      INTEGER      :: I, M, MM, MP, Q, V

      REAL(kind=dp) :: ID0, IEOF, TERM_1(N_EOFS), TERM_2(N_EOFS),  X, XSQ
      REAL(kind=dp) :: IDPLUS, IDMINUS, HI_SS, LO_SS, HI_SS_MP, LO_SS_MP, HI_SS_MM, LO_SS_MM

      REAL(kind=dp) :: LC_ID0(N_ATMOSWFS), LC_IEOF(N_ATMOSWFS)
      REAL(kind=dp) :: LC_TERM_1(N_EOFS,N_ATMOSWFS), LC_TERM_2(N_EOFS,N_ATMOSWFS)
      REAL(kind=dp) :: LC_IDPLUS, LC_IDMINUS, LC_HI_SS, LC_LO_SS, LC_X, LC_XSQ

      REAL(kind=dp), parameter ::  HALF = 0.5_dp, two = 2.0_dp

!  Zero output
!  -----------

!    PLACEHOLDER,

!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  mean value
!  ----------

!  mean value calculation

      HI_SS = INTENSITY_LD_BIN(1,V)+INTENSITY_FO_BIN(1,V)
      LO_SS = INTENSITY_2S_BIN(1,V)+INTENSITY_FO_BIN(1,V)
      ID0 = DLOG ( HI_SS / LO_SS )

!  Linearized mean value calculation

      if ( do_Jacobian ) then
         do Q = 1, n_atmoswfs
            LC_HI_SS  = LC_JACOBIANS_LD_BIN(1,V,Q) + LC_JACOBIANS_FO_BIN(1,V,Q)
            LC_LO_SS  = LC_JACOBIANS_2S_BIN(1,V,Q) + LC_JACOBIANS_FO_BIN(1,V,Q)
            LC_ID0(Q) = ( LC_HI_SS / HI_SS ) - ( LC_LO_SS / LO_SS )
         enddo
      endif

!  Plus/minus values
!  -----------------

      DO M = 1, N_EOFS

!  Indices

        MP = 2*M ; MM = MP + 1

!  Plus/Minus values. Finite Differences TERM_1 and TERM_2

        HI_SS_MP   = INTENSITY_LD_BIN(MP,V)+INTENSITY_FO_BIN(MP,V)
        LO_SS_MP   = INTENSITY_2S_BIN(MP,V)+INTENSITY_FO_BIN(MP,V)
        IDPLUS     = LOG ( HI_SS_MP / LO_SS_MP )

        HI_SS_MM   = INTENSITY_LD_BIN(MM,V)+INTENSITY_FO_BIN(MM,V)
        LO_SS_MM   = INTENSITY_2S_BIN(MM,V)+INTENSITY_FO_BIN(MM,V)
        IDMINUS    = LOG ( HI_SS_MM / LO_SS_MM )

        TERM_1(M) = ( IDPLUS - IDMINUS )             * HALF
        TERM_2(M) = ( IDPLUS + IDMINUS - TWO * ID0 ) * HALF

!  Linearizations

        if ( do_Jacobian ) then
           do Q = 1, n_atmoswfs
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MP,V,Q)+LC_JACOBIANS_FO_BIN(MP,V,Q)
              LC_LO_SS   = LC_JACOBIANS_2S_BIN(MP,V,Q)+LC_JACOBIANS_FO_BIN(MP,V,Q)
              LC_IDPLUS  = ( LC_HI_SS / HI_SS_MP ) - ( LC_LO_SS / LO_SS_MP )
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MM,V,Q)+LC_JACOBIANS_FO_BIN(MM,V,Q)
              LC_LO_SS   = LC_JACOBIANS_2S_BIN(MM,V,Q)+LC_JACOBIANS_FO_BIN(MM,V,Q)
              LC_IDMINUS = ( LC_HI_SS / HI_SS_MM ) - ( LC_LO_SS / LO_SS_MM )
              LC_TERM_1(M,Q)  = ( LC_IDPLUS - LC_IDMINUS )                * HALF
              LC_TERM_2(M,Q)  = ( LC_IDPLUS + LC_IDMINUS - TWO*LC_ID0(Q) ) * HALF
           enddo
        endif

!  End EOF loop

      ENDDO

!  set correction for all wavelengths
!  ----------------------------------

!  Including Jacobian. Return when done.

      if ( do_Jacobian ) then
         DO I = 1, npoints
            IEOF = ID0 ; LC_IEOF = LC_ID0
            DO M = 1, N_EOFS
               X   = PrinComps(M,I)   ; XSQ   = X * X
               IEOF   = IEOF   + TERM_1(M)   * X + TERM_2(M) * XSQ
               do q = 1, n_atmoswfs
                  LC_X = LC_PrinComps(M,I,Q) ; LC_XSQ = two * X * LC_X
                  LC_IEOF(Q) = LC_IEOF(Q) + LC_TERM_1(M,Q) * X   + TERM_1(M) * LC_X   &
                                          + LC_TERM_2(M,Q) * XSQ + TERM_2(M) * LC_XSQ
               enddo
            ENDDO
            IEOF = EXP(IEOF)
            INTENSITY_CORRFACS(I,V)                 = IEOF
            LC_JACOBIANS_CORRFACS(I,V,1:n_atmoswfs) = IEOF * LC_IEOF(1:n_atmoswfs)
         ENDDO
         RETURN          ! RETURN
      endif

!  No Jacobians, just intensity

      DO I = 1, npoints
         IEOF = ID0 
         DO M = 1, N_EOFS
            X    = PrinComps(M,I) ; XSQ   = X * X
            IEOF = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
         ENDDO
         IEOF = EXP(IEOF)
         INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_3M_correction_cplus

subroutine pca_3M_correction_csplus  &
          ( max_eofs_2p1, Max_geoms, max_atmoswfs, max_surfacewfs, do_Jacobian, & ! Input dimensioning 
            do_SJacobian, n_eofs, npoints, ngeoms, n_atmoswfs, n_surfacewfs,    & ! Input control
            PrinComps, LC_PrinComps, LS_PrinComps,                              & ! Input Principal Components 
            intensity_LD_bin,    intensity_2S_bin,    intensity_FO_bin,         & ! Input (Bin-Intensities)
            LC_Jacobians_LD_bin, LC_Jacobians_2S_bin, LC_Jacobians_FO_bin,      & ! Input (Linearized
            LS_Jacobians_LD_bin, LS_Jacobians_2S_bin, LS_Jacobians_FO_bin,      & ! Input (Linearized
            intensity_corrfacs, LC_Jacobians_corrfacs, LS_Jacobians_corrfacs )    ! Output (Corrections)

!  Column (bulk) Jacobians. 03 February 2012.

      IMPLICIT NONE

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms, max_atmoswfs, max_surfacewfs

!  Control

      LOGICAL      , intent(in) :: do_Jacobian, do_SJacobian
      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms, n_atmoswfs, n_surfacewfs

!  Principal components, and Jacobians thereunto. Pre-allocated.

      REAL(kind=dp), intent(in) :: PrinComps    (n_eofs,npoints)
      REAL(kind=dp), intent(in) :: LC_PrinComps (n_eofs,npoints,n_atmoswfs) 
      REAL(kind=dp), intent(in) :: LS_PrinComps (n_eofs,npoints,n_surfacewfs) 

!  Bin-Intensities 
!     (2S_BIN = 2-stream, LD_BIN = LIDORT N-stream, FO_BIN = single-scattering)

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_2S_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,Max_geoms)

!  Linearized Bin-intensities

      REAL(kind=dp), intent(in) :: LC_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LC_JACOBIANS_2S_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LC_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)

      REAL(kind=dp), intent(in) :: LS_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,max_surfacewfs)
      REAL(kind=dp), intent(in) :: LS_JACOBIANS_2S_BIN(max_eofs_2p1,Max_geoms,max_surfacewfs)
      REAL(kind=dp), intent(in) :: LS_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,max_surfacewfs)

!  outputs; Intensity corrections and their linearizations. Pre-allocated.

      REAL(kind=dp), intent(out) :: INTENSITY_CORRFACS   (npoints,ngeoms)
      REAL(kind=dp), intent(out) :: LC_JACOBIANS_CORRFACS(npoints,ngeoms,n_atmoswfs)
      REAL(kind=dp), intent(out) :: LS_JACOBIANS_CORRFACS(npoints,ngeoms,n_surfacewfs)

!  local variables
!  ---------------

      INTEGER      :: I, M, MM, MP, Q, V

      REAL(kind=dp) :: ID0, IEOF, TERM_1(N_EOFS), TERM_2(N_EOFS),  X, XSQ
      REAL(kind=dp) :: IDPLUS, IDMINUS, HI_SS, LO_SS, HI_SS_MP, LO_SS_MP, HI_SS_MM, LO_SS_MM

      REAL(kind=dp) :: LC_ID0(N_ATMOSWFS), LC_IEOF(N_ATMOSWFS)
      REAL(kind=dp) :: LC_TERM_1(N_EOFS,N_ATMOSWFS), LC_TERM_2(N_EOFS,N_ATMOSWFS)
      REAL(kind=dp) :: LC_IDPLUS, LC_IDMINUS, LC_HI_SS, LC_LO_SS, LC_X, LC_XSQ

      REAL(kind=dp) :: LS_ID0(N_SURFACEWFS), LS_IEOF(N_SURFACEWFS)
      REAL(kind=dp) :: LS_TERM_1(N_EOFS,N_SURFACEWFS), LS_TERM_2(N_EOFS,N_SURFACEWFS)
      REAL(kind=dp) :: LS_IDPLUS, LS_IDMINUS, LS_HI_SS, LS_LO_SS, LS_X, LS_XSQ

      REAL(kind=dp), parameter ::  HALF = 0.5_dp, two = 2.0_dp

!  Zero output
!  -----------

!    PLACEHOLDER - is it really necessary?

!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  mean value
!  ----------

!  mean value calculation

      HI_SS = INTENSITY_LD_BIN(1,V)+INTENSITY_FO_BIN(1,V)
      LO_SS = INTENSITY_2S_BIN(1,V)+INTENSITY_FO_BIN(1,V)
      ID0   = LOG ( HI_SS / LO_SS )

!  Linearized mean value calculation

      if ( do_Jacobian ) then
         do Q = 1, n_atmoswfs
            LC_HI_SS  = LC_JACOBIANS_LD_BIN(1,V,Q)+LC_JACOBIANS_FO_BIN(1,V,Q)
            LC_LO_SS  = LC_JACOBIANS_2S_BIN(1,V,Q)+LC_JACOBIANS_FO_BIN(1,V,Q)
            LC_ID0(Q) = ( LC_HI_SS / HI_SS ) - ( LC_LO_SS / LO_SS )
         enddo
      endif

      if ( do_SJacobian ) then
         do Q = 1, n_surfacewfs
            LS_HI_SS  = LS_JACOBIANS_LD_BIN(1,V,Q)+LS_JACOBIANS_FO_BIN(1,V,Q)
            LS_LO_SS  = LS_JACOBIANS_2S_BIN(1,V,Q)+LS_JACOBIANS_FO_BIN(1,V,Q)
            LS_ID0(Q) = ( LS_HI_SS / HI_SS ) - ( LS_LO_SS / LO_SS )
         enddo
      endif

!  Plus/minus values
!  -----------------

      DO M = 1, N_EOFS

!  Indices

        MP = 2*M ; MM = MP + 1

!  Plus/Minus values. Finite Differences TERM_1 and TERM_2

        HI_SS_MP   = INTENSITY_LD_BIN(MP,V)+INTENSITY_FO_BIN(MP,V)
        LO_SS_MP   = INTENSITY_2S_BIN(MP,V)+INTENSITY_FO_BIN(MP,V)
        IDPLUS     = LOG ( HI_SS_MP / LO_SS_MP )

        HI_SS_MM   = INTENSITY_LD_BIN(MM,V)+INTENSITY_FO_BIN(MM,V)
        LO_SS_MM   = INTENSITY_2S_BIN(MM,V)+INTENSITY_FO_BIN(MM,V)
        IDMINUS    = LOG ( HI_SS_MM / LO_SS_MM )

        TERM_1(M) = ( IDPLUS - IDMINUS )           * HALF
        TERM_2(M) = ( IDPLUS + IDMINUS - TWO*ID0 ) * HALF

!  Column Linearizations

        if ( do_Jacobian ) then
           do Q = 1, n_atmoswfs
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MP,V,Q)+LC_JACOBIANS_FO_BIN(MP,V,Q)
              LC_LO_SS   = LC_JACOBIANS_2S_BIN(MP,V,Q)+LC_JACOBIANS_FO_BIN(MP,V,Q)
              LC_IDPLUS  = ( LC_HI_SS / HI_SS_MP ) - ( LC_LO_SS / LO_SS_MP )
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MM,V,Q)+LC_JACOBIANS_FO_BIN(MM,V,Q)
              LC_LO_SS   = LC_JACOBIANS_2S_BIN(MM,V,Q)+LC_JACOBIANS_FO_BIN(MM,V,Q)
              LC_IDMINUS = ( LC_HI_SS / HI_SS_MM ) - ( LC_LO_SS / LO_SS_MM )
              LC_TERM_1(M,Q)  = ( LC_IDPLUS - LC_IDMINUS )                 * HALF
              LC_TERM_2(M,Q)  = ( LC_IDPLUS + LC_IDMINUS - TWO*LC_ID0(Q) ) * HALF
           enddo
        endif

!  surface Linearizations

        if ( do_SJacobian ) then
           do Q = 1, n_surfacewfs
              LS_HI_SS   = LS_JACOBIANS_LD_BIN(MP,V,Q)+LS_JACOBIANS_FO_BIN(MP,V,Q)
              LS_LO_SS   = LS_JACOBIANS_2S_BIN(MP,V,Q)+LS_JACOBIANS_FO_BIN(MP,V,Q)
              LS_IDPLUS  = ( LS_HI_SS / HI_SS_MP ) - ( LS_LO_SS / LO_SS_MP )
              LS_HI_SS   = LS_JACOBIANS_LD_BIN(MM,V,Q)+LS_JACOBIANS_FO_BIN(MM,V,Q)
              LS_LO_SS   = LS_JACOBIANS_2S_BIN(MM,V,Q)+LS_JACOBIANS_FO_BIN(MM,V,Q)
              LS_IDMINUS = ( LS_HI_SS / HI_SS_MM ) - ( LS_LO_SS / LO_SS_MM )
              LS_TERM_1(M,Q)  = ( LS_IDPLUS - LS_IDMINUS )                 * HALF
              LS_TERM_2(M,Q)  = ( LS_IDPLUS + LS_IDMINUS - TWO*LS_ID0(Q) ) * HALF
           enddo
        endif

!  End EOF loop

      ENDDO

!  set correction for all wavelengths
!  ----------------------------------

!  Including Jacobians (BOTH). Return when done.

      if ( do_Jacobian .and. do_SJacobian  ) then
         DO I = 1, npoints
            IEOF = ID0 ; LC_IEOF(1:n_atmoswfs) = LC_ID0(1:n_atmoswfs); LS_IEOF(1:n_surfacewfs) = LS_ID0(1:n_surfacewfs)
            DO M = 1, N_EOFS
               X   = PrinComps(M,I)   ; XSQ   = X * X
               IEOF   = IEOF   + TERM_1(M)   * X + TERM_2(M) * XSQ
               do q = 1, n_atmoswfs
                  LC_X = LC_PrinComps(M,I,Q) ; LC_XSQ = TWO * X * LC_X
                  LC_IEOF(Q) = LC_IEOF(Q) + LC_TERM_1(M,Q) * X   + TERM_1(M) * LC_X   &
                                          + LC_TERM_2(M,Q) * XSQ + TERM_2(M) * LC_XSQ
               enddo
               do q = 1, n_surfacewfs
                  LS_X = LS_PrinComps(M,I,Q) ; LS_XSQ = TWO * X * LS_X
                  LS_IEOF(Q) = LS_IEOF(Q) + LS_TERM_1(M,Q) * X   + TERM_1(M) * LS_X   &
                                          + LS_TERM_2(M,Q) * XSQ + TERM_2(M) * LS_XSQ
               enddo
            ENDDO
            IEOF = EXP(IEOF)
            INTENSITY_CORRFACS(I,V)             = IEOF
            LC_JACOBIANS_CORRFACS(I,V,1:n_atmoswfs)   = IEOF * LC_IEOF(1:n_atmoswfs)
            LS_JACOBIANS_CORRFACS(I,V,1:n_surfacewfs) = IEOF * LS_IEOF(1:n_surfacewfs)
         ENDDO
         RETURN          ! RETURN
      endif

!  No Jacobians, just intensity

      DO I = 1, npoints
         IEOF = ID0 
         DO M = 1, N_EOFS
            X    = PrinComps(M,I) ; XSQ   = X * X
            IEOF = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
         ENDDO
         IEOF = EXP(IEOF)
         INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_3M_correction_csplus

subroutine pca_2M_correction_pplus &
          ( max_eofs_2p1, Max_geoms, maxlayers, max_atmoswfs, do_Jacobian, & ! Input dimensioning 
            n_eofs, npoints, ngeoms, nlayers, n_atmoswfs, Kvary, Kpars,    & ! Input control
            PrinComps, LP_PrinComps,                                       & ! Input Principal Components 
            intensity_LD_bin,    intensity_FO_bin,                         & ! Input (Bin-Intensities)
            LP_Jacobians_LD_bin, LP_Jacobians_FO_bin,                      & ! Input (Bin-Jacobians)
            intensity_corrfacs, LP_Jacobians_corrfacs )                      ! Output (Corrections)

!  Corrections with Profile (layer) Jacobians. 06 February 2012.

      IMPLICIT NONE

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms, maxlayers, max_atmoswfs

!  Control

      LOGICAL      , intent(in) :: do_Jacobian
      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms, nlayers, n_atmoswfs

!  Profile lienarization control

      LOGICAL      , intent(in) :: Kvary(maxlayers)
      INTEGER      , intent(in) :: Kpars(maxlayers)

!  Principal components, and Jacobians thereunto. Pre-allocated.

      REAL(kind=dp), intent(in) :: PrinComps    (n_eofs,npoints)
      REAL(kind=dp), intent(in) :: LP_PrinComps (n_eofs,npoints,nlayers,n_atmoswfs) 

!  Bin-Intensities 
!     (LD_BIN = LIDORT N-stream, FO_BIN = single-scattering)

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,Max_geoms)

!  Linearized Bin_intensities

      REAL(kind=dp), intent(in) :: LP_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,maxlayers,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LP_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,maxlayers,max_atmoswfs)

!  outputs; Intensity corrections and their linearizations. Pre-allocated.

      REAL(kind=dp), intent(out) :: INTENSITY_CORRFACS   (npoints,ngeoms)
      REAL(kind=dp), intent(out) :: LP_JACOBIANS_CORRFACS(npoints,ngeoms,nlayers,n_atmoswfs)

!  local variables
!  ---------------

      INTEGER      :: I, M, MM, MP, Q, K, V

      REAL(kind=dp) :: ID0, IEOF, TERM_1(N_EOFS), TERM_2(N_EOFS), X, XSQ
      REAL(kind=dp) :: IDPLUS, IDMINUS, HI_SS, LO_SS, HI_SS_MP, LO_SS_MP, HI_SS_MM, LO_SS_MM

      REAL(kind=dp) :: LP_ID0(NLAYERS,N_ATMOSWFS), LP_IEOF(NLAYERS,N_ATMOSWFS)
      REAL(kind=dp) :: LP_TERM_1(N_EOFS,NLAYERS,N_ATMOSWFS), LP_TERM_2(N_EOFS,NLAYERS,N_ATMOSWFS)
      REAL(kind=dp) :: LP_IDPLUS, LP_IDMINUS, LP_HI_SS, LP_LO_SS, LP_X, LP_XSQ

      REAL(kind=dp), parameter ::  HALF = 0.5_dp, two = 2.0_dp

!  Zero output
!  -----------

!    PLACEHOLDER - is it really necessary?

!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  mean value
!  ----------

!  Mean-value calculation

      HI_SS = INTENSITY_LD_BIN(1,V)
      LO_SS = INTENSITY_FO_BIN(1,V)
      ID0   = LOG ( HI_SS / LO_SS )

!  Linearized mean value calculation

      if ( do_Jacobian ) then
         LP_ID0 = 0.0_dp    ! Safety-first zeroing
         do k = 1, nlayers
            if ( kvary(k) ) then
               do Q = 1, kpars(k)
                  LP_HI_SS    = LP_JACOBIANS_LD_BIN(1,V,K,Q)
                  LP_LO_SS    = LP_JACOBIANS_FO_BIN(1,V,K,Q)
                  LP_ID0(K,Q) = ( LP_HI_SS / HI_SS ) - ( LP_LO_SS / LO_SS )
               enddo
            endif
         enddo
      endif

!  Plus/minus values
!  -----------------

      DO M = 1, N_EOFS

!  Indices

         MP = 2*M ; MM = MP + 1

!  Plus/Minus values. Finite Differences TERM_1 and TERM_2

         HI_SS_MP   = INTENSITY_LD_BIN(MP,V)
         LO_SS_MP   = INTENSITY_FO_BIN(MP,V)
         IDPLUS     = LOG ( HI_SS_MP / LO_SS_MP )
         HI_SS_MM   = INTENSITY_LD_BIN(MM,V)
         LO_SS_MM   = INTENSITY_FO_BIN(MM,V)
         IDMINUS    = LOG ( HI_SS_MM / LO_SS_MM )
         TERM_1(M)  = ( IDPLUS - IDMINUS )            * HALF
         TERM_2(M)  = ( IDPLUS + IDMINUS - TWO*ID0 )  * HALF

!  Linearizations

         if ( do_Jacobian ) then
            do k = 1, nlayers
               if ( kvary(k) ) then
                  do Q = 1, kpars(k)
                     LP_HI_SS   = LP_JACOBIANS_LD_BIN(MP,V,k,Q)
                     LP_LO_SS   = LP_JACOBIANS_FO_BIN(MP,V,k,Q)
                     LP_IDPLUS  = ( LP_HI_SS / HI_SS_MP ) - ( LP_LO_SS / LO_SS_MP )
                     LP_HI_SS   = LP_JACOBIANS_LD_BIN(MM,V,k,Q)
                     LP_LO_SS   = LP_JACOBIANS_FO_BIN(MM,V,k,Q)
                     LP_IDMINUS = ( LP_HI_SS / HI_SS_MM ) - ( LP_LO_SS / LO_SS_MM )
                     LP_TERM_1(M,k,Q)  = ( LP_IDPLUS - LP_IDMINUS )                   * HALF
                     LP_TERM_2(M,k,Q)  = ( LP_IDPLUS + LP_IDMINUS - TWO*LP_ID0(k,Q) ) * HALF
                  enddo
               endif
            enddo
         endif

!  End EOF loop

      ENDDO

!  set correction for all wavelengths
!  ----------------------------------

!  Including Jacobians. Return when done.

      if ( do_Jacobian ) then
         DO I = 1, npoints
            IEOF = ID0 ; LP_IEOF = LP_ID0
            DO M = 1, N_EOFS
               X   = PrinComps(M,I) ; XSQ   = X * X
               IEOF   = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
               do k = 1, nlayers
                  if ( kvary(k) ) then
                     do Q = 1, kpars(k)
                        LP_X = LP_PrinComps(M,I,K,Q) ; LP_XSQ = TWO * X * LP_X
                        LP_IEOF(k,Q) = LP_IEOF(k,Q) + LP_TERM_1(M,k,Q) * X   + TERM_1(M) * LP_X   &
                                                    + LP_TERM_2(M,k,Q) * XSQ + TERM_2(M) * LP_XSQ
                     enddo
                  endif
               enddo
            ENDDO
            IEOF = EXP(IEOF)
            INTENSITY_CORRFACS(I,V) = IEOF
            do k = 1, nlayers
               LP_JACOBIANS_CORRFACS(I,V,k,1:n_atmoswfs) = IEOF * LP_IEOF(k,1:n_atmoswfs)
            enddo
         ENDDO
         RETURN          ! RETURN
      endif

!  No Jacobians, just intensity

      DO I = 1, npoints
         IEOF = ID0 
         DO M = 1, N_EOFS
            X    = PrinComps(M,I) ; XSQ   = X * X
            IEOF = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
         ENDDO
         IEOF = EXP(IEOF)
         INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_2M_correction_pplus


subroutine pca_2M_correction_cplus &
          ( max_eofs_2p1, Max_geoms, max_atmoswfs,             & ! Input dimensioning 
            do_Jacobian, n_eofs, npoints, ngeoms, n_atmoswfs,  & ! Input control
            PrinComps, LC_PrinComps,                           & ! Input Principal Components 
            intensity_LD_bin,    intensity_FO_bin,             & ! Input (Bin-Intensities)
            LC_Jacobians_LD_bin, LC_Jacobians_FO_bin,          & ! Input (Linearized
            intensity_corrfacs, LC_Jacobians_corrfacs )          ! Output (Corrections)

!  Column (bulk) Jacobians. 03 February 2012.

      IMPLICIT NONE

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms, max_atmoswfs

!  Control

      LOGICAL      , intent(in) :: do_Jacobian
      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms, n_atmoswfs

!  Principal components, and Jacobians thereunto. Pre-allocated

      REAL(kind=dp), intent(in) :: PrinComps    (n_eofs,npoints)
      REAL(kind=dp), intent(in) :: LC_PrinComps (n_eofs,npoints,n_atmoswfs) 

!  Bin-Intensities 
!     (LD_BIN = LIDORT N-stream, FO_BIN = single-scattering)

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,Max_geoms)

!  Linearized Bin_intensities

      REAL(kind=dp), intent(in) :: LC_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LC_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)

!  outputs; Intensity corrections and their linearizations. Pre-allocated.

      REAL(kind=dp), intent(out) :: INTENSITY_CORRFACS   (npoints,ngeoms)
      REAL(kind=dp), intent(out) :: LC_JACOBIANS_CORRFACS(npoints,ngeoms,n_atmoswfs)

!  local variables
!  ---------------

      INTEGER      :: I, M, MM, MP, Q, V

      REAL(kind=dp) :: ID0, IEOF, TERM_1(N_EOFS), TERM_2(N_EOFS),  X, XSQ
      REAL(kind=dp) :: IDPLUS, IDMINUS, HI_SS, LO_SS, HI_SS_MP, LO_SS_MP, HI_SS_MM, LO_SS_MM

      REAL(kind=dp) :: LC_ID0(N_ATMOSWFS), LC_IEOF(N_ATMOSWFS)
      REAL(kind=dp) :: LC_TERM_1(N_EOFS,N_ATMOSWFS), LC_TERM_2(N_EOFS,N_ATMOSWFS)
      REAL(kind=dp) :: LC_IDPLUS, LC_IDMINUS, LC_HI_SS, LC_LO_SS, LC_X, LC_XSQ

      REAL(kind=dp), parameter ::  HALF = 0.5_dp, two = 2.0_dp

!  Zero output
!  -----------

!    PLACEHOLDER - is it really necessary?

!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  mean value
!  ----------

!  mean value calculation

      HI_SS = INTENSITY_LD_BIN(1,V)
      LO_SS = INTENSITY_FO_BIN(1,V)
      ID0 = DLOG ( HI_SS / LO_SS )

!  Linearized mean value calculation

      if ( do_Jacobian ) then
         do Q = 1, n_atmoswfs
            LC_HI_SS  = LC_JACOBIANS_LD_BIN(1,V,Q)
            LC_LO_SS  = LC_JACOBIANS_FO_BIN(1,V,Q)
            LC_ID0(Q) = ( LC_HI_SS / HI_SS ) - ( LC_LO_SS / LO_SS )
         enddo
      endif

!  Plus/minus values
!  -----------------

      DO M = 1, N_EOFS

!  Indices

        MP = 2*M ; MM = MP + 1

!  Plus/Minus values. Finite Differences TERM_1 and TERM_2

        HI_SS_MP   = INTENSITY_LD_BIN(MP,V)
        LO_SS_MP   = INTENSITY_FO_BIN(MP,V)
        IDPLUS     = LOG ( HI_SS_MP / LO_SS_MP )

        HI_SS_MM   = INTENSITY_LD_BIN(MM,V)
        LO_SS_MM   = INTENSITY_FO_BIN(MM,V)
        IDMINUS    = LOG ( HI_SS_MM / LO_SS_MM )

        TERM_1(M) = ( IDPLUS - IDMINUS )             * HALF
        TERM_2(M) = ( IDPLUS + IDMINUS - TWO * ID0 ) * HALF

!  Linearizations

        if ( do_Jacobian ) then
           do Q = 1, n_atmoswfs
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MP,V,Q)
              LC_LO_SS   = LC_JACOBIANS_FO_BIN(MP,V,Q)
              LC_IDPLUS  = ( LC_HI_SS / HI_SS_MP ) - ( LC_LO_SS / LO_SS_MP )
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MM,V,Q)
              LC_LO_SS   = LC_JACOBIANS_FO_BIN(MM,V,Q)
              LC_IDMINUS = ( LC_HI_SS / HI_SS_MM ) - ( LC_LO_SS / LO_SS_MM )
              LC_TERM_1(M,Q)  = ( LC_IDPLUS - LC_IDMINUS )                * HALF
              LC_TERM_2(M,Q)  = ( LC_IDPLUS + LC_IDMINUS - TWO*LC_ID0(Q) ) * HALF
           enddo
        endif

!  End EOF loop

      ENDDO

!  set correction for all wavelengths
!  ----------------------------------

!  Including Jacobian. Return when done.

      if ( do_Jacobian ) then
         DO I = 1, npoints
            IEOF = ID0 ; LC_IEOF = LC_ID0
            DO M = 1, N_EOFS
               X   = PrinComps(M,I)   ; XSQ   = X * X
               IEOF   = IEOF   + TERM_1(M)   * X + TERM_2(M) * XSQ
               do q = 1, n_atmoswfs
                  LC_X = LC_PrinComps(M,I,Q) ; LC_XSQ = two * X * LC_X
                  LC_IEOF(Q) = LC_IEOF(Q) + LC_TERM_1(M,Q) * X   + TERM_1(M) * LC_X   &
                                          + LC_TERM_2(M,Q) * XSQ + TERM_2(M) * LC_XSQ
               enddo
            ENDDO
            IEOF = EXP(IEOF)
            INTENSITY_CORRFACS(I,V)                 = IEOF
            LC_JACOBIANS_CORRFACS(I,V,1:n_atmoswfs) = IEOF * LC_IEOF(1:n_atmoswfs)
         ENDDO
         RETURN          ! RETURN
      endif

!  No Jacobians, just intensity

      DO I = 1, npoints
         IEOF = ID0 
         DO M = 1, N_EOFS
            X    = PrinComps(M,I) ; XSQ   = X * X
            IEOF = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
         ENDDO
         IEOF = EXP(IEOF)
         INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_2M_correction_cplus

subroutine pca_2M_correction_csplus  &
          ( max_eofs_2p1, Max_geoms, max_atmoswfs, max_surfacewfs, do_Jacobian, & ! Input dimensioning 
            do_SJacobian, n_eofs, npoints, ngeoms, n_atmoswfs, n_surfacewfs,    & ! Input control
            PrinComps, LC_PrinComps, LS_PrinComps,                              & ! Input Principal Components 
            intensity_LD_bin,    intensity_FO_bin,                              & ! Input (Bin-Intensities)
            LC_Jacobians_LD_bin, LC_Jacobians_FO_bin,                           & ! Input (Bin-Jacobians)
            LS_Jacobians_LD_bin, LS_Jacobians_FO_bin,                           & ! Input (Bin-Jacobians)
            intensity_corrfacs, LC_Jacobians_corrfacs, LS_Jacobians_corrfacs )    ! Output (Corrections)

!  Column (bulk) Jacobians. 03 February 2012.

      IMPLICIT NONE

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms, max_atmoswfs, max_surfacewfs

!  Control

      LOGICAL      , intent(in) :: do_Jacobian, do_SJacobian
      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms, n_atmoswfs, n_surfacewfs

!  Principal components, and Jacobians thereunto. Pre-allocated.

      REAL(kind=dp), intent(in) :: PrinComps    (n_eofs,npoints)
      REAL(kind=dp), intent(in) :: LC_PrinComps (n_eofs,npoints,n_atmoswfs) 
      REAL(kind=dp), intent(in) :: LS_PrinComps (n_eofs,npoints,n_surfacewfs) 

!  Bin-Intensities 
!     (LD_BIN = LIDORT N-stream, FO_BIN = single-scattering)

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,Max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,Max_geoms)

!  Linearized Bin_intensities

      REAL(kind=dp), intent(in) :: LC_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)
      REAL(kind=dp), intent(in) :: LC_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,max_atmoswfs)

      REAL(kind=dp), intent(in) :: LS_JACOBIANS_LD_BIN(max_eofs_2p1,Max_geoms,max_surfacewfs)
      REAL(kind=dp), intent(in) :: LS_JACOBIANS_FO_BIN(max_eofs_2p1,Max_geoms,max_surfacewfs)

!  outputs; Intensity corrections and their linearizations. Pre-allocated.

      REAL(kind=dp), intent(out) :: INTENSITY_CORRFACS   (npoints,ngeoms)
      REAL(kind=dp), intent(out) :: LC_JACOBIANS_CORRFACS(npoints,ngeoms,n_atmoswfs)
      REAL(kind=dp), intent(out) :: LS_JACOBIANS_CORRFACS(npoints,ngeoms,n_surfacewfs)

!  local variables
!  ---------------

      INTEGER      :: I, M, MM, MP, Q, V

      REAL(kind=dp) :: ID0, IEOF, TERM_1(N_EOFS), TERM_2(N_EOFS),  X, XSQ
      REAL(kind=dp) :: IDPLUS, IDMINUS, HI_SS, LO_SS, HI_SS_MP, LO_SS_MP, HI_SS_MM, LO_SS_MM

      REAL(kind=dp) :: LC_ID0(N_ATMOSWFS), LC_IEOF(N_ATMOSWFS)
      REAL(kind=dp) :: LC_TERM_1(N_EOFS,N_ATMOSWFS), LC_TERM_2(N_EOFS,N_ATMOSWFS)
      REAL(kind=dp) :: LC_IDPLUS, LC_IDMINUS, LC_HI_SS, LC_LO_SS, LC_X, LC_XSQ

      REAL(kind=dp) :: LS_ID0(N_SURFACEWFS), LS_IEOF(N_SURFACEWFS)
      REAL(kind=dp) :: LS_TERM_1(N_EOFS,N_SURFACEWFS), LS_TERM_2(N_EOFS,N_SURFACEWFS)
      REAL(kind=dp) :: LS_IDPLUS, LS_IDMINUS, LS_HI_SS, LS_LO_SS, LS_X, LS_XSQ

      REAL(kind=dp), parameter ::  HALF = 0.5_dp, two = 2.0_dp

!  Zero output
!  -----------

!    PLACEHOLDER - is it really necessary?

!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  mean value
!  ----------

!  mean value calculation

      HI_SS = INTENSITY_LD_BIN(1,V)
      LO_SS = INTENSITY_FO_BIN(1,V)
      ID0   = LOG ( HI_SS / LO_SS )

!  Linearized mean value calculation

      if ( do_Jacobian ) then
         do Q = 1, n_atmoswfs
            LC_HI_SS  = LC_JACOBIANS_LD_BIN(1,V,Q)
            LC_LO_SS  = LC_JACOBIANS_FO_BIN(1,V,Q)
            LC_ID0(Q) = ( LC_HI_SS / HI_SS ) - ( LC_LO_SS / LO_SS )
         enddo
      endif

      if ( do_SJacobian ) then
         do Q = 1, n_surfacewfs
            LS_HI_SS  = LS_JACOBIANS_LD_BIN(1,V,Q)
            LS_LO_SS  = LS_JACOBIANS_FO_BIN(1,V,Q)
            LS_ID0(Q) = ( LS_HI_SS / HI_SS ) - ( LS_LO_SS / LO_SS )
         enddo
      endif

!  Plus/minus values
!  -----------------

      DO M = 1, N_EOFS

!  Indices

        MP = 2*M ; MM = MP + 1

!  Plus/Minus values. Finite Differences TERM_1 and TERM_2

        HI_SS_MP   = INTENSITY_LD_BIN(MP,V)
        LO_SS_MP   = INTENSITY_FO_BIN(MP,V)
        IDPLUS     = LOG ( HI_SS_MP / LO_SS_MP )

        HI_SS_MM   = INTENSITY_LD_BIN(MM,V)
        LO_SS_MM   = INTENSITY_FO_BIN(MM,V)
        IDMINUS    = LOG ( HI_SS_MM / LO_SS_MM )

        TERM_1(M) = ( IDPLUS - IDMINUS )           * HALF
        TERM_2(M) = ( IDPLUS + IDMINUS - TWO*ID0 ) * HALF

!  Column Linearizations

        if ( do_Jacobian ) then
           do Q = 1, n_atmoswfs
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MP,V,Q)
              LC_LO_SS   = LC_JACOBIANS_FO_BIN(MP,V,Q)
              LC_IDPLUS  = ( LC_HI_SS / HI_SS_MP ) - ( LC_LO_SS / LO_SS_MP )
              LC_HI_SS   = LC_JACOBIANS_LD_BIN(MM,V,Q)
              LC_LO_SS   = LC_JACOBIANS_FO_BIN(MM,V,Q)
              LC_IDMINUS = ( LC_HI_SS / HI_SS_MM ) - ( LC_LO_SS / LO_SS_MM )
              LC_TERM_1(M,Q)  = ( LC_IDPLUS - LC_IDMINUS )                 * HALF
              LC_TERM_2(M,Q)  = ( LC_IDPLUS + LC_IDMINUS - TWO*LC_ID0(Q) ) * HALF
           enddo
        endif

!  surface Linearizations

        if ( do_SJacobian ) then
           do Q = 1, n_surfacewfs
              LS_HI_SS   = LS_JACOBIANS_LD_BIN(MP,V,Q)
              LS_LO_SS   = LS_JACOBIANS_FO_BIN(MP,V,Q)
              LS_IDPLUS  = ( LS_HI_SS / HI_SS_MP ) - ( LS_LO_SS / LO_SS_MP )
              LS_HI_SS   = LS_JACOBIANS_LD_BIN(MM,V,Q)
              LS_LO_SS   = LS_JACOBIANS_FO_BIN(MM,V,Q)
              LS_IDMINUS = ( LS_HI_SS / HI_SS_MM ) - ( LS_LO_SS / LO_SS_MM )
              LS_TERM_1(M,Q)  = ( LS_IDPLUS - LS_IDMINUS )                 * HALF
              LS_TERM_2(M,Q)  = ( LS_IDPLUS + LS_IDMINUS - TWO*LS_ID0(Q) ) * HALF
           enddo
        endif

!  End EOF loop

      ENDDO

!  set correction for all wavelengths
!  ----------------------------------

!  Including Jacobians (BOTH). Return when done.

      if ( do_Jacobian .and. do_SJacobian  ) then
         DO I = 1, npoints
            IEOF = ID0
            LC_IEOF(1:n_atmoswfs)   = LC_ID0(1:n_atmoswfs)
            LS_IEOF(1:n_surfacewfs) = LS_ID0(1:n_surfacewfs)
            DO M = 1, N_EOFS
               X   = PrinComps(M,I)   ; XSQ   = X * X
               IEOF   = IEOF   + TERM_1(M)   * X + TERM_2(M) * XSQ
               do q = 1, n_atmoswfs
                  LC_X = LC_PrinComps(M,I,Q) ; LC_XSQ = TWO * X * LC_X
                  LC_IEOF(Q) = LC_IEOF(Q) + LC_TERM_1(M,Q) * X   + TERM_1(M) * LC_X   &
                                          + LC_TERM_2(M,Q) * XSQ + TERM_2(M) * LC_XSQ
               enddo
               do q = 1, n_surfacewfs
                  LS_X = LS_PrinComps(M,I,Q) ; LS_XSQ = TWO * X * LS_X
                  LS_IEOF(Q) = LS_IEOF(Q) + LS_TERM_1(M,Q) * X   + TERM_1(M) * LS_X   &
                                          + LS_TERM_2(M,Q) * XSQ + TERM_2(M) * LS_XSQ
               enddo
            ENDDO
            IEOF = EXP(IEOF)
            INTENSITY_CORRFACS(I,V)                   = IEOF
            LC_JACOBIANS_CORRFACS(I,V,1:n_atmoswfs)   = IEOF * LC_IEOF(1:n_atmoswfs)
            LS_JACOBIANS_CORRFACS(I,V,1:n_surfacewfs) = IEOF * LS_IEOF(1:n_surfacewfs)
         ENDDO
         RETURN          ! RETURN
      endif

!  No Jacobians, just intensity

      DO I = 1, npoints
         IEOF = ID0 
         DO M = 1, N_EOFS
            X    = PrinComps(M,I) ; XSQ   = X * X
            IEOF = IEOF + TERM_1(M) * X + TERM_2(M) * XSQ
         ENDDO
         IEOF = EXP(IEOF)
         INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_2M_correction_csplus

!  End module

End module pca_correction_plus_m

