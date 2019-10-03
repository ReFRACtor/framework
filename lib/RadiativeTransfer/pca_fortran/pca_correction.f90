module pca_correction_m

!  Derivation of Radiance corrections at all points 
!     using principal-component projections and PCA-binned intensity calculations

!  Version 1. Original code.       19 November 2011
!  Version 2. Robustness, 2M/3M.   08 November 2012
!  Version 3. Multigeometry.       05 April    2013

!  Subroutines:

!          pca_3M_correction     --> Correction with LD/2S/FO (3 Intensities)
!          pca_2M_correction     --> Correction with LD/2S    (2 Intensities)

!          pca_3M_correction_3OD --> Correction with LD/2S/FO (3 Intensities), Third-order Differencing

!  Remark:
!   Max_eofs_2p1 = 2 * Max_eofs + 1
!   Max_eofs_4p1 = 4 * Max_eofs + 1      Should be the input for the 3OD routine

use iso_c_binding

public

contains

subroutine pca_3M_correction &
          ( max_eofs_2p1, Max_geoms, n_eofs, npoints, ngeoms, PrinComps, & ! Inputs (control, PCs)
            intensity_LD_bin, intensity_2S_bin, intensity_FO_bin,        & ! Inputs (Bin-Intensities)
            Intensity_Corrfacs ) bind(C)                                   ! Output (Corrections)

      IMPLICIT NONE

!  precision parameters

!      INTEGER, PARAMETER :: sp     = selected_real_kind(6)
      INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensioning

      INTEGER(c_int)      , intent(in) :: max_eofs_2p1, Max_geoms

!  Control

      INTEGER(c_int)      , intent(in) :: n_eofs, npoints, ngeoms

!  Principal components (Pre-allocated)

      REAL(kind=c_double), intent(in) :: PrinComps(n_eofs,npoints)

!  Bin-Intensities 

      REAL(kind=c_double), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,max_geoms)
      REAL(kind=c_double), intent(in) :: INTENSITY_2S_BIN(max_eofs_2p1,max_geoms)
      REAL(kind=c_double), intent(in) :: INTENSITY_FO_BIN(max_eofs_2p1,max_geoms)

!  Correction factors (Pre-allocated)

      REAL(kind=c_double), intent(out) :: Intensity_Corrfacs(npoints,ngeoms)

!  local variables

      INTEGER       :: I, M, M2, M21, V
      REAL(kind=dp) :: ID0, IDPLUS(n_eofs), IDMINUS(n_eofs), IEOF, PC, PCSQ
      REAL(kind=dp), parameter :: half = 0.5_dp, two = 2.0_dp
 
!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  Mean value ratios

      ID0 = LOG ( ( INTENSITY_LD_BIN(1,V) + INTENSITY_FO_BIN(1,V) ) / &
                  ( INTENSITY_2S_BIN(1,V) + INTENSITY_FO_BIN(1,V) ) )

!  +/- EOF ratios

      DO M = 1, N_EOFS
        M2 = 2 * M ; M21 = M2 + 1
        IDPLUS(M)  = LOG ( ( INTENSITY_LD_BIN(M2,V)  + INTENSITY_FO_BIN(M2,V)  ) / &
                           ( INTENSITY_2S_BIN(M2,V)  + INTENSITY_FO_BIN(M2,V)  ) )
        IDMINUS(M) = LOG ( ( INTENSITY_LD_BIN(M21,V) + INTENSITY_FO_BIN(M21,V) ) / &
                           ( INTENSITY_2S_BIN(M21,V) + INTENSITY_FO_BIN(M21,V) ) )
      ENDDO

!  Generate corrections

      DO I = 1, NPOINTS
        IEOF = ID0
        DO M = 1, N_EOFS
          PC = PrinComps(M,I) ; PCSQ = PC * PC
          IEOF = IEOF + half * PC   * ( IDPLUS(M) - IDMINUS(M) )            &
                      + half * PCSQ * ( IDPLUS(M) - two*ID0 + IDMINUS(M) )
        ENDDO
        IEOF = EXP(IEOF)
        INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_3M_correction

!

subroutine pca_2M_correction &
          ( max_eofs_2p1, Max_geoms, n_eofs, npoints, ngeoms, PrinComps, & ! Inputs (control, PCs)
            intensity_LD_bin, intensity_2S_bin,                          & ! Inputs (Bin-Intensities)
            Intensity_Corrfacs )                                           ! Output (Corrections)

      IMPLICIT NONE

!  precision parameters

!      INTEGER, PARAMETER :: sp     = selected_real_kind(6)
      INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_2p1, Max_geoms

!  Control

      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms

!  Principal components (Pre-allocated)

      REAL(kind=dp), intent(in) :: PrinComps(n_eofs,npoints)

!  Bin-Intensities 

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_2p1,max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_2S_BIN(max_eofs_2p1,max_geoms)

!  Correction factors (Pre-allocated)

      REAL(kind=dp), intent(out) :: Intensity_Corrfacs(npoints,ngeoms)

!  local variables

      INTEGER       :: I, M, M2, M21, V
      REAL(kind=dp) :: ID0, IDPLUS(n_eofs), IDMINUS(n_eofs), IEOF, PC, PCSQ
      REAL(kind=dp), parameter :: half = 0.5_dp, two = 2.0_dp
 
!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  Mean value ratios

!write(0,*) intensity_2s_bin(1,v),intensity_ld_bin(1,v)
      ID0 = LOG ( INTENSITY_LD_BIN(1,V) / INTENSITY_2S_BIN(1,V) )

!  +/- EOF ratios

      DO M = 1, N_EOFS
        M2 = 2 * M ; M21 = M2 + 1
        IDPLUS(M)  = LOG ( INTENSITY_LD_BIN(M2,V)  / INTENSITY_2S_BIN(M2,V) )
        IDMINUS(M) = LOG ( INTENSITY_LD_BIN(M21,V) / INTENSITY_2S_BIN(M21,V) )
      ENDDO

!  Generate corrections

      DO I = 1, NPOINTS
        IEOF = ID0
        DO M = 1, N_EOFS
          PC = PrinComps(M,I) ; PCSQ = PC * PC
          IEOF = IEOF + half * PC   * ( IDPLUS(M) - IDMINUS(M) )            &
                      + half * PCSQ * ( IDPLUS(M) - two*ID0 + IDMINUS(M) )
        ENDDO
        IEOF = EXP(IEOF)
        INTENSITY_CORRFACS(I,V) = IEOF
      ENDDO

!  End geometry

   ENDDO

   return
end subroutine pca_2M_correction

!

subroutine pca_3M_correction_3OD &
          ( max_eofs_4p1, Max_geoms, n_eofs, npoints, ngeoms, PrinComps, & ! Inputs (control, PCs)
            intensity_LD_bin, intensity_2S_bin, intensity_FO_bin,        & ! Inputs (Bin-Intensities)
            Intensity_Corrfacs )                                           ! Output (Corrections)

      IMPLICIT NONE

!  precision parameters

!      INTEGER, PARAMETER :: sp     = selected_real_kind(6)
      INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensioning

      INTEGER      , intent(in) :: max_eofs_4p1, Max_geoms

!  Control

      INTEGER      , intent(in) :: n_eofs, npoints, ngeoms

!  Principal components (Pre-allocated)

      REAL(kind=dp), intent(in) :: PrinComps(n_eofs,npoints)

!  Bin-Intensities 

      REAL(kind=dp), intent(in) :: INTENSITY_LD_BIN(max_eofs_4p1,max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_2S_BIN(max_eofs_4p1,max_geoms)
      REAL(kind=dp), intent(in) :: INTENSITY_FO_BIN(max_eofs_4p1,max_geoms)

!  Correction factors (Pre-allocated)

      REAL(kind=dp), intent(out) :: Intensity_Corrfacs(npoints,ngeoms)

!  local variables

      INTEGER       :: I, M, M1, M2, M3, M4, V
      REAL(kind=dp) :: IDPLUS2(n_eofs), IDMINUS2(n_eofs), ID1_DIFF, ID2_DIFF, ID3_DIFF
      REAL(kind=dp) :: ID0, IDPLUS(n_eofs), IDMINUS(n_eofs), IEOF, H1, H2, H3, PC
      REAL(kind=dp), parameter :: half = 0.5_dp, two = 2.0_dp
 
!  Start Geometry loop
!  -------------------

   DO V = 1, ngeoms

!  Mean

      ID0 = LOG( (INTENSITY_LD_BIN(1,v)+INTENSITY_FO_BIN(1,v)) / &
                 (INTENSITY_2S_BIN(1,v)+INTENSITY_FO_BIN(1,v)))

!  +/- EOF, +/- 2*EOF ratios

      DO M = 1, N_EOFS
         M1 = 4*M-2 ; M2 = M1 + 1 ; M3 = M2 + 1 ; M4 = M3 + 1
 
         IDPLUS(M)   = LOG(  ( INTENSITY_LD_BIN(M2,v) + INTENSITY_FO_BIN(M2,v) ) / &
                             ( INTENSITY_2S_BIN(M2,v) + INTENSITY_FO_BIN(M2,v) )  )

         IDMINUS(M)  = LOG(  ( INTENSITY_LD_BIN(M3,v) + INTENSITY_FO_BIN(M3,v) ) / &
                             ( INTENSITY_2S_BIN(M3,v) + INTENSITY_FO_BIN(M3,v) )  )

         IDPLUS2(M)  = LOG(  ( INTENSITY_LD_BIN(M1,v) + INTENSITY_FO_BIN(M1,v) ) / &
                             ( INTENSITY_2S_BIN(M1,v) + INTENSITY_FO_BIN(M1,v) )  )

         IDMINUS2(M) = LOG(  ( INTENSITY_LD_BIN(M4,v) + INTENSITY_FO_BIN(M4,v) ) / &
                             ( INTENSITY_2S_BIN(M4,v) + INTENSITY_FO_BIN(M4,v) )  )

      ENDDO

!  Get corrections

      DO I = 1, NPOINTS
         IEOF = ID0
         DO M = 1, N_EOFS
            PC = PrinComps(M,I) ; H1 = half * PC ; H2 = H1 *PC ; H3 = H2 * PC / 6.0d0
            ID1_DIFF = IDPLUS(M)  - IDMINUS(M)
            ID2_DIFF = IDPLUS(M)  + IDMINUS(M)  - two*ID0
            ID3_DIFF = IDPLUS2(M) - IDMINUS2(M) - two*ID1_DIFF 
            IEOF = IEOF + ID1_DIFF*H1 + ID2_DIFF*H2 + ID3_DIFF*H3
         ENDDO
         INTENSITY_CORRFACS(I,V) = EXP(IEOF)
      ENDDO

!  End geometry

   ENDDO

!  finish

   return
end subroutine pca_3M_correction_3OD

!  finish module

end module pca_correction_m

