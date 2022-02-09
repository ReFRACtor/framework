module raman_sioris

use iso_c_binding

! =========
! Constants
! =========

INTEGER, PARAMETER :: N2Jmax=28, O2maxJ=53, O2max=94
INTEGER, PARAMETER :: dp = KIND(1.0D0) 

REAL (KIND=dp), DIMENSION(0:N2Jmax)      :: N2E = (/ &
    0.0000_dp, 3.9791_dp, 11.9373_dp, 23.8741_dp, 39.7892_dp, 59.6821_dp, &
    83.5521_dp, 111.3983_dp, 143.2197_dp, 179.0154_dp, 218.7839_dp, 262.524_dp, &
    310.2341_dp, 361.9126_dp, 417.5576_dp, 477.1673_dp, 540.7395_dp, 608.2722_dp, &
    679.7628_dp, 755.209_dp, 834.6081_dp, 917.9574_dp, 1005.254_dp, 1096.4948_dp, &
    1191.6766_dp, 1290.7963_dp, 1393.8503_dp, 1500.835_dp, 1611.7467_dp &
/)

REAL (KIND=dp), DIMENSION(0:2*N2Jmax-3)  :: N2b = (/ &
    0.2000_dp, 0.2571_dp, 0.2857_dp, 0.3030_dp, 0.3147_dp, 0.3231_dp, &
    0.3294_dp,  0.3344_dp, 0.3383_dp, 0.3416_dp, 0.3443_dp, 0.3467_dp, &
    0.3487_dp, 0.3504_dp, 0.3519_dp, 0.3532_dp, 0.3544_dp, 0.3555_dp, &
    0.3565_dp, 0.3573_dp, 0.3581_dp, 0.3589_dp, 0.3595_dp, 0.3601_dp, &
    0.3607_dp, 0.3612_dp, 0.3617_dp, 1.0000_dp, 0.6000_dp, 0.5143_dp, &
    0.4762_dp, 0.4545_dp, 0.4406_dp, 0.4308_dp, 0.4235_dp, 0.4180_dp, &
    0.4135_dp, 0.4099_dp, 0.4070_dp, 0.4044_dp, 0.4023_dp, 0.4004_dp, &
    0.3988_dp, 0.3974_dp, 0.3961_dp, 0.3950_dp, 0.3940_dp, 0.3931_dp, &
    0.3922_dp, 0.3915_dp, 0.3908_dp, 0.3902_dp, 0.3896_dp, 0.3890_dp  &
/)

REAL (KIND=dp), DIMENSION(0:O2maxJ)      :: O2EnZ = (/ &
    0.0000_dp, 2.0843_dp, 3.9611_dp, 16.2529_dp, 16.3876_dp, 18.3372_dp, &
    42.2001_dp, 42.224_dp, 44.2117_dp, 79.5646_dp, 79.607_dp, 81.5805_dp, &
    128.3978_dp, 128.4921_dp, 130.4376_dp, 188.7135_dp, 188.8532_dp, &
    190.7749_dp, 260.5011_dp, 260.6826_dp, 262.5829_dp, 343.7484_dp, &
    343.9697_dp, 345.85_dp, 438.4418_dp, 438.7015_dp, 440.562_dp, &
    544.5658_dp, 544.8628_dp, 546.705_dp, 662.103_dp, 662.4368_dp, &
    664.261_dp, 791.0344_dp, 791.4045_dp, 793.21_dp, 931.339_dp, 931.745_dp, &
    933.533_dp, 1082.9941_dp, 1083.4356_dp, 1085.206_dp, 1245.975_dp, 1246.4518_dp, &
    1248.204_dp, 1420.2552_dp, 1420.7672_dp, 1422.502_dp, 1605.8064_dp, 1606.3533_dp, &
    1608.071_dp, 1802.5983_dp, 1803.1802_dp, 1804.881_dp &
/)

REAL (KIND=dp), DIMENSION(0:2*O2max-7)   :: O2E = (/ &
    16.3876_dp, 546.705_dp, 440.562_dp, 3.9611_dp, 345.85_dp, 262.5829_dp, &
    190.7749_dp, 130.4376_dp, 18.3372_dp, 81.5805_dp, 44.2117_dp, 44.2117_dp, &
    81.5805_dp, 130.4376_dp, 190.7749_dp, 262.5829_dp, 18.3372_dp, 2.0843_dp, &
    345.85_dp, 440.562_dp, 546.705_dp, 16.2529_dp, 16.2529_dp, 16.3876_dp, &
    18.3372_dp, 16.2529_dp, 18.3372_dp, 42.2001_dp, 42.2001_dp, 42.224_dp, &
    44.2117_dp, 42.2001_dp, 44.2117_dp, 79.607_dp, 79.5646_dp, 81.5805_dp, &
    79.607_dp, 81.5805_dp, 42.2001_dp, 128.4921_dp, 128.3978_dp, 130.4376_dp, &
    128.4921_dp, 130.4376_dp, 188.8532_dp, 188.7135_dp, 190.7749_dp, &
    188.8532_dp, 190.7749_dp, 260.6826_dp, 260.5011_dp, 262.5829_dp, &
    260.6826_dp, 262.5829_dp, 343.9697_dp, 343.7484_dp, 345.85_dp, 343.9697_dp, &
    345.85_dp, 438.7015_dp, 438.4418_dp, 440.562_dp, 438.7015_dp, 440.562_dp, &
    544.8628_dp, 544.5658_dp, 546.705_dp, 544.8628_dp, 546.705_dp, 662.103_dp, &
    664.261_dp, 662.4368_dp, 791.0344_dp, 793.21_dp, 791.4045_dp, 931.339_dp, &
    933.533_dp, 931.745_dp, 1082.9941_dp, 1085.206_dp, 1083.4356_dp, 1245.975_dp, &
    1248.204_dp, 1246.4518_dp, 1420.2552_dp, 1422.502_dp, 1420.7672_dp, &
    1605.8064_dp, 1608.071_dp, 1606.3533_dp, 16.2529_dp, 544.8628_dp, &
    438.7015_dp, 2.0843_dp, 343.9697_dp, 260.6826_dp, 188.8532_dp, 128.4921_dp, &
    16.3876_dp, 79.607_dp, 42.224_dp, 42.2001_dp, 79.5646_dp, 128.3978_dp, &
    188.7135_dp, 260.5011_dp, 16.2529_dp, 0.0000_dp, 343.7484_dp, 438.4418_dp, &
    544.5658_dp, 3.9611_dp, 2.0843_dp, 2.0843_dp, 3.9611_dp, 0.0000_dp, &
    2.0843_dp, 18.3372_dp, 16.3876_dp, 16.3876_dp, 18.3372_dp, 16.2529_dp, &
    16.3876_dp, 44.2117_dp, 42.224_dp, 44.2117_dp, 42.2001_dp, 42.224_dp, &
    2.0843_dp, 81.5805_dp, 79.5646_dp, 81.5805_dp, 79.607_dp, 79.5646_dp, &
    130.4376_dp, 128.3978_dp, 130.4376_dp, 128.4921_dp, 128.3978_dp, 190.7749_dp, &
    188.7135_dp, 190.7749_dp, 188.8532_dp, 188.7135_dp, 262.5829_dp, 260.5011_dp, &
    262.5829_dp, 260.6826_dp, 260.5011_dp, 345.85_dp, 343.7484_dp, 345.85_dp, &
    343.9697_dp, 343.7484_dp, 440.562_dp, 438.4418_dp, 440.562_dp, 438.7015_dp, &
    438.4418_dp, 546.705_dp, 544.5658_dp, 546.705_dp, 544.8628_dp, 544.5658_dp, &
    662.103_dp, 664.261_dp, 662.4368_dp, 791.0344_dp, 793.21_dp, 791.4045_dp, &
    931.339_dp, 933.533_dp, 931.745_dp, 1082.9941_dp, 1085.206_dp, 1083.4356_dp, &
    1245.975_dp, 1248.204_dp, 1246.4518_dp, 1420.2552_dp, 1422.502_dp, 1420.7672_dp &
/)

REAL (KIND=dp), DIMENSION(0:2*O2max-7)   :: O2b = (/ &
    0.000207_dp, 0.002045_dp, 0.00255_dp, 0.2308_dp, 0.003268_dp, 0.004338_dp, &
    0.006036_dp, 0.008967_dp, 0.05132_dp, 0.0147_dp, 0.02207_dp, 0.02842_dp, &
    0.01222_dp, 0.007753_dp, 0.005353_dp, 0.003916_dp, 0.0769_dp, 0.1077_dp, &
    0.002989_dp, 0.002356_dp, 0.001905_dp, 0.1615_dp, 0.017752_dp, 0.262699_dp, &
    0.1714_dp, 0.09234_dp, 0.06596_dp, 0.04342_dp, 0.001458_dp, 0.304466_dp, &
    0.2727_dp, 0.259873_dp, 0.02613_dp, 0.01979_dp, 0.32365_dp, 0.3077_dp, &
    0.303719_dp, 0.01387_dp, 0.00103_dp, 0.01127_dp, 0.33465_dp, 0.3251_dp, &
    0.323348_dp, 0.008577_dp, 0.007271_dp, 0.341777_dp, 0.3345_dp, &
    0.334499_dp, 0.005822_dp, 0.005076_dp, 0.346767_dp, 0.3422_dp, &
    0.341691_dp, 0.004209_dp, 0.003743_dp, 0.3504_dp, 0.3471_dp, 0.346713_dp, &
    0.003184_dp, 0.002874_dp, 0.3532_dp, 0.3506_dp, 0.35042_dp, 0.002492_dp, &
    0.002275_dp, 0.3555_dp, 0.3534_dp, 0.3533_dp, 0.002079_dp, 0.3573_dp, &
    0.3556_dp, 0.3589_dp, 0.3589_dp, 0.3574_dp, 0.3574_dp, 0.3601_dp, &
    0.3589_dp, 0.3589_dp, 0.3612_dp, 0.3602_dp, 0.3602_dp, 0.3622_dp, 0.3613_dp, &
    0.3613_dp, 0.363_dp, 0.3622_dp, 0.3622_dp, 0.3637_dp, 0.363_dp, 0.363_dp, &
    0.000373_dp, 0.002156_dp, 0.002705_dp, 0.1385_dp, 0.003493_dp, 0.004685_dp, &
    0.00661_dp, 0.010022_dp, 0.03991_dp, 0.01696_dp, 0.01867_dp, 0.03474_dp, &
    0.01078_dp, 0.007014_dp, 0.004924_dp, 0.003646_dp, 0.1077_dp, 0.5383_dp, &
    0.002808_dp, 0.002229_dp, 0.001812_dp, 0.2692_dp, 0.017752_dp, 0.4729_dp, &
    0.4_dp, 0.4617_dp, 0.09234_dp, 0.05583_dp, 0.001458_dp, 0.439784_dp, &
    0.4286_dp, 0.467772_dp, 0.03193_dp, 0.02339_dp, 0.423234_dp, 0.4196_dp, &
    0.438706_dp, 0.016_dp, 0.001854_dp, 0.01278_dp, 0.413391_dp, 0.4118_dp, &
    0.422839_dp, 0.009586_dp, 0.008037_dp, 0.406877_dp, 0.406_dp, 0.413205_dp, &
    0.006377_dp, 0.005517_dp, 0.40225_dp, 0.4017_dp, 0.406775_dp, 0.004546_dp, &
    0.00402_dp, 0.398795_dp, 0.3985_dp, 0.402187_dp, 0.003403_dp, 0.003059_dp, &
    0.3961_dp, 0.3959_dp, 0.3988_dp, 0.002643_dp, 0.002405_dp, 0.393979_dp, &
    0.3938_dp, 0.3961_dp, 0.002112_dp, 0.001941_dp, 0.3922_dp, 0.3921_dp, &
    0.394_dp, 0.001726_dp, 0.3908_dp, 0.3907_dp, 0.3922_dp, 0.3895_dp, &
    0.3895_dp, 0.3908_dp, 0.3885_dp, 0.3885_dp, 0.3896_dp, 0.3876_dp, &
    0.3876_dp, 0.3885_dp, 0.3868_dp, 0.3868_dp, 0.3876_dp, 0.3861_dp, &
    0.3861_dp, 0.3868_dp &    
/)

INTEGER, DIMENSION (0:O2maxJ)            :: O2JZ = (/ &
    0, 2, 1, 2, 4, 3, 4, 6, 5, 8, 6, 7, 10, 8, 9, 12, 10, 11, 14, 12, 13, &
    16, 14, 15, 18, 16, 17, 20, 18, 19, 22, 20, 21, 24, 22, 23, 26, 24, 25, &
    28, 26, 27, 30, 28, 29, 32, 30, 31, 34, 32, 33, 36, 34, 35 &
/)

INTEGER, DIMENSION (0:2*O2max-7)         :: O2J2 = (/ &
    4, 19, 17, 1, 15, 13, 11, 9, 3, 7, 5, 5, 7, 9, 11, 13, 3, 2, 15, 17, 19, &
    2, 2, 4, 3, 2, 3, 4, 4, 6, 5, 4, 5, 6, 8, 7, 6, 7, 4, 8, 10, 9, 8, 9, 10, &
    12, 11, 10, 11, 12, 14, 13, 12, 13, 14, 16, 15, 14, 15, 16, 18, 17, 16, 17, &
    18, 20, 19, 18, 19, 22, 21, 20, 24, 23, 22, 26, 25, 24, 28, 27, 26, 30, 29, &
    28, 32, 31, 30, 34, 33, 32, 2, 18, 16, 2, 14, 12, 10, 8, 4, 6, 6, 4, 8, 10, &
    12, 14, 2, 0, 16, 18, 20, 1, 2, 2, 1, 0, 2, 3, 4, 4, 3, 2, 4, 5, 6, 5, 4, &
    6, 2, 7, 8, 7, 6, 8, 9, 10, 9, 8, 10, 11, 12, 11, 10, 12, 13, 14, 13, 12, &
    14, 15, 16, 15, 14, 16, 17, 18, 17, 16, 18, 19, 20, 19, 18, 20, 22, 21, 20, &
    24, 23, 22, 26, 25, 24, 28, 27, 26, 30, 29, 28, 32, 31, 30 &
/)

INTEGER, DIMENSION (0:2*O2max-7)         :: O2shift = (/ &
    0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, &
    -2, -2, -12, -14, -14, -14, -16, -16, -24, -26, -26, -26, -26, -28, -35, &
    -37, -37, -37, -39, -40, -47, -49, -49, -49, -51, -58, -60, -60, -60, -62, &
    -70, -72, -72, -72, -74, -81, -83, -83, -83, -85, -93, -95, -95, -95, -97, &
    -104, -106, -106, -106, -108, -118, -118, -118, -129, -129, -129, -140, &
    -140, -140, -152, -152, -152, -163, -163, -163, -174, -174, -174, -186, &
    -186, -186, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
    12, 14, 14, 14, 16, 16, 24, 26, 26, 26, 26, 28, 35, 37, 37, 37, 39, 40, 47, &
    49, 49, 49, 51, 58, 60, 60, 60, 62, 70, 72, 72, 72, 74, 81, 83, 83, 83, 85, &
    93, 95, 95, 95, 97, 104, 106, 106, 106, 108, 116, 118, 118, 118, 120, 129, &
    129, 129, 140, 140, 140, 152, 152, 152, 163, 163, 163, 174, 174, 174, 186, &
    186, 186 &
/)

INTEGER, DIMENSION (0:2*N2Jmax-3)        :: N2shift = (/ &
    -12, -20, -28, -36, -44, -52, -60, -68, -76, -84, -91, -99, -107, -115, &
    -123, -131, -139, -147, -155, -163, -171, -179, -186, -194, -202, -210, &
    -218, 12, 20, 28, 36, 44, 52, 60, 68, 76, 84, 91, 99, 107, 115, 123, &
    131, 139, 147, 155, 163, 171, 179, 186, 194, 202, 210, 218 &
/) 

contains

! This subroutine combines spline and splint function
! in Numberical Recipes by Press et al., 1997.
SUBROUTINE BSPLINE(xa, ya, n, x, y, np, errstat)

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  INTEGER, INTENT (OUT)                     :: errstat
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y

  
  REAL (KIND=dp), DIMENSION(n)              :: y2a, xacpy
  REAL (KIND=dp), DIMENSION(np)             :: xcpy
  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp)                            :: xmin, xmax, xamin, xamax
  
  errstat = 0
  IF (n < 3) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF

  xmin = MINVAL(x); xmax = MAXVAL(x)
  xamin = MINVAL(xa); xamax = MAXVAL(xa)
  IF (xmin < xamin .OR. xmax > xamax) THEN
     errstat =  -3; RETURN
  ENDIF
  
  IF (xa(1) < xa(n)) THEN
     CALL SPLINE(xa, ya, n, y2a)
     CALL SPLINT(xa, ya, y2a, n, x, y, np)
  ELSE
     xacpy = -xa; xcpy = -x
     CALL SPLINE(xacpy, ya, n, y2a)
     CALL SPLINT(xacpy, ya, y2a, n, xcpy, y, np)
  ENDIF

  RETURN
END SUBROUTINE BSPLINE

! modified to always use "natural" boundary conditions
SUBROUTINE SPLINE (x, y, n, y2)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: x, y  
  REAL (KIND=dp), DIMENSION(n), INTENT(OUT) :: y2
  
  REAL (KIND=dp), DIMENSION(n) :: u
  INTEGER       :: i, k
  REAL(KIND=dp) :: sig, p, qn, un
  
  y2 (1) = 0.0
  u (1) = 0.0
  
  DO i = 2, n - 1
     sig = (x (i) - x (i - 1)) / (x (i + 1) -x (i - 1))
     p = sig * y2 (i - 1) + 2.D0
     y2 (i) = (sig - 1.) / p
     u (i) = (6._dp * ((y (i + 1) - y (i)) / (x (i + 1) - x (i)) -  & 
          (y (i) - y (i - 1)) / (x (i) - x (i - 1))) / (x (i + 1) - &
          x (i - 1)) - sig * u (i - 1)) / p
  ENDDO
  
  qn = 0.0
  un = 0.0
  y2 (n) = (un - qn * u (n - 1)) / (qn * y2 (n - 1) + 1.0)
  DO k = n - 1, 1, -1
     y2 (k) = y2 (k) * y2 (k + 1) + u (k)
  ENDDO
  
  RETURN
END SUBROUTINE SPLINE

! This code could be optimized if x is in increasing/descreasing order
SUBROUTINE SPLINT (xa, ya, y2a, n, x, y, m)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n, m
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: xa, ya, y2a
  REAL (KIND=dp), DIMENSION(m), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(m), INTENT(OUT):: y
  
  INTEGER        :: ii, klo, khi, k 
  REAL (KIND=dp) :: h, a, b

  !klo = 1; khi = n
  DO ii = 1, m 
     klo = 1; khi = n
    
     !IF ( khi - klo == 1) THEN
     !   IF (x(ii) > xa(khi) )   THEN
     !       khi = n
     !   ENDIF
     !ENDIF

     DO WHILE (khi - klo > 1)
        k = (khi + klo) / 2
        IF (xa (k) > x(ii)) THEN
           khi = k
        ELSE
           klo = k
        ENDIF
     ENDDO
     
     h = xa (khi) - xa (klo)
     IF (h == 0.0) STOP 'Bad xa input in: splint!!!'
     a = (xa (khi) - x(ii)) / h
     b = (x(ii) - xa (klo)) / h
     
     y(ii) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + &
          (b**3 - b) * y2a (khi)) * (h**2) / 6.0
  ENDDO
  
  RETURN
END SUBROUTINE SPLINT

! ==================================================================================!
! Purpose: This routine provides an interface to Christopher E. Sioris's single     !
! scattering first-order rotational raman scattering model (Sioris and Evans, 2001) !
! Xiong Liu, 07/08/2007                                                             !
!                                                                                   !
! Inputs/Outputs:                                                                   !
! nz: number of atmospheric layers                                                  !
! nw: number of wavelengths                                                         !
! nw_out: number of output wavelengths                                                         !
! sza: solar zenith angle  (degree)                                                 !
! vza: viewing zenith angle (degree)                                                !
! sca: scattering angle (degree, backscattering: 180 )                              !
! albedo: Lambertian surface albedo (only for satellite observations)               !
! do_upwelling: true for satellite observations else: ground-based measurements     !
! ts:  temperature profile fro TOA to BOS (K)                                       !
! wave: wavelengths (nm), increasing order of sol an taus                                          !
! wave_out: wavelengths of output (nm), increasing order of sol an taus                                          !
! sol:  solar irradiance                                                            !
! rhos: number of air molecules at each layer (molecules cm^-2)                     !
! taus: optical thickness at each wavelength and each layer                         !
! rspec: Ring spectrum at wave_out (i.e., relative diff. with and without Ring effect)  !
! problems: if true, errors occur during computing the Ring spectrum                !
!                                                                                   !
! Notes:                                                                            !
! (1) Ring effect is not computed (i.e., 0) for the few nms (~3nm) at each end. So  !
!     if you want to calculate Ring effect for 310-330 nm, you need to provide      !
!     sol/taus for 307-333 nm.                                                      !
! (2) If taus and sol are at high spectral resolution, the Ring effect is also at   !
!     high spectral resolution and you need to convolve and sample it to low        !
!     resolution. On the other hand, if taus and sol are at low resolution          !
!     (e.g. OMI), so is the ring spectrum. You just need to sample/interpolate it   !
!     to the grid of you want.                                                      !
! ==================================================================================! 

SUBROUTINE GET_RAMAN (nz, nw, nw_out, maxnu, sza, vza, sca, albedo, do_upwelling, ts, rhos, &
     wave, wave_out, sol, taus, rspec, problems) bind(C)
 
  IMPLICIT NONE

  ! ========================
  ! Input/output variables
  ! ========================
  INTEGER(c_int),        INTENT(IN)  :: nz, nw, nw_out, maxnu
  LOGICAL(c_bool),       INTENT(IN)  :: do_upwelling
  LOGICAL(c_bool),       INTENT(OUT) :: problems
  REAL (KIND=c_double),  INTENT(IN)  :: sza, sca, vza, albedo
  REAL (KIND=c_double),  DIMENSION(nz), INTENT(IN)     :: ts, rhos
  REAL (KIND=c_double),  DIMENSION(nw), INTENT(IN)     :: wave, sol
  REAL (KIND=c_double),  DIMENSION(nw_out), INTENT(IN)     :: wave_out
  REAL (KIND=c_double),  DIMENSION(nw, nz), INTENT(IN) :: taus
  REAL (KIND=c_double),  DIMENSION(nw_out),    INTENT(OUT) :: rspec

  ! ========================
  ! Local Variables
  ! ========================
  INTEGER,        PARAMETER :: NedgePos = 218
  REAL (KIND=dp), PARAMETER :: pi      = 3.14159265358979_dp
  REAL (KIND=dp), PARAMETER :: deg2rad = pi / 180.0_dp

  INTEGER                              :: nuhi, nulo, nu, i, j, errstat
  REAL (KIND=dp)                       :: scl, cosvza, cossza, tmpalb
  REAL (KIND=dp), DIMENSION(nz)        :: ctau
  REAL (KIND=dp), DIMENSION(nw, nz)    :: strans, vtrans
  REAL (KIND=dp), DIMENSION(MAXNU, nz) :: st, vt
  REAL (KIND=dp), DIMENSION(MAXNU)     :: ring, ramanwav

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=9), PARAMETER :: modulename = 'GET_RAMAN'
  
  problems = .FALSE.
  IF ( .NOT. do_upwelling) THEN
     ! Ignore the surface contribution for ground-based observations.
     tmpalb = 0.0
  ELSE
     tmpalb = albedo
  ENDIF
  
  ! Get position for raman calculation
  nuhi = INT(1.0D7 / wave(1))
  nulo = INT(1.0D7 / wave(nw)) + 1
  nu = 0
  DO i = nulo, nuhi
     nu = nu + 1
     ramanwav(nu) = 1.0D7 / REAL(i, KIND=dp)
  ENDDO
  ramanwav(nu) = wave(1); ramanwav(1) = wave(nw)
  
  IF (nuhi - nulo + 1 > MAXNU) THEN
     WRITE(*, *) modulename, ': Need to increase MAXNU!!! to at least ', (nuhi - nulo + 1), ' current value: ', MAXNU
     errstat = 1; RETURN
  ELSE IF (nuhi <= nulo) THEN
     WRITE(*, *) modulename, ': nulo>=nuhi, should never happen!!!'
     errstat = 1; RETURN
  ENDIF

  cossza = COS(sza * deg2rad); cosvza = COS(vza * deg2rad)

  ! Compute optical depth
  DO i = 1, nw 
     ctau = 0.0
     DO j = 2, nz
        ctau(j) = ctau(j-1) + taus(i, j)
     ENDDO
     
     scl = sol(i) * ((wave(i) / wave(1)) ** (-4.0))

     ! Assume plane parallel. In spherical geometry, you can replace 
     ! ctau/cossza with slant optical thickness
     strans(i, 1:nz) = scl * EXP(-ctau(1:nz) / cossza)
     IF (do_upwelling) THEN
        vtrans(i, 1:nz) = EXP(-ctau(1:nz) / cosvza)
     ELSE
        vtrans(i, 1:nz) = EXP(-(ctau(nz) - ctau(1:nz)) / cosvza)
     ENDIF
  ENDDO

  ! Interpolate to raman grid in wavenumber
  DO i = 1, nz
     CALL BSPLINE(wave, strans(:, i), nw, ramanwav(1:nu), st(1:nu, i), nu, errstat)     
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
        problems = .TRUE.; RETURN
     ENDIF

     CALL BSPLINE(wave, vtrans(:, i), nw, ramanwav(1:nu), vt(1:nu, i), nu, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
        problems = .TRUE.; RETURN
     ENDIF
  ENDDO

  ! Call raman program
  CALL RAMAN(nulo, nuhi, nu, nz, sca, tmpalb, ts, rhos, st(1:nu,1:nz), vt(1:nu,1:nz), ring(1:nu))
 
  ! Interpolate calculated ring back to input radiance grids
  CALL BSPLINE(ramanwav(1:nu), ring(1:nu), nu, wave_out(1:nw_out), rspec(1:nw_out), nw_out, errstat)
    
  ! Set edge pixels to zero
  DO i = 1, nw_out
     IF (wave_out(i) >= ramanwav(nu-NedgePos) ) THEN
        EXIT
     ELSE
        rspec(i) = 0.D0
     ENDIF
  ENDDO

  DO i = nw_out, 1, -1
     IF (wave_out(i) <= ramanwav(NedgePos+1) ) THEN
        EXIT
     ELSE
        rspec(i) = 0.D0
     ENDIF
  ENDDO

  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     problems = .TRUE.; RETURN
  ENDIF
  
  RETURN
  
END SUBROUTINE GET_RAMAN

! ======================================================================================
!The model is developed by C.E. Sioris (Sioris and Evans, 2001)
! It is converted to FORTRAN 90 by Xiong Liu
! 1. Converted to FORTRAN 90
! 2. Read input parameters only once
! 3. Optimize some computation
!
! Notes from the Chris
! ======================================================================================
! This has been swtiched to a forward model, but R is for elastic and I is for inelastic
! Christopher E. Sioris, Centre for Research in Earth and Space Science
! to appear in Can.J. Phys. 2001.  
! 
! CAVEATS: 1) apply only to obs at moderate spectral res (>=1cm^-1)
! 2) valid only for 1000nm>lambda>200nm
! 3) valid only in weakly absorbing regions or where absorption
! does not exhibit fine structure (i.e. Chappuis). Thus avoid 
! Hartley/Huggins, O2 and H2O bands        (not true anymore)
! 4) assumes no aerosols, no water vapour  (not true anymore)
! 5) assumes all scattering is in line-of-sight (i.e. ignores
! phase, polarization of multiply scattered radiation)
! 6) assumes an effective temperature for the entire atmosphere
! 
! newer versions of Raman scattering model address 1-6
! 
! Using this model
! 
! 1) Change nuhi and nulo to the appropriate upper and lower 
! wavenumber boundaries.
! 2) Change TH to the number of measured spectra
! 3) Save the temperatures(T) and single scattering angles(SSA) 
! as the file tempSSA.txt in the following format:
! T1 (in Kelvin)    SSA1 (in degrees)
! T2                SSA2
! 99.9999 (or any real number to prevent EOF error)
! SEE SAMPLE tempSSA.txt
! 4) Save the measured radiances as follows:   
! spectrum1 (1 data point/line, in ascending order wrt lambda) 
! spectrum2(don't skip any lines between spectra, see sample Int.txt) 
! ======================================================================================

SUBROUTINE raman(nulo, nuhi, nline, nz, sca, albedo, T, rhos, R, tran, ring)

  IMPLICIT NONE
  INTEGER, PARAMETER        :: dp = KIND(1.0D0)
  INTEGER, PARAMETER        :: maxpos=218
  REAL (KIND=dp), PARAMETER :: pi = 3.14159265358979_dp, O2mix= 0.20949858, &
       N2mix= 0.78079469, CO2mix=0.0003668  
  REAL (KIND=dp), PARAMETER :: c1=1.438769, NL=2.686763D19

  ! ========================
  ! Input/output variables
  ! ========================
  INTEGER,        INTENT(IN) :: nulo, nuhi, nline, nz
  REAL (KIND=dp), INTENT(IN) :: sca, albedo
  REAL (KIND=dp), DIMENSION(nulo:, :), INTENT(IN) :: R, tran
  REAL (KIND=dp), DIMENSION(nz),            INTENT(IN) :: T, rhos
  REAL (KIND=dp), DIMENSION(nline),        INTENT(OUT) :: ring

  ! ========================
  ! Local variables
  ! ========================
  INTEGER        :: nu, iz, j, k, fidx, lidx
  REAL (KIND=dp) :: ZN2, ZO2, temp, temp1, temp2, phasefnc, pi3, N2so, O2so, &
       N2sumin, O2sumin
  REAL (KIND=dp), DIMENSION (nz)          :: tempz
  REAL (KIND=dp), DIMENSION (nulo:nuhi)   :: gammaN2, gammaO2
  REAL (KIND=dp), DIMENSION(0:N2Jmax)     :: N2pop
  REAL (KIND=dp), DIMENSION(0:2*N2Jmax-3) :: N2csec
  REAL (KIND=dp), DIMENSION(0:O2maxJ)     :: O2popz
  REAL (KIND=dp), DIMENSION(0:2*O2max-7)  :: O2csec, O2pop
  REAL (KIND=dp), DIMENSION(nulo+maxpos:nuhi-maxpos)     :: RaylP, nr, &
       Raylro, I_tot, R_tot, e, Raylcsec, diff
  REAL (KIND=dp), DIMENSION(nulo+maxpos:nuhi-maxpos, nz) :: I
  REAL (KIND=dp), DIMENSION(0:2 * N2Jmax-3, nulo+maxpos:nuhi-maxpos) :: &
       N2so_temp, N2gamma_temp
  REAL (KIND=dp), DIMENSION(0:2 * O2max-7, nulo+maxpos:nuhi-maxpos) :: &
       O2so_temp, O2gamma_temp
  fidx = nulo + maxpos; lidx = nuhi - maxpos
  ! calculate dynamic optical parameters
  DO nu = nulo, nuhi
     temp = (1d0 * nu) ** 2 
     gammaN2(nu) = -6.01466E-25 + 2.38557E-14 / (1.86099E10 -temp)
     gammaO2(nu) = 7.149E-26    + 4.59364E-15 / (4.82716E9 - temp)  
  ENDDO
	
  pi3 = pi ** 3
  DO nu = fidx, lidx 
     temp   = (nu / 1d4) ** 2
     nr(nu) = 1d-4 * (0.7041 + 315.9 / (157.39 - temp ) ) &
          + 8.4127d-4 / (50.429 - temp )
     
     e(nu) = O2mix * (0.096 + 0.001385 * temp + 1.448d-4 * temp ** 2) &
          + (N2mix * (0.034 + 0.000317 * temp) + CO2mix * 0.15 )
     e(nu) = e(nu) * 4.5d0 !9.0d0 /2.0d0

     ! code doesn't work below ~316 nm because nu**2 > maxint, use real instead        
     Raylcsec(nu) = 32.0d0 * (REAL(nu))**4 * pi3 * (nr(nu))**2 * &
          (1.0d0 + e(nu) / 4.5d0) / 3.d0 / NL / NL     
     Raylro(nu)   = 6.0d0 * e(nu) / (45.d0 + 7.d0 * e(nu))         
  ENDDO

  ! calculate Rayleigh and Raman scattering phase functions
  temp = COS(sca * pi / 180.d0) ** 2
  phasefnc = (13.d0 + temp) * 3.d0 / 40.d0  
  RaylP(fidx:lidx) = (1.d0 + Raylro(fidx:lidx) + (1.d0 - Raylro(fidx:lidx) ) &
       * temp) * 3.d0 / (4.d0 + 2.d0 * Raylro(fidx:lidx) )         
  
  temp = 256 * pi ** 5 / 27

  ! We precompute some values for the iz loop below for efficiency
  DO nu = fidx, lidx
    temp1 = gammaN2(nu) * gammaN2(nu); temp2 = (REAL(nu))**4
    DO j = 0, 2 * N2Jmax-3
       N2so_temp(j, nu) = (REAL(nu - N2shift(j)))**4 * temp1
       N2gamma_temp(j, nu) = gammaN2( nu + N2shift(j)) ** 2 * temp2
    ENDDO

    temp1 = gammaO2(nu) * gammaO2(nu)
    DO k = 0, 2 * O2max-7
       O2so_temp(k, nu) = (REAL(nu - O2shift(k)))**4 * temp1
       O2gamma_temp(k, nu) = gammaO2(nu + O2shift(k)) ** 2 * temp2
    ENDDO
  ENDDO

  ! The following loop could be avoided by just using an effective temperature
  ! Will have small effect on the computed Ring effect spectrum (the effect on
  ! troospheric column ozone for one orbit is (0.014+/-0.06 DU)
  ! Use effective temperature, replace T(iz) with T

  DO iz = 1, nz
     !calculate partitioning of N2  	
     !iz = 12
     ZN2=0
     DO j = 0, N2Jmax
        IF ( ifix(10* (j / 2.0 - ifix( j / 2.0) ) ) == 5) THEN
           N2pop(j) = 3 * ( 2 * j + 1.0d0) * EXP(-c1 * N2E(j) / T(iz)) 
        ELSE
           N2pop(j) = 6 * ( 2 * j + 1.0d0) * EXP(-c1 * N2E(j) / T(iz))
        ENDIF
        ZN2= ZN2 + N2pop(j)  
     ENDDO
  
     ! calculate static part of Anti-Stokes cross sections for N2
     DO j = 0, N2Jmax - 2 
        N2csec(j) = temp * N2pop(j + 2) * N2b(j) /  ZN2
     ENDDO
  
     ! calculate static part of Stokes cross sections for N2	
     DO j = N2Jmax - 1, 2 * N2Jmax - 3 
        N2csec(j) = temp * N2pop(j - N2Jmax + 1) * N2b(j) / ZN2
     ENDDO
  
     ! calculate state sum for O2 
     ZO2 = 0
     DO k = 0, O2maxJ
        O2popz(k) = ( 2 * O2JZ(k) + 1 ) * EXP(- c1 * O2EnZ(k) / T(iz) )
        ZO2 = ZO2 + O2popz(k)
     ENDDO
       
     ! O2 cross sections
     DO k = 0, 2 * O2max-7
        O2pop(k) = (2 * O2J2(k) + 1 ) * EXP(-c1 * O2E(k) / T(iz) )
        O2csec(k)= temp * O2pop(k) * O2b(k) / ZO2
     ENDDO
  
  
     ! calculate relative amounts of light shifted in/out of a given nu     
     DO nu = fidx, lidx
        ! set accumulators to zero initially
        N2so = 0.0; N2sumin = 0.0
        O2so = 0.0; O2sumin = 0.0

        DO j = 0, 2 * N2Jmax-3
           N2so = N2so + N2csec(j) * N2so_temp(j, nu)
           N2sumin = N2sumin + R(nu + N2shift(j), iz) / R(nu, iz) &
                * N2csec(j) * N2gamma_temp(j, nu)
        ENDDO
     
        DO k = 0, 2 * O2max-7
           O2so = O2so + O2csec(k) * O2so_temp(k, nu)
           O2sumin = O2sumin + R(nu + O2shift(k), iz) / R(nu, iz) &
                * O2csec(k) * O2gamma_temp(k, nu)
        ENDDO
        
        diff(nu) = N2mix * (N2sumin - N2so) + O2mix * (O2sumin - O2so)
        I(nu, iz) = R(nu, iz) * ( 1.0 + phasefnc * diff(nu) &
             / (Raylcsec(nu) * RaylP(nu)))
     ENDDO
  ENDDO                     ! end nz loop
  
  DO nu = fidx, lidx
     tempz = rhos(1:nz) * tran(nu, 1:nz) 
     ! Note for ground-based observations, surface reflection is a second-order effect
     ! The surface albedo has already been set to zero in get_raman.f90
     I_tot(nu) = SUM(I(nu, 1:nz) * tempz(1:nz)) + I(nu, nz) * albedo * tempz(nz)
     R_tot(nu) = SUM(R(nu, 1:nz) * tempz(1:nz)) + R(nu, nz) * albedo * tempz(nz)
     ring(nu - nulo + 1) = I_tot(nu) / R_tot(nu) - 1.0d0
  ENDDO
  
  ring(1:maxpos) = 0.d0
  ring(nline - maxpos + 1 : nline) = 0.d0

  RETURN
  
END SUBROUTINE raman

end module raman_sioris
