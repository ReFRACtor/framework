module UVVSWIR_Eofpc_CreateBin_m

use iso_c_binding

!  ANNOTATION
!  ==========

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   UV   HISTORY
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!    25-27 October 2011: Version 1: 8-11 bins based on Gas Absorption only
!                        Worked well for Rayleigh, fractional errors 0.001 level
!                        Not so well for Aerosols, fractional errors 0.01  level
!    27 October 2011: Version 2: *** bins based on Gas Absorption, total scattering
!    06 August 2013, Added Version V0 for the 325-335 total ozone application
!                    Also good for the Continuum case

!  VN 3/16/19. Added three more Versions, from the old UVVis code.
!  V4 is the same as V3 in old UVVis code, with more focus around tau = 1, significantly better results
!  V5 is David Crisp's method for smart.
!  V6 is Peter Somkuti's method.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   TIR   HISTORY
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  May 2015. Versions V0, V1 and V2 were used for the Thermal Rad problem
!            V0 is a 2-bin job, V1 uses Fixed limits
!            V2 uses Max/Min increment in Log(taugas), + bisection.

!    11 January 2016, Added Version V3. Set Fixed bin limits to start, then
!                     then small-occupancy bins get merged to neighbors.
!                     Object is to avois special full-calculations.

!  VN 3/16/19. Added three more Versions, from the old UVVis code.
!  V4 is the same as V3 in old UVVis code, with more focus around tau = 1, significantly better results
!  V5 is David Crisp's method for smart.
!  V6 is Peter Somkuti's method.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   XOVER   HISTORY
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!   14 January 2016, Took over the TOR binnin of 1/11/16.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

private  midval
public   create_bin_UVVSWIR_V0, create_bin_UVVSWIR_V1, create_bin_UVVSWIR_V2, create_bin_UVVSWIR_V3, &
         create_bin_UVVSWIR_V4, create_bin_UVVSWIR_V5, create_bin_UVVSWIR_V6

contains

subroutine create_bin_UVVSWIR_V4 &
            ( E_nlayers, E_ndat, E_maxbins,           &
              ndat, nlay, nbin, gasdat, taudp, omega, abs_flag, &
              ncnt, index, bin ) bind(C)

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer(c_int), intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer(c_int), intent(in) :: NLAY, NBIN, NDAT

!  Optical data

   real(kind=c_double), intent(in) :: taudp  (E_nlayers,E_ndat)
   real(kind=c_double), intent(in) :: gasdat (E_nlayers,E_ndat)
   real(kind=c_double), intent(in) :: omega  (E_nlayers,E_ndat)
   integer(c_int), intent(in)       :: abs_flag(E_ndat)

!  outputs
!  -------

   INTEGER(c_int), intent(out) ::  NCNT(0:E_Maxbins), INDEX(E_ndat)
   INTEGER(c_int), intent(out) ::  BIN(NDAT)

!  local variables (Dynamic memory here!!)
!  ---------------

   INTEGER K, KC, IGBIN, NGBIN
   INTEGER W, COUNT, TEMPINDEX
   INTEGER BINDEX(NDAT), BIN_G(NDAT), NCNT_G(0:NBIN)
   DOUBLE PRECISION TAUTOT, STAUTOT(NDAT), MINIM, MAXIM, INCREMENT,TEMPARRAY(NDAT)
   DOUBLE PRECISION GASBINLIMS(0:NBIN), HALFVALUE, MID_STAU(0:NBIN), LIST_BINS(0:(NBIN-1))

!  Initialize

     NCNT = 0 ; NCNT_G = 0; INDEX = 0

      !write(*,*)'gasdat',SUM(GASDAT(:,1))
      !maxim= DLOG(SUM(GASDAT(:,1)))     !initialize minimum and maximum values to first total OD
      !minim= DLOG(SUM(GASDAT(:,1)))
      minim=DLOG(0.01_dp)
      !DO W=1,NDAT
      !TAUTOT = SUM(GASDAT(:,W))+1e-10
      !if (DLOG(TAUTOT).lt.minim) minim=DLOG(TAUTOT)
      !if (DLOG(TAUTOT).gt.maxim) maxim=DLOG(TAUTOT)
      !ENDDO
      !write(*,*)'minim',minim,maxim
      !increment=(maxim-minim)/dble(nbin)

!Pushkar's edit, assign gas bin limits before actual wavelenth loop starts

          !NGBIN = nbin/2
!      gasbinlims (0) = DLOG(0.01_dp)
!      gasbinlims (1) = DLOG(0.1_dp)
!      gasbinlims (2) = DLOG(0.25_dp)
!      gasbinlims (3) = DLOG(0.75_dp)
!      gasbinlims (4) = DLOG(1.0_dp)
!      gasbinlims (5) = DLOG(1.5_dp)
!      gasbinlims (6) = DLOG(10.0_dp)
!      gasbinlims (7) = DLOG(200000.0_dp)
!       LIST_BINS = (/ 0.01,1.0,10.0,1000.0/) !4
!       LIST_BINS = (/ 0.01,0.1,1.0,10.0,1000.0/) !5
!       LIST_BINS = (/ 0.01,0.1,0.5,1.0,10.0,1000.0/) !6
!      LIST_BINS = (/ 0.01,0.05,0.1,0.2,0.3,0.4,0.5,1.0,10.0,200000.0 /) !10
      LIST_BINS = (/ 0.01,0.025,0.05,0.1,0.25,0.5,0.625,0.75,1.0,5.0,10.0 /) !11 but this should be bad
!      LIST_BINS = (/ 0.01,0.025,0.05,0.1,0.25,0.5,0.625,0.75,1.0,10.0,200000.0 /) !11
!      LIST_BINS = (/ 0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.625,0.75,1.0,5.0,10.0 /) !13
!LIST_BINS = (/ 0.0004,0.0012,0.0041,0.0138,0.0462,0.0716,0.111,0.172,0.266,0.413,0.64,0.99,1.53,2.38,9.59,38.65,155.7,627.3 /) !Peter bins 18
      !Do W = 1, NBIN
      Do W = 1, 11!NBIN
      gasbinlims (W-1) = DLOG(LIST_BINS(W-1))
      !gasbinlims (W-1) = minim+w*increment
      ENDDO
      !write(*,*),'gasbinlims',2.718**gasbinlims

!  start wavelength loop

      DO W = 1, NDAT

!  basic input, each wavelength

         TAUTOT = SUM(GASDAT(:,W))+1e-10
         !write(*,*),tautot

!  find Gas bins

        if (nbin.eq.11) then
        Do k=0,nbin-1
        if (dlog(tautot).lt.gasbinlims(k)) then
        IGBIN = k
                exit
        endif
        IGBIN = nbin-1
        enddo
        endif

!Now accounting for dominant absorber

        if (nbin.eq.21) then
        Do k=0,10
        if (dlog(tautot).lt.gasbinlims(k)) then
                if (k.eq.0) then
                        IGBIN = k
                        exit
                elseif (abs_flag(w).eq.1) then
                        IGBIN = k
                        !write(*,*)'dom abs'
                        exit
                else
                        IGBIN = k+10
                        !write(*,*)'sec abs'
                        exit
                endif
        endif
        if (abs_flag(w).eq.1) then      !only enters this loop if optical depth is bigger than anything specific in binlims
                IGBIN = 10
        else
                IGBIN = 10+10
        endif
        enddo
        endif

!  get indices to backmap wavelengths to locations in bins

          NCNT_G(IGBIN) = NCNT_G(IGBIN)+1
          BIN_G(W) = IGBIN

!  end wavelength loop

      ENDDO
    !write(*,*),gasbinlims

!  Find scat bins

      do w = 1, ndat 
            K = BIN_G(W)
            NCNT(K) = NCNT(K)+1
            BIN(W) = K
            BINDEX(W) = NCNT(K)
      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
       IF (BIN(W) .EQ. 0) THEN 
          COUNT = BINDEX(W)
        ELSE
          COUNT = SUM(NCNT(0:BIN(W)-1))+BINDEX(W)
        ENDIF   
        
        INDEX(COUNT) = W
      ENDDO

!  debug

!      do k = 0, nbin 
!         write(*,*)k,ncnt(k)
!      enddo 
!      pause

!  finish

   return
end subroutine create_bin_UVVSWIR_V4

subroutine create_bin_UVVSWIR_V5 &
            ( E_nlayers, E_ndat, E_maxbins,           &
              ndat, nlay, nbin, gasdat, taudp, omega, &
              ncnt,index,bin )

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs   
!  ------

!  Dimensions (only for I/O)

   integer, intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer, intent(in) :: NLAY, NBIN, NDAT

!  Optical data

   real(kind=dp), intent(in) :: taudp  (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: gasdat (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: omega  (E_nlayers,E_ndat)

!  outputs
!  -------

   INTEGER, intent(out) ::  NCNT(0:E_Maxbins), INDEX(E_ndat)

!  local variables (Dynamic memory here!!)
!  ---------------

   INTEGER K, KC, IGBIN, NGBIN, I, CURBIN
   INTEGER W, COUNT
   INTEGER BIN(NDAT), BINDEX(NDAT), BIN_G(NDAT), NCNT_G(0:NBIN)
   DOUBLE PRECISION TAUTOT, STAUTOT(NDAT), INCREMENT
   DOUBLE PRECISION GASBINLIMS(0:NBIN)
   LOGICAL FLAG 

!  Initialize

     NCNT = 0 ; NCNT_G = 0; INDEX = 0; CURBIN = 2;
     NGBIN = 0

!  Start wavelength loop

      DO W = 1, NDAT  

!  Basic input, each wavelength

      FLAG = .TRUE.
      TAUTOT = SUM(GASDAT(:,W))

!  Find Gas bins

        IF (TAUTOT.lt.(0.01)) THEN
                IGBIN = 0
        ELSEIF (TAUTOT.gt.(10.0)) THEN
                IGBIN = 1
        ELSEIF (CURBIN.eq.2) THEN
                IGBIN = 2
                CURBIN = CURBIN + 1
        ELSE    
                DO I = 1, W-1
                        IF (ALL(abs(GASDAT(:,W)-GASDAT(:,I)).le.(GASDAT(:,I)*0.1))) THEN
                                IGBIN = BIN_G(I)
                                FLAG = .FALSE. !say you found the bin, now move on
                                EXIT
                        ELSEIF (I.eq.(W-1).and.FLAG) THEN
                                IGBIN = CURBIN  !doesnt match anyprevious profile, put in a new bin
                                CURBIN = CURBIN + 1
                        ENDIF
                ENDDO
        ENDIF

!  Get indices to backmap wavelengths to locations in bins

          NCNT_G(IGBIN) = NCNT_G(IGBIN)+1
          BIN_G(W) = IGBIN

!  End wavelength loop

      ENDDO
    !write(*,*),gasbinlims

!  Find scat bins

      do w = 1, ndat
         K = BIN_G(W) ; KC=K+NGBIN!KC = 2*K
            NCNT(K) = NCNT(K)+1
            BIN(W) = K
            BINDEX(W) = NCNT(K)
      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
       IF (BIN(W) .EQ. 0) THEN
          COUNT = BINDEX(W)
        ELSE
          COUNT = SUM(NCNT(0:BIN(W)-1))+BINDEX(W)
        ENDIF        
        INDEX(COUNT) = W
      ENDDO   

!  finish

   return
end subroutine create_bin_UVVSWIR_V5

subroutine create_bin_UVVSWIR_V6 &
            ( E_nlayers, E_ndat, E_maxbins,           &
              ndat, nlay, nbin, gasdat, taudp, omega, &
              ncnt,index,bin ) 

   IMPLICIT NONE

!  precision parameters  
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------       

!  Dimensions (only for I/O)

   integer, intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer, intent(in) :: NLAY, NBIN, NDAT

!  Optical data

   real(kind=dp), intent(in) :: taudp  (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: gasdat (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: omega  (E_nlayers,E_ndat)

!  outputs
!  -------

   INTEGER, intent(out) ::  NCNT(0:E_Maxbins), INDEX(E_ndat)

!  local variables (Dynamic memory here!!)
!  ---------------

   INTEGER K, KC, IGBIN, NGBIN
   INTEGER W, COUNT, TEMPINDEX
   INTEGER BIN(NDAT), BINDEX(NDAT), BIN_G(NDAT), NCNT_G(0:NBIN)
   DOUBLE PRECISION TAUTOT, STAUTOT(NDAT), MINIM, MAXIM, INCREMENT,TEMPARRAY(NDAT)
   DOUBLE PRECISION GASBINLIMS(0:NBIN), HALFVALUE, MID_STAU(0:NBIN)

!  Initialize

      NCNT = 0 ; NCNT_G = 0; INDEX = 0

      maxim= DLOG(SUM(GASDAT(:,1)))     !initialize minimum and maximum values to first total OD
      minim= DLOG(SUM(GASDAT(:,1)))

      DO W=1,NDAT
      TAUTOT = SUM(GASDAT(:,W))
      if (DLOG(TAUTOT).lt.minim) minim=DLOG(TAUTOT)
      if (DLOG(TAUTOT).gt.maxim) maxim=DLOG(TAUTOT)
      ENDDO
      increment=(maxim-minim)*2.00_dp/dble(nbin)

!Pushkar's edit, assign gas bin limits before actual wavelenth loop starts

      NGBIN = nbin/2  
      Do W = 1, NGBIN
      gasbinlims (W-1) = minim+w*increment
      ENDDO

!  start wavelength loop

      DO W = 1, NDAT

!  basic input, each wavelength

         TAUTOT = SUM(GASDAT(:,W))
         stautot(w) = TAUTOT/SUM(TAUDP(:,W))

!  find Gas bins

        IGBIN = CEILING((DLOG(TAUTOT)-MINIM)/INCREMENT)-1
        IF (IGBIN.lt.0) IGBIN = 0   !the one case where tautot is equal to the minimum value

!  get indices to backmap wavelengths to locations in bins

          NCNT_G(IGBIN) = NCNT_G(IGBIN)+1
          BIN_G(W) = IGBIN

!  end wavelength loop

      ENDDO

    write(*,*) ncnt_g

      do k=0, ngbin-1
         TEMPARRAY=0.00_dp
         TEMPINDEX=0
         if (ncnt_g(k).eq.0) then
            MID_STAU(K)=0
         else
            do w=1,ndat
               if (K.eq.BIN_G(w)) then
                  TEMPINDEX=TEMPINDEX+1
                  TEMPARRAY(TEMPINDEX) = STAUTOT(W)
               endif
            enddo
            call MIDVAL(TEMPARRAY(1:TEMPINDEX),TEMPINDEX,HALFVALUE) !finds more or less the middlemost value in the bin
            MID_STAU(K)=HALFVALUE
         endif
      enddo

!  Find scat bins

      do w = 1, ndat
         K = BIN_G(W) ; KC=K+NGBIN!KC = 2*K
!         HALFVALUE = (MAX_STAU(K) + MIN_STAU(K) ) * 0.5d0
         HALFVALUE = MID_STAU(K)
         if ( STAUTOT(W).GT.HALFVALUE) then!KC = 2*K+1
            NCNT(KC) = NCNT(KC)+1
            BIN(W) = KC
            BINDEX(W) = NCNT(KC)
         else
            NCNT(K) = NCNT(K)+1
            BIN(W) = K
            BINDEX(W) = NCNT(K)
         endif
      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
       IF (BIN(W) .EQ. 0) THEN
          COUNT = BINDEX(W)
        ELSE    
          COUNT = SUM(NCNT(0:BIN(W)-1))+BINDEX(W)
        ENDIF

        INDEX(COUNT) = W
      ENDDO

!  debug

!      do k = 0, nbin           
!         write(*,*)k,ncnt(k)
!      enddo
!      pause

!  finish

   return
end subroutine create_bin_UVVSWIR_V6

subroutine create_bin_UVVSWIR_V0 &
            ( E_nlayers, E_ndat, E_maxbins,      &
              ndat, nlay, nbin, binlims, taudat, &
              ncnt,index )

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer, intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer, intent(in) :: NLAY, NBIN, NDAT

!  Optical data

   real(kind=dp), intent(in) :: taudat (E_nlayers,E_ndat)

!  outputs
!  -------

   INTEGER, intent(out)       :: NCNT(0:E_Maxbins-1), INDEX(E_ndat)
   real(kind=dp), intent(out) :: BINLIMS(0:E_Maxbins-1)

!  local variables (Dynamic memory here!!)
!  ---------------

   INTEGER :: IBIN, W, COUNT, K
   INTEGER :: BIN(NDAT), BINDEX(NDAT)
   real(kind=dp) :: LOGTAUTOT(ndat), mxt, mnt, halft

!  Initialize

   NCNT = 0 ; BINLIMS = 0.0_dp; INDEX = 0

!  Bin Creation (10 November 2011, 1 or 2 bins
!  ------------

!  zero the counts

    DO K = 0, nbin-1
       NCNT(K) = 0
    ENDDO

!  For Nbin = 1, Find Max and Min of LOG(TAUTOT)

   if ( nbin .eq. 2 ) then
      MNT=+9000
      MXT=-MNT
      DO W = 1, NDAT
         LOGTAUTOT(W) = DLOG(SUM(taudat(1:NLAY,W)))
         MXT = MAX(MXT,LOGTAUTOT(W))
         MNT = MIN(MNT,LOGTAUTOT(W))
      ENDDO
      HALFT = 0.5d0 * (MXT+MNT)
    ENDIF

!  start wavelength loop

    DO W = 1, NDAT

!  1 or 2 bins

       IF (nbin.EQ.1) then
          IBIN = 0 ; BINLIMS(IBIN) = 999.9d0
       ELSE IF (nbin.EQ.2) THEN
          IF (LOGTAUTOT(W).LT.HALFT) THEN
             IBIN = 0 ; BINLIMS(IBIN) = HALFT
          ELSE
             IBIN = 1 ; BINLIMS(IBIN) = MXT
          ENDIF
       ENDIF

!  get indices to backmap wavelengths to locations in bins

       NCNT(IBIN) = NCNT(IBIN)+1
       BIN(W) = IBIN
       BINDEX(W) = NCNT(IBIN)

!  debug
!       write(33,*)w,bindex(w),bin(w)

!  end wavelength loop

    ENDDO

!  get the indices to change between bin space and wavenumber space

    DO W = 1, NDAT
       IF (BIN(W) .EQ. 0) THEN
          COUNT = BINDEX(W)
       ELSE
          COUNT = SUM(NCNT(1:BIN(W)-1))+BINDEX(W)
       ENDIF
       INDEX(COUNT) = W
    ENDDO

!  finish

   return
end subroutine create_bin_UVVSWIR_V0

!

subroutine create_bin_UVVSWIR_V2 &
            ( E_nlayers, E_ndat, E_maxbins,           &
              ndat, nlay, nbin, gasdat, taudp, omega, &
              ncnt,index,bin )

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer, intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer, intent(in) :: NLAY, NBIN, NDAT

!  Optical data

   real(kind=dp), intent(in) :: taudp  (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: gasdat (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: omega  (E_nlayers,E_ndat)

!  outputs
!  -------

   INTEGER, intent(out)       :: NCNT(0:E_Maxbins-1), INDEX(E_ndat)
   INTEGER, intent(out)       :: BIN(E_ndat)

!  local variables (Dynamic memory here!!)
!  ---------------

   INTEGER K, KC, IGBIN, NGBIN
   INTEGER W, COUNT, TEMPINDEX
   INTEGER BINDEX(NDAT), BIN_G(NDAT), NCNT_G(0:NBIN)
   DOUBLE PRECISION TAUTOT, STAUTOT(NDAT), MINIM, MAXIM, INCREMENT,TEMPARRAY(NDAT)
   DOUBLE PRECISION GASBINLIMS(0:NBIN), HALFVALUE, MID_STAU(0:NBIN)

!  Initialize

     NCNT = 0 ; NCNT_G = 0; INDEX = 0

!  v1 NOTES. 25 October 2011
!  From 290 to 340, total atmosphere optical depth of Ozone absoprtion
!   goes from ~11 down to ~0.008
!   about 5 to -3

!  v2 NOTES. 27 October 2011
!    total SCATTERING OPTICAL DEPTH VARIES FROM 1.1 TO 1.73 IN INITIAL EXAMPLE
!    Split all GasTau bins (except  < 300 nm) into 2 sub-omega bins

      maxim= DLOG(SUM(GASDAT(:,1)))     !initialize minimum and maximum values to first total OD
      minim= DLOG(SUM(GASDAT(:,1)))

      DO W=1,NDAT
      TAUTOT = SUM(GASDAT(:,W))
      if (DLOG(TAUTOT).lt.minim) minim=DLOG(TAUTOT)
      if (DLOG(TAUTOT).gt.maxim) maxim=DLOG(TAUTOT)
      ENDDO
      increment=(maxim-minim)*2.00_dp/dble(nbin)
       
!Pushkar's edit, assign gas bin limits before actual wavelenth loop starts
      NGBIN = nbin/2
      Do W = 1, NGBIN
      gasbinlims (W-1) = minim+w*increment
      ENDDO

!  start wavelength loop

      DO W = 1, NDAT

!  basic input, each wavelength

         TAUTOT = SUM(GASDAT(:,W))
         stautot(w) = 0.0d0
         do k = 1, nlay
            STAUTOT(W) = STAUTOT(W) + TAUDP(k,W)*OMEGA(k,W)
         enddo

  
!  find Gas bins

        IGBIN = CEILING((DLOG(TAUTOT)-MINIM)/INCREMENT)-1
        IF (IGBIN.lt.0) IGBIN = 0   !the one case where tautot is equal to the minimum value


!  get indices to backmap wavelengths to locations in bins

          NCNT_G(IGBIN) = NCNT_G(IGBIN)+1
          BIN_G(W) = IGBIN 
!  end wavelength loop

      ENDDO

!mick debug 
!      write(*,*)
!      write(*,'(1x,a,i2)') '# of possible gas bins: ',nbin+1
!      write(*,'(1x,a,i2)') '# of actual gas bins  : ',ngbin
!      write(*,*) '# of wavenumbers in each gas bin:'

      write(*,*) ncnt_g

      do k=0, ngbin-1
         TEMPARRAY=0.00_dp
         TEMPINDEX=0
         if (ncnt_g(k).eq.0) then
            MID_STAU(K)=0
         else
            do w=1,ndat
               if (K.eq.BIN_G(w)) then
                  TEMPINDEX=TEMPINDEX+1
                  TEMPARRAY(TEMPINDEX) = STAUTOT(W) 
               endif
            enddo
            call MIDVAL(TEMPARRAY(1:TEMPINDEX),TEMPINDEX,HALFVALUE) !finds more or less the middlemost value in the bin
            MID_STAU(K)=HALFVALUE
         endif
      enddo

!  Find scat bins

      do w = 1, ndat
         K = BIN_G(W) ; KC=K+NGBIN!KC = 2*K
!         HALFVALUE = (MAX_STAU(K) + MIN_STAU(K) ) * 0.5d0
         HALFVALUE = MID_STAU(K)
         if ( STAUTOT(W).GT.HALFVALUE) then!KC = 2*K+1
            NCNT(KC) = NCNT(KC)+1
            BIN(W) = KC
            BINDEX(W) = NCNT(KC)
         else 
            NCNT(K) = NCNT(K)+1
            BIN(W) = K
            BINDEX(W) = NCNT(K)
         endif
      ENDDO
       
!  get the indices to change between bin space and wavenumber space
        
      DO W = 1, NDAT
       IF (BIN(W) .EQ. 0) THEN
          COUNT = BINDEX(W)
        ELSE
          COUNT = SUM(NCNT(0:BIN(W)-1))+BINDEX(W)
        ENDIF
        
        INDEX(COUNT) = W
      ENDDO
      
!  debug

!      do k = 0, nbin
!         write(*,*)k,ncnt(k)
!      enddo
!      pause

!  finish

!   return
end subroutine create_bin_UVVSWIR_V2

subroutine create_bin_UVVSWIR_V1 &
            ( E_nlayers, E_ndat, E_maxbins,             &
              ndat, nlay, nbin, binlims, gasdat, omega, &
              ncnt,index )

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer, intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer, intent(in) :: NLAY, NBIN, NDAT

!  Optical data

   real(kind=dp), intent(in) :: gasdat (E_nlayers,E_ndat)
   real(kind=dp), intent(in) :: omega  (E_nlayers,E_ndat)

!  outputs
!  -------

   INTEGER, intent(out)       :: NCNT(0:E_Maxbins), INDEX(E_ndat)
   real(kind=dp), intent(out) :: BINLIMS(0:E_Maxbins-1)

!  local variables (Dynamic memory here!!)
!  ---------------

      INTEGER IBIN
      INTEGER W, COUNT
      INTEGER BIN(NDAT), BINDEX(NDAT)
      DOUBLE PRECISION TAUTOT, OMEGA1W
      DOUBLE PRECISION minim, maxim, increment



!  Initialize

     NCNT = 0 ; BINLIMS = 0.0_dp; INDEX = 0

!  Error

    if ( nbin.ne.11 ) stop 'This Create-bin routine requires 11 bins'

!  NOTES. 25 October 2011
!  From 290 to 340, total atmosphere optical depth of Ozone absoprtion
!   goes from ~11 down to ~0.008
!   about 5 to -3

!Pushkar's additions, slightly more general way of choosing bins
maxim= DLOG(SUM(GASDAT(:,1)))     !initialize minimum and maximum values to first total OD
minim= DLOG(SUM(GASDAT(:,1)))
      DO W=1,NDAT
      TAUTOT = SUM(GASDAT(:,W))
      if (DLOG(TAUTOT).lt.minim) minim=DLOG(TAUTOT)
      if (DLOG(TAUTOT).gt.maxim) maxim=DLOG(TAUTOT)
      
      ENDDO

increment = (maxim-minim)/11.0_dp
!  start wavelength loop

      DO W = 1, NDAT

!  basic input, each wavelength

        TAUTOT = SUM(GASDAT(:,W))
        OMEGA1W = OMEGA(NLAY,W)

!  find bin


!CHANGED TO UNIFORM BINS
          IF (DLOG(TAUTOT) .LT. (MINIM+INCREMENT)) THEN
            IBIN = 0 ; BINLIMS(IBIN) = (MINIM+INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+2*INCREMENT)) THEN
            IBIN = 1 ; BINLIMS(IBIN) = (MINIM+2*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+3*INCREMENT)) THEN
            IBIN = 2 ; BINLIMS(IBIN) = (MINIM+3*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+4*INCREMENT)) THEN
            IBIN = 3 ; BINLIMS(IBIN) = (MINIM+4*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+5*INCREMENT)) THEN
            IBIN = 4 ; BINLIMS(IBIN) = (MINIM+5*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+6*INCREMENT)) THEN
            IBIN = 5 ; BINLIMS(IBIN) = (MINIM+6*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+7*INCREMENT)) THEN
            IBIN = 6 ; BINLIMS(IBIN) = (MINIM+7*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+8*INCREMENT)) THEN
            IBIN = 7 ; BINLIMS(IBIN) = (MINIM+8*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+9*INCREMENT)) THEN
            IBIN = 8 ; BINLIMS(IBIN) = (MINIM+9*INCREMENT)
          ELSE IF (DLOG(TAUTOT) .LT. (MINIM+10*INCREMENT)) THEN
            IBIN = 9 ; BINLIMS(IBIN) = (MINIM+10*INCREMENT)
          ELSE 
            IBIN = 10 ; BINLIMS(IBIN) = MAXIM
          ENDIF

 !         IF (DLOG(TAUTOT) .LT. -4.0D0) THEN
  !          IBIN = 0 ; BINLIMS(IBIN) = -4.0d0
  !        ELSE IF (DLOG(TAUTOT) .LT. -3.0D0) THEN
  !          IBIN = 1 ; BINLIMS(IBIN) = -3.0d0
   !       ELSE IF (DLOG(TAUTOT) .LT. -2.0D0) THEN
   !         IBIN = 2 ; BINLIMS(IBIN) = -2.0d0
    !      ELSE IF (DLOG(TAUTOT) .LT. -1.5D0) THEN
     !       IBIN = 3 ; BINLIMS(IBIN) = -1.5d0
       !   ELSE IF (DLOG(TAUTOT) .LT. -1.0D0) THEN
      !      IBIN = 4 ; BINLIMS(IBIN) = -1.0d0
        !  ELSE IF (DLOG(TAUTOT) .LT. -0.5D0) THEN
         !   IBIN = 5 ; BINLIMS(IBIN) = -0.5d0
          !ELSE IF (DLOG(TAUTOT) .LT. 0.0D0) THEN
           ! IBIN = 6 ; BINLIMS(IBIN) = 0.0d0
          !ELSE IF (DLOG(TAUTOT) .LT. 0.5D0) THEN
           ! IBIN = 7 ; BINLIMS(IBIN) = 0.5d0
          !ELSE IF (DLOG(TAUTOT) .LT. 1.0D0) THEN
           ! IBIN = 8 ; BINLIMS(IBIN) = 1.0d0
          !ELSE IF (DLOG(TAUTOT) .LT. 2.0D0) THEN
           ! IBIN = 9 ; BINLIMS(IBIN) = 2.0d0
          !ELSE 
           ! IBIN = 10 ; BINLIMS(IBIN) = 10.0d0
          !ENDIF


!  get indices to backmap wavelengths to locations in bins

        NCNT(IBIN) = NCNT(IBIN)+1
        BIN(W) = IBIN
        BINDEX(W) = NCNT(IBIN)

!  end wavelength loop

      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
        IF (BIN(W) .EQ. 0) THEN
          COUNT = BINDEX(W)
        ELSE
          COUNT = SUM(NCNT(0:BIN(W)-1))+BINDEX(W)
        ENDIF
        INDEX(COUNT) = W
      ENDDO

!  finish

   return
end subroutine create_bin_UVVSWIR_V1

subroutine create_bin_UVVSWIR_Rob_8bin_V3    &
            ( E_nlayers, E_ndat, E_maxbins,  &
              ndat, nbin, binlims, gasdat,   &
              ncnt, index, bin )

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer, intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers

   integer, intent(in) :: NBIN, NDAT

!  Optical data

   real(kind=dp), intent(in) :: gasdat (E_nlayers,E_ndat)

!  outputs
!  -------

   INTEGER, intent(out)       :: NCNT(0:E_Maxbins-1), INDEX(E_ndat), BIN(E_ndat)
   real(kind=dp), intent(out) :: BINLIMS(0:E_MaxBins-1)

!  local variables (Dynamic memory here!!)
!  ---------------

      INTEGER IBIN
      INTEGER W, COUNT, K
      INTEGER BINDEX(NDAT)
      DOUBLE PRECISION TAUTOT, DATLIMS(8), LOGDATLIMS(8)

      Data Datlims / 0.40d0, 0.60d0, 0.8d0, 1.0d0, 1.3d0, 1.7d0, 4.0d0, 20000.0d0 /

!  Initialize

      NCNT = 0 ; BINLIMS = 0.0_dp; INDEX = 0

!  Error if not the right number of bins.....

      if ( nbin.ne.8 ) stop 'This Create-bin routine requires 8 bins'
      do K = 1, nbin
         LOGDATLIMS(K) = LOG(DATLIMS(K))
!write(*,*) 'K = ',K,' LOGDATLIMS(K) = ',LOGDATLIMS(K)
      ENDDO

!  start wavelength loop

      DO W = 1, NDAT

!  basic input, each wavelength

        TAUTOT = SUM(GASDAT(:,W))

!mick fix 2/16/2015 - avoid singularity
          IF (TAUTOT .LT. 1.0D-12) TAUTOT = 1.0D-12

!  find bin

          IF (DLOG(TAUTOT) .LT. LOGDATLIMS(1) ) THEN
             IBIN = 0 ; BINLIMS(IBIN) = LOGDATLIMS(1)
          ELSE IF (DLOG(TAUTOT) .LT. LOGDATLIMS(2) ) THEN
             IBIN = 1 ; BINLIMS(IBIN) = LOGDATLIMS(2)
          ELSE IF (DLOG(TAUTOT) .LT. LOGDATLIMS(3) ) THEN
             IBIN = 2 ; BINLIMS(IBIN) = LOGDATLIMS(3)
          ELSE IF (DLOG(TAUTOT) .LT. LOGDATLIMS(4) ) THEN
             IBIN = 3 ; BINLIMS(IBIN) = LOGDATLIMS(4)
          ELSE IF (DLOG(TAUTOT) .LT. LOGDATLIMS(5) ) THEN
             IBIN = 4 ; BINLIMS(IBIN) = LOGDATLIMS(5)
          ELSE IF (DLOG(TAUTOT) .LT. LOGDATLIMS(6) ) THEN
             IBIN = 5 ; BINLIMS(IBIN) = LOGDATLIMS(6)
          ELSE IF (DLOG(TAUTOT) .LT. LOGDATLIMS(7) ) THEN
             IBIN = 6 ; BINLIMS(IBIN) = LOGDATLIMS(7)
          ELSE 
             IBIN = 7 ; BINLIMS(IBIN) = LOGDATLIMS(8)
          endif

!  get indices to backmap wavelengths to locations in bins

        NCNT(IBIN) = NCNT(IBIN)+1
        BIN(W) = IBIN
        BINDEX(W) = NCNT(IBIN)

!  end wavelength loop

      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
        IF (BIN(W) .EQ. 0) THEN
          COUNT = BINDEX(W)
        ELSE
          COUNT = SUM(NCNT(0:BIN(W)-1))+BINDEX(W)
        ENDIF
        INDEX(COUNT) = W
      ENDDO

!  finish

   return
end subroutine create_bin_UVVSWIR_Rob_8bin_V3

subroutine create_bin_UVVSWIR_V3 &
            ( E_nlayers, E_ndat, E_maxbins,      &
              ndat, nlay, nbin, binlims, gasdat, &
              ncnt_new, index_new, bin_new ) bind(C)

   IMPLICIT NONE

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer(c_int), intent(in) :: E_nlayers, E_ndat, E_maxbins

!  Numbers
!     1/11/16 Flexible number of bins, NBIN isIntent(inout)

   integer(c_int), intent(in)    :: NLAY, NDAT
   integer(c_int), intent(inout) :: NBIN

!  Optical data

   real(kind=c_double), intent(in) :: gasdat (E_nlayers,E_ndat)

!  outputs
!  -------

!  1/11/16. Changed dimensioning on BINLIMS to symbolic 

   INTEGER(c_int), intent(inout) :: BIN_NEW(NDAT)
   INTEGER(c_int), intent(out)   :: NCNT_NEW(0:E_Maxbins), INDEX_NEW(E_ndat)
   real(kind=c_double), intent(out)  :: BINLIMS(0:E_Maxbins)

!  local variables (Dynamic memory here!!)
!  ---------------

      logical :: iterating
      INTEGER :: IBIN, W, COUNT, K, NBIN_OLD, NBIN_NEW, THRESH_COUNT
      INTEGER :: NCNT_OLD(0:E_Maxbins), INDEX_OLD(E_ndat)
      INTEGER :: BIN_OLD(NDAT), BINDEX_OLD(NDAT), BINDEX_NEW(NDAT)

      DOUBLE PRECISION TAUTOT, DTAU, REAL_INITIAL_BINLIMS(8), INITIAL_BINLIMS(8)
      Data REAL_INITIAL_BINLIMS / 0.01d0, 0.1d0, 0.213d0, 0.467d0, 1.0d0, 2.13d0, 4.67d0, 10.0d0 /

!  Threshholds

!      Data THRESH_COUNT / 30 /     ! Set 1, 12 January 2016
!      Data THRESH_COUNT / 50 /      ! Set 2, 12 January 2016
      Data THRESH_COUNT / 30 /      ! Set 2, 12 January 2016

!  Initialize

      DO IBIN = 1, 8
         INITIAL_BINLIMS(IBIN) = DLOG(REAL_INITIAL_BINLIMS(IBIN))
         BINLIMS(IBIN-1) = INITIAL_BINLIMS(IBIN)
      ENDDO
      NCNT_NEW = 0 ; BINLIMS = 0.0d0 ; INDEX_NEW = 0 ; BIN_NEW = 0

!  initial value of NBIN should always be 9

      NBIN_OLD = NBIN ; NCNT_OLD = 0
      if ( NBIN .ne.9 ) stop 'bad boy, change input of NBIN, Brian'

!  NOTES. 29 May 2015
!  mick fix - changed argument passing on GASDAT from ":" to "1:nlay" at several locations
!             [e.g., GASDAT(:,1) to GASDAT(1:nlay,1)] 

!  start wavelength loop

      DO W = 1, NDAT

!  basic input, each wavelength

        TAUTOT = SUM(GASDAT(1:nlay,W)) ; DTAU = DLOG(TAUTOT)

        IF (DTAU .LT. INITIAL_BINLIMS(1)) THEN
          IBIN = 0 ; BINLIMS(IBIN) = INITIAL_BINLIMS(1)
        ELSE IF (DTAU .LT. INITIAL_BINLIMS(2)) THEN
          IBIN = 1 ; BINLIMS(IBIN) = INITIAL_BINLIMS(2)
        ELSE IF ( DTAU .LT. INITIAL_BINLIMS(3)) THEN
          IBIN = 2 ; BINLIMS(IBIN) = INITIAL_BINLIMS(3)
        ELSE IF (DTAU .LT. INITIAL_BINLIMS(4)) THEN
          IBIN = 3 ; BINLIMS(IBIN) = INITIAL_BINLIMS(4)
        ELSE IF (DTAU .LT. INITIAL_BINLIMS(5)) THEN
          IBIN = 4 ; BINLIMS(IBIN) = INITIAL_BINLIMS(5)
        ELSE IF (DTAU .LT. INITIAL_BINLIMS(6)) THEN
          IBIN = 5 ; BINLIMS(IBIN) = INITIAL_BINLIMS(6)
        ELSE IF (DTAU .LT. INITIAL_BINLIMS(7)) THEN
          IBIN = 6 ; BINLIMS(IBIN) = INITIAL_BINLIMS(7)
        ELSE IF (DTAU .LT. INITIAL_BINLIMS(8)) THEN
          IBIN = 7 ; BINLIMS(IBIN) = INITIAL_BINLIMS(8)
        ELSE 
          IBIN = 8 ; BINLIMS(IBIN) = 50000.0d0
        ENDIF
        NCNT_OLD(IBIN) = NCNT_OLD(IBIN)+1
        BIN_OLD(W) = IBIN
        BINDEX_OLD(W) = NCNT_OLD(IBIN)
      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
        IF (BIN_OLD(W) .EQ. 0) THEN
          COUNT = BINDEX_OLD(W)
        ELSE
          COUNT = SUM(NCNT_OLD(0:BIN_OLD(W)-1)) + BINDEX_OLD(W)
        ENDIF
        INDEX_OLD(COUNT) = W
      ENDDO

!  debug before
!      write(*,*)'OLD NCNT, Sum/Details = ',SUM(NCNT_OLD(0:NBIN_OLD-1)),NCNT_OLD(0:NBIN_OLD-1)
!      do w = 1, ndat
!         write(77,*)W, BINDEX_OLD(W), BIN_OLD(W), DLOG(SUM(GASDAT(1:nlay,W)))
!      enddo

!  continuation point

567   continue

!  Now reduce the number of bins by 1 (from the top), fills up penumltimate bin
!    LAST BIN ONLY
      
      K = NBIN_OLD-1 ; iterating = .false.
      if ( NCNT_OLD(K) .lt. THRESH_COUNT ) THEN
        ITERATING = .true.
        NBIN_NEW = NBIN_OLD - 1
        NCNT_NEW(0:K-2) = NCNT_OLD(0:K-2)
        NCNT_NEW(K-1) = NCNT_OLD(K-1) !+ NCNT_OLD(K)
      endif

!  Not iterating : copy Old to new and exit

      if ( .not. iterating ) then
        NBIN_NEW = NBIN_OLD
        NCNT_NEW = NCNT_OLD
        DO W = 1, NDAT
          BINDEX_NEW(W) = BINDEX_OLD(W)
          BIN_NEW(W)    = BIN_OLD(W)
          IF (BIN_NEW(W) .EQ. 0) THEN
            COUNT = BINDEX_NEW(W)
          ELSE
            COUNT = SUM(NCNT_NEW(0:BIN_NEW(W)-1)) + BINDEX_NEW(W)
          ENDIF
          INDEX_NEW(COUNT) = W
        enddo
        NBIN = NBIN_NEW
        RETURN
      endif

!  Fill up BIN_NEW

      DO W = 1, NDAT
        IF ( BIN_OLD(W).eq. NBIN_OLD - 1) then
          IBIN = NBIN_NEW - 1
          BIN_NEW(W) = IBIN
          NCNT_NEW(IBIN) = NCNT_NEW(IBIN) + 1
          BINDEX_NEW(W)  = NCNT_NEW(IBIN)
        ELSE
          BINDEX_NEW(W) = BINDEX_OLD(W)
          BIN_NEW(W)    = BIN_OLD(W)
        ENDIF
      ENDDO

!  get the indices to change between bin space and wavenumber space

      DO W = 1, NDAT
        IF (BIN_NEW(W) .EQ. 0) THEN
          COUNT = BINDEX_NEW(W)
        ELSE
          COUNT = SUM(NCNT_NEW(0:BIN_NEW(W)-1)) + BINDEX_NEW(W)
        ENDIF
        INDEX_NEW(COUNT) = W
      ENDDO

!  Re-assign output

      NBIN = NBIN_NEW
      BINLIMS(NBIN_NEW) = BINLIMS(NBIN_OLD) 
      BINLIMS(NBIN_OLD) = 0.0d0 

!  debug after

!      write(*,*)'NEW NCNT, Sum/Details = ',SUM(NCNT_NEW(0:NBIN_NEW-1)),NCNT_NEW(0:NBIN_NEW-1)
!      do w = 1, ndat
!         write(78,*)W, BINDEX_NEW(W), BIN_NEW(W), DLOG(SUM(GASDAT(1:nlay,W)))
!      enddo

!  Re-set and return to iteration

      NBIN_OLD = NBIN_NEW
      NCNT_OLD = 0
      NCNT_OLD(0:NBIN_OLD-1) = NCNT_NEW(0:NBIN_OLD-1)
      BIN_OLD(1:NDAT)    = BIN_NEW(1:NDAT)
      BINDEX_OLD(1:NDAT) = BINDEX_NEW(1:NDAT) 

      go to 567
      
!  finish

   return
end subroutine create_bin_UVVSWIR_V3


subroutine MIDVAL (array,nelements,outmid)

   IMPLICIT NONE

!precision parameters
!====================

   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions (only for I/O)

   integer, intent(in) :: nelements

   real(kind=dp), intent(in) :: array(nelements)
   
   real(kind=dp), intent(out) :: outmid

! Local parameters

   integer ::  I,J
   double precision :: locarray(nelements),temp

   locarray=array

   DO I = 1, nelements-1
      DO J = 1, nelements-1
        IF (locarray(J) > locarray(J+1)) THEN
            TEMP=locarray(J)
            locarray(J)=locarray(J+1)
            locarray(J+1)=TEMP
         END IF
      END DO
   END DO

   outmid = locarray(nint(nelements/2.0))

RETURN
end subroutine MIDVAL

!  finish module

end module UVVSWIR_Eofpc_CreateBin_m
