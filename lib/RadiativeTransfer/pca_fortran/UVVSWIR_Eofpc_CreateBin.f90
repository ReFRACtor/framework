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
! V0 and V1 are outdated. V2 is breaking up entire region into 10 uniform bins in log
! V3 is similar, with more focus around tau = 1, significantly better results
! V4 is David Crisp's method for smart.
! V5 is Peter Somkuti's method.

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   XOVER   HISTORY
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!   14 January 2016, Took over the TOR binnin of 1/11/16.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

private  midval
public   create_bin_UVVSWIR_V0, create_bin_UVVSWIR_V1, create_bin_UVVSWIR_V2, create_bin_UVVSWIR_V3, &
         create_bin_UVVSWIR_V4, create_bin_UVVSWIR_V5

contains

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
      !write(*,*),ncnt_g
   
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

!

subroutine create_bin_UVVSWIR_V1 &
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
       
end subroutine create_bin_UVVSWIR_V1

!

!attempted quick fix for large PCA errors with extra bin around tau = 1, ASSUMES 8 bins
subroutine create_bin_UVVSWIR_V3 &
            ( E_nlayers, E_ndat, E_maxbins,           &
              ndat, nlay, nbin, gasdat, taudp, omega, abs_flag, &
              ncnt,index,bin ) bind(c)

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

   integer(c_int), intent(out) ::  NCNT(0:E_Maxbins), INDEX(E_ndat), BIN(E_ndat)

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
                
!   return
end subroutine create_bin_UVVSWIR_V3

!

subroutine create_bin_UVVSWIR_V4 & 
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

end subroutine create_bin_UVVSWIR_V4

!

!VERSION 5 START HERE
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
      !write(*,*)ncnt_g

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
end subroutine create_bin_UVVSWIR_V5

!

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
