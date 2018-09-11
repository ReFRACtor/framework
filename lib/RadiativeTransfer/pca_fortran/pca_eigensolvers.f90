module pca_eigensolvers_m

!   Version 3. Robustness and Recoding, 08 November 2012
!     Maximum dimensions are given, but need not be used

!  History of PCA developments:

!   27 January  2012. Alternative EOF Scheme, using Eigenmethods.
!   06 February 2012. Proper exception handling
!   24 May      2012. Albedo Closure included
!   23 April    2013. Continuum-only PCA routines (no second bin). Renamed, 8/6/13

!  Subroutines:
!  ------------

!    PUBLIC
!          pca_eigensolver      --> PCA for atmospheric optical properties only 
!          pca_eigensolver_alb  --> PCA for atmospheric optical properties and surface albedo 
!          pca_eigensolver_Continuum      --> PCA for atmospheric optical properties Continuum only 
!          pca_eigensolver_Continuum_alb  --> PCA for atmospheric optical properties Continuum and surface albedo 
!          pca_eigensolver_aer  --> PCA for atmospheric optical properties including aerosol Continuum
!          pca_eigensolver_Continuum_aer  --> PCA for atmospheric optical properties including aerosol Continuum

!    PRIVATE
!          Prepare_Eigenmatrix  --> Prepares Autocorrelation matrix (Eigenmatrix)

!    AUXILIARY calls (from module pca_auxiliaries)
!          PCA_Ranker           --> Ordering routine from Numerical Recipes
!          PCA_ASYMTX           --> Generic Eigensolver from DISORT/LIDORT

use pca_auxiliaries_m, only : PCA_ASYMTX, PCA_Ranker
use iso_c_binding

private
public pca_eigensolver_alb, pca_eigensolver, pca_eigensolver_aer,  pca_eigensolver_Continuum_alb, pca_eigensolver_Continuum

!  precision parameters

!      INTEGER, PARAMETER :: sp     = selected_real_kind(6)
      INTEGER, PARAMETER :: dp     = selected_real_kind(15)

contains

subroutine pca_eigensolver_alb &
           ( Max_Eofs, maxpoints, maxlayers, maxlayers21,   & ! Input dimensions
             n_Eofs, npoints, nlayers, nlayers2, nlayers21, & ! Input control
             taudp, omega, albedo,                          & ! Input optical properties
             Atmosmean, Albmean, Eofs, PrinComps,           & ! Outputs
             fail, message_len, message, trace_len, trace ) bind(C) ! Exception handling

!  PCA using Eigenproblem methods

   implicit none

!  inputs
!  ------

!  Dimensioning

   integer(c_int), intent(in)       :: Max_Eofs, maxpoints, maxlayers, maxlayers21

!  Control

   integer(c_int), intent(in)       :: n_Eofs, npoints, nlayers, nlayers2, nlayers21

!  Optical properties

   real(kind=c_double), intent(in) :: taudp(maxlayers,maxpoints)
   real(kind=c_double), intent(in) :: omega(maxlayers,maxpoints)
   real(kind=c_double), intent(in) :: albedo(maxpoints)

!  outputs
!  -------

!  Mean values

   real(kind=c_double), intent(out) :: atmosmean(maxlayers,2),albmean

!  EOFS and Principal Components (Latter is allocatable)

   real(kind=c_double), intent(out) :: eofs  (Max_Eofs,maxlayers21)
   real(kind=c_double), intent(out) :: PrinComps(n_eofs,npoints)

!  Exception handling and status

   logical(c_bool), intent(out)       :: fail
   integer(c_int), intent(in) :: message_len, trace_len
   character(kind=c_char), intent(out) :: message(message_len)
   character(kind=c_char), intent(out) :: trace(trace_len)

!  local data arrays
!  -----------------

   real(kind=dp) ::  data  (npoints,nlayers21)
   real(kind=dp) ::  o3data(npoints,nlayers21)

   real(kind=dp) ::  o3_in (nlayers21,npoints)
   real(kind=dp) ::  o3_flt(nlayers21,npoints)

!  Other help arrays

   integer       :: order(nlayers21)
   real(kind=dp) :: lambda(nlayers21), evec_2(nlayers21,nlayers21)
   real(kind=dp) :: KSQ_ordered(nlayers21), ksq_abs(nlayers21)

!  tolerance input to Eigenpackage module ASYMTX

   real(kind=dp) :: tol

!  (output from Eigenpackage module ASYMTX)

   real(kind=dp) :: KSQ(nlayers21), WK(2*nlayers21)
   real(kind=dp) :: eigenmat(nlayers21,nlayers21)
   real(kind=dp) :: evec(nlayers21,nlayers21)
   INTEGER       :: IER
   LOGICAL(c_bool) :: ASYMTX_FAILURE

!  Help variables

   character*3   :: ci
   integer       :: i,k,i1,w,aa,it,nlayers1,nlayers42
   real(kind=dp) :: stdv, norm, ddim
   real(kind=dp) :: eigenmat_save(nlayers21,nlayers21)

!  initialize exception handling

   fail = .false. ; message = ' ' ; trace = ' '

!  setups

   ddim = 1.0_dp / real(npoints,dp)
   nlayers1  = nlayers + 1
   nlayers42 = 2*nlayers21

!  1. Get & process data
!     ==================

!  process data into one array, take logarithm

   do it = 1,npoints
      data(it,1:nlayers)         = taudp(1:nlayers,it)
      data(it,nlayers1:nlayers2) = omega(1:nlayers,it)
      data(it,nlayers21)         = albedo(it)
   enddo
   o3data(:,:) = log(data(:,:))

!  compute mean values for output

   do i = 1,nlayers
      atmosmean(i,1) = sum(o3data(:,i))*ddim
   enddo
   do i = nlayers+1,nlayers2
      atmosmean(i-nlayers,2) = sum(o3data(:,i))*ddim
   enddo
   albmean = sum(o3data(:,nlayers21))*ddim

!  Transpose Log-data

   o3_in = transpose(o3data)

!  Remove mean

   do i = 1,nlayers21
      o3_flt(i,:) = o3_in(i,:)-sum(o3_in(i,:))*ddim
   enddo

!  Sanity Check
!   do i = 1, nlayers21
!      write(65,'(1p10e16.8)')(o3_flt(i,k),k=1,npoints,10)
!   enddo

!  2. Prepare Eigenmatrix and Solve Eigenproblem
!     ==========================================

!  Prepare matrix (cross covariances)

   call Prepare_Eigenmatrix(nlayers21,nlayers21,npoints,o3_flt,o3_flt,eigenmat)

!  Save eigenmat in case tol needs to be changed

   eigenmat_save = eigenmat

!  Sanity Check
!   do i = 1, nlayers21
!      write(56,'(1p10e20.10)')(eigenmat(i,k),k=1,6)
!   enddo
!   pause'Eigennmat'

!  Solve eigenproblem using PCA_ASYMTX

   tol = 1.0d-6

   CALL PCA_ASYMTX ( eigenmat, nlayers21,  nlayers21,  nlayers21, nlayers42, tol, &
                     EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )

!  Change tolerance and rerun PCA_ASYMTX if eigenvalue has not converged

   DO WHILE ( IER.GT.0 )
      eigenmat = eigenmat_save
      tol = tol * 10.d0
      CALL PCA_ASYMTX ( eigenmat, nlayers21,  nlayers21,  nlayers21, nlayers42, tol, &
                        EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )
   ENDDO

!  Exception handling 1

   IF ( ASYMTX_FAILURE  ) THEN
     TRACE   = 'ASYMTX error in pca_eigensolver_Alb'
     fail = .true. ; RETURN
   ENDIF

!  Exception handling 2

   IF ( IER.GT.0 ) THEN
      WRITE(CI,'(I3)')IER
      MESSAGE = 'Eigenvalue '//CI//' has not converged'
      TRACE   = 'ASYMTX error in pca_eigensolver_Alb'
      fail = .true. ; RETURN
   ENDIF

!  Normalize vectors

   do i = 1, nlayers21
      norm   = sum( evec(:,i)*evec(:,i) )
      do k = 1, nlayers21
         evec(k,i)= evec(k,i) / norm
      enddo
   enddo

!  rank absolute(eigenvalues). [subroutine "Ranker" is same as Indexx1]

   do i = 1, nlayers21
      ksq_abs(i) = dabs(ksq(i))
   enddo
   call PCA_Ranker(nlayers21,ksq_abs,order)

!  Rank the Eigenvectors

   do i = 1, nlayers21
      i1 = nlayers21 + 1 - i
      ksq_ordered(i) = ksq_abs(order(i1))
      do k = 1, nlayers21
          evec_2(k,i)= evec(k,order(i1))
      enddo
   enddo

!  3. Prepare EOF and PC output 
!     =========================

!    *** Only perform for the first few EOFs

   do aa = 1, N_Eofs

!  set the usable eigenvalue

      lambda(aa) = ksq_ordered(aa)

!  Normalize everything according to SQRT(Lambda)

      stdv = sqrt(lambda(aa))

!  EOFs (Unnormalized) --> Transpose the Eigenvectors

      eofs(aa,:) = evec_2(:,aa)

!  Project data onto E1 basis
!    -- Set the principal components (unnormalized)

      do w = 1, npoints
         PrinComps(aa,w) = sum(eofs(aa,:)*o3_flt(:,w))
      enddo

!  Sanity Check
!      write(*,*)lambda(aa)
!      write(*,*)eofs(aa,3)
!      write(*,*)PrinComps(aa,22)

!  Final normalization of EOFs and PCs

      eofs(aa,:) = eofs(aa,:) * stdv
      PrinComps(aa,:) = PrinComps(aa,:)/stdv

!  End User-defined EOF loop

   enddo

!  Finish

   return
end subroutine pca_eigensolver_alb

subroutine pca_eigensolver &
           ( Max_Eofs, maxpoints, maxlayers,  maxlayers2, & ! Input dimensions
             n_Eofs, npoints, nlayers, nlayers2,          & ! Input control
             taudp, omega,                                & ! Input optical properties
             Atmosmean, Eofs, PrinComps,                  & ! Outputs
             fail, message_len, message, trace_len, trace) bind(C) ! Exception handling

!  PCA using Eigenproblem methods

   implicit none

!  inputs
!  ------

!  Dimensioning

   integer(c_int), intent(in)       :: Max_Eofs, maxpoints, maxlayers, maxlayers2

!  Control

   integer(c_int), intent(in)       :: n_Eofs, npoints, nlayers, nlayers2

!  Optical properties

   real(kind=c_double), intent(in) :: taudp(maxlayers,maxpoints)
   real(kind=c_double), intent(in) :: omega(maxlayers,maxpoints)

!  outputs
!  -------

!  Mean values

   real(kind=c_double), intent(out) :: atmosmean(maxlayers,2)

!  EOFS and Principal Components (Latter is allocatable)

   real(kind=c_double), intent(out) :: eofs  (Max_Eofs,maxlayers2)
   real(kind=c_double), intent(out) :: PrinComps(n_Eofs,npoints)

!  Exception handling and status

   logical(c_bool), intent(out)       :: fail
   integer(c_int), intent(in) :: message_len, trace_len
   character(kind=c_char), intent(out) :: message(message_len)
   character(kind=c_char), intent(out) :: trace(trace_len)

!  local data arrays
!  -----------------

   real(kind=dp) ::  data  (npoints,nlayers2)
   real(kind=dp) ::  o3data(npoints,nlayers2)

   real(kind=dp) ::  o3_in (nlayers2,npoints)
   real(kind=dp) ::  o3_flt(nlayers2,npoints)

!  Other help arrays

   integer       :: order(nlayers2)
   real(kind=dp) :: lambda(nlayers2), evec_2(nlayers2,nlayers2)
   real(kind=dp) :: KSQ_ordered(nlayers2), ksq_abs(nlayers2)

!  tolerance input to Eigenpackage module ASYMTX

   real(kind=dp) :: tol

!  (output from Eigenpackage module ASYMTX)

   real(kind=dp) :: KSQ(nlayers2), WK(2*nlayers2)
   real(kind=dp) :: eigenmat(nlayers2,nlayers2)
   real(kind=dp) :: evec(nlayers2,nlayers2)
   INTEGER       :: IER
   LOGICAL(c_bool) :: ASYMTX_FAILURE

!  Help variables

   character*3   :: ci
   integer       :: i,k,i1,w,aa,it,nlayers1,nlayers4
   real(kind=dp) :: stdv, norm, ddim
   real(kind=dp) :: eigenmat_save(nlayers2,nlayers2)

!  initialize exception handling

   fail = .false. ; message = ' ' ; trace = ' '

!  Setups

   ddim = 1.0_dp / real(npoints,dp)
   nlayers1 = nlayers + 1
   nlayers4 = 2*nlayers2

!  1. Get & process data
!     ==================

!  process data into one array, take logarithm

   do it = 1,npoints
      data(it,1:nlayers) = taudp(1:nlayers,it)
      data(it,nlayers1:nlayers2) = omega(1:nlayers,it)
   enddo
   o3data(:,:) = log(data(:,:))

!  compute mean value for output

   do i = 1,nlayers
      Atmosmean(i,1) = sum(o3data(:,i))*ddim
   enddo
   do i = nlayers1,nlayers2
      Atmosmean(i-nlayers,2) = sum(o3data(:,i))*ddim
   enddo

!  Transpose

   o3_in = transpose(o3data)

!  Remove time-mean

   do i = 1,nlayers2
      o3_flt(i,:) = o3_in(i,:)-sum(o3_in(i,:))*ddim
   enddo

!  Sanity Check
!   do i = 1, nlayers2
!      write(65,'(1p10e16.8)')(o3_flt(i,k),k=1,npoints,20)
!   enddo

!  2. Prepare Eigenmatrix and Solve Eigenproblem
!     ==========================================

!  Prepare matrix (cross covariances)

   call Prepare_Eigenmatrix(nlayers2,nlayers2,npoints,o3_flt,o3_flt,eigenmat)

!  Save eigenmat in case tol needs to be changed

   eigenmat_save = eigenmat

!  Debug
!    do i = 1, 26
!       do j = 1, 26
!           write(677,*),i,j,eigenmat(i,j)
!       enddo
!    enddo

!  Solve eigenproblem using PCA_ASYMTX

   tol = 1.0d-6

   CALL PCA_ASYMTX ( eigenmat, nlayers2,  nlayers2,  nlayers2, 2*nlayers2, tol, &
                     EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )

!  Change tolerance and rerun PCA_ASYMTX if eigenvalue has not converged

   DO WHILE ( IER.GT.0 )
      eigenmat = eigenmat_save
      tol = tol * 10.d0
      CALL PCA_ASYMTX ( eigenmat, nlayers2,  nlayers2,  nlayers2, 2*nlayers2, tol, &
                        EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )
   ENDDO

!  Exception handling 1

   IF ( ASYMTX_FAILURE  ) THEN
     TRACE   = 'ASYMTX error in pca_eigensolver'
     fail = .true. ; RETURN
   ENDIF

!  Exception handling 2

   IF ( IER.GT.0 ) THEN
      WRITE(CI,'(I3)')IER
      MESSAGE = 'Eigenvalue '//CI//' has not converged'
      TRACE   = 'ASYMTX error in pca_eigensolver'
      fail = .true. ; RETURN
   ENDIF

!  Normalize vectors

   do i = 1, nlayers2
      norm   = sum( evec(:,i)*evec(:,i) )
      do k = 1, nlayers2
         evec(k,i)= evec(k,i) / norm
      enddo
   enddo

!  rank absolute(eigenvalues). [subroutine "Ranker" is same as Indexx1]

   do i = 1, nlayers2
      ksq_abs(i) = dabs(ksq(i))
   enddo
   call PCA_Ranker(nlayers2,ksq_abs,order)

!  Rank the Eigenvectors

   do i = 1, nlayers2
      i1 = nlayers2 + 1 - i
      ksq_ordered(i) = ksq_abs(order(i1))
      do k = 1, nlayers2
          evec_2(k,i)= evec(k,order(i1))
      enddo
   enddo

!  3. Prepare EOF and PC output 
!     =========================

!    *** Only perform for the first few EOFs

   do aa = 1, N_Eofs

!  set the usable eigenvalue

      lambda(aa) = ksq_ordered(aa)

!  Normalize everything according to SQRT(Lambda)

      stdv = sqrt(lambda(aa))

!  EOFs (Unnormalized) --> Transpose the Eigenvectors

      eofs(aa,:) = evec_2(:,aa)

!  Project data onto E1 basis
!    -- Set the principal components (unnormalized)

      do w = 1, npoints
         PrinComps(aa,w) = sum(eofs(aa,:)*o3_flt(:,w))
      enddo

!  Sanity Check
!      write(*,*)aa,lambda(aa)
!      write(*,*)eofs(aa,3)
!      write(*,*)PrinComps(aa,22)

!  Final normalization of EOFs and PCs

      eofs(aa,:) = eofs(aa,:) * stdv
      PrinComps(aa,:) = PrinComps(aa,:)/stdv

!  End User-defined EOF loop

   enddo

!  Finish

   return
end subroutine pca_eigensolver

subroutine pca_eigensolver_aer &
           ( Max_Eofs, maxpoints, maxlayers, maxlayers2pa, maxaggregates, & ! Input dimensions
             n_Eofs, npoints, nlayers, nlayers2pa, naggregates,           & ! Input control
             taudp, omega, aod,                             & ! Input optical properties
             Atmosmean, aodmean, Eofs, PrinComps,           & ! Outputs
             fail, message_len, message, trace_len, trace ) bind(C)                          ! Exception handling

!  PCA using Eigenproblem methods

   implicit none

!  inputs
!  ------

!  Dimensioning

   integer(c_int), intent(in)       :: Max_Eofs, maxpoints, maxlayers, maxlayers2pa, maxaggregates

!  Control

   integer(c_int), intent(in)       :: n_Eofs, npoints, nlayers, nlayers2pa, naggregates

!  Optical properties

   real(kind=c_double), intent(in) :: taudp(maxlayers,maxpoints)
   real(kind=c_double), intent(in) :: omega(maxlayers,maxpoints)
   real(kind=c_double), intent(in) :: aod(maxaggregates,maxpoints)

!  outputs
!  -------

!  Mean values

   real(kind=c_double), intent(out) :: atmosmean(maxlayers,2), aodmean(maxaggregates)

!  EOFS and Principal Components (Latter is allocatable)

   real(kind=c_double), intent(out) :: eofs  (Max_Eofs,maxlayers2pa)
   real(kind=c_double), intent(out) :: PrinComps(n_Eofs,npoints)

!  Exception handling and status

   logical(c_bool), intent(out)       :: fail
   integer(c_int), intent(in) :: message_len, trace_len
   character(kind=c_char), intent(out) :: message(message_len)
   character(kind=c_char), intent(out) :: trace(trace_len)

!  local data arrays
!  -----------------

   real(kind=dp) ::  data  (npoints,nlayers2pa)
   real(kind=dp) ::  o3data(npoints,nlayers2pa)

   real(kind=dp) ::  o3_in (nlayers2pa,npoints)
   real(kind=dp) ::  o3_flt(nlayers2pa,npoints)

!  Other help arrays

   integer       :: order(nlayers2pa)
   real(kind=dp) :: lambda(nlayers2pa), evec_2(nlayers2pa,nlayers2pa)
   real(kind=dp) :: KSQ_ordered(nlayers2pa), ksq_abs(nlayers2pa)

!  tolerance input to Eigenpackage module ASYMTX

   real(kind=dp) :: tol

!  (output from Eigenpackage module ASYMTX)

   real(kind=dp) :: KSQ(nlayers2pa), WK(2*nlayers2pa)
   real(kind=dp) :: eigenmat(nlayers2pa,nlayers2pa)
   real(kind=dp) :: evec(nlayers2pa,nlayers2pa)
   INTEGER       :: IER
   LOGICAL(c_bool) :: ASYMTX_FAILURE

!  Help variables

   character*3   :: ci
   integer       :: i,k,i1,w,aa,it,nlayers1,nlayers2,nlayers21
   real(kind=dp) :: stdv, norm, ddim
   real(kind=dp) :: eigenmat_save(nlayers2pa,nlayers2pa)

!  initialize exception handling

   fail = .false. ; message = ' ' ; trace = ' '

!  Setups

   ddim = 1.0_dp / real(npoints,dp)
   nlayers1 = nlayers + 1
   nlayers2 = 2*nlayers
   nlayers21 = 2*nlayers + 1

!  1. Get & process data
!     ==================

!  process data into one array, take logarithm

   do it = 1,npoints
      data(it,1:nlayers) = taudp(1:nlayers,it)
      data(it,nlayers1:nlayers2) = omega(1:nlayers,it)
      data(it,nlayers21:nlayers2pa) = aod(1:naggregates,it)
   enddo
   o3data(:,:) = log(data(:,:))

!  compute mean value for output

   do i = 1,nlayers
      Atmosmean(i,1) = sum(o3data(:,i))*ddim
   enddo
   do i = nlayers1,nlayers2
      Atmosmean(i-nlayers,2) = sum(o3data(:,i))*ddim
   enddo
   do i = 1,naggregates
      aodmean(i) = sum(o3data(:,nlayers2+i))*ddim
   enddo

!  Transpose

   o3_in = transpose(o3data)

!  Remove time-mean

   do i = 1,nlayers2pa
      o3_flt(i,:) = o3_in(i,:)-sum(o3_in(i,:))*ddim
   enddo

!  Sanity Check
!   do i = 1, nlayers2pa
!      write(65,'(1p10e16.8)')(o3_flt(i,k),k=1,npoints,20)
!   enddo

!  2. Prepare Eigenmatrix and Solve Eigenproblem
!     ==========================================

!  Prepare matrix (cross covariances)

   call Prepare_Eigenmatrix(nlayers2pa,nlayers2pa,npoints,o3_flt,o3_flt,eigenmat)

!  Save eigenmat in case tol needs to be changed

   eigenmat_save = eigenmat

!  Debug
!    do i = 1, 26
!       do j = 1, 26
!           write(677,*),i,j,eigenmat(i,j)
!       enddo
!    enddo

!  Solve eigenproblem using PCA_ASYMTX

   tol = 1.0d-6

   CALL PCA_ASYMTX ( eigenmat, nlayers2pa,  nlayers2pa,  nlayers2pa, 2*nlayers2pa, tol, &
                     EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )

!  Change tolerance and rerun PCA_ASYMTX if eigenvalue has not converged

   DO WHILE ( IER.GT.0 )
      eigenmat = eigenmat_save
      tol = tol * 10.d0
      CALL PCA_ASYMTX ( eigenmat, nlayers2pa,  nlayers2pa,  nlayers2pa, 2*nlayers2pa, tol, &
                        EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )
   ENDDO

!  Exception handling 1

   IF ( ASYMTX_FAILURE  ) THEN
     TRACE   = 'ASYMTX error in pca_eigensolver'
     fail = .true. ; RETURN
   ENDIF

!  Exception handling 2

   IF ( IER.GT.0 ) THEN
      WRITE(CI,'(I3)')IER
      MESSAGE = 'Eigenvalue '//CI//' has not converged'
      TRACE   = 'ASYMTX error in pca_eigensolver'
      fail = .true. ; RETURN
   ENDIF

!  Normalize vectors

   do i = 1, nlayers2pa
      norm   = sum( evec(:,i)*evec(:,i) )
      do k = 1, nlayers2pa
         evec(k,i)= evec(k,i) / norm
      enddo
   enddo

!  rank absolute(eigenvalues). [subroutine "Ranker" is same as Indexx1]

   do i = 1, nlayers2pa
      ksq_abs(i) = dabs(ksq(i))
   enddo
   call PCA_Ranker(nlayers2pa,ksq_abs,order)

!  Rank the Eigenvectors

   do i = 1, nlayers2pa
      i1 = nlayers2pa + 1 - i
      ksq_ordered(i) = ksq_abs(order(i1))
      do k = 1, nlayers2pa
          evec_2(k,i)= evec(k,order(i1))
      enddo
   enddo

!  3. Prepare EOF and PC output 
!     =========================

!    *** Only perform for the first few EOFs

   do aa = 1, N_Eofs

!  set the usable eigenvalue

      lambda(aa) = ksq_ordered(aa)

!  Normalize everything according to SQRT(Lambda)

      stdv = sqrt(lambda(aa))

!  EOFs (Unnormalized) --> Transpose the Eigenvectors

      eofs(aa,:) = evec_2(:,aa)

!  Project data onto E1 basis
!    -- Set the principal components (unnormalized)

      do w = 1, npoints
         PrinComps(aa,w) = sum(eofs(aa,:)*o3_flt(:,w))
      enddo

!  Sanity Check
!      write(*,*)aa,lambda(aa)
!      write(*,*)eofs(aa,3)
!      write(*,*)PrinComps(aa,22)

!  Final normalization of EOFs and PCs

      eofs(aa,:) = eofs(aa,:) * stdv
      PrinComps(aa,:) = PrinComps(aa,:)/stdv

!  End User-defined EOF loop

   enddo

!  Finish

   return
end subroutine pca_eigensolver_aer

subroutine pca_eigensolver_Continuum_alb &
           ( Max_Eofs, maxpoints, maxlayers, maxlayers1,   & ! Input dimensions
             n_Eofs, npoints, nlayers, nlayers1,           & ! Input control
             taudp, albedo,                                & ! Input optical properties
             Atmosmean, Albmean, Eofs, PrinComps,          & ! Outputs
             fail, message_len, message, trace_len, trace ) bind(C) ! Exception handling

!  PCA using Eigenproblem methods

   implicit none

!  inputs
!  ------

!  Dimensioning

   integer(c_int), intent(in)       :: Max_Eofs, maxpoints, maxlayers, maxlayers1

!  Control

   integer(c_int), intent(in)       :: n_Eofs, npoints, nlayers, nlayers1

!  Optical properties

   real(kind=c_double), intent(in) :: taudp(maxlayers,maxpoints)
   real(kind=c_double), intent(in) :: albedo(maxpoints)

!  outputs
!  -------

!  Mean values

   real(kind=c_double), intent(out) :: atmosmean(maxlayers),albmean

!  EOFS and Principal Components (Latter is allocatable)

   real(kind=c_double), intent(out) :: eofs  (Max_Eofs,maxlayers1)
   real(kind=c_double), intent(out) :: PrinComps(n_eofs,npoints)

!  Exception handling and status

   logical(c_bool), intent(out)       :: fail
   integer(c_int), intent(in) :: message_len, trace_len
   character(kind=c_char), intent(out) :: message(message_len)
   character(kind=c_char), intent(out) :: trace(trace_len)

!  local data arrays
!  -----------------

   real(kind=dp) ::  data  (npoints,nlayers1)
   real(kind=dp) ::  o3data(npoints,nlayers1)

   real(kind=dp) ::  o3_in (nlayers1,npoints)
   real(kind=dp) ::  o3_flt(nlayers1,npoints)

!  Other help arrays

   integer       :: order(nlayers1)
   real(kind=dp) :: lambda(nlayers1), evec_2(nlayers1,nlayers1)
   real(kind=dp) :: KSQ_ordered(nlayers1), ksq_abs(nlayers1)

!  tolerance input to Eigenpackage module ASYMTX

   real(kind=dp) :: tol

!  (output from Eigenpackage module ASYMTX)

   real(kind=dp) :: KSQ(nlayers1), WK(2*nlayers1)
   real(kind=dp) :: eigenmat(nlayers1,nlayers1)
   real(kind=dp) :: evec(nlayers1,nlayers1)
   INTEGER       :: IER
   LOGICAL(c_bool) :: ASYMTX_FAILURE

!  Help variables

   character*3   :: ci
   integer       :: i,k,i1,w,aa,it,nlayers22
   real(kind=dp) :: stdv, norm, ddim
   real(kind=dp) :: eigenmat_save(nlayers1,nlayers1)

!  initialize exception handling

   fail = .false. ; message = ' ' ; trace = ' '

!  setups

   ddim = 1.0_dp / real(npoints,dp)

!  1. Get & process data
!     ==================

!  process data into one array, take logarithm

   do it = 1,npoints
      data(it,1:nlayers)         = taudp(1:nlayers,it)
      data(it,nlayers1)         = albedo(it)
   enddo
   o3data(:,:) = log(data(:,:))
   nlayers22 = 2*nlayers1

!  compute mean values for output

   do i = 1,nlayers
      atmosmean(i) = sum(o3data(:,i))*ddim
   enddo
   albmean = sum(o3data(:,nlayers1))*ddim

!  Transpose Log-data

   o3_in = transpose(o3data)

!  Remove mean

   do i = 1,nlayers1
      o3_flt(i,:) = o3_in(i,:)-sum(o3_in(i,:))*ddim
   enddo

!  Sanity Check
!   do i = 1, nlayers1
!      write(65,'(1p10e16.8)')(o3_flt(i,k),k=1,npoints,10)
!   enddo

!  2. Prepare Eigenmatrix and Solve Eigenproblem
!     ==========================================

!  Prepare matrix (cross covariances)

   call Prepare_Eigenmatrix(nlayers1,nlayers1,npoints,o3_flt,o3_flt,eigenmat)

!  Save eigenmat in case tol needs to be changed

   eigenmat_save = eigenmat

!  Sanity Check
!   do i = 1, nlayers1
!      write(56,'(1p10e20.10)')(eigenmat(i,k),k=1,6)
!   enddo
!   pause'Eigennmat'

!  Solve eigenproblem using PCA_ASYMTX

   tol = 1.0d-6

   CALL PCA_ASYMTX ( eigenmat, nlayers1,  nlayers1,  nlayers1, nlayers22, tol, &
                     EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )

!  Change tolerance and rerun PCA_ASYMTX if eigenvalue has not converged

   DO WHILE ( IER.GT.0 )
      eigenmat = eigenmat_save
      tol = tol * 10.d0
      CALL PCA_ASYMTX ( eigenmat, nlayers1,  nlayers1,  nlayers1, nlayers22, tol, &
                        EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )
   ENDDO

!  Exception handling 1

   IF ( ASYMTX_FAILURE  ) THEN
     TRACE   = 'ASYMTX error in pca_eigensolver_Alb'
     fail = .true. ; RETURN
   ENDIF

!  Exception handling 2

   IF ( IER.GT.0 ) THEN
      WRITE(CI,'(I3)')IER
      MESSAGE = 'Eigenvalue '//CI//' has not converged'
      TRACE   = 'ASYMTX error in pca_eigensolver_Alb'
      fail = .true. ; RETURN
   ENDIF

!  Normalize vectors

   do i = 1, nlayers1
      norm   = sum( evec(:,i)*evec(:,i) )
      do k = 1, nlayers1
         evec(k,i)= evec(k,i) / norm
      enddo
   enddo

!  rank absolute(eigenvalues). [subroutine "Ranker" is same as Indexx1]

   do i = 1, nlayers1
      ksq_abs(i) = dabs(ksq(i))
   enddo
   call PCA_Ranker(nlayers1,ksq_abs,order)

!  Rank the Eigenvectors

   do i = 1, nlayers1
      i1 = nlayers1 + 1 - i
      ksq_ordered(i) = ksq_abs(order(i1))
      do k = 1, nlayers1
          evec_2(k,i)= evec(k,order(i1))
      enddo
   enddo

!  3. Prepare EOF and PC output 
!     =========================

!    *** Only perform for the first few EOFs

   do aa = 1, N_Eofs

!  set the usable eigenvalue

      lambda(aa) = ksq_ordered(aa)

!  Normalize everything according to SQRT(Lambda)

      stdv = sqrt(lambda(aa))

!  EOFs (Unnormalized) --> Transpose the Eigenvectors

      eofs(aa,:) = evec_2(:,aa)

!  Project data onto E1 basis
!    -- Set the principal components (unnormalized)

      do w = 1, npoints
         PrinComps(aa,w) = sum(eofs(aa,:)*o3_flt(:,w))
      enddo

!  Sanity Check
!      write(*,*)lambda(aa)
!      write(*,*)eofs(aa,3)
!      write(*,*)PrinComps(aa,22)

!  Final normalization of EOFs and PCs

      eofs(aa,:) = eofs(aa,:) * stdv
      PrinComps(aa,:) = PrinComps(aa,:)/stdv

!  End User-defined EOF loop

   enddo

!  Finish

   return
end subroutine pca_eigensolver_Continuum_alb

subroutine pca_eigensolver_Continuum &
           ( Max_Eofs, maxpoints, maxlayers,   & ! Input dimensions
             n_Eofs, npoints, nlayers, taudp,  & ! Input control and optical
             Atmosmean, Eofs, PrinComps,       & ! Outputs
             fail, message_len, message, trace_len, trace ) ! Exception handling

!  PCA using Eigenproblem methods

   implicit none

!  inputs
!  ------

!  Dimensioning
   integer(c_int), intent(in)       :: Max_Eofs, maxpoints, maxlayers

!  Control

   integer(c_int), intent(in)       :: n_Eofs, npoints, nlayers

!  Optical properties

   real(kind=c_double), intent(in) :: taudp(maxlayers,maxpoints)


!  outputs
!  -------

!  Mean values

   real(kind=c_double), intent(out) :: atmosmean(maxlayers)

!  EOFS and Principal Components (Latter is allocatable)

   real(kind=c_double), intent(out) :: eofs  (Max_Eofs,maxlayers)
   real(kind=c_double), intent(out) :: PrinComps(n_Eofs,npoints)

!  Exception handling and status

   logical(c_bool), intent(out)       :: fail
   integer(c_int), intent(in) :: message_len, trace_len
   character(kind=c_char), intent(out) :: message(message_len)
   character(kind=c_char), intent(out) :: trace(trace_len)

!  local data arrays
!  -----------------

   real(kind=dp) ::  data  (npoints,nlayers)
   real(kind=dp) ::  o3data(npoints,nlayers)

   real(kind=dp) ::  o3_in (nlayers,npoints)
   real(kind=dp) ::  o3_flt(nlayers,npoints)

!  Other help arrays

   integer       :: order(nlayers)
   real(kind=dp) :: lambda(nlayers), evec_2(nlayers,nlayers)
   real(kind=dp) :: KSQ_ordered(nlayers), ksq_abs(nlayers)

!  tolerance input to Eigenpackage module ASYMTX

   real(kind=dp) :: tol

!  (output from Eigenpackage module ASYMTX)

   real(kind=dp) :: KSQ(nlayers), WK(2*nlayers)
   real(kind=dp) :: eigenmat(nlayers,nlayers)
   real(kind=dp) :: evec(nlayers,nlayers)
   INTEGER       :: IER
   LOGICAL(c_bool) :: ASYMTX_FAILURE

!  Help variables

   character*3   :: ci
   integer       :: i,k,i1,w,aa,it,nlayers2 
   real(kind=dp) :: stdv, norm, ddim
   real(kind=dp) :: eigenmat_save(nlayers,nlayers)

!  initialize exception handling

   fail = .false. ; message = ' ' ; trace = ' '

!  Setups

   ddim = 1.0_dp / real(npoints,dp)

   nlayers2 = 2*nlayers

!  1. Get & process data
!     ==================

!  process data into one array, take logarithm

   do it = 1,npoints
      data(it,1:nlayers) = taudp(1:nlayers,it)
   enddo
   o3data(:,:) = log(data(:,:))

!  compute mean value for output

   do i = 1,nlayers
      Atmosmean(i) = sum(o3data(:,i))*ddim
   enddo

!  Transpose

   o3_in = transpose(o3data)

!  Remove time-mean

   do i = 1,nlayers
      o3_flt(i,:) = o3_in(i,:)-sum(o3_in(i,:))*ddim
   enddo

!  Sanity Check
!   do i = 1, nlayers
!      write(65,'(1p10e16.8)')(o3_flt(i,k),k=1,npoints,20)
!   enddo

!  2. Prepare Eigenmatrix and Solve Eigenproblem
!     ==========================================

!  Prepare matrix (cross covariances)

   call Prepare_Eigenmatrix(nlayers,nlayers,npoints,o3_flt,o3_flt,eigenmat)

!  Save eigenmat in case tol needs to be changed

   eigenmat_save = eigenmat

!  Debug
!    do i = 1, 26
!       do j = 1, 26
!           write(677,*),i,j,eigenmat(i,j)
!       enddo
!    enddo

!  Solve eigenproblem using PCA_ASYMTX

   tol = 1.0d-6

   CALL PCA_ASYMTX ( eigenmat, nlayers,  nlayers,  nlayers, nlayers2, tol, &
                     EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )

!  Change tolerance and rerun PCA_ASYMTX if eigenvalue has not converged

   DO WHILE ( IER.GT.0 )
      eigenmat = eigenmat_save
      tol = tol * 10.d0
      CALL PCA_ASYMTX ( eigenmat, nlayers,  nlayers,  nlayers, nlayers2, tol, &
                        EVEC, KSQ, IER, WK, message_len, message, ASYMTX_FAILURE )
   ENDDO

!  Exception handling 1

   IF ( ASYMTX_FAILURE  ) THEN
     TRACE   = 'ASYMTX error in pca_eigensolver'
     fail = .true. ; RETURN
   ENDIF

!  Exception handling 2

   IF ( IER.GT.0 ) THEN
      WRITE(CI,'(I3)')IER
      MESSAGE = 'Eigenvalue '//CI//' has not converged'
      TRACE   = 'ASYMTX error in pca_eigensolver'
      fail = .true. ; RETURN
   ENDIF

!  Normalize vectors

   do i = 1, nlayers
      norm   = sum( evec(:,i)*evec(:,i) )
      do k = 1, nlayers
         evec(k,i)= evec(k,i) / norm
      enddo
   enddo

!  rank absolute(eigenvalues). [subroutine "Ranker" is same as Indexx1]

   do i = 1, nlayers
      ksq_abs(i) = dabs(ksq(i))
   enddo
   call PCA_Ranker(nlayers,ksq_abs,order)

!  Rank the Eigenvectors

   do i = 1, nlayers
      i1 = nlayers + 1 - i
      ksq_ordered(i) = ksq_abs(order(i1))
      do k = 1, nlayers
          evec_2(k,i)= evec(k,order(i1))
      enddo
   enddo

!  3. Prepare EOF and PC output 
!     =========================

!    *** Only perform for the first few EOFs

   do aa = 1, N_Eofs

!  set the usable eigenvalue

      lambda(aa) = ksq_ordered(aa)

!  Normalize everything according to SQRT(Lambda)

      stdv = sqrt(lambda(aa))

!  EOFs (Unnormalized) --> Transpose the Eigenvectors

      eofs(aa,:) = evec_2(:,aa)

!  Project data onto E1 basis
!    -- Set the principal components (unnormalized)

      do w = 1, npoints
         PrinComps(aa,w) = sum(eofs(aa,:)*o3_flt(:,w))
      enddo

!  Sanity Check
!      write(*,*)aa,lambda(aa)
!      write(*,*)eofs(aa,3)
!      write(*,*)PrinComps(aa,22)

!  Final normalization of EOFs and PCs

      eofs(aa,:) = eofs(aa,:) * stdv
      PrinComps(aa,:) = PrinComps(aa,:)/stdv

!  End User-defined EOF loop

   enddo

!  Finish

   return
end subroutine pca_eigensolver_Continuum

subroutine Prepare_Eigenmatrix(m,n,t,x,y,ccm) bind(C)

!  Computes the cross covariance matrix of two vector time series x,y
!  x(m,t) y(n,t) and returns the m x n matrix

   implicit none

!  inputs

   integer(c_int):: m,n,t
   real(kind=c_double) :: x(m,t),y(n,t)

!  outputs

   real(kind=c_double) ::ccm(m,n)

!  local variables

   integer       :: i,j
   real(kind=dp) :: x1(t),y1(t),x1bar,y1bar, dt
   real(kind=dp), parameter :: small = 1.0e-20_dp

   dt = 1.0_dp / real(t,dp)
   do i = 1,m
      do j = 1,n
         x1(:) = x(i,:)
         y1(:) = y(j,:)
         x1bar = sum(x1)*dt
         y1bar = sum(y1)*dt
         if ( (sum(dabs(x1)).gt.small) .and. &
              (sum(dabs(y1)).gt.small) ) then
            ccm(i,j) = sum((x1(:)-x1bar)*(y1(:)-y1bar))*dt
         else
            ccm(i,j) = 0.0_dp
         endif
      enddo
   enddo

!  Done

   return
end subroutine Prepare_Eigenmatrix

!  End module

end module pca_eigensolvers_m
