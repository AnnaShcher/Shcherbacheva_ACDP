!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mcmcrun.F90
!!! Purpose: example main program for the standalone MCMC code
!!!
!!! Marko Laine <marko.laine@fmi.fi>
!!! ------------------------------------------------------------------------
!!
!! This file provides
!!
!! 1. main program that calls mcmc_main in the mcmc library
!!
!! 2. subroutine ssfunction that calculates the minus two times the log
!! likelihood function -2*log(y|theta) aka sum-of-squares function
!!
!!  function ssfunction(theta,npar,ny) result(ss)
!!    integer*4 npar, ny
!!    real*8 theta(npar)
!!    real*8 ss(ny)
!!  end function ssfunction
!!
!! The function defined below load data from file "data.dat" and then
!! calculates sum((ydata-model(theta,xdata))**2). The model is defined
!! in function modelfunction.
!!
!! 3. checkbounds function to check for out of bounds mcmc proposals
!!
!!
!! In addtion to "data.dat" we need the following files
!!
!! 1. The run time mcmc parameters in fortran namelist "mcmcinit.nml"
!!
!! 2. The following ascii .dat files
!!
!!    mcmcpar.dat  - initial value of the parameter
!!    mcmccov.dat  - initial proposal covariance matrix
!!
!!   optionally also
!!
!!    mcmcsigma2.dat - Gaussian error variance and number of observations,
!!      optional, needed id updatesigma = 1 in mcmcinit.nml
!!    mcmcnycol.dat - optional, number of values ssfunction returns (usually=1)
!! 

!!!
!!! Main program
!!! 
program mcmcmain
 use acdc_system, only : nclust  ! number of equations
 use driver, only : acdc_driver ! driver
  implicit none
  !! call mcmc_main in the library libmcmcrun.a
  real :: start_time, stop_time

  call cpu_time(start_time)
  
  !call add_nrand_noise_param(0.05d0)
  call mcmc_main()
 
 ! call add_nrand_noise_param(0.05d0)
 ! call mcmc_main()
  !call add_nrand_noise_param(0.02d0)
  !call mcmc_main()
  call cpu_time(stop_time)

  print *, "MCMC 2 million samples, computation time:", &
      stop_time - start_time, "seconds"

end program mcmcmain
!!!
!!! -2*log(p(obs|params)) for
!!! general Gaussian likelihood, data are read from file
!!!

function ssfunction(theta,npar,ny) result(ss)

  use mcmcprec
  use matutils, only : loaddata

  implicit none
  integer(kind=ik4) :: npar, ny, i
  real(kind=dbl) :: theta(npar), modoutp(ny)
  real(kind=dbl) :: ss(ny), f(714),alpha(10),ss1(1)

  real(kind=dbl), save, pointer :: data(:,:), sig(:)
  real(kind=dbl), save, pointer :: theta0(:)
  logical, save :: first = .true.

  interface
     function modelfunction(theta,x,f,data,sig) result(y)
       use mcmcprec
       implicit none
       real(kind=dbl) :: theta(:), x(:), y(size(x)),f(714)
       real(kind=dbl), pointer :: data(:,:), sig(:)
     end function modelfunction
  end interface
  

  !! during the first call, data is loaded
  if (first) then

     call loaddata('data.dat',data)
     call loaddata('mcmcsigma2.dat',sig)
     call loaddata('mcmcpar2.dat',theta0)
     
     first = .false.
     write(*,*) 'std = ', maxval(sig)
  endif

  

 ! if (checkbounds(theta)) then
  !! calculate the actual sum of squares
  !ss = sum((data(:,2) - modelfunction(theta0,data(:,1)))**2/std**2)
  !modoutp(1:size(modoutp)) = modelfunction(theta,data(:,1),f)
  
  !ss(1) = dot_product((data(:,2) - f)**2,1/std**2)
  ! ss(1) = sum(std*(data(:,2) - modelfunction(theta0,data(:,1)))**2)
   !  ss(1) = 0
   !modoutp  = sum(((data(:,2) - modelfunction(theta,data(:,1)))/std)**2)
   ! ss(1) = sum(((data(:,2) - f)/std)**2)

  !write(*,*) 'Modoutput = ', modoutp(1)   , std(1), size(data(:,2)), size(std), size(modoutp)
  !! calculate the actual sum of squares
  !alpha = (/ 0.8d0, 0.9d0, 1.d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0 /)

  ! ss = (log10(data(:,2)) - log10(modelfunction(theta,data(:,1),f)))**2
  !ss(1) = sum((log10(data(:,2)) - log10(modelfunction(theta,data(:,1),f)))**2)
  !ss(1) = sum((data(:,2) - modelfunction(theta0,data(:,1),f))**2/(0.001*data(:,2))**2)
  modoutp(1:size(modoutp)) = modelfunction(theta,data(:,1),f,data,sig)
  
  ss(1) = 0
  do i = 1,size(data(:,2))
    if(mod(i,17) .EQ. 0) then
     ss(1) = ss(1) + (f(i) - 1)**2
    else
     ss(1) = ss(1) + (data(i,2) - f(i))**2/(0.001*data(i,2))**2
    end if
  end do
  
 
 ! write(*,*) 'SS2 = ', ss(1), size(data(:,1)),  modoutp(1), f(2)
 ! open(27,file = 'modeloutput.dat')
 ! write(27,*),f
 ! close(27)

  !else 
	 ! ss(1) = 1.d0/0.d0
  !endif

end function ssfunction

!!!
!!! model function used by ssfunction
!!!
function modelfunction(theta,x,f,data,sig) result(y)

  use mcmcprec
  use driver
  use acdc_system, only : nclust, temp  ! number of equations
  implicit none
  
  real(kind=dbl) :: theta(:), x(:), y(size(x)), f(714)
  real(kind(1.d0)), parameter :: kB = 1.3806504d-23   ! Boltzmann constant
  real(kind(1.d0)), parameter :: pres_atm = 101325d0 ! atmospheric pressure
  integer :: i,j, k, m
  logical :: solver_ok
  real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(7), Cb_vector(3) , Temperature(2)
  !real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(1), Cb_vector(1) 
  real(kind(1.d0)) :: conc_final(nclust), conv, coef(29) ,ss(1), sigma
  real :: start, finish
  real(kind=dbl), pointer :: data(:,:), sig(:)
  
   sigma = 0.001
  if (size(theta) .ne. 28) stop 'something wrong with sizes in modelfunction'

  Ca_vector = 1.d6*(/1.d7, 21544346.9003188d0,	46415888.3361278d0,	100000000d0,	215443469.003189d0,	464158883.361277d0,	5.d7/)		!mol/m3
 ! Ca_vector = 1.d6*(/1.d7 /)		!mol/m3
  Cb_vector = (/1d0, 5d0, 100d0/) !ppt
  !Cb_vector = (/1d0/) !ppt		
  Temperature = (/278.d0, 292.d0/)  
  coef(2:29) = theta(1:28)
    m = 0
	do k = 1, size(Temperature)
		do i = 1, size(Ca_vector) ! loop over all the ambient conditions
			do j = 1, size(Cb_vector)
				initial_concentrations = 0.d0
				coef(1) = Temperature(k)
				initial_concentrations(1) = Ca_vector(i)
				initial_concentrations(3) = 1.d6*Cb_vector(j)*1.d-18*pres_atm/kB/temp
				call acdc_driver(initial_concentrations,conc_final,solver_ok, coef, conv) !This is what gets timed
       ! CALL add_nrand_noise(conc_final,sigma)
				y(m+1:m+nclust) = conc_final(1:nclust)
				y(m+nclust+1) =  conv
				m = m + nclust+1
			end do
		end do
	end do		
	
  !write(*,*) 'modeloutput'
  !write(*,*) y
  f(1:714) = y(1:714)
   
  !open(27,file = 'data.dat',status='replace')
  !write (27,*) f
  !do k=1,size(f)
   ! write (27,*) k, f(k)
  !end do
  !close(27)
end function modelfunction


!!!
!!! this function returns false if any theta(i) is out of bounds
!!!
!!!
!!! function for checking parameter bounds
!!!



subroutine add_nrand_noise(y,sigma)
 implicit none
 integer, parameter :: nclust = 16
 real(kind(1.d0)) :: sigma, delta(nclust), w(nclust,2)
 real(kind(1.d0)), intent(out) :: y(nclust)
 real(kind(1.d0)), parameter :: pi = 4.d0*atan(1.d0)
 
 !write(*,*) 'Function invoked'
 !CALL init_random_seed()         ! see example of RANDOM_SEED
 CALL random_number(w)
 delta = sqrt(-2.d0*log(w(:,1)))*sin(2.d0*pi*w(:,2))*sigma*y
 y = y+delta
end subroutine add_nrand_noise	

subroutine add_nrand_noise_param(sigma)
 implicit none
 integer, parameter :: nparam = 28
 integer :: i
 real(kind(1.d0)) :: sigma, delta(nparam), w(nparam,2)
 real(kind(1.d0)) :: y(nparam)
 logical :: checkbounds_y
 real(kind(1.d0)), parameter :: pi = 4.d0*atan(1.d0)
 
 
 open(26,file = 'mcmcpar2.dat')
 read(26,*) y
 close(26)
 
 !write(*,*) 'Function invoked'
 !CALL init_random_seed()         ! see example of RANDOM_SEED
 CALL random_number(w)
 delta = sqrt(-2.d0*log(w(:,1)))*sin(2.d0*pi*w(:,2))*sigma*abs(y)
 y = y+delta
 !write(*,*) 'Function invoked'
 !write(*,*) y
    checkbounds_y  = .TRUE.

 if (any(y(1:2:27)>=0.0) .OR. any(y < -400.d0) &
   & .OR. y(3) < y(5)  .OR. y(5) < y(7) .OR. y(7) < y(11) .OR. y(11) < y(13) .OR. y(13) < y(17 ) .OR. y(17) < y(21) .OR. &
   & y(15) < y(21) .OR. y(15) < y(25) .OR. y(15) < y(23) .OR. y(25) < y(27) .OR. y(23) < y(27) .OR. y(9) <y(11)) checkbounds_y =  .FALSE.
 do while(.NOT. checkbounds_y)
  CALL random_number(w)
  delta = sqrt(-2.d0*log(w(:,1)))*sin(2.d0*pi*w(:,2))*sigma*abs(y)
  y = y+delta
  checkbounds_y  = .TRUE.
  if (any(y(1:2:27)>=0.0) .OR. any(y < -400.d0) &
     & .OR. y(3) < y(5)  .OR. y(5) < y(7) .OR. y(7) < y(11) .OR. y(11) < y(13) .OR. y(13) < y(17 ) .OR. y(17) < y(21) .OR. &
     & y(15) < y(21) .OR. y(15) < y(25) .OR. y(15) < y(23) .OR. y(25) < y(27) .OR. y(23) < y(27) .OR. y(9) <y(11)) checkbounds_y =  .FALSE.
 end do
 
 open(27,file = 'mcmcpar.dat')
 do i = 1,nparam
   write (27,*) y(i)
 end do 
 close(27)
 
end subroutine add_nrand_noise_param	

function checkbounds(theta) 
   implicit none
   
   real*8 theta(:)
   logical checkbounds
   
   checkbounds  = .TRUE.

   if (any(theta(1:2:27)>=0.0) .OR. any(theta < -400.d0) &
   & .OR. theta(3) < theta(5)  .OR. theta(5) < theta(7) .OR. theta(7) < theta(11) .OR. theta(11) < theta(13) .OR. theta(13) < theta(17 ) .OR. theta(17) < theta(21) .OR. &
   & theta(15) < theta(21) .OR. theta(15) < theta(25) .OR. theta(15) < theta(23) .OR. theta(25) < theta(27) .OR. theta(23) < theta(27) .OR. theta(9) <theta(11)) checkbounds =  .FALSE.
   return
end function checkbounds

subroutine checkbounds2(theta,checkbounds) 
   implicit none
   
   real*8 theta(:)
   logical :: checkbounds
   
   checkbounds  = .TRUE.

   if (any(theta(1:2:27)>=0.0) .OR. any(theta < -400.d0) &
   & .OR. theta(3) < theta(5)  .OR. theta(5) < theta(7) .OR. theta(7) < theta(11) .OR. theta(11) < theta(13) .OR. theta(13) < theta(17 ) .OR. theta(17) < theta(21) .OR. &
   & theta(15) < theta(21) .OR. theta(15) < theta(25) .OR. theta(15) < theta(23) .OR. theta(25) < theta(27) .OR. theta(23) < theta(27) .OR. theta(9) <theta(11)) checkbounds =  .FALSE.
   return
end subroutine  checkbounds2
