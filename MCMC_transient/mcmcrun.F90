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
 use acdc_system, only : nclust, temp  ! number of equations
 use driver, only : acdc_driver ! driver
  implicit none
  real :: start, finish
  call cpu_time(start)
  !! call mcmc_main in the library libmcmcrun.a
  !call generate_data() 
  call mcmc_main()
  call cpu_time(finish)
   print *, "Setup time:", &
      finish - start, "seconds"
end program mcmcmain


!!!
!!! -2*log(p(obs|params)) for
!!! general Gaussian likelihood, data are read from file
!!!
function ssfunction(theta,npar,ny) result(ss)

  use mcmcprec
  use matutils, only : loaddata
  use acdc_system, only : nclust, temp  ! number of equations
  
  implicit none
  integer :: i
  integer(kind=ik4) :: npar, ny
  real(kind=dbl) :: theta(npar)
  real(kind=dbl) :: ss(ny), modeloutput(ny), f(37653)

  real(kind=dbl), save, pointer :: data(:,:),theta0(:)
  logical, save :: first = .true.

  interface
     function modelfunction(theta,x,f) result(y)
       use mcmcprec
       implicit none
       real(kind=dbl) :: theta(:), x(:), f(size(x)), y(size(x))
     end function modelfunction
  end interface

  !! during the first call, data is loaded
  if (first) then

     call loaddata('data.dat',data)
     call loaddata('mcmcpar.dat',theta0)
     first = .false.

  endif
 
  !! calculate the actual sum of squares

  !ss(1) = sum((data(:,2) - modelfunction(theta,data(:,1),f))**2/((0.001**data(:,2))**2))
  modeloutput = modelfunction(theta,data(:,1),f)
  ss(1) = 0
  do i = 1,size(data(:,2))
    if(mod(i,1793 ) .EQ. 0) then
     ss(1) = ss(1) + (f(i) - 1)**2
    else
      if(data(i,2) >1.d-8)then
       ss(1) = ss(1) + (data(i,2) - f(i))**2/(0.001*data(i,2))**2
      end if 
    end if
  end do
  !print *, 'Least square sum'
  !write(*,*) 'SS = ', ss(1), ny , size(f)
end function ssfunction

!!!
!!! model function used by ssfunction
!!!
function modelfunction(theta,x,f) result(y)

  use mcmcprec
  use driver
  use acdc_system, only : nclust, temp  ! number of equations
  implicit none
  
  integer, parameter :: nt = 114
  real(kind=dbl) :: theta(:), x(:), f(37653), y(size(x))
  real(kind(1.d0)), parameter :: kB = 1.3806504d-23   ! Boltzmann constant
  real(kind(1.d0)), parameter :: pres_atm = 101325d0 ! atmospheric pressure
  integer :: i,j, k, m, Reason
  logical :: solver_ok
  real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(7), Cb_vector(3) 
  !real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(1), Cb_vector(1) 
  real(kind(1.d0)) :: conc_final(nclust*(nt-2)), conv, coef_axil
  real :: start, finish
  real(kind(1.d0)) :: tt(nt),tt_axil(nt-2)
  logical, save :: first = .true.
  
  if (first) then
     tt(1:2) = (/1.d-8, 1.d-4/)
     Reason = 0
     j = 2 
     open(25,file="time_1_5_min_1_40h.dat", status='old', access='sequential',  action='read' )
		 DO while (Reason .GE. 0 .AND. j < nt)
		 READ(25,*,IOSTAT=Reason) coef_axil 
       j = j + 1
		   tt(j) = coef_axil
		 END DO
     close(25)
     first = .false.
  endif
  if (size(theta) .ne. 39) stop 'something wrong with sizes in modelfunction'

  Ca_vector = 1.d6*(/1.d7, 21544346.9003188d0,	46415888.3361278d0,	100000000d0,	215443469.003189d0,	464158883.361277d0,	5.d7/)		!mol/m3
 ! Ca_vector = 1.d6*(/1.d7 /)		!mol/m3
  Cb_vector = (/1d0, 5d0, 100d0/) !ppt
  !Cb_vector = (/1d0/) !ppt		
  
    m = 0
 	do i = 1, size(Ca_vector) ! loop over all the ambient conditions
		do j = 1, size(Cb_vector)
			initial_concentrations = 0.d0
			initial_concentrations(1) = Ca_vector(i)
			initial_concentrations(3) = 1.d6*Cb_vector(j)*1.d-18*pres_atm/kB/temp
			call acdc_driver(initial_concentrations,conc_final,solver_ok, theta, conv,tt) 
      !print *, 'Steady-state concentration'
     ! write(*,*) (conc_final(k), k=1,16)
		    y(m+1:m+nclust*(nt-2)) = conc_final(1:nclust)
        y(m+nclust*(nt-2)+1) = conv  
		  	m = m + nclust*(nt-2)+1
		end do
	end do  
  f(1:37653) = y(1:37653)
	   !  print *, ' Model  output successfully assigned'
end function modelfunction


!!!
!!! this function returns false if any theta(i) is out of bounds
!!!
function checkbounds(theta)
  implicit none
  real*8 theta(:)
  logical checkbounds
  
!! example: all thetas must be positive
  checkbounds = .true.
  if (any(theta<=0.d0)) checkbounds = .false.
  return

end function checkbounds

subroutine generate_data() 
  use matutils, only : loaddata
  use mcmcprec
  use driver
  use acdc_system, only : nclust, temp  ! number of equations
  implicit none
  
  integer, parameter :: nt = 114
  real(kind(1.d0)), parameter :: kB = 1.3806504d-23   ! Boltzmann constant
  real(kind(1.d0)), parameter :: pres_atm = 101325d0 ! atmospheric pressure
  integer :: i,j, k, m, Reason
  logical :: solver_ok
  real(kind=dbl) ::  y(((nt-2)*nclust+1)*21)
  real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(7), Cb_vector(3) 
  !real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(1), Cb_vector(1) 
  real(kind(1.d0)) :: conc_final(nclust*(nt-2)), conv, coef_axil, sigma
  real :: start, finish
  real(kind(1.d0)) :: tt(nt),tt_axil(nt-2), theta(39)

  
  tt(1:2) = (/1.d-8, 1.d-4/)
  Reason = 0
  j = 2 
  open(25,file="time_1_5_min_1_40h.dat", status='old', access='sequential',  action='read' )
  DO while (Reason .GE. 0 .AND. j < nt)
    READ(25,*,IOSTAT=Reason) coef_axil 
    j = j + 1
    tt(j) = coef_axil
  END DO
  close(25)
  
 
  Reason = 0
  open(25,file="mcmcpar.dat", status='old', access='sequential',  action='read' )
    READ(25,*,IOSTAT=Reason) theta
  close(25)
  
  if (size(theta) .ne. 39) stop 'something wrong with sizes in modelfunction'

  Ca_vector = 1.d6*(/1.d7, 21544346.9003188d0,	46415888.3361278d0,	100000000d0,	215443469.003189d0,	464158883.361277d0,	5.d7/)		!mol/m3
 ! Ca_vector = 1.d6*(/1.d7 /)		!mol/m3
  Cb_vector = (/1d0, 5d0, 100d0/) !ppt
  !Cb_vector = (/1d0/) !ppt		
  sigma = 0.001
    m = 0
 	do i = 1, size(Ca_vector) ! loop over all the ambient conditions
		do j = 1, size(Cb_vector)
			initial_concentrations = 0.d0
			initial_concentrations(1) = Ca_vector(i)
			initial_concentrations(3) = 1.d6*Cb_vector(j)*1.d-18*pres_atm/kB/temp
			call acdc_driver(initial_concentrations,conc_final,solver_ok, theta, conv,tt) 
      call add_nrand_noise(conc_final,sigma)
      y(m+1:m+nclust*(nt-2)) = conc_final(1:nclust*(nt-2))
      y(m+nclust*(nt-2)+1) = 1   
      m = m + nclust*(nt-2)+1
		end do
	end do  
      
	  !   print *, ' Model  output'
      !     open(27,file = 'data.dat',status='replace')
     do k =1,size(y)
       write(27,*) k, y(k)
     end do
     close(27)
end subroutine generate_data

subroutine add_nrand_noise(y,sigma)
 implicit none
 
 integer, parameter :: nt = 114, nclust = 16
 real(kind(1.d0)) :: sigma, delta(nclust*(nt-2)), w(nclust*(nt-2),2)
 real(kind(1.d0)), intent(out) :: y(nclust*(nt-2))
 real(kind(1.d0)), parameter :: pi = 4.d0*atan(1.d0)
 
 !write(*,*) 'Function invoked'
 !CALL init_random_seed()         ! see example of RANDOM_SEED
 CALL random_number(w)
 delta = sqrt(-2.d0*log(w(:,1)))*sin(2.d0*pi*w(:,2))*sigma*y
 y = y+delta
end subroutine add_nrand_noise	