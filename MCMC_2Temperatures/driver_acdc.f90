module driver
use acdc_system, only : neqn => nclust					! number of equations
use acdc_system, only : temp							! temperature

implicit none

contains

subroutine acdc_driver(c0,c,ok,theta, conv)

	implicit none
	real(kind(1.d0)), intent(out) :: c(neqn)				! concentrations
	real(kind(1.d0)) :: ci(neqn)				! concentrations
	real(kind(1.d0)) :: c0(neqn), conv				! initial concentrations
	integer, parameter :: nt = 7						! number of timeranges for integration
	real(kind(1.d0)) :: t0, t(nt), tt(nt)				! start and end times and initial time step
!	real(kind(1.d0)) :: jtot(npoints)					! simulated formation rates
	integer :: ipar(4)									! parameters related to the loop over concentrations
	logical :: ok	
    real(kind(1.d0)):: coefs_dHdS(44), sources(2)	
	integer :: I_coll(37),J_coll(37),IJ_evap(37),I_evap(37),J_evap(37),N_evap, k
  real :: c_nmer_cmp(neqn), c_nmer_end(neqn)
	! .false. if solver failed
	
	! parameters for the solver
	integer, parameter :: itol = 1						! same absolute tolerance for all components
	integer, parameter :: itask = 3					! stop at or beyond t = TOUT and return
	integer, parameter :: iopt = 0						! allow using optional input in the solver
	integer, parameter :: mf = 22						! full jacobian, computed numerically
	integer, parameter :: kB = 1.3806504d-23   ! Boltzmann constant
	integer, parameter :: pres_atm = 101325.d0 ! atmospheric pressure
 	real(kind(1.d0)), parameter :: t_max = 1.08d4		! simulation length
	real(kind(1.d0)), parameter :: t_cmp = 1.02d4		! earlier time-point used for checking the convergence
	real(kind(1.d0)), parameter :: rtol = 1.d-5			! relative tolerance in the solver
	real(kind(1.d0)), parameter :: atol = 1.d-5			! absolute tolerance in the solver
	integer :: lwork, liwork
	integer :: i, j, istate, iwork(neqn+30)
	real(kind(1.d0)) :: work(22+9*neqn+2*neqn**2)
	real(kind(1.d0)) :: theta(29), coef(29) ! temperature plus 39 values for the evaporation rates
	external feval, jeval					! subroutines for equations, jacobian and formation rate
	!tt(1:18) = (/ 1.d-8, 1.d-4, 1.d-2, 1.d0, 1.d2, 3.d2, 5.d2, 7.d2, 9.d2, 1.d3, 2.d3, 4.d3,  4.5d3, 8.5d3, 1.d4, 1.2d4, 1.4d4, 5.9d4/)
   tt(1:7) = (/1.d-4, 1.d-2, 1.d0, 1.d2, 1.d3, t_cmp, t_max/)
	!tt(1:11) = (/ 1.d-8, 1.d-4, 1.d-2, 1.d0, 1.d2, 1.d3, 2.d3, 2.5d3, 3.d3,  4.d3,   4.5d3 /)
   
	coef(1) =theta(1)
	coef(2:29) = theta(2:29)
	lwork = 22+9*neqn+2*neqn**2
	liwork = neqn+30
	work(5:10) = 0.d0						! use default parameter values in the solver
	iwork(5:10) = 0							! use default parameter values in the solver
	iwork(6) = 100000000						! allow more steps in the solver (probably not needed)
	
!	write(*,*) ''
!	write(*,*) 'coef: ', coef
!	jtot=0.d0
	ipar(1:3) = 0							! sources etc. have not yet been initialized (-> 1 when they have)
 
		ipar(4) = 1							! this tells sources_and_constants at which point we are, no use otherwise
		istate = 1							! tells the solver that this is the first call
											! -> 2 in the solver after succesful integration

		! run the simulation several times until a steady state is reached
		t0 = 0.d0
		t = tt
		c = c0
		do i=1,nt-1
			ci = c												! save the previous concentrations
			call DVODE (feval,neqn,c,t0,t(i),itol,rtol,atol,&	! calling the solver
			& itask,istate,iopt,work,lwork,iwork,liwork,&
			& jeval,mf,coef,ipar)
      
      t0 = t(i)  

			if (istate .ne. 2) then								! checking everything went ok - if not, print error
				write(*,*) 'ERROR: returned istate =', istate
				write(*,*) 'coef:', coef
				if (any(isnan(coef))) then
					ok = .false.
					exit 
				end if
				write(*,*) 'c: ', c
				c = ci
				istate = 1
				cycle
			end if
			if (maxval(abs(c-ci)/max(c,1d-6))<1d-6) exit		! end the simulation when a steady state is reached
			if (i .eq. nt) then
				write(*,*) "Did not converge!"
				write(*,*) coef
!				jtot = 0.d0
			end if
		end do
      
		if (istate .ne. 2) then								! checking everything went ok - if not, print error
			write(*,*) 'ERROR: returned istate =', istate
			write(*,*) 'coef:', coefs_dHdS
			ok = .false.
		endif
   
    c_nmer_cmp = c
    ci = c
    call DVODE (feval,neqn,c,t0,t(i),itol,rtol,atol,&	! calling the solver
			& itask,istate,iopt,work,lwork,iwork,liwork,&
			& jeval,mf,coef,ipar)
    
 		if (istate .ne. 2) then								! checking everything went ok - if not, print error
			write(*,*) 'ERROR: returned istate =', istate
			write(*,*) 'coef:', coefs_dHdS
			ok = .false.
		endif  
   
   c_nmer_end = c
   
 		c(:neqn) = .5d0*(c_nmer_cmp+c_nmer_end)
		i = maxloc(abs(c_nmer_cmp-c_nmer_end)/max(c_nmer_cmp,1d-6),dim=1)
		conv = c_nmer_end(i)/c_nmer_cmp(i)
   
		ipar = 0

end subroutine acdc_driver

end module driver
