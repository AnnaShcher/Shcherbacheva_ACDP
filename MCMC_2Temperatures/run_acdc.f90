program run_acdc
 use acdc_system, only : nclust, temp  ! number of equations
 use driver, only : acdc_driver ! driver
 implicit none
 !integer :: ncoefs = 44,
	real(kind(1.d0)), parameter :: kB = 1.3806504d-23   ! Boltzmann constant
 real(kind(1.d0)), parameter :: pres_atm = 101325.d0 ! atmospheric pressure
 integer :: i,j, k, m
 logical :: solver_ok
 real(kind(1.d0)) :: theta_evap(44), coef(45), Temperature(3)
 real(kind(1.d0)) :: initial_concentrations(nclust), Ca_vector(7), Cb_vector(3) 
 real(kind(1.d0)) :: conc_final(nclust), conc_output(3*21*(nclust+1)), sigma, conv
 real :: start, finish
 Ca_vector = 1.d6*(/1.d7, 21544346.9003188d0,	46415888.3361278d0,	100000000d0,	215443469.003189d0,	464158883.361277d0,	5.d7/)		!mol/m3
 !Ca_vector = 1.d6*(/1.d7, 21544346.9003188d0,	46415888.3361278d0,	100000000d0,	215443469.003189d0,	464158883.361277d0,	1.d9/)		!mol/m3
 !Cb_vector = (/1d0, 5d0, 35d0/) !ppt		
 Cb_vector = (/1d0, 5d0, 100d0/) !ppt	
 write(*,*) 'Parameters'	
 call init_dHdS(theta_evap)
 write(*,*)  (theta_evap(k), k=1,size(theta_evap))
 call cpu_time(start) !Find the time rate  
  Temperature = (/278.d0, 282.d0, 292.d0/)  
  coef(2:45) = theta_evap(1:44)
    m = 0
	do k = 1, size(Temperature)
		do i = 1, size(Ca_vector) ! loop over all the ambient conditions
			do j = 1, size(Cb_vector)
				initial_concentrations = 0.d0
				coef(1) = Temperature(k)
				initial_concentrations(1) = Ca_vector(i)
				initial_concentrations(3) = 1.d6*Cb_vector(j)*1.d-18*pres_atm/kB/temp
				call acdc_driver(initial_concentrations,conc_final,solver_ok, coef, conv) !This is what gets timed
				conc_output(m+1:m+nclust) = conc_final(1:nclust)
				conc_output(m+nclust+1) = conv    
				m = m + nclust + 1
			end do
		end do
	end do	
	
  call cpu_time(finish) ! !Stop Timer
  print '("Time = ",f6.3," seconds.")',finish-start
  
  !write(*,*) 'Output concentrations'
	!write(*,*)  (conc_output(k), k=1,size(conc_output) )
  
      write(*,*) 'Output concentrations with gaussian noise'
	write(*,*)  (conc_output(k), k=1,size(conc_output) )
 
  sigma = 1.d-4
  CALL add_nrand_noise(conc_output,sigma)
  
 
  open(24,file = 'c_steady_state_mcmc.dat') 
  rewind(24)
  do k=1,(nclust+1)*21*3
    write (24,*) k, conc_output(k)
  end do
! write(24,*) (conc_final(j),j=1,16)
  close(24) 
 

   contains

!------------------------------------------------------------------------------------------
	
	
	
	
!	subroutine get_coefs_names(ncoefs,coef_names)
!	implicit none
!	integer ncoefs
!	character(len=11), dimension(ncoefs) :: coef_names
!
 !   coef_names(1)(:) = 'dH_2A'
	!coef_names(2)(:) = 'dS_2A/dH_2A'
!	coef_names(3)(:) = 'dH_3A'
!	coef_names(4)(:) = 'dS_3A/dH_3A'
!	coef_names(5)(:) = 'dH_4A'
!	coef_names(6)(:) = 'dS_4A/dH_4A'
!	coef_names(7)(:) = 'dH_5A'
!	coef_names(8)(:) = 'dS_5A/dH_5A'
!	coef_names(9)(:) = 'dH_1A1N'
!	coef_names(10)(:) = 'dS_1A1N/dH_1A1N'
!	coef_names(11)(:) = 'dH_2A1N'
!	coef_names(12)(:) = 'dS_2A1N/dH_2A1N'
!	coef_names(13)(:) = 'dH_3A1N/dH_2A1N'
!	coef_names(14)(:) = 'dH_4A1N/dH_2A1N'
!	coef_names(15)(:) = 'dH_5A1N/dH_2A1N'
!	coef_names(16)(:) = 'dS_3A1N/dH_2A1N'
!	coef_names(17)(:) = 'dS_4A1N/dH_2A1N'
!	coef_names(18)(:) = 'dS_5A1N/dH_2A1N'
! coef_names(19)(:) = 'dH_2N'
!	coef_names(20)(:) = 'dS_2N/dH_2N'
!	coef_names(21)(:) = 'dH_2A2N'
!	coef_names(22)(:) = 'dS_2A2N/dH_2A2N'
!	coef_names(23)(:) = 'dH_3A2N/dH_2A2N'
!	coef_names(24)(:) = 'dH_4A2N/dH_2A2N'
!	coef_names(25)(:) = 'dH_5A2N/dH_2A2N'
!	coef_names(26)(:) = 'dS_3A2N/dH_2A2N'
!	coef_names(27)(:) = 'dS_4A2N/dH_2A2N'
!	coef_names(28)(:) = 'dS_5A2N/dH_2A2N'
!	coef_names(29)(:) = 'dH_3N'
!!	coef_names(30)(:) = 'dS_3N/dH_3N'
!	coef_names(31)(:) = 'dH_3A3N'
!	coef_names(32)(:) = 'dS_3A3N/dH_3A3N'
!	coef_names(33)(:) = 'dH_4A3N/dH_3A3N'
!	coef_names(34)(:) = 'dH_5A3N/dH_3A3N'
!	coef_names(35)(:) = 'dS_4A3N/dH_3A3N'
!	coef_names(36)(:) = 'dS_5A3N/dH_3A3N'
!	coef_names(37)(:) = 'dH_4N'
!	coef_names(38)(:) = 'dS_4N/dH_4N'
!	coef_names(39)(:) = 'dH_4A4N'
!	coef_names(40)(:) = 'dS_4A4N/dH_4A4N'
!	coef_names(41)(:) = 'dH_5A4N/dH_4A4N'
!	coef_names(42)(:) = 'dS_5A4N/dH_4A4N'
!	coef_names(43)(:) = 'dH_5A5N'
!	coef_names(44)(:) = 'dS_5A5N/dH_5A5N'
!
!end subroutine get_coefs_names

subroutine init_dHdS(theta_evap)
   implicit none
   integer :: j
   integer :: Reason 
   real(kind(1.d0)) :: coef_axil, theta_evap(44) 
   
   j = 1
   Reason = 0
		! Reading IJ_evap
	!	open(23,file="evap_rates.txt", status='old', access='sequential',  action='read' )
   open(23,file="theta_init.dat", status='old', access='sequential',  action='read' )
		DO while (Reason .GE. 0 .AND. j <size(theta_evap)+1)
		 READ(23,*,IOSTAT=Reason) coef_axil 
		 theta_evap(j) = coef_axil
      
      j = j + 1
		END DO
		close(23)
end subroutine init_dHdS

subroutine add_nrand_noise(y,sigma)
 implicit none
 real(kind(1.d0)) :: sigma, delta(3*21*(nclust+1)), y(3*21*(nclust+1)), w(3*21*(nclust+1),2)
 real(kind(1.d0)), parameter :: pi = 4.d0*atan(1.d0)
 
 write(*,*) 'Function invoked'
 !CALL init_random_seed()         ! see example of RANDOM_SEED
 CALL random_number(w)
 delta = sqrt(-2.d0*log(w(:,1)))*sin(2.d0*pi*w(:,2))*sigma*y
 y = y+delta
end subroutine add_nrand_noise	

end program run_acdc 

