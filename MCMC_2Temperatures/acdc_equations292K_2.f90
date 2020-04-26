!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!
! Cluster 1: 1A
! Cluster 2: 2A
! Cluster 3: 1N
! Cluster 4: 1A1N
! Cluster 5: 2A1N
! Cluster 6: 3A1N
! Cluster 7: 2A2N
! Cluster 8: 3A2N
! Cluster 9: 4A2N
! Cluster 10: 3A3N
! Cluster 11: 4A3N
! Cluster 12: 5A3N
! Cluster 13: 4A4N
! Cluster 14: 5A4N
! Cluster 15: 4A5N
! Cluster 16: 5A5N
! 17 is for source fluxes
! 18 is for coagulation losses 
! 19 is for wall losses 
! 20 is for dilution losses 
! 21 is for recombination of positive and negative charger ions with each other
! 22 is for collisions that lead succesfully out of the system as neutrals
! 23 is for collisions that lead out of the system, but the resulting cluster is brought back
!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!

! differential equations: f = dc/dt
subroutine feval(neqn,t,c,f,coef,ipar)
	implicit none
	integer, parameter :: nclust = 16
 !integer, parameter :: nclust = 16
	integer :: neqn, ipar(4), i, j,  n
	real(kind(1.d0)) :: t, c(neqn), f(neqn), coef(29), theta(28) ! coef - temperature and 39 evaporation rates 
 	real(kind(1.d0)),save ::  sources(2)
	logical, save :: isconst(16) = .false.
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,23)=0.d0,coef_lin(23,23,nclust)=0.d0,source(16)=0.d0
	integer, save :: ind_quad_loss(23,0:34)=0,ind_quad_form(23,0:122)=0
	integer, save :: ind_quad_loss_extra(23,0:30)=0,ind_quad_form_extra(23,0:72)=0
	integer, save :: ind_lin_loss(23,0:12)=0,ind_lin_form(23,0:32)=0,fitted(16,0:16)=0
	real(kind(1.d0)), save :: n_quad_form_extra(23,0:24)=0.d0, K(nclust,nclust) = 0.d0,E(nclust,nclust) = 0.d0,cs(nclust) = 0.d0
 
	real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)

    theta(1:28) = coef(2:29)
        
	! the parameters are read at the very beginning and eg. if source terms or collision rates have changed
	if (ipar(1) .eq. 0) then
		! after this the same values are used until some other routine tells otherwise
		ipar(1) = 1
		call initialize_parameters(theta,coef_quad,coef_lin, &
            & sources, ind_quad_loss,ind_quad_form, ind_lin_loss,ind_lin_form,ind_quad_loss_extra,&
            & ind_quad_form_extra,source,isconst,fitted,coef,ipar,K,E,cs)
		n_quad_form_extra(:,1:24) = real(ind_quad_form_extra(:,3:72:3),kind=kind(1.d0))
   

	end if
! Cluster 1: 1A

  
	f = 0.d0
	do i = 1, neqn ! loop over all the clusters
		if (.not. isconst(i)) then
			! first calculate coefficients for all loss terms
			do n = 1, 2*ind_quad_loss(i,0)-1, 2 ! loop over all collisions removing this cluster
				f(i) = f(i)-coef_quad(i,ind_quad_loss(i,n),ind_quad_loss(i,n+1))*c(ind_quad_loss(i,n))
			end do
			do n = 1, 2*ind_lin_loss(i,0)-1, 2 ! loop over all evaporations + wall losses etc. removing this cluster
				f(i) = f(i)-coef_lin(ind_lin_loss(i,n),ind_lin_loss(i,n+1),i)
			end do
			f(i) = f(i)*c(i) ! multiplying loss coefficients with current concentration
			! then add all terms that form this cluster
			do n = 1, 2*ind_quad_form(i,0)-1, 2 ! loop over all collisions forming this cluster
				f(i) = f(i)+coef_quad(ind_quad_form(i,n),ind_quad_form(i,n+1),i)*&
				&      c(ind_quad_form(i,n))*c(ind_quad_form(i,n+1))
			end do
			do n = 1, 3*ind_quad_form_extra(i,0)-2, 3 ! loop over all collisions forming this cluster as an extra product
				f(i) = f(i)+coef_quad(ind_quad_form_extra(i,n),ind_quad_form_extra(i,n+1),i)*&
				&      n_quad_form_extra(i,(n+2)/3)*c(ind_quad_form_extra(i,n))*c(ind_quad_form_extra(i,n+1))
			end do
			do n = 1, 2*ind_lin_form(i,0)-1, 2 ! loop over all evaporations forming this cluster
				f(i) = f(i)+coef_lin(i,ind_lin_form(i,n),ind_lin_form(i,n+1))*c(ind_lin_form(i,n+1))
			end do
			! finally, add possible external sources
			!f(i) = f(i) + source(i)
		end if
	end do
	do n = 1, fitted(1,0) ! loop over the clusters whose concentrations are fitted (e.g. [1A]+[1A1D]=const.)
		f(fitted(n+1,0)) = 0.d0 ! don't use value from birth-death equations
		!do i = 1, fitted(n+1,1) ! loop over the clusters used for the fitting
    do i = 1,1
			f(fitted(n+1,0)) = f(fitted(n+1,0))-f(fitted(n+1,i)) ! (e.g. d[1A]/dt = -d[1A1D]/dt)
		end do
	end do                  
 
end subroutine feval

!-----------------------------------------------------------

! jacobian of the differential equations: dfdc(i,j) = df(i)/dc(j)
! not using this, since the solver works faster using a numerical jacobian&
subroutine jeval(neqn,t,c,ml,mu,dfdc,ldim,coef,ipar)
	implicit none
	integer :: ldim,neqn,ierr,ipar(4),i,j,k,n
	real(kind(1.d0)) :: t,c(neqn),dfdc(ldim,neqn),coef(40),ml,mu

end subroutine jeval


!-----------------------------------------------------------

subroutine initialize_parameters(theta,coef_quad,coef_lin,  &
            & sources, ind_quad_loss,ind_quad_form, ind_lin_loss,ind_lin_form,ind_quad_loss_extra,&
            & ind_quad_form_extra,source,isconst,fitted,coef,ipar,K,E,cs)

	!use monomer_settings, only : sources_and_constants

	implicit none
    integer, parameter :: nclust=16
	integer :: neqn=16
	logical :: isconst(16)
	real(kind(1.d0)) :: coef_quad(nclust,nclust,23),coef_lin(23,23,nclust)
	real(kind(1.d0)) :: source(16)
	integer :: ind_quad_loss(23,0:34),ind_quad_form(23,0:122),ind_quad_loss_extra(23,0:30),ind_quad_form_extra(23,0:72)
	integer :: ind_lin_loss(23,0:12),ind_lin_form(23,0:32)
	integer :: fitted(nclust,0:nclust)
	real(kind(1.d0)) :: coef(29)
	integer :: ipar(4)
	real(kind(1.d0)) :: theta(28)
	real(kind(1.d0)) :: sources(2)
	real(kind(1.d0)) :: K(nclust,nclust), E(nclust,nclust), cs(nclust)
 
	!call sources_and_constants(source,isconst,fitted,ipar)
    source = 0.d0
	!source(1) = sources(1)
	!source(3) = sources(2)
	isconst(3) = .TRUE. !ammonia monomer
	
	fitted(1,0) = 1													! 1 concentration is fitted

	fitted(2,0:1) = (/1, 4/)
																	! concentration of 1A is fitted
																	!  using n-1 other concentrations:
																	!  [1A]+... = const. 1A and 1A1N
	
	ind_quad_loss = 0
	ind_quad_form = 0
	ind_lin_loss = 0
	ind_lin_form = 0
	ind_quad_loss_extra = 0
	ind_quad_form_extra = 0

	! ind_quad_loss(i,0:2n) = (/n,j1,k1,j2,k2,&,jn,kn/) gives the cluster numbers for all the collisions
	!  i + j1 -> k1 etc. through which cluster i is lost

	! ind_quad_form(k,0:2n) = (/n,i1,j1,i2,j2,&,in,jn/) gives the cluster numbers for all the collisions
	!  i1 + j1 -> k etc. through which cluster k is formed

	! ind_lin_loss(k,0:2n) = (/n,i1,j1,i2,j2,&,lossn,lossn/) gives the cluster numbers for all the evaporations
	!  k -> i1 + j1 etc. and other losses k -> wall etc. through which cluster k is lost

	! ind_lin_form(i,0:2n) = (/n,j1,k1,j2,k2,&,jn,kn/) gives the cluster numbers for all the evaporations
	!  k1 -> i + j1 etc. through which cluster i is formed

	! ind_quad_loss_extra(i,0:3n) = (/n,j1,k1,c1,&,jn,kn,cn/) gives the cluster numbers and coefficients 
	!  i + j1 -> c1*k1 etc. for additional collision products k (e.g. monomers from the boundary)
	!  in collisions where cluster i is lost

	! ind_quad_form_extra(k,0:2n) = (/n,i1,j1,c1,&,in,jn,cn/) gives the cluster numbers and coefficients 
	!  i1 + j1 -> c1*k etc. for additional ways of forming k (e.g. as a monomer from the boundary)


	! Cluster 1: 1A
	ind_quad_loss(1,0:24) = (/ 12, 1,2, 1,2, 3,4, 4,5, 5,6, 7,8, 8,9, 10,11, 11,12&
			&, 13,14, 15,16, 16,22 /)
	ind_quad_form_extra(1,0:72) = (/ 24, 2,2,2, 2,5,1, 2,6,2, 2,8,1, 2,9,2, 2,11,1, 2,12,2, 2,13,1, 2,14,2&
			&, 4,12,1, 5,6,1, 5,9,1, 5,11,1, 5,12,2, 6,6,2, 6,8,1, 6,9,2, 6,10,1, 6,11,2&
			&, 6,12,3, 7,9,1, 8,8,1, 8,9,2, 9,9,3 /)
	ind_lin_loss(1,0:2) = (/ 1, 18,18 /)
	ind_lin_form(1,0:22) = (/ 11, 1,2, 1,2, 3,4, 4,5, 5,6, 7,8, 8,9, 10,11, 11,12&
			&, 13,14, 15,16 /)

	! Cluster 2: 2A
	ind_quad_loss(2,0:32) = (/ 16, 2,2, 2,2, 3,5, 4,6, 5,6, 6,6, 7,9, 8,9, 9,9&
			&, 10,12, 11,12, 12,12, 13,14, 14,14, 15,22, 16,22 /)
	ind_quad_loss_extra(2,0:30) = (/ 10, 2,1,2, 2,1,2, 5,1,1, 6,1,2, 8,1,1, 9,1,2, 11,1,1, 12,1,2, 13,1,1&
			&, 14,1,2 /)
	ind_quad_form(2,0:4) = (/ 2, 1,1, 2,2 /)
	ind_lin_loss(2,0:4) = (/ 2, 1,1, 18,18 /)
	ind_lin_form(2,0:8) = (/ 4, 3,5, 4,6, 7,9, 10,12 /)

	! Cluster 3: 1N
	ind_quad_loss(3,0:20) = (/ 10, 1,4, 2,5, 5,7, 6,8, 8,10, 9,11, 11,13, 12,14, 13,15&
			&, 14,16 /)
	ind_quad_form_extra(3,0:3) = (/ 1, 4,15,1 /)
	ind_lin_loss(3,0:2) = (/ 1, 18,18 /)
	ind_lin_form(3,0:20) = (/ 10, 1,4, 2,5, 5,7, 6,8, 8,10, 9,11, 11,13, 12,14, 13,15&
			&, 14,16 /)

	! Cluster 4: 1A1N
	ind_quad_loss(4,0:32) = (/ 16, 1,5, 2,6, 4,7, 4,7, 5,8, 6,9, 7,10, 8,11, 9,12&
			&, 10,13, 11,14, 12,14, 13,16, 14,22, 15,16, 16,22 /)
	ind_quad_loss_extra(4,0:6) = (/ 2, 12,1,1, 15,3,1 /)
	ind_quad_form(4,0:2) = (/ 1, 1,3 /)
	ind_lin_loss(4,0:4) = (/ 2, 1,3, 18,18 /)
	ind_lin_form(4,0:24) = (/ 12, 1,5, 2,6, 4,7, 4,7, 5,8, 6,9, 7,10, 8,11, 9,12&
			&, 10,13, 11,14, 13,16 /)

	! Cluster 5: 2A1N
	ind_quad_loss(5,0:34) = (/ 17, 1,6, 2,6, 3,7, 4,8, 5,9, 5,9, 6,9, 7,11, 8,12&
			&, 9,12, 10,14, 11,14, 12,14, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(5,0:15) = (/ 5, 2,1,1, 6,1,1, 9,1,1, 11,1,1, 12,1,2 /)
	ind_quad_form(5,0:4) = (/ 2, 1,4, 2,3 /)
	ind_lin_loss(5,0:6) = (/ 3, 1,4, 2,3, 18,18 /)
	ind_lin_form(5,0:16) = (/ 8, 1,6, 3,7, 4,8, 5,9, 5,9, 7,11, 8,12, 10,14 /)

	! Cluster 6: 3A1N
	ind_quad_loss(6,0:32) = (/ 16, 2,6, 3,8, 4,9, 5,9, 6,9, 6,9, 7,12, 8,12, 9,12&
			&, 10,14, 11,14, 12,14, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(6,0:27) = (/ 9, 2,1,2, 5,1,1, 6,1,2, 6,1,2, 8,1,1, 9,1,2, 10,1,1, 11,1,2, 12,1,3 /)
	ind_quad_form(6,0:8) = (/ 4, 1,5, 2,4, 2,5, 2,6 /)
	ind_lin_loss(6,0:6) = (/ 3, 1,5, 2,4, 18,18 /)
	ind_lin_form(6,0:6) = (/ 3, 3,8, 4,9, 7,12 /)

	! Cluster 7: 2A2N
	ind_quad_loss(7,0:32) = (/ 16, 1,8, 2,9, 4,10, 5,11, 6,12, 7,13, 7,13, 8,14, 9,14&
			&, 10,16, 11,22, 12,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(7,0:3) = (/ 1, 9,1,1 /)
	ind_quad_form(7,0:4) = (/ 2, 3,5, 4,4 /)
	ind_lin_loss(7,0:6) = (/ 3, 3,5, 4,4, 18,18 /)
	ind_lin_form(7,0:18) = (/ 9, 1,8, 2,9, 4,10, 5,11, 6,12, 7,13, 7,13, 8,14, 10,16 /)

	! Cluster 8: 3A2N
	ind_quad_loss(8,0:34) = (/ 17, 1,9, 2,9, 3,10, 4,11, 5,12, 6,12, 7,14, 8,14, 8,14&
			&, 9,14, 10,22, 11,22, 12,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(8,0:15) = (/ 5, 2,1,1, 6,1,1, 8,1,1, 8,1,1, 9,1,2 /)
	ind_quad_form(8,0:6) = (/ 3, 1,7, 3,6, 4,5 /)
	ind_lin_loss(8,0:8) = (/ 4, 1,7, 3,6, 4,5, 18,18 /)
	ind_lin_form(8,0:10) = (/ 5, 1,9, 3,10, 4,11, 5,12, 7,14 /)

	! Cluster 9: 4A2N
	ind_quad_loss(9,0:32) = (/ 16, 2,9, 3,11, 4,12, 5,12, 6,12, 7,14, 8,14, 9,14, 9,14&
			&, 10,22, 11,22, 12,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(9,0:21) = (/ 7, 2,1,2, 5,1,1, 6,1,2, 7,1,1, 8,1,2, 9,1,3, 9,1,3 /)
	ind_quad_form(9,0:16) = (/ 8, 1,8, 2,7, 2,8, 2,9, 4,6, 5,5, 5,6, 6,6 /)
	ind_lin_loss(9,0:10) = (/ 5, 1,8, 2,7, 4,6, 5,5, 18,18 /)
	ind_lin_form(9,0:4) = (/ 2, 3,11, 4,12 /)

	! Cluster 10: 3A3N
	ind_quad_loss(10,0:32) = (/ 16, 1,11, 2,12, 4,13, 5,14, 6,14, 7,16, 8,22, 9,22, 10,22&
			&, 10,22, 11,22, 12,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(10,0:3) = (/ 1, 6,1,1 /)
	ind_quad_form(10,0:4) = (/ 2, 3,8, 4,7 /)
	ind_lin_loss(10,0:6) = (/ 3, 3,8, 4,7, 18,18 /)
	ind_lin_form(10,0:10) = (/ 5, 1,11, 2,12, 4,13, 5,14, 7,16 /)

	! Cluster 11: 4A3N
	ind_quad_loss(11,0:34) = (/ 17, 1,12, 2,12, 3,13, 4,14, 5,14, 6,14, 7,22, 8,22, 9,22&
			&, 10,22, 11,22, 11,22, 12,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(11,0:9) = (/ 3, 2,1,1, 5,1,1, 6,1,2 /)
	ind_quad_form(11,0:8) = (/ 4, 1,10, 3,9, 4,8, 5,7 /)
	ind_lin_loss(11,0:10) = (/ 5, 1,10, 3,9, 4,8, 5,7, 18,18 /)
	ind_lin_form(11,0:6) = (/ 3, 1,12, 3,13, 4,14 /)

	! Cluster 12: 5A3N
	ind_quad_loss(12,0:32) = (/ 16, 2,12, 3,14, 4,14, 5,14, 6,14, 7,22, 8,22, 9,22, 10,22&
			&, 11,22, 12,22, 12,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(12,0:12) = (/ 4, 2,1,2, 4,1,1, 5,1,2, 6,1,3 /)
	ind_quad_form(12,0:20) = (/ 10, 1,11, 2,10, 2,11, 2,12, 4,9, 5,8, 5,9, 6,7, 6,8&
			&, 6,9 /)
	ind_lin_loss(12,0:12) = (/ 6, 1,11, 2,10, 4,9, 5,8, 6,7, 18,18 /)
	ind_lin_form(12,0:2) = (/ 1, 3,14 /)

	! Cluster 13: 4A4N
	ind_quad_loss(13,0:34) = (/ 17, 1,14, 2,14, 3,15, 4,16, 5,22, 6,22, 7,22, 8,22, 9,22&
			&, 10,22, 11,22, 12,22, 13,22, 13,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(13,0:3) = (/ 1, 2,1,1 /)
	ind_quad_form(13,0:6) = (/ 3, 3,11, 4,10, 7,7 /)
	ind_lin_loss(13,0:8) = (/ 4, 3,11, 4,10, 7,7, 18,18 /)
	ind_lin_form(13,0:6) = (/ 3, 1,14, 3,15, 4,16 /)

	! Cluster 14: 5A4N
	ind_quad_loss(14,0:32) = (/ 16, 2,14, 3,16, 4,22, 5,22, 6,22, 7,22, 8,22, 9,22, 10,22&
			&, 11,22, 12,22, 13,22, 14,22, 14,22, 15,22, 16,22 /)
	ind_quad_loss_extra(14,0:3) = (/ 1, 2,1,2 /)
	ind_quad_form(14,0:34) = (/ 17, 1,13, 2,13, 2,14, 3,12, 4,11, 4,12, 5,10, 5,11, 5,12&
			&, 6,10, 6,11, 6,12, 7,8, 7,9, 8,8, 8,9, 9,9 /)
	ind_lin_loss(14,0:12) = (/ 6, 1,13, 3,12, 4,11, 5,10, 7,8, 18,18 /)
	ind_lin_form(14,0:2) = (/ 1, 3,16 /)

	! Cluster 15: 4A5N
	ind_quad_loss(15,0:32) = (/ 16, 1,16, 2,22, 4,16, 5,22, 6,22, 7,22, 8,22, 9,22, 10,22&
			&, 11,22, 12,22, 13,22, 14,22, 15,22, 15,22, 16,22 /)
	ind_quad_loss_extra(15,0:3) = (/ 1, 4,3,1 /)
	ind_quad_form(15,0:2) = (/ 1, 3,13 /)
	ind_lin_loss(15,0:4) = (/ 2, 3,13, 18,18 /)
	ind_lin_form(15,0:2) = (/ 1, 1,16 /)

	! Cluster 16: 5A5N
	ind_quad_loss(16,0:32) = (/ 16, 1,22, 2,22, 4,22, 5,22, 6,22, 7,22, 8,22, 9,22, 10,22&
			&, 11,22, 12,22, 13,22, 14,22, 15,22, 16,22, 16,22 /)
	ind_quad_form(16,0:10) = (/ 5, 1,15, 3,14, 4,13, 4,15, 7,10 /)
	ind_lin_loss(16,0:10) = (/ 5, 1,15, 3,14, 4,13, 7,10, 18,18 /)

	! out_neu
	ind_quad_form(22,0:122) = (/ 61, 1,16, 2,15, 2,16, 4,14, 4,16, 5,13, 5,14, 5,15, 5,16&
			&, 6,13, 6,14, 6,15, 6,16, 7,11, 7,12, 7,13, 7,14, 7,15, 7,16&
			&, 8,10, 8,11, 8,12, 8,13, 8,14, 8,15, 8,16, 9,10, 9,11, 9,12&
			&, 9,13, 9,14, 9,15, 9,16, 10,10, 10,11, 10,12, 10,13, 10,14, 10,15&
			&, 10,16, 11,11, 11,12, 11,13, 11,14, 11,15, 11,16, 12,12, 12,13, 12,14&
			&, 12,15, 12,16, 13,13, 13,14, 13,15, 13,16, 14,14, 14,15, 14,16, 15,15&
			&, 15,16, 16,16 /)

	call get_rate_coefs(theta,coef_quad,coef_lin,coef,K,E,cs)

end subroutine initialize_parameters

!-----------------------------------------------------------

subroutine get_rate_coefs(theta,coef_quad,coef_lin,coef,K,E,cs)
	implicit none
  integer, parameter :: nclust=16
	real(kind(1.d0)) :: coef_quad(nclust,nclust,23),coef_lin(23,23,nclust)
	real(kind(1.d0)) :: theta(28),dH(22),dS(22)
	real(kind(1.d0)) :: K(nclust,nclust),E(nclust,nclust),cs(nclust)
	real(kind(1.d0)) :: coef(29), temp
	integer			 ::  i, j
  
	coef_quad = 0.d0
	coef_lin = 0.d0
	!call get_dH_dS_from_theta(dH,dS,theta)
    temp = coef(1)
	call get_coll(K,temp)
	call get_evap_mcmc(theta,E,K,temp)
	call get_losses(cs)
  
	coef_quad(1,1,2) = 0.5d0*K(1,1)	! 1A + 1A -> 2A
	coef_lin(1,1,2) = E(1,1)	! 2A -> 1A + 1A

	coef_quad(2,2,2) = 0.5d0*K(2,2)	! 2A + 2A -> boundary -> 2A
	coef_quad(2,2,1) = 0.5d0*K(2,2)	! 2A + 2A -> boundary -> 1A

	coef_quad(3,1,4) = K(1,3)	! 1N + 1A -> 1A1N
	coef_quad(1,3,4) = K(1,3)
	coef_lin(3,1,4) = E(1,3)	! 1A1N -> 1N + 1A
	coef_lin(1,3,4) = E(1,3)
	coef_quad(3,2,5) = K(2,3)	! 1N + 2A -> 2A1N
	coef_quad(2,3,5) = K(2,3)
	coef_lin(3,2,5) = E(2,3)	! 2A1N -> 1N + 2A
	coef_lin(2,3,5) = E(2,3)

	coef_quad(4,1,5) = K(1,4)	! 1A1N + 1A -> 2A1N
	coef_quad(1,4,5) = K(1,4)
	coef_lin(4,1,5) = E(1,4)	! 2A1N -> 1A1N + 1A
	coef_lin(1,4,5) = E(1,4)
	coef_quad(4,2,6) = K(2,4)	! 1A1N + 2A -> 3A1N
	coef_quad(2,4,6) = K(2,4)
	coef_lin(4,2,6) = E(2,4)	! 3A1N -> 1A1N + 2A
	coef_lin(2,4,6) = E(2,4)
	coef_quad(4,4,7) = 0.5d0*K(4,4)	! 1A1N + 1A1N -> 2A2N
	coef_lin(4,4,7) = E(4,4)	! 2A2N -> 1A1N + 1A1N

	coef_quad(5,1,6) = K(1,5)	! 2A1N + 1A -> 3A1N
	coef_quad(1,5,6) = K(1,5)
	coef_lin(5,1,6) = E(1,5)	! 3A1N -> 2A1N + 1A
	coef_lin(1,5,6) = E(1,5)
	coef_quad(5,2,6) = K(2,5)	! 2A1N + 2A -> boundary -> 3A1N
	coef_quad(2,5,6) = K(2,5)
	coef_quad(5,2,1) = K(2,5)	! 2A1N + 2A -> boundary -> 1A
	coef_quad(2,5,1) = K(2,5)
	coef_quad(5,3,7) = K(3,5)	! 2A1N + 1N -> 2A2N
	coef_quad(3,5,7) = K(3,5)
	coef_lin(5,3,7) = E(3,5)	! 2A2N -> 2A1N + 1N
	coef_lin(3,5,7) = E(3,5)
	coef_quad(5,4,8) = K(4,5)	! 2A1N + 1A1N -> 3A2N
	coef_quad(4,5,8) = K(4,5)
	coef_lin(5,4,8) = E(4,5)	! 3A2N -> 2A1N + 1A1N
	coef_lin(4,5,8) = E(4,5)
	coef_quad(5,5,9) = 0.5d0*K(5,5)	! 2A1N + 2A1N -> 4A2N
	coef_lin(5,5,9) = E(5,5)	! 4A2N -> 2A1N + 2A1N

	coef_quad(6,2,6) = K(2,6)	! 3A1N + 2A -> boundary -> 3A1N
	coef_quad(2,6,6) = K(2,6)
	coef_quad(6,2,1) = K(2,6)	! 3A1N + 2A -> boundary -> 1A
	coef_quad(2,6,1) = K(2,6)
	coef_quad(6,3,8) = K(3,6)	! 3A1N + 1N -> 3A2N
	coef_quad(3,6,8) = K(3,6)
	coef_lin(6,3,8) = E(3,6)	! 3A2N -> 3A1N + 1N
	coef_lin(3,6,8) = E(3,6)
	coef_quad(6,4,9) = K(4,6)	! 3A1N + 1A1N -> 4A2N
	coef_quad(4,6,9) = K(4,6)
	coef_lin(6,4,9) = E(4,6)	! 4A2N -> 3A1N + 1A1N
	coef_lin(4,6,9) = E(4,6)
	coef_quad(6,5,9) = K(5,6)	! 3A1N + 2A1N -> boundary -> 4A2N
	coef_quad(5,6,9) = K(5,6)
	coef_quad(6,5,1) = K(5,6)	! 3A1N + 2A1N -> boundary -> 1A
	coef_quad(5,6,1) = K(5,6)
	coef_quad(6,6,9) = 0.5d0*K(6,6)	! 3A1N + 3A1N -> boundary -> 4A2N
	coef_quad(6,6,1) = 0.5d0*K(6,6)	! 3A1N + 3A1N -> boundary -> 1A

	coef_quad(7,1,8) = K(1,7)	! 2A2N + 1A -> 3A2N
	coef_quad(1,7,8) = K(1,7)
	coef_lin(7,1,8) = E(1,7)	! 3A2N -> 2A2N + 1A
	coef_lin(1,7,8) = E(1,7)
	coef_quad(7,2,9) = K(2,7)	! 2A2N + 2A -> 4A2N
	coef_quad(2,7,9) = K(2,7)
	coef_lin(7,2,9) = E(2,7)	! 4A2N -> 2A2N + 2A
	coef_lin(2,7,9) = E(2,7)
	coef_quad(7,4,10) = K(4,7)	! 2A2N + 1A1N -> 3A3N
	coef_quad(4,7,10) = K(4,7)
	coef_lin(7,4,10) = E(4,7)	! 3A3N -> 2A2N + 1A1N
	coef_lin(4,7,10) = E(4,7)
	coef_quad(7,5,11) = K(5,7)	! 2A2N + 2A1N -> 4A3N
	coef_quad(5,7,11) = K(5,7)
	coef_lin(7,5,11) = E(5,7)	! 4A3N -> 2A2N + 2A1N
	coef_lin(5,7,11) = E(5,7)
	coef_quad(7,6,12) = K(6,7)	! 2A2N + 3A1N -> 5A3N
	coef_quad(6,7,12) = K(6,7)
	coef_lin(7,6,12) = E(6,7)	! 5A3N -> 2A2N + 3A1N
	coef_lin(6,7,12) = E(6,7)
	coef_quad(7,7,13) = 0.5d0*K(7,7)	! 2A2N + 2A2N -> 4A4N
	coef_lin(7,7,13) = E(7,7)	! 4A4N -> 2A2N + 2A2N

	coef_quad(8,1,9) = K(1,8)	! 3A2N + 1A -> 4A2N
	coef_quad(1,8,9) = K(1,8)
	coef_lin(8,1,9) = E(1,8)	! 4A2N -> 3A2N + 1A
	coef_lin(1,8,9) = E(1,8)
	coef_quad(8,2,9) = K(2,8)	! 3A2N + 2A -> boundary -> 4A2N
	coef_quad(2,8,9) = K(2,8)
	coef_quad(8,2,1) = K(2,8)	! 3A2N + 2A -> boundary -> 1A
	coef_quad(2,8,1) = K(2,8)
	coef_quad(8,3,10) = K(3,8)	! 3A2N + 1N -> 3A3N
	coef_quad(3,8,10) = K(3,8)
	coef_lin(8,3,10) = E(3,8)	! 3A3N -> 3A2N + 1N
	coef_lin(3,8,10) = E(3,8)
	coef_quad(8,4,11) = K(4,8)	! 3A2N + 1A1N -> 4A3N
	coef_quad(4,8,11) = K(4,8)
	coef_lin(8,4,11) = E(4,8)	! 4A3N -> 3A2N + 1A1N
	coef_lin(4,8,11) = E(4,8)
	coef_quad(8,5,12) = K(5,8)	! 3A2N + 2A1N -> 5A3N
	coef_quad(5,8,12) = K(5,8)
	coef_lin(8,5,12) = E(5,8)	! 5A3N -> 3A2N + 2A1N
	coef_lin(5,8,12) = E(5,8)
	coef_quad(8,6,12) = K(6,8)	! 3A2N + 3A1N -> boundary -> 5A3N
	coef_quad(6,8,12) = K(6,8)
	coef_quad(8,6,1) = K(6,8)	! 3A2N + 3A1N -> boundary -> 1A
	coef_quad(6,8,1) = K(6,8)
	coef_quad(8,7,14) = K(7,8)	! 3A2N + 2A2N -> 5A4N
	coef_quad(7,8,14) = K(7,8)
	coef_lin(8,7,14) = E(7,8)	! 5A4N -> 3A2N + 2A2N
	coef_lin(7,8,14) = E(7,8)
	coef_quad(8,8,14) = 0.5d0*K(8,8)	! 3A2N + 3A2N -> boundary -> 5A4N
	coef_quad(8,8,1) = 0.5d0*K(8,8)	! 3A2N + 3A2N -> boundary -> 1A

	coef_quad(9,2,9) = K(2,9)	! 4A2N + 2A -> boundary -> 4A2N
	coef_quad(2,9,9) = K(2,9)
	coef_quad(9,2,1) = K(2,9)	! 4A2N + 2A -> boundary -> 1A
	coef_quad(2,9,1) = K(2,9)
	coef_quad(9,3,11) = K(3,9)	! 4A2N + 1N -> 4A3N
	coef_quad(3,9,11) = K(3,9)
	coef_lin(9,3,11) = E(3,9)	! 4A3N -> 4A2N + 1N
	coef_lin(3,9,11) = E(3,9)
	coef_quad(9,4,12) = K(4,9)	! 4A2N + 1A1N -> 5A3N
	coef_quad(4,9,12) = K(4,9)
	coef_lin(9,4,12) = E(4,9)	! 5A3N -> 4A2N + 1A1N
	coef_lin(4,9,12) = E(4,9)
	coef_quad(9,5,12) = K(5,9)	! 4A2N + 2A1N -> boundary -> 5A3N
	coef_quad(5,9,12) = K(5,9)
	coef_quad(9,5,1) = K(5,9)	! 4A2N + 2A1N -> boundary -> 1A
	coef_quad(5,9,1) = K(5,9)
	coef_quad(9,6,12) = K(6,9)	! 4A2N + 3A1N -> boundary -> 5A3N
	coef_quad(6,9,12) = K(6,9)
	coef_quad(9,6,1) = K(6,9)	! 4A2N + 3A1N -> boundary -> 1A
	coef_quad(6,9,1) = K(6,9)
	coef_quad(9,7,14) = K(7,9)	! 4A2N + 2A2N -> boundary -> 5A4N
	coef_quad(7,9,14) = K(7,9)
	coef_quad(9,7,1) = K(7,9)	! 4A2N + 2A2N -> boundary -> 1A
	coef_quad(7,9,1) = K(7,9)
	coef_quad(9,8,14) = K(8,9)	! 4A2N + 3A2N -> boundary -> 5A4N
	coef_quad(8,9,14) = K(8,9)
	coef_quad(9,8,1) = K(8,9)	! 4A2N + 3A2N -> boundary -> 1A
	coef_quad(8,9,1) = K(8,9)
	coef_quad(9,9,14) = 0.5d0*K(9,9)	! 4A2N + 4A2N -> boundary -> 5A4N
	coef_quad(9,9,1) = 0.5d0*K(9,9)	! 4A2N + 4A2N -> boundary -> 1A

	coef_quad(10,1,11) = K(1,10)	! 3A3N + 1A -> 4A3N
	coef_quad(1,10,11) = K(1,10)
	coef_lin(10,1,11) = E(1,10)	! 4A3N -> 3A3N + 1A
	coef_lin(1,10,11) = E(1,10)
	coef_quad(10,2,12) = K(2,10)	! 3A3N + 2A -> 5A3N
	coef_quad(2,10,12) = K(2,10)
	coef_lin(10,2,12) = E(2,10)	! 5A3N -> 3A3N + 2A
	coef_lin(2,10,12) = E(2,10)
	coef_quad(10,4,13) = K(4,10)	! 3A3N + 1A1N -> 4A4N
	coef_quad(4,10,13) = K(4,10)
	coef_lin(10,4,13) = E(4,10)	! 4A4N -> 3A3N + 1A1N
	coef_lin(4,10,13) = E(4,10)
	coef_quad(10,5,14) = K(5,10)	! 3A3N + 2A1N -> 5A4N
	coef_quad(5,10,14) = K(5,10)
	coef_lin(10,5,14) = E(5,10)	! 5A4N -> 3A3N + 2A1N
	coef_lin(5,10,14) = E(5,10)
	coef_quad(10,6,14) = K(6,10)	! 3A3N + 3A1N -> boundary -> 5A4N
	coef_quad(6,10,14) = K(6,10)
	coef_quad(10,6,1) = K(6,10)	! 3A3N + 3A1N -> boundary -> 1A
	coef_quad(6,10,1) = K(6,10)
	coef_quad(10,7,16) = K(7,10)	! 3A3N + 2A2N -> 5A5N
	coef_quad(7,10,16) = K(7,10)
	coef_lin(10,7,16) = E(7,10)	! 5A5N -> 3A3N + 2A2N
	coef_lin(7,10,16) = E(7,10)
	coef_quad(10,8,22) = K(8,10)	! 3A3N + 3A2N -> out_neu
	coef_quad(8,10,22) = K(8,10)
	coef_quad(10,9,22) = K(9,10)	! 3A3N + 4A2N -> out_neu
	coef_quad(9,10,22) = K(9,10)
	coef_quad(10,10,22) = 0.5d0*K(10,10)	! 3A3N + 3A3N -> out_neu

	coef_quad(11,1,12) = K(1,11)	! 4A3N + 1A -> 5A3N
	coef_quad(1,11,12) = K(1,11)
	coef_lin(11,1,12) = E(1,11)	! 5A3N -> 4A3N + 1A
	coef_lin(1,11,12) = E(1,11)
	coef_quad(11,2,12) = K(2,11)	! 4A3N + 2A -> boundary -> 5A3N
	coef_quad(2,11,12) = K(2,11)
	coef_quad(11,2,1) = K(2,11)	! 4A3N + 2A -> boundary -> 1A
	coef_quad(2,11,1) = K(2,11)
	coef_quad(11,3,13) = K(3,11)	! 4A3N + 1N -> 4A4N
	coef_quad(3,11,13) = K(3,11)
	coef_lin(11,3,13) = E(3,11)	! 4A4N -> 4A3N + 1N
	coef_lin(3,11,13) = E(3,11)
	coef_quad(11,4,14) = K(4,11)	! 4A3N + 1A1N -> 5A4N
	coef_quad(4,11,14) = K(4,11)
	coef_lin(11,4,14) = E(4,11)	! 5A4N -> 4A3N + 1A1N
	coef_lin(4,11,14) = E(4,11)
	coef_quad(11,5,14) = K(5,11)	! 4A3N + 2A1N -> boundary -> 5A4N
	coef_quad(5,11,14) = K(5,11)
	coef_quad(11,5,1) = K(5,11)	! 4A3N + 2A1N -> boundary -> 1A
	coef_quad(5,11,1) = K(5,11)
	coef_quad(11,6,14) = K(6,11)	! 4A3N + 3A1N -> boundary -> 5A4N
	coef_quad(6,11,14) = K(6,11)
	coef_quad(11,6,1) = K(6,11)	! 4A3N + 3A1N -> boundary -> 1A
	coef_quad(6,11,1) = K(6,11)
	coef_quad(11,7,22) = K(7,11)	! 4A3N + 2A2N -> out_neu
	coef_quad(7,11,22) = K(7,11)
	coef_quad(11,8,22) = K(8,11)	! 4A3N + 3A2N -> out_neu
	coef_quad(8,11,22) = K(8,11)
	coef_quad(11,9,22) = K(9,11)	! 4A3N + 4A2N -> out_neu
	coef_quad(9,11,22) = K(9,11)
	coef_quad(11,10,22) = K(10,11)	! 4A3N + 3A3N -> out_neu
	coef_quad(10,11,22) = K(10,11)
	coef_quad(11,11,22) = 0.5d0*K(11,11)	! 4A3N + 4A3N -> out_neu

	coef_quad(12,2,12) = K(2,12)	! 5A3N + 2A -> boundary -> 5A3N
	coef_quad(2,12,12) = K(2,12)
	coef_quad(12,2,1) = K(2,12)	! 5A3N + 2A -> boundary -> 1A
	coef_quad(2,12,1) = K(2,12)
	coef_quad(12,3,14) = K(3,12)	! 5A3N + 1N -> 5A4N
	coef_quad(3,12,14) = K(3,12)
	coef_lin(12,3,14) = E(3,12)	! 5A4N -> 5A3N + 1N
	coef_lin(3,12,14) = E(3,12)
	coef_quad(12,4,14) = K(4,12)	! 5A3N + 1A1N -> boundary -> 5A4N
	coef_quad(4,12,14) = K(4,12)
	coef_quad(12,4,1) = K(4,12)	! 5A3N + 1A1N -> boundary -> 1A
	coef_quad(4,12,1) = K(4,12)
	coef_quad(12,5,14) = K(5,12)	! 5A3N + 2A1N -> boundary -> 5A4N
	coef_quad(5,12,14) = K(5,12)
	coef_quad(12,5,1) = K(5,12)	! 5A3N + 2A1N -> boundary -> 1A
	coef_quad(5,12,1) = K(5,12)
	coef_quad(12,6,14) = K(6,12)	! 5A3N + 3A1N -> boundary -> 5A4N
	coef_quad(6,12,14) = K(6,12)
	coef_quad(12,6,1) = K(6,12)	! 5A3N + 3A1N -> boundary -> 1A
	coef_quad(6,12,1) = K(6,12)
	coef_quad(12,7,22) = K(7,12)	! 5A3N + 2A2N -> out_neu
	coef_quad(7,12,22) = K(7,12)
	coef_quad(12,8,22) = K(8,12)	! 5A3N + 3A2N -> out_neu
	coef_quad(8,12,22) = K(8,12)
	coef_quad(12,9,22) = K(9,12)	! 5A3N + 4A2N -> out_neu
	coef_quad(9,12,22) = K(9,12)
	coef_quad(12,10,22) = K(10,12)	! 5A3N + 3A3N -> out_neu
	coef_quad(10,12,22) = K(10,12)
	coef_quad(12,11,22) = K(11,12)	! 5A3N + 4A3N -> out_neu
	coef_quad(11,12,22) = K(11,12)
	coef_quad(12,12,22) = 0.5d0*K(12,12)	! 5A3N + 5A3N -> out_neu

	coef_quad(13,1,14) = K(1,13)	! 4A4N + 1A -> 5A4N
	coef_quad(1,13,14) = K(1,13)
	coef_lin(13,1,14) = E(1,13)	! 5A4N -> 4A4N + 1A
	coef_lin(1,13,14) = E(1,13)
	coef_quad(13,2,14) = K(2,13)	! 4A4N + 2A -> boundary -> 5A4N
	coef_quad(2,13,14) = K(2,13)
	coef_quad(13,2,1) = K(2,13)	! 4A4N + 2A -> boundary -> 1A
	coef_quad(2,13,1) = K(2,13)
	coef_quad(13,3,15) = K(3,13)	! 4A4N + 1N -> 4A5N
	coef_quad(3,13,15) = K(3,13)
	coef_lin(13,3,15) = E(3,13)	! 4A5N -> 4A4N + 1N
	coef_lin(3,13,15) = E(3,13)
	coef_quad(13,4,16) = K(4,13)	! 4A4N + 1A1N -> 5A5N
	coef_quad(4,13,16) = K(4,13)
	coef_lin(13,4,16) = E(4,13)	! 5A5N -> 4A4N + 1A1N
	coef_lin(4,13,16) = E(4,13)
	coef_quad(13,5,22) = K(5,13)	! 4A4N + 2A1N -> out_neu
	coef_quad(5,13,22) = K(5,13)
	coef_quad(13,6,22) = K(6,13)	! 4A4N + 3A1N -> out_neu
	coef_quad(6,13,22) = K(6,13)
	coef_quad(13,7,22) = K(7,13)	! 4A4N + 2A2N -> out_neu
	coef_quad(7,13,22) = K(7,13)
	coef_quad(13,8,22) = K(8,13)	! 4A4N + 3A2N -> out_neu
	coef_quad(8,13,22) = K(8,13)
	coef_quad(13,9,22) = K(9,13)	! 4A4N + 4A2N -> out_neu
	coef_quad(9,13,22) = K(9,13)
	coef_quad(13,10,22) = K(10,13)	! 4A4N + 3A3N -> out_neu
	coef_quad(10,13,22) = K(10,13)
	coef_quad(13,11,22) = K(11,13)	! 4A4N + 4A3N -> out_neu
	coef_quad(11,13,22) = K(11,13)
	coef_quad(13,12,22) = K(12,13)	! 4A4N + 5A3N -> out_neu
	coef_quad(12,13,22) = K(12,13)
	coef_quad(13,13,22) = 0.5d0*K(13,13)	! 4A4N + 4A4N -> out_neu

	coef_quad(14,2,14) = K(2,14)	! 5A4N + 2A -> boundary -> 5A4N
	coef_quad(2,14,14) = K(2,14)
	coef_quad(14,2,1) = K(2,14)	! 5A4N + 2A -> boundary -> 1A
	coef_quad(2,14,1) = K(2,14)
	coef_quad(14,3,16) = K(3,14)	! 5A4N + 1N -> 5A5N
	coef_quad(3,14,16) = K(3,14)
	coef_lin(14,3,16) = E(3,14)	! 5A5N -> 5A4N + 1N
	coef_lin(3,14,16) = E(3,14)
	coef_quad(14,4,22) = K(4,14)	! 5A4N + 1A1N -> out_neu
	coef_quad(4,14,22) = K(4,14)
	coef_quad(14,5,22) = K(5,14)	! 5A4N + 2A1N -> out_neu
	coef_quad(5,14,22) = K(5,14)
	coef_quad(14,6,22) = K(6,14)	! 5A4N + 3A1N -> out_neu
	coef_quad(6,14,22) = K(6,14)
	coef_quad(14,7,22) = K(7,14)	! 5A4N + 2A2N -> out_neu
	coef_quad(7,14,22) = K(7,14)
	coef_quad(14,8,22) = K(8,14)	! 5A4N + 3A2N -> out_neu
	coef_quad(8,14,22) = K(8,14)
	coef_quad(14,9,22) = K(9,14)	! 5A4N + 4A2N -> out_neu
	coef_quad(9,14,22) = K(9,14)
	coef_quad(14,10,22) = K(10,14)	! 5A4N + 3A3N -> out_neu
	coef_quad(10,14,22) = K(10,14)
	coef_quad(14,11,22) = K(11,14)	! 5A4N + 4A3N -> out_neu
	coef_quad(11,14,22) = K(11,14)
	coef_quad(14,12,22) = K(12,14)	! 5A4N + 5A3N -> out_neu
	coef_quad(12,14,22) = K(12,14)
	coef_quad(14,13,22) = K(13,14)	! 5A4N + 4A4N -> out_neu
	coef_quad(13,14,22) = K(13,14)
	coef_quad(14,14,22) = 0.5d0*K(14,14)	! 5A4N + 5A4N -> out_neu

	coef_quad(15,1,16) = K(1,15)	! 4A5N + 1A -> 5A5N
	coef_quad(1,15,16) = K(1,15)
	coef_lin(15,1,16) = E(1,15)	! 5A5N -> 4A5N + 1A
	coef_lin(1,15,16) = E(1,15)
	coef_quad(15,2,22) = K(2,15)	! 4A5N + 2A -> out_neu
	coef_quad(2,15,22) = K(2,15)
	coef_quad(15,4,16) = K(4,15)	! 4A5N + 1A1N -> boundary -> 5A5N
	coef_quad(4,15,16) = K(4,15)
	coef_quad(15,4,3) = K(4,15)	! 4A5N + 1A1N -> boundary -> 1N
	coef_quad(4,15,3) = K(4,15)
	coef_quad(15,5,22) = K(5,15)	! 4A5N + 2A1N -> out_neu
	coef_quad(5,15,22) = K(5,15)
	coef_quad(15,6,22) = K(6,15)	! 4A5N + 3A1N -> out_neu
	coef_quad(6,15,22) = K(6,15)
	coef_quad(15,7,22) = K(7,15)	! 4A5N + 2A2N -> out_neu
	coef_quad(7,15,22) = K(7,15)
	coef_quad(15,8,22) = K(8,15)	! 4A5N + 3A2N -> out_neu
	coef_quad(8,15,22) = K(8,15)
	coef_quad(15,9,22) = K(9,15)	! 4A5N + 4A2N -> out_neu
	coef_quad(9,15,22) = K(9,15)
	coef_quad(15,10,22) = K(10,15)	! 4A5N + 3A3N -> out_neu
	coef_quad(10,15,22) = K(10,15)
	coef_quad(15,11,22) = K(11,15)	! 4A5N + 4A3N -> out_neu
	coef_quad(11,15,22) = K(11,15)
	coef_quad(15,12,22) = K(12,15)	! 4A5N + 5A3N -> out_neu
	coef_quad(12,15,22) = K(12,15)
	coef_quad(15,13,22) = K(13,15)	! 4A5N + 4A4N -> out_neu
	coef_quad(13,15,22) = K(13,15)
	coef_quad(15,14,22) = K(14,15)	! 4A5N + 5A4N -> out_neu
	coef_quad(14,15,22) = K(14,15)
	coef_quad(15,15,22) = 0.5d0*K(15,15)	! 4A5N + 4A5N -> out_neu

	coef_quad(16,1,22) = K(1,16)	! 5A5N + 1A -> out_neu
	coef_quad(1,16,22) = K(1,16)
	coef_quad(16,2,22) = K(2,16)	! 5A5N + 2A -> out_neu
	coef_quad(2,16,22) = K(2,16)
	coef_quad(16,4,22) = K(4,16)	! 5A5N + 1A1N -> out_neu
	coef_quad(4,16,22) = K(4,16)
	coef_quad(16,5,22) = K(5,16)	! 5A5N + 2A1N -> out_neu
	coef_quad(5,16,22) = K(5,16)
	coef_quad(16,6,22) = K(6,16)	! 5A5N + 3A1N -> out_neu
	coef_quad(6,16,22) = K(6,16)
	coef_quad(16,7,22) = K(7,16)	! 5A5N + 2A2N -> out_neu
	coef_quad(7,16,22) = K(7,16)
	coef_quad(16,8,22) = K(8,16)	! 5A5N + 3A2N -> out_neu
	coef_quad(8,16,22) = K(8,16)
	coef_quad(16,9,22) = K(9,16)	! 5A5N + 4A2N -> out_neu
	coef_quad(9,16,22) = K(9,16)
	coef_quad(16,10,22) = K(10,16)	! 5A5N + 3A3N -> out_neu
	coef_quad(10,16,22) = K(10,16)
	coef_quad(16,11,22) = K(11,16)	! 5A5N + 4A3N -> out_neu
	coef_quad(11,16,22) = K(11,16)
	coef_quad(16,12,22) = K(12,16)	! 5A5N + 5A3N -> out_neu
	coef_quad(12,16,22) = K(12,16)
	coef_quad(16,13,22) = K(13,16)	! 5A5N + 4A4N -> out_neu
	coef_quad(13,16,22) = K(13,16)
	coef_quad(16,14,22) = K(14,16)	! 5A5N + 5A4N -> out_neu
	coef_quad(14,16,22) = K(14,16)
	coef_quad(16,15,22) = K(15,16)	! 5A5N + 4A5N -> out_neu
	coef_quad(15,16,22) = K(15,16)
	coef_quad(16,16,22) = 0.5d0*K(16,16)	! 5A5N + 5A5N -> out_neu

	coef_lin(18,18,16) = cs(16)	! 5A5N -> coag
	coef_lin(18,18,15) = cs(15)	! 4A5N -> coag
	coef_lin(18,18,14) = cs(14)	! 5A4N -> coag
	coef_lin(18,18,13) = cs(13)	! 4A4N -> coag
	coef_lin(18,18,12) = cs(12)	! 5A3N -> coag
	coef_lin(18,18,11) = cs(11)	! 4A3N -> coag
	coef_lin(18,18,10) = cs(10)	! 3A3N -> coag
	coef_lin(18,18,9) = cs(9)	! 4A2N -> coag
	coef_lin(18,18,8) = cs(8)	! 3A2N -> coag
	coef_lin(18,18,7) = cs(7)	! 2A2N -> coag
	coef_lin(18,18,6) = cs(6)	! 3A1N -> coag
	coef_lin(18,18,5) = cs(5)	! 2A1N -> coag
	coef_lin(18,18,4) = cs(4)	! 1A1N -> coag
	coef_lin(18,18,3) = cs(3)	! 1N -> coag
	coef_lin(18,18,2) = cs(2)	! 2A -> coag
	coef_lin(18,18,1) = cs(1)	! 1A -> coag

end subroutine get_rate_coefs

!-----------------------------------------------------------

subroutine get_losses(cs)
	implicit none
	real(kind(1.d0)) :: cs(16)

	! coagulation sink

	cs(1) = 1.00000000000000d-003	! coagulation loss of 1A
	cs(2) = 6.90956439983888d-004	! coagulation loss of 2A
	cs(3) = 1.51871684105631d-003	! coagulation loss of 1N
	cs(4) = 8.18186956498848d-004	! coagulation loss of 1A1N
	cs(5) = 6.19159137757835d-004	! coagulation loss of 2A1N
	cs(6) = 5.16067750384088d-004	! coagulation loss of 3A1N
	cs(7) = 5.65331546703696d-004	! coagulation loss of 2A2N
	cs(8) = 4.83012784469347d-004	! coagulation loss of 3A2N
	cs(9) = 4.27811993608647d-004	! coagulation loss of 4A2N
	cs(10) = 4.55394606310573d-004	! coagulation loss of 3A3N
	cs(11) = 4.08002160844876d-004	! coagulation loss of 4A3N
	cs(12) = 3.72486494677588d-004	! coagulation loss of 5A3N
	cs(13) = 3.90619472920971d-004	! coagulation loss of 4A4N
	cs(14) = 3.58979709679220d-004	! coagulation loss of 5A4N
	cs(15) = 3.75208865762285d-004	! coagulation loss of 4A5N
	cs(16) = 3.46791584003977d-004	! coagulation loss of 5A5N

end subroutine get_losses

!-----------------------------------------------------------

subroutine get_coll(K,temperature)
	implicit none
	integer, parameter :: nclust = 16
	real(kind(1.d0)) :: K(nclust,nclust), temperature
 integer :: i,j

	! collision coefficients

	K = 0.d0
K(1,1) = 2.00300435981201d-017*sqrt(temperature);	! 1A + 1A
K(2,2) = 2.24829637648719d-017*sqrt(temperature);	! 2A + 2A
K(3,1) = 2.88382629050621d-017*sqrt(temperature);	! 1N + 1A
K(1,3) = K(3,1);
K(3,2) = 3.64984531663620d-017*sqrt(temperature);	! 1N + 2A
K(2,3) = K(3,2);
K(4,1) = 2.19360737069264d-017*sqrt(temperature);	! 1A1N + 1A
K(1,4) = K(4,1);
K(4,2) = 2.35870463410389d-017*sqrt(temperature);	! 1A1N + 2A
K(2,4) = K(4,2);
K(4,4) = 2.37591578600701d-017*sqrt(temperature);	! 1A1N + 1A1N
K(5,1) = 2.36147395363845d-017*sqrt(temperature);	! 2A1N + 1A
K(1,5) = K(5,1);
K(5,2) = 2.36204584489056d-017*sqrt(temperature);	! 2A1N + 2A
K(2,5) = K(5,2);
K(5,3) = 3.96575479722936d-017*sqrt(temperature);	! 2A1N + 1N
K(3,5) = K(5,3);
K(5,4) = 2.50045077951280d-017*sqrt(temperature);	! 2A1N + 1A1N
K(4,5) = K(5,4);
K(5,5) = 2.47357831269868d-017*sqrt(temperature);	! 2A1N + 2A1N
K(6,2) = 2.45627767749069d-017*sqrt(temperature);	! 3A1N + 2A
K(2,6) = K(6,2);
K(6,3) = 4.54394981584564d-017*sqrt(temperature);	! 3A1N + 1N
K(3,6) = K(6,3);
K(6,4) = 2.67742187448669d-017*sqrt(temperature);	! 3A1N + 1A1N
K(4,6) = K(6,4);
K(6,5) = 2.55235015071227d-017*sqrt(temperature);	! 3A1N + 2A1N
K(5,6) = K(6,5);
K(6,6) = 2.57047339860308d-017*sqrt(temperature);	! 3A1N + 3A1N
K(7,1) = 2.49321274842930d-017*sqrt(temperature);	! 2A2N + 1A
K(1,7) = K(7,1);
K(7,2) = 2.46224102249420d-017*sqrt(temperature);	! 2A2N + 2A
K(2,7) = K(7,2);
K(7,4) = 2.62717025408809d-017*sqrt(temperature);	! 2A2N + 1A1N
K(4,7) = K(7,4);
K(7,5) = 2.57143413276816d-017*sqrt(temperature);	! 2A2N + 2A1N
K(5,7) = K(7,5);
K(7,6) = 2.63523697119182d-017*sqrt(temperature);	! 2A2N + 3A1N
K(6,7) = K(7,6);
K(7,7) = 2.66687529977200d-017*sqrt(temperature);	! 2A2N + 2A2N
K(8,1) = 2.67745025144024d-017*sqrt(temperature);	! 3A2N + 1A
K(1,8) = K(8,1);
K(8,2) = 2.54485029577225d-017*sqrt(temperature);	! 3A2N + 2A
K(2,8) = K(8,2);
K(8,3) = 4.79536376022719d-017*sqrt(temperature);	! 3A2N + 1N
K(3,8) = K(8,3);
K(8,4) = 2.78850667019684d-017*sqrt(temperature);	! 3A2N + 1A1N
K(4,8) = K(8,4);
K(8,5) = 2.63927178942724d-017*sqrt(temperature);	! 3A2N + 2A1N
K(5,8) = K(8,5);
K(8,6) = 2.64519112843184d-017*sqrt(temperature);	! 3A2N + 3A1N
K(6,8) = K(8,6);
K(8,7) = 2.72044050159376d-017*sqrt(temperature);	! 3A2N + 2A2N
K(7,8) = K(8,7);
K(8,8) = 2.71878691622467d-017*sqrt(temperature);	! 3A2N + 3A2N
K(9,2) = 2.65066489103025d-017*sqrt(temperature);	! 4A2N + 2A
K(2,9) = K(9,2);
K(9,3) = 5.28620037657595d-017*sqrt(temperature);	! 4A2N + 1N
K(3,9) = K(9,3);
K(9,4) = 2.95750109795969d-017*sqrt(temperature);	! 4A2N + 1A1N
K(4,9) = K(9,4);
K(9,5) = 2.73516064944408d-017*sqrt(temperature);	! 4A2N + 2A1N
K(5,9) = K(9,5);
K(9,6) = 2.69726220089580d-017*sqrt(temperature);	! 4A2N + 3A1N
K(6,9) = K(9,6);
K(9,7) = 2.80666110366871d-017*sqrt(temperature);	! 4A2N + 2A2N
K(7,9) = K(9,7);
K(9,8) = 2.76318549168896d-017*sqrt(temperature);	! 4A2N + 3A2N
K(8,9) = K(9,8);
K(9,9) = 2.77649777952540d-017*sqrt(temperature);	! 4A2N + 4A2N
K(10,1) = 2.78571653410085d-017*sqrt(temperature);	! 3A3N + 1A
K(1,10) = K(10,1);
K(10,2) = 2.62739328560071d-017*sqrt(temperature);	! 3A3N + 2A
K(2,10) = K(10,2);
K(10,4) = 2.89271405648393d-017*sqrt(temperature);	! 3A3N + 1A1N
K(4,10) = K(10,4);
K(10,5) = 2.72013113222697d-017*sqrt(temperature);	! 3A3N + 2A1N
K(5,10) = K(10,5);
K(10,6) = 2.71417804563259d-017*sqrt(temperature);	! 3A3N + 3A1N
K(6,10) = K(10,6);
K(10,7) = 2.79956764945603d-017*sqrt(temperature);	! 3A3N + 2A2N
K(7,10) = K(10,7);
K(10,8) = 2.78663262833214d-017*sqrt(temperature);	! 3A3N + 3A2N
K(8,10) = K(10,8);
K(10,9) = 2.82353824429649d-017*sqrt(temperature);	! 3A3N + 4A2N
K(9,10) = K(10,9);
K(10,10) = 2.85332506980185d-017*sqrt(temperature);	! 3A3N + 3A3N
K(11,1) = 2.96170172731738d-017*sqrt(temperature);	! 4A3N + 1A
K(1,11) = K(11,1);
K(11,2) = 2.72637366871267d-017*sqrt(temperature);	! 4A3N + 2A
K(2,11) = K(11,2);
K(11,3) = 5.50323276483635d-017*sqrt(temperature);	! 4A3N + 1N
K(3,11) = K(11,3);
K(11,4) = 3.05264747323826d-017*sqrt(temperature);	! 4A3N + 1A1N
K(4,11) = K(11,4);
K(11,5) = 2.80948410170128d-017*sqrt(temperature);	! 4A3N + 2A1N
K(5,11) = K(11,5);
K(11,6) = 2.76115422785223d-017*sqrt(temperature);	! 4A3N + 3A1N
K(6,11) = K(11,6);
K(11,7) = 2.87955907067646d-017*sqrt(temperature);	! 4A3N + 2A2N
K(7,11) = K(11,7);
K(11,8) = 2.82619059510194d-017*sqrt(temperature);	! 4A3N + 3A2N
K(8,11) = K(11,8);
K(11,9) = 2.83301566570214d-017*sqrt(temperature);	! 4A3N + 4A2N
K(9,11) = K(11,9);
K(11,10) = 2.88564385947654d-017*sqrt(temperature);	! 4A3N + 3A3N
K(10,11) = K(11,10);
K(11,11) = 2.88886077054498d-017*sqrt(temperature);	! 4A3N + 4A3N
K(12,2) = 2.83315137707714d-017*sqrt(temperature);	! 5A3N + 2A
K(2,12) = K(12,2);
K(12,3) = 5.94007958257262d-017*sqrt(temperature);	! 5A3N + 1N
K(3,12) = K(12,3);
K(12,4) = 3.21189862867815d-017*sqrt(temperature);	! 5A3N + 1A1N
K(4,12) = K(12,4);
K(12,5) = 2.90894045782561d-017*sqrt(temperature);	! 5A3N + 2A1N
K(5,12) = K(12,5);
K(12,6) = 2.82561289488665d-017*sqrt(temperature);	! 5A3N + 3A1N
K(6,12) = K(12,6);
K(12,7) = 2.97186170159675d-017*sqrt(temperature);	! 5A3N + 2A2N
K(7,12) = K(12,7);
K(12,8) = 2.88511464036123d-017*sqrt(temperature);	! 5A3N + 3A2N
K(8,12) = K(12,8);
K(12,9) = 2.86745707849676d-017*sqrt(temperature);	! 5A3N + 4A2N
K(9,12) = K(12,9);
K(12,10) = 2.93914831862297d-017*sqrt(temperature);	! 5A3N + 3A3N
K(10,12) = K(12,10);
K(12,11) = 2.91867228270076d-017*sqrt(temperature);	! 5A3N + 4A3N
K(11,12) = K(12,11);
K(12,12) = 2.92938672270813d-017*sqrt(temperature);	! 5A3N + 5A3N
K(13,1) = 3.05632763782947d-017*sqrt(temperature);	! 4A4N + 1A
K(1,13) = K(13,1);
K(13,2) = 2.79853668847300d-017*sqrt(temperature);	! 4A4N + 2A
K(2,13) = K(13,2);
K(13,3) = 5.71221366785979d-017*sqrt(temperature);	! 4A4N + 1N
K(3,13) = K(13,3);
K(13,4) = 3.14368165817671d-017*sqrt(temperature);	! 4A4N + 1A1N
K(4,13) = K(13,4);
K(13,5) = 2.88025060419474d-017*sqrt(temperature);	! 4A4N + 2A1N
K(5,13) = K(13,5);
K(13,6) = 2.82172621443084d-017*sqrt(temperature);	! 4A4N + 3A1N
K(6,13) = K(13,6);
K(13,7) = 2.94889890466117d-017*sqrt(temperature);	! 4A4N + 2A2N
K(7,13) = K(13,7);
K(13,8) = 2.88586830440131d-017*sqrt(temperature);	! 4A4N + 3A2N
K(8,13) = K(13,8);
K(13,9) = 2.88633722375958d-017*sqrt(temperature);	! 4A4N + 4A2N
K(9,13) = K(13,9);
K(13,10) = 2.94441857242271d-017*sqrt(temperature);	! 4A4N + 3A3N
K(10,13) = K(13,10);
K(13,11) = 2.94150422833882d-017*sqrt(temperature);	! 4A4N + 4A3N
K(11,13) = K(13,11);
K(13,12) = 2.96677019580232d-017*sqrt(temperature);	! 4A4N + 5A3N
K(12,13) = K(13,12);
K(13,13) = 2.99346631156775d-017*sqrt(temperature);	! 4A4N + 4A4N
K(14,2) = 2.90064222637761d-017*sqrt(temperature);	! 5A4N + 2A
K(2,14) = K(14,2);
K(14,3) = 6.13484439809026d-017*sqrt(temperature);	! 5A4N + 1N
K(3,14) = K(14,3);
K(14,4) = 3.29680962444622d-017*sqrt(temperature);	! 5A4N + 1A1N
K(4,14) = K(14,4);
K(14,5) = 2.97519821385214d-017*sqrt(temperature);	! 5A4N + 2A1N
K(5,14) = K(14,5);
K(14,6) = 2.88258174632520d-017*sqrt(temperature);	! 5A4N + 3A1N
K(6,14) = K(14,6);
K(14,7) = 3.03686306118240d-017*sqrt(temperature);	! 5A4N + 2A2N
K(7,14) = K(14,7);
K(14,8) = 2.94132673785796d-017*sqrt(temperature);	! 5A4N + 3A2N
K(8,14) = K(14,8);
K(14,9) = 2.91793447623058d-017*sqrt(temperature);	! 5A4N + 4A2N
K(9,14) = K(14,9);
K(14,10) = 2.99459556920187d-017*sqrt(temperature);	! 5A4N + 3A3N
K(10,14) = K(14,10);
K(14,11) = 2.96859305407371d-017*sqrt(temperature);	! 5A4N + 4A3N
K(11,14) = K(14,11);
K(14,12) = 2.97524119773608d-017*sqrt(temperature);	! 5A4N + 5A3N
K(12,14) = K(14,12);
K(14,13) = 3.01613054079603d-017*sqrt(temperature);	! 5A4N + 4A4N
K(13,14) = K(14,13);
K(14,14) = 3.02064598732666d-017*sqrt(temperature);	! 5A4N + 5A4N
K(15,1) = 3.14733987259492d-017*sqrt(temperature);	! 4A5N + 1A
K(1,15) = K(15,1);
K(15,2) = 2.86762697093825d-017*sqrt(temperature);	! 4A5N + 2A
K(2,15) = K(15,2);
K(15,4) = 3.23112784969842d-017*sqrt(temperature);	! 4A5N + 1A1N
K(4,15) = K(15,4);
K(15,5) = 2.94794044869865d-017*sqrt(temperature);	! 4A5N + 2A1N
K(5,15) = K(15,5);
K(15,6) = 2.87944585924147d-017*sqrt(temperature);	! 4A5N + 3A1N
K(6,15) = K(15,6);
K(15,7) = 3.01516645245239d-017*sqrt(temperature);	! 4A5N + 2A2N
K(7,15) = K(15,7);
K(15,8) = 2.94269160095049d-017*sqrt(temperature);	! 4A5N + 3A2N
K(8,15) = K(15,8);
K(15,9) = 2.93693150533855d-017*sqrt(temperature);	! 4A5N + 4A2N
K(9,15) = K(15,9);
K(15,10) = 3.00033993924085d-017*sqrt(temperature);	! 4A5N + 3A3N
K(10,15) = K(15,10);
K(15,11) = 2.99141933719303d-017*sqrt(temperature);	! 4A5N + 4A3N
K(11,15) = K(15,11);
K(15,12) = 3.01222350034258d-017*sqrt(temperature);	! 4A5N + 5A3N
K(12,15) = K(15,12);
K(15,13) = 3.04270062138490d-017*sqrt(temperature);	! 4A5N + 4A4N
K(13,15) = K(15,13);
K(15,14) = 3.06102318278758d-017*sqrt(temperature);	! 4A5N + 5A4N
K(14,15) = K(15,14);
K(15,15) = 3.09125602332556d-017*sqrt(temperature);	! 4A5N + 4A5N
K(16,1) = 3.30760255827228d-017*sqrt(temperature);	! 5A5N + 1A
K(1,16) = K(16,1);
K(16,2) = 2.96574245481000d-017*sqrt(temperature);	! 5A5N + 2A
K(2,16) = K(16,2);
K(16,4) = 3.37891719764573d-017*sqrt(temperature);	! 5A5N + 1A1N
K(4,16) = K(16,4);
K(16,5) = 3.03906178400585d-017*sqrt(temperature);	! 5A5N + 2A1N
K(5,16) = K(16,5);
K(16,6) = 2.93733589882114d-017*sqrt(temperature);	! 5A5N + 3A1N
K(6,16) = K(16,6);
K(16,7) = 3.09947293401889d-017*sqrt(temperature);	! 5A5N + 2A2N
K(7,16) = K(16,7);
K(16,8) = 2.99532120462343d-017*sqrt(temperature);	! 5A5N + 3A2N
K(8,16) = K(16,8);
K(16,9) = 2.96629507065426d-017*sqrt(temperature);	! 5A5N + 4A2N
K(9,16) = K(16,9);
K(16,10) = 3.04782465891510d-017*sqrt(temperature);	! 5A5N + 3A3N
K(10,16) = K(16,10);
K(16,11) = 3.01639456993767d-017*sqrt(temperature);	! 5A5N + 4A3N
K(11,16) = K(16,11);
K(16,12) = 3.01904182963712d-017*sqrt(temperature);	! 5A5N + 5A3N
K(12,16) = K(16,12);
K(16,13) = 3.06337045226974d-017*sqrt(temperature);	! 5A5N + 4A4N
K(13,16) = K(16,13);
K(16,14) = 3.06399498472483d-017*sqrt(temperature);	! 5A5N + 5A4N
K(14,16) = K(16,14);
K(16,15) = 3.10770227047524d-017*sqrt(temperature);	! 5A5N + 4A5N
K(15,16) = K(16,15);
K(16,16) = 3.10689119145310d-017*sqrt(temperature);	! 5A5N + 5A5N

end subroutine get_coll




subroutine get_evap_mcmc(theta,E,K,temperature)
	implicit none
	real(kind(1.d0)) :: theta(28), E(16,16), K(16,16), temperature

	! evaporation coefficients
	!  Converting thermodynamics data from temperature 298.15 K to 278 K.

	E = 0.d0
	E(1,1) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((theta(1)/temperature-(theta(1)*theta(2))/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(1,1)	! 2A -> 1A + 1A
	E(3,1) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(3)/temperature-(theta(3)*theta(4))/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(3,1)	! 1A1N -> 1N + 1A
	E(1,3) = E(3,1)
	E(4,1) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(5)/temperature-(theta(5)*theta(6))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(4,1)	! 2A1N -> 1A1N + 1A
	E(1,4) = E(4,1)
	E(3,2) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(5)/temperature-(theta(5)*theta(6))/1.d3)&
			 &-0.d0&
			 &-(theta(1)/temperature-(theta(1)*theta(2))/1.d3))&
			 &*5.03218937158374d+002)*K(3,2)	! 2A1N -> 1N + 2A
	E(2,3) = E(3,2)
	E(5,1) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(7))/temperature-(theta(7)*theta(8))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(5,1)	! 3A1N -> 2A1N + 1A
	E(1,5) = E(5,1)
	E(4,2) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(7))/temperature-(theta(7)*theta(8))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3)&
			 &-(theta(1)/temperature-(theta(1)*theta(2))/1.d3))&
			 &*5.03218937158374d+002)*K(4,2)	! 3A1N -> 1A1N + 2A
	E(2,4) = E(4,2)
	E(5,3) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(5,3)	! 2A2N -> 2A1N + 1N
	E(3,5) = E(5,3)
	E(4,4) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((theta(9)/temperature-(theta(9)*theta(10))/1.d3) &
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(4,4)	! 2A2N -> 1A1N + 1A1N
	E(7,1) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(7,1)	! 3A2N -> 2A2N + 1A
	E(1,7) = E(7,1)
	E(6,3) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-((theta(7))/temperature-(theta(7)*theta(8))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(6,3)	! 3A2N -> 3A1N + 1N
	E(3,6) = E(6,3)
	E(5,4) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(5,4)	! 3A2N -> 2A1N + 1A1N
	E(4,5) = E(5,4)
	E(8,1) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(13))/temperature-(theta(13)*theta(14))/1.d3)&
			 &-((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(8,1)	! 4A2N -> 3A2N + 1A
	E(1,8) = E(8,1)
	E(7,2) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(13))/temperature-(theta(13)*theta(14))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-(theta(1)/temperature-(theta(1)*theta(2))/1.d3))&
			 &*5.03218937158374d+002)*K(7,2)	! 4A2N -> 2A2N + 2A
	E(2,7) = E(7,2)
	E(6,4) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(13))/temperature-(theta(13)*theta(14))/1.d3)&
			 &-((theta(7))/temperature-(theta(7)*theta(8))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(6,4)	! 4A2N -> 3A1N + 1A1N
	E(4,6) = E(6,4)
	E(5,5) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &(((theta(13))/temperature-(theta(13)*theta(14))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3))&
			 &*5.03218937158374d+002)*K(5,5)	! 4A2N -> 2A1N + 2A1N
	E(8,3) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(8,3)	! 3A3N -> 3A2N + 1N
	E(3,8) = E(8,3)
	E(7,4) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(7,4)	! 3A3N -> 2A2N + 1A1N
	E(4,7) = E(7,4)
	E(10,1) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-(theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(10,1)	! 4A3N -> 3A3N + 1A
	E(1,10) = E(10,1)
	E(9,3) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-((theta(13))/temperature-(theta(13)*theta(14))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(9,3)	! 4A3N -> 4A2N + 1N
	E(3,9) = E(9,3)
	E(8,4) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(8,4)	! 4A3N -> 3A2N + 1A1N
	E(4,8) = E(8,4)
	E(7,5) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3))&
			 &*5.03218937158374d+002)*K(7,5)	! 4A3N -> 2A2N + 2A1N
	E(5,7) = E(7,5)
	E(11,1) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(19))/temperature-(theta(19)*theta(20))/1.d3)&
			 &-((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(11,1)	! 5A3N -> 4A3N + 1A
	E(1,11) = E(11,1)
	E(10,2) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(19))/temperature-(theta(19)*theta(20))/1.d3)&
			 &-(theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-(theta(1)/temperature-(theta(1)*theta(2))/1.d3))&
			 &*5.03218937158374d+002)*K(10,2)	! 5A3N -> 3A3N + 2A
	E(2,10) = E(10,2)
	E(9,4) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(19))/temperature-(theta(19)*theta(20))/1.d3)&
			 &-((theta(13))/temperature-(theta(13)*theta(14))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(9,4)	! 5A3N -> 4A2N + 1A1N
	E(4,9) = E(9,4)
	E(8,5) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(19))/temperature-(theta(19)*theta(20))/1.d3)&
			 &-((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3))&
			 &*5.03218937158374d+002)*K(8,5)	! 5A3N -> 3A2N + 2A1N
	E(5,8) = E(8,5)
	E(7,6) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(19))/temperature-(theta(19)*theta(20))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-((theta(7))/temperature-(theta(7)*theta(8))/1.d3))&
			 &*5.03218937158374d+002)*K(7,6)	! 5A3N -> 2A2N + 3A1N
	E(6,7) = E(7,6)
	E(11,3) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(21)/temperature-(theta(21)*theta(22))/1.d3)&
			 &-((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(11,3)	! 4A4N -> 4A3N + 1N
	E(3,11) = E(11,3)
	E(10,4) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(21)/temperature-(theta(21)*theta(22))/1.d3)&
			 &-(theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(10,4)	! 4A4N -> 3A3N + 1A1N
	E(4,10) = E(10,4)
	E(7,7) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((theta(21)/temperature-(theta(21)*theta(22))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3))&
			 &*5.03218937158374d+002)*K(7,7)	! 4A4N -> 2A2N + 2A2N
	E(13,1) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(23))/temperature-(theta(23)*theta(24))/1.d3)&
			 &-(theta(21)/temperature-(theta(21)*theta(22))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(13,1)	! 5A4N -> 4A4N + 1A
	E(1,13) = E(13,1)
	E(12,3) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(23))/temperature-(theta(23)*theta(24))/1.d3)&
			 &-((theta(19))/temperature-(theta(19)*theta(20))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(12,3)	! 5A4N -> 5A3N + 1N
	E(3,12) = E(12,3)
	E(11,4) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(23))/temperature-(theta(23)*theta(24))/1.d3)&
			 &-((theta(17))/temperature-(theta(17)*theta(18))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(11,4)	! 5A4N -> 4A3N + 1A1N
	E(4,11) = E(11,4)
	E(10,5) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(23))/temperature-(theta(23)*theta(24))/1.d3)&
			 &-(theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-(theta(5)/temperature-(theta(5)*theta(6))/1.d3))&
			 &*5.03218937158374d+002)*K(10,5)	! 5A4N -> 3A3N + 2A1N
	E(5,10) = E(10,5)
	E(8,7) = 7.33893243358348d+027/temperature*exp(&
			 &(((theta(23))/temperature-(theta(23)*theta(24))/1.d3)&
			 &-((theta(11))/temperature-(theta(11)*theta(12))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3))&
			 &*5.03218937158374d+002)*K(8,7)	! 5A4N -> 3A2N + 2A2N
	E(7,8) = E(8,7)
	E(13,3) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(25)/temperature-(theta(25)*theta(26))/1.d3)&
			 &-(theta(21)/temperature-(theta(21)*theta(22))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(13,3)	! 4A5N -> 4A4N + 1N
	E(3,13) = E(13,3)
	E(15,1) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(27)/temperature-(theta(27)*theta(28))/1.d3)&
			 &-(theta(25)/temperature-(theta(25)*theta(26))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(15,1)	! 5A5N -> 4A5N + 1A
	E(1,15) = E(15,1)
	E(14,3) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(27)/temperature-(theta(27)*theta(28))/1.d3)&
			 &-((theta(23))/temperature-(theta(23)*theta(24))/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(14,3)	! 5A5N -> 5A4N + 1N
	E(3,14) = E(14,3)
	E(13,4) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(27)/temperature-(theta(27)*theta(28))/1.d3)&
			 &-(theta(21)/temperature-(theta(21)*theta(22))/1.d3)&
			 &-(theta(3)/temperature-(theta(3)*theta(4))/1.d3))&
			 &*5.03218937158374d+002)*K(13,4)	! 5A5N -> 4A4N + 1A1N
	E(4,13) = E(13,4)
	E(10,7) = 7.33893243358348d+027/temperature*exp(&
			 &((theta(27)/temperature-(theta(27)*theta(28))/1.d3)&
			 &-(theta(15)/temperature-(theta(15)*theta(16))/1.d3)&
			 &-(theta(9)/temperature-(theta(9)*theta(10))/1.d3))&
			 &*5.03218937158374d+002)*K(10,7)	! 5A5N -> 3A3N + 2A2N
	E(7,10) = E(10,7)
end subroutine get_evap_mcmc