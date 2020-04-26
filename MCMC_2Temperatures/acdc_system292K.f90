module acdc_system

implicit none

integer, parameter :: nclust = 16						! number of clusters, molecules and ions
integer, parameter :: neq = 16							! number of equations

real(kind(1.d0)), parameter :: temp = 278.d0				! temperature in K
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+000			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 1.00000000000000d-003

integer, parameter :: n1A = 1, n1N = 3			! cluster indices for monomers and ions

integer, parameter :: n_mol_types = 2
integer, parameter :: nmolA = 1, nmolN = 2			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 2				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 16			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 0			! number of negative molecules and clusters
integer, parameter :: n_positives = 0			! number of positive molecules and clusters

integer, parameter :: n_neutral_clusters = 14			! number of neutral clusters
integer, parameter :: n_negative_clusters = 0			! number of negative clusters
integer, parameter :: n_positive_clusters = 0			! number of positive clusters

real(kind(1.d0)), parameter :: mass_max = 575.60
real(kind(1.d0)), parameter :: diameter_max = 1.07
real(kind(1.d0)), parameter :: mob_diameter_max = 1.37

integer, parameter :: nmols_out_neutral(1, 2) = reshape((/6, 5/),(/1, 2/))			! criteria for outgrowing neutrals


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(16)

	n_A = (/1, 2, 0, 1, 2, 3, 2, 3, 4, 3, &
		&4, 5, 4, 5, 4, 5/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(2)

	cluster_numbers = (/1, 4/)

end subroutine clusters_with_1_A

subroutine cluster_arrays(neutral_clusters)
	implicit none
	integer :: neutral_clusters(14)

	neutral_clusters = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)

end subroutine cluster_arrays

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(16)

	mass = (/98.08, 196.16, 17.04, 115.12, 213.20, 311.28, 230.24, 328.32, 426.40, 345.36, &
		&443.44, 541.52, 460.48, 558.56, 477.52, 575.60/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(16)

	 diameter = (/0.55, 0.70, 0.43, 0.63, 0.75, 0.84, 0.79, 0.87, 0.94, 0.91, &
		&0.97, 1.03, 1.00, 1.05, 1.02, 1.07/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(16)

	 mob_diameter = (/0.85, 1.00, 0.73, 0.93, 1.05, 1.14, 1.09, 1.17, 1.24, 1.21, &
		&1.27, 1.33, 1.30, 1.35, 1.32, 1.37/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
	implicit none
	character(len=11), dimension(16) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '1N'
	clust(4)(:) = '1A1N'
	clust(5)(:) = '2A1N'
	clust(6)(:) = '3A1N'
	clust(7)(:) = '2A2N'
	clust(8)(:) = '3A2N'
	clust(9)(:) = '4A2N'
	clust(10)(:) = '3A3N'
	clust(11)(:) = '4A3N'
	clust(12)(:) = '5A3N'
	clust(13)(:) = '4A4N'
	clust(14)(:) = '5A4N'
	clust(15)(:) = '4A5N'
	clust(16)(:) = '5A5N'

end subroutine cluster_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(2) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'N'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(2)

	n_monomers = (/1, 3/)

end subroutine monomer_indices


end module acdc_system

