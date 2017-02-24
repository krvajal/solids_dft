module energy
    use types
    use constants
    use gvect
    implicit none
    real(dp) :: kinetic_energy
    real(dp) :: short_range_pp_enery
    real(dp) :: local_pp_energy
    real(dp) :: exc_corr_energy
    real(dp) :: hartree_energy
    real(dp) :: local_core_energy
    real(dp) :: self_energy
    real(dp) :: total_energy

    contains
    subroutine  compute_kinetic_energy(num_orbitals,filling_factor,coeffs,g_vects, length)
        implicit none

        integer :: num_orbitals
        real(dp) :: filling_factor(num_orbitals),length
        real(dp) :: coeffs(num_orbitals,num_orbitals)
        integer :: g_vects(3,num_orbitals)
        ! aux variables
        integer :: i, j
        real(dp) :: temp

        kinetic_energy = 0.0_dp
        do i =1,num_orbitals
            temp = 0.0
            do j = 1, num_orbitals
                temp = temp  +  sum(g_vects(:,j)**2) * coeffs(j,i)**2*2*pi/length
            enddo
            kinetic_energy = kinetic_energy + filling_factor(i) * temp
        enddo
        kinetic_energy = kinetic_energy * 0.5_dp

    end subroutine  compute_kinetic_energy
 end module energy