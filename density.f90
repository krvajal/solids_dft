module density
    use types
    use gvect
    use fft
    use pseudopot
    use constants
    implicit none

    
complex(dp),allocatable :: density_r(:,:,:) !density in real space
complex(dp),allocatable :: density_g(:,:,:) !density in reciprocal space
complex(dp),allocatable :: core_density_g(:,:,:), total_density_g(:,:,:)
real(dp),allocatable ::  FillingFactor(:)


contains
    subroutine init_density(Nx,Ny,Nz,num_atoms,numGVects, box_length, params)
        integer, intent(in) :: Nx,Ny,Nz,numGVects
        integer, intent(in) :: num_atoms
        real(dp) :: box_length
        type(GthPotParams),intent(in) :: params(num_atoms)
        real(dp) :: chi, two_pi_over_l
        integer :: i

     
        two_pi_over_l = 2 * pi / box_length
        allocate(density_r(0:Nx-1,0:Ny-1,0:Nz-1))
        allocate(density_g(0:Nx-1,0:Ny-1,0:Nz-1))
        allocate(core_density_g(0:Nx-1,0:Ny-1,0:Nz-1))
        allocate(total_density_g(0:Nx-1,0:Ny-1,0:Nz-1))

        density_r = 0.0_dp
        density_g = 0.0_dp
        allocate(FillingFactor(numGVects))
        FillingFactor = 0
        core_density_g = 0
        do i = 1,num_atoms
            chi = params(i)%chi 
            core_density_g = core_density_g  - params(i)%Zeff/params(i)%omega&
                             * exp(-0.5 *(g_grid_norm * chi * two_pi_over_l)**2 ) *structure_factor(i,:,:,:)
        enddo
        
        
    end subroutine init_density


function layoutKIndexForFft(Nx,Ny,Nz,num_orbitals,C) result(retval)
!-----------------------------------------------------------    
    integer :: num_orbitals, Nx,Ny,Nz
    complex(dp) :: C(num_orbitals)
    complex(dp) ::  retval(0:Nx-1,0:Ny-1,0:Nz-1)
    integer :: size , i
    integer :: toI,toJ, toK
    retval = 0.0 !init to zero


    do i = 1, num_orbitals
         toI = g_indexes(1,i)
         if(toI < 0) toI = toI + Nx
         toJ = g_indexes(2,i)
         if(toJ < 0) toJ = toJ + Ny
         toK = g_indexes(3,i)
         if(toK < 0) toK = toK + Nz
            !print *, toI,toJ,toK, C(i)
         retval(toI,toJ,toK) = C(i)
        ! print *,g_grid(:,toI,toJ,toK), g_indexes(:,i)     
    enddo
    ! ok, checkeado
end function layoutKIndexForFft
    !---------------------------------------------
    ! compute the system density
    ! in real and complex espace
    ! upon the coefficients
    ! of the expansion of the wave
    ! functions in plane waves in
    ! reciprocal space
    ! --------------------------------------------
    subroutine compute_density(N,C,omega)
        implicit none
        integer,intent(in) :: N !number of G vectors
        complex(dp) :: C(N,N)  ! the wave functions expansion coeffs, C(:,j) for the jth one
        real(dp) :: omega   ! the unit cell volumen

        real(dp),allocatable :: Nj(:,:,:) 

        ! aux variables for calculation
        real(dp) :: val
        complex(dp), allocatable :: psi_r(:,:,:)
        complex(dp), allocatable :: psi_k(:,:,:)

        integer :: i,k,l,m
        allocate(psi_r(0:Nx-1,0:Ny-1,0:Nz-1))  
        density_r = 0.0_dp
     
        do i = 1,numGVects
          
            psi_k = layoutKIndexForFft(Nx,Ny,Nz,numGVects,C(:,i))
       
            call fft_forward_3d(Nx,Ny,Nz,psi_k,psi_r)     ! exp(1)  

            psi_r = psi_r/sqrt(omega)
    
            density_r =  density_r + FillingFactor(i)*psi_r*conjg(psi_r)
        
        enddo

        call fft_backward_3d(Nx,Ny,Nz,density_r,density_g)      ! exp(-1)  

        density_g = density_g/(Nx*Ny*Nz)    
    end subroutine compute_density


!--------------------------------------
! map to linear layout, unimplemented
subroutine to_linear_layout(Nx,Ny,Nz, in, out)
!----------------------------------------
    integer, intent(in) :: Nx,Ny,Nz
    complex(dp), intent(in) :: in(Nx,Ny,Nz)
    complex(dp), intent(out) :: out(Nx,Ny,Nz)
end subroutine to_linear_layout

end module density