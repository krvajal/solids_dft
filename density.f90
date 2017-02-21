module density
    use types
    implicit none

integer :: Nx, Ny, Nz    
complex(dp),allocatable :: NtotalR(:,:,:) !density in real space
complex(dp),allocatable :: NtotalK(:,:,:) !density in reciprocal space

contains
    subroutine init_density(Nx,Ny,Nz)
        integer, intent(in) :: Nx,Ny,Nz
        allocate(NtotalR(Nx,Ny,Nz))
        allocate(NtotalK(Nx,Ny,Nz))
        NtotalR = 0.0_dp
        NtotalK = 0.0_dp
    end subroutine init_density
end module density