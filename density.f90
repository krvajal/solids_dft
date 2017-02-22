module density
    use types
    use gvect
    implicit none

    
complex(dp),allocatable :: NtotalR(:,:,:) !density in real space
complex(dp),allocatable :: NtotalK(:,:,:) !density in reciprocal space
real(dp),allocatable ::  FillingFactor(:)

contains
    subroutine init_density(Nx,Ny,Nz,numGVects)
        integer, intent(in) :: Nx,Ny,Nz,numGVects

        allocate(NtotalR(Nx,Ny,Nz))
        allocate(NtotalK(Nx,Ny,Nz))
        NtotalR = 0.0_dp
        NtotalK = 0.0_dp
        allocate(FillingFactor(numGVects))
        FillingFactor = 0
    end subroutine init_density
end module density