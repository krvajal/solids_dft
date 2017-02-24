module density
    use types
    use gvect
    use fft
    implicit none

    
complex(dp),allocatable :: density_r(:,:,:) !density in real space
complex(dp),allocatable :: density_g(:,:,:) !density in reciprocal space
real(dp),allocatable ::  FillingFactor(:)

contains
    subroutine init_density(Nx,Ny,Nz,numGVects)
        integer, intent(in) :: Nx,Ny,Nz,numGVects

        allocate(density_r(0:Nx-1,0:Ny-1,0:Nz-1))
        allocate(density_g(0:Nx-1,0:Ny-1,0:Nz-1))
        density_r = 0.0_dp
        density_g = 0.0_dp
        allocate(FillingFactor(numGVects))
        FillingFactor = 0
    end subroutine init_density


    function layoutKIndexForFft(Nx,Ny,Nz,num_orbitals,C) result(retval)
        integer :: num_orbitals, Nx,Ny,Nz
        real(dp) :: C(num_orbitals)
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
             retval(toI,toJ,toK) = complex(C(i),0.0_dp)
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
        real(dp) :: C(N,N)  ! the wave functions expansion coeffs, C(:,j) for thr jth one
        real(dp) :: omega   ! the unit cell volumen

        real(dp),allocatable :: Nj(:,:,:) 
        ! aux variables for calculation

        complex(dp), allocatable :: psi_r(:,:,:)
        complex(dp), allocatable :: psi_k(:,:,:)

        integer :: i,k,l,m
        allocate(psi_r(0:Nx-1,0:Ny-1,0:Nz-1))  
        density_r = 0.0_dp
        print *,"ok"
        do i = 1,numGVects
          
            psi_k = layoutKIndexForFft(Nx,Ny,Nz,numGVects,C(:,i))
       
            call fft_backward_3d(Nx,Ny,Nz,psi_k,psi_r)     ! exp(1)  

            psi_r = psi_r/sqrt(omega)
  
            density_r =  density_r + FillingFactor(i)*(psi_r)*conjg(psi_r)
        
        enddo

        call fft_forward_3d(Nx,Ny,Nz,density_r,density_g)      ! exp(-1)  
        density_g = density_g/(Nx*Ny*Nz)    
    end subroutine compute_density
end module density