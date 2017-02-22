! solves the problem of of a particle with the gth potential
program particle_gth
    use types
    use constants
    use linalg
    use gth_potential
    use gvect 
    use fft
    use xc
    use density
    integer,parameter :: basisDim =  7
    
    integer :: i, j, k, l
    real(dp) :: length = 5.0_dp
    real(dp) :: kLength
    real(dp) :: ecut = 1.3
    real(dp) :: hamiltMatrix(basisDim, basisDim)
    real(dp) :: energies (basisDim)
    real(dp) :: C(basisDim,basisDim)
    real(dp) :: omega

    omega = length**3
    
    kLength  = 2 * pi/length

    call init_system()
  
    do l = 1, 10

        call fillHamiltMatrix(numGVects, gVects, hamiltMatrix)
        ! solve the eigenproblem
        call eigh(hamiltMatrix, energies,C)
        call computeDensity (numGVects, C)
        call computeVxc()
        
        print *,"Energies = ",energies(1:3)
    enddo
    ! results from thijssen
    ! e1 = −0.03572203
    ! e2 =  0.68175686 
    ! e3 =  0.80555307 (3 times) 
    ! e4 =  0.83735807 (2 times)
    
    

contains


! evaluate the hamiltonian matrix with K + V
! 
subroutine fillHamiltMatrix(basisDim, KBasisSet, hamiltMatrix)
    integer,intent(in) :: basisDim
    real(dp),intent(in) :: KBasisSet(:,:)
    real(dp), intent(out) :: hamiltMatrix(basisDim, basisDim)
    type(GthPotParams) :: paramsHidrogen
    complex(dp),allocatable :: vxc_g (:,:,:), vxc_aux(:,:,:);
    real(dp) :: vhartee
    integer :: deltaG(3) ! G - G'
    ! =================================
    ! hidrogen parameters
    ! =================================
    paramsHidrogen%c1 = -4.0663326_dp
    paramsHidrogen%c2 = 0.6778322_dp
    paramsHidrogen%chi = 0.2
    paramsHidrogen%Zeff = 1 
    paramsHidrogen%omega = omega ! cell with edge length of 5.0 a.u.
    ! =================================


    allocate(vxc_g(Nx,Ny,Nz))
    allocate(vxc_aux(Nx,Ny,Nz))

    vxc_aux = CMPLX(vxc_r, kind=dp)

    call fftBackward3D(Nx,Ny,Nz, vxc_aux, vxc_g)
    vxc_g = vxc_g/(Nx*Ny*Nz)

    hamiltMatrix = 0 !init to zero the hamiltonina matrix
    
    do i = 1, basisDim
        ! sum the kinetic energy
        hamiltMatrix(i,i) = sum( (2 * pi/length * KBasisSet(:,i))**2) *  0.5_dp !kinetic energy diagonal term

        ! add the pseudopotential
        do j = 1, i-1
            hamiltMatrix(i,j) = hamiltMatrix(i,j) + &
                                gthPotential(paramsHidrogen, &
                                 KBasisSet(:,i) * 2*pi/length,&
                                  KBasisSet(:,j)* 2 * pi/length)
            

            ! add the vxc
            deltaG = KBasisSet(:,i) - KBasisSet(:,j)
            ! print *,"deltaG", deltaG
            ! print *, "Nx",Nx
            ! print *, "Ny",Ny
            ! print *, "Nz",Nz

            ! get the position on the potential
            if (deltaG(1) < 1) deltaG(1) = deltaG(1) + Nx
            if (deltaG(2) < 1) deltaG(2) = deltaG(2) + Ny
            if (deltaG(3) < 1) deltaG(3) = deltaG(3) + Nz
            hamiltMatrix(i,j) = hamiltMatrix(i,j) + vxc_g(deltaG(1),deltaG(2),deltaG(3))
            
            ! add hartree potential
            vhartee = 0
        
            vhartee = NtotalK(deltaG(1),deltaG(2),deltaG(3))*length*length/pi
            vhartee = vhartee / sum((KBasisSet(:,i) - KBasisSet(:,j))**2)
            
            hamiltMatrix(i,j) = hamiltMatrix(i,j)  + vhartee
            ! upper part of matrix
            hamiltMatrix(j,i) = hamiltMatrix(i,j)            
        enddo
    enddo


end subroutine fillHamiltMatrix

subroutine computeDensity(N,C)
    integer,intent(in) :: N !number of G vectors
    real(dp) :: C(N,N)
    real(dp),allocatable :: Nj(:,:,:)
    complex(dp), allocatable :: NjR(:,:,:)
 
    integer :: i,k,l,m
  
    allocate(NjR(Nx,Ny,Nz))  
    NtotalR = 0.0_dp

    do i = 1,numGVects
        C(:,i) = C(:,i)/sum(C(:,i)**2) ! normalize 
        NjR = CMPLX(layoutKIndexForFft(C(:,i)))
        call fftForward3D(Nx,Ny,Nz,NjR,NjR)        
        NjR = NjR/sqrt(omega)
        NtotalR = FillingFactor(i)*ABS(NjR)**2  + NtotalR
    enddo
    call fftBackward3D(Nx,Ny,Nz,NtotalR,NtotalK)        
    NtotalK = NtotalK/(Nx*Ny*Nz)
    

end subroutine computeDensity
    
    subroutine init_system()

        !TODO: read params from file 
        call ggen(length, ecut)
        call init_density(Nx,Ny,Nz,numGVects)
        FillingFactor(1) = 1
        call init_xc()
    
    end subroutine init_system
end program particle_gth
