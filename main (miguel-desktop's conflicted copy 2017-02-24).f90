! solves the problem of of a particle with the gth potential
program main
    use types
    use constants
    use linalg
    use gth_potential
    use gvect 
    use fft
    use xc
    use energy
    use density
    implicit none
    
    integer :: i, j, k, l
    real(dp) ::length = 5.0_dp
    real(dp) :: ecut = 1.3
    real(dp),allocatable :: hamiltMatrix(:, :)
    real(dp),allocatable :: energies (:)
    real(dp),allocatable :: psi_coeffs_g(:,:)
    real(dp) :: omega
    real(dp),allocatable :: vhartee_r(:,:,:)

    omega = length**3

    call init_system()
    call khon_sham_loop()
  
    ! results from thijssen
    ! e1 = −0.03572203
    ! e2 =  0.68175686 
    ! e3 =  0.80555307 (3 times) 
    ! e4 =  0.83735807 (2 times)
    

contains

subroutine init_system()

    !TODO: read params from file 
    call ggen(length, ecut)
    allocate(vhartee_r(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(hamiltMatrix(numGVects,numGVects))
    allocate(energies(numGVects))
    allocate(psi_coeffs_g(numGVects,numGVects))
    print *, "GridDim:", Nx,Ny,Nz
    print *, "num_orbitals:", numGVects
    print *, g_indexes(:,1)
    call init_density(Nx,Ny,Nz,numGVects)
    FillingFactor(1) = 1 ! only one atom
    call init_xc()
end subroutine init_system

subroutine khon_sham_loop()

      type(GthPotParams) :: paramsHidrogen
        ! hidrogen parameters
        ! =================================
        paramsHidrogen%c1 = -4.0663326_dp
        paramsHidrogen%c2 = 0.6778322_dp
        paramsHidrogen%chi = 0.2_dp
        paramsHidrogen%Zeff = 1 
        paramsHidrogen%omega = omega ! cell with edge length of 5.0 a.u.
        paramsHidrogen%box_length = length
       ! =================================
      do l = 1, 4

        call fill_hamilt_matrix(paramsHidrogen,numGVects, g_indexes, hamiltMatrix)
        ! solve the eigenproblem
        call eigh(hamiltMatrix, energies,psi_coeffs_g)
         print *,"Energies = ",energies(1:3)
        ! print *,matmul(hamiltMatrix,C(:,1))/C(:,1)
        
        call compute_density (numGVects, psi_coeffs_g,omega)
        call compute_vxc()
       ! call compute_kinetic_energy(numGVects, FillingFactor, psi_coeffs_g,g_indexes,length)    
       
        ! print *, "kinetic = ", kinetic_energy 

    enddo
end subroutine khon_sham_loop

! evaluate the hamiltonian matrix with K + V
! 
subroutine fill_hamilt_matrix(params,num_orbitals, KBasisSet, hamiltMatrix)
    implicit none
    type(GthPotParams) :: params
    integer,intent(in) :: num_orbitals
    integer,intent(in) :: KBasisSet(3,num_orbitals)
    real(dp), intent(out) :: hamiltMatrix(num_orbitals, num_orbitals)
    
    complex(dp),allocatable :: vxc_g (:,:,:), vxc_aux(:,:,:);
    real(dp) :: vhartee
    integer :: deltaG(3) ! G - G'
    real(dp) :: delta_g2 ! delta g squared
    ! =================================

    vhartee_r = 0;

    allocate(vxc_g(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(vxc_aux(0:Nx-1,0:Ny-1,0:Nz-1))

    vxc_aux = CMPLX(vxc_r, kind=dp)

    call fft_forward_3d(Nx,Ny,Nz, vxc_aux, vxc_g)
    vxc_g = vxc_g/(Nx*Ny*Nz)

    hamiltMatrix = 0 !init to zero the hamiltonina matrix
    ! density_g = density_g - 1/(omega*sqrt(4*pi))
    do i = 1, num_orbitals
        ! sum the kinetic energy, this term does not change and can be precomputed
        hamiltMatrix(i,i) = sum( (2* pi/length * KBasisSet(:,i))**2) /2 !kinetic energy diagonal term
        ! add the pseudopotential
        do j = 1,num_orbitals
            hamiltMatrix(i,j) = hamiltMatrix(i,j) + &
                                pseudo_pot_gth(params, &
                                 KBasisSet(:,i)  * 2*pi/length,&
                                  KBasisSet(:,j) * 2*pi/length)
            
            
            ! add the vxc
            deltaG = KBasisSet(:,i) - KBasisSet(:,j)
            delta_g2 = sum(deltaG**2)
            ! print *,"deltaG", deltaG
            ! print *, "Nx",Nx
            ! print *, "Ny",Ny
            ! print *, "Nz",Nz
            
            ! get the position on the potential
            if (deltaG(1) < 0) deltaG(1) =  deltaG(1) + Nx 
            if (deltaG(2) < 0) deltaG(2) =  deltaG(2) + Ny
            if (deltaG(3) < 0) deltaG(3) = deltaG(3) + Nz
            !hamiltMatrix(i,j) = hamiltMatrix(i,j) + vxc_g(deltaG(1),deltaG(2),deltaG(3))
            ! print *, "density",density_g(deltaG(1),deltaG(2),deltaG(3))

            !==================================================================
            ! Hartree term
            !================================================================
            if (i /= j) then

                vhartee = density_g(deltaG(1),deltaG(2),deltaG(3))* length**2 /( pi*delta_g2)
                vhartee_r(deltaG(1),deltaG(2),deltaG(3)) = vhartee;
                ! print *, "hartree",vhartee
                hamiltMatrix(i,j) = hamiltMatrix(i,j) + vhartee
            endif 
            !=================================================================
            ! upper part of matrix
            
        enddo
    enddo
    ! do i = 0,Nx-1
    !     print *,vhartee_r(0,0,i)
    !     ! print *,vhartee_r(0,:,0)
    !     ! print *,vhartee_r(:,0,0)
    ! enddo
end subroutine fill_hamilt_matrix


    


end program main
