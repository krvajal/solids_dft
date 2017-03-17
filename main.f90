! solves the problem of of a particle with the gth potential
program autoconsistente
    use types
    use constants
    use linalg
    use pseudopot
    use gvect 
    use fft
    use xc

    use density
    implicit none
    

    character(len = 23) :: argument
    integer :: i, j, k, l
    real(dp),parameter ::length = 5.0_dp
    real(dp) :: ecut
    complex(dp),allocatable :: hamiltMatrix(:, :)
    real(dp),allocatable :: energies (:)
    complex(dp),allocatable :: psi_coeffs_g(:,:)
    real(dp) :: omega
    complex(dp),allocatable :: vhartee_r(:,:,:)
    complex(dp),allocatable :: vxc_g (:,:,:), vxc_aux(:,:,:);
    type(GthPotParams) :: paramsHidrogen
    real(dp) :: norm
    integer :: iter
    type(gth_pp_t) :: pseudopotential


    omega = length**3
    ! hidrogen parameters
    ! =================================
    paramsHidrogen%c1 = -4.0663326_dp
    paramsHidrogen%c2 = 0.6778322_dp
    paramsHidrogen%chi = 0.2_dp
    paramsHidrogen%Zeff = 1 
    paramsHidrogen%omega = omega ! cell with edge length of 5.0 a.u.
    paramsHidrogen%box_length = length
   ! =================================

    

    call init_system()
    print  *,"Starting sc loop"
    call khon_sham_loop()
  
    ! results from thijssen
    ! e1 = −0.03572203
    ! e2 =  0.68175686 
    ! e3 =  0.80555307 (3 times) 
    ! e4 =  0.83735807 (2 times)
    

contains

subroutine init_system()

    ecut = 1.3
    
    !TODO: read params from file 
    call ggen(length, ecut)
    allocate(vhartee_r(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(hamiltMatrix(numGVects,numGVects))
    allocate(energies(numGVects))
    allocate(psi_coeffs_g(numGVects,numGVects))
    allocate(vxc_g(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(vxc_aux(0:Nx-1,0:Ny-1,0:Nz-1))

    print *, "GridDim:", Nx,Ny,Nz
    print *, "num_orbitals:", numGVects

    call init_density(Nx,Ny,Nz,numGVects,paramsHidrogen)

    FillingFactor(1) = 1 ! only one atom
    call init_xc()

end subroutine init_system

subroutine khon_sham_loop()

   
    density_g = 0.d0
    density_g(0,0,0) = 1.d0
    norm =  Omega*density_g(0,0,0)

    pseudopotential%params = paramsHidrogen

    print *, "density norm", norm
    density_g = density_g/norm

    CALL fft_forward_3d(Nx,Ny,Nz, density_g, density_r)

    do l = 1, 10
        print *, "=============================="
        print *, "Starting KS iteration no ", l
        print *, "=============================="
        call fill_hamilt_matrix(pseudopotential,numGVects, g_indexes, hamiltMatrix)
        ! solve the eigenproblem
        call eigh(hamiltMatrix, energies,psi_coeffs_g)

        ! print *, psi_coeffs_g(:,1)
        
        ! do iter = 1,7
        !     print *, hamiltMatrix(iter,:)
        ! enddo
        
         print *,"eigenvalue = ",energies(1:3)
         
         ! stop
        ! print *,matmul(hamiltMatrix,C(:,1))/C(:,1)
        
        call compute_density (numGVects, psi_coeffs_g,omega)
        call compute_total_energy(pseudopotential)
       ! call compute_kinetic_energy(numGVects, FillingFactor, psi_coeffs_g,g_indexes,length)    
       
        ! print *, "kinetic = ", kinetic_energy 

    enddo
end subroutine khon_sham_loop

! evaluate the hamiltonian matrix with K + V
! 
subroutine fill_hamilt_matrix(pseudpot,num_orbitals, KBasisSet, hamiltMatrix)
    implicit none
    type(gth_pp_t) :: pseudpot
    integer,intent(in) :: num_orbitals
    integer,intent(in) :: KBasisSet(3,num_orbitals)
    complex(dp), intent(out) :: hamiltMatrix(num_orbitals, num_orbitals)
    real(dp) :: two_pi_over_a = 2*pi/length
   
    real(dp) :: vhartee, norm_delta_G
    integer :: delta_g_idx(3) ! G - G'
    real(dp) :: delta_g2 ! delta g squared
    ! =================================
    call compute_vxc()
    vhartee_r = 0;


    vxc_aux = CMPLX(vxc_r, kind=dp)

    call fft_forward_3d(Nx,Ny,Nz, vxc_aux, vxc_g)
    vxc_g = vxc_g/(Nx*Ny*Nz)

   
    hamiltMatrix = 0 !init to zero the hamiltonian matrix
    ! density_g = density_g - 1/(omega*sqrt(4*pi))
    do i = 1, num_orbitals
        ! sum the kinetic energy, this term does not change and can be precomputed
        hamiltMatrix(i,i) = sum( (two_pi_over_a* KBasisSet(:,i))**2) /2 !kinetic energy diagonal term
        ! add the pseudopotential

        do j = 1, num_orbitals

            delta_g_idx = KBasisSet(:,i) - KBasisSet(:,j)
            norm_delta_G = norm2(delta_g_idx* two_pi_over_a)
    
            delta_g2 = norm_delta_G**2
            ! print *,"delta_g_idx", delta_g_idx
            ! print *, "Nx",Nx
            ! print *, "Ny",Ny
            ! print *, "Nz",Nz
            
            ! get the position on the potential
            if (delta_g_idx(1) < 0) delta_g_idx(1) =  delta_g_idx(1) + Nx 
            if (delta_g_idx(2) < 0) delta_g_idx(2) =  delta_g_idx(2) + Ny
            if (delta_g_idx(3) < 0) delta_g_idx(3) = delta_g_idx(3) + Nz


            ! pseudopotential contribution
            hamiltMatrix(i,j) = hamiltMatrix(i,j) + &
                                pseudpot%local(norm_delta_G) * &
                                 structure_factor(delta_g_idx(1),delta_g_idx(2), delta_g_idx(3))




            ! add the vxc

            hamiltMatrix(i,j) = hamiltMatrix(i,j) + vxc_g(delta_g_idx(1),delta_g_idx(2),delta_g_idx(3))
            ! 

            !==================================================================
            ! Hartree term
            !================================================================
            if (i /= j) then

              vhartee = density_g(delta_g_idx(1),delta_g_idx(2),delta_g_idx(3))* 4 * pi /delta_g2        
              hamiltMatrix(i,j) = hamiltMatrix(i,j) + vhartee

            endif 
            
            
        enddo
    enddo
    ! do  i = 1, numGVects
    !     do j = 1, i -1
    !         print *, hamiltMatrix(i,j) - conjg(hamiltMatrix(j,i))
    !     enddo
    ! enddo
    ! stop
    ! the hamiltonian is hermitic
    call fft_forward_3d(Nx,Ny,Nz,vhartee_r,vhartee_r)
    ! do i = 0,Nx-1
        ! print *, i,realpart(vhartee_r(0,0,i)),realpart(vhartee_r(0,i,0)),realpart(vhartee_r(i,0,0))
    !     ! print *,vhartee_r(0,:,0)
    !     ! print *,vhartee_r(:,0,0)
    ! enddo
end subroutine fill_hamilt_matrix

subroutine compute_total_energy(pseudopot)

    type(gth_pp_t),intent(in) :: pseudopot
    real(dp) :: total_energy
    real(dp) :: kinetic_energy
    real(dp) :: exc_corr_energy
    real(dp) :: pp_local_energy
    real(dp) :: overlap_energy
    real(dp) :: hartree_energy
    real(dp) :: self_energy, electrostatic_energy
    real(dp) :: local_core_energy
    integer :: j, i = 1
    real(dp) :: fact,fact2


    total_energy = 0;
    kinetic_energy  = 0;

    fact = 2*pi/pseudopot%params%box_length;
    fact2 =  pseudopot%params%omega * pseudopot%params%box_length**2/(2*pi)
    
    ! compute kinetic energy
    ! this is give by the formula
    ! E_{kin} = \sum_{j} f_j \sum_{K} | c^{j}(K)| K^2 /2
    do i = 1, numGVects
        do j = 1, numGVects  
            kinetic_energy = kinetic_energy +  FillingFactor(i)*(psi_coeffs_g(j,i)*conjg(psi_coeffs_g(j,i)) &
                             *sum(g_indexes(:,j)**2))*0.5_dp
        enddo
    enddo

    kinetic_energy = kinetic_energy * fact**2

    exc_corr_energy  = 0
    call compute_exc()
    ! compute the exchange correlation energy 
    ! in reciprocal space 
    ! given by $E_{xc} = \Omega \sum_{K} \epsilon_{xc}(K)n(K)$
    exc_corr_energy = sum(exc_g*conjg(density_g))*omega
    print '(A23 F15.8 A23)', "kinetic energy:", kinetic_energy, "ok!"

    total_energy = total_energy + kinetic_energy

    print '(A23 F15.8 A23)', "exc energy:", exc_corr_energy,"ok!"

    total_energy  = total_energy + exc_corr_energy

    ! compute the local pseudopotential energy    
    pp_local_energy = omega*sum(pseudopot%short(g_grid_norm*fact)*structure_factor*conjg(density_g))

    print '(A23 F15.8 A23)', "pp_local_short_energy:", pp_local_energy,  "ok!"

    total_energy = total_energy + pp_local_energy

    ! compute the n total
    total_density_g = core_density_g + density_g

    hartree_energy = fact2 * sum(density_g*conjg(density_g)*(g_grid_norm_inv)**2)
    print '(A23 F15.8)', "hartree_energy", hartree_energy

    
    print *, "Electrostatic energy terms"
    print  '(A23 F15.8)', "Local core energy", fact2 * sum((core_density_g)*conjg(core_density_g)*(g_grid_norm_inv)**2)


    electrostatic_energy = fact2 * sum((total_density_g)*conjg(total_density_g)*(g_grid_norm_inv)**2)
    print '(A23 F15.4)', "total_density_energy", electrostatic_energy



    ! compute the self energy
    self_energy = pseudopot%params%zeff**2/(2*sqrt(pi) * pseudopot%params%chi )
    print '(A23 F15.4 )', "self_energy", self_energy

    electrostatic_energy =  electrostatic_energy - self_energy
    !compute the overlap energy
    overlap_energy = 0

    electrostatic_energy = electrostatic_energy + overlap_energy

    print '(A23 F15.4 A23)', "electrostatic_energy", electrostatic_energy, "ok!" 
    total_energy =  total_energy + electrostatic_energy

    !===============================================================================================
    ! excludes g = 0 automatically
    local_core_energy = fact2 * sum( core_density_g*conjg(core_density_g) * (g_grid_norm_inv)**2) 
    ! print '(A23 F15.8 )', "local_core_energy", local_core_energy   
    !================================================================================================

   ! total_energy = kinetic_energy + exc_corr_energy + pp_local_energy + electrostatic_energy
    print  '(A23 F15.4)', "total_energy", total_energy
end subroutine compute_total_energy
end program autoconsistente