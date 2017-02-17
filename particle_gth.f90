! solves the problem of of a particle with the gth potential
program particle_gth
    use types
    use constants
    use linalg
    use gth_potential
    integer,parameter :: basisDim =  7
    real(dp) :: KSet(basisDim, 3)
    integer :: i, j, k
    real(dp) :: length = 5.0_dp
    real(dp) :: kLength
    real(dp):: hamiltMatrix(basisDim, basisDim)
    real(dp):: energies (basisDim)
    kLength  = 2 * pi/length
    
    KSet(1,:) = (/ 0,0,0/)
    KSet(2,:) = (/ kLength,0.0_dp,0.0_dp/)
    KSet(3,:) = (/ 0.0_dp, kLength,0.0_dp/)
    KSet(4,:) = (/ 0.0_dp, 0.0_dp, kLength/)
    KSet(5,:) = (/ -kLength, 0.0_dp, 0.0_dp/)
    KSet(6,:) = (/ 0.0_dp, -kLength, 0.0_dp/)
    KSet(7,:) = (/ 0.0_dp, 0.0_dp, -kLength/)


    call fillHamiltMatrix(basisDim, KSet, hamiltMatrix)
    energies = eigvals(hamiltMatrix)
    print *,energies
contains
subroutine fillHamiltMatrix(basisDim, KBasisSet, hamiltMatrix)
    integer,intent(in) :: basisDim
    real(dp),intent(in) :: KBasisSet(:,:)
    real(dp), intent(out) :: hamiltMatrix(basisDim, basisDim)
    type(GthPotParams) :: paramsHidrogen
    paramsHidrogen%c1 = -4.0663326_dp
    paramsHidrogen%c2 = 0.6778322_dp
    paramsHidrogen%chi = 0.2
    paramsHidrogen%Zeff = 1 
    paramsHidrogen%omega = 125.0 ! cell with edge length of 5.0 a.u.

    hamiltMatrix = 0
    do i = 1, basisDim
        hamiltMatrix(i,i) = norm2(KBasisSet(i,:))**2 *  0.5_dp
        do j = 1, basisDim
            hamiltMatrix(i,j) = hamiltMatrix(i,j) + &
                                gthPotential(paramsHidrogen, KBasisSet(i,:), KBasisSet(j,:))
        enddo
    enddo
  
end subroutine fillHamiltMatrix



end program particle_gth
