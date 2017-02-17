! solves the problem of a free particle in a box
program freeparticle
    use types
    use constants
    use linalg
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
    hamiltMatrix = 0
    do i = 1, basisDim
        hamiltMatrix(i,i) = dot_product(KBasisSet(i,:),KBasisSet(i,:))**2/ 0.5_dp
    enddo
end subroutine fillHamiltMatrix

end program freeparticle
