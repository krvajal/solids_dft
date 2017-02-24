module gvect
    use types
    use constants    
implicit none
private 
public g_indexes, numGVects, ggen, ni,nj,nk, Nx,Ny,Nz

integer :: Nx, Ny, Nz ! number of grid points in x, y and z directions
integer,allocatable :: g_indexes(:,:) ! indexes of g vectors in a linear array
integer :: numGVects  = 0 
integer :: ni,nj,nk
    

contains 
! ---------------------------------
subroutine ggen(a, eneryCutoff)
    !----------------------------------
    ! Doc here
    !----------------------------------
    real(dp), intent(in) :: a, eneryCutoff

   
    integer :: i,j,k
    real(dp) :: twoPiOverA2 
    real(dp) :: energy
    integer :: maxNumGVects
    
    numGVects = 0

    twoPiOverA2 = (2 * pi / a)**2
    ni = sqrt(2 * eneryCutoff/twoPiOverA2)
    nj = sqrt(2 * eneryCutoff/twoPiOverA2)
    nk = sqrt(2 * eneryCutoff/twoPiOverA2)
    Nx = 4*ni + 3
    Ny = 4*nj + 3
    Nz = 4*nk + 3

    maxNumGVects = 2 * ni + 1
    maxNumGVects = maxNumGVects * (2*nj + 1)
    maxNumGVects = maxNumGVects * (2*nk + 1)

    allocate(g_indexes(3, maxNumGVects))
    do i = -ni,ni
        do j = -nj,nj
            do k =  -nk, nk
                energy = (i**2 + j**2 + k**2)*twoPiOverA2/2
                if(energy < eneryCutoff) then
                    numGVects = numGVects + 1 
                    g_indexes(1,numGVects) = i   
                    g_indexes(2,numGVects) = j
                    g_indexes(3,numGVects) = k
                endif
            enddo
        enddo
    enddo
    

end subroutine ggen


end module gvect