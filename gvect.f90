module gvect
    use types
    use constants
implicit none

private 
public gVects, numGVects, ggen, ni,nj,nk

real(dp),allocatable :: gVects(:,:)
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

    maxNumGVects = 2*ni + 1;
    maxNumGVects = maxNumGVects * (2*nj + 1)
    maxNumGVects = maxNumGVects * (2*nk + 1)
    print *, maxNumGVects
    allocate(gVects(3, maxNumGVects))
    do i = -ni,ni
        do j = -nj,nj
            do k =  -nk, nk
                energy = (i**2 + j**2 + k**2)*twoPiOverA2/2
                if(energy < eneryCutoff) then
                    numGVects = numGVects + 1 
                    gVects(1,numGVects) = i   
                    gVects(2,numGVects) = j
                    gVects(3,numGVects) = k
                endif
            enddo
        enddo
    enddo
    print *, numGVects
end subroutine ggen

function layoutKIndexForFft() result(retval)
    real(dp),allocatable :: retval(:)
    integer :: size , i
    integer :: toI,toJ, toK
    size = (2*ni + 1)* (3*nj + 1)*(2*nk + 1)
    allocate(retval(size))
    retval = 0.0 !init to zero
    do i = 1, numGVects
         toI = gVects(1,i)
         toJ = gVects(2,i)
         toK = gVects(3,i)
            
    enddo


end function layoutKIndexForFft
end module gvect