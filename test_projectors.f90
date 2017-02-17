program test_projectors
use projectors
use mesh    
use types
implicit none

real(dp) :: k(1000)
integer :: i
    call init_projectors(1.0_dp, 1.0_dp)
    k = linspace(0.0_dp, 10.0_dp, 1000)

    do i = 1, 1000
        print *, k(i), p1s(k(i))
    enddo

end program test_projectors