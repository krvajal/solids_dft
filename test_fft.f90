program test_fft
    use fft
    use types

    use constants
    implicit none
    complex(C_DOUBLE_COMPLEX) :: in(2,3,2), out(2,3,2)
    complex(dp) :: result
    integer :: i
    real(dp) :: x, y
    
    in = cmplx(0.0_dp)
    in(1,1,2) = cmplx(1.0_dp)
    in(1,2,1) = cmplx(-1.0_dp)
    in(1,3,1) = cmplx(1.0_dp)
    in(2,1,1) = cmplx(1.0_dp)

    result =  2.0_dp - exp( (0,1)* 2*pi/3.0) + exp((0,1)*4*pi/3);
    call fft_backward_3d(2,3,2,in,out)
    print *, out(1,2,1), result ! ok

end program test_fft
