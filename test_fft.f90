program test_fft
    use fft
    use types
    implicit none
    complex(C_DOUBLE_COMPLEX) :: in(5000), out(5000)
    integer :: i
    real(dp) :: x, y

    do i = 1,5000
        x =  i * 0.01
        y = cos(2.0_dp *x) + sin(1.5 * x)
        in(i) = CMPLX(y )
        
    enddo

    do i = 1,5000
        print *, ABS(out(i))
    enddo
end program test_fft
