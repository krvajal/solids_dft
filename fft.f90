! some routines to compute the 
! discrete fourier transform 
! using the fftw library

module fft

use, intrinsic :: iso_c_binding ! allow c binding
include 'fftw3.f03'

contains

subroutine fftForward3D(Nx,Ny,Nz,in,out)
    integer,intent(in) :: Nx,Ny,Nz 
    complex(C_DOUBLE_COMPLEX),intent(inout) :: in(:,:,:)
    complex(C_DOUBLE_COMPLEX),intent(out) :: out(:,:,:)
    type(C_PTR) :: plan
     
    plan = fftw_plan_dft_3d(Nz,Ny,Nx,in,out,FFTW_FORWARD, FFTW_ESTIMATE)
    call fftw_execute_dft(plan,in,out)
    call fftw_destroy_plan(plan)

end subroutine fftForward3D

subroutine fftBackward3D(Nx,Ny,Nz,in,out)
    integer,intent(in) :: Nx,Ny,Nz
    complex(C_DOUBLE_COMPLEX),intent(inout) :: in(:,:,:)
    complex(C_DOUBLE_COMPLEX),intent(out) :: out(:,:,:)
    type(C_PTR) :: plan
     
    plan = fftw_plan_dft_3d(Nz,Ny,Nx,in,out,FFTW_BACKWARD, FFTW_ESTIMATE)
    call fftw_execute_dft(plan,in,out)
    call fftw_destroy_plan(plan)

end subroutine fftBackward3D

end module fft