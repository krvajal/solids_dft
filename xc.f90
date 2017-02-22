module xc

use types
use density
use gvect
use constants
implicit none
private  
public exc, compute_exc,init_xc, vxc_r, computeVxc
real(dp),allocatable :: exc(:,:,:)
real(dp),allocatable :: num1(:,:,:), num2(:,:,:),den(:,:,:)
real(dp),allocatable :: vxc_r(:,:,:) !   xc potential in real space
real(dp),allocatable :: exc_(:,:,:)
real(dp),allocatable :: rs(:,:,:)
real(dp) ::  a(4) = [0.4581652932831429, &
                     2.217058676663745,  &
                     0.7405551735357053, &
                     0.01968227878617998]
real(dp) :: b(4) = [1.0, &
                    4.504130959426697, &
                    1.110667363742916, &
                    0.02359291751427506 ]


contains 

subroutine init_xc()

    allocate(exc(Nx,Ny,Nz))
    allocate(num1(Nx,Ny,Nz))
    allocate(num2(Nx,Ny,Nz))
    allocate(den(Nx,Ny,Nz))
    allocate(rs(Nx,Ny,Nz))
    allocate(vxc_r(Nx,Ny,Nz))

end subroutine init_xc

! compute the exc using the GTH 
! parametrization of of the pseudopotential
! of Perdew and Wang
subroutine compute_exc()
    integer :: i,j   
    
    rs = (3/(4*pi*NtotalR))**(1./3) 
    num1 = 0
    den = 0
    do i = 1,4
        num1 = num1 + a(i) * rs**(i-1)
        den =  den + b(i) * rs**i
    enddo
    exc = - num1/den
end subroutine compute_exc

! compute the xc potential of pw parametrization
! vxc = exc + n * d exc /dn
subroutine compute_vxc()

    integer :: i 

    call compute_exc()
    
    rs = (3/(4*pi*NtotalR))**(1./3) 
    compute_exc
    num1 = 0
    num2 = 0
    Den = 0
    do i = 1,4
        den =  Den + b(i) * rs**i
        num1 = num1 + a(i) * (i-1)*rs**(i-2)
        num2 = num2 + b(i) * i * rs**(i-1)
    enddo

    vxc_r = - (num1 - num2 * exc)
    
    vxc_r = exc - 1/3.0_dp * rs * v_xcr

end subroutine compute_vxc

end module xc