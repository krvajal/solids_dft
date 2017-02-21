module xc

use types
use density
use gvect
use constants
implicit none
private  
public exc, compute_exc,init_xc, vxc_r, computeVxc
real(dp),allocatable :: exc(:,:,:)
real(dp),allocatable :: Num(:,:,:),Den(:,:,:)
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
        allocate(num(Nx,Ny,Nz))
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
    Num = 0
    Den = 0
    do i = 1,4
        Num = Num + a(i) * rs**(i-1)
        Den =  Den + b(i) * rs**i
    enddo
    exc = - Num/Den
end subroutine compute_exc

subroutine computeVxc()
    integer :: i,j   

    rs = (3/(4*pi*NtotalR))**(1./3) 

    Num = 0
    Den = 0
    do i = 1,4
        Num = Num + a(i) * (i-1)* rs**(i-2)
        Den =  Den + b(i) * rs**i
    enddo
    vxc_r = -Num/Den

    Num = 0
    
    do i = 1,4
        Num = Num + a(i) * rs**(i-1)
    enddo
    Num = Num * (b(1) + b(2)*2*rs + b(3)*3 * rs**2 + b(4)*4 * rs**3)
    vxc_r = vxc_r +  Num/(Den*Den)
    vxc_r = vxc_r * NtotalR + exc
    ! once it is compute in the real space, compute in the reciprocal space


end subroutine computeVxc

end module xc