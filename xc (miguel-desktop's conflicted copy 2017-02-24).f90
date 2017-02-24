module xc

use types
use density
use gvect
use constants
implicit none
private  
public exc, compute_exc,init_xc, vxc_r, compute_vxc
real(dp),allocatable :: exc(:,:,:)

real(dp),allocatable :: num(:,:,:), dnum(:,:,:),den(:,:,:),dden(:,:,:)

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

    allocate(exc(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(num(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(dnum(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(den(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(dden(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(rs(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(vxc_r(0:Nx-1,0:Ny-1,0:Nz-1))

end subroutine init_xc

! compute the exc using the GTH 
! parametrization of of the pseudopotential
! of Perdew and Wang
subroutine compute_exc()
    integer :: i,j   
    
    rs = (3.0/(4.0*pi*density_r))**(1./3) 

    num = 0
    den = 0
    do i = 1,4
        num = num + a(i) * rs**(i-1)
        den =  den + b(i) * rs**i
    enddo
    exc = - num/den
end subroutine compute_exc

! compute the xc potential of pw parametrization
! vxc = exc + n * d exc / dn
subroutine compute_vxc()

    integer :: i 
    real(dp) , allocatable :: dvxc_drs(:,:,:)

    call compute_exc()
    
    rs = (3.0/(4.0_dp*pi*density_r))**(1./3) 
    
    num = 0
    dnum = 0
    den = 0
    do i = 1,4
        den =  den + b(i) * rs**i
        num = num + a(i) * (i-1)*rs**(i-2)
        dnum = dnum + b(i) * i * rs**(i-1)
    enddo
      
    num = 0
    dnum = 0
    den = 0
    dden = 0
    do i = 1, 4
        num = num + a(i)*Rs**(i-1)
        dnum = dnum + a(i)*(i-1)*Rs**(i-2)
        den = den + b(i)*Rs**(i+3)
        dden = dden + b(i)*(i+3)*Rs**(i+2)
    enddo
    vxc_r = (den*dnum-num*dden)/den**2*Rs**4/3.0_dp
  
    ! dvxc_drs = - (num + dnum * exc)/den
    
    ! vxc_r = exc - 1.0_dp/3.0_dp * rs * dvxc_drs
    ! vxc_r =  -dvxc_drs*rs/3
end subroutine compute_vxc

end module xc