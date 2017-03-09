module xc

use types
use density
use gvect
use constants
use fft
implicit none
private  
<<<<<<< HEAD
public exc, exc_g, compute_exc,init_xc, vxc_r, compute_vxc

complex(dp),allocatable :: exc_g(:,:,:),exc(:,:,:)

real(dp),allocatable :: num(:,:,:), dnum(:,:,:),den(:,:,:),dden(:,:,:)
=======
public exc, compute_exc,init_xc, vxc_r, compute_vxc
real(dp),allocatable :: exc(:,:,:)

real(dp),allocatable :: num1(:,:,:), num2(:,:,:),den(:,:,:)
>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14

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

<<<<<<< HEAD
    allocate(exc(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(exc_g(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(num(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(dnum(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(den(0:Nx-1,0:Ny-1,0:Nz-1))

    allocate(dden(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(rs(0:Nx-1,0:Ny-1,0:Nz-1))
    allocate(vxc_r(0:Nx-1,0:Ny-1,0:Nz-1))
=======
    allocate(exc(Nx,Ny,Nz))

    allocate(num1(Nx,Ny,Nz))
    allocate(num2(Nx,Ny,Nz))

    allocate(den(Nx,Ny,Nz))
    allocate(rs(Nx,Ny,Nz))
    allocate(vxc_r(Nx,Ny,Nz))
>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14

end subroutine init_xc

! compute the exc using the GTH 
! parametrization of of the pseudopotential
! of Perdew and Wang
subroutine compute_exc()
    integer :: i,j   
    
<<<<<<< HEAD
    rs = (3.0/(4.0*pi*density_r))**(1./3) 

    num = 0
    den = 0
    do i = 1,4
        num = num + a(i) * rs**(i-1)
        den =  den + b(i) * rs**i
    enddo
    exc = - num/den
    call fft_backward_3d(Nx,Ny,Nz,exc, exc_g)
    exc_g = exc_g/(Nx*Ny*Nz)
=======
    rs = (3/(4*pi*NtotalR))**(1./3) 

    num1 = 0
    den = 0
    do i = 1,4
        num1 = num1 + a(i) * rs**(i-1)
        den =  den + b(i) * rs**i
    enddo
    exc = - num1/den
>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14

end subroutine compute_exc

! compute the xc potential of pw parametrization
<<<<<<< HEAD
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
=======
! vxc = exc + n * d exc /dn

subroutine compute_vxc()

    integer :: i 

    call compute_exc()
    
    rs = (3/(4*pi*NtotalR))**(1./3) 
    
    num1 = 0
    num2 = 0
    Den = 0
    do i = 1,4
        den =  Den + b(i) * rs**i
        num1 = num1 + a(i) * (i-1)*rs**(i-2)
        num2 = num2 + b(i) * i * rs**(i-1)
    enddo

    vxc_r = - (num1 - num2 * exc)
    
    vxc_r = exc - 1/3.0_dp * rs * vxc_r

>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14
end subroutine compute_vxc

end module xc