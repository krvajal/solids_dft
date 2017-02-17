! set of projectors 
! p_i^l(K) to be used 
! in the GTH psudopotential calculations
module projectors
    use types
    use constants
implicit none
private 
public init_projectors, p1s

!internal parameter for the module
real(dp) :: rs_
real(dp) :: omega_
real(dp) :: omegaSqrtInv_
real(dp),parameter :: pi54 = pi **(5.0/4.0) ! pi^(5/4)

contains

subroutine init_projectors(rs, omega)
    real(dp) :: rs, omega
    rs_ = rs
    omega_ = omega
    omegaSqrtInv_ = 1.0/sqrt(omega)
end subroutine init_projectors

real(dp) function p1s(k) result(retval)
    real(dp) :: k
    
    retval = omegaSqrtInv_ * 4 * rs_
    retval = retval * sqrt(2*rs_)*pi54
    retval = retval * exp(-0.5*(k*rs_)**2)
end function p1s


end module projectors