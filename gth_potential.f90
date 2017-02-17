module gth_potential
    use types
    use constants
    implicit none

    private 
    public gthPotential, GthPotParams
    type :: GthPotParams
        real(dp) :: chi
        real(dp) :: c1
        real(dp) :: c2
        real(dp) :: omega
        real(dp) :: Zeff
        real(dp) :: rloc
    end type GthPotParams

contains

real(dp) function gthPotential(params, k1, k2) result(retval)

    type(GthPotParams) :: params
    real(dp) :: k1(3), k2(3)
    real(dp) :: k
    k = norm2(k1 - k2)
    retval = 0
    if( k  >  tiny(1.0_dp)) retval =  gthPotCoreTerm(params,k) 
    retval = retval + gthPotLocalTerm(params,k)

end function gthPotential
 
real(dp) function gthPotCoreTerm(params, k) result(retval)
    type(GthPotParams) :: params
    real(dp) :: k
    retval = -4 *pi * params%Zeff / params%omega
    retval = retval * exp(- 0.5 * (k * params%chi)**2)
    retval = retval /(k*k)

end function gthPotCoreTerm


real(dp) function gthPotLocalTerm(params, k) result(retval)
    type(GthPotParams) :: params
    real(dp) :: k
    retval = (2*pi)**(3./2)* params%chi**3 /params%omega
    retval = retval * exp(-0.5 * (k * params%chi)**2)
    retval = retval * (params%c1 + params%c2 * (3 - (k*params%chi)**2))
end function gthPotLocalTerm



end module gth_potential