module gth_potential
    use types
    use constants
    implicit none

    private 
    public pseudo_pot_gth, GthPotParams, pseudo_pot_local,pseudo_pot_core

    
    type :: GthPotParams
        real(dp) :: chi
        real(dp) :: c1
        real(dp) :: c2
        real(dp) :: omega
        integer(dp) :: Zeff
        real(dp) :: rloc
        real(dp) :: box_length
    end type GthPotParams

contains

real(dp) function pseudo_pot_gth(params, k1, k2) result(retval)

    type(GthPotParams) :: params
    real(dp) :: k1(3), k2(3)
    real(dp) :: k
    k = norm2(k1 - k2)
    retval =  pseudo_pot_core(params,k) 
    retval = retval + pseudo_pot_local(params,k)

end function pseudo_pot_gth
 
real(dp) function pseudo_pot_core(params, k) result(retval)
    type(GthPotParams) :: params
    real(dp) :: k, chi
    chi = params%chi /params%box_length    
    retval = 0.0
    if(abs(k) > tiny(1.0_dp)) then
        retval = - 4 *pi * params%Zeff / params%omega
        retval = retval * exp(- 0.5_dp * (k * chi)**2)
        retval = retval /(k**2)
    endif 
end function pseudo_pot_core


real(dp) function pseudo_pot_local(params, k) result(retval)
    type(GthPotParams) :: params
    real(dp) :: k,chi

    chi =  params%chi / params%box_length
    retval = (2*pi)**(1.5_dp)* chi**3 /params%omega
    retval = retval * exp(-0.5_dp * (k * chi)**2)
    retval = retval * (params%c1 + params%c2 * (3 - (k*chi)**2))

end function pseudo_pot_local

end module gth_potential