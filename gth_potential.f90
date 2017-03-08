module gth_potential
    use types
    use constants
    implicit none

    private 
    public pseudo_pot_gth, GthPotParams, gth_pp_t

    type :: GthPotParams
        real(dp) :: chi
        real(dp) :: c1
        real(dp) :: c2
        real(dp) :: omega
        integer(dp) :: Zeff
        real(dp) :: rloc
        real(dp) :: box_length
    end type GthPotParams

    type :: gth_pp_t
        type(GthPotParams) :: params
        contains 
        procedure :: local => pseudo_pot_local
        procedure :: core => pseudo_pot_core
        procedure :: total => pseudo_pot_gth

    end type gth_pp_t

contains

elemental real(dp) function pseudo_pot_gth(this, k) result(retval)
    class(gth_pp_t),intent(in) :: this 
    real(dp),intent(in) :: k
       
    if (k /= 0) then
        retval =  this%core(k)
        retval = retval + this%local(k)
    else 
        retval = 0
    endif
end function pseudo_pot_gth
 
elemental real(dp) function pseudo_pot_core(this, k) result(retval)
    class(gth_pp_t),intent(in) :: this
    real(dp),intent(in) :: k
    real(dp) :: chi

    chi = this%params%chi /this%params%box_length    
    
    retval = 0.0

    if(abs(k) > tiny(1.0_dp)) then
        retval = - 4 *pi * this%params%Zeff / this%params%omega
        retval = retval * exp(- 0.5_dp * (k * chi)**2)
        retval = retval /(k**2)
    endif 

end function pseudo_pot_core


elemental real(dp) function pseudo_pot_local(this, k) result(retval)

    class(gth_pp_t),intent(in) :: this
    real(dp),intent(in) :: k

    type(GthPotParams) :: params
    real(dp) :: chi

    params = this%params
    chi =  params%chi / params%box_length
    
    retval = (2*pi)**(1.5_dp)* chi**3 /params%omega
    retval = retval * exp(-0.5_dp * (k * chi)**2)
    retval = retval * (params%c1 + params%c2 * (3 - (k*chi)**2))

end function pseudo_pot_local

end module gth_potential