module iondata
    use types
implicit none

integer              :: num_ions = 1
real(dp),allocatable :: ion_pos(:,:)
contains
    subroutine init_ions()
        implicit none
        allocate(ion_pos(num_ions,3))
        ion_pos(1,:) = 0.0_dp
        ion_pos(1,1) = 1 !1.0_dp
        ion_pos(1,2) = 3 !1.0_dp
        ion_pos(1,3) = 0.2 !1.0_dp
    end subroutine init_ions

end module iondata