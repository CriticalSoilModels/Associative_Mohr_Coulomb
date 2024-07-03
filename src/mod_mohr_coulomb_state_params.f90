module mod_mohr_coulomb_state_params
    implicit none
    
    type :: MohrCoulombStateParameters
        real :: cohesion
        real :: fric_angle

    end type

end module mod_mohr_coulomb_state_params