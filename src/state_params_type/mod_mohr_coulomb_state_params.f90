module mod_mohr_coulomb_state_params
    use kind_precision_module, only: dp ! Import double precision floats

    implicit none
    
    type :: MohrCoulombStateParameters
        real(kind = dp) :: cohesion
        real(kind = dp) :: fric_angle

    end type

end module mod_mohr_coulomb_state_params