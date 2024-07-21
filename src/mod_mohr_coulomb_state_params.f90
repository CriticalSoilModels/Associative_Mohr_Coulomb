module class_assoc_MC_state_params
    use degrees_utils, only: convert_deg_2_radians
    use kind_precision_module, only: dp ! Import double precision floats
    implicit none
    
    type, public :: assoc_MC_state_params

        real(kind = dp) :: phi  ! Friction angle
        real(kind = dp) :: c    ! Cohesion
        real(kind = dp) :: G    ! Shear modulus
        real(kind = dp) :: enu   ! Poissons ratio
    
    end type

    interface assoc_MC_state_params
        module procedure initialize
    end interface

contains

    function initialize(c, phi, G, enu) result(this)
        real(kind = dp), intent(in) :: c, phi, G, enu
        type(assoc_MC_state_params) :: this

        this%phi = convert_deg_2_radians(phi) ! Set the friction angle for the model and convert it to radians
        
        this%c  = c  ! Set the cohesion for the model
        this%G  = G  ! Set the shear modulus
        this%enu = enu ! Set the poisson's ratio
    end function initialize
    
end module class_assoc_MC_state_params