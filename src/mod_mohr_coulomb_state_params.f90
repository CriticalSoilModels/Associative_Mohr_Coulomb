module mod_mohr_coulomb_state_params
    use kind_precision_module, only: dp ! Import double precision floats
    use degrees_utils, only: convert_deg_2_radians
    implicit none
    
    type :: MohrCoulombStateParameters

        real(kind = dp) :: c    ! Cohesion
        real(kind = dp) :: phi  ! Friction angle
        real(kind = dp) :: G    ! Shear modulus
        real(kind = dp) :: enu   ! Poissons ratio
    
    contains
        procedure, pass(self) :: initialize

    end type

contains

    subroutine initialize(self, c, phi, G, nu)
        class(MohrCoulombStateParameters), intent(inout) :: self
        real(kind = dp), intent(in) :: c, phi, G, nu

        self%phi = convert_deg_2_radians(phi) ! Set the friction angle for the model and convert it to radians
        print *, self%phi
        
        self%c  = c  ! Set the cohesion for the model
        self%G  = G  ! Set the shear modulus
        self%enu = nu ! Set the poisson's ratio
    end subroutine initialize
    
end module mod_mohr_coulomb_state_params