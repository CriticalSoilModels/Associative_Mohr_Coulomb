module mod_mohr_coulomb_stress_params
    use kind_precision_module, only : dp
    use mod_stress_invariants, only : calc_potts_MC_stress_invariants
    implicit none
    
    type :: MohrCoulombStressParams
        real(kind = dp) :: p
        real(kind = dp) :: J
        real(kind = dp) :: lode_angle
    contains
        procedure, pass(self) :: update_stress_invariants
    end type
contains

    ! Initialize the intial values of the state parameters using the stress matrix
    subroutine initialize(self, sigma)
        class(MohrCoulombStressParams), intent(inout) :: self
        real(kind = dp), intent(in) :: sigma(6)

        call calc_potts_MC_stress_invariants(sigma, self%p, self%J, self%lode_angle)
    end subroutine

    subroutine update_stress_invariants(self, sigma)
        class(MohrCoulombStressParams), intent(inout) :: self
        real(kind = dp), intent(in) :: sigma(6)
        
        ! Evaluate and store the stress invariants
        call calc_potts_MC_stress_invariants(sigma, self%p, self%J, self%lode_angle)
        
    end subroutine update_stress_invariants
    
end module mod_mohr_coulomb_stress_params