module mod_mohr_coulomb_stress_params
    use kind_precision_module, only : dp
    use mod_stress_invariants, only : calc_stress_invariants
    implicit none
    
    type :: MohrCoulombStressParams
        real(kind = dp) :: p
        real(kind = dp) :: q
        real(kind = dp) :: theta
    contains
        procedure, pass(self) :: update_stress_invariants
    end type
contains

    subroutine update_stress_invariants(self, sigma)
        class(MohrCoulombStressParams), intent(inout) :: self
        real(kind = dp), intent(in) :: sigma(6)
        
        ! Evaluate and store the stress invariants
        call calc_stress_invariants(sigma, self%p, self%q, self%theta)
    
        
    end subroutine update_stress_invariants
    
end module mod_mohr_coulomb_stress_params