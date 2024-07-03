module mod_mohr_coulomb_stress_params
    use general_stress_equations, only : calc_stress_invariants
    implicit none
    
    type :: MohrCoulombStressParams
        real :: p
        real :: q
        real :: theta
    contains
        procedure, pass(self) :: update_stress_invariants
    end type
contains

    subroutine update_stress_invariants(self, sigma)
        class(MohrCoulombStressParams), intent(inout) :: self
        real, intent(in) :: sigma(6)
        
        ! Evaluate and store the stress invariants
        call calc_stress_invariants(sigma, self%p, self%q, self%theta)
    
        
    end subroutine update_stress_invariants
    
end module mod_mohr_coulomb_stress_params