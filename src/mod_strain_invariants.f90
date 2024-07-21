module mod_strain_invariants
    use kind_precision_module, only : dp, i32
 
    ! use mod_general_voigt, only: multiply_voigt_vectors, get_voigt_identity_vector, trace_voigt_vector
 
    implicit none
 
contains
    
    ! Calc the volumetric strain invariant
    pure function calc_volumetric_strain_invariant(strain) result(eps_p)
        real(kind = dp), intent(in) :: strain(6)
        real(kind = dp) :: eps_p

        ! Calc the valumetric strain invariant
        eps_p = sum( strain(1:3) )

    end function calc_volumetric_strain_invariant

    pure function calc_dev_strain_invariant(strain) result(eps_q)
        real(kind = dp), intent(in) :: strain(6)
        real(kind = dp) :: eps_q

        ! Local variables
        real(kind = dp) :: first_term, second_term, third_term, shear_term
        real(kind = dp), parameter :: TWO = 2.0_dp             ,&
                                      THREE = 3.0_dp           ,&
                                      ONE_THIRD = 1.0_dp/THREE

        ! Calc (eps_{yy} - eps_{zz})^{2}
        first_term = (strain(2) - strain(3))**2

        ! Calc (eps_{zz} - eps_{xx})^{2}
        second_term = (strain(3) - strain(1))**2

        ! Calc (eps_{xx} - eps_{yy})^{2}
        third_term = (strain(1) - strain(2))**2

        ! Calc the shear term, 
        shear_term = sum(strain(4:6)**2)

        ! Calc the deviatoric strain invariant
        eps_q = ONE_THIRD * sqrt(TWO * ( first_term + second_term + third_term ) + THREE * shear_term)
    end function calc_dev_strain_invariant

end module mod_strain_invariants