program test_strain_invariants
    use kind_precision_module, only: dp, i32

    use mod_strain_invariants, only: calc_volumetric_strain_invariant, calc_dev_strain_invariant

    implicit none

    ! Tolerance for the value check
    real(kind = dp), parameter :: tolerance = 1.0e-6
    
    ! Variables for checking the strain invariants
    real(kind = dp) :: strain(6)

    ! Strain invariants
    real(kind = dp) :: eps_q, eps_p
    real(kind = dp), parameter :: eps_q_correct = 5.196152422706632_dp, &
                                  eps_p_correct = 6.0_dp

    ! Init a strain vector to check the functions                                  
    strain = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]    

    ! Calc the equivalent strain invariant
    eps_q = calc_dev_strain_invariant(strain)

    ! Calc the volumetric strain invariant
    eps_p = calc_volumetric_strain_invariant(strain)

    call test_assert("Checking Deviatoric Strain Invar.: ", eps_q, eps_q_correct, tolerance)
    call test_assert("Checking Volumetric Strain Invar.: ", eps_p, eps_p_correct, tolerance)

contains
    subroutine test_assert(test_name, actual, expected, tol)
        character(len=*), intent(in) :: test_name
        real(kind = dp), intent(in) :: actual, expected, tol

        if (abs(actual - expected) > tol) then
            print *, "Test '", trim(test_name), "' Failed! Expected:", expected, ", Actual:", actual
        else
            print *, "Test '", trim(test_name), "' Passed."
        endif
    end subroutine test_assert
end program test_strain_invariants