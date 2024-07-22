

program test_calc_stress_invariants
    use kind_precision_module, only: dp, i32

    use mod_stress_invariants, only: calc_mean_stress, calc_J2_invariant, &
                                        calc_J_invariant, calc_q_invariant, &
                                        calc_s_determinant, calc_lode_angle_s, &
                                        calc_potts_dJ3_dsigma

    use mod_general_voigt, only: voigt_2_matrix

    use stdlib_linalg, only: eigh, det

    use mod_eigen_stress_invariants, only: calc_eig_mean_stress, calc_eig_J2, calc_eig_J, &
                                          calc_eig_q, calc_eig_J3, calc_eig_lode_angle_bar_s

    implicit none

    ! Constants for tolerance in floating-point comparisons
    real(kind = dp), parameter :: tolerance = 1.0e-6, &
                                  PI = atan(1.0_dp) * 4.0_dp

    ! Variables for input and output
    real(kind = dp) :: Sig(6), stress_matrix(3, 3), dev_matrix(3,3), expected_Sig(6) ! Stress tensor
    real(kind = dp) :: stress_eig_vals(3), dev_eig_vals(3)
    
    real(kind = dp) :: p    , J2    , J    , q    , J3    , lode_angle   ! Invariants
    real(kind = dp) :: p_eig, J2_eig, J_eig, q_eig, J3_eig, lode_angle_eig
    real(kind = dp) :: test_det_dev
    integer(kind = i32) :: i

    ! Derivatives
    real(kind = dp) :: dJ3_dsigma(6)


    ! Test case 1: Pure hydrostatic stress
    Sig = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    
    p = calc_mean_stress(Sig)
    call test_assert("Test 1: Mean Stress Calc: ", p, 1.0_dp, tolerance)
    
    ! Generate a random stress vector
    call random_number(Sig)

    ! Form the matrix
    stress_matrix = voigt_2_matrix( Sig )

    ! Find the eigen values
    call eigh( stress_matrix, stress_eig_vals )

    ! Calc the mean stress
    p = calc_mean_stress(Sig)
    p_eig = calc_eig_mean_stress( stress_eig_vals )

    ! Calc the J2 invariant
    J2 = calc_J2_invariant(Sig)
    J2_eig = calc_eig_J2(stress_eig_vals)

    ! Calc the J invariant
    J = calc_J_invariant(Sig)
    J_eig = calc_eig_J(stress_eig_vals)

    ! Calc the equivalent stress
    q = calc_q_invariant(Sig)
    q_eig = calc_eig_q(stress_eig_vals)

    ! Calc the J3 invariant ( det(s) )
    
    ! Store the values of the stress matrix so "s" can be formed
    dev_matrix = stress_matrix


    ! Calc s
    do i = 1, 3
        dev_matrix(i, i) = stress_matrix(i, i) - p 
    end do 
    
    ! Get the eigen values
    call eigh(dev_matrix, dev_eig_vals)

    ! Calc J3
    J3 = calc_s_determinant( Sig )
    J3_eig = calc_eig_J3( dev_eig_vals )

    ! Calc the lode angle
    lode_angle = calc_lode_angle_s( Sig )
    
    lode_angle_eig = calc_eig_lode_angle_bar_s( stress_eig_vals )

    print *, ""
    print *, "Checking stress invariants"
    call test_assert("Test 4: Mean Stress Calc: ," , p , p_eig , tolerance)
    call test_assert("Test 4: J2 Invariant Calc: ,", J2, J2_eig, tolerance)
    call test_assert("Test 4: J  Invariant Calc: ,", J , J_eig , tolerance)
    call test_assert("Test 4: Dev. Stress Calc: ," , q , q_eig , tolerance)
    call test_assert("Test 4: J3 Invariant Calc: ,", J3, J3_eig, tolerance)
    call test_assert("Test 4: Lode Angle Calc  : ,", lode_angle, lode_angle_eig, tolerance)


    print *, "Done checking stress invariants"

    ! Test the derivatives
    Sig = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp]
    expected_Sig = [100.0_dp, -1100.0_dp, 1000.0_dp, 5200.0_dp, 5800.0_dp, 4000.0_dp]

    dJ3_dsigma = calc_potts_dJ3_dsigma(Sig)

    call test_assert_vector("Test 5: dJ3/dsigma: ," , expected_Sig , dJ3_dsigma , tolerance)

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

    subroutine test_assert_vector(test_name, actual, expected, tol)
        character(len=*), intent(in) :: test_name
        real(kind = dp), intent(in) :: actual(:), expected(:), tol

        if (norm2(actual - expected) > tol) then
            print *, "Test '", trim(test_name), "Failed!"
            print *, "Expected: ", expected
            print *, "Actual  : ", actual
        else
            print *, "Test '", trim(test_name), "' Passed."
        endif
    end subroutine test_assert_vector

end program test_calc_stress_invariants