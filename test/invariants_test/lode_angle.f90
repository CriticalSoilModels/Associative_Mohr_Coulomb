program lode_angle_test

    use kind_precision_module, only: dp, i32
    use mod_stress_invariants, only: calc_lode_angle_s

    implicit none
    
    real(kind = dp), parameter :: tolerance = 1e-6, &
                                  PI = atan(1.0_dp) * 4.0_dp
                                  
    real(kind = dp) :: trx_compr_stress(6), trx_exten_stress(6), shear_stress(6)
    real(kind = dp) :: trx_compr_lode_angle, trx_exten_lode_angle, shear_lode_angle

    ! Testing that the lode angle gives the right value at for:
    ! See table: https://en.wikipedia.org/wiki/Lode_coordinates
    ! See Potts and Zdravković pg. 116
        ! Triaxial Compression: sigma_1 >= sigma_2 = sigma_3    -> \theta = -pi/6
        ! Shear               : sigma_2 = (sigma_1 + sigma_3)/2 -> \theta = 0.0 
        ! Triaxial Extension  : sigma_1 = sigma_2 >= sigma_3    -> \theta = pi/6 

    trx_compr_stress = [10.0_dp, 1.0_dp, 1.0_dp , 0.0_dp, 0.0_dp, 0.0_dp]
    trx_compr_lode_angle = calc_lode_angle_s(trx_compr_stress)
    
    shear_stress = [10.0_dp, 5.0_dp, 15.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    shear_lode_angle     = calc_lode_angle_s(shear_stress)
    
    trx_exten_stress     = [3.0_dp , 3.0_dp, 1.0_dp , 0.0_dp, 0.0_dp, 0.0_dp]
    trx_exten_lode_angle = calc_lode_angle_s(trx_exten_stress)

    call test_assert("Test 4: Lode Angle Trx Compression  : ,", trx_compr_lode_angle, -PI/6.0_dp, tolerance)
    call test_assert("Test 4: Lode Angle Shear  : ,"          , shear_lode_angle    , 0.0_dp    , tolerance)
    call test_assert("Test 4: Lode Angle Trx Extension    : ,", trx_exten_lode_angle, PI/6.0_dp , tolerance)

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
end program lode_angle_test