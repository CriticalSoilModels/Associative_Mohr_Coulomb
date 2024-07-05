module mod_eigen_stress_invariants
    use kind_precision_module, only: dp
    use integer_precision_module, only: i32
    use stdlib_sorting, only: sort

    implicit none

contains

    function voight_2_matrix(voight_vector) result(matrix)
        real(kind = dp), intent (in) :: voight_vector(6)
        real(kind = dp) :: matrix(3,3)

        ! Local variables
        integer(kind = i32) :: i

        ! The order of the voight vector is assumed to be:
        ![ \sigma_{xx}
        !  \sigma_{yy}
        !  \sigma_{zz}
        !  \sigma_{xy}
        !  \sigma_{yz}
        !  \sigma_{zx}
        ! ]

        ! Fill the diagonal
        do i = 1, 3
            matrix(i,i) = voight_vector(i)
        end do
        
        matrix(2, 3) = voight_vector(5)
        matrix(3, 2) = voight_vector(5)
        matrix(1, 3) = voight_vector(6)
        matrix(3, 1) = voight_vector(6)
        matrix(2, 1) = voight_vector(4)
        matrix(1, 2) = voight_vector(4)

    end function

    function calc_eig_mean_stress(eig_vals) result (mean_stress)
        real(kind = dp), intent(in) :: eig_vals(3)
        real(kind = dp)             :: mean_stress
        ! Eqn:
        !   p = ( \sigma_{1} + \sigma_{2} + \sigma_{3} ) / 3.0

        ! Calc the mean stress
        mean_stress = sum(eig_vals) / 3.0

    end function calc_eig_mean_stress

    function calc_eig_J2(eig_vals) result (J2)
        real(kind = dp), intent(in) :: eig_vals(3)
        real(kind = dp) :: J2

        ! Local variables
        real(kind = dp) :: inside
        
        ! The order of the eigen values doens't matter in this case
        ! Since all permutations are down and the result is squared

        inside = ( eig_vals(1) - eig_vals(2) )**2 &
               + ( eig_vals(2) - eig_vals(3) )**2 &
               + ( eig_vals(3) - eig_vals(1) )**2

        J2 = inside / 6.0

    end function calc_eig_J2

    function calc_eig_J(eig_vals) result(J)
        real(kind = dp), intent(in) :: eig_vals(3)
        real(kind = dp) :: J

        ! Local variables
        real(kind = dp) :: J2

        J2 = calc_eig_J2(eig_vals)

        J = sqrt(J2)
    end function calc_eig_J

    function calc_eig_q(eig_vals) result(q)
        real(kind = dp), intent(in) :: eig_vals(3)
        real(kind = dp) :: q

        ! Local variables
        real(kind = dp) :: J2

        J2 = calc_eig_J2(eig_vals)

        q = sqrt( 3.0 * J2 )

    end function calc_eig_q

    ! Calc the determinant of the deviatoric stress matrix using the deviatoric stress eigen values
    function calc_eig_J3(dev_eig_vals) result (det_s)
        real(kind = dp), intent(in) :: dev_eig_vals(3)
        real(kind = dp) :: det_s
        ! J3 = det(s) = 1/3 ( s_{1}^{3} + s_{2}^{3} + s_{3}^{3} )

        ! Local varaibles
        real(kind = dp) :: inside

        ! Cube each eigen value than sum them up
        inside = sum( dev_eig_vals**3 )

        det_s = inside / 3.0
    end function calc_eig_J3

    ! Calc the lode angle
    function calc_eig_lode_angle(eig_vals) result(lode_angle)
        real(kind = dp) :: eig_vals(3)
        real(kind = dp) :: lode_angle

        ! Local variables
        real(kind = dp) :: numer, denom, inside
        

        ! Need to sort the eigen values
        call sort(eig_vals, reverse = .True.)

        ! Calc the numerator
        numer = eig_vals(2) - eig_vals(3)
        
        ! Calc the denominator
        denom = eig_vals(1) - eig_vals(3)

        ! Calc the inside of the fraction
        inside = 2.0_dp * numer/denom - 1.0_dp

        ! Calc the lode angle
        lode_angle = atan( inside / sqrt( 3.0_dp ) )

    end function

end module mod_eigen_stress_invariants

