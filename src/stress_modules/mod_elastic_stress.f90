module mod_elastic_stress
    use kind_precision_module, only: dp
    use integer_precision_module, only: i32

    implicit none

contains
    !> \brief Constructs the stiffness matrix DE based on shear modulus and Poisson's ratio.
    !!
    !! This function calculates and constructs a 6x6 stiffness matrix based on the given
    !! shear modulus and Poisson's ratio.
    !!
    !! \param[in] shear_modulus Shear modulus of the material.
    !! \param[in] poisson_ratio Poisson's ratio of the material.
    !! \return Stiffness matrix DE of size (6, 6).
    !!
    function construct_stiffness_matrix( shear_modulus, poisson_ratio ) result(stiff_matrix)
        ! Construct the stiffness matrix DE
        real(kind = dp), intent(in) :: shear_modulus, poisson_ratio
        real(kind = dp) :: stiff_matrix(6, 6)

        ! Local variables
        real(kind = dp) :: F1, F2

        stiff_matrix(:, :) = 0.0_dp

        ! Calc the values needed for the stiffness matrix
        F1 = 2 * shear_modulus * ( 1 - poisson_ratio ) / ( 1- 2 * poisson_ratio )
        F2 = 2 * shear_modulus * poisson_ratio  ( 1 - 2 * poisson_ratio )

        !---- Fill the stiffness matrix ----
        
        ! Zero the matrix
        stiff_matrix = 0.0_dp

        ! Fill the upper block
        stiff_matrix(1:3, 1:3) = F2
        
        ! Loop over the first three diagonals
        do i = 1, 3
            stiff_matrix(i, i) = F1 
        end do
        ! Loop over the 4th through 6th diagonal
        do i = 4, 6
            stiff_matrix(i, i) = shear_modulus
        end do
    end function construct_stiffness_matrix

    !> \brief Calculates the increment of stress.
    !!
    !! This function calculates the increment of stress \f$ d\sigma \f$ given the stiffness 
    !! matrix \f$ DE \f$ and the strain \f$ \epsilon \f$.
    !!
    !! \param[in] strain Strain vector \f$ \epsilon \f$ of size (6).
    !! \param[in] stiff_matrix Stiffness matrix \f$ DE \f$ of size (6, 6).
    !! \return Stress increment vector \f$ d\sigma \f$ of size (6).
    !!
    function calc_elastic_stress_increment(strain, stiff_matrix ) result(stress_inc)
        real(kind = dp), intent(in) :: strain(6), stiff_matrix(6,6)
        real(kind = dp) :: stress_inc(6)

        ! Calc the stress increment
        stress_inc = matmul(stiff_matrix, strain)

    end function

    !> \brief Calculates the new stress based on the given strain, stiffness matrix, and initial stress.
    !!
    !! This function calculates the new stress state by computing the stress increment
    !! using the stiffness matrix and strain, and then adding this increment to the initial stress.
    !!
    !! \param[in] strain Strain vector \f$ \epsilon \f$ of size (6).
    !! \param[in] stiff_matrix Stiffness matrix \f$ DE \f$ of size (6, 6).
    !! \param[in] init_stress Initial stress vector \f$ \sigma_0 \f$ of size (6).
    !! \return New stress vector \f$ \sigma \f$ of size (6).
    !!
    function calc_elastic_stress(strain, stiff_matrix, init_stress) result(new_stress)
        real(kind = dp), intent(in) :: strain(6), stiff_matrix(6,6), init_stress(6)
        real(kind = dp) :: new_stress

        ! Local variables
        real(kind =dp) :: stress_inc(6) ! Stress increment

        ! Calc the stress increment
        stress_inc = calc_elastic_stress_increment(strain, stiff_matrix)

        ! Calc the new stress
        new_stress = init_stress + stress_inc
    end function

end module mod_elastic_stress