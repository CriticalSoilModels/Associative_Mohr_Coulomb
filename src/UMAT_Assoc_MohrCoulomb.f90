module mod_UMAT_assoc_MC
    use kind_precision_module   , only: dp, i32
    use class_assoc_MC_state_params , only: assoc_MC_state_params
    use class_assoc_MC_stress_params, only: assoc_MC_stress_params
    use class_assoc_MC_yield_surf   , only: assoc_MC_yield_surf
    use mod_elastic_stress          , only: construct_stiffness_matrix, calc_elastic_stress
    
    implicit none

    private
    public :: UMAT_assoc_MC
    
    
contains
    !> \brief Evaluate the Mohr-Coulomb stress-strain relationship.
    !>
    !> This subroutine evaluates the stress-strain relationship using the Mohr-Coulomb yield criterion. 
    !> It updates the stress tensor and calculates the stiffness matrix given the strain increment, state variables, 
    !> and material properties.
    !>
    !> \param[inout] stress              Stress tensor components in Voigt notation (6 components).
    !> \param[inout] strain_increment          Incremental strain tensor components.
    !> \param[inout] state_variables     State variables (array).
    !> \param[in]    props               Material properties (array).
    !> \param[out]   stiff_matrix        Stiffness matrix (6x6 components).
    subroutine UMAT_assoc_MC(stress, strain_increment, state_vars, props, stiff_matrix, max_iterations)
        use kind_precision_module   , only: dp, i32

        real(kind = dp), intent(inout) :: stress(6), state_vars(:)
        real(kind = dp), intent(in)    :: strain_increment(6), props(:)  
        real(kind = dp), intent(out)   :: stiff_matrix(6, 6)
        integer(kind = i32), intent(in):: max_iterations ! Make this an input parameter later
        
        ! Local variables
        type(assoc_MC_state_params) :: Xi
        type(assoc_MC_stress_params) :: stress_params
        type(assoc_MC_yield_surf) :: yield_surf

        real(kind = dp)     :: yield_surf_tolerance = 1e-8
        integer(kind = i32) :: integration_scheme = 1      ! Make this an input
        integer(kind = i32), parameter :: ortiz_simo = 1                                             
        logical, parameter :: DEBUG = .False.

        ! integer(kind = i32) :: max_iterations = 100 ! Make this an input parameter

        ! Flag for setting the integration scheme this can be an input later on
        integration_scheme = 1

        ! Set all of the possible integration schemes here
        ! Another one...
        ! Another one...

        ! Order is cohesion, phi (friction angle), shear modulus, poisson's ratio 
        Xi = assoc_MC_state_params( state_vars(1), state_vars(2), props(1), props(2))  ! Initialize the state paramters object

        stress_params = assoc_MC_stress_params(stress)                                  ! Initialize the stress parameters object
        
        yield_surf = assoc_MC_yield_surf(yield_surf_tolerance, max_iterations, 0.0_dp, 0.0_dp)       ! Initialize the yield surface object
        
        select case (integration_scheme) ! If yielding select the integration method
            case(ortiz_simo)
                call ortiz_simo_assoc_mohr_coulomb( yield_surf, stress, Xi, stress_params, strain_increment, stiff_matrix)
                
            case default
                print *, "The integration method selected for the Associative Mohr Coulomb surface is not valid"
                print *, "Valid inputs are 1 (ortiz-simo)"

        end select

        ! Make sure no other variables need to be updated

        ! End the subroutine
        ! Init
        if (DEBUG) then
            print *, "yield surface value", yield_surf%val
            print *, "stress", stress
        end if 
    end subroutine UMAT_assoc_MC

    subroutine ortiz_simo_assoc_mohr_coulomb(yield_surf, stress, Xi, stress_params, strain_increment, stiff_matrix)
        use kind_precision_module   , only: dp, i32

        real(kind = dp), intent(inout) :: stress(6)
        real(kind = dp), intent(in)    :: strain_increment(6)
        real(kind = dp), intent(inout)    :: stiff_matrix(6, 6)

        ! Local variables
        class(assoc_MC_state_params) :: Xi
        class(assoc_MC_stress_params) :: stress_params
        class(assoc_MC_yield_surf) :: yield_surf
        integer(kind = i32) :: counter = 0 ! Init counter of iterations to zero
        real(kind = dp) :: dF_dsigma_D_dP_dsigma, D_dp_dsigma(6), dLambda
        logical, parameter :: DEBUG = .FALSE.

        ! Form the stiffness matrix
        stiff_matrix = construct_stiffness_matrix(Xi%G, Xi%enu) 

        ! Do the elastic prediction
        stress = calc_elastic_stress(strain_increment, stiff_matrix, stress)

        if (DEBUG) print *, "elastic update: ", stress

        ! Update the stress params
        call stress_params%update_stress_invariants(stress)

        ! Update the state params
            ! Not needed for assoc MC
        
        ! Evaluate the yield surface
        call yield_surf%evaluate_surface(Xi, stress_params)
        
        ! Update the stiffness matrix
            ! Don't need to do it in this case

        if (DEBUG) then
            print *, "Init yield surf value: ", yield_surf%val
            print *, "yielding condition", yield_surf%is_yielding()
        end if

        ! If the surface isn't yielding this statement won't activate
        do while (counter <= yield_surf%max_iterations .and. yield_surf%is_yielding())
        ! If yielding is occuring
            
            ! Evaluate the n vector (in this case n = m)
            call yield_surf%evaluate_dF_dsigma(Xi, stress_params, stress)
            
            ! Evaluate the stiffness matrix
                ! This is a constant in this case

            ! Evaluate the derivative of the plastic potential wrt. stress
                ! This is the same as dF/dsigma in this case
            
            ! Evaluate Stiffness matrix : dPlastic potential_dSigma
             
            D_dp_dsigma = matmul( stiff_matrix, yield_surf%dF_dsigma )
            
            ! n:D:m, tensor product of dF_dsigma : Stiffness matrix : dPlastic potential_dSigma
            dF_dsigma_D_dP_dsigma = dot_product( yield_surf%dF_dsigma, D_dp_dsigma )
                        
            ! Calc the increment of the plastic multiplier
            dLambda = yield_surf%val / (dF_dsigma_D_dP_dsigma)

            ! Calc the update stress value
            stress = stress - dLambda * (D_dp_dsigma)

            ! Update the state params
                ! No update to do in this case

            ! Update the strain params
                ! No updates needed for this case

            ! Update the stress params
            call stress_params%update_stress_invariants(stress)

            ! Evaluate the yield surface
            call yield_surf%evaluate_surface(Xi, stress_params)

            ! increment the counter
            counter = counter + 1

        end do
        if (DEBUG) print *, "Counter", counter
    end subroutine ortiz_simo_assoc_mohr_coulomb

end module mod_UMAT_assoc_MC