module Mohr_Couloumb
    use kind_precision_module   , only: dp, i32
    use mod_mohr_coulomb_state_params , only: MohrCoulombStateParameters
    use mod_mohr_coulomb_stress_params, only: MohrCoulombStressParams
    use associative_mohr_coulomb          , only: MohrCoulombYieldSurface
    use mod_elastic_stress            , only: construct_stiffness_matrix, calc_elastic_stress

    implicit none
    
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
    !> \param[in]    material_properties Material properties (array).
    !> \param[out]   stiff_matrix        Stiffness matrix (6x6 components).
    subroutine Eval_Mohr_Couloumb(stress, stiff_matrix, strain_increment, state_vars, material_properties)

        real(kind = dp), intent(inout) :: stress(6), state_vars(:), material_properties(:)  
        real(kind = dp), intent(in)    :: strain_increment
        real(kind = dp), intent(out)   :: stiff_matrix(6, 6)
        
        ! Local variables
        class(MohrCoulombStateParameters) :: Xi
        class(MohrCoulombStressParams) :: stress_params
        class(MohrCoulombYieldSurface) :: yield_surf

        real(kind = dp)     :: yield_surf_tolerance = 1e-8 ! Make this an input parameter later
        integer(kind = i32) :: integration_scheme = 1      ! Make this an input
        integer(kind = i32), parameter :: ortiz_simo = 1                                             
        
        integer(kind = i32) :: max_iterations = 100 ! Make this an input parameter

        ! Flag for setting the integration scheme this can be an input later on
        integration_scheme = 1

        ! Set all of the possible integration schemes here
        ortiz_simo = 1
        ! Another one...
        ! Another one...

        ! Order is cohesion, phi (friction angle), shear modulus, poisson's ratio 
        Xi%initialize( state_vars(1), state_vars(2), props(1), props(2))  ! Initialize the state paramters object

        stress_params%initialize(stress)                                  ! Initialize the stress parameters object
        
        yield_surf%initialize(yield_surf_tolerance, max_iterations)       ! Initialize the yield surface object
        
        stiff_matrix = construct_stiffness_matrix(Xi%G, Xi%enu) ! Form the stiffness matrix

        select case (integration_scheme) ! If yielding select the integration method
            case(elastic)
                continue

            case(ortiz_simo)
                print *, "Do the evluation here"
                call ortiz_simo_assoc_mohr_coulomb( yield_surf, stress, Xi, stress_params, strain_increment )
                
            case default
                print *, "The integration method selected for the Associative Mohr Coulomb surface is not valid"
                print *, "Valid inputs are 1 (ortiz-simo)"

        end select

        ! Make sure no other variables need to be updated

        ! End the subroutine
        ! Init
        
    end subroutine Eval_Mohr_Couloumb

    subroutine ortiz_simo_assoc_mohr_coulomb(yield_surf, stress, Xi, stress_params, strain_increment, stiff_matrix)


        real(kind = dp), intent(inout) :: stress(6), strain_increment(6)
        real(kind = dp), intent(in)    :: stiff_matrix(6, 6)

        ! Local variables
        class(MohrCoulombStateParameters) :: Xi
        class(MohrCoulombStressParams) :: stress_params
        class(MohrCoulombYieldSurface) :: yield_surf
        integer(kind = i32) :: counter = 0 ! Init counter of iterations to zero
        real(kind = dp) :: dF_dsigma_D_dP_dsigma, D_dp_dsigma(6)
        
        ! Do the elastic prediction
        stress = calc_elastic_stress(strain, stiff_matrix, stress)

        ! Update the stress params
        stress_params%update_stress_invariants(stress)

        ! Update the state params
            ! Not needed for assoc MC

        ! Evaluate the yield surface
        yield_surf%evaluate_surface(Xi, stress_params)
        
        ! Update the stiffness matrix
            ! Don't need to do it in this case

        ! If the surface isn't yielding this statement won't activate
        do while (yield_surf%is_yielding .and. counter <= yield_surf%max_iterations)
        ! If yielding is occuring
            
            ! Evaluate the n vector (in this case n = m)
            yield_surf%evaluate_dF_dsigma(Xi, stress_params, stress)
            
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
            stress_params%update_stress_invariants(stress)

            ! Evaluate the yield surface
            yield_surf%evaluate_surface(Xi, stress_params)

            ! increment the counter
            counter = counter + 1

        end do
        
    end subroutine ortiz_simo_assoc_mohr_coulomb

end module Mohr_Couloumb