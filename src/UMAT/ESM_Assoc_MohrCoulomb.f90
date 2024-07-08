module Mohr_Couloumb

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
        use kind_precision_module, only: dp
        use integer_precision_module, only: i32
        use mod_mohr_coulomb_state_params , only: MohrCoulombStateParameters
        use mod_mohr_coulomb_stress_params, only: MohrCoulombStressParams
        use associative_mohr_coulomb      , only: MohrCoulombYieldSurface
        use mod_elastic_stress            , only: construct_stiffness_matrix, calc_elastic_stress

        real(kind = dp), intent(inout) :: stress(6), state_vars(:), material_properties(:)  
        real(kind = dp), intent(in)    :: strain_increment
        real(kind = dp), intent(out)   :: stiff_matrix(6, 6)
        
        ! Local variables
        class(MohrCoulombStateParameters) :: Xi
        class(MohrCoulombStressParams) :: stress_params
        class(MohrCoulombYieldSurface) :: yield_surf

        real(kind = dp)     :: yield_surf_tolerance = 1e-8 ! Make this an input parameter later
        integer(kind = i32) :: integration_scheme = 1      ! Make this an input
        integer(kind = i32), parameter ::   elastic = -1, & ! This flag will be used when elastic loading is happening
                                            ortiz_simo = 1
        
        integer(kind = i32) :: max_iterations = 100 ! Make this an input parameter

        ! Flag for setting the integration scheme this can be an input later on
        integration_scheme = 1

        ! Set all of the possible integration schemes here
        ortiz_simo = 1
        ! Another one...
        ! Another one...

        ! Order is cohesion, phi (friction angle), shear modulus, poisson's ratio 
        Xi%initialize( state_vars(1), state_vars(2), props(1), props(2))  ! Initialize the state paramters object

        stress_params%initialize(stress)                        ! Initialize the stress parameters object
        
        yield_surf%initialize(yield_surf_tolerance, max_iterations)             ! Initialize the yield surface object
        
        stiff_matrix = construct_stiffness_matrix(Xi%G, Xi%enu) ! Form the stiffness matrix
        
        ! Might want to move the stress into the stress params and make this function inside
        stress = calc_elastic_stress(strain_increment, stiff_matrix, stress) ! update the stress

        stress_params%update_stress_invariants(stress)    ! Update the stress variables

        yield_surf%evaluate_surface(Xi, stress_params)    ! Evaluate the yield surface

        if not (yield_surf%is_yielding) integration = -1  ! If the material is not yielding don't need to do the integration

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
        use kind_precision_module   , only: dp
        use integer_precision_module, only: i32
        use mod_mohr_coulomb_state_params , only: MohrCoulombStateParameters
        use mod_mohr_coulomb_stress_params, only: MohrCoulombStressParams
        use associative_mohr_coulomb          , only: MohrCoulombYieldSurface
        use mod_elastic_stress            , only: construct_stiffness_matrix, calc_elastic_stress

        class(MohrCoulombStateParameters), intent(in) :: Xi         ! Don't need to update the state paramters for this model
        class(MohrCoulombYieldSurface), intent(inout) :: yield_surf ! Yield surface
        class(MohrCoulombStressParams), intent(inout) :: stress_params

        real(kind = dp), intent(inout) :: stress(6), strain_increment(6)
        real(kind = dp), intent(in)    :: stiff_matrix(6, 6)

        ! Local variables
        integer(kind = i32) :: counter = 0 ! Init counter of iterations to zero

        ! Elastic prediction has already been done

        do while (yield_surf%is_yielding .and. counter <= yield_surf%max_iterations)
            ! Calc the increment of the plastic multiplier

            ! Calc the update stress value

            ! Calc the updated stress value

            ! 


        end do


        
    end subroutine ortiz_simo_assoc_mohr_coulomb
end module Mohr_Couloumb