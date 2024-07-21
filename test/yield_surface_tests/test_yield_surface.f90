program test_yield_func
    use kind_precision_module, only: dp, i32
    use class_assoc_MC_yield_surf, only: assoc_MC_yield_surf
    use class_assoc_MC_state_params, only: assoc_MC_state_params
    use class_assoc_MC_stress_params, only: assoc_MC_stress_params
  
    implicit none
    
    type(assoc_MC_yield_surf)    :: yield_surface
    type(assoc_MC_state_params) :: Xi
    type(assoc_MC_stress_params)    :: stress_var
  
    ! Yield surface values
    real(kind = dp)    , parameter :: tolerance = 1e-8
    integer(kind = i32), parameter :: max_iterations = 1000
  
    real(kind = dp) :: stress(6)
    real(kind = dp) :: phi, c, G, enu
    real(kind = dp) :: lode_angle
  
    ! Initialize state parameters
    c   = 10.0 ! cohesion
    phi = 30.0 ! friction angle, Expects value in degrees internally changes to radians
    G   = 50   ! Shear modulus
    enu = 0.0 ! Poisson's ratio
  
    ! initialize the object
    Xi = assoc_MC_state_params(c, phi, G, enu)
  
    ! print *, "State param info", Xi%c, Xi%phi, Xi%G, Xi%enu
    ! Initialize the yield surface
    yield_surface = assoc_MC_yield_surf(tolerance, max_iterations, 0.0, 0.0)
  
    ! ! Init the mock stress state
    stress = [20.0_dp, 10.0_dp, 10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  
    ! Init the stress variables
    stress_var = assoc_MC_stress_params(stress)
  
    print *, "Stress values", stress_var%J, stress_var%lode_angle, stress_var%p
  
    ! Evaluate the yield surface
    call yield_surface%evaluate_surface(Xi, stress_var)
  
    ! ! Print the result
    print *, "Yield Surface Value:", yield_surface%val
  
  end program test_yield_func
  