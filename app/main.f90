program main
  use mod_stress_invariants, only: calc_mean_stress, calc_stress_invariants
  use associative_mohr_coulomb, only: MohrCoulombYieldSurface
  use mod_mohr_coulomb_state_params, only: MohrCoulombStateParameters
  use mod_mohr_coulomb_stress_params, only: MohrCoulombStressParams

  implicit none
  
  ! real :: stress(6), mean_stress, dev_stress, lode_angle

  ! stress(:) = 1.0

  ! mean_stress = calc_mean_stress(stress)

  ! print *, mean_stress

  ! call calc_stress_invariants(stress, mean_stress, dev_stress, lode_angle)

  ! print *, "Mean Stress", mean_stress

  ! print *, "Deviatoric Stress", dev_stress

  ! print *, "Lode angle", lode_angle

  type(MohrCoulombStateParameters) :: Xi
  type(MohrCoulombStressParams) :: stress_var
  type(MohrCoulombYieldSurface) :: yield_surface

  ! Initialize parameters
  Xi%fric_angle = 30.0
  Xi%cohesion = 10.0
  stress_var%theta = 45.0
  stress_var%q = 20.0
  stress_var%p = 5.0

  ! Evaluate the yield surface
  call yield_surface%evaluate_surface(Xi, stress_var)

  ! Print the result
  print *, "Yield Surface Value:", yield_surface%val

end program main