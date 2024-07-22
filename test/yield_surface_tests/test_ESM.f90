program main
  use kind_precision_module, only: dp, i32
  use mod_UMAT_assoc_MC, only: UMAT_assoc_MC
  use mod_ESM_assoc_MC, only: UMAT

  implicit none
  
  ! Yield surface values
  ! real(kind = dp)    , parameter :: tolerance = 1e-8, max_iterations

  real(kind = dp) :: stress(6), strain_increment(6)
  real(kind = dp) :: phi, c, G, enu
  real(kind = dp) :: state_vars(2), props(2)
  real(kind = dp) :: stiff_matrix(6,6)
  real(kind = dp) :: eunloading 
  real(kind = dp)  :: plasticmultiplier, additionalvar(1)
  integer(kind = i32), parameter :: npt = 0, noel = 0, idset = 0, numberofphases = 0, &
                                    ntens = 0, naddvar = 0, nstatev = 0, nprops = 0
  character(3) :: cmname
  cmname = "Bob"

  stiff_matrix(:, :) = 0.0_dp
  
  ! Initialize state parameters
  c   = 10.0 ! cohesion
  phi = 30.0 ! friction angle, Expects value in degrees internally changes to radians
  G   = 50   ! Shear modulus
  enu = 0.0 ! Poisson's ratio

  state_vars = [c, phi]
  props      = [G, enu] 

  ! Init the mock stress state
  stress = [20.0_dp, 10.0_dp, 10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  strain_increment = [0.5_dp, 0.1_dp, 0.001_dp, 10.0_dp, 0.00_dp, 0.000_dp]
  
  call UMAT(npt, noel, idset, stress, eunloading, plasticmultiplier, strain_increment, nstatev, &
                    state_vars, naddvar, additionalvar, cmname, nprops, props, numberofphases, ntens)

end program main
