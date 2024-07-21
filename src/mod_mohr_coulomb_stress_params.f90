module class_assoc_MC_stress_params
   use kind_precision_module, only : dp
   use mod_stress_invariants, only : calc_mean_stress, calc_J_invariant, calc_lode_angle_s
   implicit none

   type :: assoc_MC_stress_params
      real(kind = dp) :: p
      real(kind = dp) :: J
      real(kind = dp) :: lode_angle
   contains
      procedure, pass(this) :: update_stress_invariants
   end type

   interface assoc_MC_stress_params
      module procedure initialize
   end interface 
contains

   ! Initialize the intial values of the state parameters using the stress matrix
   function initialize(sigma) result(this)
      real(kind = dp), intent(in)  :: sigma(6)
      type(assoc_MC_stress_params):: this

      call calc_potts_MC_stress_invariants(sigma, this%p, this%J, this%lode_angle)
   end function initialize

   subroutine update_stress_invariants(this, sigma)
      class(assoc_MC_stress_params), intent(inout) :: this
      real(kind = dp), intent(in) :: sigma(6)

      ! Evaluate and store the stress invariants
      call calc_potts_MC_stress_invariants(sigma, this%p, this%J, this%lode_angle)

   end subroutine update_stress_invariants

   ! Wrapper that calculates the stress invariants required for the
   ! Mohr-Coulomb in the form presented by Potts
   subroutine calc_potts_MC_stress_invariants(Sig, p, J, lode_angle)
      real(kind = dp), intent(in) :: Sig(6)
      real(kind = dp), intent(out) :: p, J, lode_angle

      ! Calc the mean stress
      p = calc_mean_stress(Sig)

      ! Calc the J invariants
      J = calc_J_invariant(Sig)

      ! Calc the lode angle
      ! TODO: Need to check this for bugs
      ! May need to add conditions that capture the case where the
      ! invariant calculations aren't exactly where they should be
      ! Alternatively should compare this with the results from
      ! the version using arctan and the principal stresses
      lode_angle = calc_lode_angle_s(Sig)

   end subroutine

end module class_assoc_MC_stress_params
