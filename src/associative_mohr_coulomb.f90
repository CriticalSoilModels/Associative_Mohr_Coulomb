module class_assoc_MC_yield_surf
   ! Have the imports here
   use class_assoc_MC_state_params, only: assoc_MC_state_params
   use class_assoc_MC_stress_params, only : assoc_MC_stress_params
   use kind_precision_module, only : dp, i32
   use mod_stress_invariants, only: calc_dp_dsigma, calc_dJ_dsigma, calc_dLodeAngle_bar_s_dsigma

   implicit none
   ! private
   ! public :: assoc_MC_yield_surf

   type :: assoc_MC_yield_surf
      real(kind = dp) :: tolerance           ! tolerance of the yield surface
      integer(kind = i32) :: max_iterations
      real(kind = dp) :: val                 ! Init value of the yield surface
      real(kind = dp) :: dF_dsigma(6)        ! Derivatives of the yield function wrt. the stress invariants
   contains
      procedure, pass(this) :: evaluate_surface
      procedure, pass(this) :: is_yielding
      procedure, pass(this) :: evaluate_dF_dsigma
      procedure, pass(this) :: evaluate_dF_dp
      procedure, nopass     :: calc_g_theta
      procedure, nopass     :: evaluate_dF_dJ
      procedure, nopass     :: evaluate_dF_dtheta
      ! procedure, nopass     :: eval_state_param_deriv
   end type assoc_MC_yield_surf

   interface assoc_MC_yield_surf
      module procedure initialize
   end interface

contains

   function initialize(tolerance, max_iterations, val, dF_dsigma) result(this)
      real(kind = dp)    , intent(in) :: tolerance, val, dF_dsigma
      integer(kind = i32), intent(in) :: max_iterations
      type(assoc_MC_yield_surf) :: this

      ! Set the tolerance of the yield function
      this%tolerance = tolerance

      ! Set the max number of iterations for the stress integration
      this%max_iterations = max_iterations

      this%val = val ! Set the initial value to zero
      this%dF_dsigma = dF_dsigma
   end function initialize

   function calc_g_theta(Xi, stress_var) result(g_theta)
      ! Calculates g(\theta) from Potts and Zdravkovic (Finite element analysis in geotechnical engineering. 2)

      class(assoc_MC_state_params), intent(in)    :: Xi
      class(assoc_MC_stress_params)   , intent(in)    :: stress_var
      logical :: DEBUG = .False.

      ! Local variables
      real :: g_theta_numerator, g_theta_denominator, g_theta
      
      
      ! Evaluate g(theta)
      g_theta_numerator   = sin(Xi%phi)
      g_theta_denominator = cos(stress_var%lode_angle) + sin(stress_var%lode_angle) * sin(Xi%phi) / sqrt(3.0_dp)

      ! Calc g_theta - This is the variable that is returned
      g_theta = g_theta_numerator / g_theta_denominator
      
      if (DEBUG) then
         print *, "g_theta_num: "  , g_theta_numerator
         print *, "g_theta_denom: ", g_theta_denominator
         print *, "lode angle   : ", stress_var%lode_angle
         print *, "fric_angle   : ", Xi%phi
         print *, "cos value    : ", cos(stress_var%lode_angle)
         print *, "sin terms    : ", sin(stress_var%lode_angle) * sin(Xi%phi) / sqrt(3.0_dp)
      end if
   end function

   ! Evaluate the value of the yield function

   ! Eqn:
   !   J - (c'/tan(phi') + p') g(\theta)
   !   g(theta) = \frac{ sin(phi') }{ cos(\theta) + \frac{sin(theta) sin(phi)}{ \sqrt{3} } }
   subroutine evaluate_surface(this, Xi, stress_var)
      class(assoc_MC_yield_surf)   , intent(inout) :: this
      class(assoc_MC_state_params), intent(in)    :: Xi
      class(assoc_MC_stress_params)   , intent(in)    :: stress_var

      logical :: DEBUG = .FALSE.

      ! Local variables
      real :: g_theta

      g_theta =  this%calc_g_theta(Xi, stress_var)
      
      ! Store the value of the yield function in the val param
      this%val = stress_var%J - ( Xi%c / tan(Xi%phi )  + stress_var%p ) * g_theta

      if (DEBUG) then
         print *, "g_theta: ", g_theta
         print *, "F val  : ", this%val
      end if
   end subroutine

   function is_yielding(this) result(outside_surface)
      class(assoc_MC_yield_surf), intent(in) :: this
      logical :: outside_surface
      ! Returns true if yielding (ie. F > tolerance)
      ! Returns false if elastic regime (F < tolerance)

      if ( this%val > this%tolerance ) then
         outside_surface = .True.
      end if

   end function is_yielding

   ! Evalute the yield function derivatives wrt. to stress
   ! Stores the value of the derivatives inside of the function
   subroutine evaluate_dF_dsigma(this, Xi, stress_var, stress)
      class(assoc_MC_yield_surf), intent(inout) :: this
      class(assoc_MC_state_params), intent(in) :: Xi
      class(assoc_MC_stress_params)   , intent(in) :: stress_var
      real(kind = dp), intent(in) :: stress(6)
      
      ! Local variables
      real(kind = dp) :: dF_dp, dF_dJ, dF_dtheta 
      real(kind = dp) :: dp_dsigma(6), dJ_dsigma(6), dtheta_dsigma(6)

      ! Evaluate the derivative of the yield surface wrt. the mean stress
      dF_dp = this%evaluate_dF_dp(Xi, stress_var)

      ! Evaluate the derivative of the yield surface wrt. the deviatoric stress
      dF_dJ = this%evaluate_dF_dJ(Xi, stress_var)

      ! Evaluate the derivative of the yield surface wrt. the lode angle
      dF_dtheta = this%evaluate_dF_dtheta(Xi, stress_var)

      ! Evaluate the invariant derivatives

      ! Calc the deriv of p wrt. to stress (No inputs for this function)
      dp_dsigma = calc_dp_dsigma()

      ! Calc the deriv of J wrt. stress
      dJ_dsigma = calc_dJ_dsigma(stress)

      ! Calc the deriv of the lode angle wrt. stress
      dtheta_dsigma = calc_dLodeAngle_bar_s_dsigma(stress)

      ! Evaluate the complete derivative using the product and chain rules
      this%dF_dsigma = dF_dp * dp_dsigma + dF_dJ * dJ_dsigma + dF_dtheta * dtheta_dsigma
      
   end subroutine evaluate_dF_dsigma

   ! Calc the derivative of the yield surface wrt. to the mean stress
   function evaluate_dF_dp(this, Xi, stress_var) result(dF_dp)
      class(assoc_MC_yield_surf)   , intent(in) :: this
      class(assoc_MC_state_params), intent(in) :: Xi
      class(assoc_MC_stress_params)   , intent(in) :: stress_var

      ! Local variable
      real :: dF_dp

      ! Calc the derivative of the yield function wrt. to the mean effective stress
      dF_dp = -1.0 * this%calc_g_theta(Xi, stress_var)

   end function evaluate_dF_dp

   ! Calc the derivative of the yield function wrt. to the deviatoric stres
   function evaluate_dF_dJ(Xi, stress_var) result(dF_dJ)
      ! Might not need to pass this stuff, just passing it because it's general
      class(assoc_MC_state_params), intent(in) :: Xi
      class(assoc_MC_stress_params)   , intent(in) :: stress_var

      ! Local variable
      real :: dF_dJ

      ! dF_dJ = 1.0 for MC
      dF_dJ = 1.0
   end function evaluate_dF_dJ

   ! Calculates the derivative of the yield surface with respect to the Lode angle (\theta).
   !
   ! Parameters:
   ! Xi          - Input, type(class(assoc_MC_state_params)), intent(in) :: Xi
   !               Object containing material state parameters (cohesion, friction angle).
   ! stress_var  - Input, type(class(assoc_MC_stress_params)), intent(in) :: stress_var
   !               Object containing stress parameters (pressure, Lode angle).
   !
   ! Returns:
   ! dF_dtheta   - Output, real :: dF_dtheta
   !               Derivative of the yield surface with respect to the Lode angle (\theta).
   !
   ! Local Variables:
   ! first_term         - Real :: First term in the derivative calculation.
   ! second_term_numer  - Real :: Numerator of the second term.
   ! second_term_denom  - Real :: Denominator of the second term.
   ! thrid_term         - Real :: Third term in the derivative calculation.
   ! dF_dtheta          - Real :: Resulting derivative of the yield surface.
   !
   ! Details:
   ! This function computes the derivative of the yield surface with respect to the Lode angle (\theta)
   ! based on the Mohr-Coulomb yield criterion. It utilizes material state parameters and stress variables
   ! to evaluate the derivative expression, involving cohesion, friction angle, pressure, and Lode angle.
   !
   ! Formula:
   ! dF_dtheta = [ (c' / tan( \phi' ) + p') * sin( \phi' ) / [ cos( \theta + sin( \theta ) * sin( \phi' ) / \sqrt(3.0) ) ]^2 ] *
   !             * [ sin( \theta ) - cos( \theta ) * sin( \phi' ) / sqrt(3.0) ]
   !
   ! Reference:
   ! Potts, D. M., & Zdravkovic, L. (1999). Finite element analysis in geotechnical engineering.
   ! Chapter 2: Mohr-Coulomb plasticity, yield surfaces.
   !
   function evaluate_dF_dtheta(Xi, stress_var) result(dF_dtheta)
      class(assoc_MC_state_params), intent(in) :: Xi
      class(assoc_MC_stress_params)   , intent(in) :: stress_var

      ! Local variables
      real :: first_term, second_term_numer, second_term_denom, thrid_term, dF_dtheta

      ! (c' / tan( \phi' ) + p')
      first_term = Xi%c/tan(Xi%phi) + stress_var%p

      ! sin( \phi' )
      second_term_numer = sin(Xi%phi)

      ! [ cos( \theta +  sin( \theta ) sin( \phi' ) / \sqrt(3.0) ) ]^{2}
      second_term_denom = cos(stress_var%lode_angle) + sin(stress_var%lode_angle) * sin(Xi%phi) / sqrt(3.0)

      ! sin ( \phi' ) - cos( \theta ) sin( \phi' ) / sqrt(3.0)
      thrid_term = sin(stress_var%lode_angle) - cos(stress_var%lode_angle) * sin(Xi%phi) / sqrt(3.0)

      dF_dtheta  = first_term * second_term_numer/second_term_denom**2 * thrid_term

   end function

end module class_assoc_MC_yield_surf
