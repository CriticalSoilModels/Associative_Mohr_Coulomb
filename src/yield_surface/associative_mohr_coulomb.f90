module associative_mohr_coulomb
   ! Have the imports here
   use mod_mohr_coulomb_state_params, only: MohrCoulombStateParameters
   use mod_mohr_coulomb_stress_params, only : MohrCoulombStressParams
   use kind_precision_module, only : dp
   use integer_precision_module, only : i32
   
   implicit none
   private
   public :: MohrCoulombYieldSurface

   type :: MohrCoulombYieldSurface
      real(kind = dp) :: tolerance               ! tolerance of the yield surface
      integer(kind = i32) :: max_iterations
      real(kind = dp) :: val                     ! Init value of the yield surface
      real(kind = dp) :: dF_dp, dF_dJ, dF_dtheta ! Derivatives of the yield function wrt. the stress invariants
   contains
      procedure, pass(self) :: evaluate_surface
      procedure, pass(self) :: is_yielding
      procedure, nopass     :: calc_g_theta
      procedure, pass(self) :: evaluate_stress_deriv
      procedure, pass(self) :: evaluate_dF_dp
      procedure, nopass     :: evaluate_dF_dJ
      procedure, nopass     :: evaluate_dF_dtheta
      ! procedure, nopass     :: eval_state_param_deriv
   end type MohrCoulombYieldSurface

contains

   subroutine initialize(self, tolerance)
      class(MohrCoulombYieldSurface), intent(inout) :: self
      real(kind = dp)    , intent(in) :: tolerance

      ! Set the tolerance of the yield function
      self%tolerance = tolerance
      
      self%val = 0.0_dp ! Set the initial value to zero
   end subroutine initialize

   function calc_g_theta(Xi, stress_var) result(g_theta)
      ! Calculates g(\theta) from Potts and Zdravkovic (Finite element analysis in geotechnical engineering. 2)

      class(MohrCoulombStateParameters), intent(in)    :: Xi
      class(MohrCoulombStressParams)   , intent(in)    :: stress_var

      ! Local variables
      real :: g_theta_numerator, g_theta_denominator, g_theta

      ! Evaluate g(theta)
      g_theta_numerator   = sin(Xi%phi)
      g_theta_denominator = cos(stress_var%lode_angle) + sin(stress_var%lode_angle) * sin(Xi%phi) / sqrt(3.0)

      ! Calc g_theta - This is the variable that is returned
      g_theta = g_theta_numerator / g_theta_denominator

   end function

   ! Evaluate the value of the yield function

   ! Eqn:
   !   J - (c'/tan(phi') + p') g(\theta)
   !   g(theta) = \frac{ sin(phi') }{ cos(\theta) + \frac{sin(theta) sin(phi)}{ \sqrt{3} } }
   subroutine evaluate_surface(self, Xi, stress_var)
      class(MohrCoulombYieldSurface)   , intent(inout) :: self
      class(MohrCoulombStateParameters), intent(in)    :: Xi
      class(MohrCoulombStressParams)   , intent(in)    :: stress_var

      ! Local variables
      real :: g_theta

      g_theta =  self%calc_g_theta(Xi, stress_var)

      ! Store the value of the yield function in the val param
      self%val = stress_var%J - ( Xi%c / tan(Xi%phi )  + stress_var%p ) * g_theta
   end subroutine

   function is_yielding(self) result(outside_surface)
      class(MohrCoulombYieldSurface), intent(in) :: self
      logical :: outside_surface
      ! Returns true if yielding (ie. F > tolerance)
      ! Returns false if elastic regime (F < tolerance)

      if ( self%val > self%tolerance ) then
         outside_surface = .True.
      end if

   end function is_yielding

   ! Evalute the yield function derivatives wrt. to stress
   ! Stores the value of the derivatives inside of the function
   subroutine evaluate_stress_deriv(self, Xi, stress_var)
      class(MohrCoulombYieldSurface), intent(inout) :: self
      class(MohrCoulombStateParameters), intent(in) :: Xi
      class(MohrCoulombStressParams)   , intent(in) :: stress_var

      ! Evaluate the derivative of the yield surface wrt. the mean stress
      self%dF_dp = self%evaluate_dF_dp(Xi, stress_var)

      ! Evaluate the derivative of the yield surface wrt. the deviatoric stress
      self%dF_dJ = self%evaluate_dF_dJ(Xi, stress_var)

      ! Evaluate the derivative of the yield surface wrt. the lode angle
      self%dF_dtheta = self%evaluate_dF_dtheta(Xi, stress_var)
   end subroutine evaluate_stress_deriv

   ! Calc the derivative of the yield surface wrt. to the mean stress
   function evaluate_dF_dp(self, Xi, stress_var) result(dF_dp)
      class(MohrCoulombYieldSurface)   , intent(in) :: self
      class(MohrCoulombStateParameters), intent(in) :: Xi
      class(MohrCoulombStressParams)   , intent(in) :: stress_var

      ! Local variable
      real :: dF_dp

      ! Calc the derivative of the yield function wrt. to the mean effective stress
      dF_dp = -1.0 * self%calc_g_theta(Xi, stress_var)

   end function evaluate_dF_dp

   ! Calc the derivative of the yield function wrt. to the deviatoric stres
   function evaluate_dF_dJ(Xi, stress_var) result(dF_dJ)
      ! Might not need to pass this stuff, just passing it because it's general
      class(MohrCoulombStateParameters), intent(in) :: Xi
      class(MohrCoulombStressParams)   , intent(in) :: stress_var

      ! Local variable
      real :: dF_dJ

      ! dF_dJ = 1.0 for MC
      dF_dJ = 1.0
   end function evaluate_dF_dJ

   ! Calculates the derivative of the yield surface with respect to the Lode angle (\theta).
   !
   ! Parameters:
   ! Xi          - Input, type(class(MohrCoulombStateParameters)), intent(in) :: Xi
   !               Object containing material state parameters (cohesion, friction angle).
   ! stress_var  - Input, type(class(MohrCoulombStressParams)), intent(in) :: stress_var
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
      class(MohrCoulombStateParameters), intent(in) :: Xi
      class(MohrCoulombStressParams)   , intent(in) :: stress_var

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

end module associative_mohr_coulomb
