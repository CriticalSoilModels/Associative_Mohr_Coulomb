module mod_stress_invariants
   use kind_precision_module, only : dp, i32

   use mod_general_voigt, only: multiply_voigt_vectors, get_3d_voigt_identity_vector, trace_voigt_vector

   implicit none

contains
   !> Calculates the mean stress from a given stress vector.
   pure function calc_mean_stress(stress) result(mean_stress)
      real(kind = dp), intent(in) :: stress(:)
      real(kind = dp) :: mean_stress
      ! Local variables
      real(kind = dp), parameter :: three = 3.0_dp
      
      mean_stress = sum(stress(1:3)) / three

   end function calc_mean_stress

   pure function calc_q_invariant(stress) result(q)
      real(kind = dp), intent(in) :: stress(6)
      
      real(kind = dp) :: q
      
      ! Local variables
      real(kind = dp) :: J2
      real(kind = dp), parameter :: three = 3.0_dp
 
      ! Calc the J2 invariant
      J2 = calc_J2_invariant(stress)

      q = sqrt( three * J2 )

   end function calc_q_invariant

   ! Calculates the J2 stress invariant
   pure function calc_J2_invariant(stress) result(J2)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: J2
      
      ! Local variables
      real(kind = dp) :: mean_stress
      real(kind = dp) :: stress_copy(6)
      real(kind = dp), parameter :: one_half = 0.5, two = 2.0
      
      ! Make a local copy of the stress vector
      stress_copy = stress

      mean_stress = calc_mean_stress(stress_copy)

      ! Generate the deviatoric_stress_tensor
      ! Need to operate on the diagonal of the stress tensor
      ! Which is the first 3 elements in voigt notation

      ! s_{ij} = \sigma_{ij} - \sigma_{kk}/3 \delta_{ij}
      stress_copy(1:3) = stress_copy(1:3) - mean_stress
      
      ! Square all the terms
      stress_copy = stress_copy**2
      
      ! Multiply the shear terms by two since there's two sets of them
      stress_copy(4:6) = two * stress_copy(4:6)

      ! J2 = 1/2 s_{ij} s_{ji}
      J2 = one_half *  sum(stress_copy)
      
   end function calc_J2_invariant

   pure function calc_s_determinant(stress) result(s_det)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: s_det

      ! Local variables
      real(kind = dp) :: s(6), mean_stress
      real(kind = dp) :: first_term, second_term, third_term, fourth_term, &
                         fifth_term
      real(kind = dp), parameter :: one = 1.0, two = 2.0
      real(kind = dp) :: stress_copy(6) ! Used to make a local copy of the stress vector

      ! Make a copy just to be safe
      stress_copy = stress

      ! Calc the mean stress
      mean_stress = calc_mean_stress(stress_copy)

      ! Store the stress matric in s 
      s(:) = stress_copy(:)
      
      ! Calc s by subtracting the mean stress from the diagonal
      s(1:3) = stress_copy(1:3) - mean_stress

      ! Calc (sigma_x - p) (sigma_y - p) (sigma_z - p)
      first_term = product( s(1:3) )
      
      ! Calc -1.0 * (sigma_x - p) * tau_{yz}^{2}
      second_term = -one * s(1) * s(5)**2
      
      ! Calc -1.0 * (sigma_y - p) * tau_{zx}^{2}
      third_term  = -one * s(2) * s(6)**2
      
      ! Calc (sigma_z - p) * tau_{xy}^{2}
      fourth_term = -one * s(3) * s(4)**2

      ! Calc 2 tau_{xy} tau_{yz} tau_{zx}
      fifth_term = two * product( s(4:6) )
      
      ! Sum all of the terms
      s_det = first_term + second_term + third_term + fourth_term + fifth_term

   end function

   pure function calc_J_invariant(stress) result(J)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: J

      ! Local variables
      real(kind = dp) :: J2

      ! Calc the J2 stress invariant
      J2 = calc_J2_invariant(stress)

      ! J = \sqrt{J2}
      J = sqrt(J2)
   end function calc_J_invariant

   pure function calc_lode_angle_s(stress) result(lode_angle)
      ! TODO: The lode angle should be bounded between
      ! -pi/6 <= \theta <= pi / 6
      ! -pi/6 from trx compression: sigma_1 >= sigma_2 = sigma_3
      ! pi/6 from trx extension   : sigma_1 = sigma_2 >= sigma_3
      ! 0 from shear              : sigma_2 = (sigma_1 + sigma_3)/2

      ! This equation comes from Potts and Zdravković

      ! Sin is an odd function so $ -sin(x) = sin(-x)$
      ! Eqn:
         ! $$ sin( -3 \bar{\theta}_{s} ) = \frac{ 3\sqrt{3} }{ 2 } \frac{ det(s) }{ J^{3} }$$
      ! This is the same equation as the {\theta}_{s} value on wikipedia, link:
      ! https://en.wikipedia.org/wiki/Lode_coordinates#Lode_angle_%E2%80%93_angular_coordinate_%7F'%22%60UNIQ--postMath-0000002D-QINU%60%22'%7F

      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: lode_angle

      ! Local variables
      real(kind = dp) :: J, det_s, inside
      real(kind = dp), parameter :: one = 1.0, two = 2.0, three = 3.0
      real(kind = dp), parameter :: PI = atan(1.0) * 4.0 ! Calc pi
      real(kind = dp), parameter :: tolerance = 1e-8

      ! TODO: Need to update this function to catch the edge cases

      ! Calc J
      J = calc_J_invariant(stress)
      
      ! Calc the determinant of s
      det_s = calc_s_determinant(stress)

      ! calc the inside of the paranthesis
      ! 3 \sqrt{3} / 2 * det(s) / J^{3}

      if (abs(J) < tolerance) then
         inside = 0.0_dp
      else
         ! J isn't zero
         inside = (three * sqrt(three) / two) * det_s / J**3
      end if 

      ! Catch if the inside is slightly less than negative one
      if (inside > -1.0_dp - tolerance .and. inside < -1.0_dp + tolerance) then
         inside = -1.0_dp
      else if (inside > 1.0_dp - tolerance .and. inside < 1.0_dp + tolerance) then
         ! Check if inside is slighty larger than one
         inside =  1.0_dp
      end if
      ! Calc the lode angle
      lode_angle = - one * asin(inside) / 3.0

      if (lode_angle < -PI/6.0_dp .or. lode_angle > PI/6.0_dp) then
         ! Stop the program if this thing is outside the right image
         ! See for info about the domain https://en.wikipedia.org/wiki/Lode_coordinates
         error stop
      end if
   end function calc_lode_angle_s
   
   subroutine calc_stress_invariants(Sig, p, q, theta)
      !*********************************************************************
      ! Takes the stress tensor Sig and return invariants p, q, and theta  *
      !																	 *
      !*********************************************************************
      !input variables
      real(kind = dp), dimension(6), intent(in):: Sig
      !output variables
      real(kind = dp), intent(out)::p, q, theta

      !local variables
      real(kind = dp):: dev(6), J2, J3, sin3theta

      p= calc_mean_stress(Sig) !mean stress

      dev=Sig
      dev(1)=dev(1)-p !computes deviatoric stress tensor
      dev(2)=dev(2)-p
      dev(3)=dev(3)-p

      J2 = norm2(dev)

      J2=(J2**2)/2.0 !J_2 invariant
      q=sqrt(3*J2) ! deviatoric stress

      !J3 stress invariant
      J3=dev(1)*dev(2)*dev(3)-dev(1)*dev(6)**2-dev(2)*dev(4)**2-dev(3)*dev(5)**2+2.0*dev(4)*dev(5)*dev(6)

      !sin3theta
      if (J2>0.0d0) then
         sin3theta=0.5*J3*(3.0/J2)**(1.5d0)
      else !Assume triaxial compression
         sin3theta=-1.0d0
      endif
      if (sin3theta<-1.0) sin3theta=-1.0d0
      if (sin3theta>1.0) sin3theta=1.0d0


      theta=-asin(sin3theta)/3.0d0 !Lode's angle

   end subroutine calc_stress_invariants

   ! Calc the derivative of the p invariant wrt. the stress (sigma)
   pure function calc_dp_dsigma() result(dp_dsigma)
      real(kind = dp) :: dp_dsigma(6)

      dp_dsigma = 1.0_dp / 3.0_dp * [1, 1, 1, 0, 0, 0]

   end function calc_dp_dsigma

   pure function calc_dJ_dsigma(stress) result(dJ_dsigma)
      ! Eqn:
         ! 1 / ( 2 J ) [\sigma_{x} - p, \sigma_{y} - p, \sigma_{z} - p,...
         !               2 \tau_{xy} 2 \tau_{yz} 2 \tau_{zx}]
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: dJ_dsigma(6)

      ! Local variables
      real(kind = dp) :: J, mean_stress, s(6)

      real(kind =dp), parameter :: tolerance = 1e-8

      ! Calc the J invariant
      J = calc_J_invariant(stress)
      
      ! Calc the mean stress
      mean_stress = calc_mean_stress(stress)

      ! Copy the values of the stress
      s(:) = stress(:)

      ! Calc the deviatoric stress matrix
      s(1:3) = stress(1:3) - mean_stress

      ! Scale the shear terms by two
      ! 2 \tau_{xy} 2 \tau_{yz} 2 \tau_{zx}
      s(4:6) = 2.0_dp * s(4:6)

      ! Check that J isn't close to being zero or is zero
      if (abs(J) < tolerance) J = tolerance

      dJ_dsigma = 1.0_dp / (2.0_dp * J) * s
   
   end function calc_dJ_dsigma

   pure function calc_dLodeAngle_bar_s_dsigma(stress) result(dLodeAngle_bar_s_dsigma)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: dLodeAngle_bar_s_dsigma(6)

      ! Local variables
      real(kind = dp) :: outside_term, inside_term(6)     ! Place holder variables
      real(kind = dp) :: lode_angle, J, det_s             ! Required invariants
      real(kind = dp) :: dJ_dsigma(6), dJ3_dsigma(6)      ! Required invariant derivatives
      real(kind = dp), parameter :: tolerance = 1e-8

      ! Calc the required invariants
      lode_angle = calc_lode_angle_s(stress)
      J = calc_J_invariant(stress)
      
      det_s = calc_s_determinant(stress)

      ! Calc the required derivatives
      dJ_dsigma = calc_dJ_dsigma(stress)
      dJ3_dsigma = calc_potts_dJ3_dsigma(stress)

      ! Check the denominators
      if (abs(J) < tolerance) J = tolerance
      if (abs(lode_angle) < tolerance) lode_angle = tolerance

      ! Calc the terms of the equation
      ! Term in front of the paranthesis
      outside_term = sqrt( 3.0_dp ) / (2.0_dp * cos(3.0_dp * lode_angle) * J**3)
      
      inside_term  = ( det_s / J ) * dJ_dsigma !- dJ3_dsigma

      ! Term inside the parenthesis

      ! Calc the derivative
      dLodeAngle_bar_s_dsigma = outside_term * inside_term

   end function calc_dLodeAngle_bar_s_dsigma
   
   ! Calc the derivative of the J3 Stress invariant wrt/ the stress voigt vector
   pure function calc_potts_dJ3_dsigma(stress) result(dJ3_dsigma)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: dJ3_dsigma(6)

      ! Local variables
      real(kind = dp), parameter :: ONE_NINTH = 1.0_dp/9.0_dp, &
                                    ONE_THIRD = 1.0_dp/3.0_dp, &
                                    TWO       = 2.0_dp       , &
                                    FOUR      = 4.0_dp
      real(kind = dp) :: sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_zx

      ! Unpack the stress value to make it easier
      sigma_x = stress(1)
      sigma_y = stress(2)
      sigma_z = stress(3)
      tau_xy  = stress(4)
      tau_yz  = stress(5)
      tau_zx  = stress(6)

      ! Calc the dJ3/dsigma_x term
      dJ3_dsigma(1) = ONE_NINTH * ( TWO*sigma_x**2 - TWO*sigma_x*sigma_y      &
                                   -  TWO*sigma_x*sigma_z - sigma_y**2        &
                                   + FOUR*sigma_y*sigma_z - sigma_z**2 )      &
                    + ONE_THIRD * ( tau_xy**2 - TWO*tau_yz**2 + tau_zx**2)
      
      ! Calc the dJ3/dsigma_y term
      dJ3_dsigma(2) = ONE_NINTH * ( -sigma_x**2 - TWO*sigma_x*sigma_y         &
                                    -  TWO*sigma_y*sigma_z - sigma_z**2       &
                                    + FOUR*sigma_x*sigma_z + TWO*sigma_y**2)  &
                    + ONE_THIRD * ( tau_xy**2 + tau_yz**2 - TWO*tau_zx**2)
      
      ! Calc the dJ3/dsigma_z
      dJ3_dsigma(3) = ONE_NINTH * ( -sigma_x**2 + FOUR*sigma_x*sigma_y        &
                                 - TWO*sigma_x*sigma_z - sigma_y**2           &
                                 - TWO*sigma_y*sigma_z + TWO*sigma_z**2 )     &
                    + ONE_THIRD * ( -TWO*tau_xy**2 + tau_yz**2 + tau_zx**2)

      ! Calc dJ3/tau_xy
      dJ3_dsigma(4) = TWO*tau_yz*tau_zx + TWO*tau_xy* ONE_THIRD * (sigma_x + sigma_y - TWO*sigma_z) 
      
      ! Calc dJ3/tau_yz
      dJ3_dsigma(5) = TWO*tau_xy*tau_zx + TWO*tau_yz* ONE_THIRD * (-TWO*sigma_x + sigma_y + sigma_z)

      ! Calc dJ3/tau_zx
      dJ3_dsigma(6) = TWO*tau_xy*tau_yz + TWO*tau_zx* ONE_THIRD * (sigma_x - TWO*sigma_y + sigma_z)
      
   end function calc_potts_dJ3_dsigma

end module mod_stress_invariants