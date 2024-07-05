module mod_stress_invariants
   use kind_precision_module, only : dp

   implicit none

contains
   !> Calculates the mean stress from a given stress tensor.
   !!
   !! Computes the mean stress as the average of the first three components
   !! of the input stress tensor.
   !!
   !! @param[in]  stress  Input stress tensor of length 6.
   !! @return     mean_stress  Mean stress value computed from the input tensor.
   !!

   ! q / \sqrt{3} = \sqrt{ J_{2} }

   ! From Potts and Zdravkoic -> J = \sqrt{J_{2}}

   function calc_mean_stress(stress) result(mean_stress)
      real(kind = dp), intent(in) :: stress(:)

      ! Local variables
      real(kind = dp) :: mean_stress, three = 3.0
      
      mean_stress = sum(stress(1:3)) / three

   end function calc_mean_stress

   function calc_deviatoric_stress(stress) result(deviatoric_stress)
      real(kind = dp), intent(in) :: stress(6)
      
      real(kind = dp) :: deviatoric_stress
      
      ! Local variables
      real(kind = dp) :: J2, three = 3.0

      ! Calc the J2 invariant
      J2 = calc_J2_invariant(stress)

      deviatoric_stress = sqrt( three * J2 )

   end function calc_deviatoric_stress

   ! Calculates the J2 stress invariant
   function calc_J2_invariant(stress) result(J2)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: J2
      
      ! Local variables
      real(kind = dp) :: mean_stress
      real(kind = dp) :: one_half = 0.5, two = 2.0
      real(kind = dp) :: stress_copy(6)
      
      ! Make a local copy of the stress vector
      stress_copy = stress

      mean_stress = calc_mean_stress(stress_copy)

      ! Generate the deviatoric_stress_tensor
      ! Need to operate on the diagonal of the stress tensor
      ! Which is the first 3 elements in voight notation

      ! s_{ij} = \sigma_{ij} - \sigma_{kk}/3 \delta_{ij}
      stress_copy(1:3) = stress_copy(1:3) - mean_stress
      
      ! Square all the terms
      stress_copy = stress_copy**2
      
      ! Multiply the shear terms by two since there's two sets of them
      stress_copy(4:6) = two * stress_copy(4:6)

      ! J2 = 1/2 s_{ij} s_{ji}
      J2 = one_half *  sum(stress_copy)
      
   end function calc_J2_invariant

   function calc_s_determinant(stress) result(s_det)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: s_det

      ! Local variables
      real(kind = dp) :: s(6), mean_stress
      real(kind = dp) :: first_term, second_term, third_term, fourth_term, &
                         fifth_term
      real(kind = dp), parameter :: one = 1.0, two = 2.0
      real(kind = dp) :: stress_copy(6) ! Used to make a local copy of the stress vector
      
      logical :: DEBUG = .True.


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

   function calc_J_invariant(stress) result(J)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: J
      ! Local variables
      real(kind = dp) :: J2

      ! Calc the J2 stress invariant
      J2 = calc_J2_invariant(stress)

      ! J = \sqrt{J2}
      J = sqrt(J2)
   end function calc_J_invariant

   function calc_lode_angle(stress) result(lode_angle)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: lode_angle

      ! Local variables
      real(kind = dp) :: J, det_s, inside
      real(kind = dp) :: zero = 0.0, one = 1.0, two = 2.0, three = 3.0

      ! Calc J
      J = calc_J_invariant(stress)
      
      if (J > zero) then
         ! Calc the determinant of s
         det_s = calc_s_determinant(stress)

         ! calc the inside of the paranthesis
         ! 3 \sqrt{3} / 2 * det(s) / J^{3}

         inside = (three * sqrt(three) / two) * det_s / J**3
         
         ! Calc the lode angle
         lode_angle = - one * asin(inside) / 3.0
      else
         lode_angle = 0.0
      end if
   end function

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

end module mod_stress_invariants
