module general_stress_equations
   implicit none

contains

   function calc_mean_stress(stress) result(mean_stress)
      real, intent(in) :: stress(:)
      real :: mean_stress

      ! Local variables
      integer :: i

      mean_stress = 0.0
      do i = 1, 3
         mean_stress = mean_stress + stress(i)
      end do

      mean_stress = mean_stress / 3.0
   end function calc_mean_stress

   subroutine calc_stress_invariants(Sig, p, q, theta)
      !*********************************************************************
      ! Takes the stress tensor Sig and return invariants p, q, and theta  *
      !																	 *
      !*********************************************************************
      !input variables
      real, dimension(6), intent(in):: Sig
      !output variables
      real, intent(out)::p, q, theta
      !local variables
      real:: dev(6), J2, J3, sin3theta

      p=(Sig(1)+Sig(2)+Sig(3))/3.0 !mean stress
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

end module general_stress_equations
