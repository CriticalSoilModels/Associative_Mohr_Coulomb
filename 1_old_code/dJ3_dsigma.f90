module dJ3_dsigma_wrong
    ! Something is wrong with these functions. The calculated shear terms 
    ! of the calculated dJs_dsigma value is double what it should be 
    ! Leaving them here now in case I want to debug in the future

    ! WARNING - See the above message
    implicit none
   
contains
   ! Calc \frac{ \partial (det s) }{ \partial \sigma }
   ! Derivative of the determinant of the deviatoric stress tensor wrt.
   ! the stress tensor
   ! function calc_dJ3_dsigma(stress) result(dDetS_dsigma)
   !    real(kind = dp), intent(in) :: stress(6)
   !    real(kind = dp) :: dDetS_sigma(6)
    ! pure function calc_dJ3_dsigma(stress) result(dJ3_dsigma)
    !     real(kind = dp), intent(in) :: stress(6)
    !     real(kind = dp) :: dJ3_dsigma(6)

    !     ! Local variables
    !     real(kind = dp) :: ds_dsigma(6, 6), dJ3_ds(6)

    !     ! Calc dJ3/ds
    !     dJ3_ds = calc_J3_ds(stress)

    !     ! Calc ds/dsigma
    !     ds_dsigma  = get_ds_dsigma()

    !     dJ3_dsigma = matmul(ds_dsigma, dJ3_ds)
    ! end function calc_dJ3_dsigma

    ! ! Calc the derivative of the J3 deviatoric stress invariant wrt. the deviatoric stress (s)
    ! pure function calc_J3_ds(stress) result(dJ3_ds)
    !     real(kind = dp), intent(in) :: stress(6)
    !     real(kind = dp) :: dJ3_ds(6)

    !     ! Local variables
    !     real(kind = dp) :: s(6), s_squared(6), voigt_identity(6)
    !     real(kind = dp) :: trace_s_squared
        
    !     ! Calc the deviatoric matrix
    !     s = calc_s_voigt_vector(stress)

    !     ! Calc s^{2}
    !     s_squared = multiply_voigt_vectors(s, s)

    !     ! Calc the trace of s^{2}
    !     trace_s_squared = trace_voigt_vector( s_squared )

    !     ! Get the voigt identity matrix
    !     voigt_identity = get_3d_voigt_identity_vector()

    !     ! Get the derivative
    !     ! dJ3_ds = s^{2} - 1/2 tr( s^{2} ) \bar{ 1 } 
    !     dJ3_ds = s_squared - 0.5_dp * trace_s_squared * voigt_identity

    ! end function calc_J3_ds

    ! pure function get_ds_dsigma() result(ds_dsigma)
    !     real(kind = dp) :: ds_dsigma(6, 6)
        
    !     ! Local variables
    !     integer(kind = i32) :: i
    !     ! Store the value of the matrix
        
    !     ! Fill the upper block
    !     ds_dsigma(1:3, 1:3) = -1.0_dp/3.0_dp

    !     ! Fill the first part of the diagonal
    !     do i = 1, 3
    !         ds_dsigma(i,i) = 2.0_dp / 3.0_dp
    !     end do

    !     ! Fill the second part of the diagonal
    !     do i = 4, 6
    !         ds_dsigma(i, i) = 1.0_dp
    !     end do
    ! end function get_ds_dsigma

    ! pure function calc_s_voigt_vector(stress) result(s)
    !     real(kind=  dp), intent(in) :: stress(6)
    !     real(kind = dp) :: s(6)

    !     ! Local variables
    !     real(kind = dp) :: mean_stress

    !     ! Calc the mean stress
    !     mean_stress = calc_mean_stress(stress)

    !     ! Store all of the terms
    !     s(:) = stress(:)

    !     ! Subtract off the mean stress
    !     s(1:3) = s(1:3) - mean_stress

    ! end function calc_s_voigt_vector

end module dJ3_dsigma_wrong