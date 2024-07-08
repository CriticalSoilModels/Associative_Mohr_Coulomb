! Module contains functions and subroutines that are helpful with working with voigt form
! vectors
module mod_general_voigt
   use kind_precision_module, only: dp
   use integer_precision_module, only: i32

   implicit none

contains

   function voigt_2_matrix(voigt_vector) result(matrix)
      real(kind = dp), intent (in) :: voigt_vector(6)
      real(kind = dp) :: matrix(3,3)

      ! Local variables
      integer(kind = i32) :: i

      ! The order of the voigt vector is assumed to be:
      ! vector = { 11
      !            22
      !            33
      !            12
      !            23
      !            31}

      ! Fill the diagonal
      do i = 1, 3
         matrix(i,i) = voigt_vector(i)
      end do

      matrix(2, 3) = voigt_vector(5)
      matrix(3, 2) = voigt_vector(5)
      matrix(1, 3) = voigt_vector(6)
      matrix(3, 1) = voigt_vector(6)
      matrix(2, 1) = voigt_vector(4)
      matrix(1, 2) = voigt_vector(4)

   end function voigt_2_matrix

   ! TODO: For the time being I'm converting the vectors to matrices
   ! multiplying them then converting it back to a voigt vector
   ! In the future just write out the multiplication
   function multiply_voigt_vectors(a, b) result(a_times_b)
      ! This function multiplies the vectors as the matrices they actually represent
      ! Note: that this function assumes that the values of the vector are stored in
      ! the following componenet order
      ! vector = { 11
      !            22
      !            33
      !            12
      !            23
      !            31}

      real(kind = dp), intent(in) :: a(6), b(6)
      real(kind = dp) :: a_times_b(6), a_times_b_matrix(3, 3)

      ! Local variables
      real(kind= dp) :: a_matrix(3,3), b_matrix(3,3)

      ! Convert the voigt to a matrix
      a_matrix = voigt_2_matrix(a)
      b_matrix = voigt_2_matrix(b)

      ! Calc the matrix times matrix, since a and b are symmetric the result will also be symmetric
      a_times_b_matrix = matmul(a_matrix, b_matrix)

      ! Convert the resulting matrix to voigt notation
      a_times_b = symmetric_matrix_2_voigt(a_times_b_matrix)

   end function multiply_voigt_vectors

   function symmetric_matrix_2_voigt(symmetric_matrix) result(voigt_vector)
      real(kind = dp), intent(in) :: symmetric_matrix(3, 3)
      real(kind = dp) :: voigt_vector(6)
      ! Note: that this function assumes that the values of the vector are stored in
      ! the following componenet order
      ! vector = { 11
      !            22
      !            33
      !            12
      !            23
      !            31}

      ! Local variables
      integer(kind = i32) :: i

      ! Store the diagonal
      do i = 1, 3
        voigt_vector(i) = symmetric_matrix(i, i)
      end do

      ! Set the off diagonal terms
      voigt_vector(4) = symmetric_matrix(1, 2)
      voigt_vector(5) = symmetric_matrix(2, 3)
      voigt_vector(6) = symmetric_matrix(3, 1)

   end function symmetric_matrix_2_voigt

   function get_voigt_identity_vector(vector_size) result(voigt_identity_vector)
      integer(kind = i32), intent(in) :: vector_size
      real(kind = dp), allocatable :: voigt_identity_vector(:)

      ! Local variables
      integer(kind = i32), parameter :: voigt_2_dim = 4
      integer(kind = i32), parameter:: voigt_3_dim = 6

      select case(vector_size)
       case(voigt_2_dim)
         ! Allocate the vector
         allocate(voigt_identity_vector(voigt_2_dim))

         ! Initialize to zero
         voigt_identity_vector = 0.0_dp

         ! Set the values (only first 2 elements are 1 for 2D case)
         voigt_identity_vector(1:2) = 1.0_dp

       case(voigt_3_dim)
         ! Allocate the vector
         allocate(voigt_identity_vector(voigt_3_dim))

         ! Initialize to zero
         voigt_identity_vector = 0.0_dp

         ! Set the values (only first 3 elements are 1 for 3D case)
         voigt_identity_vector(1:3) = 1.0_dp

       case default
         ! Don't allocate the vector and print an error message
         print *, "A vector size of: ", vector_size, "is not a valid input"
         print *, "Valid inputs are: ", voigt_2_dim, ",", voigt_3_dim
         ! Return an unallocated array in this case
      end select

   end function get_voigt_identity_vector

   function trace_voigt_vector(voigt_vector) result(trace)
        real(kind = dp) :: voigt_vector(6)
        real(kind = dp) :: trace

        ! Calc the trace of the voight vector
        ! Sum of the diagonal terms
        trace = sum( voigt_vector(1:3) )
        
   end function trace_voigt_vector

end module mod_general_voigt
