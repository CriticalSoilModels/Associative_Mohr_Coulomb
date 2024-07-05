module kind_precision_module
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64, qp=>real128

    implicit none
  
    private
    public :: sp, dp, qp
  
end module kind_precision_module