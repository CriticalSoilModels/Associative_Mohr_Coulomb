module integer_precision_module
  use, intrinsic :: iso_fortran_env, only: i8=>int8, i16=>int16, i32=>int32, i64=>int64
  implicit none

  private
  public :: i8, i16, i32, i64

end module integer_precision_module