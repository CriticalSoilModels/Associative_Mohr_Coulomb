module global_variables
    use kind_precision_module, only: dp

    ! Global variables - be very selective with what you include in here
    implicit none

    real(kind = dp), parameter :: PI = 3.1415927
    
end module global_variables