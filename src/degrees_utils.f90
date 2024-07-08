module degrees_utils
    implicit none
    
    private
    public :: convert_deg_2_radians

contains
    pure function convert_deg_2_radians(degrees) result(radians)
        use kind_precision_module, only: dp
        use global_variables, only: PI

        real(kind = dp), intent(in) :: degrees
        real(kind = dp) :: radians

        radians = degrees * PI / 180.0_dp 
    end function
    
end module degrees_utils