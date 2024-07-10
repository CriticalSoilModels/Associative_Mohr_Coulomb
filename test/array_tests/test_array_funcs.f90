program test_array_funcs
    use kind_precision_module, only: dp, i32

    use mod_array_helper, only: reorder_real_array

    implicit none
    
    real(kind = dp) :: input_arr(6), output_arr(6), expected_output(6)
    integer(kind = i32) :: new_order(6)

    ! Testing that the array reordering scheme is working
    input_arr = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]

    new_order = [5, 6, 4, 3, 2, 1]
    
    output_arr = reorder_real_array(input_arr, new_order)

    ! TODO: Add a check for this test

    expected_output = [5.0_dp, 6.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]


    if ( arrays_equal(output_arr, expected_output) ) then
        print *, "Algorithm to reorder an array using an input array is working"
    else
        print *, "reorder_real_array is not working as expected"
        print *, "Expected array is: ", expected_output
        print *, "Actual output is : ", output_arr
    end if

contains

    ! Check if the values of a real array are the same as another array
    logical function arrays_equal(arr1, arr2)
        real(kind = dp), intent(in) :: arr1(:), arr2(:)

        ! Check if the sizes are the same
        if (size(arr1) /= size(arr2)) then
            arrays_equal = .false.
            return
        end if

        ! Check if all elements are the same
        arrays_equal = all(arr1 == arr2)
    end function arrays_equal

end program test_array_funcs