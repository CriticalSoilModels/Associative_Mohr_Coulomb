module mod_Inc_Driver_ESM
    use kind_precision_module, only: dp, i32
    use mod_UMAT_assoc_MC, only: UMAT_assoc_MC
    implicit none
    
contains
    subroutine Inc_Driver_ESM_Wrapper(stress,statev,ddsdde,sse,spd,scd, &
        rpl,ddsddt,drplde,drpldt, &
        stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
        ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt, &
        celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)


    ! parameters
    character(len=80), intent(in) :: cmname
    integer, intent(in) :: ndi, nshr, ntens, nstatev, nprops, noel, npt, &
    layer, kspt, kstep, kinc
    real(dp), intent(in) :: dtime, temp, dtemp, celent, pnewdt

    ! arrays
    real(dp), intent(inout) :: stress(ntens), statev(nstatev), &
    ddsdde(ntens,ntens)

    real(dp), intent(in) :: ddsddt(ntens), drplde(ntens)

    real(dp), intent(in) :: stran(ntens), dstran(ntens), time(2), &
    predef(1), dpred(1), props(nprops), &
    coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)

    ! Scalars
    real(dp), intent(in) :: SSE, SPD, SCD, RPL, DRPLDT

    ! Local variables
    integer(kind = i32), parameter :: max_iterations = 10000

    ! Using the variables so the error goes away
    if (.False.) then
    print *, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, time, dtime, temp, &
    dtemp, predef, dpred, cmname, ndi, nshr, coords, drot, pnewdt, celent, &
    dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc
    end if

    call UMAT_assoc_MC(stress, dstran, statev, props, ddsdde, max_iterations)

    end subroutine Inc_Driver_ESM_Wrapper
end module mod_Inc_Driver_ESM