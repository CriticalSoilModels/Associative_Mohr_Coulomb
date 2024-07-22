module mod_ESM_assoc_MC
    use kind_precision_module, only: dp, i32
    use mod_UMAT_assoc_MC, only: UMAT_assoc_MC

    implicit none
    private
    public ESM_assoc_MC

contains
!------- Variables and what they mean -------

    !------- Variables that I don't need -------
    ! npt: Number of integration points
    ! noel: Number of elements
    ! idset: Identifier for the set of elements
    ! plasticmultiplier: I think this is lambda
    ! nstatev: Number of state variables
    ! naddvar: Number of additional variables
    ! additionalvar: Additionalvariables
    ! cmname: Constitutive model name
    ! nprops: Number of material properties
    ! numberofphases: Number of phases in the material
    ! ntens: Number of tensor components
    
    ! ------- Variables that I need -------
    ! stress: Stress tensor
    ! eunloading: ! Constrianed Modulus that is used to calculate the wave speed for the time step
    ! dstran: strain increment
    ! statev: State variables
    ! props: Material properties

    subroutine ESM_assoc_MC(npt,noel,idset,stress,eunloading, plasticmultiplier, dstran, nstatev,&
                                    statev,naddvar,additionalvar,cmname,nprops,props,numberofphases,ntens)

        ! Define variables that are being loaded in for historic puposes but won't be used
        integer(kind = i32), intent(in)    :: npt, noel, idset, numberofphases, ntens,  naddvar
        integer(kind = i32), intent(in)    :: nstatev, nprops
        real(kind = dp)    , intent(in)    :: plasticmultiplier, additionalvar(:)
        character(len =*)  , intent(in)    :: cmname
        

        ! Define variables that i need
        real(kind = dp)    , intent(inout) :: stress(6), eunloading, dstran(6), &
                                            statev(:), props(:)

        ! Local variables
        real(kind = dp) :: stiff_matrix(6, 6)
        integer(kind = i32), parameter :: max_iterations = 10000
        real(kind = dp), allocatable :: use_additionalvar(:)

        allocate(use_additionalvar(size(additionalvar)))

        ! Using the input variables so it doesn't complain

        if (.False.) then
            print *, npt, noel, idset, nstatev, naddvar, nprops, numberofphases, ntens, &
                    plasticmultiplier, cmname
        end if

        ! Zero the matrix
        stiff_matrix(:, :) = 0.0_dp

        ! Expects statev:
            ! statev(1) = cohesion
            ! statev(2) = friction angle
        ! Props:
            ! props(1) = shear modulus
            ! props(2) = poisson's ratio
        call UMAT_assoc_MC(stress, dstran, statev, props, stiff_matrix, max_iterations)

        ! Calc eunloading
        eunloading = max(stiff_matrix(1,1), stiff_matrix(2,2), stiff_matrix(3,3))

    end subroutine ESM_assoc_MC

 end module mod_ESM_assoc_MC




