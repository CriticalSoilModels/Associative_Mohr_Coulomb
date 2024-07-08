module ESM_Assoc_MohrCoulomb
    use kind_precision_module, only: dp
    use integer_precision_module, only: i32
    implicit none
    private
    public ESM_Assoc_MohrCoulomb

contains
subroutine ESM_Assoc_MohrCoulomb(npt,noel,idset,stress,eunloading, plasticmultiplier, dstran, nstatev,&
                                 statev,naddvar,additionalvar,cmname,nprops,props,numberofphases,ntens)

    ! Define variables that are being loaded in for historic puposes but won't be used
    integer(kind = i32), intent(in) :: npt, noel, idset, 
    real(kind = i32), intent(in) :: plasticmultiplier
    real(kind = dp) :: 

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




end subroutine ESM_Assoc_MohrCoulomb


! The purpose of this subroutine is to format the output of the Constitutive relation to fit the 
! expection of the ESM subroutine
subroutine UMAT_Assoc_MohrCoulomb(stress, ddsdde, dstran, statev, props)
    
    real(kind = dp), intent(inout) :: stress(6),  dstran(6), statev(:), props(:)
    real(kind = dp), intent(out)   :: ddsdde(6, 6)    ! Setting a variable to out means that it's value isn't stored in the repo

    call 
    
end subroutine UMAT_Assoc_MohrCoulomb
end module ESM_Assoc_MohrCoulomb




