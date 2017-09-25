!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module provides a smoother type for smoothing a system of equations.
!>        It is intended to support multigrid linear smoothers.
!>
!> @par Module Dependencies
!>  - @ref IntrType "IntrType": @copybrief IntrType
!>  - @ref BLAS "BLAS": @copybrief BLAS
!>  - @ref Times "Times": @copybrief Times
!>  - @ref ExceptionHandler "ExceptionHandler": @copybrief ExceptionHandler
!>  - @ref Allocs "Allocs": @copybrief Allocs
!>  - @ref ParameterLists "ParameterLists": @copybrief ParameterLists
!>  - @ref ParallelEnv "ParallelEnv": @copybrief ParallelEnv
!>  - @ref VectorTypes "VectorTypes": @copybrief VectorTypes
!>  - @ref MatrixTypes "MatrixTypes": @copybrief MatrixTypes
!>
!> @par EXAMPLES
!> @code
!>
!> @endcode
!>
!> @author Ben C. Yee
!>   @date 09/21/2017
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE SmootherTypes
  USE IntrType
  USE BLAS
  USE trilinos_interfaces
  USE Times
  USE ExceptionHandler
  USE Allocs
  USE ParameterLists
  USE ParallelEnv
  USE VectorTypes
  USE MatrixTypes
  USE PreconditionerTypes
  USE Strings
  USE IOUtil
  IMPLICIT NONE

#ifdef FUTILITY_HAVE_PETSC
#include <finclude/petsc.h>
!petscisdef.h defines the keyword IS, and it needs to be reset
#undef IS
#endif

  PRIVATE
!
! List of public members
  PUBLIC :: eSmootherType
  PUBLIC :: SmootherType_Base
  PUBLIC :: SmootherType_PETSc
  PUBLIC :: SmootherType_PETSc_CBJ
  PUBLIC :: IndexList

  !> Enumeration for smoother options
  INTEGER(SIK),PARAMETER,PUBLIC :: CBJ=0
  !> Enumeration for block solver options
  INTEGER(SIK),PARAMETER,PUBLIC :: LU=0,SOR=1,ILU=2
  !> set enumeration scheme for TPLs
  INTEGER(SIK),PARAMETER,PUBLIC :: PETSC=0,NATIVE=4

  !> @brief the base linear smoother type
  TYPE,ABSTRACT :: SmootherType_Base
    !> Initialization status
    LOGICAL(SBK) :: isInit=.FALSE.
    !> Integer flag for the solution methodology desired
    INTEGER(SIK) :: smootherMethod=-1
    !> Integer flag for the solution methodology desired
    INTEGER(SIK) :: blockMethod=-1
    !> Integer flag for the solution methodology desired
    INTEGER(SIK) :: TPLType=-1
    !> Pointer to an MPI parallel environment
    TYPE(MPI_EnvType) :: MPIparallelEnv
    !> Has initial guess?
    LOGICAL(SBK) :: hasX0=.FALSE.
    !> Current local residual norm
    REAL(SRK) :: localResidual=0.0_SRK
    !> Current global residual norm
    REAL(SRK) :: globalResidual=0.0_SRK
    !> Starting local index:
    INTEGER(SIK) :: istt=-1_SIK
    !> End local index:
    INTEGER(SIK) :: istp=-1_SIK
    !> Block size (number of unknowns per point):
    INTEGER(SIK) :: blk_size=-1_SIK

  !
  !List of Type Bound Procedures
    CONTAINS
      !> Deferred routine for clearing the linear smoother
      PROCEDURE(smoother_sub_absintfc),DEFERRED,PASS :: clear
      !> Deferred routine for applying a smoother step to the linear system
      PROCEDURE(smoother_sub_absintfc),DEFERRED,PASS :: smooth
  ENDTYPE SmootherType_Base

  TYPE,ABSTRACT,EXTENDS(SmootherType_Base) :: SmootherType_PETSc
#ifdef FUTILITY_HAVE_PETSC
    !> Pointer to PETSc KSP object, should type KSPRICHARDSON
    KSP :: ksp
    !> Pointer to PETSc PC object corresponding to ksp, should type PCSHELL
    PC :: pc
#else
    !> Dummy attribute to make sure Futility compiles when PETSc is not present
    INTEGER(SIK) :: ksp=-1_SIK
    !> Dummy attribute to make sure Futility compiles when PETSc is not present
    INTEGER(SIK) :: pc=-1_SIK
#endif
  ENDTYPE SmootherType_PETSc

  !Handy structure to store list of indices (which can vary in length from
  !  color to color)
  TYPE :: IndexList
    !Number of indices for this color:
    INTEGER(SIK) :: num_indices=-1_SIK
    !List of indices for this color:
    INTEGER(SIK),ALLOCATABLE :: index_list(:)
  ENDTYPE IndexList

  !Colored block Jacobi smoother:
  TYPE,EXTENDS(SmootherType_PETSc) :: SmootherType_PETSc_CBJ
    !> Number of colors:
    INTEGER(SIK) :: num_colors=-1_SIK
    !> List of indices for each color (1:num_colors)
    TYPE(IndexList),ALLOCATABLE :: colors(:)
    !> Whether or not each of the color index lists have been set (1:num_colors)
    LOGICAL(SBK),ALLOCATABLE :: hasColorDefined(:)
    !> Whether or not all the color index lists have been set
    LOGICAL(SBK) :: hasAllColorsDefined=.FALSE.
    !> Color ID of each point (istt:istp)
    INTEGER(SIK),ALLOCATABLE :: color_ids(:)

  !
  !List of Type Bound Procedures
    CONTAINS
      !> @copybrief SmootherTypes::init_SmootherType_PETSc_CBJ
      !> @copydetails SmootherTypes::init_SmootherType_PETSc_CBJ
      PROCEDURE,PASS :: init => init_SmootherType_PETSc_CBJ
      !> @copybrief SmootherTypes::clear_SmootherType_PETSc_CBJ
      !> @copydetails SmootherTypes::clear_SmootherType_PETSc_CBJ
      PROCEDURE,PASS :: clear => clear_SmootherType_PETSc_CBJ
      !> @copybrief SmootherTypes::smooth_SmootherType_PETSc_CBJ
      !> @copydetails SmootherTypes::smooth_SmootherType_PETSc_CBJ
      PROCEDURE,PASS :: smooth => smooth_SmootherType_PETSc_CBJ
      !> @copybrief SmootherTypes::defineColor_SmootherType_PETSc_CBJ
      !> @copydetails SmootherTypes::defineColor_SmootherType_PETSc_CBJ
      PROCEDURE,PASS :: defineColor => defineColor_SmootherType_PETSc_CBJ
      !> @copybrief SmootherTypes::defineAllColors_SmootherType_PETSc_CBJ
      !> @copydetails SmootherTypes::defineAllColors_SmootherType_PETSc_CBJ
      PROCEDURE,PASS :: defineAllColors => defineAllColors_SmootherType_PETSc_CBJ
  ENDTYPE SmootherType_PETSc_CBJ

  !> Exception Handler for use in SmootherTypes
  TYPE(ExceptionHandlerType),SAVE :: eSmootherType

  !> Explicitly defines the interface for the clear, and solve routines
  ABSTRACT INTERFACE
    SUBROUTINE smoother_sub_absintfc(smoother)
      IMPORT :: SmootherType_Base
      CLASS(SmootherType_Base),INTENT(INOUT) :: smoother
    ENDSUBROUTINE smoother_sub_absintfc
  ENDINTERFACE

  !> Explicitly defines the interface for the init routine for PETSc smoothers
  ABSTRACT INTERFACE
    SUBROUTINE smoother_petsc_sub_absintfc(smoother,ksp,params)
      IMPORT :: SmootherType_PETSc
      IMPORT :: ParamType
      CLASS(SmootherType_PETSc),INTENT(INOUT) :: smoother
      TYPE(ParamType),INTENT(IN) :: params
#ifdef FUTILITY_HAVE_PETSC
      KSP,INTENT(INOUT) :: ksp
#else
      INTEGER(SIK),INTENT(INOUT) :: ksp
#endif
    ENDSUBROUTINE smoother_petsc_sub_absintfc
  ENDINTERFACE

  !> Name of module
  CHARACTER(LEN=*),PARAMETER :: modName='SMOOTHERTYPES'
!
!------------------------------------------------------------------------------
  CONTAINS
!
!-------------------------------------------------------------------------------
!> @brief Initializes the SmootherType for an CBJ PETSc smoother
!>
!> @param smoother The smoother object to act on
!> @param ksp The PETSc ksp object to act on
!> @param params Parameter list, which must contain num_colors, istt, istp,
!>        blk_size, and MPI_Comm_ID
!>
    SUBROUTINE init_SmootherType_PETSc_CBJ(smoother,ksp,params)
      CHARACTER(LEN=*),PARAMETER :: myName='init_SmootherType_PETSc_CBJ'
      CLASS(SmootherType_PETSc_CBJ),INTENT(INOUT) :: smoother
      TYPE(ParamType),INTENT(IN) :: params
      LOGICAL(SBK) :: tmpbool
      INTEGER(SIK) :: MPI_Comm_ID
#ifdef FUTILITY_HAVE_PETSC
      KSP,INTENT(INOUT) :: ksp
      PC :: pc
      PetscErrorCode :: iperr
#else
      INTEGER(SIK),INTENT(INOUT) :: ksp
      INTEGER(SIK) :: pc

      CALL eSmootherType%raiseError(modName//"::"//myName//" - "// &
          "This type should only be used with PETSc enabled!")
#endif
      !Check parameter list:
      tmpbool=params%has('SmootherType->num_colors') .AND. &
              params%has('SmootherType->istt') .AND. &
              params%has('SmootherType->istp') .AND. &
              params%has('SmootherType->blk_size') .AND. &
              params%has('SmootherType->MPI_Comm_ID')
      IF(.NOT. tmpbool) &
        CALL eSmootherType%raiseError(modName//"::"//myName//" - "// &
            "Missing a parameter from the parameter list!")

      !Extract param list info:
      CALL params%get('SmootherType->istt',smoother%istt)
      CALL params%get('SmootherType->istp',smoother%istp)
      CALL params%get('SmootherType->num_colors',smoother%num_colors)
      CALL params%get('SmootherType->blk_size',smoother%blk_size)
      CALL params%get('SmootherType->MPI_Comm_ID',MPI_Comm_ID)

      !Create the color index lists:
      ALLOCATE(smoother%colors(smoother%num_colors))
      ALLOCATE(smoother%hasColorDefined(smoother%num_colors))
      smoother%hasColorDefined=.FALSE.
      smoother%hasAllColorsDefined=.FALSE.
      !Create the color ID list:
      ALLOCATE(smoother%color_ids(smoother%istt:smoother%istp))

      MPI_Comm_ID=-1_SIK
      IF(MPI_Comm_ID /= -1_SIK) THEN
        CALL smoother%MPIparallelEnv%init(MPI_Comm_ID)
      ELSE
        CALL smoother%MPIparallelEnv%init(MPI_COMM_WORLD)
      ENDIF

      smoother%smootherMethod=CBJ
      smoother%TPLType=PETSc
      smoother%blockMethod=LU
      IF(params%has('SmootherType->blockMethod')) &
        CALL params%get('SmootherType->blockMethod', &
                          smoother%blockMethod)
      smoother%hasX0=.FALSE.

      smoother%ksp=ksp
#ifdef FUTILITY_HAVE_PETSC
      CALL KSPSetType(smoother%ksp,KSPRICHARDSON,iperr)
      CALL KSPGetPC(smoother%ksp,smoother%pc,iperr)
      CALL PCSetType(smoother%pc,PCSHELL,iperr)
#endif
      smoother%isInit=.TRUE.

    ENDSUBROUTINE init_SmootherType_PETSc_CBJ
!
!-------------------------------------------------------------------------------
!> @brief Performs a sweep over the red cells, and then the black cells
!> @param smoother The smoother object to act on
!>
    SUBROUTINE smooth_SmootherType_PETSc_CBJ(smoother)
      CHARACTER(LEN=*),PARAMETER :: myName='smooth_SmootherType_PETSc_CBJ'
      CLASS(SmootherType_PETSc_CBJ),INTENT(INOUT) :: smoother
    ENDSUBROUTINE smooth_SmootherType_PETSc_CBJ
!
!-------------------------------------------------------------------------------
!> @brief Fill out an index list for a particular color
!>
!> @param smoother The smoother object to act on
!> @param icolor Color being defined
!> @param index_list list of indices for color icolor
!>
    SUBROUTINE defineColor_SmootherType_PETSc_CBJ(smoother,icolor,index_list)
      CHARACTER(LEN=*),PARAMETER :: myName='defineColor_SmootherType_PETSc_CBJ'
      CLASS(SmootherType_PETSc_CBJ),INTENT(INOUT) :: smoother
      INTEGER(SIK),INTENT(IN) :: icolor
      INTEGER(SIK),INTENT(IN) :: index_list(:)

      INTEGER(SIK) :: i,num_indices

      IF(.NOT. smoother%isInit) &
        CALL eSmootherType%raiseError(modName//"::"//myName//" - "// &
            "Smoother must be initialized first!")

      num_indices=SIZE(index_list)
      smoother%colors(icolor)%num_indices=num_indices
      ALLOCATE(smoother%colors(icolor)%index_list(num_indices))
      smoother%colors(icolor)%index_list=index_list
      DO i=1,num_indices
         smoother%color_ids(index_list(i))=icolor
      ENDDO
      smoother%hasColorDefined(icolor)=.TRUE.
      smoother%hasAllColorsDefined=ALL(smoother%hasColorDefined)

    ENDSUBROUTINE defineColor_SmootherType_PETSc_CBJ
!
!-------------------------------------------------------------------------------
!> @brief Provide the color_ids to fill out the index_lists for all the colors
!>
!> @param smoother The smoother object to act on
!> @param color_ids List of colors for each point, must have bounds
!>         solver%istt:solver%istp
!>
    SUBROUTINE defineAllColors_SmootherType_PETSc_CBJ(smoother,color_ids)
      CHARACTER(LEN=*),PARAMETER :: myName='defineAllColors_SmootherType_PETSc_CBJ'
      CLASS(SmootherType_PETSc_CBJ),INTENT(INOUT) :: smoother
      INTEGER(SIK),INTENT(IN) :: color_ids(:)

      INTEGER(SIK) :: icolor,i
      INTEGER(SIK),ALLOCATABLE :: tmpints(:)

      IF(.NOT. smoother%isInit) &
        CALL eSmootherType%raiseError(modName//"::"//myName//" - "// &
            "Smoother must be initialized first!")

      IF(LBOUND(color_ids,DIM=1) /= smoother%istt .OR. UBOUND(color_ids,DIM=1) /= &
          smoother%istp) CALL eSmootherType%raiseError(modName//"::"//myName &
            //" - wrong colors!")

      smoother%color_ids=color_ids

      DO icolor=1,smoother%num_colors
        smoother%colors(icolor)%num_indices=0
      ENDDO
      DO i=smoother%istt,smoother%istp
        icolor=color_ids(i)
        smoother%colors(icolor)%num_indices= &
           smoother%colors(icolor)%num_indices+1
      ENDDO
      DO icolor=1,smoother%num_colors
        ALLOCATE(smoother%colors(icolor)%index_list( &
                                       smoother%istp-smoother%istt+1))
      ENDDO
      ALLOCATE(tmpints(smoother%num_colors))
      tmpints=0_SIK
      DO i=smoother%istt,smoother%istp
        icolor=color_ids(i)
        tmpints(icolor)=tmpints(icolor)+1
        smoother%colors(icolor)%index_list(tmpints(icolor))=i
      ENDDO
      DEALLOCATE(tmpints)

    ENDSUBROUTINE defineAllColors_SmootherType_PETSc_CBJ
!
!-------------------------------------------------------------------------------
!> @brief Clears the SmootherType for an CBJ PETSc smoother
!> @param smoother The smoother object to act on
!>
    SUBROUTINE clear_SmootherType_PETSc_CBJ(smoother)
      CHARACTER(LEN=*),PARAMETER :: myName='clear_SmootherType_PETSc_CBJ'
      CLASS(SmootherType_PETSc_CBJ),INTENT(INOUT) :: smoother
      INTEGER(SIK) :: icolor

      IF(ALLOCATED(smoother%colors)) THEN
        DO icolor=1,smoother%num_colors
          IF(smoother%hasColorDefined(icolor)) THEN
            DEALLOCATE(smoother%colors(icolor)%index_list)
            smoother%hasColorDefined(icolor)=.FALSE.
          ENDIF
        ENDDO
        DEALLOCATE(smoother%colors)
        DEALLOCATE(smoother%hasColorDefined)
        smoother%hasAllColorsDefined=.FALSE.
      ENDIF
      IF(ALLOCATED(smoother%color_ids)) DEALLOCATE(smoother%color_ids)
      CALL smoother%MPIparallelEnv%clear()
      smoother%isInit=.FALSE.

    ENDSUBROUTINE clear_SmootherType_PETSc_CBJ
!
ENDMODULE SmootherTypes
