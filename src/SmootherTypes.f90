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
  PUBLIC :: SmootherType_PETSc_RBBJ

  !> set enumeration scheme for TPLs
  INTEGER(SIK),PARAMETER,PUBLIC :: PETSC=0,NATIVE=4
  !> Enumeration for smoother options
  INTEGER(SIK),PARAMETER,PUBLIC :: SOR=0

  !> @brief the base linear smoother type
  TYPE,ABSTRACT :: SmootherType_Base
    !> Initialization status
    LOGICAL(SBK) :: isInit=.FALSE.
    !> Integer flag for the solution methodology desired
    INTEGER(SIK) :: smootherMethod=-1
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

  !
  !List of Type Bound Procedures
    CONTAINS
      !> Deferred routine for initializing the linear smoother system
      PROCEDURE(smoother_petsc_sub_absintfc),DEFERRED,PASS :: init
  ENDTYPE SmootherType_PETSc

  !Red black block Jacobi smoother:
  TYPE,EXTENDS(SmootherType_PETSc) :: SmootherType_PETSc_RBBJ
    !> Number of local red blocks: (swept over first)
    INTEGER(SIK) :: num_red=-1_SIK
    !> Number of local black blocks: (swept over second)
    INTEGER(SIK) :: num_black=-1_SIK
    !> Starting local block index:
    INTEGER(SIK) :: istt=-1_SIK
    !> End local block index:
    INTEGER(SIK) :: istp=-1_SIK
    !Note that we should have num_red+num_black=istp-istt+1
    !> Logical array indicating whether an index is red (istt:istp)
    LOGICAL(SBK),ALLOCATABLE :: is_red(:)
    !> Block size:
    !>  The matrix rows owned locally should be 
    !>    (istt-1)*blk_size+1:istp*blk_size
    INTEGER(SIK) :: blk_size=-1_SIK

  !
  !List of Type Bound Procedures
    CONTAINS
      !> @copybrief SmootherTypes::init_SmootherType_PETSc_RBBJ
      !> @copydetails SmootherTypes::init_SmootherType_PETSc_RBBJ
      PROCEDURE,PASS :: init => init_SmootherType_PETSc_RBBJ
      !> @copybrief SmootherTypes::clear_SmootherType_PETSc_RBBJ
      !> @copydetails SmootherTypes::clear_SmootherType_PETSc_RBBJ
      PROCEDURE,PASS :: clear => clear_SmootherType_PETSc_RBBJ
      !> @copybrief SmootherTypes::smooth_SmootherType_PETSc_RBBJ
      !> @copydetails SmootherTypes::smooth_SmootherType_PETSc_RBBJ
      PROCEDURE,PASS :: smooth => smooth_SmootherType_PETSc_RBBJ
  ENDTYPE SmootherType_PETSc_RBBJ

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
!> @brief Initializes the SmootherType for an RBBJ PETSc smoother
!> @param smoother The smoother object to act on
!>
    SUBROUTINE init_SmootherType_PETSc_RBBJ(smoother,ksp,params)
      CHARACTER(LEN=*),PARAMETER :: myName='init_SmootherType_PETSc_RBBJ'
      CLASS(SmootherType_PETSc_RBBJ),INTENT(INOUT) :: smoother
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

      IF(.NOT. tmpbool) &
        CALL eSmootherType%raiseError(modName//"::"//myName//" - "// &
            "This type should only be used with PETSc enabled!")
#endif
      !Check parameter list:
      tmpbool=params%has('SmootherType->num_red') .AND. &
              params%has('SmootherType->num_black') .AND. &
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
      ALLOCATE(smoother%is_red(smoother%istt:smoother%istp))
      CALL params%get('SmootherType->is_red',smoother%is_red)
      CALL params%get('SmootherType->blk_size',smoother%blk_size)
      CALL params%get('SmootherType->MPI_Comm_ID',MPI_Comm_ID)

      !Count the # of red/black cells:
      smoother%num_red=COUNT(smoother%is_red)
      smoother%num_black=smoother%istp-smoother%istt-smoother%num_red+1

      MPI_Comm_ID=-1_SIK
      IF(MPI_Comm_ID /= -1_SIK) THEN
        CALL smoother%MPIparallelEnv%init(MPI_Comm_ID)
      ELSE
        CALL smoother%MPIparallelEnv%init(MPI_COMM_WORLD)
      ENDIF

      smoother%ksp=ksp
#ifdef FUTILITY_HAVE_PETSC
      CALL KSPSetType(smoother%ksp,KSPRICHARDSON,iperr)
      CALL KSPGetPC(smoother%ksp,smoother%pc,iperr)
      CALL PCSetType(smoother%pc,PCSHELL,iperr)
#endif
      smoother%isInit=.TRUE.

    ENDSUBROUTINE init_SmootherType_PETSc_RBBJ
!
!-------------------------------------------------------------------------------
!> @brief Performs a sweep over the red cells, and then the black cells
!> @param smoother The smoother object to act on
!>
    SUBROUTINE smooth_SmootherType_PETSc_RBBJ(smoother)
      CHARACTER(LEN=*),PARAMETER :: myName='smooth_SmootherType_PETSc_RBBJ'
      CLASS(SmootherType_PETSc_RBBJ),INTENT(INOUT) :: smoother
    ENDSUBROUTINE smooth_SmootherType_PETSc_RBBJ

!
!-------------------------------------------------------------------------------
!> @brief Clears the SmootehrType for an RBBJ PETSc smoother
!> @param smoother The smoother object to act on
!>
    SUBROUTINE clear_SmootherType_PETSc_RBBJ(smoother)
      CHARACTER(LEN=*),PARAMETER :: myName='clear_SmootherType_PETSc_RBBJ'
      CLASS(SmootherType_PETSc_RBBJ),INTENT(INOUT) :: smoother

      IF(ALLOCATED(smoother%is_red)) DEALLOCATE(smoother%is_red)
      CALL smoother%MPIparallelEnv%clear()
      smoother%isInit=.FALSE.

    ENDSUBROUTINE clear_SmootherType_PETSc_RBBJ
!
ENDMODULE SmootherTypes
