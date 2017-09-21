!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module provides a smoother type for smoothing a system of equations.
!>        It is intended to support multigrid linear solvers.
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

  PRIVATE
!
! List of public members
  PUBLIC :: eSmootherType
  PUBLIC :: SmootherType_Base

  !> set enumeration scheme for TPLs
  INTEGER(SIK),PARAMETER,PUBLIC :: PETSC=0,NATIVE=4
  !> Enumeration for solver options
  INTEGER(SIK),PARAMETER,PUBLIC :: SOR=0


  !> @brief the base linear solver type
  TYPE,ABSTRACT :: SmootherType_Base
    !> Initialization status
    LOGICAL(SBK) :: isInit=.FALSE.
    !> Integer flag for the solution methodology desired
    INTEGER(SIK) :: smootherMethod=-1
    !> Integer flag for the solution methodology desired
    INTEGER(SIK) :: TPLType=-1
    !> Pointer to an MPI parallel environment
    CLASS(MPI_EnvType),POINTER :: MPIparallelEnv
    !> Has initial guess?
    LOGICAL(SBK) :: hasX0=.FALSE.
    !> Current local residual norm
    REAL(SRK) :: localResidual=0.0_SRK
    !> Current global residual norm
    REAL(SRK) :: globalResidual=0.0_SRK

  !
  !List of Type Bound Procedures
    CONTAINS
      !> Deferred routine for initializing the linear solver system
      PROCEDURE(smoother_sub_absintfc),DEFERRED,PASS :: init
      !> Deferred routine for clearing the linear solver
      PROCEDURE(smoother_sub_absintfc),DEFERRED,PASS :: clear
      !> Deferred routine for applying a smoother step to the linear system
      PROCEDURE(smoother_sub_absintfc),DEFERRED,PASS :: smooth
  ENDTYPE SmootherType_Base

  !> Exception Handler for use in SmootherTypes
  TYPE(ExceptionHandlerType),SAVE :: eSmootherType

  !> Explicitly defines the interface for the clear and solve routines
  ABSTRACT INTERFACE
    SUBROUTINE smoother_sub_absintfc(solver)
      IMPORT :: SmootherType_Base
      CLASS(SmootherType_Base),INTENT(INOUT) :: solver
    ENDSUBROUTINE smoother_sub_absintfc
  ENDINTERFACE

!
ENDMODULE SmootherTypes
