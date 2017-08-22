!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PROGRAM testLinearSolver_Multigrid
#include "UnitTest.h"
  USE UnitTest
  USE IntrType
  USE ExceptionHandler
  USE ParameterLists
  USE ParallelEnv
  USE VectorTypes
  USE MatrixTypes
  USE PreconditionerTypes
  USE LinearSolverTypes
  USE LinearSolverTypes_Multigrid

  IMPLICIT NONE

  TYPE(ExceptionHandlerType),TARGET :: e
  TYPE(MPI_EnvType) :: mpiTestEnv
  TYPE(ParamType) :: pList, optListLS, optListMat, vecPList

#ifdef FUTILITY_HAVE_PETSC
#include <finclude/petscdef.h>
#undef IS
  PetscErrorCode  :: ierr

  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#endif

  !> set up default parameter list
  CALL optListLS%clear()
  CALL optListLS%add('LinearSolverType->TPLType',NATIVE)
  CALL optListLS%add('LinearSolverType->solverMethod',1_SNK) ! GE or BICGSTAB
  CALL optListLS%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
  CALL optListLS%add('LinearSolverType->numberOMP',1_SNK)
  CALL optListLS%add('LinearSolverType->timerName','LinearSolver Timer')
  CALL optListLS%add('LinearSolverType->matType',SPARSE)
  ! set parameters for matrices
  CALL optListLS%add('LinearSolverType->A->MatrixType->n',-1_SNK)
  CALL optListLS%add('LinearSolverType->A->MatrixType->nnz',-1_SNK)
  CALL optListLS%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
  ! set parameters for vectors
  CALL optListLS%add('LinearSolverType->x->VectorType->n',-1_SNK)
  CALL optListLS%add('LinearSolverType->b->VectorType->n',-1_SNK)

  ! Set up vector parameter list
  CALL vecPList%add('VectorType -> n',2)
  CALL vecPList%add('VectorType -> MPI_Comm_ID',PE_COMM_SELF)

  !Configure exception handler for test
  CALL e%setStopOnError(.FALSE.)
  CALL e%setQuietMode(.TRUE.)
  CALL eParams%addSurrogate(e)
  CALL eLinearSolverType%addSurrogate(e)
  CALL mpiTestEnv%init(PE_COMM_SELF)

  CREATE_TEST('Test Linear Solvers')

  REGISTER_SUBTEST('testClear',testClear)
  REGISTER_SUBTEST('testInit',testInit)
  REGISTER_SUBTEST('testIterativeSolve_Multigrid',testIterativeSolve_Multigrid)

  FINALIZE_TEST()

  CALL pList%clear()
  CALL vecPList%clear()
  CALL optListMat%clear()
  CALL optListLS%clear()

#ifdef FUTILITY_HAVE_PETSC
  CALL PetscFinalize(ierr)
#else
  CALL mpiTestEnv%finalize()
#endif
!
!===============================================================================
CONTAINS
!
!-------------------------------------------------------------------------------
    SUBROUTINE testClear()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      LOGICAL(SBK) :: bool
    !test Iterative
      ALLOCATE(LinearSolverType_Multigrid :: thisLS)
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Multigrid)
        !first build one by hand to test
        thisLS%isInit=.TRUE.
        thisLS%TPLType=PETSC
        thisLS%solverMethod=1
        thisLS%info=2
        thisLS%normType=2
        thisLS%maxIters=2
        thisLS%iters=2
        thisLS%convTol=2._SRK
        thisLS%residual=2._SRK
        thisLS%isDecomposed=.TRUE.
        CALL thisLS%MPIparallelEnv%init(PE_COMM_SELF)
        CALL thisLS%OMPparallelEnv%init(1)
#ifdef FUTILITY_HAVE_PETSC
        CALL KSPCreate(thisLS%MPIparallelEnv%comm,thisLS%ksp,ierr)
#endif

        ! initialize matrix A
        ALLOCATE(DenseSquareMatrixType :: thisLS%A)
        CALL pList%clear()
        CALL pList%add('MatrixType->n',2_SNK)
        CALL pList%add('MatrixType->isSym',.TRUE.)
        CALL pList%validate(pList,optListMat)
        CALL thisLS%A%init(pList) !2x2, symmetric

        ! initialize preconditioner
        ALLOCATE(ILU_PreCondType :: thisLS%PreCondType)
        CALL thisLS%setupPC() !ZZZZ

        ! initialize vector X
        CALL vecPList%set('VectorType -> n',2)
        ALLOCATE(RealVectorType :: thisLS%X)
        CALL thisLS%X%init(vecPList)

        ! initialize vector b
        ALLOCATE(RealVectorType :: thisLS%b)
        CALL thisLS%b%init(vecPList)

        ! initialize matrix M
        ALLOCATE(DenseSquareMatrixType :: thisLS%M)
        CALL pList%clear()
        CALL pList%add('MatrixType->n',10_SNK)
        CALL pList%add('MatrixType->isSym',.TRUE.)
        CALL thisLS%M%init(pList)
      ENDSELECT

      CALL thisLS%clear()

      !check results
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Multigrid)
        bool=thisLS%isInit .OR.thisLS%solverMethod == 1                  &
            .OR. ALLOCATED(thisLS%M) .OR. ASSOCIATED(thisLS%A)           &
            .OR. ASSOCIATED(thisLS%X) .OR. thisLS%info /= 0               &
            .OR. thisLS%normType == 2 .OR. thisLS%maxIters == 2          &
            .OR. thisLS%iters == 2 .OR. thisLS%isDecomposed              &
            .OR. thisLS%residual == 2._SRK .OR. thisLS%convTol == 2._SRK &
            .OR. ALLOCATED(thisLS%PreCondType)
        ASSERT(.NOT.(bool),'CALL Multigrid%clear() FAILED!')
      ENDSELECT
      CALL thisLS%clear()
      DEALLOCATE(thisLS)

    ENDSUBROUTINE testClear
!
!-------------------------------------------------------------------------------
    SUBROUTINE testInit()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      INTEGER(SIK) :: nerrors1,nerrors2
      LOGICAL(SBK) :: bool
      ALLOCATE(LinearSolverType_Multigrid :: thisLS)

      !first test a correct use case, with timer name
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',MULTIGRID)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',2_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',2_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->PCType','ILU')
      CALL pList%add('LinearSolverType->PCIters',-1)
      CALL pList%add('LinearSolverType->PCSetup',0)

      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      SELECTTYPE(thisLS); TYPE IS(LinearSolvertype_Multigrid)
        bool = (thisLS%isInit .AND. thisLS%solverMethod == MULTIGRID &
           .AND. thisLS%MPIparallelEnv%isInit() &
           .AND. thisLS%OMPparallelEnv%isInit() )
        ASSERT(bool, 'Iterative%init(...)')
      ENDSELECT
      CALL thisLS%clear()

    ENDSUBROUTINE testInit
!
!-------------------------------------------------------------------------------
    SUBROUTINE testIterativeSolve_Multigrid()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      REAL(SRK),ALLOCATABLE :: soln(:),b(:)
      REAL(SRK),ALLOCATABLE :: A_temp(:,:)
      REAL(SRK),POINTER :: x(:)
      INTEGER(SIK) :: i
      LOGICAL(SBK) :: match

      INTEGER(SIK),PARAMETER :: n=65_SNK

      ALLOCATE(LinearSolverType_Multigrid :: thisLS)

#ifdef FUTILITY_HAVE_PETSC
      !The PETSC sparse matrix type
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',MULTIGRID)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',n)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',n*3-2)
      CALL pList%add('LinearSolverType->x->VectorType->n',n)
      CALL pList%add('LinearSolverType->b->VectorType->n',n)
      CALL pList%validate(pList,optListLS)

      !Geometry dimensions for multigrid:
      CALL pList%add('LinearSolverType->Multigrid->nx',n)
      CALL pList%add('LinearSolverType->Multigrid->ny',1)
      CALL pList%add('LinearSolverType->Multigrid->nz',1)
      CALL pList%add('LinearSolverType->Multigrid->num_eqns',1)

      !TODO make this test problem a coupled system of 2 equations

      CALL thisLS%init(pList)
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Multigrid)
        CALL thisLS%setupInterpMats(pList)
      ENDSELECT

      !A is a tridiagonal system with -1 on the offdiagonals, and
      !  2.5 on the diagonals.
      ALLOCATE(A_temp(n,n))
      A_temp=0.0_SRK
      SELECTTYPE(A => thisLS%A); TYPE IS(PETScMatrixType)
      DO i=1,n
        IF(i > 1) THEN
          CALL A%set(i,i-1,-1.0_SRK)
          A_temp(i,i-1)=-1.0_SRK
        ENDIF
        CALL A%set(i,i,2.5_SRK)
        A_temp(i,i)=2.5_SRK
        IF(i < 1) THEN
          CALL A%set(i,i+1,-1.0_SRK)
          A_temp(i,i+1)=-1.0_SRK
        ENDIF
      ENDDO
      ENDSELECT

      ! Create solution:
      ALLOCATE(soln(n))
      soln=1.0_SRK
      soln(n/2)=2.0_SRK

      ALLOCATE(b(n))
      b=MATMUL(A_temp,soln)
      DEALLOCATE(A_temp)

      SELECTTYPE(LS_b => thisLS%b); TYPE IS(PETScVectorType)
        CALL LS_b%setAll_array(b)
      ENDSELECT
      DEALLOCATE(b)

      ! build x0
      ALLOCATE(x(n))
      x(1:(n-1)/2)=0.5_SRK
      x((n+1)/2:n)=1.1_SRK
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Multigrid)
        CALL thisLS%setX0(x)
      ENDSELECT

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Multigrid)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT

      !solve it
      CALL thisLS%solve()

      SELECTTYPE(LS_x => thisLS%X); TYPE IS(PETScVectorType)
        CALL LS_x%getAll(x)
      ENDSELECT
      match=ALL(ABS(x-soln) < 1.0E-6_SRK)
      ASSERT(match, 'PETScIterative%solve() - Multigrid')

      DEALLOCATE(soln)
      DEALLOCATE(x)
      CALL thisLS%A%clear()
      CALL thisLS%clear()

#endif
      DEALLOCATE(thisLS)

    ENDSUBROUTINE testIterativeSolve_Multigrid

ENDPROGRAM testLinearSolver_Multigrid
