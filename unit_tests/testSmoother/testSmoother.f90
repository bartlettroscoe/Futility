!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PROGRAM testSmoother
#include "UnitTest.h"
  USE UnitTest
  USE IntrType
  USE ExceptionHandler
  USE ParameterLists
  USE ParallelEnv
  USE VectorTypes
  USE MatrixTypes
  USE SmootherTypes

  IMPLICIT NONE

  TYPE(ExceptionHandlerType),TARGET :: e
  TYPE(MPI_EnvType) :: mpiTestEnv

#ifdef FUTILITY_HAVE_PETSC
#include <finclude/petsc.h>
#undef IS
  PetscErrorCode  :: ierr

  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#endif

#ifdef HAVE_MPI
  CALL mpiTestEnv%init(MPI_COMM_WORLD)
#endif

  !Configure exception handler for test
  CALL e%setStopOnError(.FALSE.)
  CALL e%setQuietMode(.TRUE.)
  CALL eParams%addSurrogate(e)
  CALL eSmootherType%addSurrogate(e)

  CREATE_TEST('Test SmootherType_PETSc_RBBJ')
  REGISTER_SUBTEST('testInit',testInit_PETSc_RBBJ)
  REGISTER_SUBTEST('testClear',testClear_PETSc_RBBJ)
  REGISTER_SUBTEST('testSmooth',testSmooth_PETSc_RBBJ)
  FINALIZE_TEST()

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
    SUBROUTINE testClear_PETSc_RBBJ()
      TYPE(SmootherType_PETSc_RBBJ) :: smoother

      ALLOCATE(smoother%is_red(500))
      smoother%MPIparallelEnv=mpiTestEnv

      CALL smoother%clear()

      ASSERT(.NOT.ALLOCATED(smoother%is_red),'Smoother%isred deallocated')
      ASSERT(.NOT. smoother%MPIparallelEnv%isInit(),'Smoother parenv not initialized')
      ASSERT(.NOT. smoother%isInit,'Smoother not initialized')
    ENDSUBROUTINE testClear_PETSc_RBBJ
!
!-------------------------------------------------------------------------------
    SUBROUTINE testInit_PETSc_RBBJ()
      TYPE(SmootherType_PETSc_RBBJ) :: smoother
      TYPE(ParamType) :: params
      LOGICAL(SBK) :: tmpbool
      INTEGER(SIK) :: num_red,num_black,istt,istp,blk_size
      INTEGER(SIK) :: i
      LOGICAL(SBK),ALLOCATABLE :: is_red(:)
#ifdef FUTILITY_HAVE_PETSC
      KSP :: ksp
      PC :: pc
      KSPType :: myksptype
      PCType :: mypctype
      PetscErrorCode :: iperr

      CALL KSPCreate(mpiTestEnv%comm,ksp,iperr)
      CALL KSPSetType(ksp,KSPGMRES,iperr)
      CALL KSPGetPC(ksp,pc,iperr)
      CALL PCSetType(pc,PCJACOBI,iperr)
#else
      INTEGER(SIK) :: ksp=-1_SIK
      INTEGER(SIK) :: pc
#endif
      num_red=250
      num_black=250
      istt=1
      istp=500
      blk_size=1_SIK
      ALLOCATE(is_red(istt:istp))
      is_red=.FALSE.
      DO i=istt,istp,2
        is_red(i)=.TRUE.
      ENDDO

      CALL params%add('SmootherType->is_red',is_red)
      CALL params%add('SmootherType->istt',istt)
      CALL params%add('SmootherType->istp',istp)
      CALL params%add('SmootherType->blk_size',blk_size)
      CALL params%add('SmootherType->MPI_Comm_ID',mpiTestEnv%comm)
      CALL smoother%init(ksp,params)
      CALL params%clear()

      tmpbool=ALLOCATED(smoother%is_red) .AND. &
                LBOUND(smoother%is_red,DIM=1) == istt .AND. &
                UBOUND(smoother%is_red,DIM=1) == istp
      ASSERT(tmpbool,'Smoother%isred allocated with correct bounds')

      tmpbool=smoother%num_red == num_red .AND. &
                smoother%num_black == num_black .AND. &
                smoother%istt == istt .AND. &
                smoother%istp == istp .AND. &
                smoother%blk_size == blk_size
      ASSERT(tmpbool,'smoother parameters are the correct value')

#ifdef FUTILITY_HAVE_PETSC
      CALL KSPGetType(ksp,myksptype,iperr)
      CALL PCGetType(pc,mypctype,iperr)
      tmpbool=iperr == 0_SIK .AND. &
              myksptype == KSPRICHARDSON .AND. &
              mypctype == PCSHELL
      ASSERT(tmpbool,'ksp and pc set to the correct types in PETSc')
      CALL KSPGetType(smoother%ksp,myksptype,iperr)
      CALL PCGetType(smoother%pc,mypctype,iperr)
      tmpbool=iperr == 0_SIK .AND. &
              myksptype == KSPRICHARDSON .AND. &
              mypctype == PCSHELL
      ASSERT(tmpbool,"smoother's copy of ksp and pc set to the correct types in PETSc")
#endif

      tmpbool=smoother%MPIparallelEnv%isInit()
      ASSERT(tmpbool,'smoother ParEnv initialized.')

      ASSERT(smoother%isInit,'Smoother  initialized')

      CALL smoother%clear()
      CALL KSPDestroy(ksp,iperr)

    ENDSUBROUTINE testInit_PETSc_RBBJ
!
!-------------------------------------------------------------------------------
    SUBROUTINE testInit_PETSc_RBBJ()
      TYPE(SmootherType_PETSc_RBBJ) :: smoother
      TYPE(ParamType) :: params
      LOGICAL(SBK) :: tmpbool
      INTEGER(SIK) :: num_red,num_black,istt,istp,blk_size
      INTEGER(SIK) :: i
      LOGICAL(SBK),ALLOCATABLE :: is_red(:)
#ifdef FUTILITY_HAVE_PETSC
      KSP :: ksp
      PetscErrorCode :: iperr

      CALL KSPCreate(mpiTestEnv%comm,ksp,iperr)
#else
      INTEGER(SIK) :: ksp=-1_SIK
#endif

      num_red=250
      num_black=250
      istt=1
      istp=500
      blk_size=1_SIK
      ALLOCATE(is_red(istt:istp))
      is_red=.FALSE.
      DO i=istt,istp,2
        is_red(i)=.TRUE.
      ENDDO

      CALL params%add('SmootherType->is_red',is_red)
      CALL params%add('SmootherType->istt',istt)
      CALL params%add('SmootherType->istp',istp)
      CALL params%add('SmootherType->blk_size',blk_size)
      CALL params%add('SmootherType->MPI_Comm_ID',mpiTestEnv%comm)
      CALL smoother%init(ksp,params)
      CALL params%clear()

      tmpbool=ALLOCATED(smoother%is_red) .AND. &
                LBOUND(smoother%is_red,DIM=1) == istt .AND. &
                UBOUND(smoother%is_red,DIM=1) == istp
      ASSERT(tmpbool,'Smoother%isred allocated with correct bounds')

      tmpbool=smoother%num_red == num_red .AND. &
                smoother%num_black == num_black .AND. &
                smoother%istt == istt .AND. &
                smoother%istp == istp .AND. &
                smoother%blk_size == blk_size
      ASSERT(tmpbool,'smoother parameters are the correct value')

#ifdef FUTILITY_HAVE_PETSC
      CALL KSPGetType(ksp,myksptype,iperr)
      CALL PCGetType(pc,mypctype,iperr)
      tmpbool=iperr == 0_SIK .AND. &
              myksptype == KSPRICHARDSON .AND. &
              mypctype == PCSHELL
      ASSERT(tmpbool,'ksp and pc set to the correct types in PETSc')
      CALL KSPGetType(smoother%ksp,myksptype,iperr)
      CALL PCGetType(smoother%pc,mypctype,iperr)
      tmpbool=iperr == 0_SIK .AND. &
              myksptype == KSPRICHARDSON .AND. &
              mypctype == PCSHELL
      ASSERT(tmpbool,"smoother's copy of ksp and pc set to the correct types in PETSc")
#endif

      tmpbool=smoother%MPIparallelEnv%isInit()
      ASSERT(tmpbool,'smoother ParEnv initialized.')

      ASSERT(smoother%isInit,'Smoother  initialized')

      CALL smoother%clear()
      CALL KSPDestroy(ksp,iperr)

    ENDSUBROUTINE testInit_PETSc_RBBJ

ENDPROGRAM testSmoother
