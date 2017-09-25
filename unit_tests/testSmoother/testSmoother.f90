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
  TYPE(ParamType) :: params
  INTEGER(SIK) :: istt,istp,blk_size,num_colors

#ifdef FUTILITY_HAVE_PETSC
#include <finclude/petsc.h>
#undef IS
  PetscErrorCode  :: ierr

  KSP :: ksp
  PetscErrorCode :: iperr

  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#else
  INTEGER(SIK) :: ksp=-1_SIK
#endif

#ifdef HAVE_MPI
  CALL mpiTestEnv%init(MPI_COMM_WORLD)
#endif

  !Configure exception handler for test
  CALL e%setStopOnError(.FALSE.)
  CALL e%setQuietMode(.TRUE.)
  CALL eParams%addSurrogate(e)
  CALL eSmootherType%addSurrogate(e)


  istt=1
  istp=65
  blk_size=1_SIK
  num_colors=2_SIK
  CALL params%clear()
  CALL params%add('SmootherType->istt',istt)
  CALL params%add('SmootherType->istp',istp)
  CALL params%add('SmootherType->num_colors',num_colors)
  CALL params%add('SmootherType->blk_size',blk_size)
  CALL params%add('SmootherType->MPI_Comm_ID',mpiTestEnv%comm)
#ifdef FUTILITY_HAVE_PETSC
  CALL KSPCreate(mpiTestEnv%comm,ksp,iperr)
#endif

  CREATE_TEST('Test SmootherType_PETSc_CBJ')
  REGISTER_SUBTEST('testInit',testInit_PETSc_CBJ)
  REGISTER_SUBTEST('testClear',testClear_PETSc_CBJ)
  REGISTER_SUBTEST('testDefineColor',testDefineColor_PETSc_CBJ)
  REGISTER_SUBTEST('testDefineAllColors',testDefineAllColors_PETSc_CBJ)
  REGISTER_SUBTEST('testSmooth',testSmooth_PETSc_CBJ)
  FINALIZE_TEST()

  CALL params%clear()
#ifdef FUTILITY_HAVE_PETSC
  CALL KSPDestroy(ksp,iperr)
  CALL PetscFinalize(ierr)
#else
  CALL mpiTestEnv%finalize()
#endif
!
!===============================================================================
  CONTAINS
!
!-------------------------------------------------------------------------------
    SUBROUTINE testClear_PETSc_CBJ()
      TYPE(SmootherType_PETSc_CBJ) :: smoother

      ALLOCATE(smoother%color_ids(500))
      smoother%MPIparallelEnv=mpiTestEnv

      CALL smoother%clear()

      ASSERT(.NOT.ALLOCATED(smoother%color_ids),'Smoother%isred deallocated')
      ASSERT(.NOT. smoother%MPIparallelEnv%isInit(),'Smoother parenv not initialized')
      ASSERT(.NOT. smoother%isInit,'Smoother not initialized')
    ENDSUBROUTINE testClear_PETSc_CBJ
!
!-------------------------------------------------------------------------------
    SUBROUTINE testInit_PETSc_CBJ()
      TYPE(SmootherType_PETSc_CBJ) :: smoother
      LOGICAL(SBK) :: tmpbool
#ifdef FUTILITY_HAVE_PETSC
      PC :: pc
      KSPType :: myksptype
      PCType :: mypctype

      CALL KSPSetType(ksp,KSPGMRES,iperr)
      CALL KSPGetPC(ksp,pc,iperr)
      CALL PCSetType(pc,PCJACOBI,iperr)
#endif

      CALL smoother%init(ksp,params)

      tmpbool=ALLOCATED(smoother%color_ids) .AND. &
                LBOUND(smoother%color_ids,DIM=1) == istt .AND. &
                UBOUND(smoother%color_ids,DIM=1) == istp
      ASSERT(tmpbool,'Smoother%color_ids allocated with correct bounds')

      tmpbool=smoother%num_colors == num_colors .AND. &
                smoother%istt == istt .AND. &
                smoother%istp == istp .AND. &
                smoother%blk_size == blk_size .AND. &
                smoother%TPLType == PETSc .AND. &
                smoother%smootherMethod == CBJ .AND. &
                smoother%blockMethod == LU
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

    ENDSUBROUTINE testInit_PETSc_CBJ
!
!-------------------------------------------------------------------------------
    SUBROUTINE testSmooth_PETSc_CBJ
      TYPE(SmootherType_PETSc_CBJ) :: smoother

      CALL smoother%init(ksp,params)

      CALL smoother%clear()
    ENDSUBROUTINE testSmooth_PETSc_CBJ
!
!-------------------------------------------------------------------------------
    SUBROUTINE testDefineColor_PETSc_CBJ
      TYPE(SmootherType_PETSc_CBJ) :: smoother
      INTEGER(SIK),PARAMETER :: num_indices=25_SIK
      INTEGER(SIK) :: index_list(num_indices)
      INTEGER(SIK) :: icolor,i
      LOGICAL(SBK) :: tmpbool

      CALL smoother%init(ksp,params)

      DO i=1,num_indices
        index_list(i)=2*i-1
      ENDDO
      smoother%color_ids=1_SIK

      icolor=2_SIK
      CALL smoother%defineColor(icolor,index_list)

      ASSERT(smoother%colors(icolor)%num_indices == num_indices,'Correct # of indices')
      ASSERT(ALL(smoother%colors(icolor)%index_list == index_list),'Correct index list')
      tmpbool=.TRUE.
      DO i=1,num_indices
        IF(smoother%color_ids(index_list(i)) /= icolor) THEN
          tmpbool=.FALSE.
          EXIT
        ENDIF
      ENDDO
      ASSERT(tmpbool,'Correct color ids for color 2')

      CALL smoother%clear()

    ENDSUBROUTINE testDefineColor_PETSc_CBJ
!
!-------------------------------------------------------------------------------
    SUBROUTINE testDefineAllColors_PETSc_CBJ
      TYPE(SmootherType_PETSc_CBJ) :: smoother
      LOGICAL(SBK) :: tmpbool
      INTEGER(SIK) :: i
      INTEGER(SBK),ALLOCATABLE :: color_ids(:)

      CALL smoother%init(ksp,params)

      ALLOCATE(color_ids(istt:istp))
      color_ids=2
      DO i=istt,istp,2
        color_ids(i)=1
      ENDDO
      CALL smoother%defineAllColors(color_ids)

      tmpbool=ALL(smoother%color_ids == color_ids)
      ASSERT(tmpbool,'Smoother%color_ids correctly set')
      tmpbool=smoother%colors(1)%num_indices == 33_SIK .AND. &
                smoother%colors(2)%num_indices == 32_SIK
      ASSERT(tmpbool,'Correct number of red and black indices')

      CALL smoother%clear()

    ENDSUBROUTINE testDefineAllColors_PETSc_CBJ

ENDPROGRAM testSmoother
