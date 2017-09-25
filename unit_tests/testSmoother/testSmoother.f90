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
  TYPE(ParamType) :: params,params_2proc
  INTEGER(SIK) :: istt,istp,blk_size,num_colors
  INTEGER(SIK) :: mpierr

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
  blk_size=2_SIK
  num_colors=2_SIK
  CALL params%clear()
  CALL params%add('SmootherType->istt',istt)
  CALL params%add('SmootherType->istp',istp)
  CALL params%add('SmootherType->num_colors',num_colors)
  CALL params%add('SmootherType->blk_size',blk_size)
  CALL params%add('SmootherType->MPI_Comm_ID',mpiTestEnv%comm)
  CALL params_2proc%clear()
  params_2proc=params
  IF(.NOT. mpiTestEnv%master) THEN
    CALL params_2proc%set('SmootherType->istt',istt+istp)
    CALL params_2proc%set('SmootherType->istp',istp+istp)
  ENDIF
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
!
!-------------------------------------------------------------------------------
    SUBROUTINE testSmooth_PETSc_CBJ
#ifdef FUTILITY_HAVE_PETSC
#ifdef HAVE_MPI
      !Test smoother on 2-proc,2-group problem
      TYPE(SmootherType_PETSc_CBJ) :: smoother

      INTEGER(SIK),PARAMETER :: nlocal=65_SIK
      INTEGER(SIK),PARAMETER :: num_eqns=2_SIK
      INTEGER(SIK) :: n
      REAL(SRK) :: b(nlocal*num_eqns)
      INTEGER(SIK) :: nstart,nend
      INTEGER(SIK) :: i,ieqn
      INTEGER(SIK) :: row,col
      INTEGER(SIK) :: dnnz(nlocal*num_eqns),onnz(nlocal*num_eqns)

      Mat :: A_petsc
      Vec :: b_petsc,x_petsc
      PetscErrorCode :: iperr

      n=nlocal*mpiTestEnv%nproc

      dnnz=num_eqns+2
      dnnz(1:num_eqns)=num_eqns+1
      dnnz((nlocal-1)*num_eqns+1:nlocal*num_eqns)=num_eqns+1
      onnz=0_SIK
      nstart=mpiTestEnv%rank*nlocal+1
      nend=(mpiTestEnv%rank+1)*nlocal
      IF(mpiTestEnv%rank < mpiTestEnv%nproc-1) THEN
        onnz((nend-1)*num_eqns+1:nend*num_eqns)=1_SIK
      ENDIF
      IF(mpiTestEnv%rank > 0) THEN
        onnz(1:num_eqns)=1_SIK
      ENDIF

      CALL MatCreate(MPI_COMM_WORLD,A_petsc,iperr)
      CALL VecCreate(MPI_COMM_WORLD,b_petsc,iperr)
      CALL VecCreate(MPI_COMM_WORLD,x_petsc,iperr)
      CALL MatSetSizes(A_petsc,nlocal*num_eqns,nlocal*num_eqns, &
                        n*num_eqns,n*num_eqns,iperr)
      CALL VecSetSizes(b_petsc,nlocal*num_eqns,n*num_eqns,iperr)
      CALL VecSetSizes(x_petsc,nlocal*num_eqns,n*num_eqns,iperr)
      CALL VecSetType(b_petsc,VECMPI,iperr)
      CALL MatSetType(A_petsc,MATMPIAIJ,ierr)
      CALL VecSetType(x_petsc,VECMPI,iperr)
      CALL MatMPIAIJSetPreallocation(A_petsc,0,dnnz,0,onnz,iperr)
      CALL VecSet(b_petsc,0.0_SRK,iperr)
      CALL VecSet(x_petsc,0.0_SRK,iperr)
      CALL VecAssemblyBegin(b_petsc,iperr)


      !2G homogeneous diffusion problem:
      DO i=nstart,nend
        DO ieqn=1,num_eqns
          row=(i-1)*num_eqns+ieqn
          IF(i > nstart) THEN
            col=(i-2)*num_eqns+ieqn
            CALL MatSetValue(A_petsc,row-1,col-1,-1.0_SRK,INSERT_VALUES,iperr)
          ENDIF
          IF(i < nend) THEN
            col=i*num_eqns+ieqn
            CALL MatSetValue(A_petsc,row-1,col-1,-1.0_SRK,INSERT_VALUES,iperr)
          ENDIF
        ENDDO
        row=(i-1)*num_eqns+1
        CALL MatSetValue(A_petsc,row-1,row-1,2.1_SRK,INSERT_VALUES,iperr)
        CALL MatSetValue(A_petsc,row,row,2.2_SRK,INSERT_VALUES,iperr)
        !Upscatter:
        CALL MatSetValue(A_petsc,row-1,row,-0.01_SRK,INSERT_VALUES,iperr)
        !Downscatter:
        CALL MatSetValue(A_petsc,row,row-1,-0.5_SRK,INSERT_VALUES,iperr)
      ENDDO

      CALL VecAssemblyBegin(x_petsc,iperr)
      CALL MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY,iperr)

      CALL VecAssemblyEnd(b_petsc,iperr)
      CALL VecAssemblyEnd(x_petsc,iperr)
      CALL MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY,iperr)

      CALL KSPSetOperators(ksp,A_petsc,A_petsc,iperr)

      CALL smoother%init(ksp,params_2proc)

      !CALL PCShellSetSetUp(smoother%pc,smoother%PCSetup_CBJ,iperr)

      CALL smoother%clear()

      CALL MatDestroy(A_petsc,iperr)
      CALL VecDestroy(b_petsc,iperr)
      CALL VecDestroy(x_petsc,iperr)
#endif
#endif

    ENDSUBROUTINE testSmooth_PETSc_CBJ

ENDPROGRAM testSmoother
