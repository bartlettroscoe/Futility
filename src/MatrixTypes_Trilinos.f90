!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Implementations of MatrixTypes with Trilinos support
!>
!> @par Module Dependencies
!>  - @ref IntrType "IntrType": @copybrief IntrType
!>  - @ref ExceptionHandler "ExceptionHandler": @copybrief ExceptionHandler
!>
!> @author Adam Nelson and Brendan Kochunas
!>   @date 02/14/2012
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE MatrixTypes_Trilinos
  USE IntrType
  USE ExceptionHandler
  USE ParameterLists
  USE MatrixTypes_Base
  USE trilinos_interfaces

  IMPLICIT NONE
  PRIVATE

!
! List of public members
  PUBLIC :: TrilinosMatrixType
  
  TYPE,EXTENDS(DistributedMatrixType) :: TrilinosMatrixType
#ifdef FUTILITY_HAVE_Trilinos
    INTEGER(SIK) :: A
    INTEGER(SIK) :: currow
    INTEGER(SIK) :: ncol
    INTEGER(SIK),ALLOCATABLE :: jloc(:)
    REAL(SRK),ALLOCATABLE :: aloc(:)
#endif
!
!List of Type Bound Procedures
    CONTAINS
      !> @copybrief MatrixTypes::clear_TrilinosMatrixType
      !> @copydetails MatrixTypes::clear_TrilinosMatrixType
      PROCEDURE,PASS :: clear => clear_TrilinosMatrixType
      !> @copybrief MatrixTypes::init_TrilinosMatrixType
      !> @copydetails MatrixTypes::init_TrilinosMatrixType
      PROCEDURE,PASS :: init => init_TrilinosMatrixParam
      !> @copybrief MatrixTypes::set_TrilinosMatrixType
      !> @copydetails MatrixTypes::set_TrilinosMatrixType
      PROCEDURE,PASS :: set => set_TrilinosMatrixType
      !> @copybrief MatrixTypes::set_TrilinosMatrixType
      !> @copydetails MatrixTypes::set_TrilinosMatrixType
      PROCEDURE,PASS :: setShape => setShape_TrilinosMatrixType
      !> @copybrief MatrixTypes::get_TrilinosMatrixType
      !> @copydetails MatrixTypes::get_TrilinosMatrixType
      PROCEDURE,PASS :: get => get_TrilinosMatrixType
      !> @copybrief MatrixTypes::assemble_TrilinosMatrixType
      !> @copydetails MatrixTypes::assemble_TrilinosMatrixType
      PROCEDURE,PASS :: assemble => assemble_TrilinosMatrixType
      !> @copybrief MatrixTypes::transpose_TrilinosMatrixType
      !> @copydetails MatrixTypes::transpose_TrilinosMatrixType
      PROCEDURE,PASS :: transpose => transpose_TrilinosMatrixType
  ENDTYPE TrilinosMatrixType

  !> Name of module
  CHARACTER(LEN=*),PARAMETER :: modName='MATRIXTYPES_TRILINOS'
!
!===============================================================================
  CONTAINS
!
!-------------------------------------------------------------------------------
!> @brief Initializes Trilinos Matrix Type with a Parameter List
!> @param matrix the matrix type to act on
!> @param pList the parameter list
!>
    SUBROUTINE init_TrilinosMatrixParam(matrix,Params)
      CHARACTER(LEN=*),PARAMETER :: myName='init_TrilinosMatrixParam'
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: matrix
      CLASS(ParamType),INTENT(IN) :: Params
      TYPE(ParamType) :: validParams
      INTEGER(SIK) :: n, matType, MPI_COMM_ID, nlocal, ierr, rnnz
      INTEGER(SIK),ALLOCATABLE :: dnnz(:), onnz(:)
      LOGICAL(SBK) :: isSym

#ifdef FUTILITY_HAVE_Trilinos

      !Check to set up required and optional param lists.
      IF(.NOT.MatrixType_Paramsflag) CALL MatrixTypes_Declare_ValidParams()
      !Validate against the reqParams and OptParams
      validParams=Params
      CALL validParams%validate(DistributedMatrixType_reqParams, &
          DistributedMatrixType_optParams)

      ! Pull Data From Parameter List
      CALL validParams%get('MatrixType->n',n)
      CALL validParams%get('MatrixType->isSym',isSym)
      CALL validParams%get('MatrixType->matType',matType)
      CALL validParams%get('MatrixType->MPI_COMM_ID',MPI_COMM_ID)
      CALL validParams%get('MatrixType->nlocal',nlocal)
      ALLOCATE(dnnz(nlocal))
      ALLOCATE(onnz(nlocal))
      CALL validParams%get('MatrixType->dnnz',dnnz)
      CALL validParams%get('MatrixType->onnz',onnz)
      CALL validParams%clear()

      rnnz=MAXVAL(dnnz)+MAXVAL(onnz)
      IF(rnnz==-2) rnnz=n
      IF(.NOT. matrix%isInit) THEN
        IF(n < 1) THEN
          CALL eMatrixType%raiseError('Incorrect input to '// &
            modName//'::'//myName//' - Number of rows (n) must be '// &
              'greater than 0!')
        ELSEIF(nlocal < 1) THEN
          CALL eMatrixType%raiseError('Incorrect input to '// &
            modName//'::'//myName//' - Number of local rows (nlocal) must '// &
              'be greater than 0!')
        ELSEIF(rnnz < 1) THEN
          CALL eMatrixType%raiseError('Incorrect input to '// &
            modName//'::'//myName//' - Number of non-zero elements (dnnz,onnz) '// &
              'must be greater than 0!')
        ELSEIF(isSym) THEN
          CALL eMatrixType%raiseError('Incorrect input to '// &
            modName//'::'//myName//' - Symmetric matrices are not supported.')
        ELSE
          matrix%isInit=.TRUE.
          matrix%n=n
          matrix%comm=MPI_COMM_ID
          matrix%isAssembled=.FALSE.
          matrix%nlocal=nlocal
          matrix%currow=0
          matrix%ncol=0
          ALLOCATE(matrix%jloc(rnnz))
          ALLOCATE(matrix%aloc(rnnz))
          IF(isSym) THEN
            matrix%isSymmetric=.TRUE.
          ELSE
            matrix%isSymmetric=.FALSE.
          ENDIF
          IF(.NOT.matrix%isCreated) THEN
            CALL ForPETRA_MatInit(matrix%A,n,nlocal,rnnz,matrix%comm)
            matrix%isCreated=.TRUE.
          ENDIF

          IF (matType /= SPARSE) THEN
            CALL eMatrixType%raiseError('Invalid matrix type in '// &
              modName//'::'//myName//' - Only sparse square '// &
              'matrices are available with Trilinos.')
          ENDIF
        ENDIF
      ELSE
        CALL eMatrixType%raiseError('Incorrect call to '// &
          modName//'::'//myName//' - MatrixType already initialized')
      ENDIF
#else
      CALL eMatrixType%raiseFatalError('Incorrect call to '// &
              modName//'::'//myName//' - Trilinos not enabled.  You will'// &
              'need to recompile with Trilinos enabled to use this feature.')
#endif
    ENDSUBROUTINE init_TrilinosMatrixParam
!
!-------------------------------------------------------------------------------
!> @brief Clears the Trilinos sparse matrix
!> @param matrix the matrix type to act on
!>
    SUBROUTINE clear_TrilinosMatrixType(matrix)
      CHARACTER(LEN=*),PARAMETER :: myName='clear_TrilinosMatrixType'
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: matrix
#ifdef FUTILITY_HAVE_Trilinos

      !TODO add routine to clear memory
      matrix%isInit=.FALSE.
      matrix%n=0
      matrix%isAssembled=.FALSE.
      matrix%isCreated=.FALSE.
      matrix%isSymmetric=.FALSE.
      IF(ALLOCATED(matrix%jloc)) DEALLOCATE(matrix%jloc)
      IF(ALLOCATED(matrix%aloc)) DEALLOCATE(matrix%aloc)
      matrix%currow=0
      matrix%ncol=0
      CALL ForPETRA_MatDestroy(matrix%a)
      matrix%A=-1
#else
      CALL eMatrixType%raiseFatalError('Incorrect call to '// &
              modName//'::'//myName//' - Trilinos not enabled.  You will'// &
              'need to recompile with Trilinos enabled to use this feature.')
#endif
    ENDSUBROUTINE clear_TrilinosMatrixType
!
!-------------------------------------------------------------------------------
!> @brief Sets the values in the Trilinos matrix
!> @param declares the matrix type to act on
!> @param i the ith location in the matrix
!> @param j the jth location in the matrix
!> @param setval the value to be set
!>
    SUBROUTINE set_TrilinosMatrixType(matrix,i,j,setval)
      CHARACTER(LEN=*),PARAMETER :: myName='set_TrilinosMatrixType'
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: matrix
      INTEGER(SIK),INTENT(IN) :: i
      INTEGER(SIK),INTENT(IN) :: j
      REAL(SRK),INTENT(IN) :: setval
#ifdef FUTILITY_HAVE_Trilinos
      INTEGER(SIK)  :: ierr

      IF(matrix%isInit) THEN
        IF(((j <= matrix%n) .AND. (i <= matrix%n)) &
          .AND. ((j > 0) .AND. (i > 0))) THEN
          IF(matrix%isAssembled) CALL ForPETRA_MatReset(matrix%A)
          IF(i==matrix%currow) THEN
            matrix%ncol=matrix%ncol+1
            matrix%jloc(matrix%ncol)=j
            matrix%aloc(matrix%ncol)=setval
          ELSE
            IF(matrix%currow>0) THEN
              CALL ForPETRA_MatSet(matrix%A,matrix%currow,matrix%ncol,matrix%jloc,matrix%aloc)
            ENDIF
            matrix%jloc=0
            matrix%aloc=0.0_SRK
            !Need to store index from the incomming data
            matrix%ncol=1
            matrix%jloc(1)=j
            matrix%aloc(1)=setval
            matrix%currow=i
          ENDIF
!TODO
          matrix%isAssembled=.FALSE.
        ENDIF
      ENDIF
#else
      CALL eMatrixType%raiseFatalError('Incorrect call to '// &
              modName//'::'//myName//' - Trilinos not enabled.  You will'// &
              'need to recompile with Trilinos enabled to use this feature.')
#endif
    ENDSUBROUTINE set_TrilinosMatrixType
!
!-------------------------------------------------------------------------------
!> @brief Sets the values in the Trilinos matrix
!> @param declares the matrix type to act on
!> @param i the ith location in the matrix
!> @param j the jth location in the matrix
!> @param setval the value to be set
!>
    SUBROUTINE setShape_TrilinosMatrixType(matrix,i,j,setval)
      CHARACTER(LEN=*),PARAMETER :: myName='set_TrilinosMatrixType'
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: matrix
      INTEGER(SIK),INTENT(IN) :: i
      INTEGER(SIK),INTENT(IN) :: j
      REAL(SRK),INTENT(IN) :: setval
#ifdef FUTILITY_HAVE_Trilinos
      INTEGER(SIK)  :: ierr

      IF(matrix%isInit) THEN
        IF(((j <= matrix%n) .AND. (i <= matrix%n)) &
          .AND. ((j > 0) .AND. (i > 0))) THEN
          IF(i==matrix%currow) THEN
            matrix%ncol=matrix%ncol+1
            matrix%jloc(matrix%ncol)=j
            matrix%aloc(matrix%ncol)=setval
          ELSE
            IF(matrix%currow>0) CALL ForPETRA_MatSet(matrix%A,i,matrix%ncol,matrix%jloc,matrix%aloc)
            matrix%aloc=0.0_SRK
            matrix%jloc=0
            matrix%ncol=0
            matrix%currow=i
          ENDIF
!TODO
          matrix%isAssembled=.FALSE.
        ENDIF
      ENDIF
#else
      CALL eMatrixType%raiseFatalError('Incorrect call to '// &
              modName//'::'//myName//' - Trilinos not enabled.  You will'// &
              'need to recompile with Trilinos enabled to use this feature.')
#endif
    ENDSUBROUTINE setShape_TrilinosMatrixType
!
!-------------------------------------------------------------------------------
!> @brief Gets the values in the Trilinos matrix - presently untested
!> @param declare the matrix type to act on
!> @param i the ith location in the matrix
!> @param j the jth location in the matrix
!>
!> This routine gets the values of the sparse matrix.  If the (i,j) location is
!> out of bounds, then -1051.0 (an arbitrarily chosen key) is returned.
!>
    SUBROUTINE get_TrilinosMatrixType(matrix,i,j,getval)
      CHARACTER(LEN=*),PARAMETER :: myName='get_TrilinosMatrixType'
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: matrix
      INTEGER(SIK),INTENT(IN) :: i
      INTEGER(SIK),INTENT(IN) :: j
      REAL(SRK),INTENT(INOUT) :: getval
#ifdef FUTILITY_HAVE_Trilinos
      INTEGER(SIK)  :: ierr

      getval=0.0_SRK
      IF(matrix%isInit) THEN
        ! assemble matrix if necessary
        IF (.NOT.(matrix%isAssembled)) CALL matrix%assemble()

        IF((i <= matrix%n) .AND. (j <= matrix%n) .AND. ((j > 0) .AND. (i > 0))) THEN
          CALL ForPETRA_MatGet(matrix%a,i,j,getval)
        ELSE
          getval=-1051._SRK
        ENDIF
      ENDIF
#else
      CALL eMatrixType%raiseFatalError('Incorrect call to '// &
              modName//'::'//myName//' - Trilinos not enabled.  You will'// &
              'need to recompile with Trilinos enabled to use this feature.')
#endif
    ENDSUBROUTINE get_TrilinosMatrixtype
!
!-------------------------------------------------------------------------------
    SUBROUTINE assemble_TrilinosMatrixType(thisMatrix,ierr)
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: thisMatrix
      INTEGER(SIK),INTENT(OUT),OPTIONAL :: ierr
      INTEGER(SIK) :: ierrc
#ifdef FUTILITY_HAVE_Trilinos
      INTEGER(SIK) :: iperr

      ierrc=0
      IF(.NOT.thisMatrix%isAssembled) THEN
        CALL ForPETRA_MatSet(thisMatrix%A,thisMatrix%currow,thisMatrix%ncol,thisMatrix%jloc,thisMatrix%aloc)
        thisMatrix%aloc=0.0_SRK
        thisMatrix%jloc=0
        thisMatrix%ncol=0
        thisMatrix%currow=0
        CALL ForPETRA_MatAssemble(thisMatrix%A)
        thisMatrix%isAssembled=.TRUE.
        ierrc=0
      ENDIF
      IF(PRESENT(ierr)) ierr=ierrc
#else
      CHARACTER(LEN=*),PARAMETER :: myName='assemble_TrilinosMatrixType'
      IF(PRESENT(ierr)) ierr=-1
      CALL eMatrixType%raiseFatalError('Incorrect call to '// &
         modName//'::'//myName//' - Trilinos not enabled.  You will'// &
         'need to recompile with Trilinos enabled to use this feature.')
#endif
    ENDSUBROUTINE assemble_TrilinosMatrixType
!
!-------------------------------------------------------------------------------
!> @brief tranpose the matrix
!> @param matrix declare the matrix type to act on
!>
!>
    SUBROUTINE transpose_TrilinosMatrixType(matrix)
      CHARACTER(LEN=*),PARAMETER :: myName='transpose_TrilinosMatrixType'
      CLASS(TrilinosMatrixType),INTENT(INOUT) :: matrix
      CALL eMatrixType%raiseFatalError(modName//'::'//myName// &
        ' - routine is not implemented!')
    ENDSUBROUTINE transpose_TrilinosMatrixType


ENDMODULE MatrixTypes_Trilinos
