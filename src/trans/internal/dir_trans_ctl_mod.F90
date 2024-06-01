! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIR_TRANS_CTL_MOD

USE LINKED_LIST_M

INTEGER :: NBATCHES ! NUMBER OF BATCHES
integer :: batch_sz ! BATCH SIZE IN LEVELS
integer, parameter :: MAX_ACTIVE = 3 ! MAXIMUM NUMBER OF ACTIVE BATCHES
integer :: ncomm_started=1 ! NUMBER OF OUTSTANDING COMMUNICATIONS
TYPE(LINKED_LIST) :: ACTIVE_BATCHES ! LIST OF ACTIVES BATCHES
TYPE(LINKED_LIST) :: ACTIVE_COMMS   ! LIST OF OUTSTANDING COMMUNICATIONS
INTEGER, PARAMETER :: WAITING = 1
INTEGER, PARAMETER :: READY = 2
INTEGER, PARAMETER :: PENDING=3
INTEGER, PARAMETER :: MAX_COMMS = 4  ! MAXIMUM NUMBER OF OUTSTANDING COMMUNICATIONS
INTEGER :: NLEV=8 ! NUMBER OF VERTICAL LEVELS FOR THIS TASK
INTEGER :: MAXVAR=4 ! MAXIMUM NUMBER OF VARIABLES (ALLOWS FOR SURFACE VARIABLES AT LEVEL 1)
INTEGER, ALLOCATABLE :: NVAR(:) ! NUMBER OF VARIABLES FOR EACH BATCH (ALLOWS FOR SURFACE VARIABLES AT LEVEL 1)

type AB
   integer :: stage,status,id
!   TYPE(TCOMM), POINTER :: comm_dep
     class(*), pointer               :: comm_dep
  end type AB

TYPE SR  ! send/receive structure
   INTEGER, ALLOCATABLE :: RECV_REQS(:,:),SEND_REQS(:,:)
   REAL(JPRB), ALLOCATABLE :: RECVBUF(:,:)
   LOGICAL INUSE=.FALSE.
   
   SUBROUTINE INIT
     IF(.NOT. ALLOCATED (RECV_REQS)) THEN
        ALLOCATE(RECV_REQS(MAX_RECV_NEIGHBORS,MAXVAR))
     ENDIF
     IF(.NOT. ALLOCATED (SEND_REQS)) THEN
        ALLOCATE(SEND_REQS(MAX_SEND_NEIGHBORS,MAXVAR))
     ENDIF
     if(.NOT. ALLOCATED (RECVBUF)) THEN
        ALLOCATE(RECVBUF(MAX_TOT,MAXVAR))
     ENDIF
   END SUBROUTINE INIT

   SUBROUTINE DELETE
     DEALLOCATE(RECV_REQS,SEND_REQS)
   END SUBROUTINE DELETE
END TYPE SR

TYPE(SR), ALLOCATABLE :: REQS ! ARRAY OF SEND/RECV REQUEST GROUPS, ONE GROUP PER LEVEL

type TCOMM
   INTEGER, POINTER :: RECV_REQS,SEND_REQS
   INTEGER :: nNEIGHBORS,STATUS,NCOMM
   TYPE(AB), POINTER :: MYBATCH
   REAL(JPRB), POINTER :: RECVBUF
   
   SUBROUTINE INIT
     INTEGER I
     
     DO I=1,MAX_COMMS
        IF(.NOT. REQS(I)%INUSE) THEN
           REQS(I)%INUSE = .TRUE.
           RECV_REQS => REQS(I)%RECV_REQS
           SEND_REQS => REQS(I)%SEND_REQS
           RECVBUF => REQS(I)%RECVBUF
           EXIT
        ENDIF
     ENDDO

     NCOMM = I
     IF(NCOMM .GT. MAX_COMMS) THEN
        PRINT *,'ERROR IN TCOMM INIT: EXCEEDED MAX_COMMS'
     ENDIF
     
   END SUBROUTINE INIT

   SUBROUTINE TEST_COMM(FLAG)

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) :: FLAG

    CALL MPI_TESTALL(nNEIGHBORS*NVAR(MYBATCH%ID),RECV_REQS,FLAG,MPI_STATUSES_IGNORE)
    
  END SUBROUTINE TEST_COMM

  SUBROUTINE complete_comm(batch)   ! Complete the communication call, release requests

    IMPLICIT NONE
    TYPE(AB), INTENT(OUT), POINTER :: BATCH
    INTEGER IBATCH
    
    CALL MPI_WAITALL(nNEIGHBORS*NVAR(MYBATCH%ID),RECV_REQS,MPI_STATUSES_IGNORE)
    CALL MPI_WAITALL(nNEIGHBORS*NVAR(MYBATCH%ID),SEND_REQS,MPI_STATUSES_IGNORE)
    BATCH => MYBATCH

    REQS(NCOMM)%INUSE = .FALSE.
    CALL PROCESS_BUF(MYBATCH%STAGE,MYBATCH%ID,RECVBUF,NVAR(MYBATCH%ID))
    
  END SUBROUTINE complete_comm
end type TCOMM


CONTAINS

SUBROUTINE PROCESS_BUF(STAGE,BATCHID,RECVBUF,NVARS)

END SUBROUTINE PROCESS_BUF



! Append a new batch to the end of the linked list of active batches
  ! Initiate first transpose (trgtol)
  SUBROUTINE ACTIVATE(N)

    USE LINKED_LIST_M

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: N

    TYPE(ab), POINTER :: NEWBATCH
    TYPE(TCOMM), POINTER :: NEWCOMM
    
    ALLOCATE(NEWBATCH)
    NEWBATCH%STAGE = 1
    NEWBATCH%STATUS = WAITING
    NEWBATCH%ID = N
    
    ALLOCATE(NEWCOMM)
    NEWCOMM%NNEIGHBORS = NEI1
    NEWCOMM%BATCH => NEWBATCH
    CALL NEWCOMM%INIT
    CALL ACTIVE_COMMS%APPEND(NEWCOMM)
    NEWBATCH%COMM_DEP => NEWCOMM
    
    CALL TRGTOL(BATCH=N,comm=NEWCOMM,NEWCOMM%RECVBUF)
    IF(NCOMM_STARTED .LT. MAX_COMMS) THEN
       ncomm_started = ncomm_started+1
       CALL NEWCOMM%START
       NEWBATCH%STATUS = WAITING
    ELSE
       NEWBATCH%STATUS = PENDING
    ENDIF

    CALL ACTIVE_BATCHES%APPEND(NEWBATCH)

  END SUBROUTINE ACTIVATE

     subroutine execute(BATCH)

       USE LINKED_LIST_M

       IMPLICIT NONE

       TYPE(AB), TARGET(INOUT) :: BATCH
       TYPE(TCOMM), POINTER :: COMM
       
       IF(BATCH%STATUS .NE. READY) THEN
          PRINT *,'ERROR IN EXECUTE: THIS BATCHIS NOT READY'
       ENDIF

       IF(BATCH%COMM_DEP .NE. NULL) THEN
          PRINT *,'EXECUTE FROM BATCH ',ID',: COMMUNICATION DEPENDENCY NOT MET'
       ENDIF

       IF(BATCH%STAGE .EQ. 1) THEN ! DO FFT

          CALL FTDIR(BATCH%ID,)
          CALL FOURIER_OUT(BATCH%ID)
          
          ALLOCATE(COMM)
          COMM%NNEIGHBORS = 1  ! BECAUSE IT IS ONE ALLTOALLV
          COMM%BATCH => BATCH
          CALL COMM%INIT
          CALL ACTIVE_COMMS%APPEND(COMM)
          BATCH%COMM_DEP => COMM
          CALL COMM_START_ATAV(COMM,BATCH)
          
       ELSEIF(STAGE .EQ. 2) THEN ! dO LEGENDRE

          CALL LTDIR(BATCH%ID)
          
       ENDIF

       STAGE = STAGE+1


     end subroutine execute

  SUBROUTINE DIR_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)

!**** *DIR_TRANS_CTL* - Control routine for direct spectral transform.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity
!     PSPDIV(:,:)  - spectral divergence
!     PSPSCALAR(:,:) - spectral scalarvalued fields
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     PGP(:,:,:)  - gridpoint fields

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTDIR_CTL   - control of Legendre transform
!                 FTDIR_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
!USE TPM_TRANS
!USE TPM_DISTR

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTDIR_CTL_MOD   ,ONLY : LTDIR_CTL
USE FTDIR_CTL_MOD   ,ONLY : FTDIR_CTL

USE LINKED_LIST_M
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC2(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPUV(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3A(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3B(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP2(:,:,:)

! Local variables

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_GPB

logical :: alldone,productive,flg,comm_compl
integer :: nactive,i,stage,status,c,NBATCHES,MAX_NEIGHBORS,NDONE
type(ab), pointer :: batch,IB
type(tcomm), pointer :: ic

!     ------------------------------------------------------------------

! Perform transform

IF_GPB = 2*KF_UV_G+KF_SCALARS_G

BATCH_SZ = 1
NLEV = 
NBATCHES = NLEV/BATCH_SZ
MAX_NEIGHBORS =

ALLOCATE(NVAR(NBATCHES))
MAXVAR = 1
DO I=1,NBATCHES
   NVAR(I) =
   IF(MAXVAR .LT. NVAR(I)) THEN
      MAXVAR = NVAR(I)
   ENDIF
ENDDO

ALLOCATE(REQS(MAX_COMMS)
DO I=1,MAX_COMMS
   CALL REQS(I)%INIT
ENDDO

alldone = .false.
ncomm_started = 0
nactive = 1
NDONE = 0
call activate(1)   ! Start the first batch, post communication
!active_batches(1)%stage = 1

do while(NDONE .LT. NBATCHES)

   comm_compl = .false.
   productive = .false.
   IC = ACTIVE_COMMS%HEAD
   DO WHILE(IC .NE. NULL)
      call IC%test_comm(flg)  ! Test/progress
      if(flg) then  ! If completed
         call IC%complete_comm(batch)   ! Complete the communication call, release requests
         CALL ACTIVE_COMMS%REMOVE(IC)
         comm_compl = .true.
         NCOMM_STARTED = NCOMM_STARTED -1
         exit
      endif
      IC = IC%NEXT
   enddo
   
   if(comm_compl) then    ! If a communication group has completed

      do WHILE(IB .NE. NULL) ! CHECK IF THERE ARE BATCHES PENDING COMMUNICATION; IF SO, START COMM ON THE FREED UP CHANNEL
         if(IB%status .EQ. PENDING) then
            CALL IB%COMM_DEP%START
            NCOMM_STARTED = NCOMM_STARTED +1
            EXIT
         ENDIF
      ENDDO

      productive = .true.
      call batch%execute  ! execute the next step in the algorithm for the batch that had communication completed, incl. start next communication if needed
      if(stage .eq. stage_final) then   ! If we're not done with this batch, post next communication
         call active_batches%remove(batch)
         NACTIVE = NACTIVE -1
         NDONE = NDONE+1
      endif
   else  ! if no communication is complete, see if we can do some computational work on current batches, if not we can initiate a new batch

      IB = ACTIVE_BATCHES%HEAD
      do WHILE(IB .NE. NULL)
         if(IB%status .ne. WAITING) then
            productive = .true.
            call IB%execute   ! Execute the next stage. This advances the stage marker to the next stage in the algorithm
            if(stage .EQ. stage_final) then
               call ACTIVE_BATCHES%REMOVE(IB)
               NACTIVE = NACTIVE -1
               NDONE = NDONE+1
               if(ACTIVE_BATCHES%HEAD .EQ. NULL) then
                  exit
               endif
            endif
         endif
         IB = IB%NEXT
      enddo
   
      if(.not. productive) then   ! If everyone is waiting, start a new batch
         if(nactive .lt. MAX_ACTIVE) then
            nactive = nactive+1
            call activate(nactive)  ! Start a new batch
         endif
      endif

   endif
   
enddo
   

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTL
END MODULE DIR_TRANS_CTL_MOD

