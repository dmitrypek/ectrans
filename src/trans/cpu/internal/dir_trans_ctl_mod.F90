! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIR_TRANS_CTL_MOD

USE PARKIND1,          ONLY: JPIM, JPRB
USE OVERLAP_TYPES_MOD, ONLY: BATCH, BATCHLIST, STAT_WAITING, STAT_PENDING

IMPLICIT NONE

TYPE(BATCHLIST) :: ACTIVE_BATCHES
REAL(KIND=JPRB), ALLOCATABLE :: BGTF(:,:),BIN(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IREQ_RECV(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZCOMBUFR(:,:),ZCOMBUFS(:,:)
INTEGER, PUBLIC :: MAX_COMM(2)
INTEGER, PUBLIC, PARAMETER :: MAX_ACTIVE_BATCHES = 5
INTEGER, PUBLIC, PARAMETER :: STAGE_FINAL = 3
INTEGER :: NCOMM_STARTED(2)
INTEGER :: LAST_SUBMITTED(2)

CONTAINS
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

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE TPM_GEN,         ONLY: NPROMATR, NOUT
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
USE LINKED_LIST_M,   ONLY: LINKEDLISTNODE
USE TPM_DISTR,       ONLY: D, NPROC, NPRTRNS
USE TPM_TRANS,       ONLY: NGPBLKS,FOUBUF,FOUBUF_IN
USE TRGTOL_MOD,      ONLY: TRGTOL_PROLOG
USE LTDIR_CTL_MOD, ONLY: LTDIR_CTL_COMP

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

TYPE(BATCH) :: NEW_BATCH
TYPE(BATCH), POINTER :: COMPLETE_COMM_BATCH
INTEGER(KIND=JPIM) :: IBLKS, JBLK
INTEGER(KIND=JPIM) :: NDONE, KNRECV, KRECVCOUNT, KSENDCOUNT, KNSEND
INTEGER(KIND=JPIM), ALLOCATABLE :: KSENDTOT(:), KRECVTOT(:), KSEND(:), KRECV(:), KNDOFF(:)
INTEGER(KIND=JPIM) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: KINDEX(D%NLENGTF)  
TYPE(LINKEDLISTNODE), POINTER :: IB
INTEGER(KIND=JPIM)  :: IOFFSEND, IOFFRECV, IOFFGTF
INTEGER(KIND=JPIM) :: IST,NACTIVE,IBLEN
LOGICAL :: PRODUCTIVE, COMM_COMPL
INTEGER(KIND=JPIM) :: KSENDCOUNT_GLOB

!     ------------------------------------------------------------------

! Perform transform

WRITE(NOUT,*) "KF_UV_G = ", KF_UV_G
WRITE(NOUT,*) "KF_SCALARS_G = ", KF_SCALARS_G
WRITE(NOUT,*) "KF_GP = ", KF_GP
WRITE(NOUT,*) "KF_FS = ", KF_FS
WRITE(NOUT,*) "KF_UV = ", KF_UV
WRITE(NOUT,*) "KF_SCALARS = ", KF_SCALARS

IF (NPROMATR > 0) THEN

  ! Determine number of batches
  IBLKS = (KF_GP - 1) / NPROMATR + 1

  ! ================================================================================================
  ! Initialise work arrays shared by all batches
  ! ================================================================================================

  ! Allocate grid-to-Fourier transform buffer
  IF (ALLOCATED(BGTF) .AND. SIZE(BGTF, 1) /= KF_FS .AND. SIZE(BGTF, 2) /= (D%NLENGTF)) THEN
    DEALLOCATE(BGTF)
  ENDIF
  IF (.NOT. ALLOCATED(BGTF)) THEN
    ALLOCATE(BGTF(KF_FS,D%NLENGTF))
    ALLOCATE(BIN(D%NLENGTF,KF_FS))
  ENDIF

  ALLOCATE(KSENDTOT(NPROC))
  ALLOCATE(KRECVTOT(NPROC))
  ALLOCATE(KSEND(NPROC))
  ALLOCATE(KRECV(NPROC))
  ALLOCATE(KNDOFF(NPROC))

  ! Create a combined V-set array
  IST = 1
  IF (KF_UV_G > 0) THEN
    IVSET(IST:IST+KF_UV_G-1) = KVSETUV(:)
    IST = IST+KF_UV_G
    IVSET(IST:IST+KF_UV_G-1) = KVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF (KF_SCALARS_G > 0) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = KVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF

  ! Call TRGTOL_PROLOG on "global" parameters to determine sizes of communication buffers
  CALL TRGTOL_PROLOG(KF_FS, KF_GP, IVSET, KSENDCOUNT, KRECVCOUNT, KNSEND, KNRECV, KSENDTOT, &
    &                KRECVTOT, KSEND, KRECV, KINDEX, KNDOFF, KGPTRSEND)

  ! Allocate receive request handle array
  KSENDCOUNT_GLOB = SUM(KSENDTOT)

  IF (ALLOCATED(IREQ_RECV) .AND. SIZE(IREQ_RECV,1) /= KNRECV .AND. SIZE(IREQ_RECV,2) /= IBLKS) THEN
    DEALLOCATE(IREQ_RECV)
  ENDIF
  IF (.NOT. ALLOCATED(IREQ_RECV)) THEN
    ALLOCATE(IREQ_RECV(KNRECV,IBLKS))
  ENDIF

  ! Allocate communication buffers
  IF (ALLOCATED(ZCOMBUFR) .AND. SIZE(ZCOMBUFR,2) /= KNRECV .AND. SIZE(ZCOMBUFR,1) /= KRECVCOUNT) THEN
    DEALLOCATE(ZCOMBUFR)
  ENDIF
  IF (.NOT. ALLOCATED(ZCOMBUFR)) THEN
    ALLOCATE(ZCOMBUFR(KRECVCOUNT,KNRECV))
  ENDIF
  IF (ALLOCATED(ZCOMBUFS) .AND. SIZE(ZCOMBUFS,2) /= KNSEND .AND. SIZE(ZCOMBUFS,1) /= KSENDCOUNT_GLOB) THEN
    DEALLOCATE(ZCOMBUFS)
  ENDIF
  IF (.NOT. ALLOCATED(ZCOMBUFS)) THEN
    ALLOCATE(ZCOMBUFS(KSENDCOUNT_GLOB,KNSEND))
  ENDIF

  IBLEN = D%NLENGT0B*2*KF_FS
  IF (ALLOCATED(FOUBUF)) THEN
    IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
      DEALLOCATE(FOUBUF)
      ALLOCATE(FOUBUF(MAX(1,IBLEN)))
    ENDIF
  ELSE
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
  IF (ALLOCATED(FOUBUF_IN)) THEN
    IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
      DEALLOCATE(FOUBUF_IN)
      ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
    ENDIF
  ELSE
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF

  ! ================================================================================================
  ! Begin overlap loop
  ! ================================================================================================

  IOFFSEND = 1
  IOFFRECV = 1
  IOFFGTF = 1

  LAST_SUBMITTED(1) = 0
  LAST_SUBMITTED(2) = 0

  MAX_COMM(1) = 1
  MAX_COMM(2) = 1
  JBLK = 1 ! This keeps track of the last activated batch
  NDONE = 0 ! This keeps track of the number of completed batches
  NCOMM_STARTED(1) = 0 ! This keeps track of the batches in an active communication
  NCOMM_STARTED(2) = 0 ! This keeps track of the batches in an active communication
  NACTIVE = 1 ! This keeps track of the overall number of active batches
  
  CALL ACTIVATE(JBLK, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, PGP, IOFFSEND, IOFFRECV, &
    &           IOFFGTF, KSENDCOUNT, KRECVCOUNT)

  DO WHILE (NDONE < IBLKS)
    COMM_COMPL = .FALSE.
    PRODUCTIVE = .FALSE.

    ! Check whether any active batches have a completed communication
    IB => ACTIVE_BATCHES%HEAD
    DO WHILE (ASSOCIATED(IB))
      SELECT TYPE (THISBATCH => IB%VALUE)
      TYPE IS (BATCH)
        IF (THISBATCH%STATUS == STAT_WAITING) THEN
          IF (THISBATCH%COMM_COMPLETE(IREQ_RECV(:,THISBATCH%NBLK))) THEN
            COMM_COMPL = .TRUE.
            NCOMM_STARTED(THISBATCH%STAGE) = NCOMM_STARTED(THISBATCH%STAGE) - 1
            COMPLETE_COMM_BATCH => THISBATCH ! Keep track of which batch's communication completed
            EXIT
          END IF
        ENDIF
        IB => IB%NEXT
      END SELECT
    END DO

    IF (COMM_COMPL) THEN
      ! Now that one comm has completed, we can start the comm on the next batch whose comm
      ! is pending, if there is one
      IB => ACTIVE_BATCHES%HEAD
      DO WHILE (ASSOCIATED(IB))
        SELECT TYPE (THISBATCH => IB%VALUE)
        TYPE IS (BATCH)
          IF (THISBATCH%STATUS == STAT_PENDING .AND. NCOMM_STARTED(THISBATCH%STAGE) < &
                &  MAX_COMM(THISBATCH%STAGE) .AND. THISBATCH%NBLK == &
                & LAST_SUBMITTED(THISBATCH%STAGE) + 1) THEN
            CALL THISBATCH%START_COMM(PGP, IREQ_RECV(:,THISBATCH%NBLK), BIN, ZCOMBUFS, ZCOMBUFR)
            NCOMM_STARTED(THISBATCH%STAGE) = NCOMM_STARTED(THISBATCH%STAGE) + 1
            LAST_SUBMITTED(THISBATCH%STAGE) = LAST_SUBMITTED(THISBATCH%STAGE) + 1
            EXIT
          END IF
          IB => IB%NEXT
        END SELECT
      END DO

      PRODUCTIVE = .TRUE.
      CALL COMPLETE_COMM_BATCH%EXECUTE(BIN,BGTF, IREQ_RECV(:,COMPLETE_COMM_BATCH%NBLK), ZCOMBUFR, &
                                     & PSPVOR, PSPDIV, PSPSCALAR)

      ! If this batch has finished, remove it from the active batches list
      IF (COMPLETE_COMM_BATCH%STAGE == STAGE_FINAL) THEN
        IB => ACTIVE_BATCHES%HEAD
        DO WHILE (ASSOCIATED(IB))
          SELECT TYPE (LISTBATCH => IB%VALUE)
          TYPE IS (BATCH)
            IF (LISTBATCH%NBLK == COMPLETE_COMM_BATCH%NBLK) THEN
              CALL ACTIVE_BATCHES%REMOVE(IB)
              EXIT
            END IF
            IB => IB%NEXT
          END SELECT
        END DO
        NACTIVE = NACTIVE - 1
        NDONE = NDONE + 1
      ELSEIF (NCOMM_STARTED(2) < MAX_COMM(2) .AND. COMPLETE_COMM_BATCH%NBLK == LAST_SUBMITTED(2) + 1) THEN
        ! IREQ_RECV is not used, so pass 1st element
        CALL COMPLETE_COMM_BATCH%START_COMM(PGP, IREQ_RECV(:,1), BIN, ZCOMBUFS, ZCOMBUFR)
        NCOMM_STARTED(2) = NCOMM_STARTED(2) + 1
        LAST_SUBMITTED(2) = LAST_SUBMITTED(2) + 1
      ENDIF
    ELSE
      IF (.NOT. PRODUCTIVE .AND. NACTIVE < MAX_ACTIVE_BATCHES .AND. JBLK < IBLKS) THEN
        JBLK = JBLK + 1
        NACTIVE = NACTIVE +1
        CALL ACTIVATE(JBLK, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, PGP, IOFFSEND, &
                 &    IOFFRECV, IOFFGTF, KSENDCOUNT, KRECVCOUNT)
      ENDIF
    ENDIF
  ENDDO
  
  CALL LTDIR_CTL_COMP(KF_FS, KF_UV, KF_SCALARS, &
       &     PSPVOR=PSPVOR, PSPDIV=PSPDIV, PSPSCALAR=PSPSCALAR)
  


ELSE

  ! No splitting of fields, transform done in one go
  CALL ABORT_TRANS("DIR_TRANS_CTL: NPROMATR = 0 feature disabled for overlap version")

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTL

SUBROUTINE ACTIVATE(N, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, PGP, IOFFSEND, IOFFRECV, &
  &                 IOFFGTF, SENDCNTMAX, RECVCNTMAX)

  INTEGER,                      INTENT(IN)    :: N
  INTEGER(KIND=JPIM),           INTENT(IN)    :: KF_GP
  INTEGER(KIND=JPIM),           INTENT(IN)    :: KF_SCALARS_G
  INTEGER(KIND=JPIM),           INTENT(IN)    :: KF_UV_G
  INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KVSETUV(:)
  INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KVSETSC(:)
  REAL(KIND=JPRB),    OPTIONAL, INTENT(IN)    :: PGP(:,:,:)
  INTEGER(KIND=JPIM),           INTENT(INOUT) :: IOFFSEND
  INTEGER(KIND=JPIM),           INTENT(INOUT) :: IOFFRECV
  INTEGER(KIND=JPIM),           INTENT(INOUT) :: IOFFGTF
  INTEGER(KIND=JPIM),           INTENT(IN)    :: SENDCNTMAX
  INTEGER(KIND=JPIM),           INTENT(IN)    :: RECVCNTMAX

  CLASS(BATCH), POINTER :: NEW_BATCH

  ! Add a new batch to the list
  CALL ACTIVE_BATCHES%APPEND(BATCH(N, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, IOFFSEND, &
       &                     IOFFRECV, IOFFGTF, SENDCNTMAX, RECVCNTMAX))

  SELECT TYPE (NEW_BATCH => ACTIVE_BATCHES%TAIL%VALUE)
  TYPE IS (BATCH)
    IF (NCOMM_STARTED(1) < MAX_COMM(1)) THEN
      CALL NEW_BATCH%START_COMM(PGP, IREQ_RECV(:,N), BIN, ZCOMBUFS, ZCOMBUFR)
      NCOMM_STARTED(1) = NCOMM_STARTED(1) + 1
      LAST_SUBMITTED(1) = LAST_SUBMITTED(1) + 1
    ENDIF
  END SELECT

END SUBROUTINE ACTIVATE

END MODULE DIR_TRANS_CTL_MOD
