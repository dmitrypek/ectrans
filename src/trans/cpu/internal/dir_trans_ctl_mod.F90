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
USE PROGRESS_THREAD
use mpi

IMPLICIT NONE
  
!REAL(KIND=JPRB), ALLOCATABLE :: BGTF(:,:),BIN(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IREQ_RECV(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IREQ_SEND(:)
LOGICAL :: FIRST_PASS
REAL(KIND=JPRB), ALLOCATABLE :: ZCOMBUFR(:,:),ZCOMBUFS(:,:)
INTEGER(KIND=JPIM) :: SEND_ID(2),RECV_ID(2)

CONTAINS
SUBROUTINE DIR_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2,LUSE_WAITANY)

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
USE TPM_DISTR,       ONLY: D, NPROC, NPRTRNS
USE TPM_TRANS,       ONLY: NGPBLKS
USE TRGTOL_MOD,      ONLY: TRGTOL_PROLOG!

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
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LUSE_WAITANY

! Local variables

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV
INTEGER(KIND=JPIM) :: NDONE, KNRECV, KRECVCOUNT, KSENDCOUNT, KNSEND,KSENDCOUNT_GLOB
INTEGER(KIND=JPIM), ALLOCATABLE :: KSENDTOT(:), KRECVTOT(:), KSEND(:), KRECV(:), KNDOFF(:)
INTEGER(KIND=JPIM) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM) :: IVSET(KF_GP),IST
INTEGER(KIND=JPIM) :: KINDEX(D%NLENGTF)  

INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3

!     ------------------------------------------------------------------
!print *,'kf_uv_g,kf_scalars_g,kf_gp=',kf_uv_g,kf_scalars_g,kf_gp
!ivset(412) = -1
! Perform transform

IF(NPROMATR > 0 .AND. KF_GP > NPROMATR) THEN

  ! Fields to be split into packets

  CALL SHUFFLE(KF_UV_G,KF_SCALARS_G,ISHFUV_G,IVSETUV,ISHFSC_G,IVSETSC,&
 & KVSETUV,KVSETSC)

  IBLKS=(KF_GP-1)/NPROMATR+1

  DO JBLK=1,IBLKS
  
    CALL FIELD_SPLIT(JBLK,KF_GP,KF_UV_G,IVSETUV,IVSETSC,&
     & ISTUV_G,IENUV_G,IF_UV_G,ISTSC_G,IENSC_G,IF_SCALARS_G,&
     & ISTUV,IENUV,IF_UV,ISTSC,IENSC,IF_SCALARS)

    IF_FS = 2*IF_UV + IF_SCALARS
    IF_GP = 2*IF_UV_G+IF_SCALARS_G
    DO JFLD=1,IF_UV_G
      IPTRGP(JFLD) = ISHFUV_G(ISTUV_G+JFLD-1)
      IPTRGP(JFLD+IF_UV_G) = KF_UV_G+ISHFUV_G(ISTUV_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_SCALARS_G
      IPTRGP(JFLD+2*IF_UV_G) = 2*KF_UV_G+ISHFSC_G(ISTSC_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_UV
      IPTRSPUV(JFLD) = ISTUV+JFLD-1
    ENDDO
    DO JFLD=1,IF_SCALARS
      IPTRSPSC(JFLD) = ISTSC+JFLD-1
    ENDDO

    IF(IF_UV_G > 0 .AND. IF_SCALARS_G > 0) THEN
!      CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
!       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
!       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,PGP=PGP)
    ELSEIF(IF_UV_G > 0) THEN
!      CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
!       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
!       & KPTRGP=IPTRGP,PGP=PGP)
    ELSEIF(IF_SCALARS_G > 0) THEN
!      CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
!       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,PGP=PGP)
    ENDIF
!    CALL LTDIR_CTL(IF_FS,IF_UV,IF_SCALARS, &
!     & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
!     & KFLDPTRUV=IPTRSPUV,KFLDPTRSC=IPTRSPSC)
    
  ENDDO
ELSE

  ALLOCATE(KSENDTOT(NPROC))
  ALLOCATE(KRECVTOT(NPROC))
  ALLOCATE(KSEND(NPROC))
  ALLOCATE(KRECV(NPROC))
  ALLOCATE(KNDOFF(NPROC))

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IVSETSC(:) = -1
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF(PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF(PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    DO J3=1,UBOUND(PGP3A,3)
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF(PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    DO J3=1,UBOUND(PGP3B,3)
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
ENDIF

! Create a combined V-set array
  IST = 1
  IF (KF_UV_G > 0) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
!  print *,'IST =',IST
  IF (KF_SCALARS_G > 0) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF

  ! Call TRGTOL_PROLOG on "global" parameters to determine sizes of communication buffers
  CALL TRGTOL_PROLOG(KF_FS, KF_GP, IVSET, KSENDCOUNT, KRECVCOUNT, KNSEND, KNRECV, KSENDTOT, &
    &                KRECVTOT, KSEND, KRECV, KINDEX, KNDOFF, KGPTRSEND)
  ! Allocate receive request handle array
!  KSENDCOUNT_GLOB = SUM(KSENDTOT)

  FIRST_PASS = .FALSE.
   IF(.NOT. ALLOCATED(ZCOMBUFS)) THEN
      FIRST_PASS = .TRUE.
      ALLOCATE(ZCOMBUFS(-1:KSENDCOUNT,KNSEND))
      ALLOCATE(ZCOMBUFR(-1:KRECVCOUNT,KNRECV))
      ALLOCATE(IREQ_SEND(KNSEND))
      ALLOCATE(IREQ_RECV(KNRECV))
      call init_recvs(ZCOMBUFR,KNRECV,KRECV,KRECVTOT,IREQ_RECV)
      call init_sends(ZCOMBUFS,KNSEND,KSEND,KSENDTOT,KF_FS*2,IREQ_SEND)
      RECV_ID(2) = SEND_ID(2)
   ENDIF
   
   ! No splitting of fields, transform done in one go

  CALL FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,ZCOMBUFS,ZCOMBUFR,SEND_ID,RECV_ID,&
   & KSENDCOUNT,KNSEND,KRECVCOUNT,KNRECV,KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2,LUSE_WAITANY=LUSE_WAITANY)

  CALL LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, SEND_ID, &
   &PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
   &PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2)

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTL

  SUBROUTINE INIT_SENDS( PCOMBUFS,KNSEND,KSEND,KSENDTOT,KFIELD,IREQ_SEND)

    USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
    USE TPM_DISTR       ,ONLY : D,NPRCIDS,MTAGGL,NPRTRW
    USE MPL_MODULE      ,ONLY : MPL_ALL_MS_COMM
    USE ISO_C_BINDING
    USE MPI
    USE PROGRESS_THREAD
    
    IMPLICIT NONE

    INTEGER(KIND=JPIM),        INTENT(INOUT) :: IREQ_SEND(:)
!    INTEGER(KIND=JPIM),        INTENT(INOUT)   :: SEND_ID(2)
    REAL(KIND=JPRB),           INTENT(INOUT) :: PCOMBUFS(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: KNSEND,KFIELD
    INTEGER(KIND=JPIM),  INTENT(IN) :: KSENDTOT(:), KSEND(:)

    INTEGER INS, ISEND,NREQ,STATUS,J,A2AREQ
    INTEGER(KIND=JPIM) :: IST,IEN,IERR,REQ_TMP(1),dest,FLG
    INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

    FLG = 0;
    
    DO J=1,NPRTRW
       ILENS(J) = D%NLTSGTB(J)*KFIELD
       IOFFS(J) = D%NSTAGT1B(D%MSTABF(J))*KFIELD
       ILENR(J) = D%NLTSFTB(J)*KFIELD
       IOFFR(J) = D%NSTAGT1B(J)*KFIELD
    ENDDO

!    IF(THIS%STAGE .EQ. 1) THEN
       NREQ = KNSEND
       DO INS=1,NREQ
          ISEND = KSEND(INS)
          dest = nprcids(isend) -1
          CALL MPI_SEND_INIT(PCOMBUFS(:,INS), KSENDTOT(ISEND),MPI_REAL, &
                    &                 dest, MTAGGL, MPI_COMM_WORLD, IREQ_SEND(INS),IERR)
       ENDDO
       CALL PT_REQSET_REGISTER(NREQ,IREQ_SEND,FLG,SEND_ID(1),STATUS)
       print *,'Registered send reqset ',send_id(1)
       IF(STATUS .EQ. MPI_ERR_ARG) THEN
          PRINT *,'Error in pt_reqset, INIT_SENDS'
       ENDIF
       
 !   ELSE
       NREQ = 1
       CALL MPIX_ALLTOALLV_INIT(FOUBUF_IN,ILENS,IOFFS,MPI_REAL,FOUBUF, &
            &                  ILENR,IOFFR, MPI_REAL, MPL_ALL_MS_COMM,MPI_INFO_NULL,A2AREQ,IERR )
       REQ_TMP(1) = A2AREQ
       CALL PT_REQSET_REGISTER(1,REQ_TMP,FLG,SEND_ID(2),STATUS)
 !      print *,'Registered send reqset ',send_id(2)
       IF(STATUS .EQ. MPI_ERR_ARG) THEN
          PRINT *,'Error in pt_reqset, INIT_SENDS'
       ENDIF
!       THIS%RECV_ID(2) = THIS%SEND_ID(2)
!    ENDIF


  END SUBROUTINE INIT_SENDS

    SUBROUTINE INIT_RECVS(PCOMBUFR,KNRECV,KRECV,KRECVTOT,IREQ_RECV)

    USE TPM_DISTR       ,ONLY : D,NPRCIDS,MTAGGL
    USE ISO_C_BINDING
    USE MPI
    USE PROGRESS_THREAD
    
    IMPLICIT NONE

    INTEGER(KIND=JPIM),        INTENT(INOUT) :: IREQ_RECV(:)
!    INTEGER(KIND=JPIM),        INTENT(INOUT)   :: RECV_ID
    REAL(KIND=JPRB),           INTENT(INOUT) :: PCOMBUFR(:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: KNRECV
    INTEGER(KIND=JPIM),  INTENT(IN) :: KRECVTOT(:), KRECV(:)

    INTEGER INR, IRECV,NREQ,STATUS
    INTEGER(KIND=JPIM) :: IST,IEN,IERR,SRC,FLG

    FLG = 8

    
!    IF(THIS%STAGE .EQ. 1) THEN
       NREQ = KNRECV
       DO INR=1,NREQ
          IRECV = KRECV(INR)
          src = nprcids(irecv) -1
          CALL MPI_RECV_INIT(PCOMBUFR(:,INR), &
                    &                 KRECVTOT(IRECV),MPI_REAL, &
                    &                 src, MTAGGL, MPI_COMM_WORLD, IREQ_RECV(INR),IERR)
       ENDDO
       CALL PT_REQSET_REGISTER(NREQ,IREQ_RECV,FLG,RECV_ID(1),STATUS)
!       print *,'Registered receive reqset ',recv_id(1)
       IF(STATUS .EQ. MPI_ERR_ARG) THEN
          PRINT *,'Error in pt_reqset, INIT_RECVS'
       ENDIF
!    ENDIF

  END SUBROUTINE INIT_RECVS

END MODULE DIR_TRANS_CTL_MOD
