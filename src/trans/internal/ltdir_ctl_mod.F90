! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_CTL_MOD
CONTAINS
SUBROUTINE LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC)

!**** *LTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : LALLOPERM
!USE TPM_DIM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

USE LTDIR_MOD       ,ONLY : LTDIR
USE TRLTOM_MOD      ,ONLY : TRLTOM
USE TPM_DISTR, ONLY : MYPROC
USE TPM_DIM         ,ONLY : R
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2
integer(KIND=JPIM) :: jmmaxtmp
integer ila,i,j,d1,d2
character(len=50) :: s1,s2,s3,str
!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IBLEN = D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF

CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
CALL GSTATS(153,1)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)

! Direct Legendre transform

CALL GSTATS(103,0)
ILED2 = 2*KF_FS
CALL GSTATS(1645,0)
IF(KF_FS>0) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
!   jmmaxtmp=min(D%NUMP,3)
!   DO JM=1,jmmaxtmp
   DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIR(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1645,1)

IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)
CALL GSTATS(103,1)

!#define DEBUG_COMM

#ifdef DEBUG_COMM

d1 = size(pspvor,1)
d2 = size(pspvor,2)
print *,'Array dimensions in ltdir_ctl:',d1,d2


s1 = 'new/pspvor.'
write(s2,9) myproc-1
9 format(i0)

s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')

!do i=1,R%NTMAX
!   do j=1,KF_UV
do i=1,d2
   do j=1,d1
      write(11,10) j,i,pspvor(j,i)
   enddo
enddo
10 format(i3,i8,' ',E18.10)

close(11)

s1 = 'new/pspdiv.'
write(s2,9) myproc-1
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')

do i=1,d2
   do j=1,d1
      write(11,10) j,i,pspdiv(j,i)
   enddo
enddo

close(11)

s1 = 'new/pspsc2.'
write(s2,9) myproc-1
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')

do i=1,d2
   write(11,12) i,pspsc2(1,i)
enddo
12 format(i8,E18.10)

close(11)

#endif
!     -----------------------------------------------------------------

END SUBROUTINE LTDIR_CTL
END MODULE LTDIR_CTL_MOD
