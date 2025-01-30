! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SPNORMD_MOD
CONTAINS
SUBROUTINE SPNORMD(PSPEC,KFLD,PMET,PSM,sout)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D,myproc
!

IMPLICIT NONE

REAL(KIND=JPRB)    ,INTENT(IN)  :: PSPEC(:,:)
REAL(KIND=JPRB)    ,INTENT(IN)  :: PMET(0:R%NSMAX)
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KFLD
REAL(KIND=JPRB)    ,INTENT(OUT) :: PSM(:,:)
character(*), INTENT(IN) :: sout

INTEGER(KIND=JPIM) :: JM ,JFLD ,JN ,IM ,ISP
character(len=50) :: s1,s2,s3,str

!     ------------------------------------------------------------------

s1 = 'norm.'
write(s2,9) myproc-1
9 format('.',i0)
s3 = trim(s2)
str = trim(s1) // sout // s3
open(11,file=str,form='formatted',status='unknown',action='write')

CALL GSTATS(1651,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1)  PRIVATE(JM,IM,JN,ISP,JFLD)
DO JM=1,D%NUMP
  PSM(:,JM) = 0.0_JPRB
  IM = D%MYMS(JM)
  IF(IM == 0)THEN
    DO JN=0,R%NSMAX
      ISP = D%NASM0(0)+JN*2
      DO JFLD=1,KFLD
        PSM(JFLD,JM) = PSM(JFLD,JM)+PMET(JN)*PSPEC(JFLD,ISP)**2
        write(11,10) jfld,jm,psm(jfld,jm)
      ENDDO
    ENDDO
  ELSE
    DO JN=IM,R%NSMAX
      ISP = D%NASM0(IM)+(JN-IM)*2
      DO JFLD=1,KFLD
        PSM(JFLD,JM) = PSM(JFLD,JM)+2.0_JPRB*PMET(JN)*&
         &(PSPEC(JFLD,ISP)**2+PSPEC(JFLD,ISP+1)**2)
        write(11,10) jfld,jm,psm(jfld,jm)
10 format(i3,i8,' ',E18.10)        
     ENDDO
    ENDDO
  ENDIF
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1651,1)

close(11)

!     ------------------------------------------------------------------

END SUBROUTINE SPNORMD
END MODULE SPNORMD_MOD





