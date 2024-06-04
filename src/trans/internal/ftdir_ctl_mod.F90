! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_CTL_MOD


CONTAINS
SUBROUTINE FTDIR_CTL(BATCHID)


!**** *FTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTL(..)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     PGP     -  gridpoint array
!     KVSETUV - "B" set in spectral/fourier space for
!                u and v variables
!     KVSETSC - "B" set in spectral/fourier space for
!                scalar variables
!     KPTRGP  -  pointer array to fields in gridpoint space

!     Method.
!     -------

!     Externals.  TRGTOL      - transposition routine
!     ----------  FOURIER_OUT - copy fourier data to Fourier buffer
!                 FTDIR       - fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!USE TPM_DIM
!USE TPM_GEOMETRY
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC

USE TRGTOL_MOD      ,ONLY : TRGTOL
USE FOURIER_OUT_MOD ,ONLY : FOURIER_OUT
USE FTDIR_MOD       ,ONLY : FTDIR
USE DIR_TRANS_CTL_MOD
!

IMPLICIT NONE

! Dummy arguments


! Local variables
!REAL(KIND=JPRB) :: ZGTF(KF_FS,D%NLENGTF)
INTEGER, INTENT(IN) :: BATCHID
INTEGER(KIND=JPIM) :: IST,JGL,IGL
INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space
! Fourier transform

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

CALL GSTATS(1640,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,IGL)
DO JGL=IBEG,IEND,IINC
  IGL = JGL
  IF(KF_FS>0) THEN
    CALL FTDIR(BATCHID,IGL)
  ENDIF

! Save Fourier data in FOUBUF_IN

call gstats(811,0)
    CALL FOURIER_OUT(BATCHID,IGL)
call gstats(811,1)

ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1640,1)
CALL GSTATS(106,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTL
END MODULE FTDIR_CTL_MOD



