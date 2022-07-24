!----RENAME-DIR.F90                              1/30/04 
!----updated it on 4-17-06
!
!     PROGRAM to insert user's directory into various places       
!
      PROGRAM RENAME_DIR
!
      IMPLICIT NONE
!
      INTEGER :: I, NE, NEND, NLAST, NLONG, PCOL, SCOL
!
      CHARACTER (1) :: CURDIR1(132)
      CHARACTER (1) :: LBLANK(140) = ' ' 
      CHARACTER (10) :: CDATE, CTIME
      CHARACTER (132) :: LINE, CURDIR 
      CHARACTER (140) :: LONGLINE, BLANK = ' '
!
!     EQUIVALENCE (BLANK, LBLANK)
      EQUIVALENCE (CURDIR, CURDIR1)
!
      call date_and_time (cdate, ctime)
      WRITE (*,10) CDATE(5:6), CDATE(7:8), CDATE(1:2), CDATE(3:4),&
     &          CTIME(1:2), CTIME(3:4), CTIME(5:6)
10    FORMAT (/' RENAME-DIR (version 3/18/03) for structure/density',&
     &        ' predictions (',2(A2,'-'),2A2,' at ',2(A2,':'),A2,')')

!----Get name of top directory
      OPEN(UNIT=11,FILE='dir',STATUS='OLD')
      READ(11,110)CURDIR
! 110 FORMAT(A132)
110   FORMAT (A)
      CLOSE (UNIT = 11)
      DO I=1,132
!----How many characters in name of directory
        IF( CURDIR1(I).EQ.' ' ) THEN
          NE=I
          GO TO 15
        END IF
      END DO 
   15 NE=NE-1
   20 FORMAT('setenv dir ',A)
      WRITE (*,2001) CURDIR
2001  FORMAT (/' The top program directory is ',A)
!
!----Changes to make_files.com
      OPEN(UNIT=12,FILE='UTILITIES/make-files.com',STATUS='OLD')
  999 READ(12,110,END=1003)LINE
      IF (INDEX(LINE,'PREDICTIONS') .NE. 0) THEN
         WRITE(13,20) CURDIR(1:NE)     ! new version on unit # 13
      ELSE
         NLONG = LEN_TRIM(LINE)
         WRITE(13,110)LINE(1:NLONG)
      END IF 
      GO TO 999
1003  CLOSE (UNIT = 12)
      CLOSE (UNIT = 13)
      PRINT *, '...changed make-files.com'
!
!
!----Changes to summarize_tab.com
      OPEN(UNIT=12,FILE='UTILITIES/summarize_tab.com',STATUS='OLD')
 1005 READ(12,110,END=1009)LINE
      IF(INDEX(LINE,'PREDICTIONS') .NE. 0) THEN 
         WRITE(14,20) CURDIR(1:NE)     ! new version on unit # 14
      ELSE
         NLONG = LEN_TRIM(LINE)
         WRITE(14,110)LINE(1:NLONG)
      ENDIF
      GO TO 1005
1009  CLOSE (UNIT = 12)
      CLOSE (UNIT = 14)
      PRINT *, '...changed summarize_tab.com'
!
!----Changes to make_files.f90
      OPEN(UNIT=12,FILE='UTILITIES/make_files.f90',STATUS='OLD')
      LONGLINE = BLANK
1010  LINE = BLANK(1:132)
      READ (12,110,END=1055) LINE
      PCOL = INDEX(LINE,'PREDICTIONS')     ! column # of the P 
      IF (PCOL .NE. 0) THEN
         SCOL = PCOL + 10                  ! column # of the S 
         LONGLINE(1:PCOL-1) = LINE(1:PCOL-1)
!----Figure out right most character
         DO NLAST=132,1,-1
            IF (LINE(NLAST:NLAST) .NE. ' ') GO TO 23
         ENDDO
23       LONGLINE(PCOL:PCOL+NE-1) = CURDIR(1:NE)
         NEND = NLAST - SCOL + 1
         LONGLINE(PCOL+NE:PCOL+NE+NEND) = LINE (SCOL+1:NLAST)
!----How long in the line.  Look for first right-most blank.
         NLONG = LEN_TRIM(LONGLINE)
         WRITE (15,88) LONGLINE(1:NLONG)  ! new version on unit # 15 
88       FORMAT (A)
87          FORMAT (A,'&')
89          FORMAT ('     &',A)
      ELSE
         NLONG = LEN_TRIM(LINE)
         WRITE(15,110)LINE(1:NLONG)
      END IF
      GO TO 1010
1055  CLOSE (UNIT = 12)
      CLOSE (UNIT = 15)
      PRINT *, '...changed make_files.f90'
!
!
!
!----Changes to volume-additivity-densities-4.f
      OPEN(UNIT=12,FILE='UTILITIES/volume-additivity-densities-4.f',&
     &     STATUS='OLD')
 1015 READ(12,110,END=1100)LINE
      IF (INDEX(LINE,'PREDICTIONS') .NE. 0) THEN
         WRITE(LONGLINE,30) CURDIR(1:NE)
   30    FORMAT(6X,'DIR=',1H',A,'/UTILITIES/table3',1H')
!----How long in the line.  Look for first right-most blank.
      NLONG = LEN_TRIM(LONGLINE)
         IF (NLONG .LE. 72) THEN
            WRITE (16,88) LONGLINE(1:NLONG)  ! new version on unit # 15 
         ELSE
            WRITE (16,87) LONGLINE(1:72)      ! first 72 characters
            WRITE (16,89) LONGLINE(73:NLONG)  ! 73rd character onward
         ENDIF
         GO TO 1015
      ELSE
         NLONG = LEN_TRIM(LINE)
         WRITE(16,110)LINE(1:NLONG)              ! new version on unit # 16
      END IF
      GO TO 1015
 1100 CLOSE (UNIT=12)
      CLOSE (UNIT=16)
      PRINT *, '...changed volume-additivity-densities-4.f'
!
!----Changes to atom_code-vol-additivity.f
      OPEN(UNIT=12,FILE='UTILITIES/atom_code-vol-additivity.f',&
     &     STATUS='OLD')
1070  READ (12,110,END=1200) LINE
      IF (INDEX(LINE,'PREDICTIONS') .NE. 0) THEN
         WRITE(LONGLINE,31) CURDIR(1:NE)
   31    FORMAT(6X,'DIR=',1H',A,'/UTILITIES/atom_code_volume.data',1H')
!----How long in the line.  Look for first right-most blank.
      NLONG = LEN_TRIM(LONGLINE)
         IF (NLONG .LE. 72) THEN
            WRITE (17,88) LONGLINE(1:NLONG)  ! new version on unit # 15 
         ELSE
            WRITE (17,87) LONGLINE(1:72)      ! first 72 characters
            WRITE (17,89) LONGLINE(73:NLONG)  ! 73rd character onward
         ENDIF
         GO TO 1070
      ELSE
         NLONG = LEN_TRIM(LINE)
         WRITE(17,110)LINE(1:NLONG)              ! new version on unit # 17
      END IF
      GO TO 1070
 1200 CLOSE (UNIT=12)
      CLOSE (UNIT=17)
      PRINT *, '...changed atom_code-vol-additivity.f'
!
!----Changes to resort-summarize.com
      OPEN(UNIT=12,FILE='UTILITIES/resort-summarize.com',STATUS='OLD')
 1095 READ(12,110,END=1099)LINE
      IF (INDEX(LINE,'PREDICTIONS') .NE. 0) THEN
         WRITE(18,20) CURDIR(1:NE)     ! new version on unit # 14
      ELSE
         NLONG = LEN_TRIM(LINE)
         WRITE(18,110)LINE(1:NLONG)
      ENDIF
      GO TO 1095
1099  CLOSE (UNIT = 12)
      CLOSE (UNIT = 18)
      PRINT *, '...changed resort-summarize.com'
!
!
      END PROGRAM RENAME_DIR
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
