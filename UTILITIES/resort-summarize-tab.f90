!     Take a densities.list file (created by ...com script) and 
!      produce two new files sorted on energy and density.       
!
      PROGRAM RESORT_SUMMARIZE_TAB          ! 5/31/04 
!
!     updated on 9-18-07 (Increased # geoms to 54)
!     Increased # geoms to 40                                ! 6-11-06 
      IMPLICIT NONE
!    
      INTEGER, PARAMETER :: N_VALUES = 54      ! 54 max geoms  7-18-07
      CHARACTER (130) :: LINE
      CHARACTER (45) :: LINEA(N_VALUES)
      CHARACTER (65) :: LINEB(N_VALUES)
      CHARACTER (15) :: LINEC(N_VALUES)
      CHARACTER (95) :: LINED(N_VALUES)
      CHARACTER (6) :: FTEMP
      CHARACTER (9) :: FTEMP2
      REAL :: ENG(N_VALUES), DEN(N_VALUES)
      INTEGER :: NA, I, J
!
      OPEN (UNIT=10, FILE='densities.list', STATUS='OLD')   ! densities.list
      OPEN (UNIT=11, FILE = 'E-temp', STATUS='UNKNOWN')     ! created by
      OPEN (UNIT=12, FILE = 'D-temp', STATUS='UNKNOWN')     ! ...com script
!----Find 1st valid header line id'd by ' Compound ID'  5/31/04
1     READ (10,5) LINE                                ! 5/31/04
         IF (LINE(1:12) .NE. ' Compound ID') GO TO 1  ! 5/31/04
         WRITE (11,5) LINE  ! write 1st line            5/31/04
         WRITE (12,5) LINE                            ! 5/31/04
      DO I=1,5              ! next 5 header lines       5/31/04
         READ(10,5) LINE
         WRITE(11,5) LINE
         WRITE(12,5) LINE
      END DO 
  5   FORMAT(A)
      NA = 1
 10   READ (10,5) LINE
         IF (LINE(1:5) .EQ. ' ----') GO TO 15      ! ---- means end of values
         READ (LINE,'(A45,F10.2,A65)') LINEA(NA),ENG(NA),LINEB(NA)
         READ (LINE,'(A15,F10.3,A95)') LINEC(NA),DEN(NA),LINED(NA)
         NA = NA + 1
      GO TO 10
 15   NA = NA - 1
      WRITE (11,5) LINE
      WRITE (12,5) LINE
      Loop_1 : DO I=1,NA-1             ! sort energies
         Loop_2 : DO J=I+1,NA 
                      IF (ENG(I) .GT. ENG(J)) THEN
                           CALL EXCHANGE  ( ENG(I) , ENG(J)  )
                           CALL EXACHANGE ( LINEA(I), LINEA(J) )  
                           CALL EXBCHANGE ( LINEB(I), LINEB(J) )
                      END IF
                  END DO Loop_2
               END DO Loop_1
      DO I=1,NA
          WRITE (11,25) LINEA(I),ENG(I),LINEB(I) 
 25       FORMAT (A45,F10.2,A65)
      END DO
      Loop_3 : DO I=1,NA-1                ! sort densities
         Loop_4 : DO J=I+1,NA
                      IF (DEN(I) .LT. DEN(J)) THEN
                           CALL EXCHANGE  ( DEN(I) , DEN(J)  )
                           CALL EXCCHANGE ( LINEC(I), LINEC(J) )
                           CALL EXDCHANGE ( LINED(I), LINED(J) )
                      END IF
                  END DO Loop_4
               END DO Loop_3
      DO I=1,NA
          WRITE (12,35) LINEC(I),DEN(I),LINED(I)
 35       FORMAT (A15,F10.3,A95)
      END DO
      WRITE (11,30)             ! write temp E file
      WRITE (12,30)             ! write temp density file
 30   FORMAT (1X,70('-'))
 20   READ (10,5,END=999) LINE
         WRITE (11,5) LINE
         WRITE (12,5) LINE
      GO TO 20
!----Remove all right-hand blanks in temp files and write
!     final energies and densities
999   CLOSE (UNIT = 11)
      CLOSE (UNIT = 12)      
      FTEMP = 'E-temp';  FTEMP2 = 'energies'
      CALL SHORTEN_LINES (FTEMP, FTEMP2)
      WRITE (*,*) ' Sorted file named energies created'
      FTEMP = 'D-temp';  FTEMP2 = 'densities'
      CALL SHORTEN_LINES (FTEMP, FTEMP2)   
      WRITE (*,*) ' Sorted file named densities created'
      CALL SYSTEM ('rm -f D-temp E-temp')  ! delete temp files
      STOP   
      END PROGRAM RESORT_SUMMARIZE_TAB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE EXCHANGE(X,Y)
      REAL :: X,Y
      T = X
      X = Y
      Y = T
      RETURN
      END SUBROUTINE EXCHANGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE EXACHANGE(A,B)
      CHARACTER (45) :: A, B, C
      C = A
      A = B
      B = C
      RETURN
      END SUBROUTINE EXACHANGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      SUBROUTINE EXBCHANGE(A,B)
      CHARACTER (65) :: A, B, C
       C = A
       A = B
       B = C
      RETURN
      END SUBROUTINE EXBCHANGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE EXCCHANGE(A,B)
      CHARACTER (15) :: A, B, C
      C = A
      A = B
      B = C
      RETURN
      END SUBROUTINE EXCCHANGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE EXDCHANGE(A,B)
      CHARACTER (95) A, B, C
      C = A
      A = B
      B = C
      RETURN
      END SUBROUTINE EXDCHANGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SHORTEN_LINES(FILE_IN, FILE_OUT)
!
!----Read the line of a file and remove all right 
!     hand blanks.
!
      CHARACTER (LEN=120) :: LINE, BLANK
      CHARACTER (LEN=1) :: BLANK1(120)
      CHARACTER (LEN=6) :: FILE_IN
      CHARACTER (LEN=9) :: FILE_OUT
      DATA BLANK1 /120*' '/
      EQUIVALENCE (BLANK, BLANK1)
!
      OPEN (UNIT = 20, FILE = FILE_IN, STATUS = 'OLD')      ! original
      OPEN (UNIT = 21, FILE = FILE_OUT, STATUS = 'UNKNOWN') ! shortened
!
100   LINE = BLANK
      READ (20,15,END=200) LINE
15    FORMAT (A)
      DO 190 I=120,1,-1                       ! identify the
         IF (LINE(I:I) .EQ. ' ') GO TO 190    ! right-most
         WRITE (21,15) LINE(1:I+1)            ! blank & add 1 more
         GO TO 100
190   CONTINUE
      WRITE (21,15) BLANK(1:10)   ! write out blank line
      GO TO 100                   ! with 10 characters
200   CLOSE (UNIT = 20)
      CLOSE (UNIT = 21)
      RETURN
      END SUBROUTINE SHORTEN_LINES
