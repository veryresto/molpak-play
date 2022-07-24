 program table_2                        ! 2006
!
!----Read a ..tab file and produce files sorted on density (tab.den)
!     and energy (tab.eng).  Also, write new files with just the unique
!     solutions listed.
!
 implicit none

 integer, parameter :: nset=10000

 character(110) :: line(nset),linea(nset)=' '
 character(110) :: linei, linek, line8(8), title (3)
 character (24) :: line24
 character(4) :: blank
 character(1) :: star(nset)

 integer :: ilne(nset),knum(nset),lnum(nset), nelim(nset), numln(nset)
 integer :: i,j,k,l,m,na,nc, ntot_elim, kk, kn, km, kl, ij, ii

 real :: den(nset), eng(nset), values1(9), values2(9)
!
!----climits shows the comparison defaults...
!                                              --unit cell--
!                        V     density    E    length  angle
 real :: climit(5) = (/0.11,   0.0011,  0.011, 0.011,  0.02/) 

 interface
   subroutine exachange(a,b)
     implicit none
     character(1) a,b
   end subroutine exachange

   subroutine exchange(x,y)
     implicit none
     real :: x,y
   end subroutine exchange

   subroutine exichange(m,n)
     implicit none
     integer :: m,n
   end subroutine exichange
 end interface
!
!----tab.den file will be on unit 16
!
!----Start reading the ...tab file on unit 15.  Initially the top half.
      DO I=1,8                 ! heading lines
         READ(15,11) LINE8(I)
 11         FORMAT(A110)
         WRITE(16,11) LINE8(I)
      ENDDO
!----Start reading the density and energy lines
      NA = 0
9     READ (15,10) LINE24
      IF (LINE24(1:4) == '    ') GO TO 15  ! blank signals end of this section
      NA = NA + 1
      READ(LINE24,24) ILNE(NA), DEN(NA), ENG(NA), STAR(NA) ! pull out good soln #,
24       FORMAT(I4,F10.3,F9.2,A1)                          ! density, energy and flag
!     print 24, ILNE(NA), DEN(NA), ENG(NA), STAR(NA)
      GO TO 9
!15    if (NA.GE.1000) NA=999
!15    DO I =1, NA
!      print 24, ILNE(I), DEN(I), ENG(I), STAR(I)
!     END DO
15    DO 16 I=1,NA-1
         DO 17 J=I+1,NA
             IF (DEN(I) .GE. DEN(J)) GO TO 17 
             CALL EXICHANGE(ILNE(I),ILNE(J))
             CALL EXCHANGE(DEN(I),DEN(J))
             CALL EXCHANGE(ENG(I),ENG(J))
             CALL EXACHANGE(STAR(I),STAR(J))
17       CONTINUE   
16    CONTINUE   
!     PRINT '(I4)',NA
      DO I=1,NA              ! write sorted densities
         WRITE(16,24) ILNE(I), DEN(I), ENG(I), STAR(I)
!        print 1234,I,ILNE(I), DEN(I), ENG(I), STAR(I)
1234         format (I5,1X,I4,F10.3,F9.2,A1)
26          FORMAT(I4,F10.3,F9.2,A1)
      END DO
      WRITE(16,10)
10       FORMAT (A) 
      DO I=1,3             ! heading for the second part of tab.den
         READ(15,10) LINEI
         TITLE(I) = LINEI
         WRITE(16,10) LINEI
      END DO
!----Read NA lines from the second part
      DO I=1,NA
         READ (15,10) LINEK
         READ (LINEK,27) K
27       FORMAT (3X,I4)
         LINE(K) = LINEK 
      ENDDO
!----Write density sorted lines for second part
      DO I=1,NA
         K = ILNE(I)
         WRITE (16,10) LINE(K)
      ENDDO
!----Read stuff on the bottom...place in LINEA for later use in tab.eng
      NC = 0
30    NC = NC + 1
      READ (15,10,END=35) LINEA(NC)
      GO TO 30
 35   NC = NC - 1
      DO I=1,NC
         WRITE(16,10) LINEA(I)
      ENDDO
!
!----Pass over the sorted densities and elimiate duplicates.  Use the second 
!     list as the guide.
!     
      DO I=1,NA
         NELIM(I) = 0   ! initialize to OK
      ENDDO
      NTOT_ELIM = 0	   ! count of total eliminated
      DO 220 I=1,NA-1
         K = ILNE(I)       ! pick up the line number
         IF (NELIM(K) /= 0) GO TO 220   ! non-0 value indicates rejection
         READ (LINE(K),215) VALUES1     ! 9 values...Vf, Df, Ef, a, b, c,
215         FORMAT (32X,F7.1,7X,F7.3,F8.2,3F8.4,3F7.3)   ! alp, beta, gamma 
         DO 210 J=I+1,NA
            L = ILNE(J)    ! line number of comparison value
            IF (NELIM(L) /= 0) GO TO 210
            READ (LINE(L),215) VALUES2  ! 9 values from 2nd line
            DO M=1,3
               IF (ABS(VALUES1(M) - VALUES2(M)) > CLIMIT(M)) GO TO 210
               IF (ABS(VALUES1(M+4) - VALUES2(M+4)) > CLIMIT(4)) GO TO 210
               IF (ABS(VALUES1(M+6) - VALUES2(M+6)) > CLIMIT(5)) GO TO 210
            ENDDO
            NELIM(L) = -1        ! all below tolerance, eliminate Lth
            NTOT_ELIM = NTOT_ELIM + 1
210      CONTINUE
220   CONTINUE
!
!----Summarize shortened list
      WRITE (16,240) NA, CLIMIT, NTOT_ELIM, (NA - NTOT_ELIM)
240   FORMAT (/I4,' items in original density-sorted list have been compared and'/ &
     &         ' duplicates removed based on the following tolerances...'/ &
     &           '   unit cell volume =',F5.2,' Angs**3;'/ &
     &           '   density =',F6.3,' g/cc; lattice energy =',F5.2,' kcal/mol;'/ &
     &           '   unit cell lengths =',F6.3,' Angs; unit cell angles =',F5.2,' degs'/ &
     &         I4,' duplicates eliminated;',I4,' survivors follow...')
      WRITE(16,10)                   ! headings for the list
      WRITE (16,10) (TITLE(I),I=1,3)
      DO 230 I=1,NA                  ! now the shortened list
         K = ILNE(I)
         IF (NELIM(K) /= 0) GO TO 230
         WRITE (16,10) LINE(K)
230   CONTINUE   
      CLOSE (UNIT = 16)
!
!----Resort on energy
!
      OPEN (UNIT=17, FILE='tab.eng', STATUS='UNKNOWN')
      DO 46 I=1,NA-1
         DO 47 J=I+1, NA
            IF (ENG(I) .LE. ENG(J)) GO TO 47
               CALL EXICHANGE (ILNE(I), ILNE(J))
               CALL EXCHANGE  (DEN(I), DEN(J))
               CALL EXCHANGE  (ENG(I), ENG(J))
               CALL EXACHANGE (STAR(I), STAR(J))
47       CONTINUE
46    CONTINUE   
!     re-order with the number of line  9-27-06
      KN = 1
      KK = 1
      DO I = 1,NA
       if (ENG(I+1) .EQ. ENG(I)) then
        KK = KK + 1
        numln(KN) = KK
!       print '(2I6)', KN,KK
       else
        KK = 1
        KN = KN + 1
       end if
      END DO
      KM = 0
      DO I = 1,KN-1
       if (numln(I) .EQ. 0) numln(I) = 1
       KM = KM + numln(I) 
      END DO   
!     print '(20I6)', (numln(I),I =1,KN-1)
      print '(2I6)', KN, KM
      II = 1
      DO 50 IJ = 1, KN-1
       if (IJ .EQ. 1) then
        KL = numln(IJ)
       else
        KL = numln(IJ) + II - 1
       end if
       print '(3I6,F8.2)', II,KL,KL-II+1,ENG(KL)
       DO 48 I=II, KL-1
         DO 49 J=I+1, KL
            IF (ILNE(I) .LT. ILNE(J)) GO TO 49
               CALL EXICHANGE (ILNE(I), ILNE(J))
               CALL EXCHANGE  (DEN(I), DEN(J))
               CALL EXCHANGE  (ENG(I), ENG(J))
               CALL EXACHANGE (STAR(I), STAR(J))
 49       CONTINUE
 48     CONTINUE
        II = 1 + KL
 50    CONTINUE
!     9-27-06
      DO I=1,8             ! write first 8 lines on tab.eng
         WRITE(17,11) LINE8(I)
      ENDDO
      DO I=1,NA            ! write sorted densities
         WRITE(17,26) ILNE(I), DEN(I), ENG(I), STAR(I)
      END DO
      WRITE (17,10)        ! heading for second part of tab.eng
      DO I=1,3             
         WRITE(17,10) TITLE(I)
      END DO
      DO I=1,NA            ! write energy sorted lines for second part
         K = ILNE(I)
         WRITE (17,10) LINE(K)
      ENDDO
      DO I=1,NC
         WRITE(17,10) LINEA(I)   ! write summary stuff on bottom
      ENDDO
!
!----Pass over the sorted energies and elimiate duplicates.  Use the second 
!     list as the guide.
!     
      DO I=1,NA
         NELIM(I) = 0   ! initialize to OK
      ENDDO
      NTOT_ELIM = 0	   ! count of total eliminated
      DO 320 I=1,NA-1
         K = ILNE(I)       ! pick up the line number
         IF (NELIM(K) /= 0) GO TO 320   ! non-0 value indicates rejection
         READ (LINE(K),215) VALUES1     ! 9 values...Vf, Df, Ef, a, b, c,
         DO 310 J=I+1,NA
            L = ILNE(J)    ! line number of comparison value
            IF (NELIM(L) /= 0) GO TO 310
            READ (LINE(L),215) VALUES2  ! 9 values from 2nd line
            DO M=1,3
               IF (ABS(VALUES1(M) - VALUES2(M)) > CLIMIT(M)) GO TO 310
               IF (ABS(VALUES1(M+4) - VALUES2(M+4)) > CLIMIT(4)) GO TO 310
               IF (ABS(VALUES1(M+6) - VALUES2(M+6)) > CLIMIT(5)) GO TO 310
            ENDDO
            NELIM(L) = -1        ! all below tolerance, eliminate Lth
            NTOT_ELIM = NTOT_ELIM + 1
310      CONTINUE
320   CONTINUE
!
!----Summarize shortened list
      WRITE (17,340) NA, CLIMIT, NTOT_ELIM, (NA - NTOT_ELIM)
340   FORMAT (/I4,' items in original energy-sorted list have been compared and'/ &
     &         ' duplicates removed based on the following tolerances...'/ &
     &           '   unit cell volume =',F5.2,' Angs**3;'/ &
     &           '   density =',F6.3,' g/cc; lattice energy =',F5.2,' kcal/mol;'/ &
     &           '   unit cell lengths =',F6.3,' Angs; unit cell angles =',F5.2,' degs'/ &
     &         I4,' duplicates eliminated;',I4,' survivors follow...')
      WRITE(17,10)                   ! headings for the list
      WRITE (17,10) (TITLE(I),I=1,3)
      DO 330 I=1,NA                  ! now the shortened list
         K = ILNE(I)
         IF (NELIM(K) /= 0) GO TO 330
         WRITE (17,10) LINE(K)
330   CONTINUE  
! 
      END program table_2

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 subroutine exachange(a,b)
  implicit none
  character(1) a,b,c
    c = a
    a = b
    b = c
 return
 end
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 subroutine exchange(x,y)
 implicit none
 real :: x,y,t
   t = x
   x = y
   y = t
 return
 end
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 subroutine exichange(m,n)
 implicit none
 integer :: k,m,n
   k = m
   m = n
   n = k
 return
 end
