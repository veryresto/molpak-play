!-----Read a MOLPAK output file (MOLPAK angles and unit cell volumes), 20000 max
!     values at present, and prepare a unique minimum volume output file.
!     Substantial change made to accommodate and properly handle decimal
!     rotation angles and steps. Fort.11 (becomes volume.min file) now has a 
!     sorted list of all points (angles & V's) that survive culling and 
!     averaging from smallest to largest V.
!
      PROGRAM MIN_VOL_SEARCH          ! major update on 9/7/04 at 3:45 PM.
!                                       improved identical soln filter, 5/25/07
!
      IMPLICIT NONE
!
      LOGICAL COM
      CHARACTER(1) :: JUNK(80)
      INTEGER, PARAMETER :: NSET=60000         ! 4-3-07         
!
      REAL :: EANGLE(NSET,3), VOL(NSET), XAAV(3)
      REAL :: EANGLE_SAVE(NSET,3), VOL_SAVE(NSET)
      REAL :: STEP(3), S_ORDER(3,8), VOL_MIN, ANGLETOLERANCE=0.1
      INTEGER :: N_AVG(27), N_SOLN_SEQ(NSET), LOWEST
      INTEGER :: NSOL_SAVE(NSET), N_WRITTEN=0, N_OUT, N_ON_11
!
      INTEGER :: I, IDEL, IWHICH, J, K, L, M, N, IE, IK, II, JJ, NK
      INTEGER :: NELIM=0, NMAX_SOL, NSOL, NSUM 
      REAL :: ASUM1, ASUM2, ASUM3, ASUM4, ASUM5 
      REAL :: A, VAV, WT, WTSUM
!     
      DATA S_ORDER / +1.,+1.,+1.,  &
     &               -1.,+1.,+1.,  &
     &               +1.,-1.,+1.,  &
     &               +1.,+1.,-1.,  &
     &               -1.,-1.,+1.,  &
     &               -1.,+1.,-1.,  &
     &               +1.,-1.,-1.,  &
     &               -1.,-1.,-1./
!
      WRITE(16,707)
707   FORMAT(/5X,3('***  Begin MIN-VOL-SEARCH  ***')/)
      READ (15,*)                            ! skip first 3 lines
      READ (15,*)
      READ (15,*)
      READ (15,*) NMAX_SOL, (STEP(I),I=1,3)       ! good stuff on 4th line
      WRITE(16,368) NMAX_SOL, (STEP(I),I=1,3)
368   FORMAT(/' Maximum no. of solutions from MIN-VOL-SEARCH to be',  &
     & ' processed for 2nd MOLPAK & refinement runs = ',I4/  &
     & ' MOLPAK angle increments used in initial search =',3F5.1)
      PRINT 389
389   FORMAT (' Begin unique minimum volume search')
      N=0
      READ (10,12) JUNK                  ! angles and volumes on unit # 10
12    FORMAT (80A1)
      WRITE (11,12) JUNK
      WRITE (11,'(4X,5HSTEPS,3F8.1)') (STEP(I),I=1,3)
      WRITE(16,363)
363   FORMAT(/'  Input angles and unit cell volumes from 1st MOLPAK',  &
     &        ' run'//,  &
     &'            Angle1  Angle2  Angle3   Volume '/)
10    N = N + 1
      READ (10,*,END=100) K, (EANGLE(N,I),I=1,3), VOL(N)
      WRITE(16,2) N, (EANGLE(N,I),I=1,3), VOL(N)
!1     FORMAT (9X,3F8.1,F10.3)
2     FORMAT (1X,I8,3F8.1,F10.3)
      GO TO 10
100   N = N - 1
!----The following forces the 0,0,0 point to be included in the output
!    file (11) as a test for a known crystal structure
!    DO I=1,N                                                        ! DU 9-91
!      IF(IANGLE(I,1).EQ.0.AND.IANGLE(I,2).EQ.0.AND.                 ! DU 9-91
!   +  IANGLE(I,3).EQ.0.AND.VOL(I).GT.0.0) THEN                      ! DU 9-91
!      WRITE(11,333) I,(IANGLE(I,K),K=1,3),VOL(I)                    ! DU 9-91
!      GO TO 11
!      END IF                                                        ! DU 9-91
!    END DO                                                          ! DU 9-91
!
!----The A1,90,A3 and A2,-90,-A3 points are identical.  Search for
!    these possible pairs and eliminate the A2 = -90 one.
!    The 2 volumes must agree within 0.01 also.
!11    DO 19 I=1,N-1
!        IF (VOL(I) .GE. 9998.0) GO TO 19
!        IF (IANGLE(I,2) .NE. 90) GO TO 19
!        DO 18 J=I+1,N
!            IF (VOL(J) .GE. 9998.0) GO TO 18
!            IF (IANGLE(J,2) .NE. -90) GO TO 18
!            IF (IANGLE(I,1) .EQ. IANGLE(J,1) .AND.
!    X           IANGLE(I,3) .EQ. -IANGLE(J,3) .AND.
!    X           ABS(VOL(I) - VOL(J)) .LE. 0.01) THEN
!                VOL(J) = 9999.0
!               NELIM = NELIM + 1
!               GO TO 19
!           ENDIF
!18       CONTINUE
!19     CONTINUE
!
!----Look at the A2 = 0 deg points.  Those for which A1 + A3 = N and  
!    equivalents.  Eliminate all but the first occurrence.  Assumption is made
!    that the search range does NOT go out of the -90 to +90 deg range
      COM = .false.
      DO 98 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 98
         IF (EANGLE(I,2) .NE. 0.0) GO TO 98
         ASUM1 = EANGLE(I,1) + EANGLE(I,3)
         ASUM2 = ASUM1 + 180.
         ASUM3 = ASUM1 - 360.
         ASUM4 = ASUM2 - 360.
         ASUM5 = 360. + ASUM1
         DO 27 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 27
            IF (EANGLE(J,2) .NE. 0.0) GO TO 27
            A = EANGLE(J,1) + EANGLE(J,3)
            IF (.NOT. ((A .EQ. ASUM1) .OR. (A .EQ. ASUM2) .OR. &
     &                 (A .EQ. ASUM3) .OR. (A .EQ. ASUM4) .OR. &
     &                 (A .EQ. ASUM5))) GO TO 27
            IF (ABS(VOL(J) - VOL(I)) .GE. 0.005) GO TO 27    ! V's must agree to 0.005 
            IF (COM) GO TO 97 
            write (16,26) 
26          format (/'  Minimum volume filtering of A2 = 0 deg', &
     &               ' equivalent points, # 2 eliminated....')
            COM = .true.
97          write (16,95) (eangle(i,l),l=1,3), vol(i), &
     &                    (eangle(j,l),l=1,3), vol(j)
95          format (' # 1 =',3F6.1,F8.3,'; # 2 =',3F6.1,F8.3)
            NELIM = NELIM + 1
            VOL(J) = 9999.0
27       CONTINUE
98    CONTINUE 
!
!----Look at the A1+A1' = 180, A2=-A2' and A3+A3' = 180 deg points ! 5/4/07 
!     and eliminate all but the first occurrence. 
      COM = .false.
      DO 198 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 198
         DO 127 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 127
            IF (.not. (abs(eangle(i,2) + EANGLE(J,2)) .le. 0.05)) GO TO 127
            asum1 = abs(EANGLE(i,1) - EANGLE(J,1)) - 180.0
            asum3 = abs(eangle(i,3) - eangle(j,3)) - 180.0
            IF (.NOT. (abs(asum1) .le. 0.05)) go to 127   ! sums must agree to 0.05 
            if (.not. (abs(asum3) .le. 0.05)) go to 127
            IF (ABS(VOL(J) - VOL(I)) .GE. 0.005) GO TO 127    ! V's must agree to 0.005 
            IF (COM) GO TO 197 
            write (16,126) 
126         format (/'  Minimum volume filtering of |A1i+A1j| = 180,' &
     &               ' A2i = -A2j, |A3i+A3j| = 180 equivalent points, # 2 eliminated....')
            COM = .true.
197         write (16,95) (eangle(i,l),l=1,3), vol(i), &
     &                    (eangle(j,l),l=1,3), vol(j)
            NELIM = NELIM + 1
            VOL(J) = 9999.0
127       CONTINUE
198    CONTINUE                                                  ! 5/4/07 
!
!----Look for A1+A1' = 360, A2=-A2' and A3=-A3' and              ! 5/4/07 
!             A1+A1' = 360, A2 = A2' and A3 = A3'
!     and eliminate all but the first occurrence. 
      COM = .false.
      DO 298 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 298
         DO 227 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 227
            IF ((abs(eangle(i,2) + EANGLE(J,2)) .le. 0.05) .and. &      ! A2 = -A2'
                (abs(eangle(i,3) + EANGLE(J,3)) .le. 0.05)) GO TO 221   ! A3 = -A3'
            IF ((abs(eangle(i,2) - EANGLE(J,2)) .le. 0.05) .and. &      ! A2 = A2'
                (abs(eangle(i,3) - EANGLE(J,3)) .le. 0.05)) GO TO 221   ! A3 = A3'
            go to 227
221         asum1 = abs(EANGLE(i,1) - EANGLE(J,1)) - 360.0
            IF (.NOT. (abs(asum1) .le. 0.05)) go to 227   ! sums must agree to 0.05 
            IF (ABS(VOL(J) - VOL(I)) .GE. 0.005) GO TO 227    ! V's must agree to 0.005 
            IF (COM) GO TO 297 
            write (16,226)
226         format (/'  Minimum volume filtering of |A1i+A1j| = 360, A2i=-A2j,', &
                     ' A3i=-A3j or'/ &
                     '                              |A1i+A1j| = 360, A2i=A2j,', &
                     ' A3i=A3j equivalent points, # 2 eliminated....')
            COM = .true.
297         write (16,95) (eangle(i,l),l=1,3), vol(i), &
     &                    (eangle(j,l),l=1,3), vol(j)
            NELIM = NELIM + 1
            VOL(J) = 9999.0
227       CONTINUE
298    CONTINUE                                                  ! 5/4/07 
!
!----Look for A2+A2' = 360, A1=-A1' and A3=-A3' and              ! 5/4/07 
!             A2+A2' = 360, A1= A1' and A3= A3'
!     and eliminate all but the first occurrence. 
      COM = .false.
      DO 598 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 598
         DO 527 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 527
            IF ((abs(eangle(i,1) + EANGLE(J,1)) .le. 0.05) .and. &      ! A2 = -A2'
                (abs(eangle(i,3) + EANGLE(J,3)) .le. 0.05)) GO TO 521   ! A3 = -A3'
            IF ((abs(eangle(i,1) - EANGLE(J,1)) .le. 0.05) .and. &      ! A2 = A2'
                (abs(eangle(i,3) - EANGLE(J,3)) .le. 0.05)) GO TO 521   ! A3 = A3'
            go to 527
521         asum1 = abs(EANGLE(i,2) - EANGLE(J,2)) - 360.0
            IF (.NOT. (abs(asum1) .le. 0.05)) go to 527   ! sums must agree to 0.05 
            IF (ABS(VOL(J) - VOL(I)) .GE. 0.005) GO TO 527    ! V's must agree to 0.005 
            IF (COM) GO TO 597 
            write (16,526)
526         format (/'  Minimum volume filtering of |A2i+A2j| = 360, A1i=-A1j,', &
                     ' A3i=-A3j or'/ &
                     '                              |A2i+A2j| = 360, A1i=A1j,', &
                     ' A3i=A3j equivalent points, # 2 eliminated....')
            COM = .true.
597         write (16,95) (eangle(i,l),l=1,3), vol(i), &
     &                    (eangle(j,l),l=1,3), vol(j)
            NELIM = NELIM + 1
            VOL(J) = 9999.0
527       CONTINUE
598    CONTINUE                                                  ! 5/4/07 
!
!----Look for A3+A3' = 360, A1 =-A1' and A2 =-A2' and            ! 5/4/07 
!             A3+A3' = 360, A1 = A1' and A2 = A2'
!     and eliminate all but the first occurrence. 
      COM = .false.
      DO 398 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 398
         DO 327 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 327
            IF ((abs(eangle(i,2) + EANGLE(J,2)) .le. 0.05) .and. &      ! A2 = -A2'
                (abs(eangle(i,1) + EANGLE(J,1)) .le. 0.05)) GO TO 428   ! A3 = -A3'
            IF ((abs(eangle(i,2) - EANGLE(J,2)) .le. 0.05) .and. &      ! A2 = A2'
                (abs(eangle(i,1) - EANGLE(J,1)) .le. 0.05)) GO TO 428   ! A3 = A3'
            go to 327
428         asum1 = abs(EANGLE(i,3) - EANGLE(J,3)) - 360.0
            IF (.NOT. (abs(asum1) .le. 0.05)) go to 327   ! sums must agree to 0.05 
            IF (ABS(VOL(J) - VOL(I)) .GE. 0.005) GO TO 327    ! V's must agree to 0.005 
            IF (COM) GO TO 397 
            write (16,426)
426         format (/'  Minimum volume filtering of |A3i+A3j| = 360, A1i=-A1j,', &
                     ' A2i=-A2j or'/ &
                     '                              |A3i+A3j| = 360, A1i=A1j,', &
                     ' A2i=A2j equivalent points, # 2 eliminated....')
            COM = .true.
397         write (16,95) (eangle(i,l),l=1,3), vol(i), &
     &                    (eangle(j,l),l=1,3), vol(j)
            NELIM = NELIM + 1
            VOL(J) = 9999.0
327       CONTINUE
398    CONTINUE                                                  ! 5/4/07 
!
!----Look for A1 = A1' and A2 = -A2' and A3 + A3' = 180 or              ! 5/4/07 
!             A3 = A3' and A2 = -A2' and A1 + A1' = 180               ! 5/4/07 
!             and eliminate all but the first occurrence. 
      COM = .false.
      DO 698 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 698
         DO 627 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 627   ! A2 = -A2'  
            IF (.not. (abs(eangle(i,2) + EANGLE(J,2)) .le. 0.05)) go to 627  
            if (abs(eangle(i,1) - EANGLE(J,1)) .le. 0.05 .and. &
                abs(abs(eangle(i,3) - eangle(j,3)) - 180.0) .le. 0.05) go to 648                
            if (abs(eangle(i,3) - EANGLE(J,3)) .le. 0.05 .and. &
                abs(abs(eangle(i,1) - eangle(j,1)) - 180.0) .le. 0.05) go to 648                
            go to 627
648         IF (ABS(VOL(J) - VOL(I)) .GE. 0.005) GO TO 627    ! V's must agree to 0.005 
            IF (COM) GO TO 697 
            write (16,626)
626         format (/'  Minimum volume filtering of |A1i-A1j| = 180, A2i=-A2j,', &
                     ' A3i=A3j or'/ &
                     '                               A1i=A1j, A2i=-A2j, |A3i-A3j| = 180', &
                     ' equivalent points, # 2 eliminated....')
            COM = .true.
697         write (16,95) (eangle(i,l),l=1,3), vol(i), &
     &                    (eangle(j,l),l=1,3), vol(j)
            NELIM = NELIM + 1
            VOL(J) = 9999.0
627       CONTINUE
698    CONTINUE                                                  ! 5/4/07 
!
!----Show the filtered list.
      WRITE (16,130) (N - NELIM) 
130   FORMAT (/' Remaining',I5,' points after filtering...'/  &
     &'            Angle1  Angle2  Angle3   Volume '/)
      K = 0
      DO 235 I=1,N
         IF (VOL(I) .GE. 9998.0) GO TO 235
         K = K + 1
         WRITE (16,237) K, I, (EANGLE(I,J),J=1,3), VOL(I)
237      FORMAT (I4,'/',I4,3F8.1,F10.3)
235   CONTINUE
!
!----We also filter if 2 volumes agree within 0.001 and the 
!    angles are related as A1, A2, A3 and A1, A2, 180+A3 or
!    180+A1, A2, A3. If a match is found, check the IS_ORDER array so that
!    consistent points are eliminated (this is for averaging).
      COM = .false.
      DO 29 I=1,N-1
         IF (VOL(I) .GE. 9998.0) GO TO 29
         DO 28 J=I+1,N
            IF (VOL(J) .GE. 9998.0) GO TO 28
            IF(ABS(VOL(I) - VOL(J)) .LT. 0.001) THEN
!              DO 112 K=1,3
                  IF (EANGLE(I,2) .NE. EANGLE(J,2)) GO TO 28
                  IF (((EANGLE(I,1) .EQ. EANGLE(J,1)) .AND.  &
                       (ABS(ABS(EANGLE(I,3) - EANGLE(J,3)) - 180.0) .LT. AngleTolerance)) .OR. &
                      ((EANGLE(I,3) .EQ. EANGLE(J,3)) .AND.  &
                       (ABS(ABS(EANGLE(I,1) - EANGLE(J,1)) - 180.0) .LT. AngleTolerance)))     &
                          GO TO 331
!112           CONTINUE
               GO TO 28            ! failed...no match
!----NK .ge. 1.   Which point, I or J, goes?
331          DO 51 M=1,8           ! success...a match
             DO 50 L=1,3
               IF (SIGN(1.0,EANGLE(I,L)) .NE. S_ORDER(L,M)) GO TO 51
50           CONTINUE
             II = M
             GO TO 52
51          CONTINUE
52          DO 61 M=1,8
             DO 60 L=1,3
               IF (SIGN(1.0,EANGLE(J,L)) .NE. S_ORDER(L,M)) GO TO 61
60           CONTINUE
             JJ = M
             GO TO 62
61          CONTINUE
!----Eliminate the point (II or JJ) with the larger IS_ORDER index
62             IF (II .LT. JJ) THEN
                   IE = J
                   IK = I
               ELSE
                   IE = I
                   IK = J
               ENDIF
               IF (COM) GO TO 64
!----Eliminate point in which |A1| or |A3| > 90
               write (16,326) 
326            format (/'  Minimum volume filtering of (A1,A2,A3), ',  &
                  '(A1,A2,180+A3), (180+A1,A2,A3)'/   &
                  '   pairs of equivalent points, eliminate one...')
               COM = .true.
64               IDEL = IE                  ! initialize point IE to be removed 
                 IWHICH = 2                 ! point # 2
               IF (ABS(EANGLE(IK,1)) .GT. 90.0 .OR.  &
                   ABS(EANGLE(IK,3)) .GT. 90.0) THEN 
                      IDEL = IK                    ! remove point IK
                      IWHICH = 1                   ! point # 1
                   ENDIF
               write (16,68) (eangle(ik,l),l=1,3), vol(ik), &
                          (eangle(ie,l),l=1,3), vol(ie), IWHICH
68          format (' # 1 =',3F6.1,F8.3,'; # 2 =',3F6.1,F8.3,' remove #',I2)
               NELIM = NELIM + 1      
               VOL(IDEL) = 9999.0
            ENDIF
!----Stay in the J loop to look for more of these
28       CONTINUE
29    CONTINUE
!
!----Final filtered list.
      WRITE (16,730) (N - NELIM) 
730   FORMAT (/' Remaining',I5,' points after all filtering...'/ &
     &'            Angle1  Angle2  Angle3   Volume '/)
      K = 0
      DO 735 I=1,N             
         IF (VOL(I) .GE. 9998.0) GO TO 735
         K = K + 1
         WRITE (16,237) K, I, (EANGLE(I,J),J=1,3), VOL(I)
735   CONTINUE
      WRITE(16,364)
364   FORMAT(/10X,'Solutions........',  &
     &      '   Angle1 Angle2 Angle3   Volume'/)
!
!----Special situation if N = 1
       IF (N .GT. 1) GO TO 199
!         WRITE (11,333) N, (EANGLE(1,K),K=1,3), VOL(1)
          N_WRITTEN = N_WRITTEN + 1   ! temportary save
          NSOL_SAVE(N_WRITTEN) = N
          DO K=1,3
             EANGLE_SAVE(N_WRITTEN,K) = EANGLE(1,K)  
          ENDDO
          VOL_SAVE(N_WRITTEN) = VOL(I)
          VOL(I) = 9999.0
          NSOL = 1
          WRITE(16,322) NSOL,N,(EANGLE(1,K),K=1,3),VOL(1)
       GO TO 265 
!
!----Search for points in which all three angle differ by the
!    step increment (i.e., adjacent points)
199   NSOL = 0
      DO 30 I=1,N-1
         IF(VOL(I) .GE. 9998.0) GO TO 30
         NSUM = 1
         DO 39 J= I+1,N
            IF(VOL(J) .GE. 9998.0) GO TO 39
            NK = 0
            DO 111 K= 1,3
               IF(ABS(EANGLE(I,K) - EANGLE(J,K)) .LE. STEP(K))  &
     &            NK = NK + 1
111         CONTINUE
!----Collect adjacent points for averaging
         IF(NK .NE. 3) GO TO 39
            IF (NSUM .EQ. 1) N_AVG(1) = I
            NSUM = NSUM + 1
            N_AVG(NSUM) = J
39       CONTINUE
!----NSUM = 1 means no neighbors
         IF (NSUM .EQ. 1) GO TO 65
            NSOL = NSOL + 1
            DO 122 K=1,3
122            XAAV(K) = 0
!----Calc weighted average of volumes and angles...weight by
!    1/V**2
            VAV = 0.0
            WTSUM = 0.0
            DO 116 L=1,NSUM
               M = N_AVG(L)
               WT = (1.0/VOL(M))**2        ! wt = (1/volume**2)
               DO 115 K=1,3
115               XAAV(K) = XAAV(K) + WT*EANGLE(M,K)  ! wtd angle sums
                  VAV = VAV + 1.0/VOL(M)              ! wtd volume sum 
                  WTSUM = WTSUM + WT                  ! sum of wts
116         CONTINUE
            DO 133 K=1,3
133            XAAV(K) = XAAV(K)/WTSUM
            VAV = VAV/WTSUM
!----Flag the used solutions
            DO 135 L=1,NSUM
135            VOL(N_AVG(L)) = 9999.0
               WRITE (16,423) NSOL, (XAAV(K),K=1,3), VAV
423            FORMAT('...Possible solution #',I6,':', &
     &                3F7.1,F10.3)
421            FORMAT (27X,3F7.1,F10.3/)
               IF (NSUM .LE. 14) THEN
                  WRITE(16,321) (N_AVG(L),L=1,NSUM)
321               FORMAT('   Average between sets',14I5)
               ELSE
                  WRITE (16,1321) (N_AVG(L),L=1,NSUM)
1321              FORMAT('   Average between sets',14I5/ &
     &                   '     ..................',13I5)
               ENDIF
!              WRITE(11,333) NSOL, (XAAV(K),K=1,3), VAV
333            FORMAT(1X,I8,3F8.1,F10.3)
               N_WRITTEN = N_WRITTEN + 1   ! temportary save
               NSOL_SAVE(N_WRITTEN) = NSOL
               DO K=1,3
                  EANGLE_SAVE(N_WRITTEN,K) = XAAV(K)  
               ENDDO
               VOL_SAVE(N_WRITTEN) = VAV 
               GO TO 311
!----Section for no neighbors
65          NSOL = NSOL + 1
            WRITE(16,322) NSOL, (EANGLE(I,K),K=1,3), VOL(I), I
322         FORMAT('   Possible solution #',I6,':',3F7.1,F10.3/  &
     &             '     ...there are no neighbors around',I6)
!           WRITE(11,333) NSOL, (EANGLE(I,K),K=1,3), VOL(I)
            N_WRITTEN = N_WRITTEN + 1   ! temportary save
            NSOL_SAVE(N_WRITTEN) = NSOL
            DO K=1,3
               EANGLE_SAVE(N_WRITTEN,K) = EANGLE(I,K)  
            ENDDO
            VOL_SAVE(N_WRITTEN) = VOL(I)
            VOL(I) = 9999.0
311      continue
30    CONTINUE
!
!----Sort the NSOL solutions from smallest to largest volumes
!
265   DO 1249 K=1,NSOL
1210  DO 1235 J=1,NSOL
1222     IF (VOL_SAVE(J) < 0.0) GO TO 1235     ! has Kth volume already been fingered?
         LOWEST = J
         VOL_MIN = VOL_SAVE(J) 
         IF (K == NSOL) GO TO 1236             ! store away the last one if +, no need to test 
      DO 1234 I=J+1,NSOL
         IF (VOL_SAVE(I) < 0.0) GO TO 1234         ! - volume indicates point already done
         IF (VOL_SAVE(I) .GE. VOL_MIN) GO TO 1234  ! is VOL(I) smaller? 
            VOL_MIN = VOL_SAVE(I)                  ! yes
            LOWEST = I                             ! also save soln # 
1234  CONTINUE
1236  VOL_SAVE(LOWEST) = -VOL_SAVE(LOWEST)      ! make volume - as a flag that it is accounted for
      N_SOLN_SEQ(K) = LOWEST                    
      GO TO 1249
1235  CONTINUE 
1249  CONTINUE
!----Output the sorted solutions to fort.11
      N_ON_11 = 0
      DO I=1,NSOL
         N_OUT = N_SOLN_SEQ(I)       ! pick up the solution #
         IF (VOL_SAVE(N_OUT) .GE. 0.0) CYCLE
         N_ON_11 = N_ON_11 + 1
         WRITE (11,333) NSOL_SAVE(N_OUT), (EANGLE_SAVE(N_OUT,K),K=1,3), ABS(VOL_SAVE(N_OUT))
      ENDDO 
!
      PRINT 820, N, (N - NELIM), N_ON_11 
820   FORMAT (I5,' angle/volume points from MOLPAK search'/  &
     &        I5,' points after filtering'/  &
     &        I5,' points passed to next step for elaboration') 
      WRITE(16,808)
808   FORMAT(/5X,3('*** MIN-VOL-SEARCH complete ***')/)
      STOP
      END PROGRAM MIN_VOL_SEARCH
