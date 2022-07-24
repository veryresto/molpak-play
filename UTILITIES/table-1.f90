 program table_1              ! last update on 2/23/04
!                               last update on 5-4-2010

 implicit none

 integer, parameter :: single = kind(0.0e0)
 integer, parameter :: nset=30000

 character(1) :: blank,junk,cmode,comp_id(30),spgr(2),&
               & title(80),star
 character(4) :: id
 character(65) :: line(nset)
 character(60) :: first_line
 character(86) :: second_line

 integer :: mode(nset),np(nset),ns,ns2(nset)
 integer :: ii,i,iend,j,k,n,na,nelim,nelim2,n_short, n_solns,n_eng,n_mode 
 real :: angle(3,nset),cell(6),den(2,nset),&
                  eng(nset),fpangle(3,nset),vo(nset),vol(2,nset)
 real :: E_avg, E_lsq, E_sum                  ! 5-4-10
 character(8)  :: MODE_1                      ! 5-4-10
 character(86) :: sum_line                    ! 5-4-10
 character(9)  :: E_lsq_line                  ! 5-4-10

      OPEN (UNIT=60,FILE='volume.find',STATUS='UNKNOWN')  ! save angles
!
      NELIM = 0
      NA = 0
      IEND = 0
      STAR = ' '                      ! 9-29-97
      N_SHORT = 0
      N_SOLNS = 0
      N_Eng   = 0
!
!    calculate avg E                                   ! 5-4-10
     print '(31H Total #     SUM E        Avg E)'      ! 5-4-10
     E_sum = 0.0                                       ! 5-4-10
118  READ (13,115,END=119) sum_line                    ! 5-4-10 
115  FORMAT(A86)                                       ! 5-4-10 
     READ(sum_line,117) E_lsq_line , n_mode            ! 5-4-10
117  FORMAT(68X,A,6X,I3)                               ! 5-4-10
     if (E_lsq_line == '*********'.OR.E_lsq_line == '      NaN') go to 118                  ! 5-4-10
     READ (E_lsq_line, "(F9.2)") E_lsq                 ! 5-4-10
!    READ (MODE_1,"(5X,I3)") n_mode                    ! 5-4-10
     if ( (n_mode == 1 .OR. n_mode == 5) .AND. E_lsq < -10.0 .AND. E_lsq > -50.0) then      ! 5-4-10
       E_sum = E_sum + E_lsq                           ! 5-4-10
       n_eng = n_eng + 1                               ! 5-4-10
!      print '(2I6,F9.2, F15.5)', n_mode, n_eng, E_lsq, E_sum                 ! temp
     endif 
     GO TO 118                                         ! 5-4-10
119  E_avg = E_sum / n_eng                             ! 5-4-10
     print '(I6,F15.5,1X, F9.2)', n_eng , E_sum, E_avg         ! 5-4-10  
!----Read general title info from 'search.data' file
      READ (14,116) COMP_ID, SPGR, TITLE
116   FORMAT (30A1/2A1/80A1)
      DO 14 I=1,2*NSET
19       READ (15,39,END=999) FIRST_LINE                ! line # 1
39       FORMAT (A)
            ID = FIRST_LINE(5:8)
            IF (ID == 'MOLP' .OR. ID == 'ROTP') GO TO 18 ! look for proper 1st line 
               N_SOLNS = N_SOLNS + 1
               GO TO 19
18          READ (FIRST_LINE,20) (ANGLE(J,I),J=1,3), VO(I)  ! get angles
20          FORMAT (29X,3F6.1,F7.2)               ! 6-11-02
            WRITE(60,61) I, (ANGLE(J,I),J=1,3), VO(I) ! 60 = volume.find
61          FORMAT(1X,I8,3F8.1,F10.3)             ! 6-11-02
!----Line # 1 says "MOLP" or "ROTP" in columns 5-8; line # 2 must say
!     "WMIN" or "DMA" or "PMIN".  If not abort and search for a new "MOLP" line.
         READ (15,39,END=997) SECOND_LINE               ! line # 2
         ID = SECOND_LINE(5:8)                          ! get line ID
         IF (ID.EQ.'WMIN' .OR. ID.EQ.' DMA' .OR. ID.EQ.'PMIN') GO TO 201    !7-24-07
996           NELIM = NELIM + 1
!             NP(NELIM) = N_SOLNS
              NP(NELIM) = I                              ! 7-21-07
              DO 42 J=1,3
42               FPANGLE(J,NELIM) = ANGLE(J,I)
              ANGLE(1,I) = 99999.                     ! flag it
              IF (IEND .EQ. 1) GO TO 999   ! if MOLPAK is last line..qui
                 BACKSPACE(UNIT=15)        ! 5-15-06
                 GO TO 14
997           IEND = 1
              GO TO 996
!----If failure with refinement, write message & go to next solution 
 201        READ (SECOND_LINE,25,ERR=300) VOL(1,I), DEN(1,I)
            READ (15,39,ERR=300) SECOND_LINE               ! line # 3 4-24-96
            ID = SECOND_LINE(5:8)                          ! get line ID
            IF (.not. (ID.EQ.'WMIN' .OR. ID.EQ.' DMA' .OR. ID.EQ.'PMIN')) GO TO 300  !7-24-07
            READ (SECOND_LINE,26,ERR=300) VOL(2,I), DEN(2,I),&
                         ENG(I),MODE(I)   ! 5-17-96
!           if ( ENG(I) < (1.50 *E_avg) .OR. ENG(I) > 0.0) GO TO 300     ! 5-4-10 DU 
!           if ( ENG(I) > 0.0) GO TO 300      ! 7-14-08 DU
            if ( ENG(I) < -100.0 .OR. ENG(I) > 0.0) GO TO 300            ! 5-4-10 DU
!-----Look for reduced cell parameters, check that each of next 6 lines is as expected...
         READ (15,39,END=999) SECOND_LINE(1:10)          ! 1/6     
            IF (SECOND_LINE(1:6) /= ' T 1= ') GO TO 300
         READ (15,39,END=999) SECOND_LINE(1:10)          ! 2/6     
            IF (SECOND_LINE(1:6) /= ' T 2= ') GO TO 300
         READ (15,39,END=999) SECOND_LINE(1:10)          ! 3/6     
            IF (SECOND_LINE(1:10) /= ' T 2 INV= ') GO TO 300
         READ (15,39,END=999) SECOND_LINE(1:10)          ! 4/6     
            IF (SECOND_LINE(1:6) /= '+     ') GO TO 300
         READ (15,39,END=999) SECOND_LINE(1:10)          ! 5/6     
            IF (SECOND_LINE(1:6) /= '      ') GO TO 300
         READ (15,39,END=999) SECOND_LINE(1:10)          ! 6/6     
            IF (SECOND_LINE(1:9) /= ' CELL  1=') GO TO 300
         READ (15,39,END=999) SECOND_LINE(1:76) ! 7th line has reduced cell info
            IF (SECOND_LINE(1:9) /= ' CELL  2=') GO TO 300
         LINE(I) = SECOND_LINE(12:76)  
!----Save soln number to make a short table with final densities and energies
            N_SHORT = N_SHORT + 1
            NS2(N_SHORT) = I
         NA = NA + 1
         GO TO 14
 300     NELIM = NELIM + 1
         NP (NELIM) = I
         DO J=1,3
          FPANGLE(J,NELIM) = ANGLE(J,I)
         END DO
         ANGLE(1,I) = 99999.
14    CONTINUE
  25   FORMAT(36X,F8.2,8X,F8.3)                           ! 4-24-96
  26   FORMAT(36X,8X,F8.2,8X,F8.3,F9.2,6X,I3)             ! 5-17-96
  50   FORMAT(11x,A65)

!C----Create final summary tables

999   WRITE (16,1006) COMP_ID, SPGR, TITLE
1006  FORMAT (' Compound ID: ',30a1,' sp grp: ',2a1/&
     &        ' Title: ',80a1/)
      WRITE(16,45)
   45 FORMAT(16X,'Summary of MOLPAK/REFINEMENT/Cell Reduction Results'/&
             16X,'---------------------------------------------------'/)
!C----Make the short table first
      WRITE (16,950)
950   FORMAT ('   #',7x,'Df',7x,'Ef'/)
      NELIM2 = 0                              ! 4-19-06
      DO 960 I=1,N_SHORT
         K = NS2(I)
!---------------------------------------------! 4-19-06
         IF ( ENG(K) .GE. 0.0 ) STAR='*'      ! 9-29-97
!        IF ( ENG(K) .GE. 0.0 ) THEN
!         NELIM2 = NELIM2 + 1
!         NELIM = NELIM + 1
!         NP(NELIM) = K
!         DO II = 1,3
!          FPANGLE(II,NELIM) = ANGLE(II,K)
!         END DO
!         ANGLE(1,K) = 99999.0
!         GOTO 960
!        END IF                               ! 4-19-06
!---------------------------------------------! 4-19-06
         WRITE (16,955) K, DEN(2,K), ENG(K), STAR
955      FORMAT (I4,F10.3,F9.2,A1)
         STAR=' '                             ! 9-29-97
960   CONTINUE
      WRITE(16,30)
   30 FORMAT(/,7x,'.....MOLPAK......  ',&
                '............REFINEMENT............',&
                  '  ...........Reduced Cell Parameters.........'/&
             5X,'#',1X,'angl1 angl2 angl3 ',2X,'Vi',5X,'Vf',5X,&
             'Di',5X,'Df',6X,'Ef',&
             T61,'    a       b       c     alpha   beta  gamma'/)

      NS = 0
      DO 98 J=1,NA - NELIM2
34       NS = NS + 1
         IF (ANGLE(1,NS) .GE. 999.) GO TO 34        ! bad solution
         READ (LINE(NS),55,ERR=34) CELL
55       FORMAT (F8.4,2F10.4,1X,3F10.3)
!C        added it for printing information of MODE of refinement 5/16/96
         IF (MODE(NS).EQ.1) THEN
           CMODE='L'
         ELSE
           IF (MODE(NS).EQ.3) THEN
            CMODE='R'
           ELSE
            IF (MODE(NS).EQ.5) THEN
              CMODE='D'
            END IF
           END IF
         END IF
         WRITE(16,35) NS, (ANGLE(I,NS),I=1,3), (VOL(K,NS),K=1,2),&
     &               (DEN(N,NS),N=1,2), ENG(NS), CELL, CMODE
!  35    FORMAT(1x,'>>',I3,3F6.1,2F7.1,2F7.3,F8.2,3F8.4,3F7.3,&
   35    FORMAT(1x,'>>',I4,3F6.1,2F7.1,2F7.3,F8.2,3F8.4,3F7.3,&
     &          1X,A1)
98    CONTINUE
!C----Provide column information
      WRITE (16,200)
200   FORMAT (//' MOLPAK angles 1-3 are initial model rotation angles',&
     &          ' before refinement'/&
     &          ' Vi and Di are initial unit cell volume and crystal',&
     &          ' density before refinement'/&
     &          ' Vf and Df are final unit cell volume and crystal',&
     &          ' density after refinement'/&
     &          ' Ef is the final lattice energy in kcal/mol'/&
     &          ' the a, b, c, alpha, beta, gamma given are for the',&
     &          ' final reduced cell'/&
     &          ' R = refinement with Rosenbrock search '/&
     &          ' L = refinement with Least Squares '/&
     &          ' D = refinement with DMAREL program')

!C----Have any solutions been incomplete?
      IF (NELIM .EQ. 0) STOP
      WRITE (16,495) NELIM
495     FORMAT (/' Problems are found with',I4,' solution(s)...',&   ! 2/23/04
     &          'see master.log file for more information')          ! 2/23/04
      DO 500 I=1,NELIM
         WRITE (16,505) NP(I), (FPANGLE(J,I),J=1,3)                
505      FORMAT ('   #',I3,', MOLPAK angles =',2(F7.1,','),&
     &           F7.1)
500   CONTINUE
         write(16,'(20I4,2H \)') (NP(I),I=1,NELIM)                    ! 5-4-10
      STOP
      END PROGRAM TABLE_1
