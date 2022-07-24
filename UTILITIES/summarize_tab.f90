      PROGRAM SUMMARIZE_TAB                                       ! 2008 
!
!     update on 9-18-07 
!     Added geoms au as                                           ! 12-27-06
!     Added geoms ac, an, ao, bj, bk, fb and cf for molecules on  ! 1/26/04 
!        mirrors.                                                 ! 7/5/03 
!     Added geom al for molecule in P2/c and on center.           ! 12/7/03
!
!----Read the ....tab files for a compound and prepare a density,
!     energy and space group summary
!
!
!     CHARACTER INPUT*40, LINE*33, BLANK*33, FMT1*200, FMT2*80
!      CHARACTER LINE45*45, LINE70*70, STAR(NSG)*1, 
!    2          TABSTAR(MAX_ENTRIES)*1
!      CHARACTER GEOMS(NSG)*2, SPGRPS(NSG)*7, GTYPE*2, L_OUT(NSG)*90
!     CHARACTER LINE110*110, MISSING(NSG)*2, TYPE*2
!     BYTE I_HERE(NSG) 
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: NSG = 54                             ! 6-4-07 
      INTEGER, PARAMETER :: MAX_ENTRIES = 5000
      CHARACTER (40) :: INPUT
      CHARACTER (33) :: LINE, BLANK
      CHARACTER (200) :: FMT1
      CHARACTER (80) :: FMT2 
      CHARACTER (45) :: LINE45
      CHARACTER (70) :: LINE70
      CHARACTER (1) :: STAR(NSG), TABSTAR(MAX_ENTRIES)
      CHARACTER (2) :: GEOMS(NSG), GTYPE, MISSING(NSG), TYPE
      CHARACTER (7) :: SPGRPS(NSG)
      CHARACTER (90) :: L_OUT(NSG)
      CHARACTER (110) :: LINE110
      EQUIVALENCE (LINE45,LINE70)
      INTEGER :: I_HERE(NSG), I, I_BLANK, NSKIP, NC, NL, INDEX1,  &
     &           INDEX2, NSOL(NSG), J, N_MISSING, IG, K, MNUB
      REAL :: VALUES(2,MAX_ENTRIES), DMAX, EMIN
      DATA BLANK / '                                 '/
      DATA GEOMS   /'aa', 'ab', 'ca', 'ah', 'af',  &
     &              'ai', 'ak', 'am', 'fa', 'fc',  &
     &              'da', 'db', 'dc', 'dd', 'de',  &
     &              'ap', 'ba', 'bb', 'aq', 'az',  &
     &              'ay', 'bh', 'av', 'bd', 'bf',  &
     &              'au', 'as', &                                ! 12-27-06
     &              'cc', 'cb', 'cd', 'ce', 'ac', 'ad', 'ag',  & ! 9-18-07
     &              'ae', 'an', 'ao', 'fb',  &                   ! 9-18-07
     &              'aj', 'al', 'fd', 'ar',  &                   ! 9-18-07
     &              'at', 'be', 'aw', 'bg', 'ax',  &             ! 9-18-07
     &              'bi', 'bj', 'bk', 'bc', 'df', 'dg', 'cf'/    ! 9-18-07
      DATA SPGRPS /'P1     ', 2*'P-1    ', 2*'P21    ',  &
     &             5*'P21/c  ',  &
     &             'Cc     ', 'C2     ', 3*'C2/c   ',    &
     &             3*'P21212 ',2*'P212121', 2*'Pca21  ', &
     &             5*'Pna21  ', 2*'Pbca   ',2*'Pbcn   ', &       ! 12-27-06
     &             'Pm     ',2*'Pc     ', 'P2     ','P21/m  ', 'P2/m   ',& ! 9-18-07
     &             'P2/m  ',       &                             ! 9-18-07 
     &             3*'P2/c   ',3*'Pnn2   ',    &                 ! 9-18-07
     &           2*'Pba2   ',2*'Pnc2   ','Pmn21  ', 'Pma2   ',&  ! 9-18-07
     &             'P2221  ',2*'Fdd2,C2','Pnma,Cs'/              ! 9-18-07
      DATA I_HERE /NSG*0/
!
      DO I=1,NSG
          L_OUT(I)(1:33) = BLANK
          STAR(I) = ' '             
      ENDDO
!
!----put ....tab files for each geom together in one file named "Summary.save"
         OPEN (UNIT=66, FILE='Summary.save', STATUS='UNKNOWN')     
!----Open unit # 9 which contains the names of the various ...tab files 
         OPEN (UNIT=9, FILE='TAB_FILES', STATUS='OLD', ERR=1000)
      I_BLANK = 0
      I = 0
40    I = I + 1
         READ (9,50,END=1000,ERR=1000) INPUT      ! get name of a ...tab file
         WRITE(66,51) INPUT(4:40)                 
50       FORMAT (A)                        
51       FORMAT (//'##### ',A40/)
!----How many characters in complete file name?
         IF (I_BLANK .EQ. 0) THEN
            I_BLANK = INDEX(INPUT, ' ')
         ENDIF
         OPEN (UNIT=10, FILE=INPUT, STATUS='OLD', ERR=1000)
         NSKIP = 8
         IF (I .EQ. 1) THEN
            READ (10,50) LINE45      ! for 1st tab file,
            WRITE(66,50) LINE45      
!----Remove last 3 characters (eg '_AB') from cpd id
            NC = 14
            DO I=15,45
                IF (LINE45(I:I) .EQ. ' ') GO TO 46
                NC = NC + 1
            ENDDO
            GO TO 47
46          NC = NC - 3 
47          WRITE (*,48) (LINE45(I:I),I=1,NC)
48          FORMAT (45A1)
            READ (10,50) LINE70
            WRITE(66,54) LINE70      ! 4-5-00
            WRITE (*,50) LINE70
54          FORMAT (A70/)
!----Create variable format.
            I_BLANK = I_BLANK + 8
            WRITE (6,85)
85             FORMAT (' '/15X,'[density, g/cc...energy, kcal/mol]')
            WRITE (FMT1,90) I_BLANK, (I_BLANK+30), I_BLANK, (I_BLANK+10),  &
     &                     (I_BLANK+20), (I_BLANK+30)
90          FORMAT ('( T',I2,',','''   max''',',T',I2,',','''   min''',  &
     &         '/',''' tab file name   ''',',T',I2,',','''  density''',  &
     &         ',T',I2,',','''  energy''',',T',I2,',','''  density''',  &
     &         ',T',I2,',','''  energy  sp grp''',')') 
             WRITE (6,FMT1) 
!----Create variable format for lines with numbers
            WRITE (FMT2,8) (I_BLANK-4)
8           FORMAT ('(A',I2,',F6.3,F10.2,5H  |  ,F5.3,F10.2,3X,A7)')
            NSKIP = 6
         ENDIF
         DO J=1,NSKIP
            READ (10,*)   ! skip first 6 or 8 lines
         ENDDO
         NL = 0
150      NL = NL + 1
         LINE = BLANK
         READ (10,50) LINE
         WRITE(66,50) LINE    ! 4-5-00
!155     FORMAT (A24)          ! 9-30-97
         IF (LINE .EQ. BLANK) GO TO 200  
!----Pick up data from tab file
         READ (LINE,175) (VALUES(K,NL),K=1,2), TABSTAR(NL) 
175      FORMAT (4X,F10.3,F9.3,A1)
         GO TO 150
200      NL = NL - 1
         IF (NL .LE. 0) GO TO 40  ! quit if no density data
!----Find maximum density value
         DMAX = VALUES(1,1) 
         INDEX1 = 1
         DO 250 K=2,NL
            IF (DMAX .GE. VALUES (1,K)) GO TO 250
            DMAX = VALUES(1,K)
            INDEX1 = K
250      CONTINUE  
!----Find the minimum energy
         EMIN = VALUES(2,1)
         INDEX2 = 1
         DO 350 K=2,NL
            IF (EMIN .LE. VALUES (2,K)) GO TO 350
            EMIN = VALUES(2,K)
            INDEX2 = K
350      CONTINUE                                 	
!----Which geometry is it...get space group
      GTYPE = INPUT(1:2)
      DO IG=1,NSG
         IF (GTYPE .EQ. GEOMS(IG)) GO TO 365
      ENDDO
365       WRITE (L_OUT(IG), FMT2) INPUT, DMAX, VALUES(2,INDEX1),  &
     &                            VALUES(1,INDEX2), EMIN, SPGRPS(IG)
          STAR(IG) = TABSTAR(INDEX1)          ! 9-30-97  
      NSOL(IG) = 0
31    READ (10,50,END=30) LINE110
      WRITE (66,50) LINE110
      IF (LINE110(1:24).EQ.' Problems are found with') THEN 
         READ (LINE110,32) MNUB
32       FORMAT (24X,I4)                   ! 7-14-08 DU 
         NSOL(IG) = MNUB
      END IF 
      GO TO 31
!----Close present tab file
30    CLOSE (UNIT=10)
      GO TO 40
1000  CONTINUE
!----Write the contents of L_OUT...should be ordered
      J = 0
loop_1 : DO I=1,NSG
         IF (L_OUT(I)(1:10) .EQ. '          ') CYCLE loop_1
         J = J + 1
           IF (NSOL(I).EQ.0) THEN
              WRITE (*,997) STAR(I), J, L_OUT(I)
997           FORMAT (A1,I2,'. ',A90)
           ELSE
              WRITE (*,995) STAR(I), J, L_OUT(I), NSOL(I)
995           FORMAT (A1,I2,'. ',A90, I4,' solution(s) failed')            
           END IF
         ENDDO  loop_1
      WRITE (*,998)          
998   FORMAT(1X,64(1H-))
      DO I = 1,NSG
         IF ( STAR(I) .EQ. '*' ) GOTO 410
      END DO
      GO TO 420
410   WRITE (*,400)                  
400      FORMAT(//1X,'* Warning - Lattice Energy is positive.',  &
     &' Please check tab file'/'   in appropriate directory for',  &
     &' other values of density and energy.')
!----Check to determine if any tab files are missing
420   REWIND 9
      N_MISSING = 0       ! # of missing tab files
      DO 500 I=1,NSG
         READ (9,50,END=600) TYPE
         DO 490 J=1,NSG
            IF (GEOMS(J) .NE. TYPE) GO TO 490
            I_HERE(J) = 1      ! signal that this geometry is present
            GO TO 500
490      CONTINUE   ! dropping thru 490 loop means tab file not found
500   CONTINUE
!----Have any tab files not been found
600   DO 630 I=1,NSG
         IF (I_HERE(I) .EQ. 1) GO TO 630
         N_MISSING = N_MISSING + 1
         MISSING(N_MISSING) = GEOMS(I)
630   CONTINUE
      IF (N_MISSING .NE. 0) THEN   
          PRINT 520, (MISSING(I),I=1,N_MISSING)
520       FORMAT (' One or more tab files were not found, please',  &
     &            ' check the following'/' directory(s): ',  &
     &            10(A2,1X)/,3(15X,10(A2,1X)/))
      ENDIF
      CLOSE (UNIT=9)
      STOP
      END PROGRAM SUMMARIZE_TAB
