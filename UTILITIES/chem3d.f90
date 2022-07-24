!----This program converts Chem3D Cartesian coordinate files
!     to ChemX coordinate files, MacroModel files to Chem3D
!     cartesian coordinate files, Chem3D files to molecular
!     fit format and general cell param, xyz to ChemX format. 
!
      PROGRAM CHEM3D                     ! 2/04/03
!
      IMPLICIT NONE
      LOGICAL I_SEQ2
      CHARACTER (3) :: ATYPE(500), T_ATYPE
      CHARACTER (5) :: CTYPE(500), CGNAME_2, NAME(500)
      CHARACTER (2) :: TYPE(64), ATOM_NAME(9), AT2(500)
      CHARACTER (1) :: AT(500), I_SEQ, CGNAME(5)
      CHARACTER (60) :: FMT
      CHARACTER (72) :: LINE
      CHARACTER (32) :: INFILE, OUTFIL
      INTEGER :: ICON(6,500), NATM, IPICK, NA, N_HYDROGEN, NUM, I, N_HEAVY,&
     &		 N_COUNT, N_PASSES, J, K, KM, JJ, J2, NREMOV, NINS, L, N_ADJ,&
     &           N_OXY, IS, N
      INTEGER :: NCON(500), NTYPE(500), IM(500), JM(500), ITEXT(500),&
     &           ICON_T(6)
      REAL :: XYZ(3,500), X
      REAL :: CELL(6), T(3,3), XYZ_CART(3), CHEMD(64), T_XYZ(3)
      COMMON /CONNECT_INFO/ NATM, ICON, XYZ, AT, AT2, NCON	      
      DATA ATYPE /500*'   '/
      DATA TYPE / 9*'C ',&
     & 4*'XX', 'C ', 5*'O ', 3*'XX', 'O ',&
     & 14*'N ', 2*'XX', 'N ',&
     & 5*'H ', 2*'XX', 'H ',&
     & 2*'S ', 'XX', 'S ',&
     & 'P ', 2*'XX', 'F ', 'CL', 'BR', 'I ',&
     & 'SI', 4*'XX' /
       DATA ATOM_NAME /'H', 'HE', 'LI', 'BE', 'B', 'C', 'N',&
     &                'O', 'F' /
      DATA CHEMD / 4, 2, 1, 1, 1, 1, 2, 2, 1,&                        !  1-9
     &            99,99,99,99, 1, 7, 6, 6,47, 6,99,99,99, 6,&         ! 10-23
     &            10, 9, 8, 8, 8, 9, 9,46,39,39,39,39,50,50,99,99, 8,&! 24-40
     &             5,21,23,48,28,99,99, 5,&                           ! 41-48
     &            15,15,99,15,&                                       ! 49-52
     &            25,99,99,11,12,13,14,&                              ! 53-59
     &            19,99,99,20,0 /                                     ! 60-64
      EQUIVALENCE (CGNAME, CGNAME_2)
!
!----INTRODUCTION.  ASK USER NAME OF INPUT FILE, CHEM3D OR MACMOD.
      WRITE (*,4)
4     FORMAT (' Coordinate transformations (12/26/96)...'/&
     &        '   1. Chem3D files (Cartesian1 format) into ChemX',&
     &        ' files'/&
     &        '   2. MacroModel files into Chem3D (Cartesian1 format)',&
     &        ' files'/&
     &        '   3. ChemX to Chem3D files with dummy atom ID'/&
     &        '   4. ChemX to Chem3D, ChemX line numbers as Chem3D',&
     &        ' connectivity codes'/&
     &        '   (in Chem3D select Cart Coord 1, By Position & ',&
     &        'Expect Atom Type Text Numbers)'/&
     &        '   5. Chem3D files to molecular fit format'/&
     &        '   6. SHELX ...ins file to Chem3D file'/&
     &        '   7. Cell param, general xyz to ChemX; '/&
     &        '      input line 1 = a, b, c, alpha, beta, gamma'/&
     &        '      input line 2 = atom ID, x, y, z; line format',&
     &        ' as (Ax,....) will be requested'/& 
     &        '   8. Cell param, general xyz to ChemX; '/&
     &        '      input line 1 = a, b, c, alpha, beta, gamma'/&
     &        '      input line 2 = atomic number, x, y, z; line',&
     &        ' format as (Ix,....) will be requested'/&
     &        ' Enter transformation number (1 - 8): ',$)
      READ (*,18) IPICK
18    FORMAT(I1)
307   WRITE (*,8)
8     FORMAT (' Enter name of the input coordinate file: ',$)
       READ (*,21) INFILE
21     FORMAT (A32)
!----OPEN THE INPUT FILE
      OPEN (UNIT=10,FILE=INFILE,STATUS='OLD')
!
!----ASK USER NAME OF THE OUTPUT FILE...THE CHEM3D VERSION
      WRITE (*,9)
9     FORMAT (' Enter name of the output coordinate file: ',$)
      READ (*,21) OUTFIL
!----OPEN OUTPUT FILE
      OPEN (UNIT=11,FILE=OUTFIL,STATUS='UNKNOWN')
      NA = 0  
      GO TO (99, 777, 999, 9999, 99, 8000, 8300, 8300), IPICK
!
!----CONVERT CHEM3D FILES TO CHEMX or molecular fit format
!     READ AN ATOM LINE FROM CHEM3D FILE...new Cartesian1 format
!     Discard the 1st line...number of atoms
99    READ (10,*)
         WRITE (*,101) 
101      FORMAT (' Does the Chem3D cc1 file contain atom',&
     &           ' sequence numbers [y]: ',$)
         READ (*,'(A1)') I_SEQ
         I_SEQ2 = .FALSE.
         IF (I_SEQ .EQ. ' ' .OR. I_SEQ .EQ. 'Y' .OR. I_SEQ .EQ.&
     &       'y') I_SEQ2 = .TRUE. 
      N_HYDROGEN = 0
10    NA = NA + 1
      IF (I_SEQ2) THEN
         READ (10,19,END=303) T_ATYPE, T_XYZ, NUM, ICON_T
19       FORMAT (A3,5X,3F12.6,7I5)
      ELSE 
         READ (10,15,END=303) T_ATYPE, T_XYZ, NUM, ICON_T
15       FORMAT (A3,3F12.6,7I5)
      ENDIF
27       ATYPE(NA) = T_ATYPE         ! results.  Some atom numbers could be
         IF (T_ATYPE .EQ. '  H') N_HYDROGEN = N_HYDROGEN + 1
         DO I=1,3                    ! missing.
            XYZ(I,NA) = T_XYZ(I)
         ENDDO
         DO I=1,6
            ICON(I,NA) = ICON_T(I)
         ENDDO         
      GO TO 10
303   NATM = NA - 1
!
!----Finished with Chem3D file input....ChemX or molecular fit format?
      IF (IPICK .EQ. 1) GO TO 327 
!----Write 5 header lines for the molecular fit file
      N_HEAVY = NATM - N_HYDROGEN
      WRITE (11,345) N_HEAVY, NATM, NATM 
345   FORMAT ('/u/ammon/predictions/molecular-fit.out << END'/&
     &       '    1'/'  molecular fit calculation'/&
     &       3I5,'    0    1    5    0    1'/&
     &       '    1.0000    1.0000    1.0000    90.000    90.000',&
     &       '    90.000')
      N_COUNT = 0
      DO 365 N_PASSES=1,2
         DO 375 J=1,NATM
            IF (ATYPE(J) .NE. '  H' .AND. N_PASSES .EQ. 1) GO TO 379
            IF (.NOT. (N_PASSES .EQ. 2 .AND. ATYPE(J) .EQ. '  H'))& 
     &                     GO TO 375
379            N_COUNT = N_COUNT + 1
               IF (J .LE. 9) THEN
                 WRITE (11,380) N_COUNT, ATYPE(J)(3:3), J,& 
     &                          (XYZ(K,J),K=1,3)
380              FORMAT (I4,1X,A1,'_',I1,2x,3F10.5) 
               ELSE
                 WRITE (11,381) N_COUNT, ATYPE(J)(3:3), J,& 
     &                          (XYZ(K,J),K=1,3)
381              FORMAT (I4,1X,A1,'_',I2,1x,3F10.5) 
               ENDIF
375      CONTINUE
365   CONTINUE
      WRITE (11,382)
382   FORMAT ('END'/&
     &        'mv fort.16 molecular-fit.output-file')
      STOP
!
!----PREPARE CHEMX FILE
327   WRITE(11,16) NATM
16    FORMAT(38X, 3(3X, '1.000')/21X, 3(2X, '90.000')/I4, 3X,&
     &'1'/5X, '0')
!----SQUEEZE BLANKS OUT OF ATOM NAMES
      DO 20 I = 1, NATM
         DO J = 2, 5
            CGNAME(J) = ' '
         ENDDO
         WRITE (CGNAME_2,51) ATYPE(I), I
51          FORMAT(A3,I2)
         KM = 0
42    IF (CGNAME(1) .EQ. ' ') THEN
         KM = 1
         GO TO 66
      ELSE
         KM = 2
      END IF
      IF (CGNAME(2) .NE. ' ') KM = 3
43    IF (CGNAME(KM) .NE. ' ') GO TO 997
66    DO 54 K = KM,4
54    CGNAME(K) = CGNAME(K+1)
      CGNAME(5) = ' '
      IF (KM .EQ. 1) THEN
         GO TO 42
      ELSE
         GO TO 43
      END IF
997   NAME(I) = CGNAME_2
!----WRITE AN ATOM LINE FOR CHEMX
      WRITE(11,79) I, NAME(I), (XYZ(J,I), J=1,3), (ICON(J,I),&
     & J=1,6)
79    FORMAT (1X,I3,1X,A5,3F10.6,6I4,'   0   0   0   0   1')
20    CONTINUE
      GO TO 304
!
!----PREPARING CHEM3D FILE FROM MACROMODEL FILE----
777   READ(10,*) NATM        ! number of atoms on line 1
300   DO 301 NA=1,NATM
        READ (10,22) NTYPE(NA),(ICON(J,NA),J=1,6), (XYZ(J,NA),J=1,3)
22      FORMAT(I4,6(I6,2X),3F12.6)
301   CONTINUE
      JJ = 0
      J2 = 0
      DO 3000 I = 1, NATM
      IF (TYPE (NTYPE(I)) .EQ. 'XX' ) THEN
      JJ = JJ + 1  
      IM ( JJ ) = I
      ELSE
      J2 =J2 + 1
      JM (J2) = I
      END IF
3000  CONTINUE
      DO 4000 I = 1, J2
      DO 4001 J =1, 3
4001  XYZ ( J, I ) = XYZ ( J, JM (I) )    
      DO 4002 J =1, 6
4002  ICON ( J, I ) = ICON ( J, JM(I) )
4000  CONTINUE    
      NATM = NATM - JJ
      DO 1000 I =1, NATM
      DO 2001 K =1, JJ
      DO 2000 J =1, 6
      IF ( ICON ( J, I ) .EQ. IM ( K ) ) ICON ( J, I ) = 0
2000  CONTINUE
2001  CONTINUE      
1000  CONTINUE
      DO 100 I= 1, NATM
      DO 201 J =1, 6
      NREMOV = 0
      DO 200 K =1, JJ
      IF ( ICON ( J, I ) .GE. IM ( K ) ) NREMOV = NREMOV + 1
200   CONTINUE
      IF (NREMOV .EQ. 0) GO TO 201
      ICON(J,I) = ICON(J,I) - NREMOV
201   CONTINUE
100   CONTINUE
!----How many atoms are connected to each atom?
      DO 450 I=1,NATM
        NCON(I) = 0
          DO 450 J=1,6
            IF (ICON(J,I) .NE. 0) NCON(I) = NCON(I) + 1
450   CONTINUE
      WRITE ( 11, 50 ) NATM
50    FORMAT(I3)
      DO 91 I=1,J2
        ITEXT(I)= CHEMD(NTYPE(JM(I)))
        NTYPE(I) = NTYPE(JM(I))
91    CONTINUE
!----Special fix for nitro group...any CHEM3D type 47 atom (oxygen)
!     that is connected to only 1 other atom which is a N, will be 
!     given type 7 as a nitro oxygen
      DO 470 I=1,NATM
        IF (NCON(I) .GT. 1) GO TO 470
        IF (ITEXT(I) .NE. 47) GO TO 470
        IF (TYPE(NTYPE(ICON(1,I))) .EQ. 'N ') ITEXT(I) = 7
470   CONTINUE
      DO 60  I = 1, NATM
!----WRITE AN ATOM LINE FOR CHEM3D FILE
 56   WRITE(11,58) TYPE(NTYPE(I)),I,(XYZ(J,I),J=1,3),ITEXT(I),&
     &             (ICON(J,I),J=1,6)
58    FORMAT(2X,A2,I4,3F12.6,7I5)
60    CONTINUE
      GO TO 304
!
!----PREPARING CHEM3D FILE (with dummy atoms) FROM ChemX file....
999   READ (10,*)
      READ (10,*)
      READ(10,*) NATM        ! number of atoms on line 1
      WRITE (11,74) NATM
74    FORMAT (I3)
      READ (10,*)
      DO 1301 NA=1,NATM
        READ (10,221) CTYPE(NA), (XYZ(J,NA),J=1,3), (ICON(J,NA),J=1,8)
221     FORMAT (5X,A5,3F10.5,I5,7I4) 
        WRITE(11,258) CTYPE(NA), (XYZ(J,NA),J=1,3),&
     &             (ICON(J,NA),J=1,6)
258     FORMAT (2X,A5,1X,3F12.6,'   99',6I5)
1301   CONTINUE
       GO TO 304
!
!----Prepare CHEM3D file from ChemX file with ChemX
!     line numbers as Chem3D connectivity codes or from
!     SHELX ...ins file
!
!----Get cell params from SHELX ...ins file
8000  READ (10,8005) LINE
8005  FORMAT (A72)
      IF (LINE(1:4) .NE. 'CELL') GO TO 8000
      READ (LINE(5:72),*) X, CELL
!
      NATM = 999
      NINS = 0
!----Skip down to FVAR
8004  READ (10,8005) LINE
      IF (LINE(1:4) .eq. 'FVAR' .or. line(1:4) .eq. 'SFAC') GO TO 8010
      GO TO 8004 
!----Read cell parameters from the first 2 lines...a, b, c, 3 angles
9999  READ (10,11) CELL
11    FORMAT (38X,3F8.3/21X,3F8.3)
!----Number of atoms from 3rd CHEMX line
      READ (10,14) NATM
14    FORMAT (I4)
!----Skip 1 line
      READ (10,*) 
      WRITE (11,74) NATM                        ! first CHEM3D output line
8010  DO NA=1,NATM
        IF (IPICK .EQ. 6) THEN
8011       READ (10,8005) LINE
           IF (LINE(1:4) .EQ. 'HKLF' .OR. LINE(1:4) .EQ. 'END ') GO TO 7000
           IF (.NOT. (LINE(1:1) .EQ. 'H' .OR.&
     &                LINE(1:1) .EQ. 'C' .OR.&
     &                LINE(1:1) .EQ. 'S' .OR.&
     &                LINE(1:1) .EQ. 'O' .OR.&
     &                LINE(1:1) .EQ. 'F' .OR.&
     &                LINE(1:1) .EQ. 'N' )) GO TO 8011
           NINS = NINS + 1
           READ (LINE,8013) AT2(NA), (XYZ(J,NA),J=1,3)
8013       FORMAT (A1,9X,3F13.5)
           DO J=1,6
              ICON(J,NA) = 0
           ENDDO
        ELSE          
           READ (10,721) AT2(NA), (XYZ(J,NA),J=1,3)  ! chemx format
721        FORMAT (5X,A2,3X,3F10.5)                  ! coordinates   
        ENDIF
        NTYPE(NA) = NA
        NCON(NA) = 0
!----The 2nd character can only be an "I"...as in SI.  Blank out all others
        IF (AT2(NA) .EQ. 'SI' .OR. AT2(NA) .EQ. 'Si') CYCLE 
           AT2(NA)(2:2) = ' '
       ENDDO
7000   IF (IPICK .EQ. 6) THEN 
         NATM = NINS
         WRITE (11,74) NATM
      ENDIF
!
!----Is the cell Cartesian or fractional?  Test the cell lengths
      IF (CELL(1) .LE. 1.001 .AND. CELL(2) .LE. 1.001 .AND.&
     &    CELL(3) .LE. 1.001) GO TO 1400
!----It's fractional...convert to Cartesian
      CALL MATRXT (CELL, T)                    ! fractional to cart matrix
      DO 1390 I=1,NATM
         CALL MTIMES2 (XYZ(1,I), T, XYZ_CART)  ! apply transformation 
         DO 1395 J=1,3
             XYZ(J,I) = XYZ_CART(J)
1395     CONTINUE     
1390  CONTINUE
!
!----Is it necessary to build the ICON array?  If all the ICON(1,...) 
!     are 0, then yes.
1400  DO 1370 NA=1,NATM
          IF (ICON(1,NA) .NE. 0) GO TO 1401
1370  CONTINUE  
      CALL CONNECT                          ! build ICON
1401  DO 1303 NA=1,NATM
!----How many connected?
        IF( NCON(NA) .NE. 0) GO TO 1375
        DO L=1,6
           IF (ICON(L,NA) .EQ. 0) GO TO 1375
           NCON(NA) = NCON(NA) + 1
        ENDDO
1375    IF (AT2(NA)(1:1) .EQ. 'H') THEN                    ! H
           K = 5
           GO TO 158
        ENDIF
        IF (AT2(NA)(1:1) .EQ. 'F') THEN                    ! F
           K = 11                                  
           GO TO 158
        ENDIF
        IF (AT2(NA)(1:1) .EQ. 'C') THEN
           GO TO (1001, 1001, 1003, 1004), NCON(NA)
1001       K = 4                                        ! alkyne C
           GO TO 158
1003       K = 2                                        ! alkene C
           GO TO 158
1004       K = 1                                        ! alkane C
           GO TO 158
        ENDIF
        IF (AT2(NA) .EQ. 'SI' .OR. AT2(NA) .EQ. 'Si') THEN  ! Si
           K = 19
           GO TO 158
        ENDIF
        IF (AT2(NA)(1:1) .EQ. 'S') THEN
           GO TO (1011, 1011, 1013, 1014), NCON(NA)
1011       K = 15                                       ! thio ether
           GO TO 158
1013       K = 17                                       ! suloxide S
           GO TO 158
1014       K = 18                                       ! sulfone S
        ENDIF
        IF (AT2(NA)(1:1) .EQ. 'O') THEN
           GO TO (1021, 1022), NCON(NA)
1021          N_ADJ = ICON(1,NA)           ! get ID of single connected atom
              K = 7                        ! default is carbonyl
              IF (AT2(N_ADJ)(1:1) .EQ. 'C') GO TO 158       
              IF (NCON(N_ADJ) .NE. 3) GO TO 158
                 N_OXY = 0                 ! count # O's connected to N_ADJ
                 DO I=1,3
                    IS = ICON(I,N_ADJ)
                    IF (AT2(IS)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
                 ENDDO
                 IF (N_OXY .NE. 2) GO TO 158
                    K = 47                              ! nitro O
                 GO TO 158
1022       K = 6                                        ! ether O
           GO TO 158
        ENDIF
        IF (AT2(NA)(1:1) .EQ. 'N') THEN
           GO TO (1031, 1032, 1033, 1034), NCON(NA)
1031       K = 10                                       ! nitrile N
           GO TO 158
1032       K = 40                                       ! enamine N
           GO TO 158
1033       K = 8                                        ! amine N default
              N_OXY = 0                ! count # of O's on this N
              DO I=1,3
                 N_ADJ = ICON(I,NA)                 
                 IF (AT2(N_ADJ)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
              ENDDO
              IF (N_OXY .EQ. 2) K = 46                  ! nitro N
           GO TO 158
1034       K = 39                                       ! ammonium N
           GO TO 158
        ENDIF
158     WRITE(11,159) AT2(NA), (XYZ(J,NA),J=1,3), K,&   ! chem3d
     &             (NTYPE(ICON(J,NA)),J=1,NCON(NA))     ! format
159     FORMAT (2X,A2,F11.6,2F12.6,7I5)                 ! coordinates
1303  CONTINUE
      GO TO 304              ! finished
!
!----Option # 7.  Cell param, xyz to ChemX format.
8300  READ (10,*) CELL      ! a b c alpha beta gamma
      WRITE (11,11) CELL     ! 2 lines in ChemX format
!----Read atom coord line...first obtain format
      IF (IPICK .EQ. 7) THEN
        PRINT 8320
8320    FORMAT ('Format of atom coordinates lines (id, x, y, z), eg. ',&
     &          '(A4,3F10.5)...')      
      ELSE
        PRINT 8325
8325    FORMAT ('Format of atom coordinates lines',&
     &          ' (at #, x, y, z), eg. (I4,3F10.5)...')
      ENDIF
      READ (5,'(A60)') FMT
      N = 0
8330  N = N + 1
      IF (IPICK .EQ. 7) THEN
         READ (10,FMT,END=8331) CTYPE(N), (XYZ(J,N),J=1,3)
      ELSE
         READ (10,FMT,END=8331) K, (XYZ(J,N),J=1,3)
         CTYPE(N) = ATOM_NAME(K)
      ENDIF
      GO TO 8330
8331  N = N - 1
      WRITE (11,8340) N
8340  FORMAT (I4)
      WRITE (11,8340)            ! 4th line is blank
      DO 8350 I=1,N
         WRITE (11,8360) I, CTYPE(I), (XYZ(J,I),J=1,3)
8360     FORMAT (I4,1X,A5,3F10.5)
8350  CONTINUE
!
304   PRINT 305
305   FORMAT (' Coordinate transformation complete')
      STOP
      END PROGRAM CHEM3D
!
!----CONNECT...build ICON array
      SUBROUTINE CONNECT
!
      CHARACTER AT*1, AT2*2
      COMMON /CONNECT_INFO/ NATM, ICON(6,500), XYZ(3,500), AT(500),&
     &                 AT2(500), NCON(500)
!
      DO 100 I=1,NATM
         NCON(I) = 0
         DO 90 J=1,NATM
            IF (I .EQ. J) GO TO 90
            IF (AT2(I)(1:1) .EQ. 'H' .AND. AT2(J)(1:1) .EQ. 'H')&
     &                         GO TO 90
               IF (AT2(I) .EQ. 'SI' .OR. AT2(J) .EQ. 'SI') THEN
                  DMAX = 1.90
               ELSE
                  DMAX = 1.60
               ENDIF
               D = 0.0
               DO K=1,3
                  D = D + (XYZ(K,I) - XYZ(K,J))**2
               ENDDO
               IF (SQRT(D) .LE. DMAX) THEN
                  NCON(I) = NCON(I) + 1
                  ICON(NCON(I),I) = J
               ENDIF
90       CONTINUE
100   CONTINUE
      RETURN
      END SUBROUTINE CONNECT  

!----Matrix for coordinate system transformation
      SUBROUTINE MATRXT (A,T)
      REAL :: W
      DIMENSION T(3,3),A(6)
!
      T(1,1)=A(1)
      T(1,2)=A(2)*COSD(A(6))
      T(1,3)=A(3)*COSD(A(5))

      T(2,1)=0.0
      T(2,2)=A(2)*SIND(A(6))
      T(2,3)=A(3)*(COSD(A(4))-COSD(A(5))*COSD(A(6)))/SIND(A(6))
 
      W=SQRT(1-COSD(A(4))**2-COSD(A(5))**2-COSD(A(6))**2+&
     &2*COSD(A(4))*COSD(A(5))*COSD(A(6)))
      T(3,1)=0.0
      T(3,2)=0.0
      T(3,3)=A(3)*W/SIND(A(6))

!     PRINT 100,((T(I,J),I=1,3),J=1,3)
! 100 FORMAT(3X,3F8.5)
      RETURN
      END SUBROUTINE MATRXT 

!----Matrix conversion
      SUBROUTINE MTIMES2 (XYZG,RV,XYZT) 
      DIMENSION XYZT(3),XYZG(3),RV(3,3)
      XYZT(1)=RV(1,1)*XYZG(1)+RV(1,2)*XYZG(2)&
    &+RV(1,3)*XYZG(3)
      XYZT(2)=RV(2,1)*XYZG(1)+RV(2,2)*XYZG(2)&
    &+RV(2,3)*XYZG(3)
      XYZT(3)=RV(3,1)*XYZG(1)+RV(3,2)*XYZG(2)&
    &+RV(3,3)*XYZG(3)
      RETURN
      END SUBROUTINE MTIMES2 
!
      real function SIND(x)
      SIND = SIN(x/57.29577951)
      return
      end function SIND
!
      real function cosd(x)
      if (abs(abs(x) - 90.0) .le. 0.000001) then
         cosd = 0.0
      else
        cosd = cos(x/57.29577951)
      endif
      return
      end function cosd
