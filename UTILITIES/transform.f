C----TRANSFORM.F            #38               

C    latest version version of Dec 2010                 
C    Added ability yo write a new wmin.save format file 12/19/06
C    Added params for pyridinium N+(-H), N = 85 and    7/20/06
C        H = 86                                        7/20/06
C    Added params # 83 and 84 for the N and O atoms    7/20/06
C        in NO3- (nitrate anion).                      7/20/06
C    Added diazo N...C-N=N-C (# 82)                    2/13/06
C    Added ester...C-C(=O)-O-C; # 79, 80, 81. Note that 2/3/06
C        anhydrides have the same C and O codes.        2/4/06
C    Added S (# 77) and N (# 78) of thiazole.          1/31/06
C    Altered the id of H on NH2 so it can be either    1/25/06
C        C-NH2 or N-NH2 (H = # 3)                      1/25/06
C    Added iq_aminoimine_2 and iq_aminoimine_3         1/24/06
C        subroutines ti identify the 2 N's in an       1/24/06
C        aminoimine                                    1/24/06
C    Added 1,2,3,4-tetrazole subroutines to properly    1/2/06 
C        identify atoms # 70, 71, 72, and 73.         12/31/05
C    Added triazole subroutines (iq_123triazole_2 and 12/30/05
C        iq_123triazole_3) to properly id atoms       12/30/05
C        # 26 and 27.                                 12/30/05
C    For option 10, added the ability to apply a      3/13/04
C       coordinate translation.                       3/13/04
C    Added ability to write a very simple CIF file    3/9/04
C       from a WMIN ...save file.                     3/9/04
C    Added ability to write a simple SHELX res file   3/2/04
C       from a WMIN ...save file.                     3/2/04
C    In writing the molpak.xyz file, H (#64) and      2/3/04
C       O (#65) will have atom names of "H 64"        2/3/04
C       and "O 65".                                   2/3/04
C    Added H of [Csp3]2-NH. N = 69, H = 68            1/15/04
C    Added an internal alkyne, C = 67                 12/31/03
C    Added the hydantoin moiety...                    12/3/03
C       --N(H)--C(O)--N(H)--C(O)--                    12/3/03
C       H = 64, O = 65, N = 67                        12/3/03
C    Added amide H (61), N (62) and O (63).           11/15/03
C       C(sp2)-NH-C(=O)                               11/15/03
C    Fixed format for reading molpak.xyz file under   7/3/03
C      option 8a                                      7/3/03
C    Altered format for reading ESP charge from       5/31/03
C      Gaussian log file to accomodate G03            5/31/03
C    Added separate O & H for alcohols (# 58, 59)     3/25/03
C      and ability to stretch O-H                     3/25/03
C    Added separate N's for azides (# 32,56,57)       3/5/03
C    Added separate NO2 for nitramines                2/24/03
C    Added nitrate esters (-O-NO2)                    12/14/02
C    Added F or C-F                                   1/15/03
C    Modified it for nitrocubane                      11/12/02
C
C----Carry out coordinate transformations between CHEMX, MOLPAK,
C     MOPAC, MACROMODEL, G03, CADPAC and CRYSTAL.
C
C----Altered to determine atom types for the C-SO2-N=C group.   12/11/00
C----Altered to include atom types for imides; the N, O, H of
C     C-C(=O)-NH-C(=O)-C                                        1/1/01
C
      PROGRAM TRANSFORM
      PARAMETER (NSWITCH = 34)
      CHARACTER C_SWITCH(NSWITCH)*2, AH*2 
      DIMENSION N_SWITCH(NSWITCH) 
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      DATA C_SWITCH /'1a', '1b', '3a', '4d', '4a', '5 ', '4b', '4c',
     x               '6a', '6b', '2a', '2b', '6c', '7a', '7b', '2c', 
     x               '8a', '3b', '1c', '1d', '8c', '2d',  '9', '10',
     x               '4e', '7c', '6d', '11', '8b', '7d', '7e', '12',
     x               '13', '8d'/
      DATA N_SWITCH / 1,    17,   2,    3,    4,    5,    6,    7,
     x                8,    9,    10,   11,   12,   13,   14,   15,
     x                16,   18,   19,   20,   21,   22,   23,   24,       
     x                25,   26,   27,   28,   29,   30,   31,   32,
     x                33,   34/
C
       CONNECT_FLAG = .FALSE.
       ICHARGE_FLAG = -1
       G03_CHARGE = 0.0
       XMNDO_CHARGE = 0.0
C
       PRINT 50
50    FORMAT   (' Coordinate transformations between CHEMX, CHEM3D,',
     X          ' MOPAC, MOLPAK, MACROMODEL, WMIN,'/
     x          '    CRYSTAL, G03, free and CARTESIAN formats'/
     x          '    (version of Dec 2010)....'//
     x          '    What do you want to do?.....'//
     x          '   (1a) transform a CHEMX file to a MOPAC ...dat ',
     x                  'file for geometry optimization'/
     x          '   (1b) transform a CHEMX file to MOLPAK format'/
     x          '   (1c) transform a CHEMX file to MOLPAK format & ',
     x                  'reorient molecule to ly along a C2 axis'/
     x          '   (1d) transform a CHEMX file to MOLPAK format for ',
     x                  'a molecule on a center of symmetry; the'/
     x          '        MOLPAK file will contain 1/2 molecule wrt ',
     x                  'a center at 0,0,0'/
     x          '   (2a) transform a CHEM3D cart coord 1 file to',
     x                 ' CARTESIAN COORDINATE format'/
     x          '   (2b) transform a CHEM3D cart coord 1 file to',
     x                 ' generic atom type, x, y, z file for g94'/
     x          '   (2c) transform a CHEM3D file to MOLPAK format'/
     x          '   (2d) transform a CHEM3D file to MACROMODEL format'/
     x          '   (3a) transform a MOPAC file to CHEMX format'/
     x          '   (3b) transform a MOPAC file to MACROMODEL format'/
     x          '   (4a) transform a MOPAC file to MOLPAK format'/
     x          '   (4b) transform a MOPAC file to',
     x                 ' CARTESIAN COORDINATE format'/
     x          '   (4c) transform a MOPAC coordinate file to a new ',
     x                  'MOPAC ...dat ',
     x                  'file for a MNDO-ESP calculation'/ 
     x          '   (4d) transform a MOPAC coordinate file to a new ',
     x                  'MOPAC ...dat file for geometry ',
     x                  'optimization'/
     x          '   (4e) transform a MOPAC coordinate file to a g94 ',
     x                   'input file for ESP/CHELPG charge ',
     x                   'calculation'/
     x          '    (5) transform a CARTESIAN COORDINATE file to a ',
     x                  'MOPAC ...dat file for geometry optimization'/
     x          '   (6a) transform a MACROMODEL file to CHEMX format'/
     x          '   (6b) transform a MACROMODEL file to ',
     x                  'CARTESIAN COORDINATE format'/
     x          '   (6c) transform a MACROMODEL file to MOLPAK format'/
     x          '   (6d) transform a MACROMODEL file to a G03 input ',
     x                   'file for ESP/CHELPG charge calculation'/
     x          '   (7a) transform a G03 intermediate xyz file to',
     x                 ' CARTESIAN COORDINATE format'/
     x          '   (7b) transform a free format file with at. no., x,',
     x                 ' y, z to CARTESIAN COORDINATE format'/
     x          '   (7c) transform a G03 intermediate xyz file to',
     x                 ' CHEMX format'/
     x          '   (7d) transform a final, compressed-format G03 xyz',
     x                 ' file to CHEMX format'/ 
     x          '   (7e) transform a final, compressed-format G03 xyz',
     x                 ' file to normal G94 xyz format'/ 
     x          '   (8a) transform a MOLPAK file to CHEM3D format'/
     x          '   (8b) transform a MOLPAK file for a model with Ci ',
     x                 ' symmetry to a 1/2 molecule file;'/
     x          '        average coordinates & atomic charges where ',
     x                 'necessary'/
     x          '   (8c) transform a general format file with atom ',
     x                 'name, x, y, z to CHEM3D format...input format',
     x                 ' required'/
     x          '   (8d) transform a general format file with atomic ',
     x                 'number, x, y, z to CHEM3D format...input ',
     x                 'format required'/
     x          '    (9) transform a molecular fit output file ',
     x                 '(coords of 2nd molecule) to WMIN xyz format'/
     x          '   (10) transform final WMIN fractional coordinates ',
     x                 'to CHEMX, CRYSTAL, CCDC, CIF',
     x                 ' and wmin.save formats;'/                    ! 12/19/06 
     x          '         a reduced cell matrix may be applied'/       ! 3/9/04
     x          '   (11) transform a HIN file to MACROMODEL FORMAT'/,
     x          '   (12) transform coordinates (Bohrs) in a CADPAC ',
     x                 'DMA output file to G03 xyz format'/,
     x          '   (13) transform coordinates (Angs) in a CADPAC ',
     x                 'input file to G03 xyz format'/,
     X          '    (0) exit')
      print 100
100   FORMAT (' Please select a number: ',$)
      READ (5,'(A2)') AH
C----translate the type code into a number
      DO 101 I=1,NSWITCH
         IF (AH .NE. C_SWITCH(I)) GO TO 101
         KIND = N_SWITCH(I)
         GO TO 95
101   CONTINUE
      PRINT 105, AH
105   FORMAT (' **TILT...cannot decode selection ',A2)
      STOP
C----        0  1a 3a 4d 4a 5  4b 4c 6a 6b  2a  2b  6c  7a  7b  2c  8a---
C----        0  1  2  3  4  5  6  7  8   9  10  11  12  13  14  15  16----
C
C----       1b  3b  1c  1d  8c  2d   9  10  4e  7c  6d  11  8b  7d  7e----
C----       17  18  19  20  21  22  23  24  25  26  27  28  29  30  31----
C             
C----       12  13  8d
C----       32  33  34
95    GO TO (1, 2, 3, 4, 5, 2, 3, 4, 9, 9,  10, 10,  5,  3,  3,  5,  5,
     x       6, 11, 19, 21, 5,  11, 23, 24,  4,  3,  4, 11, 29, 30, 30, 
     x      32, 32,  5) KIND + 1 
2        CALL CHEMX2MOPAC (KIND)
            GO TO 20 
3        CALL MOPAC2CHEMX (KIND)
            GO TO 20
4        CALL NEWMOPAC_OR_G03 (KIND)
            GO TO 20
5        CALL MOPAC2MOLPAK (KIND)
            GO TO 20
6        CALL CHEMX2MOLPAK
            GO TO 20
9        CALL MMOD2OTHER (KIND)
            GO TO 20
10       CALL CHEM3D2CC (KIND)
            GO TO 20
11       CALL MOPAC2MMOD (KIND)
            GO TO 20
19       CALL CHEMX_C2_MOLPAK
            GO TO 20
21       CALL CHEMX_CI_MOLPAK
            GO TO 20
23       CALL MOLFIT2WMIN
            GO TO 20
24       CALL WMIN2CHEMX
            GO TO 20
29       CALL MOLPAK_CI_MOLPAK
            GO TO 20
30       CALL G94TOCHEMX (KIND)
            GO TO 20
32       CALL CADPAC_DMAREL_TO_G94 (KIND)
C
20     PRINT *,(' Transformation Finished ')  
1      STOP
       END PROGRAM TRANSFORM    
C
C-----------------------------------------------------------------------
C 
C     Convert from CHEMX or fractional coordinates to MOPAC
C
      SUBROUTINE CHEMX2MOPAC (KIND)
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      DIMENSION IK(100),XYZ(3,100),A(6), TXYZ(3), CELL(6), T(3,3)
      DIMENSION LIST(3), XYZ_CART(3,100)
      DATA A/1.0,1.0,1.0,90.0,90.0,90.0/
      CHARACTER*72 TITLE
      CHARACTER NAME(100)*2, TNAME*2, OUTFILE2(32)*1,
     X          INFILE2(32)*1, TRIAL_FILE(32)*1, DOT_DAT(4)*1
      CHARACTER*32 INFILE,OUTFILE
      CHARACTER FORMAT1(46)*1, FORMAT2*46, FORMAT3*2, NUMBERS(9)*1
      EQUIVALENCE (OUTFILE,OUTFILE2), (INFILE, INFILE2)
      EQUIVALENCE (FORMAT1, FORMAT2), (FORMAT1(31),FORMAT3)
      DATA FORMAT2 /'(25H Mopac OUTPUT file name [,xxA1,3H]: ,$) '/      
C                    1234567890123456789012345678901234567890123456
      DATA OUTFILE2(1), TRIAL_FILE/33* ' '/
      DATA DOT_DAT/'.', 'd', 'a', 't'/
      DATA NUMBERS/'1','2','3','4','5','6','7','8','9'/
C
      IF (KIND .EQ. 1) THEN
         print 1
    1    FORMAT(' ChemX INPUT file name: ',$) 
      ELSE
         PRINT 11
11       FORMAT (' Fractional coordinate input file name: ',$)
      ENDIF
      read (5,2) INFILE
    2 FORMAT(A32)

C     OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
C----Create a possible OUTFILE name
      DO 98 I=1,32
         IF (INFILE2(I) .eq. '.' .or. infile2(i) .eq. ' ') GO TO 96
            GO TO 95
96       DO 97 J=0,3
97         TRIAL_FILE(I+J) = DOT_DAT(J+1)                     
         NC = I + 3
         GO TO 99
95       TRIAL_FILE(I) = INFILE2(I)         
98    CONTINUE
99    IF (NC .GE. 10) THEN    
         WRITE (FORMAT3,301) NC            ! create format...nnA1
301      FORMAT (I2)
      ELSE
         WRITE (FORMAT3,302) NC            ! or nA1
302      FORMAT (1x,I1)
      ENDIF
      PRINT FORMAT2, (TRIAL_FILE(J),J=1,NC)
      read (5,2) OUTFILE
      IF (OUTFILE2(1) .NE. ' ') GO TO 72
      DO 86 I=1,32
86       OUTFILE2(I) = ' '
      DO 87 I=1,NC                           ! use the default output file
87       OUTFILE2(I) = TRIAL_FILE(I)         ! name
C     OPEN THE OUTFILE
72    OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
      print 6
6     FORMAT (' TITLE: ',$)
      read (5,4) TITLE
    4 FORMAT(A72)
C----If ChemX file, read cell params from lines 1 and 2
      IF (KIND .EQ. 1) THEN
         READ (10,770) CELL
770      FORMAT (38X,3F8.3/21X,3F8.3)
         NS = 2                         ! skip 2 lines        
      ELSE
         NS = 4                         ! skip 4 lines
      ENDIF
      DO I=1,NS
         READ (10,*)           
      END DO
      NA=0
777   NA=NA+1
         READ (10,30,END=999) IK(NA), NAME(NA), (XYZ(I,NA),I=1,3)
30       FORMAT (I4,1X,A2,3X,3F10.5)
         DO L=1,9
            IF (NAME(NA)(2:2) .EQ. NUMBERS(L)) THEN
               NAME(NA)(2:2) = ' '
               GO TO 777
            ENDIF
         ENDDO
      GO TO 777
  999 NA=NA-1
C----If a ChemX file, go thru conversion from fractional to cartesian 
C     just for good measure
      IF (KIND .NE. 1) GO TO 7777
          CALL MATRXT (CELL, T)
          CALL MTIMES1 (NA, CELL, XYZ, T, XYZ_CART)
          DO 779 I=1,NA
             DO 779 J=1,3
                XYZ(J,I) = XYZ_CART(J,I)
779       CONTINUE
7777  WRITE (11,10)
   10 FORMAT(1X,' XYZ NOINTER MMOK AM1 PRECISE T=240H  ')
      WRITE (11,15) TITLE
   15 FORMAT (2X,A72)
      WRITE (11,20)
   20 FORMAT (1X,'--------------------------------------------------')
C----Before outputing coordinates, make certain that first three atoms
C     are connected to each other
      ICNT = 0
      DO 100 I=1,NA
         IF (NAME(I) .EQ. 'H') GO TO 100
         ICNT = ICNT + 1
         LIST(ICNT) = I
         DO 101 J=I+1,NA
             IF (NAME(J) .EQ. 'H') GO TO 101
             DSUM = 0.0
             DO L=1,3                                   ! calc distance between
                DSUM = DSUM + (XYZ(L,I) - XYZ(L,J))**2  ! atoms I and J
             ENDDO
             DSUM = SQRT(DSUM) 
             IF (DSUM .LE. 1.85) GO TO 127  ! Yes, if d .le. 1.85 Angs
                GO TO 101
127          ICNT = ICNT + 1
             LIST(ICNT) = J 
             IF (ICNT .GE. 3) GO TO 105   
101      CONTINUE
100   CONTINUE
C----The 3 attached atom sequence numbers are now in LIST(1-3) 
105   DO 120 I=1,3
         IF (LIST(I) .EQ. I) GO TO 120
            L = LIST(I) 
            TNAME = NAME(L)
            NAME(L) = NAME(I)
            NAME(I) = TNAME
            DO K=1,3
               TXYZ(K) = XYZ(K,L)
               XYZ(K,L) = XYZ(K,I)
               XYZ(K,I) = TXYZ(K)
            ENDDO
120   CONTINUE            
      WRITE (11,25)(NAME(K),(XYZ(I,K),I=1,3),K=1,NA)  
   25 FORMAT (1X,A2,3X,F10.5,2X,'1',F10.5,2X,'1',F10.5,2X,'1')
      WRITE (11,35)
   35 FORMAT (1X,'0'//)
      RETURN
      END
C-----------------------------------------------------------------
	BLOCK DATA ONE
        CHARACTER*5 NAME
        COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)     
        DATA CELL /1.0, 1.0, 1.0, 90.0, 90.0, 90.0/
        END
C
C------------------------------------------------------------------------
C
C     Convert from MOPAC, MACROMODEL, CHEM3D or general format to MOLPAK 
C      format and from MOLPAK to CHEM3D
C
      SUBROUTINE MOPAC2MOLPAK (KIND)
      CHARACTER NAME*5, YN*1, MMOD(64)*2, AT2*2, LINE*60
      CHARACTER*32 INFILE,OUTFILE
      CHARACTER INPUT_FORMAT*60, Q*1, ATOMNAME(9)*2
      LOGICAL ISEQ, CONNECT_FLAG
      DIMENSION IK(200)
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE  
      DATA ATOMNAME /'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ',
     X               'O ', 'F '/
      DATA MMOD /9*'C ', '  ', 2*'C ', '  ', 'C ',
     x           7*'O ', '  ', 'O ',
     x           14*'N ', 2*'  ', 'N ',
     x           5*'H ', 2*'  ', 'H ',
     x           4*'S ', 'P ', 2*'B ', 'F ', 'CL', 'BR', 'I ',
     x           'SI', 4*'  '/
C
      IF (KIND .EQ. 12) PRINT 11
11      FORMAT (' MACROMODEL input file name: ',$)
      IF (KIND .EQ. 4)  print 1
    1   FORMAT(' MOPAC input file name: ',$)  
      IF (KIND .EQ. 15) PRINT 13
13      FORMAT (' CHEM3D input file name: ',$)
      IF (KIND .EQ. 16) PRINT 14
14      FORMAT (' MOLPAK input file name: ',$)     
      IF (KIND .EQ. 21) PRINT 180               ! KIND = 21, general format
180     FORMAT (' General format input file name with ',
     x          'atom symbol, x, y, z: ',$)
      IF (KIND .EQ. 34) PRINT 181
181     FORMAT (' General format input file name with ',
     x          'atomic number, x, y, z: ',$)
      READ (5,2) INFILE
2       FORMAT(A32)
      IF (KIND .EQ. 4) THEN
        PRINT 188
188     FORMAT ('   Enter (1/2) for (initial/final) MOPAC',
     x          ' format: ',$)
        READ (5,*) IF_TYPE
      ENDIF
      IF (KIND. NE. 15) GO TO 182
         PRINT 822
822      FORMAT (' Does the CHEM3D file contain atom sequence',
     x           ' numbers? [Y]:',$)
         READ (5,824) Q
824      FORMAT (A1)
         IF (Q .EQ. 'Y' .OR. Q .EQ. 'y' .OR. Q .EQ. ' ') THEN
            ISEQ = .TRUE.
         ELSE
            ISEQ = .FALSE.
         ENDIF        
182   IF (KIND .EQ. 21) THEN
        PRINT 183
183     FORMAT ('  Format for input line as (..A2....): ',$)
        GO TO 191
      ENDIF
      IF (KIND .EQ. 34) THEN
        PRINT 184
184     FORMAT ('  Format for input line as (..I2....): ',$)
        GO TO 191
      ENDIF
      GO TO 187
191     READ (5,186) INPUT_FORMAT
186     FORMAT (A60)
        print 388, input_format
388        format ('Input file format = ',a)
C----Open the input file 
187   OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
      IF (KIND .EQ. 16 .OR. KIND .EQ. 21 .OR. KIND .EQ. 34) THEN
        PRINT 15
15      FORMAT (' CHEM3D output file (cc1) name: ',$)  
      ELSE      
        print 3
    3   FORMAT(' MOLPAK output coordinate file name: ',$)
      ENDIF
      read (5,2) OUTFILE
C     OPEN THE OUTFILE
      OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
      IF (KIND .EQ. 16 .OR. KIND .EQ. 21 .OR. KIND .EQ. 34) GO TO 104
      GO TO 100
C----MOLPAK to CHEM3D
104     NA = 0
105     NA = NA + 1
        IF (KIND .EQ. 16) THEN
           name(na) = '     '
           READ (10,99,END=200) NAME(NA)(1:2), J, (XYZ(I,NA),I=1,3)   ! 7/3/03
99            FORMAT (5X,A2,I3,1X,3F10.6)                             ! 7/3/03
           IF (.NOT. (NAME(NA)(1:2) .EQ. 'BR' .OR. NAME(NA)(1:2) .EQ. 
     X         'Br')) NAME(NA)(2:2) = ' ' 
        ELSE
           IF (KIND .EQ. 21) READ (10,INPUT_FORMAT,END=200) 
     x             NAME(NA)(1:2), (XYZ(I,NA),I=1,3) 
           IF (KIND .EQ. 34) THEN 
                READ (10,INPUT_FORMAT,END=200) J, (XYZ(I,NA),I=1,3)
                NAME(NA)(1:2) = ATOMNAME(J)
           ENDIF
        ENDIF
        GO TO 105
200     NA = NA - 1
        CALL CONNECT (NA)
        WRITE (11,'(I3)') NA              ! number of atom on first line      
C----Determine CHEM3D atom types
      NATM = NA
      DO 1303 NA=1,NATM
      AT2 = NAME(NA)(1:2) 
C----How many connected?
        IF( NCON(NA) .NE. 0) GO TO 1375
        DO L=1,6
           IF (ICON(L,NA) .EQ. 0) GO TO 1375
           NCON(NA) = NCON(NA) + 1
        ENDDO
1375    IF (AT2(1:1) .EQ. 'H') THEN                    ! H
           K = 5
           GO TO 158
        ENDIF
        IF (AT2(1:1) .EQ. 'F') THEN                    ! F
           K = 11                                  
           GO TO 158
        ENDIF
        IF (AT2(1:1) .EQ. 'C') THEN
           GO TO (1001, 1001, 1003, 1004), NCON(NA)
1001       K = 4                                        ! alkyne C
           GO TO 158
1003       K = 2                                        ! alkene C
           GO TO 158
1004       K = 1                                        ! alkane C
           GO TO 158
        ENDIF
        IF (AT2 .EQ. 'SI' .OR. AT2 .EQ. 'Si') THEN  ! Si
           K = 19
           GO TO 158
        ENDIF
        IF (AT2(1:1) .EQ. 'S') THEN
           GO TO (1011, 1011, 1013, 1014), NCON(NA)
1011       K = 15                                       ! thio ether
           GO TO 158
1013       K = 17                                       ! suloxide S
           GO TO 158
1014       K = 18                                       ! sulfone S
        ENDIF
        IF (AT2(1:1) .EQ. 'O') THEN
           GO TO (1021, 1022), NCON(NA)
1021          N_ADJ = ICON(1,NA)           ! get ID of single connected atom
              K = 7                        ! default is carbonyl
              IF (NAME(N_ADJ)(1:1) .EQ. 'C') GO TO 158       
              IF (NCON(N_ADJ) .NE. 3) GO TO 158
                 N_OXY = 0                 ! count # O's connected to N_ADJ
                 DO I=1,3
                    IS = ICON(I,N_ADJ)
                    IF (NAME(IS)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
                 ENDDO
                 IF (N_OXY .NE. 2) GO TO 158
                    K = 47                              ! nitro O
                 GO TO 158
1022       K = 6                                        ! ether O
           GO TO 158
        ENDIF
        IF (AT2(1:1) .EQ. 'N') THEN
           GO TO (1031, 1032, 1033, 1034), NCON(NA)
1031       K = 10                                       ! nitrile N
           GO TO 158
1032       K = 40                                       ! enamine N
           GO TO 158
1033       K = 8                                        ! amine N default
              N_OXY = 0                ! count # of O's on this N
              DO I=1,3
                 N_ADJ = ICON(I,NA)                 
                 IF (NAME(N_ADJ)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
              ENDDO
              IF (N_OXY .EQ. 2) K = 46                  ! nitro N
           GO TO 158
1034       K = 39                                       ! ammonium N
           GO TO 158
        ENDIF
158     WRITE(11,159) AT2, (XYZ(J,NA),J=1,3), K,
     X             (ICON(J,NA),J=1,NCON(NA))
159     FORMAT (2X,A2,F11.6,2F12.6,7I5)
1303  CONTINUE
      RETURN          ! finished converting MOLPAK to CHEM3D file
C
100   IF (KIND .EQ. 12 .OR. KIND .EQ. 15) THEN  ! MACROMODEL or CHEM3D input files
        READ (10,'(I6)') NA           ! number of atoms on 1st line
        DO 18 K=1,NA 
          IF (KIND .EQ. 12) THEN
             READ (10,17,END=997) IT, (XYZ(I,K),I=1,3)  ! MACROMODEL
   17        FORMAT (I4,48X,3F12.6)
C---------Temporary fix for B which has MM code of -5
             IF (IT .NE. -5) GO TO 2222
                NAME(K) = 'B    '
                GO TO 18
2222         NAME(K)(1:2) = MMOD(IT)    ! convert from the macromodel
             NAME(K)(3:5) = '   '       ! atom code to element symbol
          ELSE
             IF (ISEQ) THEN
                READ (10,29,END=997) NAME(K)(1:2), (XYZ(I,K),I=1,3)  ! CHEM3D
   29           FORMAT (2x,A2,4x,3F12.6)           ! w atom sequence numbers
             ELSE
                READ (10,19,END=997) NAME(K)(1:2), (XYZ(I,K),I=1,3)  ! CHEM3D
   19           FORMAT (2x,A2,F11.6,2F12.6)        ! w/o atom sequence numbers
             ENDIF
          ENDIF
18      CONTINUE
      ELSE
        NA = 0
45      NA = NA + 1
        IF (IF_TYPE .EQ. 2) THEN
           READ (10,16,END=999) IK(NA),NAME(NA),(XYZ(I,NA),I=1,3)
16         FORMAT (I6,9X,A5,10X,3F10.4)   ! final mopac format
        ELSE
           READ (10,2031,END=999) NAME(NA),(XYZ(I,NA),I=1,3)
2031       FORMAT (1X,A5,3(F10.5,3X))    ! initial mopac format
        ENDIF
        GO TO 45
  999   NA = NA - 1
      ENDIF
C-----change C-H, N-H, O-H and/or N-O                 3/25/03
997   CALL ADJUST_2
      X=0.0
      Y=0.0
      Z=0.0
      DO J=1,NA
       X=X+XYZ(1,J)
       Y=Y+XYZ(2,J)
       Z=Z+XYZ(3,J)
      END DO
      X=X/NA
      Y=Y/NA
      Z=Z/NA
      DO J=1,NA
      XYZ(1,J)=XYZ(1,J)-X
      XYZ(2,J)=XYZ(2,J)-Y
      XYZ(3,J)=XYZ(3,J)-Z
      END DO 
      DO J=1,NA
        L = NPOT(J)
        IF (L .EQ. 64 .OR. L .EQ. 65) THEN                   ! 2/3/04
           WRITE (11,10) NAME(J)(1:2), L, (XYZ(I,J),I=1,3),  ! 2/3/04
     X                   L, G03_CHARGE, XMNDO_CHARGE         ! 2/4/04
10         FORMAT('ATOM ',A2,I3,1X,3F10.6,I5,2F10.6)  !  MOLPAK format
        ELSE                                                 ! 2/4/03
           WRITE (11,10) NAME(J)(1:2), J, (XYZ(I,J),I=1,3),  ! 2/3/04
     X                   L, G03_CHARGE, XMNDO_CHARGE
        ENDIF                                                ! 2/4/04
      ENDDO
       RETURN
       END
C
C---------------------------------------------------------------
C
      SUBROUTINE CHANGE (XYZF)   ! change C-H & N-H bonds
      DIMENSION XYZ(3,200),T(3,3),MN(200),MO(200),
     1          XYZCAR(3,200),ACAR(6),NH(200),NC(200),
     2          LK(3),M(3),IK(200),A(6),XYZF(3,200),NCKIND(200)
      CHARACTER*5 NAME
      character NO*1
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZT(3,200),
     x              ICON(6,200), NCON(200)
      COMMON /COM2/XYZ 
C
      print 91
91    FORMAT (' Do you want to change the nitro N-O distance? [Y/N] ',$)
      read (5,92) NO
92    format (a1)
      DO I=1,NA
       DO J=1,3
         XYZ(J,I)=XYZF(J,I)
       END DO
      END DO
                CALL BOND2 (NUMH,NH,NC,NCKIND)  
      IF (NO .EQ. 'Y' .or. no .eq. 'y') THEN
                CALL BONDOXY (NUMO,MO,MN)
      END IF
                CALL MATRXT (CELL,T)          
                CALL MTIMES1 (NA,CELL,XYZ,T,XYZCAR)
      DO J=1,NUMH
                NN=NH(J)               
                MC=NC(J)
                NKIND=NCKIND(J)
      IF (NAME(MC)(1:1).EQ.'C') THEN
        IF (NKIND .EQ. 2) THEN
                CALL CHBOND (NN,1.098,MC,NAME,XYZCAR)
          ELSE 
              IF (NKIND .EQ. 1) THEN
                CALL CHBOND (NN,1.084,MC,NAME,XYZCAR)
                ENDIF         
          ENDIF 
      ELSE
        IF (NAME(MC)(1:1).EQ.'N') THEN
                CALL CHBOND (NN,1.013,MC,NAME,XYZCAR)
        ENDIF
      ENDIF
      END DO

      IF (NO .eq. 'Y' .or. NO .eq. 'y') THEN
         DO J=1,NUMO
                NN=MO(J)               
                MC=MN(J)
                CALL CHBOND (NN,1.220,MC,NAME,XYZCAR)
         END DO
      END IF

      DO 110 J=1,3
          ACAR(J)=1.0
  110 CONTINUE
      DO 115 J=4,6
          ACAR(J)=90.0
  115 CONTINUE

                CALL MINV (T,3,DET,LK,M)
                CALL MTIMES2 (NA,XYZCAR,T,XYZ)
      DO I=1,NA
      XYZ(1,I)=XYZ(1,I)/CELL(1)
      XYZ(2,I)=XYZ(2,I)/CELL(2)
      XYZ(3,I)=XYZ(3,I)/CELL(3)
      END DO

C-----WRITE CHEMX FILE
C     WRITE (11,80) (CELL(J),J=1,3)
C  80 FORMAT (38X,3F8.3)
C     WRITE (11,70) (CELL(J+3),J=1,3)
C  70 FORMAT (21X,3F8.3)
C     WRITE (11,60) NA
C  60 FORMAT (I4)
C     WRITE (11,40)
C  40 FORMAT (60X)
C     DO 10 J=1,NA
C     WRITE (11,50) J, NAME(J),(XYZ(I,J),I=1,3)
C  50 FORMAT (I4,1X,A4,1X,3F10.5)
C  10 CONTINUE
      DO I=1,NA
       DO J=1,3
         XYZF(J,I)=XYZ(J,I)
       END DO
      END DO
       RETURN
       END
C----------------------------------------------------------------
      SUBROUTINE BOND2 (N,NH,NNC,NCKIND) 
      DIMENSION XYZ(3,200), NH(200),NNC(200),NCKIND(200)            
      INTEGER NA,BNUM
      CHARACTER*5 NAME
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZT(3,200),
     x              ICON(6,200), NCON(200)
      COMMON /COM2/XYZ 
      DA=COSD(CELL(4))
      DB=COSD(CELL(5))
      DG=COSD(CELL(6))
        N=1
      DO 30 IA=1,NA
      DO 35 JA=IA+1,NA
      IF ((NAME(IA)(1:1).EQ.'C'.AND.NAME(JA)(1:1).EQ.'H').OR.
     1    (NAME(IA)(1:1).EQ.'N'.AND.NAME(JA)(1:1).EQ.'H').OR.
     2    (NAME(IA)(1:1).EQ.'H'.AND.NAME(JA)(1:1).EQ.'C').OR.
     3    (NAME(IA)(1:1).EQ.'H'.AND.NAME(JA)(1:1).EQ.'N')) THEN
      DX=XYZ(1,IA)-XYZ(1,JA)
      DY=XYZ(2,IA)-XYZ(2,JA)
      DZ=XYZ(3,IA)-XYZ(3,JA)
      DD=SQRT(DX**2*CELL(1)**2+
     1DY**2*CELL(2)**2+
     2DZ**2*CELL(3)**2+2*DX*DY*CELL(1)*CELL(2)*DG+2*DX*DZ*CELL(1)*
     3CELL(3)*DB+2*DY*DZ*CELL(2)*CELL(3)*DA)
         IF (DD.LE.1.3)  THEN
           IF (NAME(IA)(1:1).EQ.'H') THEN
             NH(N) =IA            
             NNC(N)=JA
             CALL BONDNUM (JA,BNUM)
             IF (BNUM.EQ.4) THEN
              NCKIND(N)=2
             ELSE
              NCKIND(N)=1
             END IF
           ELSE
             NH(N) =JA
             NNC(N)=IA
             CALL BONDNUM (IA,BNUM)
             IF (BNUM.EQ.4) THEN
              NCKIND(N)=2
             ELSE
              NCKIND(N)=1
             END IF
           ENDIF
             N=N+1
         ENDIF
       ENDIF
   35 CONTINUE
   30 CONTINUE
      N=N-1
C     PRINT *,N
C     PRINT *,(NH(J),NNC(J),J=1,N)
       RETURN
       END
C------------------------------------------------------------
      SUBROUTINE BONDOXY (N,NO,NN)                ! 9/23/93
      DIMENSION XYZ(3,200), NO(200),NN(200)            
      INTEGER NA,BNUM
      CHARACTER*5 NAME
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZT(3,200),
     x              ICON(6,200), NCON(200)
      COMMON /COM2/XYZ 
      DA=COSD(CELL(4))
      DB=COSD(CELL(5))
      DG=COSD(CELL(6))
        N=1
      DO 30 IA=1,NA
      DO 35 JA=IA+1,NA
      IF ((NAME(IA)(1:1).EQ.'O'.AND.NAME(JA)(1:1).EQ.'N').OR.
     1    (NAME(IA)(1:1).EQ.'N'.AND.NAME(JA)(1:1).EQ.'O')) THEN
      DX=XYZ(1,IA)-XYZ(1,JA)
      DY=XYZ(2,IA)-XYZ(2,JA)
      DZ=XYZ(3,IA)-XYZ(3,JA)
      DD=SQRT(DX**2*CELL(1)**2+
     1DY**2*CELL(2)**2+
     2DZ**2*CELL(3)**2+2*DX*DY*CELL(1)*CELL(2)*DG+2*DX*DZ*CELL(1)*
     3CELL(3)*DB+2*DY*DZ*CELL(2)*CELL(3)*DA)
         IF (DD .LE. 1.25) THEN       ! leave alone any N-O's > 1.25 Angs 
           IF (NAME(IA)(1:1).EQ.'O') THEN
             CALL BONDNUM (IA,BNUM)                            ! 7/10/92
             IF (BNUM.EQ.1) THEN
               NO(N) =IA            
               NN(N) =JA
               N=N+1
             ENDIF
           ELSE
             CALL BONDNUM(JA,BNUM)
             IF (BNUM.EQ.1) THEN
             NO(N) =JA
             NN(N) =IA
             N=N+1
             ENDIF
           ENDIF
         ENDIF
       ENDIF
   35 CONTINUE
   30 CONTINUE
      N=N-1
C     PRINT *,N
C     PRINT *,(NO(J),NN(J),J=1,N)
       RETURN
       END
C------------------------------------------------------------
      SUBROUTINE CHBOND (N,BLEN,M,NAME,XYZC)
      DIMENSION XYZC(3,200)
      CHARACTER*5 NAME(200)
C----N = NUMBER OF HYDROGEN, M = NUMBER OF OTHER ATOM   
      DHX=XYZC(1,N)-XYZC(1,M)
      DHY=XYZC(2,N)-XYZC(2,M)
      DHZ=XYZC(3,N)-XYZC(3,M)
      CH=SQRT(DHX**2+DHY**2+DHZ**2)
      SCALE = BLEN/CH
C----CH 1.098 if sp3, NH 1.013 
      XYZC(1,N)=DHX*SCALE + XYZC(1,M)      
      XYZC(2,N)=DHY*SCALE + XYZC(2,M)      
      XYZC(3,N)=DHZ*SCALE + XYZC(3,M)      
      DHX=XYZC(1,N)-XYZC(1,M)
      DHY=XYZC(2,N)-XYZC(2,M)
      DHZ=XYZC(3,N)-XYZC(3,M)
      CHnew=SQRT(DHX**2+DHY**2+DHZ**2)
      PRINT 1234, NAME(M), NAME(N),CH, chnew
1234  format (1x,a5,'-',a5,': original =',F6.3,', new =',F6.3)  
      RETURN
      END

      real function SIND(x)
      SIND = SIN(x/57.29577951)
      return
      end

      real function cosd(x)
      if (abs(abs(x) - 90.0) .le. 0.000001) then
         cosd = 0.0
      else
        cosd = cos(x/57.29577951)
      endif
      return
      end

      real function ACOSD(x)
      ACOSD = ACOS(x)*57.29577951
      return 
      end

      SUBROUTINE MATRXT (A,T)
      DIMENSION T(3,3),A(6)
C     print 1238, A
C1238  format (' A in matrxt =',6E12.5)
      T(1,1)=1.0     
      T(1,2)=COSD(A(6))
C     print 1239, a(6), t(1,2)
C1239  format (' a(6), t(1,2) =',2E12.5)
      T(1,3)=COSD(A(5))

      T(2,1)=0.0
      T(2,2)=SIND(A(6))
      T(2,3)=(COSD(A(4))-COSD(A(5))*COSD(A(6)))/SIND(A(6))
 
      W=SQRT(1-COSD(A(4))**2-COSD(A(5))**2-COSD(A(6))**2+
     12*COSD(A(4))*COSD(A(5))*COSD(A(6)))
      T(3,1)=0.0
      T(3,2)=0.0
      T(3,3)=W/SIND(A(6))

C     PRINT 100,((T(I,J),j=1,3),i=1,3)
C 100 FORMAT(' T matrix in MATRXT = ',3E12.5/
C    x     2('                      ',3E12.5/))
      RETURN
      END

      SUBROUTINE MTIMES1 (N,A,XYZG,RV,XYZT) 
      DIMENSION XYZT(3,100),XYZG(3,100),RV(3,3),A(6)
      DO 10 K=1,N
      XYZT(1,K)=RV(1,1)*XYZG(1,K)*A(1)+RV(1,2)*XYZG(2,K)*A(2)
     1+RV(1,3)*XYZG(3,K)*A(3)
      XYZT(2,K)=RV(2,1)*XYZG(1,K)*A(1)+RV(2,2)*XYZG(2,K)*A(2)
     1+RV(2,3)*XYZG(3,K)*A(3)
      XYZT(3,K)=RV(3,1)*XYZG(1,K)*A(1)+RV(3,2)*XYZG(2,K)*A(2)
     1+RV(3,3)*XYZG(3,K)*A(3)
   10 CONTINUE
      RETURN
      END

      SUBROUTINE MTIMES2 (N,XYZG,RV,XYZT) 
      DIMENSION XYZT(3,4),XYZG(3,4),RV(3,3)
      DO 10 K=1,N
      XYZT(1,K)=RV(1,1)*XYZG(1,K)+RV(1,2)*XYZG(2,K)
     1+RV(1,3)*XYZG(3,K)
      XYZT(2,K)=RV(2,1)*XYZG(1,K)+RV(2,2)*XYZG(2,K)
     1+RV(2,3)*XYZG(3,K)
      XYZT(3,K)=RV(3,1)*XYZG(1,K)+RV(3,2)*XYZG(2,K)
     1+RV(3,3)*XYZG(3,K)
   10 CONTINUE
      RETURN
      END

      SUBROUTINE MINV(A,N,D,L,M)                                        MINV   5
C                                                                       MINV  15
C     ..................................................................MINV  20
C                                                                       MINV  25
C        SUBROUTINE MINV                                                MINV  30
C                                                                       MINV  35
C        PURPOSE                                                        MINV  40
C           INVERT A MATRIX                                             MINV  45
C                                                                       MINV  50
C        USAGE                                                          MINV  55
C           CALL MINV(A,N,D,L,M)                                        MINV  60
C                                                                       MINV  65
C        DESCRIPTION OF PARAMETERS                                      MINV  70
C           A - INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY  MINV  70
C               RESULTANT INVERSE.                                      MINV  80
C           N - ORDER OF MATRIX A                                       MINV  85
C           D - 0.0 IF SINGULAR, 1.0 IF ALL RIGHT.                      MINV  90
C           L - WORK VECTOR OF LENGTH N                                 MINV  95
C           M - WORK VECTOR OF LENGTH N                                 MINV 100
C                                                                       MINV 105
C        REMARKS                                                        MINV 110
C           MATRIX A MUST BE A GENERAL MATRIX                           MINV 115
C                                                                       MINV 120
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MINV 125
C           NONE                                                        MINV 130
C                                                                       MINV 135
C        METHOD                                                         MINV 140
C           THE STANDARD GAUSS-JORDAN METHOD IS USED.                   MINV 145
C                                                                       MINV 150
C     ..................................................................MINV 155
C                                                                       MINV 160
      DIMENSION A(1),L(1),M(1)                                          MINV 165
C                                                                       MINV 170
C        ...............................................................MINV 175
C                                                                       MINV 180
C        SEARCH FOR LARGEST ELEMENT                                     MINV 185
C                                                                       MINV 190
      D=1.0                                                             MINV 195
      NK=-N                                                             MINV 200
      DO 590 K=1,N                                                      MINV 205
      NK=NK+N                                                           MINV 210
      L(K)=K                                                            MINV 215
      M(K)=K                                                            MINV 220
      KK=NK+K                                                           MINV 225
      BIGA=A(KK)                                                        MINV 230
      DO 510 J=K,N                                                      MINV 235
      IZ=N*(J-1)                                                        MINV 240
      DO 510 I=K,N                                                      MINV 245
      IJ=IZ+I                                                           MINV 250
  500 IF(ABS(BIGA)-ABS(A(IJ))) 505,510,510                              MINV 255
  505 BIGA=A(IJ)                                                        MINV 260
      L(K)=I                                                            MINV 265
      M(K)=J                                                            MINV 270
  510 CONTINUE                                                          MINV 270
C                                                                       MINV 280
C        INTERCHANGE ROWS                                               MINV 285
C                                                                       MINV 290
      J=L(K)                                                            MINV 295
      IF(J-K) 525,525,515                                               MINV 300
  515 KI=K-N                                                            MINV 305
      DO 520 I=1,N                                                      MINV 310
      KI=KI+N                                                           MINV 315
      HOLD=-A(KI)                                                       MINV 320
      JI=KI-K+J                                                         MINV 325
      A(KI)=A(JI)                                                       MINV 330
  520 A(JI) =HOLD                                                       MINV 335
C                                                                       MINV 340
C        INTERCHANGE COLUMNS                                            MINV 345
C                                                                       MINV 350
  525 I=M(K)                                                            MINV 355
      IF(I-K) 540,540,530                                               MINV 360
  530 JP=N*(I-1)                                                        MINV 365
      DO 535 J=1,N                                                      MINV 370
      JK=NK+J                                                           MINV 375
      JI=JP+J                                                           MINV 380
      HOLD=-A(JK)                                                       MINV 385
      A(JK)=A(JI)                                                       MINV 390
  535 A(JI) =HOLD                                                       MINV 395
C                                                                       MINV 400
C        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS        MINV 405
C        CONTAINED IN BIGA)                                             MINV 410
C                                                                       MINV 415
  540 IF(BIGA) 550,545,550                                              MINV 420
  545 D=0.0                                                             MINV 425
      RETURN                                                            MINV 430
  550 DO 560 I=1,N                                                      MINV 435
      IF(I-K) 555,560,555                                               MINV 440
  555 IK=NK+I                                                           MINV 445
      A(IK)=A(IK)/(-BIGA)                                               MINV 450
  560 CONTINUE                                                          MINV 455
C                                                                       MINV 460
C        REDUCE MATRIX                                                  MINV 465
C                                                                       MINV 470
      DO 575 I=1,N                                                      MINV 470
      IK=NK+I                                                           MINV 480
      IJ=I-N                                                            MINV 485
      DO 575 J=1,N                                                      MINV 490
      IJ=IJ+N                                                           MINV 495
      IF(I-K) 565,575,565                                               MINV 500
  565 IF(J-K) 570,575,570                                               MINV 505
  570 KJ=IJ-I+K                                                         MINV 510
      A(IJ)=A(IK)*A(KJ)+A(IJ)                                           MINV 515
  575 CONTINUE                                                          MINV 520
C                                                                       MINV 525
C        DIVIDE ROW BY PIVOT                                            MINV 530
C                                                                       MINV 535
      KJ=K-N                                                            MINV 540
      DO 585 J=1,N                                                      MINV 545
      KJ=KJ+N                                                           MINV 550
      IF(J-K) 580,585,580                                               MINV 555
  580 A(KJ)=A(KJ)/BIGA                                                  MINV 560
  585 CONTINUE                                                          MINV 565
C                                                                       MINV 570
C        REPLACE PIVOT BY RECIPROCAL                                    MINV 575
C                                                                       MINV 580
      A(KK)=1.0/BIGA                                                    MINV 585
  590 CONTINUE                                                          MINV 590
C                                                                       MINV 595
C        FINAL ROW AND COLUMN INTERCHANGE                               MINV 600
C                                                                       MINV 605
      K=N                                                               MINV 610
  595 K=(K-1)                                                           MINV 615
      IF(K) 630,630,600                                                 MINV 620
  600 I=L(K)                                                            MINV 625
      IF(I-K) 615,615,605                                               MINV 630
  605 JQ=N*(K-1)                                                        MINV 635
      JR=N*(I-1)                                                        MINV 640
      DO 610 J=1,N                                                      MINV 645
      JK=JQ+J                                                           MINV 650
      HOLD=A(JK)                                                        MINV 655
      JI=JR+J                                                           MINV 660
      A(JK)=-A(JI)                                                      MINV 665
  610 A(JI) =HOLD                                                       MINV 670
  615 J=M(K)                                                            MINV 670
      IF(J-K) 595,595,620                                               MINV 680
  620 KI=K-N                                                            MINV 685
      DO 625 I=1,N                                                      MINV 690
      KI=KI+N                                                           MINV 695
      HOLD=A(KI)                                                        MINV 700
      JI=KI-K+J                                                         MINV 705
      A(KI)=-A(JI)                                                      MINV 710
  625 A(JI) =HOLD                                                       MINV 715
      GO TO 595                                                         MINV 720
  630 RETURN                                                            MINV 725
      END                                                               MINV 730
C
C----------------------------------------------------------------------
C
C     Convert from Mopac, g94 intermediate coord and free format files
C      to CHEMX or fractional coord files
C
      SUBROUTINE MOPAC2CHEMX (KIND)
      CHARACTER*32 INFILE, OUTFILE  
      CHARACTER NAME*5, T_CHAR(200)*1, LINE*65, LINE1*25, L(65)*1,
     x          LINE3*40
      CHARACTER BLANK*25, ATOM_NAME(9)*1
      DIMENSION IK(200)  
      COMMON /COM1/ NA, A(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)  
      EQUIVALENCE (LINE3, LINE, L), (L(41), LINE1)  
      DATA ATOM_NAME /'H', 4*' ', 'C', 'N', 'O', 'F'/
      DATA BLANK /'                         '/
      DATA T_CHAR /200*' '/
C
      IF (KIND .EQ. 14) PRINT 62
62       FORMAT (' Name of free format file with at no, x, y, z: ',$)
      IF (KIND .EQ. 13 .OR. KIND .EQ. 26) THEN
         PRINT 63
63       FORMAT (' G03 intermediate coord file name: ',$)
      ELSE
        print 1
    1   FORMAT(' MOPAC input file name: ',$) 
      ENDIF
      read (5,2) INFILE
    2 FORMAT(A32)
C     OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
      IF (KIND .EQ. 2. .OR. KIND .EQ. 26) THEN
         print 3
    3    FORMAT(' CHEMX output file name: ',$)
      ELSE
         PRINT 33
33       FORMAT (' Fract coord output file name: ',$)
      ENDIF
      read (5,2) OUTFILE
C     OPEN THE OUTFILE
      OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
      IF (KIND .EQ. 13 .OR. KIND .EQ. 14 .OR. KIND .EQ. 26)  GO TO 73
      PRINT 181
181   format (' What xyz would you like to read?.... '/
     1        '   1 = read original xyz; 2 = read final xyz '/
     2        '   Enter 1 or 2: ',$)
      read (5,*) KKIND
73    NA=0
      DO 945 JN=1,100
      NA=NA+1
      IF (KIND .EQ. 2 .OR. KIND .EQ. 6) GO TO 961
         IF (KIND .EQ. 13 .OR. KIND .EQ. 26) THEN
            READ (10,955,END=999) M, (XYZ(I,JN),I=1,3)
955         FORMAT (11X,I5,7X,3F12.6)
         ELSE
            READ (10,*,END=999) M, (XYZ(I,JN),I=1,3)   ! free format 
         ENDIF
         NAME(JN) = '     '
         NAME(JN)(1:1) = ATOM_NAME(M)
         GO TO 945
961   IF (KKIND.EQ.1) THEN 
         READ (UNIT=10,END=999,FMT=15) NAME(JN),(XYZ(I,JN),I=1,3)
      ELSE
         READ (10,16,END=999) IK(JN), T_CHAR(JN), NAME(JN),
     X                        (XYZ(I,JN),I=1,3)
      END IF
945   CONTINUE
   15 FORMAT (15x,A5,3F10.5)
   16 FORMAT (I6,8X,A1,A5,10X,3F10.4)
  999 NA=NA-1
      IF (KIND .EQ. 2 .OR. KIND .EQ. 26) THEN
           WRITE (11,10) (A(J),J=1,3)
   10      FORMAT (38X,3F8.3)
           WRITE (11,20) (A(J),J=4,6)
   20      FORMAT (21X,3F8.3)
           WRITE (11,25) NA
   25      FORMAT (I4/)
      ELSE
           WRITE (11,998) NA
998        FORMAT ('REFERENCE STRUCTURE =     1   A,B,C =   1.000',
     X             '   1.000   1.000'/
     X             '  ALPHA,BETA,GAMMA =  90.000  90.000  90.000',
     X             '    SPGR =  1 P1'/
     x             I3,'   0 CODON=          SYMOPS=    1'/
     x             '  0     RFAC=   .1 ERRFLAG=0 (C-C)ESD=0')
C----For frac coord, need to establish connectivity
      CALL CONNECT (NA) 
      ENDIF
C----May need to alter atom name from MOPAC if it has 2 characters
      LINE1 = BLANK
      DO 700 JN=1,NA
         IF (KIND .EQ. 2 .OR. KIND .EQ. 26) GO TO 29 
            LINE1 = BLANK
            LL = NCON(JN)
            WRITE (LINE1,28) (ICON(M,JN),M=1,LL)
28          FORMAT (I5,5I4)
29       IF (T_CHAR(JN) .EQ. ' ') THEN 
            IF (JN .LE. 9) THEN
                 WRITE (LINE3,30) JN, NAME(JN), JN, (XYZ(I,JN),I=1,3)
30               FORMAT (I4,1X,A1,I1,3X,3F10.5)  
            ELSE
                WRITE (LINE3,330) JN, NAME(JN), JN, (XYZ(I,JN),I=1,3)
330             FORMAT (I4,1X,A1,I2,2X,3F10.5)
            ENDIF
         ELSE
            IF (JN .LE. 9) THEN
               WRITE (LINE3,332) JN, T_CHAR(JN), NAME(JN), JN,
     X                          (XYZ(I,JN),I=1,3)
332            FORMAT (I4,1X,2A1,I1,2X,3F10.5)
            ELSE 
               WRITE (LINE3,333) JN, T_CHAR(JN), NAME(JN), JN,
     X                           (XYZ(I,JN),I=1,3)
333            FORMAT (I4,1X,2A1,I2,1X,3F10.5)
            ENDIF
         ENDIF
       WRITE (11,699) LINE
699    FORMAT (A65)
700    CONTINUE
       RETURN
       END
C
C-----------------------------------------------------------------
C     Read an initial or final format mopac coordinate file and create 
C      a new mopac file for geometry optimization or mndo/esp charge 
C      calculation or a g03 file for 6-31g* esp/chelpg charge calculation
C
      SUBROUTINE NEWMOPAC_OR_G03 (KIND)
      LOGICAL CONNECT_FLAG
      CHARACTER NAME*5, TITLE*72, A_NAME(64)*2
      CHARACTER*32 INFILE,OUTFILE  
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              N_NEAR(200), N_ATOMS(6,200)    
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      DATA A_NAME /9*'C ', '  ', 2*'C ', '  ', 'C ',
     x             7*'O ', '  ', 'O ',
     x            14*'N ', 2*'  ', 'N ',
     x             5*'H ', 2*'  ', 'H ',
     x             4*'S ', 'P ', 2*'B ', 'F ', 'CL', 'BR', 'I ',
     x               'SI', 4*'  '/                                    
      DIMENSION IK(200)
      DIMENSION LIST(3)
c
      IF (KIND .EQ. 27) THEN
         PRINT 1000
1000     FORMAT (' Give name of MACROMODEL input file: ',$)
      ELSE
         PRINT 1 
    1    FORMAT(' Input MOPAC coordinate file name: ',$) 
      ENDIF
      READ (5,2) INFILE
    2 FORMAT(A32)
      IF (KIND .NE. 27) THEN
         PRINT 18                  ! what kind of MOPAC file?
18       FORMAT ('   Enter (1/2) for (initial/final) MOPAC',
     x           ' format: ',$)
         READ (5,*) IF_TYPE
      ENDIF
C     OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
C
      IF (KIND .EQ. 27) READ (10,*)    ! skip 1st line of MM file
C----Output file...mopac or g94
      IF (KIND .EQ. 25 .OR. KIND .EQ. 27) THEN
         PRINT 44                   ! file for g03 ESP/chelpg charge calcn
44       FORMAT (' Create a g03 input file named ',
     x           '631gstar_chelpg.com')
         OPEN (UNIT=11, FILE='631gstar_chelpg.com',
     x         STATUS='UNKNOWN')
      ELSE 
         print 3                                    ! mndo
    3    FORMAT(' New MOPAC file name: ',$)
         read (5,2) OUTFILE
         OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
      ENDIF
         print 6
6     FORMAT (' Title: ',$)
      read (5,4) TITLE
    4 FORMAT(A72)
      NA=0
      DO 103 JN=1,200 
        NA=NA+1
        IF (KIND .NE. 27) GO TO 45
           READ (10,43,END=999) IT, (XYZ(I,JN),I=1,3)  ! MACROMODEL input
43         FORMAT (I4,48X,3F12.6)                      ! file
           NAME(JN)(1:2) = A_NAME(IT)    ! MM atom type to element symbol            
           NAME(JN)(3:5) = '   '
           GO TO 103
45      IF (IF_TYPE .NE. 1) THEN
           READ (10,16,END=999) IK(JN), NAME(JN), (XYZ(I,JN),I=1,3)
16         FORMAT (I6,8X,A5,11X,3F10.4)              ! final MOPAC format
        ELSE
           READ (10,2031,END=999) NAME(JN), (XYZ(I,JN),I=1,3)
2031       FORMAT (1X,A5,3(F10.5,3X))              ! initial mopac format
        ENDIF
        IF (NAME(JN)(1:1) .EQ. ' ') THEN
           NAME(JN)(1:2) = NAME(JN)(2:3)
        ENDIF
103   CONTINUE 
  999 NA=NA-1 
      CALL ADJUST_2
      IF (KIND .EQ. 3) THEN                 ! AM1 geometry optimization
         WRITE (11,10)
   10    FORMAT(' XYZ NOINTER MMOK AM1 PRECISE T=240H ')
      ENDIF
      IF (KIND .EQ. 7) THEN                 ! mndo ESP charge calcn
          WRITE (11,100)
100       FORMAT(1X,' MNDO ESP 1SCF GNORM=0.1 T=240H SCFCRT=1.0D-4')
      ENDIF
      IF (KIND .EQ. 25 .OR. KIND .EQ. 27) THEN  ! g03 631g* esp/chelpg charge calcn
          WRITE (11,108)
108       FORMAT ('$RunGauss'/
     x            '%chk=631gstar_chelpg.chk'/
     x            '#p 6-31G* scf=direct pop=CHELPG geom=coord'/)
      ENDIF
      IF (KIND .EQ. 3 .OR. KIND .EQ. 7) THEN
         WRITE (11,15) TITLE
   15    FORMAT (2X,A72)
         WRITE (11,20)
   20 FORMAT (1X,'--------------------------------------------------')
      ELSE
         WRITE (11,28) TITLE
28       FORMAT (1X,A72//'   0   1')
      ENDIF
      IF (KIND .EQ. 3 .OR. KIND .EQ. 7) THEN             ! mndo
         WRITE (11,25)(NAME(K),(XYZ(I,K),I=1,3),K=1,NA)  
   25    FORMAT (1X,A2,3X,F10.5,2X,'1',F10.5,2X,'1',F10.5,2X,'1')
         WRITE (11,35)
   35    FORMAT (1X,'0'//)
      ELSE
         WRITE (11,125) (NAME(K),(XYZ(I,K),I=1,3),K=1,NA)
  125    FORMAT (3X,A5,3F12.6) 
         WRITE (11,*)
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------
C
C     Convert from MACROMODEL to CHEMX or fractional coord files
C
      SUBROUTINE MMOD2OTHER (KIND)
C
      DIMENSION XYZ(3,100), A(6), NCON(100), ICON(6,100)
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      DATA A /3*1.0, 3*90.0/
      CHARACTER NAME(64)*2, A_CHAR(100)*2, LINE*65, LINE1*25, L(65)*1,
     x          LINE3*40
      CHARACTER BLANK*25
      DATA BLANK /'                         '/
      EQUIVALENCE (LINE3, LINE, L), (L(41), LINE1)
      CHARACTER*32 INFILE,OUTFILE
      DATA NAME /9*'C ', '  ', 2*'C ', '  ', 'C ',
     x           7*'O ', '  ', 'O ',
     x           14*'N ', 2*'  ', 'N ',
     x           5*'H ', 2*'  ', 'H ',
     x           4*'S ', 'P ', 2*'B ', 'F ', 'CL', 'BR', 'I ',
     x           'SI', 4*'  '/
C
      print 1
    1 FORMAT(' MACROMODEL input file name: ',$) 
      read (5,2) INFILE
    2 FORMAT(A32)
C     OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
      IF (KIND .EQ. 8) THEN
         print 3
    3    FORMAT(' CHEMX output file name: ',$)
      ELSE
         PRINT 33
33       FORMAT (' Fract coord output file name: ',$)
      ENDIF
      read (5,2) OUTFILE
C     OPEN THE OUTFILE
      OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
C----MM file...atom type, connectivity info, x, y, z;
C     need to figure NCON
      READ (10,*)        ! skip 1st line....a title
      NA=0
      DO 995 JN=1,100
        NA=NA+1
        READ (10,15,END=999) IT, (ICON(I,JN),I=1,6), (XYZ(I,JN),I=1,3)
   15   FORMAT (I4,6(I6,2x),3F12.6)
        A_CHAR(NA) = NAME(IT)    ! MM atom type to element symbol
        NCON(NA) = 0       
        DO 990 I=1,6
           IF (ICON(I,NA) .EQ. 0) GO TO 995
           NCON(NA) = NCON(NA) + 1
990     CONTINUE  
995   CONTINUE
  999 NA=NA-1
C
      IF (KIND .EQ. 8) THEN
           WRITE (11,10) (A(J),J=1,3)
   10      FORMAT (38X,3F8.3)
           WRITE (11,20) (A(J),J=4,6)
   20      FORMAT (21X,3F8.3)
           WRITE (11,25) NA
   25      FORMAT (I4/)
      ELSE
           WRITE (11,998) NA
998        FORMAT ('REFERENCE STRUCTURE =     1   A,B,C =   1.000',
     X             '   1.000   1.000'/
     X             '  ALPHA,BETA,GAMMA =  90.000  90.000  90.000',
     X             '    SPGR =  1 P1'/
     x             I3,'   0 CODON=          SYMOPS=    1'/
     x             '  0     RFAC=   .1 ERRFLAG=0 (C-C)ESD=0')
      ENDIF
C----May need to alter atom name from MOPAC if it has 2 characters
      LINE1 = BLANK
      DO 700 JN=1,NA
         IF (KIND .EQ. 8) GO TO 29 
            LINE1 = BLANK
            LL = NCON(JN)
            WRITE (LINE1,28) (ICON(M,JN),M=1,LL)
28          FORMAT (I5,5I4)
29       IF (A_CHAR(JN)(2:2) .EQ. ' ') THEN 
            IF (JN .LE. 9) THEN
                 WRITE (LINE3,30) JN, A_CHAR(JN)(1:1), JN,
     x                            (XYZ(I,JN),I=1,3)
30               FORMAT (I4,1X,A1,I1,3X,3F10.5)  
            ELSE
                WRITE (LINE3,330) JN, A_CHAR(JN)(1:1), JN, 
     x                            (XYZ(I,JN),I=1,3)
330             FORMAT (I4,1X,A1,I2,2X,3F10.5)
            ENDIF
         ELSE
            IF (JN .LE. 9) THEN
               WRITE (LINE3,332) JN, A_CHAR(JN), JN, (XYZ(I,JN),I=1,3)
332            FORMAT (I4,1X,A2,I1,2X,3F10.5)
            ELSE 
               WRITE (LINE3,333) JN, A_CHAR(JN), JN, (XYZ(I,JN),I=1,3)
333            FORMAT (I4,1X,A2,I2,1X,3F10.5)
            ENDIF
         ENDIF
       WRITE (11,699) LINE
699    FORMAT (A65)
700    CONTINUE
       RETURN
       END
C
C----------------------------------------------------------------------
C
C     Convert from CHEM3D cart coord 1 format to
C       fractional coord file or generic atom, x, y,z file
C       for g94
C
      SUBROUTINE CHEM3D2CC (KIND)
C
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      DIMENSION XYZ(3,100), A(6), NCON(100), ICON(6,100)
      DATA A /3*1.0, 3*90.0/
      CHARACTER NAME(64)*2, A_CHAR(100)*2, LINE*65, LINE1*25, L(65)*1,
     x          LINE3*40
      CHARACTER BLANK*25, SERIAL*1
      DATA BLANK /'                         '/
      EQUIVALENCE (LINE3, LINE, L), (L(41), LINE1)
      CHARACTER*32 INFILE,OUTFILE
C
      print 1
    1 FORMAT(' CHEM3D cart coord 1 input file name: ',$)
      read (5,2) INFILE
    2 FORMAT(A32)
C     OPEN THE INFILE
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
      PRINT 5
5       FORMAT ('  Does the CHEM3D file have atom serial ',
     x          'numbers (Y/N) [Y]: ',$)
      READ (5,'(A1)') SERIAL
      IF (SERIAL .EQ. 'y' .OR. SERIAL .EQ. 'Y'
     x       .OR. SERIAL .EQ. ' ') THEN
         IS = 1         ! serial numbers
      ELSE
         IS = 2         ! no serial numbers
      ENDIF
      IF (KIND .EQ. 10) THEN
        PRINT 33
33      FORMAT (' Fract coord output file name: ',$)
      ELSE
        PRINT 34
34      FORMAT (' G03 output file name: ',$)
      ENDIF
      READ (5,2) OUTFILE
C     OPEN THE OUTFILE
      OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN')
C
C----Atoms and connectivity from CHEM3D file
      READ (10,20) NA      ! first line contains number of atoms
20    FORMAT (I3)
      IF (KIND .EQ. 10) THEN
C----Write header for output file
           WRITE (11,998) NA
998        FORMAT ('REFERENCE STRUCTURE =     1   A,B,C =   1.000',
     X             '   1.000   1.000'/
     X             '  ALPHA,BETA,GAMMA =  90.000  90.000  90.000',
     X             '    SPGR =  1 P1'/
     x             I3,'   0 CODON=          SYMOPS=    1'/
     x             '  0     RFAC=   .1 ERRFLAG=0 (C-C)ESD=0')
      ENDIF
C
      DO 30 I=1,NA
      GO TO (131, 132), IS
131     READ (10,40) A_CHAR(I), (XYZ(M,I),M=1,3), (ICON(M,I),M=1,6)
40      FORMAT (1X,A2,5X,3F12.5,5X,6I5)      ! with serial numbers
        GO TO 47
132     READ (10,140) A_CHAR(I), (XYZ(M,I),M=1,3), (ICON(M,I),M=1,6)
140     FORMAT (1X,A2,3F12.5,5X,6I5)         ! without serial numbers
47      IF (KIND .EQ. 10) GO TO 42
           IF (I .LE. 9) THEN
              WRITE (11,43) A_CHAR(I), I, (XYZ(M,I),M=1,3)
43            FORMAT (1X,A2,'_',I1,3F12.5)
           ELSE
              WRITE (11,44) A_CHAR(I), I, (XYZ(M,I),M=1,3)
44            FORMAT (1X,A2,'_',I2,3F12.5)
           ENDIF
        GO TO 30
42      IF (A_CHAR(I)(1:1) .NE. ' ') GO TO 41      ! atom id's could be ' X';
          A_CHAR(I)(1:1) = A_CHAR(I)(2:2)          ! change to 'X '
          A_CHAR(I)(2:2) = ' '
41      NCON(I) = 0   ! figure out the NCON values
        DO 45 K=1,6
           IF (ICON(K,I) .EQ. 0) GO TO 30
           NCON(I) = NCON(I) + 1
45      CONTINUE
30    CONTINUE
C----Make rest of output file
      IF (KIND .EQ. 11) RETURN
      DO 700 JN=1,NA
         LINE1 = BLANK
         LL = NCON(JN)
         WRITE (LINE1,28) (ICON(M,JN),M=1,LL)
28       FORMAT (I5,5I4)
         IF (A_CHAR(JN)(2:2) .EQ. ' ') THEN
            IF (JN .LE. 9) THEN
                 WRITE (LINE3,31) JN, A_CHAR(JN)(1:1), JN,
     x                            (XYZ(I,JN),I=1,3)
31               FORMAT (I4,1X,A1,I1,3X,3F10.5)
            ELSE
                WRITE (LINE3,330) JN, A_CHAR(JN)(1:1), JN,
     x                            (XYZ(I,JN),I=1,3)
330             FORMAT (I4,1X,A1,I2,2X,3F10.5)
            ENDIF
         ELSE
            IF (JN .LE. 9) THEN
               WRITE (LINE3,332) JN, A_CHAR(JN), JN, (XYZ(I,JN),I=1,3)
332            FORMAT (I4,1X,A2,I1,2X,3F10.5)
            ELSE
               WRITE (LINE3,333) JN, A_CHAR(JN), JN, (XYZ(I,JN),I=1,3)
333            FORMAT (I4,1X,A2,I2,1X,3F10.5)
            ENDIF
         ENDIF
       WRITE (11,699) LINE
699    FORMAT (A65)
700    CONTINUE
       RETURN
       END
C
C-----------------------------------------------------------------------
C 
C     Convert from CHEMX to MOLPAK format...can switch axes and molecule
C       is moved to the origin
C
      SUBROUTINE CHEMX2MOLPAK
      DIMENSION IK(200), T(3,3)
      DIMENSION XYZ_CART(3,200), IAXIS(3), VECTOR(3)
      CHARACTER NAME*5
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              N_NEAR(200), N_ATOMS(6,200)
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE                                          
      CHARACTER OUTFILE2(32)*1, INFILE2(32)*1, AX_ORDER*3 
      CHARACTER*32 INFILE, OUTFILE, MOLPAKXYZ, TRANS_V1
      DATA MOLPAKXYZ /'molpak.xyz                      '/
      EQUIVALENCE (OUTFILE,OUTFILE2), (INFILE, INFILE2)
      DATA OUTFILE2 /32*' '/
      DATA  IAXIS /1, 2, 3/, VECTOR /3*0.0/
C
      PRINT 1
 1    FORMAT(' CHEMX input file name: ',$) 
      READ (5,2) INFILE
    2 FORMAT(A32)
C----OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
C----Output file
      PRINT 11
11    FORMAT (' MOLPAK output coordinate file name [molpak.xyz]: ',$)
      READ (5,2) OUTFILE
      IF (OUTFILE2(1) .EQ. ' ') OUTFILE = MOLPAKXYZ     
72    OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
C----Read cell params from lines 1 and 2 or CHEMX file
      READ (10,770) CELL
770   FORMAT (38X,3F8.3/21X,3F8.3)
      READ (10,*)         ! skip 2 lines
      READ (10,*)                       
      NA=0
777   NA=NA+1
         READ (10,30,END=999) IK(NA), NAME(NA)(1:1), (XYZ(I,NA),I=1,3)
30       FORMAT (I4,1X,A1,4X,3F10.6)
      GO TO 777
  999 NA=NA-1
C----Conversion from fractional to cartesian coordinates
          CALL MATRXT (CELL, T)
          CALL MTIMES1 (NA, CELL, XYZ, T, XYZ_CART)
          DO 779 I=1,NA
             DO 779 J=1,3
                XYZ(J,I) = XYZ_CART(J,I)
779       CONTINUE
C-----Should the C-H, N-H, O-H and/or N-O vectors be adjusted?   3/25/03
997   CALL ADJUST_2
      X=0.0                     ! move center of molecule
      Y=0.0                     ! to 0, 0, 0
      Z=0.0
      DO J=1,NA
       X=X+XYZ(1,J)
       Y=Y+XYZ(2,J)
       Z=Z+XYZ(3,J)
      END DO
      X=X/NA
      Y=Y/NA
      Z=Z/NA
      DO J=1,NA
        XYZ(1,J)=XYZ(1,J)-X
        XYZ(2,J)=XYZ(2,J)-Y
        XYZ(3,J)=XYZ(3,J)-Z
      ENDDO
C----Ask user if the axes should be rearranged
      PRINT 800
800   FORMAT (' Should the axes be rearranged... possible answers ='/
     x        '   NO, no, 123 = xyz, 213 = yxz, 312 = zxy, etc]:', $)
      READ (5,'(A3)') AX_ORDER
        IF (AX_ORDER(1:1) .EQ. 'n' .OR. AX_ORDER(1:1) .EQ. 'N' .OR.
     X      AX_ORDER(1:1) .EQ. ' ') GO TO 820
        READ (AX_ORDER,'(3I1)') IAXIS
820   PRINT 806
806   FORMAT (' Cartesian translation vector wrt rearranged axes ',
     x'[possible answers ='/
     x        '   none = return and x, y, z components:', $)
      READ (5,2) TRANS_V1
      IF (TRANS_V1(1:1) .EQ. ' ') GO TO 821 
      READ (TRANS_V1,*) VECTOR
C----Output MOLPAK format and possibly rearrange axes
821   DO 810 I=1,NA
         L = NPOT(I)
         IF ( L.EQ.49) THEN     ! 11-12-02 for OX of nitrocubane
          WRITE (11,8111) NAME(I)(1:2), I,
     x         ((XYZ(IAXIS(J),I) - VECTOR(J)),J=1,3), L,
     X         G03_CHARGE, XMNDO_CHARGE
         ELSE 
          WRITE (11,811) NAME(I)(1:1), I, 
     x         ((XYZ(IAXIS(J),I) - VECTOR(J)),J=1,3), L,
     X         G03_CHARGE, XMNDO_CHARGE 
         END IF
811      FORMAT('ATOM ',A1,I4,1X,3F10.6,I5,2F10.6)     ! MOLPAK format
8111     FORMAT('ATOM ',A2,I3,1X,3F10.6,I5,2F10.6)     ! MOLPAK format
810   CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------------
      SUBROUTINE BONDNUM (IP,N)
      CHARACTER*5 NAME
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZT(3,200),
     x              ICON(6,200), NCON(200)
      COMMON /COM2/XYZ(3,200) 
C
      DA=cosd(CELL(4))
      DB=cosd(CELL(5))
      DG=cosd(CELL(6))
        N=0
      DO 30 IA=1,NA
      IF (IP .EQ. IA) GO TO 30
      DX=(XYZ(1,IP)-XYZ(1,IA))*CELL(1)
      DY=(XYZ(2,IP)-XYZ(2,IA))*CELL(2)
      DZ=(XYZ(3,IP)-XYZ(3,IA))*CELL(3)
C     PRINT 996, ip, (xyz(i,ip),i=1,3), ia, (xyz(i,ia),i=1,3) 
C 996 FORMAT (' In BONDNUM...IP, XYZ, IA, XYZ =',
C    X        2(I3,3e12.5))
C     PRINT 997, ip, ia, dx, dy, dz
C997  FORMAT (' In BONDNUM...IP, IA, DX, DY, DZ =',2I4,3e12.5)
      A = DX**2
      B = DY**2
      C = DZ**2
      D =  2.0*DX*DY*DG
      E =  2.0*DX*DZ*DB
      F =  2.0*DY*DZ*DA
      DD = SQRT(A + B + C + D + E + F)
C     print 999, ip, ia, a, b, c, d, e, f, dd
C999  format (' In BONDNUM, IP, IA, A-F, DD =',2i4,7E12.5)
       IF (DD .LE. 1.65) N = N + 1
C     PRINT 1234, IP, NAME(IP),NAME(IA),N,DD
C1234 format (' In BONDNUM...IP, NAME(IP), NAME(IA), N, DD =',
C    x        I4,2(1x,a5),i3,f8.3)
   30 CONTINUE
       RETURN
       END
C
C--------------------------------------------------------------------
C
C----CONNECT...build ICON array
      SUBROUTINE CONNECT (NATOMS)
      CHARACTER ID*5
      COMMON /COM1/ NA, CELL(6), ID(200), XYZO(3,200),
     x              ICON(6,200), NCON(200)
C
      DO 100 I=1,NATOMS
         IF (ID(I)(1:1) .EQ. 'X') GO TO 100
         NCON(I) = 0
         DO 90 J=1,NATOMS
            DMAX = 1.63                       ! 7/26/06 
            IF (ID(J)(1:1) .EQ. 'X') GO TO 90
            IF (I .EQ. J) GO TO 90
            IF (ID(I)(1:1) .EQ. 'H' .AND. ID(J)(1:1) .EQ. 'H')
     X                         GO TO 90
            IF (ID(I)(1:1) .EQ. 'F' .AND. ID(J)(1:1) .EQ. 'F')
     X                         GO TO 90
        IF ((ID(I)(1:1) .EQ. 'H' .AND. ID(J)(1:1) .EQ. 'O') .OR.
     X      (ID(I)(1:1) .EQ. 'O' .AND. ID(J)(1:1) .EQ. 'H')) 
     X         DMAX = 1.55                  ! for O...H set dmax to 1.55 Angs
        IF ((ID(I)(1:1) .EQ. 'C' .AND. ID(J)(1:1) .EQ. 'S') .OR.    ! C-S =
     X      (ID(I)(1:1) .EQ. 'S' .AND. ID(J)(1:1) .EQ. 'C'))        ! 1.91 Angs
     X         DMAX = 1.91                
        IF ((ID(I)(1:1) .EQ. 'N' .AND. ID(J)(1:1) .EQ. 'S') .OR.    ! N-S =
     X      (ID(I)(1:1) .EQ. 'S' .AND. ID(J)(1:1) .EQ. 'N'))        ! 1.75 Angs
     X         DMAX = 1.75               
        IF ((ID(I)(1:1) .EQ. 'C' .AND. ID(J)(1:2) .EQ. 'BR') .OR.   ! C-BR =
     X      (ID(I)(1:2) .EQ. 'BR' .AND. ID(J)(1:1) .EQ. 'C'))       ! 1.99 Angs
     X         DMAX = 1.99               
        IF ((ID(I)(1:1) .EQ. 'C' .AND. ID(J)(1:2) .EQ. 'I') .OR. ! C-I =    ! 11/6/03
     X      (ID(I)(1:2) .EQ. 'I' .AND. ID(J)(1:1) .EQ. 'C'))     ! 2.3 Angs ! 11/6/03
     X         DMAX = 2.3                                                   ! 11/6/03
               D = 0.0
               DO K=1,3
                  D = D + (XYZO(K,I) - XYZO(K,J))**2
               ENDDO
               IF (SQRT(D) .LE. DMAX) THEN
                  NCON(I) = NCON(I) + 1
                  ICON(NCON(I),I) = J
               ENDIF
               IF (.NOT. (ID(I)(1:1) .EQ. 'C' .AND.
     X                    ID(J)(1:2) .EQ. 'BR') .OR.
     x                   (ID(I)(1:2) .EQ. 'BR' .AND.
     X                    ID(J)(1:1) .EQ. 'C')) GO TO 90
                  IF (SQRT(D) .GE. 1.95) GO TO 90
                    NCON(I) = NCON(I) + 1
                    ICON(NCON(I),I) = J
90       CONTINUE
100   CONTINUE
      RETURN
      END  
C
C------------------------------------------------------------------
C
      SUBROUTINE ADJUST_2                     ! 3/25/03 
C
C----Adjust C-H, N-H, N-F, O-H and N-O lengths...for use with
C     MOPAC2MOLPAK subroutine in TRANSFORM2MOPAC 
C
      CHARACTER LABEL*5, AH*1   
      COMMON /COM1/ NATOMS, CELL(6), LABEL(200), O_XYZ(3,200),
     x              N_NEAR(200), N_ATOMS(6,200)
C
      DMAX = 1.65              ! establish connectivity based on 1.65 Angs
      DMAXBR = 1.99            ! for bonds involving bromine
      DMAXS = 1.87             ! for bonds involving sulfur
C
C----Should the C-H lengths be adjusted?
      PRINT 10
10    FORMAT (' Should the C-H positions be stretched [Y]: ',$)
      READ (5,12) AH
12    FORMAT (A1)
      IACH = 0
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') IACH = 1 
C
C----Should the N-H lengths be altered?
      PRINT 101
101   FORMAT (' Should the N-H lengths be stretched? [Y]: ',$)
      READ (5,12) AH
      IANH = 0
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') IANH = 1
C
C----Should the NO2 distances be altered to 1.22 Angs
      PRINT 15
15    FORMAT (' Should the nitro group N-O distances be set to 1.22',
     x        ' Angs? [Y]: ',$)
      READ (5,12) AH                                                        
      IANO = 0
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') IANO = 1 
C
C----Should the N-F distances be altered to 1.40 Angs
      PRINT 45
45    FORMAT (' Should the N-F distances be set to 1.40',
     x        ' Angs? [Y]: ',$)
      READ (5,12) AH
      IANF = 0
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') IANF = 1
C                                                                       3/25/03
C----Should the alcohol O-H distances be altered to 0.97 Angs           3/25/03
      PRINT 46                                                        ! 3/25/03  
46    FORMAT (' Should the alcohol O-H distances be set to 0.97',     ! 3/25/03
     x        ' Angs? [Y]: ',$)                                       ! 3/25/03
      READ (5,12) AH                                                  ! 3/25/03
      IAN_ALC = 0                                                     ! 3/25/03
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') IAN_ALC = 1  ! 3/25/03
C
C----Determine connectivity only if some positions are to be 
C     altered....use max distance of 1.65 Angs    
      IF ((IACH + IANH + IANO + IANF + IAN_ALC) .EQ. 0) GO TO 75      ! 3/25/03
      DO 502 I=1,NATOMS
         N_NEAR(I) = 0
         DO 490 J=1,NATOMS
            IF (I .EQ. J) GO TO 490
            IF (LABEL(I)(1:1) .EQ. 'H' .AND.
     x          LABEL(J)(1:1) .EQ. 'H') GO TO 490
       IF ((LABEL(I)(1:1) .EQ. 'H' .AND. LABEL(J)(1:1) .EQ. 'O') .OR.  
     x     (LABEL(I)(1:1) .EQ. 'O' .AND. LABEL(J)(1:1) .EQ. 'H')) THEN
               DMAX = 1.55
            ELSE
               DMAX = 1.65
            ENDIF 
            D = DISTANCE (I, J)
            IF ((LABEL(I)(1:1) .EQ. 'C' .AND. LABEL(J)(1:2)
     X           .EQ. 'BR') .OR. (LABEL(J)(1:1) .EQ. 'C'
     x           .AND. LABEL(I)(1:2) .EQ. 'BR')) THEN
               IF (D .LE. DMAXBR) GO TO 488
            ENDIF
            IF ((LABEL(I)(1:1) .EQ. 'C' .AND. LABEL(J)(1:1)
     X           .EQ. 'S') .OR. (LABEL(J)(1:1) .EQ. 'C'
     x           .AND. LABEL(I)(1:1) .EQ. 'S')) THEN
               IF (D .LE. DMAXS) GO TO 488
            ENDIF
            IF (D .GT. DMAX) GO TO 490
488         N_NEAR(I) = N_NEAR(I) + 1
            N_ATOMS(N_NEAR(I),I) = J
490      CONTINUE
502   CONTINUE            
C
C----Report atom connectivities...
      WRITE (6,2220)
      WRITE (23,2220)
2220  FORMAT (' Atom connectivity....'/
     x        ' atom ....  linked to ...')                     
      DO I=1,NATOMS
         WRITE (6,2222) LABEL(I), (LABEL(N_ATOMS(K,I)),K=1,N_NEAR(I))
         WRITE (23,2222) LABEL(I), (LABEL(N_ATOMS(K,I)),K=1,N_NEAR(I))
2222     FORMAT (1X,A5,'....  ',6(A5,2X))
      ENDDO
C
C----Should the N-H lengths be altered?
75    IF (IANH .EQ. 0) GO TO 176
      DO 2001 I=1,NATOMS            ! just stretch the X-H position
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 2001
         J = N_ATOMS(1,I)           ! J = atom to which H is linked
         IF (LABEL(J)(1:1) .NE. 'N') GO TO 2001     ! C-H bond
         CALL STRETCH (J, I, 1.013)     ! N-H extension to 1.013
2001  CONTINUE
C
C----Should the N-F lengths be altered?
176   IF (IANF .EQ. 0) GO TO 276
      DO 2101 I=1,NATOMS            ! just stretch the N-F position
         IF (LABEL(I)(1:1) .NE. 'F') GO TO 2101
         J = N_ATOMS(1,I)           ! J = atom to which F is linked
         IF (LABEL(J)(1:1) .NE. 'N') GO TO 2101     ! N-F bond
         CALL STRETCH (J, I, 1.40)     ! N-F extension to 1.40
2101  CONTINUE
C
C----Should the alcohol O-H lengths be altered?                     3/25/03
276   IF (IAN_ALC .EQ. 0) GO TO 76                                ! 3/25/03
      DO 2201 I=1,NATOMS      ! just stretch the O-H position       3/25/03
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 2201                   ! 3/25/03
         J = N_ATOMS(1,I)     ! J = atom to which O is linked       3/25/03
         IF (LABEL(J)(1:1) .NE. 'O') GO TO 2201     ! O-H bond      3/25/03
         CALL STRETCH (J, I, 0.97)     ! O-H extension to 0.97      3/25/03
2201  CONTINUE                                                    ! 3/25/03
C
C----Should the NO2 distances be altered to 1.22 Angs
76    IF (IANO .EQ. 0) GO TO 77
      DO 100 I=1,NATOMS
         IF (LABEL(I)(1:1) .NE. 'O') GO TO 100
         IF (N_NEAR(I) .NE. 1) GO TO 100     ! not an N-O type O
         J = N_ATOMS(1,I)                    ! ID of 1 atom linked to the O
         IF (LABEL(J)(1:1) .NE. 'N') GO TO 100        ! gotta be a N
C----Check if the N is connected to 2 O's, each of which has only 1 bond
         NB2O = 0
            DO 110 L=1,N_NEAR(J)   ! J is the N connected to the O (I)
               K = N_ATOMS(L,J)    ! K = # of an atom connected to N (J)
               IF (LABEL(K)(1:1) .NE. 'O') GO TO 110
               IF (N_NEAR(K) .EQ. 1) NB2O = NB2O + 1    ! count # of singly
110         CONTINUE                                    ! bonded O's, must  
            IF (NB2O .NE. 2) GO TO 100                  ! be 2 of them   
C----Set a position along the N-O (I to J) vector
         CALL STRETCH (J, I, 1.22)
100   CONTINUE    
C
C----Check the C-H situation
77    IF (IACH .EQ. 0) RETURN
      DO 200 I=1,NATOMS            
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 200
         J = N_ATOMS(1,I)           ! J = atom to which H is linked
         IF (LABEL(J)(1:1) .NE. 'C') GO TO 200     ! C-H bond ?
C----For C-H extensions, check the C hydridization
3000     GO TO (3001, 3001, 3003, 3004), N_NEAR(J)  ! how many atoms on J ?
3001     PRINT 3002, I, J 
3002     FORMAT (' Cannot adjust distance between atoms ',A5,
     x           ' and ',A5)
         GO TO 200
3003     D = 1.084                      ! sp2
         GO TO 3010
3004     D = 1.098                      ! sp3
3010     CALL STRETCH (J, I, D)         ! adjust to length D
200   CONTINUE
      RETURN                  
      END
C
C--------------------------------------------------------------
C
C----stretch.f                    10/11/93 
C
C----Set an atom along a bond vector at a particular distance
C
      SUBROUTINE STRETCH (I_FIX, I_MOVE, B_LENGTH)
C
      CHARACTER LABEL*5, AH*1
      COMMON /COM1/ NATOMS, CELL(6), LABEL(200), O_XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      D = DISTANCE (I_MOVE, I_FIX)
      RATIO = B_LENGTH/D
      DO 10 I=1,3
        O_XYZ(I,I_MOVE) = (O_XYZ(I,I_MOVE) - O_XYZ(I,I_FIX))*RATIO
     x                   + O_XYZ(I,I_FIX)
10    CONTINUE
      PRINT 20, LABEL(I_FIX), LABEL(I_MOVE), D,
     x          (DISTANCE(I_MOVE, I_FIX))
      WRITE (23,20) LABEL(I_FIX), LABEL(I_MOVE), D,
     x          (DISTANCE(I_MOVE, I_FIX))
20      FORMAT (' Bond vector for ',A5,'-',A5,' changed from',
     x          F6.3,' to',F6.3,' Angs')
      RETURN
      END
C
C---------------------------------------------------------------
C
C----distance.f                  10/5/93 
C
C----calc distance between 2 atoms
C
      FUNCTION DISTANCE (I1, I2)
      CHARACTER LABEL*5, AH*1
      COMMON /COM1/ NATOMS, CELL(6), LABEL(200), O_XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      DISTANCE = 0.0
      DO I=1,3
         DISTANCE = DISTANCE + (O_XYZ(I,I1) - O_XYZ(I,I2))**2
      ENDDO
      DISTANCE = SQRT(DISTANCE)
      RETURN
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE MOPAC2MMOD (KIND)
C
C----Convert initial and final MOPAC formats to macromodel format.
C     Also Chem3D and hin format to macromodel format.
C
        PARAMETER (MAXA=200)
        PARAMETER (MAXB=6)

C *** SYMBOL DESCRIPTION
C     I,J,K  -- LOOP COUNTER
C     TOTATM -- TOTAL ATOM NUMBER
C     TOTBND -- TOTAL BOND NUMBER
C     TOTCHN -- TOTAL CHAIN NUMBER
C     NUMBER -- ATOM NUMBER OF CHEMX
C     MMTYPE -- ATOM TYPE OF MACROMODEL FILE
C     CGTYPE -- ATOM TYPE OF CHEMX FILE
C     ATTACH -- ATTACHED ATOM NUMBER
C     ORDERS -- BOND ORDER 
C     XYZPOS -- XYZ COORDINATES OF ATOM

        INTEGER   I,J,K
        INTEGER   TOTATM,TOTBND
        INTEGER   NUMBER(MAXA), NCON(MAXA)                                        
        INTEGER   MMTYPE(MAXA)
        CHARACTER CGTYPE(MAXA)*5, CC1*1, CC2*1, SERIAL*1, LINE*70
        REAL      XYZPOS(MAXA,3)
        INTEGER   ATTACH(MAXA,MAXB)
        INTEGER   ORDERS(MAXA,MAXB)

C     TITELS -- HEADLINE OF DATA FILE, USUALLY INCLUDING STRUCTURE
C               NAME OR ENERGY INFORMATCGTPON
C     LINE   -- HEADLINES IN CHEMX FILE
C     CGTYP1 -- FIRST CHARACTER OF CHEMX ATOM TYPE
C     CGTYP2 -- FIRST TWO CHARACTERS OF CHEMX ATOM TYPE
C     ATCGTP -- ATTACHED CHEMX ATOM TYPE
C     TEMP   -- LOCAL VARIABLE
C     DOUBLE -- DOUBLE BOND
C     UNDECD -- UNDECIDED BOND 
C     BOND0  -- ANY BOND
C     BOND1  -- FIRST BOND
C     BOND2  -- SECOND BOND
C     BOND3  -- THIRD BOND
C     BOND4  -- FORTH BOND
C     NUMATT -- NUMBER OF ATTACHED ATOMS OF ONE ATOM
C     NUMORD -- NUMBER OF ORDERS OF ONE ATOM
C     ATTAC0 -- ANY ATTACHED NUMBER
C     ATTAC3 -- THIRD ATTACHED ATOM NUMBER
C     ATTAC4 -- FORTH ATTACHED ATOM NUMBER
C     STARTP -- STAR POINT FOR DECIDING BOND AT RING
C     OK     -- LOGICAL FLAG 
C     DONE   -- LOGICAL FLAG
C     ATTAC5 -- ZERO VALUE FOR OUTPUT OF FIFTH ATTACHED ATOM 
C     ATTAC6 -- ZERO VALUE FOR OUTPUT OF SIXTH ATTACHED ATOM 
C     ORDER5 -- ZERO VALUE FOR OUTPUT OF FIFTH BOND ORDER
C     ORDER6 -- ZERO VALUE FOR OUTPUT OF SIXTH BOND ORDERR
    
        CHARACTER CGTYP1*1,CGTYP2*2,ATCGTP*2, FIRST*1
        CHARACTER INPFIL*30, OUTFIL*30
        INTEGER   TEMP
        INTEGER   DOUBLE,UNDECD
        INTEGER   BOND0,BOND1,BOND2,BOND3,BOND4
        INTEGER   NUMATT,NUMORD,ATTAC0,ATTAC3,ATTAC4
        INTEGER   STARTP
        LOGICAL   OK,DONE
        INTEGER   ATTAC5,ORDER1,ATTAC6,ORDER6
        
      TOTBND = 4
C----Request file names and open them
      IF (KIND .NE. 28) GO TO 300
          PRINT 290
290       FORMAT (' HIN coordinate file name: ',$)
          READ (5,2) INPFIL
          GO TO 128
C
300   IF (KIND .EQ. 22) THEN 
          PRINT 121
121       FORMAT (' CHEM3D coordinate file name: ',$)
          READ (5,2) INPFIL
2         FORMAT (A30)
          PRINT 57
57        FORMAT ('   Does the CHEM3D file have atom serial ',
     x            'numbers (Y/N) [Y]: ',$)
          READ (5,'(A1)') SERIAL      ! serial number?
            IF (SERIAL .EQ. 'y' .OR. SERIAL .EQ. 'Y'
     x           .OR. SERIAL .EQ. ' ') THEN
              IS = 1             ! serial numbers
            ELSE
              IS = 2             ! no serial numbers
            ENDIF    
      ELSE
          PRINT 1
1         FORMAT (' MOPAC coordinate file name: ',$)
          READ (5,2) INPFIL
          PRINT 18
18        FORMAT ('   Enter (1/2) for (initial/final) MOPAC',
     x           ' format: ',$)
          READ (5,*) IF_TYPE                      
      ENDIF
128   OPEN(UNIT=10,FILE=INPFIL,STATUS='OLD')
      PRINT 3 
3       FORMAT (' Macromodel coordinate file name: ',$)
      READ (5,2) OUTFIL
      OPEN (UNIT=11, FILE=OUTFIL, STATUS='UNKNOWN')
C
C----Read Mopac, Chem3d or hin atoms...if Chem3d, skip first line
      IF (KIND .EQ. 22) READ (10,*) JUNK           ! Chem3D input file
C
      I = 0
      TOTATM = 0
5     I = I + 1
C
      IF (KIND .EQ. 28) THEN                       ! HIN input file
20       READ (10,'(A70)',END=2098) LINE
         IF (LINE(1:4) .NE. 'atom') GO TO 20
40       L = INDEX(LINE,' - ')
         TOTATM = TOTATM + 1
         CGTYPE(TOTATM) = '     '
         READ (LINE(L+3:L+3),50) CGTYPE(TOTATM)(1:1)
50       FORMAT (A1)
         L = INDEX(LINE,' - 0 ')
         READ (LINE(L+5:70),*) (XYZPOS(TOTATM,J),J=1,3)
         GO TO 20
      ENDIF
      IF (KIND .NE. 22) GO TO 148
        CGTYPE(I) = '     '
        IF (IS .EQ. 2) THEN         ! IS = 2 means no serial numbers
          READ(10,2038,END=2099) CC1, CC2, (XYZPOS(I,J),J=1,3)
2038      FORMAT (1X,2A1,3F12.5)  ! Chem3D format w/o serial numbers
        ELSE
          READ(10,2041,END=2099) CC1, CC2, (XYZPOS(I,J),J=1,3)
2041      FORMAT (1X,2A1,5X,3F12.5)  ! Chem3D format with serial numbers  
        ENDIF
          IF (CC1 .EQ. ' ') THEN
             CGTYPE(I)(1:1) = CC2
          ELSE
             CGTYPE(I)(1:1) = CC1
             CGTYPE(I)(2:2) = CC2
          ENDIF
        GO TO 158
148   IF (IF_TYPE .NE. 1) THEN                    ! Mopac input file 
        READ(10,2030,END=2099) NUMBER(I), FIRST, CGTYPE(I), 
     x                         (XYZPOS(I,J),J=1,3)
2030    FORMAT (I6,8X,A1,A5,10X,3F10.4)     ! final Mopac format
        IF (FIRST .NE. ' ') THEN
           DO L=5,2,-1
              CGTYPE(I)(L:L) = CGTYPE(I)(L-1:L-1)
           ENDDO
           CGTYPE(I)(1:1) = FIRST
        ENDIF
      ELSE
        READ (10,2031,END=2099) CGTYPE(I), (XYZPOS(I,J),J=1,3)
2031    FORMAT (1X,A5,3(F10.5,3X))       ! initial mopac format
      ENDIF
158       CGTYP2 = CGTYPE(I)
          IF ( CGTYP2 .EQ. 'CL' 
     1    .OR. CGTYP2 .EQ. 'BR'
     2    .OR. CGTYP2 .EQ. 'SI'
     3    .OR. CGTYP2 .EQ. 'Cl'
     4    .OR. CGTYP2 .EQ. 'Br'
     5    .OR. CGTYP2 .EQ. 'Si'
     6    .OR. CGTYP2 .EQ. 'Lp') THEN
            CGTYPE(I) = CGTYP2
          ELSE
            CGTYP1 = CGTYP2
            CGTYP2 = CGTYP1
            CGTYPE(I) = CGTYP2
          END IF
      GO TO 5
2099  TOTATM = I - 1
C
C----Establish connectivity info...special version of CONNECT
2098  CALL CONNECT2 (TOTATM, ATTACH, NCON, CGTYPE, XYZPOS)
      do i=1,totatm
C      print 1234, i, cgtype(i), ncon(i), (xyzpos(i,j),j=1,3)
C1234  format (' #, id, ncon, xyz =',i3,2x,a5,i4,3f10.5)
      enddo
      DO 201 I=1,TOTATM
          IF ( CGTYP2 .EQ. 'N ' ) THEN
            DO J = 1,4
              IF ( ATTACH(I,J) .EQ. 0 ) THEN
                IF ( ATTACH(I,J+1) .NE. 0 ) THEN
                  TEMP = ATTACH(I,J)
                  ATTACH(I,J) = ATTACH(I,J+1)
                  ATTACH(I,J+1) = TEMP
                END IF
              END IF
            END DO
          END IF
201    CONTINUE
2000    FORMAT(1X,I4,1X,A2,1X,4(I4,1X,I2))

C *** SET UNDECIDED BOND ORDER TO -1 SO AS TO BE PROCESSED BY THE VARIOUS
C     STEPS BELLOW 
        DO I = 1,TOTATM
          DO J = 1,NCON(I)
            IF ( ATTACH(I,J) .NE. 0 )  THEN
              ORDERS(I,J) = -1
            END IF
          END DO
        END DO

C       TYPE *,' DECIDE BOND ORDER OF H,F,CL,BR,I,O,S AND C,N WHICH
C     # ATTACHES 4 ATOMS ... '
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)

          IF ( CGTYP2 .EQ. 'H '
     1    .OR. CGTYP2 .EQ. 'F ' 
     2    .OR. CGTYP2 .EQ. 'CL' 
     3    .OR. CGTYP2 .EQ. 'BR' 
     4    .OR. CGTYP2 .EQ. 'Br'
     5    .OR. CGTYP2 .EQ. 'I '
     6    .OR. CGTYP2 .EQ. 'Lp') THEN
            ORDERS(I,1) = 1
            CALL CONNECT3(I,1,1,TOTBND,ATTACH,ORDERS)
          END IF

          IF ( ( (CGTYP2 .EQ. 'C ')
     #    .OR. (CGTYP2 .EQ. 'SI')
     #    .OR. (CGTYP2 .EQ. 'N ') )
     #    .AND. (ATTACH(I,4) .NE. 0) ) THEN
            DO J = 1,4
              ORDERS(I,J) = 1
              CALL CONNECT3(I,J,1,TOTBND,ATTACH,ORDERS)
            END DO
          END IF

          IF ( (CGTYP2 .EQ. 'O ')
     #    .OR. (CGTYP2 .EQ. 'S ') ) THEN
            IF (ATTACH(I,2) .NE. 0) THEN
              DO J = 1,2
                ORDERS(I,J) = 1
                CALL CONNECT3(I,J,1,TOTBND,ATTACH,ORDERS)
              END DO
            ELSE
              ORDERS(I,1) = 2
              CALL CONNECT3(I,1,2,TOTBND,ATTACH,ORDERS)
            END IF
          END IF
        END DO

C       TYPE *,'MODIFY BOND ORDER OF P WHICH ATTACHES OXYGEN ... ' 
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          IF ( CGTYP2 .EQ. 'P ' ) THEN
            DOUBLE = 0
            DO J = 1,4
              IF ( ORDERS(I,J) .EQ. 2 ) THEN
                DOUBLE = DOUBLE + 1
              END IF
            END DO
            IF ( DOUBLE .EQ. 0 ) THEN
              DO J = 1,4
                ORDERS(I,J) = 1
                CALL CONNECT3(I,J,1,TOTBND,ATTACH,ORDERS)
              END DO
            END IF
            IF ( DOUBLE .GE. 2 ) THEN
              J = 4 
              DONE = .FALSE.
              DO WHILE ( (.NOT. DONE) .AND. (J .GT. 1) )
                IF ( ORDERS(I,J) .EQ. 2 ) THEN
                  ORDERS(I,J) = 1
                  CALL CONNECT3(I,J,1,TOTBND,ATTACH,ORDERS)
                  DONE = .TRUE.  
                END IF
                J = J-1
              END DO
            END IF
          END IF
        END DO

C      TYPE *,'DECIDE BOND ORDER OF C,N WHICH ATTACHES TWO ATOMS ... '
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          IF ( (CGTYP2 .EQ. 'C ' .OR. CGTYP2 .EQ. 'N ')          ! 9/6/91
     #   .AND. ATTACH(I,3) .EQ. 0 ) THEN
            BOND1 = ORDERS(I,1)
            BOND2 = ORDERS(I,2)
            IF ( (BOND1 .EQ. -1) .OR. (BOND2 .EQ. -1) ) THEN
              IF ( (BOND1 .EQ. 1) .AND. (BOND2 .EQ. -1) ) THEN 
C *** BOND ORDER 3 : -C%C-  -C%N ***
                ORDERS(I,2) = 3
                CALL CONNECT3(I,2,3,TOTBND,ATTACH,ORDERS)
              END IF
              IF ( (BOND2 .EQ. 1) .AND. (BOND1 .EQ. -1) ) THEN
                ORDERS(I,1) = 3
                CALL CONNECT3(I,1,3,TOTBND,ATTACH,ORDERS)
              END IF
C *** ATTACH TWO ATOMS: -HC=C=CH-   -HC=N=CH- ***
              IF ( (BOND1 .EQ. -1) .AND. (BOND2 .EQ. -1) ) THEN
                DO J = 1,2
                  ORDERS(I,J) = 2
                  CALL CONNECT3(I,J,2,TOTBND,ATTACH,ORDERS)
                END DO
              END IF
            ENDIF
          END IF
        END DO

C       TYPE *,' DECIDE BOND ORDER N WHICH ATTACHES 3 ATOMS ... '
C *** THE BOND ORDER OF NITROGEN WHICH ATTACHES THREE ATOMS IS ASSIGNED
C     TO BE SINGLE BOND.  THIS IS SUTIBALE FOR MOST OF THE CASES. HOWEVER,
C     THERE MUST BE SOME EXCEPTIONS IN WHICH THE BOND ORDER COULD BE A
C     DOUBLE BOND WITH A POSITIVE CHARGE.  THIS PROBLEM CAN NOT BE SOLVED
C     UPTIL NOW AND LEFT FOR FURTHER IMPROVEMENT. 
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          IF ( (CGTYP2 .EQ. 'N ') 
     #    .AND. (ATTACH(I,3) .NE. 0) 
     #    .AND. (ATTACH(I,4) .EQ. 0) ) THEN
            DO  J = 1,3
              IF ( ORDERS(I,J) .EQ. -1 ) THEN
                ORDERS(I,J) = 1
                ATTAC0 = ATTACH(I,J)
                DO K = 1,TOTBND
                  IF ( (ATTACH(ATTAC0,K) .EQ. I)
     #            .AND. (ORDERS(ATTAC0,K) .EQ. -1 ) ) THEN
                     ORDERS(ATTAC0,K) = 1
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
         
C       TYPE *, 'DECIDE BOND ONDER OF C,N WHICH ATTACHES 3 ATOMS AND
C     # NOT IN THE RINGS ... '
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          DOUBLE = 0
          UNDECD = 0
          IF ( (CGTYP2 .EQ. 'C ' .OR. CGTYP2 .EQ. 'N ')
     #    .AND. ATTACH(I,3) .NE. 0
     #    .AND. ATTACH(I,4) .EQ. 0 ) THEN
            NUMATT = 3
          END IF
          IF ( (CGTYP2 .EQ. 'N ') .AND. ATTACH(I,3) .EQ. 0 ) THEN
            NUMATT = 2
          END IF
          DO J = 1,NUMATT
            IF ( ORDERS(I,J) .EQ. 2 ) DOUBLE = DOUBLE + 1
            IF ( ORDERS(I,J) .EQ. -1 ) UNDECD = UNDECD + 1
          END DO
                    
          IF ( UNDECD .EQ. 1 ) THEN
            IF ( DOUBLE .EQ. 1 ) THEN
              BOND0 = 1
            END IF
            IF ( DOUBLE .EQ. 0 ) THEN
              BOND0 = 2
            END IF
            DO J = 1,NUMATT
              IF ( ORDERS(I,J) .EQ. -1 ) THEN
                ORDERS(I,J) = BOND0
                CALL CONNECT3(I,J,BOND0,TOTBND,ATTACH,ORDERS)
              END IF
            END DO
          END IF
        END DO

C        TYPE *,'DECIDE BOND ORDER OF C,N WHICH ATTACHES 3 ATOMS AND
C     # IN RINGS ... '
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          IF ( (CGTYP2 .EQ. 'C ') 
     #    .AND. (ATTACH(I,3) .NE. 0)
     #    .AND. (ATTACH(I,4) .EQ. 0) ) THEN
            UNDECD = 0
            DO J = 1,3
              IF (ORDERS(I,J) .EQ. -1 ) UNDECD = UNDECD + 1
            END DO
            IF ( UNDECD .EQ. 2 ) THEN
              STARTP = I
C       TYPE *,'ORIGINAL START CURRPT :  ',STARTP
              CALL RING(STARTP,TOTATM,TOTBND,CGTYPE,ATTACH,ORDERS)
            ELSE 
              IF ( UNDECD .EQ. 1 ) THEN
                STARTP = I
C       TYPE *,'ORIGINAL START CURRPT :  ',STARTP
                CALL RING(STARTP,TOTATM,TOTBND,CGTYPE,ATTACH,ORDERS)
              END IF
            END IF
          END IF   
        END DO

C        TYPE *,'MODIFY BOND ORDER OF C,N WHICH ATTACHES TWO OXYGENS ...'
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          NUMORD = 0
          IF ( (CGTYP2 .EQ. 'C ') .OR. (CGTYP2 .EQ. 'N ') ) THEN
            ATTAC3 = ATTACH(I,3)
            ATTAC4 = ATTACH(I,4)
            IF ( ATTAC4 .NE. 0 )THEN 
              NUMATT = 4
            END IF
            IF ( (ATTAC3 .NE. 0) .AND. (ATTAC4 .EQ. 0) ) THEN
              NUMATT = 3
            END IF
            DO J = 1,NUMATT
              NUMORD = NUMORD + ORDERS(I,J) 
            END DO
          END IF      
          IF ( NUMORD .GT. 4 ) THEN
            DONE = .FALSE.
            J = 3
            DO WHILE ( (.NOT. DONE) .AND. (J. GT. 1) )
              ATTAC0 = ATTACH(I,J)
              IF ( CGTYPE(ATTAC0) .EQ. 'O ' ) THEN
                ORDERS(I,J) = 1
                CALL CONNECT3(I,J,1,TOTBND,ATTACH,ORDERS)
                DONE = .TRUE.
              END IF
              J = J-1 
            END DO
          END IF
        END DO     
        
C        TYPE *,'DECIDE ATOM TYPE ... '
        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          IF ( CGTYP2 .EQ. 'S ' ) MMTYPE(I) = 52
          IF ( CGTYP2 .EQ. 'S ' ) THEN
            BOND2 = ORDERS(I,2)
            IF ( BOND2 .EQ. 1 ) THEN
              MMTYPE(I) = 49
            END IF
          END IF
          IF ( CGTYP2 .EQ. 'P ' ) MMTYPE(I) = 53
          IF ( CGTYP2 .EQ. 'F ' ) MMTYPE(I) = 56
          IF ( CGTYP2 .EQ. 'CL' ) MMTYPE(I) = 57
          IF ( CGTYP2 .EQ. 'BR' ) MMTYPE(I) = 58
          IF ( CGTYP2 .EQ. 'Br' ) MMTYPE(I) = 58
          IF ( CGTYP2 .EQ. 'I ' ) MMTYPE(I) = 59
          IF ( CGTYP2 .EQ. 'SI' ) MMTYPE(I) = 60
          IF ( CGTYP2 .EQ. 'Lp' ) MMTYPE(I) = 63
          IF ( CGTYP2 .EQ. 'H ' ) THEN
            ATTAC0 = ATTACH(I,1)
            ATCGTP = CGTYPE(ATTAC0)

C *** THE CODE 48 IS ANY TYPE OF HYDROGEN IN MACROMODEL.  HOWEVER, IF A 
C     HYDROGEN WITH THIS CODE OCCURES IN A MOLECULE, THE PROCEDURE OF 
C     ENERGY MINIMIZATION CAN NOT BE DONE BECAUSE OF THE FORCE FIELD 
C     PROBLEM.  IT HAS BEEN SUCCESSFULLY TESTED TO USE CODE 41 INSTEAD 
C     OF 48 AS THE ATOM TYPE OF HYDROGEN.
C           MMTYPE(I) = 48
C            IF ( (ATCGTP .EQ. 'C ') 
C     1      .OR. (ATCGTP .EQ. 'S ') 
C     2      .OR. (ATCGTP .EQ. 'SI') ) MMTYPE(I) = 41
            MMTYPE(I) = 41
            IF ( ATCGTP .EQ. 'O ' ) MMTYPE(I) = 42
            IF ( ATCGTP .EQ. 'N ' ) THEN 
              BOND4 = ORDERS(ATTAC0,4)
              IF ( BOND4 .EQ. 0 ) THEN
                MMTYPE(I) = 43
              END IF
              IF ( BOND4 .EQ. 1 ) THEN
                MMTYPE(I) = 44
              END IF
            END IF
          END IF

          IF ( CGTYP2 .EQ. 'O ' ) THEN
            MMTYPE(I) = 23
            BOND1 = ORDERS(I,1)
            BOND2 = ORDERS(I,2)
            IF ( BOND2 .NE. 0 ) THEN
              MMTYPE(I) = 16
            ELSE 
              IF ( BOND1 .EQ. 2 ) THEN
                MMTYPE(I)= 15
              END IF
              IF ( BOND1 .EQ. 1 ) THEN
                MMTYPE(I)= 18
              END IF
            END IF
          END IF

          IF ( CGTYP2 .EQ. 'C ' ) THEN
            MMTYPE(I) = 14
            BOND1 = ORDERS(I,1)
            BOND2 = ORDERS(I,2)
            BOND3 = ORDERS(I,3)
            BOND4 = ORDERS(I,4)
            IF ( BOND3 .EQ. 0 ) THEN
              IF ( BOND1 .EQ. BOND2 ) THEN 
                MMTYPE(I) = 2
              ELSE
                MMTYPE(I) = 1
              END IF
            ELSE 
              IF ( BOND4 .EQ. 0 ) THEN
                MMTYPE(I) = 2
              ELSE 
                MMTYPE(I) = 3
              END IF
            END IF
          END IF
        END DO

        DO I = 1,TOTATM
          CGTYP2 = CGTYPE(I)
          IF ( CGTYP2 .EQ. 'N ' ) THEN
            MMTYPE(I) = 40
            BOND2 = ORDERS(I,2)
            BOND3 = ORDERS(I,3)
            BOND4 = ORDERS(I,4)
            IF ( BOND2 .EQ. 0 ) THEN
              MMTYPE(I) = 24
            ELSE
              IF ( BOND3 .EQ. 0 ) THEN
                MMTYPE(I) = 25
              ELSE
                IF ( BOND4 .EQ. 0 ) THEN
                  DOUBLE = 0
                  DO J = 1,3
                    IF ( ORDERS(I,J) .EQ. 2 ) DOUBLE = DOUBLE + 1
                  END DO
                  IF ( DOUBLE .GE. 2 ) THEN
                    PRINT *,' WRONG BOND WAS FOUND IN NITROGEN !'
                    STOP
                  END IF
                  IF ( DOUBLE .EQ. 0 ) THEN
                    MMTYPE(I) = 26
                    DO J = 1,3
                      ATTAC0 = ATTACH(I,J)
                      IF ( MMTYPE(ATTAC0) .EQ. 2 ) THEN
                        MMTYPE(I) = 25
                      END IF
                    END DO
                  END IF
                  IF ( DOUBLE .EQ. 1 ) THEN
                    MMTYPE(I) = 31
                  END IF
                ELSE
                  MMTYPE(I) = 32
                END IF
              END IF
            END IF
          END IF
        END DO 

        OK = .TRUE.
        DO I = 1,TOTATM
           OK = .TRUE.
           DO J = 1,TOTBND
             IF ( ORDERS(I,J) .EQ. -1 ) OK = .FALSE.
           END DO
           IF ( .NOT. OK ) THEN
             PRINT 1990, I, CGTYPE(I) 
1990         FORMAT (' Undetermined bond orders (-1) associated with',
     x               ' atom #',I3,' (',A5,')...')
             DO 1993 L=1,TOTBND
                IF (ATTACH(I,L) .NE. 0)   
     X          PRINT 1996, ATTACH(I,L),  ORDERS(I,L)
1996            FORMAT ('   attached to atom #',I3,', b. o. =',I3)
1993         CONTINUE
               PRINT *, '   check the output file (*) and make',
     x                  ' appropriate corrections'
           END IF
        END DO 
2040    FORMAT(I4,6(1X,I5,1X,I1),1X,3(F11.6,1X))
        WRITE(11,2020) TOTATM
2020    FORMAT (I6)
        ATTAC5 = 0
        ORDER1 = 0
        ATTAC6 = 0
        ORDER6 = 0
        DO I = 1,TOTATM
          WRITE(11,2040)
     1    MMTYPE(I),
     2   (ATTACH(I,J),ORDERS(I,J),J=1,TOTBND),
     3   ATTAC5,ORDER1,ATTAC6,ORDER6,
     3   (XYZPOS(I,J),J=1,3)
        END DO
        CLOSE (11)
        RETURN
        END
C
C--------------------------------------------------------------------
C
C *** SUBROUTINE CONNECT3 
C     IF A BOND ORDER OF ANY ATOM HAS BEEN DECIDED, THE BOND ORDER OF
C     THE ATOM WHICH ATTACHES THIS PREDECIDED ATOM CAN ALSO BE DECIDED.

        SUBROUTINE CONNECT3(CONNAN,CONNBN,CONNBO,TOTBND,ATTACH,ORDERS)

        PARAMETER (MAXA=200)
        PARAMETER (MAXB=6)

C *** SYMBOL DESCRIPTION
C     I,J,K  -- LOOP COUNTER
C     CONNAN -- CONNECTED ATOM NUMBER
C     CONNBN -- CONNECTED BOND NUMBER
C     CONNBO -- CONNECTED BOND ORDER
C     TOTBND -- TOTAL BOND NUMBER
C     ATTACH -- ATTACHED ATOM NUMBER
C     ORDERS -- BOND ORDER
C     ATTAC0 -- ANY ATTACHED ATOM NUMBER
  
        INTEGER   I,J,K
        INTEGER   CONNAN,CONNBN,CONNBO,TOTBND
        INTEGER   ATTACH(MAXA,MAXB)
        INTEGER   ORDERS(MAXA,MAXB)
        INTEGER   ATTAC0

        DO J = 1,TOTBND
          ATTAC0 = ATTACH(CONNAN,CONNBN)
          IF (ATTACH(ATTAC0,J) .EQ. CONNAN) THEN
            IF ( ORDERS(ATTAC0,J) .NE. ORDERS(CONNAN,CONNBN) ) THEN
              ORDERS(ATTAC0,J) = CONNBO
            END IF
          END IF
        END DO

        RETURN
        END
C----------------------------------------------------------------
C *** SUBROUTINE RING -- SUBPROCEDURE OF CG_MM TO SOLVE THE PROBLEM OF
C     DECIDING THE BOND ORDER OF C,N WHICH ARE AT THE RINGS OF BENZENE
C     OR OTHER MOLECULES

        SUBROUTINE RING(CURRPT,TOTATM,TOTBND,CGTYPE,ATTACH,ORDERS)

        PARAMETER (MAXA=200)
        PARAMETER (MAXB=6)

C *** SYMBOL DESCRIPTION
C     I,J,K  -- LOOP COUNTER
C     TOTATM -- TOTAL ATOM NUMBER
C     TOTBND -- TOTAL BOND NUMBER
C     CGTYPE -- ATOM TYPE OF CHEMX FILE
C     ATTACH -- ATTACHED ATOM NUMBER
C     ORDERS -- BOND ORDER 
C     ATTAC0 -- ANY ATTACHED ATOM NUMBER
C     CGTYP2 -- TEMPERARY CHEMX ATOM TYPE 
C     CURRPT -- CURRENT POINT
C     NEXTPT -- NEXT POINT
C     RESEPT -- RESERVED POINT
C     RESPTN -- RESERVED POINT NUMBER
C     RESPNT -- RESERVED POINT ARRAY        

C     DOUBLE -- DOUBLE BOND
C     UNDECD -- UNDECIDED BOND 
C     BOND0  -- ANY BOND
C     BOND1  -- FIRST BOND
C     BOND2  -- SECOND BOND
C     BOND3  -- THIRD BOND
C     BOND4  -- FORTH BOND
C     NUMATT -- NUMBER OF ATTACHED ATOMS OF ONE ATOM

        INTEGER   I,J,K
        INTEGER   TOTATM,TOTBND
        CHARACTER CGTYPE(MAXA)*5
        INTEGER   ATTACH(MAXA,MAXB)
        INTEGER   ORDERS(MAXA,MAXB)
        CHARACTER CGTYP2*2
        INTEGER   CURRPT,NEXTPT,RESEPT,RESPTN
        INTEGER   RESPNT(20)
        INTEGER   DOUBLE,UNDECD
        INTEGER   BOND0,BOND1,BOND2,BOND3
        INTEGER   NUMATT

2110     FORMAT(1X,I4,1X,A2,1X,4(I4,1X,I2))

C        TYPE *,'START CURRPT =',CURRPT
        RESPTN = 0
2115    CGTYP2 = CGTYPE(CURRPT)
C        TYPE *,'CURRENT CURRPT IS : ',CURRPT
        DOUBLE = 0
        UNDECD = 0
        IF ( (CGTYP2 .EQ. 'C ' .OR. CGTYP2 .EQ. 'N ')
     #  .AND. ATTACH(CURRPT,3) .NE. 0
     #  .AND. ATTACH(CURRPT,4) .EQ. 0 ) THEN
          NUMATT = 3
        END IF
        IF ( (CGTYP2 .EQ. 'N ' ) 
     #  .AND. ATTACH(CURRPT,3) .EQ. 0 ) THEN
          NUMATT = 2
        END IF
        DO J = 1,NUMATT
          IF ( ORDERS(CURRPT,J) .EQ. 2 ) DOUBLE = DOUBLE + 1
          IF ( ORDERS(CURRPT,J) .EQ. -1 ) UNDECD = UNDECD + 1
        END DO

        IF ( UNDECD .EQ. 0 ) THEN
C        TYPE *,'NO BOND UNDECIDED'
          RETURN
        END IF
        IF ( UNDECD .EQ. 1 ) THEN
C        TYPE *,'ONE BOND UNDECIDED'
          IF ( DOUBLE .EQ. 1 ) THEN
            BOND0 = 1
          END IF
          IF ( DOUBLE .EQ. 0 ) THEN
            BOND0 = 2
          END IF
          DO J = 1,NUMATT
            IF ( ORDERS(CURRPT,J) .EQ. -1 ) THEN
              ORDERS(CURRPT,J) = BOND0
              NEXTPT = ATTACH(CURRPT,J)
            END IF
          END DO
          DO J = 1,3
            IF ( ATTACH(NEXTPT,J) .EQ. CURRPT 
     #      .AND. ORDERS(NEXTPT,J) .EQ. -1 ) THEN
              ORDERS(NEXTPT,J) = BOND0
            END IF
          END DO
        END IF
        IF ( UNDECD .EQ. 2 ) THEN
C        TYPE *,'TWO TOTBND UNDECIDED'
          IF ( DOUBLE .EQ. 1 ) THEN
            BOND1 = 1
            BOND2 = 1
            CALL BOND(BOND1,BOND2,CURRPT,NEXTPT,RESEPT,RESPTN,
     #                ATTACH,ORDERS)
C        TYPE *,'RESERVED POINT NUMBER =',RESPTN
            RESPNT(RESPTN) = RESEPT
          END IF
          IF ( DOUBLE .EQ. 0 ) THEN
            IF ( RESPTN .EQ. 0 ) THEN 
              BOND1 = 2
              BOND2 = 1
              CALL BOND(BOND1,BOND2,CURRPT,NEXTPT,RESEPT,RESPTN,
     #                  ATTACH,ORDERS)
              RESPNT(RESPTN) = RESEPT
            ELSE
              CURRPT = RESPNT(RESPTN)
              RESPTN = RESPTN - 1
              GO TO 2115
            END IF
          END IF
        END IF
C        WRITE(*,2110)CURRPT,CGTYPE(CURRPT),(ATTACH(CURRPT,J),
C     #                     ORDERS(CURRPT,J),J=1,TOTBND)
        CURRPT = NEXTPT
C        TYPE *,'NEXT POINT IS :',NEXTPT
C        TYPE *,'  '
        GO TO 2115
        END
C--------------------------------------------------------------
C *** SUBROUTINE BOND -- SUBPROCEDURE OF SUBROUTINE RING
C     DECIDE BOND ORDER OF ATOM C,N WHICH ATTACHES THREE ATOMS AND HAS
C     TWO UNDECIDED BONDS BY ASSIGNING THE FIRST BOND IS A DOUBLE BOND
C     AND THE SECOND BOND IS A SINGLE BOND.

        SUBROUTINE BOND(BOND1,BOND2,CURRPT,NEXTPT,RESEPT,RESPTN,
     #                  ATTACH,ORDERS)

        PARAMETER (MAXA=200)
        PARAMETER (MAXB=6)

C *** SYMBOL DESCRIPTION
C     I,J,K  -- LOOP COUNTER
C     ATTACH -- ATTACHED ATOM NUMBER
C     ORDERS -- BOND ORDER 
C     BOND1  -- FIRST BOND
C     BOND2  -- SECOND BOND
C     CURRPT -- CURRENT POINT
C     NEXTPT -- NEXT POINT
C     RESEPT -- RESERVED POINT
C     RESPTN -- RESERVED POINT NUMBER
C     RESPNT -- RESERVED POINT ARRAY        
C     DONE1  -- LOGICAL FLAG
C     DONE2  -- LOGICAL FLAG

        INTEGER   I,J,K
        INTEGER   BOND1,BOND2
        INTEGER   CURRPT,NEXTPT,RESEPT,RESPTN
        INTEGER   ATTACH(MAXA,MAXB)
        INTEGER   ORDERS(MAXA,MAXB)
        LOGICAL   DONE1,DONE2

        DONE1 = .FALSE.
        DONE2 = .FALSE.
        DO J = 1,3
          IF ( ORDERS(CURRPT,J) .EQ. -1 
     #    .AND. .NOT. DONE1 ) THEN
            ORDERS(CURRPT,J) = BOND1
            NEXTPT = ATTACH(CURRPT,J)
            DONE1 = .TRUE.
          END IF
          IF ( ORDERS(CURRPT,J) .EQ. -1 
     #    .AND. .NOT. DONE2 ) THEN
            ORDERS(CURRPT,J) = BOND2
            RESPTN = RESPTN + 1
            RESEPT = ATTACH(CURRPT,J)
            DONE2 = .TRUE.
          END IF
        END DO
        DO J = 1,3
          IF ( ATTACH(NEXTPT,J) .EQ. CURRPT 
     #    .AND. ORDERS(NEXTPT,J) .EQ. -1 ) THEN
            ORDERS(NEXTPT,J)=BOND1
          END IF
        END DO
        DO J = 1,3
          IF ( ORDERS(RESEPT,J) .EQ. -1
     #    .AND. ATTACH(RESEPT,J) .EQ. CURRPT ) THEN
            ORDERS(RESEPT,J) = BOND2
          END IF
        END DO
        RETURN
        END
C
C--------------------------------------------------------------------
C
C----CONNECT2...build ATTACH array
      SUBROUTINE CONNECT2 (TOTATM, ATTACH, NCON, ID, XYZPOS)
C
      PARAMETER (MAXA=200, MAXB=6)
      CHARACTER ID(MAXA)*5
      INTEGER ATTACH(MAXA,MAXB), NCON(MAXA), TOTATM
      DIMENSION XYZPOS(MAXA,3)
C
      DO 100 I=1,TOTATM
         NCON(I) = 0
         DO 90 J=1,TOTATM
            IF (I .EQ. J) GO TO 90
            IF (ID(I)(1:1) .EQ. 'H' .AND. ID(J)(1:1) .EQ. 'H')
     X                         GO TO 90
               D = 0.0
               DO K=1,3
                  D = D + (XYZPOS(I,K) - XYZPOS(J,K))**2
               ENDDO
               IF (SQRT(D) .LE. 1.63) THEN
                  NCON(I) = NCON(I) + 1
                  ATTACH(I,NCON(I)) = J
               ENDIF
90       CONTINUE
100   CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C 
C     Convert from CHEMX to MOLPAK format...reorient molecule so that
C       2 atoms ly on a C2 axis
C
      SUBROUTINE CHEMX_C2_MOLPAK
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      DIMENSION IK(200), T(3,3), LL(3), MM(3), V1(3), V2(3),
     X          V3(3), V4(3), XM1(3,3), XYZ(3,200), C2XYZ(3), 
     X          JAX(2,3), NPOT2(200), C_G03(200), C_MNDO(200)
      CHARACTER NAME*5, S*1
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ_CART(3,200),
     x              N_NEAR(200), N_ATOMS(6,200)
      CHARACTER OUTFILE2(32)*1, INFILE2(32)*1 
      CHARACTER*32 INFILE, OUTFILE, MOLPAKXYZ
      DATA MOLPAKXYZ /'molpak.xyz                      '/
      EQUIVALENCE (OUTFILE,OUTFILE2), (INFILE, INFILE2)
      DATA OUTFILE2 /32*' '/
      DATA JAX /2,3,1,3,1,2/
C
      PRINT 1
 1    FORMAT(' CHEMX input file name: ',$) 
      READ (5,2) INFILE
    2 FORMAT(A32)
C----OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
C----Output file
      PRINT 11
11    FORMAT (' MOLPAK output coordinate file name [molpak.xyz]: ',$)
      READ (5,2) OUTFILE
      IF (OUTFILE2(1) .EQ. ' ') OUTFILE = MOLPAKXYZ     
72    OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
C----Read cell params from lines 1 and 2 or CHEMX file
      READ (10,770) CELL
770   FORMAT (38X,3F8.3/21X,3F8.3)
      READ (10,*)         ! skip 2 lines
      READ (10,*)                       
      NA=0
777   NA=NA+1
         READ (10,30,END=999) IK(NA), NAME(NA)(1:1), (XYZ(I,NA),I=1,3)
30       FORMAT (I4,1X,A1,4X,3F10.5)
      GO TO 777
  999 NA=NA-1
C----Conversion from fractional to cartesian coordinates
          CALL MATRXT (CELL, T)
          CALL MTIMES1 (NA, CELL, XYZ, T, XYZ_CART)
C-----Should the C-H, N-H and/or N-O vectors be adjusted?
997   CALL ADJUST_2
C----Get information on two-fold axis and some other direction
      PRINT 600
600   FORMAT (' The molecule will be reoriented with 2 atoms along',
     x        ' a C2 axis...'/
     x        '   give the sequence numbers of 2 atoms on the axis: ',$)
      READ (5,*) N1, N2 
      PRINT 608
608   FORMAT ('   give the sequence numbers of 2 atoms not on the ',
     x        'C2 axis to'/
     x        '      specify a vector approximately',
     x        ' normal to the first: ',$)
      READ (5,*) N3, N4
      PRINT 609
609   FORMAT ('   type 1, 2 or 3 for C2 axis should be along',
     x        ' x, y, z: ',$)
      READ (5,*) IAX1
      IAX2 = JAX(1,IAX1)
      IAX3 = JAX(2,IAX1)
      PRINT 611
611   FORMAT (' Should the C2-related atoms be averaged and ',
     x        'only a unique set output [Y]: ',$)
      READ (5,'(A1)') S
      IS = 0
      IF (S .EQ. 'Y' .OR. S .EQ. 'y' .OR. S .EQ. ' ') IS = 1
C----Form V1 (N1 - N2), normalize and form V2 (N3 - N4)
      CALL VECTOR (V1, XYZ_CART(1,N1), XYZ_CART(1,N2))
      CALL N_VECTOR (V1)
      CALL VECTOR (V2, XYZ_CART(1,N3), XYZ_CART(1,N4))
C----Vector cross-product (V3 = V1 x V2); normalize V3
      CALL V_CROSS (V3, V1, V2)
      CALL N_VECTOR (V3)
C----Vector cross-product (V4 = V3 x V1); V4 approximately parallel to V2
      CALL V_CROSS (V4, V1, V3)
      CALL N_VECTOR (V4)      ! just to be sure
C----Form 3 x 3 matrix...XM1 x new = old
      DO 612 I=1,3
         XM1(I,IAX1) = V1(I)
         XM1(I,IAX2) = V3(I)
         XM1(I,IAX3) = V4(I)
612   CONTINUE
C----Invert the sucker
      CALL MINV (XM1, 3, DETERM, LL, MM)           ! fixed, 2/18/01
C----Convert XYZ_CART to new system 
      CALL MTIMES2 (NA, XYZ_CART, XM1, XYZ)
      X=0.0                     ! move center of molecule
      Y=0.0                     ! to 0, 0, 0
      Z=0.0
      I = 0
      DO 614 J=1,NA
         IF (NAME(J)(1:1) .EQ. 'X') GO TO 614
         I = I + 1
         X=X+XYZ(1,J)
         Y=Y+XYZ(2,J)
         Z=Z+XYZ(3,J)
614   CONTINUE
      X=X/I
      Y=Y/I
      Z=Z/I
      DO J=1,NA
        XYZ(1,J)=XYZ(1,J)-X
        XYZ(2,J)=XYZ(2,J)-Y
        XYZ(3,J)=XYZ(3,J)-Z
      ENDDO                               
C----Do we want to average C2-related atoms and only produce a 
C     unique set?
      IF (IS .EQ. 0) GO TO 821
      DO 620 I=1,NA                
         IF (NAME(I)(1:1) .EQ. 'X') GO TO 620
         NPOT2(I) = NPOT(I)
         C_G03(I) = G03_CHARGE
         C_MNDO(I) = XMNDO_CHARGE                 ! any atoms for
         IF (ABS(XYZ(IAX2,I)) .LE. 0.05 .AND.    ! which both the new
     X       ABS(XYZ(IAX3,I)) .LE. 0.05) THEN    ! y and z are .le.
            XYZ(IAX2,I) = 0.0                     ! 0.05, will have
            XYZ(IAX3,I) = 0.0                     ! y = z = 0
            WRITE (11,811) NAME(I)(1:1), I, (XYZ(J,I),J=1,3), NPOT2(I),
     x                     C_G03(I), C_MNDO(I) 
            XYZ(IAX1,I) = 10000.                    ! set a flag 
         ENDIF
620   CONTINUE
C----Average 2-fold related atoms (within 0.2 Angs) and keep those
C     that are in the +/-y, + z hemisphere
      DO 650 I=1,NA-1
         IF (NAME(I)(1:1) .EQ. 'X') GO TO 650
         IF (XYZ(IAX1,I) .GE. 9999.) GO TO 650
         C2XYZ(IAX1) = XYZ(IAX1,I)
         C2XYZ(IAX2) = -XYZ(IAX2,I)
         C2XYZ(IAX3) = -XYZ(IAX3,I)
         DO 625 K=I+1,NA
            IF (XYZ(IAX1,K) .GE. 9999.0) GO TO 625
            D = 0.0
            DO 626 L=1,3
               D = D + (XYZ(L,K) - C2XYZ(L))**2
626         CONTINUE 
            IF (SQRT(D) .GT. 0.2) GO TO 625
            DO 627 J=1,3
               XYZ(J,K) = (C2XYZ(J) + XYZ(J,K))/2.0  ! av xyz's
627         CONTINUE 
            IF (XYZ(IAX3,K) .GE. 0) GO TO 630
               XYZ(IAX2,K) = -XYZ(IAX2,K)               ! place into +/- y
               XYZ(IAX3,K) = -XYZ(IAX3,K)               ! + z
630            WRITE (11,811) NAME(K)(1:1), I, (XYZ(J,K),J=1,3), 
     x                        NPOT2(K), (0.5*(C_G03(I) + C_G03(K))), 
     x                        (0.5*(C_MNDO(I) + C_MNDO(K))) 
               XYZ(IAX1,I) = 10000.
               XYZ(IAX1,K) = 10000.
625      CONTINUE
650   CONTINUE
      RETURN
C----Output MOLPAK format
821   DO 810 I=1,NA
         IF (NAME(I)(1:1) .EQ. 'X') GO TO 810
         L = NPOT(I)
         WRITE (11,811) NAME(I)(1:1), I, (XYZ(J,I),J=1,3), L,
     x                  G03_CHARGE, XMNDO_CHARGE 
811      FORMAT('ATOM ',A1,I4,1X,3F10.6,I5,2F10.6)    ! MOLPAK format
810   CONTINUE
      RETURN
      END
C
C--------------------------------------------------------------------
C----Form a vector
      SUBROUTINE VECTOR (V, A, B)
      DIMENSION V(3), A(3), B(3)
      DO I=1,3
         V(I) = A(I) - B(I)
      ENDDO
      RETURN
      END
C----Normalize a vector
      SUBROUTINE N_VECTOR (V)
      DIMENSION V(3)
      S = 0.0
      DO I=1,3
         S = S + V(I)**2
      ENDDO
      S = SQRT(S)
      DO I=1,3
         V(I) = V(I)/S
      ENDDO
      RETURN
      END
C----Vector cross product
      SUBROUTINE V_CROSS (V, A, B)               ! V = A x B
      DIMENSION V(3), A(3), B(3), IN1(3), IN2(3)
      DATA IN1 / 2,1,1/
      DATA IN2 /3,3,2/
      X = -1.0
      DO I=1,3
         X = -X
         V(I) = X*(A(IN1(I))*B(IN2(I)) - 
     X             B(IN1(I))*A(IN2(I)))
      ENDDO
      CALL N_VECTOR(V)
      RETURN
      END
C
C-----------------------------------------------------------------------
C 
C     Convert from CHEMX to MOLPAK format...take a molecule on a center,
C       determine coordinates for the middle, average center related
C       coordinates if within 0.1 Angs of each other and move middle
C       to 0,0,0
C
      SUBROUTINE CHEMX_CI_MOLPAK
      DIMENSION IK(200), XYZ(3,200), IFLAG(200)
      CHARACTER NAME*5, NAME2(100)*5
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ_CART(3,200),
     x              N_NEAR(200), N_ATOMS(6,200)
      DIMENSION T(3,3)                              ! added, 2/18/01
      CHARACTER OUTFILE2(32)*1, INFILE2(32)*1 
      CHARACTER*32 INFILE, OUTFILE, MOLPAKXYZ
      DATA MOLPAKXYZ /'molpak.xyz                      '/
      EQUIVALENCE (OUTFILE,OUTFILE2), (INFILE, INFILE2),
     x            (N_NEAR, IFLAG)
      DATA OUTFILE2 /32*' '/
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
C
      PRINT 1
 1    FORMAT(' CHEMX input file name: ',$) 
      READ (5,2) INFILE
    2 FORMAT(A32)
C----OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
C----Output file
      PRINT 11
11    FORMAT (' MOLPAK output coordinate file name [molpak.xyz]: ',$)
      READ (5,2) OUTFILE
      IF (OUTFILE2(1) .EQ. ' ') OUTFILE = MOLPAKXYZ     
72    OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
C----Read cell params from lines 1 and 2 or CHEMX file
      READ (10,770) CELL
770   FORMAT (38X,3F8.3/21X,3F8.3)
      READ (10,*)         ! skip 2 lines
      READ (10,*)                       
      NA=0
777   NA=NA+1
         READ (10,30,END=999) IK(NA), NAME(NA)(1:1), (XYZ(I,NA),I=1,3)
30       FORMAT (I4,1X,A1,4X,3F10.5)
      GO TO 777
  999 NA=NA-1
C----Conversion from fractional to cartesian coordinates
          CALL MATRXT (CELL, T)
          CALL MTIMES1 (NA, CELL, XYZ, T, XYZ_CART)
C-----Should the C-H, N-H and/or N-O vectors be adjusted?
997   CALL ADJUST_2
      X=0.0                     ! move center of molecule
      Y=0.0                     ! to 0, 0, 0
      Z=0.0
      N_ADDED = 0
      DO 990 J=1,NA
       IF (NAME(J)(1:1) .EQ. 'H') GO TO 990   ! add up all
       N_ADDED = N_ADDED + 1                  ! atoms except
       X=X+XYZ_CART(1,J)                      ! hydrogen
       Y=Y+XYZ_CART(2,J)
       Z=Z+XYZ_CART(3,J)
990   CONTINUE
      X = X/N_ADDED
      Y = Y/N_ADDED
      Z = Z/N_ADDED
      DO J=1,NA
        XYZ_CART(1,J)=XYZ_CART(1,J)-X
        XYZ_CART(2,J)=XYZ_CART(2,J)-Y
        XYZ_CART(3,J)=XYZ_CART(3,J)-Z
      ENDDO
         do n=1,na
         write (23, 30) IK(N), NAME(N)(1:1), (XYZ(I,N),I=1,3)
         enddo
C----Determine which atoms are related by the center (within
C     0.1 Angs) and average the x,y,z and -x,-y,-z coords 
      NEW = 0
      DO 620 I=1,NA-1
         IF (XYZ_CART(1,I) .GE. 9998.) GO TO 620    ! forget it flag
         DO 610 J=I+1,NA
            IF (XYZ_CART(1,J) .GE. 9998.) GO TO 610
            S = 0.0
            DO 605 K=1,3                                   ! dist between
               S = S + (XYZ_CART(K,J) + XYZ_CART(K,I))**2  ! xyz & -x-y-z
605         CONTINUE                                       ! atoms
            IF (SQRT(S) .GE. 0.1) GO TO 610
            IF (NAME(I)(1:1) .NE. NAME(J)(1:1)) GO TO 610  ! atom types must
            NEW = NEW + 1                                  ! be the same
            NAME2(NEW) = NAME(I)
            DO 606 K=1,3
               XYZ(K,NEW) = 0.5*(XYZ_CART(K,I) - XYZ_CART(K,J))
606         CONTINUE
            XYZ_CART(1,J) = 99999.0          ! forget it flags for Ith and 
            XYZ_CART(1,I) = 99999.0          ! Jth atoms
610      CONTINUE
620   CONTINUE               
C----Choose the coords (xyz or -x-y-z) which
C     makes each atom within 1.8 Angs of another
      DO 649 I=2,NEW
         IFLAG(I) = 0
649   CONTINUE
      IFLAG(1) = 1
      DO 650 I=2,NEW
         J = 1
648      J = J + 1
         IF (IFLAG(J) .EQ. 1) GO TO 648
         DO 645 L=1,NEW
            IF (J .EQ. L) GO TO 645
            IF (IFLAG(L) .EQ. 0) GO TO 645
            S = 0.0
            DO 641 M=1,3                          ! Jth atom at
               S = S + (XYZ(M,L) - XYZ(M,J))**2   ! x,y,z
641         CONTINUE
            IF (SQRT(S) .GE. 1.8) GO TO 643
            IFLAG(J) = 1
            GO TO 650
643         S = 0.0
            DO 644 M=1,3
               S = S + (XYZ(M,L) + XYZ(M,J))**2   ! Jth atom at
644         CONTINUE                              ! -x-y-z
            IF (SQRT(S) .GE. 1.8) GO TO 645
            IFLAG(J) = 1
            DO 654 M=1,3
               XYZ(M,J) = - XYZ(M,J)
654         CONTINUE
            GO TO 650
645      CONTINUE 
650   CONTINUE
C----Output MOLPAK format
821   DO 810 I=1,NEW
         L = NPOT(I)
         WRITE (11,811) NAME2(I)(1:1), I, (XYZ(J,I),J=1,3), L,
     X                  G03_CHARGE, XMNDO_CHARGE 
811      FORMAT('ATOM ',A1,I4,1X,3F10.6,I5,2F10.6)    ! MOLPAK format
810   CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C 
C     Convert from MOLPAK format for a model with Ci symmetry.  Take a 
C       molecule on a center,  determine coordinates for the middle, 
C       average center related coordinates if within 0.05 Angs of each 
C       other and move middle to 0,0,0. Output file will contain only
C       1/2 molecule and charges will be averaged.
C
      SUBROUTINE MOLPAK_CI_MOLPAK             ! 2/2/97
C
      DIMENSION IK(200), XYZ(3,200), IFLAG(200), GCHARGES(200),
     X          XMNDO_CHARGES(200), GCHAR(200), XMNDO_CHAR(200),
     x          IKNEW(200)
      CHARACTER NAME*5, NAME2(100)*5
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ_CART(3,200),
     x              N_NEAR(200), N_ATOMS(6,200)
      CHARACTER OUTFILE2(32)*1, INFILE2(32)*1 
      CHARACTER*32 INFILE, OUTFILE, MOLPAKXYZ
      DATA MOLPAKXYZ /'molpak.xyz                      '/
      EQUIVALENCE (OUTFILE,OUTFILE2), (INFILE, INFILE2),
     x            (N_NEAR, IFLAG)
      DATA OUTFILE2 /32*' '/
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
C
      PRINT 1
 1    FORMAT(' MOLPAK input file name: ',$) 
      READ (5,2) INFILE
    2 FORMAT(A32)
C----OPEN THE INFILE 
      OPEN(UNIT=10,FILE=INFILE,STATUS='OLD')
C----Output file
      PRINT 11
11    FORMAT (' MOLPAK 1/2 molecule output coordinate file name ',
     x        '[molpak.xyz]: ',$)
      READ (5,2) OUTFILE
      IF (OUTFILE2(1) .EQ. ' ') OUTFILE = MOLPAKXYZ     
72    OPEN(UNIT=11,FILE=OUTFILE,STATUS='UNKNOWN') 
C----Read MOLPAK coordinates & charges
      NA=0
777   NA=NA+1
         READ (10,30,END=999) NAME(NA), (XYZ_CART(I,NA),I=1,3), IK(NA),
     x                        GCHARGES(NA), XMNDO_CHARGES(NA)
30       FORMAT (5X,A5,1X,3F10.6,I5,2F10.6)   ! MOLPAK format
      GO TO 777
  999 NA=NA-1
C
      X=0.0                     ! move center of molecule
      Y=0.0                     ! to 0, 0, 0
      Z=0.0
      DO 990 J=1,NA
       X=X+XYZ_CART(1,J)   
       Y=Y+XYZ_CART(2,J)
       Z=Z+XYZ_CART(3,J)
990   CONTINUE
      X = X/NA
      Y = Y/NA
      Z = Z/NA
      DO J=1,NA
        XYZ_CART(1,J)=XYZ_CART(1,J)-X
        XYZ_CART(2,J)=XYZ_CART(2,J)-Y
        XYZ_CART(3,J)=XYZ_CART(3,J)-Z
      ENDDO
C----Determine which atoms are related by the center (within
C     0.1 Angs) and average the x,y,z and -x,-y,-z coords 
      NEW = 0
      DO 620 I=1,NA-1
         IF (XYZ_CART(1,I) .GE. 9998.) GO TO 620    ! forget it flag
         DO 610 J=I+1,NA
            IF (XYZ_CART(1,J) .GE. 9998.) GO TO 610
            S = 0.0
            DO 605 K=1,3                                   ! dist between
               S = S + (XYZ_CART(K,J) + XYZ_CART(K,I))**2  ! xyz & -x-y-z
605         CONTINUE                                       ! atoms
            IF (SQRT(S) .GE. 0.01) GO TO 610
            IF (NAME(I)(1:1) .NE. NAME(J)(1:1)) GO TO 610  ! atom types &
            IF (IK(I) .NE. IK(J)) GO TO 610                ! names must
            NEW = NEW + 1                                  ! be the same
            NAME2(NEW) = NAME(I)
            IKNEW(NEW) = IK(I) 
            DO 606 K=1,3
               XYZ(K,NEW) = 0.5*(XYZ_CART(K,I) - XYZ_CART(K,J))
606         CONTINUE
            GCHAR(NEW) = 0.5*(GCHARGES(I) + GCHARGES(J))
            XMNDO_CHAR(NEW) = 0.5*(XMNDO_CHARGES(I) + 
     x                              XMNDO_CHARGES(J))
            XYZ_CART(1,J) = 99999.0          ! forget it flags for Ith and 
            XYZ_CART(1,I) = 99999.0          ! Jth atoms
610      CONTINUE
620   CONTINUE               
C----Choose the coords (xyz or -x-y-z) which
C     makes each atom within 1.8 Angs of another
      DO 649 I=2,NEW
         IFLAG(I) = 0
649   CONTINUE
      IFLAG(1) = 1
      DO 650 I=2,NEW
         J = 1
648      J = J + 1
         IF (IFLAG(J) .EQ. 1) GO TO 648
         DO 645 L=1,NEW
            IF (J .EQ. L) GO TO 645
            IF (IFLAG(L) .EQ. 0) GO TO 645
            S = 0.0
            DO 641 M=1,3                          ! Jth atom at
               S = S + (XYZ(M,L) - XYZ(M,J))**2   ! x,y,z
641         CONTINUE
            IF (SQRT(S) .GE. 1.8) GO TO 643
            IFLAG(J) = 1
            GO TO 650
643         S = 0.0
            DO 644 M=1,3
               S = S + (XYZ(M,L) + XYZ(M,J))**2   ! Jth atom at
644         CONTINUE                              ! -x-y-z
            IF (SQRT(S) .GE. 1.8) GO TO 645
            IFLAG(J) = 1
            DO 654 M=1,3
               XYZ(M,J) = - XYZ(M,J)
654         CONTINUE
            GO TO 650
645      CONTINUE 
650   CONTINUE
C----Output MOLPAK format...for just 1/2 of the input molecule
821   DO 810 I=1,NEW
         WRITE (11,811) NAME2(I), (XYZ(J,I),J=1,3), IKNEW(I),
     X                  GCHAR(I), XMNDO_CHAR(I) 
811      FORMAT('ATOM ',A5,1X,3F10.6,I5,2F10.6)     ! MOLPAK format
810   CONTINUE
      RETURN
      END
C------------------------------------------------------------------------
C
C     Determine atom types for potential parameter selection for 
C      inclusion in a modified MOLPAK input file
C
      FUNCTION NPOT (NUM_ATOM)
      CHARACTER (LEN=1) :: ELEMENT                            ! combo
      CHARACTER NAME*5, YN*1, MMOD(64)*2, AT2*2, LINE*60
      CHARACTER*32 INFILE,OUTFILE
      CHARACTER  INPUT_FORMAT*60, CHARGE_FILE*30, G03_LINE*21,
     x           MOPAC_LINE*31, MOPAC1*6, CIN*1
      LOGICAL CONNECT_FLAG
      DIMENSION IK(200), NUM_O(4)
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
C
C----Does the atom connectivty need to be established...once is enough
      IF (CONNECT_FLAG) GO TO 7         
          CALL CONNECT (NA)
          CONNECT_FLAG = .TRUE.
C----Are charges to be used?
7     IF (ICHARGE_FLAG .GE. 0) GO TO 10      
          PRINT 8
8         FORMAT (' Are atomic charges to be added to the MOLPAK',
     x            ' coordinate file?'/
     x            '   0 or blank = no; 1 = from Gaussian calcns;',
     x            ' 2 = from MOPAC MNDO/ESP calcns,'/
     x            '   3 = from both Gaussian and MNDO calcns [1]:',$)
          READ 28, CIN
28        FORMAT (A1)
          IF (CIN .EQ. ' ') THEN     ! default = G94 charges
             ICHARGE_FLAG = 1
          ELSE
             READ (CIN,'(I1)') ICHARGE_FLAG
          ENDIF
      GO TO (10,701,702,703), ICHARGE_FLAG + 1
701   CALL OPEN_G03        ! open g94 charge file on unit # 17
         GO TO 10
702   CALL OPEN_MNDO       ! open MNDO charge file on unit # 18
         GO TO 10
703   CALL OPEN_G03
      CALL OPEN_MNDO
C----Determine atom types
10    AT2 = NAME(NUM_ATOM)(1:2) 
      NPOT = 0
C----look at H's first
1375    IF (AT2(1:1) .NE. 'H') GO TO 1340              ! H
           IT = ICON(1,NUM_ATOM)     ! number of atom to which H is attached
C---------Check for alcohol H                                       3/25/03
           IF (.NOT. (NAME(IT)(1:1) .EQ. 'O' .AND.                ! 3/25/03
     2                NCON(IT) .EQ. 2)) GO TO 1378                ! 3/25/03
           NPOT = IQ_ALCOHOL(NUM_ATOM)                            ! 3/25/03
              IF (NPOT .NE. 0) GO TO 1303                         ! 3/25/03
1378       IF (NAME(IT)(1:1) .EQ. 'C') THEN                       ! 3/25/03
              NC = NCON(IT)             ! it's C, how many atoms bonded to C?
              GO TO (101, 101, 103, 104), NC
103           NPOT = 1                         ! H bonded to sp2 C; param = 1
                 GO TO 1303
104           IF (IQ_CUBANE(IT) .EQ. 0) THEN  ! 4 bonded...test for cubane     
                 NPOT = 2      ! no, H bonded to normal sp3 C; param = 2
              ELSE
                 NPOT = 47     ! H bonded to cubane skeleton; param = 2 
              ENDIF            ! CAN BE CHANGED LATER; 11-12-02 changed
              GO TO 1303
           ENDIF
           IF (NAME(IT)(1:1) .NE. 'N') GO TO 101
              NC = NCON(IT)             ! H is bonded to N
              GO TO (101, 101, 1103, 101), NC
C-----------Check for pyridinium-N(+)-H atom                                7/20/06
1103             if (iq_pyridinium_n(it) .eq. 85) then                    ! 7/20/06
                    npot = 86     ! H of pyridinium-N(+)-H, param = 86      7/20/06
                    go to 1303                                            ! 7/20/06
                 endif                                                    ! 7/20/06
C-----------Check on H linked to pyrrole-type N                              1/9/05
                 npot = iq_pyrrole_nh(it) ! it = number of 3-linked N        1/9/06
                 IF (NPOT .GT. 0) GO TO 1303  ! H of pyrrole, param = 74     1/9/06
C-----------Call IQ_HN_IMIDE to determine if the H is an imide NH,          12/3/03
C               H of a benzamide N or H or a hydantoin.                     12/3/03
                 NPOT = IQ_HN_IMIDE(NUM_ATOM,IT)  ! # of H and 3-linked N   1/15/04
                 IF (NPOT .GT. 0) GO TO 1303! H of an imide, param = 39     11/15/03
C                                             H or an amide, param = 61     11/15/03
C                                             H of a hydantoin, param = 64  12/3/03
                 NPOT = AMINE(NUM_ATOM,IT)    ! # of H and 3-linked N       1/15/04
                 IF (NPOT .GT. 0) GO TO 1303  ! H of C-NH2, param = 3       1/15/04
C                                               H of [Csp3]2-NH = 68        1/15/04
C---------Error call...can't ID the atom
101        PRINT 110, NAME(NUM_ATOM), NAME(IT)
110        FORMAT (1X,A5,' bonded to ',A5,', unknown type')
           GO TO 2002
C
C----Look at C's second
1340    IF (AT2(1:1) .NE. 'C') GO TO 2340
           GO TO (1001, 1002, 1003, 1004), NCON(NUM_ATOM)
1001          PRINT 800, NAME(NUM_ATOM), NCON(NUM_ATOM)
800           FORMAT (1X,A5,' has',I2,' connections...does not compute')
              GO TO 2002
C------2 connections...possible -C#C-, =C= or -C#N
C------Check for internal alkyne, C-C#C-C, first                   12/31/03
1002       NPOT = IQ_ALKYNE(NUM_ATOM)                            ! 12/31/03
              IF (NPOT .NE. 0) GO TO 1303  ! alkyne C, param = 67  12/31/03
           DO 2701 I1=1,2
              I2 = ICON(I1,NUM_ATOM)              ! attached atom
                 IF (.NOT. (NCON(I2) .EQ. 1 .AND. 
     x                      NAME(I2)(1:1) .EQ. 'N')) THEN
                    NPOT = 34       ! nitrile C
                    GO TO 1303
                 ENDIF
2701       CONTINUE
           PRINT 2705, NAME(NUM_ATOM), NCON(NUM_ATOM)  ! not the C of a nitrile
2705       FORMAT (1X,A5,' has',I2,' connections...what is it?')
           GO TO 2002
C---------3-linked C
1003       NSO1 = 0        ! # of connected 1-linked O's          ! 2/3/06
           nso2 = 0        ! # of connected 2-linked O's          ! 2/3/06
           DO 7338 I=1,3
              IT = ICON(I,NUM_ATOM)
              IF (name(IT)(1:1) .EQ. 'O' .AND. NCON(IT) .EQ. 1) 
     2            NSO1 = NSO1 + 1
              IF (name(IT)(1:1) .EQ. 'O' .AND. NCON(IT) .EQ. 2) 
     2            NSO2 = NSO2 + 1
7338       CONTINUE
           if ((nso1 .le. 1 .and. nso2 .eq. 0) .or.               ! 2/7/06
     2         (nso1 .eq. 0 .and. nso2 .eq. 1)) then              ! 2/7/06
              NPOT = 4      ! alkene, aldehyde, ketone, or          2/7/06
              GO TO 1303    !   vinyl ether C...sp2, param = 4      2/7/06
           ENDIF
           if (nso1 + nso2 .eq. 2) then ! check for ester/anhydride C 2/3/06
              npot = iq_ester_c(num_atom)   ! param = 79          ! 2/3/06
              if (npot .gt. 0) go to 1303
           endif
           GO TO 2007        ! can't id 3-linked C                  2/3/06
C----------4-linked C
1004          IF (IQ_CUBANE(NUM_ATOM) .EQ. 0) THEN ! 4 bonded...test for cubane
                 NPOT = 5      ! no, alkane C...sp3; param = 5
              ELSE
                 NPOT = 46     ! cubane C; param = 5; 11-12-02 changed 
              ENDIF            ! CAN BE CHANGED LATER
                 GO TO 1303        
C
C----Look at O's third
2340    IF (AT2(1:1) .NE. 'O') GO TO 1024
           GO TO (1021, 1022), NCON(NUM_ATOM)
1021          IT = ICON(1,NUM_ATOM)     ! get ID of single connected atom
C----------Check if O of a NO3-                                    7/20/06
              if (name(it)(1:1) .eq. 'N' .and.                  ! 7/20/06 
     2            ncon(it) .eq. 3) then                          ! 7/20/06
                 no1 = 0  ! number of 1-linked O's               ! 7/20/06
                 do i=1,3                                        ! 7/20/06
                    it2 = icon(i,it)                             ! 7/20/06
                    if (name(it2)(1:1) .eq. 'O' .and.           ! 7/20/06
     2                  ncon(it2) .eq.  1) no1 = no1 + 1         ! 7/20/06
                 enddo                                           ! 7/20/06
                 if (no1 .eq. 3) then                            ! 7/20/06
                    npot = 84  ! param = 84, O of NO3-           ! 7/20/06
                    go to 1303                                   ! 7/20/06
                 endif                                           ! 7/20/06
              endif                                              ! 7/20/06
C----------Check if linked to ester or anhydride carbonyl C         2/3/06
              if (name(it)(1:1) .eq. 'C' .and. ncon(it) .eq. 3    ! 2/3/06
     2            .and. iq_ester_c(it) .eq. 79) then              ! 2/3/06
                 npot = 80   ! param = 80, ester carbonyl O       ! 2/3/06
                 go to 1303                                       ! 2/3/06
              endif                                               ! 2/3/06
C---------Check for O of sulfonimine
           IF (NAME(IT)(1:1) .EQ. 'S' .AND. NCON(IT) .EQ. 4) THEN
              ELEMENT = 'O'
              NPOT = IQ_SULFONIMINE(NUM_ATOM,ELEMENT) ! sulfonimine O
              IF (NPOT .GT. 0) GO TO 1303             ! param = 37
           ENDIF
C-------Check if the O is connected to a 3-linked C
         IF (NAME(IT)(1:1) .EQ. 'C' .AND. NCON(IT) .EQ. 3) THEN
            NPOT = IQ_OC_IMIDE(IT)   ! imide or amide carbonyl O?      11/15/03
            IF (NPOT .GT. 0) GO TO 1303    ! O of imide, param = 41    11/15/03
C                                            O of amide, param = 63    11/15/03
C                                            O of hydantoin, param = 65 11/15/03
         ENDIF
C------------Check for carbonyl C with only 1-linked O, like ketone
              IF (NAME(IT)(1:1) .EQ. 'C') THEN
                 NPOT = 15        ! carbonyl O; param = 15
                 GO TO 1303
              ENDIF
C------------Check for nitroso ...N=O
              IF (NAME(IT)(1:1) .EQ. 'N' .AND. NCON(IT) .EQ. 2) THEN
                 NPOT = 25      ! X-N=O
                 GO TO 1303 
              ENDIF
C------------Check for N-oxide O of furoxan
              NPOT = IQ_FUROXAN_O1(NUM_ATOM)   ! param = 23
                 IF (NPOT .NE. 0) GO TO 1303 
C
              IF (NAME(IT)(1:1) .EQ. 'N' .AND.        ! N is linked to the O
     x               NCON(IT) .EQ. 3) THEN            ! & has 3 connections
                 N_OXY = 0
                 N_CARBS = 0
                 DO I=1,3                   ! count # of O's & C's linked
                    ITT = ICON(I,IT)        ! to the 2nd atom
                    IF (NAME(ITT)(1:1) .EQ. 'C') N_CARBS = N_CARBS + 1
                    IF (NAME(ITT)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
                 ENDDO
                 IF (N_OXY .EQ. 2 .OR. N_OXY .EQ. 3) THEN ! probably a nitro O; ONO2 included 11-12-02
                   NPOT = IQ_NITRO_OX(IT)
                   IF( NPOT .EQ. 49) NAME(NUM_ATOM)(1:2) = 'OX'    ! 11-12-02 
C                   PRINT *, NUM_ATOM,  IT, NPOT, NAME(NUM_ATOM)(1:2)         ! 11-12-02
                   IF (NPOT .EQ. 0) GO TO 101 ! can't identify the kind of O
                 ELSE
                   NPOT = 17    ! 2 C's & 1 O = a pyridine N-oxide O
                 ENDIF          ! ...param = 17
                 GO TO 1303
              ENDIF
C----Check for 2-linked O or nitrate ester
1022          NPOT = IQ_OXYGEN (NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
              NPOT = IQ_ALCOHOL (NUM_ATOM)                     ! 3/25/03
                 IF (NPOT .NE. 0) GO TO 1303                   ! 3/25/03
C----Check for the ring O of a furoxan...2-linked
              NPOT = IQ_FUROXAN_O2(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303       ! param = 22 
C----Check for other 2-linked O's...C-O-C, X-O-B *********************           
              IT1 = ICON(1,NUM_ATOM)
              IT2 = ICON(2,NUM_ATOM)
              NC = 0
              NN = 0
              NB = 0
              IF (NAME(IT1)(1:1) .EQ. 'B') NB = NB + 1
              IF (NAME(IT2)(1:1) .EQ. 'B') NB = NB + 1
              IF (NAME(IT1)(1:1) .EQ. 'C') NC = NC + 1
              IF (NAME(IT2)(1:1) .EQ. 'C') NC = NC + 1
              IF (NAME(IT1)(1:1) .EQ. 'N') NN = NN + 1
              IF (NAME(IT2)(1:1) .EQ. 'N') NN = NN + 1
C----Check for "ether" O of an ester or andydride               2/3/06
              If (nc .eq. 2 .and. (iq_ester_c(it1) .eq. 79    ! 2/3/06
     2            .or. iq_ester_c(it2) .eq. 79)) then         ! 2/3/06
                 npot = 81   ! "ether" O of an ester            2/3/06
                 go to 1303                                   ! 2/3/06
              endif                                           ! 2/3/06
              IF (NC .EQ. 2 .OR. (NC .EQ. 1 .AND. NN .EQ. 1) .OR.
     x           (NC .EQ. 1 .AND. NB .EQ. 1))  THEN
                     NPOT = 14      ! ether O; C-O-C, C-O-N or C-O-B
              ELSE
                  NPOT = 999     ! unknown
              ENDIF
              GO TO 1303 
C
C----Look at N's fourth
1024    IF (AT2(1:1) .NE. 'N') GO TO 2002
           GO TO (1031, 1032, 1033, 2034), NCON(NUM_ATOM)
C------------Check for terminal azide N                          ! 3/5/03
1031          NPOT = IQ_AZIDE_N (NUM_ATOM)                       ! 3/5/03
                 IF (NPOT .NE. 0) GO TO 1303                     ! 3/5/03
C------------Possible nitrile N...check to see it there's a 2-linked C attached
              IT2 = ICON(1,NUM_ATOM)   ! IT2 = attached atom     ! 3/5/03
              IF (NCON(IT2) .EQ. 2 .AND. NAME(IT2)(1:1) .EQ. 'C') THEN
                  NPOT = 33       ! nitrile N, param = 33
              ELSE
                  NPOT = 999      ! flag as unknown...1-linked N but not a nitrile
              ENDIF
              GO TO 1303
C----N with 4 connections....like NH4+
2034          NPOT = 999         ! flag as unknown
                 GO TO 1303
C----N with 2 connections
C------------Check for 2-linked N of amino-imine                1/24/06
1032          npot = iq_aminoimine_2(num_atom)                ! 1/24/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 76   ! 1/24/06
C------------Check for 2-linked N of thiazole                   1/31/06
              npot = iq_thiazole_n(num_atom)                  ! 1/31/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 78   ! 1/31/06
C------------Check for 2-linked N (N3) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n3(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 72    ! 1/1/06
C------------Check for 2-linked N (N2) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n2(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 71    ! 1/1/06
C------------Check for 2-linked N (N4) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n4(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 73    ! 1/1/06
C------------Check for 2-linked N of 123triazole                 12/30/05
              npot = iq_123triazole_2(num_atom)                ! 12/30/05
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 27    ! 12/30/05
C------------Check for internal azide N's                        ! 3/5/03
              NPOT = IQ_AZIDE_N (NUM_ATOM)                       ! 3/5/03
                 IF (NPOT .NE. 0) GO TO 1303                     ! 3/5/03
              NPOT = IQ_N_AZAPENTALENE(NUM_ATOM)                 ! 3/5/03
                 IF (NPOT .NE. 0) GO TO 1303
              NPOT = IQ_FUROXAN_N2(NUM_ATOM)   ! param = 20
                 IF (NPOT .NE. 0) GO TO 1303
C------------Check for diazo N....C-N=N-C                         2/13/06
              NPOT = IQ_diazo_n(NUM_ATOM)      ! param = 82     ! 2/13/06
                 IF (NPOT .NE. 0) GO TO 1303                    ! 2/13/06
C                      
1037          N_OXY = 0       ! 2 connections, how many O's are attached
              N_CARBS = 0     ! # of C's attached
              N_NIT = 0       ! # of N's attached
              N_HYD = 0       ! # of H's attached
              N_F = 0         ! # of F's attached
              N_S_4 = 0       ! # of 4-linked S's attached         ! 12/11/00
              DO I=1,2
                 N_ADJ = ICON(I,NUM_ATOM)
                 IF (NAME(N_ADJ)(1:1) .EQ. 'O') THEN
                     N_OXY = N_OXY + 1
                     NUM_O(N_OXY) = N_ADJ
                 ENDIF
                 IF (NAME(N_ADJ)(1:1) .EQ. 'C') N_CARBS = N_CARBS + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'N') N_NIT = N_NIT + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'H') N_HYD = N_HYD + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'F') N_F = N_F + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'S' .AND.                     ! 12/11/00 
     X                 NCON(N_ADJ) .EQ. 4)  N_S_4 = N_S_4 + 1            ! 12/11/00
              ENDDO
C
C------------Check on a 4-linked S attached to the 2-linked N.           ! 12/11/00
              IF (N_S_4 .GT. 0) THEN 
                  ELEMENT = 'N'                                   
                  NPOT = IQ_SULFONIMINE(NUM_ATOM, ELEMENT)  ! it's a N     
                  IF (NPOT .NE. 0) GO TO 1303                            ! 12/11/00
              ENDIF                                                      ! 12/11/00
c
              IF (N_OXY .EQ. 1 .AND. NCON(NUM_O(N_OXY)) .EQ. 1) THEN
                  NPOT = 24        ! N of X-N=O
                  GO TO 1303
              ENDIF
              IF (N_OXY .EQ. 2) THEN
                  N = MIN0(NCON(NUM_O(1)), NCON(NUM_O(2)))
                  IF (N .EQ. 1) THEN
                     NPOT = 24     ! N of -O-N=O
                     GO TO 1303
                  ENDIF
              ENDIF
C                  
              IF ((N_CARBS .EQ. 2) .OR. (N_CARBS .EQ. 1 .AND.
     2             N_NIT .EQ. 1)) THEN
                 NPOT = 16    ! pyridine N...param = 16
                 GO TO 1303
              ENDIF
              GO TO 2007
C----N with 3 connections
!------------Check for 3-linked N of aminoimine                 1/24/06
1033          npot = iq_aminoimine_3(num_atom)                ! 1/24/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 76   ! 1/24/06
C------------Check for 3-linked N (N1) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n1(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 70    ! 1/1/06
C------------Check for 3-linked N of 123triazole                 12/30/05
              npot = iq_123triazole_3(num_atom)                ! 12/30/05
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 26    ! 12/30/05
C------------Check for nitrate ester N and nitramine nitro N        2/24/03
              NPOT = IQ_NITRATE_N(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
              NPOT = IQ_N_AZAPENTALENE(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
C------------Check if central N of imide, benzamide N or a hydantoin    12/3/03
              NPOT = IQ_HN_IMIDE(NUM_ATOM,NUM_ATOM)       ! same atom #'s = a N
                 IF (NPOT .GT. 0) GO TO 1303! N of an imide, param = 40 12/3/03
C                                             N or an amide, param = 62 12/3/03
C                                             N of a hydantoin, param = 67 12/3/03
              NPOT = IQ_FUROXAN_N3(NUM_ATOM)   ! param = 21
                 IF (NPOT .NE. 0) GO TO 1303
1035          N_OXY = 0       ! 3 connections, how many O's are attached
              N_CARBS = 0     ! # of C's attached
              n_carbs_3 = 0   ! # of 3-linked C's
              N_NIT = 0       ! # of N's attached
              N_HYD = 0       ! # of H's attached
              N_F = 0         ! # of F's attached
              DO I=1,3
                 N_ADJ = ICON(I,NUM_ATOM)                 
C                PRINT *, NUM_ATOM, N_ADJ, NAME(N_ADJ)(1:1)                     ! 11-12-02
                 IF (NAME(N_ADJ)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'C') THEN
                    N_CARBS = N_CARBS + 1
                    NC_ADJ = N_ADJ                ! 11-12-02
                 END IF
                 if (name(n_adj)(1:1) .eq. 'C' .and.
     2               ncon(n_adj) .eq. 3) n_carbs_3 = n_carbs_3 + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'N') N_NIT = N_NIT + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'H') N_HYD = N_HYD + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'F') N_F = N_F + 1
              ENDDO
              IF (N_OXY .EQ. 2 .OR. N_OXY .EQ. 3) THEN   ! 11-12-02 included ONO2
C                PRINT *, N_OXY, N_CARBS, NUM_ATOM, NAME(NUM_ATOM)(1:1)        ! 11-12-02
                 IF (N_CARBS .EQ. 1) THEN
                  IF ( IQ_CUBANE(NC_ADJ) .EQ. 1) THEN
                      NPOT = 48    ! nitro N of nitrocubane   11-12-02 
                      GO TO 1303
                  ELSE 
                      NPOT = 7     ! N of C-NO2 except for nitrocubane
                      GO TO 1303   
                  END IF                                                                            
                 ELSE                                                                               
                    NPOT = 7     ! nitro N...param = 7; right now NO2's in                         
                    GO TO 1303   ! C-NO2 and N-NO2 are identical  
                 END IF                                                                             
              ENDIF                                                                                 
              IF (N_CARBS .EQ. 2 .AND. N_NIT .EQ. 1) THEN                                           
                 NPOT = 8    ! nitramine N...param = 8                                              
                 GO TO 1303                                                                         
              ENDIF                                                                                 
              IF (N_CARBS .EQ. 3) THEN                                                              
                 NPOT = 9   ! C3N, tertiary N...param = 9                                           
                 GO TO 1303                                                                         
              ENDIF                                                                                 
C----------Check for pyridinium-N....[Csp2]2-N(+)-H              ! 7/20/06
              IF (N_CARBS_3 .EQ. 2 .AND. n_hyd .EQ. 1) THEN      ! 7/20/06
                 NPOT = iq_pyridinium_n(num_atom)   ! pyridinium N(+)-H 7/20/06
                 if (npot .gt. 0) GO TO 1303        ! N param = 85      7/20/06
              ENDIF                                              ! 7/20/06
C----Test for N of [Csp3]2NH                                   1/15/04
              NPOT = AMINE(NUM_ATOM, NUM_ATOM)               ! 1/15/04
                 IF (NPOT .NE. 0) GO TO 1303   ! N PARAM = 68  1/15/04
              IF (N_HYD .GE. 1 .AND. N_CARBS .LE. 2) THEN 
                 NPOT = 10  !  N of NH2 or NH, param = 10   
                 GO TO 1303                                                                         
              ENDIF                                                                                 
              IF (N_F .GE. 1) THEN                                                                  
                 NPOT = 45  ! N of NHF or NF2 param = 45                                       
                 GO TO 1303                                                                         
              END IF        
              IF ((N_CARBS .EQ. 2 .AND. N_OXY .EQ. 1) .OR.
     2            (N_CARBS .EQ. 1 .AND. N_OXY .EQ. 1 .AND.
     3             N_NIT .EQ. 1) .OR.
     4            (N_NIT .EQ.2 .AND. N_OXY .EQ. 1)) THEN         
                 NPOT = 18  ! N of pyridine N-oxide, param = 18                                   
                 GO TO 1303  
              ENDIF 
       GO TO 2007                                                                                   
C----Look at Br's 5th
2002   IF (AT2 .NE. 'BR') GO TO 2008
          NPOT = 12                  ! Br, param = 12
          GO TO 1303
C----Look at F's 6th
2008   IF (AT2(1:1) .NE. 'F') GO TO 3009 
          N_ADJ = ICON(1,NUM_ATOM)
          IF (NAME(N_ADJ)(1:1) .EQ. 'C') THEN
             NPOT = 53               ! F of C-F, param = 53
          ELSE
             NPOT = 13               ! F of N-F, param = 13 
          ENDIF 
          GO TO 1303 
C----Look at B's 7th (O2-B-C)
3009   IF (AT2(1:1) .NE. 'B') GO TO 3020
              IF (NCON(NUM_ATOM) .NE. 3) GO TO 2007
              NC = 0
              NO = 0
              DO 3015 I=1,3
                 IT = ICON(I,NUM_ATOM)
                 IF (NAME(IT)(1:1) .EQ. 'O') NO = NO + 1 
                 IF (NAME(IT)(1:1) .EQ. 'C') NC = NC + 1
3015          CONTINUE
              IF  (NC .EQ. 1 .AND. NO .EQ. 2)  THEN
                   NPOT = 35      ! boron...O2-B-C
                   GO TO 1303 
              ENDIF
C
C----Look at S's 8th (look for C-SO2-N=C)
3020   IF (AT2(1:1) .NE. 'S') GO TO 3035                  ! 11/6/03 
          if (ncon(num_atom) .ne. 2) go to 4009           ! 1/31/06
             npot = iq_thiazole_s(num_atom)               ! 1/31/06
             if (npot .ne. 0) go to 1303  ! param = 77      1/31/06
4009   IF (NCON(NUM_ATOM) .NE. 4) GO TO 3035              ! 11/6/03 
          ELEMENT = 'S'
          NPOT = IQ_SULFONIMINE(NUM_ATOM, ELEMENT) ! it's a S 
          IF (NPOT .NE. 0) GO TO 1303                     ! 12/11/00
C
C----Look at I's 9th                                      ! 11/6/03 
3035   IF (AT2(1:1) .NE. 'I') GO TO 2007                  ! 11/6/03 
       IF (NCON(NUM_ATOM) .NE. 1) GO TO 2007              ! 11/6/03 
       NPOT = 60           ! I, param = 60                ! 11/6/03
       GO TO 1303                                         ! 11/6/03 
C
2007   PRINT 2004, NAME(NUM_ATOM), NUM_ATOM
2004   FORMAT (' Atom ',A5,I2,' is unknown...type set to 999')
1034      NPOT = 999     ! 4 connections ... flag as unknown 
1303   GO TO (2000,1500,1600,1700), ICHARGE_FLAG + 1
1500      READ (17,1502) G03_CHARGE      ! G03 charge
1502      FORMAT (9X,F12.6)                                       ! 5/31/03
       GO TO 2000
1600      READ (18,1602) XMNDO_CHARGE      ! mopac charge
1602      FORMAT (33x,f12.4)
       GO TO 2000
1700      READ (17,1502) G03_CHARGE
          READ (18,1602) XMNDO_CHARGE
2000   RETURN
       END
C-
C----------------------------------------------------------
C----Function to determine S, O and N atoms of sulfonimines      
C
C       O      O = 37               
C       |      S = 36                 
C     C-S-N=C  N = 38
C       |
C       O      O = 37 
C 
C
      FUNCTION IQ_SULFONIMINE (IS, ELEMENT)  
C                   ELEMENT = O, S or N (1,2,3)
C                   IS = input atom #
C     ELEMENT is 4-linked S, 1-linked O or 2-linked N
C
      CHARACTER NAME*5
      CHARACTER (LEN=1) :: ELEMENT       
      DIMENSION N_NAMES(2), NUM_O(2), NUM_N(3), NUM_N2(3)      ! 11/3/00
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_SULFONIMINE = 0
C
C----Sulfur???   must be 4-linked
      IF (ELEMENT .EQ. 'S' .AND. NCON(IS) .EQ. 4) THEN
         ISS = IS                                 ! ISS = S #
         GO TO 20
      ENDIF
C----Oxygen???   must be 1-linked
      IF (ELEMENT .EQ. 'O' .AND. NCON(IS) .EQ. 1) THEN   
         ISS = ICON(1,IS)         ! atom linked to O should be S
         IF (NAME(ISS)(1:1) .EQ. 'S' .AND. NCON(ISS) .EQ. 4) GO TO 20
      ENDIF 
C----Nitrogen???  must be 2-linked
      IF (ELEMENT .EQ. 'N' .AND. NCON(IS) .EQ. 2) THEN   
         DO K=1,2               ! one of the 2-linked atoms
            ISS = ICON(K,IS)    ! must be a 4-linked S
            IF (NAME(ISS)(1:1) .EQ. 'S' .AND. NCON(ISS) .EQ. 4) GO TO 20
         ENDDO
      ENDIF
      RETURN       ! can't make sense of anything
C
20    NS = NCON(ISS)          ! ISS is the 4-linked S
C
C----Count # C's, # 1 and # 2-linked O's and 3-linked N's bonded
C     to the 4-linked S
      NC = 0                            ! # C's linked to S    ! 3/20/00
      NO_1 = 0                          ! # 1-linked O's       ! 3/18/00
      NN_2 = 0                          ! # 2-linked N's       ! 11/3/00
      DO 200 I=1,NS      
         I2 = ICON(I,ISS)                                       ! 3/18/00
         IF (NAME(I2)(1:1) .EQ. 'C') NC = NC + 1               ! 3/18/00
         IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 1) 
     x      NO_1 = NO_1 + 1                                    ! 3/18/00
         IF (NAME(I2)(1:1) .EQ. 'N' .AND. NCON(I2) .EQ. 2) THEN! 11/3/00
            NN_2 = NN_2 + 1         ! looking for 2-linked N's 
            NUM_N2(NN_2) = I2
         ENDIF                                                 ! 11/3/00
200   CONTINUE                                                 ! 3/18/00
C
      IF (NC .NE. 1 .OR. NO_1 .NE. 2 .OR. NN_2 .NE. 1) RETURN
C
C-------Sulfonimine...C-SO2-N=C ??    
222         KN = NUM_N2(1)       ! check imine N &    
            DO 224 J=1,2         !   look for 3-linked C 
               IC = ICON(J,KN)                                 ! 11/3/00
               IF (NAME(IC)(1:1) .EQ. 'S') GO TO 224           ! 11/3/00
               IF (NAME(IC)(1:1) .EQ. 'C' .AND.                ! 11/3/00
     X             NCON(IC) .EQ. 3) GO TO 226      ! ok        ! 11/3/00
               RETURN                              ! not ok    ! 11/3/00
224         CONTINUE                                           ! 11/3/00   
226         IF (ELEMENT .EQ. 'O') IQ_SULFONIMINE = 37  ! param = 37   
            IF (ELEMENT .EQ. 'N') IQ_SULFONIMINE = 38  ! param = 38
            IF (ELEMENT .EQ. 'S') IQ_SULFONIMINE = 36  ! param = 36
      RETURN
      END
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C 
C----IQ_ALKYNE   determine is a 2-linked C is an internal      12/31/03
C                alkyne C                                      12/31/03 
C
      FUNCTION IQ_ALKYNE (NUM_ATOM)                          ! 12/31/03
C 
C---- NUM_ATOM = a 2-linked C                                ! 12/31/03
C 
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      DIMENSION NUM_C(2)                                     ! 12/31/03
C 
      IQ_ALKYNE = 0        ! intialize to NO                   12/31/03
C 
C----Count # of 2-linked C's and total C's                     12/31/03
      NC = 0                                                 ! 12/31/03
      NC_TOTAL = 0                                           ! 12/31/03
      DO 100 I=1,2                                            ! 12/3/03
         L = ICON(I,NUM_ATOM)                                 ! 12/3/03    
         IF (NAME(L)(1:1) .EQ. 'C') NC_TOTAL = NC_TOTAL + 1  ! 12/31/03
         IF (NAME(L)(1:1) .EQ. 'C' .AND. NCON(L) .EQ. 2) THEN  ! 12/31/03
            NC = NC + 1                                      ! 12/31/03
            NUM_C(NC) = L                                    ! 12/31/03
         ENDIF                                               ! 12/31/03
100   CONTINUE                                                ! 12/3/03
C
      IF (NC .NE. 1 .OR. NC_TOTAL .NE. 2) RETURN             ! 12/31/03
C----Check the other C in NUM_C to see if it is 2-linked       12/31/03
      L = NUM_C(1)                                           ! 12/31/03
      IF (NCON(L) .NE. 2) RETURN                             ! 12/31/03
      IQ_ALKYNE = 67        ! internal alkyne, param = 67      12/31/03
      RETURN                                                 ! 12/31/03
C
      END                                                    ! 12/31/03
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C 
C----N_CARBONYL...determine if a 3-linked C is a carbonyl C     12/3/03 
C
      FUNCTION N_CARBONYL (NUM_ATOM)                          ! 12/3/03
C 
C----If NUM_ATOM = a 3-linked C, N_CARBONYL = 0/1 for no/yes    12/3/03
C       itis a carbonyl                                         12/3/03
C 
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C 
      N_CARBONYL = 0        ! intialize to NO                   12/3/03
C 
C----Count # of 3-linked and total C's on the N                 12/3/03
      DO 100 I=1,3                                            ! 12/3/03
         L = ICON(I,NUM_ATOM)                                 ! 12/3/03    
         IF (NAME(L)(1:1) .EQ. 'O' .AND. NCON(L) .EQ. 1) THEN ! 12/3/03
             N_CARBONYL = 1     ! Yes, a carbonyl             ! 12/3/03
         ENDIF                                                ! 12/3/03
100   CONTINUE                                                ! 12/3/03
C
      END                                                     ! 12/3/03
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C 
C----N_3LINK...count # of 3-linked N's attached to NUM_ATOM     12/3/03 
C
      FUNCTION N_3LINK (NUM_ATOM)                             ! 12/3/03
C 
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C 
      N_3LINK = 0                                             ! 12/3/03
C 
C----Count # of 3-linked N's                                    12/3/03
      DO 100 I=1,NCON(NUM_ATOM)                               ! 12/3/03
         L = ICON(I,NUM_ATOM)                                 ! 12/3/03    
         IF (NAME(L)(1:1) .EQ. 'N' .AND. NCON(L) .EQ. 3) THEN ! 12/3/03
             N_3LINK = N_3LINK + 1   ! yes, a 3-linked N        12/3/03
         ENDIF                                                ! 12/3/03
100   CONTINUE                                                ! 12/3/03
C
      END                                                     ! 12/3/03
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C 
C----AMINE...                                                   1/15/04
C  
C         N,C---NH2              [Csp3]2---N-H                  1/27/06
C            H = 3              N = 69     H = 68               1/27/06
C
      FUNCTION AMINE(NUM_ATOM,IT)                             ! 1/15/04
C
C----If NUM_ATOM .ne. IT, NUM_ATOM = H; IT = 3-linked N         1/15/04
C    If NUM_ATOM = IT, IT and NUM_ATOM = 3-linked N             1/15/04
C 
      CHARACTER NAME*5
      LOGICAL TYPE_ATOM                                       ! 1/15/04
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      AMINE = 0                                               ! 1/15/04 
      TYPE_ATOM = .true.        ! set signal to H               1/15/04
      IF (NUM_ATOM .EQ. IT) TYPE_ATOM = .false. ! set signal to N 1/15/04
C 
C----Count # of C's and 4-linked C's & H's on the N             1/15/04
C
      nn = 0                                                  ! 1/27/06
      NC = 0   ! # of C's                                       1/15/04
      NC4 = 0  ! # of 4-linked C's                              1/15/04
      NH = 0   ! # of H's                                       1/15/04
      DO I=1,3                                                ! 1/15/04
         L = ICON(I,IT)   ! look at connections to N            1/15/04    
         if (name(l)(1:1) .eq. 'N') nn = nn + 1               ! 1/27/06
         IF (NAME(L)(1:1) .EQ. 'C') NC = NC + 1               ! 1/15/04 
         IF (NAME(L)(1:1) .EQ. 'H') NH = NH + 1               ! 1/15/04
         IF (NAME(L)(1:1) == 'C' .AND. NCON(L) == 4) NC4 = NC4 + 1 ! 1/15/04 
      ENDDO                                                   ! 1/15/04
      if (NH == 2 .AND. (NC == 1 .or. nn .eq. 1)) THEN        ! 7/20/06
            if (type_atom) AMINE = 3        ! R-NH2, param = 3 for H  7/20/06
            if (.not. type_atom) amine = 10 ! R-NH2, param = 10 for N 7/20/06
      ELSE                                                    ! 1/15/04
         IF (NH == 1 .AND. NC4 == 2) THEN                     ! 1/15/04
            IF (TYPE_ATOM) AMINE = 68         ! for H           1/15/04
            IF (.not. TYPE_ATOM) AMINE = 69   ! for N           1/15/04
         ENDIF                                ! [Csp3]2NH       1/15/04
      ENDIF
      RETURN                                                  ! 1/15/04
      END                                                      ! 1/15/04 
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C 
C----IQ_HN_IMIDE... identify H and N of imides,                 12/3/03 
C                   benzamides and hydantoins                   12/3/03 
C
C                    40                     62                  11/15/03
C               --C---N---C--      -phenyl---N---C--            12/3/03
C                 "   |   "                  |   "              11/15/03
C                 O   H   O                  H   O              11/15/03
C                41  39  41                 61  63              11/15/03
C
C                    66      66                                 11/15/03
C               --C---N---C---N---C--                           12/3/03 
C                 "   |   "   |                                 12/3/03     
C                 O   H   O   H                                 12/3/03
C                65  64  65  64                                 12/3/03 
C
      FUNCTION IQ_HN_IMIDE(NUM_ATOM,IN)                       ! 12/3/03
C 
C----If NUM_ATOM .NE. IN, atom is a H; if NUM_ATOM = IN, atom is a N 
C 
      CHARACTER NAME*5
      LOGICAL HYDROGEN
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      DIMENSION NUM_C(3), NUM_CARB(3), NUM_N(3)               ! 12/3/03
C 
      IQ_HN_IMIDE = 0 
      HYDROGEN = .true. 
C 
      IF (NUM_ATOM .EQ. IN) HYDROGEN = .false. 
C 
C----Count # of 3-linked, total and carbonyl C's on the N       12/3/03
      NC_3 = 0                                                ! 12/3/03 
      NC_TOTAL = 0                                            ! 12/3/03
      N_CARB = 0
      DO 100 I=1,3           ! count # C's on N                 12/3/03
         L = ICON(I,IN) 
         IF (NAME(L)(1:1) .EQ. 'C') NC_TOTAL = NC_TOTAL + 1  ! 12/3/03
         IF (NAME(L)(1:1) .EQ. 'C' .AND. NCON(L) .EQ. 3) THEN 
            NC_3 = NC_3 + 1                                   ! 12/3/03
            NUM_C(NC_3) = L         ! id #'s of 3-linked C's    12/3/03
            IF (N_CARBONYL(L) .EQ. 1) THEN   ! test for         12/3/03
               N_CARB = N_CARB + 1           ! carbonyl C &     12/3/03     
               NUM_CARB(N_CARB) = L          ! record number    12/3/03
            ENDIF                                             ! 12/3/03
         ENDIF                                                ! 12/3/03 
100   CONTINUE                                                ! 12/3/03
C
      IF (N_CARB == 0) RETURN    ! not a N-C(=O)               12/3/03
C
C----This is a test for an amide...phenyl-NH-C(=O)-C              12/3/03
      IF (N_CARB .EQ. 1 .AND. N_3LINK(NUM_CARB(1)) .EQ. 1     ! 12/3/03 
     1          .AND. NC_3 .EQ. 2) THEN                       ! 12/3/03
         IF (HYDROGEN) THEN            ! it is an amide         11/15/03
            IQ_HN_IMIDE = 61           ! H param = 61           11/15/03
         ELSE                                                 ! 11/15/03 
            IQ_HN_IMIDE = 62           ! N param = 62           11/15/03 
         ENDIF                                                ! 11/15/03 
         RETURN                                               ! 11/15/03
      ENDIF                                                   ! 11/15/03
C
C----Test for imides                                            12/3/03
      IF (N_CARB .EQ. 2) THEN   ! how many 3-linked N's         12/3/03
         N_3LINKSUM = 0                                       ! 12/3/03
         DO I=1,2                                             ! 12/3/03
            N_3LINKSUM = N_3LINKSUM + N_3LINK(NUM_CARB(I))    ! 12/3/03
         ENDDO                                                ! 12/3/03
         IF (N_3LINKSUM .EQ. 2) THEN  ! an imide                12/3/03
            IF (HYDROGEN) THEN            ! it is an imide      12/3/03
               IQ_HN_IMIDE = 39           ! H param = 39        12/3/03   
            ELSE                                              ! 12/3/03
               IQ_HN_IMIDE = 40           ! N param = 40        12/3/03
            ENDIF                                             ! 12/3/03
            RETURN                                            ! 12/3/03 
         ENDIF                                                ! 12/3/03
C
C----Test for central N of hydantoin, N_3LINKSUM = 3?           12/3/03
         IF (N_3LINKSUM .EQ. 3) THEN    ! yes                   12/3/03 
            IF (HYDROGEN) THEN            ! it is an N          12/3/03 
               IQ_HN_IMIDE = 64           ! H param = 64        12/3/03 
            ELSE                                              ! 12/3/03 
               IQ_HN_IMIDE = 66           ! N param = 66        12/3/03 
            ENDIF                                             ! 12/3/03 
            RETURN                                            ! 12/3/03
         ENDIF                                                ! 12/3/03
      ENDIF                                                   ! 12/3/03
C
C----Finally test for end N of hydantoin                       12/3/03
      IF (.not. (N_CARB .EQ. 1 .AND.                         ! 12/3/03
     1      N_3LINK(NUM_CARB(1)) .EQ. 2)) RETURN  ! not end N  12/3/03  
C----Look at 2 N's on carbonyl and get N id's                  12/3/03 
      NN = 0    ! count of # 3-linked N's on carbonyl          12/3/03 
      DO 200 I=1,3                                           ! 12/3/03 
         K = ICON(I,NUM_CARB(1))                             ! 12/3/03 
         IF (NAME(K)(1:1) .NE. 'N' .OR. NCON(K) .NE. 3) GO TO 200 ! 12/3/03
            NN = NN + 1                                      ! 12/3/03
            NUM_N(NN) = K   ! id of 3-linked N                 12/3/03
200   CONTINUE                                               ! 12/3/03
      IF (NN .NE. 2) RETURN  ! must have two 3-linked N's      12/3/03
      N_CARB = 0                                             ! 12/3/03
      DO 300 I=1,NN                                          ! 12/3/03
         DO 290 L=1,NCON(NUM_N(I))                           ! 12/3/03
            M = ICON(L,NUM_N(I))                             ! 12/3/03
            IF (NAME(M)(1:1) .NE. 'C' .OR. NCON(M) .NE. 3) GO TO 290 ! 12/3/03
            N_CARB = N_CARB + N_CARBONYL(M)                  ! 12/3/03
290      CONTINUE                                            ! 12/3/03
300   CONTINUE                                               ! 12/3/03 
      IF (N_CARB .NE. 3) RETURN                                ! 12/3/03
      IF (HYDROGEN) THEN             ! it is an N              12/3/03 
          IQ_HN_IMIDE = 64           ! H param = 64            12/3/03 
      ELSE                                                   ! 12/3/03 
          IQ_HN_IMIDE = 66           ! N param = 66            12/3/03 
      ENDIF                                                  ! 12/3/03 
      RETURN                                                 ! 12/3/03
C
      END                                                    ! 12/3/03
C 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
C
C----IQ_OC_IMIDE... identify O of imide, benzamide or hydantoin 11/15/03 
C 
C                    40                     62                  11/15/03
C               --C---N---C--       -phenyl--N--C--             12/3/03
C                 "   |   "                  |  "               11/15/03
C                 O   H   O                  H  O               11/15/03
C                41  39  41                 61 63               11/15/03
C
C                    66      66                                 11/15/03
C               --C---N---C---N---C--                           12/3/03 
C                 "   |   "   |                                 12/3/03     
C                 O   H   O   H                                 12/3/03
C                65  64  65  64                                 12/3/03 
C
      FUNCTION IQ_OC_IMIDE(IN) 
C 
C----Have a 1-linked O; IN is the 3-linked C                    11/15/03
C 
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      DIMENSION NUM_N(3)
C
      IQ_OC_IMIDE = 0 
C 
C----Count # of atoms linked to the 3-linked C                  11/15/03 
      NC = 0          ! count # C's 
      NN = 0          ! count # 3-linked N's 
      NO = 0          ! count # of 1-linked O's 
      DO 100 I=1,3 
         L = ICON(I,IN) 
         IF (NAME(L)(1:1) .EQ. 'N' .AND. NCON(L) .EQ. 3) THEN 
            NN = NN + 1 
            NUM_N(NN) = L                                     ! 11/15/03
         ENDIF 
         IF (NAME(L)(1:1) .EQ. 'C') NC = NC + 1 
         IF (NAME(L)(1:1) .EQ. 'O' .AND. NCON(L) .EQ. 1) NO = NO + 1 
100   CONTINUE 
C----Check for C-C(=O)-N and N-C(=O)-N                          11/15/03
      IF (NC .EQ. 1 .AND. NN .EQ. 1 .AND.     ! atom numbers for         12/3/03
     &    NO .EQ. 1) GO TO 120                ! amide, imide, hydantoin  12/3/03
      IF (NC .EQ. 0 .AND. NN .EQ. 2 .AND.     ! atom numbers     12/3/03 
     &    NO .EQ. 1) GO TO 120                ! for hydantoin    12/3/03 
      RETURN                                                   ! 12/3/03
120   DO I=1,NN                                                ! 12/3/03 
C----Call IQ_HN_IMIDE to test if the N is an imide or amide N   11/15/03 
      L = IQ_HN_IMIDE (NUM_N(I),NUM_N(I))                     ! 11/15/03 
      IF (L .EQ. 40) THEN         ! imide N param = 40          11/15/03 
         IQ_OC_IMIDE = 41         ! imide O param = 41          11/15/03
         RETURN                                               ! 11/15/03
      ENDIF                                                   ! 11/15/03
      IF (L .EQ. 62) THEN      ! amide N param = 62           ! 11/15/03
         IQ_OC_IMIDE = 63      ! amide O param = 63           ! 11/15/03
         RETURN                                               ! 11/15/03 
      ENDIF                                                   ! 11/15/03 
      IF (L. EQ. 66) THEN      ! hydantoin N = 66               11/15/03
         IQ_OC_IMIDE = 65      ! hydantoin O = 65               11/15/03
         RETURN 
      ENDIF                                                   ! 11/15/03
      ENDDO                                                   ! 11/15/03
      END    
C 
C----------------------------------------------------------------
C----Function to determine if a single-linked oxygen is an 
C     N-oxide furoxan O
C
      FUNCTION IQ_FUROXAN_O1(I1)    ! I1 is a single-linked O
C
      CHARACTER (LEN=5) :: NAME
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_FUROXAN_O1 = 0     ! no
C
C----Must be linked to a N, which is linked to a C and O & the
C     O must be linked to another N
C                 C---C
C                 "   "
C          20 --> N   N <-- 21
C                  \ / \  
C           22 -->  O   O <-- 23
C
      N = ICON(1,I1)                         ! is the connected
         IF (NAME(N)(1:1) .NE. 'N') RETURN   ! atom a N with
         IF (NCON(N) .NE. 3) RETURN          ! NCON = 3?
      NO = 0           
      NC = 0
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,N)
         IF (NAME(I2)(1:1) .EQ. 'O') NO = NO + 1
            IF (I2 .NE. I1) IO = I2     ! IO is the "other" O 
         IF (NAME(I2)(1:1) .EQ. 'C') NC = NC + 1
100   CONTINUE
      IF (NC .NE. 1 .OR. NO .NE. 2) RETURN 
C----Test the other O...must have 2 N links
      IF (NCON(IO) .NE. 2) RETURN
      DO 120 I=1,2
         IF (NAME(ICON(I,IO))(1:1) .NE. 'N') RETURN
120   CONTINUE
      IQ_FUROXAN_O1 = 23     ! param = 23 for N-oxide O 
      RETURN
C
C----Entry point for the ring O in a furoxan
      ENTRY IQ_FUROXAN_O2(I1)   ! I1 is a double linked O
C
      IQ_FUROXAN_O2 = 0 
C
      NN = 0       ! count of # of N's linked to O
      NLINKED = 0  ! total number of things linked to the N's
      DO 200 I=1,2                
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .NE. 'N') GO TO 200
            NN = NN + 1
            NLINKED = NLINKED + NCON(I2)
            IF (NCON(I2) .EQ. 2) I_2 = I2   ! the 2-linked ring N
            IF (NCON(I2) .EQ. 3) I_3 = I2   ! the 3-linked ring N
200   CONTINUE
      IF (NN .NE. 2 .OR. NLINKED .NE .5) RETURN
      NO = 0
      NC = 0
      DO 300 I=1,2                         ! check the 2-linked N
         IF (NAME(ICON(I,I_2))(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(ICON(I,I_2))(1:1) .EQ. 'O') NO = NO + 1      
300   CONTINUE
      IF (NC. NE. 1 .OR. NO .NE. 1) RETURN
      NO = 0
      NC = 0
      DO 400 I=1,3                         ! check the 3-linked N
         IF (NAME(ICON(I,I_3))(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(ICON(I,I_3))(1:1) .EQ. 'O') NO = NO + 1      
400   CONTINUE
      IF (NC. NE. 1 .OR. NO .NE. 2) RETURN
      IQ_FUROXAN_O2 = 22     ! param = 22 for in-ring O
      RETURN
C
C----Entry point for the 2-linked ring N in a furoxan
C
      ENTRY IQ_FUROXAN_N2(I1)   ! I1 is a double linked N
C
      IQ_FUROXAN_N2 = 0
C
      NO = 0            ! count number of C's & O's
      NC = 0
      DO 500 I=1,2                ! check the 2-linked N
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'C') THEN
             NC = NC + 1
             I_C = I2
         ENDIF 
         IF (NAME(I2)(1:1) .EQ. 'O') THEN
            NO = NO + 1
            I_O = I2
         ENDIF
500   CONTINUE
      IF (NO .NE. 1 .OR. NC .NE. 1) RETURN
      IF (NCON(I_O) .NE. 2) RETURN  ! connected O should be 2-linked
      IF (NCON(I_C) .NE. 3) RETURN  ! connected C should be 3-linked
      IQ_FUROXAN_N2 = 20      ! param = 20
      RETURN
C
C----Entry point for the 3-linked ring N in a furoxan
      ENTRY IQ_FUROXAN_N3(I1)   ! I1 is a triple linked N
C
      IQ_FUROXAN_N3 = 0
C
      NO = 0            ! count number of C's & O's
      NC = 0
      NCON_O = 0
      DO 600 I=1,3                ! check the 3-linked N
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'C') THEN
             NC = NC + 1
             I_C = I2
         ENDIF
         IF (NAME(I2)(1:1) .EQ. 'O') THEN
            NO = NO + 1
            NCON_O = NCON_O + NCON(I2)
         ENDIF
600   CONTINUE
      IF (NO .NE. 2 .OR. NC .NE. 1) RETURN
      IF (NCON(I_C) .NE. 3) RETURN  ! connected C should be 3-linked
      IF (NCON_O .NE. 3) RETURN     ! connected O's must have total
                                    !  of 3 connections 
      IQ_FUROXAN_N3 = 21      ! param = 21
      RETURN
      END
C
C----------------------------------------------------------------
C----Function to determine if a 2-bond O is in nitrate ester
C
      FUNCTION IQ_OXYGEN (IT)    ! IT is the 2-linked O
C    
C                   O (#52)
C                  /
C-----      C--O--N (#51)
C           (#50)  \
C                   O
C
      CHARACTER (LEN=5) :: NAME
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_OXYGEN = 50     ! nitrate -O-
C
      NC = 0
      NN3 = 0      ! # of 3 bond N's
      DO 100 I=1,2
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1 
         IF (NAME(I1)(1:1) .EQ. 'N' .AND. NCON(I1) .EQ. 3) THEN
            NN3 = NN3 + 1
            N_NAME = I1
         ENDIF    
100   CONTINUE
      IF (NC .EQ. 1 .AND. NN3 .EQ. 1) THEN 
         NO2 = 0
         NO1 = 0
         DO 110 I=1,3       ! look at 3-linked N
            I1 = ICON(I,N_NAME)
            IF (NAME(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 2) NO2 = NO2 + 1
            IF (NAME(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 1) NO1 = NO1 + 1
110      CONTINUE
         IF (NO2 .EQ. 1 .AND. NO1 .EQ. 2) GO TO 115 
           ELSE
              IQ_OXYGEN = 0
115        ENDIF
      RETURN
      END FUNCTION IQ_OXYGEN      
C
C---------------------------------------------------------------- 
C----Function to determine if we have an alcohol O-H.                3/25/03 
C 
      FUNCTION IQ_ALCOHOL(ITT)    ! ITT is a 2-linked O or H         3/25/03
C                                   or H bonded to 2-linked O        3/25/03
C 
C              C[sp3]---O---H                                        3/25/03 
C                     #58  59                                        3/25/03
C
      CHARACTER (LEN=5) :: NAME
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_ALCOHOL = 0                                               ! 3/25/03 
C----Atom could be O or H; look at the O                           ! 3/25/03
      IF (NAME(ITT)(1:1) .EQ. 'O') THEN                            ! 3/25/03 
         INCREMENT = 0                                             ! 3/25/03  
         IT = ITT                                                  ! 3/25/03
      ELSE                                                         ! 3/25/03
         IT = ICON(1,ITT)     ! pick up attached O                   3/25/03
         IF (NAME(IT)(1:1) .EQ. 'O') INCREMENT = 1                 ! 3/25/03
      ENDIF                                                        ! 3/25/03
      NH = 0                                                       ! 3/25/03
      NC4 = 0                                                      ! 3/25/03
      DO I=1,2                                                     ! 3/25/03
         IT1 = ICON(I,IT)                                          ! 3/25/03
         IF (NAME(IT1)(1:1) .EQ. 'H') NH = NH + 1                  ! 3/25/03
         IF (NAME(IT1)(1:1) .EQ. 'C' .AND. NCON(IT1) .EQ. 4)       ! 3/25/03
     2          NC4 = NC4 + 1                                      ! 3/25/03
      ENDDO                                                        ! 3/25/03
      IF (NH .EQ. 1 .AND. NC4 .EQ. 1) THEN                         ! 3/25/03
         IQ_ALCOHOL = 58 + INCREMENT   ! O of alcohol, param = 58    3/25/03
      ENDIF                            ! or alcohol H, param = 59    3/25/03
      RETURN                                                       ! 3/25/03
      END FUNCTION IQ_ALCOHOL                                      ! 3/25/03
C     
C----------------------------------------------------------------
C----Function to determine if a 3-linked N is a nitrate ester N
C      a nitramine nitro N or a nitrate (NO3-) N.            7/20/06
C
      FUNCTION IQ_NITRATE_N (IT)    ! IT is the 3-linked N
C    
C
C                   O (#52)         O (#55)   (#84)    7/20/06
C                  /               /            |      7/20/06
C-----      C--O--N (#51)       N--N (#54)    N(O3)-   7/20/06
C           (#50)  \                \         |        7/20/06
C                   O                O      (#83)      7/20/06
C
      CHARACTER (LEN=5) :: NAME
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_NITRATE_N = 51     ! nitrate ester N
C
      NO1 = 0      ! # of 1 bond O's
      NO2 = 0      ! # of 2 bond O's
      NN3 = 0      ! # of 3 bond N's                         ! 2/24/03
      DO 100 I=1,3
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 1) NO1 = NO1 + 1  
         IF (NAME(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 2) NO2 = NO2 + 1  
         IF (NAME(I1)(1:1) .EQ. 'N' .AND. NCON(I1) .EQ. 3) NN3 = NN3 + 1   ! 2/24/03
100   CONTINUE
      if (no1 .eq. 3) then                                     ! 7/20/06
         iq_nitrate_n = 83     ! param = 83 for N of NO3-      ! 7/20/06
         return                                                ! 7/20/06
      endif                                                    ! 7/20/06
      IF (NO1 .EQ. 2 .AND. NO2 .EQ. 1) RETURN
      IF (NO1 .EQ. 2 .AND. NN3 .EQ. 1) THEN
         IQ_NITRATE_N = 54             ! param = 54 for n of N-nO2    2/24/03
      ELSE                                                          ! 2/24/03
         IQ_NITRATE_N = 0
      ENDIF                                                         ! 2/24/03
      RETURN
      END FUNCTION IQ_NITRATE_N      
C----------------------------------------------------------------------
      SUBROUTINE OPEN_G03     ! open the G03 charge file on unit # 17
C
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      CHARACTER  CHARGE_FILE1*30, G03_LINE*21,
     x           MOPAC_LINE*31, MOPAC1*6
      DATA CHARGE_FILE1/'                              '/
C
          PRINT 9
9         FORMAT ('   Name of the G03 log file', 
     x            ' that contains the charge information;'/
     x            '     this file can have a single line header of',    
     x            ' "#G03" just ahead of the'/
     x            '     first atom charge line',
     x            ' [631gstar_chelpg.log]:',$)    
          READ 27, CHARGE_FILE1
27        FORMAT (A30)
          IF (CHARGE_FILE1 .EQ. '        ') THEN
             CHARGE_FILE1 = '631gstar_chelpg.log'
          ENDIF
          OPEN (UNIT=17, FILE=CHARGE_FILE1, STATUS='OLD')
C----Position the charge log file just before the 1st atom
14          READ (17,13,END=18) G03_LINE
13          FORMAT (A21)
            IF (G03_LINE(1:4) .EQ. '#G03') GO TO 10
            IF (G03_LINE .NE. ' Charges from ESP fit') GO TO 14
            READ (17,*)       ! got it....skip the 
            READ (17,*)       ! next 2 lines
            GO TO 10
18        PRINT 19
19        FORMAT (' Cannot find appropriate header line in ',
     x            'the G03 file...abort') 
          STOP                                 ! 11-4-94
10        RETURN
          END
C
C----------------------------------------------------------------
C
      SUBROUTINE OPEN_MNDO     ! open the MNDO charge file on unit # 18
C
      LOGICAL CONNECT_FLAG
      COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G03_CHARGE,
     X              XMNDO_CHARGE
      CHARACTER  CHARGE_FILE2*30, G03_LINE*21,
     x           MOPAC_LINE*31, MOPAC1*6
          PRINT 9
9         FORMAT ('   Name of the MOPAC MOPOUT file',
     x            ' that contains the charge information;'/
     x            '     this file can have a single line header of',
     x            ' "#MOPAC" just ahead of the'/
     x            '     first atom charge line:',$)
          READ 27, CHARGE_FILE2
27        FORMAT (A30)
          OPEN (UNIT=18, FILE=CHARGE_FILE2, STATUS='OLD')
15          READ (18,21,END=25) MOPAC1, MOPAC_LINE
21          FORMAT (A6,6x,A31)
            IF (MOPAC1 .EQ. '#MOPAC') GO TO 10
            IF (MOPAC_LINE .NE. 'ELECTROSTATIC POTENTIAL CHARGES')
     x         GO TO 15
            READ (18,*)    ! skip 2 lines,
            READ (18,*)    ! then ready to read 1st charge line
            GO TO 10
25        PRINT 26
26        FORMAT (' Cannot find appropriate header line in',
     x            ' MOPAC file...abort')
          STOP
10        RETURN
          END
C------------------------------------------------------------
C----Function to determine if a C atom is a cubane C
C
      FUNCTION IQ_CUBANE(I1)   
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_CUBANE = 0   ! 0 indicates its not cubane-related
      IF (NCON(I1) .NE. 4) RETURN   ! C must be linked to 4 atom
C
C----How many C's bonded to I1...must be at least 3
      NC = 0
      DO 20 I=1,4
          I2 = ICON(I,I1)
          IF (NAME(I2)(1:1) .EQ. 'C') NC = NC + 1
20    CONTINUE
      IF (NC .LT. 3) RETURN   ! not at least 3 C's...return
C----Look for 3 C-C-C bond angles .le. 93 degs
      NC = 0
      DO 200 I=1,3
         DO 200 J=I+1,4
            I2 = ICON(I,I1)
            I3 = ICON(J,I1)
            IF (.NOT. (NAME(I2)(1:1) .EQ. 'C' .AND.        ! must both
     x                 NAME(I3)(1:1) .EQ. 'C')) GO TO 200  ! be C's
            DO 40 L=1,2           ! calc the I2 - I1 & I3 - I1
               XL(L) = 0.0        ! vectors & their lengths
               DO 38 K=1,3
                  VEC(K,L) = XYZ(K,IND(L)) - XYZ(K,I1)
                  XL(L) = XL(L) + VEC(K,L)**2
38             CONTINUE
               XL(L) = SQRT(XL(L))
40          CONTINUE
            D = 0.0
            DO 50 K=1,3
               D = D + VEC(K,1)*VEC(K,2)
50          CONTINUE
            COSINE = D/(XL(1)*XL(2))              
            IF (COSINE .GE. -0.052336) NC = NC + 1  ! cos(93) = -0.052336
200   CONTINUE
      IF (NC .NE. 3) RETURN    ! not 3 angles .le. 93 deg
         IQ_CUBANE = 1        ! 3 angles .le. 93 deg, assume its 
      RETURN                  ! a cubane C atom
      END 
C
C----------------------------------------------------------------
C----Function to determine if an O is a nitro O, if the nitro        2/24/03
C     is cubane-linked, a nitrate esters or nitramine                2/24/03
C
      FUNCTION IQ_NITRO_OX(IT)    ! IT is the N linked to the O
C
C                   O (#52)                     O (#55)              2/24/03
C                  /                           /                     2/24/03
C-----      C--O--N (#51)                  N--N (#54)                2/24/03
C           (#50)  \                           \                     2/24/03
C                   O                           O                    2/24/03
C
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      IQ_NITRO_OX = 6     ! normal nitro group O
C
      NO = 0
      NC = 0
      NO1 = 0      ! # of 1 bond O's
      NO2 = 0      ! # of 2 bond O's
      NN3 = 0      ! # of 3 bond N's
      DO 100 I=1,3
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'O') NO = NO + 1
         IF (NAME(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 1) NO1 = NO1 + 1  
         IF (NAME(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 2) NO2 = NO2 + 1  
         IF (NAME(I1)(1:1) .EQ. 'N' .AND. NCON(I1) .EQ. 3) NN3 = NN3 + 1
         IF (NAME(I1)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            IC = I1  ! IC = C attached to the N
         ENDIF
100   CONTINUE
      IF (NC .EQ. 0 .AND. NO1 .EQ. 2 .AND. NO2 .EQ. 1) THEN
         IQ_NITRO_OX = 52     ! param = 52; 1 bond nitrate ester O
         RETURN
      ENDIF
      IF (NC .EQ. 1 .AND. NO1 .EQ. 2) GO TO 200
      IF (NN3 .EQ. 1 .AND. NO1 .EQ. 2) THEN
         IQ_NITRO_OX = 55        ! param = 55 for O of N-NO2    2/24/03
      ELSE
         IQ_NITRO_OX = 0
      ENDIF
      RETURN
200   IF (IQ_CUBANE(IC) .EQ. 1) IQ_NITRO_OX = 49 ! cubane-linked nitro
      RETURN                                     ! group; O param = 49  11-12-02
      END FUNCTION IQ_NITRO_OX
C
C---------------------------------------------------------------------
C----Function to identify the 3 N's in an azide                     ! 3/5/03
C
      FUNCTION IQ_AZIDE_N (IT)    ! IT is a 1 or 2-linked N         ! 3/5/03
C
C         C---N===N===N                                             ! 3/5/03
C            32  56  57                                             ! 3/5/03
C
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
C
      NC = 0    ! number of attached C's                            ! 3/5/03
      NN1 = 0   ! number of 1-linked N's                            ! 3/5/03
      NN2 = 0   ! number of 2-linked N's                            ! 3/5/03
C
      DO I=1,NCON(IT)                                               ! 3/5/03
         IT1 = ICON(I,IT)                                           ! 3/5/03
         IF (NAME(IT1)(1:1) .EQ. 'C') NC = NC + 1                   ! 3/5/03
         IF (NAME(IT1)(1:1) .EQ. 'N' .AND. NCON(IT1) .EQ. 1)        ! 3/5/03
     2        NN1 = NN1 + 1                                         ! 3/5/03
         IF (NAME(IT1)(1:1) .EQ. 'N' .AND. NCON(IT1) .EQ. 2)        ! 3/5/03
     2        NN2 = NN2 + 1                                         ! 3/5/03
      ENDDO                                                         ! 3/5/03
C----First look for the terminal N atom                             ! 3/5/03
      IF (NCON(IT) .EQ. 1 .AND. NN2 .EQ. 1) THEN                    ! 3/5/03
         IT2 = ICON(1,IT)   ! check out the middle N to be certain  ! 3/5/03
         MN1 = 0   ! number of 1-linked N's                         ! 3/5/03
         MN2 = 0   ! number of 2-linked N's                         ! 3/5/03 
         DO I=1,2                                                   ! 3/5/03
            IT3 = ICON(I,IT2)                                       ! 3/5/03
            IF (NAME(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 1)     ! 3/5/03
     2           MN1 = MN1 + 1                                      ! 3/5/03
            IF (NAME(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 2)     ! 3/5/03
     2           MN2 = MN2 + 1                                      ! 3/5/03
         ENDDO                                                      ! 3/5/03
         IF (MN1 .EQ. 1 .AND. MN2 .EQ. 1) THEN                      ! 3/5/03
            IQ_AZIDE_N = 57     ! terminal N, param = 57            ! 3/5/03
            RETURN                                                  ! 3/5/03
         ENDIF                                                      ! 3/5/03
      ENDIF                                                         ! 3/5/03
C----Second look for the central N atom                             ! 3/5/03
      IF (NCON(IT) .EQ. 2 .AND. NN1 .EQ. 1 .AND. NN2 .EQ. 1) THEN   ! 3/5/03
          IQ_AZIDE_N = 56     ! central N, param = 56               ! 3/5/03
          RETURN                                                    ! 3/5/03
      ENDIF                                                         ! 3/5/03
C----Third look for the C-linked N atom                             ! 3/5/03
      IF (NCON(IT) .EQ. 2 .AND. NC .EQ. 1 .AND. NN2 .EQ. 1) THEN    ! 3/5/03
         DO 300 I=1,2                                               ! 3/5/03
            IT2 = ICON(I,IT)                                        ! 3/5/03
            IF (NAME(IT2)(1:1) .EQ. 'C') GO TO 300                  ! 3/5/03
            NN1 = 0    ! IT2 is now the central N                   ! 3/5/03
            NN2 = 0                                                 ! 3/5/03
            DO J=1,2                                                ! 3/5/03
               IT3 = ICON(J,IT2)                                    ! 3/5/03
               IF (NAME(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 1)  ! 3/5/03
     2              NN1 = NN1 + 1                                   ! 3/5/03
               IF (NAME(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 2)  ! 3/5/03
     2              NN2 = NN2 + 1                                   ! 3/5/03
            ENDDO                                                   ! 3/5/03
            IF (NN2 .EQ. 1 .AND. NN1 .EQ. 1) THEN                   ! 3/5/03
               IQ_AZIDE_N = 32    ! C-linked N, param = 32          ! 3/5/03
               RETURN                                               ! 3/5/03
            ENDIF                                                   ! 3/5/03
300      CONTINUE                                                   ! 3/5/03
      ENDIF                                                         ! 3/5/03
      IQ_AZIDE_N = 0                                                ! 3/5/03
      END FUNCTION IQ_AZIDE_N                                       ! 3/5/03 
C
C---------------------------------------------------------------------
      SUBROUTINE MOLFIT2WMIN
C
      CHARACTER AID*3, IN_FILE*30, LINE*89
      DIMENSION XYZ(3)
C
      PRINT 9
9     FORMAT (' Name of the molecular fit output file',
     x        ' [molecular-fit.output-file]:',$)
      READ 27, IN_FILE
27    FORMAT (A30)
      IF (IN_FILE(1:1) .NE. ' ') GO TO 29
         IN_FILE(1:25) = 'molecular-fit.output-file'
29    OPEN (UNIT=18, FILE=IN_FILE, STATUS='OLD')
      OPEN (UNIT=19, FILE='wmin.xyz', STATUS='UNKNOWN')
C----Search for 'FRACT. COORDS. OF 2ND MOLECULE'
50    READ (18,55) LINE
55    FORMAT (A89)
      IF (LINE(60:89) .NE. 'FRACT. COORDS. OF 2ND MOLECULE')
     X      GO TO 50
      READ (18,*)          !  skip 1 line
      N = 0
60    READ (18,55,END=100) LINE
      IF (LINE(1:20) .EQ. '                    ') GO TO 100
      READ (LINE,65) AID, XYZ
65    FORMAT (52X,A3,1x,3F10.5)
      N = N + 1
      WRITE (19,70) AID, N, XYZ
70    FORMAT (A3,1X,I5,18X,3F9.5)
      GO TO 60
100   RETURN
      END
C--------------------------------------------------------------
      SUBROUTINE WMIN2CHEMX
C
C----Convert WMIN output file coords to CHEMX, CRYSTAL, CCDC,       3/9/04 
C     SHELX res, CIF iand new wmin.save formats.                  12/19/06 
C
      LOGICAL OUT_CHEMX, OUT_CRYSTAL, L_MATRIX, OUT_CCDC, OUT_RES,! 3/9/04
     2        OUT_CIF, T_MATRIX, OUT_new_wmin                     ! 12/19/06 
      PARAMETER (NUMBER_ELEMENTS = 6)
      CHARACTER AID(200)*4, IN_FILE*30, OUT_FILE*30, LINE*50,
     2          OUT_FILE2*30, AT_SYMBOL(NUMBER_ELEMENTS)*1, BLANK*30,
     3          IF_MATRIX*1, OUT_FILE3*30, OUT_FILE4*30,          ! 3/2/04
     4          E_NAME(10)*1, OUT_FILE5*30, SP_GRP*10, CPD_ID*30, ! 3/9/04
     5          OUT_FILE6*30                                      ! 12/19/06
      DIMENSION XYZ(200,3), CELL(6), NUMBER_AT(NUMBER_ELEMENTS),
     2          T2_MATRIX(3,3), CELL_NEW(6), TEMP(3,3), LL(3), MM(3),
     3          XYZ_NEW(3), NUMBER_E(10), XYZ_TRANS(3), cosines(3),
     4          xyzsum(3)                                         ! 12/19/06 
      DATA AT_SYMBOL /'H', 'B', 'C', 'N', 'O', 'F'/
      DATA NUMBER_AT / 1,   5,   6,   7,   8,   9/
      DATA BLANK /'                              '/
      DATA OUT_CHEMX, OUT_CRYSTAL, OUT_CCDC, OUT_RES, OUT_CIF,    ! 3/9/04
     2     OUT_new_wmin /6*.false./                               ! 12/19/06
C
      NATOMS = 0
      L_MATRIX = .false.
      T_MATRIX = .false.                                          ! 3/14/04
C
      PRINT 9
9     FORMAT (' Name of the WMIN fractional coordinate ',
     x        'input file:',$)
      READ 27, IN_FILE
27    FORMAT (A30)
      OPEN (UNIT=18, FILE=IN_FILE, STATUS='OLD')
C----Chemx output file
      PRINT 10
10    FORMAT (' Name of the CHEMX output file [ ] = none:',$)
      READ 27, OUT_FILE
      IF (OUT_FILE .NE. BLANK) THEN
         OPEN (UNIT=19, FILE=OUT_FILE, STATUS='UNKNOWN')
         OUT_CHEMX = .true.
      ENDIF
C----CRYSTAL output file
      PRINT 101
101   FORMAT (' Name of the CRYSTAL output file [ ] = none:',$)
      READ 27, OUT_FILE2
      IF (OUT_FILE2 .NE. BLANK) THEN
         OPEN (UNIT=20, FILE=OUT_FILE2, STATUS='UNKNOWN')
         OUT_CRYSTAL = .true.
      ENDIF
C----SHELX res output file                                          3/2/04
      PRINT 105                                                   ! 3/2/04
105   FORMAT (' Name of the SHELX res output file [ ] = none:',$) ! 3/2/04
      READ 27, OUT_FILE3                                          ! 3/2/04
      IF (OUT_FILE3 .NE. BLANK) THEN                              ! 3/2/04
         OPEN (UNIT=21, FILE=OUT_FILE3, STATUS='UNKNOWN')         ! 3/2/04
         OUT_RES = .true.                                         ! 3/2/04
      ENDIF                                                       ! 3/2/04
C----CIF output file                                                3/9/04
      PRINT 106                                                   ! 3/9/04
106   FORMAT (' Name of the CIF output file [ ] = none:',$)       ! 3/9/04
      READ 27, OUT_FILE5                                          ! 3/9/04
      IF (OUT_FILE5 .NE. BLANK) THEN                              ! 3/9/04
         OPEN (UNIT=28, FILE=OUT_FILE5, STATUS='UNKNOWN')         ! 3/9/04
         OUT_CIF = .true.                                         ! 3/9/04
      ENDIF                                                       ! 3/9/04
C----CCDC output file
      PRINT 199
199   FORMAT (' Name of the CCDC output file [ ] = none:',$)
      READ 27, OUT_FILE4
      IF (OUT_FILE4 .NE. BLANK) THEN
         OPEN (UNIT=23, FILE=OUT_FILE4, STATUS='UNKNOWN')
         OUT_CCDC = .true.
      ENDIF
C----new wmin.save format file                                 ! 12/19/06
      PRINT 201                                                ! 12/19/06
201   FORMAT (' Name of the new wmin.save format output file', ! 12/19/06
     2        ' [ ] = none:',$)                                ! 12/19/06
      READ 27, OUT_FILE6                                       ! 12/19/06
      IF (OUT_FILE6 .NE. BLANK) THEN                           ! 12/19/06
         OPEN (UNIT=29, FILE=OUT_FILE6, STATUS='UNKNOWN')      ! 12/19/06
         out_new_wmin = .true.                                 ! 12/19/06 
         do l=1,3                                              ! 12/19/06
            xyzsum(l) = 0.0                                    ! 12/19/06
         enddo                                                 ! 12/19/06
      ENDIF                                                    ! 12/19/06
C----Should a translation vector be read?                        3/14/04
        PRINT 980                                              ! 3/14/04
980     FORMAT (' Read a fractional coordinate translation vector,' ! 3/14/04
     x          ,' Y/N, [ ] = N:',$)                            ! 3/14/04
        READ 1002, IF_MATRIX
1002    FORMAT (A)
        IF (IF_MATRIX .EQ. 'Y' .OR. IF_MATRIX .EQ. 'y') THEN
           T_MATRIX = .true.                                   ! 3/14/04
           PRINT 1009                                          ! 3/14/04
1009       FORMAT (' Supply matrix as del X, del Y, del Z:',$) ! 3/14/04
           READ (5,*) (XYZ_TRANS(I),I=1,3)                     ! 3/14/04
        ENDIF
C----Should a transformation matrix be read?
        PRINT 1001
1001    FORMAT (' Read a transformation matrix (eg. T2 from ',
     x          'NIST*LATTICE), Y/N, [ ] = N:',$)
        READ 1002, IF_MATRIX
        IF (IF_MATRIX .EQ. 'Y' .OR. IF_MATRIX .EQ. 'y') THEN
           L_MATRIX = .true.
           PRINT 1013
1013       FORMAT (' Supply matrix as A11, A12, A13....C31..C33:',$)
           READ (5,*) ((T2_MATRIX(I,J),J=1,3),I=1,3)
        ENDIF
C----Read cell params from line 1 of WMIN input (save) file 
      READ (18,*) CELL
C----Apply cell transformation matrix?
      IF (L_MATRIX) THEN
         DO 1005 I=1,3           ! 1st get new cell lengths
            DO J=1,3
               TEMP(I,J) = T2_MATRIX(I,J)*CELL(J)
            ENDDO
            SUM = 0.0         ! vector dot product to get cell length
            DO 1015 K=1,3
               SUM = SUM + TEMP(I,K)*TEMP(I,K)
1015        CONTINUE
            SUM = SUM + 2.0*(TEMP(I,1)*TEMP(I,2)*CELL(6) +
     2                       TEMP(I,1)*TEMP(I,3)*CELL(5) +
     3                       TEMP(I,2)*TEMP(I,3)*CELL(4))
            CELL_NEW(I) = SQRT(SUM)     ! new unit cell length
1005     CONTINUE
C----Now the new unit cell angles
         DO 1029 I=1,2
            DO 1030 M=I+1,3
               SUM = 0.0
               DO 1035 L=1,3
               DO 1035 J=1,3
                  IF (L .EQ. J) THEN
                      SUM = SUM + TEMP(I,L)*TEMP(M,J)
                  ELSE
                      SUM = SUM + TEMP(I,L)*TEMP(M,J)*CELL(9-L-J)
                  ENDIF
1035           CONTINUE
               CELL_NEW(9-I-M) = SUM/(CELL_NEW(I)*CELL_NEW(M))   ! new cosine  
1030        CONTINUE
1029     CONTINUE
         PRINT 1050, CELL, CELL_NEW
1050     FORMAT (' Original cell     =',6f10.5/
     x           ' Transformed cell  =',6f10.5)
      DO I=1,6                    ! place new cell info in cell array
         CELL(I) = CELL_NEW(I)
      ENDDO
C----Get inverse of T2_MATRIX and then its transpose
         CALL MINV (T2_MATRIX, 3, D, LL, MM)
         IF (D .GT. 0.0) GO TO 1070
            PRINT 1077
1077        FORMAT (' Transform matrix singular...stop')
            STOP
1070     DO 1080 I=1,3               ! take transpose & store in TEMP
            DO 1080 J=1,3
               TEMP(I,J) = T2_MATRIX(J,I)
1080     CONTINUE
      ENDIF
      DO I=4,6
         cosines(i-3) = cell(i)                                  ! 12/19/06
         CELL(I) = ACOS(CELL(I))*57.2957795  ! convert to angles
      ENDDO
      IF (OUT_CHEMX) WRITE (19,77) CELL    ! lines 1 & 2 of CHEMX file
77    FORMAT (38X,3F8.4/21X,3F8.4)
      IF (OUT_CRYSTAL) WRITE (20,177) CELL ! cell line of CRYSTAL file
177   FORMAT (6F10.4)
      IF (OUT_CCDC) WRITE (23,179) CELL    ! cell line for CCDC file
179   FORMAT ('CELL',6F10.4)
      IF (OUT_new_wmin) WRITE (29,279) (CELL(i),i=1,3), cosines ! cell line for 
279   FORMAT (6F9.5)                ! new wmin.save format file    ! 12/19/06
      IF (OUT_RES) THEN                    ! first lines for SHELX res file  3/2/04
         WRITE (21,181) CELL                                        ! 3/2/04
181      FORMAT ('TITL UNKNWN'/'CELL 0.71073', 3F7.3,3F7.2/         ! 3/2/04
     2           'ZERR 4 0.00   0.00   0.00   0.00   0.00   0.00'/  ! 3/2/04
     3           'LATT 1 <- ???'/'SYMM <- ???')                     ! 3/2/04
         PRINT 184                                                  ! 3/2/04
184      FORMAT ('Number of molecules in unit cell: ',$)            ! 3/2/04
         READ (*,*) N_MOLECULES                                     ! 3/2/04 
      ENDIF                                                         ! 3/2/04
      IF (OUT_CIF) THEN                                             ! 3/9/04
         PRINT 510                                                  ! 3/9/04
510      FORMAT (' Cpd ID to follow data_, eg ''ammon.III''',       ! 3/9/04
     2           ' would give ''data_ammon.III'':',$)               ! 3/9/04
         READ (*,1002) CPD_ID                                       ! 3/9/04
         WRITE (28,501) CPD_ID                                      ! 3/9/04
501      FORMAT ('data_',A/'# comments'/)   ! line 1 & comments       3/9/04
         WRITE (28,502) CELL    ! 6 lines with cell params          ! 3/9/04 
502      FORMAT ('_cell_length_a ',F7.4/                            ! 3/9/04
     2           '_cell_length_b ',F7.4/                            ! 3/9/04
     3           '_cell_length_c ',F7.4/                            ! 3/9/04
     4           '_cell_angle_alpha ',F7.3/                         ! 3/9/04
     5           '_cell_angle_beta ',F7.3/                          ! 3/9/04
     6           '_cell_angle_gamma ',F7.3)                         ! 3/9/04
         PRINT 506                                                  ! 3/9/04
506      FORMAT (' Number of formula units per unit cell:',$)       ! 3/9/04
         READ (*,*) NZ                                              ! 3/9/04
         WRITE (28,507) NZ                                          ! 3/9/04
507      FORMAT ('_cell_formula_units_Z',I2)                        ! 3/9/04
         PRINT 503                                                  ! 3/9/04
503      FORMAT (' Space group with defining ''s, eg ''P 21/c'':',$) ! 3/9/04
         READ (*,1002) SP_GRP                                       ! 3/9/04
         WRITE (28,504) SP_GRP                                      ! 3/9/04
504      FORMAT ('_symmetry_space_group_name_H-M ',A/'#'/'loop_'/   ! 3/9/04           
     2           '_atom_site_label'/'_atom_site_fract_x'/           ! 3/9/04
     3           '_atom_site_fract_y'/'_atom_site_fract_z')         ! 3/9/04
      ENDIF                                                         ! 3/9/04
C----Read atoms from WMIN file
100   NATOMS = NATOMS + 1
      READ (18,120) AID(NATOMS), (XYZ(NATOMS,I),I=1,3)  ! read WMIN
120   FORMAT (A4,23X,3F9.5)                             ! save file
      IF (AID(NATOMS) .NE. 'XTRA')  GO TO 100           !looking for XTRA
      NATOMS = NATOMS - 1
C----Finished with atom input, write CHEMX
      IF (OUT_CHEMX) WRITE (19,130) NATOMS
130   FORMAT (I4,'   1'/)                         ! 3rd & 4th CHEMX lines
      WRITE (20,131) NATOMS                       ! atom # line for CRYSTAL
131   FORMAT (I2)
      DO 200 I=1,NATOMS
C--------Apply the fract coord translation vector?                 3/14/04
         IF (T_MATRIX) THEN
             DO M=1,3                                          ! 3/14/04
                XYZ(I,M) = XYZ_TRANS(M) + XYZ(I,M)             ! 3/14/04
             ENDDO                                             ! 3/14/04
         ENDIF                                                   ! 3/14/04
C--------Apply the T2 transformation matrix?                       3/14/04
         IF (L_MATRIX) THEN
            DO 1085 L=1,3
               XYZ_NEW(L) = 0.0
               DO 1085 M=1,3
                  XYZ_NEW(L) = XYZ_NEW(L) + TEMP(L,M)*XYZ(I,M)
1085        CONTINUE
         DO L=1,3
            XYZ(I,L) = XYZ_NEW(L)
         ENDDO
         ENDIF
         IF (OUT_CHEMX) WRITE (19,300) I, AID(I), (XYZ(I,J),J=1,3) ! CHEMX atoms
300      FORMAT (I4,1X,A4,1X,3F10.6)
         IF (OUT_new_wmin) then                                      ! 12/19/06
            WRITE (29,305) AID(I), i, (XYZ(I,J),J=1,3) ! new wmin.save 
305         FORMAT (A4,i5,18x,3F9.5)                                 ! 12/29/06
            do l=1,3                                    ! collect sums   12/19/06
               xyzsum(l) = xyzsum(l) + xyz(i,l)           ! for XTRA line  12/19/06
            enddo                                         ! 12/19/06
         endif                                            ! 12/19/06
         IF (OUT_CCDC) THEN 
            IF (I .LE. 9) THEN
               WRITE (23,311) AID(I), I, (XYZ(I,J),J=1,3)  ! CCDC atoms
311            FORMAT ('ATOM  ',A1,I1,3X,3F10.5)
            ELSE
               WRITE (23,312) AID(I), I, (XYZ(I,J),J=1,3)  ! CCDC atoms
312            FORMAT ('ATOM  ',A1,I2,2X,3F10.5)
            ENDIF
         ENDIF
         IF (OUT_CIF) THEN            ! CIF format atoms     3/9/04 
            NCS = 1                                        ! 3/9/04
            IF (AID(I)(2:2) .NE. ' ') NCS = 2              ! 3/9/04 
            IF (I .LE. 9) THEN                             ! 3/9/04 
               WRITE (28,321) AID(I)(1:NCS), I, (XYZ(I,J),J=1,3) ! 3/9/04 
321            FORMAT (1X,A,I1,4X,3F10.5)                  ! 3/9/04
            ELSE                                           ! 3/9/04
               WRITE (28,322) AID(I)(1:NCS), I, (XYZ(I,J),J=1,3) ! 3/9/04 
322            FORMAT (1X,A,I2,3X,3F10.5)                  ! 3/9/04
            ENDIF                                          ! 3/9/04
         ENDIF                                             ! 3/9/04
         IF (.not.OUT_CRYSTAL) GO TO 200
         DO 250 J=1,NUMBER_ELEMENTS  ! convert atom symbol to number for CRYSTAL
            IF (AID(I)(1:1) .NE. AT_SYMBOL(J)) GO TO 250
            WRITE (20,245) NUMBER_AT(J), (XYZ(I,K),K=1,3)
245         FORMAT (I2,3X,3F10.6)
            GO TO 200
250      CONTINUE
         PRINT 1111, AID(I)(1:1)
1111     FORMAT (' Cannot identify atomic symbol ',A1,'; abort')
         STOP
200   CONTINUE
      if (out_new_wmin) then                                ! 12/19/06
         do l=1,3                                         ! 12/19/06
            xyzsum(l) = xyzsum(l)/natoms                    ! 12/19/06
         enddo                                              ! 12/19/06
         write (29,307) (natoms+1), xyzsum                  ! 12/19/06
307      FORMAT ('XTRA',i5,18x,3F9.5)                       ! 12/19/06
      endif                                                 ! 12/19/06
      IF (.not. OUT_RES) RETURN                             ! 3/2/04
         N_ELEMENTS = 0                                     ! 3/2/04
         DO 304 I=1,NATOMS                                  ! 3/2/04
            IF (N_ELEMENTS .EQ. 0) GO TO 405                ! 3/2/04
               DO 400 K=1,N_ELEMENTS                        ! 3/2/04
                  IF (AID(I)(1:1) .NE. E_NAME(K)) GO TO 400 ! 3/2/04
                  NUMBER_E(K) = NUMBER_E(K) + 1             ! 3/2/04
                  GO TO 304                                 ! 3/2/04
400            CONTINUE                                     ! 3/2/04
405            N_ELEMENTS = N_ELEMENTS + 1                  ! 3/2/04
               E_NAME(N_ELEMENTS) = AID(I)(1:1)             ! 3/2/04
               NUMBER_E(N_ELEMENTS) = 1                     ! 3/2/04 
304   CONTINUE                                              ! 3/2/04
      WRITE (21,410) (E_NAME(K),K=1,N_ELEMENTS)             ! 3/2/04
410     FORMAT ('SFAC ',10A2)                               ! 3/2/04
      WRITE (21,420) (N_MOLECULES*NUMBER_E(K),K=1,N_ELEMENTS) ! 3/2/04      
420     FORMAT ('UNIT ',10I3)                               ! 3/2/04
      WRITE (21,430)                                        ! 3/2/04
430     FORMAT ('FVAR 1.00')                                ! 3/2/04
      DO 500 I=1,NATOMS                                     ! 3/2/04
         DO 520 K=1,N_ELEMENTS                              ! 3/2/04
            IF (AID(I)(1:1) .NE. E_NAME(K)) GO TO 520       ! 3/2/04
            IF (I .LE. 9) THEN                              ! 3/2/04
               WRITE (21,540) AID(I)(1:1), I, K,            ! 3/2/04
     2                       (XYZ(I,L),L=1,3)               ! 3/2/04
540            FORMAT (A,I1,I5,3F10.6,'  1.000000 0.05000') ! 3/2/04
            ELSE                                            ! 3/2/04
               WRITE (21,541) AID(I)(1:1), I, K,            ! 3/2/04
     2                       (XYZ(I,L),L=1,3)               ! 3/2/04
541            FORMAT (A,I2,I4,3F10.6,'  1.000000 0.05000') ! 3/2/04
            ENDIF                                           ! 3/2/04
520      CONTINUE                                           ! 3/2/04
500   CONTINUE                                              ! 3/2/04
      WRITE (21,550)                                        ! 3/2/04
550     FORMAT ('END')                                      ! 3/2/04
      CLOSE (UNIT=21)                                       ! 3/2/04
      RETURN
      END
C
C------------------------------------------------------------------------
C----Function to determine if a N is one of the 4 N's is the
C     Z-tetraazapentalene moiety 
C
      FUNCTION IQ_N_AZAPENTALENE(IT) ! IT is a N 
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      DIMENSION  IND(2), XL(2), VEC(3,2), INN(2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
C
      IQ_N_AZAPENTALENE = 0   ! 0 means can't identify type
C
      IF (NCON(IT) .LT. 2 .OR. NCON(IT) .GT. 3) RETURN
      GO TO (10,20), (NCON(IT) - 1)
C----2 attachments...check for an end N
10    NN = 0
      NC = 0
      DO 100 I=1,2
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            IN = I1  ! IN = N attached to the N in question
         ENDIF
100   CONTINUE
      IF (NC .NE. 1 .OR. NN .NE. 1) RETURN 
      IF (NCON(IN) .NE. 3) RETURN 
C----Have 2 attachments, 1 C & 1 N and N has 3 attachments
      NN = 0
      NC = 0
      DO 110 I=1,3
         I1 = ICON(I,IN)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') NN = NN + 1
110   CONTINUE
      IF (NC .EQ. 1 .AND. NN .EQ. 2) THEN
         IQ_N_AZAPENTALENE = 19   ! param = 19 for Z-tetraazapentalene N
      ENDIF 
      RETURN
C----3 attachments...check for a middle N
20    NN = 0
      NC = 0
      DO 210 I=1,3
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            INN(NN) = I1      ! number of the attached N
         ENDIF
210   CONTINUE
      IF (NC .EQ. 1 .AND. NN .EQ. 2) THEN     ! attached N's must have
         I1 = MIN0(NCON(INN(1)),NCON(INN(2)))   ! 2 and 3 attachments
         I2 = MAX0(NCON(INN(1)),NCON(INN(2)))   ! respectively
         IF (I1 .EQ. 2 .AND. I2 .EQ. 3) IQ_N_AZAPENTALENE = 19 ! param = 19
      ENDIF                                   ! for Z-tetraazapentalene N
      RETURN
      END
C
C-------------------------------------------------------------------------
C----Convert final g94 coordinate summary information to CHEMX format
C                                    7/7/97
C
      SUBROUTINE G94TOCHEMX (KIND)
      CHARACTER BUFFER*4000, LINE*80, BLANK*10, ID(100)*1,
     X          IN_FILE*30, OUT_FILE*30, BSLASH*1, SIGNAL*5,
     X          PIPE*1, SIGNAL2*5
      DIMENSION XYZ(3,100)
      DATA BLANK /'          '/,
     X     BSLASH /Z"5C"/,                            ! BSLASH = \
     X     PIPE /'|'/,                                ! PIPE = |
     X     SIGNAL(1:2), SIGNAL(3:5) /Z"5C5C",'0,1'/,  ! SIGNAL = \\0,1
     X     SIGNAL2 /'||0,1'/                          ! SIGNAL2 = ||0,1    
C
      PRINT 9
9     FORMAT (' Name of the G94 compressed-format xyz ', 
     x        'input file:',$)
      READ 27, IN_FILE            ! name of input file
27    FORMAT (A30)
      OPEN (UNIT=18, FILE=IN_FILE, STATUS='OLD')
      IF (KIND .EQ. 31) THEN
         PRINT 118
118      FORMAT (' Name of the G94 xyz output file:',$)
      ELSE
         PRINT 18
18       FORMAT (' Name of the CHEMX output file:',$)
      ENDIF
      READ 27, OUT_FILE           ! name of output file
      OPEN (UNIT=19, FILE=OUT_FILE, STATUS='UNKNOWN')
C
      NATOMS = 0
      ISTART = -69
10    ISTART = ISTART + 70
      READ (18,100) LINE
100   FORMAT (1X,A70)
      IF (LINE(1:10) .EQ. BLANK) GO TO 200
C----Store LINE into BUFFER 
      BUFFER(ISTART:ISTART+69) = LINE
      GO TO 10
200   IEND = ISTART + 69    ! limits of BUFFER
C----Blank LINE encoutered...assume the end
C     Search for '\\0,1' or '||0,1' ...this signals start of atoms
      DO 300 I=1,IEND-4
         K = 0
         K = INDEX(BUFFER, SIGNAL)
         IF (K .GT. 0) GO TO 304
         K = INDEX(BUFFER, SIGNAL2)
         IF (K .GT. 0) GO TO 304
300   CONTINUE
304   J1 = K + 5    ! '\' or '|' in front of 1st atom
310   J1 = J1 + 1   ! 1st atom ID
      NATOMS = NATOMS + 1
      ID(NATOMS) = BUFFER(J1:J1)
      J1 = J1 + 2   ! skip the ',' following atom ID
2080  format (1x,a70)
C----Find the next '\' or '|'
      K = 0
      K = INDEX(BUFFER(J1:IEND), BSLASH)
      IF (K .EQ. 0) K = INDEX(BUFFER(J1:IEND), PIPE) 
      READ (BUFFER(J1:J1+K-2),*) (XYZ(L,NATOMS),L=1,3)
C----Find the next '\' or '|'
         K = 0
         K = INDEX(BUFFER(J1:IEND), BSLASH)
         IF (K .EQ. 0) K = INDEX(BUFFER(J1:IEND), PIPE)
         IF (K .EQ. 0) THEN
            PRINT 3000
3000        FORMAT (' Unrecoverable error interpreting G94 input file',
     x              '....STOP')
            STOP
         ELSE
C----Two '\' or '|' in a row ('\\' or '||') signals the end
            IF (BUFFER(J1+K:J1+K) .EQ. BSLASH .OR.
     x          BUFFER(J1+K:J1+K) .EQ. PIPE) GO TO 600   ! make CHEMX file
         ENDIF
500   J1 = J1 + K - 1
      GO TO 310
C
C----Write CHEMX or G94 xyz files on unit # 19
600   IF (KIND .EQ. 31) GO TO 608
      WRITE (19,605) NATOMS
605   FORMAT (38X,3('   1.000')/21X,3('   90.00')/I4,'   0'/)
608   DO 620 I=1,NATOMS
         IF (KIND .EQ. 31) THEN
           WRITE (19,625) ID(I),  (XYZ(J,I),J=1,3)       ! G94 xyz format
625        FORMAT (1X,A1,3X,3F20.11) 
         ELSE
           WRITE (19,615) I, ID(I), I, (XYZ(J,I),J=1,3)  ! CHEMX format
615        FORMAT (I4,1X,A1,I3,1X,3F10.6)
         ENDIF
620   CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C----Transform cadpac DMA coordinate (Bohr) and cadpac input
C     coordinate (Angs) files to G94 xyz format
C
       SUBROUTINE CADPAC_DMAREL_TO_G94 (KIND)
C
       CHARACTER LINE*70, ENG_NAME(5)*3, US_NAME(5)*2,
     x           IN_FILE*30, OUT_FILE*30
       DIMENSION XYZ(3)     
       DATA ENG_NAME /'CAR', 'OXY', 'HPD', 'NIT', '???'/
       DATA US_NAME  /'C ',  'O ',  'H ',  'N ',  '? '/
C
      PRINT 9
9     FORMAT (' Name of the cadpac DMA output or cadpac input ', 
     x        'file for transformation:',$)
      READ 27, IN_FILE            ! name of input file
27    FORMAT (A30)
      OPEN (UNIT=18, FILE=IN_FILE, STATUS='OLD')
      PRINT 18
18    FORMAT (' Name of the G94 xyz output file:',$)
      READ 27, OUT_FILE           ! name of output file
      OPEN (UNIT=19, FILE=OUT_FILE, STATUS='UNKNOWN')
C
100   READ (18,'(A70)',END=1000) LINE
      DO 200 I=1,5      ! look for cadpac/dmarel atom names
         DO 190 J=1,68
            IF (INDEX(LINE(J:J+2), ENG_NAME(I)) .EQ. 0) GO TO 190
            IF (KIND .EQ. 32) THEN
               READ (LINE,180) XYZ    ! got one...get xyz's from cadpac DMA file
180            FORMAT (18X,3F10.6)
               XM = 0.529177249       ! conversion from Bohrs to Angs
            ELSE
               READ (LINE,181) XYZ    ! xyz's from cadpac input file
181            FORMAT (10X,3F13.6)
               XM = 1.0
            ENDIF
            DO L=1,3
                XYZ(L) = XYZ(L)*XM
            ENDDO
            WRITE (19,185) US_NAME(I), XYZ  ! G94 output format, replace 
185         FORMAT (2X,A2,3X,3F12.6)        !  Engl name with US name
            GO TO 100
190      CONTINUE
200   CONTINUE
      GO TO 100
C
1000  RETURN
      END
C----------------------------------------------------------------------
C----Function to determine if a 2-linked nitrogen is           12/31/05
C     N3 of a 1,2,3,4-tetrazole.                               12/31/05
C
      FUNCTION IQ_1234tetrazole_n3(I1) ! I1 is a 2-linked N  ! 12/31/05
C
C
      CHARACTER (LEN=5) :: label
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_1234tetrazole_n3, n, nc, nn, c3id, c3id_2 ! 12/31/05
      integer :: n3id, nnew, i1, i2, i3, i, j, n2id(2)
C
      IQ_1234tetrazole_n3 = 0     ! no
C
C----Must be linked to two 2-linked N's
C
C     (N3) 72 --> N---N <-- 73 (N4)
C                 "   "
C     (N2) 71 --> N   C-- 
C                  \ /  
C      (N1) 70 -->  N   
C                   |
C
      nn = 0        ! count # of 2-linked N's
      DO 100 I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
             NN = NN + 1
             n2id(nn) = i2        ! id of 2-linked N
         endif
100   CONTINUE
      IF (NN .ne. 2) RETURN       ! not a # 72 N
C----Look at the tw0 2-linked N's..is one N2 and one N4?
      n3id = 0
      c3id = 0
      do 200 j=1,2
         i2 = n2id(j)      ! a 2-linked N
             do i=1,2
                i3 = icon(i,i2)
                if (i3 .eq. i1) cycle         ! the original
                if (label(i3)(1:1) .eq. 'N' .and. ncon(i3) .eq. 3)  
     2              n3id = i3
                if (label(i3)(1:1) .eq. 'C' .and. ncon(i3) .eq. 3)  
     2              c3id = i3
             enddo
200   continue
C----Are both n3id and c3id non-0?
      if (.not. (n3id .ne. 0 .and. c3id .ne. 0)) return  
C----Check to determine that n3id and c3id are bonded
      do i=1,3
         if (icon(i,n3id) .eq. c3id) then
            iq_1234tetrazole_n3 = 72      ! success, param = 72
            return
         endif
      enddo
C
      return                ! no match...not N3
C
      END FUNCTION IQ_1234tetrazole_N3                       ! 1/1/06
C
C----------------------------------------------------------------
C----Function to determine if a 2-linked nitrogen is           1/1/06
C     N2 of a 1,2,3,4-tetrazole.                               1/1/06
C
      FUNCTION IQ_1234tetrazole_n2(I1) ! I1 is a 2-linked N  ! 1/1/06
C
C
      CHARACTER (LEN=5) :: label
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_1234tetrazole_n2, n3, n2                 ! 1/1/06
      integer :: iq_1234tetrazole_n3                         ! 1/1/06
      integer :: i1, i2, i, n2id
C
      IQ_1234tetrazole_n2 = 0     ! no
C
C----Must be linked to two 2-linked N's
C
C     (N3) 72 --> N---N <-- 73 (N4)
C                 "   "
C     (N2) 71 --> N   C-- 
C                  \ /  
C      (N1) 70 -->  N   
C                   |
C
      n2 = 0        ! count # of 2-linked N's
      n3 = 0        ! count # of 3-linked N's
      DO 100 I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 3) N3 = N3 + 1
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
            N2 = N2 + 1
            n2id = i2                   ! id of 2-linked N
         endif
100   CONTINUE
      IF (.not. (N2 .eq. 1 .and. n3 .eq. 1)) RETURN     ! not a # 71 N
C----Look at the 2-linked N..is it an N3 with param = 72 
      if (iq_1234tetrazole_n3(n2id) .eq. 72) then
         iq_1234tetrazole_n2 = 71          !yes, param = 71
      endif
C
      return
C
      END FUNCTION IQ_1234tetrazole_N2                       ! 1/1/06
C
C----------------------------------------------------------------
C----Function to determine if a 3-linked nitrogen is           1/3/06
C     N1 of a 1,2,3,4-tetrazole.                               1/3/06
C
      FUNCTION IQ_1234tetrazole_n1(I1) ! I1 is a 3-linked N  ! 1/1/06
C
      CHARACTER (LEN=5) :: label 
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_1234tetrazole_n2, c3, n2                 ! 1/1/06
      integer :: iq_1234tetrazole_n1                         ! 1/1/06
      integer :: i1, i2, i, n2id(2), ntotal                  ! 1/3/06
C
      IQ_1234tetrazole_n1 = 0     ! no
C
C----Must be linked to two 2-linked N's
C
C     (N3) 72 --> N---N <-- 73 (N4)
C                 "   "
C     (N2) 71 --> N   C-- 
C                  \ /  
C      (N1) 70 -->  N   
C                   |
C                   X <-- can be C, H or N                    1/3/06
C
      n2 = 0        ! count # of 2-linked N's
      c3 = 0        ! count # of 3-linked C's
      ntotal = 0
      do i=1,2
         n2id(i) = 0   ! id of 2-linked N
      enddo
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,I1)
         if (LABEL(I2)(1:1) .EQ. 'C' .or. 
     2       LABEL(I2)(1:1) .EQ. 'N' .or.
     3       LABEL(I2)(1:1) .EQ. 'H') ntotal = ntotal + 1
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) c3 = c3 + 1
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
            N2 = N2 + 1
            n2id(n2) = i2                   ! id of 2-linked N
         endif
100   CONTINUE
      IF (.not. ((N2 .eq. 1 .or. n2 .eq. 2) .and.  
     2           (c3 .eq. 1 .or. c3 .eq. 2) .and.
     3           (ntotal .eq. 3))) RETURN     ! not a # 70 N
C----Look at the 2-linked N's..if > 1, one must be a N2 with param = 71 
      do i=1,n2
         if (iq_1234tetrazole_n2(n2id(i)) .eq. 71) then
            iq_1234tetrazole_n1 = 70          !yes, param = 70
            return
         endif
      enddo
C
      return
C
      END FUNCTION IQ_1234tetrazole_N1                       ! 1/1/06
C
C----------------------------------------------------------------
C----Function to determine if a 2-linked nitrogen is           1/1/06
C     N4 of a 1,2,3,4-tetrazole.                               1/1/06
C
      FUNCTION IQ_1234tetrazole_n4(I1) ! I1 is a 2-linked N  ! 1/1/06
C
      CHARACTER (LEN=5) :: label
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_1234tetrazole_n4, c3, n2                 ! 1/1/06
      integer :: iq_1234tetrazole_n3                         ! 1/1/06
      integer :: i1, i2, i, n2id
C
      IQ_1234tetrazole_n4 = 0     ! no
C
C----Must be linked to two 2-linked N's
C
C     (N3) 72 --> N---N <-- 73 (N4)
C                 "   "
C     (N2) 71 --> N   C-- 
C                  \ /  
C      (N1) 70 -->  N   
C                   |
C
      n2 = 0        ! count # of 2-linked N's
      c3 = 0        ! count # of 3-linked C's
      DO 100 I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) c3 = c3 + 1
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
            N2 = N2 + 1
            n2id = i2                   ! id of 2-linked N
         endif
100   CONTINUE
      IF (.not. (N2 .eq. 1 .and. c3 .eq. 1)) RETURN     ! not a # 73 N
C----Look at the 2-linked N..is it an N3 with param = 72 
      if (iq_1234tetrazole_n3(n2id) .eq. 72) then
         iq_1234tetrazole_n4 = 73          !yes, param = 73
      endif
C
      return
C
      END FUNCTION IQ_1234tetrazole_N4                       ! 1/1/06
C
C----------------------------------------------------------------
C----Function to determine if a 2-linked nitrogen is           12/30/05
C     part of a 1,2,3-triazole.                                12/30/05
C
      FUNCTION IQ_123triazole_2(I1) ! I1 is a 2-linked N     ! 12/30/05
C
      CHARACTER (LEN=5) :: label
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_123triazole_2, n, nc, nn, c3id_1, c3id_2
      integer :: n3id, nnew, i1, i2, i, n2id(2)
C
      IQ_123triazole_2 = 0     ! no
C
C----Must be linked to a N, which is linked to 3-linked  C and N
C
C                 C---C
C                 "   "
C          27 --> N   N <-- 27
C                  \ / 
C           26 -->  N   
C                   |
C
      NC = 0        ! count # of 3-linked C's
      nn = 0        ! count # of 3-linked N's
      DO 100 I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then 
             NC = NC + 1
             c3id_1 = i2      ! id of first 3-linked C to be found
         endif
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 3) then
             NN = NN + 1
             n3id = i2        ! id of 3-linked N
         endif
100   CONTINUE
      IF (.not.(NC .eq. 1 .and. NN .eq. 1)) RETURN   ! not a # 27 N
C----Look at 3-linked N and count # of 2-linked N's 
      nn = 0               ! # of 2-linked N's
      nnew = 0             ! id of the 'new' 2-linked N
      DO 120 I=1,3
         i2 = icon(i,n3id) 
         IF (LABEL(I2)(1:1) .eq. 'N' .and. ncon(i2) .eq. 2) then
             nn = nn + 1
             if (i2 .ne. i1) nnew = i2
         endif
120   CONTINUE
      if (nn .ne. 2 .or. nnew .eq. 0) return ! 3-linked N not bonded to two 2-linked N's 
C----Look at the other 2-linked N (nnew), must be bonded to a 3-linked C and a N
      nc = 0
      nn = 0
      do 140 i=1,2
         I2 = ICON(i,nnew)
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then
             NC = NC + 1
             c3id_2 = i2     ! id of second 3-linked C 
         endif
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 3) NN = NN + 1
140   CONTINUE
      IF (NC .NE. 1 .OR. NN .NE. 1) RETURN     ! not a # 27 N
C----Finally check the c3id_1 and c3id_2 are bonded
      do 150 i=1,3
         i2 = icon(i,c3id_1)
         if (i2 .ne. c3id_2) go to 150    ! compare with id of second C
            iq_123triazole_2 = 27         ! param = 27
            return
150   continue
C
      return      ! return means a failure for iq_123triazole_2 
      END FUNCTION IQ_123triazole_2                     ! 12/30/05
C
C-----------------------------------------------------------------
C    Function to determine if 3-linked N is # 26 in triazole ring     12/30/05
C
      function iq_123triazole_3(i1)   ! i1 is a 3-linked N            12/30/05
C
      CHARACTER (LEN=5) :: label
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_123triazole_2, iq_123triazole_3, nn
      integer :: i1, i2, i, n2id(2)
C
      iq_123triazole_3 = 0
C
      nn = 0
      do 20 i=1,3
         i2 = icon(i,i1)      ! connected to 3-linked N
         if (label(i2)(1:1) .eq. 'N' .and. ncon(i2) .eq. 2) then
             nn = nn + 1
             n2id(nn) = i2
         endif
20    continue
      if (nn .ne. 2) return     ! not connected to two 2-linked N's
C----check that they are param = 27's
      do 40 i=1,2
         if (iq_123triazole_2(n2id(i)) .ne. 27) return
40    continue
C----Dropping thru loop means both are 27's
      iq_123triazole_3 = 26             ! param = 26
C
      end function iq_123triazole_3                 ! 12/30/05
C----------------------------------------------------------
C----Function to determine if a H bonded to a 3-linked           1/9/06
C     N that is part of a pyrrole ring.                          1/9/06
C
      FUNCTION IQ_pyrrole_nh(I1) ! I1 is a 3-linked N          ! 1/9/06
C
      CHARACTER (LEN=5) :: label
      LOGICAL CONNECT_FLAG
      COMMON /COM1/ NA, CELL(6), label(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_pyrrole_nh, n, nc, nn                      ! 1/9/06
      integer :: i1, i2, i3, i4, i, j, atom2id(2)              ! 1/9/06
C
      IQ_pyrrole_nh = 0     ! no
C
C----Must be linked to a N, which is linked to 3-linked  C and N
C
C                 X---X <-- C or N
C                 "   "
C      C or N --> X   X <-- C or N
C                  \ / 
C                   N   
C                   |
C                   H <-- 74                    

      NC = 0        ! count # of 3-linked C's
      nn = 0        ! count # of 3-linked N's
      n = 0
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,I1)
         if (label(i2)(1:1) .eq. 'H') go to 100
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then 
             NC = NC + 1
             n = n + 1
             atom2id(n) = i2      ! id of atom linked to N1
         endif
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. (ncon(i2) .eq. 3 .or.  
     2      ncon(i2) .eq. 2)) then
             NN = NN + 1
             n = n + 1
             atom2id(n) = i2      ! id of atom linked to N1
         endif
100   CONTINUE
      IF (.not.(NC .le. 2 .and. NN .le. 2 .and.  
     2   n .eq. 2)) RETURN           ! not proper neighbors
C----Look at two atoms in atom2id...do they have connected atoms? 
      i2 = atom2id(1)      ! start with the first of the two atoms
      if (.not. (ncon(i2) .eq. 2 .or. ncon(i2) .eq. 3)) return   ! not pi-bonded
      i4 = atom2id(2)
      if (.not. (ncon(i4) .eq. 2 .or. ncon(i4) .eq. 3)) return   ! not pi-bonded
      DO 120 I=1,ncon(i2)        
         i3 = icon(i,i2) 
         do 115 j=1,ncon(i4)
         IF (icon(j,i4) .eq. i3) then   ! connected..ring complete
            iq_pyrrole_nh = 74
            return
         endif
115      continue
120   CONTINUE
C
      return      ! return means a failure for iq_pyrrole_nh 
C
      END FUNCTION IQ_pyrrole_nh                          ! 1/3/06
C
C----------------------------------------------------------------
C----Function to determine if a 2-linked N is part              1/20/06
C     of an amino-imine                                         1/20/06
C
      FUNCTION IQ_aminoimine_2(I1) ! I1 is a 2-linked N       ! 1/23/06
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_aminoimine_2, nc, nn, nh                  ! 1/27/06
      integer :: i1, i2, i, j, atomid_n                       ! 1/27/06
C
      IQ_aminoimine_2 = 0     ! no
C
C----The functionality...
C          H or C                                             ! 1/27/06
C                \ 
C                 N---N===C (sp2)
C                /^   ^
C               C |   | 
C                 75  76 
C----2-linked N
      nc = 0
      nn = 0
       do i=1,2
          i2 = icon(i,i1)
          IF (name(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) nc = nc + 1 
          IF (name(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 3) then   ! 1/27/06
             nn = nn + 1                                           ! 1/27/06
             atomid_n = i2                                         ! 1/27/06
          endif
       enddo
       if (nc .ne. 1) return                                       ! 1/27/06
       if (nn .ne. 1) return                                       ! 1/27/06
       nc = 0                                                      ! 1/27/06
       nh = 0                                                      ! 1/27/06
       do i=1,3
          i2 = icon(i,atomid_n)                                    ! 1/27/06
          if (name(i2)(1:1) .eq. 'C') nc = nc + 1                  ! 1/27/06
          if (name(i2)(1:1) .eq. 'H') nh = nh + 1                  ! 1/27/06
       enddo
       if ((nc + nh) .ne. 2) return                                ! 1/27/06
       iq_aminoimine_2 = 76   ! param = 76, 2-linked N
       return    ! return with either a 0 or 76
      end function iq_aminoimine_2                                 ! 1/23/06
C
C--------------------------------------------------------------------------
C          
C----Function to determine if a 3-linked N is part              1/23/06
C     of an amino-imine                                         1/20/06
C
      FUNCTION IQ_aminoimine_3(I1) ! I1 is a 3-linked N       ! 1/23/06
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_aminoimine_2, iq_aminoimine_3, nn         ! 1/23/06
      integer :: i1, i2, i, j, atom2id(3)                     ! 1/23/06
C
C----The functionality...
C                \ 
C                 N---N===C (sp2)
C                /^   ^
C                 |   | 
C                 75  76 
C
       iq_aminoimine_3 = 0                                         ! 1/23/06
       nn = 0
         DO I=1,3                           ! examine the N
            I2 = ICON(I,I1)
            IF (name(I2)(1:1) .ne. 'N' .or. ncon(i2) .ne. 2) cycle
               nn = nn + 1 
               atom2id(nn) = i2
         enddo
         if (nn .ne. 1) return   ! not linked to one possible imino N
         if (iq_aminoimine_2(atom2id(1)) .ne. 76) return  ! NO
         iq_aminoimine_3 = 75     ! param = 75, 3-linked N       1/23/06
         return
C
      END FUNCTION IQ_aminoimine_3                                ! 1/23/06
C
C--------------------------------------------------------------
C
C----Function to determine if a 2-linked nitrogen is           1/31/06
C     part of a 1,3-thiazole.                                  1/31/06 
C
      FUNCTION IQ_thiazole_n(I1) ! I1 is a 2-linked N         ! 2/2/06
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_thiazole_n, nc, i, j, k, j1, j2, k1, k2 
      integer :: i1, i2, id_n(2), kk1, kk2, l3
C
C          78 --> N---C
C                 "   "
C                 C   C 
C                  \ / 
C           77 -->  S   
C                   
      IQ_thiazole_n = 0     ! no
C
      NC = 0                             ! count # of 3-linked C's
      DO I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (name(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then 
             NC = NC + 1
             id_n(nc) = i2               ! id of 3-linked C 
         endif
      enddo
      IF (NC .ne. 2) RETURN              ! not a possible # 78 N 
C----Look at two 3-linked C's...one should be bonded to 3-linked
C     C and the other to a 2-linked S 
      DO 120 i=1,3                ! examine first of 3-linked C's
         j1 = icon(i,id_n(1))     ! id's of atoms linked to first C 
         do 130 j=1,3             ! examine second of 3-linked C's
            j2 = icon(j,id_n(2))  ! id's of atoms linked to second C
            do 140 k1=1,ncon(j1)
               kk1 = icon(k1,j1)
               do 150 k2=1,ncon(j2)
                  kk2 = icon(k2,j2)
                  IF (name(kk1)(1:1) .eq. 'S' .or.      ! one of these must 
     2                name(kk2)(1:1) .eq. 'S') then     ! be S to keep going
                  do 160 l3=1,ncon(kk2)
                     if (kk1 .ne. icon(l3,kk2)) go to 160
                        iq_thiazole_n = 78      ! param = 78
                        return  
160               continue
                  endif
150            continue
140         continue
130      continue
120   CONTINUE
C 
      return           ! returning here means no luck
C
      end function iq_thiazole_n                             ! 1/31/06
C
C----------------------------------------------------------------
C
C----Function to determine if a 2-linked sulfur is             1/31/06
C     part of a thiazole.                                      1/31/06 
C
      FUNCTION IQ_thiazole_s(I1) ! I1 is a 2-linked S        ! 1/31/06
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_thiazole_s, iq_thiazole_n, nc, i, j
      integer :: i1, i2, i3, id_n(2)
C
      IQ_thiazole_s = 0     ! no
C
C          78 --> N---C
C                 "   "
C                 C   C 
C                  \ / 
C           77 -->  S   
C                   
      NC = 0                             ! count # of 3-linked C's
      DO I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (name(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then 
             NC = NC + 1
             id_n(nc) = i2               ! id of 3-linked C 
         endif
      enddo
      IF (NC .ne. 2) RETURN              ! not a possible # 79 S 
C----Look at two 3-linked C's...one should be bonded to a # 78 N
      DO 120 I=1,2
         i2 = id_n(i)      ! id of 3-linked C 
         do 130 j=1,3
            i3 = icon(j,i2)
            IF (name(I3)(1:1) .eq. 'N' .and. ncon(i3) .eq. 2) then  
              if (iq_thiazole_n(i3) .eq. 78) then
                 iq_thiazole_s = 77     ! param = 77
                 return      ! found the S
              endif
            endif
130      continue
120   CONTINUE
C 
      return           ! returning here means no luck
C
      end function iq_thiazole_s                             ! 1/31/06
C
C-----------------------------------------------------------------
C----Function to determine if a 3-linked C is                   2/3/06
C     carbonyl C of an ester or ahhydride                       2/3/06 
c
      FUNCTION IQ_ester_c(I1) ! I1 is a 3-linked C            ! 2/3/06
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_ester_c, nc, i, j, k, j1, j2, k1, k2 
      integer :: i1, i2, no1, no2, ido2(2)
C
C         80 -->  O               O
C                 "               " 
C         79 -->  C   C or i-->   C 
C                / \ /             \
C               C   O <--81         C
C                   
      IQ_ester_c = 0     ! no
C
      no1 = 0                           ! count # of 1-linked O's
      no2 = 0                           ! count # of 2-linked O's
      nc = 0                            ! count number of C's
      DO I=1,3                    
         I2 = ICON(I,I1)
         IF (name(I2)(1:1) .EQ. 'C') NC = NC + 1
         if (name(i2)(1:1) .eq. 'O' .and. ncon(i2) .eq. 1) no1 = no1 + 1
         if (name(i2)(1:1) .eq. 'O' .and. ncon(i2) .eq. 2) then
            no2 = no2 + 1
            ido2(no2) = i2
         endif
      enddo 
      IF (.not. (nc .eq. 1 .and. no1 .eq. 1 .and. no2 .eq. 1)) RETURN ! no 
C----Look at 2-linked O...must if bonded to 2 C's
      nc = 0
      DO i=1,2                
         i2 = icon(i,ido2(1))     ! id's of atoms linked to 2-linked O 
         if (name(i2)(1:1) .eq. 'C') nc = nc + 1
      enddo
C
      if (nc .ne. 2) return     ! no
C
      iq_ester_c = 79      ! param = 79
C
      end function iq_ester_c                       ! 2/3/06
C------------------------------------------------------------------
C----Function to determine if a 2-linked N of                  2/13/06
C     a diazo linkage                                          2/13/06 
C
      FUNCTION IQ_diazo_n(I1) ! I1 is a 2-linked N           ! 2/13/06
C
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_diazo_n, nc, i, nn_2                     ! 2/10/06
      integer :: idn2, i1, i2                                ! 2/10/06
      real :: distance                                       ! 2/10/06
C
C             C--N==N--C                                       2/10/06
C              (82)                                            2/10/06
C                                                  
      iq_diazo_n = 0
C
      nc = 0                       ! number of C's
      nn_2 = 0                     ! number of 2-linked n's
      DO I=1,2                    
         I2 = ICON(I,I1)
         IF (name(I2)(1:1) .EQ. 'C') NC = NC + 1
         if (name(i2)(1:1) .eq. 'N' .and. ncon(i2) .eq. 2) then
            nn_2 = nn_2 + 1
            idn2 = i2
         endif
      enddo 
      IF (.not. (nc .eq. 1 .and. nn_2 .eq. 1)) RETURN ! no 
C----Look at the other 2-linked N...must if bonded to a C
      nc = 0
      DO i=1,2                
         i2 = icon(i,idn2)  ! id's of atoms linked to other 2-linked N 
         if (name(i2)(1:1) .eq. 'C') nc = nc + 1
      enddo
      if (nc .ne. 1) return     ! no
C
C----One final check...N=N distance must be .le. 1.286 Angs
      if (distance(i1,idn2) .gt. 1.286) return ! too BIG
C
      iq_diazo_n = 82      ! param = 82
C
      end function iq_diazo_n                       ! 2/13/06
C------------------------------------------------------------------
C
C------------------------------------------------------------------
C----Function to determine if a 3-linked nitrogen is a
C     pyridinium-N;   -N(+)-H                                7/20/06 
C
      FUNCTION iq_pyridinium_n(I1)    ! I1 is a 3-linked N
C                                       with 1 bond to H  
C
      CHARACTER NAME*5
      COMMON /COM1/ NA, CELL(6), NAME(200), XYZ(3,200),
     x              ICON(6,200), NCON(200)
      integer :: iq_pyridinium_n, c3, n2, nh          
      integer :: i1, i2, i, idH, id(3), id2(6), id3(6)
      integer :: j, j2, nsave, nsave2, nsave3    
C
      iq_pyridinium_n = 0     ! no
C
C----Must be linked to two 2-linked N's
C                   Z
C                  / \
C                 Z   Z 
C                 |   |
C                 Z   Z      Z = 3-linked C or 2-linked N
C                  \ /
C           85 -->  N
C                   |
C           86 -->  H 
C
      n2 = 0        ! count # of 2-linked N's
      c3 = 0        ! count # of 3-linked C's
      nh = 0
      nsave = 0
      do i=1,3
         id(i) = 0      ! id of 2-linked N or 3-linked C
      enddo
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,I1)
         if (name(I2)(1:1) .EQ. 'H') then
            nh = nh + 1
            idH = i2
         endif
         IF (name(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then
             c3 = c3 + 1
             nsave = nsave + 1
             id(nsave) = i2
         endif
         IF (name(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
            N2 = N2 + 1
            nsave = nsave + 1
            id(nsave) = i2                   ! id of 2-linked N
         endif
100   CONTINUE
      if (nH .ne. 1) return       ! must be 1 H on N
      IF (nsave .ne. 2) return    ! attached N's + C's must = 2
C----id contains the C and/or N #'s of the 2 non-H atoms linked
C     to the pyridinium N
      nsave2 = 0    ! nsave2 will count 2nd atoms out from pyridinium N
      do i=1,2      ! look at attached non-H atoms, eliminate pyridinium N 
         i2 = id(i)
            do j=1,ncon(i2)
               j2 = icon(j,i2)
               if (name(j2)(1:1) .eq. 'H') cycle ! don't need H
               if (j2 .eq. i1) cycle              ! i1 = ID of pyridinium N
               if (name(j2)(1:1) .eq. 'N' .and.  
     2             ncon(j2) .ne. 2) cycle
               if (name(j2)(1:1) .eq. 'C' .and.  
     2             ncon(j2) .ne. 3) cycle 
               nsave2 = nsave2 + 1
               id2(nsave2) = j2
            enddo
      enddo
C----Make a list of ID's of atoms attached to 2nd atoms out...if 2 are same
C     then have a 6-membered ring.  These would be 3rd atoms out.
      nsave3 = 0
      do i=1,nsave2
         i2 = id2(i)
         do j=1,ncon(i2)
            j2 = icon(j,i2)
            if (name(j2)(1:1) .eq. 'H') cycle
            if (name(j2)(1:1) .eq. 'N' .and.  
     2          ncon(j2) .ne. 2) cycle
            if (name(j2)(1:1) .eq. 'C' .and.  
     2          ncon(j2) .ne. 3) cycle
            nsave3 = nsave3 + 1
            id3(nsave3) = j2  ! ID's of all C and N atoms connected
         enddo                !  to 2nd atoms out
      enddo
      if (nsave3 .lt. 2) return
C----To establish 6-membered ring, two of atoms in id3 must be the same
      do i=1,nsave3-1
         do j=i+1,nsave3
            if (id3(i) .eq. id3(j)) go to 110
         enddo
      enddo
C----dropping thru the loop means we've failed to find a ring
      return
110   iq_pyridinium_n = 85   ! pyridinium N+, param = 85
      return
C
      end 
C--------------------------------------------------------------------
