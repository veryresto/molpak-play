C----VOLUME-ADDITIVITY-DENSITIES.F                                       7/15/05
C
      PARAMETER (NPOTS=166, MAX_PARAM=108, NTYPES_ATOMS=11)             ! 6/7/05
      CHARACTER*1 LINE(80), RF2(8), OUTNAM(15), DATNAM(7), LINE2(80)
      CHARACTER LINE80*80, ATSYM*1, COMMENT*120,                       ! 2/26/00
     x          ATSYM2*2, MWOUT*36, DESCRIPTION(NPOTS)*60,
     x          DES_TEMP*60
      EQUIVALENCE (LINE, LINE80, LINE2), (LLINE, LINE500)
      CHARACTER LINE3(80)*1, LLINE(500)*1, LINE4(80)*1, LINE500*80,
     x          DIR*80
      CHARACTER*7 DATNM2
      CHARACTER*15 OUTFIL
      CHARACTER*48 INFILE
      CHARACTER*8 REFCOD, RF, REFCOD_2
      CHARACTER*8 SPGRP
      CHARACTER*5 LABEL(500), NEW_LABEL
      CHARACTER*2 LABEL_OUT                                            ! 2/26/00
      CHARACTER*2 ELEMENTS_OUT(40)                                     ! 2/28/00
      INTEGER*2 IP(6), IS_OUT(500)
      EQUIVALENCE (RF2, REFCOD), (OUTNAM, OUTFIL), (DATNAM, DATNM2)
      EQUIVALENCE (LINE(2), RF)
      DIMENSION N_AT_TYPE(NPOTS)
      DIMENSION L_ORDER(NPOTS), DATA_TABLE(NPOTS,4), DATA_TEMP(4)
      DIMENSION ICELL(6), NVALUE(37), IXYZ(3,500), DEN(2),           ! 2/26/00
     x          AT_WTS(NTYPES_ATOMS)                                 ! 2/26/00
      DIMENSION XYZ(3,500), R(3,3), O_XYZ(3), IND(NPOTS),
     x          NATTYPES(NPOTS), NSIZES(4)
      EQUIVALENCE (IXYZ, XYZ)
      COMMON /COM1/ NCON(500), ICON(5,500), NEW_LABEL(500),
     x                 OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      COMMON /CONTENTS/ NCON_FDAT(NTYPES_ATOMS),                     ! 2/26/00
     x                  NCON_FBIB(NTYPES_ATOMS),ATSYM(NTYPES_ATOMS)  ! 2/26/00
C----N_PER_UNIT counts # hits as function of molec/asym unit         ! 3/6/00
      COMMON /ELEMENTS/ N_AT_DIFF, N_PER_UNIT(8), ATSYM2(NTYPES_ATOMS) ! 3/6/00
      COMMON /PROBLEMS/ JPROB, COMMENT, REFCOD_2, ICSP3, JCSP3
      DATA DEN /10000., 100000./
      DATA AT_WTS /12.011,1.0079,18.9984,14.0067,15.9994,79.904,     ! 6/7/05
     x             35.4527,32.066,30.9738,10.811,126.904/            ! 6/7/05
C                     C     H       F       N       O      Br        ! 6/7/05
C                     Cl    S       P       B       I                ! 6/7/05
C----Br will be recoded as atom ID "X" to use only 1 character       ! 2/26/00
C     And CL is recorded as Z.
C
C----File # 12 = file that contains molpak...input names
      OPEN (UNIT = 12, FILE = 'MOLPAK.NAME', STATUS='OLD')
      READ (12,'(A40)') INFILE
2     CLOSE (UNIT=12)
C
C----Open the molpak xyz input file (unit # 10)
      OPEN (UNIT=10, FILE=INFILE, STATUS='OLD')
C
C----3rd line of unit 9...IGROUPS, ICC1, ITWO_ALL, ICOORD, JPROB
C 1.  IGROUPS =  (0/1) for (separate atoms)/(groups atom into
C                            functional group units.....
C                               nitro (NO2)
C                               cubane nitro
C                               nitramine (N-NO2)
C                               nitrile (C#N)
C                               diazo (C=N(+)=N(-)
C                               carboxylic acid (CO2H)
C                               carbonyl (C=O)
C                               amide, imide, urea carbonyl
C                               hydroxyl (OH), sp2 & sp3 linked
C                               nitroso (N=O)
C                               C-N(-)-N(C)(+)=N-C
C 2.  ICC1 =     (0/1) for (do/do not) produce *.cc1 files
C 3.  ITWO_ALL = (0/N) for (do not/do) accept structures with 2-N molecules  ! 3/6/00
C                             per symmetric unit.                            ! 3/6/00
C 4.  ICOORD =   (0/1/2) for (do not/
C                             print the connectivity list for problem structures/
C                             print list for all structures)
C 5.  JPROB =    (0/1) for (do not/do) write the PROBLEM refcodes on unit # 15
C 6.  ICSP3 =    (0/1) for (do not/do) separate Csp3 with 2 bonds in rings into
C                             3, 4, 5, 6, 7, 8 and > 8 size rings
C 7.  JCSP3 =    (0/1) FOR (do not/do) calculate separate increments for CH3,
C                             CH2 (not in ring) and CH (not in ring)
C 8.  JCALCD =   (0/1) for (do not/do) calculate the linear and nonlinear
C                             volume volumes and densities...requires the table2
C                             file
      IGROUPS = 1
      ICC1 = 0
      ITWO_ALL = 1
      ICOORD = 0
      JPROB = 0
      ICSP3 = 1
      JCSP3 = 1
      JCALCD = 1
C
C----Pick up volumes and weights from table3 & store in DATA_TABLE?        ! 12-18-96
      DIR='/PREDICTIONS/UTILITIES/table3'
C     DIR='/export/software/predictions-100cyc04/UTILITIES/table3'         ! 6/12/05
771   OPEN (UNIT=18, FILE=DIR, STATUS='OLD')
         READ (18,*)                                                       ! 12/16/96
         READ (18,*)                                                       ! 12/16/96
         DO I=1,MAX_PARAM                                                  ! 12/18/96
            READ (18,785) M, DATA_TEMP, DES_TEMP
785         FORMAT (6X,I3,1X,F7.3,6X,2(F8.3,6X),F9.4,9X,A60)               ! 6//8/05
            L_ORDER(I) = M
            DO L=1,4
	       DATA_TABLE(M,L) = DATA_TEMP(L)
            ENDDO
            DESCRIPTION(M) = DES_TEMP
         ENDDO                                                        ! 12/16/96
      CLOSE (UNIT=18)                                                 ! 12/16/96
C----Report results of DATA_TABLE construction                        ! 12-18-96
      OPEN (UNIT=45,FILE='additivity_data',STATUS='UNKNOWN')
      WRITE (45,777)
777   FORMAT (' #  old#      lin V    nonl V    lin d       wt',
     x        '     atom/group'/
     x        '--- ----     ------    ------    ------    ------',
     x        '   ----------')
      J = 0                                                           ! 12-18-96
      DO 790 I=1,MAX_PARAM                                            ! 12-18-96
         M = L_ORDER(I)
         WRITE (45,789) I, M, (DATA_TABLE(M,K),K=1,4),
     x              DESCRIPTION(M)
789      FORMAT (I3,'(#',i3,')',3F10.3,F10.4,2X,A60)
790   CONTINUE                                                        ! 12-18-96
      CLOSE (UNIT=45)
C
      IF (JPROB .NE. 0) OPEN (UNIT = 15, FILE = 'problems.gcd', STATUS =
     X                        'UNKNOWN', FORM = 'FORMATTED')
C
      LCOORD = ICOORD
C
C----Begin reading molpak xyz file.  Look for 'ATOM' lines.
120   NATOMS = 0
122   READ (10,'(A)',END=128) LINE80
      IF (LINE80(1:4) .NE. 'ATOM') GO TO 122
      NATOMS = NATOMS + 1
      READ (LINE80,125) LABEL(NATOMS)(1:2), (XYZ(J,NATOMS),J=1,3)
125   FORMAT (5X,A2,4X,3F10.6)
      GO TO 122
C
C----Search fdat atoms labels for BR & replace as X             ! 3/4/00
C     and CL & replace as Z                                     ! 3/4/00
128    DO 124 I=1,NATOMS                                       ! 2/26/00
           IF (LABEL(I)(1:2) .NE. 'BR') GO TO 655               ! 2/28/00
               LABEL(I)(1:2) = 'X '                             ! 2/28/00
               GO TO 124                                        ! 2/28/00
655        IF (LABEL(I)(1:2) .NE. 'CL') GO TO 124               ! 3/4/00
               LABEL(I)(1:2) = 'Z '                             ! 3/4/00
124     CONTINUE                                                ! 2/28/00
C
C----Sort atom list...C's 1st, O's 2nd, H's last                ! 3/20/00
      NA = 0
      DO 161 IPASS=1,4
      DO 160 I=1,NATOMS
        GO TO (181,187,182,183), IPASS
181     IF (LABEL(I)(1:1) .NE. 'C') GO TO 160  ! C's at top
           GO TO 184
187     IF (LABEL(I)(1:1) .NE. 'O') GO TO 160  ! O's 2nd
           GO TO 184
182     IF (LABEL(I)(1:1) .EQ. 'H' .OR.
     x      LABEL(I)(1:1) .EQ. 'C' .OR.
     x      LABEL(I)(1:1) .EQ. 'O') GO TO 160
           GO TO 184
183     IF (LABEL(I)(1:1) .NE. 'H') GO TO 160  ! H's at bottom
184     DO 152 L=1,3
             O_XYZ(L) = XYZ(L,I)
152     CONTINUE
        NA = NA + 1
        NEW_LABEL(NA) = LABEL(I)
        DO L = 1,3
           OR_XYZ(L,NA) = O_XYZ(L)
        ENDDO
160   CONTINUE
161   CONTINUE
      DO 162 I = 1,NTYPES_ATOMS                             ! 2/26/00
         NCON_FDAT(I) = 0
162   CONTINUE
      DO 700 I = 1, NA
         DO 699, J=1,NTYPES_ATOMS                           ! 2/26/00
            IF(NEW_LABEL(I)(1:1) .NE. ATSYM(J))  GO TO 699
            NCON_FDAT(J) = NCON_FDAT(J) + 1
            GO TO 700
699      CONTINUE
700   CONTINUE
C----Calc molecular weight,
701   WT = 0.0
      DO 175 I = 1,NTYPES_ATOMS                             ! 2/26/00
         WT = WT + NCON_FDAT(I)*AT_WTS(I)
175   CONTINUE
C----Calc H fraction
      FRACT_H = NCON_FDAT(2)*AT_WTS(2)/WT
C----Build molecular formula
      IPOSITION = 0
      DO 222 I=1,NTYPES_ATOMS
          IF (NCON_FDAT(I) .EQ. 0) GO TO 222
          IPOSITION = IPOSITION + 1
          WRITE (MWOUT(IPOSITION:IPOSITION+3),218)
     x                    ATSYM2(I), NCON_FDAT(I)
218       FORMAT (A2,I2)
          IPOSITION = IPOSITION + 3
222   CONTINUE
      LAST = IPOSITION
      DO 235 I=1,IPOSITION
         K = INDEX(MWOUT(1:LAST),' ')
         IF (K .EQ. 0) GO TO 238
         IF (K .LE. LAST) GO TO 237
            LAST = LAST - 1
            GO TO 238
237      LAST = LAST -1
         DO 236 L=K,LAST
            MWOUT(L:L) = MWOUT(L+1:L+1)
236      CONTINUE
235   CONTINUE
C----Molecular formula and weight output
238   PRINT 195, MWOUT(1:LAST), WT, FRACT_H
195   FORMAT (' Predicted densities from volume additivity data base...'
     x        /' Molecular formula: ',A/
     x        ' Molecular mass:',F8.3,' g/mol; H fraction:',F6.3)
C
      ICOORD = LCOORD
      CALL CONNECTION(ITEST,ICOORD)  ! check for more than 4 connections
      IF (ITEST .EQ. 0) GO TO 282
         N = NCON(ITEST)
         WRITE (COMMENT, 283) N, ITEST, NEW_LABEL(ITEST)
         WRITE (6,283) N, ITEST, NEW_LABEL(ITEST)
283      FORMAT (' $$$$ PROBLEMS:',I3,' connections found for atom',
     x           I3,' (',A5,')')
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
         IF (N .GE. 1) THEN
            WRITE (6,285) (NEW_LABEL(ICON(I,ITEST)),I=1,N)
285         FORMAT ('      connected atoms = ',5(A5,1X))
         ENDIF
         WRITE (6,284)
284      FORMAT ('Volume additivity calculation terminated')        ! 6/7/05
         STOP
C
282   IF (ICC1 .EQ. 0) WRITE (11,199) NA
199   FORMAT (I3)
      DO I=1,NA              ! initialize IT_TABLE....0 means
         IT_TABLE(I) = 0     ! atom has not been typed
         I3D(I) = 0
         KEEP_NPOT(I) = 0    ! initialize KEEP_NPOT
      ENDDO
      DO I=1,NPOTS
         N_AT_TYPE(I) = 0
      ENDDO
      DO 400 I=1,NA
         IF (IT_TABLE(I) .NE. 1) THEN    ! already been done?
           CALL ICHEM3D(I3D(I),I,KNPOT)  ! no
           KEEP_NPOT(I) = KNPOT
           IF(KNPOT .EQ. 999) GO TO 1100
         ENDIF
400   CONTINUE
C----Count up the atom type information
      DO 405 I=1,NA
         K = KEEP_NPOT(I)
         N_AT_TYPE(K) = N_AT_TYPE(K) + 1
405   CONTINUE
C----Write coordinates to Chem3D file
      IF (ICC1 .EQ. 0) THEN
      DO I=1,NA
         LABEL_OUT(2:2) = ' '                                     ! 2/26/00
         LABEL_OUT(1:1) = NEW_LABEL(I)(1:1)                       ! 2/26/00
         IF (LABEL_OUT(1:1) .EQ. 'X') LABEL_OUT = 'BR'            ! 2/26/00
         IF (LABEL_OUT(1:1) .EQ. 'Z') LABEL_OUT = 'CL'            ! 3/4/00
         WRITE (11,320) LABEL_OUT,                                ! 2/26/00
     x      (OR_XYZ(L,I),L=1,3), I3D(I), (ICON(L,I), L=1,NCON(I))
320      FORMAT (2X,A2,F11.6,2F12.6,5I5)                          ! 2/26/00
      ENDDO
      ENDIF
C
C----Count up number of atoms of each type (or groups as appropriate)
C
C----Count number of C & S for C=S                                ! 3/11/00
      IF (N_AT_TYPE(115) .EQ. N_AT_TYPE(116)) THEN                ! 3/11/00
         N_AT_TYPE(117) = N_AT_TYPE(115)  ! count groups          ! 3/11/00
         N_AT_TYPE(115) = 0      ! zero C's                       ! 3/11/00
         N_AT_TYPE(116) = 0      ! zero S's                       ! 3/11/00
      ELSE                                                        ! 3/12/00
         WRITE (COMMENT,718) N_AT_TYPE(115), N_AT_TYPE(116)       ! 3/12/00
718      FORMAT (' $$$$ PROBLEMS: unbalanced C=S atoms, # C''s =',! 3/12/00
     x           I2,', # S''s =',I2)                              ! 3/12/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/12/00
         GO TO 77                                                 ! 3/12/00
      ENDIF                                                       ! 3/12/00
C
C----Count number of S & O sulfoxides                             ! 3/19/00
      IF (N_AT_TYPE(120) .EQ. N_AT_TYPE(119)) THEN                ! 3/19/00
         N_AT_TYPE(121) = N_AT_TYPE(120)  ! count groups          ! 3/19/00
         N_AT_TYPE(120) = 0      ! zero O's                       ! 3/19/00
         N_AT_TYPE(119) = 0      ! zero S's                       ! 3/19/00
      ELSE                                                        ! 3/19/00
         WRITE (COMMENT,818) N_AT_TYPE(119), N_AT_TYPE(120)       ! 3/19/00
818      FORMAT (' $$$$ PROBLEMS: unbalanced S & O atoms in ',    ! 3/19/00
     x           'sulfoxides, # S''s =',I2,', # O''s =',I2)       ! 3/19/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/19/00
         GO TO 77                                                 ! 3/19/00
      ENDIF                                                       ! 3/19/00
C
C----Count number of S & O sulfones                               ! 3/19/00
      IF (2*N_AT_TYPE(122) .EQ. N_AT_TYPE(123)) THEN              ! 3/19/00
         N_AT_TYPE(124) = N_AT_TYPE(122)  ! count groups          ! 3/19/00
         N_AT_TYPE(123) = 0      ! zero O's                       ! 3/19/00
         N_AT_TYPE(122) = 0      ! zero S's                       ! 3/19/00
      ELSE                                                        ! 3/19/00
         WRITE (COMMENT,820) N_AT_TYPE(119), N_AT_TYPE(120)       ! 3/19/00
820      FORMAT (' $$$$ PROBLEMS: unbalanced S & O atoms in ',    ! 3/19/00
     x           'sulfones, # S''s =',I2,', # O''s =',I2)         ! 3/19/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/19/00
         GO TO 77                                                 ! 3/19/00
      ENDIF                                                       ! 3/19/00
C
C----Count number of S's & O's in sulfonates                      ! 3/20/00
      IF ((N_AT_TYPE(125) .NE. N_AT_TYPE(127)) .OR.               ! 3/20/00
     x    (2*N_AT_TYPE(125) .NE. N_AT_TYPE(126))) THEN            ! 3/20/00
         WRITE (COMMENT,822) N_AT_TYPE(125), N_AT_TYPE(126),      ! 3/20/00
     x                       N_AT_TYPE(127)                       ! 3/20/00
822      FORMAT (' $$$$ PROBLEMS: unbalanced S & O atoms in ',    ! 3/19/00
     x           'sulfonates, # S''s =',I2, ', # O2''s =',I2,     ! 3/20/00
     x           ', # O''s =',I2)                                 ! 3/20/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/19/00
         GO TO 77                                                 ! 3/19/00
      ELSE                                                        ! 3/20/00
         N_AT_TYPE(128) = N_AT_TYPE(125)                          ! 3/20/00
         N_AT_TYPE(125) = 0          ! zero S's                   ! 3/20/00
         N_AT_TYPE(126) = 0          ! zero O's                   ! 3/20/00
         N_AT_TYPE(127) = 0          ! zero O's                   ! 3/20/00
      ENDIF                                                       ! 3/20/00
C
C----Count number of S's, O's and N's in sulfonamides             ! 3/29/00
      IF ((N_AT_TYPE(140) .NE. N_AT_TYPE(141)) .OR.               ! 3/29/00
     x    (2*N_AT_TYPE(140) .NE. N_AT_TYPE(139))) THEN            ! 3/29/00
         WRITE (COMMENT,824) N_AT_TYPE(140), N_AT_TYPE(139),      ! 3/29/00
     x                       N_AT_TYPE(141)                       ! 3/29/00
         PRINT 824, N_AT_TYPE(140), N_AT_TYPE(139),               ! 3/29/00
     x                       N_AT_TYPE(141)                       ! 3/29/00
824      FORMAT (' $$$$ PROBLEMS: unbalanced S, O and N atoms in',! 3/29/00
     x           ' sulfonamides, # S''s =',I2, ', # O''s =',I2,   ! 3/29/00
     x           ', # N''s =',I2)                                 ! 3/29/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/19/00
         GO TO 77                                                 ! 3/19/00
      ELSE                                                        ! 3/20/00
         N_AT_TYPE(142) = N_AT_TYPE(140)                          ! 3/20/00
         N_AT_TYPE(140) = 0          ! zero S's                   ! 3/29/00
         N_AT_TYPE(139) = 0          ! zero O's                   ! 3/29/00
         N_AT_TYPE(141) = 0          ! zero N's                   ! 3/29/00
      ENDIF                                                       ! 3/20/00
C
C----Count number of S's, O's, N's iand C's in sulfonimines       ! 11/3/00
      IF ((N_AT_TYPE(153) .NE. N_AT_TYPE(154)) .OR.               ! 11/3/00
     x    (2*N_AT_TYPE(153) .NE. N_AT_TYPE(152)) .OR.             ! 11/3/00
     x    (N_AT_TYPE(153) .NE. N_AT_TYPE(155))) THEN              ! 11/3/00
         WRITE (COMMENT,2824) N_AT_TYPE(153), N_AT_TYPE(152),     ! 11/3/00
     x                       N_AT_TYPE(154), N_AT_TYPE(155)       ! 11/3/00
         PRINT 2824, N_AT_TYPE(153), N_AT_TYPE(152),              ! 11/3/00
     x                       N_AT_TYPE(154), N_AT_TYPE(155)       ! 11/3/00
2824     FORMAT (' $$$$ PROBLEMS: unbalanced S, O and N atoms in',! 11/3/00
     x           ' sulfonimines, # S''s =',I2, ', # O''s =',I2,   ! 11/3/00
     x           ', # N''s =',I2,', # C''s =',I2)                 ! 11/3/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 11/3/00
         GO TO 77                                                 ! 11/3/00
      ELSE                                                        ! 11/3/00
         N_AT_TYPE(151) = N_AT_TYPE(153)                          ! 11/3/00
         N_AT_TYPE(153) = 0          ! zero S's                   ! 11/3/00
         N_AT_TYPE(152) = 0          ! zero O's                   ! 11/3/00
         N_AT_TYPE(154) = 0          ! zero N's                   ! 11/3/00
         N_AT_TYPE(155) = 0          ! zero C's                   ! 11/3/00
      ENDIF                                                       ! 11/3/00
C
C----Count number of S's & H's in thiols                          ! 3/22/00
      IF (N_AT_TYPE(131) .NE. N_AT_TYPE(130)) THEN                ! 3/22/00
         WRITE (COMMENT,924) N_AT_TYPE(131), N_AT_TYPE(130)       ! 3/22/00
924      FORMAT (' $$$$ PROBLEMS: unbalanced S & H atoms in ',    ! 3/22/00
     x           'thiols, # S''s =',I2, ', # H''s =',I2)          ! 3/22/00
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/19/00
         GO TO 77                                                 ! 3/19/00
      ELSE                                                        ! 3/20/00
         N_AT_TYPE(132) = N_AT_TYPE(131)                          ! 3/22/00
         N_AT_TYPE(131) = 0          ! zero S's                   ! 3/22/00
         N_AT_TYPE(130) = 0          ! zero H's                   ! 3/22/00
      ENDIF                                                       ! 3/20/00
C
C----Count number of imides with C-C(=O)-NH-C(=O)-C structure
C     identified by a central N with N_AT_TYPE = 150
      IF (N_AT_TYPE(150) .EQ. 0) GO TO 925
      IF ((N_AT_TYPE(150) .LE. N_AT_TYPE(42)) .AND.     ! compare N with H on HN
     x    (2*N_AT_TYPE(150) .LE. N_AT_TYPE(79)) .AND.   ! compare 2*N with C of C=O
     X    (2*N_AT_TYPE(150) .LE. N_AT_TYPE(80))) THEN   ! compare 2*N with O or C=O
         N_AT_TYPE(42) = N_AT_TYPE(42) - N_AT_TYPE(150)   ! # H on HN
         N_AT_TYPE(79) = N_AT_TYPE(79) - 2*N_AT_TYPE(150) ! # C in C=O
         N_AT_TYPE(80) = N_AT_TYPE(80) - 2*N_AT_TYPE(150) ! # O in C=O
      ELSE
         WRITE (COMMENT,926) N_AT_TYPE(150), N_AT_TYPE(43),
     x                       N_AT_TYPE(79), N_AT_TYPE(80)
926      FORMAT (' $$$$ PROBLEMS: # of imide N''s =',I2,
     x           ' # amide H''s =',I2/
     x           '                # carbony C''s =',I2,
     x           ' # carbonyl O''s =',I2)
         IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                  ! 3/19/00
         GO TO 77                                                 ! 3/19/00
      ENDIF
C
C----Divide # azide N's by 3 if multiple of 3
925   IF (MOD(N_AT_TYPE(39),3) .EQ. 0) THEN
        IF (N_AT_TYPE(39) .GE. 3) N_AT_TYPE(39) = N_AT_TYPE(39)/3
      ELSE
        WRITE (COMMENT, 715) N_AT_TYPE(39)
        PRINT 715, N_AT_TYPE(39)
715     FORMAT (' $$$$ PROBLEMS: number of azide N''s (',I2,
     x          ') is not a multiple of 3')
        IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
        GO TO 77
      ENDIF
C----count nitriles (C#N)
         IF (N_AT_TYPE(45) .EQ. N_AT_TYPE(46)) THEN  ! # C's must = # N's
            N_AT_TYPE(77) = N_AT_TYPE(46)
            N_AT_TYPE(46) = 0     ! N
            N_AT_TYPE(45) = 0     ! C
         ELSE
            WRITE (COMMENT, 315) N_AT_TYPE(45), N_AT_TYPE(46)
            PRINT 315, N_AT_TYPE(45), N_AT_TYPE(46)
315         FORMAT (' $$$$ PROBLEMS: unbalanced C#N atoms, # C''s =',I2,
     x              ', # N''s =',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C----count diazo groups (C=N=N)
         N_AT_TYPE(72) = N_AT_TYPE(71)    ! # C=N=N = # C's
         N_AT_TYPE(71) = 0                ! # C's set to 0
         N_AT_TYPE(69) = N_AT_TYPE(69) - 2*N_AT_TYPE(72)
C----count nitro groups....C(sp3)-linked NO2's
         IF (N_AT_TYPE(7) .NE. 0) THEN
            IF (N_AT_TYPE(6)/2 .LT. N_AT_TYPE(7)) THEN                 ! 10-21-96
               WRITE (COMMENT, 740) N_AT_TYPE(7), N_AT_TYPE(6)
               PRINT 740, N_AT_TYPE(7), N_AT_TYPE(6)
740            FORMAT (' $$$$ PROBLEMS: unbalanced normal NO2 atoms, ',
     x                 '# N''s =',I2,', # O''s =',I2)
               IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
               GO TO 77
            ENDIF
            N_AT_TYPE(6) = N_AT_TYPE(6) - 2*N_AT_TYPE(7)              ! 10-21-96
            N_AT_TYPE(26) = N_AT_TYPE(7)
            N_AT_TYPE(7) = 0          ! #NO2's = # nitro N's
         ENDIF
C----count C(sp2)-linked nitro groups                                 ! 10-21-96
         IF (N_AT_TYPE(110) .NE. 0) THEN                              ! 10-21-96
            IF (N_AT_TYPE(110) .LE. 2*N_AT_TYPE(6)) THEN              ! 10-21-96
               N_AT_TYPE(111) = N_AT_TYPE(110)                        ! 10-21-96
               N_AT_TYPE(110) = 0                                     ! 10-21-96
               N_AT_TYPE(6) = N_AT_TYPE(6) - 2*N_AT_TYPE(111)         ! 10-21-96
            ELSE                                                      ! 10-21-96
               WRITE (COMMENT,1740) N_AT_TYPE(110), N_AT_TYPE(6)      ! 10-21-96
               PRINT 740, N_AT_TYPE(110), N_AT_TYPE(6)                ! 10-21-96
1740           FORMAT (' $$$$ PROBLEMS: unbalanced C(sp2)-linked ',   ! 10-21-96
     x                 'NO2 atoms, # N''s =',I2,', # O''s =',I2)      ! 10-21-96
               IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                ! 10-21-96
               GO TO 77                                               ! 10-21-96
            ENDIF                                                     ! 10-21-96
         ENDIF                                                        ! 10-21-96
C----count nitro groups....cubane-linked
         IF (N_AT_TYPE(37) .NE. 0) THEN
            IF (N_AT_TYPE(11) .NE. 2*N_AT_TYPE(37)) THEN
               WRITE (COMMENT, 742) N_AT_TYPE(37), N_AT_TYPE(11)
               PRINT 742, N_AT_TYPE(37), N_AT_TYPE(11)
742            FORMAT (' $$$$ PROBLEMS: unbalanced cubane NO2 atoms, ',
     x                 '# N''s =',I2,', # O''s =',I2)
               IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
               GO TO 77
            ENDIF
            N_AT_TYPE(11) = 0
            N_AT_TYPE(38) = N_AT_TYPE(37)
            N_AT_TYPE(37) = 0          ! #NO2's = # nitro N's
         ENDIF
C----count nitramines
         IF (N_AT_TYPE(8) .NE. 0) THEN
            N_AT_TYPE(28) = N_AT_TYPE(8)
            N_AT_TYPE(26) = N_AT_TYPE(26) - N_AT_TYPE(28)
            N_AT_TYPE(8) = 0
         ENDIF
         IF (N_AT_TYPE(27) .NE. 0) THEN
            N_AT_TYPE(29) = N_AT_TYPE(27)
            N_AT_TYPE(26) = N_AT_TYPE(26) - N_AT_TYPE(29)
            N_AT_TYPE(27) = 0
         ENDIF
C----count carboxylic acid groups...atom type 40 is the number of CO2H H's
         IF (N_AT_TYPE(40) .GE. 1) THEN
            N_AT_TYPE(41) = N_AT_TYPE(40)  ! 40 = # CO2H H's
            N_AT_TYPE(15) = N_AT_TYPE(15) - N_AT_TYPE(41)  ! 41 = # CO2H groups
            N_AT_TYPE(31) = N_AT_TYPE(31) - N_AT_TYPE(41)
            N_AT_TYPE(32) = N_AT_TYPE(32) - N_AT_TYPE(41)
            N_AT_TYPE(40) = 0
         ENDIF
C----Count number of carbonyl groups...C=O; done after counting CO2H's;
C     # O's must equal # not in ring (32) + # in ring (96) C's
         IF (N_AT_TYPE(15) .EQ. (N_AT_TYPE(32) +
     x                           N_AT_TYPE(96))) THEN
            N_AT_TYPE(73) = N_AT_TYPE(32)  ! C=O not in ring
            N_AT_TYPE(97) = N_AT_TYPE(96)  ! C=O in ring
            N_AT_TYPE(32) = 0
            N_AT_TYPE(96) = 0
            N_AT_TYPE(15) = 0
         ELSE
            WRITE (COMMENT, 340) N_AT_TYPE(32), N_AT_TYPE(96),
     x                           N_AT_TYPE(15)
            PRINT 340, N_AT_TYPE(32), N_AT_TYPE(96), N_AT_TYPE(15)
340         FORMAT (' $$$$ PROBLEMS: unbalanced C=O atoms, ',
     x              'C''s not in ring =',I2,', C''s in ring =',I2,
     x              ', O''s =',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C----Count number of amide, imide and urea C=O groups
         IF (N_AT_TYPE(79) .EQ. N_AT_TYPE(80)) THEN  ! C & O for the carbonyl
            N_AT_TYPE(78) = N_AT_TYPE(79)   ! amide C=O = # amide C's
            N_AT_TYPE(79) = 0          ! # C's = 0
            N_AT_TYPE(80) = 0          ! # O's = 0
         ELSE
            WRITE (COMMENT, 346) N_AT_TYPE(79), N_AT_TYPE(80)
            PRINT 346, N_AT_TYPE(79), N_AT_TYPE(80)
346         FORMAT (' $$$$ PROBLEMS: unbalanced amide, urea, imide ',
     x              'C=O atoms, # C''s =',I2,', # O''s =',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C----Count number of O-H groups linked to Csp3
         IF (N_AT_TYPE(30) .EQ. N_AT_TYPE(31)) THEN
            IF (N_AT_TYPE(30) .LE. 0) GO TO 1342                    ! 10/24/96
            N_AT_TYPE(74) = N_AT_TYPE(30)   ! # O-H's (74) = # H's (30)
            N_AT_TYPE(31) = N_AT_TYPE(31) - N_AT_TYPE(30) ! correct # O's
            N_AT_TYPE(30) = 0      ! zero # H's
         ELSE
            WRITE (COMMENT, 341) N_AT_TYPE(31), N_AT_TYPE(30)
            PRINT 341, N_AT_TYPE(31), N_AT_TYPE(30)
341         FORMAT (' $$$$ PROBLEMS: unbalanced O-H atoms, # O''s =',I2,
     x              ', # H''s =',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C----Count number of O-H groups linked to Csp2
1342     IF (N_AT_TYPE(68) .EQ. N_AT_TYPE(67)) THEN                 ! 10/24/96
            N_AT_TYPE(75) = N_AT_TYPE(68)   ! # O-H's (75) = # H's (68)
            N_AT_TYPE(67) = N_AT_TYPE(67) - N_AT_TYPE(68) ! correct # O's
            N_AT_TYPE(68) = 0      ! zero # H's
         ELSE
            WRITE (COMMENT, 342) N_AT_TYPE(67), N_AT_TYPE(68)	
            PRINT 342, N_AT_TYPE(67), N_AT_TYPE(68)
342         FORMAT (' $$$$ PROBLEMS: unbalanced O-H atoms, # O''s =',I2,
     x              ', # H''s =',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C----count number of C-N(-)-N(+)(C)=N-C groups
         IF (N_AT_TYPE(83) .EQ. 2*N_AT_TYPE(82)) THEN
            N_AT_TYPE(84) = N_AT_TYPE(82)
            N_AT_TYPE(82) = 0
            N_AT_TYPE(83) = 0
         ELSE
            WRITE (COMMENT, 748) N_AT_TYPE(83), N_AT_TYPE(82)
            PRINT 748, N_AT_TYPE(83), N_AT_TYPE(82)
748         FORMAT (' $$$$ PROBLEMS: for the C-N(-)-N(+)(C)=N-C group;',
     x       ' # of end N''s =',I2,', # of middle N''s =',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C----count nitroso groups (-N=O)
         IF (N_AT_TYPE(61) .NE. N_AT_TYPE(62)) THEN                    ! 12/28/95
            WRITE (COMMENT, 749) N_AT_TYPE(62), N_AT_TYPE(61)
            PRINT 749, N_AT_TYPE(62), N_AT_TYPE(61)                    ! 12/28/95
749         FORMAT (' $$$$ PROBLEMS: unbalanced -N=O atoms, # N =',I2,  ! 12/28/95
     x                              ', # O =',I2)                      ! 12/28/95
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)                    ! 12/28/95
            GO TO 77                                                   ! 12/28/95
         ELSE                                                          ! 12/28/95
            N_AT_TYPE(63) = N_AT_TYPE(61)   ! # N=O = # O              ! 12/28/95
            N_AT_TYPE(61) = 0    ! zero O count                        ! 12/28/95
            N_AT_TYPE(62) = 0    ! zero N count                        ! 12/28/95
         ENDIF                                                         ! 12/28/95
C----count azoxy groups, N(+)-O(-), from -N(O)=N- & -N(O)=N(O)-
         IF (N_AT_TYPE(56) .NE. N_AT_TYPE(57)) THEN
            WRITE (COMMENT, 849) N_AT_TYPE(56), N_AT_TYPE(57)
            PRINT 849, N_AT_TYPE(56), N_AT_TYPE(57)
849         FORMAT (' $$$$ PROBLEMS: unbalanced N & O atoms, # N =',I2,
     x              ', # O =',I2,', from -N(O)=N- or -N(O)=N(O)-')
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ELSE
            N_AT_TYPE(90) = N_AT_TYPE(57)   ! # N=O = # N
            N_AT_TYPE(56) = 0    ! zero O count
            N_AT_TYPE(57) = 0    ! zero N count
         ENDIF
C----count pyridine N-oxide (C=N(+,C)-O(-)) N-O groups,
C     amine N-oxide (R3N(+)-O(-)) N-O groups &
C     amine oxide (R2N-O) groups
      ISUM = N_AT_TYPE(18) + N_AT_TYPE(70) + N_AT_TYPE(89)
      IF (N_AT_TYPE(17) .NE. ISUM) THEN
            WRITE (COMMENT, 949) N_AT_TYPE(18), N_AT_TYPE(70),
     x                           N_AT_TYPE(89), ISUM, N_AT_TYPE(17)
            PRINT 949, N_AT_TYPE(18), N_AT_TYPE(70), N_AT_TYPE(89),
     x                           ISUM, N_AT_TYPE(17)
949      FORMAT (' $$$$ PROBLEMS: N sum from 18 (',I1,') + 70 (',I1,')',
     x        ' + 89 (',I1,') of',I2,', does not equal # O (17) of',I2)
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
      ELSE
            N_AT_TYPE(91) = N_AT_TYPE(18)   ! pyridine N-oxide (C=N(+,C)-O(-)) N-O groups
            N_AT_TYPE(18) = 0
            N_AT_TYPE(92) = N_AT_TYPE(70)   ! amine N-oxide (R3N(+)-O(-)) N-O groups
            N_AT_TYPE(70) = 0
            N_AT_TYPE(93) = N_AT_TYPE(89)   ! amine oxide (R2N-O) groups
            N_AT_TYPE(89) = 0
            N_AT_TYPE (17) = 0              ! zero O sum
      ENDIF
C----Count number of -O-NO2 groups
      IF (N_AT_TYPE(60) .GT. 0) THEN       ! # of -O-             ! 10/18/96
         IF (N_AT_TYPE(26) .LT. N_AT_TYPE(60)) THEN               ! 10/18/96
            WRITE (COMMENT, 951) N_AT_TYPE(26), N_AT_TYPE(60)     ! 10/18/96
            PRINT 951, N_AT_TYPE(26), N_AT_TYPE(60)               ! 10/18/96
951      FORMAT (' $$$$ PROBLEMS:',I2,' normal NO2 groups <',I2,  ! 10/18/96
     x           ' -O- in -O-NO2')                                ! 10/18/96
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)               ! 10/18/96
            GO TO 77                                              ! 10/18/96
         ENDIF                                                    ! 10/18/96
         N_AT_TYPE(109) = N_AT_TYPE(60)   ! -O-NO2, param = 109   ! 10/18/96
         N_AT_TYPE(60) = 0                                        ! 10/18/96
         N_AT_TYPE(26) = N_AT_TYPE(26) - N_AT_TYPE(109)           ! 10/18/96
      ENDIF                                                       ! 10/18/96
C----Are separate non-ring CH, CH2 and CH3's being counted
      IF (JCSP3 .EQ. 0) GO TO 408
         N_HYDROGENS = N_AT_TYPE(144) + 2*N_AT_TYPE(145) +        ! 4/11/00
     x                 3*(N_AT_TYPE(146) + N_AT_TYPE(147) +       ! 4/11/00
     X                    N_AT_TYPE(148) + N_AT_TYPE(149))        ! 4/11/00
         IF (N_AT_TYPE(2) .GE. N_HYDROGENS) THEN
            N_AT_TYPE(2) = N_AT_TYPE(2) - N_HYDROGENS
C-------Now code CH, CH2 and various CH3 types                    ! 4/11/00
        DO I=1,3                                                  ! 4/11/00
           N_AT_TYPE(104+I) = N_AT_TYPE(143+I)  ! CH, CH2, CH3    ! 4/11/00
           N_AT_TYPE(143+I) = 0                                   ! 4/11/00
           N_AT_TYPE(133+I) = N_AT_TYPE(146+I)  ! CH3-N, CH3-O    ! 4/11/00
           N_AT_TYPE(146+I) = 0                 ! & CH3-S         ! 4/11/00
        ENDDO                                                     ! 4/11/00
         ELSE
            WRITE (COMMENT, 960) (N_AT_TYPE(I+104),I=1,3),        ! 3/24/00
     x                 (N_AT_TYPE(I+133),I=1,3),                  ! 3/24/00
     x                 N_HYDROGENS, N_AT_TYPE(2)
            PRINT 960, (N_AT_TYPE(I+143),I=1,3),                  ! 4/11/00
     x                 (N_AT_TYPE(I+146),I=1,3),                  ! 4/11/00
     x                 N_HYDROGENS, N_AT_TYPE(2)
960         FORMAT (' $$$$ PROBLEMS: # hydrogens associated with',! 3/24/00
     x       I3,' CH,',I3,' CH2,',I3,' CH3-C,'/                   ! 3/24/00
     x       15X,I3,' CH3-N,',I3,' CH3-O &',I3,' CH3-S =',I3,     ! 3/24/00
     x       '; # Csp3 H''s =',I3)                                ! 3/24/00
            IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
            GO TO 77
         ENDIF
C
C----Show which atom types have been determined
408   I1 = 0
      DO 415 I2=1,NPOTS
         IF (N_AT_TYPE(I2) .EQ. 0) GO TO 415
         I1 = I1 + 1
         IND(I1) = I2                  ! atom/group type
         NATTYPES(I1) = N_AT_TYPE(I2)  ! number of this atom/group type
415   CONTINUE
C
      WRITE (6,418) NATTYPES(1), IND(1),
     x              DESCRIPTION(IND(1))
418   FORMAT (' Atom types:',I3,' [#',I3,'] of ',A60)
      IF (I1 .EQ. 1) GO TO 7777
      DO I2=2,I1
         WRITE (6,918) NATTYPES(I2), IND(I2),
     x                 DESCRIPTION(IND(I2))
918      FORMAT ('            ',I3,' [#',I3,'] of ',A60)
      ENDDO
C
C----Calc volumes and densities from DATA_TABLE (table3) information
7777  SUMLV = 0.0                                            ! 12/16/96
      SUMNLV = 0.0                                           ! 12/16/96
      SUMDEN = 0.0
      SUMWT = 0.0                                            ! 12/16/96
      DO 4100 I=1,I1                                         ! 12/16/96
         L = IND(I)
           SUMLV = SUMLV + DATA_TABLE(L,1)*NATTYPES(I)              ! 12/16/96
           SUMNLV = SUMNLV + DATA_TABLE(L,2)*NATTYPES(I)            ! 12/16/96
           SUMDEN = SUMDEN + DATA_TABLE(L,3)*NATTYPES(I)
           SUMWT = SUMWT + DATA_TABLE(L,4)*NATTYPES(I)              ! 12/16/96
4100  CONTINUE                                                      ! 12/16/96
      DL = SUMWT/(.6022*SUMLV)                                      ! 12/16/96
      DN = SUMWT/(.6022*SUMNLV)                                     ! 12/16/96
      SUMDEN = (1.0/SUMWT)*SUMDEN
      WRITE (6,4120) DL, DN, SUMDEN
4120  FORMAT (' Calcd densities: from linear volumes =',F6.3,' g/cc'/
     x        '               from nonlinear volumes =',F6.3,' g/cc'/
     x        '                from linear densities =',F6.3,' g/cc')
C
C----Finished
404   IF (ICC1 .EQ. 0) CLOSE (UNIT=11)
      STOP
77    CLOSE (UNIT = 12)
1100  WRITE (6,1105) REFCOD, I, NEW_LABEL(I)
      WRITE (COMMENT, 1105) REFCOD, I, NEW_LABEL(I)
1105  FORMAT (' ---> PROBLEMS with atom types for ',A8,
     x        ', atom # =',I3,' (',A5,')')
      IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
      GO TO 404
      END
C--------------------------------------------------------
      FUNCTION SIND(X)
      SIND = SIN(X/57.29577951)
      RETURN
      END
C--------------------------------------------------------
      FUNCTION COSD(X)
      COSD = COS(X/57.29577951)
      RETURN
      END
C--------------------------------------------------------
      SUBROUTINE CONNECTION(ITEST,ICOORD)               ! 7/15/05 		
C
C----Determine connectivity
C
      CHARACTER LABEL*5, NAME(500)*5
      COMMON /COM1/ NCON(500), ICON(5,500), LABEL(500),
     x                 OR_XYZ(3,500), CELL(6), NATOMS,
     x       IT_TABLE(500), I3D(500), KEEP_NPOT(500)

      EQUIVALENCE (LABEL, NAME)
      DIMENSION D_H(2), N_H(2), dlist(8)
C
      ITEST = 0
C
C----Determine connectivity....use max distance of DMAX for most
      DO 502 I=1,NATOMS
         NCON(I) = 0
         DO 490 J=1,NATOMS
            DMAX = 1.65
            IF (I .EQ. J) GO TO 490
            IF (LABEL(I)(1:1) .EQ. 'H' .AND.
     x          LABEL(J)(1:1) .EQ. 'H') GO TO 490
            D = DISTANCE (I, J)
            IF (LABEL(I)(1:1) .EQ. 'S' .OR.          ! S          ! 3/13/00
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 1.95               ! 3/28/00
            IF (LABEL(I)(1:1) .EQ. 'P' .OR.          ! P          ! 4/5/00
     X          LABEL(J)(1:1) .EQ. 'P') DMAX = 1.90               ! 4/5/00
            IF (LABEL(I)(1:1) .EQ. 'S' .AND.         ! S-S        ! 3/28/00
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.15               ! 3/28/00
            IF (LABEL(I)(1:1) .EQ. 'X' .OR.          ! X = BR     ! 2/26/00
     X          LABEL(J)(1:1) .EQ. 'X') DMAX = 2.04  ! X = BR     ! 2/26/00
C
            IF (LABEL(I)(1:1) .EQ. 'C' .AND.         ! C-C = 1.78 ! 7/14/05 
     X          LABEL(J)(1:1) .EQ. 'C') DMAX = 1.78               ! 7/14/05
C
            IF (LABEL(I)(1:1) .EQ. 'Z' .OR.          ! Z = CL     ! 3/4/00
     X          LABEL(J)(1:1) .EQ. 'Z') DMAX = 1.90  ! Z = CL     ! 3/15/00
            IF (LABEL(I)(1:1) .EQ. 'I' .OR.          ! Z = I      ! 6/7/05
     X          LABEL(J)(1:1) .EQ. 'I') DMAX = 2.25  ! Z = I      ! 6/7/05
            IF (LABEL(I)(1:1) .EQ. 'H' .OR.
     X          LABEL(J)(1:1) .EQ. 'H') DMAX = 1.30
            IF (D .GT. DMAX) GO TO 490
            NCON(I) = NCON(I) + 1    ! no more than 4 connections
            ICON(NCON(I),I) = J
            IF (NCON(I) .LE. 4) GO TO 490
               ITEST = I
               GO TO 700
490      CONTINUE
502   CONTINUE
C----Try to correct any H's with NCON > 2...quit if NCON = 0
      DO 600 I=1,NATOMS     ! look only at H's
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 600
         IF (NCON(I) .EQ. 1) GO TO 600
         IF (NCON(I) .LT. 1) THEN   ! can't deal with 0's
            ITEST = I        ! quit
            GO TO 700
         ENDIF
	 L = 0               ! eliminate any H...H connection first
         DO 510 J=NCON(I),1,-1
            K = ICON(J,I)
            IF (NAME(K)(1:1) .NE. 'H') GO TO 510
            L = L + 1   ! # of H's found & eliminated
            ICON(J,I) = 0
510      CONTINUE
         NCON(I) = NCON(I) - L
C
         IF (NCON(I) .GT. 2) THEN          ! # > 2, bail out
            ITEST = I
            GO TO 700
         ENDIF
C
         DO 560 J=1,NCON(I)    ! only handling cases with 2 links
            K = ICON(J,I)
            N_H(J) = K
            D_H(J) = DISTANCE(K,I)   ! H...? distance
560      CONTINUE
         LARGEST = N_H(1)      ! keep the ? with smallest distance to H
         ISMALLEST = N_H(2)
         IF (D_H(1) .LE. D_H(2)) THEN
            LARGEST = N_H(2)
            ISMALLEST = N_H(1)
         ENDIF    ! LARGEST = # of atom with H contact to be eliminated
C----Fix up the connectivity list
         NCON(I) = 1
         ICON(1,I) = ISMALLEST   ! so much for the H with smallest dist to ?
C
         N = NCON(LARGEST)        ! LARGEST = atom with extra link to H
         DO 580 K=1,N
            IF (ICON(K,LARGEST) .NE. I) GO TO 580
            IF (K .EQ. N) THEN
               ICON(K,LARGEST) = 0
            ELSE
               DO 585 M=K,N-1
                  ICON(M,LARGEST) = ICON(M+1,LARGEST)
585            CONTINUE
            ENDIF
            NCON(LARGEST) = N - 1
            GO TO 600
580     CONTINUE
600   CONTINUE
      IF (ICOORD .EQ. 2) GO TO 701
      RETURN
C
700   if (ICOORD .eq. 0) return
701   print *, ' +++++++++++ Connection Table +++++++++++'
       if (itest .ne. 0) then
          print 1265, itest
1265      format (' partial connectivity list through atom #',i3)
          n = itest
       else
          n = natoms
       endif
       do i=1,n
         do j=1,ncon(i)
            m = icon(j,i)
            dlist(j) = distance(i,m)
         enddo
         print 1251, i, label(i), (icon(j,i),dlist(j),j=1,ncon(i))
1251     format (i3,2x,a5,5(i4,' /',f6.3))
       enddo
      RETURN
      END
C
C------------------------------------------------------------
C
C----Determine atom types for volume additivity calculations
C
      FUNCTION NPOT (NUM_ATOM)
      CHARACTER NAME*5, YN*1, MMOD(64)*2, AT2*2, LINE*60
      CHARACTER COMMENT*120, REFCOD_2*8
      DIMENSION IK(500), NSIZES(10), NUM_O(3), DB(3)               ! 10-23-96
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500),KEEP_NPOT(500)
      COMMON/T1/IAT_USED(500),NAROMATIC
      COMMON /GROUPS/ ICOORD
      COMMON /PROBLEMS/ JPROB, COMMENT, REFCOD_2, ICSP3, JCSP3
C
10    AT2 = NAME(NUM_ATOM)(1:2)
      NPOT = 0
c     print 2211, at2, num_atom, it_table(num_atom)   ! $$$$$$$$$$$$$$$$$$
c2211  format ('top of npot: at2,num_atom,it_table =',a5,2i4) ! $$$$$$$$$$$$$$$$$
C----Perhaps this atom is already known?                      ! 3/20/00
      IF (IT_TABLE(NUM_ATOM) .NE. 0) THEN                     ! 3/20/00
         NPOT = KEEP_NPOT(NUM_ATOM)                           ! 3/20/00
         RETURN                                               ! 3/20/00
      ENDIF                                                   ! 3/20/00
C
C----look at H's first....................................
C
1375    IF (AT2(1:1) .NE. 'H') GO TO 1340              ! H
C
           NPOT = IQ_CUBANE_H(NUM_ATOM)     ! cubane-linked H?
              IF (NPOT .NE. 0) GO TO 1303   ! cubane-H, param = 36
C
           IT = ICON(1,NUM_ATOM)     ! number of atom to which H is attached
C
           IF (NAME(IT)(1:1) .EQ. 'O') THEN   ! H-O?
              NPOT = 30                       ! hydroxyl H, param = 30
              GO TO 1303
           ENDIF
C-----Check for thiol H                                   ! 3/22/00
           IF (NAME(IT)(1:1) .EQ. 'S') THEN               ! 3/22/00
              NPOT = 130      ! H of S-H, param = 130     ! 3/22/00
              GO TO 1303                                  ! 3/22/00
           ENDIF                                          ! 3/22/00
C
           IF (NAME(IT)(1:1) .EQ. 'C') THEN   ! H-C
              NC = NCON(IT)                   ! how many atoms bonded to C?
              GO TO (101, 108, 103, 104), NC
C
108           NPOT = 54           ! H bonded to sp C; param = 54
                 GO TO 1303
C
C----H bonded to 3-linked C...is it an aomatic C?
103           CALL BOND_CHECK(IT, N_INRINGS, NSIZES, NBONDS)
              GO TO (1000,3005,9002,3005,3005), NBONDS + 1
C                     0    ?    2    ?    ?    <--- # bonds in rings
1000          NPOT = 1         ! Csp2 not in any rings, H param = 1
                GO TO 1303
C----Csp2 with 2 bonds in rings...aromatic or non-aromatic?
9002          IF (NEG_COUNT(NSIZES) .EQ. 1) THEN  ! how many -5's and -6's?
                NPOT = 108   ! aromatic C with 2 bonds in rings, H param = 108
              ELSE
                NPOT = 1     ! 2 bonds in rings...C not aromatic, H param = 1
              ENDIF
                GO TO 1303
C
104           NPOT = 2            ! H bonded to normal sp3 C; param = 2
                 GO TO 1303
           ELSE
              IF (NAME(IT)(1:1) .NE. 'N') GO TO 101   ! H-N ?
                 NC = NCON(IT)                      ! H is bonded to N
                 IF (NC .NE. 3 .OR. NC .NE. 2) THEN  ! R2NH, RNH2 or =N-H
                    NPOT = 3         ! H of R2NH, RNH2 or =N-H; param = 3
                    GO TO 1303
                 ENDIF
                 GO TO 101
           ENDIF
           GO TO 1303
101        PRINT 110, NAME(NUM_ATOM), NUM_ATOM, NAME(IT), IT
110        FORMAT (1X,A5,' (#',I2,') bonded to ',A5,' (#',
     x             I2,')...unknown type')
           GO TO 2002
C
C----Look at C's second.......................................
C
1340    IF (AT2(1:1) .NE. 'C') GO TO 2340
           GO TO (1001, 1002, 1003, 1004, 1001), NCON(NUM_ATOM)
C
C---- 1-linked C...check for isocyanides
1001       NPOT = IQ_ISOCYANIDE(NUM_ATOM)
              IF (NPOT .NE. 0) GO TO 1303
              GO TO 1007
C
C---- 2-linked C...check for nitrile or alkyne C
1002       NPOT = IQ_CN(NUM_ATOM,1)  ! params = 46 & 45 (N#C); 53 (C#C)
              IF (NPOT .NE. 0) GO TO 1303    ! 55 (C#C in a ring)
C
C----check for central allene-type C...like sp
           NPOT = IQ_ALLENE_2(NUM_ATOM,1) ! central allene C, C=C=C ?
              IF (NPOT .NE. 0) GO TO 1303    ! PARAM = 65
C
C----can't handle the problem
1007          PRINT 800, NAME(NUM_ATOM), NUM_ATOM, NCON(NUM_ATOM),
     x                   (ICON(LL,NUM_ATOM),LL=1,NCON(NUM_ATOM))
              WRITE (COMMENT, 800) NAME(NUM_ATOM), NUM_ATOM,
     X    NCON(NUM_ATOM), (ICON(LL,NUM_ATOM),LL=1,NCON(NUM_ATOM))
800           FORMAT (' $$$$ PROBLEM: ',A5,' (atom #',I3,') has',I2,
     x            ' connections...does not compute...icon =',10i4)
              IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
              GO TO 2002
C
C---- 3-linked is it -CO2H group ?
C     x            num_atom, name(num_atom)
1003  NPOT = IQ_CO2H(NUM_ATOM)
         IF (NPOT .NE. 0) GO TO 1303    ! carbonyl C param = 32
C
C----Check for carbonyl C of C=S; C can be connected to           ! 3/28/00
C     C, S, O, N                                                  ! 3/28/00
      NS_1 = 0                   ! look for 1-linked S            ! 3/28/00
      DO 2015 I=1,3                                               ! 3/11/00
         IT = ICON(I,NUM_ATOM)                                    ! 3/11/00
         IF (NAME(IT)(1:1) .EQ. 'S' .AND. NCON(IT) .EQ. 1)        ! 3/11/00
     X       NS_1 = NS_1 + 1                                      ! 3/11/00
2015  CONTINUE                                                    ! 3/11/00
      IF (NS_1 .EQ. 1) THEN                                       ! 3/28/00
         NPOT = 115              ! param = 115 for C of C=S       ! 3/11/00
         GO TO 1303                                               ! 3/11/00
      ENDIF                                                       ! 3/11/00
C
C----Check for carbonyl C
      NPOT = IQ_CARBONYL(NUM_ATOM,1)
         IF (NPOT .GT. 0) GO TO 1303    ! carbonyl C; param = 32, 96 or 79
C
      NPOT = IQ_ALLENE_END(NUM_ATOM)    ! terminal C=C=C
         IF (NPOT .NE. 0) go to 1303    ! param = 64
C
C---- 3-linked C, must be sp2
      CALL BOND_CHECK(NUM_ATOM, N_INRINGS, NSIZES, NBONDS)
      GO TO (2020,3005,2022,2023,3005), NBONDS + 1
C             0    ?    2    3    ?    <--- # bonds in rings
2020  NPOT = 4                           ! Csp2 not in any rings, param = 4
      GO TO 1303
C----Csp2 with 2 bonds in rings...aromatic or non-aromatic?
2022     IF (NEG_COUNT(NSIZES) .EQ. 1) THEN  ! how many -5's and -6's?
            NPOT = 33        ! aromatic with 2 bonds in rings, param = 33
            GO TO 1303                                                    ! 10-23-96
         ELSE
C----------2 bonds in a non-aromatic ring.  Is C=C in the ring?           ! 10-23-96
            NPOT = 24            ! initialize to C=C not in a ring        ! 10-23-96
            DO 7002 I=1,3                                                 ! 10-23-96
               J = ICON(I,NUM_ATOM)                                       ! 10-23-96
               IF (NAME(J)(1:1) .EQ. 'H') THEN                            ! 10-23-96
                  NPOT = 112                      ! C=C in the ring       ! 10-23-96
                  GO TO 1303                                              ! 10-23-96
               ENDIF                                                      ! 10-23-96
               DB(I) = DISTANCE(NUM_ATOM, J)                              ! 10-23-96
7002        CONTINUE                                                      ! 10-23-96
C---------Find shortest bond to the 3 atoms connected to NUM_ATOM         ! 10-23-96
            ISHORT = 1                                                    ! 10-23-96
            DO 7004 K=2,3                                                 ! 10-23-96
               IF (DB(ISHORT) .LE. DB(K)) GO TO 7004                      ! 10-23-96
               ISHORT = K                                                 ! 10-23-96
7004        CONTINUE                                                      ! 10-23-96
            ISHORT = ICON(ISHORT,NUM_ATOM)                                ! 10-23-96
C-------------Is the atom (ISHORT) linked to NUM_ATOM by the shortest     ! 10-23-96
C              bond length in a ring?                                     ! 10-23-96
              CALL BOND_CHECK(ISHORT, N_INRINGS, NSIZES, NBONDS)          ! 10-23-96
              GO TO (7000,7000,7022,7022,7022), NBONDS + 1                ! 10-23-96
C                      0    ?    2    3    ?    <--- # bonds in rings     ! 10-23-96
7022          NPOT = 112                ! C=C in the ring                 ! 10-23-96
              GO TO 1303                                                  ! 10-23-96
         ENDIF
7000     GO TO 1303                                                       ! 10-23-96
C----Csp2 with 3 bonds in rings...all aromtaic ?
2023     GO TO (5010, 5012, 5013), NEG_COUNT(NSIZES) + 1
C                 0     2     3  <--- # aromatic bonds
5010         NPOT = 25        ! Csp2 with 3 bonds in rings, param = 25
           GO TO 1303
5012         NPOT = 34        ! 3 bonds in rings, 2 aromatic, param = 34
           GO TO 1303
5013         NPOT = 35        ! 3 bonds in rings, 3 aromatic, param = 35
           GO TO 1303
C
C----Problem section...........................................
C
3005    WRITE(6,2010) NUM_ATOM, NAME(NUM_ATOM), NBONDS,
     x                NEG_COUNT(NSIZES), (NSIZES(L),L=1,N_INRINGS)
        WRITE (COMMENT, 2010) NUM_ATOM, NAME(NUM_ATOM), NBONDS,
     x                NEG_COUNT(NSIZES), (NSIZES(L),L=1,N_INRINGS)
2010    FORMAT (' $$$$ PROBLEMS: atom number =',I3,', name = ',A5,
     x          ', number of bonds in rings =',I2,', neg_count =',I3,
     x          ', ring sizes =',4I3)
        IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
        NPOT = 999
        GO TO 1303
C
C---- 4-linked C..................
C
C----Is Csp3 a cubane C?
1004  IF (IQ_CUBANE(NUM_ATOM) .NE. 0) THEN
        NPOT = 5                         ! cubane C; param = 5
        GO TO 1303
      ELSE
C----Is it part of a ring?
      CALL BOND_CHECK(NUM_ATOM, NRINGS, NSIZES, NBONDS)
        GO TO (3000,3005,3002,2003,3004), NBONDS + 1
C                0    ?    2    3    4     <--- # bonds in rings
C------Should the 2 bonds in a ring Csp3 be grouped by rings?
3002    IF (ICSP3 .EQ. 0) THEN
           NPOT = 20                 ! no....Csp3 with 2 bonds in a ring
        ELSE
           IF (NSIZES(1) .GE. 9) THEN
              NPOT = 104             ! 9-ring or larger
           ELSE
              NPOT = 95 + NSIZES(1)  ! 98, 99, 100, 101, 102, 103 for
           ENDIF                     !  3   4    5    6    7    8-membered rings
        ENDIF
        GO TO 1303
2003    NPOT = 21                          ! Csp3 with 3 bonds in rings
        GO TO 1303
3004    NPOT = 22                          ! Csp3 with 4 bonds in rings
        GO TO 1303
C------Should separate increments for non-ring CH3, CH2 and CH be calculated
3000    IF (JCSP3 .EQ. 0) THEN
           NPOT = 23                          ! Csp3 that's not in a ring
        ELSE
           NPOT = IQ_CH3CH2CH(NUM_ATOM)       ! make them separate, CH3, CH2
        ENDIF                                 !  and CH not in rings
        GO TO 1303
      ENDIF
C
C----Look at O's third.................................................
C
2340    IF (AT2(1:1) .NE. 'O') GO TO 1024
           GO TO (1021, 1022), NCON(NUM_ATOM)
C
C---- 1-linked O...............
C
1021          IT = ICON(1,NUM_ATOM)     ! get ID of single connected atom
C
C------------Check for O to S
              IF (NAME(IT)(1:1) .EQ. 'S' .AND. NCON(IT) .GE. 3) THEN   ! 3/18/00
                 NPOT = IQ_S_O(NUM_ATOM,IT)   ! #'s of O and S         ! 3/18/00
                 IF (NPOT .NE. 0) GO TO 1303                           ! 3/18/00
              ENDIF                                                    ! 3/18/00
C
              IF (NAME(IT)(1:1) .EQ. 'N' .AND. NCON(IT) .EQ. 2) THEN
                 NPOT = IQ_NITROSO(NUM_ATOM)   ! -N=O O ?
                 IF (NPOT .NE. 0) GO TO 1303
              ENDIF
C
              IF (NAME(IT)(1:1) .EQ. 'N' .AND. NCON(IT) .EQ. 3) THEN
                 NPOT = IQ_AZOXY(NUM_ATOM)      ! -N(O)=N- or -N(O)=N(O)- O ?
                 IF (NPOT .NE. 0) GO TO 1303
              ENDIF
C
C-------------R3N(+)-O(-); pyridine-N(+)-O(-); R2N-O; C=N(C,+)-O(-)
C     params    N = 70       N = 18            N = 89    N = 18
C             param for O = 17
              NPOT = IQ_AMINEOXIDE(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
C
              IF (NAME(IT)(1:1) .EQ. 'N' .AND. NCON(IT) .EQ. 3) THEN
                 NPOT = IQ_AZOXY(NUM_ATOM)      ! -N(O)=N- or -N(O)=N(O)- O ?
                 IF (NPOT .NE. 0) GO TO 1303
              ENDIF
C
              IF (NAME(IT)(1:1) .NE. 'C') GO TO 1019
              NPOT = KEEP_NPOT(NUM_ATOM)
              IF (IT_TABLE(NUM_ATOM) .EQ. 1 .AND.
     X            NPOT .NE. 0) GO TO 1303
              GO TO 101   ! problems...C-O, O should already be known
C
1019          IF (NAME(IT)(1:1) .EQ. 'N' .AND.        ! N is linked to the O
     x               NCON(IT) .EQ. 3) THEN            ! & has 3 connections
                 N_OXY = 0
                 N_CARBS = 0
                 N_1 = 0
                 DO I=1,3                   ! count # of O's & C's linked
                    ITT = ICON(I,IT)        ! to the 2nd atom
                    IF (NAME(ITT)(1:1) .EQ. 'C') N_CARBS = N_CARBS + 1
                    IF (NAME(ITT)(1:1) .EQ. 'O') THEN
                        N_OXY = N_OXY + 1
                        IF (NCON(ITT) .EQ. 1) N_1 = N_1 + 1  ! # of 1-linked O's
                    ENDIF
                 ENDDO
                 IF (N_1 .EQ. 2) THEN    ! probably a nitro O (C-NO2 or O-NO2)
                    NPOT = IQ_NITRO_OX(IT)    ! param = 6 for normal NO2 O atom
                    IF (NPOT .NE. 0) GO TO 1303 ! can't identify the kind of O
                 ENDIF
C
                 IF (N_1 .EQ. 1) THEN
                    NPOT = IQ_PYRIDINE_OXIDE(NUM_ATOM)
                    IF (NPOT .NE. 0) GO TO 1303  ! pyridine or aza-pyrrole N-oxide O
                 ENDIF                   ! N param = 18; O param = 17
              ENDIF
              GO TO 2007
C
C---- 2-linked O.........................................................
C
1022          NPOT = IQ_OH_O(NUM_ATOM)         ! is it an hydroxyl O?
                 IF (NPOT .NE. 0) GO TO 1303   ! O of C-O-H, N-O-H, O-O-H; param = 31 or 67
C                                                Csp2 OH = 67; Csp3, N-OH, O-OH = 31
C
C------------Check for 2-linked O of a sulfonate                        ! 3/20/00
              DO I=1,2                                                  ! 3/20/00
                 IS = ICON(I,NUM_ATOM)                                  ! 3/20/00
                 IF (NAME(IS)(1:1) .EQ. 'S' .AND. NCON(IS) .EQ. 4) THEN ! 3/20/00
                     NPOT = IQ_S_O(NUM_ATOM,IS)                         ! 3/20/00
                     IF (NPOT .NE. 0) GO TO 1303                        ! 3/20/00
                 ENDIF                                                  ! 3/20/00
              ENDDO                                                     ! 3/20/00
C
              NPOT = IQ_EPOXIDE(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303   ! epoxide O, param = 76
C
              NPOT = IQ_ESTER(NUM_ATOM)        ! is it an ester or anhydride O?
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 94/95 for yes/no in a ring
C
              NPOT = IQ_ONO2(NUM_ATOM)         ! is it a nitrate ester O ?
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 60
C
              NPOT = IQ_PEROXY(NUM_ATOM)       ! is it a peroxide O ?
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 66; -O-O-
C
C----check for C(N)-O-C(N) where the C's aren't carbonyls
1023          IT1 = ICON(1,NUM_ATOM)
              IT2 = ICON(2,NUM_ATOM)
              IF ((NAME(IT1)(1:1) .EQ. 'C' .OR.
     x             NAME(IT1)(1:1) .EQ. 'N') .AND.
     X            (NAME(IT2)(1:1) .EQ. 'C' .OR.
     x             NAME(IT2)(1:1) .EQ. 'N')) THEN
                  NPOT = 14      ! "ether" O; C(N)-O-C(N)
                  CALL BOND_CHECK(NUM_ATOM, N_INRINGS, NSIZES, NBONDS)
                  GO TO (1025,1026,1027), NBONDS + 1
C                          0    ?    2           <--- # bonds in rings
1025                 NPOT = 14     ! ether O; C-O-C; not in a ring
                   GO TO 1303
1026                 NPOT = 999     ! unknown
                   GO TO 1303
1027               IF (NEG_COUNT(NSIZES) .NE. 0) THEN
                     NPOT = 52  ! furan O
                   ELSE
                     NPOT = 44     ! ether O; C-O-C in a ring
                   ENDIF
                   GO TO 1303
              ENDIF
C
C----Look at N's fourth................................................
C
1024    IF (AT2(1:1) .NE. 'N') GO TO 2002
c      print 2209, at2(1:2), num_atom, ncon(num_atom)   ! $$$$$$$$$$$$$$$$$$$$$$$
c2209   format ('# N''s: at2, num_atom, ncon ='a5,1x,2i3)  ! $$$$$$$$$$$
           GO TO (1031, 1032, 1033, 1038), NCON(NUM_ATOM)
C
C---- 1-linked N.................
C
1031          NPOT = IQ_AZIDE_N1(NUM_ATOM)    ! terminal azide N [N=N=N]
C                                             ! or diazo [C=N=N] ?
                 IF (NPOT .NE. 0) GO TO 1303  ! azide N, param = 39
C                                             ! diazo N, param = 69
              NPOT = 999      ! unknown 1-linked N, already have
                 GO TO 2007   ! handled nitriles
C
C---- 2-linked N....................
C
1032  NPOT = IQ_N_AZAPENTALENE(NUM_ATOM)  ! azapentalene N, param = 19
             IF (NPOT .NE. 0) GO TO 1303
c     print 2281, num_atom       ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c2281  format ('after iq_n_azapentalene: num_atom =',i2) ! $$$$$$$$$$$$$$
C
C----Check if pyridine or aza-pyrrole type N
      CALL BOND_CHECK(NUM_ATOM, N_INRINGS, NSIZES, NBONDS)
         IF (NBONDS .NE. 2) GO TO 1040
         IF (NEG_COUNT(NSIZES) .NE. 1) GO TO 1040
            NPOT = 16    ! pyridine or aza-pyrrole N; param = 16
            GO TO 1303
C
1040  NPOT = IQ_AZIDE_N2(NUM_ATOM)     ! middle or central azide N ?
         IF (NPOT .NE. 0) GO TO 1303   ! azide N, param = 39
C
      NPOT = IQ_NNN_END(NUM_ATOM)      ! end N of C-N(-)-N(+)=N-C   param = 83
         IF (NPOT .NE. 0) GO TO 1303   !                 |
C                                                        C
      NPOT = IQ_AZO_N(NUM_ATOM)        ! diazo N, param = 51
ccc        print 2277, npot                   ! $$$$$$$$$$$$$$$$$$$$$$$$$$$
ccc2277    format ('iq_azo_n return: npot =',i4)  ! $$$$$$$$$$$$$$$$$$$$$$$$$$
         IF (NPOT .GT. 0) GO TO 1303                           ! 10/17/00
C
C----Is it an sp2 N in C=N- ?
      NPOT = IQ_NDOUBLEC(NUM_ATOM)   ! param = 47 (no ring) or 50 (1 ring)
         IF (NPOT .NE. 0) GO TO 1303
         GO TO 2007
C
C---- 3-linked N........................
C
1033          NPOT = IQ_AMIDE_H(NUM_ATOM)  ! look for amides & imides ! 10/17/00
c             print 2222, npot                    ! $$$$$$$$$$$$$$$$$$$$$$$$$$
c2222          format ('->>> npot =',i3)           ! $$$$$$$$$$$$$$$$$$$$$$$$$
                 IF (NPOT .NE. 0) GO TO 1303
C
              NPOT = IQ_N_AZAPENTALENE(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
C
1035          NPOT = IQ_PYRROLE_N(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303   ! pyrrole N; param = 48 or 49
C
              NPOT = IQ_NNN_MIDDLE(NUM_ATOM)   ! end N of C-N(-)-N(+)=N-C  param = 82
                 IF (NPOT .NE. 0) GO TO 1303   !                 |
C                                                                C
C
              N_OXY = 0       ! # of O's attached
              N_CARBS = 0     ! # of C's attached
              N_NIT = 0       ! # of N's attached
              N_HYD = 0       ! # of H's attached
              N_F = 0         ! # of F's attached
              DO I=1,3
                 N_ADJ = ICON(I,NUM_ATOM)
                 IF (NAME(N_ADJ)(1:1) .EQ. 'O') THEN
                     N_OXY = N_OXY + 1
                     NUM_O(N_OXY) = N_ADJ
                 ENDIF
                 IF (NAME(N_ADJ)(1:1) .EQ. 'C') N_CARBS = N_CARBS + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'H') N_HYD = N_HYD + 1
                 IF (NAME(N_ADJ)(1:1) .EQ. 'N') THEN
                     N_NIT = N_NIT + 1
                     NUM_N = N_ADJ       ! keep N atom #
                 ENDIF
                 IF (NAME(N_ADJ)(1:1) .EQ. 'F') N_F = N_F + 1
              ENDDO
c       print 2288, n_oxy,n_carbs,n_nit,n_hyd,n_f  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$
c2288   format ('new place: n_oxy,n_carbs,n_nit,n_hyd,n_f =',5i3)   ! $$$$$$$$$$$$$$$$
              IF (N_OXY .EQ. 2 .OR. N_OXY .EQ. 3) THEN
                 NPOT = IQ_NITRO_N(NUM_ATOM)   ! cubane-linked nitro N?
                 IF (NPOT .NE. 0) GO TO 1303  ! cubane nitro N, param = 37 or
              ENDIF             ! other nitro N (C-NO2 or O-NO2) ...param = 7
C
              IF (N_CARBS .EQ. 2 .AND. N_NIT .EQ. 1) THEN   ! nitramine N ?
                 IF (IQ_ANY_NO2(NUM_N) .EQ. 0) GO TO 3506 ! is attached N a NO2 ?
                 CALL BOND_CHECK(NUM_ATOM,N_INRINGS,NSIZES,NBONDS)
                 GO TO (3500,3005,3502,3005), NBONDS + 1
C                        0    ?    2    ?    <--- # bond in rings
3502             NPOT = 27             ! nitramine N with 2 bonds in a ring...
                 GO TO 1303            ! param = 27
3500             NPOT = 8    ! nitramine N with no bonds in rings...param = 8
                 GO TO 1303
              ENDIF
C
C----Check out tertiary N's...C3N, C2-N-N, C2-N-O,C-N-N2, C-N(-O-,N),
C                             C-N(OR)2
C
3506          IF ((N_CARBS .EQ. 3) .OR.
     x           (N_CARBS .EQ. 2 .AND. N_OXY .EQ. 1) .OR.
     x           (N_CARBS .EQ. 1 .AND. N_OXY .EQ. 2) .OR.
     x           (N_CARBS .EQ. 1 .AND. N_NIT .EQ. 2) .OR.
     x           (N_CARBS .EQ. 2 .AND. N_NIT .EQ. 1) .OR.
     x           (N_CARBS .EQ. 1 .AND. N_NIT .EQ. 1 .AND. N_OXY .EQ. 1
     x               .AND. NCON(NUM_O(1)) .EQ. 2)) THEN
                 CALL BOND_CHECK(NUM_ATOM,N_INRINGS,NSIZES,NBONDS)
                 GO TO (4500,3005,4502,4504), NBONDS + 1
C                        0    ?    2    3    <--- # bond in rings
C----Is the tertiary N linked to a C(sp2)?
4500              DO 4501 I=1,3
                    J = ICON(I,NUM_ATOM)
                    IF (NAME(J)(1:1) .EQ. 'C' .AND. NCON(J) .EQ. 3) THEN
                       NPOT = 86     ! C(sp2)-linked tertiary N
                       GO TO 1303    ! param = 86
                    ENDIF
4501             CONTINUE
                 NPOT = 9       ! C3N, tertiary N not linked to C(sp2)
                    GO TO 1303  ! param = 9
4502             NPOT = 58  ! C3N, tert N in 1 ring; param = 58
                    GO TO 1303
4504             NPOT = 59  ! C3N, tert N in 2 rings; param = 59
                    GO TO 1303
              ENDIF
C----Is the N linked to a C(sp2)?
              IF (N_HYD .GE. 1 .OR. N_F .GE. 1) THEN
                 DO 4505 I=1,3
                    J = ICON(I,NUM_ATOM)
                    IF (NAME(J)(1:1) .EQ. 'C' .AND. NCON(J) .EQ. 3) THEN
                       NPOT = 85     ! C(sp2)-linked NH2, NH, NHF, NF2
                       GO TO 1303    ! param = 85
                    ENDIF
4505             CONTINUE
                 NPOT = 10      ! not C(sp2)-linked NH2, NH, NHF, NF2
                   GO TO 1303   ! param = 10
              ENDIF
              IF (N_CARBS .EQ. 2 .AND. N_OXY .EQ. 1) THEN
                 NPOT = 18      ! pyridine N-oxide N...param = 18
                 GO TO 1303
              ENDIF
       GO TO 2007
C
C---- 4-linked N...must be an amine N-oxide or other 4-lined N  ! 6/9/05
C
1038  NO = 0    ! count # of O's; must = 1
      DO I=1,4
         J = ICON(I,NUM_ATOM)
         IF (NAME(J)(1:1) .EQ. 'O') NO = NO + 1                ! 3/5/00
      ENDDO
      IF (NO .EQ. 1) THEN                                      ! 6/9/05
         NPOT = 70         ! amine N-oxide N...param = 70      ! 6/9/05
      ELSE                                                     ! 6/9/05
         NPOT = 158        ! 4-linked N+, not amine oxide      ! 6/9/05
      ENDIF                ! param = 158                       ! 6/9/05
         GO TO 1303                                            ! 6/9/05

C
C----Look at (bromine) Br's 5th..............................
2002   IF (AT2(1:1) .NE. 'X') GO TO 2004     ! X = BR          ! 3/4/00
          if (ncon(num_atom) .ne. 0) go to 7003                ! 5/22/05
             npot = 165              ! BR-                       5/22/05
             go to 1303                                        ! 5/22/05
7003       J = ICON(1,NUM_ATOM)                                 ! 3/17/00
           IF (.NOT. NAME(J)(1:1) .EQ. 'C' .AND. J .GE. 3)     ! 3/17/00
     x             GO TO 2007                                  ! 3/17/00
             IF (J .EQ. 4) THEN                                ! 3/17/00
                NPOT = 12            ! BR-C(sp3), param = 12   ! 3/17/00
             ELSE
                NPOT = 118           ! BR-C(sp2), param = 118  ! 3/17/00
             ENDIF                                             ! 3/17/00
             GO TO 1303
C
C----Look at (chlorine) Cl's 6th............................   ! 3/4/00
2004    IF (AT2(1:1) .NE. 'Z') GO TO 2006    ! Z = CL          ! 3/4/00
           if (ncon(num_atom) .ne. 0) go to 7005               ! 5/22/05
             npot = 164              ! CL-                       5/22/05
             go to 1303                                        ! 5/22/05
C-------Only accept CL to C(sp3) & C(sp2)                      ! 3/5/00
7005       J = ICON(1,NUM_ATOM)                                ! 3/5/00
           IF (NAME(J)(1:1) .NE. 'C') GO TO 2007               ! 3/5/00
              IF (NCON(J) .EQ. 4) THEN                         ! 3/5/00
                 NPOT = 113         ! CL-C(sp3) param = 113    ! 3/5/00
                 GO TO  1303                                   ! 3/5/00
              ENDIF                                            ! 3/5/00
              IF (NCON(J) .EQ. 3) THEN                         ! 3/5/00
                 NPOT = 114         ! CL-C(sp2) param = 114    ! 3/5/00
                 GO TO 1303                                    ! 3/5/00
              ENDIF                                            ! 3/5/00
C--------Discard a few CL-C(sp) structures                     ! 3/5/00
         GO TO 2007                                            ! 3/5/00
C
C----Look at (iodine) I's                                      ! 2/7/04
2006    IF (AT2(1:1) .NE. 'I') GO TO 2008                      ! 2/7/04
           if (ncon(num_atom) .ne. 0) go to 7007               ! 5/22/05
              npot = 166              ! I-                       5/22/05
              go to 1303                                       ! 5/22/05
7007    IF (NCON(NUM_ATOM) .GT. 1) GO TO 2007  ! only 1-linked I  2/7/04
C-------Only accept I to C(sp3) & C(sp2)                       ! 2/7/04
           J = ICON(1,NUM_ATOM)                                ! 3/5/00
           IF (NAME(J)(1:1) .NE. 'C') GO TO 2007               ! 3/5/00
              IF (NCON(J) .EQ. 4) THEN                         ! 3/5/00
                 NPOT = 160         ! I-C(sp3) param = 160     ! 2/7/04
                 GO TO  1303                                   ! 3/5/00
              ENDIF                                            ! 3/5/00
              IF (NCON(J) .EQ. 3) THEN                         ! 3/5/00
                 NPOT = 159         ! I-C(sp2) param = 159     ! 2/7/04
                 GO TO 1303                                    ! 3/5/00
              ENDIF                                            ! 3/5/00
C--------Discard any I-C(sp) structures                        ! 2/7/04
         GO TO 2007                                            ! 3/5/00
C
C----Look at F's 7th........................................
2008   IF (AT2(1:1) .NE. 'F') GO TO 3007                       ! 3/11/00
          I1 = ICON(1,NUM_ATOM)
          IF (NAME(I1)(1:1) .EQ. 'N') THEN
             NPOT = 13    ! N-F, param = 13
             GO TO 1303
          ENDIF
                         IF (NAME(I1)(1:1) .NE. 'C') GO TO 2007    ! F-???
             NPOT = 0
                               IF (NCON(I1) .EQ. 3) NPOT = 87         ! F-C(sp2)
          IF (NCON(I1) .EQ. 4) NPOT = 88         ! F-C(sp3)
          IF (NPOT .EQ. 0) GO TO 2007
             GO TO 1303
C
C----Look at S's 8th........................................    ! 3/11/00
C
C----Check for S of C=S                                         ! 3/20/00
3007    IF (AT2(1:1) .NE. 'S') GO TO 3030                       ! 4/4/00
        IF (AT2(1:1) .EQ. 'S' .AND. NCON(NUM_ATOM) .EQ. 1) THEN ! 4/4/00
           NPOT = 116      ! S of thio ald or ketone            ! 3/11/00
           GO TO 1303      ! param = 116                        ! 3/11/00
        ENDIF                                                   ! 3/11/00
C
C----Check for S of sulfoxide, sulfone or sulfonate             ! 3/20/00
        IF (AT2(1:1) .EQ. 'S' .AND. NCON(NUM_ATOM) .GE. 3) THEN ! 3/20/00
           NPOT = IQ_S_O(NUM_ATOM,NUM_ATOM)                     ! 3/20/00
           GO TO 1303                                           ! 3/20/00
        ENDIF                                                   ! 3/20/00
C
C----Check for 2-linked S...thiophene, sulfides, thiols and     ! 3/27/00
C     central S of C-S-S and S-S-S (S's must be 2-linked).      ! 3/27/00
C     Note that thiol H has been accounted for in H section.    ! 3/27/00
        IF (AT2(1:1) .EQ. 'S' .AND. NCON(NUM_ATOM) .EQ. 2) THEN ! 3/22/00
           NC = 0                      ! # C's linked           ! 3/27/00
           NS_2 = 0                    ! # 2-linked S's         ! 3/27/00
           DO I=1,2                                             ! 3/22/00
              IT = ICON(I,NUM_ATOM)                             ! 3/22/00
              IF (NAME(IT)(1:1) .EQ. 'H') THEN                  ! 3/22/00
                 NPOT = 131       ! S of S-H, param = 131       ! 3/22/00
                 GO TO 1303                                     ! 3/22/00
              ELSE                                              ! 3/27/00
                 IF (NAME(IT)(1:1) .EQ. 'C') NC = NC + 1        ! 3/27/00
                 IF (NAME(IT)(1:1) .EQ. 'S' .AND.               ! 3/27/00
     X               NCON(IT) .EQ. 2) NS_2 = NS_2 + 1           ! 3/27/00
              ENDIF                                             ! 3/27/00
           ENDDO                                                ! 3/22/00
C-----Check of central S on C-S-S and S-S-S                     ! 3/27/00
        IF (NC .EQ. 1 .AND. NS_2 .EQ. 1) THEN                   ! 3/27/00
            NPOT = 137         ! S of C-S-S, param = 137        ! 3/27/00
            GO TO 1303                                          ! 3/27/00
        ENDIF                                                   ! 3/27/00
        IF (NS_2 .EQ. 2) THEN                                   ! 3/27/00
           NPOT = 138         ! central S of S-S-S, param = 138 ! 3/27/00
           GO TO 1303                                           ! 3/27/00
        ENDIF                                                   ! 3/27/00
C---------Rings??                                               ! 3/22/00
           CALL BOND_CHECK(NUM_ATOM, N_INRINGS, NSIZES, NBONDS) ! 3/22/00
           GO TO (750, 751, 752), NBONDS + 1                    ! 3/22/00
C                  0    1    2    <---# bonds in rings          ! 3/22/00
751           NPOT = 999          ! unknown                     ! 3/22/00
                 GO TO 1303                                     ! 3/22/00
750           NPOT = 133          ! sulfide, C-S-C, param = 133 ! 3/22/00
                 GO TO 1303                                     ! 3/22/00
C------------2 bonds in rings                                   ! 3/22/00
752           IF (NEG_COUNT(NSIZES) .NE. 0) THEN                ! 3/22/00
                 NPOT = 129     ! thiophene S                   ! 3/22/00
              ELSE                                              ! 3/22/00
                 NPOT = 133     ! regular sulfide, param = 133  ! 3/22/00
              ENDIF
              GO TO 1303                                        ! 3/22/00
        ENDIF                                                   ! 3/22/00
C
C----Look at P's 9th........................................    ! 4/4/00
C
3030    IF (AT2(1:1) .EQ. 'P' .AND. NCON(NUM_ATOM) .EQ. 3) THEN ! 4/4/00
           NC = 0                                               ! 4/4/00
           DO I=1,3                                             ! 4/4/00
              IT = ICON(I,NUM_ATOM)                             ! 4/4/00
              IF (NAME(IT)(1:1) .EQ. 'C') NC = NC + 1           ! 4/4/00
           ENDDO                                                ! 4/4/00
           IF (NC. EQ. 3) THEN                                  ! 4/4/00
              NPOT = 143      ! P of C3P                        ! 4/4/00
              GO TO 1303      ! param = 143                     ! 4/4/00
           ENDIF                                                ! 4/4/00
        ENDIF                                                   ! 4/4/00
C
2007   PRINT 2011, NAME(NUM_ATOM)
2011   FORMAT (' Atom ',A5,' is unknown...type set to 999')
1034      NPOT = 999     ! 4 connections ... flag as unknown
1303   continue
       RETURN
       END
C---------------------------------------------------------------
C----Function to determine if a H atom is a cubane-linked H
C
      FUNCTION IQ_CUBANE_H(I1)
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_CUBANE_H = 0   ! 0 indicates its not cubane-related
C
      I2 = ICON(1,I1)
      IF (NAME(I2)(1:1) .NE. 'C' .OR. NCON(I2) .NE. 4) RETURN
C
C----Is the C a cubane C?
      IF (IQ_CUBANE(I2) .EQ. 0) RETURN
      IQ_CUBANE_H = 36
C
      RETURN
      END
C---------------------------------------------------------------
C----Function to determine if a C atom is a cubane C
C
      FUNCTION IQ_CUBANE(I1)   ! I1 = 4-linked C
      logical phelp, pphelp
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      DIMENSION NUM_C(4), NUM_ALL(12), NUM_ALL2(30),
     x          NCON2(500)
C
      phelp = .false.
      pphelp = .false.
      IQ_CUBANE = 0
      IF (NCON(I1) .NE. 4) RETURN   ! C must be linked to 4 atom
C
C----Find the last C atom & adjust NCON2 to count only C's
      DO 10 I=1,NA
         IF (NAME(I)(1:1) .EQ. 'C') THEN
            NCON2(I) = 0
         ELSE
            LAST = I - 1
            GO TO 15
         ENDIF
         DO 8 J=1,NCON(I)
            I2 = ICON(J,I)
            IF (NAME(I2)(1:1) .NE. 'C') GO TO 10
            NCON2(I) = NCON2(I) + 1
8        CONTINUE
10    CONTINUE
15    continue
      if (phelp) then
         print 1232, i1, name(i1)
1232     format (' in iq_cubane: i1, name(i1) =',i3,1x,a5)
         do i=1,last
            n = ncon2(i)
            print 1234, i, name(i), n, (icon(j,i),j=1,n)
1234        format (' i, name, ncon2, icon =',i3,1x,a5,i4,4i3)
         enddo
      endif
C
      IF (NCON2(I1) .LT. 3) RETURN
C
C----How many C's bonded to I1...must be at least 3
      N_3 = 0   ! # of linked C's with 3 C links
      DO 300 I=1,NCON2(I1)
          I2 = ICON(I,I1)                 ! I2-I1 link
          IF (NAME(I2)(1:1) .EQ. 'C') THEN
             IF (NCON2(I2) .LT. 3) GO TO 300
             N_3 = N_3 + 1
             NUM_C(N_3) = I2
          ENDIF
300   CONTINUE
      if (pphelp) print 1228, n_3, (num_c(i),i=1,n_3)
1228  format (' n_3, num_c =',5i3)
      IF (N_3 .LT. 3) RETURN  ! must have at least 3 3-linked C's
C----Look at the 3 or 4-linked C's attached to I1
      NT = 0
      DO 400 I=1,N_3
         I2 = NUM_C(I)
         N2 = NCON2(I2)
         DO 400 J=1,N2
            I3 = (ICON(J,I2))    ! I3...I2...I1
            DO 340 L=1,N_3
               IF (I3 .EQ. NUM_C(L)) GO TO 400
340         CONTINUE
            IF (I3 .EQ. I1) GO TO 400
            IF (NCON2(I3) .LT. 3) GO TO 400
         NT = NT + 1
         NUM_ALL(NT) = I3
400   CONTINUE
      if (pphelp) print 1239, nt, (num_all(i),i=1,nt)
1239  format (' nt, num_all =',13i3)
      IF (NT .LT. 6) RETURN
C----Must be 3 pairs of identical corners in NUM_ALL
      NPAIRS = 0
      DO 410 I=1,NT-1
         IF (NUM_ALL(I) .LT. 0) GO TO 410
         DO 409 J=I+1,NT
            IF (NUM_ALL(J) .LT. 0) GO TO 409
            IF (NUM_ALL(I) .EQ. NUM_ALL(J)) THEN
               NPAIRS = NPAIRS + 1
               NUM_ALL(J) = -1
            ENDIF
409      CONTINUE
410      CONTINUE
      IF (NPAIRS .NE. 3) RETURN
C----The I3 atoms in NUM_ALL should have a common connection (I4)
      NT2 = 0
      DO 500 I=1,NT
         I3 = NUM_ALL(I)
         IF (I3 .LT. 0) GO TO 500
         DO 450 J=1,NCON2(I3)
            I4 = ICON(J,I3)
            IF (I4 .EQ. I3) GO TO 450
            IF (NCON2(I4) .GE. 3) THEN
               DO 440 L=1,NCON2(I1)
                  IF (I4 .EQ. ICON(L,I1)) GO TO 450
440            CONTINUE
               DO 460 L=1,NT
                  IF (I4 .EQ. NUM_ALL(L)) GO TO 450
460            CONTINUE
               NT2 = NT2 + 1
               NUM_ALL2(NT2) = I4
            ENDIF
450      continue
500   CONTINUE
      if (pphelp) print 1299, nt2, (num_all2(i),i=1,nt2)
1299  format (' nt2, num_all2 =',31i3)
      IF (NT2 .LT. 3) RETURN
C----How many of the same I4's are in NUM_ALL2
      DO 600 I=1,NT2-1
         N = 1
         DO 550 J=I+1,NT2
            IF (NUM_ALL2(I) .EQ. NUM_ALL2(J)) N = N + 1
550      CONTINUE
         IF (N .EQ. 3) THEN
            IQ_CUBANE = 5
            RETURN
         ENDIF
600   CONTINUE
      RETURN
      END
C----------------------------------------------------------------
C----Function to determine if an O is a nitro O and if the nitro
C     is cubane-linked; also does the same of the N.
C
      FUNCTION IQ_NITRO_OX(IT)    ! IT is the N linked to the O
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_NITRO_OX = 6     ! normal nitro group O
C
      NO = 0           ! count total O's
      NO_1 = 0         ! count 1-linked O's
      NC = 0
      DO 100 I=1,3
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'O') THEN
            NO = NO + 1
            IF(NCON(I1) .EQ. 1) NO_1 = NO_1 + 1
         ENDIF
         IF (NAME(I1)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            IC = I1  ! IC = C attached to the N
         ENDIF
100   CONTINUE
C
      IF (NO .EQ. 3) RETURN   ! -O-NO2 nitro O ?
C
      IF (NC .EQ. 1 .AND. NO_1 .EQ. 2) GO TO 200
         RETURN
200   IF (IQ_CUBANE(IC) .NE. 0) IQ_NITRO_OX = 11 ! cubane-linked nitro
      RETURN                                     ! group; O param = 11
C
C----Is is a cubane-linked nitro N?
C
      ENTRY IQ_NITRO_N(IT)
C
      IQ_NITRO_N = 0     ! NO!
C
      NO = 0         ! count total O's
      NO_1 = 0       ! count 1-linked O's
      NC = 0
      NN = 0
      DO 300 I=1,3
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'O') THEN
             NO = NO + 1
             IF (NCON(I1) .EQ. 1) NO_1 = NO_1 + 1
         ENDIF
         IF (NAME(I1)(1:1) .EQ. 'N') NN = NN + 1
         IF (NAME(I1)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            IC = I1  ! IC = C attached to the N
         ENDIF
300   CONTINUE
C----If C-NO2, check to see if the C is 3-linked.  If yes code it as
C     C(sp2)-linked N....param = 110
      IF (NO_1 .EQ. 2 .AND. NC .EQ. 1 .AND. NCON(IC) .EQ. 3) THEN      ! 10-21-96
         IQ_NITRO_N = 110      ! C(sp2)-linked nitro N                 ! 10-21-96
         RETURN                                                        ! 10-21-96
      ENDIF                                                            ! 10-21-96
      IF (NO_1 .NE. 2) RETURN
      IF (NO .EQ. 3 .OR. NN .EQ. 1) THEN
         IQ_NITRO_N = 7  ! param = 7; code nitro N for N-NO2 or O-NO2
         RETURN          ! as normal nitro N
      ENDIF
      IF (NC .EQ. 1 .AND. NO .EQ. 2) GO TO 400
         RETURN
400   IF (IQ_CUBANE(IC) .NE. 0) THEN
         IQ_NITRO_N = 37  ! cubane-linked nitro group; N param = 37
      ELSE
         IQ_NITRO_N = 7       ! C-NO2, not cubane linked, param = 7
      ENDIF
      RETURN
      END
C----------------------------------------------------------------
C----Function to determine if 1 and 2-linked N's are azide N's
C     (-N=N=N) or diazo N's (C=N=N)
C
      FUNCTION IQ_AZIDE_N1(I1)    ! I1 is a 1-linked N
      CHARACTER NAME*5
      DIMENSION N_NAMES(2), NUM_C(2)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
C
      IQ_AZIDE_N1 = 0
C
C----What is the 1-linked N attached to?
      I2 = ICON(1,I1)  ! atom linked to I1
      IF (NCON(I2) .NE. 2 .OR. NAME(I2)(1:1) .NE. 'N') RETURN
C----N(I1)-N(I2)-?    What is N(I2) attached to?
      NN = 0
      NOTHER = 0
      DO I=1,2
         I3 = ICON(I,I2)
         IF (NAME(I3)(1:1) .EQ. 'N') THEN
            NN = NN + 1
         ELSE
            NOTHER = NOTHER + 1
            NO = I3
         ENDIF
      ENDDO
         IF (NN .EQ. 2) THEN
            IQ_AZIDE_N1 = 39        ! azide N, param = 39
            RETURN
         ENDIF
         IF (NN .EQ. 1 .AND. NAME(NO)(1:1) .EQ. 'C') THEN
            D = DISTANCE(I1,I2)
            IF (D .GE. 1.16) RETURN        ! N-N must be < 1.16 Angs
            IQ_AZIDE_N1 = 69               ! end N of a C=N=N
         ENDIF
         RETURN
C
C---- 2-linked N..is it a central (N2) or C-linked (N1) azide N
C      or central N (N2) in diazo group ?
C
      ENTRY IQ_AZIDE_N2(N1)    ! N1 or N2 of C-N1-N2-N3
C                              ! or N2 or C-N2-N1
      IQ_AZIDE_N2 = 0
C
      NN = 0
      NOTHER = 0
      NC = 0
      DO 100 I=1,2
         J1 = ICON(I,N1)
         IF (NAME(J1)(1:1) .EQ. 'N') THEN
            NN = NN + 1         ! # N's attached to the N
            N_NAMES(NN) = J1    ! N's attached to N1
            GO TO 100
         ENDIF
         IF (NAME(J1)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            NUM_C(NC) = J1    ! NUM_C = C attached to N1
            GO TO 100
         ENDIF
            NOTHER = NOTHER + 1  ! something else
100   CONTINUE
C
      IF (NOTHER .NE. 0) RETURN  ! something else...quit
C
      IF (NC .EQ. 1 .AND. NN .EQ. 1) THEN  ! inner most...-C-N1 <--?
         N2 = N_NAMES(1)
         N_N2 = NCON(N2)   ! # links to N2...-C-N1-N2
         IF (N_N2 .EQ. 2) GO TO 105   ! not the inner azide N
C----Is it the inner N of a diazo group ?
         IF (N_N2 .EQ. 1) THEN
            IQ_AZIDE_N2 = 69    ! n of C=n=N
C
            IT_TABLE(N1) = 1          ! n of C=n=N
            CALL ICHEM3D (I3D(N1),-1,69) ! type = 69
            KEEP_NPOT(N1) = 69
C
            IT_TABLE(N2) = 1          ! n of C=N=n
            CALL ICHEM3D (I3D(N2),-1,69) ! type = 69
            KEEP_NPOT(N2) = 69
C
            NC = NUM_C(1)
            IT_TABLE(NC) = 1          ! C of C=N=N
            CALL ICHEM3D (I3D(NC),-1,71) ! type = 71
            KEEP_NPOT(NC) = 71
         ENDIF
         RETURN
c
105      DO 110 I=1,2
            N3 = ICON(I,N2)
            IF (N3 .EQ. N1) GO TO 110
110      CONTINUE
         IF (NCON(N3) .EQ. 1 .AND.              ! is N3 1-linked & a N ?,
     X       NAME(N3)(1:1) .EQ. 'N') IQ_AZIDE_N2 = 39  ! yes...it's an azide
         RETURN
      ENDIF
C
      IF (NN .NE. 2) RETURN     ! not the middle N
         M1 = N_NAMES(1)
         M2 = N_NAMES(2)
         IF (NCON(M1) .EQ. 1 .OR. NCON(M2) .EQ. 1) THEN
            IQ_AZIDE_N2 = 39    ! its the central azide N
         ENDIF
      RETURN
      END
C----------------------------------------------------------------
C----Function to determine S and O of sulfoxides, sulfones,    ! 3/29/00
C      sulfonates, sulfonamides and sulfonimines               ! 11/7/00
C
C       O    O = 120       O2   O = 123         O2   O = 126   ! 3/20/00
C       |    S = 119       |    S = 122         |    S = 125   ! 3/20/00
C     C-S-C              C-S-C                C-S-O-C  O = 127 ! 3/20/00
C     [S-O =121]         [S-O2 = 124]         [SO3 = 128]      ! 3/20/00
C                                                              ! 3/18/00
C         O2   O = 139     O2     O = 152                      ! 11/3/00
C         |    S = 140     |      S = 153                      ! 11/3/00
C       C-S-N  N = 141   C-S-N=C  N = 154 & C = 155            ! 11/3/00
C       [SO2-N = 142]    [SO2-N = 151]                         ! 11/3/00
C
      FUNCTION IQ_S_O (IO, IS)  ! if IO .ne. IS, IO is input O ! 11/3/00
C                               !   and IS is attached S       ! 11/3/00
C                               ! if IO = IS, input atom is a  ! 11/3/00
      CHARACTER NAME*5          !  a 3 or 4-linked S           ! 11/3/00
      LOGICAL S_ATOM                                           ! 3/20/00
      DIMENSION N_NAMES(2), NUM_O(2), NUM_N(3), NUM_N2(3)      ! 11/3/00
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
C
      S_ATOM = .false.                                         ! 3/20/00
C
C----If IO = IS, the atom is a S                               ! 3/20/00
      IF (IO .EQ. IS) S_ATOM = .true.                          ! 3/20/00
      IQ_S_O = 0                                               ! 3/18/00
      NS = NCON(IS)                     ! # links to S         ! 3/20/00
      IF (NS .GT. 4) RETURN                                    ! 3/20/00
C
C----Count # C's, # 1 and # 2-linked O's and 3-linked N's      ! 3/29/00
C     to the S                                                 ! 3/29/00
      NC = 0                            ! # C's linked to S    ! 3/20/00
      NO_1 = 0                          ! # 1-linked O's       ! 3/18/00
      NO_2 = 0                          ! # 2-linked O's       ! 3/20/00
      NN_2 = 0                          ! # 2-linked N's       ! 11/3/00
      NN_3 = 0                          ! # 3-linked N's       ! 3/29/00
      DO 200 I=1,NS                                            ! 3/18/00
         I2 = ICON(I,IS)                                       ! 3/18/00
         IF (NAME(I2)(1:1) .EQ. 'C') NC = NC + 1               ! 3/18/00
         IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 1) THEN! 3/18/00
            NO_1 = NO_1 + 1                                    ! 3/18/00
            NUM_O(NO_1) = I2                                   ! 3/18/00
         ENDIF                                                 ! 3/18/00
         IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 2) THEN! 3/20/00
            NO_2 = NO_2 + 1                                    ! 3/20/00
            NUM_O_2 = I2                                       ! 3/20/00
         ENDIF                                                 ! 3/20/00
         IF (NAME(I2)(1:1) .EQ. 'N' .AND. NCON(I2) .EQ. 3) THEN! 3/29/00
            NN_3 = NN_3 + 1        ! looking for               ! 3/29/00
            NUM_N(NN_3) = I2       ! sulfonamides              ! 3/29/00
         ENDIF                                                 ! 3/29/00
         IF (NAME(I2)(1:1) .EQ. 'N' .AND. NCON(I2) .EQ. 2) THEN! 11/3/00
            NN_2 = NN_2 + 1         ! looking for              ! 11/3/00
            NUM_N2(NN_2) = I2       ! sulfonimines             ! 11/3/00
         ENDIF                                                 ! 11/3/00
200   CONTINUE                                                 ! 3/18/00
C
      GO TO (210, 210, 213, 214), NS                           ! 3/18/00
210      RETURN                                                ! 3/18/00
C                                                              ! 3/18/00
C-------Sulfoxide?                                             ! 3/18/00
213      IF (NC .EQ. 2 .AND. NO_1 .EQ. 1) THEN                 ! 3/18/00
            IQ_S_O = 120            ! sulfoxide O param = 120  ! 3/18/00
            IT_TABLE(IS) = 1              ! S of sulfoxide     ! 3/18/00
            CALL ICHEM3D (I3D(IS),-1,119) ! param = 119        ! 3/18/00
            KEEP_NPOT(IS) = 119                                ! 3/18/00
            IF (S_ATOM) IQ_S_O = 119   ! it's a S              ! 3/20/00
         ENDIF                                                 ! 3/18/00
         RETURN                                                ! 3/18/00
C
C-------Sulfone, sulfonate or sulfonamide?                     ! 3/29/00
214      IF (NC .EQ. 2 .AND. NO_1 .EQ. 2 .AND. NO_2 .EQ. 0)    ! 3/20/00
     x           GO TO 216      ! sulfone                      ! 11/3/00
         IF (NC .EQ. 1 .AND. NO_1 .EQ. 2 .AND. NO_2 .EQ. 1)    ! 3/20/00
     x           GO TO 218      ! sulfonate                    ! 11/3/00
         IF (NC .EQ. 1 .AND. NO_1 .EQ. 2 .AND. NN_3 .EQ. 1)    ! 3/29/00
     x           GO TO 220      ! sulfonamide                  ! 3/29/00
         IF (NC .EQ. 1 .AND. NO_1 .EQ. 2 .AND. NN_2 .EQ. 1)    ! 11/3/00
     x           GO TO 222      ! sulfonimine ?                ! 11/3/00
         RETURN                         ! none of these        ! 3/29/00
C
C-------Sulfone!                                               ! 3/20/00
216         IQ_S_O = 123            ! sulfone O, param = 123   ! 3/18/00
            IT_TABLE(IS) = 1                 ! S of sulfone    ! 3/18/00
            CALL ICHEM3D (I3D(IS),-1,122)    ! param = 122     ! 3/18/00
            KEEP_NPOT(IS) = 122                                ! 3/18/00
            IF (S_ATOM) IQ_S_O = 122    ! it's a S             ! 3/20/00
C----------Code the 2nd sulfone O
            DO 250 I=1,2                                       ! 3/18/00
               J = NUM_O(I)             ! look for other O     ! 3/18/00
               IF (J .EQ. IO) GO TO 250                        ! 3/18/00
               IT_TABLE(J) = 1                                 ! 3/18/00
               CALL ICHEM3D (I3D(J),-1,123)                    ! 3/18/00
               KEEP_NPOT(J) = 123                              ! 3/18/00
250         CONTINUE                                           ! 3/18/00
      RETURN                                                   ! 3/20/00
C
C--------Sulfonate!                                            ! 3/20/00
218         IF (S_ATOM) THEN    ! S of a sulfonate             ! 3/20/00
               IQ_S_O = 125     ! param = 125                  ! 3/20/00
               RETURN                                          ! 3/20/00
            ENDIF                                              ! 3/20/00
C
            IF (NCON(IO) .EQ. 1) THEN    ! 1-linked O          ! 3/20/00
               IQ_S_O = 126              ! param = 126         ! 3/20/00
            ELSE                                               ! 3/20/00
               IQ_S_O = 127              ! 2-linked O          ! 3/20/00
            ENDIF                                              ! 3/20/00
      RETURN                                                   ! 3/18/00
C
C-------Sulfonamide!                                           ! 3/29/00
220         IQ_S_O = 139        ! sulfonamide O, param = 139   ! 3/29/00
            IT_TABLE(IS) = 1               ! S of sulfonamide  ! 3/29/00
            CALL ICHEM3D (I3D(IS),-1,140)    ! param = 140     ! 3/29/00
            KEEP_NPOT(IS) = 140                                ! 3/29/00
            IF (S_ATOM) IQ_S_O = 140    ! it's a S             ! 3/29/00
C----------Code the 2nd sulfonamide O                          ! 3/29/00
            DO 253 I=1,2                                       ! 3/29/00
               J = NUM_O(I)             ! look for other O     ! 3/29/00
               IF (J .EQ. IO) GO TO 253                        ! 3/29/00
               IT_TABLE(J) = 1                                 ! 3/29/00
               CALL ICHEM3D (I3D(J),-1,139)                    ! 3/29/00
               KEEP_NPOT(J) = 139                              ! 3/29/00
253         CONTINUE                                           ! 3/29/00
C----------Code the sulfonamide N (N's??)                      ! 3/29/00
            DO 256 I=1,NN_3                                    ! 3/29/00
               J = NUM_N(I)                                    ! 3/29/00
               IT_TABLE(J) = 1                                 ! 3/29/00
               CALL ICHEM3D (I3D(J),-1,141)                    ! 3/29/00
               KEEP_NPOT(J) = 141                              ! 3/29/00
256         CONTINUE                                           ! 3/29/00
      RETURN                                                   ! 3/29/00
C
C-------Sulfonimine...C-SO2-N=C                                ! 11/3/00
222         KN = NUM_N2(1)       ! check imine N               ! 11/3/00
            DO 224 J=1,2         ! look for 3-linked C         ! 11/3/00
               IC = ICON(J,KN)                                 ! 11/3/00
               IF (NAME(IC)(1:1) .EQ. 'S') GO TO 224           ! 11/3/00
               IF (NAME(IC)(1:1) .EQ. 'C' .AND.                ! 11/3/00
     X             NCON(IC) .EQ. 3) GO TO 226      ! ok        ! 11/3/00
               RETURN                              ! not ok    ! 11/3/00
224         CONTINUE                                           ! 11/3/00
226         IQ_S_O = 152        ! sulfonimine O, param = 152   ! 11/3/00
            IT_TABLE(IS) = 1               ! S of sulfonimine  ! 11/3/00
            CALL ICHEM3D (I3D(IS),-1,153)    ! param = 153     ! 11/3/00
            KEEP_NPOT(IS) = 153                                ! 11/3/00
            IF (S_ATOM) IQ_S_O = 153    ! it's a S             ! 11/3/00
C----------Code the 2nd sulfonimine O                          ! 11/3/00
            DO 1253 I=1,2                                      ! 11/3/00
               J = NUM_O(I)             ! look for other O     ! 11/3/00
               IF (J .EQ. IO) GO TO 1253                       ! 11/3/00
               IT_TABLE(J) = 1                                 ! 11/3/00
               CALL ICHEM3D (I3D(J),-1,152)                    ! 11/3/00
               KEEP_NPOT(J) = 152                              ! 11/3/00
1253        CONTINUE                                           ! 11/3/00
C----------Code the sulfonimine N                              ! 11/3/00
            DO 1256 I=1,NN_2                                   ! 11/3/00
               J = NUM_N2(I)                                   ! 11/3/00
               IT_TABLE(J) = 1                                 ! 11/3/00
               CALL ICHEM3D (I3D(J),-1,154)                    ! 11/3/00
               KEEP_NPOT(J) = 154                              ! 11/3/00
1256         CONTINUE                                          ! 11/3/00
C----------Code the sulfonimine C                              ! 11/3/00
               IT_TABLE(IC) = 1                                ! 11/3/00
               CALL ICHEM3D (I3D(IC),-1,155)                   ! 11/3/00
               KEEP_NPOT(IC) = 155                             ! 11/3/00
      RETURN                                                   ! 11/3/00
      END
C-------------------------------------------------------------------
C----Function to determine if a N is one of the 4 N's in the
C     Z-tetraazapentalene moiety   -N(-)-N--N(+)=N-C
C                                        |  |
C                                        C  C
C
      FUNCTION IQ_N_AZAPENTALENE(ITT) ! ITT is a N
C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      DIMENSION  IND(2), XL(2), VEC(3,2), INN(2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3), (IT,IN1)
C
      IQ_N_AZAPENTALENE = 0   ! 0 means can't identify type
C
      IT = ITT
      IF (NCON(IT) .LT. 2 .OR. NCON(IT) .GT. 3) RETURN
      GO TO (10,20), (NCON(IT) - 1)
C----IT has 2 attachments...check for an end N....N1 of N1-N2-N3-N4
10    NN = 0
      NC = 0
      DO 100 I=1,2
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            IF (NCON(I1) .EQ. 3) IN2 = I1  ! a possible N2
         ENDIF
100   CONTINUE
      IF (NC .NE. 1 .OR. NN .NE. 1) RETURN
C----IN2 has 3 attachments...possible N2 of N1-N2-N3-N4
      NN = 0
      NC = 0
      DO 110 I=1,3
         I1 = ICON(I,IN2)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            IF (I1 .NE. IT) IN3 = I1   ! a possible N3
         ENDIF
110   CONTINUE
      IF (NC .NE. 1 .OR. NN .NE. 2) RETURN
C----IN3 is a possible N3 of N1-N2-N3-N4
      NN = 0
      NC = 0
      DO 210 I=1,3
         I1 = ICON(I,IN3)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            IF (I1 .NE. IN2 .AND. NCON(I1) .EQ. 2) IN4 = I1  ! a possible N4
         ENDIF
210   CONTINUE
      IF (NC .NE. 1 .OR. NN .NE. 2) RETURN
C----IN4 is a possible N4 of N1-N2-N3-N4
      NN = 0
      NC = 0
      DO 300 I=1,2
         I1 = ICON(I,IN4)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') NN = NN + 1
300   CONTINUE
      IF (.NOT. (NC .EQ. 1 .AND. NN .EQ. 1)) RETURN
C--------Have atom numbers for N1-N2-N3-N4.  Confirm that N1 & N3 are part
C          are joined by 2 C's.  Ditto for N2 and N4.
         DO 400 I=1,NCON(IN1)            ! start with N1
            J1 = ICON(I,IN1)
            IF (NAME(J1)(1:1) .NE. 'C') GO TO 400
               DO 395 J=1,NCON(J1)       ! J1 is a C linked to N1
                  J3 = ICON(J,J1)
                  IF (NAME(J3)(1:1) .NE. 'C') GO TO 395
                  DO 390 K=1,NCON(J3)    ! J3 is a C linked to J1
                     K1 = ICON(K,J3)     ! is J3 linked to IN3
                     IF (K1 .NE. IN3) GO TO 390
                     GO TO 410
390               CONTINUE
395            CONTINUE
400      CONTINUE
C--------Now check the N2...N4 linkage
410      DO 500 I=1,NCON(IN4)            ! start with N4
            J1 = ICON(I,IN4)
            IF (NAME(J1)(1:1) .NE. 'C') GO TO 500
               DO 495 J=1,NCON(J1)       ! J1 is a C linked to N4
                  J3 = ICON(J,J1)
                  IF (NAME(J3)(1:1) .NE. 'C') GO TO 495
                  DO 490 K=1,NCON(J3)    ! J3 is a C linked to J1
                     K1 = ICON(K,J3)     ! is J3 linked to IN2
                     IF (K1 .NE. IN2) GO TO 490
                     GO TO 510
490               CONTINUE
495            CONTINUE
500      CONTINUE
         RETURN
510      IQ_N_AZAPENTALENE = 19   ! param = 19
      RETURN                                       ! for Z-tetraazapentalene N
C----3 attachments...check for a middle N, N2 of N1-N2-N3-N4
20    NN = 0
      NC = 0
      DO 220 I=1,3
         I1 = ICON(I,ITT)
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            INN(NN) = I1      ! number of the attached N
         ENDIF
220   CONTINUE
      IF (.NOT. (NC .EQ. 1 .AND. NN .EQ. 2)) RETURN     ! attached N's must have
         I1 = MIN0(NCON(INN(1)),NCON(INN(2)))      ! NCON's must be 2 & 3,
         I2 = MAX0(NCON(INN(1)),NCON(INN(2)))      ! respectively
         IF (.NOT. (I1 .EQ. 2 .AND. I2 .EQ. 3)) RETURN
            IT = INN(1)
            IF (NCON(INN(2)) .EQ. 2) IT = INN(2)
            GO TO 10
      END
C---------------------------------------------------------------
C----Function to determine if a 2-linked O is an hydroxyl oxygen
C        -C-O-H or N-O-H or -O-O-H
C
      FUNCTION IQ_OH_O(I1)      ! I1 = 2-linked O
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
       DIMENSION NUM_C(2), NUM_H(2), NUM_N(2)
C
      IQ_OH_O = 0   ! 0 means NO!
C
C----Must be a least one C/N bonded to I1 and one H.
      NC = 0
      NH = 0
      NN = 0
      NO = 0
      DO 20 I=1,2
          I2 = ICON(I,I1)    ! I2 is linked to I1
          IF (NAME(I2)(1:1) .EQ. 'C') THEN
              NC = NC + 1
              NUM_C(NC) = I2
          ENDIF
          IF (NAME(I2)(1:1) .EQ. 'H') THEN
              NH = NH + 1
              NUM_H(NH) = I2
           ENDIF
          IF (NAME(I2)(1:1) .EQ. 'N') THEN
              NN = NN + 1
              NUM_N(NN) = I2
          ENDIF
          IF (NAME(I2)(1:1) .EQ. 'O') NO = NO + 1
20    CONTINUE
      IF (NH .NE. 1) RETURN
      IF (NO .EQ. 1) THEN
          IQ_OH_O = 31             ! hydroxyl O in -O-O-H
          IH = 30
          GO TO 60
      ENDIF
      IF (NN .EQ. 1 .AND. NCON(NUM_N(1)) .EQ. 3) THEN
          IQ_OH_O = 31             ! hydroxyl O in -NH-OH or R2N-OH
          IH = 30
          GO TO 60
      ENDIF
      IF (NN .EQ. 1 .AND. NCON(NUM_N(1)) .EQ. 2) THEN
          IQ_OH_O = 67  ! hydroxyl O...like in C=N-OH
          IH = 68
          GO TO 60
      ENDIF
      IF (NC .NE. 1) RETURN
         N = NUM_C(1)
         IF (NCON(N) .EQ. 3) THEN
            IQ_OH_O = 67  ! Csp2-O-H
            IH = 68
            GO TO 60
         ENDIF
         IF (NCON(N) .EQ. 4) THEN
            IQ_OH_O = 31  ! Csp3-O-H
            IH = 30
         ENDIF
C----Take care of H & O
60    NH = NUM_H(1)
      IT_TABLE(NH) = 1          ! H on O atom
      CALL ICHEM3D (I3D(NH),-1,IH) ! type = 30 or 68
      KEEP_NPOT(NH) = IH
C
      IT_TABLE(I1) = 1          ! O atom
      CALL ICHEM3D (I3D(I1),-1,IQ_OH_O) ! type = 31 or 67
      KEEP_NPOT(I1) = IQ_OH_O
C
      RETURN
      END
C
C---------------------------------------------------------------
C----Function to determine if a 3-linked C is a carbonyl carbon
C
      FUNCTION IQ_CARBONYL(I1, ICHANGE)     ! I1 = 3-linked C
C
      logical phelp
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x              XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x              I3D(500), KEEP_NPOT(500)
      COMMON /CARBONYL/ N_O_CO
      DIMENSION NUM_O(3), NSIZES(10)                            ! 10/17/00
C
      phelp = .false.
C
      IQ_CARBONYL = 0   ! NO...initially
C
C----Find the carbonyl oxygen...1-linked
      NO = 0
      NN = 0
      DO 20 I=1,3
          I2 = ICON(I,I1)
          IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 1) THEN
             NO = NO + 1
             NUM_O(NO) = I2
          ENDIF
          IF (NAME(I2)(1:1) .EQ. 'N') NN = NN + 1
20    CONTINUE
ccc      print 1122, i1, name(i1), no, (num_o(jj),name(num_o(jj)),jj=1,no)  ! #########
ccc1122  format ('in iq_carbonyl: i1, name, no, num_o =',                   ! #########
ccc     x        i3,1x,a5,i4,3(i4,1x,a5))                                     ! #########
      IF (NO .NE. 1) RETURN    ! not carbonyl
C----Regular carbonyl C or N-linked?
      IF (NN .EQ. 0) THEN
           CALL BOND_CHECK (I1, NRINGS, NSIZES, NBONDS)
ccc           print 1111, nn, i1, name(i1), nbonds, nrings,    ! ###################
ccc     x                 no,(num_o(jj),jj=1,no)               ! ###################
ccc1111       format ('**> nn, i1, name, nbonds nrings =',2i3,1x,a5,2i3/  ! ###################
ccc     x             '    no, num_o =',4i3)                              ! ###################
           IF (NBONDS .EQ. 0 .OR. NBONDS .EQ. 2) GO TO 21  ! return if bonds .ne.
                RETURN                                     !   0 or 2
21         IF (NBONDS .EQ. 0) THEN
              J = 32   ! regular carbonyl C not in a ring
           ELSE
              J = 96   ! regular carbonyl C in a ring
           ENDIF
      ELSE
         J = 79     ! N-linked
      ENDIF
C----Calculate distance...must be between 1.13 & 1.31 Angs
22    I3 = NUM_O(1)
ccc      print 1198, i3, num_o(1)                         ! ####################
ccc1198  format ('before distance: i3, num_o(1) =',2i4)   ! ####################
      D = DISTANCE(I3,I1)
ccc      print 1199, i3, i1, d                            ! ####################
ccc1199  format ('after distance: i3, i1, d =',2i4,f6.2)  ! ####################
      IF (D .GE. 1.13 .AND. D .LE. 1.311) THEN
         IQ_CARBONYL = J     ! carbonyl C; param = 32, 96 or 79
         N_O_CO = I3
C
         IT_TABLE(I1) = 1             ! carbonyl C
         CALL ICHEM3D (I3D(I1),-1,J)  ! type = 32, 96 or 79
         KEEP_NPOT(I1) = J
C
         IT_TABLE(I3) = 1     ! carbonyl O, regular or N-associated
         IF (NN .EQ. 0) THEN
           J = 15     ! regular
         ELSE
           J = 80     ! N-associated
         ENDIF
         IT_TABLE(I3) = 1
         CALL ICHEM3D (I3D(I3),-1,J)       ! type = 15 or 80
         KEEP_NPOT(I3) = J
      ENDIF
      if (phelp) print 1211, i1, name(i1), i3, name(i3),
     x           keep_npot(i1), keep_npot(i3), d, iq_carbonyl
1211  format ('   end of iq_carbonyl: i1, name(i1), i3, name(i3), ',
     x        'keep_npot(i1), keep_npot(i3), d, iq_carbonyl ='/
     x         2(i3,1x,a5,1x),2i3,f6.3,i3)
      RETURN
      END
C
C-------------------------------------------------------------
C----Function to determine if a 3-linked C is a carbonyl carbon;
C     special routine to be called by BOND_CHECK
C
      FUNCTION IQ_CARBONYL2(I1)     ! I1 = 3-linked C
C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      COMMON /CARBONYL/ N_O_CO
      DIMENSION NUM_O(3)
C
      IQ_CARBONYL2 = 0   ! NO...initially
C
C----Find the carbonyl oxygen...1-linked
      NO = 0
      NN = 0
      DO 20 I=1,3
          I2 = ICON(I,I1)
          IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 1) THEN
             NO = NO + 1
             NUM_O(NO) = I2
          ENDIF
          IF (NAME(I2)(1:1) .EQ. 'N') NN = NN + 1
20    CONTINUE
      IF (NO .NE. 1) RETURN    ! not carbonyl
C----Regular carbonyl C or N-linked?
      IF (NN .EQ. 0) THEN
         J = 32   ! regular carbonyl C...ring or no ring
      ELSE
         J = 79   ! N-linked
      ENDIF
C----Calculate distance...must be between 1.13 & 1.31 Angs
22    I3 = NUM_O(1)
      D = DISTANCE(I3,I1)
      IF (D .GE. 1.13 .AND. D .LE. 1.311) THEN
         IQ_CARBONYL2 = J     ! carbonyl C; param = 32 or 79
         N_O_CO = I3
      ENDIF
      RETURN
      END
C
C---------------------------------------------------------------
C----Function to reset a carbonyl carbon from the amide
C     to regular designation.  This is when the N of an amide
C     is also the N of a nitramine.  Nitramine N takes
C     precedence.
C
      FUNCTION IQ_CARBONYL3(I1)     ! I1 = 3-linked C
C
      logical phelp
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x              XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x              I3D(500), KEEP_NPOT(500)
      COMMON /CARBONYL/ N_O_CO
      DIMENSION NUM_O(3), NSIZES(10)
C
      phelp = .false.
C
      IQ_CARBONYL3= 0   ! NO...initially
C
C----Find the carbonyl oxygen...1-linked
      NO = 0
      NN = 0
      DO 20 I=1,3
          I2 = ICON(I,I1)
          IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 1) THEN
             NO = NO + 1
             NUM_O(NO) = I2
          ENDIF
          IF (NAME(I2)(1:1) .EQ. 'N') NN = NN + 1
20    CONTINUE
      IF (NO .NE. 1) RETURN    ! not carbonyl
C----What kind of carbonyl C...don't test for an N-link
           CALL BOND_CHECK (I1, NRINGS, NSIZES, NBONDS)
           IF (NBONDS .EQ. 0 .OR. NBONDS .EQ. 2) GO TO 21  ! return if bonds .ne.
                RETURN                                     !   0 or 2
21         IF (NBONDS .EQ. 0) THEN
              J = 32   ! regular carbonyl C not in a ring
           ELSE
              J = 96   ! regular carbonyl C in a ring
           ENDIF
C----Calculate distance...must be between 1.13 & 1.31 Angs
22    I3 = NUM_O(1)
      D = DISTANCE(I3,I1)
      IF (D .GE. 1.13 .AND. D .LE. 1.311) THEN
         IQ_CARBONYL3 = J     ! carbonyl C; param = 32 or 96
         N_O_CO = I3
C
         IT_TABLE(I1) = 1             ! carbonyl C
         CALL ICHEM3D (I3D(I1),-1,J)  ! type = 32 or 96
         KEEP_NPOT(I1) = J
C
         IT_TABLE(I3) = 1             ! carbonyl O
         CALL ICHEM3D (I3D(I3),-1,15)  ! type = 15
         KEEP_NPOT(I3) = 15
      ENDIF
      if (phelp) print 1211, i1, name(i1), i3, name(i3),
     x           keep_npot(i1), keep_npot(i3), d, iq_carbonyl3
1211  format ('   end of iq_carbonyl: i1, name(i1), i3, name(i3), ',
     x        'keep_npot(i1), keep_npot(i3), d, iq_carbonyl3 ='/
     x         2(i3,1x,a5,1x),2i3,f6.3,i3)
      RETURN
      END
C
C---------------------------------------------------------------
C----Function to determine if a 1-linked C is an isocyanide carbon
C          R-N#C:
C            + -
C
      FUNCTION IQ_ISOCYANIDE(I1)     ! I1 = 1-linked C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_ISOCYANIDE = 0
C
      I2 = ICON(1,I1)
      IF (NAME(I2)(1:1) .EQ. 'N' .AND. NCON(I2) .EQ. 2) THEN
C
         IQ_ISOCYANIDE = 45
C
         IT_TABLE(I1) = 1        ! isocyanide C, same as nitrile
         CALL ICHEM3D (I3D(I1),-1,45) ! type = 45
         KEEP_NPOT(I1) = 45
C
         IT_TABLE(I2) = 1        ! isocyanide N, same as nitrile
         CALL ICHEM3D (I3D(I2),-1,46) ! type = 46
         KEEP_NPOT(I2) = 46
      ENDIF
      RETURN
      END
C
C---------------------------------------------------------------
C----Function to determine if a 1-linked O attached to a 2-linked
C     N is a nitroso group
C
      FUNCTION IQ_NITROSO(NUM_ATOM)
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_NITROSO = 0   ! NO!
C
      IN = ICON(1,NUM_ATOM)   ! # of the attached N
C
      IF (NCON(IN) .NE. 2) RETURN    ! not nitroso
C
      IT_TABLE(NUM_ATOM) = 1             ! nitroso O atom
      CALL ICHEM3D (I3D(NUM_ATOM),-1,61) ! type = 61
      KEEP_NPOT(NUM_ATOM) = 61
C
      IT_TABLE(IN) = 1             ! nitroso N atom
      CALL ICHEM3D (I3D(IN),-1,62) ! type = 62
      KEEP_NPOT(IN) = 62
      IQ_NITROSO = 61
C
      RETURN
      END
C---------------------------------------------------------------
C----Function to determine if a 1-linked O is an amine or
C     N(+)-oxide O
C                                               C
C                                               |
C     R3N(+)-O(-); pyridine-N(+)-O(-); R2N-O; C=N(+)-O(-);
C       N = 70       N = 18            N = 89    N = 18
C     C=N(+)-O-
C       |                all O's = 17
C       O(-)
C     N = 18
C
      FUNCTION IQ_AMINEOXIDE(NUM_ATOM)   ! NUM_ATOM = 1-linked oxygen
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      DIMENSION NUM_C(3), NUM_O(3), NSIZES(10)
C
      IQ_AMINEOXIDE = 0   ! NO!
C
      IN = ICON(1,NUM_ATOM)   ! IN = # of the attached N
      IF (NAME(IN)(1:1) .NE. 'N') RETURN
      IF (NCON(IN) .LT. 3) RETURN    ! not amine oxide
C
      NO = 0                ! # of connected O's
      NC = 0                ! # of connected C's
      NO_1 = 0              ! # of connected O's with NCON = 1
      DO 20 I=1,NCON(IN)
         I1 = ICON(I,IN)
         IF (NAME(I1)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            NUM_C(NC) = I1
         ENDIF
         IF (NAME(I1)(1:1) .EQ. 'O') THEN
            NO = NO + 1
            NUM_O(NO) = I1
            IF (NCON(I1) .EQ. 1) NO_1 = NO_1 + 1
         ENDIF
20    CONTINUE
      IF (NO_1 .GT. 1) RETURN   ! eliminate possible NO2 groups
C----R3N(+)-O(-) ?
      IF (NC .EQ. 3 .AND. NO_1 .EQ. 1) THEN
          IT_TABLE(IN) = 1             ! amine N-oxide N atom
          CALL ICHEM3D (I3D(IN),-1,70) ! type = 70
          KEEP_NPOT(IN) = 70
          GO TO 100                    ! success, O param = 17
      ENDIF
C----C=N(+)-O- ?
C      |
C      O(-)
      IF (NC .EQ. 1 .AND. NO .EQ. 2) THEN     ! 1 O must be 2 -linked
         DO 25 I=1,2
            IF (NUM_O(I) .EQ. NUM_ATOM) GO TO 25
            IF (NCON(NUM_O(I)) .EQ. 2) GO TO 41
25       CONTINUE
         RETURN
41       IT_TABLE(IN) = 1             ! N(+)-O(-) N atom
         CALL ICHEM3D (I3D(IN),-1,18) ! type = 18
         KEEP_NPOT(IN) = 18
         GO TO 100
      ENDIF
      IF (.NOT. (NC .EQ. 2 .AND. NO_1 .EQ. 1)) RETURN
C----pyridine-N(+)-O(-) ?
      CALL BOND_CHECK(IN, N_INRINGS, NSIZES, NBONDS)
      IF (NBONDS .EQ. 0) GO TO 45
      IF (NBONDS .EQ. 2 .AND. NEG_COUNT(NSIZES) .EQ. 1) THEN
         IT_TABLE(IN) = 1             ! pyridine N-oxide N atom
         CALL ICHEM3D (I3D(IN),-1,18) ! type = 18
         KEEP_NPOT(IN) = 18
         GO TO 100
      ENDIF
C----not a pyridine N-oxide N,   C
C                                |
C                              C=N(+)-O(-) ?
C     look for C=N (d .le. 1.36 Angs)
45    DO 40 I=1,2
         IC = NUM_C(I)
         IF (DISTANCE(IN,IC) .GT. 1.36) GO TO 40
            IT_TABLE(IN) = 1             ! C=N(+,C)-O(-)
            CALL ICHEM3D (I3D(IN),-1,18) ! type = 18
            KEEP_NPOT(IN) = 18
            GO TO 100
40    CONTINUE
C----If reach here, all possibilities exhausted except R2N-O
      IT_TABLE(IN) = 1
      CALL ICHEM3D (I3D(IN),-1,89) ! param for N = 89
      KEEP_NPOT(IN) = 89
C
100   IQ_AMINEOXIDE = 17
C
      RETURN
      END
C
C---------------------------------------------------------------
C----Function to determine if a 2-linked N is an sp2 N in -C=N-
C
      FUNCTION  IQ_NDOUBLEC(I1)   ! I1 = 2-linked N
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      DIMENSION NUM_C(2), NSIZES(10)
C
      IQ_NDOUBLEC = 0   ! NO!
C
C----Count the number of  2 & 3-linked C's
      NC = 0      ! look for R2C=N- & C=C=N-
      NC_3 = 0
      DO 20 I=1,2
          I2 = ICON(I,I1)
          IF (NAME(I2)(1:1) .EQ. 'C') THEN
             NC = NC + 1
             IF (NCON(I2) .EQ. 3 .OR. NCON(I2) .EQ. 2) THEN
                NC_3 = NC_3 + 1
                NUM_C(NC_3) = I2
             ENDIF
          ENDIF
20    CONTINUE
      IF (NC .LT. 1) RETURN ! must be at least 1 C on the N
      IF (NC_3 .EQ. 0) RETURN
C----Calculate C-N distances...one must be < 1.38 Angs
      DO I=1,NC_3
         I3 = NUM_C(I)
         D = DISTANCE(I3,I1)
         IF (D .LE. 1.38) THEN
            IQ_NDOUBLEC = 47     ! param = 47
         ENDIF
      ENDDO
      IF (IQ_NDOUBLEC .EQ. 0) RETURN ! not C=N
C----Now check to see if it's in a ring
      CALL BOND_CHECK (I1, NRINGS, NSIZES, NBONDS)
         IF (NBONDS .EQ. 2) IQ_NDOUBLEC = 50   ! param = 50
      RETURN
      END
C-----------------------------------------------------------------
C----Function to determine if a 3-linked C is is a carboxylic acid
C     group.
C
      FUNCTION IQ_CO2H(I1)      ! I1 = 3-linked C
      CHARACTER NAME*5
      logical phelp
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x        XYZ(3,500), CELL(6), NA, IT_TABLE(500),I3D(500),
     x        KEEP_NPOT(500)
      DIMENSION N_OXY(3)
C
      phelp = .false.
C
      IQ_CO2H = 0
C
      NO = 0
      N = IQ_CARBONYL(I1,0)   ! is it a carbonyl C ?
      if (phelp) print 1233, i1, name(i1), n
1233  format (' -------> in iq_co2h: i1, name, iq_carbonyl =',
     x        i3,1x,a5,i3)
      IF (N .EQ. 32) GO TO 20
         RETURN    ! return...not a carbonyl C that's not in a ring
20    DO 101 I=1,3
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'O') THEN
            NO = NO + 1
            N_OXY(NO) = I2
         ENDIF
101   CONTINUE
      if (.not. phelp) go to 1244
         print 1211, no, (n_oxy(l),l=1,no)
1211     format (' in iq_co2h: no of O''s, at #''s =',4i3)
         print 1212, (name(n_oxy(l)),l=1,no)
1212     format ('                    atom names =',3(1x,a5))
1244  continue
C
      IF (NO .LT. 2) RETURN   ! not 2 O's linked to C
C----Examine the C-linked O's
      I_OK = 0      ! must = 2 for both O's to have been id'd
      I_OH = 0
      I_CO = 0
      DO 102 J = 1,NO
         N = NCON(N_OXY(J))
      if (phelp) then
         l = n_oxy(j)
         print 1256, j, no, ncon(l), l, name(l)
1256     format (' in iq_co2h: j, # O''s, ncon, # atom, name =',
     x           4i3,1x,a5)
      endif
         IF (N .EQ. 2) THEN
            L = N_OXY(J)
            IF (IQ_OH_O(L) .EQ. 31 .OR. IQ_OH_O(L) .EQ. 67) THEN
               I_OH = N_OXY(J)    ! # of OH oxygen
               I_OK = I_OK + 1
            ENDIF
         ENDIF
105      IF (N .EQ. 1) THEN
            I_CO = N_OXY(J)     ! # of C=O oxygen
            I_OK = I_OK + 1
         ENDIF
102   CONTINUE
      if (phelp) print 1255, i_ok, i_oh, name(i_oh), i_co, name(i_co)
1255  format (' in iq_co2h: i_ok, i_oh, name, i_co, name =',
     x        2i3,1x,a5,i3,1x,a5)
      IF (I_OK .NE. 2) RETURN
C----Success...record the appropriate info in IT_TABLE
106   IT_TABLE(I1) = 1   ! carbonyl C
      CALL ICHEM3D (I3D(I1),-1,32)   ! C=O carbon atom type = 32
      KEEP_NPOT(I1) = 32
C
      IT_TABLE(I_OH) = 1   ! O-H oxygen
      CALL ICHEM3D (I3D(I_OH),-1,31)   ! OH oxygen atom type = 31
      KEEP_NPOT(I_OH) = 31
C
      IT_TABLE(I_CO) = 1         ! O of carbonyl
      CALL ICHEM3D(I3D(I_CO),-1,15)  ! C=O oxygen = 15 atom type = 15
      KEEP_NPOT(I_CO) = 15
C
      IQ_CO2H = 32    ! param for CO2H C atom
C
      DO 200 I = 1,2
         N = ICON(I,I_OH)
         IF (NAME(N)(1:1) .NE. 'H') GO TO 200
            IT_TABLE(N) = 1
            CALL ICHEM3D(I3D(N),-1,40)
            KEEP_NPOT(N) = 40
      if (phelp) print 1288, (j, it_table(j),j=1,na)
1288  format (' j, it_table(j): ',20i3)
      if (phelp) print 1289, (j, keep_npot(j),j=1,na)
1289  format (' j, keep_npot(j): ',20i3)
      if (phelp) print 1290, IQ_CO2H
1290  format (' IQ_CO2H =',i3)
            RETURN
200   CONTINUE
      RETURN
      END
C------------------------------------------------------------------------
C----Function to locate amide and imide N's and H's and to          ! 10/17/00
C     identify an imide functional group with the structure.....    ! 10/17/00
C                C--C--N--C--C                                      ! 10/17/00
C                   "  |  "                                         ! 10/17/00
C                   O  H  O                                         ! 10/17/00
C
      FUNCTION IQ_AMIDE_H(I1)    ! I1 = 3-linked N                  ! 10/17/00
      CHARACTER NAME*5
      LOGICAL IMIDE                                                 ! 10/17/00
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x        XYZ(3,500), CELL(6), NA, IT_TABLE(500), I3D(500),
     x        KEEP_NPOT(500)
      COMMON /CARBONYL/ N_O_CO
      DIMENSION NUM_H(2), NUM_C(3), NUM_N(3)
C
      IMIDE = .false.                                               ! 10/17/00
      IQ_AMIDE_H = 0
C
      NH = 0
      NC = 0
      NN = 0          ! could be a nitramine N
      DO 10 I = 1,3
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'H') THEN
            NH = NH + 1
            NUM_H(NH) = I2
         ENDIF
         IF (NAME(I2)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            NUM_C(NC) = I2
         ENDIF
         IF (NAME(I2)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            NUM_N(NN) = I2
         ENDIF
10    CONTINUE
c     print 2244, nh, nc, nn   ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c2244  format ('near top of iq_amide_h: nh, nc, nn =',3i3)   !$$$$$$$$$$$$$$$$$$$$$$
C----Nitramine N takes precedence over amide N
      IF (NN .EQ. 1) THEN
         I3 = NUM_N(1)
         NNN = 0
         NNC = 0
         NNO = 0
         DO J=1,NCON(I3)
            L = ICON(J,I3)
            IF (NAME(L)(1:1) .EQ. 'C') NNC = NNC + 1
            IF (NAME(L)(1:1) .EQ. 'O') NNO = NNO + 1
            IF (NAME(L)(1:1) .EQ. 'N') NNN = NNN + 1
         ENDDO
         IF (NCON(I3) .EQ. 3 .AND. NNO .EQ. 2) THEN ! nitramine N, don't
            DO 150 I=1,NC                           !  code as amide...
              N = IQ_CARBONYL3(NUM_C(I))            !  change it back;
150         CONTINUE                                !  done by IQ_CARBONYL3
            RETURN
         ENDIF
      ENDIF
C
      IF ((NH + NC + NN) .NE. 3) RETURN
C
C----Determine if C-C(=O)-NH-C(=O)-C is present                     ! 10/17/00
      IF ((NC + NH) .NE. 3) GO TO 152    ! imide would have 2 C's + 1 H on N ! 10/17/00
      IF (NC .NE. 2) GO TO 152           ! imide has 2 C's on N              ! 10/17/00
      N_CARBONYLS = 0                    ! count C=O's on the N              ! 10/17/00
      DO 154 I=1,NC                                                 ! 10/17/00
         N = IQ_CARBONYL(NUM_C(I),0)                                ! 10/17/00
         IF (N .NE. 0) N_CARBONYLS = N_CARBONYLS + 1                ! 10/17/00
154   CONTINUE                                                      ! 10/17/00
c     print 2222, nn, nc, nh, n_carbonyls              ! $$$$$$$$$$$$$$$$$$$$$$$
c2222  format ('in iq_amide_h: nn,nc,nh,n_carbonyls='4i3) !$$$$$$$$$$$$$$$$$$$$$$$$
      IF (N_CARBONYLS .NE. 2) GO TO 152  ! must be 2 C=O's on the N ! 10/17/0
C-------For this imide func group, must be 1 C, 1 O and 1 N attached to    ! 10/17/00
C        each carbonyl C                                                   ! 10/17/00
      DO 156 I=1,2                                                  ! 10/17/00
          J = NUM_C(I)                                              ! 10/17/00
c         print 2223, i, j, ncon(j)           ! $$$$$$$$$$$$$$$$$$$$$
c2223      format ('in iq_amide_h: i,j,ncon(j)=',4i3) ! $$$$$$$$$$$$$$$$
          IF (NCON(J) .NE. 3) GO TO 152                             ! 10/17/00
          NC1 = 0                                                   ! 10/17/00
          NO1 = 0                                                   ! 10/17/00
          NN1 = 0                                                   ! 10/17/00
          DO K=1,3                                                  ! 10/17/00
             L = ICON(K,J)                                          ! 10/17/00
             IF (NAME(L)(1:1) .EQ. 'C') NC1 = NC1 + 1               ! 10/17/00
             IF (NAME(L)(1:1) .EQ. 'O') NO1 = NO1 + 1               ! 10/17/00
             IF (NAME(L)(1:1) .EQ. 'N') NN1 = NN1 + 1               ! 10/17/00
          ENDDO                                                     ! 10/17/00
          IF ((NC1 .NE. 1) .OR. (NO1 .NE. 1) .OR.                   ! 10/17/00
     X        (NN1 .NE. 1)) GO TO 152                               ! 10/17/00
156   CONTINUE  ! dropping thru 156 loop means have two C-C(=O)-N groups   ! 10/17/00
      IMIDE = .true.                                                ! 10/17/00
152   DO 15 I = 1,NC
         N = IQ_CARBONYL(NUM_C(I),0)
         IF (N .EQ. 0) GO TO 15
           IT_TABLE(NUM_C(I)) = 1
           CALL ICHEM3D(I3D(NUM_C(I)),-1,79)   ! amide carbonyl C, param = 79
           KEEP_NPOT(NUM_C(I)) = 79
C
           IT_TABLE(N_O_CO) = 1
           CALL ICHEM3D(I3D(N_O_CO),-1,80)  ! amide carbonyl O, param = 80
           KEEP_NPOT(N_O_CO) = 80
C
           IF (IMIDE) THEN                                          ! 10/17/00
              IQ_AMIDE_H = 150      ! imide N in C-C(=O)-NH-C(=O)-C ! 10/17/00
           ELSE                                                     ! 10/17/00
              IQ_AMIDE_H = 43       ! amide N...param = 43          ! 10/17/00
           ENDIF                                                    ! 10/17/00
           GO TO 20
15    CONTINUE
      RETURN
C
20    IF (NH .EQ. 0) RETURN    ! any amide H's
      DO 200 I=1,NH
         N = NUM_H(I)
         IT_TABLE(N) = 1
         CALL ICHEM3D(I3D(N),-1,42)   ! amide H, param = 42
         KEEP_NPOT(N) = 42
200   CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C
C----Function to determine the presence of C#N or C#C groups
C
      FUNCTION IQ_CN(I1,ICHANGE)        ! I1 = 2-linked C
      CHARACTER NAME*5, COMMENT*120, REFCOD_2*8
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x        XYZ(3,500), CELL(6), NA, IT_TABLE(500), I3D(500),
     x        KEEP_NPOT(500)
      DIMENSION NUM_N(2), NUM_C(2), NSIZES(10)
      COMMON /PROBLEMS/ JPROB, COMMENT, REFCOD_2, ICSP3, JCSP3
C
      IQ_CN = 0
C
      NN = 0
      NC = 0
C----Is there a N or C attached?  Could be C-C#N or N-C#N.
      DO 101 I = 1,2
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            NUM_N(NN) = I2
         ENDIF
         IF (NAME(I2)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            NUM_C(NC) = I2
         ENDIF
101   CONTINUE
C
      IF (NN .EQ. 0) GO TO 110
C----Check the N's for no. of attachments...N-C#N, C-C#N or N-C#C
         DO 102 I = 1,NN
            I3 = NCON(NUM_N(I))
            IF (I3 .EQ. 1) GO TO 103
102      CONTINUE
         GO TO 110    ! check out C#C, no C#N here
C----It's a CN group...encode it
103      IQ_CN = 45
         IF (ICHANGE .NE. 0) THEN
           IT_TABLE(I1) = 1
           CALL ICHEM3D(I3D(I1),-1,45)      ! C#N nitrile C
           KEEP_NPOT(I1) = 45               ! atom type = 45
C
           IT_TABLE(NUM_N(I)) = 1
           CALL ICHEM3D(I3D(NUM_N(I)),-1,46)! C#N nitrile N
           KEEP_NPOT(NUM_N(I)) = 46         ! atom type = 46
         ENDIF
      RETURN
C----Not nitrile...one of attached C's must be 2-linked
110   DO 120 I=1,NC
         I3 = NUM_C(I)
         N3 = NCON(I3)
         IF (N3 .NE. 2) GO TO 120
C----It's an alkyne...perhaps.  Could be an allene. First determine if part of a ring.
130      D = DISTANCE(I1, I3)         ! calc C-C distance
C      print 7711, refcod_2, i1, i3, name(i1), name(i3), d
C7711  format (' in IQ_CN, refcod_2, i1, i3, name(i1),',
C    x        ' name(i3), d = ',a8,1x,2i3,2(1x,a5),f8.3)
         IF (D .GE. 1.24) GO TO 120   ! C#C must be .le. 1.24 Angs
         CALL BOND_CHECK (I1, NRINGS, NSIZES, NBONDS)
         IF (NBONDS .EQ. 2) THEN
           IB = 55     ! alkyne in a ring
         ELSE
           IB = 53     ! alkyne not in a ring
         ENDIF
C
      IQ_CN = IB
      IF (ICHANGE .NE. 0) THEN
         IT_TABLE(I1) = 1
         CALL ICHEM3D(I3D(I1),-1,IB)     ! C#C alkyne
         KEEP_NPOT(I1) = IB              ! atom type = 53 or 55
C
         IT_TABLE(I3) = 1
         CALL ICHEM3D(I3D(I3),-1,IB)   ! C#C
         KEEP_NPOT(I3) = IB
      ENDIF
      RETURN
120   CONTINUE
      RETURN
      END
C
C------------------------------------------------------------
C    Look for pyrrole-type N's
C
      FUNCTION IQ_PYRROLE_N(NUM_ATOM)  ! NUM_ATOM = 3-linked N
C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x        XYZ(3,500), CELL(6), NA, IT_TABLE(500), I3D(500),
     x        KEEP_NPOT(500)
      DIMENSION NSIZES(10)
C
      IQ_PYRROLE_N = 0
C
      CALL BOND_CHECK(NUM_ATOM, N_INRINGS, NSIZES, NBONDS)
         IF (NEG_COUNT(NSIZES) .EQ. 0) RETURN
      IF (NBONDS .EQ. 2)  IQ_PYRROLE_N = 48 ! pyrrole N in 1 ring; param = 48
      IF (NBONDS .EQ. 3)  IQ_PYRROLE_N = 49 ! pyrrole N in 2 rings; param = 49
      RETURN
      END
C-----------------------------------------------------------------
C----Function to determine the presence of an azo group from a two
C      linked N
C
      FUNCTION IQ_AZO_N(I1)    ! I1 = 2-linked N
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500),KEEP_NPOT(500)
      DIMENSION NUM_C(2), NUM_N(2), NUM_O(2), NSIZES(10)
C
      IQ_AZO_N = 0
      NC = 0
      NN = 0
      NO = 0
C----Enter with a two-linked N
      DO 10 I = 1,2
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'C') THEN   ! check for C connection
            NC = NC + 1
            NUM_C(NC) = I2
         ENDIF
         IF (NAME(I2)(1:1) .EQ. 'N' .AND. NCON(I2) .EQ. 2) THEN   ! check for N connection ! 10/17/00
            NN = NN + 1                                          !  linked to only 2 atoms ! 10/17/00
            NUM_N(NN) = I2
         ENDIF
         IF (NAME(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 2) THEN
            NO = NO + 1                     ! check for O connection,
            NUM_O(NO) = I2                  ! in which O is 2-linked
         ENDIF
10    CONTINUE
c     print 2299, i1, nc, nn, no  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c2299  format ('near top of iq_azo_n: i1,nc,nn,no =',4i2)   ! $$$$$$$$$$$
C----Eliminate if not C-N=N-C or N-N=N-X
      IF (NN .LT. 1) RETURN
C      print 1234, nn, no, nc
C1234  format (' in iq_azo_n: nn, no, nc =',3i3)
C----Azo N must be connected to N and to (C or O)
      IF (.NOT. (NN .EQ. 1) .AND. ((NO + NC) .EQ. 1)) RETURN
      IF (NN .EQ. 1 .AND. KEEP_NPOT(NUM_N(1)) .NE. 0) RETURN  ! 12/28/95
C----Calculate bond distances between N's...must be =< 1.33   ! 10/17/00
      DO 35 I=1,NN
         I3 = NUM_N(I)
         D = DISTANCE(I3,I1)
         IF (D .GT. 1.33) GO TO 35                            ! 10/17/00
            GO TO 40
35    CONTINUE
C     IQ_AZO_N = -999              ! too long N=N             ! 10/17/00
      RETURN
C
40    CALL BOND_CHECK(I1, N_INRINGS, NSIZES, NBONDS)
      IF (NBONDS .EQ. 2 .AND. NEG_COUNT(NSIZES) .NE. 0) RETURN
C
      IQ_AZO_N = 51          ! param = 51
C
      IT_TABLE(I1) = 1               ! diazo N (-N=N)
      CALL ICHEM3D (I3D(I1),-1,51)   ! atom type = 51
      KEEP_NPOT(I1) = 51
C
      IT_TABLE(I3) = 1               ! diazo N (-N=N)
      CALL ICHEM3D (I3D(I3),-1,51)   ! atom type = 51
      KEEP_NPOT(I3) = 51
C
      RETURN
      END
C----------------------------------------------------------------------
C
C----Identity a pyridine N-oxide O or aza-pyrrole N-oxide O
C
      FUNCTION IQ_PYRIDINE_OXIDE(NUM_ATOM) ! NUM_ATOM = oxygen
C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500),KEEP_NPOT(500)
      DIMENSION NSIZES(10)
C
      IQ_PYRIDINE_OXIDE = 0
C
      IN = ICON(1,NUM_ATOM)    ! number of attached N
C
      CALL BOND_CHECK(IN, N_INRINGS, NSIZES, NBONDS) ! check N
      IF (NBONDS .NE. 2) RETURN
      IF (NSIZES(1) .EQ. -6 .OR. NSIZES(1) .EQ. -5) THEN
C
         IT_TABLE(IN) = 1   ! pyridine N-oxide or aza-pyrrole type N
         CALL ICHEM3D (I3D(IN),-1,18)   ! atom type = 18
         KEEP_NPOT(IN) = 18
C
         IT_TABLE(NUM_ATOM) = 1   ! pyridine N-oxide O
         CALL ICHEM3D (I3D(NUM_ATOM),-1,17)  ! atom type = 17
         KEEP_NPOT(NUM_ATOM) = 17
C
         IQ_PYRIDINE_OXIDE = 17          ! param = 17 for O
      ENDIF
      RETURN
      END
C------------------------------------------------------------
C----Function to determine the presence of an azoxy group
C     staring with a single-linked O.
C     -N(O)=N- or -N(O)=N(O)- or even -O-N=N(O)-
C
      FUNCTION IQ_AZOXY(I1)   ! I1 = 1-linked O
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500), KEEP_NPOT(500)
      DIMENSION NUM_C(3), NSIZES(10), NUM_N(3), NUM_O(3)
C
      IQ_AZOXY = 0
C
C----Enter with a single-linked O...I1; N-O = I2-I1
      I2 = ICON(1,I1)    ! is it linked to a N?
      IF (NAME(I2)(1:1) .NE. 'N') RETURN
c
C----Check N (I2) for connections to other N and C...O(I1)-N(I2)
15    NN = 0
      NC = 0
      DO 20 I = 1,NCON(I2)
         I3 = ICON(I,I2)
         IF (NAME(I3)(1:1) .EQ. 'C') THEN
            NC = NC + 1
            NUM_C(NC) = I3
         ENDIF
         IF (NAME(I3)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            NUM_N(NN) = I3
         ENDIF
20    CONTINUE
C----Eliminate if N not attached to one N and one C
      IF (NN .NE. 1) RETURN
      IF (NC .NE. 1) RETURN
C
C----Return if linked N is part of an aromatic ring.
      CALL BOND_CHECK (I2, N_INRINGS, NSIZES, NBONDS)
      IF (NBONDS .EQ. 2 .AND. NEG_COUNT(NSIZES) .NE. 0) RETURN
C
      NN2 = 0
      NC = 0
      NO_1LINK = 0
      NO_2LINK = 0
C----Is 2nd N linked to 1 C & 1 N, 1 N & 1 O or 2 N's ?
C     -C(O)-N=N(O)-, -C-N(O)=N(O)-, -N-N=N(O)-
      DO 25 I = 1,NCON(NUM_N(1))
        I4 = ICON(I,NUM_N(1))
        IF (NAME(I4)(1:1) .EQ. 'C') NC = NC + 1
        IF (NAME(I4)(1:1) .EQ. 'N') NN2 = NN2 + 1
        IF (NAME(I4)(1:1) .EQ. 'O' .AND. NCON(I4) .EQ. 1) THEN
              NO_1LINK = NO_1LINK + 1     ! count # of 1-linked O's
              NUM_O(NO_1LINK) = I4        ! on 2nd N
        ENDIF
        IF (NAME(I4)(1:1) .EQ. 'O' .AND. NCON(I4) .EQ. 2)
     X        NO_2LINK = NO_2LINK + 1  ! count # of 2-linked O's on 2nd N
25    CONTINUE
C----Return if second N not attached to 1 N, 1 C, 1 O or 2 N's
      IF (NN2 .LT. 1) RETURN
      IF (.NOT. (NC .EQ. 1 .OR. NO_2LINK .EQ. 1 .OR. NN2 .EQ. 2)) RETURN
C
      IQ_AZOXY = 56     ! O param = 56
C
      IT_TABLE(I1) = 1               ! azoxy O
      CALL ICHEM3D (I3D(I1),-1,56)   ! atom type = 56
      KEEP_NPOT(I1) = 56
C
      IT_TABLE(I2) = 1               ! N linked to the O
      CALL ICHEM3D (I3D(I2),-1,57)   ! atom type = 57
      KEEP_NPOT(I2) = 57
C
C----Handle 2nd N and its O (if there is one)
      IF (NO_1LINK .EQ. 1) THEN          ! -N(O)=N(O)-
         N = NUM_O(1)
         IT_TABLE(N) = 1               ! azoxy O
         CALL ICHEM3D (I3D(N),-1,56)   ! atom type = 56
         KEEP_NPOT(N) = 56
C
         N = NUM_N(1)
         IT_TABLE(N) = 1               ! N linked to the O
         CALL ICHEM3D (I3D(N),-1,57)   ! atom type = 57
         KEEP_NPOT(N) = 57
      ELSE
         N = NUM_N(1)                 !-N(O)=N-
         IT_TABLE(N) = 1               ! N not linked to an O
         CALL ICHEM3D (I3D(N),-1,51)   ! atom type = 51
         KEEP_NPOT(N) = 51
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
C
      FUNCTION IQ_ONO2(NUM_ATOM)  ! nitrate O; O-no2
C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500), KEEP_NPOT(500)
      DIMENSION NUM_N(3)
C
      IQ_ONO2 = 0   ! initialize to NO !
C
      NC = 0
      NN = 0
      DO 100 I=1,NCON(NUM_ATOM)
          I1 = ICON(I,NUM_ATOM)
          IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
          IF (NAME(I1)(1:1) .EQ. 'N') THEN
              NN = NN + 1
              NUM_N(NN) = I1
          ENDIF
100   CONTINUE
C
      IF (.NOT. (NC .EQ. 1 .AND. NN .EQ. 1)) RETURN
C
      NN = NCON(NUM_N(1))   ! # of N attached to NUM_ATOM
      IF (NN. NE. 3) RETURN
      NO = 0
      DO 200 I=1,3
          I1 = ICON(I,NUM_N(1))
          IF (NAME(I1)(1:1) .EQ. 'O') NO = NO + 1
200   CONTINUE
C
      IF (NO .EQ. 3) IQ_ONO2 = 60  ! O-NO2; param = 60
C
      RETURN
      END
C---------------------------------------------------------------
C----Function to determine if a 2-linked C is a central
C     allene C (C=C=C) or N=C=O or C=C=N
C
      FUNCTION IQ_ALLENE_2(NUM_ATOM, ICHANGE)
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2), NUM_N(2), NUM_O(2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
C
      IQ_ALLENE_2 = 0
      ICENTRAL = NUM_ATOM
C
C----How many C's bonded to NUM_ATOM...must be 2 or 1 C & 1 N
      NN = 0           ! look for C=C=C, C=C=N or N=C=O
      NC = 0
      NO = 0
      DO 20 I=1,2
          I4 = ICON(I,ICENTRAL)
          IF (NAME(I4)(1:1) .NE. 'C') GO TO 19
          IF (KEEP_NPOT(I4) .EQ. 53 .OR.
     X        KEEP_NPOT(I4) .EQ. 55) GO TO 20
C---------Check to determine if I4 is an alkyne C
             J = IQ_CN(I4,0)
             IF (J .EQ. 53 .OR. J .EQ. 55) RETURN
             NC = NC + 1
             IND(NC+NN+NO) = I4
19        IF (NAME(I4)(1:1) .EQ. 'N' .AND. NCON(I4) .EQ. 2) THEN
             NN = NN + 1        ! if N, must be 2-linked
             IND(NC+NN+NO) = I4
             NUM_N(NN) = I4
          ENDIF
          IF (NAME(I4)(1:1) .EQ. 'O' .AND. NCON(I4) .EQ. 1) THEN
             NO = NO + 1        ! if O, must be 1-linked
             IND(NC+NN+NO) = I4
             NUM_O(NO) = I4
          ENDIF
20    CONTINUE
      IF ((NC + NN + NO) .LT. 2) RETURN
C----Look for C-C-C bond angle .gt. 160 deg
            DO 40 L=1,2           ! calc the I2 - I1 & I3 - I1
               XL(L) = 0.0        ! vectors & their lengths
               DO 38 K=1,3
                  VEC(K,L) = XYZ(K,IND(L)) - XYZ(K,ICENTRAL)
                  XL(L) = XL(L) + VEC(K,L)**2
38             CONTINUE
               XL(L) = SQRT(XL(L))
40          CONTINUE
            D = 0.0
            DO 50 K=1,3
               D = D + VEC(K,1)*VEC(K,2)
50          CONTINUE
            COSINE = D/(XL(1)*XL(2))
            IF (COSINE .GE. -0.939693) RETURN  ! cos(160) = -0.939693
C
            IQ_ALLENE_2 = 65
            IF (ICHANGE .EQ. 0) RETURN
            IF (IT_TABLE(ICENTRAL) .NE. 0) GO TO 60
            IT_TABLE(ICENTRAL) = 1   ! central allene-like C (sp)
            CALL ICHEM3D (I3D(ICENTRAL),-1,65)   ! atom type = 65
            KEEP_NPOT(ICENTRAL) = 65
C
60          IF (NO .EQ. 0) RETURN     ! =C=O ?
            N = NUM_O(1)
            IT_TABLE(N) = 1                ! O or C=C=O or C=C=O
            CALL ICHEM3D (I3D(N),-1,81)    ! atom type = 81
            KEEP_NPOT(N) = 81
      RETURN
      END
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C----An end allene C ?
C
      FUNCTION IQ_ALLENE_END(IEND)  ! IEND = # of 3-linked C
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2), NUM_N(2), NUM_O(2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
C
      IQ_ALLENE_END = 0
C
C----Is there a 2-linked C connected to IEND
      IND(1) = IEND
      DO 200 I=1,3
         ICENTRAL = ICON(I,IEND)     ! look for a =C=
         IF (ICENTRAL .EQ. IEND) GO TO 200
         IF (NCON(ICENTRAL) .NE. 2) GO TO 200
         IF (NAME(ICENTRAL)(1:1) .NE. 'C') GO TO 200
         IF (KEEP_NPOT(ICENTRAL) .EQ. 53) GO TO 200
         IF (KEEP_NPOT(ICENTRAL) .EQ. 55) GO TO 200
         GO TO 210 ! passed all tests...not an alkyne
200   CONTINUE
      RETURN
C
210   NTEST = IQ_ALLENE_2(ICENTRAL,0)
      IF (NTEST .NE. 65) RETURN
            IQ_ALLENE_END = 64
            IT_TABLE(IEND) = 1   ! terminal allene R2C= (R2C=C=CR2)
            CALL ICHEM3D (I3D(IEND),-1,64)   ! atom type = 64
            KEEP_NPOT(IEND) = 64
      RETURN
      END
C
C------------------------------------------------------------------
C----Function to determine if a 2-linked N is an end N in a
C    C-N(-)-N(+)=N-C
C           |
C           C
C
      FUNCTION IQ_NNN_END(I1)   ! I1 = 2-linked N
      CHARACTER NAME*5
      DIMENSION  NUM_N(2), NUM_N2(3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
C
C                       Z-N(-)-N(+)=N-Z   Z = C or (R2)N
C                              |
C                              C
C
      IQ_NNN_END = 0
C
C----How many C's & N's bonded to I1...must be 1 C & 1 N or
C     2 N's (one 2-linked and one 3-linked)
5     NN = 0         ! count number of 3-linked N's
      NC = 0
      NTOTAL = 0
      DO 20 I=1,2
          I4 = ICON(I,I1)
          IF (NAME(I4)(1:1) .EQ. 'C') NC = NC + 1
          IF (NAME(I4)(1:1) .EQ. 'N') NTOTAL = NTOTAL + 1
          IF (NAME(I4)(1:1) .EQ. 'N' .AND. NCON(I4) .EQ. 3) THEN
             NN = NN + 1        ! if N, must be 3-linked
             NUM_N(NN) = I4
          ENDIF
20    CONTINUE
      IF (NC .EQ. 1 .AND. NN .EQ. 1) GO TO 24
      IF (NTOTAL .EQ. 2 .AND. NN .GE. 1) GO TO 24
         RETURN
C----Now C-I1-I4 or N-I1-I4...may need to check 2 possible I4's
24    DO 28 L=1,NN
         I4 = NUM_N(L)
         NTEST = IQ_NNN_MIDDLE(I4)
         IF (NTEST .EQ. 82) GO TO 29
28    CONTINUE
      RETURN      ! dropping thru loop = can't properly id I4
C----Is there another 2-linked N (other than I1) bonded to I4
29    NN = 0         ! N must be 2-linked
      DO 35 I=1,3
         I3 = ICON(I,I4)
         IF (I3 .EQ. I1) GO TO 35  ! we know about I1...1st N
         IF (NAME(I3)(1:1) .EQ. 'N' .AND. NCON(I3) .EQ. 2) THEN
            NN = NN + 1
            NUM_N(NN) = I3
         ENDIF
35    CONTINUE
      IF (NN .EQ. 1) IQ_NNN_END = 83
      RETURN
      END
C------------------------------------------------------------
C----Function to determine if a 3-linked N is a middle N in a
C
C    Z-N(-)-N(+)=N-Z      Z = C or (R2)N
C           |
C           C
C
      FUNCTION IQ_NNN_MIDDLE(I5)    ! I5 = 3-linked N
      CHARACTER NAME*5
      DIMENSION  NUM_N(2), NUM_N2(3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500),KEEP_NPOT(500)
C
      IQ_NNN_MIDDLE = 0
C
55    NC2 = 0
      NN2 = 0
      DO 40 I=1,3
         N = ICON(I,I5)
         IF (NAME(N)(1:1) .EQ. 'C') NC2 = NC2 + 1
         IF (NAME(N)(1:1) .EQ. 'N' .AND. NCON(N) .EQ. 2) THEN
            NN2 = NN2 + 1
            NUM_N2(NN2) = N
         ENDIF
40    CONTINUE
      IF (NC2 .EQ. 1 .AND. NN2 .EQ. 2) GO TO 42
         RETURN
C----Check the two connected N's...
C     must both be 2-linked & have 1 link to C or N
42    DO 100 I=1,2
         N = NUM_N2(I)
         IF (NAME(N)(1:1) .EQ. 'N' .AND. NCON(N) .EQ. 2) GO TO 44
            RETURN
44          MC = 0     ! check the 2 N's linked to I5
            N3N = 0    ! number of 3-linked N's to the end N
            DO 95 J=1,2
               L = ICON(J,N)
               IF (NAME(L)(1:1) .EQ. 'C') MC = MC + 1  ! count # C's
               IF (NAME(L)(1:1) .EQ. 'N' .AND.         ! count # of
     x             NCON(L) .EQ. 3) N3N = N3N + 1       !  3-linked N's
95          CONTINUE
            IF ((MC .EQ. 1 .AND. N3N .EQ. 1) .OR.  ! must be linked to...
     x          (N3N .EQ. 2)) GO TO 100            ! (1) 1 C & one 3-linked N or
               RETURN                              ! (2) two 3-linked N's
100   CONTINUE
C
      IQ_NNN_MIDDLE = 82
      RETURN
      END
C-----------------------------------------------------------------
C----Function to determine the presence of a peroxide group from
C             a two-linked O
C
      FUNCTION IQ_PEROXY(I1)
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_PEROXY = 0   ! 0 means not a peroxide
C
C----How many O's bonded to I1...must be at least 1
      NC = 0
      NO = 0
      DO 20 I=1,2
          I2 = ICON(I,I1)
          IF (NAME(I2)(1:1) .EQ. 'C') NC = NC + 1
          IF (NAME(I2)(1:1) .EQ. 'O') NO = NO + 1
20    CONTINUE
      IF (NO .NE. 1) RETURN
      IQ_PEROXY = 66  ! peroxide O
      RETURN
      END
C--------------------------------------------------------------
C----Function to determine a N is part of any kind of
C     nitro group
C     is cubane-linked; also does the same of the N.
C
      FUNCTION IQ_ANY_NO2(IT)    ! IT is the 3-linked N
      CHARACTER NAME*5
      DIMENSION  IND(2), XL(2), VEC(3,2)
      EQUIVALENCE (IND(1),I2), (IND(2),I3)
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_ANY_NO2 = 0     ! NO
C
      NO = 0
      NC = 0
      NN = 0
      DO 100 I=1,3
         I1 = ICON(I,IT)
         IF (NAME(I1)(1:1) .EQ. 'O') NO = NO + 1
         IF (NAME(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (NAME(I1)(1:1) .EQ. 'N') NN = NN + 1
100   CONTINUE
C
      IF (NO .LT. 2) RETURN   ! doesn't link to at least 2 O's
      IF (NO .EQ. 3) IQ_ANY_NO2 = 1   ! -O-NO2
      IF (NC .EQ. 1 .AND. NO .EQ. 2) IQ_ANY_NO2 = 2  ! C-NO2
      IF (NN .EQ. 1 .AND. NO .EQ. 2) IQ_ANY_NO2 = 3  ! N-NO2
      RETURN
      END
C--------------------------------------------------------------
C----Function to determine if a 2-linked O is an epoxide
C
      FUNCTION IQ_EPOXIDE(I1)   ! I1 is a 2-linked O
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      DIMENSION NUM_C(2)
C
      IQ_EPOXIDE = 0     ! NO
C
      DO 20 I=1,2
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .NE. 'C') RETURN
20    CONTINUE
C----The 2 C's linked to O must be linked to each other
      I3 = ICON(1,I1)
      N3 = NCON(I3)
      I4 = ICON(2,I1)
      DO 40 I=1,N3
         IF (I4 .EQ. ICON(I,I3)) GO TO 60
40    CONTINUE
      RETURN
60    IQ_EPOXIDE = 76
      RETURN
      END
C--------------------------------------------------------------
C----Function to determine if a 2-linked O is linked to 2 C's
C     one of which is a carbonyl...like an ester or anhydride
C
      FUNCTION IQ_ESTER(I1)   ! I1 is a 2-linked O
      CHARACTER NAME*5
      LOGICAL CARBONYL
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      DIMENSION NSIZES(10)
C
      IQ_ESTER = 0     ! NO
C
      CARBONYL = .false.
      DO 20 I=1,2
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .NE. 'C') RETURN
         IF (IQ_CARBONYL(I2,0) .NE. 0) CARBONYL = .true.  ! test for adjacent C=O
20    CONTINUE
      IF (.NOT. CARBONYL) RETURN       ! no adjacent C=O's?...return
C----The 2 C's linked to O...at least one is a carbonly.  Do ring test
      CALL BOND_CHECK(I1, N_INRINGS, NSIZES, NBONDS)
      GO TO (30,30,40), NBONDS + 1
C             0  ?  2           <--- # bonds in rings
40       IQ_ESTER = 94     ! carbonyl linked O in a ring
         RETURN
30       IQ_ESTER = 95     ! carbonyl linked O not in a ring
         RETURN
      END
C
C--------------------------------------------------------------
C----Function to determine if a 4-linked C that is not in a
C     rings is CH3, CH2 or CH.  If CH3, is it to C, to N, to O    ! 3/24/00
C     or to 2-linked S?                                           ! 3/24/00
C
      FUNCTION IQ_CH3CH2CH(I1)   ! I1 is a 4-linked C
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      IQ_CH3CH2CH = 0     ! NO
C
      NH = 0
      I_NOS = 0                                                   ! 3/24/00
      DO 20 I=1,4
         I2 = ICON(I,I1)
         IF (NAME(I2)(1:1) .EQ. 'H') NH = NH + 1
         IF (NAME(I2)(1:1) .EQ. 'O') I_NOS = 2                    ! 3/24/00
         IF (NAME(I2)(1:1) .EQ. 'N') I_NOS = 1                    ! 3/24/00
         IF (NAME(I2)(1:1) .EQ. 'S' .AND. NCON(I2) .EQ. 2)        ! 3/24/00
     X              I_NOS = 3                                     ! 3/24/00
20    CONTINUE
      IF (NH .EQ. 0) THEN
30       IQ_CH3CH2CH = 23  ! quaternary C...no H's
         RETURN
      ENDIF
      IF (NH .EQ. 3 .AND. I_NOS .NE. 0)                           ! 3/24/00
     x       GO TO (31,32,33),I_NOS                               ! 3/24/00
C------C of non-ring CH, CH2 or CH3-C                             ! 4/11/00
        IQ_CH3CH2CH = 143 + NH  ! params = 144, 145, 146 for C's  ! 4/11/00
           RETURN               !   with 1, 2 & 3 H's             ! 4/11/00
C------CH3-N                                                      ! 3/24/00
31      IQ_CH3CH2CH = 147       ! C of CH3-N                      ! 4/11/00
           RETURN                                                 ! 3/24/00
C------CH3-O                                                      ! 3/24/00
32      IQ_CH3CH2CH = 148       ! C of CH3-O                      ! 4/11/00
           RETURN                                                 ! 3/24/00
C------CH3-S                                                      ! 3/24/00
33      IQ_CH3CH2CH = 149       ! C of CH3-S                      ! 4/11/00
           RETURN                                                 ! 3/24/00
      END
C
C--------------------------------------------------------------------
C----Calculate distances between 2 atoms
C
      FUNCTION DISTANCE (I1, I2)
      CHARACTER LABEL*5
C
      COMMON /COM1/ NCON(500), ICON(5,500), LABEL(500),
     x                 OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      DISTANCE = 0.0
      DO I=1,3
         DISTANCE = DISTANCE + (OR_XYZ(I,I1) - OR_XYZ(I,I2))**2
      ENDDO
      DISTANCE = SQRT(DISTANCE)
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE ICHEM3D(I3D,NUMAT,J)
C
      PARAMETER (N = 132)                                         ! 6/9/05
      CHARACTER COMMENT*120, REFCOD_2*8
      DIMENSION ITYPE_3D(N), IPOT(N)
C----IPOT = atom type from subprogram NPOT.
C    ITYPE = Chem3D atom type
      DIMENSION IT1(8), IT2(8), IP1(8), IP2(8), IT3(6), IP3(6),
     x          IT4(7), IP4(7), IT5(8), IP5(8), IT6(6), IP6(6),
     x          IT7(7), IP7(7), IT8(7), IP8(7), IT9(7), IP9(7),
     x          IT10(7), IP10(7), IT11(7), IP11(7),
     x          IT12(7), IP12(7), IT13(7), IP13(7),               ! 3/5/00
     x          IT14(7), IP14(7), IT15(7), IP15(7),               ! 3/22/00
     x          IT16(7), IP16(7), IT17(4), IP17(4),               ! 4/11/00
     x          IT18(7), IP18(7), IT19(8), IP19(8)                ! 6/9/05
      EQUIVALENCE (IT1,ITYPE_3D), (IT2,ITYPE_3D(9)),
     x            (IT3,ITYPE_3D(17)), (IT4,ITYPE_3D(23)),
     x            (IT5,ITYPE_3D(30)), (IT6,ITYPE_3D(38)),
     x            (IT7,ITYPE_3D(44)), (IT8,ITYPE_3D(51)),
     x            (IT9,ITYPE_3D(58)), (IT10,ITYPE_3D(65)),
     x            (IT11,ITYPE_3D(72)), (IT12,ITYPE_3D(79)),
     x            (IT13,ITYPE_3D(86)), (IT14,ITYPE_3D(93)),       ! 3/4/00
     x            (IT15,ITYPE_3D(100)), (IT16,ITYPE_3D(107)),     ! 3/22/00
     x            (IT17,ITYPE_3D(114)), (IT18,ITYPE_3D(118)),     ! 4/11/00
     x            (IT19,ITYPE_3D(125)),                           ! 11/3/00
     x            (IP1,IPOT),(IP2,IPOT(9)),
     x            (IP3,IPOT(17)), (IP4,IPOT(23)),
     x            (IP5,IPOT(30)), (IP6,IPOT(38)),
     X            (IP7,IPOT(44)), (IP8,IPOT(51)),
     x            (IP9,IPOT(58)), (IP10,IPOT(65)),
     x            (IP11,IPOT(72)), (IP12,IPOT(79)),
     x            (IP13,IPOT(86)), (IP14,IPOT(93)),               ! 3/4/00
     x            (IP15,IPOT(100)), (IP16,IPOT(107)),             ! 3/22/00
     x            (IP17,IPOT(114)), (IP18,IPOT(118)),             ! 11/3/00
     x            (IP19,IPOT(125))                                ! 11/3/00
C
      COMMON /PROBLEMS/ JPROB, COMMENT, REFCOD_2, ICSP3, JCSP3
      DATA IT1 / 1,    2,     5,    5,    1,      1,    1,    2/
      DATA IP1 / 5,    4,     1,    2,   20,     21,   22,   23/
C----           cub-  Csp2  H to    H   Csp3    Csp3   Csp3  Csp3
C               ane   0 in  sp2 C       2 in    3 in   4 in   no
C                C    ring  nonarom     ring    ring   ring  ring
C
      DATA IT2 / 2,    2,    46,    47,     8,     8,     8,     7 /
      DATA IP2 /24,   25,     7,     6,     8,    27,     9,    15 /
C              C=C   Csp2   nitro  nitro  N-NO2  N-NO2  amine carbonyl    ! 10-23-96
C             C(sp2) 3 in     N      O      N      N      N      O        ! 10-23-96
C            not in  ring                not in  ring   no rng            ! 10-23-96
C             ring                        ring                            ! 10-23-96
C
      DATA IT3 /21,    6,     3,       2,        2,         2/
      DATA IP3 /30,   31,    32,      33,       34,        35/
C               H     O   carbonyl aromatic  aromatic  aromatic
C               of    of   C not    2 bonds   3 rg bnds 3 rg bnds
C               O-H   O-H  in rg    in rings  2 aromat  3 aromat
      DATA IT4 /  5,     47,      8,    23,     6,    46,    45 /   ! all azides N's
      DATA IP4 / 36,     11,     10,     3,    14,    37,    39 /   ! have same Chem3d
C               H to    O of   amine   amine  ether nitro  azide    ! atom type
C               cubane  cubane   N       H      O   N to    N's
C               ring    nitro          (N-H)        cubane
C
      DATA IT5 / 24,   28,   9,     6,    37,     40,    2,    10  /
      DATA IP5 / 40,   42,  43,    44,    16,     47,    45,   46  /
C              -CO2H amide amide ether  pyrid or sp2 N  nit-  nit-
C                H     H     N    O in  aza-pyr   in    rile  rile
C                                 ring    N      -C=N-   N     C
C
      DATA IT6 /  40,     40,      11,     40,      37,    41 /
      DATA IP6 /  48,     49,      13,     50,      51,    52 /
C               pyrrole pyrrole    F    sp2 N in   azo    furan
C                N in    N in     in    -C=N- in    N       O
C               1 ring  2 rings   N-F    1 ring
C
      DATA IT7 /   4,      5,      4,     81,      50,     37,   81   /
      DATA IP7 /  53,     54,     55,     17,      18,     56,   57  /
C                sp C    H of    sp C   pyridine pyridine    azoxy
C               alkyne  alkyne  alkyne  N-oxide  N-oxide    O     N
C               no rng          in rng     O        N      -N(O)=N-
C
      DATA IT8 /   8,      8,     99,     7,    40,    2,    4  /
      DATA IP8 /  58,     59,     60,    61,    62,   64,   65  /
C                 tertiary N     O of     nitroso      allene
C                1 rng   2 rngs  c-O-no2   O=N-         C=C=C
C                                         O      N    sp2   sp
C
      DATA IT9 /    6,     6,    21,    45,    23,     2,    49 /
      DATA IP9 /   66,    67,    68,    69,    70,    71,    76 /
C               peroxide   O      H   diazo  amine  diazo  epoxide
C                -O-O-     of Csp2    C=N=N  oxide  C=N=N     O
C                  O         O-H        N's    N      C
C
      DATA IT10 /   3,    7,     7,    50,      40,      8,      8 /
      DATA IP10 /  79,   80,    81,    82,      83,     85,     86 /
C                  C     O    O in   middle & end N's  amine    C3N
C                 of amide,   C=C=O  C-N(-)-N(+)=N,    bonded   bonded
C                urea, imide    &           |           to       to
C                   C=O       N=C=O         C          C(sp2)   C(sp2)
C
      DATA IT11 /   40,      11,     11,      8,     6,     6,      3/
      DATA IP11 /   19,      87,     88,     89,    94,    95,     96/
C                tetraaza    F to    F to   amine     ester or   carbonyl
C                pentalene   C(sp2)  C(sp3)   N     anhydride O   C in a
C                                                   in     not in  ring
C                                                  ring     ring
C
      DATA IT12 /   1,     1,     1,     1,     1,     1,     1/
      DATA IP12 /  98,    99,   100,   101,   102,   103,   104/
C                   3      4      5      6      7      8 .ge. 9
C                 ...Csp3 with 2 bonds in above sized rings...
C
      DATA IT13 /   1,     1,     1,     5,      46,     2,   13/     ! 2/26/00
      DATA IP13 / 105,   106,   107,   108,     110,   112,   12/     ! 2/26/00
C                  CH    CH2    CH3  H linked   N of   C=C    Br      ! 2/26/00
C                  not in rings C's  to arom C  C(sp2) (sp2)  to      ! 3/17/00
C                                               nitro  in     C(sp3)  ! 3/17/00
C                                                      ring           ! 10-23-96
      DATA IT14 /   12,     12,     3,  161,   12,   17,  81/         ! 3/19/00
      DATA IP14 /  113,    114,   115,  116,  118,  119, 120/         ! 3/19/00
C                   CL      CL     C     S     Br     S    O          ! 3/19/00
C                   to      to        of       to       in            ! 3/19/00
C                 C(sp3)  C(sp2)     C=S     C(sp2)  sulfoxide        ! 3/19/00
      DATA IT15 /   18,  81,  18,   7,   6,   44,    15/              ! 3/22/00
      DATA IP15 /  122, 123, 125, 126, 127,  130,   131/              ! 3/22/00
C                   S    O    S    O2   O     H      S                ! 3/22/00
C                     in           in            in                   ! 3/22/00
C                  sulfone    sulfonates       thiol                  ! 3/22/00
C
      DATA IT16 /   42,    15,      1,     1,     1,   15,   15/      ! 3/27/00
      DATA IP16 /  129,   133,    134,   135,   136,  137,  138/      ! 3/27/00
C                   S      S      CH3    CH3    CH3    C     S        ! 3/27/00
C                 thio-   sul-     |      |      |     S*    S*       ! 3/27/00
C                 phene   fide     N      O      S     S     S        ! 3/27/00
C
      DATA IT17 /   18,    81,      8,     25/                        ! 4/4/00
      DATA IP17 /  140,   139,    141,    143/                        ! 4/4/00
C                   S      O       N       P                          ! 4/4/00
C                   of     of      of      of                         ! 4/4/00
C                     sulfonamide         C3P                         ! 4/4/00
C
      DATA IT18 /    1,   1,   1,   1,   1,   1,   9/                ! 11/3/00
      DATA IP18 /  144, 145, 146, 147, 148, 149, 150/                ! 11/3/00
C                   C    C    C    C    C    C    N                  ! 10/17/00
C                  of   of   of   of   of   of   of                  ! 10/17/00
C                  CH   CH2  CH3  CH3  CH3  CH3 imide                ! 10/17/00
C                   not in    |    |    |    |                       ! 10/17/00
C                   ring      C    N    O    S                       ! 10/17/00
      DATA IT19 /   18,   81,   40,    2,  170,  350,  530,   39/    ! 6/9/05
      DATA IP19 /  152,  153,  154,  155,  164,  165,  166,  158/    ! 6/9/05
C                   O     S     N     C    Cl-   Br-   I-     N+     ! 6/9/05
C                  of    of    of    of                       of     ! 6/9/05
C                      sulfonimine                         ammonium  ! 6/9/05
C                       SO2-N=C                             cation   ! 6/9/05
      I3D = 0
      IF (NUMAT .LE. 0) GO TO 90
C      write (36,2288) numat                                ! 1$$$$$$$$$$$
C2288  format ('top of ichem3d before npot: numat =',i2)    ! 1$$$$$$$$$$$
      jj = numat
      j = npot(jj)
c     J = NPOT(NUMAT)
c      print 2266, numat, j                                  ! 1$$$$$$$$$$$$$$$$
c2266  format ('top of ichem3d after npot: numat, j =',2i4)  ! 1$$$$$$$$$$$$$$$
      IF (J .EQ. 999) RETURN
90    DO 100 I = 1,N
         IF (J .NE. IPOT(I)) GO TO 100
         I3D = ITYPE_3D(I)
         RETURN
100   CONTINUE
      WRITE (COMMENT, 110) J
      PRINT 110, J
110   FORMAT (' $$$$ PROBLEM: No Chem3D type corresponding to atom',
     x        ' type',I3,' in table')
      IF (JPROB .NE. 0) CALL WRFILE(REFCOD_2)
      RETURN
      END
C----------------------------------------------------------------
      BLOCK DATA STUFF
      PARAMETER (NTYPES_ATOMS=11)                                 ! 6/7/05
      CHARACTER ATSYM*1, ATSYM2*2                                 ! 2/26/00
      COMMON /CONTENTS/ NCON_FDAT(NTYPES_ATOMS),                  ! 2/26/00
     x                  NCON_FBIB(NTYPES_ATOMS),                  ! 2/26/00
     x                  ATSYM(NTYPES_ATOMS)                       ! 2/26/00
      DATA ATSYM /'C','H','F','N','O','X','Z','S','P','B','I'/    ! 6/7/05
C                                      BR  CL                     ! 3/4/00
      COMMON /ELEMENTS/ N_AT_DIFF, N_PER_UNIT(8), ATSYM2(NTYPES_ATOMS) ! 3/6/00
      DATA N_PER_UNIT, N_AT_DIFF /9*0/                            ! 3/6/00
      DATA ATSYM2 /'C ', 'H ', 'F ', 'N ', 'O ', 'BR', 'CL', 'S ',! 3/11/00
     x             'P ', 'B ', 'I '/                              ! 6/7/05
      END
C---------------------------------------------------------------------
      SUBROUTINE BOND_CHECK (NUM_ATOM, NRINGS, NSIZES, NBONDS)
C
      PARAMETER (NENTRIES = 1500, MAXSTEPS=33, NUM_LEVELS=MAXSTEPS+2)
      CHARACTER NAME*5, COMMENT*120, SPACER*8
      logical phelp, pphelp
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      COMMON /PROBLEMS/ JPROB, COMMENT, SPACER, ICSP3, JCSP3
      DIMENSION KEEP(3,NENTRIES,NUM_LEVELS), NKEEP(NENTRIES),
     x          NSIZES(10), NTEMP(60,60), NICT(60), NTRYS(500)
      DIMENSION NCON2(500), ICON2(4,500), ILAST(4), N_IN_NTEMP(20)
C
      phelp = .false.
      pphelp = .false.
C
      if (pphelp) then
         print *, ' +++++++++top of BOND_CHECK++++++++++++++'
         do i=1,na
            print 1251, i, name(i), (icon(j,i),j=1,ncon(i))
1251        format (i3,2x,a5,4i4)
         enddo
      endif
C
      N_PATHS = 0        ! initialize # of possible ring paths
      DO I=1,20
         NICT(I) = 0     ! # of atoms in each ring path
      ENDDO
      DO I=1,10
         NSIZES(I) = 0    ! size of ring for each bond to NUM_ATOM
      ENDDO
      NRINGS = 0          ! # of rings examined for NUM_ATOM
C
      if (phelp) then
         print 1240, num_atom, name(num_atom),
     x               (icon(i,num_atom),i=1,ncon(num_atom))
1240     format (' ---------------------------------------------'/
     x           ' START: num_atom =',i3,' (',a5,')',' icon =',7i3)
      endif
C
C----Find the last atom before the H's
      DO 10 I=1,NA
         IF (NAME(I)(1:1) .EQ. 'H') GO TO 15
10    CONTINUE
      LAST = I
      GO TO 20    ! falling thru loop means no H's
15    LAST = I - 1
C
C----Create working arrays (NCON2 & ICON2) w/o H's
20    DO 200 I=1,LAST
         N = NCON(I)
         NCON2(I) = 0
         DO 190 J=1,N
            IF (ICON(J,I) .GT. LAST) GO TO 200
            NCON2(I) = NCON2(I) + 1         ! transfer info to NCON2 &
            ICON2(NCON2(I),I) = ICON(J,I)   ! ICON2 arrays
190      CONTINUE
200   CONTINUE
C
      if (phelp) then
         print 1241, (icon2(i,num_atom),i=1,ncon2(num_atom))
1241     format (28X,' icon2 =',7i3)
      endif
C
      N1 = NCON2(NUM_ATOM)         ! # of atoms linked to NUM_ATOM
      IF (N1 .LT. 2) GO TO 5000    ! < 2, cannot be in a ring
C
      DO I=1,500
         NTRYS(I) = 0
      ENDDO
      NSTEPS = 1
      NKEEP(NSTEPS) = 0
C
C----Place the atoms connected to NUM_ATOM into the KEEP arrays
      DO 30 L1=1,N1
         IB = ICON2(L1,NUM_ATOM)  ! start with IB...linked to NUM_ATOM
	 IF (NCON2(IB) .EQ. 1) GO TO 30          ! IB can't be in a ring
         NKEEP(NSTEPS) = NKEEP(NSTEPS) + 1
         KEEP(1,NKEEP(NSTEPS),NSTEPS) = IB        ! IB is first atom in the pyramid...
         KEEP(2,NKEEP(NSTEPS),NSTEPS) = NUM_ATOM  ! IB came from NUM_ATOM?
         KEEP(3,NKEEP(NSTEPS),NSTEPS) = 0         ! position in NSTEP - 1 that gave this atom
30    CONTINUE
C
      if (phelp) then
         print *, ' after 30 continue'
         is = nkeep(nsteps)
         print 1211, num_atom, nsteps, is, (keep(1,ll,nsteps),
     x               ll=1,is)
         print 1212, (keep(2,ll,nsteps),ll=1,is)
         print 1213, (keep(3,ll,nsteps),ll=1,is)
      endif
C
40       MSTEPS = NSTEPS
         NSTEPS = NSTEPS + 1
         NKEEP(NSTEPS) = 0
         IS = 0
C
         DO 1000 J=1,NKEEP(MSTEPS)
            N = KEEP(1,J,MSTEPS)
            IF (N .LT. 0) GO TO 1000   ! -; don't use to propagate
            DO 100 I1=1,NCON2(N)
               I2 = ICON2(I1,N)
               IF (NCON2(I2) .EQ. 1) GO TO 100      ! can't be in ring
               IF (I2 .EQ. KEEP(2,J,MSTEPS)) GO TO 100
               IF (IS + 1 .GT. NENTRIES) THEN
                  WRITE (COMMENT, 64) NAME(NUM_ATOM), NUM_ATOM,
     x                                NAME(IB), IB, IS + 1
                  PRINT 64, NAME(NUM_ATOM), NUM_ATOM, NAME(IB), IB,
     x                      IS + 1
64                FORMAT (' $$$$ PROBLEMS: failure in BOND_CHECK for ',
     x                   A5,' (#',I2,') - ',A5,' (#',I2,')',
     x                   ', IS + 1 =',I4,'...stopping')
                  STOP
               ENDIF
               IS = IS + 1                ! got one to put away
               KEEP(1,IS,NSTEPS) = I2
               KEEP(2,IS,NSTEPS) = N    ! I2 came from N in NSTEPS - 1
               KEEP(3,IS,NSTEPS) = J    ! I2's position in NSTEPS -1
               NKEEP(NSTEPS) = IS       ! update NKEEP
100         CONTINUE
1000     CONTINUE
C
      if (phelp) then
         print 1204, is, nkeep(nsteps)
1204     format (' after 1000 continue: is, nkeep(nsteps) =',2i3,'...')
         print 1211, num_atom, nsteps, is, (keep(1,ll,nsteps),
     x               ll=1,is)
1211     format (' line done:  num_atom =',i3,', nsteps =',i3,
     x           ', is =',i3/
     x           '             keep(1...) =',2(30I3))
         print 1212, (keep(2,ll,nsteps),ll=1,is)
1212     format ('             keep(2...) =',2(30i3))
         print 1213, (keep(3,ll,nsteps),ll=1,is)
1213     format ('             keep(3...) =',2(30i3))
      endif
C
      IF (NSTEPS .EQ. 1) GO TO 40
      IF (NKEEP(NSTEPS) .LE. 0) GO TO 4000
C
C----Is an atom in KEEP the same as an atom in the generating path?
      DO 500 I2=1,NKEEP(NSTEPS)
         M = KEEP(1,I2,NSTEPS)
         IF (M .LE. 0) GO TO 500
         IPLACE = KEEP(3,I2,NSTEPS)
         DO 480 L=NSTEPS-1,1,-1
            I4 = KEEP(1,IPLACE,L)
            IF (M .EQ. I4) THEN
               KEEP(1,I2,NSTEPS) = -M
               GO TO 500
            ENDIF
            IPLACE = KEEP(3,IPLACE,L)
480      CONTINUE
500   CONTINUE
C
C----Is an atom in KEEP the same as NUM_ATOM?
      DO 700 I3=1,NKEEP(NSTEPS)
         M = KEEP(1,I3,NSTEPS)
         IF (M .LE. 0) GO TO 700    ! prevents going around a ring
         IF (M .NE. NUM_ATOM) GO TO 700
C----Look at the atom just before M.  The number of times this atom
C     is accessed should not be greater than NCON2 - 1.
         N = KEEP(2,I3,NSTEPS)
         NTRYS(N) = NTRYS(N) + 1
C      print 2133, nsteps, i3, m, n, ntrys(n), ncon2(n)
C2133  format ('##nsteps, i3, m, n, ntrys(n), ncon2(n) =',6i3)
         IF (NTRYS(N) .GE. NCON2(N)) THEN
            KEEP(1,I3,NSTEPS) = -M
         ENDIF
C----Have a possible new ring.  Check the atom path to see if any
C     non-adjacent atoms are connected.  If yes, this path should
C     be rejected.  Don't bother if NSTEPS .le. 3.  NTEMP will contain
C     the path atoms, starting with the last one (same as NUM_ATOM).
C--------Get path information
         IPLACE = KEEP(3,I3,NSTEPS)
         N_PATHS = N_PATHS + 1  ! increment # of possible ring paths
         ICT = 1            ! counter for atoms in the path
         NTEMP(N_PATHS,ICT) = M     ! record the last atom
         DO 580 L=NSTEPS-1,1,-1
            ICT = ICT + 1                    ! get atom path
            NTEMP(N_PATHS,ICT) = KEEP(1,IPLACE,L)    ! into NTEMP
            IPLACE = KEEP(3,IPLACE,L)
580      CONTINUE
      if (phelp) then
         print 2121, n_paths, ict, (ntemp(n_paths,i),i=1,ict)
2121     format (' path info: n_paths, ict, ntemp(n_paths,i) =',32i3)
      endif
         NICT(N_PATHS) = ICT
C--------If NSTEPS .le. 3, don't check non-adjacent atoms
         IF (NSTEPS .LE. 3) GO TO 680
         DO 590 I=1,ICT-2
            I1 = NTEMP(N_PATHS,I)       ! look at positions no closer than
            DO 595 J=I+2,ICT            ! 1 & 3 and don't do 1st and last
               J1 = NTEMP(N_PATHS,J)    ! because they are linked
               IF (I .EQ. 1 .AND. J .EQ. ICT) GO TO 590
               DO 598 K=1,NCON2(J1)
                  IF (ICON2(K,J1) .NE. I1) GO TO 598
       if (phelp) then
          print 2122, i1, j1
2122      format (' intermediate ring found between',i3,' &',i3)
       endif
                  KEEP(1,I3,NSTEPS) = -M   ! flag...recognized as ring former
                  NICT(N_PATHS) = -ICT    ! bad path, flag as negative
                  GO TO 700
598            CONTINUE
595         CONTINUE
590      CONTINUE
680         NRINGS = NRINGS + 1
            NSIZES(NRINGS) = NSTEPS
            KEEP(1,I3,NSTEPS) = -M   ! flag...recognized as ring former
      if (phelp) then
          print 2244, num_atom, m, nsteps, nrings, nsizes(nrings)
2244      format (' ring found: num_atom, m, nsteps, nrings, ',
     x            'nsizes(nrings) =',5i3)
      endif
700   CONTINUE
C
C----If there are identical atoms in KEEP(1...), which have the same
C     precusor atom [KEEP(2...)], eliminate the second one
701   DO 520 I2=1,NKEEP(NSTEPS)-1
        M1 = KEEP(1,I2,NSTEPS)
        IF (M1 .LE. 0) GO TO 520
        DO 510 I3=I2+1,NKEEP(NSTEPS)
           M2 = KEEP(1,I3,NSTEPS)
           IF (M2 .LE. 0) GO TO 510
           IF (M1 .NE. M2) GO TO 510
           IF (KEEP(2,I2,NSTEPS) .NE. KEEP(2,I3,NSTEPS)) GO TO 510
              KEEP(1,I3,NSTEPS) = -M2  ! match...eliminate 2nd one
510      CONTINUE
520   CONTINUE
C
3995     IF (NKEEP(NSTEPS) .EQ. 0) GO TO 4000
         IF (NRINGS .GE. N1 + MOD(N1,2)) GO TO 4000  ! N1 = NCON2(NUM_ATOM)
         IF (NSTEPS .LE. MAXSTEPS) GO TO 40
C
4000  CONTINUE
C
      if (phelp) then
         print 2265, n_paths
2265     format (' after 4000 continue, n_paths =',I3)
         do i=1,n_paths
            print 2277, i, nict(i), (ntemp(i,j),j=1,nict(i))
2277        format (' rings before modifications: i, nict(i), ',
     x              '(ntemp(i,j)j=1,nict(i)) =',30i3)
         enddo
      endif
C
C----Check the possible ring paths...remove any duplicates
      IF (N_PATHS .LE. 0) GO TO 5000
      NRINGS = 0                ! initialize...will count unique rings
      IF (N_PATHS .EQ. 1) THEN
         NRINGS = NRINGS + 1
         NSIZES(NRINGS) = NICT(1)
         N_IN_NTEMP(NRINGS) = 1          ! record line in NTEMP array
         GO TO 4101
      ENDIF
      DO 4100 I=1,N_PATHS-1
         IF (NICT(I) .LE. 0) GO TO 4100
         DO 4090 J=I+1,N_PATHS
            IF (NICT(J) .LE. 0) GO TO 4090
            IF (NICT(I) .NE. NICT(J)) GO TO 4090        ! ring sizes must be equal
            IF (NTEMP(I,1) .NE. NTEMP(J,1)) GO TO 4090  ! starting atoms must be the same
               DO 4080 K1=2,NICT(I)
                  K2 = NICT(I) + 2 - K1                 ! check the ring atom paths
                  IF (NTEMP(I,K1) .NE. NTEMP(J,K2)) GO TO 4090
4080           CONTINUE
C--------------Got a match...flag NICT(J) and save info
               NRINGS = NRINGS + 1
               NSIZES(NRINGS) = NICT(I)
               N_IN_NTEMP(NRINGS) = I    ! record line in NTEMP array
               NICT(J) = -NICT(J)
               GO TO 4100
4090     CONTINUE
C--------Dropping thru 4090 loop means only 1 path for this ring
         NRINGS = NRINGS + 1
         NSIZES(NRINGS) = NICT(I)
         N_IN_NTEMP(NRINGS) = I          ! record line in NTEMP array
4100  CONTINUE
4101  CONTINUE
      if (phelp) then
         print 2299, nrings, (nsizes(i),i=1,nrings)
2299     format (' unique rings: nrings, nsizes =',5i3)
      endif
C
C----Start of section to look for aromatic 5 & 6-rings
C
      IF (NRINGS .EQ. 0) GO TO 5000
      DO 4900 I=1,NRINGS
         J = N_IN_NTEMP(I)
         ISIZE = NICT(J)
C----Can it be a 5-atom aromatic ring...ISIZE = 5?
         IF (ISIZE .NE. 5) GO TO 4200   ! check for aromatic
         NSP2 = 0   ! must have 5 sp2 atoms in a ring for aromatic
         IF (NAME(NUM_ATOM)(1:1) .EQ. 'C' .AND.
     X       NCON(NUM_ATOM) .EQ. 3) THEN       ! if C, must be 3-linked
             IF (IQ_CARBONYL2(NUM_ATOM) .EQ. 0) GO TO 4010 ! reject if C=O
         ENDIF
         IF (NAME(NUM_ATOM)(1:1) .EQ. 'N' .AND. NCON(NUM_ATOM) .EQ. 2)
     X       GO TO 4010                   ! if N, must be 2 or 3-linked
         IF (NAME(NUM_ATOM)(1:1) .EQ. 'N' .AND. NCON(NUM_ATOM) .EQ. 3
     X       .AND. ANGLES_N(NUM_ATOM) .GE. 350.0) GO TO 4010
         IF (NAME(NUM_ATOM)(1:1) .EQ. 'O' .AND.
     X       NCON(NUM_ATOM) .EQ. 2) GO TO 4010 ! if O, must be 2-linked
         IF (NAME(NUM_ATOM)(1:1) .EQ. 'S' .AND.
     X       NCON(NUM_ATOM) .EQ. 2) GO TO 4010 ! if S, must be 2-linked   ! 3/22/00
         GO TO 4900
4010     NSP2 = NSP2 + 1
C----For 5-ring, count number of O & N and S in the ring         ! 3/22/00
      N_S_5RING = 0                                              ! 3/22/00
      N_O_5RING = 0
      IF (NAME(NUM_ATOM)(1:1) .EQ. 'O' .OR.
     X    NAME(NUM_ATOM)(1:1) .EQ. 'N') N_O_5RING = 1
      IF (NAME(NUM_ATOM)(1:1) .EQ. 'S') N_S_5RING = 1
      DO 4901 I1=2,ISIZE     ! first atom (NUM_ATOM) already processed
C--------Get ring atom info from Jth line of NTEMP
         M = NTEMP(J,I1)
         IF (NAME(M)(1:1) .EQ. 'O' .OR.
     X       NAME(M)(1:1) .EQ. 'N') N_O_5RING = N_O_5RING + 1
         IF (NAME(M)(1:1) .EQ. 'S') N_S_5RING = N_S_5RING + 1    ! 3/22/00
C--------Look for sp2 C and N.  If N, can be 2 or 3-linked.
         IF (NAME(M)(1:1) .EQ. 'C' .AND. NCON(M) .EQ. 3) THEN ! if C, must be 3-linked
            IF (IQ_CARBONYL2(M) .EQ. 0) GO TO 4890 ! reject if C=O
         ENDIF
         IF (NAME(M)(1:1) .EQ. 'N' .AND. NCON(M) .EQ. 2) GO TO 4890   ! if N, must be 2 or 3-linked
         IF (NAME(M)(1:1) .EQ. 'N' .AND. NCON(M) .EQ. 3 .AND.   ! if 3-linked N, angle sum
     x         ANGLES_N(M) .GE. 350.0)  GO TO 4890              ! must be .ge. 350.0 deg
         IF (NAME(M)(1:1) .EQ. 'O' .AND. NCON(M) .EQ. 2) GO TO 4890 ! if O, must be 2-linked
         IF (NAME(M)(1:1) .EQ. 'S' .AND. NCON(M) .EQ. 2) GO TO 4890 ! S must be 2-linked  ! 3/22/00
         GO TO 4900                 ! not acceptable
4890        NSP2 = NSP2 + 1
4901     CONTINUE
         IF (NSP2 .EQ. 5) THEN    ! 5 sp2 atoms in a ring...now check C-C bond lengths
            NSIZES(I) = -NSIZES(I)   ! 5 sp2 atoms (C, N, O, S) in a ring     ! 3/22/00
C-----------Check C-C distances in the ring path.  If any > 1.445 Angs,
C            code ring as non-aromatic by making NSIZES +
            DO 4970 I1=1,ISIZE
               M = NTEMP(J,I1)
               IF (NAME(M)(1:1) .NE. 'C') GO TO 4970
               L = NTEMP(J,I1+1)
               IF (I1 + 1 .GT. ISIZE) L = NTEMP(J,1)
               IF (NAME(L)(1:1) .NE. 'C') GO TO 4970
               D = DISTANCE(M,L)
C------------If C-S or N-S, make DRMAX = 1.79                             ! 3/22/00
              IF (NAME(M)(1:1) .EQ. 'S' .OR. NAME(L)(1:1) .EQ. 'S') THEN  ! 3/22/00
                 DRMAX = 1.79                                             ! 3/22/00
                 GO TO 4930                                               ! 3/22/00
              ENDIF                                                       ! 3/22/00
C      print 3388, num_atom, l, m, d
C3388  format (' num_atom, l, m, distance =',3i3,f10.3)
               IF (N_O_5RING .GE. 1 .OR. N_S_5RING .GE. 1) THEN  ! if 5-ring contains O or N  ! 3/22/00
                  DRMAX = 1.48               ! increase permissable d max
               ELSE                          ! to 1.48 Angs
                  DRMAX = 1.445
               ENDIF
4930           IF (D .GT. DRMAX) THEN                                      ! 3/22/00
                  NSIZES(I) = IABS(NSIZES(I))  ! make it +
                  GO TO 4975
               ENDIF
4970        CONTINUE
        ENDIF
4975    CONTINUE
      if (phelp) then
         print 1280, nrings, (nsizes(ll),ll=1,nrings)
1280     format (' 5-ring found: nrings, nsizes =',6i3)
      endif
      GO TO 4900
C
C----Can it be a 6-atom aromatic ring...ISIZE = 6?  C's must be 3-linked
C     and N's must be 2 or 3-linked
4200  IF (ISIZE .NE. 6) GO TO 4900
      NSP2 = 0   ! must have 6 sp2 atoms in a ring for aromatic
      IF (NCON(NUM_ATOM) .EQ. 3 .AND. NAME(NUM_ATOM)(1:1) .EQ. 'C')
     x       GO TO 4020
      IF (NCON(NUM_ATOM) .EQ. 3 .AND. NAME(NUM_ATOM)(1:1) .EQ.'N'
     x       .AND. ANGLES_N(NUM_ATOM) .GE. 350.0) GO TO 4020
      IF (NCON(NUM_ATOM) .EQ. 2 .AND. NAME(NUM_ATOM)(1:1) .EQ. 'N')
     x       GO TO 4020
      GO TO 4900       ! not acceptable
4020  NSP2 = NSP2 + 1
      DO 4910 I1=2,ISIZE
         M = NTEMP(J,I1)
C--------Look for sp2 C and N.  If N, can be 2 or 3-linked.
         IF (NCON(M) .EQ. 3 .AND. NAME(M)(1:1) .EQ. 'C') THEN
            IF (IQ_CARBONYL2(M) .EQ. 0) GO TO 3890
         ENDIF
         IF (NCON(M) .EQ. 3 .AND. NAME(M)(1:1) .EQ. 'N' .AND.
     x           ANGLES_N(M) .GE. 350.0) GO TO 3890
         IF (NCON(M) .EQ. 2 .AND. NAME(M)(1:1) .EQ. 'N') GO TO 3890
         GO TO 4900           ! not acceptable
3890     NSP2 = NSP2 + 1
4910     CONTINUE
C----C-C bond length check section
         IF (NSP2 .EQ. 6) THEN    ! 6 sp2 atoms in a ring...now check C-C bond lengths
            NSIZES(I) = -NSIZES(I)   ! 6 sp2 atoms (C, N, O) in a ring
C-----------Check C-C distances in the ring path.  If any > 1.44 Angs,
C            code ring as non-aromatic by making NSIZES +
            DO 5970 I1=1,ISIZE
               M = NTEMP(J,I1)
               IF (NAME(M)(1:1) .NE. 'C') GO TO 5970
               L = NTEMP(J,I1+1)
               IF (I1 + 1 .GT. ISIZE) L = NTEMP(J,1)
               IF (NAME(L)(1:1) .NE. 'C') GO TO 5970
               D = DISTANCE(M,L)
C       print 3388, num_atom, l, m, d
               IF (D .GT. 1.44) THEN
                  NSIZES(I) = IABS(NSIZES(I))  ! make it +
                  GO TO 5975
               ENDIF
5970        CONTINUE
        ENDIF
5975    CONTINUE
      if (phelp) then
         print 1238, nrings, (nsizes(ll),ll=1,nrings)
1238     format (' 6-ring found: nrings, nsizes =',6i3)
      endif
4900  CONTINUE
C----Finished
5000  IF (NRINGS .EQ. 0) GO TO 5100
C----Check the contents of the ring paths to determine how many of
C     the atoms linked to NUM_ATOM are in rings
      DO I=1,4            ! put a 1 into ILAST for each occurence
         ILAST(I) = 0     !  of a linked atom in the NTEMP array
      ENDDO
      DO 5050 I=1,N1  ! look at NCON2(NUM_ATOM) atoms
         K = 0
         DO 5040 J=1,NRINGS
5035        K = K + 1
            IF (NICT(K) .LE. 0) GO TO 5035
            DO 5025 L=1,NICT(K)
               IF (ICON(I,NUM_ATOM) .EQ. NTEMP(K,L)) ILAST(I) = 1
5025        CONTINUE
5040     CONTINUE
5050  CONTINUE
      NBONDS = 0
      DO J=1,4
         IF (ILAST(J) .GE. 1) NBONDS = NBONDS + 1
      ENDDO
5100  IF (NRINGS .EQ. 0) NBONDS = 0
      if (phelp) then
         print 2291, num_atom, nbonds, nrings, (nsizes(i),i=1,nrings)
2291     format (' final summary: num_atom, nbonds in rings, nrings,',
     x           ' nsizes =',
     x           7i3)
      endif
      RETURN
      END
C------------------------------------------------------------
C    Count the number of ring-sizes that are negative
      FUNCTION NEG_COUNT (NSIZES)
C
      DIMENSION NSIZES(10)
C
      NEG_COUNT = 0
C
      DO 100 I=1,10
         IF (NSIZES(I) .LT. 0) NEG_COUNT = NEG_COUNT + 1
100   CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------
      SUBROUTINE WRFILE(REFCOD_2)
C
      COMMON /PROBLEMS/ JPROB, COMMENT, SPACER, ICSP3, JCSP3
      CHARACTER REFCOD_2*8, COMMENT*120, SPACER*8
      IF (REFCOD_2(1:6) .EQ. '      ') RETURN
      WRITE (15, '(A8,3X,A120)') REFCOD_2, COMMENT
      RETURN
      END
C
C-------------------------------------------------------------
C----Calculate bond angle     1-2-3
      FUNCTION ANGLE(XYZ1, XYZ2, XYZ3)
C
      DIMENSION XYZ1(3), XYZ2(3), XYZ3(3), XL(2), VEC(3,2)
C
      DO 38 K=1,3
         VEC(K,1) = XYZ1(K) - XYZ2(K)
         VEC(K,2) = XYZ3(K) - XYZ2(K)
38    CONTINUE
C
      DO 45 L=1,2
         XL(L) = 0.0
         DO 44 K=1,3
            XL(L) = XL(L) + VEC(K,L)**2
44       CONTINUE
         XL(L) = SQRT(XL(L))
45    CONTINUE
C
      ANGLE = 0.0
      DO 50 K=1,3
         ANGLE = ANGLE + VEC(K,1)*VEC(K,2)
50    CONTINUE
      ANGLE = 57.296*ACOS(ANGLE/(XL(1)*XL(2)))
C
      RETURN
      END
C
C--------------------------------------------------------------
      FUNCTION ANGLES_N(NUM_ATOM)    ! calc sum of angles at a
C                                      3-linked N
      CHARACTER NAME*5
      COMMON /COM1/ NCON(500), ICON(5,500), NAME(500),
     x                 OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
C
      ANGLES_N = 0.0
      DO 100 I=1,NCON(NUM_ATOM)-1
         I1 = ICON(I,NUM_ATOM)
         DO 100 J=I+1,NCON(NUM_ATOM)
            I2 = ICON(J,NUM_ATOM)
            ANGLES_N = ANGLES_N + ANGLE(OR_XYZ(1,I1),
     x                 OR_XYZ(1,NUM_ATOM), OR_XYZ(1,I2))
100   CONTINUE
      RETURN
      END
C------------------------------------------------------------------
