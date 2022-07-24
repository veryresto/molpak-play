C----ATOM_CODE-VOL-ADDITIVITY.F                       7/15/05
C
C----Program to use the molpak.xyz coordinates and table of atom codes, volumes
C     and sigmas to calculate a density.  Codes have been extended by including
C     multiple fragment structures in the initial database, halides anions
C     (no F-), ClO4- and 1-6 valent C + 3-6 valent B (to handle carboranes).  2/28/05
C
C----Build a code for each atom in a molecule that contains the following
C     information...
C          position     function
C             1         element type: H, C, N, O, F, Z (CL), X (BR), S, P, I, B 2/22/05
C             2         # 1-linked C's
C             3           2
C             4           3
C             5           4
C             6           5
C             7           6
C             8         # 1-linked N's
C             9           2
C            10           3
C            11           4
C            12         # 1-linked O's
C            13         # 2
C            14         # 1-linked S's
C            15           2
C            16           3
C            17           4
C            18           5
C            19           6
C            20         # 1-linked P's
C            21           2
C            22           3
C            23           4
C            24           5
C            25           6
C            26         # 3-linked B's (only accept 3-6
C            27           4             linked B's)
C            28           5
C            29           6
C            30         # H's
C            31         # F's
C            32         # Cl's (Z)
C            33         # Br's (X)
C            34         # I's (only accept 1-linked I's)
C            35         # Cl-'s (Z)
C            36         # Br-'s (X)
C            37         # I-'s
C            38         # Cl in ClO4
C
      PARAMETER (NPOTS=3000, NTYPES_ATOMS=11)                         ! 2/28/05
      CHARACTER (LEN=120) :: COMMENT
      CHARACTER (LEN=8) :: REFCOD_2
      CHARACTER ALLTYPES(NPOTS)*38, ATYPE*38, BLANK*38                ! 2/28/05
      character interpretation(npots)*28                              ! 3/8/05
      COMMON /TYPES/ ALLTYPES, n_types, blank                         ! 2/28/05
      DIMENSION KINDS(NPOTS), AT_VOL(NPOTS), AT_SIGMA(NPOTS),
     X          NUM_VOLS(NPOTS)
      CHARACTER ATSYM*1, ATSYM2*2, MWOUT*36, DIR*120
      CHARACTER*7 DATNM2
      CHARACTER*48 INFILE
      CHARACTER*5 LABEL(500), NEW_LABEL
      CHARACTER*2 LABEL_OUT                                            ! 2/26/00
      EQUIVALENCE (OUTNAM, OUTFIL), (DATNAM, DATNM2)
      DIMENSION N_AT_TYPE(NPOTS)
      DIMENSION ICELL(6), NVALUE(37), XYZ(3,500), DEN(2),           ! 2/26/00
     x          AT_WTS(NTYPES_ATOMS)                                 ! 2/26/00
      DIMENSION O_XYZ(3), IND(NPOTS),
     x          NATTYPES(NPOTS), NSIZES(4)
      COMMON /COM1/ NCON(500), ICON(7,500), NEW_LABEL(500),          ! 5/8/00
     x                 OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x                 I3D(500), KEEP_NPOT(500)
      COMMON /CONTENTS/ NCON_FDAT(NTYPES_ATOMS),                  ! 2/26/00
     x                  NCON_FBIB(NTYPES_ATOMS),                  ! 2/26/00
     x                  MAX_CONN(NTYPES_ATOMS),                   ! 5/8/00
     x                  ATSYM(NTYPES_ATOMS),                      ! 2/28/05
     x                  n_zxi, n_clo4, n_Z, n_X, n_I              ! 2/28/05
      COMMON /ELEMENTS/ N_AT_DIFF(11), ATSYM2(NTYPES_ATOMS)       ! 2/28/05
      COMMON /PROBLEMS/ JPROB, COMMENT, REFCOD_2, ICSP3, JCSP3
      DATA DEN /10000., 100000./
      DATA AT_WTS /12.011,1.0079,18.9984,14.0067,15.9994,79.904,     ! 3/4/00
     x             35.4527,32.066,30.9738,126.9045,10.811/          ! 2/28/05
C                     C     H       F       N       O      BR       ! 2/28/05
C                     CL    S       P       I       B               ! 2/28/05
C----Br will be recoded as atom ID "X" and CL as "Z" to use only
C     1 character
C
C----File # 10 = MOLPAK.NAME file
      OPEN (UNIT=10,FILE='MOLPAK.NAME',STATUS='OLD')
      READ (10,'(A)') INFILE             ! name of xyz file
      CLOSE (UNIT=10)
C
C----Open xyz file
      OPEN (UNIT=9,FILE=INFILE,STATUS='OLD')
C
C----File # 24 = atom code and volume data
      DIR='/PREDICTIONS/UTILITIES/atom_code_volume.data'
!     DIR='/export/software/predictions-100cyc04/UTILITIES/atom_code_vol&
!    &ume.data'
         OPEN (UNIT=24, FILE=DIR, STATUS='OLD')
         READ (24,*)              ! skip first, title line
         N_TYPES = 0
753      N_TYPES = N_TYPES + 1
         READ (24,754,END=750) I, AT_VOL(N_TYPES), AT_SIGMA(N_TYPES),
     2                       NUM_VOLS(N_TYPES), ALLTYPES(N_TYPES),      ! 3/8/05
     3                       interpretation(n_types)                    ! 3/8/05
         write (26,755) i, at_vol(n_types), at_sigma(n_types),
     2                  num_vols(n_types), alltypes(n_types),           ! 3/8/05
     3                  interpretation(n_types)                         ! 3/8/05
755      FORMAT (I5,2F8.3,I8,' /',A38,'/  ',a)                          ! 3/8/05
754      FORMAT (I5,1X,2F9.3,I8,6X,A38,1x,a)                            ! 3/8/05
         IF (I .EQ. N_TYPES) GO TO 753                                     ! 5/5/00
            WRITE (6,741) I, N_TYPES, ALLTYPES(N_TYPES)                    ! 5/5/00
741         FORMAT ('Problem reading atom_codes file, file line # =',
     2              I3,', file # =',I3,';'/
     3              'atom type code = /',A38,'/')                ! 2/28/05
            STOP                                                 ! 5/5/00
750      N_TYPES = N_TYPES - 1                                   ! 5/1/00
         CLOSE (UNIT=24)                                         ! 5/1/00
         WRITE (6,745) N_TYPES
745      FORMAT (' Predicted density from atom_code volume ',
     2           'additivity data base...'/
     3           I5,' atom_codes read and copied to fort.26')
C
C----Read atom id's and xyz's from molpak.xyz file
        NATOMS = 0
118     NATOMS = NATOMS + 1
        READ (9,125,END=119) LABEL(NATOMS), (XYZ(J,NATOMS),J=1,3)
125     FORMAT (5X,A5,1X,3F10.6)
        GO TO 118
119     NATOMS = NATOMS - 1
C----Search atoms labels for BR & replace as X and CL & replace as Z
        DO 124 I=1,NATOMS                                       ! 2/26/00
           IF (LABEL(I)(1:2) .NE. 'BR') GO TO 655               ! 2/28/00
               LABEL(I)(1:2) = 'X '                             ! 2/28/00
               GO TO 124                                        ! 2/28/00
655        IF (LABEL(I)(1:2) .NE. 'CL') GO TO 124               ! 3/4/00
               LABEL(I)(1:2) = 'Z '
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
184     NA = NA + 1
        NEW_LABEL(NA) = LABEL(I)
        DO L = 1,3
           OR_XYZ(L,NA) = XYZ(L,I)
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
C----Calc molecular weight
701   WT = 0.0
      DO 175 I = 1,NTYPES_ATOMS
         WT = WT + NCON_FDAT(I)*AT_WTS(I)
175   CONTINUE
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
238   PRINT 195, MWOUT(1:LAST), WT, FRACT_H
195   FORMAT (' Molecular formula: ',A/
     x        ' Molecular mass:',F8.3,' g/mol',
     x        '; H fraction:',F6.3)
C
C----Check for more than 6 connections for C, 4 for N, 2 for O, ! 3/8/05
C     1 for F, 6 for S, 6 for P, etc                            ! 3/8/05
      CALL CONNECTION(ITEST,ICOORD)  ! check for connections
      IF (ITEST .EQ. 0) GO TO 282
         N = NCON(ITEST)
         WRITE (6,283) N, ITEST, NEW_LABEL(ITEST)
283      FORMAT (' ---> PROBLEMS:',I3,' connections found for atom',
     x           I3,' (',A5,')...quit')
         STOP
C
C----Loop through the atoms to identify types
282   V_MOLEC = 0.0             ! molecular volume
      DO 295 L=1,NA
         NUM_ATOM = L                                          ! 4/20/00
         J = NPOT(NUM_ATOM, ATYPE)                             ! 4/20/00
         IF (J .NE. 0) THEN           ! J .eq. 0 = success     ! 4/24/00
            WRITE (6,285) NUM_ATOM                             ! 4/20/00
285         FORMAT (' $$$$ PROBLEMS with determining code for',
     x              ' atom #',I3,'...quit')
            STOP
         ENDIF                                                 ! 4/20/00
         DO 289 M=1,N_TYPES
            IF (ATYPE .NE. ALLTYPES(M)) GO TO 289
               KINDS(L) = M
C----Check for - NUM_VOLS
               IF (NUM_VOLS(M) .LE. 0) THEN
                  WRITE (6,444) ALLTYPES(M)
444               FORMAT ('--> No volume data for atom_code /',
     x                    A38,'/  ....quitting')               ! 2/28/05
                  STOP
                  ENDIF
               V_MOLEC = V_MOLEC + AT_VOL(M)   ! add up atom volumes
               GO TO 295
289      CONTINUE
         WRITE (6,1133) L, NEW_LABEL(L), ATYPE
1133     format ('Cannot locate tabulated atom code for atom #',I3,
     x            ' (',A5,')'/' code = /',A38,'/',                  ! 3/8/05
     x            '  ....quitting')
         STOP
295   CONTINUE                                                ! 4/20/00
C
C----Record number of atoms of each type
      DO L=1,NPOTS                                             ! 4/20/00
          N_AT_TYPE(L) = 0                                     ! 4/20/00
      ENDDO                                                    ! 4/20/00
      DO 380 L=1,NA            ! count kinds of atoms          ! 4/20/00
          K = KINDS(L)                                         ! 4/20/00
          N_AT_TYPE(K) = N_AT_TYPE(K) + 1  ! local molecule total
380   CONTINUE                                                 ! 4/20/00
C
C----Show which atom types have been determined
408   I1 = 0
      DO 415 I2=1,NPOTS
         IF (N_AT_TYPE(I2) .EQ. 0) GO TO 415
         I1 = I1 + 1
         IND(I1) = I2
         NATTYPES(I1) = N_AT_TYPE(I2)
415   CONTINUE
C
      vol_sigma = 0.0                                          ! 3/8/05
C----Summarize molecule atom types and sum info for vol sigma  ! 3/8/05
        WRITE (6,418)                                          ! 3/8/05
418     FORMAT ('  # atoms, atom type, volume, vsigma,'        ! 3/8/05
     2          ' atom_code, interpretation...')               ! 4/8/05
      DO I2=1,I1
        WRITE (6,419) NATTYPES(I2), IND(I2),                         ! 3/8/05
     2       AT_VOL(IND(I2)), at_sigma(ind(i2)), alltypes(ind(i2)),  ! 3/8/05
     3       interpretation(ind(i2))                                 ! 3/8/05
419     format (5x,I3,',',I4,',',f7.3,',',f7.3,', /',A38,'/',1x,a,'/') ! 4/8/05
      vol_sigma = vol_sigma + (nattypes(I2)*at_sigma(i2))**2   ! 3/8/05
      ENDDO
C
      DENSITY = WT/(.6022*V_MOLEC)
      density_sigma = (density/v_molec)*sqrt(vol_sigma)        ! 3/8/05
      WRITE (6,4121) V_MOLEC, DENSITY, density_sigma           ! 3/8/05
4121  FORMAT (' Calcd molecular volume =',F9.3,' Angs**3, ',
     2        ' calcd density =',F6.3,' g/cc,',                ! 3/8/05
     3        ' sigma (density) =',f6.3,' g/cc')               ! 3/8/05
C----Finished...phew!
      STOP
      END
C
C------------SUBROUTINE CONNECTION-----------------------
C
      SUBROUTINE CONNECTION(ITEST,ICOORD)                       ! 7/15/05
C
C----Determine connectivity
C
      PARAMETER (NTYPES_ATOMS=11)                              ! 2/28/05
      CHARACTER LABEL*5, NAME(500)*5
      CHARACTER ATSYM*1
      COMMON /COM1/ NCON(500), ICON(7,500), LABEL(500),
     x              OR_XYZ(3,500), CELL(6), NATOMS,
     x              IT_TABLE(500), I3D(500), KEEP_NPOT(500)
      EQUIVALENCE (LABEL, NAME)
      DIMENSION D_H(2), N_H(2), dlist(8)
      COMMON /CONTENTS/ NCON_FDAT(NTYPES_ATOMS),               ! 2/26/00
     x                  NCON_FBIB(NTYPES_ATOMS),               ! 2/26/00
     x                  MAX_CONN(NTYPES_ATOMS),                ! 5/8/00
     x                  ATSYM(NTYPES_ATOMS),                   ! 2/28/05
     x                  n_zxi, n_clo4, n_Z, n_X, n_I           ! 2/18/05
C
      ITEST = 0
      n_good = 0                                               ! 2/28/05
      n_gZ = 0                                                 ! 2/28/05
      n_gX = 0                                                 ! 2/28/05
      n_gI = 0                                                 ! 2/28/05
C
C----Determine connectivity....use max distance of DMAX for most
      DO 502 I=1,NATOMS
         NCON(I) = 0
         DO 490 J=1,NATOMS
            DMAX = 1.70                                           ! 4/28/00
            IF (I .EQ. J) GO TO 490
            IF (LABEL(I)(1:1) .EQ. 'H' .AND.
     x          LABEL(J)(1:1) .EQ. 'H') GO TO 490
            D = DISTANCE (I, J)
C
            IF (LABEL(I)(1:1) .EQ. 'B' .OR.          ! B          ! 2/22/05
     X          LABEL(J)(1:1) .EQ. 'B') DMAX = 1.90               ! 2/22/05
C
            IF (LABEL(I)(1:1) .EQ. 'S' .OR.          ! S          ! 3/13/00
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 1.95               ! 3/28/00
C
            IF (LABEL(I)(1:1) .EQ. 'P' .OR.          ! P          ! 4/5/00
     X          LABEL(J)(1:1) .EQ. 'P') DMAX = 1.90               ! 4/5/00
C
            IF (LABEL(I)(1:1) .EQ. 'X' .OR.          ! X = BR     ! 2/26/00
     X          LABEL(J)(1:1) .EQ. 'X') DMAX = 2.16  ! X = BR     ! 5/3/00
C
            IF (LABEL(I)(1:1) .EQ. 'X' .AND.         ! X = BR     ! 5/2/00
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.25  ! S          ! 5/2/00
            IF (LABEL(J)(1:1) .EQ. 'X' .AND.         ! X = BR     ! 5/2/00
     X          LABEL(I)(1:1) .EQ. 'S') DMAX = 2.25  ! S          ! 5/2/00
C
            IF (LABEL(I)(1:1) .EQ. 'Z' .OR.          ! Z = CL     ! 3/4/00
     X          LABEL(J)(1:1) .EQ. 'Z') DMAX = 2.10  ! Z = CL     ! 2/4/05
C
            IF (LABEL(I)(1:1) .EQ. 'I' .OR.          ! I          ! 4/3/04
     X          LABEL(J)(1:1) .EQ. 'I') DMAX = 2.25  ! I          ! 4/3/04
C
            IF (LABEL(I)(1:1) .EQ. 'S' .AND.         ! S-S        ! 3/28/00
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.15               ! 3/28/00
C
            IF (LABEL(I)(1:1) .EQ. 'S' .AND.         ! S-P        ! 1/17/05
     X          LABEL(J)(1:1) .EQ. 'P') DMAX = 2.20               ! 1/17/05
            IF (LABEL(I)(1:1) .EQ. 'P' .AND.         ! S-P        ! 1/17/05
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.20               ! 1/17/05
C
            IF (LABEL(I)(1:1) .EQ. 'Z' .AND.         ! P-Z (Cl)   ! 2/12/05
     X          LABEL(J)(1:1) .EQ. 'P') DMAX = 2.20               ! 2/12/05
            IF (LABEL(I)(1:1) .EQ. 'P' .AND.         ! P-Z (Cl)   ! 2/12/05
     X          LABEL(J)(1:1) .EQ. 'Z') DMAX = 2.20               ! 2/12/05
C
            IF (LABEL(I)(1:1) .EQ. 'X' .AND.         ! X = BR     ! 5/2/00
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.25  ! S          ! 5/2/00
            IF (LABEL(J)(1:1) .EQ. 'X' .AND.         ! X = BR     ! 5/2/00
     X          LABEL(I)(1:1) .EQ. 'S') DMAX = 2.25  ! S          ! 5/2/00
C
            IF (LABEL(I)(1:1) .EQ. 'Z' .AND.         ! Z = Cl     ! 2/14/05
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.40  ! S          ! 2/14/05
            IF (LABEL(J)(1:1) .EQ. 'Z' .AND.         ! Z = Cl     ! 2/14/05
     X          LABEL(I)(1:1) .EQ. 'S') DMAX = 2.40  ! S          ! 2/14/05
C
            IF (LABEL(I)(1:1) .EQ. 'B' .AND.         ! B-S = 2.1  ! 2/23/05
     X          LABEL(J)(1:1) .EQ. 'S') DMAX = 2.10  !            ! 2/14/05
            IF (LABEL(J)(1:1) .EQ. 'B' .AND.         !            ! 2/23/05
     X          LABEL(I)(1:1) .EQ. 'S') DMAX = 2.10  !
C
            IF (LABEL(I)(1:1) .EQ. 'C' .AND.         ! C-C = 1.78 ! 7/15/05
     X          LABEL(J)(1:1) .EQ. 'C') DMAX = 1.78  !            ! 7/17/05
C
            IF (LABEL(I)(1:1) .EQ. 'H' .OR.
     X          LABEL(J)(1:1) .EQ. 'H') DMAX = 1.30
C
            IF (LABEL(I)(1:1) .EQ. 'B' .AND.        ! B-H = 1.4; ! 2/24/05
     X          LABEL(J)(1:1) .EQ. 'H') DMAX = 1.40 ! 1.4 to try ! 2/14/05
            IF (LABEL(J)(1:1) .EQ. 'B' .AND.        ! to eliminate ! 2/24/05
     X          LABEL(I)(1:1) .EQ. 'H') DMAX = 1.40 ! H..B..H bridges ! 2/24/04
C
            IF (LABEL(I)(1:1) .EQ. 'I' .AND.         ! eliminate I...I
     X          LABEL(J)(1:1) .EQ. 'I') DMAX = 3.5   ! anions, 2/8/05
            IF (LABEL(I)(1:1) .EQ. 'X' .AND.         ! eliminate Br...Br
     X          LABEL(J)(1:1) .EQ. 'X') DMAX = 3.1   ! anions, 2/9/05
            IF (LABEL(I)(1:1) .EQ. 'Z' .AND.         ! eliminate Cl...Cl
     X          LABEL(J)(1:1) .EQ. 'Z') DMAX = 2.5   ! anions, 2/8/05
            IF ((LABEL(I)(1:1) .EQ. 'I' .AND.
     X           LABEL(J)(1:1) .EQ. 'X') .OR.        ! eliminate I...Br
     X          (LABEL(I)(1:1) .EQ. 'X' .AND.
     X           LABEL(J)(1:1) .EQ. 'I')) DMAX = 3.3 ! anions, 2/8/05
            IF ((LABEL(I)(1:1) .EQ. 'Z' .AND.
     X           LABEL(J)(1:1) .EQ. 'X') .OR.        ! eliminate Cl...Br
     X          (LABEL(I)(1:1) .EQ. 'X' .AND.
     X           LABEL(J)(1:1) .EQ. 'Z')) DMAX = 2.7 ! anions, 2/8/05
            IF (D .GT. DMAX) GO TO 490
            NCON(I) = NCON(I) + 1
            ICON(NCON(I),I) = J
            IF (D .GT. DMAX) GO TO 490
            DO 350 L=1,NTYPES_ATOMS
               IF(LABEL(I)(1:1) .NE. ATSYM(L)) GO TO 350
               IF (NCON(I) .LE. MAX_CONN(L)) GO TO 490
C              if (ncon(i) .le. 4 .and. label(i)(1:1) .eq. 'Z')    ! 2/11/05
C    2             go to 490                                       ! 2/11/05
               ITEST = I
               GO TO 700
350         CONTINUE
490      CONTINUE
C----For boron, eliminate any with ncon .lt. 3, ie 1 or 2         ! 2/22/05
         if (label(i)(1:1) .eq. 'B' .and.                         ! 2/22/05
     2       ncon(i) .lt. 3) then                                 ! 2/22/05
            itest = i                                             ! 2/22/05
            return                                                ! 2/22/05
         endif                                                    ! 2/22/05
C----For halogen (Z, X, I), eliminated any ncon = 1 with another halogen   ! 2/16/05
         if (label(i)(1:1) .eq. 'X' .or.                          ! 2/16/05
     x       label(i)(1:1) .eq. 'Z' .or.                          ! 2/16/05
     x       label(i)(1:1) .eq. 'I') then                         ! 2/16/05
             if (ncon(i) .eq. 1) then                             ! 2/16/05
                 j = icon(1,i)                                    ! 2/16/05
                 if (label(j)(1:1) .eq. 'X' .or.                  ! 2/16/05
     2               label(j)(1:1) .eq. 'Z' .or.                  ! 2/16/05
     3               label(j)(1:1) .eq. 'I') then                 ! 2/26/05
                    itest = i                                     ! 2/16/05
                    return                                        ! 2/16/05
                 endif                                            ! 2/16/05
             endif
             if (ncon(i) .eq. 0) then         ! assume a halide  ! 2/16/05
                  print 1122, i, label(i), ncon(i)
1122              format ('>>>>Halide, atom #',i3,', label = |',
     2                    a,'|, ncon =',i2)
             endif
         endif                                                    ! 2/16/05
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
      go to 715                                                 ! 2/28/05
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
715   DO 721 I=1,NATOMS
C----Only allow ncon = 0, 1, 4 for Cl, and ncon = 0 or 1 for Br and I
C----Eliminate any Cl's with ncon = 2 or 3                          ! 2/14/05
       if (label(i)(1:1) .eq. 'Z' .and. (ncon(i) .eq. 2 .or.        ! 2/14/05
     2     ncon(i) .eq. 3)) go to 788                               ! 2/14/05
C----Eliminate any Br's with ncon .gt.1                             ! 2/14/05
       if (label(i)(1:1) .eq. 'X' .and. ncon(i) .gt. 1) go to 788   ! 2/14/05
C----Eliminate any I's with ncon .gt.1                              ! 2/14/05
       if (label(i)(1:1) .eq. 'I' .and. ncon(i) .gt. 1) go to 788   ! 2/14/05
C----Perchlorate, CLO4-
       if (ncon(i) .eq. 4 .and. label(i)(1:1) .eq. 'Z') then ! 2/11/05
          n_ox = 0                    ! # 1-linked O's on CL     2/14/05
          do l=1,4                                             ! 2/14/05
             j = icon(l,i)                                     ! 2/14/05
             if (ncon(j) .eq. 1 .and. label(j)(1:1) .eq. 'O')  ! 2/14/05
     2              n_ox = n_ox + 1                            ! 2/14/05
          enddo
          if (n_ox .ne. 4) go to 788  ! 4 1-linked O's?        ! 2/14/05
          n_clo4 = n_clo4 + 1                               ! 2/11/05
          write (46,748) label(i)(1:1), n_clo4              ! 2/28/05
748          format (10x,'CLO4- is ',a,', #',i3)            ! 2/28/05
          print 749, i, label(i), ncon(i)                   ! 2/22/05
749          format ('>>>>CLO4-, Cl atom #',i3,            ! 2/11/05
     2                ', label = |',a,'|, ncon =',i2)       ! 2/11/05
          go to 721                                         ! 2/11/05
       endif
C----Non-halogens
        if ((label(i)(1:1) .ne. 'X' .or.     ! eliminate ncon = 0   ! 2/4/05
     x         label(i)(1:1) .ne. 'Z' .or.     ! non-halogens           2/4/05
     x         label(i)(1:1) .ne. 'I') .and.   ! 2/4/05
     x        (ncon(i) .ne. 0)) go to 721                         ! 2/13/05
C----Halides anions...Cl-, Br-, I-
          if ((label(i)(1:1) .eq. 'X' .or.     ! if atom = X, Z   ! 2/4/05
     x         label(i)(1:1) .eq. 'Z' .or.     ! or I and         ! 2/4/05
     x         label(i)(1:1) .eq. 'I') .and.   ! ncon = 0         ! 2/4/05
     x        (ncon(i) .eq. 0)) then           ! assume a halide  ! 2/6/05
               n_good = n_good + 1                                ! 2/7/05
            write (45,730) label(i)(1:1), n_good                  ! 2/28/05
730             format (10x,'halide is ',a,', #',i2)              ! 2/28/05
            if (label(i)(1:1) .eq. 'Z') n_gZ = 1                  ! 2/14/05
            if (label(i)(1:1) .eq. 'X') n_gX = 1                  ! 2/14/05
            if (label(i)(1:1) .eq. 'I') n_gI = 1                  ! 2/14/05
             go to 721                                            ! 2/6/05
          endif                                                   ! 2/11/05
          if ((label(i)(1:1) .eq. 'X' .or.       ! single linked   2/4/05
     x         label(i)(1:1) .eq. 'Z' .or.       ! halogen         2/4/05
     x         label(i)(1:1) .eq. 'I') .and.     ! 2/4/05
     x        (ncon(i) .eq. 1)) then             ! 2/6/05
             go to 721                                            ! 2/11/05
          endif                                                   ! 2/11/05
788       ITEST = I                                               ! 2/11/05
          RETURN
721    CONTINUE
      if (n_good .gt. 0) n_zxi = n_zxi + 1                        ! 2/7/05
      if (n_gZ .ne. 0) n_Z = n_Z + 1                            ! 2/14/05
      if (n_gX .ne. 0) n_X = n_X + 1                            ! 2/14/05
      if (n_gI .ne. 0) n_I = n_I + 1                            ! 2/14/05
      RETURN
      END
C
C-------------FUNCTION NPOT------------------------------------
C
C----Determine atom types for volume additivity calculations
C
      FUNCTION NPOT (NUM_ATOM,ATYPE)                              ! 4/20/00
      CHARACTER*38 ATYPE, BLANK                                   ! 2/28/05
      CHARACTER NAME*5, YN*1, MMOD(64)*2, AT2*2, LINE*60
      CHARACTER COMMENT*120, REFCOD_2*8
      DIMENSION IK(500), NSIZES(4), NUM_O(3), DB(3)           ! 10-23-96
      COMMON /ZERO_COUNTS/ NCA(6), NNL(4), NOL(2),            ! 2/28/05
     2                     nsu(6), nph(6), nb(4),
     3                     n_hy,   nfl,    ncl,
     4                     nbr,    ni,     nclm,
     5                     nbrm,   nim,    nclo4
C                          C       N       O
C                          S       P       B
C                          H       F       CL
C                          BR      I       CL-
C                          BR-     I-      CLO4               ! 2/28/05
      COMMON /COM1/ NCON(500), ICON(7,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500),KEEP_NPOT(500)
      COMMON/T1/IAT_USED(500),NAROMATIC
      COMMON /GROUPS/ ICOORD
      COMMON /PROBLEMS/ JPROB, COMMENT, REFCOD_2, ICSP3, JCSP3
      DATA BLANK /'                                      '/       ! 2/28/05
C
10    AT2 = NAME(NUM_ATOM)(1:2)
      NPOT = 0
      ATYPE = BLANK                                               ! 4/20/00
C
C----look at H's first....................................
C
      IF (AT2(1:1) .NE. 'H') GO TO 1340     ! H
C
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 4/20/00
         GO TO 1303
C
C----Look at C's second.......................................
C
1340  IF (AT2(1:1) .NE. 'C') GO TO 2340
C
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 4/20/00
         GO TO 1303
C
C----Look at O's third.................................................
C
2340  IF (AT2(1:1) .NE. 'O') GO TO 1024
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 4/20/00
         GO TO 1303
C
C----Look at N's fourth...............................................
C
1024  IF (AT2(1:1) .NE. 'N') GO TO 2008
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 4/20/00
         GO TO 1303
C
C----Look at F's fifth........................................
2008  IF (AT2(1:1) .NE. 'F') GO TO 3007                      ! 4/20/00
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 5/1/00
         GO TO 1303
C
C----Look at CL's (CL = Z) sixth..............................
3007  IF (AT2(1:1) .NE. 'Z') GO TO 4007                      ! 4/20/00
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 5/1/00
         GO TO 1303
C
C----Look at BR's (BR = X) seventh...........................
4007  IF (AT2(1:1) .NE. 'X') GO TO 5007                      ! 4/20/00
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 5/1/00
         GO TO 1303
C
C----Look at S's eigth......................................
5007  IF (AT2(1:1) .NE. 'S') GO TO 6007                      ! 4/20/00
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
             IF (NPOT .EQ. 999) GO TO 2007                   ! 5/1/00
         GO TO 1303
C
C----Look at P's ninth......................................
6007  IF (AT2(1:1) .NE. 'P') GO TO 7007                      ! 4/3/04
          NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/20/00
              IF (NPOT .EQ. 999) GO TO 2007                   ! 5/1/00
          GO TO 1303
C
C----Look at I's tenth......................................   4/3/04
7007  IF (AT2(1:1) .NE. 'I') GO TO 8007                      ! 4/3/04
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/3/04
             IF (NPOT .EQ. 999) GO TO 2007                   ! 4/3/04
         GO TO 1303                                          ! 4/3/04
C
C----Look at B's 11th......................................   2/28/05
8007  IF (AT2(1:1) .NE. 'B') GO TO 2007                      ! 2/22/05
         NPOT = I_ATOM_COUNT (NUM_ATOM, ATYPE)               ! 4/3/04
             IF (NPOT .EQ. 999) GO TO 2007                   ! 4/3/04
         GO TO 1303                                          ! 4/3/04
C
2007   PRINT 2011, NAME(NUM_ATOM)
2011   FORMAT (' Atom ',A5,' is unknown...type set to 999')
1034      NPOT = 999     ! 4 connections ... flag as unknown
c1303   print 12345, num_atom, name(num_atom), atype
c12345  format ('-->',i3,2x,a5,2x,a38)                       ! 2/28/05
1303   continue
       RETURN
       END
C
C--------------FUNCTION I_ATOM_COUNT----------------------------
C
      FUNCTION I_ATOM_COUNT (NUM_ATOM, ATYPE)                 ! 6/12/05
      CHARACTER ATYPE*38                                      ! 2/28/05
      CHARACTER AT2*2
      CHARACTER NAME*5
      COMMON /ZERO_COUNTS/ NCA(6), NNL(4), NOL(2),            ! 2/28/05
     2                     nsu(6), nph(6), nb(4),
     3                     n_hy,   nfl,    ncl,
     4                     nbr,    ni,     nclm,
     5                     nbrm,   nim,    nclo4
C                          C       N       O
C                          S       P       B
C                          H       F       CL
C                          BR      I       CL-
C                          BR-     I-      CLO4-
      COMMON /COM1/ NCON(500), ICON(7,500), NAME(500),
     x          OR_XYZ(3,500), CELL(6), NA, IT_TABLE(500),
     x          I3D(500),KEEP_NPOT(500)
C
C----What is NUM_ATOM, what atoms (and links) are connected?     ! 4/20/00
C
      I_ATOM_COUNT = 0                                   ! 4/20/00
      AT2 = NAME(NUM_ATOM)                               ! 4/20/00
      ATYPE(1:1) = AT2(1:1)                              ! 4/20/00
      IN = NCON(NUM_ATOM)                                ! 4/20/00
      CALL INITIALIZE  ! initialize NCA, NNL, NOL,       ! 2/28/05
C                        NFL, N_HY, NCL, NBR, NSU,       ! 2/4/05
C                        NPH, NI, NCLM, NBRM, NIM        ! 2/4/05
C                        nclo4, NB                       ! 2/28/05

      k = num_atom                                       ! 6/12/05
      if (in .gt. 0) go to 1298                          ! 6/12/05
         if (name(k)(1:1) .eq. 'Z') nclm = nclm + 1 ! CL-! 2/4/05
         if (name(k)(1:1) .eq. 'X') nbrm = nbrm + 1 ! BR-! 2/4/05
         if (name(k)(1:1) .eq. 'I') nim = nim + 1   ! I- ! 2/4/05
         go to 1297                                      ! 6/12/05
1298  if (in .eq. 4 .and. name(k)(1:1) .eq. 'Z') then    ! 6/12/05
         nclo4 = nclo4 + 1                               ! 2/10/05
         go to 1296                                      ! 6/12/05
      endif                                              ! 2/10/05
1300       DO 1345 I=1,IN                                ! 2/4/05
              K = ICON(I,NUM_ATOM)                       ! 4/20/00
              M = I_ATOM(NAME(K)(1:2))                   ! 4/20/00
                 IF (M .EQ. 999) GO TO 102               ! 4/20/00
              N = NCON(K)                                ! 4/20/00
              IF (NAME(K)(1:1) .EQ. 'C') THEN            ! 4/20/00
                 NCA(N) = NCA(N) + 1                     ! 5/1/00
                 GO TO 1345                              ! 4/20/00
              ENDIF                                      ! 4/20/00
              IF (NAME(K)(1:1) .EQ. 'N') THEN            ! 4/20/00
                 NNL(N) = NNL(N) + 1                     ! 4/20/00
                 GO TO 1345                              ! 4/20/00
              ENDIF                                      ! 4/20/00
              IF (NAME(K)(1:1) .EQ. 'O') THEN            ! 4/20/00
                 NOL(N) = NOL(N) + 1                     ! 4/20/00
                 GO TO 1345                              ! 4/20/00
              ENDIF
              IF (NAME(K)(1:1) .EQ. 'S') THEN            ! 5/2/00
                 NSU(N) = NSU(N) + 1                     ! 5/2/00
                 GO TO 1345                              ! 5/2/00
              ENDIF                                      ! 5/2/00
              IF (NAME(K)(1:1) .EQ. 'P') THEN            ! 5/2/00
                 NPH(N) = NPH(N) + 1                     ! 5/2/00
                 GO TO 1345                              ! 5/2/00
              ENDIF                                      ! 5/2/00
              if (name(k)(1:1) .eq. 'B') then            ! 2/22/05
                 nb(n-2) = nb(n-2) + 1  ! save ncon's of ! 2/22/05
                 go to 1345             ! 3-6 in nb(1) to! 2/22/05
              endif                     ! nb(4)          ! 2/22/05
              IF (NAME(K)(1:1) .EQ. 'H') N_HY = N_HY + 1 ! 4/20/00
              IF (NAME(K)(1:1) .EQ. 'F') NFL = NFL + 1   ! 4/20/00
              IF (NAME(K)(1:1) .EQ. 'Z') NCL = NCL + 1   ! 5/1/00
              IF (NAME(K)(1:1) .EQ. 'X') NBR = NBR + 1   ! 5/1/00
              IF (NAME(K)(1:1) .EQ. 'I') NI = NI + 1     ! 4/3/04
1345       CONTINUE                                      ! 4/20/00
C----Build ATYPE contents                                ! 4/20/00
1346    IPOS = 1                 ! C                     ! 2/10/05
           DO 1351 I=1,6                                 ! 2/22/05
              IF (NCA(I) .EQ. 0) GO TO 1351              ! 4/20/00
              J = IPOS + I     ! position in ATYPE       ! 4/20/00
              ATYPE(J:J) = ACHAR(O'060' + NCA(I))        ! 4/20/00
1351       CONTINUE                                      ! 4/20/00
        IPOS = 7                 ! N                     ! 2/22/05
           DO 1352 I=1,4                                 ! 4/20/00
              IF (NNL(I) .EQ. 0) GO TO 1352              ! 4/20/00
              J = IPOS + I     ! position in ATYPE       ! 4/20/00
              ATYPE(J:J) = ACHAR(O'060' + NNL(I))        ! 4/20/00
1352       CONTINUE                                      ! 4/20/00
        IPOS = 11                ! O                     ! 2/22/05
           DO 1353 I=1,2                                 ! 4/20/00
              IF (NOL(I) .EQ. 0) GO TO 1353              ! 4/20/00
              J = IPOS + I     ! position in ATYPE       ! 4/20/00
              ATYPE(J:J) = ACHAR(O'060' + NOL(I))        ! 4/20/00
1353       CONTINUE                                      ! 4/20/00
        IPOS = 13                ! S                     ! 2/22/05
           DO 1354 I=1,6                                 ! 5/3/00
              IF (NSU(I) .EQ. 0) GO TO 1354              ! 5/3/00
              J = IPOS + I     ! position in ATYPE       ! 5/3/00
              ATYPE(J:J) = ACHAR(O'060' + NSU(I))        ! 5/3/00
1354       CONTINUE                                      ! 5/3/00
        IPOS = 19                ! P                     ! 2/22/05
           DO 1356 I=1,6                                 ! 2/22/05
              IF (NPH(I) .EQ. 0) GO TO 1356              ! 5/3/00
              J = IPOS + I     ! position in ATYPE       ! 5/3/00
              ATYPE(J:J) = ACHAR(O'060' + NPH(I))        ! 5/3/00
1356       CONTINUE                                      ! 5/3/00
        IPOS = 25                ! B                     ! 2/22/05
           DO 1358 I=1,4  ! 1-4 = 3-6 in terms of connectivity ! 2/22/05
              IF (NB(I) .EQ. 0) GO TO 1358               ! 2/23/05
              J = IPOS + I     ! position in ATYPE       ! 5/3/00
              ATYPE(J:J) = ACHAR(O'060' + NB(I))         ! 2/23/05
1358       CONTINUE                                      ! 2/22/05
        IF (N_HY .NE. 0) ATYPE(30:30) = ACHAR(O'060' + N_HY) ! 2/25/05
        IF (NFL .NE. 0) ATYPE(31:31) = ACHAR(O'060' + NFL)   ! 2/25/05
        IF (NCL .NE. 0) ATYPE(32:32) = ACHAR(O'060' + NCL)   ! 2/25/05
        IF (NBR .NE. 0) ATYPE(33:33) = ACHAR(O'060' + NBR)   ! 2/25/05
        IF (NI .NE. 0) ATYPE(34:34) = ACHAR(O'060' + NI)     ! 2/25/05
1297    if (nclm .ne. 0) atype(35:35) = achar(O'060' + nclm) ! 6/12/05
        if (nbrm .ne. 0) atype(36:36) = achar(O'060' + nbrm) ! 2/25/05
        if (nim .ne. 0) atype(37:37) = achar(O'060' + nim)   ! 2/25/05
1296    if (nclo4 .ne. 0) atype(38:38) = achar(O'060' + nclo4) ! 6/12/05
      RETURN                                             ! 4/20/00
C----Problem section
102   WRITE (6,1400) K, NAME(K), NUM_ATOM, NAME(NUM_ATOM) ! 4/20/00
1400  FORMAT ('****Atom #',I3,' (',A5,') linked to atom #', ! 4/20/00
     2     I3,' (',A5,') is unknown')                    ! 4/20/00
1333  I_ATOM_COUNT = 999                                 ! 4/20/00
      RETURN                                             ! 4/20/00
      END                                                ! 4/20/00
C
C-------------FUNCTION I_ATOM-------------------------------------
C
      FUNCTION I_ATOM(AT)
      PARAMETER (NTYPES_ATOM = 11)                               ! 2/21/05
      CHARACTER A_ORDER(NTYPES_ATOM)*2, AT*2                     ! 2/2/05
      DIMENSION IPOS(NTYPES_ATOM)                                ! 2/2/05
      DATA A_ORDER /'C ', 'N ', 'O ', 'S ', 'P ',
     2              'B ', 'H ', 'F ', 'Z ', 'X ', 'I '/          ! 2/23/05
C                                      CL    BR
      DATA IPOS /1, 7, 11, 13, 19, 25, 29,                       ! 2/23/05
C                C  N  O   S   P   B   H                         ! 2/23/05
     2          30, 31, 32, 33/                         ! position 2/23/05
C                F  CL  BR  I
C----Use A_ORDER array to assign a number (1, 2 or 3) to          ! 4/20/00
C     atom type (C, N, O)
C                                                                 ! 4/20/00
C     print 1122, at               !##################
C1122  format ('***at =','|',a2,'|') !############
      I_ATOM = 0                                                  ! 4/24/00
      DO 20 J=1,NTYPES_ATOM                                       ! 2/2/05
C       print 1133, at, j, a_order(j), ipos(j)    ! ##########
C1133   format ('***at, j, a_order(j), ipos(j) = ',2('|',a2,'|',i3))  ! #########
         IF (AT(1:1) .EQ. A_ORDER(J)(1:1)) THEN                             ! 4/20/00
            I_ATOM = IPOS(J)                                      ! 4/20/00
            RETURN                                                ! 4/20/00
         ENDIF                                                    ! 4/20/00
20    CONTINUE                                                    ! 4/20/00
      I_ATOM = 999                                                ! 4/24/00
      RETURN                                                      ! 2/2/05
      END                                                         ! 4/20/00
C
C------------FUNCTION INITIALIZE-----------------------------
C
       SUBROUTINE INITIALIZE                              ! 4/20/00
       DIMENSION NALL(37)                                 ! 2/22/05
      COMMON /ZERO_COUNTS/ nall                           ! 2/22/05
C                          C     N     O     F    H
C                          Cl    Br    S     P    I
C                          Cl-   Br-   I-   CLO4-  B         ! 2/22/05
C
       DO I=1,37                                          ! 2/22/05
          NALL(I) = 0                                     ! 4/20/00
       ENDDO                                              ! 4/20/00
       RETURN                                             ! 4/20/00
       END                                                ! 4/20/00
C
C---------------------------------------------------------------------
      BLOCK DATA STUFF
      PARAMETER (NTYPES_ATOMS=11, NPOTS=3000)                     ! 2/21/05
      CHARACTER ALLTYPES(NPOTS)*(38), ATYPE*(38),                 ! 2/22/05
     x          BLANK*(38)                                        ! 2/22/05
      COMMON /TYPES/ ALLTYPES, N_TYPES, BLANK                     ! 4/3/04
      CHARACTER ATSYM*1, ATSYM2*2                                 ! 2/26/00
      COMMON /CONTENTS/ NCON_FDAT(NTYPES_ATOMS),                  ! 2/26/00
     x                  NCON_FBIB(NTYPES_ATOMS),                  ! 2/26/00
     x                  MAX_CONN(NTYPES_ATOMS),                   ! 5/8/00
     x                  ATSYM(NTYPES_ATOMS), n_zxi, n_clo4,       ! 2/11/05
     X                  n_Z, n_X, n_I                             ! 2/14/05
      data n_zxi, n_Z, n_X, n_I, n_clo4 /5*0/                     ! 2/14/05
      DATA ATSYM    /'C','H','F','N','O','X','Z','S','P','I','B'/ ! 2/21/05
C                                         BR  CL                  ! 3/4/00
      DATA MAX_CONN / 6,  1,  1,  4,  2,  1,  4,  6,  6,  1,  6/  ! max allowed connections, will only keep
C                     C   H   F   N   O   Br  Cl  S   P   I   B   ! 0, 1 and 4 for Cl; 3-6 for B
      COMMON /ELEMENTS/ N_AT_DIFF, N_PER_UNIT(10), ATSYM2(NTYPES_ATOMS) ! 12/27/04
      DATA N_PER_UNIT, N_AT_DIFF /11*0/                           ! 12/27/04
      DATA ATSYM2 /'C ', 'H ', 'F ', 'N ', 'O ', 'BR', 'CL', 'S ',! 3/11/00
     x             'P ', 'I ', 'B '/                              ! 2/21/05
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
      COMMON /COM1/ NCON(500), ICON(7,500), NAME(500),
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
C
C-----------FUNCTION DISTANCE----------------------------------------
C
C----Calculate distances between 2 atoms
C
      FUNCTION DISTANCE (I1, I2)
      CHARACTER LABEL*5
C
      COMMON /COM1/ NCON(500), ICON(7,500), LABEL(500),
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
