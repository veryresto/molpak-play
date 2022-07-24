      PROGRAM MAKE_FILES
!     The latest version 09-04-2014
!     updated it for repalcing DMAREL with DMACRYS       1-29-09    DU
!     updated on 7-9-2008 for adding PMIN                           DU
!     54 geoms codes added                               9-18-07    DU
!     updated...added "limit cputime for running wmin"   2-13-07    DU
!     Updated...add ability to read cross-term coef information for
!        PMIN refinement in molpak-pmin calcns                       6/24/09   HLA
!
!----Create all files necessary of MOLPAK + WMIN or PMIN and MOLPAK +
!           DMACRYS structure predictions.
!    For WMIN refinement, procedure requires only model coordinate
!           files with charges.  Use of DMACRYS is sufficently more
!           complicated and requires four additional input files.
!
!
      USE make_filesCommonMod
!
      IMPLICIT NONE
!
!----Set initial values for symmetry flags
      INVER2 = 0
      IMR = 0
      IAX = 0
!----NSOL is the solution #
      DO I=1,26
         BLNK(I) = BLANK
      ENDDO
!----Heading information
      WRITE (*,60)
60    FORMAT (63(1H#)/'#',61X,'#'/'#Create UNIX script for crystal', &
     &       ' structure/density predictions #'/ &
     &       '#',18X,'version of Sept. 2014',19X,'#'/ &
     &       '#',61X,'#')
      WRITE (*,61)
61    FORMAT ('#',10X,'Space groups/MOLPAK coordination geometries',8X,'#'/ &
     & '#',10X,'that have been extensively tested...       ',8X,'#'/ &
     &  '#',9X,'       Sg      cd geom   cd geom   cd geom ',9X,'#'/ &
     &  '#',9X,' tric  P1      AA                          ',9X,'#'/ &
     &  '#',9X,' tric  P-1     AB        CA                ',9X,'#'/ &
     &  '#',9X,' mono  P2      AE                          ',9X,'#'/ &
     &  '#',9X,' mono  Pm      AC                          ',9X,'#'/ &
     &  '#',9X,' mono  Pc      AD        AG                ',9X,'#'/ &
     &  '#',9X,' mono  P21     AH        AF                ',9X,'#'/ &
     &  '#',9X,' mono  P2/c    AJ        AL        FD      ',9X,'#'/ &
     &  '#',9X,' mono  P21/m   AN                          ',9X,'#'/ &
     &  '#',9X,' mono  P2/m    AO        FB                ',9X,'#'/ &
     &  '#',9X,' mono  P21/c   AI        AK        AM      ',9X,'#'/ &
     &  '#',9X,' mono  P21/c   FA        FC                ',9X,'#'/ &
     &  '#',9X,' mono  Cc      DA                          ',9X,'#'/ &
     &  '#',9X,' mono  C2      DB                          ',9X,'#'/ &
     &  '#',9X,' mono  C2/c    DC        DD        DE      ',9X,'#'/ &
     &  '#',9X,' orth  Pnn2    AR        AT        BE      ',9X,'#'/ &
     &  '#',9X,' orth  Pba2    AW        BG                ',9X,'#'/ &
     &  '#',9X,' orth  Pnc2    AX        BI                ',9X,'#'/ &
     &  '#',9X,' orth  P221    BC                          ',9X,'#'/ &
     &  '#',9X,' orth  Pmn21   BJ                          ',9X,'#'/ &
     &  '#',9X,' orth  Pma2    BK                          ',9X,'#'/ &
     &  '#',9X,' orth  P21212  AP        BA        BB      ',9X,'#'/ &
     &  '#',9X,' orth  P212121 AQ        AZ                ',9x,'#'/ &
     &  '#',9X,' orth  Pca21   AY        BH                ',9X,'#'/ &
     &  '#',9X,' orth  Pna21   AV        BD        BF      ',9X,'#'/ &
     &  '#',9X,' orth  Pna21   AU        AS                ',9X,'#'/ &
     &  '#',9X,' orth  Pnma    CF (on mirror)              ',9X,'#'/ &
     &  '#',9X,' orth  Fdd2    DF        DG  (on 2 fold)   ',9X,'#'/ &
     &  '#',9X,' orth  Pbcn    CD        CE                ',9X,'#'/ &
     &  '#',9X,' orth  Pbca    CC        CB                ',9X,'#'/ &
     &  '#',61X,'#'/ &
     &  '#',' MOLPAK assumes that a whole molecule is supplied', &
     &  ' regardless #'/'#',' of any internal molecular symmetry ', &
     &   25X,'#'/63(1H#))
!----Request 1st line...title for output
      WRITE (*,62)
62    FORMAT (/'Title line: ',$)
      READ (*,337) TITLE
337   FORMAT (A)
!----Read 2nd line...compound ID
      WRITE (*,33)
33    FORMAT('Compound ID (8 char max): ',$)
      READ (*,338) COMPDID
338   FORMAT(A8)                             ! 1-29-09
!----How many characters in compound ID?
      NID = 0
      DO I=1,8
         IF(COMPD_ID(I) .EQ. BLANK) GO TO 1332
         NID = NID + 1
         MASTER_NAME(I) = COMPD_ID(I)
      ENDDO
!----Read 3rd line...2 character space group designation
1332  WRITE (*,70 )
70    FORMAT('Give MOLPAK coord geom symbol(s), separate symbols with'/  &
     &       '   spaces, commas or nothing (40 chars max); symbols/31/54 [default=31]: ',$)     ! 9-18-07
      IALL = 0
      READ (5,21) SPGR
21    FORMAT(40A1)
!----Some coord geoms can be eliminated.
      IF (SPGR(1) .EQ. ' ' .OR. SPGR(1) .EQ. '3'.OR. SPGR(1) .EQ. '5') THEN
         WRITE (*,210)
210      FORMAT ('Do you want to eliminate some coord geoms? [N]: ',$)
         READ (5,337) Y_N   ! some coord geoms can be eliminated
      END IF
      IF (Y_N .EQ. 'Y' .OR. Y_N .EQ. 'y') THEN
         WRITE (*,213)
213      FORMAT(' Give coord geom code(s) eliminate'/ &
     &         '   (example =  DA DB CB CC; spaces are required):')
         READ (5,212) ELM_SPGR
212      FORMAT (50(A2,1X))
      END IF
!----Select program for refinement (WMIN or PMIN or DMACRYS)
3300  WRITE (*,3302)
3302     FORMAT ('Program for lattice refinement...1 = WMIN,', &
     &          ' 2 = PMIN, 3 = DMACRYS [1]: ',$)
      READ (*,337) R_KIND
      IF (R_KIND .EQ. '1' .OR. R_KIND .EQ. ' ') THEN
         KIND = 1
         GO TO 3155
      ELSE
       IF (R_KIND .EQ. '2') THEN                                         ! 7-9-08
         KIND = 2          ! not blank or 1 or 3  assume PMIN refinement ! 7-9-08
!        WRITE (*,3305)
!3305        FORMAT (' *** NOTE: The use of PMIN requires the following', &
!    &              ' 1 file ***'/ &
!    &              ' PMIN.inp')
!        OPEN (UNIT=12, FILE='PMIN.inp', STATUS='OLD', ERR=3203)
!        CLOSE(UNIT=12)
!        GO TO 3401
!3203     WRITE (*,3204)
!3204        FORMAT('*****PMIN file error, file availability problem - QUIT*****'/ &
!    &          '   PMIN.inp file indicated above was not found')
!        STOP
3401     CALL REF_PMIN
         STOP
       ELSE
         KIND = 3          ! not blank or 1 or 2  assume DMACRYS refinement
         WRITE (*,3304)
3304        FORMAT (' NOTE: The use of DMACRYS requires the following', &
     &              ' 4 files...'/ &
     &              ' pote.dat or will01.pots, cadpac.charges, dmarel.axis & bondlengths.')
         OPEN (UNIT=12, FILE='pote.dat', STATUS='OLD', ERR=3202)
         CLOSE (UNIT=12)
         GO TO 3205
3202     CLOSE (UNIT=12)
         OPEN (UNIT=12, FILE='will01.pots', STATUS='OLD', ERR=3200)
         CLOSE (UNIT=12)
3205     OPEN (UNIT=12, FILE='dmarel.axis', STATUS='OLD', ERR=3200)
         CLOSE (UNIT=12)
         OPEN (UNIT=12, FILE='bondlengths', STATUS='OLD', ERR=3200)
         CLOSE (UNIT=12)
         OPEN(UNIT=12, FILE='cadpac.charges',STATUS='OLD',ERR=3200)
         GO TO 3400
3200     WRITE (*,3201)
3201        FORMAT('**DMACRYS file error, file availability problem - QUIT**'/ &
     &          '   one or more of the 4 files indicated above was not found')
         STOP
3400     WRITE (*,1555)
1555     FORMAT ('The potential parameter file is...'/&
     &           '  1 = UMD PC, 2 = FIT  PC, 3 = Williams : ',$)
         READ(*,'(I1)') I_POTE
         CALL REF_DMACRYS
         STOP
       END IF
      ENDIF
!----Set IDCHG = 2...hold-over from original program.  Only one type
!     of charge is input with atom coordinates.
3155  IDCHG = 2
      POT_FILE = ' '
      WRITE (*,155)
155      FORMAT ('The potential parameter file is...'/&
     &           '  1 = UMD_potential, 2 = ARL_potential,',&
     &           ' 3 = other [1]: ',$)
      READ (*,337) POT_FILE
         SELECT CASE (POT_FILE)
            CASE (' ','1')
               POTENTIAL_FILE = "UMD_potential"    ! 7-9-08
            CASE ("2")
               POTENTIAL_FILE = "ARL_potential"
            CASE ("3")
               WRITE (*,160)
160            FORMAT ('Provide name of potential parameter file ',&
     &                 '(40 chars max) :',$)
               READ (*,337) POTENTIAL_FILE
         END SELECT
!     IF (GROUP .EQ. '  ' .OR. GROUP .EQ. 'AL' .OR. GROUP .EQ. 'al') &
      IF (GROUP .EQ. '  ' .OR. GROUP .EQ. '31' .OR. GROUP .EQ. '54') &        ! 9-18-07
     &      IALL = 1
!----4th line...question about COMPLETE output file
      WRITE (*,487)
487      FORMAT (' Should a complete output file be created and saved?',  &
     &        ' [N]: ', $)
      I_COMP = .true.             ! about a HUGE, complete output file
      I_COMPLETE = ' '
      READ (*,337) I_COMPLETE
      IF (I_COMPLETE .EQ. 'N' .OR. I_COMPLETE .EQ. 'n' .OR. &
     &       I_COMPLETE .EQ. ' ') I_COMP = .false.
!----5th line...question about the ....save file
      WRITE (*,588)
588      FORMAT (' Save the refined coordinate and cell parameter',&
     &           ' file (...save)? [N]: ', $)
      I_SAVE = .true.    ! true means keep the ...save file
      I_SAVED = ' '
      READ (*,337) I_SAVED
      IF (I_SAVED .EQ. 'N' .OR. I_SAVED .EQ. 'n' .OR. &
     &       I_SAVED .EQ. ' ') I_SAVE = .false.
!----Read 6th line...name of MOLPAK coordinate input file.  Try to
!     pick up file name automatically
      OPEN (UNIT=14, FILE='MOLPAK.NAME', STATUS='OLD')  ! look for files
      I = 0
45    READ (14,337,END=47) MOLPAK_NAME                   ! that have
      I = INDEX(MOLPAK_NAME,'molpak.xyz')               ! molpak.xyz
      IF (I .NE. 0) THEN
         L = INDEX(MOLPAK_NAME,'.xyz ') + 3
         GO TO 48
      ENDIF
      GO TO 45                                          ! at the end
47    IF (I .EQ. 0) THEN
         MOLPAK_NAME = 'molpak.xyz'
         L = 10
      ENDIF
48    WRITE (*,50) MOLPAK_NAME(1:L)
50    FORMAT('Name of the MOLPAK atom coord file', &
     &       ' [',A,']: ',$)
      READ (*,337) WHAT_LINE
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,337) ATOM_LIST
51       FORMAT (30A1)
      ELSE
         ATOM1_LIST = MOLPAK_NAME
      ENDIF
!----Inquire about possible 2-fold symmetry in model
      WRITE (*,486)
 486  FORMAT ('The following questions are about molecules with', &
     &        ' C2, Cs and Ci symmetry'/ &
     &       '   and will restrict searches to space groups with', &
     &       ' appropriate symmetries;'/ &
     &       '   reply N/RETURN for a general search that ignores', &
     &       ' molecular symmetry...'/ &
     &       '    Does the molecule contain two-fold symmetry?', &
     &       ' (Y/N [N]): ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. ' ' .OR. Y_N .EQ. 'N' .OR. Y_N .EQ. 'n')  THEN
         WHAT_LINE(1:1) = ' '
      ELSE
         WRITE (*,4861)
 4861    FORMAT('    Which axis is parallel to two-fold axis?  '/ &
     &          '      1 = X, 2 = Y, 3 = Z: ',$)
         READ (*,337) WHAT_LINE
      END IF
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,3303) IAX
3303        FORMAT (I1)
         IF (IAX .GE. 1 .AND. IAX .LE. 3) THEN
            ELMSPGR(1:40) = 'AAABCAAHAFAIAKAMFAFCDA&
                           &AQAZAYBHAVBDBFCCCB'
         ELSE
            IAX = 0
         ENDIF
      ENDIF
      IF (IAX .NE. 0) GO TO 505 ! if axis, forget mirror & inversion question
!----Inquire about possible mirror plane symmetry in model
      WRITE (*,4862)
4862     FORMAT('    Does the molecule contain mirror symmetry? (Y/N', &
     &       ' [N]): ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. ' ' .OR. Y_N .EQ. 'N' .OR. Y_N .EQ. 'n')  THEN
        WHAT_LINE(1:1) = ' '
      ELSE
         WRITE (*,4863)
 4863    FORMAT('    Which axis is perpendicular to the mirror plane?  '/ &
     &          '      1 = X, 2 = Y, 3 = Z: ',$)
         READ (*,337) WHAT_LINE(1:1)
      END IF
      IMR = 0
      IF (WHAT_LINE .NE. INPUT_LINE)  READ (WHAT_LINE,3303) IMR
      IF (IMR .NE. 0) GO TO 505  ! if mirror, forget inversion question
!----Inquire about possible inversion symmetry in model
      WRITE (*,496)
496      FORMAT('    Does the molecule contain a center of symmetry', &
     &          ' at the origin? (Y/N [N]): ',$)
      READ (5,337) INVER1
      IF (INVER1 .EQ. ' ' .OR. INVER1 .EQ. 'N' .OR. &
     &    INVER1 .EQ. 'n') THEN
          INVER2 = 0
      ELSE
         ELMSPGR(1:30) = 'AAAHAFDADBAPBABBAQAZAYBHAVBDBF'
         INVER2 = 1
      ENDIF
!----Open atom coordinate file
505   OPEN (UNIT=15, FILE=ATOM1_LIST, STATUS='OLD')
      WRITE (*,109)
109      FORMAT (4X,'atom id     x         y         z      Id    ', &
     &        '6-31g*')
      N = 0
111   N = N + 1
        READ (15,1,END=100) (ATOM_TY(N,J),J=1,2), (ATOM(N,I),I=1,4), &
     &                    X(N), Y(N), Z(N), IDATOM(N), G92CHARGE(N)
1          FORMAT (5X,2A1,4A1,3F10.6,I5,F10.6)
        WRITE (*,112) (ATOM_TY(N,J),J=1,2), (ATOM(N,I),I=1,4), &
     &                X(N), Y(N), Z(N), IDATOM(N), G92CHARGE(N)
112        FORMAT(4X,2A1,4A1,3F10.6,I5,F10.6)
      GO TO 111
100   N = N -1
!
!----Question about using an existing solution
      NSOL = 0             ! solution #
      Y_N = ' '
      WRITE (*,10)
10       FORMAT ('Do you want to rerun a solution from an existing', &
     &           ' volume.min file [N]: ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. 'Y' .OR. Y_N .EQ. 'y') THEN
          WRITE (*,12)
12           FORMAT (' Enter solution number: ',$)
          READ (*,*) NSOL
          WHAT_LINE = '1'
          OPEN (UNIT=30, FILE='volume.min', STATUS="OLD", ERR=800)
          GO TO 801
800       WRITE (*,804)
804          FORMAT (' The required file volume.min is not exist!')
          STOP
801       READ (30,*)
          READ (30,'(9X,3F8.1)') (ZSTEP(I), I=1,3)
          DO L=1,3
             ZMINA(L) = -90.
             ZMAXA(L) =  90.
          END DO
      ELSE
!----Read 7th line...angle search limits for initial molpak search
          WRITE (*,80)
80            FORMAT('Minimum, maximum and step for 3 angles for 1st ', &
     &               'MOLPAK map'/ &
     &               '  [-90,90,10,-90,90,10,-90,90,10,]: ',$)
          READ (*,337) WHAT_LINE
          IF (WHAT_LINE .NE. INPUT_LINE) THEN
              READ (WHAT_LINE,*) (ZMINA(I),ZMAXA(I),ZSTEP(I),I=1,3)
          ELSE
!----Set up defaults...
              DO L=1,3
                 ZMINA(L) = -90.
                 ZMAXA(L) =  90.
                 ZSTEP(L) =  10.
              ENDDO
          ENDIF
!----Read number of minimum volume points to be saved in MOLPAK search
             WRITE (*,703)
703             FORMAT ('Number of minimum volumes from initial MOLPAK', &
     &                  ' search to pass to next step,'/  &
     &                  '   suggest 7000 for 10 deg step search and', &
     & ' 60000 for 5 deg step search or range(-180 to 180) [7000]: ',$)
             WHAT_LINE = INPUT_LINE
             READ (*,337) WHAT_LINE
             IF (WHAT_LINE .EQ. INPUT_LINE) THEN
                N_VOLS_SAVED = 7000
             ELSE
                READ (WHAT_LINE,*) N_VOLS_SAVED
             ENDIF
!----Read number of minimum volume points to be elaborated
             WRITE (*,803)
803             FORMAT('Maximum number of 2nd MOLPAK calcns desired [500]: ',$)
             WHAT_LINE = INPUT_LINE
             READ (*,337) WHAT_LINE
      END IF
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,*) NMOL
      ELSE
         NMOL = 500       ! default -> elaborate 500 solns
      ENDIF
!----How many coord geoms have been specified?  Each 2 letter code
!     should be separated from its neighbors by a space or comma.
      IF (IALL .NE. 0) GO TO 29
      N_SPGRPS = 0
      NSCHARS = 0
      DO 22 I=1,40
         IF (SPGR(I) .EQ. ' ' .OR. SPGR(I) .EQ. ',') GO TO 22
         IF (NSCHARS .NE. 0) THEN
            IF (NSCHARS .EQ. 1) GO TO 24
            NSCHARS = 0
         ENDIF
23       N_SPGRPS = N_SPGRPS + 1
24       NSCHARS = NSCHARS + 1
         MOL_SPGRP(NSCHARS,N_SPGRPS) = SPGR(I)
22    CONTINUE
!----All or just some space groups
29    IF (IALL .EQ. 0) THEN
         N_GO = N_SPGRPS
      ELSE
       IF (GROUP .EQ. '  ' .OR. GROUP .EQ. '31') THEN    ! 9-18-07
         N_GO = N_SPGR-23                                ! 9-18-07
       ELSE
         N_GO = N_SPGR
       END IF
      ENDIF
!
!----Big LOOP for creating the N_GO sets of files for execution
!
Big_loop : DO LGRPS=1,N_GO
               IF (IALL .EQ. 0) THEN
                  SPGR(1) = MOL_SPGRP(1,LGRPS)
                  SPGR(2) = MOL_SPGRP(2,LGRPS)
               ELSE
                  DO J=1,N_GO
                     IF (ALL_SPGR(LGRPS) .EQ. ELM_SPGR(J)) CYCLE Big_loop
                  END DO
                  SPGR(1) = ALL_SPGR(LGRPS)(1:1)
                  SPGR(2) = ALL_SPGR(LGRPS)(2:2)
               ENDIF
!----Go through the N_SPGRPS coord geoms & make appropriate files.
!     Tack 3 characters (_SPGR) onto end of COMPD_ID and MASTER_NAME.
      COMPD_ID(NID+1) = '_'
      COMPD_ID(NID+2) = OLD(1)
      COMPD_ID(NID+3) = OLD(2)
      COMPD_ID(NID+4) = OLD(3)
      COMPD_ID(NID+5) = '_'
      COMPD_ID(NID+6) = SPGR(1)
      COMPD_ID(NID+7) = SPGR(2)
      MASTER_NAME(NID+1) = '_'
      MASTER_NAME(NID+2) = SPGR(1)
      MASTER_NAME(NID+3) = SPGR(2)
!----Construct name for G98 file
      DO I=1,NID
         G92FILE(I)=COMPD_ID(I)
      END DO
      G92FILE(NID+1) = '_'
      G92FILE(NID+2) = G92(1)
      G92FILE(NID+3) = G92(2)
      G92FILE(NID+4) = G92(3)
      G92FILE(NID+5) = '_'
      G92FILE(NID+6) = SPGR(1)
      G92FILE(NID+7) = SPGR(2)
!----Construct name for MASTER.COM file that will have compound ID 1st,
!     space group symbol 2nd and "master.com" 3rd.
      MID = NID + 3
      DO I=1,11
         MASTER_NAME(MID+I) = MASTERCOM(I)
      ENDDO
      MOLPAK_FILE(1:7) = 'molpak_'      ! create name of the
      MOLPAK_FILE(8:9) = GROUP          ! molpak input file
      MOLPAK_FILE(10:15) = '.input'
      SEARCH_FILE(1:7) = 'search_'      ! create name for search
      SEARCH_FILE(8:9) = GROUP          ! data file
      SEARCH_FILE(10:14) = '.data'
!----Open files
      OPEN (UNIT=20, FILE=MASTER_NAME2, STATUS='UNKNOWN')
      OPEN (UNIT=25, FILE=MOLPAK_FILE,  STATUS='UNKNOWN')
!----Start to write master.com... files
      WRITE (20,39) (NMOL+1)
39    FORMAT ('#!/bin/csh'/'pwd'/'#'/'time'/'date'/'#'/ &
     &   'setenv dir /c/Users/ACER/Downloads/PREDICTIONS'/ &
     &   '#mkdir WMININP_LS WMININP_MPK'/&
     &   '#'/'#****initial MOLPAK search****'/'#'/ &
     &   '# if input.* exists, restart WMIN'/ &
     &   'set k=1'/'while ($k != ',I4,')'/ &
     &   '  if ( -e input.$k ) then'/'    goto runwmin'/'  endif'/ &
     &   '  @ k = ($k + 1)'/'end'/ &
     &   '# if volume.min  is existing, skip first MOLPAK search'/ &
     &   'if ( -e volume.min  ) then'/'   cp -f volume.min  fort.11'/ &
     &   '   goto mkcom'/'else'/'   if ( -e volume.find ) goto mkcom'/ &
     &   'endif'/ &
     &   '# if molpak.vol exists, skip MOLPAK search'/&
     &   'if ( -e molpak.vol ) cp -f molpak.vol molpak.min'/'#')
      WRITE (20,41) MOLPAK_FILE
41    FORMAT ('cp -f ',A15,'  fort.15       # MOLPAK input file'/ &
     &       '$dir/MOLPAK/molpak_12.exe',' # run MOLPAK'/ &
     &       'if ( -e molpak.min ) rm -f molpak.min'/ &
     &       'rm -f molpak.vol'/ &
     &       'time'/'date')
      NDIF = 30 - MID
!----Write title on MOLPAK input file + rest of it
      WRITE (25,245) COMPDID(1:NID), TITLE
245   FORMAT ('HEAD  ',A,8X,A)
!     using default  2-27-04
!     WRITE (25,146)                      ! line # 2...NSEG values
146   FORMAT ('NSEG  128 128   8')
!----Place atom coordinates into molpak input file
      DO L=1,N
         WRITE(25,483) (ATOM_TY(L,J),J=1,2), (ATOM(L,I),I=1,4), &
     &                 X(L), Y(L), Z(L), IDATOM(L), G92CHARGE(L)
483         FORMAT('ATOM ',2A1,4A1,3F10.6,I5,F10.6)
      ENDDO
      IF (INVER2 .EQ. 0) GO TO 7631 ! INVER2 = 0  means no center
      WRITE (25,764)                ! write "CENT" line into
764       FORMAT ('CENT',30X,' ')   ! MOLPAK input file
      GO TO 488
 7631 IF (IMR .EQ. 0) GO TO 763    ! IMR = 0 means no mirror
      WRITE (25,7632) IMR          ! write "PLAN" line into
 7632    FORMAT ('PLAN',I2,28X,' ') ! MOLPAK input
      GO TO 488
763   IF (IAX .EQ. 0) GO TO 488    ! IAX = 0 means no 2-fold axis
      WRITE (25,787) IAX           ! write "AXIS" line into
787      FORMAT ('AXIS',I2,28X,' ')! MOLPAK input file
488   IF (I_COMP)  WRITE (20,491) (COMPD_ID(I),I=1,MID), COMPLT
491                FORMAT ('mv -f fort.16 ',46A1)
      WRITE(25,246) GROUP
246   FORMAT('INCL ',A2)
      IF (IAX .EQ. 0) THEN
         WRITE (25,147) N_VOLS_SAVED  ! # min vols for 1st MOLPAK
147         FORMAT ('VOLS ',I5,'    8')
         WRITE(25,81) (ZMINA(I),ZMAXA(I),ZSTEP(I),I=1,3) ! angle search info
81          FORMAT('SEEK 1',9F6.1)
      ELSE
         WRITE (25,148)
 148        FORMAT ('VOLS    19    8'/'SEEK') ! special for molec on 2-fold axis
      ENDIF
      WRITE (25,91) POTENTIAL_FILE
91       FORMAT('WMIN',11X,A/'FINI')
      WRITE (*,1111) MASTER_NAME
1111     FORMAT (' >>>>>>Name of the new COM file: ',41A1)
      CLOSE (UNIT=25)
!
!----Finished writing MOLPAK input file on unit # 25
!
!----Data for searches on unit # 25
      OPEN (UNIT=25, FILE=SEARCH_FILE, STATUS='UNKNOWN')
      WRITE(25,231) G92FILE2(1:NID+7), GROUP, TITLE
231      FORMAT (A/A2/A73)
      WRITE (25,'(I5,3F5.1,I5)') NMOL, ZSTEP, NSOL ! and soln to be elaborated
      CLOSE (UNIT=25)
!----Section for minimum volume search
      WRITE (20,442)
442   FORMAT('#'/'# ****section for minimum volume search****'/'#'/ &
!            'echo "....distances before wmin...." > wmin.dis'/ &
!            'echo "....distances after wmin...." > wmin-lsq.dis'/ &
             'mv -f fort.8 fort.10  # min vol file on unit # 10')
      WRITE (20,448) SEARCH_FILE
448   FORMAT ('rm -f fort.15'/'cp -f ',A14,'  fort.15'/ &
     &        '$dir/UTILITIES/min-vol-search.exe', &
     &        '   # make unique minimum volume file')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
493                FORMAT (' cat fort.16 >> ',46A1)
      WRITE (20,551) MOLPAK_FILE, SEARCH_FILE
551   FORMAT ('#'/'# ****create MOLPAK command files****'/'#'/ &
     &        'cp -f fort.11 volume.min '/'mkcom:'/ &
     &        'mv -f fort.11 fort.8  # min vol file now on unit # 8'/ &
     &        'rm -f fort.15 fort.10'/ &
     &        'cp -f ',A15,'  fort.10'/ &
     &        'cp -f ',A14,'  fort.15  # unit 15 contains angle search data'/ &
     &        '#          and number of initial volumes to elaborate'/ &
     &        '$dir/UTILITIES/make-molpak-com-find.exe')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,552)
552      FORMAT ('#'/'#****create files for MOLPAK/WMIN/NBSLATTICE****'/'#')
!
!----FINISHED
!
!***************** LOOP THROUGH NMOL RUN******************************
!
!----Assign files for running MOLPAK
      WRITE (20,600) NMOL+1
600      FORMAT ('runwmin:'/'set i=1'/'while ($i != ',I4,')')
      WRITE (20,456)
456      FORMAT('if ( -e input.$i ) then'/'rm -f fort.15'/ &
     &          'cp -f input.$i fort.15    # MOLPAK input file')
      WRITE (20,458)
458   FORMAT ('#...xyz file on unit # 9'/ &
     &        '#...sum file on unit # 21'/ &
     &        'echo "------------------------------------------------------------"'/ &
     &        'echo "solution # $i .... "'/ &
     &        '$dir/MOLPAK/molpak_12.exe         # run MOLPAK...')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,1951)
1951     FORMAT (' if ($i == 1) then')
      WRITE (20,1952) G92FILE2(1:NID+7), SUM, G92FILE2(1:NID+7), SUM
1952     FORMAT('   echo "solution # $i .... " > ',A,6A1/ &
     &'   cat fort.21 >> ',A,6A1/'   rm -f fort.21')
      WRITE (20,952) G92FILE2(1:NID+7), SUM, G92FILE2(1:NID+7), SUM
952      FORMAT(' else'/'   echo "solution # $i .... " >> ',A,6A1/ &
     &'   cat fort.21 >> ',A,6A1/ &
     &       '   rm -f fort.21')
      WRITE(20,955)
955      FORMAT (' endif'/'rm -f fort.15'/ &
     & ' if ($i == 1) then'/'  cp -f fort.49 molpak.cell'/' else'/ &
     &           '  cat fort.49 >> molpak.cell'/' endif'/ &
     &           'cp -f fort.9 fort.15   # input file for WMIN')
!    & 'cp fort.15 WMININP_MPK/wmininput.$i')
      IF(NSOL .NE. 0) WRITE (20,555)
555                   FORMAT ('cp -f fort.9 wmin.input')
      WRITE (20,957)
957      FORMAT ('#Files for WMIN .....')
!      NECS = NEC(3)
      WRITE  (20,6331)
6331     FORMAT ('#....sum file on unit # 21'/ &
     &           '#...save file on unit # 22'/ &
     &           '#...WMIN cell on unit # 11'/ &
     &           'if ($i == 1) cp fort.15 wmininput.$i'/ &
!    &           'echo " run WMIN cycle 1 to 8; NCYS = final cycle "'/ &
     &'(limit cputime 600 ; $dir/WMIN/wmin_100.exe)     #  run WMIN')  ! replaced wmin.exe 11-4-04
      WRITE (20,1665)
1665  FORMAT('if (-e core) rm -f core' / &
     & 'rm -f fort.16')
!    & ' then'/' rm -f core'/'else'/ &
!    & ' cp fort.22 save.$i'/' mv fort.16 out.$i'/ &
!    & ' mv fort.29 dis.$i' /' mv fort.99 list.$i'/'endif')
!    &'echo " ### soln $i ### " >> wmin.dis'/&
!    &'cat fort.29 >> wmin.dis')
      WRITE (20,500) G92FILE2(1:NID+7), SUM , G92FILE2(1:NID+7), G92FILE2(1:NID+7), &
     &               G92FILE2(1:NID+7)                      ! 9-4-14
!     chenaged them on 11-4-04
      WRITE (20,1666)
 1666 FORMAT('#make new wmin input for LSQ'/ &
     & 'if ($i == 1) cp -f fort.22 wmin_RSS.save' /&
     & 'mv -f fort.22 wmin.save'/ &
     & '$dir/UTILITIES/rewrite-wmininp-100.exe'/ &
     & 'if (-e wminlsq.inp) then'/ &
     & 'mv -f  wminlsq.inp fort.15'/ &
!    & 'echo " refinement with least squares"'/ &
     & '#refinement with least squares'/ &
     & '#cp fort.15 WMININP_LS/wmininput.$i'/ &
     & '(limit cputime 600 ; $dir/WMIN/wmin_100.exe)'/ &
     & 'if ($i == 1) cp -f fort.22 wmin_LSQ.save'/ &
!    & 'cat  fort.22 >> wmin.save'/&                     ! 9-24-07
     & 'if (-e core) rm -f core'/'endif'/ &
     & 'if(-e wminnew.inp) rm -f wminnew.inp' )
!    & ' then'/' rm -f core'/'else'/ &
!    & ' cp fort.22 save_lsq.$i'/' mv fort.16 out_lsq.$i'/ &
!    & ' mv fort.29 dis_lsq.$i' /' mv fort.99 list_lsq.$i'/'endif')
!     WRITE (20,1666)
!1666 FORMAT(' set n=1'/' while ($n != 20)'/ &
!    &          ' if ( -e fort.22 ) then'/ &
!    &          '  mv -f fort.22 wmin.save'/ &
!    &          '  $dir/UTILITIES/rewrite-wmininp.exe'/ &
!    &          '  if ( -e wminnew.inp ) then'/ &
!    &          '    @ cycle1 = $n * 8'/'    @ cycle2 = $cycle1 + 1'/ &
!    &'    @ cycle3 = $cycle2 + 7'/'    mv -f wminnew.inp fort.15'/ &
!    &'    echo " continue running WMIN, cycle $cycle2 to $cycle3;', &
!    &' NCYS = final cycle "'/'    $dir/WMIN/wmin.exe'/'  else'/ &
!    &'   if ( -e wminlsq.inp ) then'/ &
!    &'    mv -f wminlsq.inp fort.15 '/ &
!    &'    echo " refinement with least squares"'/&
!    &'    $dir/WMIN/wmin.exe   # refinement with least squares'/&
!    &       '    echo " ### soln $i ### " >> wmin-lsq.dis'/&
!    &       '    cat fort.29 >> wmin-lsq.dis'/&
!    &       '    mv -f fort.22 wmin.save'/ &
!    &       '    goto wminend'/'   else'/ &
!    &       '    goto wminend'/'   endif'/'  endif'/' endif'/ &
!    &       '   @ n = ($n + 1)'/' end'/'wminend:'/&
!    &       'cp fort.15 WMININP_LS/wmininput.$i')
      WRITE(20,500) G92FILE2(1:NID+7), SUM, G92FILE2(1:NID+7), G92FILE2(1:NID+7), &
     &              G92FILE2(1:NID+7)                      ! 9-4-14
500   FORMAT ('if(-e fort.21) cat fort.21 >> ',A,6A1 / &
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >> ', A, '.sum'/'endif')              ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
      WRITE (20,1955)
1955     FORMAT ('if ( -e fort.11 ) then')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
!     WRITE(20,399) G92FILE2(1:NID+7), SAVE
399      FORMAT (' if ($i == 1) then'/ &
     &           '   mv -f wmin.save ',A,5A1/ &
     &           '   cp -f fort.15 wmin.inp')
!     WRITE(20,395) G92FILE2(1:NID+7), SAVE
395      FORMAT (' else'/'   cat wmin.save >> ',A,5A1/' endif')
      WRITE (20,398)
398      FORMAT ('# ..... cell reduction # .....')
      WRITE (20,361)
361   FORMAT (' rm -f fort.15'/ &
     &        ' cp -f fort.11 fort.15  # establish WMIN cell as input'/ &
     &        ' $dir/UTILITIES/nbslattice.exe  # cell reduction')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE(20,501) G92FILE2(1:NID+7), SUM
501      FORMAT (' cat fort.11 >> ',A,6A1/' rm -f fort.11')
      WRITE (20,1956)
1956  FORMAT ('else '/ '  if ($i == 1) then'/'   cp -f fort.15 wmin.inp'/'  else'/ &
     &        '   cat fort.15 >> wmin.inp'/'  endif'/ &
     &        ' echo " failure in WMIN, No = $i " '/'endif')
      WRITE (20,889)
889   FORMAT ('endif'/' rm -f input.$i wmin.save fort.*'/' @ i = ($i + 1)'/'end')
!----Finished producing the IMOL fine search packages
      WRITE (20,708)
708      FORMAT ('#'/'#****create summary tables****'/'#')
      WRITE (20,709) G92FILE2(1:NID+7)                      ! 5-4-10
709      FORMAT ('rm -f fort.*'/'cp -f ',A,'.sum fort.13')  ! 5-4-10
      WRITE(LINE_OUT,710) G92FILE2(1:NID+7), SUM
710      FORMAT ('cp -f ',A,6A1)
      LINE_OUT2 = FORT15
      WRITE (20,337) LINE_OUT
      WRITE (20,712) SEARCH_FILE
712      FORMAT ('cp -f ',A14,'  fort.14'/&
     &           '$dir/UTILITIES/table-1.exe   # make final tables'/&
     &           '  cp -f fort.16 temp.tab')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE(20,711) (SPGR(I),I=1,2),TAB
 711     FORMAT ('  mv -f fort.16 ',2A1,4A1/'  mv -f temp.tab fort.15'/&
     &           '$dir/UTILITIES/table-2.exe   # resort table'/&
     &           '  mv -f fort.16 tab.den')
!----Finish up by deleting various files
      WRITE (20,716)
716   FORMAT ('#****clean up files****')
      IF (.NOT. I_SAVE) WRITE (20,718)
718      FORMAT ('  rm -f *.save')
719   WRITE (20,720)
720      FORMAT ('  rm -f fort.* '/'  rm -f wmin.inp *.wmin'/ &
!    &'gzip wmin.dis'/'gzip wmin-lsq.dis'/&
!    &'cd WMININP_MPK'/&
!    &'tar -cvf WMININP.tar wmininput.*'/&
!    &'gzip WMININP.tar'/'rm wmininput.*'/&
     &'#cd WMININP_LS'/&
     &'#tar -cvf WMININP.tar wmininput.*'/&
     &'#gzip WMININP.tar'/'#rm wmininput.*'/'#cd ..'/&
     &           'time'/'date'/'#'/'#**** FINISHED****')
!
      CLOSE (UNIT = 15)
      CLOSE (UNIT = 20)
      CLOSE (UNIT = 25)
      ENDDO Big_loop                 ! finished
!
!----If N_GO > 2 , then CALL WRUNJUB to write a command file
!     to make new directories, move various files to the new
!     directories and submit one job for each coord geom.
3100  IF (N_GO .GT. 1) CALL WRUNJOB
      STOP
      END PROGRAM MAKE_FILES

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE WRUNJOB
!
      USE make_filesCommonMod, only: COMPD_ID, COMPDID, IDCHG, I, J,&
     &                               KIND, NID, POTENTIAL_FILE

      IMPLICIT NONE

      INTEGER, PARAMETER ::  N_SPGR = 54
      CHARACTER (2) ::  DIR(N_SPGR) = (/'aa','db','ah','ab','da','af',&
     &          'bh','ba','bb','av','bf','fc','dc','ai','bd','ap','aq',&
     &               'ay','fa','az','ak','am','de','dd','ca','ce','cc',&
     &               'cd','cb','au','as',&
     &             'ae', 'ac', 'ad', 'ag', 'aj', 'al', 'fd', 'an', 'ao',&
     &             'fb', 'ar', 'at', 'be', 'aw', 'bg', 'ax', 'bi', 'bc',&
     &             'bj', 'bk', 'cf', 'df', 'dg'/)
      CHARACTER (2) :: FILE(N_SPGR) = (/'AA','DB','AH','AB','DA','AF',&
     &          'BH','BA','BB','AV','BF','FC','DC','AI','BD','AP','AQ',&
     &               'AY','FA','AZ','AK','AM','DE','DD','CA','CE','CC',&
     &               'CD','CB','AU','AS', &                                ! default 31
     &              'AE', 'AC', 'AD', 'AG', 'AJ', 'AL', 'FD', 'AN', 'AO',&
     &              'FB', 'AR', 'AT', 'BE', 'AW', 'BG', 'AX', 'BI', 'BC',&
     &              'BJ', 'BK', 'CF', 'DF', 'DG'/)                         ! additional 23

!
      OPEN (UNIT=30, FILE='runjobs', STATUS='UNKNOWN')
!----Make new directories & move *.com to new directories
      WRITE (30,10) (FILE(I),I=1,29), (FILE(I),I=30,N_SPGR), &
     &              (DIR(I), I=1,29), (DIR(I), I=30,N_SPGR), &
     &              N_SPGR+1,(COMPDID(1:NID),J=1,2), &
     &              POTENTIAL_FILE                        ! 26feb03
10        FORMAT ('#!/bin/csh'/ &
     &   'setenv dir /c/Users/ACER/Downloads/PREDICTIONS'/ &
     &            '# write new directories & move *.com file'/ &
     &            'set n=1'/'set f=',1H(,29(A2,1X),1H\,/25(A2,1X),1H)/ &
     &            'set d=',1H(,29(A2,1X),1H\,/25(A2,1X),1H)/ &
     &            'set dir = ($d)'/ &
     &            'set fl  = ($f)'/'while ($n !=',I3,')'/ &
     &            ' if (-e ',A,'_$fl[$n]_master.com) then'/ &
     &            '  if ( -d $dir[$n] ) then'/ &
     &            '   goto skip'/'  else'/'   mkdir $dir[$n]'/ &
     &            '  endif'/'skip:'/ &
     &            '   if (-e volume.min) mv -f volume.min $dir[$n]'/ &
     &            '   mv -f ',A,'_$fl[$n]_master.com $dir[$n]'/ &
     &            '   mv -f ','molpak_$fl[$n].input $dir[$n]'/ &
     &            '   mv -f ','search_$fl[$n].data $dir[$n]'/ &
     &            '   cp -f ',A,' $dir[$n]/.'/ &
     &            ' endif'/'@ n = ($n + 1)'/'end'/ &                     ! 26feb03
     &            'echo " directories created & files moved "')
!----Submmit one job at a time
      WRITE (30,15) N_SPGR+1, COMPDID(1:NID)
15       FORMAT ('#****one job at a time****'/ &
     &       'date'/'time'/'set n=1'/'while ($n !=',I3,')'/ &
     &       ' if ( -d $dir[$n] ) then'/'   cd $dir[$n]'/ &
     &       '   echo "current job number: $n " '/ &
     &       '   echo "current directory : $dir[$n]" '/ &
     &       '   if ( -e ',A,'_$fl[$n]_master.com ) then'/ &
     &       '    if (-e master.log) then'/ &
     &       '     echo "the code $fl[$n] is already running" '/ &
     &       '    else'/ &
     &       '     ./*.com > master.log '/'    endif'/ &
     &       '    date'/'    time'/'   endif'/'   cd ..'/' endif'/ &
     &       '@ n = ($n + 1)'/'end'/'exit')
!----Run several summary files
        WRITE (30,20)
20         FORMAT ('#create a file named densities.list and sorted files'/&
     &             '#   named densities and energies'/&
     &             '$dir1/UTILITIES/resort-summarize.com'/&
     &             '#****FINISHED****')
      RETURN
      END SUBROUTINE WRUNJOB

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE REF_PMIN                     ! 7-9-08
!
      USE make_filesCommonMod
!
      IMPLICIT NONE
!
!----Set initial values for symmetry flags
      INVER2 = 0
      IMR    = 0
      IAX    = 0
      IDCHG  = 1
!----How many characters in compound ID?
      NID = 0
      DO I=1,8
         IF(COMPD_ID(I) .EQ. BLANK) GO TO 1332
         NID = NID + 1
         MASTER_NAME(I) = COMPD_ID(I)
      ENDDO
!
1332  IALL = 0
!     IF (GROUP .EQ. '  ' .OR. GROUP .EQ. 'AL' .OR. GROUP .EQ. 'al') &
      IF (GROUP .EQ. '  ' .OR. GROUP .EQ. '31' .OR. GROUP .EQ. '54') &
     &      IALL = 1
!----Set IDCHG = 2...hold-over from original program.  Only one type
!     of charge is input with atom coordinates.
      IDCHG = 2
      POT_FILE = ' '
      WRITE (*,155)
155      FORMAT ('The potential parameter file is...'/&
     &           '  1 = UMD_potential, 2 = ARL_potential,',&
     &           ' 3 = other [1]: ',$)
      READ (*,337) POT_FILE
         SELECT CASE (POT_FILE)
            CASE (' ','1')
               POTENTIAL_FILE = "UMD_potential"    ! 7-9-08
            CASE ("2")
               POTENTIAL_FILE = "ARL_potential"
            CASE ("3")
               WRITE (*,160)
160            FORMAT ('Provide name of potential parameter file ',&
     &                 '(40 chars max) :',$)
               READ (*,337) POTENTIAL_FILE
         END SELECT
!----Should information be read for pmin cross-term interactions?      6/24/09
      write (*,156)                                                  ! 6/24/09
156      format ('Add PMIN cross-term coefficient data?...'/ &       ! 6/24/09
            '   0/blank = none, 1-5 = number of interactions: ',$)   ! 6/24/09
      read (*,337) cross                                             ! 6/24/09
         select case (cross)                                         ! 6/24/09
            case (' ','0')                                           ! 6/24/09
               go to 157
            case ('1':'5')                                           ! 6/24/09
               read (cross,'(i1)') i_cross                           ! 6/24/09
               write (*,259) i_cross                                 ! 6/24/09
259               format ('Read',i2,' set(s) of cross-term coefficients...') ! 6/24/09
               go to 159                                             ! 6/24/09
         end select                                                  ! 6/24/09
!----Read i_cross sets of cross-term potentials...1 per line         ! 6/24/09
159      do i=1,i_cross                                              ! 6/24/09
            read (*,337) cross_terms(i)                              ! 6/24/09
         enddo                                                       ! 6/24/09
!
157   I_COMP = .false.
      I_SAVE = .false.
!----Read 6th line...name of MOLPAK coordinate input file.  Try to
!     pick up file name automatically
      OPEN (UNIT=14, FILE='MOLPAK.NAME', STATUS='OLD')  ! look for files
      I = 0
45    READ (14,337,END=47) MOLPAK_NAME                   ! that have
337   FORMAT (A)
      I = INDEX(MOLPAK_NAME,'molpak.xyz')               ! molpak.xyz
      IF (I .NE. 0) THEN
         L = INDEX(MOLPAK_NAME,'.xyz ') + 3
         GO TO 48
      ENDIF
      GO TO 45                                          ! at the end
47    IF (I .EQ. 0) THEN
         MOLPAK_NAME = 'molpak.xyz'
         L = 10
      ENDIF
48    WRITE (*,50) MOLPAK_NAME(1:L)
50    FORMAT('Name of the MOLPAK atom coord file', &
     &       ' [',A,']: ',$)
      READ (*,337) WHAT_LINE
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,337) ATOM_LIST
51       FORMAT (30A1)
      ELSE
         ATOM1_LIST = MOLPAK_NAME
      ENDIF
!----Inquire about possible 2-fold symmetry in model
      WRITE (*,486)
 486  FORMAT ('The following questions are about molecules with', &
     &        ' C2, Cs and Ci symmetry'/ &
     &       '   and will restrict searches to space groups with', &
     &       ' appropriate symmetries;'/ &
     &       '   reply N/RETURN for a general search that ignores', &
     &       ' molecular symmetry...'/ &
     &       '    Does the molecule contain two-fold symmetry?', &
     &       ' (Y/N [N]): ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. ' ' .OR. Y_N .EQ. 'N' .OR. Y_N .EQ. 'n')  THEN
         WHAT_LINE(1:1) = ' '
      ELSE
         WRITE (*,4861)
 4861    FORMAT('    Which axis is parallel to two-fold axis?  '/ &
     &          '      1 = X, 2 = Y, 3 = Z: ',$)
         READ (*,337) WHAT_LINE
      END IF
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,3303) IAX
3303        FORMAT (I1)
         IF (IAX .GE. 1 .AND. IAX .LE. 3) THEN
            ELMSPGR(1:40) = 'AAABCAAHAFAIAKAMFAFCDA&
                           &AQAZAYBHAVBDBFCCCB'
         ELSE
            IAX = 0
         ENDIF
      ENDIF
      IF (IAX .NE. 0) GO TO 505 ! if axis, forget mirror & inversion question
!----Inquire about possible mirror plane symmetry in model
      WRITE (*,4862)
4862     FORMAT('    Does the molecule contain mirror symmetry? (Y/N', &
     &       ' [N]): ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. ' ' .OR. Y_N .EQ. 'N' .OR. Y_N .EQ. 'n')  THEN
        WHAT_LINE(1:1) = ' '
      ELSE
         WRITE (*,4863)
 4863    FORMAT('    Which axis is perpendicular to the mirror plane?  '/ &
     &          '      1 = X, 2 = Y, 3 = Z: ',$)
         READ (*,337) WHAT_LINE(1:1)
      END IF
      IMR = 0
      IF (WHAT_LINE .NE. INPUT_LINE)  READ (WHAT_LINE,3303) IMR
      IF (IMR .NE. 0) GO TO 505  ! if mirror, forget inversion question
!----Inquire about possible inversion symmetry in model
      WRITE (*,496)
496      FORMAT('    Does the molecule contain a center of symmetry', &
     &          ' at the origin? (Y/N [N]): ',$)
      READ (5,337) INVER1
      IF (INVER1 .EQ. ' ' .OR. INVER1 .EQ. 'N' .OR. &
     &    INVER1 .EQ. 'n') THEN
          INVER2 = 0
      ELSE
         ELMSPGR(1:30) = 'AAAHAFDADBAPBABBAQAZAYBHAVBDBF'
         INVER2 = 1
      ENDIF
!----Open atom coordinate file
505   OPEN (UNIT=15, FILE=ATOM1_LIST, STATUS='OLD')
      WRITE (*,109)
109      FORMAT (4X,'atom id     x         y         z      Id    ', &
     &        '6-31g*')
      N = 0
111   N = N + 1
        READ (15,1,END=100) (ATOM_TY(N,J),J=1,2), (ATOM(N,I),I=1,4), &
     &                    X(N), Y(N), Z(N), IDATOM(N), G92CHARGE(N)
1          FORMAT (5X,2A1,4A1,3F10.6,I5,F10.6)
        WRITE (*,112) (ATOM_TY(N,J),J=1,2), (ATOM(N,I),I=1,4), &
     &                X(N), Y(N), Z(N), IDATOM(N), G92CHARGE(N)
112        FORMAT(4X,2A1,4A1,3F10.6,I5,F10.6)
      GO TO 111
100   N = N -1
!
!----Question about using an existing solution
      NSOL = 0             ! solution #
      Y_N = ' '
      WRITE (*,10)
10       FORMAT ('Do you want to rerun a solution from an existing', &
     &           ' volume.min file [N]: ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. 'Y' .OR. Y_N .EQ. 'y') THEN
          WRITE (*,12)
12           FORMAT (' Enter solution number: ',$)
          READ (*,*) NSOL
          WHAT_LINE = '1'
          OPEN (UNIT=30, FILE='volume.min', STATUS="OLD", ERR=800)
          GO TO 801
800       WRITE (*,804)
804          FORMAT (' The required file volume.min is not exist!')
          STOP
801       READ (30,*)
          READ (30,'(9X,3F8.1)') (ZSTEP(I), I=1,3)
          DO L=1,3
             ZMINA(L) = -90.
             ZMAXA(L) =  90.
          END DO
      ELSE
!----Read 7th line...angle search limits for initial molpak search
          WRITE (*,80)
80            FORMAT('Minimum, maximum and step for 3 angles for 1st ', &
     &               'MOLPAK map'/ &
     &               '  [-90,90,10,-90,90,10,-90,90,10,]: ',$)
          READ (*,337) WHAT_LINE
          IF (WHAT_LINE .NE. INPUT_LINE) THEN
              READ (WHAT_LINE,*) (ZMINA(I),ZMAXA(I),ZSTEP(I),I=1,3)
          ELSE
!----Set up defaults...
              DO L=1,3
                 ZMINA(L) = -90.
                 ZMAXA(L) =  90.
                 ZSTEP(L) =  10.
              ENDDO
          ENDIF
!----Read number of minimum volume points to be saved in MOLPAK search
             WRITE (*,703)
703             FORMAT ('Number of minimum volumes from initial MOLPAK', &
     &                  ' search to pass to next step,'/  &
     &                  '   suggest 7000 for 10 deg step search and', &
     & ' 60000 for 5 deg step search or range(-180 to 180) [7000]: ',$)
             WHAT_LINE = INPUT_LINE
             READ (*,337) WHAT_LINE
             IF (WHAT_LINE .EQ. INPUT_LINE) THEN
                N_VOLS_SAVED = 7000
             ELSE
                READ (WHAT_LINE,*) N_VOLS_SAVED
             ENDIF
!----Read number of minimum volume points to be elaborated
             WRITE (*,803)
803             FORMAT('Maximum number of 2nd MOLPAK calcns desired [500]: ',$)
             WHAT_LINE = INPUT_LINE
             READ (*,337) WHAT_LINE
      END IF
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,*) NMOL
      ELSE
         NMOL = 500       ! default -> elaborate 500 solns
      ENDIF
!----How many coord geoms have been specified?  Each 2 letter code
!     should be separated from its neighbors by a space or comma.
      IF (IALL .NE. 0) GO TO 29
      N_SPGRPS = 0
      NSCHARS = 0
      DO 22 I=1,40
         IF (SPGR(I) .EQ. ' ' .OR. SPGR(I) .EQ. ',') GO TO 22
         IF (NSCHARS .NE. 0) THEN
            IF (NSCHARS .EQ. 1) GO TO 24
            NSCHARS = 0
         ENDIF
23       N_SPGRPS = N_SPGRPS + 1
24       NSCHARS = NSCHARS + 1
         MOL_SPGRP(NSCHARS,N_SPGRPS) = SPGR(I)
22    CONTINUE
!----All or just some space groups
29    IF (IALL .EQ. 0) THEN
         N_GO = N_SPGRPS
      ELSE
       IF (GROUP .EQ. '  ' .OR. GROUP .EQ. '31') THEN    ! 9-18-07
         N_GO = N_SPGR-23                                ! 9-18-07
       ELSE
         N_GO = N_SPGR
       END IF
      ENDIF
!
!----Big LOOP for creating the N_GO sets of files for execution
!
Big_loop : DO LGRPS=1,N_GO
               IF (IALL .EQ. 0) THEN
                  SPGR(1) = MOL_SPGRP(1,LGRPS)
                  SPGR(2) = MOL_SPGRP(2,LGRPS)
               ELSE
                  DO J=1,N_GO
                     IF (ALL_SPGR(LGRPS) .EQ. ELM_SPGR(J)) CYCLE Big_loop
                  END DO
                  SPGR(1) = ALL_SPGR(LGRPS)(1:1)
                  SPGR(2) = ALL_SPGR(LGRPS)(2:2)
               ENDIF
!----Go through the N_SPGRPS coord geoms & make appropriate files.
!     Tack 3 characters (_SPGR) onto end of COMPD_ID and MASTER_NAME.
      COMPD_ID(NID+1) = '_'
      COMPD_ID(NID+2) = OLD(1)
      COMPD_ID(NID+3) = OLD(2)
      COMPD_ID(NID+4) = OLD(3)
      COMPD_ID(NID+5) = '_'
      COMPD_ID(NID+6) = SPGR(1)
      COMPD_ID(NID+7) = SPGR(2)
      MASTER_NAME(NID+1) = '_'
      MASTER_NAME(NID+2) = SPGR(1)
      MASTER_NAME(NID+3) = SPGR(2)
!----Construct name for G98 file
      DO I=1,NID
         G92FILE(I)=COMPD_ID(I)
      END DO
      G92FILE(NID+1) = '_'
      G92FILE(NID+2) = G92(1)
      G92FILE(NID+3) = G92(2)
      G92FILE(NID+4) = G92(3)
      G92FILE(NID+5) = '_'
      G92FILE(NID+6) = SPGR(1)
      G92FILE(NID+7) = SPGR(2)
!----Construct name for MASTER.COM file that will have compound ID 1st,
!     space group symbol 2nd and "master.com" 3rd.
      MID = NID + 3
      DO I=1,11
         MASTER_NAME(MID+I) = MASTERCOM(I)
      ENDDO
      MOLPAK_FILE(1:7) = 'molpak_'      ! create name of the
      MOLPAK_FILE(8:9) = GROUP          ! molpak input file
      MOLPAK_FILE(10:15) = '.input'
      SEARCH_FILE(1:7) = 'search_'      ! create name for search
      SEARCH_FILE(8:9) = GROUP          ! data file
      SEARCH_FILE(10:14) = '.data'
!----Open files
      OPEN (UNIT=20, FILE=MASTER_NAME2, STATUS='UNKNOWN')
      OPEN (UNIT=25, FILE=MOLPAK_FILE,  STATUS='UNKNOWN')
!----Start to write master.com... files
      WRITE (20,39) (NMOL+1)
39    FORMAT ('#!/bin/csh'/'pwd'/'#'/'time'/'date'/'#'/ &
     &   'setenv dir /c/Users/ACER/Downloads/PREDICTIONS'/&
     &   '#'/'#****initial MOLPAK search****'/'#'/ &
     &   '# if input.* exists, restart PMIN'/ &
     &   'set k=1'/'while ($k != ',I4,')'/ &
     &   '  if ( -e input.$k ) then'/'    goto runpmin'/'  endif'/ &
     &   '  @ k = ($k + 1)'/'end'/ &
     &   '# if volume.min  is existing, skip first MOLPAK search'/ &
     &   'if ( -e volume.min  ) then'/'   cp -f volume.min  fort.11'/ &
     &   '   goto mkcom'/'else'/'   if ( -e volume.find ) goto mkcom'/ &
     &   'endif'/ &
     &   '# if molpak.vol exists, skip MOLPAK search'/&
     &   'if ( -e molpak.vol ) cp -f molpak.vol molpak.min'/'#')
      WRITE (20,41) MOLPAK_FILE
41    FORMAT ('cp -f ',A15,'  fort.15       # MOLPAK input file'/ &
     &       '$dir/MOLPAK/molpak_12.exe',' # run MOLPAK'/ &
     &       'if ( -e molpak.min ) rm -f molpak.min'/ &
     &       'rm -f molpak.vol'/ &
     &       'time'/'date')
      NDIF = 30 - MID
!----Write title on MOLPAK input file + rest of it
      WRITE (25,245) COMPDID(1:NID), TITLE
245   FORMAT ('HEAD  ',A,8X,A)
!     using default  2-27-04
!     WRITE (25,146)                      ! line # 2...NSEG values
146   FORMAT ('NSEG  128 128   8')
!----Place atom coordinates into molpak input file
      DO L=1,N
         WRITE(25,483) (ATOM_TY(L,J),J=1,2), (ATOM(L,I),I=1,4), &
     &                 X(L), Y(L), Z(L), IDATOM(L), G92CHARGE(L)
483         FORMAT('ATOM ',2A1,4A1,3F10.6,I5,F10.6)
      ENDDO
      IF (INVER2 .EQ. 0) GO TO 7631 ! INVER2 = 0  means no center
      WRITE (25,764)                ! write "CENT" line into
764       FORMAT ('CENT',30X,' ')   ! MOLPAK input file
      GO TO 488
 7631 IF (IMR .EQ. 0) GO TO 763    ! IMR = 0 means no mirror
      WRITE (25,7632) IMR          ! write "PLAN" line into
 7632    FORMAT ('PLAN',I2,28X,' ') ! MOLPAK input
      GO TO 488
763   IF (IAX .EQ. 0) GO TO 488    ! IAX = 0 means no 2-fold axis
      WRITE (25,787) IAX           ! write "AXIS" line into
787      FORMAT ('AXIS',I2,28X,' ')! MOLPAK input file
488   IF (I_COMP)  WRITE (20,491) (COMPD_ID(I),I=1,MID), COMPLT
491                FORMAT ('mv -f fort.16 ',46A1)
      WRITE(25,246) GROUP
246   FORMAT('INCL ',A2)
      IF (IAX .EQ. 0) THEN
         WRITE (25,147) N_VOLS_SAVED  ! # min vols for 1st MOLPAK
147         FORMAT ('VOLS ',I5,'    8')
         WRITE(25,81) (ZMINA(I),ZMAXA(I),ZSTEP(I),I=1,3) ! angle search info
81          FORMAT('SEEK 1',9F6.1)
      ELSE
         WRITE (25,148)
 148        FORMAT ('VOLS    19    8'/'SEEK') ! special for molec on 2-fold axis
      ENDIF
      WRITE (25,91) POTENTIAL_FILE, i_cross                      ! 6/24/09
91       format ('PMIN',11x,a40,i4)                              ! 6/24/09
      if (i_cross .ge. 1) then                                   ! 6/24/09
         write (25,93) (cross_terms(i),i=1,i_cross)              ! 6/24/09
93          format ('CROS'6x,a40)                                ! 6/24/09
      endif                                                      ! 6/24/09
      write (25,94)                                              ! 6/24/09
94       format ('FINI')                                         ! 6/24/09
      WRITE (*,1111) MASTER_NAME
1111     FORMAT (' >>>>>>Name of the new COM file: ',41A1)
      CLOSE (UNIT=25)
!
!----Finished writing MOLPAK input file on unit # 25
!
!----Data for searches on unit # 25
      OPEN (UNIT=25, FILE=SEARCH_FILE, STATUS='UNKNOWN')
      WRITE(25,231) G92FILE2(1:NID+7), GROUP, TITLE
231      FORMAT (A/A2/A73)
      WRITE (25,'(I5,3F5.1,I5)') NMOL, ZSTEP, NSOL ! and soln to be elaborated
      CLOSE (UNIT=25)
!----Section for minimum volume search
      WRITE (20,442)
442   FORMAT('#'/'# ****section for minimum volume search****'/'#'/ &
!            'echo "....distances before pmin...." > pmin.dis'/ &
!            'echo "....distances after pmin...." > pmin-lsq.dis'/ &
             'mv -f fort.8 fort.10  # min vol file on unit # 10')
      WRITE (20,448) SEARCH_FILE
448   FORMAT ('rm -f fort.15'/'cp -f ',A14,'  fort.15'/ &
     &        '$dir/UTILITIES/min-vol-search.exe', &
     &        '   # make unique minimum volume file')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
493                FORMAT (' cat fort.16 >> ',46A1)
      WRITE (20,551) MOLPAK_FILE, SEARCH_FILE
551   FORMAT ('#'/'# ****create MOLPAK command files****'/'#'/ &
     &        'cp -f fort.11 volume.min '/'mkcom:'/ &
     &        'mv -f fort.11 fort.8  # min vol file now on unit # 8'/ &
     &        'rm -f fort.15 fort.10'/ &
     &        'cp -f ',A15,'  fort.10'/ &
     &        'cp -f ',A14,'  fort.15  # unit 15 contains angle search data'/ &
     &        '#          and number of initial volumes to elaborate'/ &
     &        '$dir/UTILITIES/make-molpak-com-find.exe')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,552)
552      FORMAT ('#'/'#****create files for MOLPAK/PMIN/NBSLATTICE****'/'#')
!
!----FINISHED
!
!***************** LOOP THROUGH NMOL RUN******************************
!
!----Assign files for running MOLPAK
!     IF (N_GO .LE. 1) THEN                                                          ! 6-9-08 DU
!     WRITE(20,601)                                                                  ! 6-9-08 DU
!601  FORMAT('cp -f $dir/PMIN/PMIN.inp  .  # if need, to modify PMIN.inp')          ! 6-9-08 DU
!     END IF                                                                         ! 6-9-08 DU
      WRITE (20,600) NMOL+1
600      FORMAT ('runpmin:'/'set i=1'/'while ($i != ',I4,')')
      WRITE (20,456)
456      FORMAT('if ( -e input.$i ) then'/'rm -f fort.15'/ &
     &          'cp -f input.$i fort.15    # MOLPAK input file')
      WRITE (20,458)
458   FORMAT ('#...xyz file on unit # 9'/ &
     &        '#...sum file on unit # 21'/ &
     &        'echo "------------------------------------------------------------"'/ &
     &        'echo "solution # $i .... "'/ &
     &        '$dir/MOLPAK/molpak_12.exe         # run MOLPAK...')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,1951)
1951     FORMAT (' if ($i == 1) then')
      WRITE (20,1952) G92FILE2(1:NID+7), SUM, G92FILE2(1:NID+7), SUM
1952     FORMAT('   echo "solution # $i .... " > ',A,6A1/ &
     &'   cat fort.21 >> ',A,6A1/'   rm -f fort.21')
      WRITE (20,952) G92FILE2(1:NID+7), SUM, G92FILE2(1:NID+7), SUM
952      FORMAT(' else'/'   echo "solution # $i .... " >> ',A,6A1/ &
     &'   cat fort.21 >> ',A,6A1/ &
     &       '   rm -f fort.21')
      WRITE(20,955)
955      FORMAT (' endif'/'rm -f fort.15'/ &
     & ' if ($i == 1) then'/'  cp -f fort.49 molpak.cell'/' else'/ &
     &           '  cat fort.49 >> molpak.cell'/' endif'/ &
     &           'cp -f fort.9 fort.15   ')
!     IF(NSOL .NE. 0) WRITE (20,555)
555                   FORMAT ('cp -f fort.9 wmin.input')
!      NECS = NEC(3)
      WRITE  (20,6331)
6331     FORMAT ('######### RUN PMIN #########'/ &
!    &'if (-e wmin.input) rm -f wmin.input'/ &
!    &           'cp -f fort.15 wmin.input'/ &           ! 11-29-07
     &'(limit cputime 600 ; $dir/PMIN/pmin.exe)     #  run PMIN' / &
!    &'$dir/PMIN/pmin.exe                           #  run PMIN' / &
     &'###################################################')  ! replaced  11-29-07
      WRITE (20,1665)
1665  FORMAT('if (-e core) rm -f core')
      WRITE (20,500) G92FILE2(1:NID+7), SUM, G92FILE2(1:NID+7), G92FILE2(1:NID+7), G92FILE2(1:NID+7)
500   FORMAT ('if(-e fort.21) cat fort.21 >> ',A,6A1 /&
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >>  ', A, '.sum'/'endif')             ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
      WRITE (20,1955)
1955     FORMAT ('if (-e fort.31) cp -f fort.31 fort.11'/'if ( -e fort.11 ) then')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,398)
398      FORMAT ('# ..... cell reduction # .....')
      WRITE (20,361)
361   FORMAT (' rm -f fort.15'/ &
     &        ' cp -f fort.11 fort.15  # establish PMIN cell as input'/ &
     &        ' $dir/UTILITIES/nbslattice.exe  # cell reduction')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE(20,501) G92FILE2(1:NID+7), SUM
501      FORMAT (' cat fort.11 >> ',A,6A1/' rm -f fort.11'/'endif')
      WRITE (20,889)
889   FORMAT ('endif'/' rm -f input.$i'/' @ i = ($i + 1)'/'end')
!----Finished producing the IMOL fine search packages
      WRITE (20,708)
708      FORMAT ('#'/'#****create summary tables****'/'#')
      WRITE (20,709) G92FILE2(1:NID+7)                      ! 5-4-10
709      FORMAT ('rm -f fort.*'/'cp -f ',A,'.sum fort.13')  ! 5-4-10
      WRITE(LINE_OUT,710) G92FILE2(1:NID+7), SUM
710      FORMAT ('cp -f ',A,6A1)
      LINE_OUT2 = FORT15
      WRITE (20,337) LINE_OUT
      WRITE (20,712) SEARCH_FILE
712      FORMAT ('cp -f ',A14,'  fort.14'/&
     &           '$dir/UTILITIES/table-1.exe   # make final tables'/&
     &           '  cp -f fort.16 temp.tab')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE(20,711) (SPGR(I),I=1,2),TAB
 711     FORMAT ('  mv -f fort.16 ',2A1,4A1/'  mv -f temp.tab fort.15'/&
     &           '$dir/UTILITIES/table-2.exe   # resort table'/&
     &           '  mv -f fort.16 tab.den')
!----Finish up by deleting various files
      WRITE (20,716)
716   FORMAT ('#****clean up files****')
!     IF (.NOT. I_SAVE) WRITE (20,718)
      if (NSOL .EQ. 0) write(20,718)       ! 7-9-08 DU
718      FORMAT ('  rm -f *.save *.out *last ')
719   WRITE (20,720)
720      FORMAT ('  rm -f fort.* '/ &
     &           'time'/'date'/'#'/'#**** FINISHED****')
!
      CLOSE (UNIT = 15)
      CLOSE (UNIT = 20)
      CLOSE (UNIT = 25)
      ENDDO Big_loop                 ! finished
!
!----If N_GO > 2 , then CALL WRUNJUB to write a command file
!     to make new directories, move various files to the new
!     directories and submit one job for each coord geom.
3100  IF (N_GO .GT. 1) CALL PMIN_WRUNJOB
      RETURN
      END SUBROUTINE REF_PMIN
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE PMIN_WRUNJOB              ! 7-9-08
!
      USE make_filesCommonMod, only: COMPD_ID, COMPDID, IDCHG, I, J,&
     &                               KIND, NID, POTENTIAL_FILE

      IMPLICIT NONE

      INTEGER, PARAMETER ::  N_SPGR = 54
      CHARACTER (2) ::  DIR(N_SPGR) = (/'aa','db','ah','ab','da','af',&
     &          'bh','ba','bb','av','bf','fc','dc','ai','bd','ap','aq',&
     &               'ay','fa','az','ak','am','de','dd','ca','ce','cc',&
     &               'cd','cb','au','as',&
     &             'ae', 'ac', 'ad', 'ag', 'aj', 'al', 'fd', 'an', 'ao',&
     &             'fb', 'ar', 'at', 'be', 'aw', 'bg', 'ax', 'bi', 'bc',&
     &             'bj', 'bk', 'cf', 'df', 'dg'/)
      CHARACTER (2) :: FILE(N_SPGR) = (/'AA','DB','AH','AB','DA','AF',&
     &          'BH','BA','BB','AV','BF','FC','DC','AI','BD','AP','AQ',&
     &               'AY','FA','AZ','AK','AM','DE','DD','CA','CE','CC',&
     &               'CD','CB','AU','AS', &                                ! default 31
     &              'AE', 'AC', 'AD', 'AG', 'AJ', 'AL', 'FD', 'AN', 'AO',&
     &              'FB', 'AR', 'AT', 'BE', 'AW', 'BG', 'AX', 'BI', 'BC',&
     &              'BJ', 'BK', 'CF', 'DF', 'DG'/)                         ! additional 23

!
      OPEN (UNIT=30, FILE='runjobs', STATUS='UNKNOWN')
!----Make new directories & move *.com to new directories
      WRITE (30,10) (FILE(I),I=1,29), (FILE(I),I=30,N_SPGR), &
     &              (DIR(I), I=1,29), (DIR(I), I=30,N_SPGR), &
     &              N_SPGR+1,(COMPDID(1:NID),J=1,2), &
     &              POTENTIAL_FILE                        ! 26feb03
10        FORMAT ('#!/bin/csh'/ &
     &   'setenv dir /c/Users/ACER/Downloads/PREDICTIONS'/&
!    &   'cp -f $dir/PMIN/PMIN.inp  .  # if need, to modify PMIN.inp'/&
     &            '# write new directories & move *.com file'/ &
     &            'set n=1'/'set f=',1H(,29(A2,1X),1H\,/25(A2,1X),1H)/ &
     &            'set d=',1H(,29(A2,1X),1H\,/25(A2,1X),1H)/ &
     &            'set dir = ($d)'/ &
     &            'set fl  = ($f)'/'while ($n !=',I3,')'/ &
     &            ' if (-e ',A,'_$fl[$n]_master.com) then'/ &
     &            '  if ( -d $dir[$n] ) then'/ &
     &            '   goto skip'/'  else'/'   mkdir $dir[$n]'/ &
     &            '  endif'/'skip:'/ &
     &            '   if (-e volume.min) mv -f volume.min $dir[$n]'/ &
     &            '   mv -f ',A,'_$fl[$n]_master.com $dir[$n]'/ &
     &            '   mv -f ','molpak_$fl[$n].input $dir[$n]'/ &
     &            '   mv -f ','search_$fl[$n].data $dir[$n]'/ &
     &            '   cp -f ',A,' $dir[$n]/.'/ &
!    &            '   cp -f PMIN.inp $dir[$n]/.'/ &
     &            ' endif'/'@ n = ($n + 1)'/'end'/ &                     ! 26feb03
     &            'echo " directories created & files moved "')
!----Submmit one job at a time
      WRITE (30,15) N_SPGR+1, COMPDID(1:NID)
15       FORMAT ('#****one job at a time****'/ &
     &       'date'/'time'/'set n=1'/'while ($n !=',I3,')'/ &
     &       ' if ( -d $dir[$n] ) then'/'   cd $dir[$n]'/ &
     &       '   echo "current job number: $n " '/ &
     &       '   echo "current directory : $dir[$n]" '/ &
     &       '   if ( -e ',A,'_$fl[$n]_master.com ) then'/ &
     &       '    if (-e master.log) then'/ &
     &       '     echo "the code $fl[$n] is already running" '/ &
     &       '    else'/ &
     &       '     ./*.com > master.log '/'    endif'/ &
     &       '    date'/'    time'/'   endif'/'   cd ..'/' endif'/ &
     &       '@ n = ($n + 1)'/'end'/'exit')
!----Run several summary files
        WRITE (30,20)
20         FORMAT ('#create a file named densities.list and sorted files'/&
     &             '#   named densities and energies'/&
     &             '$dir1/UTILITIES/resort-summarize.com'/&
     &             '#****FINISHED****')
      RETURN
      END SUBROUTINE PMIN_WRUNJOB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!----Produce a shell script to run molpak, dmacrys, minimum volume
!     search, nbslattice and table automatically.
!----An input atom coordinate file (molpak format) with following
!     format is required...(4A1,3F8.4) for ATOM name, X, Y, Z.
!----Also need a multipole charges file (cadpac.charges) from
!     CADPAC program.
!
      SUBROUTINE REF_DMACRYS
!
      USE make_filesCommonMod

      IMPLICIT NONE
!
!----Set initial values for symmetry flags
      INVER2 = 0
      IMR = 0
      IAX = 0
      IDCHG = 1
!----How many characters in compound ID?
      NID = 0
      DO I=1,8
         IF (COMPD_ID(I) .EQ. BLANK) GO TO 1332
         NID = NID + 1
         MASTER_NAME(I) = COMPD_ID(I)
      ENDDO
!
1332  IALL = 0
!     IF (GROUP .EQ. '  ' .OR. GROUP .EQ. 'AL' .OR. GROUP .EQ. 'al') &
      IF (GROUP .EQ. '  ' .OR. GROUP .EQ. '31' .OR. GROUP .EQ. '54') &        ! 9-18-07
     &      IALL = 1
! ################## The following lines removed #########################
!----4th line...question about COMPLETE output file
!     WRITE (*,487)
!487      FORMAT (' Should a complete output file be created and saved?', &
!    &        ' [N]: ', $)
!     I_COMP = .true.              ! should complete (HUGE) output
!     I_COMPLETE = ' '             ! file be created and saved?
!     READ (*,21) I_COMPLETE
!21   FORMAT(40A1)
!     IF (I_COMPLETE .EQ. 'N' .OR. I_COMPLETE .EQ. 'n' .OR. &
!    &       I_COMPLETE .EQ. ' ') I_COMP = .false.
!----5th line...question about ....save output file
!     WRITE (*,588)
!588      FORMAT (' Save the refined coordinate and cell parameter ',/ &
!    &           '    file (...save)? [N]: ', $)
!     I_SAVE = .true.    ! true means keep the ...save file
!     I_SAVED = ' '
!     READ (*,21) I_SAVED
!     IF (I_SAVED .EQ. 'N' .OR. I_SAVED .EQ. 'n' .OR. &
!    &       I_SAVED .EQ. ' ') I_SAVE = .false.
! #########################################################################
      I_COMP = .false.
      I_SAVE = .false.
!----Read 6th line...name of MOLPAK coordinate input file.  Try to
!     pick up file name automatically
      OPEN (UNIT=14, FILE='MOLPAK.NAME', STATUS='OLD')  ! look for files
      I = 0
45    READ (14,337,END=47) MOLPAK_NAME             ! that have
337      FORMAT (A)
      I = INDEX(MOLPAK_NAME,'molpak.xyz')          ! molpak.xyz
      IF (I .NE. 0) THEN
         L = INDEX(MOLPAK_NAME,'.xyz ') + 3
         GO TO 48
      ENDIF
      GO TO 45                                     ! at the end
47    IF (I .EQ. 0) THEN
         MOLPAK_NAME = 'molpak.xyz'
         L = 10
      ENDIF
48    WRITE (*,50) MOLPAK_NAME(1:L)
50    FORMAT('Name of the MOLPAK atom coord file', &
     &       ' [',A,']: ',$)
      READ (*,337) WHAT_LINE                              ! 7/10/03
      IF (WHAT_LINE .NE. ' ') THEN                        ! 7/10/03
         READ (WHAT_LINE,51) ATOM_LIST
51       FORMAT (30A1)
      ELSE
         ATOM1_LIST = MOLPAK_NAME
      ENDIF
!----Inquire about possible 2-fold symmetry in model
      WRITE (*,486)
 486  FORMAT ('The following questions are about molecules with', &
     &        ' C2, Cs and Ci symmetry'/ &
     &       '   and will restrict searches to space groups with', &
     &       ' appropriate symmetries;'/ &
     &       '   reply N/RETURN for a general search that ignores', &
     &       ' molecular symmetry...'/ &
     &       '    Does the molecule contain two-fold symmetry?', &
     &       ' (Y/N [N]): ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. ' ' .OR. Y_N .EQ. 'N' .OR. Y_N .EQ. 'n')  THEN
         WHAT_LINE(1:1) = ' '
      ELSE
         WRITE (*,4861)
 4861    FORMAT('    Which axis is parallel to two-fold axis?  '/ &
     &          '      1 = X, 2 = Y, 3 = Z: ',$)
         READ (*,337) WHAT_LINE
      END IF
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,4863) IAX
4863        FORMAT (I1)
         IF (IAX .GE. 1 .AND. IAX .LE. 3) THEN
            ELMSPGR(1:40) = 'AAABCAAHAFAIAKAMFAFCDA&
                           &AQAZAYBHAVBDBFCCCB'
         ELSE
            IAX = 0
         ENDIF
      ENDIF
      IF (IAX .NE. 0) GO TO 505 ! if axis, forget mirror & inversion question
!----Inquire about possible mirror plane symmetry in model
      WRITE (*,4862)
4862     FORMAT('    Does the molecule contain mirror symmetry? (Y/N', &
     &       ' [N]): ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. ' ' .OR. Y_N .EQ. 'N' .OR. Y_N .EQ. 'n')  THEN
        WHAT_LINE(1:1) = ' '
      ELSE
        WRITE (*,4868)
 4868   FORMAT('    Which axis is perpendicular to the mirror plane?  '/ &
     &         '      1 = X, 2 = Y, 3 = Z: ',$)
        READ (*,337) WHAT_LINE(1:1)
      END IF
      IMR = 0
      IF (WHAT_LINE .NE. INPUT_LINE)  READ (WHAT_LINE,4863) IMR
      IF (IMR .NE. 0) GO TO 505  ! if mirror, forget inversion question
!----Inquire about possible inversion symmetry in model
      WRITE (*,496)
496      FORMAT('    Does the molecule contain a center of symmetry', &
     &          ' at the origin? (Y/N [N]): ',$)
      READ (*,337) INVER1
      IF (INVER1 .EQ. ' ' .OR. INVER1 .EQ. 'N' .OR. &
     &    INVER1 .EQ. 'n') THEN
          INVER2 = 0
      ELSE
         ELMSPGR(1:30) = 'AAAHAFDADBAPBABBAQAZAYBHAVBDBF'
         INVER2 = 1
      ENDIF
!----Open atom coordinate file
505   OPEN (UNIT=15, FILE=ATOM1_LIST, STATUS='OLD')
      WRITE (*,109)
109      FORMAT (4X,'atom id     x         y         z      Id    ', &
     &        '6-31g*')
      N = 0
111   N = N + 1
        READ (15,1,END=100) (ATOM_TY(N,J),J=1,2), (ATOM(N,I),I=1,4), &
     &                    X(N), Y(N), Z(N), IDATOM(N), G92CHARGE(N)
1          FORMAT (5X,2A1,4A1,3F10.6,I5,F10.6)
        WRITE (*,112) (ATOM_TY(N,J),J=1,2), (ATOM(N,I),I=1,4), &
     &                X(N), Y(N), Z(N), IDATOM(N), G92CHARGE(N)
112        FORMAT(4X,2A1,4A1,3F10.6,I5,F10.6)
      GO TO 111
100   N = N -1
!----Question about using an existing solution
      NSOL = 0             ! solution #
      Y_N = ' '
      WRITE (*,10)
10       FORMAT ('Do you want to rerun a solution from an existing', &
     &           ' volume.min file [N]: ',$)
      READ (*,337) Y_N
      IF (Y_N .EQ. 'Y' .OR. Y_N .EQ. 'y') THEN
          WRITE (*,12)
12           FORMAT (' Enter solution number: ',$)
          READ (*,*) NSOL
          WHAT_LINE = '1'
          OPEN (UNIT=30, FILE='volume.min', STATUS="OLD", ERR=800)
          GO TO 801
800       WRITE (*,804)
804          FORMAT (' The required volume.min file does not exist!')
          STOP
801       READ (30,*)                 ! skip first line
          READ (30,'(9X,3F8.1)') (ZSTEP(I),I=1,3)
          DO L=1,3
             ZMINA(L) = -90.
             ZMAXA(L) =  90.
          END DO
      ELSE
!----Read 7th line...angle search limits for initial molpak search
          WRITE (*,80)
80            FORMAT('Minimum, maximum and step for 3 angles for 1st ', &
     &               'MOLPAK map'/ &
     &               '  [-90,90,10,-90,90,10,-90,90,10,]: ',$)
          READ (*,337) WHAT_LINE
          IF (WHAT_LINE .NE. INPUT_LINE) THEN
              READ (WHAT_LINE,*) (ZMINA(I),ZMAXA(I),ZSTEP(I),I=1,3)
          ELSE
!----Set up defaults...
              DO L=1,3
                 ZMINA(L) = -90.
                 ZMAXA(L) =  90.
                 ZSTEP(L) =  10.
              ENDDO
          ENDIF
!----Read number of minimum volume points to be saved in MOLPAK search
             WRITE (*,703)
703             FORMAT ('Number of minimum volumes from initial MOLPAK', &
     &                  ' search to pass to next step,'/  &
     &                  '   suggest 7000 for 10 deg step search and', &
     & ' 60000 for 5 deg step search or range(-180 to 180) [7000]: ',$)
             WHAT_LINE = INPUT_LINE
             READ (*,337) WHAT_LINE
             IF (WHAT_LINE .EQ. INPUT_LINE) THEN
                N_VOLS_SAVED = 7000
             ELSE
                READ (WHAT_LINE,*) N_VOLS_SAVED
             ENDIF
!----Read number of minimum volume points to be elaborated
             WRITE (*,803)
803             FORMAT('Maximum number of 2nd MOLPAK calcns desired [500]: ',$)
             WHAT_LINE = INPUT_LINE
             READ (*,337) WHAT_LINE
      END IF
      IF (WHAT_LINE .NE. INPUT_LINE) THEN
         READ (WHAT_LINE,*) NMOL
      ELSE
         NMOL = 500       ! default -> elaborate 500 solns
      ENDIF
!----How many coord geoms have been specified?  Each 2 letter code
!     should be separated from its neighbors by a space or comma.
      IF (IALL .NE. 0) GO TO 29
      N_SPGRPS = 0
      NSCHARS = 0
      DO 22 I=1,40
         IF (SPGR(I) .EQ. ' ' .OR. SPGR(I) .EQ. ',') GO TO 22
         IF (NSCHARS .NE. 0) THEN
            IF (NSCHARS .EQ. 1) GO TO 24
            NSCHARS = 0
         ENDIF
23       N_SPGRPS = N_SPGRPS + 1
24       NSCHARS = NSCHARS + 1
         MOL_SPGRP(NSCHARS,N_SPGRPS) = SPGR(I)
22    CONTINUE
!----All or just some space groups
29    IF (IALL .EQ. 0) THEN
         N_GO = N_SPGRPS
      ELSE
       IF (GROUP .EQ. '  ' .OR. GROUP .EQ. '31') THEN    ! 9-18-07
         N_GO = N_SPGR-23                                ! 9-18-07
       ELSE
         N_GO = N_SPGR
       END IF
      ENDIF
!
!----Big LOOP for creating the N_GO sets of files for execution
!
Big_loop : DO LGRPS=1,N_GO
               IF (IALL .EQ. 0) THEN
                  SPGR(1) = MOL_SPGRP(1,LGRPS)
                  SPGR(2) = MOL_SPGRP(2,LGRPS)
               ELSE
                  DO J=1,N_GO
                     IF (ALL_SPGR(LGRPS) .EQ. ELM_SPGR(J)) CYCLE
                  END DO
                  SPGR(1) = ALL_SPGR(LGRPS)(1:1)
                  SPGR(2) = ALL_SPGR(LGRPS)(2:2)
               ENDIF
!----Go through the N_SPGRPS coord geoms & make appropriate files.
!     Tack 3 characters (_SPGR) onto end of COMPD_ID and MASTER_NAME.
      COMPD_ID(NID+1) = '_'
      COMPD_ID(NID+2) = OLD(1)
      COMPD_ID(NID+3) = OLD(2)
      COMPD_ID(NID+4) = OLD(3)
      COMPD_ID(NID+5) = '_'
      COMPD_ID(NID+6) = SPGR(1)
      COMPD_ID(NID+7) = SPGR(2)
      MASTER_NAME(NID+1) = '_'
      MASTER_NAME(NID+2) = SPGR(1)
      MASTER_NAME(NID+3) = SPGR(2)
!----Construct name for G98 file
      DO I=1,NID
         NEWFILE(I)=COMPD_ID(I)
      END DO
      NEWFILE(NID+1) = '_'
      NEWFILE(NID+2) = NEW(1)
      NEWFILE(NID+3) = NEW(2)
      NEWFILE(NID+4) = NEW(3)
      NEWFILE(NID+5) = '_'
      NEWFILE(NID+6) = SPGR(1)
      NEWFILE(NID+7) = SPGR(2)
!----Construct name for MASTER.COM file that will have compound ID 1st,
!     space group symbol 2nd and "master.com" 3rd.
      MID = NID + 3
      DO I=1,11
         MASTER_NAME(MID+I) = MASTERCOM(I)
      ENDDO
      MOLPAK_FILE(1:7) = 'molpak_'      ! create name of the
      MOLPAK_FILE(8:9) = GROUP          ! molpak input file
      MOLPAK_FILE(10:15) = '.input'
      SEARCH_FILE(1:7) = 'search_'      ! create name for search
      SEARCH_FILE(8:9) = GROUP          ! data file
      SEARCH_FILE(10:14) = '.data'
!----Open files
      OPEN (UNIT=20, FILE=MASTER_NAME2, STATUS='UNKNOWN')
      OPEN (UNIT=25, FILE=MOLPAK_FILE,  STATUS='UNKNOWN')
!----Start to write master.com... files
      WRITE (20,39)(NMOL+1)
39    FORMAT ('#!/bin/csh'/'pwd'/'#'/'time'/'date'/'#'/ &
     &   'setenv dir /c/Users/ACER/Downloads/PREDICTIONS'/ &
     &       'mkdir FDAT-files'/ &
     &       '# ********** initial MOLPAK search **********'/'#'/ &
     &       '# if input.* is existing, restart DMACRYS'/ &
     &       'set k=1'/'while ($k != ',I4,')'/ &
     &       ' if ( -e input.$k ) then'/'  goto rundmacrys'/' endif'/ &
     &       ' @ k = ($k + 1)'/'end'/ &
     &       '# if volume.min  is existing, skip first MOLPAK search'/ &
     &       'if ( -e volume.min  ) then'/' cp volume.min  fort.11'/ &
     &       ' goto mkcom'/'else'/' if ( -e volume.find ) goto mkcom'/ &
     &       'endif'/'# if molpak.vol exists skip MOLPAK search'/ &
     &       'if ( -e molpak.vol ) then'/ &
     &       ' cp molpak.vol molpak.min'/'endif'/)
      WRITE (20,41) MOLPAK_FILE
41    FORMAT('cp ',A15,'  fort.15  # MOLPAK input file to reader'/ &
     &       '#                   minimum vol output on unit # 8'/ &
     &       '$dir/MOLPAK/molpak_12.exe',' # run MOLPAK'/ &
     &       'if ( -e molpak.min ) then'/ &
     &       '  rm -f molpak.min'/'endif'/'rm -f molpak.vol'/ &
     &       'time'/'date'/'#')
      NDIF = 30 - MID
!----Write title on MOLPAK input file + rest of it
      WRITE (25,245) COMPDID(1:NID), TITLE
245   FORMAT ('HEAD  ',A,8X,A)
!     using default  2-27-04
!     WRITE (25,146)                      ! line # 2...NSEG values
146   FORMAT ('NSEG  128 128   8')
!----Place atom coordinates into molpak input file
      DO L=1,N
         WRITE(25,483) (ATOM_TY(L,J),J=1,2), (ATOM(L,I),I=1,4), &
     &                 X(L),Y(L),Z(L),IDATOM(L),G92CHARGE(L)
483         FORMAT('ATOM ',2A1,4A1,3F10.6,I5,F10.6)
      ENDDO
      IF (INVER2 .EQ. 0) GO TO 7631 ! INVER2 = 0  means no center
      WRITE (25,764)                ! write "CENT" line into
764       FORMAT ('CENT')           ! MOLPAK input file
      GO TO 488
 7631 IF (IMR .EQ. 0) GO TO 763    ! IMR = 0 means no mirror
      WRITE (25,7632) IMR          ! write "PLAN" line into
 7632    FORMAT ('PLAN',I2)        ! MOLPAK input
      GO TO 488
763   IF (IAX .EQ. 0) GO TO 488    ! IAX = 0 means no 2-fold axis
      WRITE (25,787) IAX           ! write "AXIS" line into
787      FORMAT ('AXIS',I2)        ! MOLPAK input file
488   IF (I_COMP)  WRITE (20,491) (COMPD_ID(I),I=1,MID), COMPLT
491                FORMAT ('mv -f fort.16 ',46A1)
      WRITE (25,246) GROUP
246   FORMAT('INCL ',A2)
      IF (IAX .EQ. 0) THEN
         WRITE (25,147) N_VOLS_SAVED  ! # min vols for 1st MOLPAK
147         FORMAT ('VOLS ',I5,'    8')
         WRITE(25,81) (ZMINA(I),ZMAXA(I),ZSTEP(I),I=1,3) ! angle search info
81          FORMAT('SEEK 1',9F6.1)
      ELSE
         WRITE (25,148)
 148        FORMAT ('VOLS    19    8'/'SEEK') ! special for molec on 2-fold axis
      ENDIF
      WRITE (25,91)
91       FORMAT('WMIN',11X/'FINI')
      WRITE (*,1111) MASTER_NAME
1111     FORMAT (' >>>>>>Name of the new COM file: ',41A1)
      CLOSE (UNIT=25)
!
!----Finished writing MOLPAK input file on unit # 25
!
!----Data for searches on unit # 25
      OPEN (UNIT=25, FILE=SEARCH_FILE, STATUS='UNKNOWN')
      IF (IDCHG .EQ. 1)&
     &    WRITE (25,231) NEWFILE2(1:NID+7), GROUP, TITLE
231          FORMAT (A/A2/A73)
      WRITE (25,'(I5,3F5.1,I5)') NMOL, ZSTEP, NSOL ! and soln to be elaborated
      CLOSE (UNIT=25)
!----Section for minimum volume search
      WRITE (20,442)
442   FORMAT('#'/ &
     &       '#****minimum volume search****'/'#'/&
     &       'mv -f fort.8 fort.10  # min vol file is now unit # 10')
      WRITE (20,448) SEARCH_FILE
448   FORMAT ('rm -f fort.15'/'cp -f ',A14,'  fort.15'/ &
     &        '$dir/UTILITIES/min-vol-search.exe', &
     &        '   # make unique min V list')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
493                FORMAT ('cat fort.16 >> ',46A1)
      WRITE (20,551) MOLPAK_FILE, SEARCH_FILE
551   FORMAT ('#'/'#****make MOLPAK command files****'/ &
     &        '#'/'cp fort.11 volume.min '/'mkcom:'/ &
     &        'if ( -e fort.11 ) mv fort.11 fort.8  '/ &
     &        'rm -f fort.10 fort.15' / &
     &        'cp ',A15,'  fort.10'/ &
     &        'cp ',A14,'  fort.15  # small angle search'/ &
     &        '#          and number of initial volumes to elaborate'/ &
     &        '$dir/UTILITIES/make-molpak-com-find.exe'/ &
     &        'rm -f fort.*'/ &
     &        '# make the MOLPAK/DMARCRYS/NBSLATTICE files')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
!----FINISHED
!***************** LOOP THROUGH NMOL RUN******************************
!----Assign files for running MOLPAK
      WRITE (20,600) NMOL+1
600      FORMAT ('#'/'# ****MOLPAK/DMACRYS/CELL REDUCTION****'/'#'/ &
     &           'rundmacrys:'/'set i=1'/'while ($i != ',I4,')')
      WRITE (20,456)
456      FORMAT ('if ( -e input.$i ) then'/ &
     &           'cp input.$i fort.15    # MOLPAK input file')
      WRITE (20,458)
458      FORMAT ('#                 ...xyz file on unit # 9'/ &
     &           '#                 ...sum file on unit # 21'/ &
     &           ' echo " MOLPAK NO: $i "'/ &
     &           '$dir/MOLPAK/molpak_12.exe     # run MOLPAK')
      IF (I_COMP) WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,1951)
1951     FORMAT (' if ($i == 1) then')
      IF (IDCHG.EQ.1) THEN
         WRITE (20,1952) NEWFILE2(1:NID+7), SUM,NEWFILE2(1:NID+7), SUM
         WRITE (20,952) NEWFILE2(1:NID+7), SUM,NEWFILE2(1:NID+7), SUM
      END IF
1952     FORMAT('   echo "solution # $i .... " > ',A,6A1/ &
     &'   cat fort.21 >> ',A,6A1/'   rm -f fort.21')
952      FORMAT(' else'/'   echo "solution # $i .... " >> ',A,6A1/ &
     &'   cat fort.21 >> ',A,6A1/ &
     &       '   rm -f fort.21')
      WRITE (20,955)
955      FORMAT (' endif'/ &
     &           'if (-e molpak.cell) then'/ &
     &           'cat fort.49 >> molpak.cell'/'else'/ &
     &           'mv fort.49 molpak.cell'/'endif'/ &
     &           '#****run DMACRYS****')
      GO TO (7775, 7777, 7776), I_POTE        ! 3-19-09 DU
!     IF(I_POTE.EQ.2) GOTO 7777               ! 10-9-07
7775  IF (NSOL .EQ. 0) THEN
      WRITE (20,957)SEARCH_FILE,SEARCH_FILE,NEWFILE2(1:NID+7),SUM,NEWFILE2(1:NID+7), NEWFILE2(1:NID+7), &
     &              NEWFILE2(1:NID+7)                       ! 9-4-14
957   FORMAT ('# fort.9 from MOLPAK as an input of wmin2fdat'/ &
     &'#cp fort.9 wmininput.$i'/ &                          ! temp
     &'rm -f fort.15'/ &
     &'cp -f ',A14,' fort.15'/ &
     &'$dir/UTILITIES/wmin2fdat-umd.exe'/'cp fort.11 $i.fdat'/ &
     &'mv fort.11 FDAT-files/$i.fdat'/ &
     &'# 2nd neighbour (read cadpac.charges file)'/ &
     &'$dir/DMACRYS/NEIGHCRYS_2.0.4/neighcrys.out << EOF'/'I'/ &
     &'$i.fdat'/'bondlengths'/'4.0'/'n'/'n'/'n'/'n'/'0'/'y'/ &  ! H be not adjusted
     &'cadpac.charges'/'y'/'dmarel.axis'/'n'/'y'/ &
     &'pote.dat'/'EOF'/'mv fort.20 symm.20.$i'/ &
     &'rm -f *.mac *.ccl *.nem *.nnl fort.21'/'mv symm.20.$i fort.20'/ &
     &'echo " 2nd neighbour finished " '/'# 2nd run DMACRYS'/ &
     &'mv $i.dmain $i-2nd.dmain'/ &
     &'$dir/DMACRYS/DMACRYS_2.0.4/dmacrys.out <$i-2nd.dmain> $i.dmacrys-out'/ &
     &'#echo " DMACRYS No: $i " >> dmacrys.message'/ &
     &'#cp $i.dmacrys-out fort.25'/'#$dir/UTILITIES/read-out.exe'/ &
     &'#cat fort.26 >> dmacrys.message'/'#rm -f fort.25 fort.26'/ &
     &'#compress $i.dmacrys-out'/&
     &'rm -f $i.dmacrys-out'/ &
     &'#mv fort.15 FDAT-files/$i.fdat'/'mv fort.16 $i.SHELX'/&
     &'mv fort.12 $i.summary'/'mv fort.13 $i.save'/ &
     &'rm -f graphic $i-2nd.dmain'/'echo "2nd DMACRYS finished"'/ &
     &'if (-e fort.21) rm -f fort.21'/ &
     &'rm -f fort.15'/'cp -f ',A14,' fort.15'/ &
     &'if (-e fort.5) rm -f fort.5'/'cp -f $i.summary   fort.5'/ &
     &'$dir/UTILITIES/read-sum-4moldma.exe < fort.5'/ &
     &'rm -f fort.15'/ &
     &'if (-e fort.21) then'/' cat fort.21 >> ',A,6A1/ &
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >>  ', A, '.sum'/'endif')             ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
ELSE
      WRITE (20,958)SEARCH_FILE,SEARCH_FILE,NEWFILE2(1:NID+7),SUM, NEWFILE2(1:NID+7),NEWFILE2(1:NID+7), &
     &              NEWFILE2(1:NID+7)                      ! 9-4-14  run special solution
958   FORMAT ('# fort.9 from MOLPAK as an input of wmin2fdat'/ &
     &'#cp fort.9 wmininput.$i'/ &                          ! temp
     &'rm -f fort.15'/ &
     &'cp -f ',A14,' fort.15'/ &
     &'$dir/UTILITIES/wmin2fdat-umd.exe'/'cp fort.11 $i.fdat'/ &
     &'mv fort.11 FDAT-files/$i.fdat'/ &
     &'# 2nd neighbour (read cadpac.charges file)'/ &
     &'$dir/DMACRYS/NEIGHCRYS_2.0.4/neighcrys.out << EOF'/'I'/ &
     &'$i.fdat'/'bondlengths'/'4.0'/'n'/'n'/'n'/'n'/'0'/'y'/ &  ! H be not adjusted
     &'cadpac.charges'/'y'/'dmarel.axis'/'n'/'y'/'pote.dat'/'EOF'/ &
     &'mv fort.20 temp.20'/'rm -f  *.mac *.ccl *.nem *.nnl fort.21'/ &
     &'mv temp.20 fort.20'/ &
     &'echo " 2nd neighbour finished " '/'# 2nd run DMACRYS'/ &
     &'mv $i.dmain $i-2nd.dmain'/ &
     &'$dir/DMACRYS/DMACRYS_2.0.4/dmacrys.out <$i-2nd.dmain> $i.dmacrys-out'/ &
     &'#echo " DMACRYS No: $i " >> dmacrys.message'/ &
     &'#cp $i.dmacrys-out fort.25'/'#$dir/UTILITIES/read-out.exe'/ &
     &'#cat fort.26 >> dmacrys.message'/'#rm -f fort.25 fort.26'/ &
     &'#compress $i.dmacrys-out'/ &
     &'rm -f $i.dmacrys-out'/ &
     &'#mv fort.15 $i.fdat'/'mv fort.16 $i.SHELX'/&
     &'mv fort.12 $i.summary'/'mv fort.13 $i.save'/ &
     &'rm -f graphic $i-2nd.dmain'/'echo "2nd DMACRYS finished"'/ &
     &'if (-e fort.21) rm -f fort.21'/ &
     &'rm -f fort.15'/'cp -f ',A14,' fort.15'/ &
     &'if (-e fort.5) rm -f fort.5'/'cp -f $i.summary   fort.5'/ &
     &'$dir/UTILITIES/read-sum-4moldma.exe < fort.5'/ &
     &'rm -f fort.15'/ &
     &'if (-e fort.21) then'/' cat fort.21 >> ',A,6A1/  &
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >> ', A, '.sum'/'endif')              ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
      END IF
      GO TO 7778
7776  IF (NSOL .EQ. 0) THEN
      WRITE (20,9570)SEARCH_FILE,SEARCH_FILE,NEWFILE2(1:NID+7),SUM, NEWFILE2(1:NID+7),NEWFILE2(1:NID+7), &
     &               NEWFILE2(1:NID+7)                     ! 9-4-14
9570  FORMAT ('# fort.9 from MOLPAK as an input of wmin2fdat'/ &
     &'#cp fort.9 wmininput.$i'/ &                          ! temp
     &'rm -f fort.15'/ &
     &'cp -f ',A14,' fort.15'/ &
     &'$dir/UTILITIES/wmin2fdat-umd.exe'/'cp fort.11 $i.fdat'/ &
     &'mv fort.11 FDAT-files/$i.fdat'/ &
     &'# 2nd neighbour (read cadpac.charges file)'/ &
     &'$dir/DMACRYS/NEIGHCRYS_2.0.4/neighcrys.out << EOF'/'I'/ &
     &'$i.fdat'/'bondlengths'/'4.0'/'n'/'y'/'w'/'y'/'0'/'y'/ &  ! H be not adjusted
     &'cadpac.charges'/'y'/'dmarel.axis'/'n'/'y'/ &
     &'will01.pots'/'EOF'/'mv fort.20 symm.20.$i'/ &
     &'rm -f *.mac *.ccl *.nem *.nnl fort.21'/'mv symm.20.$i fort.20'/ &
     &'echo " 2nd neighbour finished " '/'# 2nd run DMACRYS'/ &
     &'mv $i.dmain $i-2nd.dmain'/ &
     &'$dir/DMACRYS/DMACRYS_2.0.4/dmacrys.out <$i-2nd.dmain> $i.dmacrys-out'/ &
     &'#echo " DMACRYS No: $i " >> dmacrys.message'/ &
     &'#cp $i.dmacrys-out fort.25'/'#$dir/UTILITIES/read-out.exe'/ &
     &'#cat fort.26 >> dmacrys.message'/'#rm -f fort.25 fort.26'/ &
     &'#compress $i.dmacrys-out'/&
     &'rm -f $i.dmacrys-out'/ &
     &'#mv fort.15 FDAT-files/$i.fdat'/'mv fort.16 $i.SHELX'/&
     &'mv fort.12 $i.summary'/'mv fort.13 $i.save'/ &
     &'rm -f graphic $i-2nd.dmain'/'echo "2nd DMACRYS finished"'/ &
     &'if (-e fort.21) rm -f fort.21'/ &
     &'rm -f fort.15'/'cp -f ',A14,' fort.15'/ &
     &'if (-e fort.5) rm -f fort.5'/'cp -f $i.summary   fort.5'/ &
     &'$dir/UTILITIES/read-sum-4moldma.exe < fort.5'/ &
     &'rm -f fort.15'/ &
     &'if (-e fort.21) then'/' cat fort.21 >> ',A,6A1 /&
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >> ', A, '.sum'/'endif')              ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
ELSE
      WRITE (20,9580)SEARCH_FILE,SEARCH_FILE,NEWFILE2(1:NID+7),SUM, NEWFILE2(1:NID+7), NEWFILE2(1:NID+7), &
     &               NEWFILE2(1:NID+7)                     ! 9-4-14 run special solution
9580  FORMAT ('# fort.9 from MOLPAK as an input of wmin2fdat'/ &
     &'#cp fort.9 wmininput.$i'/ &                          ! temp
     &'rm -f fort.15'/ &
     &'cp -f ',A14,' fort.15'/ &
     &'$dir/UTILITIES/wmin2fdat-umd.exe'/'cp fort.11 $i.fdat'/ &
     &'mv fort.11 FDAT-files/$i.fdat'/ &
     &'# 2nd neighbour (read cadpac.charges file)'/ &
     &'$dir/DMACRYS/NEIGHCRYS_2.0.4/neighcrys.out << EOF'/'I'/ &
     &'$i.fdat'/'bondlengths'/'4.0'/'n'/'y'/'w'/'y'/'0'/'y'/ &  ! H be not adjusted
     &'cadpac.charges'/'y'/'dmarel.axis'/'n'/'y'/'will01.pots'/'EOF'/ &
     &'mv fort.20 temp.20'/'rm -f  *.mac *.ccl *.nem *.nnl fort.21'/ &
     &'mv temp.20 fort.20'/ &
     &'echo " 2nd neighbour finished " '/'# 2nd run DMACRYS'/ &
     &'mv $i.dmain $i-2nd.dmain'/ &
     &'$dir/DMACRYS/DMACRYS_2.0.4/dmacrys.out <$i-2nd.dmain> $i.dmacrys-out'/ &
     &'#echo " DMACRYS No: $i " >> dmacrys.message'/ &
     &'#cp $i.dmacrys-out fort.25'/'#$dir/UTILITIES/read-out.exe'/ &
     &'#cat fort.26 >> dmacrys.message'/'#rm -f fort.25 fort.26'/ &
     &'#compress $i.dmacrys-out'/ &
     &'rm -f $i.dmacrys-out'/ &
     &'#mv fort.15 $i.fdat'/'mv fort.16 $i.SHELX'/&
     &'mv fort.12 $i.summary'/'mv fort.13 $i.save'/ &
     &'rm -f graphic $i-2nd.dmain'/'echo "2nd DMACRYS finished"'/ &
     &'if (-e fort.21) rm -f fort.21'/ &
     &'rm -f fort.15'/'cp -f ',A14,' fort.15'/ &
     &'if (-e fort.5) rm -f fort.5'/'cp -f $i.summary   fort.5'/ &
     &'$dir/UTILITIES/read-sum-4moldma.exe < fort.5'/ &
     &'rm -f fort.15'/ &
     &'if (-e fort.21) then'/' cat fort.21 >> ',A,6A1 / &
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >> ', A, '.sum'/'endif')              ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
     END IF
     GO TO 7778
7777 IF (NSOL .EQ. 0) THEN                                                ! 10-9-07
      WRITE (20,9577)SEARCH_FILE,SEARCH_FILE,NEWFILE2(1:NID+7),SUM, NEWFILE2(1:NID+7),NEWFILE2(1:NID+7), &
     &               NEWFILE2(1:NID+7)                     ! 9-4-14
9577  FORMAT ('# fort.9 from MOLPAK as an input of wmin2fdat'/ &
     &'#cp fort.9 wmininput.$i'/ &                          ! temp
     &'rm -f fort.15'/ &
     &'cp -f ',A14,' fort.15'/ &
     &'$dir/UTILITIES/wmin2fdat-umd.exe'/'cp fort.11 $i.fdat'/ &
     &'mv fort.11 FDAT-files/$i.fdat'/ &
     &'# 2nd neighbour (read cadpac.charges file)'/ &
     &'$dir/DMACRYS/NEIGHCRYS_2.0.4/neighcrys.out << EOF'/'I'/ &
     &'$i.fdat'/'bondlengths'/'4.0'/'n'/'y'/'n'/'n'/'0'/'y'/ &  ! H be adjusted   ! 10-9-07
     &'cadpac.charges'/'y'/'dmarel.axis'/'n'/'y'/ &
     &'pote.dat'/'EOF'/'mv fort.20 symm.20.$i'/ &
     &'rm -f *.mac *.ccl *.nem *.nnl fort.21'/'mv symm.20.$i fort.20'/ &
     &'echo " 2nd neighbour finished " '/'# 2nd run DMACRYS'/ &
     &'mv $i.dmain $i-2nd.dmain'/ &
     &'$dir/DMACRYS/DMACRYS_2.0.4/dmacrys.out <$i-2nd.dmain> $i.dmacrys-out'/ &
     &'#echo " DMACRYS No: $i " >> dmacrys.message'/ &
     &'#cp $i.dmacrys-out fort.25'/'#$dir/UTILITIES/read-out.exe'/ &
     &'#cat fort.26 >> dmacrys.message'/'#rm -f fort.25 fort.26'/ &
     &'#compress $i.dmacrys-out'/&
     &'rm -f $i.dmacrys-out'/ &
     &'#mv fort.15 FDAT-files/$i.fdat'/'mv fort.16 $i.SHELX'/&
     &'mv fort.12 $i.summary'/'mv fort.13 $i.save'/ &
     &'rm -f graphic $i-2nd.dmain'/'echo "2nd DMACRYS finished"'/ &
     &'if (-e fort.21) rm -f fort.21'/ &
     &'rm -f fort.15'/'cp -f ',A14,' fort.15'/ &
     &'if (-e fort.5) rm -f fort.5'/'cp -f $i.summary   fort.5'/ &
     &'$dir/UTILITIES/read-sum-4moldma.exe < fort.5'/ &
     &'rm -f fort.15'/ &
     &'if (-e fort.21) then'/' cat fort.21 >> ',A,6A1/ &
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >> ', A, '.sum'/'endif')              ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
      ELSE
      WRITE (20,9587)SEARCH_FILE,SEARCH_FILE,NEWFILE2(1:NID+7),SUM, NEWFILE2(1:NID+7), NEWFILE2(1:NID+7), &
     &               NEWFILE2(1:NID+7)                     ! 9-4-14  run special solution
9587  FORMAT ('# fort.9 from MOLPAK as an input of wmin2fdat'/ &
     &'#cp fort.9 wmininput.$i'/ &                          ! temp
     &'rm -f fort.15'/ &
     &'cp -f ',A14,' fort.15'/ &
     &'$dir/UTILITIES/wmin2fdat-umd.exe'/'cp fort.11 $i.fdat'/ &
     &'mv fort.11 FDAT-files/$i.fdat'/ &
     &'# 2nd neighbour (read cadpac.charges file)'/ &
     &'$dir/DMACRYS/NEIGHCRYS_2.0.4/neighcrys.out << EOF'/'I'/ &
     &'$i.fdat'/'bondlengths'/'4.0'/'n'/'y'/'n'/'n'/'0'/'y'/ &      ! H be adjusted   ! 10-9-07
     &'cadpac.charges'/'y'/'dmarel.axis'/'n'/'y'/'pote.dat'/'EOF'/ &
     &'mv fort.20 temp.20'/'rm -f  *.mac *.ccl *.nem *.nnl fort.21'/ &
     &'mv temp.20 fort.20'/ &
     &'echo " 2nd neighbour finished " '/'# 2nd run DMACRYS'/ &
     &'mv $i.dmain $i-2nd.dmain'/ &
     &'$dir/DMACRYS/DMACRYS_2.0.4/dmacrys.out <$i-2nd.dmain> $i.dmacrys-out'/ &
     &'#echo " DMACRYS No: $i " >> dmacrys.message'/ &
     &'#cp $i.dmacrys-out fort.25'/'#$dir/UTILITIES/read-out.exe'/ &
     &'#cat fort.26 >> dmacrys.message'/'#rm -f fort.25 fort.26'/ &
     &'#compress $i.dmacrys-out'/ &
     &'rm -f $i.dmacrys-out'/ &
     &'#mv fort.15 $i.fdat'/'mv fort.16 $i.SHELX'/&
     &'mv fort.12 $i.summary'/'mv fort.13 $i.save'/ &
     &'rm -f graphic $i-2nd.dmain'/'echo "2nd DMACRYS finished"'/ &
     &'if (-e fort.21) rm -f fort.21'/ &
     &'rm -f fort.15'/'cp -f ',A14,' fort.15'/ &
     &'if (-e fort.5) rm -f fort.5'/'cp -f $i.summary   fort.5'/ &
     &'$dir/UTILITIES/read-sum-4moldma.exe < fort.5'/ &
     &'rm -f fort.15'/ &
     &'if (-e fort.21) then'/' cat fort.21 >> ',A,6A1/ &
!    & 'if ($i == 1) then' / &                             ! 5-4-10
!    & 'cp -f fort.21 ',A,'.sum' / &                       ! 5-4-10
!    & 'else' / &                                          ! 5-4-10
!    & 'cat  fort.21 >> ', A, '.sum'/'endif')              ! 5-4-10
     & 'if (-e ',A,'.sum) then' / &                        ! 9-4-14
     & ' cat fort.21 >> ', A, '.sum'/'else'/ &             ! 9-4-14
     & ' cp -f fort.21 ',  A, '.sum'/'endif')              ! 9-4-14
      END IF
7778  WRITE (20,361)
361      FORMAT (' mv fort.11 fort.15  # establish DMACRYS cell as input'/ &
     &           ' $dir/UTILITIES/nbslattice.exe  # cell reduction')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      IF (IDCHG .EQ. 1) THEN
         WRITE (20,501) NEWFILE2(1:NID+7), SUM
501         FORMAT (' cat fort.11 >> ',A,6A1/' rm -f fort.* '/ &
     &              'endif'/'endif')
         IF (NSOL .EQ. 0) WRITE (20,502)
502      FORMAT('rm -f $i.summary $i.SHELX $i.save')
         WRITE (20,503)
503      FORMAT('rm -f input.$i $i.fdat'/ &
     &              'if (-e core) rm -f core'/'@ i = ($i + 1)'/'end')
      END IF
!----Finished producing the IMOL fine search packages
      WRITE (20,708)
708   FORMAT ('#'/'#****summary table****'/'#')
      WRITE (20,709) NEWFILE2(1:NID+7)                      ! 5-4-10
709      FORMAT ('rm -f fort.*'/'cp -f ',A,'.sum fort.13')  ! 5-4-10
      IF (IDCHG .EQ. 1) THEN
         WRITE (LINE_OUT,710) NEWFILE2(1:NID+7), SUM
      END IF
710   FORMAT ('cp ',A,6A1)
      LINE_OUT2 = FORT15
      WRITE (20,337) LINE_OUT
      WRITE (20,712) SEARCH_FILE
712      FORMAT ('cp ',A14,'  fort.14'/ &
     &           '$dir/UTILITIES/table-1.exe       # make final tables'/ &
     &           'cp fort.16 temp.tab')
      IF (I_COMP)  WRITE (20,493) (COMPD_ID(I),I=1,MID), COMPLT
      WRITE (20,711) (SPGR(I),I=1,2),TAB
 711     FORMAT ('  mv fort.16 ',2A1,4A1/ &
     &          'mv -f temp.tab fort.15'/ &
     &          '$dir/UTILITIES/table-2.exe     # resort table'/ &
     &          'mv fort.16 tab.den')
!----Finish up by deleting various files
!     IF (.NOT. I_SAVE) WRITE (20,718)
!718                    FORMAT ('rm -f *.save')
719   WRITE (20,720)
720      FORMAT ('rm -f fort.* '/ &
     &           '#rm -f *.summary *.fdat *.SHELX'/ &
     &           'cd FDAT-files'/'tar -cvf FDAT.tar *.fdat'/&
     &           'gzip FDAT.tar'/'rm -f *.fdat'/'cd ../'/ &
     &           'time'/'date'/'#  FINISHED')
!
      CLOSE (UNIT = 15)
      CLOSE (UNIT = 20)
      CLOSE (UNIT = 25)
      ENDDO Big_loop                ! finished
!
!-----If N_GO > 2 , CALL DMA_WRUNJOB to create script to make
!      new directories, move files, etc.
      IF (N_GO .GE. 2) CALL DMA_WRUNJOB
      RETURN
      END SUBROUTINE REF_DMACRYS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!----Write a script to make new directories, move files to new
!     directories and submit jobs.
!     special for running MOLPAK-DMACRYS w/dmacrys
!
      SUBROUTINE DMA_WRUNJOB

      USE make_filesCommonMod, only: COMPD_ID, COMPDID, I, IDCHG, J, KIND,&
     &                               NID, N_SPGR, POTENTIAL_FILE, I_POTE

      IMPLICIT NONE

      CHARACTER (2) :: DIR(N_SPGR) = (/'af','bh','bd','bf','ah','ay',&
     &             'az','av','ab','aq','ai','cc','dc','ak','ca','am',&
     &             'fa','fc','da','db','dd','de','cb','aa','ap','ba', &
     &             'bb','cd','ce','au','as',&
     &             'ae', 'ac', 'ad', 'ag', 'aj', 'al', 'fd', 'an', 'ao',&
     &             'fb', 'ar', 'at', 'be', 'aw', 'bg', 'ax', 'bi', 'bc',&
     &             'bj', 'bk', 'cf', 'df', 'dg'/)
      CHARACTER (2) :: FILE(N_SPGR) = (/'AF','BH','BD','BF','AH','AY',&
     &              'AZ','AV','AB','AQ','AI','CC','DC','AK','CA','AM',&
     &              'FA','FC','DA','DB','DD','DE','CB','AA','AP','BA',&
     &              'BB','CD','CE','AU','AS', &                            ! default 31
     &              'AE', 'AC', 'AD', 'AG', 'AJ', 'AL', 'FD', 'AN', 'AO',&
     &              'FB', 'AR', 'AT', 'BE', 'AW', 'BG', 'AX', 'BI', 'BC',&
     &              'BJ', 'BK', 'CF', 'DF', 'DG'/)                         ! additional 23
      CHARACTER (2) :: NUMBER_FORMAT

!----The following are for MOLPAK-DMACRYS
      CHARACTER (2) ::  DIRD(N_SPGR) = (/'aa','db','ah','ab','da','af',&
     &          'bh','ba','bb','av','bf','fc','dc','ai','bd','ap','aq',&
     &               'ay','fa','az','ak','am','de','dd','ca','ce','cc',&
     &               'cd','cb','au','as',&
     &             'ae', 'ac', 'ad', 'ag', 'aj', 'al', 'fd', 'an', 'ao',&
     &             'fb', 'ar', 'at', 'be', 'aw', 'bg', 'ax', 'bi', 'bc',&
     &             'bj', 'bk', 'cf', 'df', 'dg'/)
      CHARACTER (2) :: FILED(N_SPGR) = (/'AA','DB','AH','AB','DA','AF',&
     &          'BH','BA','BB','AV','BF','FC','DC','AI','BD','AP','AQ',&
     &               'AY','FA','AZ','AK','AM','DE','DD','CA','CE','CC',&
     &               'CD','CB','AU','AS', &                                ! default 31
     &              'AE', 'AC', 'AD', 'AG', 'AJ', 'AL', 'FD', 'AN', 'AO',&
     &              'FB', 'AR', 'AT', 'BE', 'AW', 'BG', 'AX', 'BI', 'BC',&
     &              'BJ', 'BK', 'CF', 'DF', 'DG'/)                         ! additional 23
      CHARACTER (500) :: V_FORMAT = '(''#!/bin/csh''/&
     &       ''# write new directories & move *.com file''//&
!    &       ''set f=('',99(A2,1X),1H)/&
!    &       ''set d=('',99(A2,1X),1H)/&
     &       ''set f=('',54(A2,1X),1H)/&
     &       ''set d=('',54(A2,1X),1H)/&
     &       ''set dir1 = ($d)''/&
     &       ''set fl  = ($f)'')'
      INTEGER :: NFORMAT
!
      OPEN (UNIT=30, FILE='runjobs', STATUS='UNKNOWN')
!----Make new directories & move *.com to new directories
      IF (KIND .EQ. 2) THEN
!        NFORMAT = INDEX (V_FORMAT, '99')
         NFORMAT = INDEX (V_FORMAT, '54')         ! 9-18-07
         WRITE (NUMBER_FORMAT,'(I2)') N_SPGR
         V_FORMAT(NFORMAT:NFORMAT+1) = NUMBER_FORMAT
         WRITE(30,V_FORMAT) (FILED(I),I=1,N_SPGR), (DIRD(I),I=1,N_SPGR)
      ELSE
         WRITE(30,V_FORMAT) (FILE(I),I=1,N_SPGR), (DIR(I),I=1,N_SPGR)
      END IF
      if (I_POTE == 3) then
         WRITE (30,20) N_SPGR+1, (COMPDID(1:NID),J=1,2)
   20    FORMAT('set n=1'/'while ($n !=',I3,')'/ &
     &          ' if (-e ',A,'_$fl[$n]_master.com) then'/ &
     &          '  if ( -d $dir1[$n] ) then'/ &
     &          '   goto skip'/'  else'/'   mkdir $dir1[$n]'/ &
     &          '  endif'/'skip:'/ &
     &          '   if (-e volume.min) mv volume.min $dir1[$n]'/ &
     &          '   mv ',A,'_$fl[$n]_master.com $dir1[$n]'/ &
     &          '   mv ','molpak_$fl[$n].input $dir1[$n]'/ &
     &          '   mv ','search_$fl[$n].data $dir1[$n]'/ &
     &          '   cp dmarel.axis $dir1[$n]'/ &
     &          '   cp cadpac.charges $dir1[$n]'/ &
     &          '   cp bondlengths $dir1[$n]'/ &
     &          '   cp will01.pots $dir1[$n]'/ &
     &          ' endif'/'@ n = ($n + 1)'/'end'/ &
     &          'echo " directories created & files moved "')
      else
         WRITE (30,10) N_SPGR+1, (COMPDID(1:NID),J=1,2)
!  10    FORMAT(48Hset r=(rshell1 rshell2 rshell3 rshell4 rshell5 \,/ &
!    &          'rshell6 rshell7 rshell8 rshell9 rshell10)',/ &
!    &          'set rs  = ($r)'/ &
!    &          'cp $hd/rshell* .'/'cp $hd/rnode.com .'/ &
!    &          'cp $hd/cp-tar.com .'/'cp $hd/tar.com .'/ &
   10    FORMAT('set n=1'/'while ($n !=',I3,')'/ &
     &          ' if (-e ',A,'_$fl[$n]_master.com) then'/ &
     &          '  if ( -d $dir1[$n] ) then'/ &
     &          '   goto skip'/'  else'/'   mkdir $dir1[$n]'/ &
     &          '  endif'/'skip:'/ &
     &          '   if (-e volume.min) mv volume.min $dir1[$n]'/ &
     &          '   mv ',A,'_$fl[$n]_master.com $dir1[$n]'/ &
     &          '   mv ','molpak_$fl[$n].input $dir1[$n]'/ &
     &          '   mv ','search_$fl[$n].data $dir1[$n]'/ &
     &          '   cp dmarel.axis $dir1[$n]'/ &
     &          '   cp cadpac.charges $dir1[$n]'/ &
     &          '   cp bondlengths $dir1[$n]'/ &
     &          '   cp pote.dat $dir1[$n]'/ &
     &          ' endif'/'@ n = ($n + 1)'/'end'/ &
     &          'echo " directories created & files moved "')
      end if
!----Submmit one job at a time
      WRITE (30,15) N_SPGR+1, COMPDID(1:NID)
15       FORMAT ('#****one job at a time****'/ &
     &       'date'/'time'/'set n=1'/'while ($n !=',I3,')'/ &
     &       ' if ( -d $dir1[$n] ) then'/'   cd $dir1[$n]'/ &
     &       '   echo "current job number: $n " '/ &
     &       '   echo "current directory : $dir1[$n]" '/ &
     &       '   if ( -e ',A,'_$fl[$n]_master.com ) then'/ &
     &       '    if (-e master.log) then'/ &
     &       '     echo "the code $fl[$n] is already running" '/ &
     &       '    else'/ &
     &       '     ./*.com > master.log '/'    endif'/ &
     &       '    date'/'    time'/'   endif'/'   cd ..'/' endif'/ &
     &       '@ n = ($n + 1)'/'end'/'exit')
!----Submit jobs w/parallel procedure
!     WRITE (30,15)
!15      FORMAT ('echo " copy rshell* and re-write rshell* "'/ &
!    &           'pwd > dire'/'set i=1'/'while ($i != 11)'/ &
!    &           ' cp $rs[$i] fort.15'/' $hd/rewrshell.exe'/ &
!    &           ' mv fort.16 $rs[$i]'/' chmod 755 $rs[$i]'/ &
!    &           ' @ i = ($i + 1)'/'end'/'rm -f fort.15'/ &
!    &           'echo " copied rshell* "'/ &
!    &           'echo " parallel process to run MOLPAK "'/ &
!    &           '#nice -10 rnode.com > log &'/'exit'/)
      RETURN
      END SUBROUTINE DMA_WRUNJOB




