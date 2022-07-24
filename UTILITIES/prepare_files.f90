      PROGRAM PREPARE_FILES
!
!----PREPARE-FILES.F90             #57                10/30/06
!
!
!    Added cubic space group Pa3                      10/30/06
!    Added space group Pcmn, permutation of Pnma (# 62). 8/28/06
!    Added space group Cm (# 8)                        8/27/06
!    Added space group Cmca (# 64)                     7/24/06
!    Added params for pyridinium N+(-H), N = 85 and    7/20/06
!        H = 86                                        7/20/06
!    Added params # 83 and 84 for the N and O atoms    7/20/06
!        in NO3- (nitrate anion).                      7/20/06
!    Added sp grps Amam and Cmcm                       5/28/06
!    Added diazo N...C-N=N-C (# 82)                    2/10/06
!    Added ester...C-C(=O)-O-C; # 79, 80, 81. Note that 2/3/06
!        anhydrides have the same C and O codes.        2/4/06
!    Added S (# 77) and N (# 78) of thiazole.          1/31/06
!    Altered the id of H on NH2 so it can be either    1/25/06
!     C-NH2 or N-NH2 (H = # 3)                         1/25/06
!    Added param # 75 and 76...n--N==C; n is 3-linked  1/20/06
!                             75 76                    1/20/06
!    Added param # 74...H linked to N1 of pyrrole-      1/3/06
!        like ring.                                     1/3/06
!    Added tetrazole subroutines to properly identify 12/31/05
!        atoms # 70, 71, 72, and 73.                  12/31/05
!    Added triazole subroutines (iq_123triazole_2 and 12/30/05
!        iq_123triazole_3) to properly id atoms       12/30/05
!        # 26 and 27.                                 12/30/05
!    Changed N_PARAMS to 74                             1/4/06
!    Added [(Csp3)2-NH].  N = # 69, H = # 68          1/15/04
!    Added sp grp C2/m; sp posn is mirror             1/12/04
!    Added an internal alkyne, C = 67                 12/31/03
!    Added the hydantoin moiety...                    12/3/03
!       --N(H)--C(O)--N(H)--C(O)--                    12/3/03
!       H = 64, O = 65, N = 66                        12/3/03
!    Added amide H (61), N (62) and O (63).           11/15/03
!       C(sp2)-NH-C(=O)                               11/15/03
!    Changed number of recognizable atom types        11/6/03
!       (N_ATYPES) from 8 to 9.  And added "I".       11/6/03
!       Also changed number of atom coefficients      11/15/03
!       (N_PARAMS) from 60 to 63.                     11/15/03
!    Changed C-C-C angle criterion in IQ_CUBANE to    8/19/03
!       handle a molecule with an angle of 93.04 deg. 8/18/03
!       Now all 3 angles must be .le. 93.1 deg for    8/19/03
!       carbon to be id'd as cubane carbon.           8/19/03
!    Added Pc and corrected pyridine N-oxide N        8/11/03
!         determination                               8/11/03
!    Added P1121 (unique axis is c)                   7/23/03
!    Altered ISB code for C2/c from 20 to 19          7/1/03
!    Altered read of G charge file to correct problem 6/19/03
!      of atom symbol comparison with WMIN file       6/19/03
!    Altered identify of pyridine-N-oxide nitrogen    6/19/03
!    Altered format to write correct chemx file with  6/1/03
!      generated H's                                  6/1/03
!    Read WMIN coefficients from a file...not stored  5/3/03
!      in the program.                                5/3/03
!    Added alcohol OH (C[sp3]-O-H) (# 58, 59)         3/25/03
!     and ability to stretch O-H                      3/25/03
!    Added separate N's for azides (# 32,56,57)       3/5/03
!    Added separate NO2 for nitramines                2/24/03
!    Added nitrate esters (-O-NO2)                    12/14/02
!    Added F or C-F                                   1/15/03
!    Modified for nitrocubanes                        11-12-2002
!    Added sulfur                                     11/2/2000
!    Added option # 12.  Like # 4 except CHEMX file
!      is created.                                    11/16/00
!    The NPOT atom ID codes were added to the wmin input file
!       atom charge lines in cols 5-9                 1/4/01
!    Updated WMIN coefs for atom types 1, 2, 4, 5, 12, 36-41    2/6/01
!    When reading a CSD file (option # 10), allowing a         4/22/06
!     transformation matrix to be input for xyz modification.  4/22/06
!
!----Program to....
!
!     (1) read a CSD FDAT file; write the 1st half of a MOLECULAR
!         FIT  input file, an input file for MOPAC and a file with
!         observed cell constants and coordinates
!     (2) convert a MNDO/ESP output file to molecular fit format and
!         create a file of the ESP charges
!     (3) read CSD FDAT file, molecular fit output file and MNDO/ESP
!         charge data file and prepare a WMIN file with the charge data
!         and so forth included
!     (4) read a CSD FDAT file, adjust nitro N-O and N-F distances,
!         adjust/idealize X-H vectors; write input files for MOPAC,
!         CHEM3D, MOLECULAR FIT (1st half) and a file with obsd cell
!         params and fractional coordinates
!     (5) convert a CSD FDAT file to CHEM3D format (this allows one to
!         view the molecule to confirm that things are alright)
!     (6) read a CHEM3D format file (from option # 4 or 5) and write
!         an input file for a G98 calcn (6-31g*, no geom optimization,
!         atom charges from ESP/CHELPG procedure)
!     (7) read G98 ...log and wmin.input files; extract ESP/CHELPG
!         charges from G98 file; replace existing charges in the
!         wmin.input file; new updated file is wmin.input2
!     (8) read a CSD FDAT file, adjust X-H and nitro N-O bond vectors
!         and prepare PCKME ...csf and ...msf files
!     (9) read an existing WMIN input file and create PCKME ...csf
!         and ...msf files
!    (10) read a CSD FDAT file, adjust nitro N-O and N-F distances, and
!         adjust/idealize X-H vectors; write an input file for MNDO
!         or G98 ESP charge calculation and a WMIN input file
!    (11) read MNDO or G98 ESP charge calculation output file and a
!         WMIN input file; create new WMIN input file that
!         incorporates the charge data & file (run.wmin) to execute
!         WMIN
!    (12) read a CSD FDAT file, adjust nitro N-O and N-F distances,
!         adjust/idealize X-H vectors; write a CHEMX output file
!
!----Files used....
!       #10 (input) -- CSD FDAT file
!       #11 (output) -- MOPAC input files (XXXXXX-mopac.dat or
!                       XXXXXX-esp.dat) or G98 input file
!                       (XXXXXX-G98.com)
!       #12 (output) -- molecular fit input file (XXXXXX-molfit.inp)
!       #13 (output) -- cartesian coordinates
!           (input) -- cart coord. for MOPAC file creation
!       #14 (output) -- CHEM3D (cc1; with serial nos.) format for options
!                       # 4 and 5; CHEMX for option # 12
!       #15 (output) -- adjusted fractional coordinates in a CHEMX file
!                       named fort.15
!       #16 (output) -- PCKME ...csf format file
!       #17 (output) -- PCKME ...msf format file
!       #18 (output) -- create MOPNAMES file for running MOPAC
!       #19 (output) -- adjusted fractional coordinates for making
!                       WMIN file with option # 10
!       #20 (output) -- name of CHEM3D file as XXXXXX.chem3d
!       #21 (output) -- name file
!          line 1 = name of molecular fit input file (XXXXXX-molfit.inp)
!          line 2 = name of mopac esp output file (XXXXXX-esp.out)
!          line 3 = name of wmin input file (the goal) (XXXXXX-wmin.input)
!       #22 -- scratch file
!          line 1 = cell parameters
!          line 2 = space group  symbol, Z, (0/1) for molecule (does not/
!                   does) occupy a special position
!          line 3 = no. of sym ops, no. of atoms
!          line 4 = symops (1 per line)
!       #23 (output) -- temporary output file for summary purposes
!       #24 (output) -- input file for G98 charge calcn,
!                       name = 631gstar_chelpg.log
!       #25 (input) -- existing wmin.input file
!       #26 (output) -- updated wmin.input file (name = wmin.input2)
!       #27 (input) -- G98 log file (name = 631gstar_chelpg.log) or
!                      MNDO MOPOUT file
!       #28 (output) -- wmin.com...file for running wmin refinment
!       #29 (output) -- chemx file from option # 12
!       #30 (input) -- file containing WMIN energy coefficients
!
      USE mod_preppot
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: N_SPGRPS = 70, N_ATYPES = 10, N_PARAMS = 86 ! 10/30/06
!
      LOGICAL :: SERIAL, CH_TYPE, ltmatrix=.false.
!CONNECT_FLAG transferred to module
!
      CHARACTER (LEN=1) :: LINE(80), tmatrix=' '         ! 4/22/06
!ID*1, AT*1, LINE4(80)*1, LINE2(80)*1, RF2(8)*1 removed, not in use
!RF2 is replaced with REFCOD, DATNAM(11)*1 is replaced with DATNM2*11
      CHARACTER (LEN=1) :: SERIAL_NUM, ATOM_ID(500)
!OUTNAM(19)*1 is replaced with OUTFIL*19
!BLANK2 is removed and data is transferred to BLANK
      CHARACTER (LEN=1) :: LLINE(1000), A_ID(80)                              ! 7/22/05
!LINE3(80)*1 removed, not in use
      CHARACTER (LEN=2) :: A_TYPE(N_ATYPES), ATOMID, ANAME2
      CHARACTER (LEN=3) :: NSPGRP, END
!  --------------------------------------   END is a fortran function
!A_LABEL*3 removed, not in use
      CHARACTER (LEN=4) :: A_NAME
!A_NAME set on 1076 but never used
      CHARACTER (LEN=5) :: TEMP, TLABEL(500), ANAME3
!NAME(500)*5 removed, not in use
      CHARACTER (LEN=8) :: REFCOD, RF, SPGRP, SG_TYPES(N_SPGRPS)
      CHARACTER (LEN=9) :: C_WORD
      CHARACTER (LEN=10) :: CDATE, CTIME
!CZONE*10 removed, not in use
      CHARACTER (LEN=12) :: WMIN_REF_CODES(2,N_SPGRPS)
      CHARACTER (LEN=27) :: LINE27
      CHARACTER (LEN=30) :: COEF_FILE, COEF_DESCRIPTION(N_PARAMS)                ! 5/3/03
      CHARACTER (LEN=31) :: LINE31
!LINE35*35 removed, not in use
      CHARACTER (LEN=36) :: LINE36
      CHARACTER (LEN=54) :: LINE54
      CHARACTER (LEN=60) :: BLANK, TITLE, INFILE
!      CHARACTER (LEN=68) :: CSF_TITLE
      CHARACTER (LEN=72) :: LINE72, MF_TITLE, SYM_OP
      CHARACTER (LEN=80) :: WMIN_COM(35), LINE80
      CHARACTER (LEN=81) :: LINE81                                      ! 5/3/03
      CHARACTER (LEN=1000) :: LLINE500                                  ! 7/22/05
!
      INTEGER (KIND=2) :: IP(6)
!IS_OUT(500), CDATE_TIME(8) removed, not in use
      INTEGER :: ISYM(3,3,48)                                           ! 7/22/05
      INTEGER :: ICELL(6), NVALUE(37), IXYZ(3,500)
      INTEGER :: IJK_WMIN(500), LATTICE(N_SPGRPS), ISB(N_SPGRPS)
      INTEGER :: IPCKME_POT(N_ATYPES)
!
      REAL :: TRANS(3,48), SYM(3,3), R2(3,3), R3(3,3)                   ! 7/22/05
      REAL :: DEN(2)
      REAL :: CHARGE(500), SUMXYZ(3), W_XYZ(3)
      REAL :: W_TYPE(N_ATYPES), ABC_COEFS(3,N_PARAMS)          ! 5/3/03
      REAL :: AT_NUMBER(N_ATYPES)
!
!      DIMENSION ISYM(3,3,16), TRANS(3,16), SYM(3,3), R2(3,3), R3(3,3)
!
      INTEGER :: I, J, K, L, LM, M, N, NA, NL, NS, NZ
      INTEGER :: I1, I2, I3, I4, IZ, I_ONE, I_WHAT
      INTEGER :: I601, I_CENTER, IEND, IPASS, I_PLUS, ISTART
      INTEGER :: K_TYPE
      INTEGER :: LATSYM
      INTEGER :: NC, NKA, N_OPS, NPOT, NSYM_OPS, N_WMIN
      INTEGER :: MULT, MMULT
!NAT, NSAT removed, not in use
      REAL    :: AC, BC, CC
      REAL    :: OLDQ, SCALED_Q, WT_MOL, tr(3,3), tt(3), new_xyz(3)
!
      EQUIVALENCE (LINE, LINE80, LINE31, LINE27)
!LINE2 is taken out of equivalence statement
!      EQUIVALENCE (MF_TITLE, SYM_OP)
!MF_TITLE read 658; line689 read st and line973 write st, no more
!SYM_OP   ; line1357 read st and line 1359 write st, no more
!      EQUIVALENCE (LINE3, A_ID)  ! LINE3 is not in use
!
!      EQUIVALENCE (MTITLE, CSF_TITLE)
! MTITLE is common, CSF_TITLE ; 578 r, 579 w - 994 r, 996 w
!CSF_TITLE*68 is changed to MTITLE*72
      EQUIVALENCE (LINE36, LINE72, LINE54)
!LINE35 is taken out of equivalence statement
!      EQUIVALENCE (LABEL, NAME)
!      EQUIVALENCE (RF2, REFCOD), (OUTNAM, OUTFIL),   &
!     &            (DATNAM, DATNM2)
!(RF2, REFCOD) removed, RF2 is replaced with REFCOD
      EQUIVALENCE (LINE(2), RF), (LLINE, LLINE500)
!LLINE(I) can be replaced w/LLINE500(I:I)
!      EQUIVALENCE (IXYZ, OR_XYZ)
!OR_XYZ is real, IXYZ is int ...!
!      EQUIVALENCE (NVALUE(11), NAT), (NVALUE(12), NSAT),   &
!     &            (NVALUE(9), N_OPS)
      EQUIVALENCE (NVALUE(9), N_OPS)
!
      DATA TLABEL /500*'     '/

!      EQUIVALENCE (BLANK, BLANK2)
!BLANK2 is removed and data is transferred to BLANK
!      DATA BLANK/ 60*' ' /
      DATA BLANK/ '                                        &
     &                    ' /
      DATA R3 /1.0, 0.0, 0.0,  0.0, 0.0, 1.0,  0.0, -1.0, 0.0/
      DATA DEN /10000., 100000./
      DATA AT_NUMBER / 1,   6,   7,   8,   9,   35, 16, 5, 53, 17/,  &  ! 6/10/05
     &    IPCKME_POT / 2,   4,   0,   0,   0,    0,  0, 0,  0,  0/      ! 6/10/05
!
      DATA WMIN_COM /'#!/bin/csh', 'pwd', '#', 'time', 'date', '#',  &  ! 5-26-05
       'setenv hd1 /export/software/predictions-july07', &
     & '# **************** Run WMIN *************************',      &
     & 'cp wmin.input2 fort.15', '#',                                &
     & '#              WMIN input file on unit # 15',                &
     & '#                  ...sum file on unit # 21',                &
     & '#                 ...save file on unit # 22',                &
     & '#        cell for NIST*LATTICE on unit # 11', '#',           &
     & '$hd1/WMIN/wmin_100.exe',                                     &
     & ' cat fort.22 >> wmin.save-total',                            &
     & ' mv fort.16 wmin.step-output',                               &
     & ' mv fort.99 step.list ',                                     &
     & ' if ( -e fort.22 )  mv fort.22 wmin.save',                   &
     & '  $hd1/UTILITIES/rewrite-wmininp-100.exe',                   &
     & '   if ( -e wminlsq.inp ) then',                              &
     & '    echo "refinement with LSQ"',                             &
     & '    mv wminlsq.inp fort.15',                                 &
     & '$hd1/WMIN/wmin_100.exe',                                     &
     & '    cat fort.22 >> wmin.save-total',                         &
     & '    mv fort.16  wmin.lsq-output',                            &
     & '    mv fort.99 lsq.list ',                                   &
     & '   endif',                                                   &
     & ' cat fort.11 >> wmin.cell',                                  &
     & ' cp wmin.cell fort.15',                                      &
     & ' $hd1/UTILITIES/nbslattice.exe',                             &
     & ' cat fort.11 >> wmin.cell',                                  &
     & '# ************End of cell reduction********************',    &
     & '#  FINISHED'/
!
!----LATTICE = 1 for triclinic,
!              2 for monoclinic,
!              3 for orthorhombic,
!              4 for tetragonal with 4/mmm Laue symmetry,
!              5 for tetragonal with 4/m Laue symmetry
!              6 for rhombohedral
!              7 for cubic
!              8 for trigonal
!     ISB is the WMIN subroutine switch;
!     N_SPGRPS = no. of space groups for which data are stored.
!     WMIN_REF_CODES = WMIN parameter refinement codes (1 = normal,
!                      2 = special position)
!----New space groups require these arrays to be updated.
      DATA SG_TYPES /'Pbn21   ', 'P-1     ', 'P21/n   ', 'P21/a   ', &
     &               'P212121 ', 'P21/c   ', 'Pbca    ', 'P21     ', &
     &               'C2/c    ', 'C2      ', 'I2/c    ', 'Pmn21   ', &
     &               'Pn21m   ', 'Pcab    ', 'Ccca    ', 'Pna21   ', &
     &               'P21nb   ', 'A2/a    ', 'Pbcn    ', 'P41212  ', &
     &               'P1      ', 'P-421c  ', 'Cc      ', 'P21/m   ', &
     &               'I41/a   ', 'Fdd2    ', 'P1121/a ', 'Pc21n   ', &
     &               'C22/c   ', 'Pa      ', 'C2cb    ', 'Aba2    ', &
     &               'P21cn   ', 'P21212  ', 'Pnma    ', 'Pca21   ', &
     &               'Pn      ', 'Cmc21   ', 'Pn21a   ', 'Pnam    ', &
     &               'Pbc21   ', 'P61     ', 'P-421c  ', 'Pcnb    ', &
     &               'P21ab   ', 'P1121/b ', 'A112/a  ', 'B112/b  ', &
     &               'Ibam    ', 'C2cb    ', 'Pccn    ', 'Iba2    ', &
     &               'R-3     ', 'P41     ', 'P22121  ', 'P1121   ', &
     &               'Pc      ', 'I2/a    ', 'Pbcm    ', 'C2/m    ', &
     &               'R-3r    ', 'I-43m   ', 'P32     ', 'P42     ', &
     &               'Amam    ', 'Cmcm    ', 'Cmca    ', 'Cm      ', &
                     'Pcmn    ', 'Pa3     '/                  ! 10/30/06
      DATA LATTICE  /    3,           1,         2,          2,   &
     &                   3,           2,         3,          2,   &
     &                   2,           2,         2,          3,   &
     &                   3,           3,         3,          3,   &
     &                   3,           2,         3,          4,   &
     &                   1,           4,         2,          2,   &
     &                   5,           3,         2,          3,   &
     &                   2,           2,         3,          3,   &
     &                   3,           3,         3,          3,   &
     &                   2,           3,         3,          3,   &
     &                   3,           1,         1,          3,   &
     &                   3,           2,         2,          2,   &
     &                   3,           3,         3,          3,   &
     &                   1,           1,         3,          3,   &
     &                   2,           2,         3,          2,   & ! 1/12/04
     &                   5,           6,         8,          5,   &
     &                   3,           3,         3,          2,   &
                         3,           7/                          ! 10/30/06
!----ISB is the ISB variable for structure refinement...all set = 19
      DATA ISB  / n_spgrps*19/                                      ! 6/9/06
      DATA WMIN_REF_CODES /  '111000111110', '111000000000',  & ! Pbn21
     &                       '111111111111', '111111111000',  & ! P-1; sp posn = center
     &                       '111010111111', '111010111000',  & ! P21/n; sp posn = center
     &                       '111010111111', '111010111000',  & ! P21/a; sp posn = center
     &                       '111000111111', '111000111111',  & ! P212121
     &                       '111010111111', '111010111000',  & ! P21/c; sp posn = center
     &                       '111000111111', '111000111000',  & ! Pbca; sp posn = center
     &                       '111010111101', '111010111101',  & ! P21
     &                       '111010111111', '111010111000',  & ! C2/c; sp posn = center
     &                       '111010111101', '111010111101',  & ! C2
     &                       '111010111111', '111010111000',  & ! I2/c; sp posn = center
     &                       '111000111110', '111000100010',  & ! Pmn21; sp posn = m
     &                       '111000111101', '111000001100',  & ! Pn21m; sp posn = m
     &                       '111000111111', '111000111000',  & ! Pcab; sp posn = center
     &                       '111000111111', '111000000000',  & ! Ccca; sp posn = 222
     &                       '111000111110', '111000111110',  & ! Pna21
     &                       '111000111011', '111000111011',  & ! P21nb
     &                       '111010111111', '111010111000',  & ! A2/a; sp posn = center
!    &                       '111000111111', '111000111000',  & ! Pbcn; spposn = center
     &                       '111000111111', '111000010010',  & ! Pbcn; sp posn = 2-fold
     &                       '101000111111', '101000001001',  & ! P41212; sp posn = 2-fold
     &                       '111111111000', '111111111000',  & ! P1
     &                       '101000111111', '101000001000',  & ! P-421c; sp posn = -4
     &                       '111010111010', '111010111010',  & ! Cc
     &                       '111010111111', '111010010101',  & ! P21/m; sp posn = m
     &                       '101000111111', '101000111111',  & ! I41/a
     &                       '111000111110', '111000001000',  & ! Fdd2; sp posn = 2-fold
     &                       '111001111111', '111001111000',  & ! P1121/a [c = unique axis]
     &                       '111000111111', '111000111111',  & ! Pc21n
     &                       '111010010010', '111010010010',  & ! C22/c [C2/c, sp posn = 2-fold]
     &                       '111010111111', '111010111111',  & ! Pa
     &                       '111000111111', '111000100100',  & ! C2cb [sp posn = 2-fold]
     &                       '111000111111', '111000001001',  & ! Aba2 [sp posn = 2-fold]
     &                       '111000111011', '111000111011',  & ! P21cn
     &                       '111000111111', '111000001001',  & ! P21212 [sp posn = 2-fold]
     &                       '111000111111', '111000010101',  & ! Pnma [sp posn = m]
     &                       '111000111110', '111000111110',  & ! Pca21
     &                       '111010111111', '111010111111',  & ! Pn
     &                       '111000111110', '111000100010',  & ! Cmc21 [sp posn = m]
     &                       '111000111101', '111000111101',  & ! Pn21a
     &                       '111000111111', '111000001110',  & ! Pnam  [sp posn = m]
     &                       '111000111110', '111000111110',  & ! Pbc21
     &                       '101000001000', '101000001000',  & ! P61
     &                       '101000100000', '101000000000',  & ! P-421c [sp posn = -4]
     &                       '111000111111', '111000001001',  & ! Pcnb  [sp posn = 2-fold]
     &                       '111000111011', '111000111011',  & ! P21ab
     &                       '111010111111', '111010111000',  & ! P1121/b [c = unique axis]
     &                       '111010111111', '111010111000',  & ! A112/a
     &                       '111010111111', '111010111000',  & ! B112/a
     &                       '111000111111', '111000001110',  & ! Ibam [sp posn = m]
     &                       '111000111111', '111000100100',  & ! C2cb [sp posn = 2-fold]
     &                       '111000111111', '111000111111',  & ! Pccn
     &                       '111000111110', '111000001000',  & ! Iba2 [sp posn = 2-fold]
     &                       '100100000000', '100100000000',  & ! R-3
     &                       '101000110111', '101000110111',  & ! P41
     &                       '111000111111', '111000100100',  & ! P22121 [sp posn = 2-fold]
     &                       '111001111111', '111001111111',  & ! P1121        ! 7/23/03
     &                       '111010111110', '111010111110',  & ! Pc    ! 8/13/03
     &                       '111010111111', '111010111111',  & ! I2/a  ! 8/13/03
     &                       '111000111111', '111000001110',  & ! Pbcm [sp posn = m] ! 9/20/03
     &                       '111010111111', '111010010101',  & ! C2/m [sp posn = m] ! 1/12/04
     &                       '100100111111', '100100100000',  & ! R-3r [sp posn = ?] ! 7/22/05
                             '100000111111', '100000100000',  & ! I43m       "            "
                             '101000111111', '101000111111',  & ! P32
                             '101000111110', '101000111110',  & ! P42                  7/25/05
                             '111000111111', '111000000010',  & ! Amam [sp posn = mm]
                             '111000111111', '111000000010',  & ! Cmcm [sp posn = mm]  5/30/06
                             '111000111111', '111000100101',  & ! Cmca [sp posn = m]   7/24/06
                             '111010111111', '111010010101',  & ! Cm   [sp posn = m]   8/26/06
                             '111000111111', '111000010101',  & ! Pcmn [sp posn = m]   8/28/06
                             '100000111111', '100000100100'/    ! P32  [sp posn = ??] 10/30/06
!
      DATA REFCOD /'        '/
      DATA A_TYPE /'C ',   'N ',    'O ',    'H ',   'F ',   &
     &             'BR',   'B ',    'S ',    'I ',   'CL'/           ! 6/10/05
      DATA W_TYPE /12.011, 14.0067, 15.9994, 1.0079, 18.9984,  &
     &             79.904, 10.811,  32.066,126.9045,35.4527/         ! 6/10/05
!
! End of data decleration -------------------------------------
!
      MOPAC_TYPE = 2
      CONNECT_FLAG = .false.
      REWIND 23
!
      call date_and_time (cdate, ctime)
      PRINT 10, CDATE(5:6), CDATE(7:8), CDATE(1:2), CDATE(3:4),          &
     &          CTIME(1:2), CTIME(3:4), CTIME(5:6)
10    FORMAT (' PREPARE-FILES (version 7/20/06) for energy parameter',  &
     &        ' refinement (',2(A2,'-'),2A2,' at ',2(A2,':'),A2,')',     &
     &        '....'/)
      WRITE (23,10) CDATE(5:6), CDATE(7:8), CDATE(1:2), CDATE(3:4),      &
     &              CTIME(1:2), CTIME(3:4), CTIME(5:6)
      PRINT 12
12    FORMAT (' 1 = read a CSD FDAT file; write the 1st half of a',      &
     &        ' MOLECULAR FIT  input file,'/                             &
     &        '     an input file for MOPAC and a file with',            &
     &        ' observed cell constants and coordinates'/                &
     &        ' 2 = read a MNDO/ESP output file; put the optimized',     &
     &        ' MOPAC coordinates'/                                      &
     &        '     into the 2nd half of a',                             &
     &        ' MOLECULAR FIT  input file and create a file of',         &
     &        ' ESP charges'/                                            &
     &        ' 3 = read MOLECULAR FIT output, ESP charge data',         &
     &        ' and a scratch file;'/                                    &
     &        '     prepare a complete WMIN file with the charge',       &
     &        ' data and so forth included'/                             &
     &        ' 4 = read a CSD FDAT file, adjust nitro N-O distances',   &
     &        ' and adjust/idealize X-H vectors;'/                       &
     &        '     write input files for MOPAC, CHEM3D, MOLECULAR FIT', &
     &        ' (1st half)'/                                             &
     &        '     and a file with obsd cell params and fractional',    &
     &        ' coordinates'/                                            &
     &        ' 5 = read a CSD FDAT file; write an input file for',      &
     &        ' CHEM3D'/                                                 &
     &        '     (cc1 format with serial numbers)'/                   &
     &        ' 6 = read a CHEM3D file; write input file for G98',       &
     &        ' charge calculation'/                                     &
     &        '     (6-31g*, no geom optim, ESP/CHELPG charge calcn)'/   &
     &        ' 7 = read G98 ...log and complete wmin.input files;',     &
     &        ' extract G98 charges and'/                                &
     &        '     replace old charges with G98 values'/                &
     &        ' 8 = read a CSD FDAT file, adjust X-H and nitro N-O',     &
     &        ' bond vectors'/                                           &
     &        '     and prepare PCKME ...csf and ...msf files'/          &
     &        ' 9 = read an existing WMIN input file and create',        &
     &        ' PCKME ...csf and ...msf files'/                          &
     &        '10 = read a CSD FDAT file, adjust nitro N-O and N-F ',    &
     &        'distances and adjust/idealize X-H vectors;'/              &
     &        '     write an input file for MNDO or G98 ESP charge',     &
     &        ' calculation and a WMIN input file;'/                     &
     &        '     allows xyz transformation matrix to be applied'/     &
     &        '11 = read MNDO or G98 ESP charge calculation output',     &
     &        ' files and a WMIN input file;'/                           &
     &        '     create a new WMIN input file that incorporates the', &
     &        ' charge data & file (run.wmin) to execute WMIN;'/         &
     &        '12 = read a CSD FDAT file, adjust nitro N-O and N-F',     &
     &             ' distances,'/                                        &
     &       '     adjust/idealize X-H vectors and write a CHEMX file'/)
!
      PRINT 1001
1001  FORMAT (/' Conversion required: ', $)
      READ *, I_WHAT
      PRINT '()'
      WRITE (23,2001) I_WHAT
2001  FORMAT (' Conversion performed:',I2)
!
!        ----  1     2     3     4     5     6     7     8  -----
      GO TO (1002, 1010, 2000, 1002, 1003, 1005, 3007, 1002,   &
!        ----  9    10    11    12 ----
     &       3070, 1002, 6030, 1002), I_WHAT
!
!----Read existing WMIN file and create ...msf and ...mcf files
!     for PCKME
!
3070  PRINT 1014
1014  FORMAT (' Give the name of the WMIN input file',   &
     &        ' [wmin.input]: ',$)
      INFILE = BLANK
      READ 2, INFILE
      IF (INFILE .EQ. BLANK) INFILE = 'wmin.input'
!----Open the WMIN input file
      OPEN (UNIT=25,FILE=INFILE,STATUS='OLD')
      WRITE (23,2015) INFILE
2015  FORMAT (' WMIN input file: ',A60)
      GO TO 77
!
!----Use # 22 as a scratch file to hold data for communicating
!     between various calculations
1002  REWIND 22
!
!********* preliminary input for reading CSD file ************
!
!----Request name of CSD FDAT input file
1003  PRINT 11
11    FORMAT (' Give the name of the FDAT input file: ',$)
      READ 2, INFILE
2     FORMAT (A60)
!----OPEN THE INPUT FILE
      OPEN (UNIT=10,FILE=INFILE,STATUS='OLD')
      WRITE (23,1007) INFILE
1007  FORMAT (' FDAT input file: ',A60)
!
!----Read REFCOD to identify the compound to be transformed.
77    PRINT 3
      OUTFIL = '            '
3     FORMAT (' Give CSD ref code: ',$)
      READ 4, REFCOD
4     FORMAT (A8)
      IF (REFCOD .EQ. 'END' .OR. REFCOD .EQ. '        ') STOP
!
!----Open the output file...but first figure out how many characters
!     are in the REFCOD name
      WRITE (23,1004) REFCOD
1004  FORMAT (' REFCOD: ',A8)
      DO 20 I=1,8
         NCC = I                ! # of non-blank chars in name
         IF (REFCOD(I:I) .EQ. ' ') GO TO 29
20    CONTINUE
29    OUTFIL(1:I) = REFCOD(1:I)
      NC = NCC
!----Request xyz transformation matrix
      print 48
48    format ('x,y,z transformation matrix as r11,r12,r13,t1,r21..', &
              ' Y/N, [ ] = no: ',$)
      read (*,'(a1)') tmatrix
      if (tmatrix .eq. 'N' .or. tmatrix .eq. 'n' .or. tmatrix .eq. ' ') go to 23
         ltmatrix = .true.
         print 49
49          format ('   matrix elements: ',$)
         read *, ((tr(i,j),j=1,3),tt(i),i=1,3)
!
23    GO TO (777, 777, 777, 777, 42, 777, 777, 778, 778,   &
     &      1077, 9988, 774), I_WHAT
9988  WRITE (6,9987) I_WHAT
9987  FORMAT ('**Inappropriate transfer based on optin #',I3/  &
     &        '**Program stop**')
!
!************ preliminary input for option # 10 *****************
!
1077  PRINT 1079
1079  FORMAT (' Type of energy coefficients for WMIN input file,'/  &
     &        '   1 = no coulombic term, 2 = MNDO/ESP,',        &     ! 5/3/03
     &        ' 3 = Gaussian/ESP: ',$)                              ! 5/3/03
      READ *, I_CH_TYPE
      WRITE (21,*)                 ! line #1 on unit # 21 is blank    5/3/03
      IF (I_CH_TYPE .EQ. 0) I_CH_TYPE = 3                           ! 5/3/03
      IF (I_CH_TYPE .LE. 1) THEN                                    ! 5/3/03
         WRITE (21,*)              ! line # 2 of unit 21 is blank   ! 5/3/03
         GO TO 1083
      ENDIF
      IF (I_CH_TYPE .EQ. 2) GO TO 212         ! MNDO, create file name
!----Read WMIN coefficients from file # 30; first get file name     ! 5/3/03
      WRITE (*,1074)                                                ! 5/3/03
1074  FORMAT (' Name of file with WMIN energy coefficients: ',$)    ! 5/3/03
      READ (*,1075) COEF_FILE                                       ! 5/3/03
      DO I=1,N_PARAMS       ! initialize coef array                   5/3/03
         DO J=1,3                                                   ! 5/3/03
            ABC_COEFS(J,I) = 0.0                                    ! 5/3/03
         ENDDO                                                      ! 5/3/03
         COEF_DESCRIPTION(I) = BLANK(1:30)                          ! 5/3/03
      ENDDO                                                         ! 5/3/03
      OPEN (UNIT = 30, FILE = COEF_FILE, STATUS = 'OLD')            ! 5/3/03
         READ (30,1075) LLINE500(1:81)                              ! 5/3/03
1075     FORMAT (A)                                                 ! 5/3/03
1076     READ (30,1075,END=977) LINE81                              ! 5/3/03
         IF (LINE81(1:3) .EQ. 'END') GO TO 977                      ! 5/3/03
         READ (LINE81,970) N, (ABC_COEFS(I,N),I=1,3),     &          ! 5/3/03
     &                      COEF_DESCRIPTION(N)                     ! 5/3/03
970      FORMAT (I3,3X,3F15.4,A)                                    ! 5/3/03
         GO TO 1076                                                 ! 5/3/03
977    OPEN (UNIT = 31, FILE = 'wmin_coefficients', STATUS = 'UNKNOWN') ! 5/3/03
         WRITE(31,1075) LLINE500(1:81)                              ! 5/3/03
         WRITE(31,971) (I, (ABC_COEFS(J,I),J=1,3),          &        ! 5/3/03
     &                  COEF_DESCRIPTION(I),I=1,N_PARAMS)           ! 5/3/03
971      FORMAT (I3,F13.4,F16.4,F14.2,2X,A)                         ! 5/3/03
      CLOSE (UNIT=31)                                               ! 5/3/03
      PRINT 1081
1081  FORMAT (' A G03 input file will be created with the',   &
     &        ' name of 631gstar_chelpg.com')
      OPEN (UNIT=11, FILE='631gstar_chelpg.com',   &    ! open G98 input file
     &      STATUS='UNKNOWN', FORM='FORMATTED')     ! on unit # 11
      WRITE (11, 1082) CDATE(5:6), CDATE(7:8), CDATE(1:2),  &
     &                 CDATE(3:4), (OUTFIL(1:NC))
1082  FORMAT ('$RunGauss'/                                           &  ! G98 input file
     &        '%chk=631gstar_chelpg.chk'/                            &
     &        '#p 6-31G* scf=direct pop=CHELPG geom=coord'//         &
     &        1X,A2,'-',A2,'-',2A2,', ',A8,', charges from 6-31g*',  &
     &        '  ESP/CHELPG,',                                       &
     &        ' cart coord input'//'   0   1')
      WRITE (21,2083)                  ! place '631gstar_chelpg.com'  	
2083  FORMAT ('631gstar_chelpg.log')    ! on line # 2 of file # 21
      GO TO 1083                       ! go to make WMIN input file name
!
!**********************************************************************
!
!----Option # 8...create PCKME ...CSF & MCF files
!
778   DATNM2 = '.csf      '           ! ...csf file
!         DO J=1,10
!            OUTFIL(J+NC:J+NC) = DATNM2(J:J)
!         ENDDO
!lines replaced with
            OUTFIL(NC:NC+10) = DATNM2
         PRINT 785, OUTFIL
785      FORMAT (' A PCKME ...csf file will be created with',   &
     &           ' the name of ',A18)
       WRITE (23,785) OUTFIL
       OPEN (UNIT=16,FILE=OUTFIL,STATUS='UNKNOWN',  & ! ...csf on 16
     &   FORM='FORMATTED')
      DATNM2 = '.mcf      '           ! ...mcf file
!         DO J=1,10
!            OUTFIL(J+NC:J+NC) = DATNM2(J:J)
!         ENDDO
!lines replaced with
            OUTFIL(NC:NC+10) = DATNM2
         PRINT 786, OUTFIL
786      FORMAT (' A PCKME ...mcf file will be created with',   &
     &           ' the name of ',A18)
       WRITE (23,786) OUTFIL
       OPEN (UNIT=17,FILE=OUTFIL,STATUS='UNKNOWN',   &   ! ...mcf on 17
     &   FORM='FORMATTED')
       IF (I_WHAT .NE. 9) GO TO 42
!
!----Prepare PCKME ...csf file; input from existing WMIN file
!
      PRINT 810
      READ (5, 815) MTITLE                   ! 68 chars max
      WRITE (16,1029) MTITLE                  ! csf card # 1...
1029  FORMAT (A68,'  -1')
!----Read WMIN input file, looking for proper REFCOD
2030  READ (25,1032,END=1040) RF, LINE54
1032     FORMAT (3X,A8,A54)
         IF (RF .EQ. REFCOD) GO TO 1034
         GO TO 2030
1040     PRINT 1042, INFILE
1042     FORMAT (' Cannot locate ref code ',A8,' on file ',A60)
         STOP
1034  WRITE (16,1036) REFCOD, LINE54
1036  FORMAT (A8,1X,A54,8X,'0')        ! 2nd title line from WMIN title info
      PRINT 1037, LINE54
1037  FORMAT (' WMIN title info: ',A54)
      READ (25,1046) NATOMS, NKA, NSYM_OPS  ! NATOMS includes XTRA
1046     FORMAT (3I3)
         NATOMS = NATOMS - 1
      WRITE (16,822) NSYM_OPS, NATOMS           ! csf card # 2
      WRITE (16,824)                            ! csf cards # 3-5
      DO I=1,4             ! skip next 4
         READ (25,*)       ! lines on WMIN
      ENDDO                ! input file
      READ (25,1048) CELL      ! cell parameters
1048  FORMAT (6F9.4)           ! angles are cosines
      DO I=4,6
         CELL(I) = ACOS(CELL(I))*57.29578  ! cosines to degrees
      ENDDO
!----Easier to ask for LATSYM then figure it out
      PRINT 1057
1057    FORMAT (' Give csf file LATSYM code (1 = tric, mono,',   &
     &          ' orth; 2 = tetr, hexa;'/                        &
     &          '   3 = trig; 4 = cubic: ',$)
      READ (5,*) LATSYM
      WRITE (16,826) CELL, LATSYM               ! csf card # 6
      READ (25,*)    ! skip WMIN esd line
      DO K=1,NSYM_OPS
         READ (25,3111) (TRANS(J,1),(SYM(I,J),I=1,3),J=1,3) ! sym
3111     FORMAT (3(4X,F11.2,3F3.0))                         ! ops
         WRITE (16,1210) (TRANS(J,1),J=1,3),       &     ! csf cards # 7
     &                   ((SYM(I,J),I=1,3),J=1,3)       ! ...symmetry
      ENDDO
!----Finished making PCKME ...csf file; now make ...mcf file
      GO TO 839
!
!----Option # 1...create molecular fit & mopac input files
!
!----Mopac input file on # 11 'XXXXX-mopac.dat' or 'XXXXX-esp.dat'
!                              ---------------
777   PRINT 201
201      FORMAT (' Type of MOPAC file, ',   &
     &           '1 = AM1 optimization, 2 = MNDO/ESP: ',$)
       READ 205, MOPAC_TYPE
205       FORMAT (I1)
       GO TO (211, 212), MOPAC_TYPE
211      DATNM2 = '-mopac.dat'           ! AM1 optimization
!         DO 25 J=1,10
!25          OUTFIL(J+NC:J+NC) = DATNM2(J:J)
!       GO TO 220
!lines replaced with ..
            OUTFIL(NC:NC+10) = DATNM2
212      DATNM2 = '-esp.dat'           ! MNDO/ESP calcns
!         DO 215 J=1,8
!215         OUTFIL(J+NC:J+NC) = DATNM2(J:)
 !lines replaced with ..
            OUTFIL(NC:NC+8) = DATNM2
220    PRINT 30, OUTFIL
30       FORMAT (' A MOPAC input file will be created with',  &
     &           ' the name of ',A18)
!----Create the MOPNAMES file for use with MOPAC ver 7
       OPEN (UNIT=18, FILE='MOPNAMES', STATUS='UNKNOWN')
       WRITE (18,31) OUTFIL            ! place name of MOPAC input
31     FORMAT ('./',A18)               ! file into MOPNAMES file
       CLOSE (UNIT=18)
!----Open MOPAC input file on unit # 11
       WRITE (23,30) OUTFIL
       OPEN (UNIT=11,FILE=OUTFIL,STATUS='UNKNOWN',  &  ! mopac on 11
     &   FORM='FORMATTED')
!----Request mopac title
          print 6
6         FORMAT ('   Mopac title: ',$)
          read (5,34) MTITLE
34        FORMAT(A72)
!
774   IF (I_WHAT .EQ. 10) GO TO 1085
      IF (I_WHAT .NE. 12) GO TO 42
!
!----If option # 12, create CHEMX file named XXXXXX.chemx
!                                            ------------
      DATNM2 = '.chemx'
!      DO J=1,6
!         OUTFIL(J+NC:J+NC) = DATNM2(J:J)
!      ENDDO
!lines replaced with ..
            OUTFIL(NC:NC+6) = DATNM2
      GO TO 38
!
!----If option # 4, 5 or 8, create name for CHEM3D (cc1) file named
!     XXXXXX.cc1
!     -------------
!
42    IF (I_WHAT .EQ. 5 .OR. I_WHAT .EQ. 4 .OR. I_WHAT .EQ. 8) THEN
      CALL OPEN_CHEM3D
      ENDIF
      IF (I_WHAT .EQ. 5 .OR. I_WHAT .EQ. 8) GO TO 38  ! skip to CSD input
!
!----Molecular fit input file on # 12 as 'XXXXX-molfit.inp'
!                                         ----------------
         DATNM2 = '-molfit.inp'
!         DO 26 J=1,11
!26         OUTFIL(J+NC:J+NC) = DATNM2(J:)
!lines replaced with ..
            OUTFIL(NC:NC+10) = DATNM2
         OPEN (UNIT=12,FILE=OUTFIL,STATUS='UNKNOWN', & ! molecular fit on
     &         FORM='FORMATTED')                     ! unit # 12
         PRINT 60, OUTFIL
60         FORMAT (' An initial molecular fit input file will',   &
     &             ' be created as ',A19)
         WRITE (23,60) OUTFIL
!----Request molecular fit title
          print 61
61        FORMAT ('   Molecular-fit title: ',$)
          read (5,34) MF_TITLE
!----Put name of molecular fit input file into unit # 21
24       WRITE (21,28) OUTFIL      ! If I_WHAT = 4, write unit 21 to
28	 FORMAT (A19)              ! keep spacing intact
!
!----Create name for mopac esp output file...XXXXXX-esp.out
1085     DATNM2 = '-esp.MOPOUT'
!         DO 27 J=1,11
!27         OUTFIL(J+NC:J+NC) = DATNM2(J:J)
!lines replaced with ..
            OUTFIL(NC:NC+11) = DATNM2
!----Put name of esp output file on line # 2 of unit # 21
	 WRITE (21,28) OUTFIL
!
!----Create name for WMIN input file
1083     DATNM2 = '-wmin.input'
!         DO 227 J=1,11
!227        OUTFIL(J+NC:J+NC) = DATNM2(J:J)
!lines replaced with ..
            OUTFIL(NC:NC+11) = DATNM2
!----Put name of wmin input file on line # 3 of unit # 21
         WRITE (21,28) OUTFIL
         CLOSE (UNIT = 21)     ! close name file
!
!----Begin reading RDAT file...search for REFCOD
38    READ (10,41,END=1000) LINE80
41    FORMAT (A80)
40    FORMAT (80A1)
!----Look for the # as the start of a compound block
      IF (LINE(1) .NE.'#') GO TO 38
!----Found one...now check the ID
      IF (RF .NE. REFCOD) GO TO 38
!----Found the one we are looking for...DECODE the line
      READ (LINE80,50) NVALUE
50    FORMAT (9X,2I1,I6,6X,11I3,22I1,I2)
      I_CENTER = NVALUE(18)
!----If there's no cell card...quit
      IF (NVALUE(15) .EQ. 0 .OR. I_CENTER .EQ. 0) THEN
         PRINT 55
55       FORMAT (' Either no cell or sym op information',   &
     &           '...transformation aborted')
         GO TO 77
      ENDIF
!----Read line # 2
      READ (10,41,END=1000) LINE80
!----Decode to get cell params, space group and
      READ (LINE80,59) ICELL, IP, NSPGRP, SPGRP, IZ
59    FORMAT (6I6,6I1,T61,A3,A8,I3)
!----Construct cell parameters
      DO 63 I=1,6
63    CELL(I) = FLOAT(ICELL(I))/(10**IP(I))
      PRINT 65, CELL, NSPGRP, SPGRP, IZ
      WRITE (23,65) CELL, NSPGRP, SPGRP, IZ
65    FORMAT (' Information from the CSD file...'/                 &
     &        '   Cell parameters...'/                             &
     & 5X,'a =',F7.3,', b =',F7.3,', c =',F7.3, ' Angs'/           &
     & 5X,'alpha =',F7.2,', beta =',F7.2,' gamma =',F7.2,' degs'/  &
     & '   Space group...','No.',A3,', ',A8,', Z =',I3)
      IF (I_CENTER .NE. 1) GO TO 64
         PRINT 67
67       FORMAT ('   Center of symmetry at the origin')
         WRITE (23,67)
!----Does molecule lie on special position...is NSAT = 0?
!64    IF (NSAT .NE. 0) THEN ! changed b/c equivalence
64    IF (NVALUE(12) .NE. 0) THEN
         ISPECIAL = 1                       ! yes
         PRINT 66
         WRITE (23,66)
66       FORMAT ('   Molecule occupies special position')
      ELSE
         ISPECIAL = 0
      ENDIF
!
!----On unit # 22: line 1 = cell constants
!                  line 2 = space group, Z, special position code,
!                           refcod, # chars in refcod and center
!                           of symmetry code (1 if center at origin)
      IF (I_WHAT .NE. 5)                                        &
     &  WRITE (22,805) CELL, SPGRP, IZ, ISPECIAL, REFCOD, NC,   &
     &                 I_CENTER
805        FORMAT (6F10.4/A8,2X,2I5,2X,A8,2I5)
!
!----Read R-factor, REMARK, DISORD and ERROR information starting
!     with line # 3
910   NC = NVALUE(5) + NVALUE(6) + NVALUE(7) + NVALUE(8)
      IF (NC .EQ. 0) GO TO 100                 ! if NC = 0, no such info
      READ (10,40,END=1000) (LLINE(I),I=1,NC)
!----Report R-factor information
      IEND = NVALUE(5)
      IF (IEND .EQ. 0) GO TO 71
      PRINT 70, (LLINE(I),I=1,IEND)
      WRITE (23,70) (LLINE(I),I=1,IEND)
70    FORMAT ('   R-factor information...',80A1)
!----Look for REMARK information
71    IF (NVALUE(6) .EQ. 0)  GO TO 80
         ISTART = IEND+ 1
         IEND = NVALUE(5) + NVALUE(6)
         PRINT 75, (LLINE(I),I=ISTART,IEND)
         WRITE (23,75) (LLINE(I),I=ISTART,IEND)
75       FORMAT ('   Remark comments...'/(5X,60A1))
!----Look for DISORD information
80    IF (NVALUE(7) .EQ. 0) GO TO 90
        ISTART = IEND + 1
        IEND = IEND + NVALUE(7)
        PRINT 85, (LLINE(I),I=ISTART,IEND)
        WRITE (23,85) (LLINE(I),I=ISTART,IEND)
85      FORMAT ('   Disorder comments...'/(5X,60A1))
!----Look for ERROR information
90    IF (NVALUE(8) .EQ. 0) GO TO 100
        ISTART = IEND + 1
        PRINT 86, (LLINE(I),I=ISTART,NC)
        WRITE (23,86) (LLINE(I),I=ISTART,NC)
86      FORMAT ('   Error comments...'/(5X,60A1))
!----Read symmetry operations...are there any?  NVALUE(9)
!     [sym op count] should not be zero.
100   IF (NVALUE(9) .EQ. 0) THEN
         PRINT 101
101      FORMAT (' *** No sym op data on CSD FDAT file...STOP')
         STOP
      ELSE
         NC = NVALUE(9)*15
         READ (10,401,END=1000) (LLINE(I),I=1,NC)
401      FORMAT (75A1,5X)
      ENDIF
!----Process sym op information
      IEND = 0
110   ISTART =  1
      IEND = NVALUE(9)
         READ (LLINE500,112) (((ISYM(I,J,K),I=1,3), TRANS(J,K),J=1,3),   &
     &                         K=ISTART,IEND)
112      FORMAT (150(3I1,F2.0))
!----Read and skip atomic radii information...are there any?
111   IF (NVALUE(10) .EQ. 0) GO TO 120
        NC = NVALUE(10)*5
        READ (10,40,END=1000) (LLINE(I),I=1,NC)
!----Is there atom information?
120   NATOMS = NVALUE(11)+ NVALUE(12)   ! atoms + symmetry related atoms
      IF (NATOMS .EQ. 0) GO TO 990
!----Which atom format?
      IF (NVALUE(17) .EQ. 0) GO TO 980
      IF (NVALUE(17) .EQ. 1) THEN
        READ (10,125) (TLABEL(I), (IXYZ(J,I),J=1,3), I=1,NATOMS)
125     FORMAT (4(A5,3I5))
      ELSE
        READ (10,128) (TLABEL(I), (IXYZ(J,I),J=1,3), I=1,NATOMS)
128     FORMAT (A5,3I7,1X,A5,3I7,1X,A5,3I7)
      ENDIF
!----Construct a new TLABEL like C_10 ... assume that the
!     element symbol has only 1 character or that it's BR
      DO 600 I=1,NATOMS
        TEMP = TLABEL(I)
          IF (TLABEL(I)(1:2) .EQ. 'BR') THEN
             TEMP(3:3) = '_'
             I601 = 4
          ELSE
             TEMP(2:2) = '_'
             I601 = 3
          ENDIF
        DO 601 L=I601,5
           LM = L - 1
           TEMP(L:L) = TLABEL(I)(LM:LM)
601     CONTINUE
        TLABEL(I) = TEMP
600   CONTINUE
!
!----Convert to fractional coordinates and simultaneously put H and O atoms
!     at the bottom of the list.  If atom_order = 1, then do not change atom
!     sequence.                                            3/6/07
1045  NA = 0
      N_OX = 0
      N_HY = 0
      print 121
121   format ('Should atom order be altered, Y/N/blank = Y: ',$)   ! 3/6/07
      read (*,'(a1)') char_atom_order                              ! 3/6/07
         if (char_atom_order .eq. 'Y' .or. char_atom_order .eq. 'y' .or. &
             char_atom_order .eq. ' ') go to 2222
             do i=1,natoms
                NA = NA + 1
                LABEL(NA) = TLABEL(I) ! LABEL and XYZ now contain label &
                DO J=1,3              ! frac coord info with H's & O's at bottom
                   XYZ(J,NA) = FLOAT(IXYZ(J,I))/DEN(NVALUE(17))
                enddo
             enddo
          go to 2226                                               ! 3/6/07
2222  DO 131 IPASS=1,3
      DO 132 I=1,NATOMS
         GO TO (134, 136, 135), IPASS
134      IF (TLABEL(I)(1:1) .EQ. 'C' .OR. TLABEL(I)(1:1) .EQ. 'F'  &
     &    .OR. TLABEL(I)(1:2) .EQ. 'BR'                            &
     &    .OR. TLABEL(I)(1:1) .EQ. 'S'                             &  ! 11/6/03
     &    .OR. TLABEL(I)(1:1) .EQ. 'I') GO TO 133                     ! 11/6/03
         IF (TLABEL(I)(1:1) .NE. 'N') GO TO 132
133          NA = NA + 1
             LABEL(NA) = TLABEL(I)     ! LABEL and XYZ now contain label &
             DO 130 J=1,3              ! frac coord info with H's & O's at bottom
130             XYZ(J,NA) = FLOAT(IXYZ(J,I))/DEN(NVALUE(17))
         GO TO 132
135      IF (TLABEL(I)(1:1) .EQ. 'H') THEN  ! H's must go at bottom
             N_HY = N_HY + 1
             GO TO 133
         ENDIF
         GO TO 132
136      IF (TLABEL(I)(1:1) .EQ. 'O') THEN  ! O's go after C & N and
             N_OX = N_OX + 1                ! before H
             GO TO 133
         ENDIF
132   CONTINUE
131   CONTINUE
!
!----Transform coordinates?
2226     if (ltmatrix) then      ! apply xyz transformation matrix
            do i=1,na
                do k=1,3
                    new_xyz(k) = tt(k)
                    do j=1,3
                       new_xyz(k) = new_xyz(k) + tr(k,j)*xyz(k,i)
                    enddo
                 enddo
                 forall (k=1:3)
                    xyz(k,i) = new_xyz(k)
                 end forall
             enddo
          endif
      N_OTHER = NATOMS - N_HY - N_OX   ! N_OTHER = # C, N, S, F's
!
!----Convert fractional to orthogonal angstrom coordinates.
!----A, B, C = new orthoginal axes; a, b, c = original crystal axes...
!     A parallel to a
!     B in the a,b plane perpendicular to a
!     C perpendicular to A and B
      CALL MATRXT (CELL, R2)
!----Now multiply by the 3x3 R3 matrix so that...
!     A is parallel to a
!     B is in the a,c plane perpendicular to a
!     C is perpendicular to A and B
      DO 1201 I=1,3
         DO 1201 J=1,3
            R(I,J) = 0.0
            DO 1201 K=1,3
               R(I,J) = R(I,J) + R3(I,K)*R2(K,J)
1201  CONTINUE
!
139   IF (I_WHAT .EQ. 5) GO TO 143
!
!----Begin to write a temporary file on unit # 13
      REWIND 13                        ! cartesian coords...no cell required
!
!----Convert from fractional (XYZ) to cartesian (OR_XYZ) coords
143   DO 148 M=1,NATOMS
         DO 152 L=1,3
           OR_XYZ(L,M) = 0.0
           DO 152 J=1,3
              OR_XYZ(L,M) = OR_XYZ(L,M) + R(L,J)*XYZ(J,M)
152     CONTINUE
148   CONTINUE
!
!----If I_WHAT = 12, make CHEMX file and quit
160   IF (I_WHAT .EQ. 12) THEN
        OPEN (UNIT = 29, FILE = OUTFIL, STATUS = 'UNKNOWN')
        WRITE (29,2911) CELL                                  ! 6/1/03
2911    FORMAT (38X,3F8.3/21X,3F8.3)                          ! 6/1/03
        CALL ADJUST
        WRITE (29,2908) NATOMS                                ! 6/1/03
2908    FORMAT (I4,'   O'/)                                   ! 6/1/03
        DO J=1,NATOMS
           K = INDEX(LABEL(J)(1:5),'_')
           IF (K .NE. 0) LABEL(J)(K:K) = ' '
           WRITE (29,2912) J, LABEL(J), (XYZ(I,J),I=1,3)
2912       FORMAT (I4,1X,A5,3F10.5)
        ENDDO
        STOP
      ENDIF
!
!----If I_WHAT = 5, make CHEM3D file and quit
      IF (I_WHAT .EQ. 5) THEN
         CALL MAKE_CHEM3D
         PRINT 163, CHEM3D_NAME
163        FORMAT (' CHEM3D file created is ',A19)
         GO TO 178
      ENDIF
!
!----If I_WHAT = 10, adjust H and other appropriate atom positions
      IF (I_WHAT .NE. 10) GO TO 1095
         CALL ADJUST
         CALL OPEN_CHEM3D
         CALL MAKE_CHEM3D
         PRINT 163, CHEM3D_NAME
         OPEN (UNIT=19, STATUS='UNKNOWN')
         IF (I_CH_TYPE .EQ. 2) OPEN (UNIT=13, STATUS='UNKNOWN')
!----Write adjusted fractional coordinates on file # 19 & orthog coords
!     on # 13 for MNDO & # 11 for G98
         DO J=1,NATOMS
            WRITE (19,1096) LABEL(J)(1:4), (XYZ(I,J),I=1,3)
1096        FORMAT (A4,2X,3F10.5)
            IF (I_CH_TYPE .EQ. 2) WRITE (13,2099) J, LABEL(J),   &
     &                            (OR_XYZ(I,J),I=1,3)   ! MNDO on 13
2099        FORMAT (I4,4X,A5,3F10.5)
            IF (I_CH_TYPE .EQ. 3) WRITE (11,3011) LABEL(J)(1:4),   &
     &                            (OR_XYZ(I,J),I=1,3)   ! G98 on 11
3011        FORMAT (3X,A4,1X,3F12.6)
         ENDDO
         IF (I_CH_TYPE .EQ. 3) THEN
            WRITE (11,*)
            CLOSE (UNIT=11)
         ENDIF
         CLOSE (UNIT=19)
      IF (I_CH_TYPE .EQ. 2) CALL CONVERT2MOPAC
         GO TO 760      ! go to write sym op info on unit # 22
!
!----If I_WHAT = 4 or 8, adjust H and nitro O positions
!
1095  IF (I_WHAT .EQ. 4 .OR. I_WHAT .EQ. 8) CALL ADJUST
      CALL MAKE_CHEM3D                                       ! CHEM3D
      print 8888, CHEM3D_NAME                                ! file
8888  format (' Finished adjusting H and/or O positions...',   & ! and
     &        'CHEM3D file name is ',A19)                    ! name
      IF (I_WHAT .EQ. 8) GO TO 760
!
!----Write info for molecular fit input file on unit # 12
         WRITE (12,146) MF_TITLE, N_OTHER, NATOMS, NATOMS    ! top 3 lines...
146      FORMAT ('/mm/predictions/molecular-fit.out << END'/  &
     &           '    1'/2X,A72/                        &    ! mol fit line # 1
     &           3I5,'    0    1    5    0    1')            ! mol fit line # 2
         WRITE (12,141) CELL                                 ! line # 3 = cell line
141      FORMAT (6F10.3)
!
!----Write no. of sym ops and no. of atoms on unit # 22
760   MULT = 1
!----If ISPECIAL = 0 and I_CENTER = 1, make center-related sym ops
      IF (ISPECIAL .EQ. 0 .AND. I_CENTER .EQ. 1) MULT = 2  ! double sym ops
      WRITE (22,'(2I5)') (MULT*N_OPS), NATOMS
      IF (I_WHAT .EQ. 4 .OR. I_WHAT .EQ.10) GO TO 750
!
!----Prepare PCKME ...csf file
      PRINT 810
810   FORMAT (' Title for PCKME csf file: ',$)  ! 68 characters
      READ (5, 815) MTITLE                   ! max
815   FORMAT (A68)
      WRITE (16,820) MTITLE                  ! csf card # 1...
820   FORMAT (A68,'  -1')                       ! 1st title line
      WRITE (16,821) REFCOD, SPGRP, IZ          ! card # 2...
821   FORMAT (1X,A8,2X,A8,'  Z =',I2,T72,'0')   ! 2nd title line
      WRITE (16,822) (MULT*N_OPS), NATOMS                    ! csf
822   FORMAT (2I5,'    1    0    1    3   20    0  0.2000')  ! card # 2
      WRITE (16,824)                            ! csf cards # 3-5
824   FORMAT ('   1.   1.   1.   1.   1.   1.   1.   1.   1.   1.'/  &
     &        '   1.   1.   1.   1.   1.   1.   1.   1.   1.   1.'/  &
     &        '     1.0     0.0     0.0     0.0     1.0     0.0',    &
     &        '     0.0     0.0     1.0'/                            &
     &        '       0.0       0.0       0.0       0.0       0.0',  &
     &        '       0.0')
!----LATSYM only correct for trigonal - tetragonal
      IF (LATTICE(K_TYPE) .LE. 3) THEN
         LATSYM = 1      ! triclinic. monoclinic, orthorhombic
      ELSE
         LATSYM = 2      ! tetragonal, hexagonal
      ENDIF
      WRITE (16,826) CELL, LATSYM               ! csf card # 6...
826   FORMAT (6F10.3,I5)                        ! lattice constants
!
!----Subtract 1 from rotation part and divide translation by 12
750     DO 113 K=1,N_OPS
            DO 114 J=1,3
               DO 115 I=1,3
                  ISYM(I,J,K) = ISYM(I,J,K) - 1
115            CONTINUE
               TRANS(J,K) = TRANS(J,K)/12.0
114         CONTINUE
!----Write symop info on unit # 22
            L = -1
            DO 116 MMULT=1,MULT
               L = -L
               IF (I_WHAT .EQ. 4 .OR. I_WHAT .EQ. 10) THEN
                  WRITE (22,1200) (TRANS(J,K),   &
     &                            (L*(ISYM(I,J,K)),I=1,3),J=1,3)
1200              FORMAT (3(F15.2,3I3))
               ELSE
                  WRITE (16,1210) (TRANS(J,K),J=1,3),       &     ! csf cards # 7
     &                    ((FLOAT(L*(ISYM(I,J,K))),I=1,3),J=1,3) ! symmetry
1210              FORMAT (3F10.4,9F4.0)
               ENDIF
116         CONTINUE
113      CONTINUE
!----Write 'END' on unit # 22
      IF (.NOT. (I_WHAT .EQ. 4 .OR. I_WHAT .EQ. 10)) GO TO 793
         WRITE (22,173)
         CLOSE (UNIT=22)
793   IF (I_WHAT .EQ. 10) GO TO 2000
      IF (I_WHAT .EQ. 4) GO TO 761
!
!----Finish PCKME ...csf file
839   WRITE (16,835)
835   FORMAT (' 7 7 7 7 7 7 7 7 7 7 7 7'/                            &
     &        ' 4 4 4 4 4 4 4 4 4 4 4 4'/                            &
     &        ' 0 0 0 0 0 0 0 0 0 0 0 0'/                            &
     &        '   0.00020   0.00050         1       -1.         0',  &
     &        '     0.000   1.0'//)
!
!----Finished creating PCKME ...csf file
!
!----Make PCKME ...mcf file on unit # 17
!
      IF (I_WHAT .NE. 9) GO TO 859
         DO I=1,NKA           ! WMIN input file; pick up atom
            READ (25,858) LABEL(I)(1:4), CHARGE(I)  ! atom charges
858         FORMAT (A4,5X,F9.4)
         ENDDO
859   WRITE (17,860) NATOMS       ! mcf line # 1
860   FORMAT (I4,'   1   2')      ! set up to use qudrupoles
      DO 890 I=1,NATOMS
         IF (I_WHAT .EQ. 8) THEN
            ATOMID = LABEL(I)(1:2)      ! input from CSD
         ELSE
            READ (25,3041) LINE54       ! input from WMIN
            ATOMID = LINE54(1:2)
         ENDIF
         DO 863 J=1,N_ATYPES            ! find atomic number
            IF (ATOMID .NE. A_TYPE(J)) GO TO 863
            NZ = AT_NUMBER(J)
            L = IPCKME_POT(J)
            GO TO 868
863      CONTINUE
         PRINT 866, ATOMID         ! if here, cannot find atom symbol match
866      FORMAT (/' Atom symbol ',A2,' is unknown...STOP')
         STOP
868      IF (I_WHAT .EQ. 8) THEN
            WRITE (17,870) NZ, (XYZ(L,I),L=1,3), L, LABEL(I)(1:4)   ! CSD input
870         FORMAT (I4,8X,3F16.9,'  1  0  0',I3,A4)
         ELSE
            READ (LINE54,872) A_NAME, W_XYZ     ! fractional coordinates from
872         FORMAT (A4,23X,3F9.5)               ! existing WMIN input file
            IF (LABEL(I)(1:4) .EQ. LINE54(1:4)) GO TO 871
               PRINT 879, LABEL(I)(1:4), LINE54(1:4)
879            FORMAT (' Atom label of ',A4,' from WMIN charge line'/   &
     &                 '  differs from label of ',A4,' found on',       &
     &                 '  coordinate line...STOP')
               STOP
871         WRITE (17,873) NZ, CHARGE(I), W_XYZ, L, LINE54(1:4)  ! put charges
873         FORMAT (I4,F7.4,1X,3F16.9,'  1  0  0',I3,A4)   ! from WMIN file into
         ENDIF                                             ! cols 5-11 of mcf line
890   CONTINUE
      STOP
!
!----Finished with both ...csf and ...mcf
!
!----Write files with fractional (molec fit) and cartesian (mopac) coords
761   DO 1160 NA=1,NATOMS
         WRITE (12,151) NA, LABEL(NA), (XYZ(J,NA),J=1,3)  ! fractional coords
         WRITE (13,151) NA, LABEL(NA), (OR_XYZ(J,NA),J=1,3)   ! cartesian coords
151      FORMAT (I4,1X,A5,3F10.5)
1160  CONTINUE
      CLOSE (UNIT=13)
!----Final line (for now) on molecular fit file
      WRITE (12,173)
173   FORMAT ('END')
      CLOSE (UNIT=12)
!
!----Convert file # 13 to mopac format
      CALL CONVERT2MOPAC
!
!----Finished
178   PRINT 170
170   FORMAT (' Transformation complete...happy trails')
      CLOSE (UNIT = 11)
      CLOSE (UNIT = 12)
      GO TO 77
980   PRINT 985
985   FORMAT (' No atom format information available...quit')
990   PRINT 995
995   FORMAT (' No atom information available...quit')
1000  STOP
!
!************************************************************
!
!----Decode MNDO/ESP file...pick up possible name from unit # 21
1010  REWIND 21
      READ (21,28)
      READ (21,28) OUTFIL       ! 2nd line is name of esp output file
!
      PRINT 1015, OUTFIL
1015  FORMAT (' Name of the MNDO/ESP input file to supply charges',  &
     &        ' and'/'    mopac-optimized coordinates: [',A19,']',$)
!----Input file
      READ 2, INFILE
      IF (INFILE(1:1) .NE. ' ') GO TO 1099     ! if blank, use default
         INFILE = OUTFIL                       ! name from unit # 21
1099  OPEN (UNIT=10,FILE=INFILE,STATUS='OLD')
!----Molecular fit output file...get the file name from unit # 21
      REWIND 21
      READ (21,'(A19)') OUTFIL
      PRINT 1017
1017  FORMAT (' The mopac-optimized coordinates will to added to',   &
     &            ' the end of the molecular fit input file')
      OPEN (UNIT=11,FILE=OUTFIL,STATUS='OLD',FORM='FORMATTED')
!
!----Search for 'POINT-CHG.' line...this points to the cartesian coordinates
1016  READ (10,41,END=1030) LINE80
      IF (LINE80(1:11) .NE. ' POINT-CHG.') GO TO 1016
!----Found it...skip the next 8 lines
      DO I=1,8
         READ (10,*,END=1030)
      ENDDO
!----9th line...write a cartesian cell line on # 11 for molecular fit
      REWIND 11
!----Search for the 'END' mark point
2011  READ (11,'(A3)',END=2022) END
         IF (END .NE. 'END') GO TO 2011
         BACKSPACE 11
         GO TO 2012
2022     PRINT 2013
2013     FORMAT (/' ***Cannot locate END point on molecular fit file')
         STOP
2012  WRITE (11,2039)
2039     FORMAT (3(7X,'1.0'),3(6x,'90.0'))
!----Read the cartesian coordinates
      NA = 0
2040  NA = NA + 1
         READ (10,1025,END=1030) LABEL(NA), (XYZ(J,NA),J=1,3)
1025     FORMAT (14X,A5,11X,3F10.4)
         IF (LABEL(NA) .EQ. '     ') GO TO 1050 ! blank line signals end of atoms
         IF (LABEL(NA)(1:1) .EQ. ' ') THEN
            DO L=1,4                              ! atom symbol has 1 character,
               LABEL(NA)(L:L) = LABEL(NA)(L+1:L+1)    ! shifted left 1 character
            ENDDO
            LABEL(NA)(5:5) = ' '
         ELSE
            IF (LABEL(NA)(1:2) .EQ. 'Br') LABEL(NA)(1:2) = 'BR'   ! make Br upper case
         ENDIF
         WRITE (11,151) NA, LABEL(NA), (XYZ(J,NA),J=1,3) ! for molecular fit
      GO TO 2040
1030  PRINT 1031
1031     FORMAT (' Unexpected end-of-file')
         STOP
1050  NA = NA - 1
!----Replace the 'END' at the bottom of the mol fit input file
      WRITE (11,173)
      WRITE (11,174)
174   FORMAT ('mv fort.16 molecular-fit.output-file')
!
!----Now get the scaled ESP charges
      CLOSE (UNIT=11)
      OPEN (UNIT=13,FILE='esp.charges',STATUS='UNKNOWN',   &
     & FORM='FORMATTED')
!----Look for the string 'SCALED CHARGE'
1060  READ (10,41,END=1030) LINE80
      IF (LINE80(37:49) .NE. 'SCALED CHARGE') GO TO 1060
      DO 1070 I=1,NA
         READ (10,1080) ANAME2, SCALED_Q
1080     FORMAT (20X,A2,16X,F7.4)
           IF (ANAME2 .EQ. 'Br') THEN
              ANAME2 = 'BR'
           ELSE
              ANAME2(1:1) = ANAME2(2:2)
              ANAME2(2:2) = ' '
           ENDIF
         WRITE (13,1090) ANAME2, SCALED_Q
1090     FORMAT (1X,A2,7X,F10.4)
1070  CONTINUE
!
      STOP
!
!********** create complete WMIN input file **************
!
!----Open various files...
2000  REWIND 21
      READ (21,*)
      READ (21,*)
      READ (21,'(A19)') OUTFIL  ! 3rd line is XXXXXX-wmin.input
      IF (I_WHAT .EQ. 10) OPEN (UNIT=19, STATUS='OLD')
!----Open unit # 14...for the complete WMIN input file
         OPEN (UNIT=14, FILE=OUTFIL, STATUS='UNKNOWN',   &
     &         FORM='FORMATTED')
      IF (I_WHAT .EQ. 10) GO TO 2010       ! I_WHAT = 10, skip file
!----Open unit # 15...with ESP charges     ! # 15 & # 10
         OPEN (UNIT=15,FILE='esp.charges',STATUS='OLD',   &
     &         FORM='FORMATTED')
      I4 = 1
!----Open unit # 10...molecular fit output file
         OPEN (UNIT=10,FILE='molecular-fit.output-file',   &
     &         STATUS='OLD',FORM='FORMATTED')
!----Open unit # 22...scratch file
         OPEN (UNIT=22,STATUS='OLD')
!----Read the molecular fit output file to extract the mopac
!     coordinates that have been fit to the CSD unit cell fractional
4050  READ (10,41,END=2090) LINE80
      IF (LINE80(64:79) .NE. 'IN CELL OF FIRST') GO TO 4050
2010  NA = 0
      WT_MOL = 0                 ! molecular weight
      DO 2059 I=1,3
2059     SUMXYZ(I) = 0.0
2060  NA = NA + 1
      IF (I_WHAT .EQ. 10) THEN
         READ (19,1096,END=2090) LABEL(NA)(1:4), (XYZ(J,NA),J=1,3)
      ELSE
         READ (10,3040,END=2090) LABEL(NA)(1:4), (XYZ(J,NA),J=1,3)
3040     FORMAT (11X,A4,42X,3F10.6)
      ENDIF
         ATOMID = LABEL(NA)(1:2)
            IF (ATOMID(2:2) .EQ. '_') ATOMID(2:2) = ' '
!----Get the atom type for WMIN
      DO 2045 K=1,N_ATYPES
         IF (ATOMID .NE. A_TYPE(K)) GO TO 2045
         WT_MOL = WT_MOL + W_TYPE(K)
         IJK_WMIN(NA) = K
         GO TO 2046
2045  CONTINUE
!----Sum the xyz's for a molecular center
2046  DO 2055 I=1,3
         SUMXYZ(I) = SUMXYZ(I) + XYZ(I,NA)
2055  CONTINUE
      GO TO 2060
!----Finished...calc xyz averages and write XTRA line
2090  DO 2070 I=1,3
         XYZ(I,NA) = SUMXYZ(I)/(NA - 1)   ! there are NA -1 atoms
2070  CONTINUE                            ! NAth line of the XTRA
      LABEL(NA)(1:4) = 'XTRA'             ! the molecular center
      IJK_WMIN(NA) = 0                    ! identifies XTRA
!
!----Ask for title line for WMIN file
         PRINT 915
915      FORMAT (' Title for WMIN file: ',$)
         READ 916, TITLE
916      FORMAT (A60)
!----Create lines 1 & 2 for WMIN
      N_WMIN = 2
      WRITE (14,905)
905   FORMAT ('        0      100       55       16',   &
     &        '        1        0        1        1')
      WRITE (14,906)
906   FORMAT ('        1      700     5500       30',   &
     &        '       68       20       40       20')
!----Get cell, space group and Z from scratch file # 22
       READ (22,'(6F10.3)',END=9999) (CELL(I),I=1,6)
       READ (22,'(A8,2X,2I5,15X,I5)',END=9999) SPGRP, IZ,   &
     &            ISPECIAL, I_CENTER
!----Create 3rd WMIN line
      WRITE (14,218) TITLE, SPGRP
218   FORMAT (A60,2X,A8)
      N_WMIN = 3
!----Create 4th WMIN line
      WRITE (14,718)
718   FORMAT ('  3 99  0  0  1  0  0  5  0   0.0001',   &
     &        ' 0.000001      0.1     0.01     10.0')
      N_WMIN = 4
!----Read number of symops  and atoms from scratch file # 22
      READ (22,'(2I5)',END=9999) N_OPS, NATOMS
!----Pick up the space group array number
      K_TYPE = 0
      DO 1058 I=1,N_SPGRPS
         IF (SPGRP .NE. SG_TYPES(I)) GO TO 1058
         K_TYPE = I
         GO TO 300
1058  CONTINUE
         PRINT 3060
3060     FORMAT (' Cannot identify space group')
         STOP
!----Use K_TYPE to set the rep lattice directions
300   GO TO (301, 302, 303, 304, 305, 306, 307, 308), LATTICE(K_TYPE)
301      I1 = 1                ! triclinic
         I2 = 1
         I3 = 0
         GO TO 1059
302      I1 = 1                ! monoclinic
         I2 = 0
         I3 = 0
         GO TO 1059
303      I1 = 0                ! orthorhombic
         I2 = 0
         I3 = 0
         GO TO 1059
304      I1 = 0                ! tetragonal, 4/mmm Laue symm
         I2 = 0
         I3 = 1                ! a = b
         GO TO 1059
305      I1 = 1                ! tetragonal, 4/m Laue symmetry
         I2 = 0
         I3 = 1                ! a = b
         go to 1059
306      I1 = 1                ! rhombohedral
         I2 = 0
         I3 = 3                ! a = b = c, alpha = beta = gamma
         go to 1059
307      i1 = 0                ! cubic
         i2 = 0
         i3 = 2                ! a = b = c
         go to 1059
308      i1 = 1                ! trigonal/hexagonal
         i2 = 1
         i3 = 1                ! a = b
!----Create 5th WMIN line
1059  WRITE (14,919) ISB(K_TYPE), TITLE
919   FORMAT (I3,A60)
      N_WMIN = 5
!----Create 6th and 7th WMIN lines
      IF (I_WHAT .NE. 10) GO TO 2037
      IF (I_CH_TYPE .LE. 1) THEN
          I4 = 0                 ! no coulombic term
      ELSE
          I4 = 1                 ! coulombic terms
      ENDIF
2037  WRITE (14,1038) (NATOMS+1), NATOMS, N_OPS, I1, I2, I4, I3
1038  FORMAT (3I3,'  2',2I3,'  0  1    -2',I3,   &
     &        '  1  2              1  0',I3/     &
     &        '     0.25  2  2')
      N_WMIN = 7
!----Create line # 8
      WRITE (14,'(F9.2,I3)') WT_MOL, IZ
      N_WMIN = 8
!----Create lines # 9-12
      WRITE (14,4100) (CELL(I),I=1,3), (COSD(CELL(I)),I=4,6)
4100  FORMAT ('      0.4      0.5'/'      5.5      6.0'/   &
     &        3F9.3,3F9.5/6('    0.001'))
      N_WMIN = 12
!----Write sym ops to WMIN input file
      DO I=1,N_OPS
         READ (22,34,END=9999) SYM_OP
         N_WMIN = N_WMIN + 1
         WRITE (14,34) SYM_OP
      ENDDO
      IF (I_WHAT .EQ. 10) GO TO 4600
!----Make the atom parameter lines with coefficients and esp
!     charges from # 15.  Pick up other atom info from the
!     LABEL, IJK and XYZ arrays.
      DO 4300 I=1,NATOMS
          READ (15,'(1X,A1,8X,F10.4)',END=4310) A_ID(I), CHARGE(I)
          IF (LABEL(I)(1:1) .EQ. A_ID(I)) THEN
             K = IJK_WMIN(I)
             WRITE (14,4350) LABEL(I)(1:4), CHARGE(I),       &     ! 5/3/03
     &                  (ABC_COEFS(M,J),M=1,3),  W_TYPE(K)        ! 5/3/03
4350         FORMAT (A4,5X,5F9.4)
             N_WMIN = N_WMIN + 1
          ELSE
             PRINT 4351
4351         FORMAT (' Atom label mismatch between esp.charges and',  &
     &               ' mol fit output files')
             STOP
          ENDIF
4300  CONTINUE
      GO TO 4700
!----I_WHAT = 10; if I_CH_TYPE = 0 or 1, then WMIN file will be complete
!     otherwise "_CHARGE__" will be placed in the atomic charge location
4600  DO I=1,NATOMS                                               ! 5/3/03
          K = IJK_WMIN(I)
          J = NPOT(I)      ! determine atom potential number
          IF (J .EQ. 999) THEN    ! J = 999
             AC = 0.0             !  indicates that
             BC = 0.0             !  atom type is
             CC = 0.0             !  unknown
          ELSE
             AC = ABC_COEFS(1,J)                                  ! 5/3/03
             BC = ABC_COEFS(2,J)                                  ! 5/3/03
             CC = ABC_COEFS(3,J)                                  ! 5/3/03
          ENDIF
          IF (I_CH_TYPE .GE. 2) THEN                              ! 5/3/03
             C_WORD = "_CHARGE__"
          ELSE
             C_WORD = "         "
          ENDIF
          WRITE (14,4610) LABEL(I)(1:4), J, C_WORD, AC, BC, CC,    &
     &                    W_TYPE(K)
4610      FORMAT (A4,I5,A9,3F9.3,F9.4)                           ! 3/27/03
      ENDDO
!----Add the NA atom lines (+ XTRA)
4700  DO 4370 I=1,NA
          N_WMIN = N_WMIN + 1
          K = I
          IF (I .EQ. NA) K = 0
          WRITE (14,4371) LABEL(I)(1:4), K, (XYZ(J,I),J=1,3)
4371      FORMAT (A4,I5,18X,3F9.5)
4370  CONTINUE
!----Create rigid body identifier line
      I_ONE = 1
      WRITE (14,4380) (I_ONE,I=1,NA)
4380  FORMAT (24I3)
      N_WMIN = N_WMIN + 1
!----Create the last line
      WRITE (14,4385) NA, WMIN_REF_CODES(ISPECIAL+1,K_TYPE)
4385  FORMAT (I3,15X,'     0.05     0.01'/   &
     &        ' '/A12/'000000000000'/)
      N_WMIN = N_WMIN + 1
!----Finished!!!!
      PRINT 4389, N_WMIN
4389  FORMAT (' Preparation of WMIN input file complete...',I4,   &
     &        ' lines placed into input file')
      STOP             ! whoopie
4310  PRINT 4311
4311  FORMAT (' Unexpected end encountered on esp.charges file')
      STOP
9999  PRINT 9998
9998  FORMAT (' Unexpected end encountered on scratch file # 22')
!
!**************# 6...CHEM3D input; G98 output********************
!
!----Read CHEM3D file and make input file for G98 charge calcn
1005  READ (20,1139,ERR=1135,END=1135) CHEM3D_NAME
      GO TO 1137
1135  PRINT 1138
1138  FORMAT (' CHEM3D name file not found; give name of CHEM3D ',   &
     &        'coordinate file: ',$)
      READ 1139, CHEM3D_NAME
1139  FORMAT (A19)
1137  PRINT 1171
1171  FORMAT ('   Does the CHEM3D file have atom serial',   &
     &        ' numbers? [N]: ',$)
      READ '(A1)', SERIAL_NUM
        IF (SERIAL_NUM .EQ. 'N' .OR. SERIAL_NUM .EQ. 'n'   &
     &      .OR. SERIAL_NUM .EQ. ' ') THEN
            SERIAL = .false.
        ELSE
            SERIAL = .true.
        ENDIF
      OPEN (UNIT=30, FILE=CHEM3D_NAME, STATUS='OLD')
      OPEN (UNIT=24, FILE='631gstar_chelpg.com', STATUS='UNKNOWN')
      PRINT 1151
1151  FORMAT (' Name of G98 input file is 631gstar_chelpg.com')
!----Extract compound name from the chem3d file name...look for '.'
      NC = 0
      DO 1152 I=1,8
         IF (CHEM3D_NAME(I:I) .EQ. '.') GO TO 1153
         NC = NC + 1
1152  CONTINUE
!----Write header info for G98 input file on unit # 24
1153  WRITE (24,1141) CHEM3D_NAME(1:NC)
1141  FORMAT ('$RunGauss'/                                               &
     &        '%chk=6-31Gstar_chelpg.chk'/                               &
     &        '#p 6-31G* scf=direct pop=CHELPG geom=coord'//             &
     &        A8,', 6-31G*, charges from ESP/CHELPG, cart coord input,', &
     &        ' no optim'//                                              &
     &        '   0   1')
      READ (30,*) N          ! number of atoms from CHEM3D file
      DO 1143 I=1,N
         IF (SERIAL) THEN
            READ (30,1145) ANAME2, LINE36       ! with serial numbers
1145        FORMAT (1X,A2,5X,A36)
	 ELSE
            READ (30,1146) ANAME2, LINE36       ! w/out serial numbers
1146        FORMAT (1X,A2,A36)
         ENDIF
         ANAME3 = '     '
         ANAME3(1:1) = ANAME2(1:1)
         J = 2
         IF (ANAME2(2:2) .EQ. ' ') GO TO 2138
            ANAME3(2:2) = ANAME2(2:2)
            J = 3
2138     ANAME3(J:J) = '_'
         J = J + 1
         IF (I .LE. 9) THEN
            WRITE (ANAME3(J:J+1),1147) I
1147        FORMAT (I1)
         ELSE
            WRITE (ANAME3(J:J+1),1148) I
1148        FORMAT (I2)
         ENDIF
         WRITE (24,1149) ANAME3, LINE36
1149     FORMAT (2X,A5,1X,A36)
1143  CONTINUE
      WRITE (24,*)
      CLOSE (UNIT=24)
      STOP
!
!**********# 7...G98 & wmin.input; updated wmin.input created***********
!
3007  OPEN (UNIT=25, FILE='wmin.input', STATUS='OLD')          ! old wmin input file
      OPEN (UNIT=26, FILE='wmin.input2', STATUS='UNKNOWN')     ! updated wmin file
      OPEN (UNIT=27, FILE='631gstar_chelpg.log', STATUS='OLD') ! G98/ESP charges
!----1st wmin line...# 11 from wmin documentation
      READ (25,3010) LINE72
3010  FORMAT (A72)
      WRITE (26,3010) LINE72
!----2nd wmin line...# 12; pick up number of atoms & sym ops
      READ (25,3010) LINE72
      READ (LINE72,3012) NA, NS
3012  FORMAT (I3,3X,I3)
      WRITE (26,3010) LINE72
!----3-8th & NS sym op wmin lines
      DO 3014 I=1,6+NS
         READ (25,3010) LINE72
         WRITE (26,3010) LINE72
3014  CONTINUE
!----Find charge area in G98 ...log file
3020  READ (27,3010) LINE72
      IF (LINE72(1:23) .NE. ' Charges from ESP fit,') GO TO 3020
      READ (27,3010) LINE72
      IF (LINE72(1:8) .NE. ' Charge=') GO TO 3020
      READ (27,*)                                  ! skip a blank line
!----Should find NA - 1 lines before a '-------' line
      NL = 0
3025  READ (27,3026) LINE36
3026  FORMAT (A36)
      IF (LINE36(2:12) .NE. '-----------') GO TO 3035
         IF (NL .EQ. (NA - 1)) GO TO 3042
         PRINT 3027
3027     FORMAT (' Problems reading G98 charge file...stop')
         STOP
3035  NL = NL + 1
      READ (LINE36,3039) ATOM_ID(NL), CHARGE(NL)
3039  FORMAT (5X,A1,F12.6)
      GO TO 3025
!----Have all G98 charges...read old wmin.input and replace charges;
!     show user what the charge differences are
3042  PRINT 3047
3047  FORMAT (' G98 (new) and old (replaced) charges...'/    &
     &        '  #  G98ID  q(G98)   oldID  q(old)    diffq')
      DO 3050 I=1,NA-1
         READ (25,3041) LINE54
3041     FORMAT (A54)
         IF (LINE54(1:1) .NE. ATOM_ID(I)) THEN
            PRINT 3043, LINE54(1:5), ATOM_ID(I)
3043        FORMAT (' Atom ID mismatch between wmin.input [',A5,   &
     &              '] and charge files [',A1,']...stop')
            STOP
         ENDIF
         READ (LINE54(11:18),'(F8.4)') OLDQ
         PRINT 3045, I, ATOM_ID(I), CHARGE(I), LINE54(1:5),   &
     &               OLDQ, (CHARGE(I) - OLDQ)
3045     FORMAT (I3,4X,A1,2X,F8.4,3X,A5,F8.4,F9.4)
         WRITE (LINE54(10:18),'(F9.4)') CHARGE(I)
         WRITE (26,3041) LINE54
3050  CONTINUE
!----Finish up
      NL = (NA + 23)/24
      DO 3069 I=1,NA+4+NL
         READ (25,3010) LINE72
         WRITE (26,3010) LINE72
3069  CONTINUE
      CLOSE (UNIT=25)
      CLOSE (UNIT=26)
      CLOSE (UNIT=27)
      STOP
!
!************** option # 11 - add charge data to a WMIN file **********
!
6030  OPEN (UNIT=21, STATUS='OLD') ! open names file
      READ (21,*)                  ! skip 1st line
      READ (21,28) OUTFIL          ! line # 2 = name of file with charges
      IF (OUTFIL .EQ. '631gstar_chelpg.log') THEN
         CH_TYPE = .true.          ! G98 631g* charges
      ELSE
         CH_TYPE = .false.         ! MNDO charges by default
      ENDIF
      OPEN (UNIT=27, FILE=OUTFIL, STATUS='OLD')
      READ (21,28) OUTFIL          ! line # 3 = name of WMIN input file
      OPEN (UNIT=25, FILE=OUTFIL, STATUS='OLD')
      OPEN (UNIT=26, FILE='wmin.input2', STATUS='UNKNOWN') ! new WMIN file
      OPEN (UNIT=28, FILE='run.wmin', STATUS='UNKNOWN')  ! for running wmin calcns
      WRITE (28,41) WMIN_COM
      CLOSE (UNIT=28)
      CALL SYSTEM ('chmod 755 run.wmin')
!----Find charge section in the MNDO or G98 file
      IF (CH_TYPE) THEN
6088     READ (27,'(A27)',END=6400) LINE27  ! G98 631g* file
         IF (LINE27 .NE. ' Fitting point charges to e') GO TO 6088   ! 12/3/03
         READ (27,*)      ! skip one line                              12/3/03
      ELSE
6090     READ (27,6100,END=6400) LINE31     ! MNDO file
6100     FORMAT (12X,A31)
         IF (LINE31 .NE. 'ELECTROSTATIC POTENTIAL CHARGES') GO TO 6090
      ENDIF
         READ (27,*)      ! skip next 2 lines
         READ (27,*)
!----Position the WMIN input file
6110  READ (25,41) LINE80
      IF (LINE80(10:18) .EQ. '_CHARGE__') GO TO 6116
         WRITE (26,41) LINE80    ! not a charge line...copy to output file
         GO TO 6110
!----Found the 1st charge line; now read both WMIN and charge files
6114  READ (25,41) LINE80
6116  IF (LINE80(10:18) .NE. '_CHARGE__') GO TO 6120
6115  READ (27,34) LINE72         ! a line from the charge file
      IF (CH_TYPE) THEN
         IF (.not.(LINE72(6:6) .EQ. LINE80(1:1) .OR.        &    ! 6/19/03
     &       LINE72(9:9) .EQ. LINE80(1:1))) GO TO 6300           ! 6/19/03
             LINE80(10:18) = LINE72(11:19)              ! 12/3/03
      ELSE
         IF (LINE72(22:22) .NE. LINE80(1:1)) GO TO 6300
            LINE80(10:18) = LINE72(37:45)
      ENDIF
      WRITE (26,41) LINE80    ! write the charge-corrected WMIN line
      GO TO 6114
!----No more charge lines; finish up the WMIN file
6118  READ (25,41,END=6320) LINE80
6120  WRITE (26,41) LINE80
      GO TO 6118
!----Finished
6320  PRINT 6321
6321  FORMAT ('Creation of charge-corrected WMIN file (wmin.input2) ',  &
     &        'complete'/' run.wmin can be used to execute WMIN')
      STOP
6300  PRINT 6301, LINE72, LINE80
6301  FORMAT (' Atom symbol mismatch...'/   &
     &        '  charge file line = ',A72/  &
     &        '  WMIN file line = ',A80)
      STOP
6400  PRINT 6401
6401  FORMAT (' End of charge data file reached without',   &
     &        ' finding the charge information')
      STOP
      END PROGRAM PREPARE_FILES
!
!--------------------------------------------------------------------
      FUNCTION SIND(X)
        REAL :: SIND, X
        SIND = SIN(X/57.29577951)
        RETURN
      END
      FUNCTION COSD(X)
        REAL :: COSD, X
        COSD = COS(X/57.29577951)
        RETURN
      END
!
!--------------------------------------------------------------------------
!----CONVERT2MOPAC
!
      SUBROUTINE CONVERT2MOPAC
!
     USE mod_preppot, only : MTITLE, MOPAC_TYPE, LABEL
!!I_CH_TYPE type is not in use
      IMPLICIT NONE
!
!      CHARACTER (LEN=1)  :: TEMP(6)
!TEMP(6)*1 and TEMP2*6 is removed and replaced with TEMP*6
!      CHARACTER (LEN=5)  :: NAME(100) ! changed to LABEL
      CHARACTER (LEN=6)  :: TEMP
!DOT_DAT*4, TNAME*1, TRIAL_FILE(32)*1 , OUTFILE2(32)*1, INFILE2(32)*1
!FORMAT1(46)*1, FORMAT2*46, FORMAT3*2 removed, not in use
      INTEGER :: I, IK(100), K, L, NA
!LIST(3) removed, not in use
      REAL    :: XYZ(3,100)
!A(6), CELL(6), TXYZ(3), T(3,3), XYZ_CART(3,100) removed, not in use
!
!      EQUIVALENCE (TEMP, TEMP2)
!---Set up temporary coordinate file
      REWIND 13
!----Mopac output file on unit # 11
!
!----Cartesian coords...no cell params needed
      NA=0
777   NA=NA+1
      READ (13,30,END=999) IK(NA), LABEL(NA), (XYZ(I,NA),I=1,3)
30    FORMAT (I4,4X,A5,3F10.5)
      GO TO 777
  999 NA=NA-1
      REWIND 11
7777  GO TO (771, 772), MOPAC_TYPE
771     WRITE (11,10)
10        FORMAT(' XYZ NOINTER MMOK AM1 PRECISE T=240H')
        GO TO 775
772     WRITE (11,773)
773       FORMAT(' MNDO ESP 1SCF GNORM=0.1 T=240H SCFCRT=1.0D-4')
775   WRITE (11,15) MTITLE
   15 FORMAT (2X,A72)
      WRITE (11,20)
   20 FORMAT (1X,'--------------------------------------------------')
      DO 50 K=1,NA
         TEMP = '      '
            IF (LABEL(K)(1:2) .EQ. 'BR') THEN
               TEMP(1:2) = LABEL(K)(1:2)
               L = 3
            ELSE
               TEMP(1:1) = LABEL(K)(1:1)
               L = 2
            ENDIF
            TEMP(L:L) = '('
         DO 51 I=L+1,5
            IF (LABEL(K)(I:I) .EQ. ' ') GO TO 52
            TEMP(I:I) = LABEL(K)(I:I)
51       CONTINUE
         IF (I .EQ. 5) I = 6
52       TEMP(I:I) = ')'
!      WRITE (11,25)(LABEL(K),(XYZ(I,K),I=1,3),K=1,NA)
!  25 FORMAT (1X,A5,F10.5,2X,'1',F10.5,2X,'1',F10.5,2X,'1')
      WRITE (11,25) TEMP, (XYZ(I,K),I=1,3)
   25 FORMAT (1X,6A1,F10.5,2X,'1',F10.5,2X,'1',F10.5,2X,'1')
50    CONTINUE
      WRITE (11,35)
   35 FORMAT (1X,'0'//)
      CLOSE (UNIT=11)
      RETURN
      END SUBROUTINE CONVERT2MOPAC
!------------------------------------------------------------------
!     subroutine mtime(t1,t2)
!
!  subroutine to return cpu and wall clock time in seconds
!  current revision date: 3-jun-1991 by mha
!
!  -----------------------------------------------------------------
!  variables in call list
!
!    t1:     on return:  contains cpu clock time in seconds
!    t2:     on return:  contains wall clock time in seconds
!
!  -----------------------------------------------------------------
!       implicit double precision(a-h,o-z)
!     double precision cpu,sys,tio
!     call timing(cpu,sys,tio)
!     t1=cpu
!     t2=sys+tio+cpu
!     return
!     end
!....should return time, date and machine type
!     subroutine timdat(tim,dat,mach)
!     character *(*) tim,dat,mach
!     character*9 a(7)
!     integer uname
!     call timec(tim)
!     call datec(dat)
!     l=len(tim)
!     if(l.gt.8) tim(9:l)=' '
!     l=len(dat)
!     if(l.gt.9) dat(10:l)=' '
!     i=uname(a)
!     mach=' '
!     mach(1:5)=a(1)(1:5)
!     mach(7:9)=a(3)(1:3)
!     mach(11:18)=a(5)(1:8)
!     return
!     end
!
!------------------------------------------------------------
!
!     Determine atom types for potential parameter selection for
!      inclusion in a WMIN input file
!
      FUNCTION NPOT (NUM_ATOM)
!
      USE mod_preppot, ONLY : ICON, NCON, CONNECT_FLAG
      USE mod_preppot, ONLY : LABEL
!
      IMPLICIT NONE
      INTEGER :: NPOT, NUM_ATOM
!
      CHARACTER (LEN=1) :: ELEMENT                    ! combo
      CHARACTER (LEN=2) :: AT2
      INTEGER :: I, I1, I2, IT, IT1, IT2, ITT, AMINE, no1              ! 7/20/06
      INTEGER :: IQ_ALCOHOL, IQ_ALKYNE, IQ_CUBANE, IQ_HN_IMIDE         ! 12/31/03
      INTEGER :: IQ_SULFONIMINE                                        ! 12/31/03
      INTEGER :: IQ_OC_IMIDE, IQ_FUROXAN_O1, IQ_NITRO_OX, IQ_OXYGEN
      INTEGER :: IQ_FUROXAN_O2, IQ_AZIDE_N, IQ_N_AZAPENTALENE
      INTEGER :: IQ_FUROXAN_N2, IQ_NITRATE_N, IQ_FUROXAN_N3
      integer :: iq_123triazole_2, iq_123triazole_3                    ! 12/30/05
      integer :: iq_1234tetrazole_n1, iq_1234tetrazole_n2              ! 12/31/05
      integer :: iq_1234tetrazole_n3, iq_1234tetrazole_n4              ! 12/31/05
      integer :: iq_pyrrole_nh, iq_thiazole_n, iq_thiazole_s           !  1/31/06
      integer :: iq_aminoimine_2, iq_aminoimine_3                       ! 1/23/06
      integer :: iq_pyridinium_n                                       ! 7/20/06
      integer :: iq_ester_c, nso1, nso2, iq_diazo_n                     ! 2/10/06
      INTEGER :: N, NB, NC, NN, NO, NSO, N_OXY, N_OXY_1, N_CARBS, N_CARBS_3
      INTEGER :: N_NIT, N_HYD, N_F, N_S_4, N_ADJ, NC_ADJ, NUM_O(4)
!
!----Does the atom connectivty need to be established...once is enough
      IF (CONNECT_FLAG) GO TO 10
          CALL CONNECT2
          CONNECT_FLAG = .TRUE.
!
!----Determine atom types
10    AT2 = LABEL(NUM_ATOM)(1:2)
      NPOT = 0
!----look at H's first
1375    IF (AT2(1:1) .NE. 'H') GO TO 1340              ! H
           IT = ICON(1,NUM_ATOM)     ! number of atom to which H is attached
!---------Check for alcohol H                                       3/25/03
           IF (.NOT. (LABEL(IT)(1:1) .EQ. 'O' .AND.      &        ! 3/25/03
     &                NCON(IT) .EQ. 2)) GO TO 1378                ! 3/25/03
           NPOT = IQ_ALCOHOL(NUM_ATOM)                            ! 3/25/03
              IF (NPOT .NE. 0) GO TO 1303                         ! 3/25/03
1378       IF (LABEL(IT)(1:1) .EQ. 'C') THEN                       ! 3/25/03
              NC = NCON(IT)             ! it's C, how many atoms bonded to C?
              GO TO (101, 101, 103, 104), NC
103           NPOT = 1                         ! H bonded to sp2 C; param = 1
                 GO TO 1303
104           IF (IQ_CUBANE(IT) .EQ. 0) THEN  ! 4 bonded...test for cubane
                 NPOT = 2      ! no, H bonded to normal sp3 C; param = 2
              ELSE
                 NPOT = 47     ! H bonded to cubane skeleton; param = 2
              ENDIF            ! CAN BE CHANGED LATER; 11-12-02 changed = 47
                 GO TO 1303
           ELSE
              IF (LABEL(IT)(1:1) .NE. 'N') GO TO 101
                 NC = NCON(IT)             ! H is bonded to N
                 GO TO (101, 101, 1103, 1104), NC
!-----------H linked to ammonium N?
1104             npot = 157            ! H-N+                               6/10/05
                    go to 1303
!-----------Check for pyridinium-N(+)-H atom                                7/20/06
1103             if (iq_pyridinium_n(it) .eq. 85) then                    ! 7/20/06
                    npot = 86     ! H of pyridinium-N(+)-H, param = 86      7/20/06
                    go to 1303                                            ! 7/20/06
                 endif                                                    ! 7/20/06
!-----------Check on H linked to pyrrole-type N                              1/3/05
                 npot = iq_pyrrole_nh(it) ! it = number of 3-linked N       7/20/06
                 IF (NPOT .GT. 0) GO TO 1303  ! H of pyrrole, param = 74     1/3/06
!-----------Call IQ_HN_IMIDE to determine if the H is an imide NH,          12/3/03
!               H of a benzamide N or H or a hydantoin.                     12/3/03
                NPOT = IQ_HN_IMIDE(NUM_ATOM,IT)  ! # of H and 3-linked N   1/15/04
                 IF (NPOT .GT. 0) GO TO 1303! H of an imide, param = 39     11/15/03
!                                             H or an amide, param = 61     11/15/03
!                                             H of a hydantoin, param = 64  12/3/03
                 NPOT = AMINE(NUM_ATOM,IT)    ! # of H and 3-linked N       1/15/04
                 IF (NPOT .GT. 0) GO TO 1303  ! H of C-NH2, param = 3       1/15/04
           ENDIF
!---------Error call...can't ID the atom
101        PRINT 110, LABEL(NUM_ATOM), LABEL(IT)
110        FORMAT (1X,A,' bonded to ',A,', is an unknown type')
           GO TO 2002
!----Look at C's second
1340    IF (AT2(1:1) .NE. 'C') GO TO 2340
           GO TO (1001, 1002, 1003, 1004), NCON(NUM_ATOM)
1001          PRINT 800, LABEL(NUM_ATOM), NCON(NUM_ATOM)
800           FORMAT (1X,A5,' has',I2,' connections...does not compute')
              GO TO 2007
!------2 connections...possible -C#C-, =C= or -C#N
!------Check for internal alkyne, C-C#C-C, first                   12/31/03
1002       NPOT = IQ_ALKYNE(NUM_ATOM)                            ! 12/31/03
              IF (NPOT /= 0) GO TO 1303  ! alkyne C, param = 67    12/31/03
           DO 2701 I1=1,2
              I2 = ICON(I1,NUM_ATOM)              ! attached atom
                 IF (.NOT. (NCON(I2) .EQ. 1 .AND.   &
     &                      LABEL(I2)(1:1) .EQ. 'N')) THEN
                    NPOT = 34       ! nitrile C
                    GO TO 1303
                 ENDIF
2701       CONTINUE
           PRINT 2705, LABEL(NUM_ATOM), NCON(NUM_ATOM)  ! not the C of a nitrile
2705       FORMAT (1X,A5,' has',I2,' connections...what is it?')
           GO TO 2007
!---------3-linked C
1003       NSO1 = 0        ! # of connected 1-linked O's          ! 2/3/06
           nso2 = 0        ! # of connected 2-linked O's          ! 2/3/06
           DO 7338 I=1,3
              IT = ICON(I,NUM_ATOM)
              IF (LABEL(IT)(1:1) .EQ. 'O' .AND. NCON(IT) .EQ. 1) NSO1 = NSO1 + 1
              IF (LABEL(IT)(1:1) .EQ. 'O' .AND. NCON(IT) .EQ. 2) NSO2 = NSO2 + 1
7338       CONTINUE
           if ((nso1 .le. 1 .and. nso2 .eq. 0) .or.  &            ! 2/7/06
               (nso1 .eq. 0 .and. nso2 .eq. 1)) then              ! 2/7/06
              NPOT = 4      ! alkene, aldehyde, ketone, or          2/7/06
              GO TO 1303    !   vinyl ether C...sp2, param = 4      2/7/06
           ENDIF
           if (nso1 + nso2 .eq. 2) then  ! check for ester/anhydride C 2/3/06
              npot = iq_ester_c(num_atom)   ! param = 79          ! 2/3/06
              if (npot .gt. 0) go to 1303
           endif
           GO TO 2007        ! can't id 3-linked C                  2/3/06
!---------4-linked C
1004       IF (IQ_CUBANE(NUM_ATOM) .EQ. 0) THEN ! 4 bonded...test for cubane
                 NPOT = 5      ! no, alkane C...sp3; param = 5
              ELSE
                 NPOT = 46     ! cubane C; param = 5;
              ENDIF            ! CAN BE CHANGED LATER; 11-12-02 changed = 46
                 GO TO 1303
!
!----Look at O's third
2340    IF (AT2(1:1) .NE. 'O') GO TO 1024
           GO TO (1021, 1022), NCON(NUM_ATOM)
1021          IT = ICON(1,NUM_ATOM)     ! get ID of single connected atom
!----------Check if O of a NO3-                                    7/20/06
              if (label(it)(1:1) .eq. 'N' .and. ncon(it) .eq. 3) then  ! 7/20/06
                 no1 = 0  ! number of 1-linked O's               ! 7/20/06
                 do i=1,3                                        ! 7/20/06
                    it2 = icon(i,it)                             ! 7/20/06
                    if (label(it2)(1:1) .eq. 'O' .and. &         ! 7/20/06
                        ncon(it2) .eq.  1) no1 = no1 + 1         ! 7/20/06
                 enddo                                           ! 7/20/06
                 if (no1 .eq. 3) then                            ! 7/20/06
                    npot = 84  ! param = 84, O of NO3-           ! 7/20/06
                    go to 1303                                   ! 7/20/06
                 endif                                           ! 7/20/06
              endif                                              ! 7/20/06
!----------Check if linked to ester or anhydride carbonyl C         2/3/06
              if (label(it)(1:1) .eq. 'C' .and. ncon(it) .eq. 3 & ! 2/3/06
                  .and. iq_ester_c(it) .eq. 79) then              ! 2/3/06
                 npot = 80   ! param = 80, ester carbonyl O       ! 2/3/06
                 go to 1303                                       ! 2/3/06
              endif                                               ! 2/3/06
!----------Check for O of sulfonimine
           IF (LABEL(IT)(1:1) .EQ. 'S' .AND. NCON(IT) .EQ. 4) THEN
              ELEMENT = 'O'
              NPOT = IQ_SULFONIMINE(NUM_ATOM,ELEMENT) ! sulfonimine O
              IF (NPOT .GT. 0) GO TO 1303             ! param = 37
           ENDIF
!-------Check if the O is connected to a 3-linked C
         IF (LABEL(IT)(1:1) .EQ. 'C' .AND. NCON(IT) .EQ. 3) THEN
            NPOT = IQ_OC_IMIDE(IT)   ! imide or amide carbonyl O?      11/15/03
            IF (NPOT .GT. 0) GO TO 1303   ! O of imide, param = 41    11/15/03
!                                           O of amide, param = 63    11/15/03
!                                           O of hydantoin, param = 65 11/15/03
         ENDIF
!------------Check for carbonyl C with only 1-linked O, like ketone
              IF (LABEL(IT)(1:1) .NE. 'C') GO TO 7221
              NC = NCON(IT)
              IF (NC .NE. 3) GO TO 7221
              NSO = 0   ! # of 1-linked O's to the C
              DO 7224 I=1,3
                 IT2 = ICON(I,IT)
                 IF (NCON(IT2) .EQ. 1 .AND. LABEL(IT2)(1:1) .EQ. 'O')   &
     &                NSO = NSO + 1
7224          CONTINUE
              GO TO (7225, 2007, 101), NSO  ! to 2007 means 2 O's not recognized
7225              NPOT = 15        ! carbonyl O; param = 15
              GO TO 1303
!------------Check for nitroso ...N=O
7221          IF (LABEL(IT)(1:1) .EQ. 'N' .AND. NCON(IT) .EQ. 2) THEN
                 NPOT = 25      ! X-N=O
                 GO TO 1303
              ENDIF
!------------Check for N-oxide O of furoxan
              NPOT = IQ_FUROXAN_O1(NUM_ATOM)   ! param = 23
                 IF (NPOT .NE. 0) GO TO 1303
!
              IF (LABEL(IT)(1:1) .EQ. 'N' .AND.    &   ! N is linked to the O
     &               NCON(IT) .EQ. 3) THEN            ! & has 3 connections
                 N_OXY = 0
                 N_CARBS = 0
                 DO I=1,3                   ! count # of O's & C's linked
                    ITT = ICON(I,IT)        ! to the 2nd atom
                    IF (LABEL(ITT)(1:1) .EQ. 'C') N_CARBS = N_CARBS + 1
                    IF (LABEL(ITT)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
                 ENDDO
                 IF (N_OXY .EQ. 2 .OR. N_OXY .EQ. 3) THEN ! probably nitro O;
!                                                           C-NO2 or O-NO2
                   NPOT = IQ_NITRO_OX(IT)
                   IF( NPOT .EQ. 49) LABEL(NUM_ATOM)(1:2) = 'OX'    ! 11-12-02
                   IF (NPOT .EQ. 0) GO TO 101 ! can't identify the kind of O
                 ELSE
                   NPOT = 17    ! 2 C's & 1 O = a pyridine N-oxide O
                 ENDIF          ! ...param = 17
                 GO TO 1303
              ENDIF
!----Check for 2-linked O or nitrate ester
1022          NPOT = IQ_OXYGEN (NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
              NPOT = IQ_ALCOHOL (NUM_ATOM)                   ! 3/25/03
                 IF (NPOT .NE. 0) GO TO 1303                 ! 3/25/03
!----Check for the ring O of a furoxan...2-linked
              NPOT = IQ_FUROXAN_O2(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303       ! param = 22
!----Check for other 2-linked O's...C-O-C, X-O-B ***********************
              IT1 = ICON(1,NUM_ATOM)
              IT2 = ICON(2,NUM_ATOM)
              NC = 0
              NN = 0
              NB = 0
              IF (LABEL(IT1)(1:1) .EQ. 'B') NB = NB + 1
              IF (LABEL(IT2)(1:1) .EQ. 'B') NB = NB + 1
              IF (LABEL(IT1)(1:1) .EQ. 'C') NC = NC + 1
              IF (LABEL(IT2)(1:1) .EQ. 'C') NC = NC + 1
              IF (LABEL(IT1)(1:1) .EQ. 'N') NN = NN + 1
              IF (LABEL(IT2)(1:1) .EQ. 'N') NN = NN + 1
!----Check for "ether" O of an ester or anhydride               2/3/06
              If (nc .eq. 2 .and. (iq_ester_c(it1) .eq. 79 &  ! 2/3/06
                  .or. iq_ester_c(it2) .eq. 79)) then         ! 2/3/06
                 npot = 81   ! "ether" O of an ester            2/3/06
                 go to 1303                                   ! 2/3/06
              endif                                           ! 2/3/06
              IF (NC .EQ. 2 .OR. (NC .EQ. 1 .AND. NN .EQ. 1) .OR.   &
     &           (NC .EQ. 1 .AND. NB .EQ. 1))  THEN
                     NPOT = 14      ! ether O; C-O-C, C-O-N or C-O-B
              ELSE
                  NPOT = 999        ! unknown
              ENDIF
              GO TO 1303
!
!----Look at N's fourth
1024    IF (AT2(1:1) .NE. 'N') GO TO 2002
           GO TO (1031, 1032, 1033, 2034), NCON(NUM_ATOM)
!------------Check for terminal azide N                          ! 3/5/03
1031          NPOT = IQ_AZIDE_N (NUM_ATOM)                       ! 3/5/03
                 IF (NPOT .NE. 0) GO TO 1303                     ! 3/5/03
!------------Possible nitrile N...check to see it there's a 2-linked C attached
              IT2 = ICON(1,NUM_ATOM)   ! IT2 = attached atom     ! 3/5/03
              IF (NCON(IT2) .EQ. 2 .AND. LABEL(IT2)(1:1) .EQ. 'C') THEN
                 NPOT = 33       ! nitrile N, param = 33
              ELSE
                 NPOT = 999      ! flag as unknown...1-linked N but not a nitrile
              ENDIF
              GO TO 1303
!----N with 4 connections....like NH4+
2034          NPOT = 999         ! flag as unknown
                 GO TO 1303
!----N with 2 connections
!------------Check for 2-linked N of amino-imine                1/20/06
1032          npot = iq_aminoimine_2(num_atom)                ! 1/23/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 76   ! 1/20/06
!------------Check for 2-linked N of thiazole                   1/31/06
              npot = iq_thiazole_n(num_atom)                  ! 1/31/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 78   ! 1/31/06
!------------Check for 2-linked N (N3) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n3(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 72    ! 1/1/06
!------------Check for 2-linked N (N2) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n2(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 71    ! 1/1/06
!------------Check for 2-linked N (N4) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n4(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 73    ! 1/1/06
!------------Check for 2-linked N of 123triazole                 12/30/05
              npot = iq_123triazole_2(num_atom)                ! 12/30/05
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 27    ! 12/30/05
!------------Check for internal azide N's                        ! 3/5/03
              NPOT = IQ_AZIDE_N (NUM_ATOM)                       ! 3/5/03
                 IF (NPOT .NE. 0) GO TO 1303                     ! 3/5/03
              NPOT = IQ_N_AZAPENTALENE(NUM_ATOM)                 ! 3/5/03
                 IF (NPOT .NE. 0) GO TO 1303
              NPOT = IQ_FUROXAN_N2(NUM_ATOM)   ! param = 20
                 IF (NPOT .NE. 0) GO TO 1303
!------------Check for diazo N....C-N=N-C                         2/10/06
              NPOT = IQ_diazo_n(NUM_ATOM)      ! param = 82     ! 2/10/06
                 IF (NPOT .NE. 0) GO TO 1303                    ! 2/10/06
!
1037          N_OXY = 0       ! 2 connections, how many O's are attached
              N_CARBS = 0     ! # of C's attached
              N_NIT = 0       ! # of N's attached
              N_HYD = 0       ! # of H's attached
              N_F = 0         ! # of F's attached
              N_S_4 = 0       ! # of 4-linked S's attached         ! 12/11/00
              DO I=1,2
                 N_ADJ = ICON(I,NUM_ATOM)
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'O') THEN
                     N_OXY = N_OXY + 1
                     NUM_O(N_OXY) = N_ADJ
                 ENDIF
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'C') N_CARBS = N_CARBS + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'N') N_NIT = N_NIT + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'H') N_HYD = N_HYD + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'F') N_F = N_F + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'S' .AND.        &            ! 12/11/00
     &                 NCON(N_ADJ) .EQ. 4)  N_S_4 = N_S_4 + 1            ! 12/11/00
              ENDDO
!------------Check on a 4-linked S attached to the 2-linked N.           ! 12/11/00
              IF (N_S_4 .GT. 0) THEN
                  ELEMENT = 'N'
                  NPOT = IQ_SULFONIMINE(NUM_ATOM, ELEMENT)  ! it's a N
                  IF (NPOT .NE. 0) GO TO 1303                            ! 12/11/00
              ENDIF                                                      ! 12/11/00
              IF (N_OXY .EQ. 1 .AND. NCON(NUM_O(N_OXY)) .EQ. 1) THEN
                  NPOT = 24        ! N of X-N=O
                  GO TO 1303
              ENDIF
!
              IF (N_OXY .EQ. 2) THEN
                  N = MIN0(NCON(NUM_O(1)), NCON(NUM_O(2)))
                  IF (N .EQ. 1) THEN
                     NPOT = 24     ! N of -O-N=O
                     GO TO 1303
                  ENDIF
              ENDIF
!
              IF ((N_CARBS .EQ. 2) .OR. (N_CARBS .EQ. 1 .AND.   &
     &             N_NIT .EQ. 1)) THEN
                 NPOT = 16    ! pyridine N...param = 16
                 GO TO 1303
              ENDIF
              GO TO 2007
!----N with 3 connections
!------------Check for 3-linked N of aminoimine                 1/20/06
1033          npot = iq_aminoimine_3(num_atom)                ! 1/23/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 76   ! 1/20/06
!------------Check for 3-linked N (N1) of 1234-tetrazole         1/1/06
              npot = iq_1234tetrazole_n1(num_atom)             ! 1/1/06
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 70    ! 1/1/06
!------------Check for 3-linked N of 123triazole                 12/30/05
              npot = iq_123triazole_3(num_atom)                ! 12/30/05
                 IF (NPOT .NE. 0) GO TO 1303   ! param = 26    ! 12/30/05
!------------Check for nitrate ester N and nitramine nitro N       2/24/03
             NPOT = IQ_NITRATE_N (NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
             NPOT = IQ_N_AZAPENTALENE(NUM_ATOM)
                 IF (NPOT .NE. 0) GO TO 1303
!------------Check if central N of imide, benzamide N or a hydantoin    12/3/03
              NPOT = IQ_HN_IMIDE(NUM_ATOM,NUM_ATOM)       ! same atom #'s = a N
                 IF (NPOT .GT. 0) GO TO 1303! N of an imide, param = 40 12/3/03
!                                             N or an amide, param = 62 12/3/03
!                                             N of a hydantoin, param = 67 12/3/03
              NPOT = IQ_FUROXAN_N3(NUM_ATOM)   ! param = 21
                 IF (NPOT .NE. 0) GO TO 1303
1035          N_OXY = 0       ! 3 connections, how many O's are attached
              N_OXY_1 = 0     ! # of 1-linked O's attached        ! 6/19/03
              N_CARBS = 0     ! # of C's attached
              N_CARBS_3 = 0   ! # of 3-linked C's attached        ! 6/19/03
              N_NIT = 0       ! # of N's attached
              N_HYD = 0       ! # of H's attached
              N_F = 0         ! # of F's attached
              DO I=1,3
                 N_ADJ = ICON(I,NUM_ATOM)
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'O') N_OXY = N_OXY + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'O' .AND.        &     ! 6/19/03
     &               NCON(N_ADJ) .EQ. 1) N_OXY_1 = N_OXY_1 + 1   ! 6/19/03
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'C') THEN
                    N_CARBS = N_CARBS + 1
                    NC_ADJ = N_ADJ                ! 11-12-02
                 END IF
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'C' .AND.     &          ! 6/19/03
     &               NCON(N_ADJ) .EQ. 3) N_CARBS_3 = N_CARBS_3 + 1  ! 6/19/03
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'N') N_NIT = N_NIT + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'H') N_HYD = N_HYD + 1
                 IF (LABEL(N_ADJ)(1:1) .EQ. 'F') N_F = N_F + 1
              ENDDO
              IF (N_OXY .EQ. 2 .OR. N_OXY .EQ. 3) THEN   ! 11-12-02 included ONO2
                 IF (N_CARBS .EQ. 1) THEN
                    IF ( IQ_CUBANE(NC_ADJ) .EQ. 1) THEN
                       NPOT = 48    ! nitro N of nitrocubane   11-12-02
                       GO TO 1303
                    ELSE
                       NPOT = 7     ! nitro N...param = 7 CNO2 except nitrocubane 11-12-02
                       GO TO 1303
                    END IF
                 ELSE
                    NPOT = 7     ! nitro N...param = 7 NNO2 or ONO2 11-12-02
                    GO TO 1303
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
!----------Check for pyridinium-N....[Csp2]2-N(+)-H              ! 7/20/06
              IF (N_CARBS_3 .EQ. 2 .AND. n_hyd .EQ. 1) THEN      ! 7/20/06
                 NPOT = iq_pyridinium_n(num_atom)   ! pyridinium N(+)-H 7/20/06
                 if (npot .gt. 0) GO TO 1303        ! N param = 85      7/20/06
              ENDIF                                              ! 7/20/06
!----Test for N of [Csp3]2NH                                   1/15/04
              NPOT = AMINE(NUM_ATOM, NUM_ATOM)               ! 1/15/04
                 IF (NPOT .NE. 0) GO TO 1303   ! N PARAM = 68  1/15/04
              IF (N_HYD .GE. 1 .and. n_carbs .ge. 1) THEN  ! must have one
                 NPOT = 10  ! NH2 or NH param = 10         !   C-N bond
                 GO TO 1303
              ENDIF
              IF (N_F .GE. 1) THEN
                 NPOT = 45  ! NHF or NF2 param = 45, 11-12-02
                 GO TO 1303
              END IF
              IF (N_CARBS_3 .EQ. 2 .AND. N_OXY_1 .EQ. 1) THEN    ! 8/11/03
                 NPOT = 18      ! pyridine N-oxide N...param = 18
                 GO TO 1303
              ENDIF                                              ! 6/19/03
       GO TO 2007
!----Look at Br's 5th
2002   IF (AT2 .NE. 'BR') GO TO 2008
          NPOT = 12                  ! Br, param = 12
          GO TO 1303
!----Look at F's 6th
2008   IF (AT2(1:1) .NE. 'F') GO TO 3009
          N_ADJ = ICON(1,NUM_ATOM)
          IF (LABEL(N_ADJ)(1:1) .EQ. 'C') THEN
             NPOT = 53               ! F of C-F, param = 53
          ELSE
             NPOT = 13               ! F of N-F, param = 13
          ENDIF
          GO TO 1303
!----Look at B's 7th (O2-B-C)
3009   IF (AT2(1:1) .NE. 'B') GO TO 4007
              IF (NCON(NUM_ATOM) .NE. 3) GO TO 4007
              NC = 0
              NO = 0
              DO 3015 I=1,3
                 IT = ICON(I,NUM_ATOM)
                 IF (LABEL(IT)(1:1) .EQ. 'O') NO = NO + 1
                 IF (LABEL(IT)(1:1) .EQ. 'C') NC = NC + 1
3015          CONTINUE
              IF (NC .EQ. 1 .AND. NO .EQ. 2)  THEN
                  NPOT = 35      ! boron...O2-B-C, param = 35
                  GO TO 1303
              ENDIF
!----Look at S's 8th
4007   IF (AT2(1:1) .NE. 'S') GO TO 3035                  ! 11/6/03
          if (ncon(num_atom) .ne. 2) go to 4009           ! 1/31/06
             npot = iq_thiazole_s(num_atom)               ! 1/31/06
             if (npot .ne. 0) go to 1303  ! param = 77      1/31/06
4009      IF (NCON(NUM_ATOM) .NE. 4) GO TO 3035           ! 1/31/06
          ELEMENT = 'S'
          NPOT = IQ_SULFONIMINE (NUM_ATOM,ELEMENT) ! sulfonimine S, param = 38
          IF (NPOT .NE. 0) GO TO 1303                     ! 11/15/03
!----Look at I's 9th                                      ! 11/6/03
3035   IF (AT2(1:1) .NE. 'I') GO TO 2007                  ! 11/6/03
          NPOT = 60           ! I, param = 60             ! 11/6/03
       GO TO 1303                                         ! 11/6/03
!
2007   PRINT 2004, LABEL(NUM_ATOM)
2004   FORMAT (/' Atom ',A5,' is unknown...type set to 999')
1034      NPOT = 999     ! flag as unknown
!
1303   RETURN
       END FUNCTION NPOT
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!
!----IQ_SULFONIMINE...identify S, O and N of C-S(O2)-N=C (sulfonimine)
!
      FUNCTION IQ_SULFONIMINE (IS,ELEMENT)
!
!----ELEMENT = 4-linked S, 1-linked O or 2-linked N
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_SULFONIMINE, IS
      CHARACTER (LEN=1) :: ELEMENT

!      DIMENSION NUM_N2(3)
! NAME is changed to LABEL
!
      INTEGER :: I, I2, IC, J, KN, NC, NO_1
      INTEGER :: NN_2, NS, K, ISS, NUM_N2(3)
!
      IQ_SULFONIMINE = 0
!
!----Sulfur???   must be 4-linked
      IF (ELEMENT .EQ. 'S' .AND. NCON(IS) .EQ. 4) THEN
         ISS = IS                                 ! ISS = S #
         GO TO 20
      ENDIF
!----Oxygen???   must be 1-linked
      IF (ELEMENT .EQ. 'O' .AND. NCON(IS) .EQ. 1) THEN
         ISS = ICON(1,IS)         ! atom linked to O should be S
         IF (LABEL(ISS)(1:1) .EQ. 'S' .AND. NCON(ISS) .EQ. 4) GO TO 20
      ENDIF
!----Nitrogen???  must be 2-linked
      IF (ELEMENT .EQ. 'N' .AND. NCON(IS) .EQ. 2) THEN
         DO K=1,2               ! one of the 2-linked atoms
            ISS = ICON(K,IS)    ! must be a 4-linked S
            IF (LABEL(ISS)(1:1) .EQ. 'S' .AND. NCON(ISS) .EQ. 4) GO TO 20
         ENDDO
      ENDIF
      RETURN       ! can't make sense of anything
!
20    NS = NCON(ISS)          ! ISS is the 4-linked S
!
!----Count # C's, # 1 and # 2-linked O's and 3-linked N's bonded
!     to the 4-linked S
      NC = 0                            ! # C's linked to S    ! 3/20/00
      NO_1 = 0                          ! # 1-linked O's       ! 3/18/00
      NN_2 = 0                          ! # 2-linked N's       ! 11/3/00
      DO 200 I=1,NS
         I2 = ICON(I,ISS)                                       ! 3/18/00
         IF (LABEL(I2)(1:1) .EQ. 'C') NC = NC + 1               ! 3/18/00
         IF (LABEL(I2)(1:1) .EQ. 'O' .AND. NCON(I2) .EQ. 1)    &
     &      NO_1 = NO_1 + 1                                    ! 3/18/00
         IF (LABEL(I2)(1:1) .EQ. 'N' .AND. NCON(I2) .EQ. 2) THEN! 11/3/00
            NN_2 = NN_2 + 1         ! looking for 2-linked N's
            NUM_N2(NN_2) = I2
         ENDIF                                                 ! 11/3/00
200   CONTINUE                                                 ! 3/18/00
!
      IF (NC .NE. 1 .OR. NO_1 .NE. 2 .OR. NN_2 .NE. 1) RETURN
!
!-------Sulfonimine...C-SO2-N=C ??
222         KN = NUM_N2(1)       ! check imine N &
            DO 224 J=1,2         !   look for 3-linked C
               IC = ICON(J,KN)                                 ! 11/3/00
               IF (LABEL(IC)(1:1) .EQ. 'S') GO TO 224           ! 11/3/00
               IF (LABEL(IC)(1:1) .EQ. 'C' .AND.      &         ! 11/3/00
     &             NCON(IC) .EQ. 3) GO TO 226      ! ok        ! 11/3/00
               RETURN                              ! not ok    ! 11/3/00
224         CONTINUE                                           ! 11/3/00
226         IF (ELEMENT .EQ. 'O') IQ_SULFONIMINE = 37  ! param = 37
            IF (ELEMENT .EQ. 'N') IQ_SULFONIMINE = 38  ! param = 38
            IF (ELEMENT .EQ. 'S') IQ_SULFONIMINE = 36  ! param = 36
      RETURN
      END FUNCTION IQ_SULFONIMINE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!----IQ_ALKYNE   determine is a 2-linked C is an internal      12/31/03
!                alkyne C                                      12/31/03
!
      FUNCTION IQ_ALKYNE (NUM_ATOM)                          ! 12/31/03
!
!---- NUM_ATOM = a 2-linked C                                ! 12/31/03
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
!
      INTEGER :: I, L, NC, NC_TOTAL, NUM_ATOM, NUM_C(2)      ! 12/31/03
      INTEGER :: IQ_ALKYNE                                   ! 12/31/03
!
      IQ_ALKYNE = 0        ! intialize to NO                   12/31/03
!
!----Count # of 2-linked C's and total C's                     12/31/03
      NC = 0                                                 ! 12/31/03
      NC_TOTAL = 0                                           ! 12/31/03
      DO 100 I=1,2                                            ! 12/3/03
         L = ICON(I,NUM_ATOM)                                 ! 12/3/03
         IF (LABEL(L)(1:1) .EQ. 'C') NC_TOTAL = NC_TOTAL + 1 ! 12/31/03
         IF (LABEL(L)(1:1) == 'C' .AND. NCON(L) == 2) THEN   ! 12/31/03
            NC = NC + 1                                      ! 12/31/03
            NUM_C(NC) = L                                    ! 12/31/03
         ENDIF                                               ! 12/31/03
100   CONTINUE                                                ! 12/3/03
!
      IF (NC /= 1 .OR. NC_TOTAL /= 2) RETURN                 ! 12/31/03
!----Check the other C in NUM_C to see if it is 2-linked       12/31/03
      L = NUM_C(1)                                           ! 12/31/03
      IF (NCON(L) /= 2) RETURN                               ! 12/31/03
      IQ_ALKYNE = 67        ! internal alkyne, param = 67      12/31/03
      RETURN                                                 ! 12/31/03
!
     END FUNCTION IQ_ALKYNE                                  ! 12/3/03

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!----N_CARBONYL...determine if a 3-linked C is a carbonyl C     12/3/03
!
      FUNCTION N_CARBONYL (NUM_ATOM)                          ! 12/3/03
!
!----If NUM_ATOM = a 3-linked C, N_CARBONYL = 0/1 for no/yes    12/3/03
!       itis a carbonyl                                         12/3/03
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
!
      INTEGER :: I, L, NUM_ATOM, N_CARBONYL                   ! 12/3/03
!
      N_CARBONYL = 0        ! intialize to NO                   12/3/03
!
!----Count # of 3-linked and total C's on the N                 12/3/03
      DO 100 I=1,3                                            ! 12/3/03
         L = ICON(I,NUM_ATOM)                                 ! 12/3/03
         IF (LABEL(L)(1:1) .EQ. 'O' .AND. NCON(L) .EQ. 1) THEN ! 12/3/03
             N_CARBONYL = 1     ! Yes, a carbonyl             ! 12/3/03
         ENDIF                                                ! 12/3/03
100   CONTINUE                                                ! 12/3/03
!
     END FUNCTION N_CARBONYL                                  ! 12/3/03
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!----N_3LINK...count # of 3-linked N's attached to NUM_ATOM     12/3/03
!
      FUNCTION N_3LINK (NUM_ATOM)                             ! 12/3/03
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
!
      INTEGER :: I, L, NUM_ATOM, N_3LINK                      ! 12/3/03
!
      N_3LINK = 0                                             ! 12/3/03
!
!----Count # of 3-linked N's                                    12/3/03
      DO 100 I=1,NCON(NUM_ATOM)                               ! 12/3/03
         L = ICON(I,NUM_ATOM)                                 ! 12/3/03
         IF (LABEL(L)(1:1) .EQ. 'N' .AND. NCON(L) .EQ. 3) THEN ! 12/3/03
             N_3LINK = N_3LINK + 1   ! yes, a 3-linked N        12/3/03
         ENDIF                                                ! 12/3/03
100   CONTINUE                                                ! 12/3/03
!
     END FUNCTION N_3LINK                                     ! 12/3/03
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!----AMINE...                                                   1/15/04
!
!        N,C---NH2               [Csp3]2---N-H                  1/15/04
!            H = 3              N = 69     H = 68               1/15/04
!
      FUNCTION AMINE(NUM_ATOM,IT)                             ! 1/15/04
!
!----If NUM_ATOM .ne. IT, NUM_ATOM = H; IT = 3-linked N         1/15/04
!    If NUM_ATOM = IT, IT and NUM_ATOM = 3-linked N             1/15/04
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
!
      LOGICAL :: TYPE_ATOM                                    ! 1/15/04
      INTEGER :: I, IT, L, NC, NC4, NH, NUM_ATOM, AMINE       ! 1/15/04
      integer :: nn                                           ! 1/25/06
!
      AMINE = 0                                               ! 1/15/04
      TYPE_ATOM = .true.        ! set signal to H               1/15/04
      IF (NUM_ATOM == IT) TYPE_ATOM = .false. ! set signal to N 1/15/04
!----Count # of C's and 4-linked C's & H's on the N             1/15/04
!
      nn = 0                                                  ! 1/26/06
      NC = 0   ! # of C's                                       1/15/04
      NC4 = 0  ! # of 4-linked C's                              1/15/04
      NH = 0   ! # of H's                                       1/15/04
      DO I=1,3                                                ! 1/15/04
         L = ICON(I,IT)   ! look at connections to N            1/15/04
         if (label(l)(1:1) .eq. 'N') nn = nn + 1              ! 1/25/06
         IF (LABEL(L)(1:1) .EQ. 'C') NC = NC + 1              ! 1/15/04
         IF (LABEL(L)(1:1) .EQ. 'H') NH = NH + 1              ! 1/15/04
         IF (LABEL(L)(1:1) == 'C' .AND. NCON(L) == 4) NC4 = NC4 + 1 ! 1/15/04
      ENDDO                                                   ! 1/15/04
      if (NH == 2 .AND. (NC == 1 .or. nn .eq. 1)) THEN        ! 7/20/06
            if (type_atom) AMINE = 3        ! R-NH2, param = 3 for H  7/20/06
            if (.not. type_atom) amine = 10 ! R-NH2, param = 10 for N 7/20/06
      ELSE                                                    ! 1/15/04
         IF (NH == 1 .AND. NC4 == 2) THEN                     ! 1/15/04
            IF (TYPE_ATOM) AMINE = 68         ! for H           1/15/04
            IF (.not. TYPE_ATOM) AMINE = 69   ! for N           1/15/04
         ENDIF                                ! [Csp3]2NH       1/15/04
      ENDIF                                                   ! 1/15/04
      RETURN                                                  ! 1/15/04
     END FUNCTION AMINE                                       ! 1/15/04
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!----IQ_HN_IMIDE... identify H and N of imides,                 12/3/03
!                   benzamides and hydantoins                   12/3/03
!
!                    40                     62                  11/15/03
!               --C---N---C--      -phenyl---N---C--            12/3/03
!                 "   |   "                  |   "              11/15/03
!                 O   H   O                  H   O              11/15/03
!                41  39  41                 61  63              11/15/03
!
!                    66      66                                 11/15/03
!               --C---N---C---N---C--                           12/3/03
!                 "   |   "   |                                 12/3/03
!                 O   H   O   H                                 12/3/03
!                65  64  65  64                                 12/3/03
!
      FUNCTION IQ_HN_IMIDE(NUM_ATOM,IN)                       ! 12/3/03
!
!----If NUM_ATOM .NE. IN, atom is a H; if NUM_ATOM = IN, atom is a N
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
!
      IMPLICIT NONE
      LOGICAL :: HYDROGEN
      INTEGER :: IQ_HN_IMIDE, IQ_OC_IMIDE, NUM_ATOM, IN       ! 12/3/03
      INTEGER :: I, K, L, M, NC_3, NC_TOTAL, NN, NO           ! 12/3/03
      INTEGER :: N_CARB, N_CARBONYL, N_3LINK, N_3LINKSUM      ! 12/3/03
      INTEGER :: NUM_C(3), NUM_CARB(3), NUM_N(3)              ! 12/3/03
!
      IQ_HN_IMIDE = 0
      HYDROGEN = .true.
!
      IF (NUM_ATOM .EQ. IN) HYDROGEN = .false.
!
!----Count # of 3-linked, total and carbonyl C's on the N       12/3/03
      NC_3 = 0                                                ! 12/3/03
      NC_TOTAL = 0                                            ! 12/3/03
      N_CARB = 0
      DO 100 I=1,3           ! count # C's on N                 12/3/03
         L = ICON(I,IN)
         IF (LABEL(L)(1:1) .EQ. 'C') NC_TOTAL = NC_TOTAL + 1  ! 12/3/03
         IF (LABEL(L)(1:1) .EQ. 'C' .AND. NCON(L) .EQ. 3) THEN
            NC_3 = NC_3 + 1                                   ! 12/3/03
            NUM_C(NC_3) = L         ! id #'s of 3-linked C's    12/3/03
            IF (N_CARBONYL(L) == 1) THEN   ! test for           12/3/03
               N_CARB = N_CARB + 1         ! carbonyl C &       12/3/03
               NUM_CARB(N_CARB) = L        ! record number      12/3/03
            ENDIF                                             ! 12/3/03
         ENDIF                                                ! 12/3/03
100   CONTINUE                                                ! 12/3/03
!
      IF (N_CARB == 0) RETURN    ! not a N-C(=O)               12/3/03
!
!----This is a test for an amide...phenyl-NH-C(=O)-C              12/3/03
      IF (N_CARB == 1 .AND. N_3LINK(NUM_CARB(1)) == 1  &      ! 12/3/03
                .AND. NC_3 == 2) THEN                         ! 12/3/03
         IF (HYDROGEN) THEN            ! it is an amide         11/15/03
            IQ_HN_IMIDE = 61           ! H param = 61           11/15/03
         ELSE                                                 ! 11/15/03
            IQ_HN_IMIDE = 62           ! N param = 62           11/15/03
         ENDIF                                                ! 11/15/03
         RETURN                                               ! 11/15/03
      ENDIF                                                   ! 11/15/03
!
!----Test for imides                                            12/3/03
      IF (N_CARB == 2) THEN   ! how many 3-linked N's           12/3/03
         N_3LINKSUM = 0                                       ! 12/3/03
         DO I=1,2                                             ! 12/3/03
            N_3LINKSUM = N_3LINKSUM + N_3LINK(NUM_CARB(I))    ! 12/3/03
         ENDDO                                                ! 12/3/03
         IF (N_3LINKSUM == 2) THEN  ! an imide                ! 12/3/03
            IF (HYDROGEN) THEN            ! it is an imide      12/3/03
               IQ_HN_IMIDE = 39           ! H param = 39        12/3/03
            ELSE                                              ! 12/3/03
               IQ_HN_IMIDE = 40           ! N param = 40        12/3/03
            ENDIF                                             ! 12/3/03
            RETURN                                            ! 12/3/03
         ENDIF                                                ! 12/3/03
!
!----Test for central N of hydantoin, N_3LINKSUM = 3?           12/3/03
         IF (N_3LINKSUM == 3) THEN    ! yes                     12/3/03
            IF (HYDROGEN) THEN            ! it is an N          12/3/03
               IQ_HN_IMIDE = 64           ! H param = 64        12/3/03
            ELSE                                              ! 12/3/03
               IQ_HN_IMIDE = 66           ! N param = 66        12/3/03
            ENDIF                                             ! 12/3/03
            RETURN                                            ! 12/3/03
         ENDIF                                                ! 12/3/03
      ENDIF                                                   ! 12/3/03
!
!----Finally test for end N of hydantoin                       12/3/03
      IF (.not. (N_CARB == 1 .AND.  &                        ! 12/3/03
              N_3LINK(NUM_CARB(1)) == 2)) RETURN  ! not end N  12/3/03
!----Look at 2 N's on carbonyl and get N id's                  12/3/03
      NN = 0    ! count of # 3-linked N's on carbonyl          12/3/03
      DO 200 I=1,3                                           ! 12/3/03
         K = ICON(I,NUM_CARB(1))                             ! 12/3/03
         IF (LABEL(K)(1:1) /= 'N' .OR. NCON(K) /= 3) GO TO 200 ! 12/3/03
            NN = NN + 1                                      ! 12/3/03
            NUM_N(NN) = K   ! id of 3-linked N                 12/3/03
200   CONTINUE                                               ! 12/3/03
      IF (NN /= 2) RETURN    ! must have two 3-linked N's      12/3/03
      N_CARB = 0                                             ! 12/3/03
      DO 300 I=1,NN                                          ! 12/3/03
         DO 290 L=1,NCON(NUM_N(I))                           ! 12/3/03
            M = ICON(L,NUM_N(I))                             ! 12/3/03
            IF (LABEL(M)(1:1) /= 'C' .OR. NCON(M) /= 3) GO TO 290 ! 12/3/03
            N_CARB = N_CARB + N_CARBONYL(M)                  ! 12/3/03
290      CONTINUE                                            ! 12/3/03
300   CONTINUE                                               ! 12/3/03
      IF (N_CARB /= 3) RETURN                                ! 12/3/03
      IF (HYDROGEN) THEN             ! it is an N              12/3/03
          IQ_HN_IMIDE = 64           ! H param = 64            12/3/03
      ELSE                                                   ! 12/3/03
          IQ_HN_IMIDE = 66           ! N param = 66            12/3/03
      ENDIF                                                  ! 12/3/03
      RETURN                                                 ! 12/3/03
!
      END FUNCTION IQ_HN_IMIDE                               ! 12/3/03
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!----IQ_OC_IMIDE... identify O of imide, benzamide or hydantoin 11/15/03
!
!                    40                     62                  11/15/03
!               --C---N---C--       -phenyl--N--C--             12/3/03
!                 "   |   "                  |  "               11/15/03
!                 O   H   O                  H  O               11/15/03
!                41  39  41                 61 63               11/15/03
!
!                    66      66                                 11/15/03
!               --C---N---C---N---C--                           12/3/03
!                 "   |   "   |                                 12/3/03
!                 O   H   O   H                                 12/3/03
!                65  64  65  64                                 12/3/03
!
      FUNCTION IQ_OC_IMIDE(IN)
!
!----Have a 1-linked O; IN is the 3-linked C                    11/15/03
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_OC_IMIDE, IN
!
!
      INTEGER :: I, L, NUM_N(3), NC, NN, NO, IQ_HN_IMIDE      ! 11/15/03
!
      IQ_OC_IMIDE = 0
!
!----Count # of atoms linked to the 3-linked C                  11/15/03
      NC = 0          ! count # C's
      NN = 0          ! count # 3-linked N's
      NO = 0          ! count # of 1-linked O's
      DO 100 I=1,3
         L = ICON(I,IN)
         IF (LABEL(L)(1:1) .EQ. 'N' .AND. NCON(L) .EQ. 3) THEN
            NN = NN + 1
            NUM_N(NN) = L                                     ! 11/15/03
         ENDIF
         IF (LABEL(L)(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(L)(1:1) .EQ. 'O' .AND. NCON(L) .EQ. 1) NO = NO + 1
100   CONTINUE
!----Check for C-C(=O)-N and N-C(=O)-N                          11/15/03
      IF (NC .EQ. 1 .AND. NN .EQ. 1 .AND.  &  ! atom numbers for         12/3/03
     &    NO .EQ. 1) GO TO 120                ! amide, imide, hydantoin  12/3/03
      IF (NC .EQ. 0 .AND. NN .EQ. 2 .AND.  &  ! atom numbers     12/3/03
     &    NO .EQ. 1) GO TO 120                ! for hydantoin    12/3/03
      RETURN                                                   ! 12/3/03
120   DO I=1,NN                                                ! 12/3/03
!----Call IQ_HN_IMIDE to test if the N is an imide or amide N   11/15/03
      L = IQ_HN_IMIDE (NUM_N(I),NUM_N(I))                     ! 11/15/03
      IF (L .EQ. 40) THEN         ! imide N param = 40          11/15/03
         IQ_OC_IMIDE = 41         ! imide O param = 41          11/15/03
         RETURN                                               ! 11/15/03
      ENDIF                                                   ! 11/15/03
      IF (L .EQ. 62) THEN      ! amide N param = 62           ! 11/15/03
         IQ_OC_IMIDE = 63      ! amide O param = 63           ! 11/15/03
         RETURN                                               ! 11/15/03
      ENDIF                                                   ! 11/15/03
      IF (L .EQ. 66) THEN      ! hydantoin N = 66               11/15/03
         IQ_OC_IMIDE = 65      ! hydantoin O = 65               11/15/03
         RETURN
      ENDIF                                                   ! 11/15/03
      ENDDO                                                   ! 11/15/03
      END FUNCTION IQ_OC_IMIDE
!
!----------------------------------------------------------------
!----Function to determine if a single-linked oxygen is an
!     N-oxide furoxan O
!
      FUNCTION IQ_FUROXAN_O1(I1)    ! I1 is a single-linked O
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_FUROXAN_O1, I1
!
!NAME is changed to LABEL
      INTEGER :: IQ_FUROXAN_N3, IQ_FUROXAN_N2, IQ_FUROXAN_O2
      INTEGER :: I, I2, IO, I_2, I_3, I_C, I_O
      INTEGER :: N, NN, NO, NC, NLINKED, NCON_O
!
      IQ_FUROXAN_O1 = 0     ! no
!
!----Must be linked to a N, which is linked to a C and O & the
!     O must be linked to another N
!
!                 C---C
!                 "   "
!          20 --> N   N <-- 21
!                  \ / \
!           22 -->  O   O <-- 23
!
      N = ICON(1,I1)                         ! is the connected
         IF (LABEL(N)(1:1) .NE. 'N') RETURN   ! atom a N with
         IF (NCON(N) .NE. 3) RETURN          ! NCON = 3?
      NO = 0
      NC = 0
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,N)
         IF (LABEL(I2)(1:1) .EQ. 'O') NO = NO + 1
            IF (I2 .NE. I1) IO = I2     ! IO is the "other" O
         IF (LABEL(I2)(1:1) .EQ. 'C') NC = NC + 1
100   CONTINUE
      IF (NC .NE. 1 .OR. NO .NE. 2) RETURN
!----Test the other O...must have 2 N links
      IF (NCON(IO) .NE. 2) RETURN
      DO 120 I=1,2
         IF (LABEL(ICON(I,IO))(1:1) .NE. 'N') RETURN
120   CONTINUE
      IQ_FUROXAN_O1 = 23     ! param = 23 for N-oxide O
      RETURN
!
!----Entry point for the ring O in a furoxan
      ENTRY IQ_FUROXAN_O2(I1)   ! I1 is a double linked O
!
      IQ_FUROXAN_O2 = 0
!
      NN = 0       ! count of # of N's linked to O
      NLINKED = 0  ! total number of things linked to the N's
      DO 200 I=1,2
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .NE. 'N') GO TO 200
            NN = NN + 1
            NLINKED = NLINKED + NCON(I2)
            IF (NCON(I2) .EQ. 2) I_2 = I2   ! the 2-linked ring N
            IF (NCON(I2) .EQ. 3) I_3 = I2   ! the 3-linked ring N
200   CONTINUE
      IF (NN .NE. 2 .OR. NLINKED .NE.5) RETURN
      NO = 0
      NC = 0
      DO 300 I=1,2                         ! check the 2-linked N
         IF (LABEL(ICON(I,I_2))(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(ICON(I,I_2))(1:1) .EQ. 'O') NO = NO + 1
300   CONTINUE
      IF (NC.NE. 1 .OR. NO .NE. 1) RETURN
      NO = 0
      NC = 0
      DO 400 I=1,3                         ! check the 3-linked N
         IF (LABEL(ICON(I,I_3))(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(ICON(I,I_3))(1:1) .EQ. 'O') NO = NO + 1
400   CONTINUE
      IF (NC.NE. 1 .OR. NO .NE. 2) RETURN
      IQ_FUROXAN_O2 = 22     ! param = 22 for in-ring O
      RETURN
!
!----Entry point for the 2-linked ring N in a furoxan
!
      ENTRY IQ_FUROXAN_N2(I1)   ! I1 is a double linked N
!
      IQ_FUROXAN_N2 = 0
!
      NO = 0            ! count number of C's & O's
      NC = 0
      DO 500 I=1,2                ! check the 2-linked N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C') THEN
             NC = NC + 1
             I_C = I2
         ENDIF
         IF (LABEL(I2)(1:1) .EQ. 'O') THEN
            NO = NO + 1
            I_O = I2
         ENDIF
500   CONTINUE
      IF (NO .NE. 1 .OR. NC .NE. 1) RETURN
      IF (NCON(I_O) .NE. 2) RETURN  ! connected O should be 2-linked
      IF (NCON(I_C) .NE. 3) RETURN  ! connected C should be 3-linked
      IQ_FUROXAN_N2 = 20      ! param = 20
      RETURN
!
!----Entry point for the 3-linked ring N in a furoxan
      ENTRY IQ_FUROXAN_N3(I1)   ! I1 is a triple linked N
!
      IQ_FUROXAN_N3 = 0
!
      NO = 0            ! count number of C's & O's
      NC = 0
      NCON_O = 0
      DO 600 I=1,3                ! check the 3-linked N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C') THEN
             NC = NC + 1
             I_C = I2
         ENDIF
         IF (LABEL(I2)(1:1) .EQ. 'O') THEN
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
      END FUNCTION IQ_FUROXAN_O1
!
!----------------------------------------------------------------
!----Function to determine if a 2-linked nitrogen is           12/31/05
!     N3 of a 1,2,3,4-tetrazole.                               12/31/05
!
      FUNCTION IQ_1234tetrazole_n3(I1) ! I1 is a 2-linked N  ! 12/31/05
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_1234tetrazole_n3, n, nc, nn, c3id, c3id_2 ! 12/31/05
      integer :: n3id, nnew, i1, i2, i3, i, j, n2id(2)
!
      IQ_1234tetrazole_n3 = 0     ! no
!
!----Must be linked to two 2-linked N's
!
!     (N3) 72 --> N---N <-- 73 (N4)
!                 "   "
!     (N2) 71 --> N   C--
!                  \ /
!      (N1) 70 -->  N
!                   |
!
      nn = 0        ! count # of 2-linked N's
      DO 100 I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
             NN = NN + 1
             n2id(nn) = i2        ! id of 2-linked N
         endif
100   CONTINUE
      IF (NN .ne. 2) RETURN       ! not a # 72 N
!----Look at the tw0 2-linked N's..is one N2 and one N4?
      n3id = 0
      c3id = 0
      do 200 j=1,2
         i2 = n2id(j)      ! a 2-linked N
             do i=1,2
                i3 = icon(i,i2)
                if (i3 .eq. i1) cycle         ! the original
                if (label(i3)(1:1) .eq. 'N' .and. ncon(i3) .eq. 3) n3id = i3
                if (label(i3)(1:1) .eq. 'C' .and. ncon(i3) .eq. 3) c3id = i3
             enddo
200   continue
!----Are both n3id and c3id non-0?
      if (.not. (n3id .ne. 0 .and. c3id .ne. 0)) return
!----Check to determine that n3id and c3id are bonded
      do i=1,3
         if (icon(i,n3id) .eq. c3id) then
            iq_1234tetrazole_n3 = 72      ! success, param = 72
            return
         endif
      enddo
!
      return                ! no match...not N3
!
      END FUNCTION IQ_1234tetrazole_N3                       ! 1/1/06
!
!----------------------------------------------------------------
!----Function to determine if a 2-linked nitrogen is           1/1/06
!     N2 of a 1,2,3,4-tetrazole.                               1/1/06
!
      FUNCTION IQ_1234tetrazole_n2(I1) ! I1 is a 2-linked N  ! 1/1/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_1234tetrazole_n2, n3, n2                 ! 1/1/06
      integer :: iq_1234tetrazole_n3                         ! 1/1/06
      integer :: i1, i2, i, n2id
!
      IQ_1234tetrazole_n2 = 0     ! no
!
!----Must be linked to two 2-linked N's
!
!     (N3) 72 --> N---N <-- 73 (N4)
!                 "   "
!     (N2) 71 --> N   C--
!                  \ /
!      (N1) 70 -->  N
!                   |
!
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
!----Look at the 2-linked N..is it an N3 with param = 72
      if (iq_1234tetrazole_n3(n2id) .eq. 72) then
         iq_1234tetrazole_n2 = 71          !yes, param = 71
      endif
!
      return
!
      END FUNCTION IQ_1234tetrazole_N2                       ! 1/1/06
!
!----------------------------------------------------------------
!----Function to determine if a 3-linked nitrogen is           1/3/06
!     N1 of a 1,2,3,4-tetrazole.                               1/1/06
!
      FUNCTION IQ_1234tetrazole_n1(I1) ! I1 is a 3-linked N  ! 1/1/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_1234tetrazole_n2, c3, n2                 ! 1/1/06
      integer :: iq_1234tetrazole_n1                         ! 1/1/06
      integer :: i1, i2, i, n2id(2), ntotal                  ! 1/3/06
!
      IQ_1234tetrazole_n1 = 0     ! no
!
!----Must be linked to two 2-linked N's
!
!     (N3) 72 --> N---N <-- 73 (N4)
!                 "   "
!     (N2) 71 --> N   C--
!                  \ /
!      (N1) 70 -->  N
!                   |
!                   X <-- can be C, H, N                     1/3/06
!
      n2 = 0        ! count # of 2-linked N's
      c3 = 0        ! count # of 3-linked C's
      ntotal = 0
      do i=1,2
         n2id = 0      ! id of 2-linked N
      enddo
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,I1)
         if (LABEL(I2)(1:1) .EQ. 'C' .or.  &
             LABEL(I2)(1:1) .EQ. 'N' .or.  &
             LABEL(I2)(1:1) .EQ. 'H') ntotal = ntotal + 1
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) c3 = c3 + 1
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
            N2 = N2 + 1
            n2id(n2) = i2                   ! id of 2-linked N
         endif
100   CONTINUE
      IF (.not. ((N2 .eq. 1 .or. n2 .eq. 2) .and.  &
                 (c3 .eq. 1 .or. c3 .eq. 2) .and.  &
                 (ntotal .eq. 3))) RETURN     ! not a # 70 N
!----Look at the 2-linked N..is it an N2 with param = 71
      do i=1,n2
      if (iq_1234tetrazole_n2(n2id(i)) .eq. 71) then
         iq_1234tetrazole_n1 = 70          !yes, param = 70
         return
      endif
      enddo
!
      return
!
      END FUNCTION IQ_1234tetrazole_N1                       ! 1/1/06
!
!----------------------------------------------------------------
!----Function to determine if a 2-linked nitrogen is           1/1/06
!     N4 of a 1,2,3,4-tetrazole.                               1/1/06
!
      FUNCTION IQ_1234tetrazole_n4(I1) ! I1 is a 2-linked N  ! 1/1/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_1234tetrazole_n4, c3, n2                 ! 1/1/06
      integer :: iq_1234tetrazole_n3                         ! 1/1/06
      integer :: i1, i2, i, n2id
!
      IQ_1234tetrazole_n4 = 0     ! no
!
!----Must be linked to two 2-linked N's
!
!     (N3) 72 --> N---N <-- 73 (N4)
!                 "   "
!     (N2) 71 --> N   C--
!                  \ /
!      (N1) 70 -->  N
!                   |
!
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
!----Look at the 2-linked N..is it an N3 with param = 72
      if (iq_1234tetrazole_n3(n2id) .eq. 72) then
         iq_1234tetrazole_n4 = 73          !yes, param = 73
      endif
!
      return
!
      END FUNCTION IQ_1234tetrazole_N4                       ! 1/1/06
!
!----------------------------------------------------------------
!----Function to determine if a 2-linked nitrogen is           12/30/05
!     part of a 1,2,3-triazole.                                12/30/05
!
      FUNCTION IQ_123triazole_2(I1) ! I1 is a 2-linked N     ! 12/30/05
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_123triazole_2, n, nc, nn, c3id_1, c3id_2
      integer :: n3id, nnew, i1, i2, i, n2id(2)
!
      IQ_123triazole_2 = 0     ! no
!
!----Must be linked to a N, which is linked to 3-linked  C and N
!
!                 C---C
!                 "   "
!          27 --> N   N <-- 27
!                  \ /
!           26 -->  N
!                   |
!
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
!----Look at 3-linked N and count # of 2-linked N's
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
!----Look at the other 2-linked N (nnew), must be bonded to a 3-linked C and a N
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
!----Finally check the c3id_1 and c3id_2 are bonded
      do 150 i=1,3
         i2 = icon(i,c3id_1)
         if (i2 .ne. c3id_2) go to 150    ! compare with id of second C
            iq_123triazole_2 = 27         ! param = 27
            return
150   continue
!
      return      ! return means a failure for iq_123triazole_2
!
      END FUNCTION IQ_123triazole_2                     ! 12/30/05
!
!----------------------------------------------------------------
!
!    Function to determine if 3-linked N is # 26 in triazole ring     12/30/05
!
      function iq_123triazole_3(i1)   ! i1 is a 3-linked N            12/30/05
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_123triazole_2, iq_123triazole_3, nn
      integer :: i1, i2, i, n2id(2)
!
      iq_123triazole_3 = 0
!
      nn = 0
      do 20 i=1,3
         i2 = icon(i,i1)      ! connected to 3-linked N
         if (label(i2)(1:1) .eq. 'N' .and. ncon(i2) .eq. 2) then
             nn = nn + 1
             n2id(nn) = i2
         endif
20    continue
      if (nn .ne. 2) return     ! not connected to two 2-linked N's
!----check that they are param = 27's
      do 40 i=1,2
         if (iq_123triazole_2(n2id(i)) .ne. 27) return
40    continue
!----Dropping thru loop means both are 27's
      iq_123triazole_3 = 26             ! param = 26
!
      end function iq_123triazole_3                 ! 12/30/05
!----------------------------------------------------------------
!----Function to determine if a H bonded to a 3-linked           1/3/06
!     N that is part of a pyrrole ring.                          1/3/06
!
      FUNCTION IQ_pyrrole_nh(I1) ! I1 is a 3-linked N          ! 1/3/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_pyrrole_nh, n, nc, nn                      ! 1/3/06
      integer :: i1, i2, i3, i4, i, j, atom2id(2)              ! 1/3/06
!
      IQ_pyrrole_nh = 0     ! no
!
!----Must be linked to a N, which is linked to 3-linked  C and N
!
!                 X---X <-- C or N
!                 "   "
!      C or N --> X   X <-- C or N
!                  \ /
!                   N
!                   |
!                   H <-- 74

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
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. (ncon(i2) .eq. 3 .or.  &
             ncon(i2) .eq. 2)) then
             NN = NN + 1
             n = n + 1
             atom2id(n) = i2      ! id of atom linked to N1
         endif
100   CONTINUE
      IF (.not.(NC .le. 2 .and. NN .le. 2 .and.  &
          n .eq. 2)) RETURN           ! not proper neighbors
!----Look at two atoms in atom2id...do they have connected atoms?
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
!
      return      ! return means a failure for iq_pyrrole_nh
!
      END FUNCTION IQ_pyrrole_nh                          ! 1/3/06
!
!----------------------------------------------------------------
!----Function to determine if a 2-linked N is part              1/25/06
!     of an amino-imine                                         1/25/06
!
      FUNCTION IQ_aminoimine_2(I1) ! I1 is a 2-linked N       ! 1/23/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_aminoimine_2, nc, nn, nh                  ! 1/27/06
      integer :: i1, i2, i, j, atomid_n                       ! 1/26/06
!
      IQ_aminoimine_2 = 0     ! no
!
!----The functionality...
!          H or C                                             ! 1/27/06
!                \
!                 N---N===C (sp2)
!                /^   ^
!               C |   |
!                75  76
!----2-linked N
      nc = 0
      nn = 0
         do i=1,2
            i2 = icon(i,i1)
            IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) nc = nc + 1
            IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 3) then  ! 1/25/06
               nn = nn + 1                                           ! 1/25/06
               atomid_n = i2                                         ! 1/25/06
            endif
         enddo
         if (nc .ne. 1) return                                       ! 1/25/06
         if (nn .ne. 1) return                                       ! 1/25/06
         nc = 0                                                      ! 1/25/06
         nh = 0                                                      ! 1/27/06
         do i=1,3
            i2 = icon(i,atomid_n)                                    ! 1/25/06
            if (label(i2)(1:1) .eq. 'C') nc = nc + 1                 ! 1/25/06
            if (label(i2)(1:1) .eq. 'H') nh = nh + 1                 ! 1/27/06
         enddo
         if ((nc + nh) .ne. 2) return                                ! 1/27/06
         iq_aminoimine_2 = 76   ! param = 76, 2-linked N
         return    ! return with either a 0 or 76
      end function iq_aminoimine_2                                 ! 1/23/06
!
!--------------------------------------------------------------------------
!
!----Function to determine if a 3-linked N is part              1/23/06
!     of an amino-imine                                         1/20/06
!
      FUNCTION IQ_aminoimine_3(I1) ! I1 is a 3-linked N       ! 1/23/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_aminoimine_2, iq_aminoimine_3, nn         ! 1/23/06
      integer :: i1, i2, i, j, atom2id(3)                     ! 1/23/06
!
!----The functionality...
!                \
!                 N---N===C (sp2)
!                /^   ^
!                 |   |
!                 75  76
!
       iq_aminoimine_3 = 0                                         ! 1/23/06
       nn = 0
         DO I=1,3                           ! examine the N
            I2 = ICON(I,I1)
            IF (LABEL(I2)(1:1) .ne. 'N' .or. ncon(i2) .ne. 2) cycle
               nn = nn + 1
               atom2id(nn) = i2
         enddo
         if (nn .ne. 1) return   ! not linked to one possible imino N
         if (iq_aminoimine_2(atom2id(1)) .ne. 76) return  ! NO
         iq_aminoimine_3 = 75     ! param = 75, 3-linked N       1/23/06
         return
!
      END FUNCTION IQ_aminoimine_3                                ! 1/23/06
!
!----------------------------------------------------------------
!----Function to determine if a 2-bond O is in nitrate ester
!
      FUNCTION IQ_OXYGEN (IT)    ! IT is the 2-linked O
!
!                   O (#52)
!                  /
!-----      C--O--N (#51)
!           (#50)  \
!                   O
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_OXYGEN, IT
!
      INTEGER :: I, I1
      INTEGER :: NC, NN3, N_NAME, NO1, NO2
!IND(2), I2, I3, XL(2), VEC(3,2) are removed because never used
!NAME is changed to LABEL

      IQ_OXYGEN = 50     ! nitrate -O-
!
      NC = 0
      NN3 = 0      ! # of 3 bond N's
      DO 100 I=1,2
         I1 = ICON(I,IT)
         IF (LABEL(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(I1)(1:1) .EQ. 'N' .AND. NCON(I1) .EQ. 3) THEN
            NN3 = NN3 + 1
            N_NAME = I1
         ENDIF
100   CONTINUE
      IF (NC .EQ. 1 .AND. NN3 .EQ. 1) THEN
         NO2 = 0
         NO1 = 0
         DO 110 I=1,3       ! look at 3-linked N
            I1 = ICON(I,N_NAME)
            IF (LABEL(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 2) NO2 = NO2 + 1
            IF (LABEL(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 1) NO1 = NO1 + 1
110      CONTINUE
            IF (NO2 .EQ. 1 .AND. NO1 .EQ. 2) GO TO 115
      ELSE
         IQ_OXYGEN = 0
115   ENDIF
      RETURN
      END FUNCTION IQ_OXYGEN
!
!----------------------------------------------------------------
!----Function to determine if a 3-linked N is a nitrate ester N,
!      a nitramine nitro N or a nitrate (NO3-) N.            7/20/06
!
      FUNCTION IQ_NITRATE_N (IT)    ! IT is the 3-linked N
!
!                   O (#52)         O (#55)   (#84)    7/20/06
!                  /               /           	|      7/20/06
!-----      C--O--N (#51)       N--N (#54)    N(O3)-   7/20/06
!           (#50)  \                \         |        7/20/06
!                   O                O      (#83)      7/20/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_NITRATE_N, IT
!
!NAME is changed to LABEL
!
      INTEGER :: I, I1, NO1, NO2, NN3
!IND(2), I2, I3, XL(2), VEC(3,2) removed, not in use
!
      IQ_NITRATE_N = 51     ! nitrate ester N
!
      NO1 = 0      ! # of 1 bond O's
      NO2 = 0      ! # of 2 bond O's
      NN3 = 0      ! # of 3 bond N's                           ! 2/24/03
      DO 100 I=1,3
         I1 = ICON(I,IT)
         IF (LABEL(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 1) NO1 = NO1 + 1
         IF (LABEL(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 2) NO2 = NO2 + 1
         IF (LABEL(I1)(1:1) .EQ. 'N' .AND. NCON(I1) .EQ. 3) NN3 = NN3 + 1   ! 2/24/03
100   CONTINUE
      if (no1 .eq. 3) then                                     ! 7/20/06
         iq_nitrate_n = 83     ! param = 83 for N of NO3-      ! 7/20/06
         return                                                ! 7/20/06
      endif                                                    ! 7/20/06
      IF (NO1 .EQ. 2 .AND. NO2 .EQ. 1) RETURN
      IF (NO1 .EQ. 2 .AND. NN3 .EQ. 1) THEN                         ! 2/24/03
         IQ_NITRATE_N = 54             ! param = 54 for n of N-nO2    2/24/03
      ELSE                                                          ! 2/24/03
         IQ_NITRATE_N = 0
      ENDIF                                                         ! 2/24/03
      RETURN
      END FUNCTION IQ_NITRATE_N
!
!--------------------------------------------------------------
!
      SUBROUTINE CONNECT2
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL, NATOMS, OR_XYZ
      IMPLICIT NONE
!
!ID is changed to LABEL
!
      INTEGER :: I, J, K
      REAL    :: D, DMAX
!
      DO 100 I=1,NATOMS
         NCON(I) = 0
         DO 90 J=1,NATOMS
            IF (I .EQ. J) GO TO 90
            IF (LABEL(I)(1:1) .EQ. 'H' .AND. LABEL(J)(1:1) .EQ. 'H')   &
     &                         GO TO 90
               DMAX = 1.62
               D = 0.0
               DO K=1,3
                  D = D + (OR_XYZ(K,I) - OR_XYZ(K,J))**2
               ENDDO
!----Check C-Br, C-S and N-S and C-I distances
               IF ((LABEL(I)(1:1) .EQ. 'C' .AND.       & ! C-BR
     &              LABEL(J)(1:2) .EQ. 'BR') .OR.      &
     &             (LABEL(I)(1:2) .EQ. 'BR' .AND.      &
     &              LABEL(J)(1:1) .EQ. 'C')) THEN
                    DMAX = 2.2
                    GO TO 92
               ENDIF
!
91             IF ((LABEL(I)(1:1) .EQ. 'C' .AND.       & !C-S
     &              LABEL(J)(1:1) .EQ. 'S') .OR.       &
     &             (LABEL(I)(1:1) .EQ. 'S' .AND.       &
     &              LABEL(J)(1:1) .EQ. 'C')) THEN
                    DMAX = 1.85
                    GO TO 92
               ENDIF
!
            IF ((LABEL(I)(1:1) .EQ. 'N' .AND.      & ! S-N
     &              LABEL(J)(1:1) .EQ. 'S') .OR.     &
     &             (LABEL(I)(1:1) .EQ. 'S' .AND.     &
     &              LABEL(J)(1:1) .EQ. 'N')) THEN
                    DMAX = 1.69
                    GO TO 92
               ENDIF
!
93            IF ((LABEL(I)(1:1) .EQ. 'C' .AND.    & ! C-I   11/6/03
     &              LABEL(J)(1:1) .EQ. 'I') .OR.   & !       11/6/03
     &             (LABEL(I)(1:1) .EQ. 'C' .AND.   & !       11/6/03
     &              LABEL(J)(1:1) .EQ. 'I')) THEN    !       11/6/03
                    DMAX = 2.2                       !       11/6/03
               ENDIF                                 !       11/6/03
!
92             IF (SQRT(D) .GE. DMAX) GO TO 90
                   NCON(I) = NCON(I) + 1
                   ICON(NCON(I),I) = J
90       CONTINUE
100   CONTINUE
      RETURN
      END SUBROUTINE CONNECT2
!
!------------------------------------------------------------
!
!----Function to determine if a C atom is a cubane C
!
      FUNCTION IQ_CUBANE(I1)
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL, OR_XYZ
      IMPLICIT NONE
      INTEGER :: IQ_CUBANE, I1
!
! NAME is changed to LABEL
!
      INTEGER :: IND(2), I, I2, I3, J, K, L, NC
      REAL    :: XL(2), VEC(3,2), D, COSINE
!      EQUIVALENCE (IND(1),I2), (IND(2),I3)
!
!
      IQ_CUBANE = 0   ! 0 indicates its not cubane-related
!
!----How many C's bonded to I1...must be at least 3
      NC = 0
      DO 20 I=1,4
          I2 = ICON(I,I1)
          IF (LABEL(I2)(1:1) .EQ. 'C') NC = NC + 1
20    CONTINUE
      IF (NC .LT. 3) RETURN   ! not at least 3 C's...return
!----Look for 3 C-C-C bond angles .le. 93 degs
      NC = 0
      DO 200 I=1,3
         DO 200 J=I+1,4
            I2 = ICON(I,I1)
            I3 = ICON(J,I1)
            IF (.NOT. (LABEL(I2)(1:1) .EQ. 'C' .AND.       &  ! must both
     &                 LABEL(I3)(1:1) .EQ. 'C')) GO TO 200  ! be C's
!
            IND(1) = I2     ! these lines added
            IND(2) = I3     ! instead of equivalence
!
            DO 40 L=1,2           ! calc the I2 - I1 & I3 - I1
               XL(L) = 0.0        ! vectors & their lengths
               DO 38 K=1,3
                  VEC(K,L) = OR_XYZ(K,IND(L)) - OR_XYZ(K,I1)
                  XL(L) = XL(L) + VEC(K,L)**2
38             CONTINUE
               XL(L) = SQRT(XL(L))
40          CONTINUE
            D = 0.0
            DO 50 K=1,3
               D = D + VEC(K,1)*VEC(K,2)
50          CONTINUE
            COSINE = D/(XL(1)*XL(2))
!            IF (COSINE .GE. -0.052336) NC = NC + 1  ! cos(93) = -0.052336
            IF (COSINE .GE. -0.0540788) NC = NC + 1  ! cos(93) = -0.0540788 ! 8/19/03
200   CONTINUE
!      IF (NC .NE. 3) RETURN    ! not 3 angles .le. 93 deg
!         IQ_CUBANE = 1        ! 3 angles .le. 93 deg, assume its
      IF (NC .NE. 3) RETURN   ! all 3 angles not .le. 93.1 deg              ! 8/19/03
         IQ_CUBANE = 1        ! all 3 angles .le. 93.1 deg, assume it's     ! 8/19/03
      RETURN                  ! a cubane C atom
      END FUNCTION IQ_CUBANE
!
!----------------------------------------------------------------
!----Function to determine if an O is a nitro O, if the nitro        2/24/03
!     is cubane-linked, a nitrate esters or nitramine                2/24/03
!
      FUNCTION IQ_NITRO_OX(IT)    ! IT is the N linked to the O
!
!                   O (#52)                     O (#55)              2/24/03
!                  /                           /                     2/24/03
!-----      C--O--N (#51)                  N--N (#54)                2/24/03
!           (#50)  \                           \                     2/24/03
!                   O                           O                    2/24/03
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_NITRO_OX, IT
!
! NAME is changed to LABEL
!
      INTEGER :: I, I1, IC, IQ_CUBANE
      INTEGER :: NC, NO, NO1, NO2, NN3
!      EQUIVALENCE (IND(1),I2), (IND(2),I3)
!I2, I3, IND(2) removed, not in use
!  XL(2), VEC(3,2) removed, not in use
!
!
      IQ_NITRO_OX = 6     ! normal nitro group O
!
      NO = 0
      NC = 0
      NO1 = 0      ! # of 1 bond O's
      NO2 = 0      ! # of 2 bond O's
      NN3 = 0      ! # of 3 bond N's
      DO 100 I=1,3
         I1 = ICON(I,IT)
         IF (LABEL(I1)(1:1) .EQ. 'O') NO = NO + 1
         IF (LABEL(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 1) NO1 = NO1 + 1
         IF (LABEL(I1)(1:1) .EQ. 'O' .AND. NCON(I1) .EQ. 2) NO2 = NO2 + 1
         IF (LABEL(I1)(1:1) .EQ. 'N' .AND. NCON(I1) .EQ. 3) NN3 = NN3 + 1
         IF (LABEL(I1)(1:1) .EQ. 'C') THEN
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
!
!----------------------------------------------------------------
!----Function to determine if we have an alcohol O-H.                3/25/03
!
      FUNCTION IQ_ALCOHOL(ITT)    ! ITT is a 2-linked O or H         3/25/03
!                                   or H bonded to 2-linked O        3/25/03
!
!              C[sp3]---O---H                                        3/25/03
!                     #58  59                                        3/25/03
!
      USE mod_preppot, ONLY :  ICON, NCON
      USE mod_preppot, ONLY :  LABEL
      IMPLICIT NONE
      INTEGER :: IQ_ALCOHOL, ITT
!
! NAME changed to LABEL
!
      INTEGER :: INCREMENT, I, IT, IT1, NH, NC4
!
      IQ_ALCOHOL = 0                                               ! 3/25/03
!----Atom could be O or H; look at the O                           ! 3/25/03
      IF (LABEL(ITT)(1:1) .EQ. 'O') THEN                            ! 3/25/03
         INCREMENT = 0                                             ! 3/25/03
         IT = ITT                                                  ! 3/25/03
      ELSE                                                         ! 3/25/03
         IT = ICON(1,ITT)     ! pick up attached O                   3/25/03
         IF (LABEL(IT)(1:1) .EQ. 'O') INCREMENT = 1                 ! 3/25/03
      ENDIF                                                        ! 3/25/03
      NH = 0                                                       ! 3/25/03
      NC4 = 0                                                      ! 3/25/03
      DO I=1,2                                                     ! 3/25/03
         IT1 = ICON(I,IT)                                          ! 3/25/03
         IF (LABEL(IT1)(1:1) .EQ. 'H') NH = NH + 1                  ! 3/25/03
         IF (LABEL(IT1)(1:1) .EQ. 'C' .AND. NCON(IT1) .EQ. 4)  &    ! 3/25/03
     &          NC4 = NC4 + 1                                      ! 3/25/03
      ENDDO                                                        ! 3/25/03
      IF (NH .EQ. 1 .AND. NC4 .EQ. 1) THEN                         ! 3/25/03
         IQ_ALCOHOL = 58 + INCREMENT   ! O of alcohol, param = 58    3/25/03
      ENDIF                            ! or alcohol H, param = 59    3/25/03
      RETURN                                                       ! 3/25/03
      END FUNCTION IQ_ALCOHOL                                      ! 3/25/03
!
!----------------------------------------------------------------
!----Function to identify the 3 N's in an azide                     ! 3/5/03
!
      FUNCTION IQ_AZIDE_N (IT)    ! IT is a 1 or 2-linked N         ! 3/5/03
!
!         C---N===N===N                                             ! 3/5/03
!            32  56  57                                             ! 3/5/03
!
      USE mod_preppot, ONLY : ICON, NCON, CONNECT_FLAG
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_AZIDE_N, IT
!
! NAME changed to LABEL
!
      INTEGER :: I, IT1, IT2, IT3, J, MN1, MN2, NC, NN1, NN2
!
      NC = 0    ! number of attached C's                            ! 3/5/03
      NN1 = 0   ! number of 1-linked N's                            ! 3/5/03
      NN2 = 0   ! number of 2-linked N's                            ! 3/5/03
!
      DO I=1,NCON(IT)                                               ! 3/5/03
         IT1 = ICON(I,IT)                                           ! 3/5/03
         IF (LABEL(IT1)(1:1) .EQ. 'C') NC = NC + 1                   ! 3/5/03
         IF (LABEL(IT1)(1:1) .EQ. 'N' .AND. NCON(IT1) .EQ. 1)   &    ! 3/5/03
     &        NN1 = NN1 + 1                                         ! 3/5/03
         IF (LABEL(IT1)(1:1) .EQ. 'N' .AND. NCON(IT1) .EQ. 2)   &    ! 3/5/03
     &        NN2 = NN2 + 1                                         ! 3/5/03
      ENDDO                                                         ! 3/5/03
!----First look for the terminal N atom                             ! 3/5/03
      IF (NCON(IT) .EQ. 1 .AND. NN2 .EQ. 1) THEN                    ! 3/5/03
         IT2 = ICON(1,IT)   ! check out the middle N to be certain  ! 3/5/03
         MN1 = 0   ! number of 1-linked N's                         ! 3/5/03
         MN2 = 0   ! number of 2-linked N's                         ! 3/5/03
         DO I=1,2                                                   ! 3/5/03
            IT3 = ICON(I,IT2)                                       ! 3/5/03
            IF (LABEL(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 1)  &  ! 3/5/03
     &           MN1 = MN1 + 1                                      ! 3/5/03
            IF (LABEL(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 2)  &  ! 3/5/03
     &           MN2 = MN2 + 1                                      ! 3/5/03
         ENDDO                                                      ! 3/5/03
         IF (MN1 .EQ. 1 .AND. MN2 .EQ. 1) THEN                      ! 3/5/03
            IQ_AZIDE_N = 57     ! terminal N, param = 57            ! 3/5/03
            RETURN                                                  ! 3/5/03
         ENDIF                                                      ! 3/5/03
      ENDIF                                                         ! 3/5/03
!----Second look for the central N atom                             ! 3/5/03
      IF (NCON(IT) .EQ. 2 .AND. NN1 .EQ. 1 .AND. NN2 .EQ. 1) THEN   ! 3/5/03
          IQ_AZIDE_N = 56     ! central N, param = 56               ! 3/5/03
          RETURN                                                    ! 3/5/03
      ENDIF                                                         ! 3/5/03
!----Third look for the C-linked N atom                             ! 3/5/03
      IF (NCON(IT) .EQ. 2 .AND. NC .EQ. 1 .AND. NN2 .EQ. 1) THEN    ! 3/5/03
         DO 300 I=1,2                                               ! 3/5/03
            IT2 = ICON(I,IT)                                        ! 3/5/03
            IF (LABEL(IT2)(1:1) .EQ. 'C') GO TO 300                  ! 3/5/03
            NN1 = 0    ! IT2 is now the central N                   ! 3/5/03
            NN2 = 0                                                 ! 3/5/03
            DO J=1,2                                                ! 3/5/03
               IT3 = ICON(J,IT2)                                    ! 3/5/03
               IF (LABEL(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 1) & ! 3/5/03
     &              NN1 = NN1 + 1                                   ! 3/5/03
               IF (LABEL(IT3)(1:1) .EQ. 'N' .AND. NCON(IT3) .EQ. 2) & ! 3/5/03
     &              NN2 = NN2 + 1                                   ! 3/5/03
            ENDDO                                                   ! 3/5/03
            IF (NN2 .EQ. 1 .AND. NN1 .EQ. 1) THEN                   ! 3/5/03
               IQ_AZIDE_N = 32    ! C-linked N, param = 32          ! 3/5/03
               RETURN                                               ! 3/5/03
            ENDIF                                                   ! 3/5/03
300      CONTINUE                                                   ! 3/5/03
      ENDIF                                                         ! 3/5/03
      IQ_AZIDE_N = 0                                                ! 3/5/03
      END FUNCTION IQ_AZIDE_N                                       ! 3/5/03
!
!---------------------------------------------------------------------
!----Create name and open chem3d file
!
      SUBROUTINE OPEN_CHEM3D
      USE mod_preppot, ONLY :  NCC, OUTFIL, CHEM3D_NAME
      IMPLICIT NONE
!DATNM2 is removed from use
!      CHARACTER (LEN=1)  :: OUTNAM(19), DATNAM(11)
!      INTEGER :: J
!
!      EQUIVALENCE (OUTFIL, OUTNAM), (DATNM2, DATNAM)
!
!      DATNM2 = '.cc1       '     ! create name for CHEM3D (cc1) file
!      DO 29 J=1,11               ! as XXXXXX.cc1
!29       OUTNAM(J+NCC) = DATNAM(J)
!
      OUTFIL = OUTFIL(1:NCC-1) // '.cc1'
      CHEM3D_NAME = OUTFIL           ! save name of CHEM3D file
      WRITE (20,'(A19)') CHEM3D_NAME   ! on unit # 20
      OPEN (UNIT=14,FILE=OUTFIL,STATUS='UNKNOWN',     & ! CHEM3D on # 14
     &         FORM='FORMATTED')
      RETURN
      END SUBROUTINE OPEN_CHEM3D
!
!----------------------------------------------------------------
!----Function to determine if a N is one of the 4 N's is the
!     Z-tetraazapentalene moiety
!
      FUNCTION IQ_N_AZAPENTALENE(IT) ! IT is a N number
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      INTEGER :: IQ_N_AZAPENTALENE, IT
!
! NAME is changed to LABEL
!
      INTEGER ::  I, I1, I2, IN, INN(2), NN, NC
!      EQUIVALENCE (IND(1),I2), (IND(2),I3)   ! IND(2) is not in use
!
!XL(2), VEC(3,2), I3 removed, not in use
!
!
      IQ_N_AZAPENTALENE = 0   ! 0 means can't identify type
!
      IF (NCON(IT) .LT. 2 .OR. NCON(IT) .GT. 3) RETURN
      GO TO (10,20), (NCON(IT) - 1)
!----2 attachments...check for an end N
10    NN = 0
      NC = 0
      DO 100 I=1,2
         I1 = ICON(I,IT)
         IF (LABEL(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(I1)(1:1) .EQ. 'N') THEN
            NN = NN + 1
            IN = I1  ! IN = N attached to the N in question
         ENDIF
100   CONTINUE
      IF (NC .NE. 1 .OR. NN .NE. 1) RETURN
      IF (NCON(IN) .NE. 3) RETURN
!----Have 2 attachments, 1 C & 1 N and N has 3 attachments
      NN = 0
      NC = 0
      DO 110 I=1,3
         I1 = ICON(I,IN)
         IF (LABEL(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(I1)(1:1) .EQ. 'N') NN = NN + 1
110   CONTINUE
      IF (NC .EQ. 1 .AND. NN .EQ. 2) THEN
         IQ_N_AZAPENTALENE = 19   ! param = 19 for Z-tetraazapentalene N
      ENDIF
      RETURN
!----3 attachments...check for a middle N
20    NN = 0
      NC = 0
      DO 210 I=1,3
         I1 = ICON(I,IT)
         IF (LABEL(I1)(1:1) .EQ. 'C') NC = NC + 1
         IF (LABEL(I1)(1:1) .EQ. 'N') THEN
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
      END FUNCTION IQ_N_AZAPENTALENE
!
!------------------------------------------------------------
      SUBROUTINE ADJUST      ! version 2               1/25/01
!
!----Adjust C-H, N-H, N-F and N-O lengths
!    added sulfur    11-2-2000
!
      USE mod_preppot, ONLY : CELL, LABEL, XYZ,        &
     &                OR_XYZ, R, NATOMS,               &
     &                N_OTHER, N_OX, ISPECIAL,         &
     &                NEW, XYZT, N_2, N_CENTRAL, XYZR, &
     &                XMAT2, C1
! O_XYZ is changed with OR_XYZ
      USE mod_preppot, ONLY : T, N_NEAR, N_ATOMS, J1, J2
      IMPLICIT NONE
!
      LOGICAL :: READFILE
      CHARACTER (LEN=1) :: AH, GEN, RF
      CHARACTER (LEN=5) :: ID2
!      CHARACTER (LEN=5) :: ID(500)
!      EQUIVALENCE (LABEL, ID)
!ID is changed to LABEL
! I_ALTERED(500) removed, not in use
      REAL :: RS(3,3), TR(3), DSUM(2)
      REAL :: C2(3), C3(3), C4(3), C(3,500)
      REAL :: VAVG(3), V1(3), V2(3), V3(3), V4(3), VH(3), VX(3)
      REAL :: X12(3), X23(3), X31(3)
      REAL :: XMAT1(3,3), XMAT3(3,3), XYZR4(3,3), XYZR3(3,3)
!
! RINV(3,3) removed, not in use
!
!      EQUIVALENCE (C, OR_XYZ)
! OR_XYZ is set later ..
      REAL ::  XYZ_NEW(3,500), XYZH(3,500)
!
      INTEGER :: IAH, IGEN, IH_TYPE
      INTEGER :: I, J, K, L, L1, L2, L3, M, N, N_3
      INTEGER :: KNH, NB2O, NEW_NATOMS, NUMBER
      REAL    :: ANGLE, BEST_ANG, D, DET, DIR, DISTANCE, DMAX, ROTN
      REAL    :: S, V_DOT, XHL, XHOR, XVER
!
      READFILE = .false.
      NEW = 0
      DMAX = 1.65
!
!----Calc matrix to go from orthogonal to fractional coordinates
!
3510  DO 3511 I=1,3                 ! copy frac to cart transformation matrix
         DO 3511 J=1,3              ! [R] to T in preparation for forming
            T(I,J) = R(I,J)         ! inverse
3511  CONTINUE
      CALL MINV (T, 3, DET, J1, J2)  ! cart to fractional, T replaced by inverse
!
!----Determine connectivity....use max distance of 1.65 Angs
      CALL CONNECTION (DMAX)
!
!----Should the C-H lengths be adjusted (lengths changed) or idealized?
      PRINT 10
10    FORMAT (' Should the C-H positions be idealized (I) or',   &
     &        ' stretched (S) [I]: ',$)
      READ (5,12) AH
12    FORMAT (A1)
      IF (AH .EQ. ' ' .OR. AH .EQ. 'I' .OR. AH .EQ. 'i') THEN
         IAH = 1               ! idealize
         PRINT 16
         WRITE (23,16)
16       FORMAT (' C-hydrogens will be idealized...all input',   &
     &           ' C-H''s have been deleted')
         NEW_NATOMS = N_OTHER + N_OX    ! drop H's from list
!----Keep the H's linked to N
         KNH = 0
         DO 2000 I=NEW_NATOMS,NATOMS
             IF (LABEL(I)(1:1) .NE. 'H') GO TO 2000
             L = N_ATOMS(1,I)   ! # of atom linked to H...is it N
             IF (LABEL(L)(1:1) .NE. 'N') GO TO 2000
             KNH = KNH + 1    ! an N-H, move the H atom
             M = NEW_NATOMS + KNH
             LABEL(M) = LABEL(I)
             DO 2005 J=1,3
                OR_XYZ(J,M) = OR_XYZ(J,I)    ! move I to M
                XYZ(J,M) = XYZ(J,I)
2005         CONTINUE
2000     CONTINUE
         NATOMS = NEW_NATOMS + KNH  ! update atom # for any N-H's kept
         IF (KNH .EQ. 0) GO TO 2010
            PRINT 2015, KNH
2015        FORMAT (' Positions of',I3,' hydrogen linked to N have',   &
     &              ' kept and moved in the atom storage list')
2010     IGEN = 0
         IF (ISPECIAL .EQ. 1) THEN
            PRINT 17
17          FORMAT (' Molecules occupies a special position...'/        &
     &              '    should the C-hydrogens also be generated by',  &
     &              ' symmetry? [Y]: ',$)
            READ (5,12) GEN
            IF (GEN .EQ. 'Y' .OR. GEN .EQ. ' ' .OR. GEN .EQ. 'y')  &
     &              IGEN = 1
            IF (IGEN .EQ. 0) GO TO 19
            PRINT 18
18          FORMAT (' Give the symmetry operation for H generation ', &
     &        'as r11, r12, r13, t1, r21, ... t2, ... r33, t3;'/      &
     &        '    the -x, 1/2-y, -z sym op should be entered ',      &
     &        'as -1, 0, 0, 0, 0, -1, 0, 0.5, 0, 0, -1, 0...'/7X,$)
            READ (5,*) ((RS(I,J),J=1,3), TR(I), I=1,3)
19       ENDIF
      ELSE
         IAH = 0                    ! just stretch (or shrink) the H position
      ENDIF
!
!----Determine connectivity again...just in case number of atoms changed
      CALL CONNECTION (DMAX)
!
!----Should the N-H lengths be altered?
      PRINT 101
101   FORMAT (' Should the N-H lengths be stretched? [Y]: ',$)
      READ (5,12) AH
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') THEN
      DO 2001 I=1,NATOMS            ! just stretch the X-H position
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 2001
         J = N_ATOMS(1,I)           ! J = atom to which H is linked
         IF (LABEL(J)(1:1) .NE. 'N') GO TO 2001     ! C-H bond
!----N-H extension to 1.013 or N(+)-H to 1.031 Angs
         IF (N_NEAR(J) .EQ. 4) CALL STRETCH (J, I, 1.031)
         IF (N_NEAR(J) .LE. 3) CALL STRETCH (J, I, 1.013)
2001  CONTINUE
      ENDIF
!                                                                3/25/03
!----Should the O-H lengths be altered?                          3/25/03
      PRINT 501                                                ! 3/25/03
501   FORMAT (' Should the O-H lengths be stretched? [Y]: ',$) ! 3/25/03
      READ (5,12) AH                                           ! 3/25/03
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') THEN  ! 3/25/03
      DO 2011 I=1,NATOMS        ! stretch the O-H position       3/25/03
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 2011                ! 3/25/03
         J = N_ATOMS(1,I)       ! J = atom to which H is linked  3/25/03
         IF (LABEL(J)(1:1) .NE. 'O') GO TO 2011     ! O-H bond   3/25/03
!----O-H extension to 0.97 Angs                                  3/25/03
         CALL STRETCH (J, I, 0.97)                             ! 3/25/03
2011  CONTINUE                                                 ! 3/25/03
      ENDIF                                                    ! 3/25/03
!
!----Should the N-F lengths be altered?
      PRINT 311
311   FORMAT (' Should the N-F lengths be stretched? [Y]: ',$)
      READ (5,12) AH
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') THEN
      DO 2071 I=1,NATOMS            ! just stretch the N-F position
         IF (LABEL(I)(1:1) .NE. 'F') GO TO 2071
         J = N_ATOMS(1,I)           ! J = atom to which F is linked
         IF (LABEL(J)(1:1) .NE. 'N') GO TO 2071     ! N-F bond
!----N-F extension to 1.40
         CALL STRETCH (J, I, 1.40)
2071  CONTINUE
      ENDIF
!
!----Should the NO2 distances be altered to 1.22 Angs
      PRINT 15
15    FORMAT (' Should the nitro group N-O distances be set to 1.22',  &
     &        ' Angs? [Y]: ',$)
      READ (5,12) AH
!
!----Report atom connectivities...
      WRITE (6,2220)
      WRITE (23,2220)
2220  FORMAT (' Atom connectivity....'/    &
     &        ' atom ....  linked to ...')
      DO I=1,NATOMS
         IF (LABEL(I)(1:1) .EQ. 'C' .AND.   & ! if a C linked to 2 other
     &       N_NEAR(I) .EQ. 2) THEN         ! atoms, show bond lengths
            DO 2224 L1=1,2
               L3 = N_ATOMS(L1,I)
               D = 0.0
                 DO 2226 L2=1,3
                   D = D + (OR_XYZ(L2,I) - OR_XYZ(L2,L3))**2
2226             CONTINUE
               DSUM(L1) = SQRT(D)
2224        CONTINUE
            WRITE (6,2232) LABEL(I), (LABEL(N_ATOMS(K,I)),K=1,2),  &
     &            (LABEL(N_ATOMS(K,I)),DSUM(K),K=1,2)
            WRITE (23,2232) LABEL(I), (LABEL(N_ATOMS(K,I)),K=1,2), &
     &            (LABEL(N_ATOMS(K,I)),DSUM(K),K=1,2)
2232        FORMAT (1X,A5,'....  ',2(A5,2X),13X,': d to ',A5,'=',F5.2,  &
     &              '; d to ',A5,'=',F5.2)
         ELSE
           WRITE (6,2222) LABEL(I), (LABEL(N_ATOMS(K,I)),K=1,N_NEAR(I))
           WRITE (23,2222) LABEL(I), (LABEL(N_ATOMS(K,I)),K=1,N_NEAR(I))
2222       FORMAT (1X,A5,'....  ',6(A5,2X))
         ENDIF
      ENDDO
!
!----Fix the N-O's ?
      IF (AH .EQ. ' ' .OR. AH .EQ. 'Y' .OR. AH .EQ. 'y') THEN
      DO 100 I=1,NATOMS
         IF (LABEL(I)(1:1) .NE. 'O') GO TO 100
         IF (N_NEAR(I) .NE. 1) GO TO 100     ! not an N-O type O
         J = N_ATOMS(1,I)                    ! ID of 1 atom linked to the O
         IF (LABEL(J)(1:1) .NE. 'N') GO TO 100        ! gotta be a N
!----Check if the N is connected to 2 O's, each of which has only 1 bond
         NB2O = 0
            DO 110 L=1,N_NEAR(J)   ! J is the N connected to the O (I)
               K = N_ATOMS(L,J)    ! K = # of an atom connected to N (J)
               IF (LABEL(K)(1:1) .NE. 'O') GO TO 110
               IF (N_NEAR(K) .EQ. 1) NB2O = NB2O + 1    ! count # of singly
110         CONTINUE                                    ! bonded O's, must
            IF (NB2O .NE. 2) GO TO 100                  ! be 2 of them
!----Set a position along the N-O (I to J) vector
         CALL STRETCH (J, I, 1.22)
100   CONTINUE
      ENDIF
!
!----Check the C-H situation
      IF (IAH .EQ. 1) GO TO 444
      DO 200 I=1,NATOMS            ! just stretch the X-H position
         IF (LABEL(I)(1:1) .NE. 'H') GO TO 200
         J = N_ATOMS(1,I)           ! J = atom to which H is linked
         IF (LABEL(J)(1:1) .NE. 'C') GO TO 200     ! C-H bond ?
!----For C-H extensions, check the C hydridization
3000     GO TO (3001, 3001, 3003, 3004), N_NEAR(J)  ! how many atoms on J ?
3001     PRINT 3002, I, J
3002     FORMAT (' Cannot adjust distance between atoms ',A5,   &
     &           ' and ',A5)
         GO TO 200
3003     D = 1.084                      ! sp2
         GO TO 3010
3004     D = 1.098                      ! sp3
3010     CALL STRETCH (J, I, D)         ! adjust XMATto length D
200   CONTINUE
      RETURN                  ! $$$$$$$$$$$$$$$$$$$$$$$
!
!----Calc ideal positions for sp, sp2 and sp3 H's linked to C
!
!----Atom names for H atom generation....
!     central atom ID,
!     1/2/3 for sp/sp2/sp3,
!     C-H distance
!    If this is a CH3-C(sp2), the generated methyl will have 1 H in the
!     sp2 plane.  Use -3 to rotate the methyl by 90 deg.
!    If a CH3-C(sp3), generated methyl H's will be eclipsed.  Use -4 to
!     rotate methyl by 60 deg.
444    WRITE (*,833)                                                 ! 3/6/03
833    FORMAT ('Read H information dynamically or from file ?'/      & ! 3/6/03
     &           ' file name = h_info'/                              & ! 3/6/03
     &           ' information = central atom name'/                 & ! 3/6/03
     &           '               1/2/3 for sp/sp2/sp3 hybridization'/& ! 3/6/03
     &           '               C-H distance'/                      & ! 3/6/03
     &           ' format = (A5,I2,F8.2) [F]: ',$)                   ! 3/6/03
       READ (5,303) RF
303    FORMAT (A1)
       IF (RF .EQ. 'f' .OR. RF .EQ. 'F' .OR. RF .EQ. ' ') THEN
          READFILE = .true.      ! get H generation info from a file
          OPEN (UNIT = 38, FILE='h_info', STATUS = 'OLD')
       ENDIF
300    ID2 = '     '
       IF (READFILE) THEN
          READ (38,302,END=1000) ID2, IH_TYPE, XHL    ! from file
       ELSE
          PRINT 301                                   ! dynamically
301       FORMAT (' Give central atom name, 1/2/3 for # H''s, C-H dist', &
     &            ' (A5,I2,F8.2): ',$)
          READ (5,302) ID2, IH_TYPE, XHL
302       FORMAT (A5,I2,F8.2)
       ENDIF
       IF (ID2(1:1) .EQ. ' ') GO TO 1000
!----Is an atom id like C_23 or a sequence number being used?
       DO 701 I=1,NATOMS                 ! search for id like
          IF (ID2 .EQ. LABEL(I)) GO TO 715  ! C_23
701    CONTINUE
       READ (ID2,*) NUMBER
       ID2 = LABEL(NUMBER)
       PRINT 705, NUMBER, ID2
705       FORMAT (' Atom sequence number',I3,' encountered...',   &
     &            'assumed to refer to CSD atom ',A5)
!300   READ (10,1234,END=1000) LINE
!1234  FORMAT (A40)
!300   READ (10,302,END=1000) ID2, IH_TYPE, XHL
!      READ (LINE,302) ID2, IH_TYPE, XHL
!302   FORMAT (A5,I2,F8.2)
715    ROTN = 0.0                  ! initialize methyl rotn to 0 degs
       IF (IH_TYPE .EQ. -3) THEN
         ROTN = 90.0               ! typical for methyl-C(sp2)
         IH_TYPE = 3
       ENDIF
       IF (IH_TYPE .EQ. -4) THEN
         ROTN = 60.0               ! typical for methyl-C(sp3)
         IH_TYPE = 3
       ENDIF
!----for sp3, this H-C-C-other dihedral angle
308   GO TO (305, 310, 315, 315), IH_TYPE
305   WRITE (6,306)
306   FORMAT (' Do not presently handle sp')
      GO TO 300
!----Sp2...find sequence numbers of central and surrounding atoms
310   DO 600 I=1,NATOMS
         IF (ID2 .NE. LABEL(I)) GO TO 600
         DO J=1,3
!           C1(J) = C(J,I)                     ! 10/9/03, C replaced
            C1(J) = OR_XYZ(J,I)                ! with OR_XYZ
         ENDDO
         N_CENTRAL = I
      GO TO 610
600   CONTINUE
604   WRITE (6,620) ID2, (LABEL(L),L=1,NATOMS)
620   FORMAT (' Cannot find sp2 atom ',A5,' in list; the atoms are...'/  &
     &        (5X,10(1X,A5)))
      GO TO 300
!----Who are the surrounding atoms?
610   IF (N_NEAR(N_CENTRAL) .EQ. 2) THEN
         DO J=1,3
!           C2(J) = C(J,N_ATOMS(1,N_CENTRAL))            ! 10/9/03, C replaced
!           C3(J) = C(J,N_ATOMS(2,N_CENTRAL))            ! with OR_XYZ
            C2(J) = OR_XYZ(J,N_ATOMS(1,N_CENTRAL))
            C3(J) = OR_XYZ(J,N_ATOMS(2,N_CENTRAL))
         ENDDO
      ELSE
         WRITE (6,630) ID2, N_NEAR(N_CENTRAL)
630      FORMAT (' Atom ',A5,' indicated to be sp2 has ',I2,  &
     &           ' neighbors')
         GO TO 300
      ENDIF
!----On with the sp2 case...C2-C1-C3
      CALL VECTOR (V1, C1, C2)     ! form vector from C2 to C1
      CALL VECTOR (V2, C1, C3)     ! from C3 to C1
      CALL N_VECTOR (V1)                    ! normalize
      CALL N_VECTOR (V2)
!----Average V1 and V2 to get mid-point
      DO I=1,3
         VAVG(I) = ((V1(I) + V2(I))/2.0)
      ENDDO
      CALL N_VECTOR(VAVG)  ! VAVG = unit vector along C-H bond
!----Get new H position
      DO I=1,3
         VH(I) = C1(I) + XHL*VAVG(I)
      ENDDO
!----Convert to fractional coordinates
      CALL MTIMES2 (VH, T, XYZH(1,1))
!----Create H label and and add orthogonal coordinates to main list
      NEW = NEW + 1
      L = NATOMS + NEW
      CALL MAKE_H_LABEL(L)                            ! 10/9/03
      DO I=1,3
         XYZ_NEW(I,NEW) = XYZH(I,1)  ! fractional
         OR_XYZ(I,L) = VH(I)          ! orthogonal into main list
         XYZ(I,L) = XYZH(I,1)        ! fractional into main list
      ENDDO
      WRITE (6,350) LABEL(L), ID2, (XYZH(I,1),I=1,3), VH
350   FORMAT (1X,A5,' attached to sp2 atom ',A5,'...'/  &
     &        '   fract coord:',3F7.4,'; orthog coord:',3F7.3)
      GO TO 300
!
!----Sp3....find seq no. of central atom and number of connections
315   DO 320 I=1,NATOMS
         IF (ID2 .NE. LABEL(I)) GO TO 320
         DO J=1,3
!           C1(J) = C(J,I)             ! 10/9/03, C replaced
            C1(J) = OR_XYZ(J,I)             ! with OR_XYZ
         ENDDO
         N_CENTRAL = I
      GO TO 330
320   CONTINUE
      WRITE (6,321) ID2, (LABEL(L),L=1,NATOMS)
321   FORMAT (' Cannot find sp3 atom ',A5,' in list; the atoms are...'/  &
     &        (5X,10(1X,A5)))
      GO TO 300
!----# of attached atoms tells us what to do
330   GO TO (510, 820, 830, 840), N_NEAR(N_CENTRAL)
840   WRITE (6,841) ID2
841   FORMAT (' Atom ',A5,' has 4 attachments already')
      GO TO 300
!----2 attached atoms....make 2 H's
820   DO J=1,3
!        C2(J) = C(J,N_ATOMS(1,N_CENTRAL))       ! 10/9/03, C replaced
!        C2(J) = C(J,N_ATOMS(1,N_CENTRAL))       ! with OR_XYZ
         C2(J) = OR_XYZ(J,N_ATOMS(1,N_CENTRAL))  ! corrected on
         C3(J) = OR_XYZ(J,N_ATOMS(2,N_CENTRAL))  !   4/22/06
      ENDDO
!----On with the sp3 case...C2-C1-C3
      CALL VECTOR (V1, C2, C1)              ! form vectors
      CALL VECTOR (V2, C3, C1)
      CALL N_VECTOR (V1)                    ! normalize
      CALL N_VECTOR (V2)
!----Average V1 and V2 to get mid-point
      DO I=1,3
         VAVG(I) = (V1(I) + V2(I))/2.0
      ENDDO
      CALL N_VECTOR(VAVG)
!----Vector cross product of V1 and V2
      CALL V_CROSS (VX, V1, V2)     ! VX is perp to V1 and V2
      CALL N_VECTOR (VX)
      XHOR = XHL*COSD(109.5/2.0)
      XVER = XHL*SIND(109.5/2.0)
!----Form 2 H positions
      DIR = -1.0
      DO 800 I=1,2
         DIR = -DIR
         DO 802 J=1,3
            XYZR4(J,I) = C1(J) - XHOR*VAVG(J) + DIR*XVER*VX(J)
802      CONTINUE
      CALL MTIMES2 (XYZR4(1,I), T, XYZH(1,I))   ! convert to frac
800   CONTINUE
      WRITE (6,813) ID2
813   FORMAT (' Two H atoms attached to sp3 atom ',A5,' are at...')
      DO 814 L=1,2
         NEW = NEW + 1
         L2 = NATOMS + NEW
         CALL MAKE_H_LABEL(L2)                        ! 10/9/03
         DO 816 I=1,3
            XYZ_NEW(I,NEW) = XYZH(I,L)   ! fractional
            OR_XYZ(I,L2) = XYZR4(I,L)     ! orthog into main list
            XYZ(I,L2) = XYZH(I,L)        ! fractional into main list
816      CONTINUE
      WRITE (6,909) LABEL(L2), (XYZ(I,L2),I=1,3),   &
     &                         (XYZR4(I,L),I=1,3)
814   CONTINUE
      GO TO 300
!
!----3 atoms attached to 1 sp3 atom...get 1 H
830   DO J=1,3                    ! C1 =
!        C1(J) = C(J,N_CENTRAL)        ! central atom       ! 10/9/03, C replaced
         C1(J) = OR_XYZ(J,N_CENTRAL)   ! central atom       ! with OR_XYZ
      ENDDO
!----Who are the surrounding atoms?
410   DO J=1,3
!        C2(J) = C(J,N_ATOMS(1,N_CENTRAL))     ! 10/9/03, replace C
!        C3(J) = C(J,N_ATOMS(2,N_CENTRAL))     ! with OR_XYZ in
!        C4(J) = C(J,N_ATOMS(3,N_CENTRAL))     ! these 3 lines
         C2(J) = OR_XYZ(J,N_ATOMS(1,N_CENTRAL))
         C3(J) = OR_XYZ(J,N_ATOMS(2,N_CENTRAL))
         C4(J) = OR_XYZ(J,N_ATOMS(3,N_CENTRAL))
      ENDDO
      print 3535, i, n_central, n_near(i), (n_atoms(j,i),j=1,n_near(i))
3535  format ('  i, n_central =',8i5)
!----On with the sp3 case...C1 is bonded to C2, C3 and C4
      CALL VECTOR (V1, C2, C1)              ! form vectors
      CALL VECTOR (V2, C3, C1)
      CALL VECTOR (V3, C4, C1)
      CALL N_VECTOR (V1)                    ! normalize
      CALL N_VECTOR (V2)
      CALL N_VECTOR (V3)
!----Form the three vector cross products
      CALL V_CROSS (X12, V1, V2)
      CALL V_CROSS (X23, V2, V3)
      CALL V_CROSS (X31, V3, V1)
!----Average the 3 vector cross products and normalize
      DO I=1,3
         VAVG(I) = X12(I) + X23(I) + X31(I)
      ENDDO
      CALL N_VECTOR(VAVG)
!----Take dot product of V1 & VAVG to determine direction
      V_DOT = 0.0
      DO I=1,3
         V_DOT = V1(I)*VAVG(I) + V_DOT
      ENDDO
      S = +1.0
      IF (V_DOT .GE. 0.0) S = -1.0
!----Get new H position
      DO I=1,3
         VH(I) = C1(I) + S*XHL*VAVG(I)
      ENDDO
!----Average V1, V2 and V3 to get mid-point
!     DO I=1,3
!        VAVG(I) = (V1(I) + V2(I) + V3(I))/3.0
!     ENDDO
!     CALL N_VECTOR(VAVG)
!----Get new H position
!     DO I=1,3
!        VH(I) = C1(I) - XHL*VAVG(I)
!     ENDDO
!----Convert to fractional coordinates
      CALL MTIMES2 (VH, T, XYZH(1,1))
      WRITE (6,401) ID2
401   FORMAT (' One H atom attached to sp3 atom ',A5,' is at...')
      NEW = NEW + 1
      L2 = NATOMS + NEW
      CALL MAKE_H_LABEL(L2)                        ! 10/9/03
      DO I=1,3
         XYZ_NEW(I,NEW) = XYZH(I,1)   ! fractional
         OR_XYZ(I,L2) = VH(I)          ! orthog into main list
         XYZ(I,L2) = XYZH(I,1)        ! fractional into main list
      ENDDO
      WRITE (6,909) LABEL(L2), (XYZ(I,L2),I=1,3),   &
     &                         (OR_XYZ(I,L2),I=1,3)
      GO TO 300
!----1 attached atom...make a methyl group.  Who is attached?
510   N_2 = N_ATOMS(1,N_CENTRAL)         ! N_2 is the # of the atom bonded
      DO J=1,3                           ! to the methyL C ATOM
!        C2(J) = C(J,N_2)                ! C2 = coords of the bonded atom
         C2(J) = OR_XYZ(J,N_2)    ! 10/9/03, replace C with OR_XYZ in above line
      ENDDO
!----Find just one atom in the C3 position....C1-C2-C3
      DO 512 I=1,N_NEAR(N_2)
         N_3 = N_ATOMS(I,N_2)               ! N_3 is # of C atom in 3rd position
         IF (N_3 .EQ. N_CENTRAL) GO TO 512  ! N_3 must not equal N_CENTRAL
         DO 513 J=1,3
!           C3(J) = C(J,N_3)                ! C3 = coords of the 3rd atom
            C3(J) = OR_XYZ(J,N_3)  ! 10/09/03, replace C with OR_XYZ in above line
513      CONTINUE
         GO TO 519
512   CONTINUE
!----VX is in the plane of V1 and V4...generate an initial H location in the
!     more-or-less -V1 direction
519   CALL VECTOR (V1, C1, C2)           ! V1 = C1 - C2
      CALL N_VECTOR (V1)
      CALL VECTOR (VX, C3, C2)           ! VX = C3 - C2
      CALL N_VECTOR (VX)
      CALL V_CROSS (V3, V1, VX)
      CALL N_VECTOR (V3)                 ! V3 is normal to V1
      CALL V_CROSS (V4, V1, V3)
      CALL N_VECTOR (V4)                 ! V4 is normal to V1 and V3
!----Form matrix so that V1 (C2 - C1) can be rotation axis
      DO 850 J=1,3
         XMAT1(J,1) = V1(J)
         XMAT1(J,2) = V3(J)
         XMAT1(J,3) = V4(J)
850   CONTINUE
!----Make a copy and invert it
      DO 851 I=1,3
         DO 851 J=1,3
            XMAT2(I,J) = XMAT1(I,J)
851   CONTINUE
      CALL MINV (XMAT1, 3, DET, J1, J2)
      XHOR = XHL*SIND(109.5 - 90.)
      XVER = XHL*COSD(109.5 - 90.)
!----Produce a H location
      XYZT(1) = XHOR                    ! along V1
      XYZT(2) = 0.0
      XYZT(3) = -XVER                   ! along -V4
      IF (IH_TYPE .NE. 4) GO TO 852
         CALL OPTIM_CH3 (BEST_ANG)
         GO TO 903
!----Rotate XYZT to form the three H positions
852   XMAT3(1,1) = 1.0
      ANGLE = -120.0 + ROTN
      DO 900 I=1,3
         ANGLE = ANGLE + 120.
         XMAT3(2,2) = COSD(ANGLE)
         XMAT3(2,3) = -SIND(ANGLE)
         XMAT3(3,2) = SIND(ANGLE)
         XMAT3(3,3) = COSD(ANGLE)
         CALL MTIMES2 (XYZT, XMAT3, XYZR(1,I))
900   CONTINUE
!----Convert back to original system
903   DO 920 J=1,3
         CALL MTIMES2 (XYZR(1,J), XMAT2, XYZR3(1,J))   ! in cart system
         DO K=1,3
            XYZR4(K,J) = C1(K) + XYZR3(K,J)             ! cart + C1
         ENDDO
920   CONTINUE
!----Convert the H position back to frac coords
      DO 950 I=1,3
         CALL MTIMES2 (XYZR4(1,I), T, XYZH(1,I))
950   CONTINUE
      WRITE (6,960) ID2
960   FORMAT (' Three H atoms linked to sp3 ',A5,' are at...')
      DO 914 L2=1,3
         NEW = NEW + 1
         L = NATOMS + NEW
         CALL MAKE_H_LABEL(L)                           ! 10/9/03
         DO 915 I=1,3
            XYZ_NEW(I,NEW) = XYZH(I,L2)        ! fractional
            OR_XYZ(I,L) = XYZR4(I,L2)           ! cartesian
            XYZ(I,L) = XYZH(I,L2)              ! fractional
915      CONTINUE
         WRITE (6,909) LABEL(L), (XYZH(I,L2),I=1,3),   &
     &                           (XYZR4(I,L2),I=1,3)
909      FORMAT (2X,A5,' fract coord:',3F7.4,'; orthog coord:',3F7.3)
914   CONTINUE
      IF (IH_TYPE .EQ. 4) WRITE (6,910) BEST_ANG
910   FORMAT (4X,'Coordinates have been optimized by a',F6.1,  &
     &          ' deg rotation from initial orientation')
      GO TO 300
1000  IF (IGEN .NE. 1) GO TO 1001
      PRINT 1002
1002    FORMAT (' The following C-hydrogens have been generated by',   &
     &     ' a sym op from previously idealized atoms...')
      N = NEW
      DO 1914 L2 = 1,N
         NEW = NEW + 1
         L = NATOMS + NEW
         CALL MAKE_H_LABEL(L)                           ! 10/9/03
         DO 1916 I=1,3
            XYZ_NEW(I,NEW) = TR(I)
            DO 1915 J=1,3
               XYZ_NEW(I,NEW) = XYZ_NEW(I,NEW) +  RS(I,J)*XYZ_NEW(J,L2)   ! fractional
1915        CONTINUE
            XYZ(I,L) = XYZ_NEW(I,NEW)                      ! fractional
1916     CONTINUE
         DO 1913 I=1,3                                      ! cartesian from
            OR_XYZ(I,L) = 0.0                                ! just made
            DO 1913 J=1,3                                   ! fractional
               OR_XYZ(I,L) = OR_XYZ(I,L) + R(I,J)*XYZ(J,L)    ! coordinates
1913     CONTINUE
!----Check to see if a newly generated H is .le. 0.1 Ang from one that
!     already exists...like on a 2-fold
         DO 2600 M=1,L-1
            IF (DISTANCE(M, L) .GT. 0.1) GO TO 2600
            NEW = NEW -1                     ! eliminate NEWth atom,
            GO TO 1914                       ! it's a duplicate by symmetry
2600     CONTINUE
         WRITE (6,909) LABEL(L), (XYZ(I,L),I=1,3),  &
     &                           (OR_XYZ(I,L),I=1,3)
1914  CONTINUE
!
!----Write CHEMX format file...fractional coordinate file
1001  PRINT 995
995   FORMAT (' Fractional coordinate CHEMX file with adjusted',  &
     &        ' coordinates written as fort.15')
      OPEN (UNIT=15, STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE (15,11) CELL
11    FORMAT (38X,3F8.3/21X,3F8.3)
      NATOMS = NATOMS + NEW
      WRITE (15,14) NATOMS
14    FORMAT (I4,'   0')
!----Skip 1 line
      WRITE (15,*)
      DO 918 J=1,NATOMS
         WRITE (15,916) J, LABEL(J)(1:1), (XYZ(I,J),I=1,3)
916      FORMAT (I4,1X,A1,4X,3F10.5)
918   CONTINUE
      CLOSE (UNIT = 15)
      RETURN
      END SUBROUTINE ADJUST
!-----------------------------------------------------------------
!----Form a vector
      SUBROUTINE VECTOR (V, A, B)
      REAL    :: V(3), A(3), B(3)
      INTEGER :: I
      DO I=1,3
         V(I) = A(I) - B(I)
      ENDDO
      RETURN
      END SUBROUTINE VECTOR
!----Normalize a vector
      SUBROUTINE N_VECTOR (V)
      REAL    :: V(3), S
      INTEGER :: I
      S = 0.0
      DO I=1,3
         S = S + V(I)**2
      ENDDO
      S = SQRT(S)
      DO I=1,3
         V(I) = V(I)/S
      ENDDO
      RETURN
      END SUBROUTINE N_VECTOR
!----Vector cross product
      SUBROUTINE V_CROSS (V, A, B)               ! V = A x B
      REAL    :: V(3), A(3), B(3), X
      INTEGER :: IN1(3), IN2(3), I
      DATA IN1 / 2,1,1/
      DATA IN2 /3,3,2/
      X = -1.0
      DO I=1,3
         X = -X
         V(I) = X*(A(IN1(I))*B(IN2(I)) -   &
     &             B(IN1(I))*A(IN2(I)))
      ENDDO
      CALL N_VECTOR(V)
      RETURN
      END SUBROUTINE V_CROSS
!----Matrix for coordinate system transformation
      SUBROUTINE MATRXT (A,T)
      REAL    :: T(3,3), A(6) , W
      T(1,1)=A(1)
      T(1,2)=A(2)*COSD(A(6))
      T(1,3)=A(3)*COSD(A(5))
      T(2,1)=0.0
      T(2,2)=A(2)*SIND(A(6))
      T(2,3)=A(3)*(COSD(A(4))-COSD(A(5))*COSD(A(6)))/SIND(A(6))
      W=SQRT(1-COSD(A(4))**2-COSD(A(5))**2-COSD(A(6))**2+  &
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
      REAL ::  XYZT(3),XYZG(3),RV(3,3)
      XYZT(1)=RV(1,1)*XYZG(1)+RV(1,2)*XYZG(2)   &
     &+RV(1,3)*XYZG(3)
      XYZT(2)=RV(2,1)*XYZG(1)+RV(2,2)*XYZG(2)   &
     &+RV(2,3)*XYZG(3)
      XYZT(3)=RV(3,1)*XYZG(1)+RV(3,2)*XYZG(2)   &
     &+RV(3,3)*XYZG(3)
      RETURN
      END SUBROUTINE MTIMES2
!----------
      SUBROUTINE MINV (A,N,D,L,M)                                       !MINV   5
        IMPLICIT NONE
!                                                                       !MINV  15
!     ..................................................................!MINV  20
!                                                                       !MINV  25
!        SUBROUTINE MINV                                                !MINV  30
!                                                                       !MINV  35
!        PURPOSE                                                        !MINV  40
!           INVERT A MATRIX                                             !MINV  45
!                                                                       !MINV  50
!        USAGE                                                          !MINV  55
!           CALL MINV(A,N,D,L,M)                                        !MINV  60
!                                                                       !MINV  65
!        DESCRIPTION OF PARAMETERS                                      !MINV  70
!           A - INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY  !MINV  70
!               RESULTANT INVERSE.                                      !MINV  80
!           N - ORDER OF MATRIX A                                       !MINV  85
!           D - 0.0 IF SINGULAR, 1.0 IF ALL RIGHT.                      !MINV  90
!           L - WORK VECTOR OF LENGTH N                                 !MINV  95
!           M - WORK VECTOR OF LENGTH N                                 !MINV 100
!                                                                       !MINV 105
!        REMARKS                                                        !MINV 110
!           MATRIX A MUST BE A GENERAL MATRIX                           !MINV 115
!                                                                       !MINV 120
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  !MINV 125
!           NONE                                                        !MINV 130
!                                                                       !MINV 135
!        METHOD                                                         !MINV 140
!           THE STANDARD GAUSS-JORDAN METHOD IS USED.                   !MINV 145
!                                                                       !MINV 150
!     ..................................................................!MINV 155
!                                                                       !MINV 160
      REAL    :: A(1), D
      INTEGER :: L(1),M(1), N                                           !MINV 165
!                                                                       !MINV 170
!        ...............................................................!MINV 175
!                                                                       !MINV 180
!        SEARCH FOR LARGEST ELEMENT                                     !MINV 185
!                                                                       !MINV 190
      INTEGER :: I, IJ, IK, IZ, J, JI, JK, JP, JQ, JR
      INTEGER :: K, KI, KJ, KK, NK
      REAL    :: BIGA, HOLD
!
      D=1.0                                                             !MINV 195
      NK=-N                                                             !MINV 200
      DO 590 K=1,N                                                      !MINV 205
      NK=NK+N                                                           !MINV 210
      L(K)=K                                                            !MINV 215
      M(K)=K                                                            !MINV 220
      KK=NK+K                                                           !MINV 225
      BIGA=A(KK)                                                        !MINV 230
      DO 510 J=K,N                                                      !MINV 235
      IZ=N*(J-1)                                                        !MINV 240
      DO 510 I=K,N                                                      !MINV 245
      IJ=IZ+I                                                           !MINV 250
  500 IF(ABS(BIGA)-ABS(A(IJ))) 505,510,510                              !MINV 255
  505 BIGA=A(IJ)                                                        !MINV 260
      L(K)=I                                                            !MINV 265
      M(K)=J                                                            !MINV 270
  510 CONTINUE                                                          !MINV 270
!                                                                       !MINV 280
!        INTERCHANGE ROWS                                               !MINV 285
!                                                                       !MINV 290
      J=L(K)                                                            !MINV 295
      IF(J-K) 525,525,515                                               !MINV 300
  515 KI=K-N                                                            !MINV 305
      DO 520 I=1,N                                                      !MINV 310
      KI=KI+N                                                           !MINV 315
      HOLD=-A(KI)                                                       !MINV 320
      JI=KI-K+J                                                         !MINV 325
      A(KI)=A(JI)                                                       !MINV 330
  520 A(JI) =HOLD                                                       !MINV 335
!                                                                       !MINV 340
!        INTERCHANGE COLUMNS                                            !MINV 345
!                                                                       !MINV 350
  525 I=M(K)                                                            !MINV 355
      IF(I-K) 540,540,530                                               !MINV 360
  530 JP=N*(I-1)                                                        !MINV 365
      DO 535 J=1,N                                                      !MINV 370
      JK=NK+J                                                           !MINV 375
      JI=JP+J                                                           !MINV 380
      HOLD=-A(JK)                                                       !MINV 385
      A(JK)=A(JI)                                                       !MINV 390
  535 A(JI) =HOLD                                                       !MINV 395
!                                                                       !MINV 400
!        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS        !MINV 405
!        CONTAINED IN BIGA)                                             !MINV 410
!                                                                       !MINV 415
  540 IF(BIGA) 550,545,550                                              !MINV 420
  545 D=0.0                                                             !MINV 425
      RETURN                                                            !MINV 430
  550 DO 560 I=1,N                                                      !MINV 435
      IF(I-K) 555,560,555                                               !MINV 440
  555 IK=NK+I                                                           !MINV 445
      A(IK)=A(IK)/(-BIGA)                                               !MINV 450
  560 CONTINUE                                                          !MINV 455
!                                                                       !MINV 460
!        REDUCE MATRIX                                                  !MINV 465
!                                                                       !MINV 470
      DO 575 I=1,N                                                      !MINV 470
      IK=NK+I                                                           !MINV 480
      IJ=I-N                                                            !MINV 485
      DO 575 J=1,N                                                      !MINV 490
      IJ=IJ+N                                                           !MINV 495
      IF(I-K) 565,575,565                                               !MINV 500
  565 IF(J-K) 570,575,570                                               !MINV 505
  570 KJ=IJ-I+K                                                         !MINV 510
      A(IJ)=A(IK)*A(KJ)+A(IJ)                                           !MINV 515
  575 CONTINUE                                                          !MINV 520
!                                                                       !MINV 525
!        DIVIDE ROW BY PIVOT                                            !MINV 530
!                                                                       !MINV 535
      KJ=K-N                                                            !MINV 540
      DO 585 J=1,N                                                      !MINV 545
      KJ=KJ+N                                                           !MINV 550
      IF(J-K) 580,585,580                                               !MINV 555
  580 A(KJ)=A(KJ)/BIGA                                                  !MINV 560
  585 CONTINUE                                                          !MINV 565
!                                                                       !MINV 570
!        REPLACE PIVOT BY RECIPROCAL                                    !MINV 575
!                                                                       !MINV 580
      A(KK)=1.0/BIGA                                                    !MINV 585
  590 CONTINUE                                                          !MINV 590
!                                                                       !MINV 595
!        FINAL ROW AND COLUMN INTERCHANGE                               !MINV 600
!                                                                       !MINV 605
      K=N                                                               !MINV 610
  595 K=(K-1)                                                           !MINV 615
      IF(K) 630,630,600                                                 !MINV 620
  600 I=L(K)                                                            !MINV 625
      IF(I-K) 615,615,605                                               !MINV 630
  605 JQ=N*(K-1)                                                        !MINV 635
      JR=N*(I-1)                                                        !MINV 640
      DO 610 J=1,N                                                      !MINV 645
      JK=JQ+J                                                           !MINV 650
      HOLD=A(JK)                                                        !MINV 655
      JI=JR+J                                                           !MINV 660
      A(JK)=-A(JI)                                                      !MINV 665
  610 A(JI) =HOLD                                                       !MINV 670
  615 J=M(K)                                                            !MINV 670
      IF(J-K) 595,595,620                                               !MINV 680
  620 KI=K-N                                                            !MINV 685
      DO 625 I=1,N                                                      !MINV 690
      KI=KI+N                                                           !MINV 695
      HOLD=A(KI)                                                        !MINV 700
      JI=KI-K+J                                                         !MINV 705
      A(KI)=-A(JI)                                                      !MINV 710
  625 A(JI) =HOLD                                                       !MINV 715
      GO TO 595                                                         !MINV 720
  630 RETURN                                                            !MINV 725
      END SUBROUTINE MINV                                               !MINV 730
!------------------------------------------------------------------
      SUBROUTINE OPTIM_CH3 (ANG_MIN)              ! 3/2/94
!
!----Rotate a newly formed methyl group about the C-C bond and
!     find minimum repulsion with other atoms with 5 Angs.
!
      USE mod_preppot, ONLY : LABEL,  OR_XYZ, NATOMS, &
     &                NEW, XYZT, N_2, N_CENTRAL, XYZR, &
     &                XMAT2, C1
!O_XYZ --> OR_XYZ and ID --> LABEL
      IMPLICIT NONE
      REAL :: ANG_MIN
!
! CHARACTER (LEN=5) :: ID2, ID(500), NAME(500) removed, not in use
!I_ALTERED(500), RS(3,3), TR(3), N_NEAR(100), N_ATOMS(4,100), J1(3), J2(3),
!C(3,500), C2(3), C3(3), C4(3) removed, not in used

      REAL    ::  CC(2), B(2), XMAT3(3,3), XYZR4(3,3), XYZR3(3,3)
!VAVG(3), V1(3), V2(3), V3(3), V4(3), VH(3), VX(3), XMAT1(3,3) and
!RINV(3,3), T(3,3) removed, not in use

!
!      REAL :: C not in use
!XYZ_NEW(3,500), XYZH(3,500) removed, not in use
      DATA B / 55.45, 269.20 /                ! B term for H, C
      DATA CC / 1.87, 1.80   /                ! C term for H, C
!
      INTEGER :: I, IA, IT, J, K, L, NEWANDOLD
      REAL    :: A, ANGLE, D, E, E_MIN, SUME
!                                  repl E = B1*B2*exp(-(C1+C2)*R)
!
      NEWANDOLD = NATOMS + NEW
      E_MIN = 9999999.
      ANG_MIN = 0.0
!
!----Rotate XYZT to form the three H positions
      XMAT3(1,1) = 1.0
      DO 100 IA=0,119
         ANGLE = IA
         SUME = 0.0
         DO 9 I=1,3
            A = ANGLE + (I - 1)*120.
            XMAT3(2,2) = COSD(A)
            XMAT3(2,3) = -SIND(A)
            XMAT3(3,2) = SIND(A)
            XMAT3(3,3) = COSD(A)
            CALL MTIMES2 (XYZT, XMAT3, XYZR(1,I))
9        CONTINUE
!----Convert back to original system
         DO 20 J=1,3
            CALL MTIMES2 (XYZR(1,J), XMAT2, XYZR3(1,J))   ! in cart system
            DO K=1,3
               XYZR4(K,J) = C1(K) + XYZR3(K,J)             ! cart + C1
            ENDDO
20       CONTINUE
!----Look at interaction between the 3 new H's and all other atoms within
!     5 Angs.  Discard the two C's that define the C-CH3 bond.
         DO 90 J=1,NEWANDOLD
            IF (J .EQ. N_2 .OR. J .EQ. N_CENTRAL) GO TO 90
            DO 80 K=1,3
               D = 0.0                                  ! distance from methyl
               DO 75 L=1,3                              ! H to another
                  D = D + (XYZR4(L,K) -  OR_XYZ(L,J))**2 ! atom
75             CONTINUE
               D = SQRT(D)
               IF (D .GT. 6.0) GO TO 80             ! if d > 6.0 Angs...quit
               IF (LABEL(J)(1:1) .EQ. 'H') THEN
                  IT = 1                         ! other atom is a H
               ELSE
                  IT = 2                         ! other atom is a C
               ENDIF
               E = B(1)*B(IT)*EXP(-(CC(1) + CC(IT))*D)
               SUME = SUME + E
80          CONTINUE
90       CONTINUE
         IF (SUME .GE. E_MIN) GO TO 100
            E_MIN = SUME
            ANG_MIN = ANGLE
100   CONTINUE
!----Pass back the optimized CH3 H's
         DO 200 I=1,3
            A = ANG_MIN + (I-1)*120.0
            XMAT3(2,2) = COSD(A)
            XMAT3(2,3) = -SIND(A)
            XMAT3(3,2) = SIND(A)
            XMAT3(3,3) = COSD(A)
            CALL MTIMES2 (XYZT, XMAT3, XYZR(1,I))
200      CONTINUE
      RETURN
      END SUBROUTINE OPTIM_CH3
!------------------------------------------------------------------
      SUBROUTINE CONNECTION (DMAX)             ! 8/9/94
!        Add S and I                            11/6/03
!
!----Determine connectivity
!
      USE mod_preppot, ONLY : LABEL, NATOMS
      USE mod_preppot, ONLY : N_NEAR, N_ATOMS
      IMPLICIT NONE
      REAL :: DMAX
!
!AH*1, GEN*1,ID2*5, ID(500)*5, NAME(500)*5 removed, not in use
!
      INTEGER :: I, J
      REAL    :: D, DISTANCE
!----Determine connectivity....use max distance of DMAX
      DO 502 I=1,NATOMS
         N_NEAR(I) = 0
         DO 490 J=1,NATOMS
            DMAX = 1.65                                  ! 11/6/03
            IF (I .EQ. J) GO TO 490
            IF (LABEL(I)(1:1) .EQ. 'H' .AND.   &
     &          LABEL(J)(1:1) .EQ. 'H') GO TO 490
            D = DISTANCE (I, J)
            IF ((LABEL(I)(1:1) .EQ. 'C' .AND.  &         ! 11/6/03
     &           LABEL(J)(1:1) .EQ. 'I') .OR.  &         ! 11/6/03
                (LABEL(I)(1:1) .EQ. 'I' .AND.  &         ! 11/6/03
     &           LABEL(J)(1:1) .EQ. 'C')) DMAX = 2.2     ! 11/6/03
            IF (LABEL(I)(1:2) .EQ. 'BR' .OR.   &
     &          LABEL(J)(1:2) .EQ. 'BR') DMAX = 2.04     ! 11/6/03
            IF (LABEL(I)(1:1) .EQ. 'S' .OR.   &
     &          LABEL(J)(1:1) .EQ. 'S') DMAX = 1.85      ! 11/6/03
            IF (D .GT. DMAX) GO TO 490
470         N_NEAR(I) = N_NEAR(I) + 1                    ! 11/6/03
            N_ATOMS(N_NEAR(I),I) = J
490      CONTINUE
502   CONTINUE
      RETURN
      END SUBROUTINE CONNECTION
!
!------------------------------------------------------------
      SUBROUTINE MAKE_CHEM3D              ! 6/10/05
!      added sulfur                         11/2/00
!      added chlorine                       6/10/05
!
!----Prepare CHEM3D file from CSD entry with line numbers
!     as Chem3D connectivity codes...output on unit # 14
!
      USE mod_preppot, ONLY : LABEL, OR_XYZ, NATOMS
      USE mod_preppot, ONLY : ICON, NCON, AT2
      IMPLICIT NONE
!
      CHARACTER (LEN=21) :: ATOMS_OUT(500)
! NAME removed, not in use
! NCON dimension is 200 in module
!
      INTEGER :: I, IS, J, K, L, NA, N_ADJ, N_OXY
!
      DO I=1,NATOMS       ! it is assumed that element symbol
        NCON(NA) = 0      ! has only 1 character in 1st position
      ENDDO               ! in the atom name...like C(10)
!
!----Build the build the ICON array
      CALL CCONNECT                          ! build ICON
!
      WRITE (14,10) NATOMS    ! # atoms as 1st line
10    FORMAT (I3)
!
!----Show the equivalence between csd and chem3d atom names
      WRITE (23,12)
12    FORMAT (' Equivalence of CHEM3D and CSD atom names...'/  &
     &        ' Chem3d name  Csd name '/                       &
     &        ' -----------  -------- ')
      DO 1303 NA=1,NATOMS
      k = 0
!----How many connected?
        IF( NCON(NA) .NE. 0) GO TO 1375
        DO L=1,6
           IF (ICON(L,NA) .EQ. 0) GO TO 1375
           NCON(NA) = NCON(NA) + 1
        ENDDO
1375    IF (AT2(NA) .EQ. 'H') THEN                    ! H
           l = icon(na,1)
           if (label(l)(1:1) .eq. 'N' .and. ncon(l) .eq. 4) then
              k = 48                 ! H lined to N+
              go to 158
           endif
           l = icon(na,1)
           if (label(l)(1:1) .eq. 'O' .and. ncon(l) .eq. 2) then
              k = 21                 ! H lined to -O-
              go to 158
           endif
              K = 5
              GO TO 158
        ENDIF
!
        IF (LABEL(NA)(1:2) .EQ. 'BR')  THEN       ! Bromine
           if (ncon(na) .ne. 0) then
              K = 13                 ! C-Br
           else
              k = 350                ! Br- (anion)
           endif
           go to 158
        endif
!
        IF (LABEL(NA)(1:2) .EQ. 'CL')  THEN       ! chlorine
           if (ncon(na) .ne. 0) then
              K = 12                 ! C-Cl
           else
              k = 170                ! Cl- (anion)
           endif
        go to 158
        endif
!
        IF (LABEL(NA)(1:1) .EQ. 'I')  THEN       ! iodine
           if (ncon(na) .ne. 0) then
              K = 14                 ! C-I
           else
              k = 530                ! I- (anion)
           endif
        go to 158
        endif
!
        IF (LABEL(NA)(1:1) .EQ. 'S')  THEN            ! S
           K = 18
           GO TO 158
        ENDIF
        IF (AT2(NA) .EQ. 'F') THEN                    ! F
           K = 11
           GO TO 158
        ENDIF
        IF (at2(NA) .EQ. 'C') THEN                      ! 6/10/05
           GO TO (1001, 1001, 1003, 1004), NCON(NA)
1001       K = 4                                        ! alkyne C
           GO TO 158
1003       K = 2                                        ! alkene C
           GO TO 158
1004       K = 1                                        ! alkane C
           GO TO 158
        ENDIF
        IF (AT2(NA) .EQ. 'S') THEN
           GO TO (1011, 1011, 1013, 1014), NCON(NA)
1011       K = 15                                       ! thio ether
           GO TO 158
1013       K = 17                                       ! suloxide S
           GO TO 158
1014       K = 18                                       ! sulfone S
        ENDIF
        IF (AT2(NA) .EQ. 'O') THEN
           GO TO (1021, 1022), NCON(NA)
1021          N_ADJ = ICON(1,NA)           ! get ID of single connected atom
              K = 7                        ! default is carbonyl
              IF (AT2(N_ADJ) .EQ. 'C') GO TO 158
              IF (NCON(N_ADJ) .NE. 3) GO TO 158
                 N_OXY = 0                 ! count # O's connected to N_ADJ
                 DO I=1,3
                    IS = ICON(I,N_ADJ)
                    IF (AT2(IS) .EQ. 'O') N_OXY = N_OXY + 1
                 ENDDO
                 IF (N_OXY .NE. 2) GO TO 158
                    K = 47                              ! nitro O
                 GO TO 158
1022       K = 6                                        ! ether O
           GO TO 158
        ENDIF
        IF (AT2(NA) .EQ. 'N') THEN
           GO TO (1031, 1032, 1033, 1034), NCON(NA)
1031       K = 10                                       ! nitrile N
           GO TO 158
1032       K = 40                                       ! enamine N
           GO TO 158
1033       K = 8                                        ! amine N default
              N_OXY = 0                    ! count # of O's on this N
              DO I=1,3
                 N_ADJ = ICON(I,NA)
                 IF (AT2(N_ADJ) .EQ. 'O') N_OXY = N_OXY + 1
              ENDDO
              IF (N_OXY .EQ. 2) K = 46                  ! nitro N
           GO TO 158
1034       K = 39                                       ! ammonium N
           GO TO 158
        ENDIF
!----Chem3d format file (cc1, with serial numbers on unit # 14)
158   if (label(na)(1:2) .eq. 'CL' .or. label(na)(1:2) .eq. 'BR') then
         WRITE(14,1169) LABEL(NA)(1:2), NA, (OR_XYZ(J,NA),J=1,3), K,  &
     &             (ICON(J,NA),J=1,NCON(NA))
1169     FORMAT (2X,A2,I4,3F12.6,7I5)
         WRITE (ATOMS_OUT(NA), 2161) label(na)(1:2), NA, LABEL(NA)
         WRITE (23,2161) label(na)(1:2), NA, LABEL(NA)
2161     FORMAT (4x,A2,I2,7X,A5,1X)
       go to 1303
       endif
168   IF (LABEL(NA)(2:2) .NE. '_') THEN
         WRITE(14,1159) LABEL(NA)(1:2), NA, (OR_XYZ(J,NA),J=1,3), K,  &
     &             (ICON(J,NA),J=1,NCON(NA))
1159     FORMAT (1X,A2,I5,3F12.6,7I5)
         WRITE (ATOMS_OUT(NA), 1161) LABEL(NA)(1:2), NA, LABEL(NA)
         WRITE (23,1161) LABEL(NA)(1:2), NA, LABEL(NA)
1161     FORMAT (4x,A2,I2,7X,A5,1X)
      ELSE
         WRITE(14,159) AT2(NA), NA, (OR_XYZ(J,NA),J=1,3), K,  &
     &             (ICON(J,NA),J=1,NCON(NA))
159      FORMAT (2X,A1,I5,3F12.6,7I5)
         WRITE (ATOMS_OUT(NA), 161) AT2(NA), NA, LABEL(NA)
         WRITE (23,161) AT2(NA), NA, LABEL(NA)
161      FORMAT (4x,A1,I3,7X,A5,1X)
      ENDIF
1303  CONTINUE
!
!----Show list of original atom names and chem3d names
      PRINT 1308
1308  FORMAT (' Equivalence of CHEM3D and CSD atom names...'/  &
     &        3(' Chem3d name  Csd name ')/                    &
     &        3(' -----------  -------- '))
      PRINT 1312, (ATOMS_OUT(L),L=1,NA)
1312  FORMAT (3(1X,A21,1X))
      CLOSE (UNIT=14)
      RETURN
      END SUBROUTINE MAKE_CHEM3D
!---------------------------------------------------------
!
!----CONNECT...build ICON array
      SUBROUTINE CCONNECT
!
      USE mod_preppot, ONLY : LABEL, NATOMS
      USE mod_preppot, ONLY : ICON, NCON, AT2
      IMPLICIT NONE
!
! NCON dimention is set 200 in module
!
      INTEGER :: I, J
      REAL    :: D, DISTANCE
!
      DO I=1,NATOMS
         AT2(I) = LABEL(I)(1:1)
      ENDDO
!
      DO 100 I=1,NATOMS
         NCON(I) = 0
         DO 90 J=1,NATOMS
            IF (I .EQ. J) GO TO 90
            IF (AT2(I) .EQ. 'H' .AND. AT2(J) .EQ. 'H')   &
     &                         GO TO 90
            D = DISTANCE(I, J)
            IF (LABEL(I)(1:2) .EQ. 'BR' .OR.             &
     &          LABEL(J)(1:2) .EQ. 'BR') THEN
               IF (D .LE. 2.0) GO TO 89
            ELSE
            IF (LABEL(I)(1:1) .EQ. 'S' .OR.              &
     &          LABEL(J)(1:1) .EQ. 'S') THEN
               IF (D .LE. 1.85) GO TO 89
            ELSE
!              IF (D .LE. 1.9) GO TO 89
               IF (D .LE. 1.7) GO TO 89       ! changed on 10-30-97
            END IF
            END IF
            GO TO 90
89          NCON(I) = NCON(I) + 1
            ICON(NCON(I),I) = J
90       CONTINUE
100   CONTINUE
      RETURN
      END SUBROUTINE CCONNECT
!
!--------------------------------------------------------------
!----Set an atom along a bond vector at a particular distance
!
      SUBROUTINE STRETCH (I_FIX, I_MOVE, B_LENGTH)
!
!
      USE mod_preppot, ONLY : LABEL, XYZ, OR_XYZ
      USE mod_preppot, ONLY : T
      IMPLICIT NONE
!
      INTEGER :: I, I_FIX, I_MOVE
      REAL    :: B_LENGTH, D, DISTANCE, RATIO
!
      D = DISTANCE (I_MOVE, I_FIX)
      RATIO = B_LENGTH/D
      DO 10 I=1,3
        OR_XYZ(I,I_MOVE) = (OR_XYZ(I,I_MOVE) - OR_XYZ(I,I_FIX))*RATIO  &
     &                   + OR_XYZ(I,I_FIX)
10    CONTINUE
      PRINT 20, LABEL(I_FIX), LABEL(I_MOVE), D,    &
     &          (DISTANCE(I_MOVE, I_FIX))
      WRITE (23,20) LABEL(I_FIX), LABEL(I_MOVE), D, &
     &          (DISTANCE(I_MOVE, I_FIX))
20      FORMAT (' Bond vector for ',A5,'-',A5,' changed from',  &
     &          F6.3,' to',F6.3,' Angs')
!----Update the fractional coordinates with the moved atom
      CALL MTIMES2 (OR_XYZ(1,I_MOVE), T, XYZ(1,I_MOVE))
      RETURN
      END SUBROUTINE STRETCH
!
!---------------------------------------------------------------
!
!----calc distance between 2 atoms
!
      FUNCTION DISTANCE (I1, I2)
!
      USE mod_preppot, ONLY : OR_XYZ
      IMPLICIT NONE
      INTEGER :: I1, I2
      REAL    :: DISTANCE
!
      INTEGER :: I
!
      DISTANCE = 0.0
      DO I=1,3
         DISTANCE = DISTANCE + (OR_XYZ(I,I1) - OR_XYZ(I,I2))**2
      ENDDO
      DISTANCE = SQRT(DISTANCE)
      RETURN
      END FUNCTION DISTANCE
!--------------------------------------------------------------
!
!----generate new H atom label
!
      SUBROUTINE MAKE_H_LABEL (L)
!
      USE mod_preppot, only : LABEL
      IMPLICIT NONE
!
      INTEGER :: L
!
!----Check for size of L                              ! 10/9/03
      IF (L .GE. 10) THEN                             ! 10/9/03
         WRITE (LABEL(L),345) L                       ! 10/9/03
345      FORMAT ('H_',I2,' ')                         ! 10/9/03
      ELSE                                            ! 10/9/03
         WRITE (LABEL(L),346) L                       ! 10/9/03
346      FORMAT ('H_',I1,'  ')                        ! 10/9/03
      ENDIF                                           ! 10/9/03
      RETURN                                          ! 10/9/03
      END SUBROUTINE MAKE_H_LABEL                     ! 10/9/03
!
!--------------------------------------------------------------
!
!----Function to determine if a 2-linked sulfur is             1/31/06
!     part of a thiazole.                                      1/31/06
!
      FUNCTION IQ_thiazole_s(I1) ! I1 is a 2-linked S        ! 1/31/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_thiazole_s, iq_thiazole_n, nc, i, j
      integer :: i1, i2, i3, id_n(2)
!
      IQ_thiazole_s = 0     ! no
!
!          78 --> N---C
!                 "   "
!                 C   C
!                  \ /
!           77 -->  S
!
      NC = 0                             ! count # of 3-linked C's
      DO I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then
             NC = NC + 1
             id_n(nc) = i2               ! id of 3-linked C
         endif
      enddo
      IF (NC .ne. 2) RETURN              ! not a possible # 79 S
!----Look at two 3-linked C's...one should be bonded to a # 78 N
      DO 120 I=1,2
         i2 = id_n(i)      ! id of 3-linked C
         do 130 j=1,3
            i3 = icon(j,i2)
            IF (LABEL(I3)(1:1) .eq. 'N' .and. ncon(i3) .eq. 2) then
              if (iq_thiazole_n(i3) .eq. 78) then
                 iq_thiazole_s = 77     ! param = 77
                 return      ! found the S
              endif
            endif
130      continue
120   CONTINUE
!
      return           ! returning here means no luck
!
      end function iq_thiazole_s                             ! 1/31/06
!
!-----------------------------------------------------------------
!
!----Function to determine if a 2-linked nitrogen is           1/31/06
!     part of a 1,3-thiazole.                                  1/31/06
!
      FUNCTION IQ_thiazole_n(I1) ! I1 is a 2-linked N         ! 2/2/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_thiazole_n, nc, i, j, k, j1, j2, k1, k2
      integer :: i1, i2, id_n(2), kk1, kk2, l3
!
!          78 --> N---C
!                 "   "
!                 C   C
!                  \ /
!           77 -->  S
!
      IQ_thiazole_n = 0     ! no
!
      NC = 0                             ! count # of 3-linked C's
      DO I=1,2                           ! examine the N
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then
             NC = NC + 1
             id_n(nc) = i2               ! id of 3-linked C
         endif
      enddo
      IF (NC .ne. 2) RETURN              ! not a possible # 78 N
!----Look at two 3-linked C's...one should be bonded to 3-linked
!     C and the other to a 2-linked S
      DO 120 i=1,3                ! examine first of 3-linked C's
         j1 = icon(i,id_n(1))     ! id's of atoms linked to first C
         do 130 j=1,3             ! examine second of 3-linked C's
            j2 = icon(j,id_n(2))  ! id's of atoms linked to second C
            do 140 k1=1,ncon(j1)
               kk1 = icon(k1,j1)
               do 150 k2=1,ncon(j2)
                  kk2 = icon(k2,j2)
                  IF (LABEL(kk1)(1:1) .eq. 'S' .or. &    ! one of these must
                      label(kk2)(1:1) .eq. 'S') then     ! be S to keep going
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
!
      return           ! returning here means no luck
!
      end function iq_thiazole_n                             ! 1/31/06
!
!-----------------------------------------------------------------
!----Function to determine if a 3-linked C is                   2/3/06
!     carbonyl C of an ester or ahhydride                       2/3/06
!
      FUNCTION IQ_ester_c(I1) ! I1 is a 3-linked C            ! 2/3/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_ester_c, nc, i, j, k, j1, j2, k1, k2
      integer :: i1, i2, no1, no2, ido2(2)
!
!         80 -->  O               O
!                 "               "
!         79 -->  C   C or i-->   C
!                / \ /             \
!               C   O <--81         C
!
      IQ_ester_c = 0     ! no
!
      no1 = 0                           ! count # of 1-linked O's
      no2 = 0                           ! count # of 2-linked O's
      nc = 0                            ! count number of C's
      DO I=1,3
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C') NC = NC + 1
         if (label(i2)(1:1) .eq. 'O' .and. ncon(i2) .eq. 1) no1 = no1 + 1
         if (label(i2)(1:1) .eq. 'O' .and. ncon(i2) .eq. 2) then
            no2 = no2 + 1
            ido2(no2) = i2
         endif
      enddo
      IF (.not. (nc .eq. 1 .and. no1 .eq. 1 .and. no2 .eq. 1)) RETURN ! no
!----Look at 2-linked O...must if bonded to 2 C's
      nc = 0
      DO i=1,2
         i2 = icon(i,ido2(1))     ! id's of atoms linked to 2-linked O
         if (label(i2)(1:1) .eq. 'C') nc = nc + 1
      enddo
!
      if (nc .ne. 2) return     ! no
!
      iq_ester_c = 79      ! param = 79
!
      end function iq_ester_c                       ! 2/3/06
!------------------------------------------------------------------
!----Function to determine if a 2-linked N of                  2/10/06
!     a diazo linkage                                          2/10/06
!
      FUNCTION IQ_diazo_n(I1) ! I1 is a 2-linked N           ! 2/10/06
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_diazo_n, nc, i, nn_2                     ! 2/10/06
      integer :: idn2, i1, i2                                ! 2/10/06
      real :: distance                                       ! 2/10/06
!
!             C--N==N--C                                       2/10/06
!              (82)                                            2/10/06
!
      iq_diazo_n = 0
!
      nc = 0                       ! number of C's
      nn_2 = 0                     ! number of 2-linked n's
      DO I=1,2
         I2 = ICON(I,I1)
         IF (LABEL(I2)(1:1) .EQ. 'C') NC = NC + 1
         if (label(i2)(1:1) .eq. 'N' .and. ncon(i2) .eq. 2) then
            nn_2 = nn_2 + 1
            idn2 = i2
         endif
      enddo
      IF (nc .eq. 1 .and. nn_2 .eq. 1) go to 10
         RETURN ! no
!----Look at the other 2-linked N...must if bonded to a C
10    nc = 0
      DO i=1,2
         i2 = icon(i,idn2)  ! id's of atoms linked to other 2-linked N
         if (label(i2)(1:1) .eq. 'C') nc = nc + 1
      enddo
      if (nc .ne. 1) return     ! no
!
!----One final check...N=N distance must be .le. 1.286 Angs
      if (distance(i1,idn2) .gt. 1.286) return ! too BIG
!
      iq_diazo_n = 82      ! param = 82
!
      end function iq_diazo_n                       ! 2/10/06
!------------------------------------------------------------------
!----Function to determine if a 3-linked nitrogen is a
!     pyridinium-N;   -N(+)-H                                7/20/06
!
      FUNCTION iq_pyridinium_n(I1)    ! I1 is a 3-linked N
!                                       with 1 bond to H
!
      USE mod_preppot, ONLY : ICON, NCON
      USE mod_preppot, ONLY : LABEL
      IMPLICIT NONE
      integer :: iq_pyridinium_n, c3, n2, nh
      integer :: i1, i2, i, idH, id(3), id2(6), id3(6)
      integer :: j, j2, nsave, nsave2, nsave3
!
      iq_pyridinium_n = 0     ! no
!
!----Must be linked to two 2-linked N's
!                   Z
!                  / \
!                 Z   Z
!                 |   |
!                 Z   Z      Z = 3-linked C or 2-linked N
!                  \ /
!           85 -->  N
!                   |
!           86 -->  H
!
      n2 = 0        ! count # of 2-linked N's
      c3 = 0        ! count # of 3-linked C's
      nh = 0
      nsave = 0
      do i=1,3
         id(i) = 0      ! id of 2-linked N or 3-linked C
      enddo
      DO 100 I=1,3                           ! examine the N
         I2 = ICON(I,I1)
         if (LABEL(I2)(1:1) .EQ. 'H') then
            nh = nh + 1
            idH = i2
         endif
         IF (LABEL(I2)(1:1) .EQ. 'C' .and. ncon(i2) .eq. 3) then
             c3 = c3 + 1
             nsave = nsave + 1
             id(nsave) = i2
         endif
         IF (LABEL(I2)(1:1) .EQ. 'N' .and. ncon(i2) .eq. 2) then
            N2 = N2 + 1
            nsave = nsave + 1
            id(nsave) = i2                   ! id of 2-linked N
         endif
100   CONTINUE
      if (nH .ne. 1) return       ! must be 1 H on N
      IF (nsave .ne. 2) return    ! attached N's + C's must = 2
!----id contains the C and/or N #'s of the 2 non-H atoms linked
!     to the pyridinium N
      nsave2 = 0    ! nsave2 will count 2nd atoms out from pyridinium N
      do i=1,2      ! look at attached non-H atoms, eliminate pyridinium N
         i2 = id(i)
            do j=1,ncon(i2)
               j2 = icon(j,i2)
               if (label(j2)(1:1) .eq. 'H') cycle ! don't need H
               if (j2 .eq. i1) cycle              ! i1 = ID of pyridinium N
               if (label(j2)(1:1) .eq. 'N' .and. &
                   ncon(j2) .ne. 2) cycle
               if (label(j2)(1:1) .eq. 'C' .and. &
                   ncon(j2) .ne. 3) cycle
               nsave2 = nsave2 + 1
               id2(nsave2) = j2
            enddo
      enddo
!----Make a list of ID's of atoms attached to 2nd atoms out...if 2 are same
!     then have a 6-membered ring.  These would be 3rd atoms out.
      nsave3 = 0
      do i=1,nsave2
         i2 = id2(i)
         do j=1,ncon(i2)
            j2 = icon(j,i2)
            if (label(j2)(1:1) .eq. 'H') cycle
            if (label(j2)(1:1) .eq. 'N' .and. &
                ncon(j2) .ne. 2) cycle
            if (label(j2)(1:1) .eq. 'C' .and. &
                ncon(j2) .ne. 3) cycle
            nsave3 = nsave3 + 1
            id3(nsave3) = j2  ! ID's of all C and N atoms connected
         enddo                !  to 2nd atoms out
      enddo
      if (nsave3 .lt. 2) return
!----To establish 6-membered ring, two of atoms in id3 must be the same
      do i=1,nsave3-1
         do j=i+1,nsave3
            if (id3(i) .eq. id3(j)) go to 110
         enddo
      enddo
!----dropping thru the loop means we've failed to find a ring
      return
110   iq_pyridinium_n = 85   ! pyridinium N+, param = 85
      return
!
      end function iq_pyridinium_n
!---------------------------------------------------------------------
