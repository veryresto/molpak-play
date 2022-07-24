  module make_filesCommonMod
! updated on 9-18-07 
! updated it on 1-1-2007 for adding new S.G codes (AS,AU)
  implicit none

! For /id/...
  character(1) ::  compd_id(40)
  character(40) :: compdid
  character(40) :: potential_file
  integer :: idchg, kind, nid
  equivalence (compd_id, compdid)
  
! For /moldml/...
  character(1) :: spgr(40), elmspgr(100)
  character(2) :: elm_spgr(50), group
  character(73) :: title
  equivalence (group, spgr), (elm_spgr, elmspgr)
  
! Big common for main pgm and ref_dmarel

      INTEGER, PARAMETER :: N_SPGR = 54        ! 9-18-07 
      
      LOGICAL I_COMP, I_SAVE
       
      CHARACTER (1) :: ATOM_LIST(40)
      CHARACTER (1) :: BLANK = ' '
      CHARACTER (1) :: BLNK(26) 
      CHARACTER (1) :: COMPLT(16) = (/'_','C','O','M','P','L','E','T','E',&
     &                                '.','O','U','T','P','U','T'/), FMT(40)
      CHARACTER (1) :: G92(3) = (/'g','9','4'/)
      CHARACTER (1) :: G92FILE(37), I_COMPLETE
      CHARACTER (1) :: INP(4) = (/'.','I','N','P'/)
      CHARACTER (1), DIMENSION (40) :: INPUT_LINE2 = ' '
      CHARACTER (1) :: INVER1, I_SAVED, LINE_OUT1(80), NEWFILE(37) 
      CHARACTER (1) :: MASTERCOM(11)
      CHARACTER (1), DIMENSION (41) ::  MASTER_NAME = ' '
      CHARACTER (1) :: MIN(4) = (/'.','M','I','N'/)
      CHARACTER (1) :: MIN_1(4), MOL_SPGRP(2,20) 
      CHARACTER (1) :: NEW(3) = (/'n','e','w'/), OLD(3) = (/'o','l','d'/)
      CHARACTER (1) :: OUTPUT_MIN(20), OUTPUT_XYZ(20), POT_FILE, R_KIND, cross  ! 6/24/09
      CHARACTER (1) :: SAVE(5) = (/'.','s','a','v','e'/)
      CHARACTER (1) :: SUM(6) = (/'.','S','U','M','R','Y'/)
      CHARACTER (1) :: TAB(4) = (/'.','t','a','b'/) 
      CHARACTER (1) :: WHAT_LINE2(40), XYZ(4) = (/'.','X','Y','Z'/), Y_N, YN
      CHARACTER (2) :: ALL_SPGR(N_SPGR) =  (/'AB', 'CA', 'AH', 'AF', 'AI',&
                           'AK', 'AM', 'FA', 'FC', 'DA', 'DB', 'AQ', 'AZ', 'AY',& 
     &                     'BH', 'AV', 'BD', 'BF', 'CC', 'CB', 'DC', 'DD', 'DE',&
     &                     'AA', 'AP', 'BA', 'BB', 'CD', 'CE', 'AU', 'AS',&       ! default 31
     &                     'AE', 'AC', 'AD', 'AG', 'AJ', 'AL', 'FD', 'AN', 'AO',&
     &                     'FB', 'AR', 'AT', 'BE', 'AW', 'BG', 'AX', 'BI', 'BC',&
     &                     'BJ', 'BK', 'CF', 'DF', 'DG'/)                         ! additional 23
      CHARACTER (10) :: FORT15 = '   fort.15' 
      CHARACTER (11) :: MASTERCOM2 = '_master.com' 
      CHARACTER (14) :: SEARCH_FILE
      CHARACTER (15) :: MOLPAK_FILE  
      CHARACTER (37) :: G92FILE2, NEWFILE2
      CHARACTER (40) :: ATOM1_LIST, INPUT_LINE, LINE_OUT2, MOLPAK_NAME 
      CHARACTER (40) :: WHAT_LINE, cross_terms(5) 
      CHARACTER (41) :: MASTER_NAME2
      CHARACTER (80) :: LINE_OUT
      CHARACTER (90) :: SPGLINE
      
      EQUIVALENCE (G92FILE, G92FILE2), (NEWFILE, NEWFILE2)
      EQUIVALENCE (MASTER_NAME, MASTER_NAME2), (MASTERCOM, MASTERCOM2)
      EQUIVALENCE (ATOM_LIST, ATOM1_LIST)
      EQUIVALENCE (INPUT_LINE, INPUT_LINE2) 
      EQUIVALENCE (WHAT_LINE, WHAT_LINE2)
      EQUIVALENCE (LINE_OUT, LINE_OUT1), (LINE_OUT1(41), LINE_OUT2)
      
      INTEGER :: I, IALL, IAX, IDATOM(100), IMR, INVER2, J, L, LGRPS, MID 
      INTEGER :: N, NDIF, N_GO, NMOL, NSCHARS, NSOL, N_SPGRPS, N_VOLS_SAVED 
      INTEGER :: IRB, NRB(5)
      INTEGER :: I_POTE, i_cross                      ! 6/24/09 
      
      REAL :: ATOM_TY(100,2), ATOM(100,4), G92CHARGE(100) 
      REAL :: X(100), Y(100), Z(100), ZMINA(3), ZMAXA(3), ZSTEP(3) 
 
  end module make_filesCommonMod 
