!     ************** MOLPAK program ***************                      TEXT
!                    --------------                                      TEXT
!     ****** The latest version 2009 ******                              2009
!     Feb. 9, 2001 added sulfur                                          DU-01
!     File MOLPAK.DUN                                                    TEXT
!                                                                        TEXT
!     modified on 13 March 1995 by Du, Seeing mark " DU-95"              DU-95
!     add the parts modified for using NEC (new energy coefficients) &   DU-95
!     atom with MOPAC MNDO/ESP or Gaussian 92 6-31g*/CHELPG charges      DU-95
!                                                                        TEXT
      PROGRAM MOLPAK                                                     !MLPK

      USE molpakCommonMod

      IMPLICIT NONE
!
!      CHARACTER*80 BUF                                                   !CCC
!      CHARACTER*60 HEAD                                                  !CCC
!      CHARACTER*4 AN,CDIJ,CHT1                                           !CCC
!      CHARACTER*2 AT,CODE                                                !CCC
! -- BUF, HEAD, AN, CDIJ, CHT1, AT and CODE are already decleared in common
!
!*****Transferred to module molpakCommonMod********************************** 
!-----General program communication and transfer variables               !CCC
!     Only 27 JTR and 83 ITR locations are currently in use.             !JH0702
!     COMMON BUF,KILL,NPAGE,LINE,NEW,JTR(30),ITR(90)                     !JH-95
!-----Energy constants                                                   !CCC
!     COMMON NSEP,CN(1200),ERM,FERM,ER,ERMT,ERMF,NV,NE,NNE               !CCC 
!-----Compound input data                                                !CCC
!      COMMON NMOD,NCTR,IA(200),IAA(200),X(3,200),AT(200),AN(200)        !CCC
!-----Compound rotation and structure parameters                         !CCC
!      COMMON CS2(37),SN2(37),CS3(19),SN3(19),V(3,100),W(3,100)          !CCC
!      COMMON NLOW,MLOW,NP2,KCT                                          !CCC
!      COMMON H1M,H2M,H3M,H1N,H2N,H3N,A1,A2,A3,V3                        !CCC
!-----Output flags and storage                                           !CCC
!      COMMON NQ,NR,NS,NWM,NSGT,IPR,HEAD,DIJ(500),HT1(14)                !CCC
!      DIMENSION IJD(500),CDIJ(500),CHT1(14)                             !CCC
!      EQUIVALENCE (DIJ(1),IJD(1)),(DIJ(1),CDIJ(1)),(HT1(1),CHT1(1))     !CCC
!     COMMON ICOUNT,MCOUNT,IORDER(500)                                   !CCC
!      COMMON ICOUNT,MCOUNT,IORDER(5000)                                 !DU-02
!     COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500)   !CCC
!      COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000),&
!    &       CODE(5000)
!-----Linear storage array and markers                                   !CCC
!      COMMON MAXT,MINT,LIMIT                                            !CCC
!      COMMON MARK,NEND,NSTP,NARK                                        !CCC
!      COMMON MARK1,NEND1,NSTP1                                          !CCC
!      COMMON MARK2,NEND2,NSTP2                                          !CCC
!      COMMON MARK3,NEND3,NSTP3                                          !CCC
!      COMMON MARK4,NEND4,NSTP4                                          !CCC
!      COMMON TI(100000)                                                 !DU-95
!      COMMON /NEW/IDCHG,IDATOM(200),G92CHARGE(200),ESPCHARGE(200)       !DU-95
!      DIMENSION IT(100000)                                              !DU-95
!      EQUIVALENCE (TI(1),IT(1))                                         !CCC
!      COMMON /ATMWT/ATWT(10),JATYPE(10)                                 !JH-01
! *****    transfer ended ***********
!
!      CHARACTER*4 FLST,BLANK,CARD,DD                                    !MLPK
!      CHARACTER*2 SYM,SYMP,SYMA,SYMC,CSYM,TDA,TDAA,AT1,AT2,FLAG,EXPD    !MLPK
!      CHARACTER*1 CC                                                    !MLPK
!      DIMENSION EC1(400),EC2(400),TDA(10),C(25),SYM(83),SYMP(12)        !JH0702
!     &,SYMA(36),SYMC(25),JTP(27),CSYM(25)                               !JH0702
!      DIMENSION FLST(19),CC(4),ORIG(3)                                  !JH-95
!      DIMENSION JORDER(500),ANGL1(500),ANGL2(500),ANGL3(500)            !JH0702
!      DIMENSION TDAA(20),JSEP(200),ATWTT(20),IATYPE(20)                 !JH-01
!*****    transfer ended point ***********
!
!
      CHARACTER(2) :: BLANK='  '
      CHARACTER(1) :: CC(4)
!     CHARACTER(2) :: CSYM(25)                                           ! 6-9-07
      CHARACTER(4) :: CARD, DD 
      CHARACTER(2) :: TDA(10),AT1,AT2,FLAG,EXPD
      CHARACTER(4) :: FLST(21) = (/'ATOM','HEAD','SEPC','ENGY',     &
     & 'INCL','SEEK','FINI','LIST','VOLS','FILE','INPT','WMIN',     &
     & 'CENT','AXIS','PLAN','KILL','FIND','NSEG','FNDS','PMIN',     & ! 6/24/09
     & 'CROS'/)                                                       ! 6/24/09 
      CHARACTER(2) :: SYM(83) =(/'AA','AB','AC','AD','AH','AE',     &
     &'AF','AG','CA','AI','AJ','AK','AL','AM','AN','AO','AP','AQ',  &
     &'AR','AS','AT','AU','AV','AW','AX','AY','AZ','BA','BB','BC',  &
     &'BD','BE','BF','BG','BH','BI','BJ','BK','CB','CC','CD','CE',  &
     &'CF','DA','DB','DC','DD','DE','DF','DG','EA','EB','EC','ED',  &
     &'EE','EF','EG','SA','SB','SC','SD','SE','SF','SG','SH','SI',  &
     &'SJ','SK','SL','SM','SN','FA','FB','FC','FD','SO','SP','SQ',  &
     &'SR','SS','ST','SU','SV'/) 
      CHARACTER(2) :: SYMP(12) =(/'AC','AN','AO','BJ','BK','CF',    &
     &'FB','SA','SL','SM','SO','SP'/) 
      CHARACTER(2) :: SYMA(36) =(/'AE','AJ','AL','AO','AP','AR',    &
     &'AT','AW','AX','BA','BB','BC','BE','BG','BI','BK','CD','CE',  &
     &'DC','DD','DE','DF','DG','FD','SA','SF','SG','SH','SI','SJ',  &
     &'SK','SQ','SR','SS','ST','SU'/)
      CHARACTER(2) :: SYMC(25)=(/'AB','CA','AI','AJ','AK','AL','AM',&
     &'AN','CB','CC','DC','DD','DE','FA','FB','FC','FD','SB','SC',  &
     &'SD','SE','SN','CD','CE','SV'/)
      CHARACTER(2) :: TDAA(20)=(/'C ','H ','HB','N ','O ','OX',     &
     &'F ','CL','BR','B ','S ','I ','NB','  ','  ','  ','  ','  ',  &
     &'  ','  '/)
!
      INTEGER :: I, II, IJ, INCL, IA1, IA2, IA3, ISEP, IEGY
      INTEGER :: J, JJ, JI, JII, JNCL, JCYC, JCOUNT  
      INTEGER :: JORDER(500), JSEP(200)
      INTEGER :: K, cross_count=0                            ! 6/24/09
      INTEGER :: L
      INTEGER :: M, MSEP
      INTEGER :: N, NFNDS, NINCL, NCOUNT
      INTEGER :: NP, NPS, NQS, NT 
      INTEGER, PARAMETER :: NFLST = 21                            ! 6/24/09 
      INTEGER, PARAMETER :: JTP(27) = (/1,6,9,10,14,15,17,19,27,39,& 
     &41,44,49,51,58,61,63,65,68,69,71,72,74,76,78,81,83/)
      INTEGER, PARAMETER :: IATYPE(20) = (/1,4,4,2,3,3,5,6,7,8,&
     &9,10,11,12,13,14,15,16,17,18/)
!
      REAL    :: ANGL1(500), ANGL2(500), ANGL3(500)   
      REAL    :: C(25), C1, C2, C3  
      REAL    :: D, DDD  
      REAL    :: E, E1, E2, EC1(400), EC2(400)  
      REAL    :: F  
      REAL    :: ORIG(3)  
      REAL    :: STEP1, STEP2, STEP3 
      REAL    :: ATWTT(20) = (/12.0111,1.0080,1.0080,14.0067,&
     &15.9994,15.9994,18.9984,35.4530,79.9040,10.811,32.066,126.904,14.0067,&
     &0.0,0.0,0.0,0.0,0.0,0.0,0.0/)    
!
!      EQUIVALENCE (DIJ(1),IJD(1)),(DIJ(1),CDIJ(1)),(HT1(1),CHT1(1))      !CCC
! 
!     The JTP array consists of sequence numbers in the ITR array where  TEXT 
!      the referenced space group (subgroup) is the first one treated by TEXT 
!      the next subroutine - ZTRTA, then ZTRTB, then ZTRTC, etc.         TEXT 
!-----The FLST array contains the names of program instruction lines     TEXT 
!      DATA FLST/'ATOM','HEAD','SEPC','ENGY','INCL','SEEK','FINI',&       JH-95 
!     &'LIST','VOLS','FILE','INPT','WMIN','CENT','AXIS','PLAN','KILL',&   MLPK    
!     &'FIND','NSEG','FNDS'/                                              MLPK    
!-----Space group (subgroup) codes                                       TEXT 
!      DATA SYM/'AA','AB','AC','AD','AH','AE','AF','AG','CA','AI','AJ',&  MLPK 
!     &'AK','AL','AM','AN','AO','AP','AQ','AR','AS','AT','AU','AV','AW',& MLPK 
!     &'AX','AY','AZ','BA','BB','BC','BD','BE','BF','BG','BH','BI','BJ',& MLPK 
!     &'BK','CB','CC','CD','CE','CF','DA','DB','DC','DD','DE','DF','DG',& MLPK 
!     &'EA','EB','EC','ED','EE','EF','EG','SA','SB','SC','SD','SE','SF',& MLPK 
!     &'SG','SH','SI','SJ','SK','SL','SM','SN','FA','FB','FC','FD','SO',& JH-95 
!     &'SP','SQ','SR','SS','ST','SU','SV'/                                JH0702 
!     Possible space group (subgroup) codes for molecules on mirrors     TEXT 
!      DATA SYMP/'AC','AN','AO','BJ','BK','CF','FB','SA','SL','SM','SO',& JH-96 
!     &'SP'/                                                              JH-96 
!-----Possible space group (subgroup) codes for molecules on axes         TEXT 
!      DATA SYMA/'AE','AJ','AL','AO','AP','AR','AT','AW','AX','BA','BB',& !MLPK 
!     &'BC','BE','BG','BI','BK','CD','CE','DC','DD','DE','DF','DG','FD',& !JH-95 
!     &'SA','SF','SG','SH','SI','SJ','SK','SQ','SR','SS','ST','SU'/       !JH-96 
!-----Space group (subgroup) codes programmed for centric molecules       !TEXT 
!      DATA SYMC/'AB','CA','AI','AJ','AK','AL','AM','AN','CB','CC','DC',& !MLPK 
!     &'DD','DE','FA','FB','FC','FD','SB','SC','SD','SE','SN','CD','CE',& !JH0702 
!     &'SV'/                                                              !JH0702 
!-----Default atom type symbols                                           !TEXT 
!      DATA TDAA/'C ','H ','HB','N ','O ','OX','F ','CL','BR','B ',&      !DU-99 
!     &          'S ','  ','  ','  ','  ','  ','  ','  ','  ','  '/       !DU-01 
!-----Table to convert program atom types to WMIN atom types              TEXT 
!      DATA IATYPE/1,4,4,2,3,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18/      !JH-01  
!-----Atomic weights                                                      TEXT 
!      DATA ATWTT/12.0111,1.0080,1.0080,14.0067,15.9994,15.9994,18.9984  &!DU-01 
!     &,35.4530,79.9040,10.811,32.066,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,   &!JH-01 
!     &0.0/                                                               !JH-01  
      DATA EC1/                                                          &!DU-04 
     &-.3272E+03,-.1250E+03,-.1250E+03,-.3400E+03,-.3423E+03,-.3423E+03, &
     &-.3423E+03,-.6840E+03,-.6720E+03,-.3272E+03,-.1125E+04,-.6720E+03, &
     &-.3403E+03,         0,         0,         0,         0,         0, &
     &         0,         0,-.1250E+03,-.4920E+02,-.4920E+02,-.1320E+03, &
     &-.1327E+03,-.1327E+03,-.1327E+03,-.2652E+03,-.2606E+03,-.1250E+03, &
     &-.4364E+03,-.2606E+03,-.1320E+03,         0,         0,         0, &
     &         0,         0,         0,         0,-.1250E+03,-.4920E+02, &
     &-.4920E+02,-.1320E+03,-.1327E+03,-.1327E+03,-.1327E+03,-.2652E+03, &
     &-.2606E+03,-.1250E+03,-.4364E+03,-.2606E+03,-.1320E+03,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &-.3400E+03,-.1320E+03,-.1320E+03,-.3540E+03,-.3560E+03,-.3560E+03, &
     &-.3560E+03,-.7115E+03,-.7115E+03,-.3400E+03,-.1171E+04,-.6989E+03, &
     &-.3540E+03,         0,         0,         0,         0,         0, &
     &         0,         0,-.3423E+03,-.1327E+03,-.1327E+03,-.3560E+03, &
     &-.3585E+03,-.3585E+03,-.3585E+03,-.7155E+03,-.7155E+03,-.3423E+03, &
     &-.1178E+04,-.7034E+03,-.3562E+03,         0,         0,         0, &
     &         0,         0,         0,         0,-.3423E+03,-.1327E+03, &
     &-.1327E+03,-.3560E+03,-.3585E+03,-.3585E+03,-.3585E+03,-.7155E+03, &
     &-.7155E+03,-.3423E+03,-.1178E+04,-.7034E+03,-.3562E+03,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &-.3423E+03,-.1327E+03,-.1327E+03,-.3560E+03,-.3585E+03,-.3585E+03, &
     &-.3585E+03,-.7155E+03,-.7155E+03,-.3423E+03,-.1178E+04,-.7034E+03, &
     &-.3562E+03,         0,         0,         0,         0,         0, &
     &         0,         0,-.6840E+03,-.2652E+03,-.2652E+03,-.7115E+03, &
     &-.7155E+03,-.7155E+03,-.7155E+03,-.1430E+04,-.1430E+04,-.6840E+03, &
     &-.2353E+04,-.1405E+04,-.7115E+03,         0,         0,         0, &
     &         0,         0,         0,         0,-.6720E+03,-.2606E+03, &
     &-.2606E+03,-.7115E+03,-.7155E+03,-.7155E+03,-.7155E+03,-.1430E+04, &
     &-.1380E+04,-.6720E+03,-.2311E+04,-.1380E+04,-.6989E+03,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &-.3272E+03,-.1250E+03,-.1250E+03,-.3400E+03,-.3423E+03,-.3423E+03, &
     &-.3423E+03,-.6840E+03,-.6720E+03,-.3272E+03,-.1125E+04,-.6720E+03, &
     &-.3403E+03,         0,         0,         0,         0,         0, &
     &         0,         0,-.1125E+04,-.4364E+03,-.4364E+03,-.1171E+04, &
     &-.1178E+04,-.1178E+04,-.1178E+04,-.2353E+04,-.2311E+04,-.1125E+04, &
     &-.3870E+04,-.2311E+04,-.1170E+04,         0,         0,         0, &
     &         0,         0,         0,         0,-.6720E+03,-.2606E+03, &
     &-.2606E+03,-.6989E+03,-.7034E+03,-.7034E+03,-.7034E+03,-.1405E+04, &
     &-.1380E+04,-.6720E+03,-.2311E+04,-.1380E+04,-.6989E+03,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &-.3403E+03,-.1320E+03,-.1320E+03,-.3540E+03,-.3562E+03,-.3562E+03, &
     &-.3562E+03,-.7115E+03,-.6989E+03,-.3403E+03,-.1170E+04,-.6989E+03, &
     &-.3540E+03,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0/
      DATA EC2/                                                          &!DU-04 
     &0.4500E+06,0.7747E+05,0.7747E+05,0.3400E+06,0.4234E+06,0.2790E+06, &
     &0.4234E+06,0.7370E+06,0.1037E+07,0.4500E+06,0.1873E+07,0.1174E+07, &
     &0.4173E+06,         0,         0,         0,         0,         0, &
     &         0,         0,0.7747E+05,0.1737E+05,0.1737E+05,0.9010E+05, &
     &0.2092E+05,0.2092E+05,0.2092E+05,0.1393E+06,0.4693E+06,0.7747E+05, &
     &0.3681E+06,0.2307E+06,0.8199E+05,         0,         0,         0, &
     &         0,         0,         0,         0,0.7747E+05,0.1737E+05, &
     &0.1737E+05,0.1491E+04,0.1491E+04,0.1491E+04,0.2092E+05,0.1393E+06, &
     &0.4693E+06,0.7747E+05,0.3681E+06,0.2307E+06,0.8199E+05,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &0.3400E+06,0.9010E+05,0.1491E+04,0.3870E+06,0.2622E+06,0.3160E+06, &
     &0.2622E+06,0.8328E+06,0.8328E+06,0.3400E+06,0.1737E+07,0.1089E+07, &
     &0.3870E+06,         0,         0,         0,         0,         0, &
     &         0,         0,0.4234E+06,0.2092E+05,0.1491E+04,0.2622E+06, &
     &0.3656E+06,0.3656E+06,0.3656E+06,0.6747E+06,0.6747E+06,0.4234E+06, &
     &0.1689E+07,0.1058E+07,0.3761E+06,         0,         0,         0, &
     &         0,         0,         0,         0,0.2790E+06,0.2092E+05, &
     &0.1491E+04,0.3160E+06,0.3656E+06,0.3656E+06,0.3656E+06,0.6747E+06, &
     &0.6747E+06,0.2790E+06,0.1689E+07,0.1058E+07,0.3761E+06,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &0.4234E+06,0.2092E+05,0.2092E+05,0.2622E+06,0.3656E+06,0.3656E+06, &
     &0.3656E+06,0.6747E+06,0.6747E+06,0.4234E+06,0.1689E+07,0.1058E+07, &
     &0.3761E+06,         0,         0,         0,         0,         0, &
     &         0,         0,0.7370E+06,0.1393E+06,0.1393E+06,0.8328E+06, &
     &0.6747E+06,0.6747E+06,0.6747E+06,0.2063E+07,0.2063E+07,0.7370E+06, &
     &0.4011E+07,0.2514E+07,0.8935E+06,         0,         0,         0, &
     &         0,         0,         0,         0,0.1037E+07,0.4693E+06, &
     &0.4693E+06,0.8328E+06,0.6747E+06,0.6747E+06,0.6747E+06,0.2063E+07, &
     &0.3063E+07,0.1037E+07,0.4888E+07,0.3063E+07,0.1089E+07,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &0.4500E+06,0.7747E+05,0.7747E+05,0.3400E+06,0.4234E+06,0.2790E+06, &
     &0.4234E+06,0.7370E+06,0.1037E+07,0.4500E+06,0.1873E+07,0.1174E+07, &
     &0.4173E+06,         0,         0,         0,         0,         0, &
     &         0,         0,0.1873E+07,0.3681E+06,0.3681E+06,0.1737E+07, &
     &0.1689E+07,0.1689E+07,0.1689E+07,0.4011E+07,0.4888E+07,0.1873E+07, &
     &0.7799E+07,0.4888E+07,0.1737E+07,         0,         0,         0, &
     &         0,         0,         0,         0,0.1174E+07,0.2307E+06, &
     &0.2307E+06,0.1089E+07,0.1058E+07,0.1058E+07,0.1058E+07,0.2514E+07, &
     &0.3063E+07,0.1174E+07,0.4888E+07,0.3063E+07,0.1089E+07,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &0.4173E+06,0.8199E+05,0.8199E+05,0.3870E+06,0.3761E+06,0.3761E+06, &
     &0.3761E+06,0.8935E+06,0.1089E+07,0.4173E+06,0.1737E+07,0.1089E+07, &
     &0.3870E+06,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0,         0,         0, &
     &         0,         0,         0,         0/
!-----NPMIN = 0  write out wmin.input NPMIN = 1 write out PMIN.inp       ! 1-29-09
      NPMIN = 0
!-----Set standard input and output file numbers                         TEXT 
!     NPS=5                  ! 5 for VMS, 15 for UNIX                    !MLPK 
!     NQS=6                  ! 6 for VMS, 16 for UNIX                    !MLPK 
       NPS=15                 ! 5 for VMS, 15 for UNIX                   !DU 
       NQS=16                 ! 6 for VMS, 16 for UNIX                   !DU 
      NP=NPS                                                             !MLPK  
      NQ=NQS                                                             !MLPK 
!-----Set flag that FNDS instruction has NOT been read.                  TEXT 
      NFNDS=0                                                            !RTPK 
!-----Set flag that INCL instruction has NOT been read.                  !JH-95 
      NINCL=0                                                            !JH-95 
!-----Set WMIN (NR) and MOLCON (NS) input file numbers to default        MLPK 
      NR=9                                                               !MLPK 
      NS=NQ                                                              !MLPK 
      NPAGE=0                                                            !MLPK    
      LINE=10                                                            !MLPK    
      KILL=0                                                             !MLPK    
      IPR=0                                                              !MLPK    
!-----Set standard values of number of segments into which lengths are   TEXT 
!      divided during the three stages of minimization searches          TEXT 
      NV=8                                                               !MLPK 
      NE=128                                                             !MLPK 
      NNE=1024                                                           !MLPK 
!-----Set default step sizes for FIND routine - modified by SEEK.        TEXT 
      DIJ(498) = 10.0                                                    !JH0502 
      DIJ(499) = 10.0                                                    !JH0502 
      DIJ(500) = 10.0                                                    !JH0502 
!-----Set default values for saving minimum volume structures            TEXT 
      MCOUNT=10                                                          !MLPK 
      ICOUNT=1                                                           !MLPK 
      NT=0                                                               !MLPK 
      IORDER(1)=1                                                        !MLPK 
      VOLUME(1)=9999.                                                    !MLPK 
      ANGLE1(1)=999.                                                     !MLPK 
      ANGLE2(1)=999.                                                     !MLPK 
      ANGLE3(1)=999.                                                     !MLPK 
      CODE(1)=BLANK                                                      !MLPK 
      NSGT=0                                                             !MLPK 
      NMOD=0                                                             !MLPK    
      NCTR=0                                                             !MLPK    
      NWM=0                                                              !MLPK    
!     Number of different atom types in stored tables                    DU-95 
      MSEP=13                                                            !DU-04 
!-----Number of atom types which have been read in ATOM instructions 
      NSEP=0  
      ISEP=0                                                             !MLPK    
      IEGY=1                                                             !MLPK    
      MINT=1                                                             !MLPK    
      MAXT=1000000                                                       !DU-03-15-06 
      MARK1=MINT                                                         !MLPK    
      LIMIT=MAXT-10                                                      !MLPK 
      ERM=.5                                                             !MLPK    
      FERM=.1                                                            !MLPK    
      OPEN(12, FILE='fort.12')                                           !MLPK 
      OPEN(17, FILE='fort.17')                                           !MLPK 
      OPEN(21, FILE='fort.21')                                           !MLPK 
      WRITE(NQ,5)                                                        !MLPK    
    5 FORMAT(4X,4(' **Begin MOLPAK** '))                                 !MLPK 
      HEAD=BLANK                                                         !MLPK    
   50 READ(NP,51,END=52) BUF                                             !MLPK 
   51 FORMAT(A80)                                                        !MLPK    
      GO TO 53                                                           !MLPK    
   52 IF(NP.EQ.NPS)GO TO 1000                                            !MLPK    
      NP=NPS                                                             !MLPK    
      GO TO 50                                                           !MLPK    
   53 IF(IPR.LE.0)GO TO 75                                               !MLPK    
      WRITE(NQ,67) BUF                                                   !MLPK 
   67 FORMAT(1H ,3X,1H*,A80)                                             !MLPK 
   75 READ (BUF,76) CARD                                                 !MLPK 
   76 FORMAT(A4)                                                         !MLPK    
      DO 77 I=1,NFLST                                                    !MLPK    
      IF(CARD.NE.FLST(I))GO TO 77                                        !MLPK    
      if (CARD == 'PMIN' .OR. I == 20) NPMIN = 1                         ! 1-29-09 DU
      GO TO(80,135,140,150,155,330,110,100,120,125,130,115,400,  &       !JH-96 
     &410,420,105,332,335,490,115,305),I                            ! 6/24/09    
   77 CONTINUE                                                           !MLPK    
      WRITE(NQ,78)CARD                                                   !MLPK    
      WRITE(17,78)CARD                                                   !MLPK    
   78 FORMAT(5H0***(,A4,') INSTRUCTION LINE FOUND')                      !MLPK 
      KILL=-1                                                            !MLPK    
      GO TO 990                                                          !MLPK    
!-----ATOM instruction                                                   TEXT 
   80 NMOD=NMOD+1                                                        !MLPK    
!     Read atom coordinates and charges                                  !DU-95 
      READ (BUF,81) AT(NMOD),AN(NMOD),(X(I,NMOD),I=1,3),IDATOM(NMOD),& 
     &              G92CHARGE(NMOD) 
   81    FORMAT (5X,A2,A4,3F10.6,I5,F10.6)                   
      DO 82 I=1,MSEP                                                     !MLPK    
      IF(AT(NMOD).EQ.TDAA(I))GO TO 84                                    !MLPK    
   82 CONTINUE                                                           !MLPK    
      WRITE(NQ,83)AT(NMOD),AN(NMOD)                                      !MLPK    
      WRITE(17,83)AT(NMOD),AN(NMOD)                                      !MLPK    
   83 FORMAT(6H0*** (,A2,1H),A4,' - ( ) NOT RECOGNIZED ATOM TYPE')       !MLPK    
      KILL=1                                                             !MLPK    
      NMOD=NMOD-1                                                        !MLPK    
      GO TO 50                                                           !MLPK    
   84 IF(NMOD.EQ.1)GO TO 86  
      DO 85 J=1,NSEP  
	IF(AT(NMOD).EQ.TDA(J))GO TO 90  
   85 CONTINUE   
   86 NSEP=NSEP+1  
      IF(NSEP.LE.10)GO TO 89                                              
	WRITE(NQ,87)AT(NMOD),AN(NMOD)  
	WRITE(17,87)AT(NMOD),AN(NMOD)  
   87 FORMAT(' *** (',A2,A4,') MAKES TOO MANY ATOM TYPES IN USE')  
      KILL=1  
	NMOD=NMOD-1  
	GO TO 50  
   89 TDA(NSEP)=TDAA(I)  
      ATWT(NSEP)=ATWTT(I)                                                !JH-01 
      JATYPE(NSEP)=IATYPE(I)                                             !JH-01 
      JSEP(NSEP)=I  
      J=NSEP  
   90 IA(NMOD)=12*(J-1)+1                                                !MLPK    
      IAA(NMOD)=120*(J-1)                                                !MLPK    
      WRITE(NQ,91)NMOD,AT(NMOD),AN(NMOD),(X(I,NMOD),I=1,3),IA(NMOD),  &  !MLPK    
     &IAA(NMOD)                                                          !MLPK    
   91 FORMAT(10X,4HATOM,I4,5X,A2,A4,3F10.6,5X,2I5)                       !DU-98  
      GO TO 50                                                           !MLPK    
!-----LIST instruction                                                   !TEXT 
  100 READ (BUF,101) IPR                                                 !TEXT 
  101 FORMAT(4X,I2)                                                      !MLPK    
      GO TO 50                                                           !MLPK    
!-----KILL instruction                                                   TEXT 
  105 KILL=0                                                             !MLPK    
      GO TO 50                                                           !MLPK    
!-----FINI instruction                                                   TEXT 
110   call wmin                                                         ! 6/26/09
         GO TO 500                                                      ! 6/26/09 
!-----WMIN or PMIN instruction; potential_file = name of file with  potential 
  115 READ (BUF,116) DD,NWM,NR,NS,POTENTIAL_FILe, n_cross                ! 6/24/09      
  116 FORMAT(4X,A2,3I3,a40,i4)                                           ! 6/24/09 
      NSGT=0                                                             !MLPK 
      IF(DD.NE.BLANK)NSGT=1                                              !MLPK 
      IF(NR.LE.0)NR=9                                                    !MLPK 
      IF(NS.LE.0)NS=NQ                                                   !MLPK 
      OPEN(8, FILE='fort.8')                                             !MLPK 
      OPEN(9, FILE='fort.9')                                             !MLPK 
!     CALL WMIN                                                    ! 6/26/09 
      GO TO 50                                                     ! 6/26/09
!----CROS line...info on PMIN cross-term coefs                        6/24/09
305   cross_count = cross_count + 1                                 ! 6/24/09
      if (cross_count .gt. 5) go to 50   ! no more than 5           ! 6/24/09
      read (buf,'(10x,a)') cross_terms(cross_count)                 ! 6/24/09
!     print 5656, n_cross, cross_count, cross_terms(cross_count)    ! temp ###########
5656     format ('**from molpak..n_cross, cross_count,cross_terms =',2i3,'|',a40) ! temp ###########
      go to 50                                                      ! 6/24/09
!-----VOLS instruction                                                   TEXT 
  120 READ (BUF,121)MCOUNT,NT                                            !MLPK 
! 121 FORMAT(5X,2I3)                                                     !MLPK    
  121 FORMAT(5X,2I5)                                                     !DU-07 
!     MCOUNT=MIN0(MCOUNT,500)                                            !MLPK 
      MCOUNT=MIN0(MCOUNT,60000)                                          !DU-07 
      IF(NT.EQ.8)OPEN(8, FILE='fort.8')                                  !RTPK 
      GO TO 50                                                           !MLPK    
!-----FILE instruction                                                   TEXT 
  125 READ (BUF,131) NQ                                                  !MLPK 
      IF(NQ.LE.0)NQ=NQS                                                  !MLPK 
      IF(NQ.EQ.22)OPEN(22, FILE='fort.22')                               !MLPK 
      IF(NQ.EQ.23)OPEN(23, FILE='fort.23')                               !MLPK 
      GO TO 50                                                           !MLPK    
!-----INPT instruction                                                   TEXT 
  130 READ (BUF,131) NP                                                  !MLPK 
  131 FORMAT(4X,I3)                                                      !MLPK 
      IF(NP.LE.0)NP=NPS                                                  !MLPK 
      IF(NP.EQ.31)OPEN(31, FILE='fort.31')                               !MLPK 
      IF(NP.EQ.32)OPEN(32, FILE='fort.32')                               !MLPK 
      IF(NP.EQ.33)OPEN(33, FILE='fort.33')                               !MLPK 
      IF(NP.EQ.34)OPEN(34, FILE='fort.34')                               !MLPK 
      IF(NP.EQ.35)OPEN(35, FILE='fort.35')                               !MLPK 
      IF(NP.EQ.36)OPEN(36, FILE='fort.36')                               !MLPK 
      IF(NP.EQ.37)OPEN(37, FILE='fort.37')                               !MLPK 
      IF(NP.EQ.38)OPEN(38, FILE='fort.38')                               !MLPK 
      IF(NP.EQ.39)OPEN(39, FILE='fort.39')                               !MLPK 
      GO TO 50                                                           !MLPK    
!-----HEAD instruction                                                   TEXT 
  135 READ (BUF,136) HEAD                                                !MLPK 
  136 FORMAT(4X,A60)                                                     !MLPK    
      GO TO 50                                                           !MLPK    
!-----SEPC instruction                                                   TEXT 
  140 READ (BUF,141) AT1,AT2,C1,C2,C3                                    !JH-01 
  141 FORMAT(5X,2A2,1X,2E10.0,F10.0)                                     !JH-01    
      IF(AT1.EQ.BLANK)AT1=AT2                                            !MLPK    
      IF(AT2.EQ.BLANK)AT2=AT1                                            !MLPK    
      IF(MSEP.LE.0)GO TO 143                                             !MLPK    
      DO 142 I=1,MSEP                                                    !MLPK    
      IF(AT1.EQ.TDAA(I))GO TO 144                                        !MLPK    
  142 CONTINUE                                                           !MLPK    
  143 IF(MSEP.GE.20)GO TO 146     
      MSEP=MSEP+1                                                        !MLPK    
      I=MSEP                                                             !MLPK    
      TDAA(I)=AT1                                                        !MLPK    
      ATWTT(I)=C3                                                        !JH-01 
  144 DO 145 J=1,MSEP                                                    !MLPK    
      IF(AT2.EQ.TDAA(J))GO TO 149                                        !MLPK    
  145 CONTINUE                                                           !MLPK    
      IF(MSEP.LT.20)GO TO 148      
  146 WRITE(NQ,147)AT1,AT2                                               !MLPK    
      WRITE(17,147)AT1,AT2                                               !MLPK    
  147 FORMAT(4H0***,2A2,' MAKES TOO MANY ATOM TYPES')                    !MLPK    
      KILL=1                                                             !MLPK    
      GO TO 50                                                           !MLPK    
  148 MSEP=MSEP+1                                                        !MLPK    
      J=MSEP                                                             !MLPK    
      TDAA(J)=AT2                                                        !MLPK    
      ATWT(J)=C3                                                         !JH-01 
  149 IJ=20*(I-1)+J                                                      !MLPK    
      JI=20*(J-1)+I                                                      !MLPK    
      EC1(IJ)=C1                                                         !MLPK    
      EC2(IJ)=C2                                                         !MLPK    
      EC1(JI)=C1                                                         !MLPK    
      EC2(JI)=C2                                                         !MLPK    
      ISEP=1                                                             !MLPK    
      IEGY=1                                                             !MLPK    
      GO TO 50                                                           !MLPK    
!-----ENGY instruction                                                   TEXT 
  150 READ (BUF,151) ERM,FERM                                            !MLPK 
  151 FORMAT(4X,2F4.0)                                                   !MLPK    
      IEGY=1                                                             !MLPK    
      GO TO 50                                                           !MLPK    
!-----INCL instruction                                                   TEXT 
  155 READ (BUF,156) CSYM                                                !JH-96 
  156 FORMAT(5X,25(A2,1X))                                               !MLPK    
!-----Set flag that INCL instruction has been read.                      !JH-95 
      NINCL=1                                                            !JH-95 
      DO 157 I=1,83                                                      !JH0702 
  157 ITR(I)=0                                                           !JH-96 
!-----Are plane, axis, or centric molecule flags turned on?              TEXT 
      IF(NCTR.NE.0)GO TO 180                                             !MLPK 
      DO 179 I=1,25                                                      !JH-95 
      IF(CSYM(I).EQ.BLANK)GO TO 179                                      !MLPK    
      DO 158 J=1,57                                                      !JH-96 
      IF(CSYM(I).NE.SYM(J))GO TO 158                                     !MLPK    
      IF(J.EQ.43)GO TO 176                                               !JH-95 
      ITR(J)=1                                                           !MLPK    
      GO TO 179                                                          !MLPK    
  158 CONTINUE                                                           !MLPK    
      DO 159 J=58,83                                                     !JH0702 
      IF(CSYM(I).NE.SYM(J))GO TO 159                                     !JH-96 
      K=J-57                                                             !JH-96 
      GO TO(161,163,163,163,163,165,165,165,165,165,165,167,167,163  &   !JH-96 
     &,169,169,169,169,167,167,165,165,165,165,165,163),K                !JH0702 
  159 CONTINUE                                                           !JH-96 
      WRITE(NQ,160)CSYM(I)                                               !MLPK    
      WRITE(17,160)CSYM(I)                                               !MLPK    
  160 FORMAT('0BAD SYMBOL (',A2,') IN INCL INSTRUCTION LINE')            !MLPK 
      GO TO 178                                                          !JH-95 
  161 WRITE(NQ,162)CSYM(I)                                               !JH-96 
      WRITE(17,162)CSYM(I)                                               !JH-96 
  162 FORMAT('0Code ',A2,' requires molecule on two-fold axis or mirror &
     &plane')                                                            !JH-96 
      GO TO 178                                                          !JH-96 
  163 WRITE(NQ,164)CSYM(I)                                               !JH-96 
      WRITE(17,164)CSYM(I)                                               !JH-96 
  164 FORMAT('0Code ',A2,' requires molecule on center of symmetry')     !JH-96 
      GO TO 178                                                          !JH-96 
  165 WRITE(NQ,166)CSYM(I)                                               !JH-96 
      WRITE(17,166)CSYM(I)                                               !JH-96 
  166 FORMAT('0Code ',A2,' requires molecule on two-fold axis')          !JH-96 
      GO TO 178                                                          !JH-96 
  167 WRITE(NQ,168)CSYM(I)                                               !JH-96 
      WRITE(17,168)CSYM(I)                                               !JH-96 
  168 FORMAT('0Code ',A2,' requires molecule on mirror plane')           !JH-96 
      GO TO 178                                                          !JH-96 
  169 ITR(J)=1                                                           !JH-96 
      GO TO 179                                                          !JH-96 
  176 WRITE(NQ,177)                                                      !JH-95 
      WRITE(17,177)                                                      !JH-95 
  177 FORMAT('0Code CF (Pnma with no molecular symmetry) not programmed')!JH-95 
  178 KILL=1                                                             !MLPK    
  179 CONTINUE                                                           !MLPK    
      GO TO 325                                                          !MLPK    
  180 GO TO(181,181,181,200,200,200,250),NCTR                            !MLPK 
!-----Molecule contains a mirror plane                                   TEXT 
  181 DO 199 I=1,25                                                      !MLPK 
      IF(CSYM(I).EQ.BLANK)GO TO 199                                      !MLPK    
      DO 197 J=1,12                                                      !JH-96 
      IF(CSYM(I).NE.SYMP(J))GO TO 197                                    !MLPK 
!            AC  AN  AO  BJ  BK  CF  FB  SA  SL  SM  SO  SP              TEXT 
      GO TO(182,183,184,185,186,187,184,188,189,190,191,192),J           !JH-96   
  182 ITR(58)=1                                                          !MLPK 
      GO TO 199                                                          !MLPK 
  183 ITR(5)=1                                                           !MLPK 
      ITR(7)=1                                                           !MLPK 
      GO TO 199                                                          !MLPK 
  184 ITR(6)=1                                                           !MLPK 
      GO TO 199                                                          !MLPK 
  185 ITR(69)=1                                                          !MLPK 
      GO TO 199                                                          !MLPK 
  186 ITR(70)=1                                                          !MLPK 
      GO TO 199                                                          !MLPK 
  187 ITR(76)=1                                                          !JH-95 
      ITR(77)=1                                                          !JH-95 
      GO TO 199                                                          !JH-95 
  188 ITR(58)=1                                                          !JH-96 
      GO TO 199                                                          !JH-96 
  189 ITR(69)=1                                                          !JH-96 
      GO TO 199                                                          !JH-96 
  190 ITR(70)=1                                                          !JH-96 
      GO TO 199                                                          !JH-96 
  191 ITR(76)=1                                                          !JH-96 
      GO TO 199                                                          !JH-96 
  192 ITR(77)=1                                                          !JH-96 
      GO TO 199                                                          !JH-96 
  197 CONTINUE                                                           !MLPK 
      WRITE(NQ,198)CSYM(I)                                               !MLPK 
      WRITE(17,198)CSYM(I)                                               !MLPK 
  198 FORMAT('0Code ',A2,' not compatible with molecules on mirrors')    !MLPK 
      KILL=1                                                             !MLPK 
  199 CONTINUE                                                           !MLPK 
      GO TO 325                                                          !MLPK 
!-----Molecule contains a two-fold axis                                  !TEXT 
  200 DO 240 I=1,25                                                      !MLPK    
      IF(CSYM(I).EQ.BLANK)GO TO 240                                      !MLPK    
      DO 238 J=1,36                                                      !JH-96 
      IF(CSYM(I).NE.SYMA(J))GO TO 238                                    !MLPK    
!            AE  AJ  AL  AO  AP  AR  AT  AW  AX  BA  BB  BC  BE  BG  BI  TEXT 
      GO TO(201,202,202,203,204,205,205,206,207,204,204,208,205,206,207& !MLPK 
!       BK  CD  CE  DC  DD  DE  DF  DG  FD  SA  SF  SG  SH  SI  SJ  SK   !JH-96 
     &,209,220,218,217,217,217,219,219,202,221,222,223,224,225,226,227&  !JH-96 
!       SQ  SR  SS  ST  SU                                               TEXT 
     &,228,229,230,231,232),J                                            !JH-96 
  201 ITR(58)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  202 ITR(8)=1                                                           !MLPK 
      GO TO 240                                                          !MLPK 
  203 ITR(3)=1                                                           !MLPK 
      GO TO 240                                                          !MLPK 
  204 ITR(63)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  205 ITR(65)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  206 ITR(66)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  207 ITR(68)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  208 ITR(64)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  209 ITR(67)=1                                                          !MLPK 
      GO TO 240                                                          !MLPK 
  217 ITR(44)=1                                                          !MLPK    
      GO TO 240                                                          !MLPK    
  218 ITR(78)=1                                                          !JH-95 
      GO TO 240                                                          !JH-95 
  219 ITR(81)=1                                                          !JH-95 
      ITR(82)=1                                                          !JH-96 
      GO TO 240                                                          !JH-95  
  220 ITR(78)=1                                                          !JH-95 
      ITR(79)=1                                                          !JH-95 
      ITR(80)=1                                                          !JH-95 
      GO TO 240                                                          !JH-95 
  221 ITR(58)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  222 ITR(63)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  223 ITR(64)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  224 ITR(65)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  225 ITR(66)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  226 ITR(67)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  227 ITR(68)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  228 ITR(78)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  229 ITR(79)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  230 ITR(80)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  231 ITR(81)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  232 ITR(82)=1                                                          !JH-96 
      GO TO 240                                                          !JH-96 
  238 CONTINUE                                                           !MLPK    
      WRITE(NQ,239)CSYM(I)                                               !MLPK    
      WRITE(17,239)CSYM(I)                                               !MLPK    
  239 FORMAT(' Code ',A2,' not compatible with molecules on axes')       !MLPK    
      KILL=1                                                             !MLPK    
  240 CONTINUE                                                           !MLPK    
      GO TO 325                                                          !MLPK    
!-----Molecule contains a center of symmetry                             TEXT 
  250 DO 269 I=1,25                                                      !MLPK    
      IF(CSYM(I).EQ.BLANK)GO TO 269                                      !MLPK    
      DO 267 J=1,25                                                      !JH0702 
      IF(CSYM(I).NE.SYMC(J))GO TO 267                                    !MLPK    
!            AB  CA  AI  AJ  AK  AL  AM  AN  CB  CC  DC  DD  DE  FA  FB  TEXT 
      GO TO(251,251,252,253,252,253,254,255,257,257,258,258,258,254,255,&
!            FC  FD  SB  SC  SD  SE  SN  CD  CE  SV                      TEXT 
     &      254,253,259,260,261,262,263,264,264,264),J                   !JH0702 
  251 ITR(1)=1                                                           !MLPK 
      GO TO 269                                                          !MLPK 
  252 ITR(61)=1                                                          !MLPK 
      GO TO 269                                                          !MLPK    
  253 ITR(62)=1                                                          !MLPK 
      GO TO 269                                                          !MLPK 
  254 ITR(60)=1                                                          !MLPK 
      ITR(61)=1                                                          !MLPK 
      GO TO 269                                                          !MLPK 
  255 ITR(59)=1                                                          !MLPK 
      GO TO 269                                                          !MLPK 
  257 ITR(71)=1                                                          !MLPK 
      GO TO 269                                                          !MLPK    
  258 ITR(45)=1                                                          !MLPK    
      GO TO 269                                                          !MLPK    
  259 ITR(59)=1                                                          !JH-96 
      GO TO 269                                                          !JH-96 
  260 ITR(60)=1                                                          !JH-96 
      GO TO 269                                                          !JH-96 
  261 ITR(61)=1                                                          !JH-96 
      GO TO 269                                                          !JH-96 
  262 ITR(62)=1                                                          !JH-96 
      GO TO 269                                                          !JH-96 
  263 ITR(71)=1                                                          !JH-96 
      GO TO 269                                                          !JH-96 
  264 ITR(83)=1                                                          !JH0702 
      GO TO 269                                                          !JH0702 
  267 CONTINUE                                                           !MLPK    
      WRITE(NQ,268)CSYM(I)                                               !MLPK    
      WRITE(17,268)CSYM(I)                                               !MLPK    
  268 FORMAT(' Code ',A2,' not compatible with centric molecules')       !MLPK 
      KILL=1                                                             !MLPK    
  269 CONTINUE                                                           !MLPK    
!-----Set ITR and JTR arrays to transfer to codes which have been        TEXT 
!     activated.                                                         TEXT 
  325 JNCL=1                                                             !MLPK    
      INCL=1                                                             !MLPK    
!      PRINT 1325,(ITR(I),I=1,83)  
! 1325 FORMAT(' ITR =',30I2/6X,30I2)  
      K=27                                                               !JH0702 
      J=83                                                               !JH0702 
  326 IF(ITR(J).NE.1)GO TO 327                                           !MLPK    
      INCL=0                                                             !MLPK    
      JNCL=0                                                             !MLPK    
  327 INCL=INCL+1                                                        !MLPK    
      ITR(J)=INCL                                                        !MLPK    
      IF(J.NE.JTP(K))GO TO 328                                           !MLPK    
      JNCL=JNCL+1                                                        !MLPK    
      JTR(K)=JNCL                                                        !MLPK    
      INCL=1                                                             !MLPK    
      K=K-1                                                              !MLPK    
  328 J=J-1                                                              !MLPK    
      IF(J.GT.0)GO TO 326                                                !MLPK    
      IF(IPR.LT.1)GO TO 50                                               !MLPK    
      WRITE(NQ,329)(JTR(I),I=1,27),(ITR(J),J=1,83)                       !JH0702    
  329 FORMAT(' SEEK POINTERS, JTR =',27I2/5X,'ITR =',33I2/10X,33I2/&     !MLPK 
     &10X,33I2)                                                          !MLPK 
      GO TO 50                                                           !MLPK    
!-----SEEK instruction                                                   TEXT 
  330 IF(KILL.NE.0)GO TO 990                                             !MLPK    
      IF(NINCL.EQ.0)GO TO 905                                            !JH-95 
      IF(ISEP.NE.0)GO TO 340                                             !MLPK    
      IF(IEGY.NE.0)GO TO 350                                             !MLPK    
!-----Can standard SEEK routine be called?                               TEXT 
      IF((NCTR.EQ.0).OR.(NCTR.EQ.7))GO TO 331                            !MLPK 
      CALL SEEKA                                                         !MLPK 
      IF(KILL.NE.0)GO TO 990                                             !MLPK    
      GO TO 50                                                           !MLPK    
  331 CALL SEEK                                                          !MLPK    
      IF(KILL.NE.0)GO TO 990                                             !MLPK    
      GO TO 50                                                           !MLPK    
!-----FIND instruction                                                   TEXT 
  332 IF(KILL.NE.0)GO TO 990                                             !MLPK    
      IF(NINCL.EQ.0)GO TO 905                                            !JH-95 
      IF(ISEP.NE.0)GO TO 340                                             !MLPK    
      IF(IEGY.NE.0)GO TO 350                                             !MLPK    
!-----Can standard FIND routine be called?                               TEXT 
      IF((NCTR.EQ.0).OR.(NCTR.EQ.7))GO TO 333                            !MLPK 
      CALL FINDA                                                         !MLPK 
      IF(KILL.NE.0)GO TO 990                                             !MLPK    
      GO TO 50                                                           !MLPK    
  333 CALL FIND                                                          !MLPK    
      IF(KILL.NE.0)GO TO 990                                             !MLPK    
      GO TO 50                                                           !MLPK    
!-----NSEG instruction - specify non-standard number of search segments  TEXT 
  335 READ(BUF,336)NNE,NE,NV                                             !MLPK 
  336 FORMAT(5X,3I4)                                                     !MLPK 
      IF(NNE.LE.0)NNE=1024                                               !MLPK 
      IF(NE.LE.0)NE=MIN0(128,NNE)                                        !MLPK 
      IF(NV.LE.0)NV=MIN0(8,NE)                                           !MLPK 
      GO TO 50                                                           !MLPK 
!-----FILL IN 'OFF-DIAGONAL' TERMS                                       TEXT 
  340 WRITE(NQ,341)                                                      !MLPK    
  341 FORMAT(' INTERACTION ENERGY PARAMETERS BETWEEN ATOMS OF DIFFERENT &
     & DITYPES'/' AT1 TO AT2  =  C1/D**6  + C2/D**12')          
      DO 349 I=1,MSEP-1                                                  !MLPK    
      II=20*(I-1)+I                                                      !MLPK    
      DO 348 J=I+1,MSEP                                                  !MLPK    
      IJ=20*(I-1)+J                                                      !MLPK    
      IF(EC1(IJ).NE.0.)GO TO 344                                         !MLPK    
      JJ=20*(J-1)+J                                                      !MLPK    
      D=EC1(II)+EC1(JJ)                                                  !MLPK    
      IF(D.EQ.0.)GO TO 342                                               !MLPK    
      EC1(IJ)=2.*EC1(II)*EC1(JJ)/D                                       !MLPK    
  342 D=EC2(II)+EC2(JJ)                                                  !MLPK 
      IF(D.EQ.0.)GO TO 343                                               !MLPK    
      EC2(IJ)=2.*EC2(II)*EC2(JJ)/D                                       !MLPK    
      JI=20*(J-1)+I                                                      !MLPK    
  343 EC1(JI)=EC1(IJ)                                                    !MLPK    
      EC2(JI)=EC2(IJ)                                                    !MLPK    
      GO TO 346                                                          !JH-01 
  344 WRITE(NQ,345)TDAA(I),TDAA(J),EC1(IJ),EC2(IJ)                       !JH-01 
  345 FORMAT(2X,A2,5X,A2,2E12.4,' - stored')                             !JH-01 
      GO TO 348                                                          !JH-01 
  346 WRITE(NQ,347)TDAA(I),TDAA(J),EC1(IJ),EC2(IJ)                       !MLPK    
  347 FORMAT(2X,A2,5X,A2,2E12.4,' - calculated')                         !JH-01    
  348 CONTINUE                                                           !MLPK    
  349 CONTINUE                                                           !MLPK    
      ISEP=0                                                             !MLPK    
      GO TO 75                                                           !MLPK    
!-----CALCULATE ENERGY VALUES AT INTERVALS OF D**2                       !MLPK    
  350 E1=4.*ERM                                                          !MLPK    
      E2=ERM*FERM                                                        !MLPK    
      IF(IPR.LT.5)GO TO 352                                              !MLPK    
      WRITE(NQ,351)E1,E2                                                 !MLPK    
  351 FORMAT('0INTERACTION ENERGIES WITH E(MAX) =',F6.2,' (4*ERM), &  
     &E(MINN) =',F6.2,'(E(MAX)*FERM)'/' AT1 TO AT2 D(MAX)',21X, & 
     &'ENERGIES',26X,'D**2  DELTA D**2')                                 !MLPK    
  352 DO 370 I=1,NSEP                                                    !MLPK    
      II=JSEP(I)  
      DO 369 J=I,NSEP                                                    !MLPK    
      JJ=JSEP(J)  
      JII=20*(JJ-1)+II  
      IJ=120*(I-1)+12*(J-1)+1                                            !MLPK    
      JI=10*(J-1)+I                                                      !MLPK    
      D=2.*EC2(JII)/(2.*SQRT(E2*EC2(JII))-EC1(JII))    
      CN(IJ+10)=D**.333333333                                            !MLPK    
      CN(IJ)=SQRT(CN(IJ+10))                                             !MLPK    
      D=2.*EC2(JII)/(2.*SQRT(E1*EC2(JII))-EC1(JII))     
      DDD=(CN(IJ+10)-D**.333333333)/9.                                   !MLPK    
      CN(IJ+11)=DDD                                                      !MLPK    
      D=CN(IJ+10)-.5*DDD                                                 !MLPK    
      E=EC1(JII)**2/(4.*EC2(JII))    
      DO 355 K=IJ+1,IJ+9                                                 !MLPK    
      F=D**3                                                             !MLPK    
      CN(K)=EC1(JII)/F+EC2(JII)/F**2+E    
  355 D=D-DDD                                                            !MLPK    
      IF(I.EQ.J)GO TO 360                                                !MLPK    
      JI=120*(J-1)+12*(I-1)                                              !MLPK    
      DO 356 K=IJ,IJ+11                                                  !MLPK    
      JI=JI+1                                                            !MLPK    
  356 CN(JI)=CN(K)                                                       !MLPK    
  360 IF(IPR.LT.5)GO TO 369                                              !MLPK    
      WRITE(NQ,361)TDA(I),TDA(J),(CN(K),K=IJ,IJ+11)                      !MLPK    
  361 FORMAT(2X,A2,5X,A2,F7.3,3X,9F5.2,3X,2F8.3)                         !MLPK    
  369 CONTINUE                                                           !MLPK    
  370 CONTINUE                                                           !MLPK    
      IEGY=0                                                             !MLPK    
      GO TO 75                                                           !MLPK    
!-----CENT instruction                                                   TEXT 
  400 READ(BUF,401)FLAG,EXPD,ORIG                                        !MLPK 
  401 FORMAT(7X,2A2,3F8.0)                                               !MLPK 
      DO 405 L=1,3                                                       !MLPK    
  405 C(L)=-1.                                                           !MLPK    
!-----Set centric molecule flag unless voided                            TEXT 
      IF(FLAG.EQ.BLANK)NCTR=7                                            !MLPK 
      GO TO 422                                                          !MLPK    
!-----AXIS instruction                                                   TEXT 
  410 READ (BUF,411) M,FLAG,EXPD,ORIG                                    !MLPK 
  411 FORMAT(5X,I1,1X,2A2,3F8.0)                                         !MLPK 
      IF((M.LE.0).OR.(M.GT.3))GO TO 901                                  !MLPK    
      DO 412 L=1,3                                                       !MLPK    
      C(L)=-1.                                                           !MLPK    
      IF(M.EQ.L)C(L)=1.                                                  !MLPK    
  412 CONTINUE                                                           !MLPK    
!-----Set molecule contains axis flag unless voided                      TEXT 
      IF(FLAG.EQ.BLANK)NCTR=M+3                                          !MLPK 
      GO TO 422                                                          !MLPK    
!-----PLAN instruction                                                   TEXT 
  420 READ (BUF,411) M,FLAG,EXPD,ORIG                                    !MLPK 
      IF((M.LE.0).OR.(M.GT.3))GO TO 901                                  !MLPK    
      DO 421 L=1,3                                                       !MLPK    
      C(L)=1.                                                            !MLPK    
      IF(M.EQ.L)C(L)=-1.                                                 !MLPK    
  421 CONTINUE                                                           !MLPK    
!-----Set molecule on mirror plane flag unless voided                    TEXT 
      IF(FLAG.EQ.BLANK)NCTR=M                                            !MLPK 
  422 IF(NMOD.LE.0)GO TO 903                                             !MLPK    
!-----Relocate origin                                                    TEXT 
      WRITE(NQ,423)                                                      !MLPK 
  423 FORMAT(' Relocated and/or Expanded Molecule')                      !MLPK 
      DO 425 I=1,NMOD                                                    !MLPK 
      DO 424 J=1,3                                                       !MLPK 
      X(J,I)=X(J,I)-ORIG(J)                                              !MLPK 
  424 CONTINUE                                                           !MLPK 
!     WRITE(NQ,91)I,AT(I),AN(I),(X(L,I),L=1,3),IA(I),IAA(I)              !MLPK 
      WRITE(NQ,92) AT(I), AN(I), (X(L,I),L=1,3),&           
     &           IDATOM(I),G92CHARGE(I) 
  425 CONTINUE                                                           !MLPK 
   92 FORMAT(' ATOM ',A2,A4,3F10.6,I5,F10.6)        
!-----Do not expand atom list if operation voided                        TEXT 
      IF(EXPD.NE.BLANK)GO TO 50                                          !MLPK 
      K=NMOD                                                             !MLPK    
      DO 435 I=1,K                                                       !MLPK    
      DO 427 J=1,K                                                       !MLPK    
      DO 426 L=1,3                                                       !MLPK    
      IF(ABS(X(L,I)-C(L)*X(L,J)).GE..1)GO TO 427                         !MLPK    
  426 CONTINUE                                                           !MLPK    
      GO TO 435                                                          !MLPK    
  427 CONTINUE                                                           !MLPK    
      NMOD=NMOD+1                                                        !MLPK    
      DO 428 L=1,3                                                       !MLPK    
  428 X(L,NMOD)=C(L)*X(L,I)                                              !MLPK    
      AT(NMOD)=AT(I)                                                     !MLPK    
      IA(NMOD)=IA(I)                                                     !MLPK    
      IAA(NMOD)=IAA(I)                                                   !MLPK  
      IDATOM(NMOD)=IDATOM(I)                                             !DU-95   
      G92CHARGE(NMOD)=G92CHARGE(I)                                       !DU-95 
      READ (AN(I),429) (CC(L),L=1,4)                                     !MLPK 
  429 FORMAT(4A1)                                                        !MLPK    
      L=3                                                                !MLPK    
  430 L=L-1                                                              !MLPK    
      IF(L.LT.1)GO TO 431                                                !MLPK    
      IF(CC(L).NE.BLANK)GO TO 430                                        !MLPK    
      CC(L)=CARD                                                         !MLPK    
  431 WRITE (AN(NMOD),429) (CC(L),L=1,4)                                 !MLPK 
!     WRITE(NQ,91)NMOD,AT(NMOD),AN(NMOD),(X(L,NMOD),L=1,3),IA(NMOD),     !MLPK    
!    1IAA(NMOD)                                                          !MLPK    
      WRITE (NQ,92) AT(NMOD), AN(NMOD), (X(L,NMOD),L=1,3),&           
     &              IDATOM(NMOD),G92CHARGE(NMOD)     
  435 CONTINUE                                                           !MLPK    
      GO TO 50                                                           !MLPK    
!-----FNDS Instruction                                                   TEXT 
  490 IF(KILL.NE.0)GO TO 990                                             !RTPK    
      READ(BUF,491)JCYC,NCOUNT,STEP1,STEP2,STEP3                         !JH0702 
      NFNDS=1                                                            !RTPK 
      OPEN(18, FILE='fort.18')                                           !RTPK
 
  491 FORMAT(4X,I2,I4,3F6.1)                                             !JH0702 
      IF((NCOUNT.LE.0).OR.(NCOUNT.GT.ICOUNT))NCOUNT=ICOUNT               !RTPK 
      IF(JCYC.LE.0)JCYC=2                                                !RTPK 
      IF(STEP1.LE..1)STEP1=DIJ(498)                                      !JH0702 
      IF(STEP2.LE..1)STEP2=DIJ(499)                                      !JH0702 
      IF(STEP3.LE..1)STEP3=DIJ(500)                                      !JH0702 
      JCOUNT=ICOUNT                                                      !RTPK 
      DO 492 I=1,ICOUNT                                                  !RTPK 
      JORDER(I)=IORDER(I)                                                !RTPK 
      ANGL1(I)=ANGLE1(I)                                                 !JH0702 
      ANGL2(I)=ANGLE2(I)                                                 !JH0702 
      ANGL3(I)=ANGLE3(I)                                                 !JH0702 
  492 CONTINUE                                                           !RTPK 
!-----Reset starting values for saving minimum volume structures         TEXT 
      ICOUNT=1                                                           !RTPK 
      IORDER(1)=1                                                        !RTPK 
      VOLUME(1)=9999.                                                    !RTPK 
      ANGLE1(1)=999.                                                     !RTPK 
      ANGLE2(1)=999.                                                     !RTPK 
      ANGLE3(1)=999.                                                     !RTPK 
      CODE(1)=BLANK                                                      !RTPK 
      DO 495 I=1,NCOUNT                                                  !RTPK 
      J=JORDER(I)  
      WRITE(BUF,493)JCYC,ANGL1(J),ANGL2(J),ANGL3(J),STEP1,STEP2,STEP3    !JH0702 
  493 FORMAT('FIND',I2,6F6.1)                                            !JH0702 
!-----Can standard FIND routine be called?                               TEXT 
      IF((NCTR.EQ.0).OR.(NCTR.EQ.7))GO TO 494                            !RTPK 
      CALL FINDA                                                         !RTPK 
      IF(KILL.NE.0)GO TO 990                                             !RTPK    
      GO TO 495                                                          !RTPK    
  494 CALL FIND                                                          !RTPK    
      IF(KILL.NE.0)GO TO 990                                             !RTPK    
  495 CONTINUE                                                           !RTPK 
      GO TO 50                                                           !RTPK    
  500 IF(NFNDS.EQ.1)WRITE(18,499)HEAD                                    !RTPK 
  499 FORMAT(' MOLPAK output for ',A60)                                  !RTPK 
!-----Write out smallest cell volumes if stored                          TEXT 
      WRITE(NQ,501)MCOUNT                                                !MLPK 
      IF(NFNDS.EQ.1)WRITE(18,501)MCOUNT                                  !RTPK 
  501 FORMAT(4H0THE,I4,' SMALLEST VOLUMES/MOLECULE WERE REQUESTED')      !MLPK 
      WRITE(NQ,502)                                                      !MLPK 
      IF(NFNDS.EQ.1)WRITE(18,502)                                        !RTPK 
  502 FORMAT(1H0,'NO.  CD      A1       A2       A3       V/MOL')        !MLPK 
      IF(NT.LE.0)GO TO 505                                               !7/93DU 
      WRITE(NT,504)'    NO','ANGLE1','ANGLE2','ANGLE3','VOLUME'          !7/93DU 
  504 FORMAT(1X,4(1X,A7),2X,A8)                                          !7/93DU 
  505 DO 509 I=1,ICOUNT                                                  !MLPK 
      J=IORDER(I)                                                        !MLPK 
      IF(ANGLE1(J).GT.998.)GO TO 1000                                    !MLPK 
      IF(NT.LE.0)GO TO 507                                               !7/93DU 
      IA1 = ANGLE1(J)                                                    !7/93DU 
      IA2 = ANGLE2(J)                                                    !7/93DU 
      IA3 = ANGLE3(J)                                                    !7/93DU 
!     WRITE(NT,506)I,IA1,IA2,IA3,VOLUME(J)                               !7/93DU 
! 506 FORMAT(1X,4I8,F10.3)                                               !7/93DU 
      WRITE(NT,506)I, ANGLE1(J),ANGLE2(J),ANGLE3(J),VOLUME(J)            !DU02  
  506 FORMAT(1X,I8,3F8.1,F10.3)                                          !DU02 
  507 WRITE(NQ,508)I,CODE(J),ANGLE1(J),ANGLE2(J),ANGLE3(J),VOLUME(J)     !MLPK 
      IF(NFNDS.EQ.1)WRITE(18,508)I,CODE(J),ANGLE1(J),ANGLE2(J),ANGLE3(J) & !MLPK 
     &,VOLUME(J)                                                         !MLPK 
  508 FORMAT(1H ,I6,2X,A2,3F9.2,F11.3)                                   !MLPK 
  509 CONTINUE                                                           !MLPK 
      GO TO 1000                                                         !MLPK 
  901 WRITE(NQ,902)CARD                                                  !MLPK    
      WRITE(17,902)CARD                                                  !MLPK    
  902 FORMAT(6H0BAD (,A4,6H) LINE)                                       !MLPK    
      KILL=1                                                             !MLPK    
      GO TO 50                                                           !MLPK    
  903 WRITE(NQ,904)CARD                                                  !MLPK    
      WRITE(17,904)CARD                                                  !MLPK    
  904 FORMAT('0MISSPLACED (',A4,') LINE')                                !MLPK    
      KILL=1                                                             !MLPK    
      GO TO 50                                                           !MLPK    
  905 WRITE(NQ,906)CARD                                                  !JH-95 
      WRITE(17,906)CARD                                                  !JH-95 
  906 FORMAT('0MISSING (INCL) INSTRUCTION BEFORE (',A4,') LINE')         !JH-95 
  990 WRITE(NQ,991)                                                      !MLPK    
      WRITE(17,991)                                                      !MLPK    
  991 FORMAT(22H0****MOLPAK ERROR STOP)                                  !MLPK    
 1000 CONTINUE                                                           !DU02 
!1000 PRINT *,(' MOLPAK is finished')                                    !MLPK 
      WRITE (NQ,1001)                                                    !MLPK 
 1001 FORMAT (' **MOLPAK end**MOLPAK end**MOLPAK end**MOLPAK end**', &   !MLPK 
     &        'MOLPAK end**MOLPAK end**MOLPAK end**MOLPAK end**')        !MLPK 
      CLOSE(22)                                                          !MLPK 
      CLOSE(23)                                                          !MLPK 
!     print coefficients of 6-12 potential funtions                      !DU-95 
      IF (IPR.EQ.4) THEN                                                 !DU-95 
      DO I=1,10                                                          !DU-95  
       DO J=1,10                                                         !DU-95 
        IJ= 10*(I-1)+J                                                   !DU-95 
        JI= 10*(J-1)+I                                                   !DU-95 
        PRINT 6,TDA(I),TDA(J),EC1(IJ),EC2(IJ),EC1(JI),EC2(JI),I,J,IJ,JI  !DU-95 
  6     FORMAT (1X,A2,'--',A2,2X,4E12.4,2X,4I4)                          !DU-95 
       END DO  
      END DO                                                             !DU-95 
      END IF   
      STOP                                                               !MLPK
      END PROGRAM MOLPAK                                                 !MLPK
