      module molpakCommonMod
      
      IMPLICIT NONE
      
      CHARACTER (80) :: BUF                                                       
      CHARACTER (60) :: HEAD                                                  
      CHARACTER (4)  :: AN(200), CDIJ(500), CHT1(14)          
      CHARACTER (2)  :: AT(200), CODE(60000) 
                   
      INTEGER :: IAA(200), IA(200), ICOUNT, IJD(500), IORDER(60000), IT(1000000)
      INTEGER :: IPR, ITR(90), JTR(30), KCT, KILL, LIMIT, LINE   
      INTEGER :: MARK, MARK1, MARK2, MARK3, MARK4, MAXT, MINT, MLOW, MCOUNT
      INTEGER :: NARK, NCTR, NE, NEW, NEND, NEND1, NEND2, NEND3, NEND4
      INTEGER :: NLOW, NMOD, NNE, NPAGE, NP2, NQ, NR, NS, NSEP, NSGT
      INTEGER :: NSTP, NSTP1, NSTP2, NSTP3, NSTP4, NV, NWM                                    
      INTEGER :: NRB , NRB_AT(5)  ! the numberof rigid bodies  7-20-06
      INTEGER :: NPMIN            ! 1-29-09  DU
      integer :: n_cross   ! # sets of cross-term lines for PMIN (max = 5)  6/24/09
           
      REAL :: A1, A2, A3, ANGLE1(60000),ANGLE2(60000),ANGLE3(60000)
      REAL :: CS2(37), CS3(19), CN(1200), DIJ(500), ERM, ER, ERMT, ERMF                                               
      REAL :: FERM, H1M, H2M, H3M, H1N, H2N, H3N, HT1(14)                                
      REAL :: SN2(37), SN3(19), TI(1000000), V3, V(3,200), VOLUME(60000) 
      REAL :: W(3,200), X(3,200)                                           
      
!     EQUIVALENCE (TI, IT)
!     EQUIVALENCE (DIJ, IJD, CDIJ)
!     EQUIVALENCE (HT1, CHT1)
      
!----For /NEW/...
      CHARACTER(40) :: POTENTIAL_FILE       ! file with potential params
      character (len=40), dimension(5) :: cross_terms    ! 5 sets of cross-term coefs for PMIN  6/24/09
      INTEGER :: IDATOM(200)
      REAL :: G92CHARGE(200) 
      
!---For /ATMWT/...   Du 04
      INTEGER :: JATYPE(13)
      REAL :: ATWT(13)
!   6-9-07 
      CHARACTER(2) :: CSYM(25)
      
      END module molpakCommonMod                                               
