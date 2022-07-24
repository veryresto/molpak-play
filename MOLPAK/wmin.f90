      SUBROUTINE WMIN       
!
!----Modified to write PMIN cross-term sets on the ....input file   6/24/09 
!      
      USE molpakCommonMod
      
      IMPLICIT NONE

!      CHARACTER*80 BUF          
!      CHARACTER*60 HEAD         
!      CHARACTER*4 AN,CHT1         
!      CHARACTER*2 AT,CODE      
!-----General program communication and transfer variables                
!      COMMON BUF,KILL,NPAGE,LINE,NEW,JTR(30),ITR(90)                      
!-----Energy constants                                                        
!      COMMON NSEP,CN(1200),ERM,FERM,ER,ERMT,ERMF,NV,NE,NNE       
!-----Compound input data                                                 
!      COMMON NMOD,NCTR,IA(200),IAA(200),X(3,200),AT(200),AN(200)            
!-----Compound rotation and structure parameters                          
!      COMMON CS2(37),SN2(37),CS3(19),SN3(19),V(3,100),W(3,100)            
!      COMMON NLOW,MLOW,NP2,KCT                                            
!      COMMON H1M,H2M,H3M,H1N,H2N,H3N,A1,A2,A3,V3                          
!-----Output flags and storage                                            
!      COMMON NQ,NR,NS,NWM,NSGT,IPR,HEAD,DIJ(500),HT1(14)                  
!      DIMENSION IJD(500),CHT1(14)                                         
!      EQUIVALENCE (DIJ(1),IJD(1)),(HT1(1),CHT1(1))                        
!     COMMON ICOUNT,MCOUNT,IORDER(500)                                    
!      COMMON ICOUNT,MCOUNT,IORDER(5000)                                   
!     COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500)    
!      COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000), &
!     &       CODE(5000)
!-----Linear storage array and markers                                    
!      COMMON MAXT,MINT,LIMIT                                              
!      COMMON MARK,NEND,NSTP,NARK                                          
!      COMMON MARK1,NEND1,NSTP1                                                
!      COMMON MARK2,NEND2,NSTP2                                                
!      COMMON MARK3,NEND3,NSTP3                                                
!      COMMON MARK4,NEND4,NSTP4                                                
!      COMMON TI(100000)                                                     
!      COMMON /NEW/IDCHG,IDATOM(200),G92CHARGE(200),ESPCHARGE(200)         
!      DIMENSION IT(100000)                                                
!      EQUIVALENCE (TI(1),IT(1))                                          
!      COMMON /ATMWT/ATWT(10),JATYPE(10)

      CHARACTER(1)  :: P(6)                                                                                                          
      CHARACTER(2)  :: ATN
      CHARACTER(60) :: HHEAD 
      
      INTEGER  :: I,IA1, IA2, IA3, ICENT, IHB, IHXNEG, IHYNEG, IHZNEG
      INTEGER  :: IRIGID, IREP, IOX, IS(3,8), J, JJ, K, KK, L, LL
      INTEGER  :: M, NCN, NKIND, NMODD, NNKIND, NOP, NTA  
      integer :: idchg = 2     ! temporary
      integer :: LH,ISML            ! 5-9-05 Du 
      INTEGER  :: I64, I65, I61, I63     ! for CCSP_1 & CCSP_3 
!
      REAL     :: SML               ! 5-9-05 Du
      REAL     :: A(3),C(3), AA(3), CC(3), CS
      REAL     :: D, D1, D2, D3, DSA(3), DSA4, DSA5, DSA6
      REAL     :: E, F, FXX, FXY
      REAL     :: H1I, H2I, H21I, H3I, H31I, H32I
      REAL     :: SI(3,8), SINE1, SINE2, SINE3, SN
      REAL     :: VXX, VYY, VZZ, WMOL 
      REAL     :: PC_A(200), PC_B(200), PC_C(200), PC_Q(200), PC_AW(200)                 ! 1-29-09 DU
      REAL     :: PM_X(200), PM_Y(200), PM_Z(200)                                        ! 1-29-09 DU
      INTEGER  :: PC_CODE(200),I_flag(200)                                               ! 1-29-09 DU
      CHARACTER(4)  :: PM_AT(200)                                                        ! 1-29-09 DU
      CHARACTER(72) :: PM_RB                                                             ! 1-29-09 DU
      
!-----Space group (subgroup) codes and corresponding SBN numbers          
      CHARACTER(2)  :: SYM(83)=(/'AA','AB','AC','AD','AH','AE','AF',   &
     & 'AG','CA','AI','AJ','AK','AL','AM','AN','AO','AP','AQ','AR','AS',&
     & 'AT','AU','AV','AW','AX','AY','AZ','BA','BB','BC','BD','BE','BF',&
     & 'BG','BH','BI','BJ','BK','CB','CC','CD','CE','CF','DA','DB','DC',&
     & 'DD','DE','DF','DG','EA','EB','EC','ED','EE','EF','EG','SA','SB',&
     & 'SC','SD','SE','SF','SG','SH','SI','SJ','SK','SL','SM','SN','FA',&
     & 'FB','FC','FD','SO','SP','SQ','SR','SS','ST','SU','SV'/)                     
!-----Signs in symops for nonprimative orthorhombics based on MLOW-50     
      INTEGER :: IW(7)= (/1,1,1,1,-1,-1,-1/)      
      INTEGER :: JW(7)= (/1,-1,1,-1,1,-1,-1/)          
      INTEGER :: KW(7)= (/1,1,-1,-1,-1,1,-1/)           
!-----Set axis conversion variables depending on MLOW                     
      INTEGER  :: I1(83)=(/ 1,1,1,2,1, 1,1,1,1,2, 2,2,2,1,1,  &
    & 1,1,1,1,2, 3,2,1,1,3, 2,1,1,1,3, 1,1,2,1,1, 2,1,1,1,1, 1,1,1,1,1, & 
    & 1,1,1,1,1, 1,1,1,1,1, 1,1,1,2,2, 1,1,1,1,1, 1,2,2,3,3, 2,2,2,1,1, &
    & 2,2,2,2,2, 1,1,3/)                                         
      INTEGER  :: I2(83)=(/ 2,2,3,3,3, 3,3,3,2,3, 3,3,3,3,3,  &
    & 3,2,2,2,1, 2,3,2,2,2, 3,2,2,3,2, 2,2,1,2,2, 1,2,2,2,2, 2,2,2,3,3, &
    & 3,3,3,3,3, 2,2,2,2,2, 2,2,2,3,3, 3,3,2,3,2, 2,1,1,2,2, 3,3,3,3,3, &
    & 3,3,3,3,3, 2,2,2/)                                          
      INTEGER  :: I3(83)=(/ 3,3,2,1,2, 2,2,2,3,1, 1,1,1,2,2,  &
    & 2,3,3,3,3, 1,1,3,3,1, 1,3,3,2,1, 3,3,3,3,3, 3,3,3,3,3, 3,3,3,2,2, &
    & 2,2,2,2,2, 3,3,3,3,3, 3,3,3,1,1, 2,2,3,2,3, 3,3,3,1,1, 1,1,1,2,2, &
    & 1,1,1,1,1, 3,3,1/)                                                      
!-----SBN numbers corresponding to values of MLOW !! COMPLETED ON 10-12-02 
      INTEGER :: ISB(83)=(/    19,  19,  19,  19,  19,  19,  &
     & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  &
     & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  &
     & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  &
     & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  &
     & 19,  19 , 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  &
     & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19/)           
!-----Atomic weights                                                      
!     DATA WATM/12.0111,1.0080,1.0080,14.0067,15.9994,15.9994,18.9984   &
!    &,35.4530,79.9040,10.811,32.066/                                      
!-----Table to convert program atom types to WMIN atom types              
!     DATA IATYPE/1,4,4,2,3,3,5,6,0,0,0/                                   
!-----Quit if no structure has been found                                 
      IF(NLOW.LE.0)GO TO 900                                                 
!-----NR for WMIN input with the old energy coefficients (no charges)     
!-----Set the number of new WMIN input files with NEC & charges            
!     seeing below                                                         
!-----19 for WMIN input with new energy coefficients (no charges)          
!-----29 for WMIN input with new energy coefficients & Gaussian 92        
!            6-31g*/CHELPG atomic charges                                 
!-----39 for WMIN input with new energy coefficients & MOPAC MNDO/ESP     
!            atomic charges                                               
!     write PMIN.inp if NPMIN = 1                     ! 1-29-09 DU      
      if (NPMIN == 1) then      
         OPEN (UNIT=59, file='PMIN.inp', STATUS='UNKNOWN')
         write(59,1001) n_cross            ! # PMIN cross_terms       ! 6/24/09
1001        format('PMIN.inp'/'2 2 0 0 0 1 0.5 0.25 6.0 50',i3,' 0 0'/ & ! 6/24/09
     &      '50 100 0.00001 0.00001 0.00001 0.00001 1.0 1.0 1.0 1.0')
         if (n_cross .ne. 0) then    ! write cross-terms in n_cross .le. 1  6/24/09
            write (59,1010) (cross_terms(i),i=1,n_cross)               ! 6/24/09
1010           format (a40)                                            ! 6/24/09
         endif                                                         ! 6/24/09
      print 9191, n_cross
9191     format ('### n_cross =',i3)                                   ! temp
      print 9194, (cross_terms(i),i=1,n_cross)                         ! temp
9194     format ('### |',a40)                                          !temp
      end if                                          ! 1-29-09 DU
!     Set c glide to n glide flag to no change                            
      NCN=0                                                               
!     Set triclinic angle adjustment flag to no change                    
      NTA=0                                                        
!-----Rotate molecule to angles of minimum volume structure               
!     PRINT '(7HIDCHG= ,I4)', IDCHG
      D=HT1(10)/57.296                                                       
      CS=COS(D)                                                              
      SN=SIN(D)                                                              
      DO 10 I=1,NMOD                                                         
      V(1,I)=X(1,I)*CS-X(2,I)*SN                                          
      W(2,I)=X(1,I)*SN+X(2,I)*CS                                          
   10 CONTINUE                                                            
      D=HT1(11)/57.296                                                       
      CS=COS(D)                                                              
      SN=SIN(D)                                                              
      DO 11 I=1,NMOD                                                         
      W(1,I)=V(1,I)*CS+X(3,I)*SN                                          
      V(3,I)=-V(1,I)*SN+X(3,I)*CS                                         
   11 CONTINUE                                                            
      D=HT1(12)/57.296                                                       
      CS=COS(D)                                                              
      SN=SIN(D)                                                              
      DO 12 I=1,NMOD                                                         
      V(1,I)=W(1,I)*CS-W(2,I)*SN                                          
      V(2,I)=W(1,I)*SN+W(2,I)*CS                                          
   12 CONTINUE                                                            
!-----Transform from orthogonal Angstrom to unit cell space               
      H1I=HT1(1)                                                          
      H21I=HT1(2)                                                         
      H2I=HT1(3)                                                          
      H3I=HT1(4)                                                          
      H31I=HT1(13)                                                        
      H32I=HT1(14)                                                        
!-----Move to new origin                                                  
      D1=-HT1(5)/2.                                                       
      D2=-HT1(6)/2.                                                       
      D3=-HT1(7)/2.                                                       
      DO 13 I=1,NMOD                                                      
      V(1,I)=V(1,I)+D1                                                    
      V(2,I)=V(2,I)+D2                                                    
   13 V(3,I)=V(3,I)+D3                                                    
!-----Is "no space group change" flag turned on?                          
      IF(NSGT.NE.0)GO TO 23                                               
!-----Should the structure be adjusted toward 90 degree angles?           
      GO TO(15,15,21,21,21, 21,15,21,21,21, 21,23,23,23,23,    &           
     &      23,23,23,23,23, 23,21,21,21,23, 23,23,23,23,21,    &           
     &      21,23,23,23,23, 23),NLOW                                          
!-----Triclinic structures                                                
   15 IF(H21I.LE.H1I/2.)GO TO 16                                          
      H21I=H21I-H1I                                                       
      NTA=1                                                               
      GO TO 17                                                            
   16 IF(H21I.GT.-H1I/2.)GO TO 17                                         
      NTA=1                                                               
      H21I=H21I+H1I                                                       
   17 IF(H32I.LE.H2I/2.)GO TO 18                                          
      NTA=NTA+2                                                           
      H32I=H32I-H2I                                                       
      H31I=H31I-H21I                                                      
      GO TO 19                                                            
   18 IF(H32I.GT.-H2I/2.)GO TO 19                                         
      NTA=NTA+2                                                           
      H32I=H32I+H2I                                                       
      H31I=H31I+H21I                                                      
   19 IF(H32I.LE.H2I/2.)GO TO 20                                          
      NTA=NTA+4                                                           
      H31I=H31I-H1I                                                       
      GO TO 23                                                            
   20 IF(H31I.GT.-H1I/2.)GO TO 23                                         
      NTA=NTA+4                                                           
      H31I=H31I+H1I                                                       
      GO TO 23                                                            
!-----Make monoclinic unit cell angle near 90 degrees                     
   21 IF(H21I.LT.H1I/2.)GO TO 22                                          
      NCN=-1                                                              
      H21I=H21I-H1I                                                       
      GO TO 23                                                            
   22 IF(H21I.GT.-H1I/2.)GO TO 23                                         
      NCN=1                                                               
      H21I=H21I+H1I                                                       
!-----Correct shape of molecule to that input to MOLPAK when required.    
   23 GO TO(26,28,26,26,24, 26,24,26,28,28, 26,28,28,24,26, &              
     &      26,26,28,26,26, 26,26,28,26,26, 28,28,26,24,26, &              
     &      28,26,24,26,28, 26,26,26,28,28, 26,26,26,24,24, &              
     &      24,24,24,26,26, 26,26,26,26,26, 26,26,26,26,26, &              
     &      26,26,26,26,26, 26,26,26,26,26, 26,26,26,24,26, &              
     &      26,26,26,26,26, 26,26,26),MLOW                                
!     Change the "hand" of the output molecule.                           
   24 DO 25 I=1,NMOD                                                      
   25 V(3,I)=-V(3,I)                                                      
      H31I=-H31I                                                          
      H32I=-H32I                                                          
      GO TO 28                                                            
!     The "hand" of the output molecule has not been checked.             
   26 WRITE(NQ,27)                                                        
   27 FORMAT(' The "hand" of the output molecule has not been checked.')  
!-----Calculate space group unit cell dimensions                          
!     Cell axis lengths                                                   
   28 A(1)=H1I                                                            
      A(2)=SQRT(H2I**2+H21I**2)                                           
      A(3)=SQRT(H3I**2+H31I**2+H32I**2)                                   
!     Cosines of cell angles                                              
      C(1)=(H21I*H31I+H2I*H32I)/(A(2)*A(3))                               
      C(2)=H31I/A(3)                                                      
      C(3)=H21I/A(2)                                                      
!-----Transform atoms from orthogonal to cell coordinates                 
      SINE2=SQRT(1.-C(2)**2)                                              
      SINE3=SQRT(1.-C(3)**2)                                              
      E=(C(2)*C(3)-C(1))/(SINE2*SINE3)                                       
      F=SQRT(1.-E**2)                                                     
      FXX=SINE2*F                                                            
      FXY=SINE2*E/SINE3                                                      
      DO 29 I=1,NMOD                                                      
!-----Transform first to Angstrom coordinates                             
      V(3,I)=V(3,I)/FXX                                                   
      V(2,I)=V(2,I)/SINE3+V(3,I)*FXY                                      
      V(1,I)=V(1,I)-V(3,I)*C(2)-V(2,I)*C(3)                               
!-----Then to fractional coordinates                                      
      V(1,I)=V(1,I)/A(1)                                                  
      V(2,I)=V(2,I)/A(2)                                                  
      V(3,I)=V(3,I)/A(3)                                                  
   29 CONTINUE                                                            
!-----Assemble space group dependent data                                 
!     Set symop variables to tentative values                             
   37 NOP=4                                                               
      DO 39 J=1,8                                                         
      DO 38 I=1,3                                                            
      IS(I,J)=1                                                              
   38 SI(I,J)=0.                                                             
   39 CONTINUE                                                            
      NMODD=NMOD+1                                                        
!-----Tentatively set WMIN refinement flags for orthorhombic space groups 
      IHXNEG=0                                                            
      IHYNEG=0                                                            
      IHZNEG=0                                                            
      DSA(1)=0.                                                           
      DSA(2)=0.                                                           
      DSA(3)=0.                                                           
      GOTO(100,101,119,110,114 ,118,114,111,101,121 ,122,121,122,124,126 & 
     &    ,129,132,131,134,135 ,136,137,138,139,140 ,141,131,132,143,144 & 
     &    ,145,146,147,148,149 ,150,151,152,190,190 ,189,189,902,180,181 & 
     &    ,183,183,183,902,902 ,160,160,160,160,160 ,160,160,191,193,192 & 
     &    ,194,195,196,197,198 ,199,200,201,202,203 ,204,121,126,124,123 & 
     &    ,205,205,206,206,206 ,207,207,208,902),MLOW                       
!-----Triclinic space groups                                              
!     SG P1                                                               
  100 NOP=1                                                               
      IHXNEG=1                                                            
      IHYNEG=1                                                            
      DSA(1)=0.001                                                        
      DSA(2)=0.001                                                        
      DSA(3)=0.001                                                        
      GO TO 250                                                           
!     SG P1bar - Subgroups 1 and 2                                        
  101 IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IHXNEG=1                                                            
      IHYNEG=1                                                            
      DSA(1)=0.001                                                        
      DSA(2)=0.001                                                        
      DSA(3)=0.001                                                        
      GO TO 249                                                           
!     SG Pc - Subgroup 1                                                  
  110 SI(1,2)=.5                                                          
      GO TO 119                                                           
!     SG Pc - Subgroup 2                                                  
  111 IF(NCN.NE.0)GO TO 211                                               
      SI(2,2)=.5                                                          
      GO TO 119                                                           
!     SG P21 - Subgroups 1 and 2                                          
  114 SI(3,2)=.5                                                          
!     SG P2 - Subgroups 1 and 2                                           
  118 IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 249                                                           
  119 IS(3,2)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 249                                                           
!-----Monoclinic space groups with Z=4                                    
!     SG P21/c - Subgroups 1, 2, and 4:                                   
  121 SI(3,3)=.5                                                          
      SI(3,4)=.5                                                          
!     SG P2/c - Subgroups 1 and 2                                         
  122 SI(1,3)=.5                                                          
      SI(1,4)=.5                                                          
      GO TO 129                                                           
!     SG P2/c - Subgroup 3                                                
  123 IF(NCN.NE.0)GO TO 215                                               
      SI(2,3)=.5                                                          
      SI(2,4)=.5                                                          
      GO TO 129                                                           
!     SG P21/c - Subgroups 3 and 5                                        
  124 IF(NCN.NE.0)GO TO 216                                               
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 129                                                           
!     SG P21/m - Subgroups 1 and 2                                        
  126 SI(3,3)=.5                                                          
      SI(3,4)=.5                                                          
!     SG P2/m                                                             
  129 IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(2,3)=-1                                                          
      IS(3,4)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=.001                                                         
      GO TO 250                                                           
!-----Orthorhombic space groups with Z=4                                  
!     SG P21212 - Subgroups 1 and 2: SG P212121 - Subgroups 1 and 2       
  131 SI(3,3)=.5                                                          
      SI(3,4)=.5                                                          
  132 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(2,3)=.5                                                          
      SI(1,4)=.5                                                          
      GO TO 158                                                           
!     SG Pnn2 - Subgroup 1                                                
  134 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      GO TO 159                                                           
!     SG Pna21 - Subgroup 1                                               
  135 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 159                                                           
!     SG Pnn2 - Subgroup 2                                                  
  136 SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 157                                                           
!     SG Pna21 - Subgroup 2                                               
  137 SI(1,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 157                                                           
!     SG Pna21 - Subgroup 3                                               
  138 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 159                                                           
!     SG Pba2 - Subgroup 1                                                
  139 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      GO TO 159                                                           
!     SG Pnc2 - Subgroup 1                                                
  140 SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(2,4)=.5                                                          
      GO TO 157                                                           
!     SG Pca21 - Subgroup 1                                               
  141 SI(1,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(2,4)=.5                                                          
      GO TO 157                                                           
!     SG P21212 - Subgroup 3                                              
  143 SI(1,2)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 158                                                           
!     SG P2221                                                            
  144 SI(1,2)=.5                                                          
      SI(1,4)=.5                                                          
      GO TO 158                                                           
!     SG Pna21 - Subgroup 2                                               
  145 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 159                                                           
!     SG Pnn2 - Subgroup 2                                                
  146 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      GO TO 159                                                           
!     SG Pna21 - Subgroup 2                                               
  147 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 159                                                           
!     SG Pba2 -Subgroup 2                                                 
  148 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      GO TO 159                                                           
!     SG Pca21 - Subgroup 2                                               
  149 SI(1,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 159                                                           
!     SG Pnc2 Subgroup 2                                                  
  150 SI(1,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(3,3)=.5                                                          
      GO TO 159                                                           
!     SG Pmn21 - Subgroup 2                                               
  151 SI(1,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(3,4)=.5                                                          
      GO TO 159                                                           
!     SG Pma2 - Subgroup 2                                                
  152 SI(1,2)=.5                                                          
      SI(1,3)=.5                                                          
      GO TO 159                                                           
  157 IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IS(2,3)=-1                                                          
      IS(3,4)=-1                                                          
      GO TO 250                                                           
  158 IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(3,3)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
      GO TO 250                                                           
  159 IS(1,2)=-1                                                          
      IS(2,3)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
      GO TO 250                                                           
!-----Non-primative, face centered orthorhombic space groups              
  160 DO 161 J=5,8                                                           
      IS(1,J)=IW(MLOW-50)                                                    
      IS(2,J)=JW(MLOW-50)                                                    
  161 IS(3,J)=KW(MLOW-50)                                                    
  168 SI(2,2)=.5                                                             
      SI(3,2)=.5                                                             
      SI(1,3)=.5                                                             
      SI(3,3)=.5                                                             
      SI(1,4)=.5                                                             
      SI(2,4)=.5                                                             
      DO 169 J=5,8                                                           
      SI(1,J)=SI(1,J-4)                                                   
      SI(2,J)=SI(2,J-4)                                                   
  169 SI(3,J)=SI(3,J-4)                                                   
      GO TO 248                                                           
!-----Non-primative, C centered monoclinic space groups                   
!     SG Cc - Code DA                                                     
  180 IS(3,2)=-1                                                          
      IS(3,4)=-1                                                          
      SI(2,2)=.5                                                          
      SI(2,4)=.5                                                          
      GO TO 182                                                           
!     SG C2 - Code DB                                                     
  181 IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
  182 SI(1,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(3,4)=.5                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 250                                                           
!     SG C2/c - Codes DC, DD, and DE                                      
  183 DO 184 I=1,3                                                        
      IS(I,2)=-1                                                          
  184 IS(I,6)=-1                                                          
      IS(1,3)=-1                                                          
      IS(2,3)=-1                                                          
      IS(1,7)=-1                                                          
      IS(2,7)=-1                                                          
      IS(3,4)=-1                                                          
      IS(3,8)=-1                                                          
      SI(2,3)=.5                                                          
      SI(2,4)=.5                                                          
      SI(2,7)=.5                                                          
      SI(2,8)=.5                                                          
      DO 185 I=5,8                                                        
      SI(1,I)=.5                                                          
  185 SI(3,I)=.5                                                          
      IHXNEG=1                                                            
      DSA(3)=.001                                                         
      GO TO 248                                                           
!-----Z=8 Orthorhombic space groups                                       
!     SG Pbcn                                                             
  189 IS(1,2)=-1                                                          
      IS(2,3)=-1                                                          
      IS(3,4)=-1                                                          
      IS(2,5)=-1                                                          
      IS(3,5)=-1                                                          
      IS(1,6)=-1                                                          
      IS(3,6)=-1                                                          
      IS(1,7)=-1                                                          
      IS(2,7)=-1                                                          
      IS(1,8)=-1                                                          
      IS(2,8)=-1                                                          
      IS(3,8)=-1                                                          
      SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
      SI(1,5)=.5                                                          
      SI(2,5)=.5                                                          
      SI(3,6)=.5                                                          
      SI(1,7)=.5                                                          
      SI(2,7)=.5                                                          
      SI(3,7)=.5                                                          
      GO TO 248                                                           
!     SG Pbca                                                             
  190 IS(1,2)=-1                                                          
      IS(2,3)=-1                                                          
      IS(3,4)=-1                                                          
      IS(2,5)=-1                                                          
      IS(3,5)=-1                                                          
      IS(1,6)=-1                                                          
      IS(3,6)=-1                                                          
      IS(1,7)=-1                                                          
      IS(2,7)=-1                                                          
      IS(1,8)=-1                                                          
      IS(2,8)=-1                                                          
      IS(3,8)=-1                                                          
      SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(3,4)=.5                                                          
      SI(1,5)=.5                                                          
      SI(2,5)=.5                                                          
      SI(2,6)=.5                                                          
      SI(3,6)=.5                                                          
      SI(1,7)=.5                                                          
      SI(3,7)=.5                                                          
      GO TO 248                                                           
!----Special monoclinic space groups for symmetric molecules              
!     SG A-P2 or P-Pm - Code SA                                           
  191 NOP=1                                                               
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 250                                                           
!     SG  C- P21/c - Code SC                                              
  192 SI(1,2)=.5                                                          
!     SG C-P21/m - Code SB                                                
  193 SI(3,2)=.5                                                          
      IS(3,2)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 249                                                           
!     SG C-P21/c - Code SD                                                
  194 IF(NCN.NE.0)GO TO 218                                               
      SI(3,2)=.5                                                          
!     SG C-P2/c - Code SE                                                 
  195 IF(NCN.NE.0)GO TO 219                                               
      SI(2,2)=.5                                                          
      IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 249                                                           
!-----Special orthorhombic space groups for symmetric molecules           
!     SG A-P21212 - Code SF                                               
  196 SI(1,2)=.5                                                          
!     SG A-P2221 - Code SG                                                
  197 SI(2,2)=.5                                                          
      IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      GO TO 249                                                           
!     SG A-Pnn2 - Code SH                                                 
  198 SI(3,2)=.5                                                          
!     SG A-Pba2 - Code SI                                                 
  199 SI(1,2)=.5                                                          
!     SG A-Pma2 - Code SJ                                                 
  200 SI(2,2)=.5                                                          
      IS(1,2)=-1                                                          
      GO TO 249                                                           
!     SG A-Pnc2 - Code SK                                                 
  201 SI(1,2)=.5                                                          
      SI(3,2)=.5                                                          
      IS(1,2)=-1                                                          
      GO TO 249                                                           
!     SG P-Pmn21 - Code SL                                                
  202 SI(1,2)=.5                                                          
!     SG P-Pma2 - Code SM                                                 
  203 SI(3,2)=.5                                                          
      IS(2,2)=-1                                                          
      GO TO 249                                                           
!     SG C-Pbca - Code SN                                                 
  204 SI(1,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
      IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(3,3)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
      GO TO 250                                                           
!     SG P-Pnma - Subgroups 1 & 2 - Codes SO & SP                         
  205 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      IS(2,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(3,3)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
      IS(3,4)=-1                                                          
      GO TO 250                                                           
!     SG A-Pbcn - Subgroups 1, 2 & 3 - Codes SQ, SR & SS                  
  206 SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      IS(2,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(3,3)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
      IS(3,4)=-1                                                          
      GO TO 250                                                           
  207 SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(2,4)=.5                                                          
      SI(1,5)=.25                                                         
      SI(2,5)=.25                                                         
      SI(3,5)=.25                                                         
      SI(1,6)=.25                                                         
      SI(2,6)=.75                                                         
      SI(3,6)=.75                                                         
      SI(1,7)=.75                                                         
      SI(2,7)=.25                                                         
      SI(3,7)=.75                                                         
      SI(1,8)=.75                                                         
      SI(2,8)=.75                                                         
      SI(3,8)=.25                                                         
      IS(1,5)=-1                                                          
      IS(1,6)=-1                                                          
      IS(1,7)=-1                                                          
      IS(1,8)=-1                                                          
      GO TO 248                                                           
!     SG C-Pbcn - Code SV                                                 
  208 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      SI(3,2)=.5                                                          
      SI(1,3)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
      IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(3,3)=-1                                                          
      IS(1,4)=-1                                                          
      IS(2,4)=-1                                                          
      GO TO 250                                                           
!-----Modified monoclinic space groups - n glide                          
!     SG Pn - Subgroup 2                                                  
  211 SI(2,2)=.5                                                          
      SI(1,2)=.5                                                          
      IS(3,2)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 249                                                           
!-----Monoclinic space groups with Z=4                                    
!     SG P2/n - Subgroup 3                                                
  215 SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(2,4)=.5                                                          
      GO TO 217                                                           
!     SG P21/n - Subgroups 3 and 5                                        
  216 SI(1,3)=.5                                                          
      SI(2,3)=.5                                                          
      SI(3,3)=.5                                                          
      SI(1,4)=.5                                                          
      SI(2,4)=.5                                                          
      SI(3,4)=.5                                                          
  217 IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IS(3,2)=-1                                                          
      IS(1,3)=-1                                                          
      IS(2,3)=-1                                                          
      IS(3,4)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=.001                                                         
      GO TO 250                                                           
!     SG C-P21/n - Code SD                                                
  218 SI(3,2)=.5                                                          
!     SG C-P2/n - Code SE                                                 
  219 SI(1,2)=.5                                                          
      SI(2,2)=.5                                                          
      IS(1,2)=-1                                                          
      IS(2,2)=-1                                                          
      IHXNEG=1                                                            
      DSA(3)=0.001                                                        
      GO TO 249                                                           
  248 NOP=8                                                               
      GO TO 250                                                           
  249 NOP=2                                                               
!-----Should axes be interchanged to International Tables standard?       
  250 IF(NSGT.EQ.0)GO TO 253                                              
!     Output as determined by structure search routines                   
      J=1                                                                 
      K=2                                                                 
      L=3                                                                 
      GO TO 254                                                           
!-----Interchange axes to standard space group nomenclature               
  253 J=I1(MLOW)                                                          
      K=I2(MLOW)                                                          
      L=I3(MLOW)                                                          
  254 AA(1)=A(J)                                                          
      AA(2)=A(K)                                                          
      AA(3)=A(L)                                                          
      CC(1)=C(J)                                                          
      CC(2)=C(K)                                                          
      CC(3)=C(L)                                                          
!-----Interchange angle refinement parameters                             
      DSA4=DSA(J)                                                         
      DSA5=DSA(K)                                                         
      DSA6=DSA(L)                                                         
!     Tranform atomic coordinates and calculate molecular weight          
      WMOL=0.                                                             
      DO 255 I=1,NMOD                                                     
      W(1,I)=V(J,I)                                                       
      W(2,I)=V(K,I)                                                       
      W(3,I)=V(L,I)                                                       
      M=IAA(I)/120+1                                                      
      WMOL=WMOL+ATWT(M)                                                   
  255 CONTINUE                                                            
      CALL PAGE(3,3,0)                                                    
      WRITE(NQ,256)WMOL,NCN,NTA                                           
  256 FORMAT(1H0,'MOLECULAR WEIGHT =',F8.3,', NCN =',I3,', NTA =',I2)     
      D=WMOL/(HT1(8)*0.602252)                                            
      WRITE(NQ,257)D                                                      
  257 FORMAT('          DENSITY =',F8.3)                                  
!----Symmetry operations                                                  
      DO 258 I=1,NOP                                                      
      JJ=IS(J,I)                                                          
      KK=IS(K,I)                                                          
      LL=IS(L,I)                                                          
      D=SI(J,I)                                                           
      E=SI(K,I)                                                           
      F=SI(L,I)                                                           
      IS(1,I)=JJ                                                          
      IS(2,I)=KK                                                          
      IS(3,I)=LL                                                          
      SI(1,I)=D                                                           
      SI(2,I)=E                                                           
  258 SI(3,I)=F                                                           
!-----Put code and rotation angles in HEAD line to be written             
      IA1=HT1(10)                                                         
      IA2=HT1(11)                                                         
      IA3=HT1(12)                                                         
      WRITE(HHEAD,259)HEAD,CHT1(9),IA1,IA2,IA3                            
  259 FORMAT(A46,A2,3I4)                                                  
!-----Write input lines for MOLCON or WMIN                                
  260 IF(NWM.LT.2)GO TO 300                                               
!     Write lines in MOLCON format                                        
!-----WRITE TITLE LINE                                                       
      WRITE(NS,262)HHEAD                                                     
  262 FORMAT(6HTITLE ,A60)                                                   
!-----Write MOLCON CELL line                                              
      WRITE(NS,263)(AA(I),I=1,3),(CC(I),I=1,3)                            
  263 FORMAT(4HCELL,9X,3F6.3,3F6.4)                                          
!-----WRITE SYMOP LINES                                                      
      DO 265 I=1,NOP                                                         
      J=2*IS(2,I)                                                            
      K=3*IS(3,I)                                                            
      WRITE(NS,264)IS(1,I),SI(1,I),J,SI(2,I),K,SI(3,I)                       
  264 FORMAT(5HSYMOP,8X,3(I6,F6.2))                                          
  265 CONTINUE                                                               
!-----Write RANGE  DIST line                                              
      WRITE(NS,266)                                                          
  266 FORMAT('RANGE  DIST',6X,'.5')                                          
!-----WRITE ATOM LINES                                                       
      DO 271 I=1,NMOD                                                        
      WRITE(NS,270)AT(I),AN(I),(W(J,I),J=1,3)                                
  270 FORMAT(7HATOMC  ,A2,A4,3F10.6)                                         
  271 CONTINUE                                                               
      IF(NWM.GT.3)GO TO 276                                               
!-----Write INPUT line - transfers control back to input stream           
      WRITE(NS,275)                                                          
  275 FORMAT(5HINPUT)                                                        
      IF(NWM.EQ.2)GO TO 1000                                                 
!-----Write END line instead of INPUT line                                
  276 WRITE(NS,278)                                                       
  278 FORMAT('END')                                                       
      IF(NWM.EQ.4)GO TO 1000                                              
!                                                                         
!     Write input lines for WMIN  program                                 
!-----WRITE WMIN LINE #1 & #2                                             
! 300 WRITE(NR,301)                                                       
! 301 FORMAT(8X,'0',T16,'100',T27,'5',T35,'16',T45,'1',T54,'0',           
!    1T63,'1',T72,'1'/8X,'1',T16,'500',T24,'3000',T35,'30',T44,'68',      
!    2T53,'20',T62,'40',T71,'20')                                         
!     First find out if both HB and OX atom types are present             
300   IHB = 0                                                             
      IOX = 0                                                             
!     for CCSP_1 (hydantoin) 2-3-04
      I64 = 0
      I65 = 0
!     for CCSP_3 (amides) 4-6-04
      I61 = 0
      I63 = 0
      DO 301 I=1,NMOD                                                     
         IF (AT(I) .EQ. 'HB') IHB = IHB + 1                               
         IF (AT(I) .EQ. 'OX') IOX = IOX + 1                               
!     OC for O of nitrocubane using spacial O---O                         
         IF (AT(I) .EQ. 'OC') IOX = IOX + 1                               
!     for CCSP_1 (hydantoin) 2-3-04
         IF (IDATOM(I) .EQ. 64) I64 = I64 + 1
         IF (IDATOM(I) .EQ. 65) I65 = I65 + 1
!     for CCSP_3 (amides) 4-8-04
         IF (IDATOM(I) .EQ. 61) I61 = I61 + 1
         IF (IDATOM(I) .EQ. 63) I63 = I63 + 1
301   CONTINUE                                                            
      IREP = 2                                                            
      IF (IHB .NE. 0 .OR. IOX .NE. 0) IREP = 4                             
      IF (I64 .NE. 0 .OR. I65 .NE. 0) IREP = 4  ! for CCSP_1 (hydantoin) 2-3-04
      IF (I61 .NE. 0 .OR. I63 .NE. 0) IREP = 4  ! for CCSP_3 (amides) 4-8-04
!     add the O...H bonding energy to calculate E                         
      WRITE(NR,302)                                                       
!     create three new WMIN input files with NEC & charges                
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,302)                                                      
       WRITE(29,302)                                                      
       WRITE(39,302)                                                      
      END IF                                                              
  302 FORMAT(8X,'0',T16,'200',T25,'200',T35,'16',T45,'1',T54,'0',   &       
     &T63,'1',T72,'1'/8X,'1',T15,'1999',T23,'19999',T35,'30',T44,'68',&      
     &T53,'20',T62,'40',T71,'20') ! 03-15-06                                        
!-----WRITE LINE #3....TITLE                                              
      WRITE(NR,305)HHEAD                                                     
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,305)HHEAD                                                 
       WRITE(29,305)HHEAD                                                 
       WRITE(39,305)HHEAD                                                 
      END IF                                                              
  305 FORMAT(4H   1,A60)                                                     
!-----WRITE LINE #4....set up for 8 cycles of mode 3 refinement           
      WRITE(NR,306)                                                       
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,306)                                                      
       WRITE(29,306)                                                      
       WRITE(39,306)                                                      
      END IF                                                              
  306 FORMAT('  3 99  0  0  1  0  0  5  0',T31,'0.0001',T38,'0.000001' &  ! replaced 8 w/99 11-4-04 
     &,T52,'0.1',T60,'0.01',T69,'10.0')                                   
!-----WRITE LINE #5.....ISB  title again                                  
!     DO 308 I=1,48                                                       
!     IF(CHT1(9).EQ.SYM(I))GO TO 312                                      
! 308 CONTINUE                                                            
      IF(ISB(MLOW).NE.0)GO TO 312                                         
  309 WRITE(NR,310)HHEAD                                                  
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,310)HHEAD                                                 
       WRITE(29,310)HHEAD                                                 
       WRITE(39,310)HHEAD                                                 
      END IF                                                              
  310 FORMAT(1X,'??',1X,A60)                                              
      GO TO 315                                                           
! 312 IF(ISB(I).LE.0)GO TO 309                                            
  312 WRITE(NR,313) ISB(MLOW),SYM(MLOW), HHEAD                            
      if (NPMIN == 1) write(59,1002) HHEAD(1:8), SYM(MLOW), IA1,IA2,IA3   ! 1-29-09 DU
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,313)ISB(MLOW),SYM(MLOW), HHEAD                            
       WRITE(29,313)ISB(MLOW),SYM(MLOW), HHEAD                            
       WRITE(39,313)ISB(MLOW),SYM(MLOW), HHEAD                            
      END IF                                                              
  313 FORMAT(1X,I2,1X,A2,A60)
1002  FORMAT(A8,2X,A2,2X,3I4)                                             ! 1-29-09 DU                                              
  315 ICENT=2                                                             
!-----Set ICENT equal to 1 if space group is P1bar                        
      IF((MLOW.EQ.2).OR.(MLOW.EQ.9))ICENT=1                               
!-----If fluorine atoms are present, set NKIND equal to 5                 
      NKIND=4                                                             
      DO 316 I=1,NMOD                                                     
      K=IAA(I)/120+1                                                      
      IF(K.NE.7)GO TO 316                                                 
      NKIND=5                                                             
      GO TO 320                                                           
  316 CONTINUE                                                            
!-----WRITE LINE #6 control card for NA,NKA,...                           
! 320 WRITE(NR,321) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG            
! 321 FORMAT(7I3,2X,'1',T29,'-9',T33,'0',T36,'1',T39,'2',T54,'1')         
  320 IF (IDCHG.EQ.4) THEN                                                
       WRITE(NR,321) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG, IREP     
       NNKIND=NKIND                                                        
       NKIND=NMOD                                                         
       WRITE(19,321) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG, IREP     
       WRITE(29,322) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG, IREP     
       WRITE(39,322) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG, IREP     
       GO TO  333                                                         
      END IF                                                               
      IF (IDCHG.GT.0) THEN                                                
        NKIND=NMOD                                                        
      END IF                                                              
      IF (IDCHG.GT.1) THEN                                                
       WRITE(NR,322) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG, IREP     
      ELSE                                                                
       WRITE(NR,321) NMODD,NKIND,NOP,ICENT,IHXNEG,IHYNEG,IHZNEG, IREP     
      END IF                                                              
  321 FORMAT(7I3,2X,'1',T29,'-9',T33,'0',T36,'1',T39,I1,T54,'1')           
  322 FORMAT(7I3,2X,'1',T29,'-9',T33,'1',T36,'1',T39,I1,T54,'1')           
!-----WRITE LINES # 7-9                                                   
!     WRITE(NR,325)                                                       
  333 WRITE(NR,325)                                                       
  325 FORMAT(5X,'0.25',T12,'2',T15,'2')                                   
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,325)                                                      
       WRITE(29,325)                                                      
       WRITE(39,325)                                                      
      END IF                                                              
!     Calculate the molecular weight                                      
      WRITE(NR,326) WMOL,NOP                                              
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,326) WMOL,NOP                                             
       WRITE(29,326) WMOL,NOP                                             
       WRITE(39,326) WMOL,NOP                                             
      END IF                                                              
  326 FORMAT (1X,F8.2,I3)                                                 
      WRITE(NR,330)                                                       
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,330)                                                      
       WRITE(29,330)                                                      
       WRITE(39,330)                                                      
      END IF                                                              
  330 FORMAT(6X,'0.4',T16,'0.5'/6X,'5.5',T16,'6.0')                       
!-----WRITE CELL DIMENSIONS LINE - WMIN LINE #10                          
      if(NPMIN == 1) then                             ! 1-29-09 DU
       WRITE(59,331)(AA(I),I=1,3),(CC(I),I=1,3)       ! 1-29-09 DU
       WRITE(59,1003) NKIND, NOP, WMOL                ! 1-29-09 DU
      end if                                          ! 1-29-09 DU 
1003  format(2I4,F10.4)                               ! 1-29-09 DU       
      WRITE(NR,331)(AA(I),I=1,3),(CC(I),I=1,3)                            
      WRITE(49,331)(AA(I),I=1,3),(CC(I),I=1,3)                           
      IF (IDCHG.EQ.4) THEN                                                
      WRITE(19,331)(AA(I),I=1,3),(CC(I),I=1,3)                            
      WRITE(29,331)(AA(I),I=1,3),(CC(I),I=1,3)                            
      WRITE(39,331)(AA(I),I=1,3),(CC(I),I=1,3)                            
      END IF                                                              
  331 FORMAT(3F9.5,3F9.5)                                                  
!-----WRITE WMIN LINE #11 - standard errors of cell                       
      WRITE(NR,332)DSA4,DSA5,DSA6                                         
      IF (IDCHG.EQ.4) THEN                                                
      WRITE(19,332)DSA4,DSA5,DSA6                                         
      WRITE(29,332)DSA4,DSA5,DSA6                                         
      WRITE(39,332)DSA4,DSA5,DSA6                                         
      END IF                                                              
  332 FORMAT(4X,'0.001',T14,'0.001',T23,'0.001',3F9.3)                    
!-----Write symmetry operation lines                                      
      DO 335 I=1,NOP 
      WRITE(59,334)(SI(J,I),IS(J,I),J=1,3)        ! 1-29-09 DU                                                       
      WRITE(NR,334)(SI(J,I),IS(J,I),J=1,3)                                   
      IF (IDCHG.EQ.4) THEN                                                
      WRITE(19,334)(SI(J,I),IS(J,I),J=1,3)                                  
      WRITE(29,334)(SI(J,I),IS(J,I),J=1,3)                                  
      WRITE(39,334)(SI(J,I),IS(J,I),J=1,3)                                  
      END IF                                                              
  334 FORMAT(10X,F5.2,I3,'  0  0',10X,F5.2,'  0',I3,'  0',10X,F5.2, &          
     &'  0  0',I3)                                                           
  335 CONTINUE                                                               
!-----WRITE ENERGY PARAMETERS - WMIN LINE 19                              
!-----If atom with charges CALL SETCHARGE                                 
      IF (IDCHG .GT. 0.AND. IDCHG.LT.4) THEN                              
        CALL SETCHARGE                                                    
        GO TO 345                                                         
      END IF                                                              
!     IF(NKIND.NE.5)GO TO 340                                             
      IF(NNKIND.NE.5)GO TO 340                                            
      WRITE(NR,336)                                                       
  336 FORMAT(3X,'C',T20,'19.01246',T29,'290.0000',T43,'1.8',T51,'12.0'/&   
     &3X,'N',T20,'20.72228',T29,'150.0000',T41,'1.805',T51,'14.0'/     &    
     &3X,'O',T20,'28.86570',T29,'390.0000',T42,'1.98',T51,'16.0'/      &     
     &3X,'H',T20,'6.70970',T29,' 43.4896',T42,'1.87',T52,'1.0'/        &    
     &3X,'F',T20,'17.60590',T29,'296.5007',T42,'2.08',T51,'19.0')         
!     GO TO 345                                                           
      GO TO 3451                                                          
  340 WRITE(NR,341)                                                       
  341 FORMAT(3X,'C',T20,'19.01246',T29,'290.0000',T43,'1.8',T51,'12.0'/&   
     &3X,'N',T20,'20.72228',T29,'150.0000',T41,'1.805',T51,'14.0'/     &   
     &3X,'O',T20,'28.86570',T29,'390.0000',T42,'1.98',T51,'16.0'/      &   
     &3X,'H',T20,'6.70970',T29,' 43.4896',T42,'1.87',T52,'1.0')           
 3451 IF (IDCHG.EQ.4) THEN                                                
!        CALL ALLWAYS       !*****remove later************************
      END IF                                                              
!     calculate center of molecule
  345 DO 365 J=1,NMOD
      VXX=VXX+W(1,J)
      VYY=VYY+W(2,J)
      VZZ=VZZ+W(3,J)
  365 CONTINUE
      VXX=VXX/NMOD
      VYY=VYY/NMOD
      VZZ=VZZ/NMOD
!     PRINT '(3E20.14)',VXX,VYY,VZZ
!  check "enantio" 5-9-05 
      CALL HAND(NMOD,AA,CC,X,W,AT,AN,LH)      
!  if LH=1 the molpak molecule is enantiomorph of input molecule
      IF(LH .EQ. 1) THEN
        DO I = 1, NMOD
         W(1,I) = -W(1,I)
         W(2,I) = -W(2,I)
         W(3,I) = -W(3,I) 
        END DO
        VXX = -VXX
        VYY = -VYY
        VZZ = -VZZ
      END IF 
!-----Write atom parameter lines                                          
      DO 360 I=1,NMOD
      READ(AT(I),'(A2)')ATN                                               
      K=IAA(I)/120+1                                                      
      READ(AT(I),350)P(1),P(2)                                               
      READ(AN(I),350)(P(J),J=3,6)                                            
  350 FORMAT(4A1)                                                            
      DO 352 J=1,4                                                           
      IF (P(J) .NE. ' ') GO TO 352                                             
      DO 351 L=J,5                                                           
  351 P(L)=P(L+1)                                                            
  352 CONTINUE                                                               
! 355 WRITE(NR,356)(P(J),J=1,4),IATYPE(K),(W(J,I),J=1,3)                     
  355 IF (IDCHG.GT.0.AND.IDCHG.LT.4) THEN                                  
!      WRITE(NR,357)ATN,I,I,(W(J,I),J=1,3)     
       WRITE(NR,356)(P(J),J=1,4),I,(W(J,I),J=1,3)      ! 2-3-04 CCSP DU
      ELSE                                                                
       WRITE(NR,356)(P(J),J=1,4),JATYPE(K),(W(J,I),J=1,3)                
      END IF                                                              
      IF (IDCHG.EQ.4) THEN                                                
       WRITE(19,357)ATN,IDATOM(I),I,(W(J,I),J=1,3)                        
       WRITE(29,357)ATN,IDATOM(I),I,(W(J,I),J=1,3)                        
       WRITE(39,357)ATN,IDATOM(I),I,(W(J,I),J=1,3)                       
      END IF                                                              
  356 FORMAT(4A1,I5,18X,3F9.5)                                             
  357 FORMAT(A2,I2,I5,18X,3F9.5)                                          
  360 CONTINUE                                                               
      WRITE(NR,366)VXX,VYY,VZZ                                            
      IF (IDCHG.EQ.4) THEN                                                
      WRITE(19,366)VXX,VYY,VZZ                                            
      WRITE(29,366)VXX,VYY,VZZ                                            
      WRITE(39,366)VXX,VYY,VZZ                                            
      END IF                                                              
  366 FORMAT('XTRA',23X,3F9.4)                                            
!-----WRITE ATOM RIGID-BODY IDENTIFIERS                                   
!     IRIGID=1                                                            
!     J=NMODD                                                             
! 370 IF(J.LE.24)GO TO 372                                                
!     J=J-24                                                              
!     WRITE(NR,371)(IRIGID,I=1,24)                                        
! 371 FORMAT(24I3)                                                        
!     GO TO 370                                                           
! 372 WRITE(NR,371)(IRIGID,I=1,J)                                         
! 374 WRITE(NR,375)NMODD                                                  
      IRIGID=1                                                            
      J=NMODD                                                             
  370 IF(J.LE.24)GO TO 372                                                
      J=J-24                                                              
      WRITE(NR,371)(IRIGID,I=1,24)                                        
      IF (IDCHG.EQ.4) THEN                                                
      WRITE(19,371)(IRIGID,I=1,24)                                        
      WRITE(29,371)(IRIGID,I=1,24)                                        
      WRITE(39,371)(IRIGID,I=1,24)                                        
      END IF                                                              
  371 FORMAT(24I3)                                                        
      GO TO 370                                                           
  372 WRITE(NR,371)(IRIGID,I=1,J)                                         
      IF (IDCHG.EQ.4) THEN                                                
      WRITE(19,371)(IRIGID,I=1,J)                                         
      WRITE(29,371)(IRIGID,I=1,J)                                         
      WRITE(39,371)(IRIGID,I=1,J)                                         
      END IF                                                              
  374 WRITE(NR,375)NMODD                                                  
      IF(IDCHG.EQ.4) THEN                                                 
       WRITE(19,375)NMODD                                                 
       WRITE(29,375)NMODD                                                 
       WRITE(39,375)NMODD                                                 
      END IF                                                               
  375 FORMAT(I3,T25,'.05',T34,'.01'/)                                     
!-----Write parameter refinement flags                                    
      GO TO                                                            &   
     &(378,380,382,384,384  ,382,380,386,386,386  ,386,388,390,392,388 &   
     &,382,384,386,390,390  ,388,398,400,400,402  ,402,402,406,408,386 &   
     &,386,388,410,404,402  ,408),NLOW                                          
  378 WRITE(NR,379)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,379)                                                     
        WRITE(29,379)                                                     
        WRITE(39,379)                                                     
      END IF                                                              
!     the flags for triclinic P-1(AA,AB),P1 (AA) & molecule with centre   
  379 FORMAT('111111111000'/'000000000000'/)                              
      GO TO 1000                                                          
  380 WRITE(NR,381)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,381)                                                     
        WRITE(29,381)                                                     
        WRITE(39,381)                                                     
      END IF                                                              
!     the flags for triclinic P-1(AA,AB),P1 (AA)                          
  381 FORMAT('111111111111'/'000000000000'/)                              
      GO TO 1000                                                          
  382 IF((NCTR.GE.4).AND.(NCTR.LE.6))GO TO 394                            
      WRITE(NR,383)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,383)                                                     
        WRITE(29,383)                                                     
        WRITE(39,383)                                                     
      END IF                                                              
!     for monoclinic Pm(AC),Pc(AD,AG),Cc(DA)                              
  383 FORMAT('111010111010'/'000000000000'/)                              
      GO TO 1000                                                          
  384 IF(NCTR.GE.7)GO TO 400                                              
      IF((NCTR.GE.1).AND.(NCTR.LE.3))GO TO 396                            
      WRITE(NR,385)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,385)                                                     
        WRITE(29,385)                                                     
        WRITE(39,385)                                                     
      END IF                                                              
!     for monoclinic P2(AE),P21(AF,AH),C2(DB)                             
  385 FORMAT('111010111101'/'000000000000'/)                              
      GO TO 1000                                                          
  386 WRITE(NR,387)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,387)                                                     
        WRITE(29,387)                                                     
        WRITE(39,387)                                                     
      END IF                                                              
!     for monoclinic P21/c(AI,AK,AM), P2/c(AJ,AL), P21/m(AN),P2/m(AO)     
!                    C2/c(DC,DD,DE)                                       
  387 FORMAT('111010111111'/'000000000000'/)                              
      GO TO 1000                                                          
  388 WRITE(NR,389)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,389)                                                     
        WRITE(29,389)                                                     
        WRITE(39,389)                                                     
      END IF                                                              
!     for orthorhombic P212121,Pbca,Pbcn                                   
  389 FORMAT('111000111111'/'000000000000'/)                              
      GO TO 1000                                                          
  390 WRITE(NR,391)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,391)                                                     
        WRITE(29,391)                                                     
        WRITE(39,391)                                                     
      END IF                                                              
!     for orthorhombic P21212,Pnn2,Pna21,Pba2,Pnc2,Pca21,Pba2,Pmn21,Pma2  
!                      Fdd2                                             
  391 FORMAT('111000111110'/'000000000000'/)                              
      GO TO 1000                                                          
  392 IF(MLOW.LE.30)GO TO 388                                             
      GO TO 390                                                           
  394 WRITE(NR,395)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,395)                                                     
        WRITE(29,395)                                                     
        WRITE(39,395)                                                     
      END IF                                                              
!     for monoclinic & molecule with 2 fold symmetry                      
  395 FORMAT('111010010010'/'000000000000'/)                              
      GO TO 1000                                                          
  396 WRITE(NR,397)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,397)                                                     
        WRITE(29,397)                                                     
        WRITE(39,397)                                                     
      END IF                                                              
!     for monoclinic & molecule with mirror symmetry                      
  397 FORMAT('111010010101'/'000000000000'/)                              
      GO TO 1000 
  398 DO I = 1, 25                                                         
      IF (CSYM(I).EQ.'AC') THEN      ! 6-9-07
      WRITE(NR,3999)
 3999 FORMAT('111001001110'/'000000000000'/)
      GO TO 1000
      END IF
      END DO                         ! 6-9-07
      WRITE(NR,399)
! 398 WRITE(NR,399)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,399)                                                     
        WRITE(29,399)                                                     
        WRITE(39,399)                                                     
      END IF                                                              
!     for pseudo triclinic & molecule with two-fold symmetry              
  399 FORMAT('111010010000'/'000000000000'/)                              
      GO TO 1000                                                          
  400 WRITE(NR,401)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,401)                                                     
        WRITE(29,401)                                                     
        WRITE(39,401)                                                     
      END IF                                                              
  401 FORMAT('111010111000'/'000000000000'/)                              
      GO TO 1000                                                          
  402 IF(MLOW.EQ.64)GO TO 404                                             
      WRITE(NR,403)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,403)                                                     
        WRITE(29,403)                                                     
        WRITE(39,403)                                                     
      END IF                                                              
!     for orthorhombic & molecule with two-fold on axis-3                 
  403 FORMAT('111000001001'/'000000000000'/)                              
      GO TO 1000                                                          
  404 WRITE(NR,405)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,405)                                                     
        WRITE(29,405)                                                     
        WRITE(39,405)                                                     
      END IF                                                              
!     for orthorhombic & molecule with two-fold on axis-2                 
  405 FORMAT('111000010010'/'000000000000'/)                              
      GO TO 1000                                                          
  406 WRITE(NR,407)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,407)                                                     
        WRITE(29,407)                                                     
        WRITE(39,407)                                                     
      END IF                                                              
!     for orthorhombic & molecule with mirror on axis-1                   
! 407 FORMAT('111000100100'/'000000000000'/)                              
  407 FORMAT('111000100011'/'000000000000'/)       ! 6-9-07
      GO TO 1000                                                          
  408 WRITE(NR,409)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,409)                                                     
        WRITE(29,409)                                                     
        WRITE(39,409)                                                     
      END IF                                                              
!     for orthorhombic & molecule with center of symmetry                 
  409 FORMAT('111000111000'/'000000000000'/)                              
      GO TO 1000                                                          
  410 WRITE(NR,411)                                                       
      IF(IDCHG.EQ.4) THEN                                                 
        WRITE(19,411)                                                     
        WRITE(29,411)                                                     
        WRITE(39,411)                                                     
      END IF                                                              
!     for orthorhombic & molecule with mirror on axis-2                   
  411 FORMAT('111000010101'/'000000000000'/)                              
      GO TO 1000                                                          
  900 CALL PAGE(2,2,0)                                                       
      WRITE(NQ,901)                                                          
  901 FORMAT('0WMIN CALLED WITH NO DATA')                                    
      GO TO 1000                                                          
  902 CALL PAGE(2,2,0)                                                    
      WRITE(NQ,903)                                                       
  903 FORMAT(1H0,'UNPROGRAMMED VALUE OF MLOW')                            
  999 KILL=-1 
!------------- write out PMIN.inp 1-29-09 -------------------------------------------
 1000 if (NPMIN == 1) then                  ! 1-29-09 DU
       CLOSE(UNIT=NR)                       ! 1-29-09 DU
       OPEN (UNIT= 9)                       ! 1-29-09 DU
       DO I = 1, 12+NOP                     ! 1-29-09 DU
        READ(9,*)                           ! 1-29-09 DU
       END DO                               ! 1-29-09 DU
       DO I = 1, NMOD                       ! 1-29-09 DU
        READ(9,1004) PC_CODE(I),PC_Q(I),PC_A(I),PC_B(I),PC_C(I),PC_AW(I)   ! 1-29-09 DU
       END DO
       DO I = 1, NMOD
        READ(9,1005) PM_AT(I),PM_X(I),PM_Y(I),PM_Z(I)
       END DO
       READ(9,*)
       READ(9,1007) (I_flag(I), I= 1,NMOD+1)
       DO I = 1,2
        READ(9,*)
       END DO
       READ(9,1008) PM_RB
       DO I = 1, NMOD
        WRITE(59,1006) PM_AT(I),PC_CODE(I),PM_X(I),PM_Y(I),PM_Z(I), & 
     &  PC_Q(I),PC_A(I),PC_B(I),PC_C(I),PC_AW(I)
       END DO
       WRITE(59,1008) PM_RB
       WRITE(59,1007) (I_flag(I), I= 1,NMOD)
      end if
1004  FORMAT(2X,I4,3X,5F9.4)
1005  FORMAT(A4,23X,3F9.5)
1006  FORMAT(5x,a4,2X,i5,8F10.5)
1007  FORMAT(24I3)
1008  FORMAT(A72)
! ----------------------- 1-29-09 ---------------------------------------
      RETURN 
      END SUBROUTINE WMIN 
