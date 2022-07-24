      SUBROUTINE ZSRTB                                                      
!                                                                                                                                
!-----This routine finds structures of centric molecules in space group  
!     P21/c.  The coordination sphere contains 6 I molecules in          
!     plane-1,3.  For SD and SE.
!
      USE molpakCommonMod
      
      IMPLICIT NONE
!                                                      
!      CHARACTER*80 BUF                                                       
!      CHARACTER*60 HEAD                                                  
!      CHARACTER*4 AN,CDIJ,CHT1                                           
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
!      DIMENSION IJD(500),CDIJ(500),CHT1(14)                              
!      EQUIVALENCE (DIJ(1),IJD(1)),(DIJ(1),CDIJ(1)),(HT1(1),CHT1(1))      
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
!
!      CHARACTER*2 ALPHA                                                  
!      DIMENSION A(6),C(6),ALPHA(2)                                       
!      DATA ALPHA/'SD','SE'/ 
!                                                
      CHARACTER(2) :: ALPHA(2) = (/'SD','SE'/)
!
      INTEGER :: I, ICT, JCT, JE, K, KE, L, LCT, MCT, N
!
      REAL :: A(6), C(6)
      REAL :: D, D1, DY1A3
      REAL :: H1A3, H2A3, H2I, H21I, H3A3, HL
      REAL :: VX
      REAL :: X1A3, X2A3, X2I
      REAL :: Y1A3, Y2A3
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)ITR(61)                                                 
    5 FORMAT(1X,'ZSRTB called',3I3)                                      
   10 ER=ERMF                                                            
!-----Approach of a grid of X molecules.                                 
!     Go to the first designated space group                             
      IF(ITR(61)-2)100,200,1000                                          
!     Space group Pc, P21 - Code SD                                      
!-----Set to use distances to atoms displaced to H3N - screw axis        
  100 ICT=-1                                                             

      LCT=1                                                              
      H3A3=H3N                                                           
  105 JCT=0                                                              
      X2A3=100.                                                          
      HL=H2M                                                             
      Y2A3=H2N                                                           
  110 KE=NE                                                              
      JE=1                                                               
      MCT=0                                                              
      MARK=MARK2                                                            
      DY1A3=H1M/NV                                                          
      Y1A3=-H1N                                                          
      MCT=0                                                                 
  115 N=MARK                                                                
      DO 119 I=MARK1,NEND1,NSTP1                                         
!-----Skip distances to atoms either offset or not offset by H3N.        
      K=ICT*IT(I)                                                           
      IF(K.LE.0)GO TO 119                                                
      D1=TI(I+5)+Y1A3                                                       
  117 D1=D1+H1M                                                             
      IF(D1.LE.CN(K))GO TO 117                                              
  118 D1=D1-H1M                                                             
      IF(D1.GT.CN(K))GO TO 118                                              
      IF(D1.LT.-CN(K))GO TO 119                                             
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 118                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+6)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 118                                                             
  119 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      CALL MINHI(Y2A3,0.,X2A3,HL,KE)                                        
      IF(KCT)904,151,152                                                 
  151 X1A3=Y1A3                                                          
      X2A3=Y2A3                                                          
      IF(MCT)156,153,156                                                 
  152 IF(MCT)156,153,155                                                 
  153 IF(JE.GE.NV)GO TO 154                                              
      JE=JE+1                                                            
      Y1A3=Y1A3+DY1A3                                                    
      GO TO 115                                                          
  154 MCT=1                                                              
      JE=2.*JE                                                           
      DY1A3=H1M/JE                                                       
      Y1A3=X1A3+DY1A3                                                    
      GO TO 115                                                          
  155 MCT=-1                                                             
      Y1A3=X1A3-DY1A3                                                    
      IF(Y1A3.LT.-H1N)Y1A3=Y1A3+H1M                              
      GO TO 115                                                          
  156 IF(JE.LT.KE)GO TO 154                                              
      IF(JE.GE.NNE)GO TO 158                                             
      KE=NNE                                                             
      GO TO 154                                                          
  158 CONTINUE                                                           
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
      H2I=2.*H2A3                                                  
      H21I=2.*H1A3                                                  
  160 VX=H1N*H2I*H3M                                                     
      IF(VX.GE.V3)GO TO 191                                                 
!-----Check derived distances of I molecules                             
      MARK=MARK2                                                         
      N=MARK                                                                
      DO 169 I=MARK1,NEND1,NSTP1                                         
      K=IT(I)                                                               
!-----Skip distances to offset atoms                                     
      IF(K.LE.0)GO TO 169                                                
      D1=TI(I+5)+H21I                                                       
  167 D1=D1+H1M                                                             
      IF(D1.LE.CN(K))GO TO 167                                              
  168 D1=D1-H1M                                                             
      IF(D1.GT.CN(K))GO TO 168                                              
      IF(D1.LT.-CN(K))GO TO 169                                             
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 168                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 168                                                             
  169 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 190                                             
      NEND=N-1                                                              
      X2I=H2M                                                            
      CALL MINHI(X2I,H2I,100.,H2M,KE)                                       
      IF(KCT)190,170,906                                                 
!-----I molecules too close - adjust axis-2 parameters                   
  170 H2A3=H2A3+(X2I-H2I)/2.                                             
      H2I=X2I                                                            
      VX=H1N*H2I*H3M                                                     
      IF(VX.GT.V3)GO TO 191                                              
  190 V3=VX                                                              
      DIJ(NP2)=V3                                                        
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 191                                             
      HT1(1)=H1M                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                         
      HT1(5)=0.                                                             
      HT1(6)=0.                                                             
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=24                                                              
      MLOW=60+LCT                                                        
  191 IF(IPR.LT.2)GO TO 197                                              
      A(1)=H1M                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3M                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      A(6)=H3A3                                                          
      A(5)=H2A3/SQRT(1.-C(3)**2)                                         
      A(4)=H1A3-A(5)*C(3)                                                
!-----Then to fractional coordinates                                     
      D=SQRT(H1A3**2+H2A3**2+H3A3**2)                                    
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,195)ALPHA(LCT),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)      
  195 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3/1X,'  COSINES'  &    
     &,3F7.4,'  ANGLES',3F7.2)                                           
      WRITE(NQ,196)H1A3,H2A3,H3A3,(A(I),I=4,6),D                         
  196 FORMAT(1X,'A3 AT',3F7.3,';  ',3F7.4,';  D =',F6.2)               
  197 IF(LCT.GE.2)GO TO 1000                                             
!-----Is space group Pc, P2 to be calculated - Code SE                   
      IF(ITR(62).GT.1)GO TO 1000                                         
!     Set to pick up atoms not displaced along axis-3                    
  200 ICT=1                                                              
      LCT=2                                                              
      H3A3=0.                                                            
      GO TO 105                                                          
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTB,5X,3(I9,I3),I9)                  
      GO TO 999                                                             
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)                                                      
  905 FORMAT('0PARAMETER TOO SMALL')                                     
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT('0PARAMETER TOO LARGE')                                     
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)                                                         
  921 FORMAT(16H0NO INTERACTIONS)                                           
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                             
      END SUBROUTINE ZSRTB
