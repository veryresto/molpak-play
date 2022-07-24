      SUBROUTINE ZSRTF                                                      
!                                                                     
!     This subroutine finds structures in space groups Pmn21 and Pma2    
!      when the molecule contains a mirror plane.  The coordination      
!      sphere contains 6 I molecules in plane-1,3. For SL and SM.                      
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
!      DIMENSION ALPHA(2)                                                    
!      DIMENSION A(6),C(6)                                                
!      DATA ALPHA/'SL','SM'/                                                 
!
      CHARACTER(2) :: ALPHA (2) = (/'SL','SM'/)
!
      INTEGER :: I, J, K, KE, L, LCT, N
!
      REAL :: A(6), C(6)
      REAL :: D, D1, D3
      REAL :: H1P2, H2I, H2P2, H3P2, HL
      REAL :: VX, X2I, X2P2, Y2P2
!
      NARK=0                                                             
      MARK=MARK2                                                         
      KE=NNE                                                             
      ER=ERMF                                                            
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=69,70)                                        
    5 FORMAT(1X,'ZSRTF called',2I3)                                      
!-----Go to first designated ZSRTF space group                           
   10 IF(ITR(69)-2)100,380,1000                                          
!-----Code SL - P2(1/2,1/2,z)                                            
  100 LCT=1                                                              
      H1P2=H1N                                                           
      H3P2=H3N                                                           
!     Approach of grids of P2 molecules from both positive and negative  
!     axis-2 directions                                                  
!-----Set up for positive direction                                      
  105 X2P2=100.                                                          
      HL=H2M                                                             
      Y2P2=H2N                                                           
  315 N=MARK                                                                
      DO 339 I=1,NMOD                                                       
      DO 338 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=H1P2+W(1,J)-W(1,I)                                                 
      D3=H3P2+W(3,J)-W(3,I)                                                 
  331 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 331                                              
  332 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 332                                              
      IF(D3.LT.-CN(K))GO TO 338                                             
  333 D1=D1+H1M                                                             
      IF(D1.LE.CN(K))GO TO 333                                              
  334 D1=D1-H1M                                                             
      IF(D1.GT.CN(K))GO TO 334                                              
      IF(D1.LT.-CN(K))GO TO 332                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 334                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=-W(2,J)-W(2,I)                                                
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 334                                                             
  338 CONTINUE                                                              
  339 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      CALL MINHI(Y2P2,0.,X2P2,HL,KE)                                        
      IF(KCT)904,340,906                                                 
  340 IF(HL.LT.0.)GO TO 350                                              
      H2P2=Y2P2                                                          
!-----Setup for negative axis-2 direction                                
      HL=-H2M                                                            
      X2P2=-100.                                                         
      GO TO 315                                                          
  350 H2I=H2P2-Y2P2                                                      
!-----Check I molecule distances                                         
      N=MARK                                                             
      DO 359 I=1,NMOD                                                       
      DO 358 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=W(1,J)-W(1,I)                                                      
      D3=W(3,J)-W(3,I)                                                      
  351 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 351                                              
  352 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 352                                              
      IF(D3.LT.-CN(K))GO TO 358                                             
  353 D1=D1+H1M                                                             
      IF(D1.LE.CN(K))GO TO 353                                              
  354 D1=D1-H1M                                                             
      IF(D1.GT.CN(K))GO TO 354                                              
      IF(D1.LT.-CN(K))GO TO 352                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 354                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)-W(2,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 354                                                             
  358 CONTINUE                                                              
  359 CONTINUE                                                              
      NEND=N-1                                                           
      CALL MINHI(X2I,0.,H2I,H2M,KE)                                      
      IF(KCT)904,361,360                                                 
!-----I molecule distances are too small                                 
  360 H2P2=H2P2+(X2I-H2I)/2.                                             
      H2I=X2I                                                            
  361 VX=H1N*H2I*H3M                                                     
      IF(V3-VX.LT..01)GO TO 362                                             
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 362                                             
      HT1(1)=H1M                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                         
      HT1(5)=0.                                                             
      HT1(6)=H2P2                                                           
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=28                                                               
      MLOW=LCT+68                                                        
  362 IF(IPR.LT.2)GO TO 370                                              
      A(1)=H1M                                                           
      A(2)=H2I                                                           
      A(3)=H3M                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=0.                                                            
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      A(6)=H3P2                                                          
      A(5)=H2P2                                                          
      A(4)=H1P2                                                          
!-----Then to fractional coordinates                                     
      D=SQRT(H1P2**2+H2P2**2+H3P2**2)                                    
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=90.                                                           
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,365)ALPHA(LCT),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)      
  365 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3/1X,'  COSINES'  &    
     &,3F7.4,'  ANGLES',3F7.2)                                           
      WRITE(NQ,366)H1P2,H2P2,H3P2,(A(I),I=4,6),D                         
  366 FORMAT(1X,'P2 AT',3F7.3,';  ',3F7.4,';  D =',F6.2)                  
  370 IF(LCT.GE.2)GO TO 1000                                             
!-----Is Code SM to be calculated                                        
      IF(ITR(70).GE.2)GO TO 1000                                         
!-----Code SM - P2(0,1/2,z)                                              
  380 LCT=2                                                              
      H1P2=0.                                                            
      H3P2=H3N                                                           
      GO TO 105                                                             
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZTRTB,5X,3(I9,I3),I9)                  
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
      END SUBROUTINE ZSRTF

