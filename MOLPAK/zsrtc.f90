      SUBROUTINE ZSRTC                                                      
!                                                                                                                          
!-----This subroutine finds structures of two-fold axis molecules        
!     in space groups P21212 and P2221. The coordination sphere contains 
!     6 I molecules in plane-1,2.  For SF and SG.                                      
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
!      DATA ALPHA/'SF','SG'/                                                 
 
      CHARACTER(2) :: ALPHA(2) = (/'SF','SG'/)
!
      INTEGER :: I, J, JE, K, KE, L, LCT, N
!
      REAL :: A(6), C(6)
      REAL :: D, D1, D2
      REAL :: H1A1, H1I, H1IH, H1IT, H2I, H2A1, H3A1, H3I
      REAL :: X3A1, X3I
      REAL :: VX

!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=63,64)                                        
    5 FORMAT(1X,'ZSRTC called',2I3)                                      
!-----Calculate I molecule separation along axis-1                       
   10 N=MARK2                                                            
      DO 22 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 22                                                 
      IF(IT(I+7).NE.0)GO TO 22                                           
      D2=TI(I+3)                                                            
      IF(ABS(D2).GT.CN(K))GO TO 22                                          
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 22                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   22 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      H1I=H1M                                                               
      CALL MINHI(H1I,0.,100.,H1M,NNE)                                      
      IF(KCT)904,23,906                                                  
   23 H1IT=2.*H1I                                                        
      H1IH=H1I/2.                                                        
!-----Calculate axis-2 separation of line of I molecules                 
      ER=ERMT                                                               
      H2I=H2M                                                            
      KE=NE                                                                 
      JE=1                                                                  
   31 N=MARK                                                                
      DO 38 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 38                                                 
      IF(IT(I+7).NE.0)GO TO 38                                           
      D1=TI(I+2)                                                            
   32 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 32                                               
   33 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 33                                               
      IF(D1.LT.-CN(K))GO TO 38                                              
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 33                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 33                                                              
   38 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(H2I,0.,100.,H2M,KE)                                       
      IF(KCT)904,40,906                                                 
   40 H2A1=H2I/2.                                                        
      IF(IPR.LT.3)GO TO 60                                                  
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,51)H2I,H1I                                                   
   51 FORMAT(6H H2I =,F7.3,8H,  H1I =,F7.3)                                 
!-----Go to the first designated space group in ZSRTC                    
   60 IF(ITR(63)-2)300,400,1000                                          
!-----Place layer of A1 molecules along both plus and minus axis-3       
!     Code SF - A1(1/2,Y,1/2)                                            
  300 LCT=1                                                                 
      H1A1=H1IH                                                          
  301 N=MARK2                                                            
      ER=ERMF                                                            
      MARK=N                                                             
      DO 319 I=1,NMOD                                                       
      DO 318 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=H2A1-W(2,J)-W(2,I)                                                  
      D1=H1A1+W(1,J)-W(1,I)                                                 
  311 D2=D2+H2I                                                          
      IF(D2.LT.CN(K))GO TO 311                                           
  312 D2=D2-H2I                                                          
      IF(D2.GT.CN(K))GO TO 312                                           
      IF(D2.LT.-CN(K))GO TO 318                                          
  313 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 313                                           
  314 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 314                                           
      IF(D1.LT.-CN(K))GO TO 312                                          
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 314                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=-W(3,J)-W(3,I)                                                
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 314                                                             
  318 CONTINUE                                                              
  319 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      H3A1=H3N                                                           
      CALL MINHI(H3A1,0.,100.,H3M,NNE)                                      
      X3A1=-H3N                                                          
      CALL MINHI(X3A1,0.,-100.,-H3M,NNE)                                    
      H3I=H3A1-X3A1                                                      
      VX=H1IH*H2I*H3I                                                    
      IF(VX.GT.V3)GO TO 360                                              
!-----Check derived I molecule distances                                 
      N=MARK                                                             
      DO 329 I=1,NMOD                                                       
      DO 328 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=W(2,J)-W(2,I)                                                      
      D1=W(1,J)-W(1,I)                                                      
  321 D2=D2+H2I                                                          
      IF(D2.LT.CN(K))GO TO 321                                           
  322 D2=D2-H2I                                                          
      IF(D2.GT.CN(K))GO TO 322                                           
      IF(D2.LT.-CN(K))GO TO 328                                          
  323 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 323                                           
  324 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 324                                           
      IF(D1.LT.-CN(K))GO TO 322                                          
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 324                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=-W(3,J)-W(3,I)                                                
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 324                                                             
  328 CONTINUE                                                              
  329 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      X3I=H3I                                                            
      CALL MINHI(X3I,H3I,100.,H3I,NNE)                                      
      IF(KCT)340,330,906                                                 
!-----I molecules too close, adjust H3I and H3P3                         
  330 H3A1=H3A1+(X3I-H3I)/2.                                             
      H3I=X3I                                                            
      VX=H1IH*H2I*H3I                                                    
      IF(VX.GT.V3)GO TO 360                                              
  340 V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 360                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=0.                                                             
      HT1(6)=0.                                                             
      HT1(7)=H3A1                                                           
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=25                                                               
      MLOW=LCT+62                                                        
  360 IF(IPR.LT.2)GO TO 370                                              
      A(1)=H1I                                                           
      A(2)=H2I                                                           
      A(3)=H3I                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=0.                                                            
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      A(6)=H3A1                                                          
      A(5)=H2A1                                                          
      A(4)=H1A1                                                          
!-----Then to fractional coordinates                                     
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=90.                                                           
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,361)ALPHA(LCT),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)      
  361 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3/1X,'  COSINES' &     
     &,3F7.4,'  ANGLES',3F7.2)                                           
      D=SQRT(H1A1**2+H2A1**2+H3A1**2)                                    
      WRITE(NQ,362)H1A1,H2A1,H3A1,(A(I),I=4,6),D                         
  362 FORMAT(1X,'A1 AT',3F8.3,';  ',3F8.4,'  D =',F8.2)                  
  370 IF(LCT.GE.2)GO TO 1000                                             
!-----Code SF completed, is Code SG to be calculated?                    
      IF(ITR(64).GE.2)GO TO 1000                                         
!-----Code SG - A1(0,Y,1/2)                                              
  400 LCT=2                                                              
      H1A1=0.                                                            
      GO TO 301                                                          
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT               
  903 FORMAT(26H0STORAGE EXCEEDED BY ZTRTA,5X,3(I9,I3),I9)                  
      GO TO 999                                                          
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)                                                      
  905 FORMAT(20H0PARAMETER TOO SMALL)                                    
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT(20H0PARAMETER TOO LARGE)                                    
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,921)                                                      
  921 FORMAT(21H0NO ATOM INTERACTIONS)                                   
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZSRTC      
