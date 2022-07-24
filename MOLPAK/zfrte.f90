      SUBROUTINE ZFRTE                                                      
!                                                                                                                                                                                                
!-----This subroutine finds structures in the space group P21/c in which 
!     the coordination sphere contains 6 I molecules in plane-1,3, 4 C,  
!     and 4 A3 molecules. For FC and FD.
!                                               
      USE molpakCommonMod
!
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
!      CHARACTER*2 ALPHA(2),SYM(3)                                        
!      DIMENSION A(12),C(6),SM(3)                                         
!      DATA ALPHA/'FC','FD'/                                                 
!      DATA SYM/'C ','A3','P3'/                                           
!
      CHARACTER(2) :: ALPHA(2) = (/'FC','FD'/)
      CHARACTER(2) :: SYM(3) = (/'C ','A3','P3'/)
!
      INTEGER :: I, ICT, J, JCT, JE, K, KE, L, LCT, LE
      INTEGER :: MCT, N, NCT, NERR
!
      REAL :: A(12),C(6)
      REAL :: D, D1, D2, D3, DX1P3, DX3P3, DY1A3
      REAL :: H1A3, H1C, H1P3, HH1P3 
      REAL :: H2A3, H2A3T, H2C, H2I,  H21I, H2P3
      REAL :: H3A3, H3C, H3P3, HH3P3, HL
      REAL :: SM(3), SINE3, VX
      REAL :: X1C, X1A3, X1P3, X2C, X2A3, X21I, X2P3, X3P3
      REAL :: Y1A3, Y2A3, Y2P3 
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)                                                        
    5 FORMAT(1X,'ZFRTE called')                                          
   10 ER=ERMF                                                            
!     Approach of grid of A molecules from negative axis-2 direction     
!-----Space group P21/c or P2/c                                          
      IF(ITR(74)-2)100,105,1000                                          
!-----Set signal to pick up distances to atoms offset by H3N; that       
!     is, a screw axis.                                                  
  100 ICT=-1                                                             
      LCT=1                                                              
      H3A3=H3N                                                           
      GO TO 110                                                          
!-----Set signal to pick up distances to atoms NOT offset by H3N; that   
!     is, a two-fold axis.                                               
  105 ICT=1                                                              
      LCT=2                                                              
      H3A3=0.                                                            
  110 X2A3=-100.                                                         
      HL=-H2M                                                            
      Y2A3=-H2N                                                          
      KE=NE                                                              
      JE=1                                                               
      MCT=0                                                              
      MARK=MARK2                                                            
      DY1A3=H1M/NV                                                          
      Y1A3=-H1N                                                          
      MCT=0                                                                 
  115 N=MARK                                                                
      DO 119 I=MARK1,NEND1,NSTP1                                         
      K=ICT*IT(I)                                                           
!-----Skip distances to either offset or non-offset atoms                
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
      IF(KCT)1904,151,152                                                
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
      H2A3T=2.*H2A3
!     Calculate standoff of grid of P3 and C molecules in along axis-2   
  300 KE=NE                                                              
      MARK=MARK2                                                         
      JE=1                                                                  
      LE=1                                                                  
      DX1P3=H1M/NV                                                          
      DX3P3=H3M/NV                                                       
      X3P3=-H3N                                                             
      X1P3=-H1N                                                             
      NCT=1                                                                 
      MCT=1                                                                 
      H2P3=100.                                                             
      X2P3=H2N                                                              
      ER=ERMF                                                               
  315 N=MARK                                                                
      IF(IPR.LT.4)GO TO 320                                              
      WRITE(NQ,316)X2P3,X1P3,X3P3,NCT,MCT,JE,LE                          
  316 FORMAT(1X,3F7.3,4I3)                                               
  320 CONTINUE                                                           
      X1C=X1P3+H1A3
      X21I=2.*X1P3
!-----Collect I molecule distances for double standoff				   
	DO 329 I=1,NMOD
	DO 328 J=1,NMOD
	K=IA(I)+IAA(J)
	L=K+10
	D1=X21I+W(1,J)-W(1,I)
	D2=W(2,J)-W(2,I)
	D3=W(3,J)-W(3,I)
  321 D3=D3+H3M
      IF(D3.LE.CN(K))GO TO 321
  322 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 322
	IF(D3.LT.-CN(K))GO TO 328
  323 D1=D1+H1M
      IF(D1.LE.CN(K))GO TO 323
  324 D1=D1-H1M
      IF(D1.GT.CN(K))GO TO 324
	IF(D1.LT.-CN(K))GO TO 322
	D=D1**2+D3**2
	IF(D.GT.CN(L))GO TO 324
	IT(N)=K
	TI(N+1)=D
	TI(N+2)=D2
	N=N+NSTP
	IF(N.GT.LIMIT)GO TO 902
	GO TO 324
  328 CONTINUE
  329 CONTINUE
      NARK=N
      DO 339 I=1,NMOD                                                       
      DO 338 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      JCT=0
      D1=X1P3+W(1,J)-W(1,I)                                                 
      D2=W(2,J)-W(2,I)
      D3=X3P3-W(3,J)-W(3,I)                                                 
  331 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 331                                              
  332 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 332                                              
      IF(D3.LT.-CN(K))GO TO 335                                             
  333 D1=D1+H1M                                                             
      IF(D1.LE.CN(K))GO TO 333                                              
  334 D1=D1-H1M                                                             
      IF(D1.GT.CN(K))GO TO 334                                              
      IF(D1.LT.-CN(K))GO TO 332                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 334                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2  
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 334                                                             
!-----Have C molecules been included                                     
  335 IF(JCT.GT.0)GO TO 338
      JCT=1
	D1=X1C-W(1,J)-W(1,I)
	D2=H2A3-W(2,J)-W(2,I)
	D3=D3+H3A3
	GO TO 331
  338 CONTINUE                                                              
  339 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      CALL MINHI(X2P3,0.,H2P3,H2M,KE)                                       
      NARK=0
      IF(KCT)2904,342,341                                                
  341 IF(MCT.LE.1)GO TO 344                                              
      J=MCT-NCT                                                          
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 357                                 
      GO TO 343                                                          
  342 H3P3=X3P3                                                          
      H1P3=X1P3                                                          
      H2P3=X2P3                                                          
      NCT=MCT                                                            
  343 GO TO(344,349,350,351,352,353,354,355,348),MCT                     
  344 IF(LE.GE.NV)GO TO 345                                              
      LE=LE+1                                                            
      X1P3=X1P3+DX1P3                                                    
      GO TO 315                                                          
  345 IF(JE.GE.NV)GO TO 346                                              
      JE=JE+1                                                            
      X3P3=X3P3+DX3P3                                                    
      LE=1                                                               
      X1P3=-H1N                                                          
      GO TO 315                                                          
  346 NCT=6                                                              
  347 JE=2.*JE                                                           
      DX3P3=H3M/JE                                                       
      DX1P3=H1M/JE                                                       
      HH3P3=H3P3                                                         
      HH1P3=H1P3                                                         
      NCT=NCT+1                                                          
      IF(NCT.GT.9)NCT=NCT-8                                              
      GO TO(908,350,351,352,353,354,355,348,349),NCT                     
  348 MCT=2                                                              
      X3P3=HH3P3+DX3P3                                                   
      X1P3=HH1P3                                                         
      GO TO 315                                                          
  349 MCT=3                                                              
      X3P3=HH3P3+DX3P3                                                   
      X1P3=HH1P3+DX1P3                                                   
      GO TO 315                                                          
  350 MCT=4                                                              
      X3P3=HH3P3                                                         
      X1P3=HH1P3+DX1P3                                                   
      GO TO 315                                                          
  351 MCT=5                                                              
      X3P3=HH3P3-DX3P3                                                   
      IF(X3P3.LT.-H3N)X3P3=X3P3+H3M                                      
      X1P3=HH1P3+DX1P3                                                   
      GO TO 315                                                          
  352 MCT=6                                                              
      X3P3=HH3P3-DX3P3                                                   
      IF(X3P3.LT.-H3N)X3P3=X3P3+H3M                                      
      X1P3=HH1P3                                                         
      GO TO 315                                                          
  353 MCT=7                                                              
      X3P3=HH3P3-DX3P3                                                   
      IF(X3P3.LT.-H3N)X3P3=X3P3+H3M                                      
      X1P3=HH1P3-DX1P3                                                   
      IF(X1P3.LT.-H1N)X1P3=X1P3+H1M                                      
      GO TO 315                                                          
  354 MCT=8                                                              
      X3P3=HH3P3                                                         
      X1P3=HH1P3-DX1P3                                                   
      IF(X1P3.LT.-H1N)X1P3=X1P3+H1M                                      
      GO TO 315                                                          
  355 MCT=9                                                              
      X3P3=HH3P3+DX3P3                                                   
      X1P3=HH1P3-DX1P3                                                   
      IF(X1P3.LT.-H1N)X1P3=X1P3+H1M                                      
      GO TO 315                                                          
  357 IF(JE.LT.KE)GO TO 347                                              
      IF(JE.GE.NNE)GO TO 360                                             
      KE=NNE                                                             
      GO TO 347                                                          
  360 H21I=2.*H1P3                                                       
      H2I=2.*H2P3                                                        
      H1C=H1A3+H1P3
	H2C=H2A3+H2P3
	H3C=H3A3+H3P3
      VX=H1N*H2I*H3N                                                     
      IF(VX.GT.V3)GO TO 400                                                 
!-----Check P and C molecules on negative side                           
      X1P3=-H21I+H1P3
      X1C=-H1P3+H1A3
      X2C=-H2I+H2A3
      N=MARK                                                             
      DO 369 I=1,NMOD                                                       
      DO 368 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      JCT=0
      D3=H3P3-W(3,J)-W(3,I)                                                
      D2=W(2,J)-W(2,I)
      D1=-H1P3+W(1,J)-W(1,I)                                                
  361 D3=D3+H3M                                                          
      IF(D3.LT.CN(K))GO TO 361                                           
  362 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 362                                           
      IF(D3.LT.-CN(K))GO TO 365                                          
  363 D1=D1+H1M                                                          
      IF(D1.LT.CN(K))GO TO 363                                           
  364 D1=D1-H1M                                                          
      IF(D1.GT.CN(K))GO TO 364                                           
      IF(D1.LT.-CN(K))GO TO 362                                          
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 364                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 364                                                             
  365 IF(JCT.GT.0)GO TO 368
      D1=X1C-W(1,J)-W(1,I)
	D2=H2A3-W(2,J)-W(2,I)
	D3=H3C-W(3,J)-W(3,I)
	JCT=1
	GO TO 361
  368 CONTINUE                                                              
  369 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      X2P3=-H2P3 
      Y2P3=-H2P3                                                         
      CALL MINHI(Y2P3,X2P3,-100.,-H2M,NNE)                                  
!     WRITE(NQ,367)X2P3,Y2P3,X1C,H2A3,H3C 
! 367 FORMAT(' TEST',5F8.3)
      IF(KCT)375,370,906                                                 
!-----Adjust the P and C molecule positions                              
  370 CONTINUE
!     WRITE(NQ,371)X2A3,H2A3
! 371 FORMAT(' X2A3, H2A3 =',2F8.3)
      H2P3=-Y2P3
      H2C=H2P3-H2A3
      H2I=2.*H2P3
      VX=H1N*H2I*H3N                                                     
  375 IF(VX.GT.V3)GO TO 400                                              
  376 V3=VX                                                              
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 400                                             
      HT1(1)=H1M                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                         
      HT1(5)=H1C                                                            
      HT1(6)=H2C                                                            
      HT1(7)=H3C                                                            
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=31                                                               
      MLOW=73+LCT                                                        
  400 IF(IPR.LT.2)GO TO 500                                              
      A(1)=H1M                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3M                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      SINE3=SQRT(1.-C(3)**2)                                             
      A(6)=H3C                                                           
      A(5)=H2C/SINE3                                                     
      A(4)=H1C-A(5)*C(3)                                                 
      A(9)=H3P3                                                          
      A(8)=H2P3/SINE3                                                    
      A(7)=H1P3-A(8)*C(3)                                                
      A(12)=H3A3                                                         
      A(11)=H2A3/SINE3                                                   
      A(10)=H1A3-A(11)*C(3)                                              
!-----Then to fractional coordinates                                     
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      A(7)=A(7)/A(1)                                                     
      A(8)=A(8)/A(2)                                                     
      A(9)=A(9)/A(3)                                                     
      A(10)=A(10)/A(1)                                                   
      A(11)=A(11)/A(2)                                                   
      A(12)=A(12)/A(3)                                                   
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,410)ALPHA,A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)           
  410 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3       & 
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
  411 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)                                 
      SM(1)=A(4)-A(4)                                                    
      SM(2)=A(5)-A(5)                                                    
      SM(3)=A(6)-A(6)                                                    
      WRITE(NQ,411)SYM(1),(A(I),I=4,6),SM                                
      SM(1)=A(10)-A(4)                                                   
      SM(2)=A(11)-A(5)                                                   
      SM(3)=A(12)                                                        
      WRITE(NQ,411)SYM(2),(A(I),I=10,12),SM                              
      SM(1)=A(7)                                                         
      SM(2)=A(8)                                                         
      SM(3)=A(9)-A(6)                                                    
      WRITE(NQ,411)SYM(3),(A(I),I=7,9),SM                                
  500 IF(LCT.GE.2)GO TO 1000                                                
      IF(ITR(75)-1)105,105,1000                                          
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZFRTE,5X,3(I9,I3),I9)                  
      GO TO 999                                                             
 1904 NERR=151
      GO TO 904
 2904 NERR=341
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)NERR                                                  
  905 FORMAT('0PARAMETER TOO SMALL NEAR', I4)                            
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT('0PARAMETER TOO LARGE')                                     
      GO TO 999                                                          
  908 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,909)                                                      
  909 FORMAT('0DTWO failure')                                            
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)                                                         
  921 FORMAT(16H0NO INTERACTIONS)                                           
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                             
      END SUBROUTINE ZFRTE  
