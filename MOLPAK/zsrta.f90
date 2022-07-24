      SUBROUTINE ZSRTA 
!                                                                                                                        
!-----This subroutine finds structures in space group P2 containing     
!      two-fold axis molecules, Pm with mirror plane molecules, and     
!      P21/m or P21/c with centric molecules.  The coordination spheres 
!      produced contain 6 I molecules in plane-1,2.  For SA, SB and SC.                   
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
!      CHARACTER*2 ALPHA                                                  
!      DIMENSION ALPHA(3)                                                    
!      DIMENSION A(6),C(6)                                                
!      DATA ALPHA/'SA','SB','SC'/                                            
!
      CHARACTER(2) :: ALPHA(3) = (/'SA','SB','SC'/)
!
      INTEGER :: I, J, JE, K, KE, L, LCT, MCT, N
!
      REAL :: A(6), C(6)
      REAL :: D, D1, D2
      REAL :: DX21I
      REAL :: H1I, H1IH, H1P3, H1IT 
      REAL :: H2I, H21I
      REAL :: H3I, H3P3
      REAL :: VX
      REAL :: X2I, X21I
      REAL :: X3I, X3P3
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=1,5)                                          
    5 FORMAT(1X,'ZSRTA called',5I3)                                      
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
      H2I=100.                                                           
      DX21I=H1I/NV                                                       
      X21I=-H1IH                                                         
      X2I=H2M                                                            
      KE=NE                                                                 
      JE=1                                                                  
      MCT=0                                                                 
   31 N=MARK                                                                
      DO 38 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 38                                                 
      IF(IT(I+7).NE.0)GO TO 38                                           
      D1=TI(I+2)+X21I                                                       
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
      CALL MINHI(X2I,0.,H2I,H2M,KE)                                         
      IF(KCT)904,41,42                                                   
   41 H21I=X21I                                                          
      H2I=X2I                                                            
      IF(MCT)46,43,46                                                    
   42 IF(MCT)46,43,45                                                    
   43 IF(JE.GE.NV)GO TO 44                                               
      JE=JE+1                                                            
      X21I=X21I+DX21I                                                    
      GO TO 31                                                           
   44 MCT=1                                                              
      JE=2.*JE                                                           
      DX21I=H1I/JE                                                       
      X21I=H21I+DX21I                                                    
      GO TO 31                                                           
   45 MCT=-1                                                             
      X21I=H21I-DX21I                                                    
      IF(X21I.LT.H1IH)X21I=X21I+H1I      
      GO TO 31                                                           
   46 IF(JE.LT.KE)GO TO 44                                               
      IF(JE.GE.NNE)GO TO 48                                              
      KE=NNE                                                             
      GO TO 44                                                           
   48 CONTINUE                                                           
      IF(IPR.LT.3)GO TO 60                                                  
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,51)H2I,H1I,H21I                                              
   51 FORMAT(6H H2I =,F7.3,8H,  H1I =,F7.3,8H, H21I =,F7.3)                 
!-----Go to the first designated space group in ZTRTA                    
   60 GO TO(100,300,400,1000),ITR(58)                                    
!-----Calculate axis-3 separation of layer of I molecules -             
!     Effective space group P1 with alpha and gamma = 0 - Code SA       
  100 ER=ERMF                                                               
      MARK=MARK2                                                            
      H3I=H3M                                                            
  113 KE=NE                                                                 
  116 N=MARK                                                                
      DO 129 I=1,NMOD                                                       
      DO 128 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=W(1,J)-W(1,I)                                                      
      D2=W(2,J)-W(2,I)                                                      
  122 D2=D2+H2I                                                             
      D1=D1+H21I                                                            
      IF(D2.LE.CN(K))GO TO 122                                              
  123 D2=D2-H2I                                                             
      D1=D1-H21I                                                            
      IF(D2.GT.CN(K))GO TO 123                                              
      IF(D2.LT.-CN(K))GO TO 128                                             
  124 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 124                                              
  125 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 125                                              
      IF(D1.LT.-CN(K))GO TO 123                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 125                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(3,J)-W(3,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 125                                                             
  128 CONTINUE                                                              
  129 CONTINUE                                                              
      NEND=N-1                                                              
      CALL MINHI(H3I,H3M,100.,H3M,KE)                                       
      IF(KCT)142,142,906                
!-----Calculate the cell volume per molecule                            
  142 VX=H1I*H2I*H3I                                                        
      IF(V3-VX.LT..01)GO TO 161                                             
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(1)                                                 
      IF(V3.GE.HT1(8))GO TO 161                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=0.                                                             
      HT1(6)=0.                                                             
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(1)                                                      
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=22                                                               
      MLOW=58                                                            
  161 IF(IPR.LT.2)GO TO 170                                              
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
!     Cosines of cell angles                                            
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
      C(4)=0.                                                            
      C(5)=0.                                                            
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,162)ALPHA(1),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)        
  162 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3/1X,'  COSINES' &     
     &,3F7.4,'  ANGLES',3F7.2)                                           
!-----Go to the next designated space group in ZTRTA                    
  170 GO TO(300,400,1000),ITR(59)                                        
!-----Place layer of X molecules at 1/2 on axis-3 and either 0 or 1/2   
!     on axis-1.                                                        
!     Effective space groups Pm, P21 and/or Pc, P21 - Codes SB and SC   
!-----Space group Pm, P21 - Code SB                                     
  300 LCT=2                                                                 
      H1P3=0.                                                            
  301 N=MARK2                                                            
      ER=ERMF                                                            
      MARK=N                                                             
      DO 319 I=1,NMOD                                                       
      DO 318 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=W(2,J)-W(2,I)                                                      
      D1=W(1,J)+H1P3-W(1,I)                                                 
  311 D2=D2+H2I                                                          
      D1=D1+H21I                                                         
      IF(D2.LT.CN(K))GO TO 311                                           
  312 D2=D2-H2I                                                          
      D1=D1-H21I                                                         
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
      H3P3=H3N                                                           
      CALL MINHI(H3P3,0.,100.,H3M,NNE)                                      
      IF(KCT)904,320,906                                                 
  320 X3P3=-H3N                                                          
      CALL MINHI(X3P3,0.,-100.,-H3M,NNE)                                    
      IF(-X3P3.GT.H3P3)H3P3=-X3P3                                        
      H3I=2.*H3P3                                                        
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
      D1=D1+H21I                                                         
      IF(D2.LT.CN(K))GO TO 321                                           
  322 D2=D2-H2I                                                          
      D1=D1-H21I                                                         
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
  330 H3P3=H3P3+(X3I-H3I)/2.                                             
      H3I=X3I                                                            
      VX=H1IH*H2I*H3I                                                    
      IF(VX.GT.V3)GO TO 360                                              
  340 V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 360                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
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
      NLOW=23                                                               
      MLOW=LCT+57                                                        
  360 IF(IPR.LT.2)GO TO 370                                              
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
!     Cosines of cell angles                                            
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates 
!-----Transform first to Angstrom coordinates                           
      A(6)=H3P3                                                          
      A(5)=0.                                                            
      A(4)=H1P3-A(5)*C(3)                                                
!-----Then to fractional coordinates                                    
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,162)ALPHA(LCT),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)      
      D=SQRT(H1P3**2+H3P3**2)                                            
      WRITE(NQ,362)H1P3,0.,H3P3,(A(I),I=4,6),D                           
  362 FORMAT(1X,'P3 AT',3F8.3,';  ',3F8.4,'  D =',F8.2)                  
  370 IF(LCT.GE.3)GO TO 1000                                             
!-----Code SB completed, go to next designated space group              
      IF(ITR(60).GT.1)GO TO 1000                                         
!-----Space group Pc, P21 - Code SC - Calculated as Pa                  
  400 LCT=3                                                              
      H1P3=H1IH                                                          
      GO TO 301                                                          
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT               
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTA,5X,3(I9,I3),I9)                  
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
      END SUBROUTINE ZSRTA
