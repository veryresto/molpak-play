      SUBROUTINE ZTRTA                        
                                                     
!-----This subroutine finds structures in the space groups P1, P1bar,    
!     Pm, Pc, and P21 in which the coordination sphere contains 6 I      
!     molecules in plane-1,2. For AA, AB, AC, AD, AH. 
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
!-----Linear storage array and markers                                   
!     COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000),
!     &       CODE(5000)
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
!      DIMENSION A(6),C(6)                                                
      CHARACTER (2) :: ALPHA(5) = (/'AA','AB','AC','AD','AH'/)     
!
      INTEGER :: I, J, JCT, JE, K, KE, L, LCT, LE, MCT, N, NCT
!
      REAL    :: A(6), C(6)
      REAL    :: D, D1, D2, DX1A3, DX2A3
      REAL    :: DX21I, DX31I, DX32I
      REAL    :: DY1, DY1C, DY2, DY2C
      REAL    :: E, F, FXX, FXY
      REAL    :: H1A3, H1C, H1I, H1IH, H1IT, H1P3
      REAL    :: H2A3, H2C, H2I, H2IH, H21I, H2P3
      REAL    :: H3A3, H3C, H3I, H31I, H32I, H3P3
      REAL    :: HH1A3, HH2A3, HH31I, HH32I
      REAL    :: HL   
      REAL    :: SINE2, SINE3
      REAL    :: X1A3, X2A3, X3A3, X3P3, XX1, XX2
      REAL    :: X2I, X21I, X21IT, X3I, X31I, X32I, X31IT, X32IT
      REAL    :: X1C, X2C, X3C, XX1C, XX2C
      REAL    :: VX, Y1C, Y2C, Y3C                  

      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=1,5)                                          
    5 FORMAT(1X,'ZTRTA called',5I3)                                      
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
   25 N=MARK                                                                
!-----Collect contacts for double standoff                               
      X21IT=2.*X21I                                                      
      DO 29 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 29                                                 
      IF(IT(I+7).NE.0)GO TO 29                                           
      D1=TI(I+2)+X21IT                                                      
   26 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 26                                               
   27 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 27                                               
      IF(D1.LT.-CN(K))GO TO 29                                              
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 27                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 27                                                              
   29 CONTINUE                                                              
      NARK=N                                                             
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
      GO TO 25                                                           
   44 MCT=1                                                              
      JE=2.*JE                                                           
      DX21I=H1I/JE                                                       
      X21I=H21I+DX21I                                                    
      GO TO 25                                                           
   45 MCT=-1                                                             
!     X21I=MAX(H21I-DX21I,-H1IH)                                        
      X21I=H21I-DX21I                                                    
!     IF(X21I.LT.-H1IH)X21I=X21I+H1I                                    
      GO TO 25                                                           
   46 IF(JE.LT.KE)GO TO 44                                               
      IF(JE.GE.NNE)GO TO 48                                              
      KE=NNE                                                             
      GO TO 44                                                           
   48 CONTINUE                                                           
      H2IH=H2I/2.                                                       
      NARK=0                                                        
      IF(IPR.LT.3)GO TO 60                                                  
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,51)H2I,H1I,H21I                                              
   51 FORMAT(6H H2I =,F7.3,8H,  H1I =,F7.3,8H, H21I =,F7.3)                 
!-----Go to the first designated space group in ZTRTA                    
   60 GO TO(100,200,300,400,500,1000),ITR(1)                             
!-----Calculate axis-3 separation of layer of I molecules -              
!     Space group P1 - Code AA                                           
  100 ER=ERMF                                                               
      MARK=MARK2                                                            
      H3I=100.                                                           
      X3I=H3N                                                            
      JE=1                                                                  
      LE=1                                                                  
      KE=NE                                                                 
      DX31I=H1I/NV                                                          
      DX32I=H2I/NV                                                          
      X31I=-H1IH                                                         
      X32I=-H2IH                                                        
      NCT=1                                                                 
      MCT=1                                                                 
  101 N=MARK                                                                
!-----Collect second layer for double standoff.                          
      X31IT=2.*X31I                                                    
      X32IT=2.*X32I                                                    
      DO 119 I=1,NMOD                                                       
      DO 118 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=X31IT+W(1,J)-W(1,I)                                                
      D2=X32IT+W(2,J)-W(2,I)                                                
  112 D2=D2+H2I                                                             
      D1=D1+H21I                                                            
      IF(D2.LE.CN(K))GO TO 112                                              
  113 D2=D2-H2I                                                             
      D1=D1-H21I                                                            
      IF(D2.GT.CN(K))GO TO 113                                              
      IF(D2.LT.-CN(K))GO TO 118                                             
  114 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 114                                              
  115 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 115                                              
      IF(D1.LT.-CN(K))GO TO 113                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 115                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(3,J)-W(3,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 115                                                             
  118 CONTINUE                                                              
  119 CONTINUE                                                              
      NARK=N                                                           
!-----Collect first layer for single standoff                            
      DO 129 I=1,NMOD                                                       
      DO 128 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=X31I+W(1,J)-W(1,I)                                                 
      D2=X32I+W(2,J)-W(2,I)                                                 
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
      CALL MINHI(X3I,0.,H3I,H3M,KE)                                         
      IF(KCT)904,142,141                                              
  141 IF(MCT.LE.1)GO TO 144                                           
      J=MCT-NCT                                                       
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 157                              
      GO TO 143                                                       
  142 H31I=X31I                                                       
      H32I=X32I                                                       
      H3I=X3I                                                         
      NCT=MCT                                                         
  143 GO TO(144,149,150,151,152,153,154,155,148),MCT                  
  144 IF(LE.GE.NV)GO TO 145                                           
      LE=LE+1                                                         
      X32I=X32I+DX32I                                                 
      GO TO 101                                                       
  145 IF(JE.GE.NV)GO TO 146                                           
      JE=JE+1                                                         
      X31I=X31I+DX31I                                                 
      LE=1                                                            
      X32I=-H2IH                                                        
      GO TO 101                                                       
  146 NCT=6                                                           
  147 JE=2.*JE                                                        
      DX31I=H1I/JE                                                    
      DX32I=H2I/JE                                                    
      HH31I=H31I                                                      
      HH32I=H32I                                                      
      NCT=NCT+1                                                       
      IF(NCT.GT.9)NCT=NCT-8                                           
      GO TO(908,150,151,152,153,154,155,148,149),NCT                  
  148 MCT=2                                                           
      X31I=HH31I+DX31I                                                
      X32I=HH32I                                                      
      GO TO 101                                                       
  149 MCT=3                                                           
      X31I=HH31I+DX31I                                                
      X32I=HH32I+DX32I                                                
      GO TO 101                                                       
  150 MCT=4                                                           
      X31I=HH31I                                                      
      X32I=HH32I+DX32I                                                
      GO TO 101                                                       
  151 MCT=5                                                           
!     X31I=MAX(HH31I-DX31I,-H1IH)                                       
      X31I=HH31I-DX31I                                                
      IF(X31I.LT.-H1IH)X31I=X31I+H1I                                    
      X32I=HH32I+DX32I                                                
      GO TO 101                                                       
  152 MCT=6                                                           
!     X31I=MAX(HH31I-DX31I,-H1IH)                                       
      X31I=HH31I-DX31I                                                
      IF(X31I.LT.-H1IH)X31I=X31I+H1I                                    
      X32I=HH32I                                                      
      GO TO 101                                                       
  153 MCT=7                                                           
!     X31I=MAX(HH31I-DX31I,-H1IH)                                       
      X31I=HH31I-DX31I                                                
      IF(X31I.LT.-H1IH)X31I=X31I+H1I                                    
!     X32I=MAX(HH32I-DX32I,-H2IH)                                       
      X32I=HH32I-DX32I                                                
      IF(X32I.LT.-H2IH)X32I=X32I+H2I                                    
      GO TO 101                                                       
  154 MCT=8                                                           
      X31I=HH31I                                                      
!     X32I=MAX(HH32I-DX32I,-H2IH)                                       
      X32I=HH32I-DX32I                                                
      IF(X32I.LT.-H2IH)X32I=X32I+H2I                                    
      GO TO 101                                                       
  155 MCT=9                                                           
      X31I=HH31I+DX31I                                                
!     X32I=MAX(HH32I-DX32I,-H2IH)                                       
      X32I=HH32I-DX32I                                                
      IF(X32I.LT.-H2IH)X32I=X32I+H2I                                    
      GO TO 101                                                       
  157 IF(JE.LT.KE)GO TO 147                                           
      IF(JE.GE.NNE)GO TO 158                                          
      KE=NNE                                                          
      GO TO 147                                                       
  158 CONTINUE                                                        
      NARK=0                                                           
!-----Calculate the cell volume per molecule                             
      VX=H1I*H2I*H3I                                                        
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
      HT1(13)=H31I                                                          
      HT1(14)=H32I                                                          
      NLOW=1                                                                
      MLOW=1                                                             
  161 IF(IPR.LT.2)GO TO 170                                              
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=SQRT(H3I**2+H31I**2+H32I**2)                                  
!     Cosines of cell angles                                             
      C(1)=(H21I*H31I+H2I*H32I)/(A(2)*A(3))                              
      C(2)=H31I/A(3)                                                     
      C(3)=H21I/A(2)                                                     
      C(4)=57.296*ACOS(C(1))                                             
      C(5)=57.296*ACOS(C(2))                                             
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,162)ALPHA(1),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)        
  162 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3/1X,'  COSINES' &      
     &,3F7.4,'  ANGLES',3F7.2)                                           
!-----Go to the next designated space group in ZTRTA                     
  170 GO TO(200,300,400,500,1000),ITR(2)                                 
!-----Calculate the separation for a layer of C molecules in both the    
!     positive and negative axis-3 directions -                          
!     Space group P1bar - Code AB                                        
  200 MARK=MARK2                                                            
      ER=ERMF                                                               
!-----Set up for positive direction                                      
      JCT=0                                                              
      X3C=100.                                                           
      Y3C=H3N                                                               
      HL=H3M                                                                
      GO TO 213                                                          
!-----Set up for negative direction                                      
  205 JCT=1                                                              
      X3C=-100.                                                          
      Y3C=-H3N                                                              
      HL=-H3M                                                               
  213 JE=1                                                                  
      LE=1                                                                  
      KE=NE                                                                 
      DY1C=H1I/NV                                                           
      DY2C=H2I/NV                                                           
      Y1C=-H1IH                                                          
      Y2C=-H2IH                                                         
      NCT=1                                                                 
      MCT=1                                                                 
  216 N=MARK                                                                
      DO 229 I=1,NMOD                                                       
      DO 228 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=Y1C-W(1,J)-W(1,I)                                                  
      D2=Y2C-W(2,J)-W(2,I)                                                  
  222 D2=D2+H2I                                                             
      D1=D1+H21I                                                            
      IF(D2.LE.CN(K))GO TO 222                                              
  223 D2=D2-H2I                                                             
      D1=D1-H21I                                                            
      IF(D2.GT.CN(K))GO TO 223                                              
      IF(D2.LT.-CN(K))GO TO 228                                             
  224 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 224                                              
  225 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 225                                              
      IF(D1.LT.-CN(K))GO TO 223                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 225                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=-W(3,J)-W(3,I)                                               
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 225                                                             
  228 CONTINUE                                                              
  229 CONTINUE                                                              
      NEND=N-1                                                              
      CALL MINHI(Y3C,0.,X3C,HL,KE)                                          
      IF(KCT)904,242,241                                              
  241 IF(MCT.LE.1)GO TO 244                                           
      J=MCT-NCT                                                       
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 257                              
      GO TO 243                                                       
  242 X1C=Y1C                                                         
      X2C=Y2C                                                         
      X3C=Y3C                                                         
      NCT=MCT                                                         
  243 GO TO(244,249,250,251,252,253,254,255,248),MCT                  
  244 IF(LE.GE.NV)GO TO 245                                           
      LE=LE+1                                                         
      Y2C=Y2C+DY2C                                                    
      GO TO 216                                                       
  245 IF(JE.GE.NV)GO TO 246                                           
      JE=JE+1                                                         
      Y1C=Y1C+DY1C                                                    
      LE=1                                                            
      Y2C=-H2IH                                                         
      GO TO 216                                                       
  246 NCT=6                                                           
  247 JE=2.*JE                                                        
      DY1C=H1I/JE                                                     
      DY2C=H2I/JE                                                     
      XX1C=X1C                                                        
      XX2C=X2C                                                        
      NCT=NCT+1                                                       
      IF(NCT.GT.9)NCT=NCT-8                                           
      GO TO(908,250,251,252,253,254,255,248,249),NCT                  
  248 MCT=2                                                           
      Y1C=XX1C+DY1C                                                   
      Y2C=XX2C                                                        
      GO TO 216                                                       
  249 MCT=3                                                           
      Y1C=XX1C+DY1C                                                   
      Y2C=XX2C+DY2C                                                   
      GO TO 216                                                       
  250 MCT=4                                                           
      Y1C=XX1C                                                        
      Y2C=XX2C+DY2C                                                   
      GO TO 216                                                       
  251 MCT=5                                                           
      Y1C=XX1C-DY1C                                                   
      IF(Y1C.LT.-H1IH)Y1C=Y1C+H1I                                       
      Y2C=XX2C+DY2C                                                   
      GO TO 216                                                       
  252 MCT=6                                                           
      Y1C=XX1C-DY1C                                                   
      IF(Y1C.LT.-H1IH)Y1C=Y1C+H1I                                       
      Y2C=XX2C                                                        
      GO TO 216                                                       
  253 MCT=7                                                           
      Y1C=XX1C-DY1C                                                   
      IF(Y1C.LT.-H1IH)Y1C=Y1C+H1I                                       
      Y2C=XX2C-DY2C                                                   
      IF(Y2C.LT.-H2IH)Y2C=Y2C+H2I                                       
      GO TO 216                                                       
  254 MCT=8                                                           
      Y1C=XX1C                                                        
      Y2C=XX2C-DY2C                                                   
      IF(Y2C.LT.-H2IH)Y2C=Y2C+H2I                                       
      GO TO 216                                                       
  255 MCT=9                                                           
      Y1C=XX1C+DY1C                                                   
      Y2C=XX2C-DY2C                                                   
      IF(Y2C.LT.-H2IH)Y2C=Y2C+H2I                                       
      GO TO 216                                                       
  257 IF(JE.LT.KE)GO TO 247                                           
      IF(JE.GE.NNE)GO TO 258                                          
      KE=NNE                                                          
      GO TO 247                                                       
  258 CONTINUE                                                        
!-----Has negative axis-3 direction been calculated                      
      IF(JCT.NE.0)GO TO 259                                              
      H1C=X1C                                                            
      H2C=X2C                                                            
      H3C=X3C                                                            
      GO TO 205                                                          
  259 H3I=H3C-X3C                                                        
      VX=H1IH*H2I*H3I                                                    
      IF(VX.GT.V3)GO TO 290                                              
!-----Check distances to derived positions of I molecules                
      X31I=H1C-X1C                                                       
      X32I=H2C-X2C                                                       
      N=MARK                                                             
      DO 269 I=1,NMOD                                                       
      DO 268 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=X31I+W(1,J)-W(1,I)                                                 
      D2=X32I+W(2,J)-W(2,I)                                                 
  262 D2=D2+H2I                                                             
      D1=D1+H21I                                                            
      IF(D2.LE.CN(K))GO TO 262                                              
  263 D2=D2-H2I                                                             
      D1=D1-H21I                                                            
      IF(D2.GT.CN(K))GO TO 263                                              
      IF(D2.LT.-CN(K))GO TO 268                                             
  264 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 264                                              
  265 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 265                                              
      IF(D1.LT.-CN(K))GO TO 263                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 265                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(3,J)-W(3,I)                                                
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 265                                                             
  268 CONTINUE                                                              
  269 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 280                                             
      NEND=N-1                                                              
      X3I=H3M                                                            
      CALL MINHI(X3I,H3I,100.,H3M,KE)                                       
      IF(KCT)280,270,906                                                 
!-----I molecules are too close, adjust H3I and H3C                      
  270 H3C=H3C+(X3I-H3I)/2.                                               
      H3I=X3I                                                            
      VX=H1IH*H2I*H3I                                                    
      IF(VX.GT.V3)GO TO 290                                              
  280 V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(2)                                                 
      IF(V3.GE.HT1(8))GO TO 290                                             
      H31I=H1C-X1C                                                       
      H32I=H2C-X2C                                                       
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=H1C                                                            
      HT1(6)=H2C                                                            
      HT1(7)=H3C                                                            
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(2)                                                      
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=H31I                                                          
      HT1(14)=H32I                                                          
      NLOW=2                                                                
      MLOW=2                                                             
  290 IF(IPR.LT.2)GO TO 295                                              
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=SQRT(H3I**2+H31I**2+H32I**2)                                  
!     Cosines of cell angles                                             
      C(1)=(H21I*H31I+H2I*H32I)/(A(2)*A(3))                              
      C(2)=H31I/A(3)                                                     
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates  
      SINE2=SQRT(1.-C(2)**2)                                             
      SINE3=SQRT(1.-C(3)**2)                                             
      E=(C(2)*C(3)-C(1))/(SINE2*SINE3)                                      
      F=SQRT(1.-E**2)                                                    
      FXX=SINE2*F                                                           
      FXY=SINE2*E/SINE3                                                     
!-----Transform first to Angstrom coordinates                            
      A(6)=H3C/FXX                                                       
      A(5)=H2C/SINE3+A(6)*FXY                                            
      A(4)=H1C-A(6)*C(2)-A(5)*C(3)                                       
!-----Then to fractional coordinates                                     
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      C(4)=57.296*ACOS(C(1))                                             
      C(5)=57.296*ACOS(C(2))                                             
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,162)ALPHA(2),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)        
      D=SQRT(H1C**2+H2C**2+H3C**2)                                       
      WRITE(NQ,293)H1C,H2C,H3C,(A(I),I=4,6),D                            
  293 FORMAT(1X,'C AT',3F8.3,';  ',3F8.4,'  D =',F8.2)                   
!-----Go to next designated space group in ZTRTA                         
  295 GO TO(300,400,500,1000),ITR(3)                                     
!-----Place layer of P molecules along both plus and minus axis-3        
!     Space groups Pm and/or Pc - Codes AC and AD                        
!-----Space group Pm - Code AC                                           
  300 LCT=3                                                                 
      H1P3=0.                                                            
      H2P3=0.                                                            
  301 N=MARK2                                                            
      ER=ERMF                                                            
      MARK=N                                                             
      DO 319 I=1,NMOD                                                       
      DO 318 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=W(2,J)+H2P3-W(2,I)                                                 
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
      X3P3=-H3N                                                          
      CALL MINHI(X3P3,0.,-100.,-H3M,NNE)                                    
      H3I=H3P3-X3P3                                                      
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
      HT1(7)=H3P3                                                           
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=3                                                                
      MLOW=LCT                                                           
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
      A(5)=H2P3/SQRT(1.-C(3)**2)                                         
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
      D=SQRT(H1P3**2+H2P3**2+H3P3**2)                                    
      WRITE(NQ,362)H1P3,H2P3,H3P3,(A(I),I=4,6),D                         
  362 FORMAT(1X,'P3 AT',3F8.3,';  ',3F8.4,'  D =',F8.2)                  
  370 IF(LCT.GE.4)GO TO 470                                              
!-----Code AC completed, go to next designated space group               
      IF(ITR(4)-2)400,500,1000                                           
!-----Space group Pc - Code AD                                           
  400 LCT=4                                                              
      H1P3=H1IH                                                          
      H2P3=0.                                                            
      GO TO 301                                                          
!-----Code AD completed, go to next designated space group               
  470 IF(ITR(5).GT.1)GO TO 1000                                          
!-----Calculate axis-3 height of layer of A molecules -                  
!     Space group P21 - Code AH                                          
  500 MARK=MARK2                                                         
      ER=ERMF                                                               
  513 JE=1                                                                  
      LE=1                                                                  
      KE=NE                                                                 
      DX1A3=H1I/NV                                                          
      DX2A3=H2I/NV                                                          
      X1A3=-H1IH                                                         
      X2A3=-H2IH                                                        
      X3A3=H3N                                                           
      H1A3=0.                                                               
      H2A3=0.                                                               
      H3A3=100.                                                             
      NCT=1                                                                 
      MCT=1                                                                 
  516 N=MARK                                                                
      DO 529 I=1,NMOD                                                       
      DO 528 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=X1A3-W(1,J)-W(1,I)                                                 
      D2=X2A3-W(2,J)-W(2,I)                                                 
  522 D2=D2+H2I                                                             
      D1=D1+H21I                                                            
      IF(D2.LE.CN(K))GO TO 522                                              
  523 D2=D2-H2I                                                             
      D1=D1-H21I                                                            
      IF(D2.GT.CN(K))GO TO 523                                              
      IF(D2.LT.-CN(K))GO TO 528                                             
  524 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 524                                              
  525 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 525                                              
      IF(D1.LT.-CN(K))GO TO 523                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 525                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(3,J)-W(3,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 525                                                             
  528 CONTINUE                                                              
  529 CONTINUE                                                              
      NEND=N-1                                                              
      CALL MINHI(X3A3,0.,H3A3,H3M,KE)                                      
      IF(KCT)904,542,541                                              
  541 IF(MCT.LE.1)GO TO 544                                           
      J=MCT-NCT                                                       
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 557                              
      GO TO 543                                                       
  542 H1A3=X1A3                                                       
      H2A3=X2A3                                                       
      H3A3=X3A3                                                       
      NCT=MCT                                                         
  543 GO TO(544,549,550,551,552,553,554,555,548),MCT                  
  544 IF(LE.GE.NV)GO TO 545                                           
      LE=LE+1                                                         
      X2A3=X2A3+DX2A3                                                 
      GO TO 516                                                       
  545 IF(JE.GE.NV)GO TO 546                                           
      JE=JE+1                                                         
      X1A3=X1A3+DX1A3                                                 
      LE=1                                                            
      X2A3=-H2IH                                                        
      GO TO 516                                                       
  546 NCT=6                                                           
  547 JE=2.*JE                                                        
      DX1A3=H1I/JE                                                    
      DX2A3=H2I/JE                                                    
      HH1A3=H1A3                                                      
      HH2A3=H2A3                                                      
      NCT=NCT+1                                                       
      IF(NCT.GT.9)NCT=NCT-8                                           
      GO TO(908,550,551,552,553,554,555,548,549),NCT                  
  548 MCT=2                                                           
      X1A3=HH1A3+DX1A3                                                
      X2A3=HH2A3                                                      
      GO TO 516                                                       
  549 MCT=3                                                           
      X1A3=HH1A3+DX1A3                                                
      X2A3=HH2A3+DX2A3                                                
      GO TO 516                                                       
  550 MCT=4                                                           
      X1A3=HH1A3                                                      
      X2A3=HH2A3+DX2A3                                                
      GO TO 516                                                       
  551 MCT=5                                                           
!     X1A3=MAX(HH1A3-DX1A3,-H1IH)                                       
      X1A3=HH1A3-DX1A3                                                
      IF(X1A3.LT.-H1IH)X1A3=X1A3+H1I                                    
      X2A3=HH2A3+DX2A3                                                
      GO TO 516                                                       
  552 MCT=6                                                           
!     X1A3=MAX(HH1A3-DX1A3,-H1IH)                                       
      X1A3=HH1A3-DX1A3                                                
      IF(X1A3.LT.-H1IH)X1A3=X1A3+H1I                                    
      X2A3=HH2A3                                                      
      GO TO 516                                                       
  553 MCT=7                                                           
!     X1A3=MAX(HH1A3-DX1A3,-H1IH)                                       
      X1A3=HH1A3-DX1A3                                                
      IF(X1A3.LT.-H1IH)X1A3=X1A3+H1I                                    
!     X2A3=MAX(HH2A3-DX2A3,-H2IH)                                       
      X2A3=HH2A3-DX2A3                                                
      IF(X2A3.LT.-H2IH)X2A3=X2A3+H2I                                    
      GO TO 516                                                       
  554 MCT=8                                                           
      X1A3=HH1A3                                                      
!     X2A3=MAX(HH2A3-DX2A3,-H2IH)                                       
      X2A3=HH2A3-DX2A3                                                
      IF(X2A3.LT.-H2IH)X2A3=X2A3+H2I                                    
      GO TO 516                                                       
  555 MCT=9                                                           
      X1A3=HH1A3+DX1A3                                                
!     X2A3=MAX(HH2A3-DX2A3,-H2IH)                                       
      X2A3=HH2A3-DX2A3                                                
      IF(X2A3.LT.-H2IH)X2A3=X2A3+H2I                                    
      GO TO 516                                                       
  557 IF(JE.LT.KE)GO TO 547                                           
      IF(JE.GE.NNE)GO TO 558                                          
      KE=NNE                                                          
      GO TO 547                                                       
  558 CONTINUE                                                        
!-----Calculate the cell volume per molecule                             
      VX=H1I*H2I*H3A3                                                       
      IF(VX.GE.V3)GO TO 1000                                                
!-----Check derived I molecule distances                                 
      N=MARK                                                             
      DO 569 I=1,NMOD                                                       
      DO 568 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=W(2,J)-W(2,I)                                                      
      D1=W(1,J)-W(1,I)                                                      
  561 D2=D2+H2I                                                          
      D1=D1+H21I                                                         
      IF(D2.LT.CN(K))GO TO 561                                           
  562 D2=D2-H2I                                                          
      D1=D1-H21I                                                         
      IF(D2.GT.CN(K))GO TO 562                                           
      IF(D2.LT.-CN(K))GO TO 568                                          
  563 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 563                                           
  564 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 564                                           
      IF(D1.LT.-CN(K))GO TO 562                                          
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 564                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=-W(3,J)-W(3,I)                                                
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 564                                                             
  568 CONTINUE                                                              
  569 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      H3I=2.*H3A3                                                        
      X3I=H3I                                                            
      CALL MINHI(X3I,H3I,100.,H3I,NNE)                                      
      IF(KCT)580,570,906                                                 
!-----I molecules too close, adjust H3I and H3A3                         
  570 H3A3=H3A3+(X3I-H3I)/2.                                             
      H3I=X3I                                                            
      VX=H1IH*H2I*H3I                                                    
      IF(VX.GT.V3)GO TO 600                                              
  580 V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(5)                                                 
      IF(V3.GE.HT1(8))GO TO 600                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=H1A3                                                           
      HT1(6)=H2A3                                                           
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(5)                                                      
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=4                                                                
      MLOW=5                                                             
  600 IF(IPR.LT.2)GO TO 1000                                             
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
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
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,162)ALPHA(5),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)        
      D=SQRT(H1A3**2+H2A3**2+H3A3**2)                                    
      WRITE(NQ,606)H1A3,H2A3,H3A3,(A(I),I=4,6),D                         
  606 FORMAT(1X,'A3 AT',3F8.3,';  ',3F8.4,'  D =',F8.2)                  
      GO TO 1000                                                            
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
  908 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,909)                                                      
  909 FORMAT('0DTWO failure')                                            
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,921)                                                      
  921 FORMAT(21H0NO ATOM INTERACTIONS)                                   
  999 KILL=1                                                                
!1000 WRITE(NQ,'(8HHT1 8-12,4F8.3)') HT1(8),(HT1(I),I=10,12) ! TEST-03 # 1
!     WRITE(NQ,'(10HExit ZTRTA,3F8.3)') A1,A2,A3             ! TEST-03 # 1
1000  RETURN                                                                
      END SUBROUTINE ZTRTA 

