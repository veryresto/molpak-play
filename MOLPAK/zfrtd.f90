      SUBROUTINE ZFRTD                                                    
                                                                      
!-----This subroutine finds structures in the space group P21/c in which 
!     the coordination sphere contains 6 I molecules in plane-1,2, 4 C,  
!     and 4 P3 molecules.  For FA and FB.                                               
!
      USE molpakCommonMod
!
      IMPLICIT NONE
!
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
!      DATA ALPHA/'FA','FB'/                                                 
!      DATA SYM/'C ','A3','P3'/                                           
!
      CHARACTER(2) :: ALPHA(2) = (/'FA','FB'/)
      CHARACTER(2) :: SYM(3) = (/'C ','A3','P3'/)
!
      INTEGER :: I, J, JE, K, KE, L, LE, LCT, MCT, N, NCT
!
      REAL :: A(12), C(6)
      REAL :: D, D1, D2, DX21I, DY1C, DY2C
      REAL :: H1A3, H1I, H1C, H1IH, H1IT, H1P3, HL
      REAL :: H2A3, H2C, H21I, H2I, H2IH, H2P3
      REAL :: H3A3, H3C, H3I,  H3P3
      REAL :: SM(3), SINE3, VX
      REAL :: X1C 
      REAL :: X2C,  X2I, X21I, X21IT
      REAL :: X3C, X3A3, XX1C, XX2C
      REAL :: Y1C, Y2C, Y3C
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                             
      WRITE(NQ,5)                                                      
    5 FORMAT(1X,'ZFRD called')                                         
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
  121 N=MARK                                                           
!-----Collect distances for second layer                  
      X21IT=2.*X21I                                      
      DO 128 I=MARK1,NEND1,NSTP1                                       
      K=IT(I)                                                          
      IF(K.LE.0)GO TO 128                                              
      IF(IT(I+7).NE.0)GO TO 128                                        
      D1=TI(I+2)+X21IT                                                 
  122 D1=D1+H1I                                                        
      IF(D1.LE.CN(K))GO TO 122                                         
  123 D1=D1-H1I                                                        
      IF(D1.GT.CN(K))GO TO 123                                         
      IF(D1.LT.-CN(K))GO TO 128                                        
      D=TI(I+1)+D1**2                                                  
      L=IT(I+4)                                                        
      IF(D.GT.CN(L))GO TO 123                                          
      IT(N)=K                                                          
      TI(N+1)=D                                                        
      TI(N+2)=TI(I+3)                                                  
      N=N+NSTP                                                         
      IF(N.GT.LIMIT)GO TO 902                                          
      GO TO 123                                                        
  128 CONTINUE                                                         
      NARK=N                                 
!-----Collect distances for first layer         
      DO 138 I=MARK1,NEND1,NSTP1                                       
      K=IT(I)                                                          
      IF(K.LE.0)GO TO 138                                              
      IF(IT(I+7).NE.0)GO TO 138                                        
      D1=TI(I+2)+X21I                                                  
  132 D1=D1+H1I                                                        
      IF(D1.LE.CN(K))GO TO 132                                         
  133 D1=D1-H1I                                                        
      IF(D1.GT.CN(K))GO TO 133                                         
      IF(D1.LT.-CN(K))GO TO 138                                        
      D=TI(I+1)+D1**2                                                  
      L=IT(I+4)                                                        
      IF(D.GT.CN(L))GO TO 133                                          
      IT(N)=K                                                          
      TI(N+1)=D                                                        
      TI(N+2)=TI(I+3)                                                  
      N=N+NSTP                                                         
      IF(N.GT.LIMIT)GO TO 902                                          
      GO TO 133                                                        
  138 CONTINUE                                                         
      IF(N.LE.MARK)GO TO 920                                           
      NEND=N-1                                                         
      CALL MINHI(X2I,0.,H2I,H2M,KE)                                    
      IF(KCT)904,141,142                                                 
  141 H21I=X21I                                                          
      H2I=X2I                                                            
      IF(MCT)146,143,146                                                 
  142 IF(MCT)146,143,145                                                 
  143 IF(JE.GE.NV)GO TO 144                                              
      JE=JE+1                                                            
      X21I=X21I+DX21I                                                    
      GO TO 121                                                          
  144 MCT=1                                                              
      JE=2.*JE                                                           
      DX21I=H1I/JE                                                       
      X21I=H21I+DX21I                                                    
      GO TO 121                                                          
  145 MCT=-1                                                             
      X21I=H21I-DX21I                                                    
      IF(X21I.LT.-H1IH)X21I=X21I+H1I      
      GO TO 121                                                          
  146 IF(JE.LT.KE)GO TO 144                                              
      IF(JE.GE.NNE)GO TO 148                                             
      KE=NNE                                                             
      GO TO 144                                                          
  148 H2IH=H2I/2.                                                        
      NARK=0                                                         
      IF(IPR.LT.3)GO TO 200                                            
      CALL PAGE(1,1,0)                                                 
      WRITE(NQ,151)H2I,H1I,H21I                                        
  151 FORMAT(6H H2I =,F7.3,8H,  H1I =,F7.3,8H, H21I =,F7.3)            
!-----Calculate the separation for a layer of C molecules in the         
!     positive axis-3 direction.                                         
!     Space group P21/c - Code CD                                        
  200 MARK=MARK2                                                       
      ER=ERMF                                                          
!-----Set up for positive direction                                      
      X3C=100.                                                         
      Y3C=H3N                                                          
      HL=H3M                                                           
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
      H1C=X1C                                                          
      H2C=X2C                                                          
      H3C=X3C                                                          
      H2P3=0.                                                          
!-----Space group P21/c or P21/m                                         
      IF(ITR(72)-2)270,275,1000                                        
  270 LCT=1                                                            
      H1P3=H1IH                                                        
      GO TO 301                                                        
  275 LCT=2                                                            
      H1P3=0.                                                          
!-----Place layer of P molecules along minus axis-3 direction.           
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
      H3P3=-H3N                                                        
      CALL MINHI(H3P3,0.,-100.,-H3M,NNE)                               
      H1A3=H1C-H1P3                                                    
      H2A3=H2C-H2P3                                                    
      H3A3=H3C-H3P3                                                    
      H3I=2.*H3A3                                                      
      VX=H1IH*H2IH*H3I                                                 
      IF(VX.GT.V3)GO TO 390                                            
!-----Check derived I and A molecule distances                           
      N=MARK                                                           
      DO 349 I=1,NMOD                                                  
      DO 348 J=1,NMOD                                                  
      K=IA(I)+IAA(J)                                                   
      L=K+10                                                           
      D2=W(2,J)-W(2,I)                                                   
      D1=W(1,J)-W(1,I)                                                 
  341 D2=D2+H2I                                                        
      D1=D1+H21I                                                       
      IF(D2.LT.CN(K))GO TO 341                                         
  342 D2=D2-H2I                                                        
      D1=D1-H21I                                                       
      IF(D2.GT.CN(K))GO TO 342                                         
      IF(D2.LT.-CN(K))GO TO 348                                        
  343 D1=D1+H1I                                                        
      IF(D1.LT.CN(K))GO TO 343                                         
  344 D1=D1-H1I                                                        
      IF(D1.GT.CN(K))GO TO 344                                         
      IF(D1.LT.-CN(K))GO TO 342                                        
      D=D1**2+D2**2                                                    
      IF(D.GT.CN(L))GO TO 344                                          
      IT(N)=K                                                          
      TI(N+1)=D                                                        
      TI(N+2)=W(3,J)-W(3,I)                                            
      N=N+NSTP                                                         
      IF(N.GT.LIMIT)GO TO 902                                          
      GO TO 344                                                        
  348 CONTINUE                                                         
  349 CONTINUE                                                         
      NARK=N                                               
      DO 359 I=1,NMOD                                                  
      DO 358 J=1,NMOD                                                  
      K=IA(I)+IAA(J)                                                   
      L=K+10                                                           
      D2=H2A3-W(2,J)-W(2,I)                                              
      D1=H1A3-W(1,J)-W(1,I)                                            
  351 D2=D2+H2I                                                        
      D1=D1+H21I                                                       
      IF(D2.LT.CN(K))GO TO 351                                         
  352 D2=D2-H2I                                                        
      D1=D1-H21I                                                       
      IF(D2.GT.CN(K))GO TO 352                                         
      IF(D2.LT.-CN(K))GO TO 358                                        
  353 D1=D1+H1I                                                        
      IF(D1.LT.CN(K))GO TO 353                                         
  354 D1=D1-H1I                                                        
      IF(D1.GT.CN(K))GO TO 354                                         
      IF(D1.LT.-CN(K))GO TO 352                                        
      D=D1**2+D2**2                                                    
      IF(D.GT.CN(L))GO TO 354                                          
      IT(N)=K                                                          
      TI(N+1)=D                                                        
      TI(N+2)=-W(3,J)-W(3,I)                                           
      N=N+NSTP                                                         
      IF(N.GT.LIMIT)GO TO 902                                          
      GO TO 354                                                        
  358 CONTINUE                                                         
  359 CONTINUE                                                         
      IF(N.LE.MARK)GO TO 920                                           
      NEND=N-1                                                         
      X3A3=H3A3                                                        
      CALL MINHI(X3A3,H3A3,100.,H3I,NNE)                               
      IF(KCT)375,370,906                                               
!-----A molecules too close in +axis-3, adjust H3C and H3P3              
  370 H3P3=H3P3-(X3A3-H3A3)/2.                                         
      H3C=X3C+(X3A3-H3A3)/2.                                           
      H3A3=X3A3                                                        
      H3I=2.*H3A3                                                      
      VX=H1IH*H2IH*H3I                                                 
      IF(VX.GT.V3)GO TO 390                                            
  375 X3A3=-H3A3                                                       
      CALL MINHI(X3A3,-H3A3,-100.,-H3I,NNE)                            
      IF(KCT)380,376,906                                               
!-----A molecules too close in -axis-3, adjust H3C and H3P3              
  376 H3P3=H3P3+(X3A3+H3A3)/2.                                         
      H3C=X3C-(X3A3+H3A3)/2.                                           
      H3A3=-X3A3                                                       
      H3I=2.*H3A3                                                      
      VX=H1IH*H2IH*H3I                                                 
      IF(VX.GT.V3)GO TO 390                                            
  380 V3=VX                                                            
      DIJ(NP2)=V3                                                      
      CDIJ(NP2+20)=ALPHA(LCT)                                          
      IF(V3.GE.HT1(8))GO TO 390                                        
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
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
      NLOW=30                                                               
      MLOW=71+LCT                                                        
  390 IF(IPR.LT.2)GO TO 400                                            
      A(1)=H1I                                                         
      A(2)=SQRT(H21I**2+H2I**2)                                        
      A(3)=H3I                                                         
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
      WRITE(NQ,395)ALPHA(LCT),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)      
  395 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3  &                   
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
  396 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)                               
      SM(1)=A(4)-A(4)                                                  
      SM(2)=A(5)-A(5)                                                  
      SM(3)=A(6)-A(6)                                                  
      WRITE(NQ,396)SYM(1),(A(I),I=4,6),SM                              
      SM(1)=A(10)-A(4)                                                 
      SM(2)=A(11)-A(5)                                                 
      SM(3)=A(12)                                                      
      WRITE(NQ,396)SYM(2),(A(I),I=10,12),SM                            
      SM(1)=A(7)                                                       
      SM(2)=A(8)                                                       
      SM(3)=A(9)-A(6)                                                  
      WRITE(NQ,396)SYM(3),(A(I),I=7,9),SM                              
!-----Has P21/m been calculated?              
  400 IF(LCT.GE.2)GO TO 1000                                           
      IF(ITR(73).GT.1)GO TO 1000                                       
      H3C=X3C                                                          
      GO TO 275                                                        
  902 CALL PAGE(2,2,0)                                                 
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT             
  903 FORMAT(25H0STORAGE EXCEEDED BY ZFRD,5X,3(I9,I3),I9)              
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
 1000 RETURN                                                           
      END SUBROUTINE ZFRTD           
