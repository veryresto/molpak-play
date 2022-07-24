      SUBROUTINE ZSRTI                                                      
!
!     For SQ, SR and SS.
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
!      DIMENSION SM(3),A(12)                                              
!      CHARACTER*2 ALPHA(3),SYM(3)                                        
!      DATA ALPHA/'SQ','SR','SS'/                                         
!      DATA SYM/'P2','A2','C '/ 
!                                          
      CHARACTER(2) :: ALPHA(3) = (/'SQ','SR','SS'/)
      CHARACTER(2) :: SYM(3) = (/'P2','A2','C '/)
!
      INTEGER :: I, ICT, J, JE, K, KE, KTR, L, LCT, LCTM
      INTEGER :: MCT, N, NCT 
!
      REAL :: A(12)
      REAL :: D, D1, D2, D3, DX3A2
      REAL :: H1A2, H1C, H1I, H1P2, HL
      REAL :: H2A2, H2C, H2I, H2P2
      REAL :: H3A2, H3C, H3I, H3P2
      REAL :: SM(3), VH, VX
      REAL :: X1I, X1A2, X3A2, X3C, XCT
      REAL :: Y1A2
!
      NARK=0                                                             
      VH=1000000.                                                        
      KTR=ITR(78)                                                        
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=78,80)                                        
    5 FORMAT(1X,'ZSRTI called',3I3)                                      
   10 H1P2=0.                                                            
      H2C=0.                                                             
      IF(KTR-2)11,100,400                                                
!-----Subgroup 1 - I molecules are in contact along axis-2.              
!     Calculate axis-2 by I molecule separation.                         
   11 LCT=1                                                              
      MARK=MARK2                                                         
      NARK=0                                                             
      N=MARK                                                             
      DO 22 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 22                                                 
      IF(IT(I+7).NE.0)GO TO 22                                           
      D1=TI(I+2)                                                            
      IF(ABS(D1).GT.CN(K))GO TO 22                                          
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 22                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   22 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      H2I=H2M                                                               
      CALL MINHI(H2I,0.,100.,H2M,NNE)                                       
      IF(KCT)904,25,906                                                  
   25 IF(IPR.LT.3)GO TO 30                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,26)H2I                                                    
   26 FORMAT(1X,'H2I =',F7.3)                                            
   30 H2A2=H2I/2.                                                        
      H2P2=H2A2                                                          
!-----Determine H3P2 and H3I by placing a line of P2 molecules parallel  
!     to axis-2 - H1P2=0 - H3I=2*H3P2                                    
      MARK=MARK2                                                         
      N=MARK                                                             
      DO 51 I=1,NMOD                                                        
      DO 50 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                             
      D1=W(1,J)-W(1,I)                                                   
      IF(ABS(D1).GT.CN(K))GO TO 50                                       
      D2=H2P2-W(2,J)-W(2,I)                                              
   47 D2=D2+H2I                                                             
      IF(D2.LE.CN(K))GO TO 47                                               
   48 D2=D2-H2I                                                             
      IF(D2.GT.CN(K))GO TO 48                                               
      IF(D2.LT.-CN(K))GO TO 50                                              
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L))GO TO 48                                             
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(3,J)-W(3,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 48                                                              
   50 CONTINUE                                                              
   51 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      H3P2=H3N                                                           
      ER=ERMT                                                            
      CALL MINHI(H3P2,0.,100.,H3M,NNE)                                   
      IF(KCT)904,55,906                                                  
   55 H3I=2.*H3P2                                                        
      GO TO 300                                                          
!-----Subgroup 2 - Identity molecules are in contact along axis-3        
  100 LCT=2                                                              
      H3I=H3M                                                            
      H3P2=H3N                                                           
!-----Determine H2P2 and H2I by placing a row of P2 molecules parallel   
!     to axis-3 - H1P2=0 - H2I=2*H2P2                                    
      MARK=MARK2                                                         
      NARK=0                                                             
      N=MARK                                                             
      DO 109 I=MARK1,NEND1,NSTP1                                         
      K=-IT(I)                                                           
      IF(K.LE.0)GO TO 109                                                
      D1=TI(I+2)                                                         
      IF(ABS(D1).GT.CN(K))GO TO 109                                      
      D=TI(I+1)+D1**2                                                    
      L=IT(I+4)                                                          
      IF(D.GT.CN(L))GO TO 109                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+6)                                                    
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
  109 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      H2P2=H2N                                                              
      CALL MINHI(H2P2,0.,100.,H2M,NNE)                                      
      IF(KCT)904,110,906                                                 
  110 H2I=2.*H2P2                                                        
      IF(IPR.LT.3)GO TO 300                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,111)H2P2,H2I                                              
  111 FORMAT(1X,'H2P2 =',F7.3,',  H2I =',F7.3)                           
!-----Determine H1A2 and H1I by placing a plane of A2 and C molecules    
!-----Collect I and P2 molecule positions for axis-1 double standoff     
  300 N=MARK                                                             
      DO 309 I=1,NMOD                                                    
      DO 308 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
!-----Set flag for I molecules                                           
      ICT=1                                                              
      D3=W(3,J)-W(3,I)                                                   
      D2=W(2,J)-W(2,I)                                                   
  301 D3=D3+H3I                                                          
      IF(D3.LT.CN(K))GO TO 301                                           
  302 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 305                                          
  303 D2=D2+H2I                                                          
      IF(D2.LT.CN(K))GO TO 303                                           
  304 D2=D2-H2I                                                          
      IF(D2.LT.-CN(K))GO TO 302                                          
      D=D3**2+D2**2                                                      
      IF(D.GT.CN(K+10))GO TO 304                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 304                                                          
!-----Have P2 molecules been added?                                      
  305 IF(ICT.LT.0)GO TO 308                                              
      ICT=-1                                                             
      D2=H2P2-W(2,J)-W(2,I)                                              
      D3=D3+H3P2                                                         
      GO TO 301                                                          
  308 CONTINUE                                                           
  309 CONTINUE                                                           
      NARK=N                                                             
      ER=ERMF                                                            
      MCT=0                                                              
      H1I=100.                                                           
      H2A2=H2P2                                                          
      HL=H1M                                                             
      KE=NNE                                                             
      X3A2=-H3P2                                                         
      DX3A2=H3I/NV                                                       
      JE=1                                                               
  316 N=NARK                                                             
      DO 339 I=1,NMOD                                                    
      DO 338 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
!-----Set flag for A2 molecules                                          
      NCT=1                                                              
      D3=X3A2-W(3,J)-W(3,I)                                              
      D2=H2A2+W(2,J)-W(2,I)                                              
  321 D2=D2+H2I                                                          
      IF(D2.LT.CN(K))GO TO 321                                           
  322 D2=D2-H2I                                                          
      IF(D2.LT.-CN(K))GO TO 330                                          
  323 D3=D3+H3I                                                          
      IF(D3.LT.CN(K))GO TO 323                                           
  324 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 322                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 324                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=-W(1,J)-W(1,I)                                             
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 324                                                          
!-----Have C molecule positions been added?                              
  330 IF(NCT.LT.0)GO TO 338                                              
      D3=D3+H3P2                                                         
      D2=H2P2+H2A2-W(2,J)-W(2,I)                                         
      NCT=-1                                                             
      GO TO 321                                                          
  338 CONTINUE                                                           
  339 CONTINUE                                                           
      IF(N.LE.NARK)GO TO 920                                             
      NEND=N-1                                                           
      X1A2=H1N                                                           
      CALL MINHI(X1A2,0.,100.,HL,KE)                                     
      IF(KCT)904,340,906                                                 
  340 Y1A2=-H1N                                                          
      CALL MINHI(Y1A2,0.,-100.,-HL,KE)                                   
      IF(KCT)904,341,906                                                 
  341 X1I=X1A2-Y1A2                                                      
      IF(X1I.GE.H1I)GO TO 352                                            
      H1A2=X1A2                                                          
      H3A2=X3A2                                                          
      H1I=X1I                                                            
      IF(MCT)356,353,356                                                 
  352 IF(MCT)356,353,355                                                 
  353 IF(JE.GE.NV)GO TO 354                                              
      JE=JE+1                                                            
      X3A2=X3A2+DX3A2                                                    
      GO TO 316                                                          
  354 MCT=1                                                              
      JE=2*JE                                                            
      DX3A2=H1I/JE                                                       
      X3A2=H3A2+DX3A2                                                    
      GO TO 316                                                          
  355 MCT=-1                                                             
      X3A2=H3A2-DX3A2                                                    
      GO TO 316                                                          
  356 IF(JE.LT.KE)GO TO 354                                              
      IF(JE.GE.NNE)GO TO 358                                             
      KE=NNE                                                             
      GO TO 354                                                          
  358 H1C=H1A2                                                           
      H3C=H3P2+H3A2                                                      
      IF(IPR.LE.3)GO TO 360                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,359)H1A2,H3A2,H1I                                         
  359 FORMAT('H1A2 =',F7.3,',  H3A2 =',F7.3,',  H1I =',F7.3)             
  360 VX=H1I*H2I*H3I/4.                                                     
      IF(IPR.LT.2)GO TO 361                                              
      IF(VX.GT.VH)GO TO 361                                              
      VH=VX                                                              
      LCTM=LCT                                                           
      A(1)=H1I                                                           
      A(2)=H2I                                                           
      A(3)=H3I                                                           
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      A(6)=H3P2                                                          
      A(5)=H2P2                                                          
      A(4)=H1P2                                                          
      A(9)=H3A2                                                          
      A(8)=H2A2                                                          
      A(7)=H1A2                                                          
      A(12)=H3C                                                          
      A(11)=H2C                                                          
      A(10)=H1C                                                          
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
  361 IF(VX.GE.V3)GO TO 365                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 365                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
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
      NLOW=34                                                               
      MLOW=LCT+77                                                        
  365 IF(LCT.NE.1)GO TO 399                                                 
      IF(ITR(79)-2)100,400,510                                           
  399 IF(ITR(80).NE.1)GO TO 510                                          
!-----Subgroup 3 - I molecules are in contact along axis-1.              
!     Calculate axis-1 by I molecule separation.                         
  400 LCT=3                                                              
      MARK=MARK2                                                         
      NARK=0                                                             
      N=MARK                                                             
      DO 422 I=MARK1,NEND1,NSTP1                                            
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 422                                                
      IF(IT(I+7).NE.0)GO TO 422                                          
      D2=TI(I+3)                                                            
      IF(ABS(D2).GT.CN(K))GO TO 422                                         
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 422                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
  422 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      H1I=H1M                                                               
      CALL MINHI(H1I,0.,100.,H1M,NNE)                                       
      IF(KCT)904,425,906                                                 
  425 IF(IPR.LT.3)GO TO 430                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,426)H1I                                                   
  426 FORMAT(1X,'H1I =',F7.3)                                            
  430 H1A2=H1I/2.                                                        
      H1C=H1A2                                                           
!-----Calculate H3C, X3C, then H3I by placing rows of C molecules,       
!     parallel to axis-1, on either side of the origin along axis-3      
      MARK=MARK2                                                         
      N=MARK                                                             
      DO 439 I=1,NMOD                                                    
      DO 438 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      D1=H1C-W(1,J)-W(1,I)                                               
      D2=-W(2,J)-W(2,I)                                                  
  431 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 431                                           
  432 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 438                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(K+10))GO TO 432                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=-W(3,J)-W(3,I)                                             
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 432                                                          
  438 CONTINUE                                                           
  439 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      H3C=H3N                                                               
      CALL MINHI(H3C,0.,100.,H3M,NNE)                                       
      IF(KCT)904,440,906                                                 
  440 X3C=-H3N                                                           
      CALL MINHI(X3C,0.,-100.,-H3M,NNE)                                  
      IF(KCT)904,441,906                                                 
  441 H3I=H3C-X3C                                                        
      H3P2=H3I/2.                                                        
      H3A2=H3C-H3P2                                                      
      IF(IPR.LT.3)GO TO 450                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,445)H3C,H3I                                               
  445 FORMAT(1X,'H3C =',F7.3,',  H3I =',F7.3)                            
!-----Determine H2I by locating a plane of P2 and A2 molecules (at 1/2)  
!     while checking C and I molecules at double standoff.               
!     Generate C and I molecule positions.                               
  450 N=MARK2                                                            
      MARK=N                                                             
      ER=ERMF                                                            
      DO 459 I=1,NMOD                                                    
      DO 458 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
!-----Set multiplier for C molecules.                                    
      XCT=-1.                                                            
      D1=H1C-W(1,J)-W(1,I)                                               
      D3=H3C-W(3,J)-W(3,I)                                               
  451 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 451                                           
  452 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 455                                          
  453 D3=D3+H3I                                                          
      IF(D3.LE.CN(K))GO TO 453                                           
  454 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 452                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(K+10))GO TO 454                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=XCT*W(2,J)-W(2,I)                                          
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 454                                                          
!-----Have I molecules been placed?                                      
  455 IF(XCT.GT.0.)GO TO 458                                             
      XCT=1.                                                             
      D1=W(1,J)-W(1,I)                                                   
      D3=W(3,J)-W(3,I)                                                   
      GO TO 451                                                          
  458 CONTINUE                                                           
  459 CONTINUE                                                           
      NARK=N                                                             
!-----Generate A2 and P2 positions                                       
      DO 469 I=1,NMOD                                                    
      DO 468 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
!-----Set multiplier for P2 molecules                                    
      XCT=-1.                                                            
      D1=W(1,J)-W(1,I)                                                   
      D3=H3P2+W(3,J)-W(3,I)                                              
  461 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 461                                           
  462 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 465                                          
  463 D3=D3+H3I                                                          
      IF(D3.LE.CN(K))GO TO 463                                           
  464 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 462                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(K+10))GO TO 464                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=XCT*W(2,J)-W(2,I)                                          
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 464                                                          
!-----Have A2 molecules been placed?                                     
  465 IF(XCT.GT.0.)GO TO 468                                             
      XCT=1.                                                             
      D1=H1A2-W(1,J)-W(1,I)                                              
      D3=H3A2-W(3,J)-W(3,I)                                              
      GO TO 461                                                          
  468 CONTINUE                                                           
  469 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      H2P2=H2N                                                           
      CALL MINHI(H2P2,0.,100.,H2M,NNE)                                   
      IF(KCT)904,470,906                                                 
  470 H2A2=H2P2                                                          
      H2I=2.*H2P2                                                        
      H2A2=H2P2                                                          
      IF(IPR.LE.3)GO TO 480                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,475)H2P2,H2A2,H2I                                         
  475 FORMAT('H2P2 =',F7.3,',  H2A2 =',F7.3,',  H2I =',F7.3)             
  480 VX=H1I*H2I*H3I/4.                                                     
      IF(VX.GE.V3)GO TO 500                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 500                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
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
      NLOW=34                                                               
      MLOW=LCT+77                                                        
  500 IF(IPR.LT.2)GO TO 1000                                             
      IF(VX.GT.VH)GO TO 520                                              
      VH=VX                                                              
      LCTM=LCT                                                           
      A(1)=H1I                                                           
      A(2)=H2I                                                           
      A(3)=H3I                                                           
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      A(6)=H3P2                                                          
      A(5)=H2P2                                                          
      A(4)=H1P2                                                          
      A(9)=H3A2                                                          
      A(8)=H2A2                                                          
      A(7)=H1A2                                                          
      A(12)=H3C                                                          
      A(11)=H2C                                                          
      A(10)=H1C                                                          
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
      GO TO 520                                                          
  510 IF(IPR.LT.2)GO TO 1000                                             
  520 CALL PAGE(4,4,0)                                                   
      WRITE(NQ,525)ALPHA(LCTM),A1,A2,A3,VH,(A(I),I=1,3)                  
  525 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3)                    
  526 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)                                 
      SM(1)=A(4)-A(10)/2.                                                
      SM(2)=A(5)-A(11)/2.                                                
      SM(3)=A(6)-A(12)/2.                                                
      WRITE(NQ,526)SYM(1),(A(I),I=4,6),SM                                
      SM(1)=A(7)-A(10)/2.                                                
      SM(2)=A(8)-A(11)/2.                                                
      SM(3)=A(9)-A(12)/2.                                                
      WRITE(NQ,526)SYM(2),(A(I),I=7,9),SM                                
      SM(1)=A(10)/2.                                                     
      SM(2)=A(11)/2.                                                     
      SM(3)=A(12)/2.                                                     
      WRITE(NQ,526)SYM(3),(A(I),I=10,12),SM                              
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTI,5X,3(I9,I3),I9)                  
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
      WRITE(NQ,921)KCT,MARK1,MARK2,NEND2,MARK,NARK,N                        
  921 FORMAT(16H0NO INTERACTIONS,I3,6I7)                                    
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZSRTI 
