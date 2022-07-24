      SUBROUTINE ZSRTH                                          
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
!      DIMENSION SM(3),A(12)                                              
!      CHARACTER*2 ALPHA(2),SYM(3)                                        
!      DATA ALPHA/'SO','SP'/                                              
!      DATA SYM/'P2','A2','C '/                                           
!
      CHARACTER(2) :: ALPHA(2) = (/'SO','SP'/)
      CHARACTER(2) :: SYM(3) = (/'P2','A2','C '/)     
!
      INTEGER :: I, J, JE, K, KE, KTR, L, LCT, LCTM, MCT, N
!
      REAL :: D, D1, D2, D3, DX1A2, DX2P2
      REAL :: H1A2, H1C, H1I, H1P2      
      REAL :: H2A2, H2C, H2I, H2P2 
      real :: H3C, H3I, H3A2, H3P2
      REAL :: A(12), SM(3), VX, VH
      REAL :: X1A2, X2A2, X2P2, X3P2, XCT
      REAL :: Y3P2
!
      VH=1000000.                                                        
      KTR=ITR(76)                                                        
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=76,77)                                        
    5 FORMAT(1X,'ZSRTH called',2I3)                                      
   10 H3A2=0.                                                             
      IF(KTR-2)11,100,1000                                               
!-----Subgroup 1 - I molecules are in contact along axis-2.              
!     Calculate axis-2 (H2I) by I molecule separation.                   
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
!-----Determine H1A2 and H1I by placing a line of A2 molecules on either 
!     side of axis-1 - H3A2 is zero.                                     
      MARK=MARK2                                                         
      N=MARK                                                             
      DO 51 I=1,NMOD                                                        
      DO 50 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                             
      D3=-W(3,J)-W(3,I)                                                  
      IF(ABS(D3).GT.CN(K))GO TO 50                                       
      D2=H2A2+W(2,J)-W(2,I)                                              
   47 D2=D2+H2I                                                             
      IF(D2.LE.CN(K))GO TO 47                                               
   48 D2=D2-H2I                                                             
      IF(D2.LT.-CN(K))GO TO 50                                              
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 48                                             
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=-W(1,J)-W(1,I)                                                
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 48                                                              
   50 CONTINUE                                                              
   51 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      H1A2=H1N                                                           
      ER=ERMT                                                            
      CALL MINHI(H1A2,0.,100.,H1M,NNE)                                   
      IF(KCT)904,55,906                                                  
   55 X1A2=-H1N                                                          
      CALL MINHI(X1A2,0.,-100.,-H1M,NNE)                                 
      IF(KCT)904,56,906                                                  
   56 H1I=H1A2-X1A2                                                      
      H1P2=H1I/2.                                                        
      IF(IPR.LT.3)GO TO 200                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,59)H1A2,H1I                                               
   59 FORMAT(' H1A2 =',F7.3,',  H1I =',F7.3)                             
      GO TO 200                                                          
!-----Subgroup 2 - Identity molecules are in contact along axis-1        
  100 LCT=2                                                              
!     Calculate axis-1 (H1I) by I molecule separation.                   
      MARK=MARK2                                                         
      NARK=0                                                             
      N=MARK                                                             
      DO 102 I=MARK1,NEND1,NSTP1                                            
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 102                                                
      IF(IT(I+7).NE.0)GO TO 102                                          
      D2=TI(I+3)                                                            
      IF(ABS(D2).GT.CN(K))GO TO 102                                         
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 102                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
  102 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      H1I=H1M                                                               
      CALL MINHI(H1I,0.,100.,H1M,NNE)                                       
      IF(KCT)904,105,906                                                 
  105 H1P2=H1I/2.                                                        
      IF(IPR.LT.3)GO TO 110                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,106)H1I                                                   
  106 FORMAT(1X,'H1I =',F7.3)                                            
!-----Determine H2A2, hence H2I, by approach of line of A2 molecules     
!     parallel to axis-1.                                                
  110 ER=ERMT                                                            
      KE=NE                                                              
      MCT=0                                                              
      JE=1                                                               
      X1A2=-H1I/2.                                                       
      DX1A2=H1I/NV                                                       
      H2A2=100.                                                          
      X2A2=H2N                                                           
      MARK=MARK2                                                         
      NARK=0                                                             
  116 N=MARK                                                             
      DO 130 I=1,NMOD                                                    
      DO 129 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D3=-W(3,J)-W(3,I)                                                  
      IF(ABS(D3).GT.CN(K))GO TO 129                                      
      D1=X1A2-W(1,J)-W(1,I)                                              
  120 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 120                                           
  121 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 129                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 121                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 121                                                          
  129 CONTINUE                                                           
  130 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      CALL MINHI(X2A2,0.,H2A2,H2M,KE)                                    
      IF(KCT)904,151,152                                                 
  151 H2A2=X2A2                                                          
      H1A2=X1A2                                                          
      IF(MCT)156,153,156                                                 
  152 IF(MCT)156,153,155                                                 
  153 IF(JE.GE.NV)GO TO 154                                              
      JE=JE+1                                                            
      X1A2=X1A2+DX1A2                                                    
      GO TO 116                                                          
  154 MCT=1                                                              
      JE=2*JE                                                            
      DX1A2=H1I/JE                                                       
      X1A2=H1A2+DX1A2                                                    
      GO TO 116                                                          
  155 MCT=-1                                                             
      X1A2=H1A2-DX1A2                                                    
      GO TO 116                                                          
  156 IF(JE.LT.KE)GO TO 154                                              
      IF(JE.GE.NNE)GO TO 158                                             
      KE=NNE                                                             
      GO TO 154                                                          
  158 CONTINUE                                                           
      H2I=2.*H2A2                                                        
      IF(IPR.LT.3)GO TO 200                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,160)H1A2,H2A2                                             
  160 FORMAT(' H1A2 =',F7.3,',  H2A2 =',F7.3)                            
!-----Collect I and A2 molecule positions for axis-3 double standoff    
  200 N=MARK                                                             
      DO 209 I=1,NMOD                                                    
      DO 208 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
!-----Set multiplier for I molecules                                     
      XCT=1.                                                             
      D1=W(1,J)-W(1,I)                                                   
      D2=W(2,J)-W(2,I)                                                   
  201 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 201                                           
  202 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 205                                          
  203 D2=D2+H2I                                                          
      IF(D2.LE.CN(K))GO TO 203                                           
  204 D2=D2-H2I                                                          
      IF(D2.LT.-CN(K))GO TO 202                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L))GO TO 204                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=XCT*W(3,J)-W(3,I)                                          
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 204                                                          
!-----Have A2 molecules been added?                                      
  205 IF(XCT.LE.0.)GO TO 208                                             
      XCT=-1.                                                            
      D1=H1A2-W(1,J)-W(1,I)                                              
      D2=D2+H2A2                                                         
      GO TO 201                                                          
  208 CONTINUE                                                           
  209 CONTINUE                                                           
      NARK=N                                                             
!-----Determine H3P2 and H3I by placing a plane of P2 and C molecules    
      ER=ERMF                                                            
      MCT=0                                                              
      H3P2=100.                                                          
      KE=NE                                                              
      X2P2=-H2A2                                                         
      DX2P2=H2I/NV                                                       
      JE=1                                                               
  216 N=NARK                                                             
      DO 249 I=1,NMOD                                                    
      DO 248 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
!-----Set multiplier for P2 molecules                                    
      XCT=1.                                                             
      D2=X2P2-W(2,J)-W(2,I)                                              
      D1=H1P2+W(1,J)-W(1,I)                                              
  221 D2=D2+H2I                                                          
      IF(D2.LE.CN(K))GO TO 221                                           
  222 D2=D2-H2I                                                          
      IF(D2.LT.-CN(K))GO TO 230                                          
  223 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 223                                           
  224 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 222                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L))GO TO 224                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=XCT*W(3,J)-W(3,I)                                          
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 224                                                          
!-----Have C molecule positions been added?                              
  230 IF(XCT.LT.0.)GO TO 248                                             
      D2=D2+H2A2                                                         
      D1=H1P2+H1A2-W(1,J)-W(1,I)                                         
      XCT=-1.                                                            
      GO TO 221                                                          
  248 CONTINUE                                                           
  249 CONTINUE                                                           
      IF(N.LE.NARK)GO TO 920                                             
      NEND=N-1                                                           
      X3P2=H3N                                                           
      CALL MINHI(X3P2,0.,H3P2,H3M,KE)                                    
      IF(KCT)904,250,252                                                 
  250 Y3P2=-H3N                                                          
      CALL MINHI(Y3P2,0.,-H3P2,-H3M,KE)                                  
      IF(KCT)904,251,252                                                 
  251 H2P2=X2P2                                                          
      H3P2=X3P2                                                          
      IF(-Y3P2.GT.X3P2)H3P2=-Y3P2                                        
      IF(MCT)256,253,256                                                 
  252 IF(MCT)256,253,255                                                 
  253 IF(JE.GE.NV)GO TO 254                                              
      JE=JE+1                                                            
      X2P2=X2P2+DX2P2                                                    
      GO TO 216                                                          
  254 MCT=1                                                              
      JE=2*JE                                                            
      DX2P2=H2I/JE                                                       
      X2P2=H2P2+DX2P2                                                    
      GO TO 216                                                          
  255 MCT=-1                                                             
      X2P2=H2P2-DX2P2                                                    
      GO TO 216                                                          
  256 IF(JE.LT.KE)GO TO 254                                              
      IF(JE.GE.NNE)GO TO 258                                             
      KE=NNE                                                             
      GO TO 254                                                          
  258 CONTINUE                                                           
      IF(IPR.LT.3)GO TO 260                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,259)H2P2,H3P2                                             
  259 FORMAT('H2P2 =',F7.3,',  H3P2 =',F7.3)                             
  260 H3I=2.*H3P2                                                        
      VX=H1I*H2I*H3I/4.                                                     
      H1C=H1P2+H1A2                                                      
      H3C=H3P2                                                           
      H2C=H2P2+H2A2                                                      
      IF(H1C.GT.H1I/2.)H1C=H1C-H1I                                       
      IF(H1C.LT.-H1I/2.)H1C=H1C+H1I                                      
      IF(H2C.GT.H2I/2.)H2C=H2C-H2I                                       
      IF(H2C.LT.-H2I/2.)H2C=H2C+H2I                                      
      IF(VX.GE.V3)GO TO 265                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 265                                             
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
      NLOW=33                                                               
      MLOW=LCT+75                                                        
  265 IF(IPR.LT.2)GO TO 270                                              
      IF(VX.GT.VH)GO TO 270                                              
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
  270 IF(LCT.NE.1)GO TO 400                                                 
      IF(ITR(77).EQ.1)GO TO 100                                             
  400 IF(IPR.LT.2)GO TO 1000                                             
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,405)ALPHA(LCTM),A1,A2,A3,VH,(A(I),I=1,3)                  
  405 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3)                    
  406 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)                                 
      SM(1)=A(4)-A(10)/2.                                                
      SM(2)=A(5)-A(11)/2.                                                
      SM(3)=A(6)-A(12)/2.                                                
      WRITE(NQ,406)SYM(1),(A(I),I=4,6),SM                                
      SM(1)=A(7)-A(10)/2.                                                
      SM(2)=A(8)-A(11)/2.                                                
      SM(3)=A(9)-A(12)/2.                                                
      WRITE(NQ,406)SYM(2),(A(I),I=7,9),SM                                
      SM(1)=A(10)/2.                                                     
      SM(2)=A(11)/2.                                                     
      SM(3)=A(12)/2.                                                     
      WRITE(NQ,406)SYM(3),(A(I),I=10,12),SM                              
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTH,5X,3(I9,I3),I9)                  
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
      WRITE(NQ,921)KCT,MARK1,MARK2,NEND2,MARK                               
  921 FORMAT(16H0NO INTERACTIONS,I3,4I9)                                    
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZSRTH 
