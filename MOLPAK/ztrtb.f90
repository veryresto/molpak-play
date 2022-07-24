       SUBROUTINE ZTRTB     
!      
!      For AE, AF, AG.
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
!     COMMON TI(100000)                                                
!      COMMON /NEW/IDCHG,IDATOM(200),G92CHARGE(200),ESPCHARGE(200)      
!      DIMENSION IT(100000)                                             
!      EQUIVALENCE (TI(1),IT(1))                                       
!
      CHARACTER(2) :: ALPHA(3) = (/'AE','AF','AG'/)
!
!      DIMENSION ALPHA(3) 
!      DIMENSION A(6),C(6)                                             
!      DATA ALPHA/'AE','AF','AG'/                                         
!
      INTEGER :: I, ICT, J, JCT, JE, K, KE, L, LCT, LE
      INTEGER :: MCT, MPL, N, NCT
!
      REAL  :: A(6), C(6)
      REAL  :: D, D1, D3, DX1P3, DX3P3, DY1A3
      REAL  :: H1A3, H1P3
      REAL  :: H2A3, H2I, H21I, H2P3
      REAL  :: H3A3, H3P3 
      REAL  :: HH1P3, HH3P3, HL, HN
      REAL  :: X1A3, X1I, X1P3, X2A3, X2I, X2P3, X3P3
      REAL  :: VX 
      REAl  :: Y1A3, Y2A3
!
      NARK=0                                                          
      IF(IPR.LT.2)GO TO 10                                            
      WRITE(NQ,5)(ITR(I),I=6,8)                                       
    5 FORMAT(1X,'ZTRTB called',3I3)                                   
!-----Go to first designated ZTRTB space group                      
   10 IF(ITR(6)-2)100,200,300                                         
!-----Space group P2 - Code AE - A molecules at zero on axis-3      
  100 LCT=1                                                           
!     Approach of grids of A molecules from both positive and negative
!     axis-2 directions                                             
!-----Set signal to pick up distances to atoms not offset by H3N; that
!     is, a two-fold axis, not a screw axis.                        
      ICT=1                                                           
      ER=ERMF                                                         
      H3A3=0.                                                         
!-----Set up for positive direction                                 
  105 JCT=0                                                           
      X2A3=100.                                                       
      HL=H2M                                                          
      HN=H2N                                                          
  110 KE=NE                                                           
      JE=1                                                            
      MCT=0                                                           
      MARK=MARK2                                                         
      DY1A3=H1M/NV                                                       
      Y1A3=-H1N                                                       
      MCT=0                                                              
  115 N=MARK                                                             
      Y2A3=HN                             
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
      IF(KCT)1903,151,152                                              
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
      GO TO 115                                                        
  156 IF(JE.LT.KE)GO TO 154                                            
      IF(JE.GE.NNE)GO TO 158                                           
      KE=NNE                                                           
      GO TO 154                                                        
  158 CONTINUE                                                          
!-----Has negative direction been calculated?                       
      IF(JCT.GT.0)GO TO 159                                           
      H1A3=X1A3                                                       
      H2A3=X2A3                                                       
      X2A3=-100.                                                      
      HN=-H2N                                                         
      HL=-H2M                                                         
      JCT=1                                                           
!-----Is the molecule centric?                                      
      IF(NCTR.NE.7)GO TO 110                                             
      H2I=2.*H2A3                                                        
      H21I=2.*H1A3                                                       
      GO TO 160                                                          
  159 H2I=H2A3-X2A3                                                   
      H21I=H1A3-X1A3                                                  
  160 VX=H1N*H2I*H3M                                                  
      IF(VX.GE.V3)GO TO 191                                              
!-----Check derived distances to I molecules                        
      MARK=MARK2                                                      
      N=MARK                                                             
      DO 169 I=1,NMOD                                                    
      DO 168 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D1=H21I+W(1,J)-W(1,I)                                              
      D3=W(3,J)-W(3,I)                                                   
  161 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 161                                           
  162 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 162                                           
      IF(D3.LT.-CN(K))GO TO 168                                          
  163 D1=D1+H1M                                                          
      IF(D1.LE.CN(K))GO TO 163                                           
  164 D1=D1-H1M                                                          
      IF(D1.GT.CN(K))GO TO 164                                           
      IF(D1.LT.-CN(K))GO TO 162                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 162                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 164                                                          
  168 CONTINUE                                                           
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
      IF(VX.GE.V3)GO TO 191                                           
  190 V3=VX                                                           
      DIJ(NP2)=V3                                                     
      CDIJ(NP2+20)=ALPHA(LCT)                                            
      IF(V3.GE.HT1(8))GO TO 191                                          
      IF(H21I.GT.H1N)H21I=H21I-H1M                                    
      IF(H21I.LT.-H1N)H21I=H21I+H1M                                   
      HT1(1)=H1M                                                         
      HT1(2)=H21I                                                        
      HT1(3)=H2I                                                         
      HT1(4)=H3M                                                      
      HT1(5)=H1A3                                                        
      HT1(6)=H2A3                                                        
      HT1(7)=0.                                                          
      HT1(8)=V3                                                          
      CHT1(9)=ALPHA(LCT)                                                 
      HT1(10)=A1                                                         
      HT1(11)=A2                                                         
      HT1(12)=A3                                                         
      HT1(13)=0.                                                         
      HT1(14)=0.                                                         
      NLOW=5                                                             
      MLOW=LCT+5                                                      
  191 IF(IPR.LT.2)GO TO 198                                           
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
!-----Which space group has been completed?                         
  198 IF(LCT.GT.1)GO TO 210                                           
!-----P2 completed, go to next designated space group in ZTRTB      
      IF(ITR(7)-2)200,300,1000                                        
!-----Space group P21 - Code AF                                     
  200 LCT=2                                                           
!-----Set to use distances to atoms displaced to H3N - screw axis   
      ICT=-1                                                          
      ER=ERMF                                                         
      H3A3=H3N                                                        
      GO TO 105                                                       
!-----P21 completed, go to next designated space group              
  210 IF(ITR(8).GT.1)GO TO 1000                                       
!-----Space group Pc - Code AG                                      
!     Calculate standoff of grid of P molecules in axis-2 direction 
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
!-----Collect I molecule positions for double standoff              
      X1I=2.*X1P3                                                     
      DO 329 I=1,NMOD                                                    
      DO 328 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D1=X1I+W(1,J)-W(1,I)                                               
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
      TI(N+2)=W(2,J)-W(2,I)                                              
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
      D1=X1P3+W(1,J)-W(1,I)                                              
      D3=X3P3-W(3,J)-W(3,I)                                              
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
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 334                                                          
  338 CONTINUE                                                           
  339 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                          
      NEND=N-1                                                           
      CALL MINHI(X2P3,0.,H2P3,H2M,KE)                                    
      IF(KCT)2903,342,341                                              
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
      X1P3=HH1P3+DX1P3                                                 
      GO TO 315                                                        
  352 MCT=6                                                            
      X3P3=HH3P3-DX3P3                                                 
      X1P3=HH1P3                                                       
      GO TO 315                                                        
  353 MCT=7                                                            
      X3P3=HH3P3-DX3P3                                                 
      X1P3=HH1P3-DX1P3                                                 
      GO TO 315                                                        
  354 MCT=8                                                            
      X3P3=HH3P3                                                       
      X1P3=HH1P3-DX1P3                                                 
      GO TO 315                                                        
  355 MCT=9                                                            
      X3P3=HH3P3+DX3P3                                                 
      X1P3=HH1P3-DX1P3                                                 
      GO TO 315                                                        
  357 IF(JE.LT.KE)GO TO 347                                            
      IF(JE.GE.NNE)GO TO 358                                           
      KE=NNE                                                           
      GO TO 347                                                        
  358 CONTINUE                                                         
      H21I=2.*H1P3                                                    
      H2I=2.*H2P3                                                     
      VX=H1N*H2I*H3M                                                  
      IF(V3-VX.LT..01)GO TO 360                                          
      V3=VX                                                              
      DIJ(NP2)=V3                                                        
      CDIJ(NP2+20)=ALPHA(3)                                              
      IF(V3.GE.HT1(8))GO TO 360                                          
      HT1(1)=H1M                                                         
      HT1(2)=H21I                                                        
      HT1(3)=H2I                                                         
      HT1(4)=H3M                                                      
      HT1(5)=0.                                                          
      HT1(6)=0.                                                          
      HT1(7)=H3P3                                                        
      HT1(8)=V3                                                          
      CHT1(9)=ALPHA(3)                                                   
      HT1(10)=A1                                                         
      HT1(11)=A2                                                         
      HT1(12)=A3                                                         
      HT1(13)=0.                                                         
      HT1(14)=0.                                                         
      NLOW=6                                                             
      MLOW=8                                                          
  360 IF(IPR.LT.2)GO TO 1000                                          
      A(1)=H1M                                                        
      A(2)=SQRT(H21I**2+H2I**2)                                       
      A(3)=H3M                                                        
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
      D=SQRT(H1P3**2+H2P3**2+H3P3**2)                                 
      C(4)=90.                                                        
      C(5)=90.                                                        
      C(6)=57.296*ACOS(C(3))                                          
      CALL PAGE(2,2,0)                                                
      WRITE(NQ,195)ALPHA(3),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)     
      WRITE(NQ,366)H1P3,H2P3,H3P3,(A(I),I=4,6),D                      
  366 FORMAT(1X,'P3 AT',3F7.3,';  ',3F7.4,';  D =',F6.2)                
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT               
  903 FORMAT(26H0STORAGE EXCEEDED BY ZTRTB,5X,3(I9,I3),I9)               
      GO TO 999                                                          
 1903 MPL=1                                    
      GO TO 904                          
 2903 MPL=2                        
  904 CALL PAGE(2,2,0)                                                
      WRITE(NQ,905)MPL                                                
  905 FORMAT('0PARAMETER TOO SMALL - MPL =',I2)                       
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
      END SUBROUTINE ZTRTB     

