      SUBROUTINE ZFRTC    
!                                                                                                                         
!     For AN and AO.
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
      CHARACTER(2) :: ALPHA(2) = (/'AN','AO'/)
      CHARACTER(2) :: SYM(3) = (/'C ','A3','P3'/)
!                                  
!      DATA ALPHA/'AN','AO'/                                                 
!      DATA SYM/'C ','A3','P3'/
!
      INTEGER :: I, ICT, J, JCT, JE, K, KE, LCT, LCTM, MCT, N
!
      REAL :: A(12), C(6)
      REAL :: D, D1, D2, D3, DY1A3
      REAL :: H1A3, H1C, H1I, H1IH 
      REAL :: H2A3, H2C, H2I, H21I
      REAL :: H3A3, H3C, H3I, H3IH, H3P3, HL 
      REAL :: SM(3), SINE3
      REAL :: VX, V3S
      REAL :: X1A3, X2A3, X3P3
      REAL :: Y1A3, Y2A3
!
!-----IDENTITIES IN CONTACT ALONG AXIS-1                                    
!-----PLANES IN CONTACT WITH CENTRAL MOLECULE ALONG AXIS-3                  
!                                                                           
      IF(IPR.LT.2)GO TO 10                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)ITR(15),ITR(16)                                         
    5 FORMAT(1X,'ZFRTC called',2I3)                                      
!-----Set starting value of V3                                           
   10 V3S=V3                                                             
!-----DETERMINE PLANE-1,3 OF PLANE AND IDENTITY MOLECULES                   
      KE=NNE                                                                
      NARK=0                                                             
      N=MARK2                                                               
      NSTP2=5                                                            
!-----Collect close distances along axis-2                               
      DO 21 I=1,NMOD                                                        
      DO 20 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      D2=W(2,J)-W(2,I)                                                      
      IF(ABS(D2).GT.CN(K))GO TO 20                                          
      IT(N)=K                                                               
      TI(N+1)=D2**2                                                         
      TI(N+2)=W(1,J)-W(1,I)                                                 
      TI(N+3)=W(3,J)-W(3,I)                                                 
      TI(N+4)=-W(3,J)-W(3,I)                                                
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
   20 CONTINUE                                                              
   21 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
!-----Calculate P molecule offsets along positive and negative axis-3    
!-----CALC. + AND - PLANE OFFSETS ALONG B AXIS                              
      DO 25 I=MARK2,NEND2,NSTP2                                             
      K=IT(I)                                                               
      IF(ABS(TI(I+2)).GT.CN(K))GO TO 25                                     
      D=TI(I+1)+TI(I+2)**2                                                  
      IF(D.GT.CN(K+10))GO TO 25                                             
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+4)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   25 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      H3P3=H3M                                                              
      ER=ERM                                                                
      CALL MINHI(H3P3,0.,100.,H3M,KE)                                       
      IF(KCT)904,26,906                                                  
   26 X3P3=-H3M                                                             
      CALL MINHI(X3P3,0.,-100.,-H3M,KE)                                     
      IF(KCT)904,27,906                                                  
   27 H3I=H3P3-X3P3                                                         
      H3IH=H3I/2.                                                        
      IF(IPR.LT.3)GO TO 30                                                  
      CALL PAGE(2,2,0)                                                      
      WRITE(NQ,28)H3P3,X3P3                                                 
   28 FORMAT(7H0H3P3 =,F7.3,8H, X3P3 =,F7.3)                                
!-----CALC. AXIS-1 LENGTH WITH BOTH IDENTITY AND PLANE CONTACTS             
   30 N=MARK                                                                
      DO 36 I=MARK2,NEND2,NSTP2                                             
      K=IT(I)                                                               
      D3=TI(I+3)                                                            
   31 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 31                                               
   32 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 32                                               
      IF(D3.LT.-CN(K))GO TO 33                                              
      D=TI(I+1)+D3**2                                                       
      IF(D.GT.CN(K+10))GO TO 32                                             
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 32                                                              
   33 D3=TI(I+4)+H3P3                                                       
   34 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 34                                               
   35 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 35                                               
      IF(D3.LT.-CN(K))GO TO 36                                              
      D=TI(I+1)+D3**2                                                       
      IF(D.GT.CN(K+10))GO TO 35                                             
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 35                                                              
   36 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      H1I=7.                                                                
      CALL MINHI(H1I,0.,100.,7.,KE)                                         
      IF(KCT)904,38,906                                                  
   38 H1IH=H1I/2.                                                        
      IF(IPR.LT.3)GO TO 40                                                  
      WRITE(NQ,39)H1I                                                       
   39 FORMAT(' H1I =',F7.3)                                                 
!-----Calculate axis-2 altitude and angle from offset of grid of A and   
!     C molecules in positive and negative directions                    
!-----CALC. AXIS-2 HALF-HIGHTH FROM OFFSET OF AXIS-CENTER PLANE             
   40 IF(ITR(15)-2)41,42,1000                                               
!-----Screw axis - space group P21/m                                     
   41 H3A3=H3IH                                                          
      LCT=1                                                                 
      GO TO 45                                                              
!-----Twofold axis - space group P2/m                                    
   42 H3A3=0.                                                               
      LCT=2                                                                 
!-----Collect short distances along axis-3                               
   45 N=MARK2                                                               
      DO 51 I=1,NMOD                                                        
      DO 50 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
!-----Collect A molecule distances first                                 
      ICT=0                                                                 
      D3=H3A3+W(3,J)-W(3,I)                                                 
   47 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 47                                               
   48 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 48                                               
      IF(D3.LT.-CN(K))GO TO 49                                              
      IT(N)=K                                                               
      TI(N+1)=D3**2                                                         
      TI(N+2)=-W(1,J)-W(1,I)                                                
      TI(N+3)=-W(2,J)-W(2,I)                                                
      IT(N+4)=K+10                                                          
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 48                                                              
!-----Have C molecule distances been collected?                          
   49 IF(ICT.NE.0)GO TO 50                                                  
!-----Switch to C molecules - the vector from an A molecule to a C       
!     molecule is that from the origin I to a P molecule                 
      D3=H3A3+H3P3-W(3,J)-W(3,I)                                            
      ICT=1                                                                 
      GO TO 47                                                              
   50 CONTINUE                                                              
   51 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
!-----Set up for calculation of offset in positive direction             
      H2A3=H2N                                                              
      HL=H2M                                                                
      X2A3=100.                                                             
      JCT=0                                                              
      ER=ERMF                                                               
  110 JE=1                                                                  
      DY1A3=H1I/NV                                                          
      Y1A3=-H1IH                                                            
      MCT=0                                                                 
      KE=NE                                                                 
  116 N=MARK                                                                
      DO 129 I=MARK2,NEND2,NSTP2                                            
      K=IT(I)                                                               
      D1=TI(I+2)+Y1A3                                                       
  121 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 121                                              
  122 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 122                                              
      IF(D1.LT.-CN(K))GO TO 129                                             
      D=TI(I+1)+D1**2                                                       
      IF(D.GT.CN(K+10))GO TO 122                                            
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 122                                                             
  129 CONTINUE                                                              
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
      GO TO 116                                                          
  154 MCT=1                                                              
      JE=2*JE                                                            
      DY1A3=H1I/JE                                                       
      Y1A3=X1A3+DY1A3                                                    
      GO TO 116                                                          
  155 MCT=-1                                                             
      Y1A3=X1A3-DY1A3                                                    
      IF(Y1A3.LT.-H1IH)Y1A3=Y1A3+H1I                                     
      GO TO 116                                                          
  156 IF(JE.LT.KE)GO TO 154                                              
      IF(JE.GE.NNE)GO TO 158                                             
      KE=NNE                                                             
      GO TO 154                                                          
  158 CONTINUE                                                               
!-----Has negative standoff been calculated?                             
      IF(JCT.NE.0)GO TO 159                                              
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
!     Set up for negative standoff calculation                           
      X2A3=-100.                                                         
      Y2A3=-H2N                                                          
      HL=-H2M                                                            
      JCT=1                                                              
      GO TO 110                                                          
  159 H21I=H1A3-X1A3                                                     
      H2I=H2A3-X2A3                                                      
      VX=H1I*H2I*H3I/4.                                                     
      H1C=H1A3                                                           
      H2C=H2A3                                                           
      H3C=H3A3+H3P3                                                      
      IF(H3C.GT.H3IH)H3C=H3C-H3I                                         
      IF(H3C.LT.-H3IH)H3C=H3C+H3I                                        
      IF(VX.GE.V3)GO TO 170                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 160                                             
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
      NLOW=11                                                               
      MLOW=LCT+14                                                        
  160 IF(IPR.LT.2)GO TO 170                                              
      LCTM=LCT                                                           
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
      A(8)=0.                                                            
      A(7)=0.                                                            
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
  170 IF(LCT.NE.1)GO TO 400                                                 
      IF(ITR(16).EQ.1)GO TO 42                                              
  400 IF(IPR.LT.2)GO TO 1000                                             
      IF(V3.GE.V3S)GO TO 1000                                            
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,405)ALPHA(LCTM),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)     
  405 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3    &                 
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
  406 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)
      SM(1)=A(4)-A(4)
      SM(2)=A(5)-A(5)
      SM(3)=A(6)-A(6)
      WRITE(NQ,406)SYM(1),(A(I),I=4,6),SM
      SM(1)=A(10)-A(4)
      SM(2)=A(11)-A(5)
      SM(3)=A(12)
      WRITE(NQ,406)SYM(2),(A(I),I=10,12),SM
      SM(1)=A(7)
      SM(2)=A(8)
      SM(3)=A(9)-A(6)
      WRITE(NQ,406)SYM(3),(A(I),I=7,9),SM
      GO TO 1000
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZFRTC,5X,3(I9,I3),I9)                  
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
  921 FORMAT(16H0NO INTERACTIONS)                                           
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZFRTC  
