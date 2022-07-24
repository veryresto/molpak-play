      SUBROUTINE ZORTB  
!                                                                        
!     For AR, AS, AT, AU, AV, AW, AX, AY                                                  
!                                                                        
      USE molpakCommonMod
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
!      DATA ALPHA/'AR','AS','AT','AU','AV','AW','AX','AY'/                  
!      DATA SYM/'P1','A1','P2','A3','P3'/ 
!
      CHARACTER(2) :: ALPHA(8) = (/'AR','AS','AT','AU','AV','AW',&
     &                             'AX','AY'/)
      CHARACTER(2) :: SYM(5) = (/'P1','A1','P2','A3','P3'/)
!
      INTEGER :: I, ICT, IE, J, JCT, JE, JPT, K, KE, KK, KPT, L, LCT, LCTM
      INTEGER :: MCT, MZ2P1, M3P2, N, N3P1, N3P2, NXZG 
!
      REAL :: D, D1, D2, D3, DY2P2, DZ1P1, DZ3P3, E
      REAL :: H1I, H1P1, H1P2, H1P3, H2I, H2P1, H2P2, H2P3, H2P2R, H2P2S 
      REAL :: H3A3, H3P1, H3P2, H3P3, HL
      REAL :: ORG1, ORG2, ORG3 
      REAL :: SM(6), VX, VY, V3S, W2P1, W2P3
      REAL :: X1I, X1P1, X1P2, X1P3, X2I, X2P1, X2P2, X2P3 
      REAL :: X3A3, X3P1, X3P2, X3P3 
      REAL :: Y1I, Y1P1, Y1P2, Y1P3, Y2I, Y2P1, Y2P2, Y2P3, Y3P3
      REAL :: Z1P1, Z1P2, Z2P1, Z2P3, Z3P3
!                               
      IF(IPR.LT.2)GO TO 15                                               
      WRITE(NQ,10)(ITR(I),I=19,26)                                       
   10 FORMAT(' ZORTB called',8I3)                                        
!-----TWO PLANE ONE AXIS SPACE GROUPS                                    
!     Only space groups which contain a P2 with an axis-1 glide and a P1 
!      or P3 with an axis-2 glide are included                           
!                                                                        
   15 NARK=0                                                             
!-----Set starting value of V3                                           
      V3S=V3                                                             
   40 LCT=0                                                                 
      M3P2=0                                                                
   41 LCT=LCT+1                                                             
!-----TRANSFER TO VARIOUS SUB-CLASSES                                       
      GO TO(45,50,55,60,65,70,75,80,300),LCT                             
!-----P2 - AXIS-3 GLIDE, P1 - AXIS-3 GLIDE, A3 TWOFOLD - PNN2            
   45 GO TO(46,51,56,61,66,71,76,81,300),ITR(19)                         
!     Signal collection of axis-3 distances displaced by H3N - axis-3    
!      glide of P2 - first four space groups                             
   46 N3P2=-1                                                               
      X3P2=H3N                                                              
!     Signal axis-3 distances for axis-3 glide of P1                     
      N3P1=-1                                                               
      X3P1=H3N                                                              
      X3A3=0.                                                            
      GO TO 100                                                             
!-----P2 - AXIS-3 GLIDE, P1 - NO AXIS-3 GLIDE, A3 SCREW - PNA2(1)        
   50 GO TO(51,56,61,66,71,76,81,300),ITR(20)                            
   51 LCT=2                                                                 
      N3P2=-1                                                               
      X3P2=H3N                                                              
!     Signal axis-3 distances for no axis-3 glide of P1                  
      N3P1=1                                                                
      X3P1=0.                                                               
      X3A3=H3N                                                           
      GO TO 100                                                             
!-----P2 - AXIS-3 GLIDE, P3 - AXIS-1 GLIDE, A1 TWOFOLD - PNN2            
   55 GO TO(56,61,66,71,76,81,300),ITR(21)                               
   56 LCT=3                                                                 
      N3P2=-1                                                               
      X3P2=H3N                                                              
      NXZG=0                                                                
!     Multiplication factor to produce axis-1 glide of P3                
      F1P3=1.                                                               
      F1A1=0.                                                            
      GO TO 100                                                             
!-----P2 - AXIS-3 GLIDE, P3 - NO AXIS-1 GLIDE, A1 SCREW - PNA2(1)        
   60 GO TO(61,66,71,76,81,300),ITR(22)                                  
   61 LCT=4                                                                 
      N3P2=-1                                                               
      X3P2=H3N                                                              
      NXZG=0                                                                
!     Multiplication factor for no axis-1 glide of P3                    
      F1P3=0.                                                               
      F1A1=1.                                                            
      GO TO 100                                                             
!-----P2 - NO AXIS-3 GLIDE, P1 - AXIS-3 GLIDE, A3 SCREW - PNA2(1)        
   65 GO TO(66,71,76,81,300),ITR(23)                                     
   66 LCT=5                                                                 
!     Signal collection of axis-3 distances NOT displaced by H3N - NO    
!      axis-3 glide of P2 - last four space groups                       
      N3P2=1                                                                
      X3P2=0.                                                               
!     Signal axis-3 distances for axis-3 glide of P1                     
      N3P1=-1                                                               
      X3P1=H3N                                                              
      X3A3=H3N                                                           
      GO TO 100                                                             
!-----P2 - NO AXIS-3 GLIDE, P1 - NO AXIS-3 GLIDE, A3 TWOFOLD - PBA2      
   70 GO TO(71,76,81,300),ITR(24)                                        
   71 LCT=6                                                                 
!     Signal axis-3 distances for NO axis-3 glide of P1                  
      N3P2=1                                                                
      X3P2=0.                                                               
      N3P1=1                                                                
      X3P1=0.                                                               
      X3A3=0.                                                            
      GO TO 100                                                             
!-----P2 - NO AXIS-3 GLIDE, P3 - AXIS-1 GLIDE, A1 TWOFOLD - PNC2         
   75 IF(ITR(25)-2)76,81,300                                             
   76 LCT=7                                                                 
      N3P2=1                                                                
      X3P2=0.                                                               
      NXZG=0                                                                
!     Multiplication factor to produce axis-1 glide of P3                
      F1P3=1.                                                               
      F1A1=0.                                                            
      GO TO 100                                                             
!-----P2 - NO AXIS-3 GLIDE, P3 - NO AXIS-1 GLIDE, A1 SCREW - PCA2(1)     
   80 IF(ITR(26)-2)81,300,300                                            
   81 LCT=8                                                                 
      N3P2=1                                                                
      X3P2=0.                                                               
      NXZG=0                                                                
!     Multiplication factor for no axis-1 glide of P3                    
      F1P3=0.                                                               
      F1A1=1.                                                            
!     Have the axis-2 offset limits for a P2 molecule with this axis-3   
!      been calculated and stored?                                       
  100 IF(M3P2.EQ.N3P2)GO TO 110                                             
!     Calculate axis-2 offset limits for P2 molecule by determination of 
!     the positions at which X1I would have its minimum value, H1M       
      M3P2=N3P2                                                             
      KE=NNE                                                                
      ER=ERMT                                                               
      MARK=MARK2                                                            
  105 N=MARK                                                                
      DO 106 I=MARK1,NEND1,NSTP1                                            
      K=N3P2*IT(I)                                                          
      IF(K.LE.0)GO TO 106                                                   
      D1=TI(I+2)+H1N                                                        
      IF(ABS(D1).GT.CN(K))GO TO 106                                         
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 106                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+6)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
  106 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      H2P2R=.75*H2M                                                         
      CALL MINHI(H2P2R,0.,100.,H2M,KE)                                      
      IF(KCT)904,107,906                                                 
  107 H2P2S=-H2P2R                                                       
      CALL MINHI(H2P2S,0.,-100.,-H2M,KE)                                 
      IF(KCT)904,108,906                                                 
  108 H2P2R=H2P2R-H2P2S                                                  
      IF(IPR.LT.3)GO TO 110                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,109)H2P2R,H2P2S                                              
  109 FORMAT(15X,7HH2P2R =,F7.3,9H, H2P2S =,F7.3)                           
!-----TOP OF P2 PLACEMENT LOOP                                              
  110 VX=1000000.         ! 03-15-06                                                   
      KE=NE                                                                 
      ICT=0                                                                 
      IE=1                                                                  
      DY2P2=H2P2R/NV                                                        
      Y2P2=H2P2S                                                            
  111 MARK=MARK2                                                            
      ER=ERMT                                                               
      N=MARK                                                                
      DO 113 I=MARK1,NEND1,NSTP1                                            
      K=N3P2*IT(I)                                                          
      IF(K.LE.0)GO TO 113                                                   
      D2=TI(I+6)+Y2P2                                                       
      IF(ABS(D2).GT.CN(K))GO TO 113                                         
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 113                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
  113 CONTINUE                                                              
      IF(N.GT.MARK)GO TO 114                                                
      Y1P2=H1N                                                              
      GO TO 115                                                             
  114 NEND=N-1                                                              
      CALL MINHI(Y1P2,H1N,100.,H1M,KE)                                      
  115 Y1I=2.*Y1P2                                                           
      IF(IPR.LT.3)GO TO 118                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,117)Y1P2,Y2P2                                             
  117 FORMAT(1X,'Y1P2 =',F8.3,';  Y2P2 =',F8.3)                          
  118 GO TO(200,200,119,119,200,200,119,119),LCT                            
!-----TOP OF P3 PLACEMENT LOOP                                           
  119 Y1P3=F1P3*Y1P2                                                        
!     Set phony value of X1P1 to be used in the HT1 array                
!      There is no P1 in this space group                                
      X1P1=0.                                                            
      N=MARK2                                                               
      NSTP2=5                                                            
      DO 125 I=1,NMOD                                                       
      DO 124 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
!-----Set up for P3 molecule interactions                                
      KK=K                                                                  
      D=1.                                                                  
      E=0.                                                                  
      D1=Y1P3+W(1,J)-W(1,I)                                                 
  120 D1=D1+Y1I                                                             
      IF(D1.LE.CN(K))GO TO 120                                              
  121 D1=D1-Y1I                                                             
      IF(D1.GT.CN(K))GO TO 121                                              
      IF(D1.LT.-CN(K))GO TO 122                                             
      IT(N)=KK                                                              
      TI(N+1)=D1**2                                                         
      TI(N+2)=D*W(2,J)-W(2,I)+E                                             
      TI(N+3)=-W(3,J)-W(3,I)                                                
      IT(N+4)=L                                                             
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 121                                                             
!-----Have A1 molecule interactions been calculated?                     
  122 IF(KK.LE.0)GO TO 124                                                  
      KK=-K                                                                 
      D=-1.                                                                 
!     The vector from a P3 molecule to an A1 molecule is equal to that   
!      from the origin to a P2 molecule                                  
      E=Y2P2                                                                
      D1=D1+Y1P2                                                            
      GO TO 120                                                             
  124 CONTINUE                                                              
  125 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
!-----TOP OF P3 PLACEMENT LOOP                                           
      DZ3P3=H3M/NV                                                          
      Z3P3=-H3N                                                             
      Y2P3=100.                                                             
      Z2P3=1.5*H2N                                                          
      HL=H2M                                                                
      ER=ERMF                                                               
      JE=1                                                                  
      MCT=0                                                                 
  130 N=MARK                                                                
      DO 135 I=MARK2,NEND2,NSTP2                                            
      K=IABS(IT(I))                                                         
      D3=TI(I+3)+Z3P3                                                       
!-----A negative value means that this is an A1 interaction, apply       
!      vector from origin to the P2 molecule                             
      IF(IT(I).LT.0)D3=D3+X3P2                                              
  131 D3=D3+H3M                                                             
      IF(D3.LE.CN(K))GO TO 131                                              
  132 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 132                                              
      IF(D3.LT.-CN(K))GO TO 135                                             
      D=TI(I+1)+D3**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 132                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 132                                                             
  135 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 140                                                
      NEND=N-1                                                              
      CALL MINHI(Z2P3,H2N,Y2P3,H2M,KE)                                      
      IF(KCT)139,136,152                                                    
!-----Check negative axis-2 direction because of A1 molecules            
  136 W2P3=-Z2P3                                                         
      CALL MINHI(W2P3,-Z2P3,-Y2P3,-H3M,KE)                               
      IF(KCT)151,137,152                                                 
  137 Z2P3=-W2P3                                                         
      GO TO 151                                                          
  139 W2P3=-H2N                                                          
      CALL MINHI(W2P3,-H2N,-Y2P3,-H3M,KE)                                
      IF(KCT)140,137,152                                                 
!-----Y2P3 has its minimum possible value, H2N, terminate search         
  140 Y2P3=H2N                                                              
      Y3P3=Z3P3                                                             
      GO TO 158                                                             
  151 Y3P3=Z3P3                                                          
      Y2P3=Z2P3                                                          
      Y2I=2.*Y2P3
      IF(MCT)156,153,156                                                 
  152 IF(MCT)156,153,155                                                 
  153 IF(JE.GE.NV)GO TO 154                                              
      JE=JE+1                                                            
      Z3P3=Z3P3+DZ3P3                                                    
      GO TO 130                                                          
  154 MCT=1                                                              
      JE=2.*JE                                                           
      DZ3P3=H3M/JE                                                       
      Z3P3=Y3P3+DZ3P3                                                    
      GO TO 130                                                          
  155 MCT=-1                                                             
      Z3P3=Y3P3-DZ3P3                                                    
      GO TO 130                                                          
  156 IF(JE.LT.KE)GO TO 154                                              
  158 GO TO 259                                                          
!-----TOP OF P1 PLACEMENT LOOP                                           
  200 DZ1P1=Y1I/NV                                                       
      Z1P1=-Y1P2                                                         
      Y1P1=100.                                                          
!     Set phony value X3P3 to be used in the HT1 array.                  
!      There is no P3 in this space group.                               
      X3P3=0.                                                            
      Z2P1=AMAX1(H2M,H2M-Y2P2)                                           
      Y2P1=100.                                                          
      ER=ERMT                                                               
      JE=1                                                                  
      MCT=0                                                                 
      MARK=MARK2                                                            
  216 N=MARK                                                                
      DO 225 I=MARK1,NEND1,NSTP1                                            
      JCT=0                                                              
      K=N3P1*IT(I)                                                          
!-----Use axis-3 glide or non-glide P1 molecules as called for           
      IF(K.LE.0)GO TO 223                                                   
      D1=TI(I+5)+Z1P1                                                       
      D2=TI(I+3)                                                            
  221 D1=D1+Y1I                                                             
      IF(D1.LE.CN(K))GO TO 221                                              
  222 D1=D1-Y1I                                                             
      IF(D1.GT.CN(K))GO TO 222                                              
      IF(D1.LT.-CN(K))GO TO 223                                             
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 222                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 222                                                             
!-----Have A3 interactions been considered?                              
  223 IF(JCT.NE.0)GO TO 225                                                 
      K=N3P2*K                                                              
!     Use either screw or twofold A3 molecules as called for - an axis-3 
!      glide by either P2 or P1 calls for a twofold A3 - a glide by      
!      neither or both calls for a screw A3                              
      IF(K.LE.0)GO TO 225                                                    
!     The vector from a P1 to an A3 is equal to that from the origin to  
!      a P2 molecule                                                     
      D1=D1+Y1P2                                                            
      D2=TI(I+6)+Y2P2                                                       
      JCT=1                                                              
      GO TO 221                                                             
  225 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 240                                                
      NEND=N-1                                                              
      ER=ERMT
      CALL MINHI(Z2P1,H2N,Y2P1,H2M,KE)                                      
      MZ2P1=0
      IF(KCT)229,230,252                                                    
!-----In the positive axis-2 direction, Z2P1 has its minimum value       
  229 MZ2P1=1
!-----Check negative axis-2 direction because of A3 molecules            
  230 W2P1=AMIN1(-H2M,-H2M-Y2P2)                                         
      CALL MINHI(W2P1,-Z2P1,-Y2P1,-H2M,KE)                               
      IF(KCT)238,239,252                                                 
  238 IF(MZ2P1)251,251,240
  239 Z2P1=-W2P1                                                         
      GO TO 251                                                          
!-----Z2P1 has its minimum possible value, H2N, terminate search         
  240 Y2P1=H2N                                                              
      Y1P1=Z1P1                                                             
      Y2I=2.*Y2P1             
      GO TO 258                                                             
  251 Y1P1=Z1P1                                                          
      Y2P1=Z2P1                                                          
      Y2I=2.*Y2P1
      IF(MCT)256,253,256                                                 
  252 IF(MCT)256,253,255                                                 
  253 IF(JE.GE.NV)GO TO 254                                              
      JE=JE+1                                                            
      Z1P1=Z1P1+DZ1P1                                                    
      GO TO 216                                                          
  254 MCT=1                                                              
      JE=2.*JE                                                           
      DZ1P1=Y1I/JE                                                       
      Z1P1=Y1P1+DZ1P1                                                    
      GO TO 216                                                          
  255 MCT=-1                                                             
      Z1P1=Y1P1-DZ1P1                                                    
      IF(Z1P1+Y1P2)256,256,216                                           
  256 IF(JE.LT.KE)GO TO 254                                              
  258 CONTINUE                                  
      VY=Y1P2*Y2P1*H3M                                                   
      IF(VY.GE.VX)GO TO 285                                              
      JPT=0
      GO TO 260
  259 VY=Y1P2*Y2P3*H3M                                                      
      IF(VY.GE.VX)GO TO 285                                              
      JPT=1
!-----Check all locations of molecule P2                                 
  260 MARK=MARK2
	N=MARK
      DO 269 I=1,NMOD
      DO 268 J=1,NMOD
	K=IA(I)+IAA(J)
      L=K+10
	D2=Y2P2-W(2,J)-W(2,I)
	D3=X3P2+W(3,J)-W(3,I)
  261 D2=D2+Y2I
      IF(D2.LT.CN(K))GO TO 261
  262 D2=D2-Y2I
      IF(D2.GT.CN(K))GO TO 262
	IF(D2.LT.-CN(K))GO TO 268
  263 D3=D3+H3M
      IF(D3.LT.CN(K))GO TO 263
  264 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 264
	IF(D3.LT.-CN(K))GO TO 262
	D=D2**2+D3**2
	IF(D.GT.CN(L))GO TO 264
	IT(N)=K
	TI(N+1)=D
	TI(N+2)=W(1,J)-W(1,I)
	N=N+NSTP
	IF(N.GT.LIMIT)GO TO 902
	GO TO 264
  268 CONTINUE
  269 CONTINUE
      IF(N.LE.MARK)GO TO 920
      NEND=N-1
	Z1P2=Y1P2
	CALL MINHI(Z1P2,Y1P2,100.,H1M,NNE)
      KPT=KCT
      IF(KCT)271,270,906
  270 Y1P2=Z1P2
      Y1I=2.*Y1P2
  271 Z1P2=-Y1P2
      CALL MINHI(Z1P2,-Y1P2,-100.,-H1M,NNE)
	IF(KCT)272,273,906
  272 IF(KPT.EQ.0)GO TO 274
      IF(JPT)281,276,281
  273 Y1P2=-Z1P2
      Y1I=2.*Y1P2
  274 IF(JPT)280,275,280
  275 VY=Y1P2*Y2P1*H3M
      IF(VY.GE.VX)GO TO 285 
  276 X1P1=Y1P1                                                          
      X2P1=Y2P1                                                          
      X2I=Y2I                                                            
      GO TO 287                                                             
  280 VY=Y1P2*Y2P3*H3M
      IF(VY.GE.VX)GO TO 285
  281 X1P3=Y1P3                                                          
      X2P3=Y2P3                                                          
      X3P3=Y3P3                                                          
      X2I=Y2I                                                            
      GO TO 287                                                          
  285 IF(ICT)291,293,286                                                    
  286 ICT=-1                                                                
      Y2P2=X2P2-DY2P2                                                       
      IF(Y2P2-H2P2S)291,291,111                                          
  287 X1P2=Y1P2                                                          
      X2P2=Y2P2                                                          
      X1I=Y1I                                                            
      VX=VY                                                                 
      IF(IPR.LT.3)GO TO 290                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,289)X1P1,X2P1,X1P2,X2P2,X1P3,X2P3,VX,NEND                 
  289 FORMAT(1X,6F8.3,' VX =',F9.3,' NEND =',I6)                         
  290 IF(ICT.EQ.0)GO TO 293                                                 
  291 IF(IE.GE.KE)GO TO 294                                                 
  292 IE=2*IE                                                               
      DY2P2=H2P2R/IE                                                        
      ICT=1                                                                 
      Y2P2=X2P2+DY2P2                                                       
      GO TO 111                                                             
  293 IF(IE.GE.NV)GO TO 292                                                 
      IE=IE+1                                                               
      Y2P2=Y2P2+DY2P2                                                       
      GO TO 111                                                             
  294 IF(KE.GE.NNE)GO TO 295                                                
      KE=NNE                                                                
      GO TO 293                                                             
  295 IF(VX.GE.V3)GO TO 41                                                  
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H2I=X2I                                                            
      H1P1=X1P1                                                          
      H2P1=X2P1                                                          
      H3P1=X3P1                                                          
      H1P2=X1P2                                                          
      H2P2=X2P2                                                          
      H3P2=X3P2                                                          
      H1P3=X1P3                                                          
      H2P3=X2P3                                                          
      H3P3=X3P3                                                          
      H3A3=X3A3                                                          
      GO TO(298,298,297,297,298,298,297,297),LCT                         
  297 ORG1=H1P1                                                          
      ORG2=H2P2+H2P3                                                     
      ORG3=H3P2+H3P3                                                     
      GO TO 299                                                          
  298 ORG1=H1P1+H1P2                                                     
      ORG2=H2P1+H2P2                                                     
      ORG3=H3P3                                                          
  299 IF(V3.GE.HT1(8))GO TO 41                                              
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                            
      HT1(5)=ORG1                                                           
      HT1(6)=ORG2                                                           
      HT1(7)=ORG3                                                           
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=13                                                               
      MLOW=LCT+18                                                           
      GO TO 41                                                           
  300 IF(IPR.LT.2)GO TO 1000                                                
      IF(V3.GE.V3S)GO TO 1000                                            
      CALL PAGE(4,4,0)                                                      
      WRITE(NQ,301)ALPHA(LCTM),A1,A2,A3,VX,H1I,H2I,H3M                   
  301 FORMAT(1X,A2,3F7.1,'  V =',F8.2,'  AXES',3F8.3)                    
      DO 319 I=1,5                                                       
      GO TO (302,304,306,308,310),I                                      
  302 GO TO (303,303,319,319,303,303,319,319),LCTM                       
  303 SM(1)=H1P1                                                         
      SM(2)=H2P1                                                         
      SM(3)=H3P1                                                         
      SM(4)=SM(1)-ORG1                                                   
      SM(5)=SM(2)                                                        
      SM(6)=SM(3)                                                        
      GO TO 312                                                          
  304 GO TO (319,319,305,305,319,319,305,305),LCTM                       
  305 SM(1)=H1P2+H1P3                                                    
      SM(2)=H2P2+H2P3                                                    
      SM(3)=H3P2+H3P3                                                    
      SM(4)=SM(1)                                                        
      SM(5)=SM(2)-ORG2                                                   
      SM(6)=SM(3)-ORG3                                                   
      GO TO 312                                                          
  306 SM(1)=H1P2                                                         
      SM(2)=H2P2                                                         
      SM(3)=H3P2                                                         
      SM(4)=SM(1)                                                        
      SM(5)=SM(2)-ORG2                                                   
      SM(6)=SM(3)                                                        
      GO TO 312                                                          
  308 GO TO (309,309,319,319,309,309,319,319),LCTM                       
  309 SM(1)=H1P1+H1P2                                                    
      SM(2)=H2P1+H2P2                                                    
      SM(3)=H3P1+H3P2                                                    
      SM(4)=SM(1)-ORG1                                                   
      SM(5)=SM(2)-ORG2                                                   
      SM(6)=SM(3)                                                        
      GO TO 312                                                          
  310 GO TO (319,319,311,311,319,319,311,311),LCTM                       
  311 SM(1)=H1P3                                                         
      SM(2)=H2P3                                                         
      SM(3)=H3P3                                                         
      SM(4)=SM(1)                                                        
      SM(5)=SM(2)                                                        
      SM(6)=SM(3)-ORG3                                                   
  312 SM(4)=SM(4)/H1I                                                    
      SM(5)=SM(5)/H2I                                                    
      SM(6)=SM(6)/H3M                                                    
      IF(SM(4).GT..5)SM(4)=SM(4)-1.                                      
      IF(SM(4).LT.-.5)SM(4)=SM(4)+1.                                     
      IF(SM(5).GT..5)SM(5)=SM(5)-1.                                      
      IF(SM(5).LT.-.5)SM(5)=SM(5)+1.                                     
      IF(SM(6).GT..5)SM(6)=SM(6)-1.                                      
      IF(SM(6).LT.-.5)SM(6)=SM(6)+1.                                     
      D=SQRT(SM(1)**2+SM(2)**2+SM(3)**2)                                 
      WRITE(NQ,315)SYM(I),SM,D                                           
  315 FORMAT(1X,A2,' AT',3X,3F8.3,';  ',3F8.4,'  D =',F6.2)              
  319 CONTINUE                                                           
      SM(1)=ORG1/2.                                                      
      SM(2)=ORG2/2.                                                      
      SM(3)=ORG3/2.                                                      
      SM(4)=SM(1)/H1I                                                    
      SM(5)=SM(2)/H2I                                                    
      SM(6)=SM(3)/H3M                                                    
      WRITE(NQ,320)SM                                                    
  320 FORMAT(1X,'ORIGIN',2X,3F8.3,';  ',3F8.4)                           
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,N,LIMIT,MAXT           
  903 FORMAT(26H0STORAGE EXCEEDED BY ZORTB,5X,3(I9,I3),3I9)                 
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
      WRITE(NQ,921)KCT,MCT,N,MARK3,NARK,MARK,NEND2                          
  921 FORMAT(16H0NO INTERACTIONS,2I3,5I5)                                   
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZORTB  
