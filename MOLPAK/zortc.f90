      SUBROUTINE ZORTC      
!                                                                        
!     For AZ, BA, BB, BC, BD, BE, BF, BG, BH, BI, BJ, BK
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
!      COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000),
!     +       CODE(5000)
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
      CHARACTER(2) :: ALPHA(12) = (/'AZ','BA','BB','BC','BD','BE', &
     &                              'BF','BG','BH','BI','BJ','BK'/) 
      CHARACTER(2) :: SYM(5) = (/'A1','A2','P1','P2','A3'/)    
!
!      DIMENSION ALPHA(12),SM(6),SYM(5)                                      
!      DATA ALPHA/'AZ','BA','BB','BC','BD','BE','BF','BG','BH','BI', &        
!     &'BJ','BK'/                                                            
!      DATA SYM/'A1','A2','P1','P2','A3'/                                 
!
      INTEGER :: I, J, JCT, JE, K, KE, L, LCT, LCTM 
      INTEGER :: MA2, MA3, MCT, MP1, MP2
      INTEGER :: N, NA2, NNP1, N3A3, N3P1 
!
      REAL :: D, D2, D3, DY2A3, DY3A2
      REAL :: H1A1, H1A2, H1A2P, H1A2N, H1A3, H1A3SN, H1A3SP, H1A3TN, H1A3TP, H1I 
      REAL :: H1P1, H1P2, H1P1N, H1P1P, H1Q 
      REAL :: H2A2, H2A3, H2A3TN, H2A3SN, H2A3SP, H2A3TP, H2P1, HN  
      REAL :: H3A2, H3A2N, H3A2P, H3A3, H3A3S, H3A3T, HL, H3P1 
      REAL :: SM(6)
      REAL :: VX, VHT, V3S
      REAL :: X1I, X1A1, X1A2, X1A3, X1P2
      REAL :: X2A1, X2A3, X2P2
      REAL :: X3A1, X3A2, X3P2
      REAL :: Y1A1, Y1A2, Y1A3, Y1P2, Y2A3, Y3A2  
!
      IF(IPR.LT.2)GO TO 15                                               
      WRITE(NQ,10)(ITR(I),I=27,38)                                       
   10 FORMAT(' ZORTC called',12I3)                                       
!-----This subroutine finds coordination spheres containing six I        
!      molecules in plane-23 and either four A2 plus four A3 molecules   
!      or four P1 and four A3 molecules.                                 
!                                                                        
!-----Set the starting values of V3 AND HT1(8)                           
   15 V3S=V3                                                             
      VHT=HT1(8)                                                    
      NARK=0                                                             
      LCTM=0                                                             
      MP2=0
   38 LCT=0                                                                 
      H1Q=H1M/4.                                                         
!     Set signal that an A2 molecule has not been placed                 
      MA2=0                                                                 
!     Set signal that a P1 molecule has not been placed                  
      MP1=0                                                                 
!     Set signal that an A3 molecule has not been placed                 
      MA3=0                                                                 
   39 LCT=LCT+1                                                             
!-----SETUP DESIGNATED SUB-CLASS                                            
      GO TO(40,45,50,55,60,65,70,75,80,85,90,95,700),LCT                 
!-----A2 - screw, A3 - screw - AZ - P212121                              
   40 GO TO(41,46,51,56,61,66,71,76,81,86,91,96,700),ITR(27)             
   41 N3P1=0                                                                
      NA2=-1                                                             
      H2A2=H2N                                                              
      N3A3=-1                                                               
      H3A3S=H3N                                                          
      GO TO 100                                                          
!-----A2 - screw, A3 - twofold -BA - P21212                              
   45 GO TO(46,51,56,61,66,71,76,81,86,91,96,700),ITR(28)                
   46 LCT=2                                                                 
      N3P1=0                                                                
      NA2=-1                                                             
      H2A2=H2N                                                              
      N3A3=1                                                                
      H3A3T=0.                                                           
      IF(MA2)250,100,100                                                 
!-----A2 - twofold, A3 - screw - BB - P21212                             
   50 GO TO(51,56,61,66,71,76,81,86,91,96,700),ITR(29)                   
   51 LCT=3                                                                 
      N3P1=0                                                                
      NA2=1                                                              
      H2A2=0.                                                               
      N3A3=-1                                                               
      H3A3S=H3N                                                          
      GO TO 100                                                          
!-----A2 - twofold, A3 - twofold - BC - P2221                            
   55 GO TO(56,61,66,71,76,81,86,91,96,700),ITR(30)                      
   56 LCT=4                                                              
      N3P1=0                                                                
      NA2=1                                                              
      H2A2=0.                                                               
      N3A3=1                                                                
      H3A3T=0.                                                           
      IF(MA2)100,100,250                                                 
!-----P1 -    axis-2 glide -    axis-3 glide, A3 - screw -BD - Pna21     
   60 GO TO(61,66,71,76,81,86,91,96,700),ITR(31)                         
   61 LCT=5                                                                 
      NNP1=1                                                             
      H2P1=H2N                                                              
      N3P1=-1                                                               
      H3P1=H3N                                                           
      N3A3=-1                                                               
      H3A3S=H3N                                                          
      GO TO 200                                                             
!-----P1 -    axis-2 glide -    axis-3 glide, A3 - twofold - BE - Pnn2   
   65 GO TO(66,71,76,81,86,91,96,700),ITR(32)                            
   66 LCT=6                                                                 
      NNP1=1                                                             
      H2P1=H2N                                                              
      N3P1=-1                                                               
      H3P1=H3N                                                           
      N3A3=1                                                                
      H3A3T=0.                                                           
      IF(MP1-1)200,250,200                                               
!-----P1 -    axis-2 glide - no axis-3 glide, A3 - screw - BF - Pna21    
   70 GO TO(71,76,81,86,91,96,700),ITR(33)                               
   71 LCT=7                                                                 
      NNP1=2                                                             
      H2P1=H2N                                                              
      N3P1=1                                                                
      H3P1=0.                                                            
      N3A3=-1                                                               
      H3A3S=H3N                                                          
      GO TO 200                                                             
!-----P1 -    axis-2 glide - no axis-3 glide, A3 - twofold - BG - Pba2   
   75 GO TO(76,81,86,91,96,700),ITR(34)                                  
   76 LCT=8                                                                 
      NNP1=2                                         
      H2P1=H2N                                                              
      N3P1=1                                                                
      H3P1=0.                                                            
      N3A3=1                                                                
      H3A3T=0.                                                           
      IF(MP1-2)200,250,200                                               
!-----P1 - no axis-2 glide -    axis-3 glide, A3 - screw - BH - Pca21    
   80 GO TO(81,86,91,96,700),ITR(35)                                     
   81 LCT=9                                                                 
      NNP1=3                                                             
      H2P1=0.                                                               
      N3P1=-1                                                               
      H3P1=H3N                                                           
      N3A3=-1                                                               
      H3A3S=H3N                                                          
      GO TO 200                                                             
!-----P1 - no axis-2 glide -    axis-3 glide, A3 - twofold - BI - Pnc2   
   85 GO TO(86,91,96,700),ITR(36)                                        
   86 LCT=10                                                                
      NNP1=3                                                             
      H2P1=0.                                                               
      N3P1=-1                                                               
      H3P1=H3N                                                           
      N3A3=1                                                                
      H3A3T=0.                                                           
      IF(MP1-3)200,250,200                                               
!-----P1 - no axis-2 glide - no axis-3 glide, A3 - screw - BJ - Pmn21    
   90 IF(ITR(37)-2)91,96,700                                             
   91 LCT=11                                                                
      NNP1=4                                                             
      H2P1=0.                                                               
      N3P1=1                                                                
      H3P1=0.                                                            
      N3A3=-1                                                               
      H3A3S=H3N
      GO TO 200                                                             
!-----P1 - no axis-2 glide - no axis-3 glide, A3 - twofold - BK - Pma2   
   95 IF(ITR(38).GT.1)GO TO 700                                          
   96 LCT=12                                                                
      NNP1=4                                                             
      H2P1=0.                                                               
      N3P1=1                                                                
      H3P1=0.                                                            
      N3A3=1                                                                
      H3A3T=0.                                                           
      IF(MP1-4)200,250,200                                               
!-----A2 placement                                                       
  100 MA2=NA2                                                            
      N=MARK2                                                               
      NSTP2=5                                                            
      DO 104 I=1,NMOD                                                       
      DO 103 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=H2A2+W(2,J)-W(2,I)                                                 
  101 D2=D2+H2M                                                             
      IF(D2.LE.CN(K))GO TO 101                                              
  102 D2=D2-H2M                                                             
      IF(D2.GT.CN(K))GO TO 102                                              
      IF(D2.LT.-CN(K))GO TO 103                                             
      IT(N)=K                                                               
      TI(N+1)=D2**2                                                         
      TI(N+2)=-W(1,J)-W(1,I)                                                
      TI(N+3)=-W(3,J)-W(3,I)                                                
      IT(N+4)=L                                                             
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 102                                                             
  103 CONTINUE                                                              
  104 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
      ER=ERMF                                                               
!-----Set up for negative axis-1 direction                               
!-----NEGATIVE X SIDE                                                       
      HL=-H1M                                                               
      X1A2=-100.                                                            
      HN=-H1Q                                                            
      Y1A2=-1.5*H1N                                                         
      JCT=0                                                              
      GO TO 110                                                             
!-----POSITIVE X SIDE                                                       
  109 CONTINUE                                                              
      HL=H1M                                                                
      X1A2=100.                                                             
      HN=H1Q                                                             
      Y1A2=1.5*H1N                                                          
  110 KE=NE                                                                 
      DY3A2=H3M/NV                                                          
      Y3A2=-H3N                                                             
      JE=1                                                                  
      MCT=0                                                                 
  116 N=MARK                                                                
      DO 128 I=MARK2,NEND2,NSTP2                                            
      K=IT(I)                                                               
      D3=Y3A2+TI(I+3)                                                       
  121 D3=D3+H3M                                                             
      IF(D3.LT.CN(K))GO TO 121                                              
  122 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 122                                              
      IF(D3.LT.-CN(K))GO TO 128                                             
      D=TI(I+1)+D3**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 122                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 122                                                             
  128 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 130                                                
      NEND=N-1                                                              
      CALL MINHI(Y1A2,HN,X1A2,HL,KE)                                     
      IF(KCT)130,151,152                                                    
  130 X1A2=HN                                                            
      X3A2=Y3A2                                                             
      GO TO 158                                                          
  151 X3A2=Y3A2                                                          
      X1A2=Y1A2                                                          
      IF(MCT)156,153,156                                                 
  152 IF(MCT)156,153,155                                                 
  153 IF(JE.GE.NV)GO TO 154                                              
      JE=JE+1                                                            
      Y3A2=Y3A2+DY3A2                                                    
      GO TO 116                                                          
  154 MCT=1                                                              
      JE=2.*JE                                                           
      DY3A2=H3M/JE                                                       
      Y3A2=X3A2+DY3A2                                                    
      GO TO 116                                                          
  155 MCT=-1                                                             
      Y3A2=X3A2-DY3A2                                                    
      GO TO 116                                                          
  156 IF(JE.LT.KE)GO TO 154                                              
      IF(JE.GE.NNE)GO TO 158                                             
      KE=NNE                                                             
      GO TO 154                                                          
  158 IF(HL.GE.0.)GO TO 160                                                 
      H1A2N=X1A2                                                            
      H3A2N=X3A2                                                            
      GO TO 109                                                             
  160 H1A2P=X1A2                                                         
      H3A2P=X3A2                                                         
      GO TO 255                                                             
!-----P1 placement                                                       
  200 MP1=NNP1                                                           
      MARK=MARK2                                                            
      N=MARK                                                                
      DO 205 I=MARK1,NEND1,NSTP1                                            
      K=N3P1*IT(I)                                                          
      IF(K.LE.0)GO TO 205                                                   
      D2=TI(I+3)+H2P1                                                       
  201 D2=D2+H2M                                                             
      IF(D2.LT.CN(K))GO TO 201                                              
  202 D2=D2-H2M                                                             
      IF(D2.GT.CN(K))GO TO 202                                              
      IF(D2.LT.-CN(K))GO TO 205                                             
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 202                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+5)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 202                                                             
  205 CONTINUE                                                              
      IF(N.GT.MARK)GO TO 206                                                
      H1P1P=H1N                                                             
      H1P1N=-H1N                                                         
      GO TO 255                                                             
  206 NEND=N-1                                                              
      KE=NNE                                                                
      H1P1P=1.5*H1N                                                         
      ER=ERMF                                                               
      CALL MINHI(H1P1P,H1N,100.,H1M,NNE)                                    
      H1P1N=-1.5*H1N                                                        
      CALL MINHI(H1P1N,-H1N,-100.,-H1M,NNE)                                 
  250 IF(N3A3)251,999,255                                                
  251 IF((MA3.EQ.1).OR.(MA3.EQ.3))GO TO 370     
      MA3=MA3+1                                                          
      GO TO 300                                                          
  255 IF(MA3.GT.1)GO TO 370                    
      MA3=2                                                              
!-----Set up for calculation of A3 placement on positive axis-1 side     
  300 HL=H1M                                                                
      X1A3=100.                                                             
      HN=H1Q                                                             
      Y1A3=1.5*H1N                                                          
      GO TO 310                                                             
!-----Set up for calculation of A3 placement on negative axis-1 side     
  305 HL=-H1M                                                               
      X1A3=-100.                                                            
      HN=-H1Q                                                            
      Y1A3=-1.5*H1N                                                         
  310 KE=NE                                                                 
      MCT=0                                                                 
      JE=1                                                                  
      DY2A3=H2M/NV                                                          
      Y2A3=-H2N                                                             
      ER=ERMF                                                               
      MARK=MARK2                                                            
  316 N=MARK                                                                
      DO 329 I=MARK1,NEND1,NSTP1                                            
      K=N3A3*IT(I)                                                          
      IF(K.LE.0)GO TO 329                                                   
      D2=TI(I+6)+Y2A3                                                       
  321 D2=D2+H2M                                                             
      IF(D2.LE.CN(K))GO TO 321                                              
  322 D2=D2-H2M                                                             
      IF(D2.GT.CN(K))GO TO 322                                              
      IF(D2.LT.-CN(K))GO TO 329                                             
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 322                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+5)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 322                                                             
  329 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 330                                                
      NEND=N-1                                                              
      CALL MINHI(Y1A3,HN,X1A3,HL,KE)                                     
      IF(KCT)330,351,352                                                    
  330 X1A3=HN                                                            
      X2A3=Y2A3                                                             
      GO TO 358                                                          
  351 X2A3=Y2A3                                                          
      X1A3=Y1A3                                                          
      IF(MCT)356,353,356                                                 
  352 IF(MCT)356,353,355                                                 
  353 IF(JE.GE.NV)GO TO 354                                              
      JE=JE+1                                                            
      Y2A3=Y2A3+DY2A3                                                    
      GO TO 316                                                          
  354 MCT=1                                                              
      JE=2.*JE                                                           
      DY2A3=H2M/JE                                                       
      Y2A3=X2A3+DY2A3                                                    
      GO TO 316                                                          
  355 MCT=-1                                                             
      Y2A3=X2A3-DY2A3                                                    
      GO TO 316                                                          
  356 IF(JE.LT.KE)GO TO 354                                              
      IF(JE.GE.NNE)GO TO 358                                             
      KE=NNE                                                             
      GO TO 354                                                          
  358 CONTINUE                       
      IF(HL.LE.0.)GO TO 365                                                 
      IF(N3A3)360,999,361                                                   
  360 H1A3SP=X1A3                                                           
      H2A3SP=X2A3                                                           
      GO TO 305                                                             
  361 H1A3TP=X1A3                                                           
      H2A3TP=X2A3                                                           
      GO TO 305                                                             
  365 IF(N3A3)366,999,367                                                   
  366 H1A3SN=X1A3                                                           
      H2A3SN=X2A3                                                           
      GO TO 370                                                          
  367 H1A3TN=X1A3                                                           
      H2A3TN=X2A3                                                           
  370 IF(N3P1.NE.0)GO TO 450                                             
      IF(N3A3)371,999,410                                                
  371 X1A1=H1A2P-H1A3SN                                                  
      X1A3=H1A3SN                                                        
      VX=X1A1*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 390                                              
!-----Check distances to derived A1 location                             
      X2A1=H2A2+H2A3SN                                                   
      X3A1=H3A2P+H3A3
      Y1A1=X1A1                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 379 I=1,NMOD                                                    
      DO 378 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2A1-W(2,J)-W(2,I)                                              
      D3=X3A1-W(3,J)-W(3,I)                                              
  372 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 372                                           
  373 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 373                                           
      IF(D2.LT.-CN(K))GO TO 378                                          
  374 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 374                                           
  375 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 375                                           
      IF(D3.LT.-CN(K))GO TO 373                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 375                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 375                                                          
  378 CONTINUE                                                           
  379 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 381                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1A1,X1A1,100.,H1M,NNE)                                 
      IF(KCT)381,380,380                                                 
  380 X1A1=Y1A1                                                          
      X1A3=H1A2P-X1A1                                                    
      VX=X1A1*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 390                                              
  381 X1I=2.*X1A1                                                        
      X1A3=H1A2P-X1A1                                                    
      X1I=2.*X1A1                                                        
      IF(IPR.LT.3)GO TO 385                                              
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
  383 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3)                    
      WRITE(NQ,384)H1A2P,H2A2,H3A2P,H1A3SN,H2A3SN,H3A3S                  
  384 FORMAT(5X,'A2 AT',3F7.3,'  A3 AT',3F7.3)                           
  385 V3=VX                                                                 
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1A1=X1A1                                                          
      H1A2=H1A2P                                                         
      H3A2=H3A2P                                                         
      H1A3=X1A3                                                          
      H2A3=H2A3SN                                                        
      H3A3=H3A3S                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  390 X1A1=H1A3SP-H1A2N                                                  
      X1A2=H1A2N
      VX=X1A1*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 408                                              
!-----Check distances to derived A1 location                             
      X2A1=H2A2+H2A3SP                                                   
      X3A1=H3A2N+H3A3
      Y1A1=X1A1                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 399 I=1,NMOD                                                    
      DO 398 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2A1-W(2,J)-W(2,I)                                              
      D3=X3A1-W(3,J)-W(3,I)                                              
  392 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 392                                           
  393 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 393                                           
      IF(D2.LT.-CN(K))GO TO 398                                          
  394 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 394                                           
  395 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 395                                           
      IF(D3.LT.-CN(K))GO TO 393                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 395                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 395                                                          
  398 CONTINUE                                                           
  399 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 401                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1A1,X1A1,100.,H1M,NNE)                                 
      IF(KCT)401,400,400                                                 
  400 X1A1=Y1A1                                                          
      X1A2=H1A3SP-X1A1                                                   
      VX=X1A1*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 408                                              
  401 X1I=2.*X1A1                                                        
      IF(IPR.LT.3)GO TO 407                                              
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,384)H1A2N,H2A2,H3A2N,H1A3SP,H2A3SP,H3A3S                  
  407 V3=VX                                                                 
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1A1=X1A1                                                          
      H1A2=X1A2                                                          
      H3A2=H3A2N                                                         
      H1A3=H1A3SP                                                        
      H2A3=H2A3SP                                                        
      H3A3=H3A3S                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  408 IF(V3.GE.VHT)GO TO 39                                                 
      VHT=V3                        
      HT1(1)=H1I                                                            
      HT1(5)=H1A2                                                           
      HT1(6)=H2A3                                                           
      HT1(7)=H3A2+H3A3                                                      
      GO TO 39                                                           
  410 X1A1=H1A2P-H1A3TN                                                  
      VX=X1A1*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 429                                              
!-----Check distances to derived A1 location                             
      X2A1=H2A2+H2A3TN                                                   
      X3A1=H3A2P+H3A3                                                    
      Y1A1=X1A1                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 419 I=1,NMOD                                                    
      DO 418 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2A1-W(2,J)-W(2,I)                                              
      D3=X3A1-W(3,J)-W(3,I)                                              
  412 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 412                                           
  413 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 413                                           
      IF(D2.LT.-CN(K))GO TO 418                                          
  414 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 414                                           
  415 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 415                                           
      IF(D3.LT.-CN(K))GO TO 413                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 415                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 920                                            
      GO TO 415                                                          
  418 CONTINUE                                                           
  419 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 421                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1A1,X1A1,100.,H1M,NNE)                                 
      IF(KCT)421,420,906                                                 
  420 X1A1=Y1A1                                                          
  421 VX=X1A1*H2M*H3N                                                    
      X1I=2.*X1A1                                                        
      X1A3=H1A2P-X1A1                                                    
  425 IF(VX.GE.V3)GO TO 429                                              
      IF(IPR.LT.3)GO TO 426                                              
      CALL PAGE (4,4,0)                                                  
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,384)H1A2P,H2A2,H3A2P,H1A3TN,H2A3TN,H3A3T                  
  426 V3=VX                                                                 
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1A1=X1A1                                                          
      H1A2=H1A2P                                                         
      H3A2=H3A2P                                                         
      H1A3=H1A3TN                                                        
      H2A3=H2A3TN                                                        
      H3A3=H3A3T                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  429 X1A1=H1A3TP-H1A2N                                                  
      VX=X1A1*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 444                                              
!-----Check distances to derived A1 location                             
      X2A1=H2A2+H2A3TP                                                   
      X3A1=H3A2N+H3A3                                                    
      Y1A1=X1A1                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 439 I=1,NMOD                                                    
      DO 438 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2A1-W(2,J)-W(2,I)                                              
      D3=X3A1-W(3,J)-W(3,I)                                              
  432 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 432                                           
  433 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 433                                           
      IF(D2.LT.-CN(K))GO TO 438                                          
  434 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 434                                           
  435 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 435                                           
      IF(D3.LT.-CN(K))GO TO 433                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 435                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 920                                            
      GO TO 435                                                          
  438 CONTINUE                                                           
  439 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 441                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1A1,X1A1,100.,H1M,NNE)                                 
      IF(KCT)441,440,906                                                 
  440 X1A1=Y1A1                                                          
  441 VX=X1A1*H2M*H3N                                                    
      X1I=2.*X1A1                                                        
      X1A2=H1A3TP-X1A1                                                   
      IF(IPR.LT.3)GO TO 442                                              
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,384)H1A2N,H2A2,H3A2N,H1A3TP,H2A3TP,H3A3                   
  442 IF(VX.GE.V3)GO TO 444                                              
      V3=VX                                                                 
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1A1=X1A1                                                          
      H1A2=H1A2N                                                         
      H3A2=H3A2N                                                         
      H1A3=H1A3TP                                                        
      H2A3=H2A3TP                                                        
      H3A3=H3A3T                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  444 IF(V3.GE.VHT)GO TO 39                                                 
      VHT=V3                                                             
      HT1(1)=H1I                                                            
      HT1(5)=H1A2                                                           
      HT1(6)=H2A3                                                           
      HT1(7)=H3A2+H3A3                                                      
      GO TO 39                                                           
  450 IF(N3A3)451,999,500                                                
  451 X1P2=H1P1P-H1A3SN                                                  
      VX=X1P2*H2M*H3N                                                    
      X1I=2.*X1P2                                                        
      IF(IPR.LT.3)GO TO 455                                              
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,452)H1P1P,H2P1,H3P1,H1A3SN,H2A3SN,H3A3S                   
  452 FORMAT(1X,'P1 AT',3F7.3,'  A3 AT',3F7.3)                           
  455 IF(VX.GE.V3)GO TO 475                                              
!-----Check distances to derived P2 location                             
      X2P2=H2P1+H2A3SN                                                   
      X3P2=H3P1+H3A3                                                     
      Y1P2=X1P2                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 469 I=1,NMOD                                                    
      DO 468 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2P2-W(2,J)-W(2,I)                                              
      D3=X3P2+W(3,J)-W(3,I)                                              
  462 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 462                                           
  463 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 463                                           
      IF(D2.LT.-CN(K))GO TO 468                                          
  464 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 464                                           
  465 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 465                                           
      IF(D3.LT.-CN(K))GO TO 463                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 465                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 920                                            
      GO TO 465                                                          
  468 CONTINUE                                                           
  469 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 471                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1P2,X1P2,100.,H1M,NNE)                                 
      IF(KCT)471,470,906                                                 
  470 X1P2=Y1P2                                                          
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 475                                              
  471 V3=VX                                                                 
      X1I=2.*X1P2
      MP2=1
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1P2=X1P2                                                          
      H1P1=H1P1P                                                         
      H1A3=H1A3SN                                                        
      H2A3=H2A3SN                                                        
      H3A3=H3A3S                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  475 X1P2=H1A3SP-H1P1N                                                  
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 495                                              
!-----Check distances to derived P2 location                             
      X2P2=H2P1+H2A3SP                                                   
      X3P2=H3P1+H3A3                                                     
      Y1P2=X1P2                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 489 I=1,NMOD                                                    
      DO 488 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2P2-W(2,J)-W(2,I)                                              
      D3=X3P2+W(3,J)-W(3,I)                                              
  482 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 482                                           
  483 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 483                                           
      IF(D2.LT.-CN(K))GO TO 488                                          
  484 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 484                                           
  485 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 485                                           
      IF(D3.LT.-CN(K))GO TO 483                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 485                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 920                                            
      GO TO 485                                                          
  488 CONTINUE                                                           
  489 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 491                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1P2,X1P2,100.,H1M,NNE)                                 
      IF(KCT)491,490,906                                                 
  490 X1P2=Y1P2                                                          
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 495                                              
  491 V3=VX                                                                 
      X1I=2.*X1P2                                                        
      IF(IPR.LT.3)GO TO 492                                              
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,452)H1P1N,H2P1,H3P1,H1A3SP,H2A3SP,H3A3S                   
  492 MP2=2
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1P2=X1P2                                                          
      H1P1=H1P1N                                                         
      H1A3=H1A3SP                                                        
      H2A3=H2A3SP                                                        
      H3A3=H3A3S                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  495 IF(V3.GE.VHT)GO TO 39                                                 
      VHT=V3                                                             
      HT1(1)=H1I                                                            
      HT1(5)=H1A3                                                           
      HT1(6)=H2A3                                                           
      HT1(7)=0.                                                             
      GO TO 39                                                           
  500 X1P2=H1P1P-H1A3TN                                                  
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 515
!-----Check distances to derived P2 location                             
      X2P2=H2P1+H2A3TN                                                   
      X3P2=H3P1+H3A3                                                     
      Y1P2=X1P2                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 509 I=1,NMOD                                                    
      DO 508 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2P2-W(2,J)-W(2,I)                                              
      D3=X3P2+W(3,J)-W(3,I)                                              
  502 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 502                                           
  503 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 503                                           
      IF(D2.LT.-CN(K))GO TO 508                                          
  504 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 504                                           
  505 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 505                                           
      IF(D3.LT.-CN(K))GO TO 503                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 505                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 920                                            
      GO TO 505                                                          
  508 CONTINUE                                                           
  509 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 511                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1P2,X1P2,100.,H1M,NNE)                                 
      IF(KCT)511,510,906                                                 
  510 X1P2=Y1P2                                                          
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 515                                              
  511 V3=VX                                                                 
      X1I=2.*X1P2                                                        
      IF(IPR.LT.3)GO TO 512                                              
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,452)H1P1P,H2P1,H3P1,H1A3TN,H2A3TN,H3A3T                   
  512 MP2=3
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1P2=X1P2                                                          
      H1P1=H1P1P                                                         
      H1A3=H1A3TN                                                        
      H2A3=H2A3TN                                                        
      H3A3=H3A3T                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  515 X1P2=H1A3TP-H1P1N                                                  
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 535
!-----Check distances to derived P2 location                             
      X2P2=H2P1+H2A3TP                                                   
      X3P2=H3P1+H3A3                                                     
      Y1P2=X1P2                                                          
      N=MARK                                                             
      NARK=0                                                             
      DO 529 I=1,NMOD                                                    
      DO 528 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2P2-W(2,J)-W(2,I)                                              
      D3=X3P2+W(3,J)-W(3,I)                                              
  522 D2=D2+H2M                                                          
      IF(D2.LE.CN(K))GO TO 522                                           
  523 D2=D2-H2M                                                          
      IF(D2.GT.CN(K))GO TO 523                                           
      IF(D2.LT.-CN(K))GO TO 528                                          
  524 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 524                                           
  525 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 525                                           
      IF(D3.LT.-CN(K))GO TO 523                                          
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 525                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 920                                            
      GO TO 525                                                          
  528 CONTINUE                                                           
  529 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 531                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y1P2,X1P2,100.,H1M,NNE)                                 
      IF(KCT)531,530,906                                                 
  530 X1P2=Y1P2                                   
      VX=X1P2*H2M*H3N                                                    
      IF(VX.GE.V3)GO TO 535                                              
  531 V3=VX                                                                 
      X1I=2.*X1P2                                                        
      IF(IPR.LT.3)GO TO 532                                              
      WRITE(NQ,383)ALPHA(LCT),A1,A2,A3,VX,X1I,H2M,H3M                    
      WRITE(NQ,452)H1P1N,H2P1,H3P1,H1A3TP,H2A3TP,H3A3T                   
  532 MP2=4
      LCTM=LCT                                                           
      H1I=X1I                                                            
      H1P2=X1P2                                                          
      H1P1=H1P1N                                                         
      H1A3=H1A3TP                                                        
      H2A3=H2A3TP                                                        
      H3A3=H3A3T                                                         
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
  535 IF(V3.GE.VHT)GO TO 39                                                 
      VHT=V3                                                             
      HT1(1)=H1I                                                            
      HT1(5)=H1A3                                                           
      HT1(6)=H2A3                                                           
      HT1(7)=0.                                                             
      GO TO 39                                                           
  700 IF(VHT.GE.HT1(8))GO TO 800                                         
      HT1(2)=0.                                                          
      HT1(3)=H2M                                                            
      HT1(4)=H3M                                                            
      HT1(8)=VHT                                                            
      CHT1(9)=ALPHA(LCTM)                                                   
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=14                                                               
      MLOW=LCTM+26                                                       
  800 IF(IPR.LT.2)GO TO 1000                                             
      IF(V3.GE.V3S)GO TO 1000                                            
      WRITE(NQ,801)ALPHA(LCTM),A1,A2,A3,V3,H1I,H2M,H3M                   
  801 FORMAT(1X,A2,3F7.1,'  V =',F8.2,'  AXES',3F8.3)                    
      DO 819 I=1,3                                                       
      IF(LCTM.LE.4)GO TO(806,808,809),I                                  
      IF(I-2)805,807,809                                                 
  805 SM(1)=H1P1                                                         
      SM(2)=H2P1                                                         
      SM(3)=H3P1                                                         
      SM(4)=SM(1)-HT1(5)                                                 
      SM(5)=SM(2)                                                        
      SM(6)=SM(3)                                                        
      N=3                                                                
      GO TO 812                                                          
  806 SM(1)=H1A1                                                         
      SM(2)=H2A2+H2A3                                                    
      SM(3)=H3A2+H3A3                                                    
      SM(4)=SM(1)                                                        
      SM(5)=SM(2)-HT1(6)                                                 
      SM(6)=SM(3)-HT1(7)                                                 
      N=1                                                                
      GO TO 812                                                          
  807 SM(1)=H1P2                                                         
      SM(2)=H2P1+H2A3                                                    
      SM(3)=H3P1+H3A3                                                    
      SM(4)=SM(1)                                                        
      SM(5)=SM(2)-HT1(6)                                                 
      SM(6)=SM(3)                                                        
      N=4                                                                
      GO TO 812                                                          
  808 SM(1)=H1A2                                                         
      SM(2)=H2A2                                                         
      SM(3)=H3A2                                                         
      SM(4)=SM(1)-HT1(5)                                                 
      SM(5)=SM(2)                                                        
      SM(6)=SM(3)-HT1(7)                                                 
      N=2                                                                
      GO TO 812                                                          
  809 SM(1)=H1A3                                                         
      SM(2)=H2A3                                                         
      SM(3)=H3A3                                                         
      SM(4)=SM(1)-HT1(5)                                                 
      SM(5)=SM(2)-HT1(6)                                                 
      SM(6)=SM(3)                                                        
      N=5                                                                
  812 SM(4)=SM(4)/H1I                                                    
      SM(5)=SM(5)/H2M                                                    
      SM(6)=SM(6)/H3M                                                    
      IF(SM(4).GT..5)SM(4)=SM(4)-1.                                      
      IF(SM(4).LE.-.5)SM(4)=SM(4)+1.                                     
      IF(SM(5).GT..5)SM(5)=SM(5)-1.                                      
      IF(SM(5).LE.-.5)SM(5)=SM(5)+1.                                     
      IF(SM(6).GT..5)SM(6)=SM(6)-1.                                      
      IF(SM(6).LE.-.5)SM(6)=SM(6)+1.                                     
      D=SQRT(SM(1)**2+SM(2)**2+SM(3)**2)                                 
      WRITE(NQ,815)SYM(N),SM,D                                           
  815 FORMAT(1X,A2,' AT',3X,3F8.3,';  ',3F8.4,'  D =',F6.2)              
  819 CONTINUE                                                           
      SM(1)=HT1(5)/2.                                                    
      SM(2)=HT1(6)/2.                                                    
      SM(3)=HT1(7)/2.                                                    
      SM(4)=SM(1)/H1I                                                    
      SM(5)=SM(2)/H2M                                                    
      SM(6)=SM(3)/H3M                                                    
      WRITE(NQ,820)SM                                                    
  820 FORMAT(1X,'ORIGIN',2X,3F8.3,';  ',3F8.4)                        
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,NARK,NSTP,LIMIT        
  903 FORMAT(26H0STORAGE EXCEEDED BY ZORTC,5X,4(I9,I3),I9)                  
      GO TO 999                                                          
  906 CALL PAGE (2,2,0)                                                  
      WRITE(NQ,907)                                                      
  907 FORMAT(' PARAMETER TOO LARGE')                                     
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)KCT,MCT,N,MARK3,NARK,MARK,NEND2                          
  921 FORMAT(16H0NO INTERACTIONS,2I3,5I5)                                   
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 IF(IPR.GE.2)PRINT 1001,MP2
 1001 FORMAT(' MP2 =',I2)
      RETURN                                                                
      END SUBROUTINE ZORTC 
