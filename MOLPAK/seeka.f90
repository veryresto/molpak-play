      SUBROUTINE SEEKA
!
      USE molpakCommonMod
      
      IMPLICIT NONE

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
!      DIMENSION IJD(500),CHT1(14),CDIJ(500)                              
!      EQUIVALENCE (DIJ(1),IJD(1)),(HT1(1),CHT1(1)),(DIJ(1),CDIJ(1))      
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
!     COMMON TI(70000)                                                       
!     DIMENSION IT(70000)                                                
!      COMMON TI(100000)                                                    
!      COMMON /NEW/IDCHG,IDATOM(200),G92CHARGE(200),ESPCHARGE(200)        
!     DIMENSION IT(100000)                                               
!     EQUIVALENCE (TI(1),IT(1))                                          
!
      CHARACTER(4)  :: BLANK= '    '
      CHARACTER(5) :: ALPHA(6) = (/'A    ','B    ','C    ','ALPHA',&
     &                             'BETA ','GAMMA'/)
     
      INTEGER :: I, IRG, J, K, KK, L, MCT, MINV 
      INTEGER :: N, N3, NA3, NAP, NAVE, NH3, NMV
!
      REAL :: A, AM3, ANG3(200), AS3, AVE 
      REAL :: D, DA3,  D1, D2, D3, E
      REAL :: H1, H2, H3, VMIN
!    
!-----Does the molecule have an internal mirror plane or two-fold axis?  
!       If not, quit                                                      
      IF((NCTR.LE.0).OR.(NCTR.GT.6))GO TO 922                            
      MCT=0                                                                 
      NLOW=0                                                                
      MLOW=0                                                                
      ERMT=2.*ERM                                                           
      ERMF=4.*ERM                                                           
      AVE=0.                                                                
      NAVE=0                                                                
      HT1(8)=1000000.      ! 03-15-06                                               
      WRITE (NQ,10) BUF                                                     
   10 FORMAT (1X,A80)                                                    
!-----Does the SEEK line specify a full non-redundant angular search at  
!      the angular intervals specified in columns 7-9 ?                  
      READ(BUF,11)IRG,DA3                                                
   11 FORMAT(4X,I2,F6.1)                                                  
      IF(IRG.NE.0)GO TO 15                                                  
!     If the angular interval is not specified, set it to 10 degrees.    
      IF(DA3.LT..1)DA3=10.                                                  
      AS3=-90.                                                              
      AM3=90.                                                               
      GO TO 24                                                              
!-----The SEEK line specifies the an angular range to be searched        
!     READ COMPLETE RANGE AND INTERVAL INSTRUCTIONS                      
   15 READ(BUF,16)IRG,AS3,AM3,DA3                                        
!  16 FORMAT(4X,I2,24X,3F4.0)                                            
   16 FORMAT(4X,I2,36X,3F6.1)                                            
!-----PRINT RANGE AND INTERVAL INSTRUCTIONS                              
      WRITE(NQ,17)AS3,AM3,DA3                                            
   17 FORMAT (4X,' SEKA BEGINNING : ', 2X, 3F6.1 )                       
!-----SET UP ORIENTATION SCAN OR SEARCH                                     
      IF(DA3.LT..1)DA3=10.                                                  
!-----CALCULATE SIN AND COS TABLES                                          
   24 CONTINUE                                                              
      DIJ(500)=DA3                                                       
   30 E=DA3/57.296                                                          
      D=AS3/57.296                                                          
      A3=AS3                                                                
      NA3=1                                                                 
      NAP=1                                                                 
   31 CS3(NA3)=COS(D)                                                       
      SN3(NA3)=SIN(D)                                                       
      ANG3(NAP)=A3                                                          
      IF(A3.GE.AM3)GO TO 35                                                 
      IF(NA3.GE.200)GO TO 900                                               
      D=D+E                                                                 
      A3=A3+DA3                                                             
      NA3=NA3+1                                                             
      NAP=NAP+1                                                             
      GO TO 31                                                              
   35 CONTINUE                                                              
!-----Orient molecule so that mirror plane is perpendicular to axis-3 or 
!     two-fold axis is parallel to axis-3                                
      GO TO(40,45,50,40,45,50),NCTR                                      
   40 A1=0.                                                              
      A2=90.                                                             
      DO 41 I=1,NMOD                                                     
      V(1,I)=X(3,I)                                                      
      V(2,I)=X(2,I)                                                      
      W(3,I)=-X(1,I)                                                     
   41 CONTINUE                                                           
      GO TO 52                                                           
   45 A1=90.                                                             
      A2=90.                                                             
      DO 46 I=1,NMOD                                                     
      V(1,I)=X(3,I)                                                      
      V(2,I)=X(1,I)                                                      
      W(3,I)=X(2,I)                                                      
   46 CONTINUE                                                           
      GO TO 52                                                           
   50 A1=0.                                                              
      A2=0.                                                              
      DO 51 I=1,NMOD                                                     
      V(1,I)=X(1,I)                                                      
      V(2,I)=X(2,I)                                                      
      W(3,I)=X(3,I)                                                      
   51 CONTINUE                                                           
!-----Calculate minimum axis-3 length                                    
   52 N=MARK1                                                               
      MARK=N                                                                
      NSTP=3                                                             
      NARK=0                                                             
      DO 56 I=1,NMOD                                                        
      DO 55 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=V(1,J)-V(1,I)                                                      
      IF(ABS(D1).GT.CN(K))GO TO 55                                          
      D2=V(2,J)-V(2,I)                                                      
      IF(ABS(D2).GT.CN(K))GO TO 55                                          
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 55                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(3,J)-W(3,I)                                              
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   55 CONTINUE                                                              
   56 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      H3M=7.                                                                
      ER=ERM                                                                
      CALL MINHI(H3M,0.,100.,7.,NNE)                                        
      IF(KCT)904,57,906                                                  
   57 H3N=.5*H3M                                                            
      IF(IPR.LT.3)GO TO 69                                                  
      WRITE(NQ,58)H3M                                                    
   58 FORMAT(1X,'H3M =',F7.3)                                            
   69 N3=1                                                                  
      A3=AS3                                                                
      NP2=1                                                                 
   70 DIJ(NP2)=.001                                                         
      CDIJ(NP2+230)=BLANK                                                   
      V3=10000.                                                             
      DO 71 I=1,NMOD                                                        
      W(1,I)=V(1,I)*CS3(N3)-V(2,I)*SN3(N3)                               
      W(2,I)=V(1,I)*SN3(N3)+V(2,I)*CS3(N3)                               
   71 CONTINUE                                                           
      IF(IPR.LT.4)GO TO 74                                               
      WRITE(NQ,72)                                                       
   72 FORMAT(1H0,'POSITIONS AFTER ROTATION')                             
      DO 73 I=1,NMOD                                                     
      WRITE(NQ,47)AT(I),AN(I),W(1,I),W(2,I),W(3,I)                       
   47 FORMAT(1X,A2,A4,3F8.3)                                             
   73 CONTINUE                                                           
!-----Separate various close axis-3 distances                            
   74 N=MARK1                                                               
      NSTP1=8                                                            
      NARK=0                                                             
      DO 79 I=1,NMOD                                                        
      DO 78 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      KK=K                                                                  
      D3=W(3,J)-W(3,I)                                                      
!-----NH3 is the number of times H3M has been added to D3                
      NH3=0                                                              
   75 D3=D3+H3M                                                             
      NH3=NH3+1                                                          
      IF(D3.LE.CN(K))GO TO 75                                               
   76 D3=D3-H3M                                                             
      NH3=NH3-1                                                          
      IF(D3.GT.CN(K))GO TO 76                                               
      IF(D3.LT.-CN(K))GO TO 77                                              
      IT(N)=KK                                                              
      TI(N+1)=D3**2                                                         
      TI(N+2)=W(1,J)-W(1,I)                                                 
      TI(N+3)=W(2,J)-W(2,I)                                                 
      IT(N+4)=L                                                             
      TI(N+5)=-W(1,J)-W(1,I)                                                
      TI(N+6)=-W(2,J)-W(2,I)                                                
      IT(N+7)=NH3                                                        
      N=N+NSTP1                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 76                                                              
!-----Have values offset by H3N been calculated?                         
   77 IF(KK.LE.0)GO TO 78                                                   
!-----A negative value signals a distance to an offset atom              
      KK=-K                                                                 
      D3=D3+H3N                                                             
      GO TO 75                                                              
   78 CONTINUE                                                              
   79 CONTINUE                                                              
      IF(N.LE.MARK1)GO TO 920                                               
      NEND1=N-1                                                             
      MARK2=N                                                            
!-----Calculate separation along axis-1 of line of I molecules along     
!     axis-3, H1M - Establish plane-1,3 I grid                           
      MARK=MARK2                                                            
      DO 80 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
!-----Skip distances to atoms offset by H3N                              
      IF(K.LE.0)GO TO 80                                                    
      IF(ABS(TI(I+3)).GT.CN(K))GO TO 80                                     
      D=TI(I+1)+TI(I+3)**2                                                  
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 80                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   80 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      ER=ERM                                                             
      H1M=7.                                                             
      CALL MINHI(H1M,0.,100.,7.,NNE)                                        
      IF(KCT)904,81,906                                                  
   81 H1N=H1M/2.                                                         
!-----Calculate separation along axis-2 of a line of I molecules along   
!     axis-3, H2M                                                        
      N=MARK2                                                            
      DO 82 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
!-----Skip distances to atoms offset by H3N                              
      IF(K.LE.0)GO TO 82                                                    
      IF(ABS(TI(I+2)).GT.CN(K))GO TO 82                                     
      D=TI(I+1)+TI(I+2)**2                                                  
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 82                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   82 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      ER=ERM                                                             
      H2M=7.                                                             
      CALL MINHI(H2M,0.,100.,7.,NNE)                                        
      IF(KCT)904,83,906                                                  
   83 H2N=H2M/2.                                                         
      IF(IPR.LT.3)GO TO 90                                               
      WRITE(NQ,85)H1M,H2M                                                
   85 FORMAT(1H ,'H1M =',F8.3,', H2M =',F8.3)                            
!-----Transfer to specified coordination sphere finding subroutines      
   90 GO TO(91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107  & 
     &,108,109,110,111,112,113,114,115,116,117,999),JTR(1)               
   91 CALL ZTRTA                                                            
      GO TO(92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108 & 
     &,109,110,111,112,113,114,115,116,117,198),JTR(2)                   
   92 CALL ZTRTB                                                            
      GO TO(93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109& 
     &,110,111,112,113,114,115,116,117,198),JTR(3)                       
   93 CALL ZTRTC                                                            
      GO TO(94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109   & 
     &,110,111,112,113,114,115,116,117,198),JTR(4)                       
   94 CALL ZFRTA                                                            
      GO TO(95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110  & 
     &,111,112,113,114,115,116,117,198),JTR(5)                           
   95 CALL ZFRTB                                                            
      GO TO(96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111 & 
     &,112,113,114,115,116,117,198),JTR(6)                                 
   96 CALL ZFRTC                                                            
      GO TO(97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112& 
     &,113,114,115,116,117,198),JTR(7)                                   
   97 CALL ZORTA                                                            
      GO TO(98,99,100,101,102,103,104,105,106,107,108,109,110,111,112   & 
     &,113,114,115,116,117,198),JTR(8)                                   
   98 CALL ZORTB                                                            
      GO TO(99,100,101,102,103,104,105,106,107,108,109,110,111          & 
     &,112,113,114,115,116,117,198),JTR(9)                               
   99 CALL ZORTC                                                            
      GO TO(100,101,102,103,104,105,106,107,108,109,110,111             & 
     &,112,113,114,115,116,117,198),JTR(10)                              
  100 CALL ZERTA                                                         
      GO TO(101,102,103,104,105,106,107,108,109,110,111                 & 
     &,112,113,114,115,116,117,198),JTR(11)                              
  101 CALL ZERTB                                                         
      GO TO(102,103,104,105,106,107,108,109,110,111,112,113,114,115,116 & 
     &,117,198),JTR(12)                                                  
  102 CALL ZNRTA                                                         
      GO TO(103,104,105,106,107,108,109,110,111,112,113,114,115,116     & 
     &,117,198),JTR(13)                                                  
  103 CALL ZNRTC                                                         
      GO TO(104,105,106,107,108,109,110,111,112,113,114,115,116         &  
     &,117,198),JTR(14)                                                  
  104 CALL ZNRTB                                                         
      GO TO(105,106,107,108,109,110,111,112,113,114,115,116             & 
     &,117,198),JTR(15)                                                  
  105 CALL ZSRTA                                                         
      GO TO(106,107,108,109,110,111,112,113,114,115,116,117,198),JTR(16) 
  106 CALL ZSRTB                                                         
      GO TO(107,108,109,110,111,112,113,114,115,116,117,198),JTR(17)     
  107 CALL ZSRTC                                                         
      GO TO(108,109,110,111,112,113,114,115,116,117,198),JTR(18)         
  108 CALL ZSRTD                                                         
      GO TO(109,110,111,112,113,114,115,116,117,198),JTR(19)             
  109 CALL ZSRTE                                                         
      GO TO(110,111,112,113,114,115,116,117,198),JTR(20)                 
  110 CALL ZSRTF                                                         
      GO TO(111,112,113,114,115,116,117,198),JTR(21)                     
  111 CALL ZSRTG                                                         
      GO TO(112,113,114,115,116,117,198),JTR(22)                         
  112 CALL ZFRTD                                                         
      GO TO(113,114,115,116,117,198),JTR(23)                             
  113 CALL ZFRTE                                                         
      GO TO(114,115,116,117,198),JTR(24)                                 
  114 CALL ZSRTH                                                         
      GO TO(115,116,117,198),JTR(25)                                     
  115 CALL ZSRTI                                                         
      IF(JTR(26)-2)116,117,198                                           
  116 CALL ZSRTJ                                                         
      IF(JTR(27).GT.1)GO TO 198                                          
  117	CALL ZSRTK                                                         
  198 CDIJ(NP2+230)=CDIJ(NP2+20)                                          
      IF(A3.GE.AM3)GO TO 213                                                
      N3=N3+1                                                               
      NP2=NP2+1                                                             
      A3=A3+DA3                                                             
      GO TO 70                                                              
  213 WRITE(NQ,214)                                                         
  214 FORMAT(1H0,'UNIT CELL VOLUME PER MOL. VS. MOLECULE ROTATION ABOUT &  
     &AXIS-3')                                                           
      WRITE(NQ,215)(ANG3(I),I=1,NAP)                                        
  215 FORMAT(1H0,6H  A3  ,6H MIN V,200F6.1)                              
      IF(MCOUNT.LE.0)GO TO 230                                           
!-----Add cell volumes/molecule to MCOUNT collection of smallest ones    
!     PRINT '(4HNP2=,1X,I4)',NP2                                         
!     DO 229 I=201,NP2                                                   
      DO 229 I = 1, NP2                                                   
      IF((ICOUNT.GE.MCOUNT).OR.(I.EQ.NP2))GO TO 221                  
      ICOUNT =ICOUNT+1                                                   
      IORDER(ICOUNT)=ICOUNT                                              
      J=ICOUNT                                                           
      GO TO 222                                                          
  221 J=IORDER(ICOUNT)                                                   
      IF(DIJ(I).GE.VOLUME(J))GO TO 229                                   
  222 ANGLE1(J) = A1                                                     
      ANGLE2(J) = A2                                                     
      ANGLE3(J) = ANG3(I)                                                
      VOLUME(J) = DIJ(I)                                                 
      CODE(J)=CDIJ(I+230)                                                
!     PRINT *, ANGLE1(J),ANGLE2(J),ANGLE2(J),VOLUME(J),CODE(J)           
      K=ICOUNT                                                           
  223 K=K-1                                                              
      IF(K.LE.0)GO TO 224                                                
      L=IORDER(K)                                                        
      IF(VOLUME(L).LE.DIJ(I))GO TO 224                                   
      IORDER(K+1)=L                                                      
      GO TO 223                                                          
  224 K=K+1                                                              
      IORDER(K)=J                                                        
  229 CONTINUE                                                           
!-----Contribute to sums used to calculate average value of volume,      
!      determine minimum volume and fix                                  
      VMIN=1000000.      ! 03-15-06                                                     
  230 DO 231 I=1,NP2                                                        
      NAVE=NAVE+1                                                           
      AVE=AVE+DIJ(I)                                                        
      IF(DIJ(I).GT.VMIN)GO TO 231                                        
      VMIN=DIJ(I)                                                        
      MINV=DIJ(I)+.5                                                     
      NMV=I                                                              
  231 IJD(I)=DIJ(I)+.5                                                      
!-----Print line of volumes/molecule                                     
      WRITE(NQ,235)ANG3(NMV),MINV,DIJ(NMV+230),(IJD(I),CDIJ(I+230),  &       
     &I=1,NP2)                                                           
  235 FORMAT(1H ,F6.1,I4,A2,200(I4,A2))                                  
!-----Calculate and print average volume and ratio to smallest volume    
      IF(NAVE.EQ.0)NAVE=1                                                
      AVE=AVE/NAVE                                                          
      IF(HT1(8).EQ.0.)HT1(8)=1.                                          
      D=AVE/HT1(8)                                                          
      WRITE(NQ,238)AVE,D                                                    
  238 FORMAT(10H0AVERAGE =,F9.3,12H, AVE/BEST =,F6.3)                       
!-----Print out parameters of smallest volume cell found by this  
!      instruction                                                      
  370 IF(HT1(8).GT.10000.)GO TO 1000                                        
      WRITE(NQ,371)(HT1(I),I=10,12),NLOW,CHT1(9),HT1(8),ERM                 
  371 FORMAT(24H0SMALLEST CELL AT   A1 =,F6.1,7H,  A2 =,F6.1,7H,  A3 =,&     
     &F6.1,7H   CASE,I3,1X,A,8H,  V/M =,F7.2,6H,  E =,F6.3/1H )             
      E=90.                                                                 
      D=0.                                                                  
      WRITE(21,375)(HT1(I),I=10,12),HT1(8),CHT1(9)                       
  375 FORMAT(/4X,'MOLPAK ANGLES & VOLUME = ',3F6.1,F7.2,2X,A2)            
      GO TO(380,380,385,385,385, 385,380,390,390,390, 390,395,395,395, &  
     & 395, 390,390,390,390,390, 395,385,385,385,395, 395,395,395,395, &   
     & 390, 390,395,395,395,395, 395),NLOW                                  
!-----ZTRT result for triclinic structures                               
  380 D=SQRT(HT1(2)**2+HT1(3)**2)                                        
      IF(HT1(14).LE.HT1(3)-HT1(14))GO TO 381                             
      HT1(14)=HT1(14)-HT1(3)                                             
      HT1(13)=HT1(13)-HT1(2)                                             
  381 IF(HT1(13).GE.HT1(1)-HT1(13))HT1(13)=HT1(13)-HT1(1)                
      D3=SQRT(HT1(4)**2+HT1(13)**2+HT1(14)**2)                           
      A1=57.96*ACOS((HT1(2)*HT1(13)+HT1(3)*HT1(14))/(D*D3))              
      A=57.96*ACOS(HT1(13)/D3)                                           
      A3=57.96*ACOS(HT1(2)/D)                                            
      WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),A1,HT1(1),D,D                
      WRITE(NQ,386)ALPHA(2),D,ALPHA(5),A,HT1(2),HT1(3),D                 
      WRITE(NQ,386)ALPHA(3),D3,ALPHA(6),A3,HT1(13),HT1(14),HT1(4)        
      IF(NLOW.EQ.1)GO TO 1000               
      GO TO 396                                                          
!-----ZTRT or ZSRT result for Z=2 monoclinic structures                  
  385 D=SQRT(HT1(2)**2+HT1(3)**2)                                        
      A3=57.96*ACOS(HT1(2)/D)                                            
      WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),E,HT1(1),D,D                 
      WRITE(NQ,386)ALPHA(2),D,ALPHA(5),E,HT1(2),HT1(3),D                 
      WRITE(NQ,386)ALPHA(3),HT1(4),ALPHA(6),A3,D,D,HT1(4)                
  386 FORMAT(1H ,3X,A1,2H =,F7.3,3X,A5,2H =,F7.1,8X,3F8.3)               
      IF(NLOW.EQ.22)GO TO 1000                                           
      GO TO 396                                                          
!-----ZFRT or ZNRT result for Z=4 monoclinic structures                  
  390 D=SQRT(HT1(2)**2+HT1(3)**2)                                        
      A3=57.96*ACOS(HT1(2)/D)                                            
      WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),E,HT1(1),D,D                 
      WRITE(NQ,386)ALPHA(2),D,ALPHA(5),E,HT1(2),HT1(3),D                 
      WRITE(NQ,386)ALPHA(3),HT1(4),ALPHA(6),A3,D,D,HT1(4)                
      GO TO 396                                                             
!-----ZORT, ZERT, ZNRT, or ZSRT result for orthorhombic structures        
  395 WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),E,HT1(1),D,D                 
      WRITE(NQ,386)ALPHA(2),HT1(3),ALPHA(5),E,D,HT1(3),D                 
      WRITE(NQ,386)ALPHA(3),HT1(4),ALPHA(6),E,D,D,HT1(4)                 
  396 H1=HT1(5)/2.                                                       
      H2=HT1(6)/2.                                                       
      H3=HT1(7)/2.                                                       
      WRITE(NQ,397)H1,H2,H3                                              
  397 FORMAT(1H0,30X,8HORIGIN  ,3F8.3)                                      
      GO TO 1000                                                         
  900 WRITE(NQ,901)                                                         
      WRITE(17,901)                                                         
  901 FORMAT(16H0TOO MANY ANGLES)                                           
      GO TO 999                                                             
  902 WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT               
      WRITE(17,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT               
  903 FORMAT('0STORAGE EXCEEDED BY MOLPAK',5X,4(I9,I3),I9)               
      GO TO 999                                                             
  904 WRITE(NQ,905)                                                      
  905 FORMAT(20H0PARAMETER TOO SMALL)                                    
      GO TO 999                                                          
  906 WRITE(NQ,907)                                                      
      WRITE(17,907)                                                      
  907 FORMAT(20H0PARAMETER TOO LARGE)                                    
      GO TO 999                                                          
  920 WRITE(NQ,921)KCT,MCT,N,MARK1,MARK,MARK,NARK                           
      WRITE(17,921)KCT,MCT,N,MARK1,MARK,MARK,NARK                           
  921 FORMAT(16H0NO INTERACTIONS,I3,5I9)                                    
      GO TO 999                                                             
  922 WRITE(NQ,923)                                                      
      WRITE(17,923)                                                      
  923 FORMAT('0SEEKA CALLED INSTEAD OF SEEK')                            
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE SEEKA

