      SUBROUTINE SEEK                                                    
                                                                    
      USE molpakCommonMod
      
      IMPLICIT NONE
      
!      CHARACTER*80 BUF                                                   CCC    
!     CHARACTER*60 HEAD                                                  CCC
!      CHARACTER*4 AN,CDIJ,CHT1                                           CCC
!      CHARACTER*2 AT,CODE                                                CCC 
!----General program communication and transfer variables               CCC
!     COMMON BUF,KILL,NPAGE,LINE,NEW,JTR(30),ITR(90)                     JH-95
!-----Energy constants                                                   CCC    
!     COMMON NSEP,CN(1200),ERM,FERM,ER,ERMT,ERMF,NV,NE,NNE               CCC    
!-----Compound input data                                                CCC
!     COMMON NMOD,NCTR,IA(200),IAA(200),X(3,200),AT(200),AN(200)         CCC  
!-----Compound rotation and structure parameters                         CCC
!     COMMON CS2(37),SN2(37),CS3(19),SN3(19),V(3,100),W(3,100)           CCC
!     COMMON NLOW,MLOW,NP2,KCT                                           CCC
!     COMMON H1M,H2M,H3M,H1N,H2N,H3N,A1,A2,A3,V3                         CCC
!-----Output flags and storage                                           CCC
!     COMMON NQ,NR,NS,NWM,NSGT,IPR,HEAD,DIJ(500),HT1(14)                 CCC
!     DIMENSION IJD(500),CHT1(14),CDIJ(500)                              CCC
!     EQUIVALENCE (DIJ(1),IJD(1)),(HT1(1),CHT1(1)),(DIJ(1),CDIJ(1))      CCC
!     COMMON ICOUNT,MCOUNT,IORDER(500)                                   CCC
!     COMMON ICOUNT,MCOUNT,IORDER(5000)                                  DU-02
!     COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500)   CCC
!     COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000),
!    +       CODE(5000)
!-----Linear storage array and markers                                   CCC
!     COMMON MAXT,MINT,LIMIT                                             CCC
!     COMMON MARK,NEND,NSTP,NARK                                         CCC
!     COMMON MARK1,NEND1,NSTP1                                           CCC    
!     COMMON MARK2,NEND2,NSTP2                                           CCC    
!     COMMON MARK3,NEND3,NSTP3                                           CCC    
!     COMMON MARK4,NEND4,NSTP4                                           CCC    
!     COMMON TI(100000)                                                  DU-95  
!     COMMON /NEW/IDCHG,IDATOM(200),G92CHARGE(200),ESPCHARGE(200)        DU-95
!     DIMENSION IT(100000)                                               DU-95
!     EQUIVALENCE (TI(1),IT(1)) 
   
      CHARACTER(4) :: BLANK = '    '
      CHARACTER(5) :: ALPHA(6) = (/'A    ','B    ','C    ','ALPHA',&
     &                           'BETA ','GAMMA'/)
      
      INTEGER :: I, IRG, J, K, KK, KIC, KID, L, MINV 
      INTEGER :: N, NA1, NA2, NA3, NAP, NAVE, N1, N2, N3
      INTEGER :: NERR, NH3, MCT, NMV

      REAL    ::  A1NEW, ANG3(200), AM1, AM2, AM3, AS1, AS2, AS3, AVE
      REAL    ::  C2(200), C3(200), CS
      REAL    ::  D, DAC, DA1, DA2, DA3, D1, D2, D3, E
      REAL    ::  H1, H2, H3 
      REAL    ::  S2(200), S3(200), SN, VMIN, WW(100)              
                                             
!     if restart job, read minimun volume saved                        
      OPEN (UNIT=7, FILE='molpak.min', STATUS='OLD', ERR=9999)         
      KID=1                                                            
      KIC=1                                                            
      READ(7,9993)A1NEW                                                
 9993 FORMAT(5X,F8.2)                                                  
 9992 READ(7,9991,END=9996)ANGLE1(KIC),ANGLE2(KIC),ANGLE3(KIC),VOLUME(KIC)                                                      
 9991 FORMAT(9X,3F8.2,F10.3)                                           
      ICOUNT=KIC                                                       
      KIC=KIC+1                                                        
      GO TO 9992                                                       
 9996 DO I=1,ICOUNT-1                                                  
       DO J=I+1,ICOUNT                                                 
        IF(VOLUME(I).GT.VOLUME(J)) THEN                                
           CALL EXCHANGE(VOLUME(I),VOLUME(J))                          
           CALL EXCHANGE(ANGLE1(I),ANGLE1(J))                   
           CALL EXCHANGE(ANGLE2(I),ANGLE2(J))        
           CALL EXCHANGE(ANGLE3(I),ANGLE3(J))    
        END IF                                                         
       END DO                                                          
      END DO                                                           
      DO I=1,ICOUNT                                                    
        IORDER(I)=I                                                    
      END DO                                                           
 9999 CLOSE (UNIT=7)                                                   
      MCT=0                                                             
      NLOW=0                                                            
      MLOW=0                                                            
      ERMT=2.*ERM                                                       
      ERMF=4.*ERM                                                       
      AVE=0.                                                            
      NAVE=0                                                            
      HT1(8)=1000000.       ! 03-15-06                                                 
   10 FORMAT (1X,A80)                                                
!-----Does the SEEK line specify a full non-redundant angular search at  
!     the angular intervals specified in columns 7-12 ?                 
      READ(BUF,11)IRG,DAC                                            
   11 FORMAT(4X,I2,F6.1)                                             
      IF(IRG.NE.0)GO TO 15                                              
!    If the angular interval is not specified, set it to 10 degrees
      IF(DAC.LT..1)DAC=10.                                              
      DA1=DAC                                                           
      DA2=DAC                                                           
      DA3=DAC                                                           
      AS1=-90.                                                          
      AS2=-90.                                                          
      AS3=-90.                                                          
      AM1=90.                                                           
      AM2=90.                                                           
      AM3=90.                                                           
      GO TO 20                                                          
!-----The SEEK line specifies the angular ranges to be searched         
!     READ COMPLETE RANGE AND INTERVAL INSTRUCTIONS                     
   15 READ(BUF,16)IRG,AS1,AM1,DA1,AS2,AM2,DA2,AS3,AM3,DA3                                                       
   16 FORMAT(4X,I2,9F6.1)                                               
      IF(DA1.LT..1)DA1=10.                                              
      IF(DA2.LT..1)DA2=10.                                              
      IF(DA3.LT..1)DA3=10.                                              
   20 CONTINUE                                                       
      DIJ(498)=DA1                                                      
      DIJ(499)=DA2                                                      
      DIJ(500)=DA3                                                      
!     Restart a job at the section AS1,where the job stoped before     
      IF (KID.EQ.1) AS1=A1NEW+DA1                                      
      NA1 = NINT((AM1 - AS1)/DA1) + 1                                   
!-----PRINT COMPLETE RANGE AND INTERVAL INSTRUCTIONS                 
      WRITE(NQ,21)AS1,AM1,DA1,AS2,AM2,DA2,AS3,AM3,DA3       
   21 FORMAT (4X,' SEEK BEGINNING : ', 2X, 9F6.1 )                      
!-----CALCULATE SIN AND COS TABLES                                      
      E=DA2/57.296                                                      
      D=AS2/57.296                                                      
      NA2=1                                                             
      A2=AS2                                                            
   25 C2(NA2)=COS(D)                                                      
      S2(NA2)=SIN(D)                                                     
      IF(A2.GE.AM2)GO TO 30                                             
      IF(NA2.GE.200)GO TO 900                                             
      D=D+E                                                             
      A2=A2+DA2                                                         
      NA2=NA2+1                                                         
      GO TO 25                                                          
   30 E=DA3/57.296                                                      
      D=AS3/57.296                                                      
      A3=AS3                                                            
      NA3=1                                                                
   31 C3(NA3)=COS(D)                                                     
      S3(NA3)=SIN(D)                                                         
      ANG3(NA3)=A3                                                      
      IF(A3.GE.AM3)GO TO 35                                             
      IF(NA3.GE.200)GO TO 900                                              
      D=D+E                                                             
      A3=A3+DA3                                                         
      NA3=NA3+1                                                           
      GO TO 31                                                          
   35 A1=AS1                                                            
!
!-----TOP OF MOLECULE ORIENTATION =LOOP=                               
!      Rotation # 1..............
      DO 237 N1=1,NA1
   40 D=A1/57.296                                                       
      CS=COS(D)                                                         
      SN=SIN(D)                                                         
      DO 45 I=1,NMOD                                                    
         WW(I) =X(1,I)*CS-X(2,I)*SN                                  
         V(2,I)=X(1,I)*SN+X(2,I)*CS                                  
   45 CONTINUE                                                       
      IF(IPR.LT.4)GO TO 49                                           
      WRITE(NQ,46)                                                   
   46 FORMAT(1H0,'POSITIONS AFTER FIRST ROTATION')                   
      DO 48 I=1,NMOD                                                 
      WRITE(NQ,47)AT(I),AN(I),WW(I),V(2,I),X(3,I)                       
   47 FORMAT(1X,A2,A4,3F8.3)                                         
   48 CONTINUE                                                       
49    A2=AS2                                                            
!----Rotation # 2 ..................
      DO 236 N2=1,NA2                                                   
   50 DO 51 I=1,NMOD                                                    
      V(1,I)=WW(I)*C2(N2)+X(3,I)*S2(N2)                                  
      W(3,I)=-WW(I)*S2(N2)+X(3,I)*C2(N2)                                
   51 CONTINUE                                                       
      IF(IPR.LT.4)GO TO 54                                           
      WRITE(NQ,52)                                                   
   52 FORMAT(1H0,'POSITIONS AFTER SECOND ROTATION')                  
      DO 53 I=1,NMOD                                                 
      WRITE(NQ,47)AT(I),AN(I),V(1,I),V(2,I),W(3,I)                   
   53 CONTINUE                                                       
!-----Calculate minimum axis-3 length                                   
   54 N=MARK1                                                           
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
      IF(N.LE.MARK)GO TO 930                                            
      NEND=N-1                                                          
      H3M=7.                                                            
      ER=ERM                                                            
      CALL MINHI(H3M,0.,100.,7.,NNE)                                    
      IF(KCT)911,57,921                                              
   57 H3N=.5*H3M                                                        
      IF(IPR.LT.3)GO TO 69                                              
      WRITE(NQ,58)A1,A2,H3M                                             
   58 FORMAT(5H A1 =,F5.1,6H, A2 =,F5.1,7H, H3M =,F7.3)                 
   69 A3=AS3                                                            
      NP2=1                                                                 
!----Rotation # 3 ...................
      DO 199 N3=1,NA3                                                    
   70 DIJ(NP2)=.001                                                     
      CDIJ(NP2+230)=BLANK                                                 
      V3=10000.                                                         
      DO 71 I=1,NMOD                                                    
         W(1,I)=V(1,I)*C3(N3)-V(2,I)*S3(N3)                         
         W(2,I)=V(1,I)*S3(N3)+V(2,I)*C3(N3)                           
   71 CONTINUE                                                       
      IF(IPR.LT.4)GO TO 74                                           
      WRITE(NQ,72)                                                   
   72 FORMAT(1H0,'POSITIONS AFTER THIRD ROTATION')                   
      DO 73 I=1,NMOD                                                 
      WRITE(NQ,47)AT(I),AN(I),W(1,I),W(2,I),W(3,I)                   
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
      NH3 = 0
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
      IF(N.LE.MARK1)GO TO 930                                           
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
      IF(N.LE.MARK)GO TO 930                                            
      NEND=N-1                                                          
      ER=ERM                                                         
      H1M=7.                                                         
      CALL MINHI(H1M,0.,100.,7.,NNE)                                    
      IF(KCT)912,81,922                                              
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
      IF(N.LE.MARK)GO TO 930                                            
      NEND=N-1                                                          
      ER=ERM                                                         
      H2M=7.                                                         
      CALL MINHI(H2M,0.,100.,7.,NNE)                                    
      IF(KCT)913,83,923                                              
   83 H2N=H2M/2.                                                     
      IF(IPR.LT.3)GO TO 90                                           
      WRITE(NQ,85)H1M,H2M                                            
   85 FORMAT(1H ,'H1M =',F8.3,', H2M =',F8.3)                        
!-----Transfer to specified coordination sphere finding subroutines      
   90 GO TO(91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,&   
     & 108,109,110,111,112,113,114,115,116,117,999),JTR(1)               
   91 CALL ZTRTA                                                        
      GO TO(92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,&  
     & 109,110,111,112,113,114,115,116,117,198),JTR(2)                   
   92 CALL ZTRTB                                                        
      GO TO(93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,& 
     & 110,111,112,113,114,115,116,117,198),JTR(3)                       
   93 CALL ZTRTC                                                        
      GO TO(94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,&    
     & 110,111,112,113,114,115,116,117,198),JTR(4)                       
   94 CALL ZFRTA                                                        
      GO TO(95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,&   
     & 111,112,113,114,115,116,117,198),JTR(5)                           
   95 CALL ZFRTB                                                        
      GO TO(96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,&  
     & 112,113,114,115,116,117,198),JTR(6)                               
   96 CALL ZFRTC                                                        
      GO TO(97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,& 
     & 113,114,115,116,117,198),JTR(7)                                   
   97 CALL ZORTA                                                        
      GO TO(98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,&    
     & 113,114,115,116,117,198),JTR(8)                                   
   98 CALL ZORTB                                                        
      GO TO(99,100,101,102,103,104,105,106,107,108,109,110,111,&           
     & 112,113,114,115,116,117,198),JTR(9)                           
   99 CALL ZORTC                                                        
      GO TO(100,101,102,103,104,105,106,107,108,109,110,111,&      
     & 112,113,114,115,116,117,198),JTR(10)                          
  100 CALL ZERTA                                                     
      GO TO(101,102,103,104,105,106,107,108,109,110,111,&               
     & 112,113,114,115,116,117,198),JTR(11)                              
  101 CALL ZERTB                                                         
      GO TO(102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,&  
     & 117,198),JTR(12)                                                  
  102 CALL ZNRTA                                                     
      GO TO(103,104,105,106,107,108,109,110,111,112,113,114,115,116,&      
     & 117,198),JTR(13)                                                  
  103 CALL ZNRTC                                                     
      GO TO(104,105,106,107,108,109,110,111,112,113,114,115,116,&          
     & 117,198),JTR(14)                                                  
  104 CALL ZNRTB                                                         
      GO TO(105,106,107,108,109,110,111,112,113,114,115,116,&             
     & 117,198),JTR(15)                                                  
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
      IF (N3 .LT. NA3) NP2=NP2+1                                        
      A3=A3+DA3                                                         
  199 CONTINUE                ! bottom of rotn # 3
  213 IF(A2.GT.AS2)GO TO 220                                            
      IF(A1.GT.AS1)GO TO 215                                            
      WRITE(NQ,214)                                                     
  214 FORMAT('UNIT CELL VOLUME PER MOL. VS. MOLECULE ROTATION ABOUT',&    
     &       ' AXIS 2 DOWN, AXIS 3 ACROSS')                          
  215 WRITE(NQ,217)A1                                                   
  217 FORMAT (1H0,30X,4HA1 =,F6.1)                                      
      PRINT 218, A1                                                  
  218 FORMAT (' Start MOLPAK section A1 =',F6.1)                     
      WRITE(NQ,219)(ANG3(I),I=1,NA3)                                       
  219 FORMAT(1H0,6H  A2  ,6H  A3  ,6H MIN V,200F6.1)                   
  220 IF(NEW.NE.0)WRITE(NQ,219)(ANG3(I),I=1,NA3)                         
      IF(MCOUNT.LE.0)GO TO 230                                       
!-----Add cell volumes/molecule found at this A2 to collection of MCOUNT 
!      smallest volumes                                                 
      DO 229 I=1,NP2                                                     
!     if restart job, eliminate the same parts                         
      IF(KID.EQ.1) THEN                                                
       DO K=1,ICOUNT                                                   
        IF (A1.EQ.ANGLE1(K).AND.A2.EQ.ANGLE2(K).AND.ANG3(I).EQ.&        
     &      ANGLE3(K)) GO TO 229                                       
       END DO                                                          
      END IF                                                           
      IF(ICOUNT.GE.MCOUNT)GO TO 221                                  
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
!      determine minimum volume for this value of A2, and fix          
      VMIN=1000000.          ! 03-15-06                                         
  230 DO 231 I=1,NP2                                                     
      NAVE=NAVE+1                                                       
      AVE=AVE+DIJ(I)                                                    
      IF(DIJ(I).GT.VMIN)GO TO 231                                    
      VMIN=DIJ(I)                                                    
      MINV=DIJ(I)+.5                                                 
      NMV=I                                                            
  231 IJD(I)=DIJ(I)+.5                                                  
!-----Print line of volumes/molecule found at THIS value of A2          
      WRITE(NQ,235)A2,ANG3(NMV),MINV,CDIJ(NMV+230),(IJD(I),&
     &             CDIJ(I+230),I=1,NP2)                                                          
  235 FORMAT(1H ,2F6.1,I4,A2,200(I4,A2))                                
      A2=A2+DA2                                                         
  236 CONTINUE                 ! bottom of rotn # 2                     
!     after a section for searching minimun volume finished, save      
!     minimun volume to 'molpak.vol' file (unit=18)                    
      OPEN (UNIT=18, FILE='molpak.vol', STATUS='UNKNOWN')              
      WRITE(18,227)A1                                                  
      DO J=1,MCOUNT                                                    
        WRITE(18,225)J,ANGLE1(J),ANGLE2(J),ANGLE3(J),VOLUME(J)         
      END DO                                                           
      CLOSE (UNIT=18)                                                  
  225 FORMAT(1X,I8,3F8.2,F10.3)                                        
  227 FORMAT(1X,'A1 =',F8.2)                                           
      A1 = A1 + DA1                                                     
  237 CONTINUE                      ! bottom of rotn # 1                
!-----Calculate and print average volume and ratio to smallest volume   
      IF(NAVE.EQ.0)NAVE=1                                            
      AVE=AVE/NAVE                                                      
      IF(HT1(8).EQ.0.)HT1(8)=1.                                      
      D=AVE/HT1(8)                                                      
      WRITE(NQ,238)AVE,D                                                
  238 FORMAT(10H0AVERAGE =,F9.3,12H, AVE/BEST =,F6.3)                   
!-----Print out parameters of smallest volume cell found by this SEEK   
!      instruction                                                      
  370 IF(HT1(8).GT.10000.)GO TO 1000                                    
      WRITE(NQ,371)(HT1(I),I=10,12),NLOW,CHT1(9),HT1(8),ERM             
  371 FORMAT(24H0SMALLEST CELL AT   A1 =,F6.1,7H,  A2 =,F6.1,7H,  A3 =,& 
     & F6.1,7H   CASE,I3,1X,A2,8H,  V/M =,F7.2,6H,  E =,F6.3/1H )        
      E=90.                                                             
      D=0.                                                              
      WRITE(21,375)(HT1(I),I=10,12),HT1(8),CHT1(9)                      
  375 FORMAT(/4X,'MOLPAK ANGLES & VOLUME = ',3F6.1,F7.2,2X,A2)       
      GOTO(380,380,385,385,385, 385,380,390,390,390, 390,395,395,395,395,& 
     &     390,390,390,390,390, 395,385,385,385,395, 395,395,395,395,390,&
     &     390,395,395,395,395, 395),NLOW                                
!-----ZTRT result for triclinic structures                               
  380 D2=SQRT(HT1(2)**2+HT1(3)**2)                                   
      IF(HT1(14).LE.HT1(3)-HT1(14))GO TO 381                         
      HT1(14)=HT1(14)-HT1(3)                                         
      HT1(13)=HT1(13)-HT1(2)                                         
  381 IF(HT1(13).GE.HT1(1)-HT1(13))HT1(13)=HT1(13)-HT1(1)            
      D3=SQRT(HT1(4)**2+HT1(13)**2+HT1(14)**2)                       
      A1=57.296*ACOS((HT1(2)*HT1(13)+HT1(3)*HT1(14))/(D2*D3))        
      A2=57.296*ACOS(HT1(13)/D3)                                     
      A3=57.296*ACOS(HT1(2)/D2)                                      
      WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),A1,HT1(1),D,D            
      WRITE(NQ,386)ALPHA(2),D2,ALPHA(5),A2,HT1(2),HT1(3),D           
      WRITE(NQ,386)ALPHA(3),D3,ALPHA(6),A3,HT1(13),HT1(14),HT1(4)    
      IF(NLOW.EQ.1)GO TO 1000                                            
      GO TO 396                                                      
!-----ZTRT or ZSRT result for Z=2 monoclinic structures             
  385 D2=SQRT(HT1(2)**2+HT1(3)**2)                                   
      A3=57.296*ACOS(HT1(2)/D2)                                      
      WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),E,HT1(1),D,D             
      WRITE(NQ,386)ALPHA(2),D2,ALPHA(5),E,HT1(2),HT1(3),D            
      WRITE(NQ,386)ALPHA(3),HT1(4),ALPHA(6),A3,D,D,HT1(4)            
  386 FORMAT(1H ,3X,A1,2H =,F7.3,3X,A5,2H =,F7.1,8X,3F8.3)           
      IF(NLOW.EQ.22)GO TO 1000                                           
      GO TO 396                                                      
!-----ZFRT or ZNRT result for Z=4 monoclinic structures                 
  390 D2=SQRT(HT1(2)**2+HT1(3)**2)                                   
      A3=57.296*ACOS(HT1(2)/D2)                                      
      WRITE(NQ,386)ALPHA(1),HT1(1),ALPHA(4),E,HT1(1),D,D             
      WRITE(NQ,386)ALPHA(2),D2,ALPHA(5),E,HT1(2),HT1(3),D            
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
  903 FORMAT(27H0STORAGE EXCEEDED BY MOLPAK,5X,4(I9,I3),I9)          
      GO TO 999                                                         
  911 NERR=1                                                         
      GO TO 918                                                      
  912 NERR=2                                                         
      GO TO 918                                                      
  913 NERR=3                                                         
  918 WRITE(NQ,919)NERR                                              
      WRITE(17,919)NERR                                              
  919 FORMAT(' PARAMETER TOO SMALL - SEEK Location',I2)              
      GO TO 999                                                      
  921 NERR=1                                                         
      GO TO 928                                                      
  922 NERR=2                                                         
      GO TO 928                                                      
  923 NERR=3                                                         
  928 WRITE(NQ,929)NERR                                              
      WRITE(17,929)NERR                                              
  929 FORMAT(' PARAMETER TOO LARGE - SEEK Location',I2)              
      GO TO 999                                                      
  930 WRITE(NQ,931)KCT,MCT,N,MARK1,MARK2,MARK,NARK                      
      WRITE(17,931)KCT,MCT,N,MARK1,MARK2,MARK,NARK                      
  931 FORMAT(16H0NO INTERACTIONS,2I3,5I9)                               
      GO TO 999                                                         
  999 KILL=1                                                            
 1000 RETURN                                                            
      END SUBROUTINE SEEK                                                            
