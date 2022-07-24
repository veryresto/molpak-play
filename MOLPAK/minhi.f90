      SUBROUTINE MINHI (HX,HMIN,HMAX,RANGE,IE)
 
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
!      DIMENSION IJD(500),CDIJ(500),CHT1(14)                              
!      EQUIVALENCE (DIJ(1),IJD(1)),(DIJ(1),CDIJ(1)),(HT1(1),CHT1(1))      
!     COMMON ICOUNT,MCOUNT,IORDER(500)                                   
!      COMMON ICOUNT,MCOUNT,IORDER(5000)                                  
!     COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500)   
!      COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000), &
!     &       CODE(5000)
!----Linear storage array and markers                                   
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

      INTEGER :: I, IE, J, JE, K, KARK, LARK, LEND, NCT, MCT, NERR

      REAL :: D, DA,  DXX, EGY, HL, HMIN, HMAX, HT, HX, RANGE, TXX, XX 
!
      IF(IPR.LT.6)GO TO 3                                                
      WRITE(NQ,1)HX,HMIN,HMAX,RANGE,IE,MARK,NARK,NEND                    
    1 FORMAT(1X,'MINHI called',4F8.3,4I6)                                
!     If there are no interactions in the designated table, set HX to    
!      its minimum possible value, HMIN, and exit.                       
    3 IF(NEND.LE.MARK)GO TO 43                                              
!----Assume no atoms at double standoff are present                     
      LARK=0                                                             
      KARK=MARK                                                          
!----If NARK is set within the table, all interactions before NARK are  
!      calculated at double standoff.                                    
      IF((NARK.LE.MARK).OR.(NARK.GE.NEND))GO TO 5                        
      LARK=1                                                             
      LEND=NARK-1                                                        
      KARK=NARK                                                          
!     Place the standoff variable, HX, in the allowed range - negative   
!      values are treated as absolutes.                                  
    5 HL=SIGN(1.,RANGE)                                                     
      IF(HL*(HMAX-HX).LT.0.)HX=HMAX                                         
      IF(HL*(HX-HMIN).LT.0.)HX=HMIN                                         
      XX=HX                                                                 
!----Set flag to indicate that the original value of HX is being tested 
      NCT=0                                                                 
!----Set flag to indicate that the original scan of possible HX values  
!      is in progress - an interaction energy near the designated        
!      threshold value, ER, has not yet been found.                      
      MCT=0                                                                 
      JE=4                                                                  
      DXX=RANGE/JE                                                          
   10 EGY=0.                                                                
      IF(LARK.LE.0)GO TO 16                                              
      TXX=2.*XX                                                             
      DO 15 I=MARK,LEND,NSTP                                                
      K=IT(I)                                                               
      DA=TXX+TI(I+2)                                                        
      IF(ABS(DA).GT.CN(K))GO TO 15                                          
      D=TI(I+1)+DA**2                                                       
      IF(D.GT.CN(K+10))GO TO 15                                             
      IF(CN(K+11).EQ.0.)GO TO 900                                        
      J=IFIX((CN(K+10)-D)/CN(K+11))+1                                       
      IF(J.GT.9)GO TO 27                                                    
      EGY=EGY+CN(K+J)                                                       
      IF(EGY.GT.ER)GO TO 27                                                 
   15 CONTINUE                                                              
   16 DO 19 I=KARK,NEND,NSTP                                                
      K=IT(I)                                                               
      DA=XX+TI(I+2)                                                         
      IF(ABS(DA).GT.CN(K))GO TO 19                                          
      D=TI(I+1)+DA**2                                                       
      IF(D.GT.CN(K+10))GO TO 19                                             
      IF(CN(K+11).EQ.0.)GO TO 901                                        
      J=IFIX((CN(K+10)-D)/CN(K+11))+1                                       
      IF(J.GT.9)GO TO 27                                                    
      EGY=EGY+CN(K+J)                                                       
      IF(EGY.GT.ER)GO TO 27                                                 
   19 CONTINUE                                                              
      IF(IPR.LT.7)GO TO 24                                               
      WRITE(NQ,20)XX,HX,EGY,ER,JE,IE                                     
   20 FORMAT(1X,'too large',4F8.3,2I6)                                   
!----The total interaction energy is below the designated threshold,    
!      set the best standoff value, HX, to the trial value, XX.          
   24 HX=XX                                                                 
      IF(MCT.EQ.1)GO TO 34                                                  
      IF(NCT)26,25,33                                                       
!----Original value of HX is too large, set flag to step inward.        
   25 NCT=-1                                                                
   26 XX=HX-DXX                                                             
      IF(HL*(XX-HMIN).GE.0.)GO TO 10                                        
      IF(JE.GE.IE)GO TO 43                                                  
      JE=2*JE                                                               
      DXX=RANGE/JE                                                          
      GO TO 26                                                              
   27 IF(IPR.LT.7)GO TO 29                                               
      WRITE(NQ,28)XX,HX,EGY,ER,JE,IE                                     
   28 FORMAT(1X,'too small',4F8.3,2I6)                                   
!----The test value of HX, XX, is too small - the energy is too large   
   29 IF(MCT.EQ.1)GO TO 34                                                  
      IF(NCT)33,30,31                                                       
!----Original value of HX is too small, set flag to step outward        
   30 NCT=1                                                                 
!----Set HT equal to the last value of XX for which energy was calcd.   
   31 HT=XX                                                                 
   32 XX=HT+DXX                                                             
      IF(HL*(HMAX-XX).GE.0.)GO TO 10                                        
      IF(JE.GE.IE)GO TO 45                                                  
      JE=2*JE                                                               
      DXX=RANGE/JE                                                          
      GO TO 32                                                              
!----Set flag to indicate that the trial standoff , XX, has produced an 
!      interaction energy which crossed the threshold value, ER.         
   33 MCT=1                                                                 
!     The following steps are the final refinement of the standoff, HX   
   34 IF(JE.GE.IE)GO TO 49                                                  
      JE=2*JE                                                               
      DXX=RANGE/JE                                                          
      XX=HX-DXX                                                             
      GO TO 10                                                              
   43 HX=HMIN                                                               
      KCT=-1                                                                
      GO TO 50                                                              
   45 HX=HMAX                                                               
      KCT=1                                                              
      GO TO 50                                                              
!----Check whether best value is significantly different from limits    
   49 IF(ABS(HX-HMAX).LE.ABS(DXX))GO TO 45                               
      IF(ABS(HX-HMIN).LE.ABS(DXX))GO TO 43                               
!----Set flag to show that an intermediate HX value has been found      
      KCT=0                                                                 
      GO TO 50                                                           
  900 NERR=1                                                             
      GO TO 905                                                          
  901 NERR=2                                                             
  905 WRITE(NQ,906)NERR,MARK,NEND,MARK1,MARK2,NARK                       
  906 FORMAT(' MINHI divide error',I2,5I6)                               
      STOP                                                               
   50 RETURN                                                                
      END SUBROUTINE MINHI  
