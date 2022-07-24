      SUBROUTINE ZNRTB                                                      
!                                                                        
!     For EA, EB, EC, ED, EE, EF, EG                                                   
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
      CHARACTER(2) :: ALPHA(7) = (/'EA','EB','EC','ED','EE','EF','EG'/)
!
!      DIMENSION ALPHA(7),WI(7),WJ(7),WK(7),WL(2)                            
!      DATA ALPHA/'EA','EB','EC','ED','EE','EF','EG'/                     
!      DATA WI/-1.,1.,1.,1.,-1.,-1.,-1./                                     
!      DATA WJ/1.,-1.,1.,-1.,1.,-1.,-1./                                     
!      DATA WK/1.,1.,-1.,-1.,-1.,1.,-1./                                     
!
      INTEGER :: I, J, JE, K, KE, L, LCT, LCTM
      INTEGER :: MCT, N, NP, NPM
!
      REAL :: AREA
      REAL :: D, D1, D2, D3
      REAL :: DX3IH
      REAL :: H1I, H1IH
      REAL :: H2IH, H2IQ, H2X
      REAL :: H3I, H3IH
      REAL :: VX
      REAL :: WI(7) =(/-1.,1.,1.,1.,-1.,-1.,-1./) 
      REAL :: WJ(7) =(/1.,-1.,1.,-1.,1.,-1.,-1./)
      REAL :: WK(7) =(/1.,1.,-1.,-1.,-1.,1.,-1./)
      REAL :: WL(7)
      REAL :: X1IH, X2X
      REAL :: X3I, X3IH
!
!-----This subroutine treats face-centered orthorhombic space groups 
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)(ITR(I),I=51,57)                                        
    5 FORMAT(1X,'ZNRTB called',7I3)                                      
!-----Find minimum area F-centered plane-1,3 by determining H1I and H3I  
   10 N=MARK2                                                               
      NSTP2=4                                                            
      DO 19 I=1,NMOD                                                        
      DO 18 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      D2=W(2,J)-W(2,I)                                                      
      IF(ABS(D2).GT.CN(K))GO TO 18                                          
      IT(N)=K                                                            
      TI(N+1)=D2**2                                                      
      TI(N+2)=W(1,J)-W(1,I)                                              
      TI(N+3)=W(3,J)-W(3,I)                                              
      N=N+NSTP2                                                          
      IF(N.GT.LIMIT)GO TO 902                                            
   18 CONTINUE                                                           
   19 CONTINUE                                                           
      IF(N.LE.MARK2)GO TO 920                                            
      NEND2=N-1                                                          
      MARK=N                                                             
      KE=NE                                                                 
      JE=NV/2                                                               
      MCT=0                                                              
      X3I=H3M                                                               
      X3IH=H3N                                                           
      X1IH=1.5*H1N                                                          
      DX3IH=H3M/NV                                                       
      ER=ERMT                                                               
!-----Set dummy starting face half area                                  
      AREA=10000.                                                           
   20 N=MARK                                                             
      DO 25 I=MARK2,NEND2,NSTP2                                          
      K=IT(I)                                                            
      L=K+10                                                             
      D3=TI(I+3)                                                         
   21 D3=D3+X3I                                                          
      IF(D3.LE.CN(K))GO TO 21                                            
   22 D3=D3-X3I                                                          
      IF(D3.GT.CN(K))GO TO 22                                            
      IF(D3.LT.-CN(K))GO TO 25                                           
      D=TI(I+1)+D3**2                                                    
      IF(D.GT.CN(L))GO TO 22                                             
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 22                                                           
   25 CONTINUE                                                           
      NARK=N                                                             
      DO 29 I=MARK2,NEND2,NSTP2                                          
      K=IT(I)                                                            
      L=K+10                                                             
      D3=X3IH+TI(I+3)                                                       
   26 D3=D3+X3I                                                             
      IF(D3.LE.CN(K))GO TO 26                                               
   27 D3=D3-X3I                                                             
      IF(D3.GT.CN(K))GO TO 27                                               
      IF(D3.LT.-CN(K))GO TO 29                                              
      D=TI(I+1)+D3**2                                                       
      IF(D.GT.CN(L))GO TO 27                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 27                                                              
   29 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(X1IH,0.,100.,H1M,KE)                                       
      D=X1IH*X3I                                                         
      IF(D.GE.AREA)GO TO 31                                              
      AREA=D                                                             
      H1IH=X1IH                                                          
      H3IH=X3IH                                                          
      H1I=2.*H1IH                                                        
      H3I=2.*H3IH                                                        
      IF(MCT)33,32,33                                                    
   31 IF(MCT)33,32,35                                                    
   32 IF(JE.GE.NV)GO TO 34                                               
      JE=JE+1                                                            
      X3IH=X3IH+DX3IH                                                    
      X3I=2.*X3IH                                                        
      GO TO 20                                                           
   33 IF(JE.GE.KE)GO TO 36                                               
   34 JE=2*JE                                                            
      DX3IH=H3M/JE                                                       
      X3IH=H3IH+DX3IH                                                    
      X3I=2.*X3IH                                                        
      MCT=1                                                              
      GO TO 20                                                           
   35 X3IH=H3IH-DX3IH                                                    
      IF(X3IH.LT.H3N)GO TO 33                                            
      X3I=2.*X3IH                                                        
      MCT=-1                                                             
      GO TO 20                                                           
   36 IF(JE.GE.NNE)GO TO 37                                              
      KE=NNE                                                             
      GO TO 34                                                           
   37 WL(1)=.25*H3I                                                         
      WL(2)=.75*H3I                                                         
      IF(IPR.LT.3)GO TO 40                                                  
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,38)A3,H1I,H3I                                                
   38 FORMAT(12X,4HA3 =,F5.1,7H, H1I =,F7.3,7H, H3I =,F7.3)                 
!-----Calculate H2IH                                                     
   40 N=MARK                                                                
      DO 46 I=1,NMOD                                                        
      DO 45 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=W(1,J)-W(1,I)                                                      
      D3=W(3,J)-W(3,I)                                                   
   41 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 41                                            
   42 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 42                                            
      IF(D1.LT.-CN(K))GO TO 45                                           
   43 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 43                                               
   44 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 44                                               
      IF(D3.LT.-CN(K))GO TO 42                                              
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 44                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)-W(2,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 44                                                              
   45 CONTINUE                                                              
   46 CONTINUE                                                              
      NARK=N                                                             
      DO 52 I=1,NMOD                                                        
      DO 51 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=W(1,J)-W(1,I)                                                      
      D3=H3IH+W(3,J)-W(3,I)                                                 
   47 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 47                                            
   48 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 48                                               
      IF(D1.LT.-CN(K))GO TO 51                                           
   49 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 49                                               
   50 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 50                                               
      IF(D3.LT.-CN(K))GO TO 48                                              
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 50                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)-W(2,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 50                                                              
   51 CONTINUE                                                              
   52 CONTINUE                                                              
      DO 58 I=1,NMOD                                                        
      DO 57 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=H1IH+W(1,J)-W(1,I)                                                 
      D3=W(3,J)-W(3,I)                                                      
   53 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 53                                            
   54 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 54                                               
      IF(D1.LT.-CN(K))GO TO 57                                           
   55 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 55                                               
   56 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 56                                               
      IF(D3.LT.-CN(K))GO TO 54                                              
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 56                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)-W(2,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 56                                                              
   57 CONTINUE                                                              
   58 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      H2IH=1.5*H2N                                                          
      CALL MINHI(H2IH,H2N,100.,H2M,KE)                                      
      H2IQ=H2IH/2.                                                          
      IF(IPR.LT.3)GO TO 60                                                  
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,59)H2IH                                                      
   59 FORMAT(23X,6HH2IH =,F7.3)                                             
   60 LCT=0                                                              
      LCTM=0                                                             
      H2X=100.                                                           
   61 LCT=LCT+1                                                          
      NP=1                                                               
      GO TO(65,70,75,80,85,90,95,115),LCT                                
   65 GO TO(66,71,76,81,86,91,96,115),ITR(51)                            
!-----P1 molecule                                                        
   66 LCT=1                                                              
      GO TO 100                                                          
   70 GO TO(71,76,81,86,91,96,115),ITR(52)                               
!-----P2 molecule                                                        
   71 LCT=2                                                              
      GO TO 100                                                          
   75 GO TO(76,81,86,91,96,115),ITR(53)                                  
!-----P3 molecule                                                        
   76 LCT=3                                                              
      GO TO 100                                                          
   80 GO TO(81,86,91,96,115),ITR(54)                                     
!-----A1 molecule                                                        
   81 LCT=4                                                              
      GO TO 100                                                          
   85 GO TO(86,91,96,115),ITR(55)                                        
!-----A2 molecule                                                        
   86 LCT=5                                                              
      GO TO 100                                                          
   90 IF(ITR(56)-2)91,96,115                                             
!-----A3 molecule                                                        
   91 LCT=6                                                              
      GO TO 100                                                          
   95 IF(ITR(57).NE.1)GO TO 115                                          
!-----C molecule                                                         
   96 LCT=7                                                              
!-----Calculate axis-2 one-fourth length                                 
  100 IF(LCTM.LT.0)GO TO 61                                              
      NP=1                                                                  
      X2X=1.5*H2IQ                                                          
      ER=ERMF                                                               
      MARK=MARK2                                                            
      NARK=0                                                             
  101 N=MARK2                                                               
      DO 109 I=1,NMOD                                                       
      DO 108 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=H1IH+W(1,J)*WI(LCT)-W(1,I)                                         
      D3=WL(NP)+W(3,J)*WK(LCT)-W(3,I)                                       
  102 D1=D1+H1I                                                             
      D3=D3+H3IH                                                            
      IF(D1.LE.CN(K))GO TO 102                                              
  103 D1=D1-H1I                                                             
      D3=D3-H3IH                                                            
      IF(D1.GT.CN(K))GO TO 103                                              
      IF(D1.LT.-CN(K))GO TO 108                                             
  104 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 104                                              
  105 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 105                                              
      IF(D3.LT.-CN(K))GO TO 103                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 105                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)*WJ(LCT)-W(2,I)                                         
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 105                                                             
  108 CONTINUE                                                              
  109 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 110                                                
      NEND=N-1                                                              
      CALL MINHI(X2X,H2IQ,H2X,H2IH,KE)                                      
      IF(KCT)110,111,112                                                    
  110 H2X=H2IQ                                                              
      LCTM=-LCT                                                             
      NPM=NP                                                                
      GO TO 120                                                             
  111 H2X=X2X                                                               
      LCTM=LCT                                                              
      NPM=NP                                                                
  112 IF(NP.GE.2)GO TO 61                                                   
      NP=2                                                                  
      GO TO 101                                                             
  115 IF(LCTM.EQ.0)GO TO 1000                                            
  120 VX=H1IH*H2X*H3I                                                       
      LCTM=IABS(LCTM)                                                    
      IF(IPR.LT.2)GO TO 125                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,121)A3,ALPHA(LCTM),H1IH,H2X,H3I,VX,WL(NPM)                   
  121 FORMAT(8X,F6.1,1X,A2,3F9.3,F9.2,5X,3F9.3,5X,3F9.3)                    
  125 IF(V3-VX.LT..1)GO TO 1000                                             
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCTM)                                              
      IF(V3.GE.HT1(8))GO TO 1000                                            
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=4.*H2X                                                         
      HT1(4)=H3I                                                            
      HT1(5)=H1IH                                                           
      HT1(6)=H2X                                                            
      HT1(7)=WL(NPM)                                                        
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCTM)                                                   
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=21                                                               
      MLOW=LCTM+50                                                          
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZNRTB,5X,3(I9,I3),I9)                  
      GO TO 999                                                             
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)KCT,MARK1,MARK2,NEND2,MARK                               
  921 FORMAT(16H0NO INTERACTIONS,I3,4I9)                                    
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZNRTB           
