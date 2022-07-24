      SUBROUTINE ZSRTD                                                      
!                                                                                                                                   
!     This subroutine finds structures in space groups Pnn2, Pba2, and   
!      Pma2 when the molecule contains a two-fold axis. The coordination 
!      sphere contains 6 I molecules in plane-1,3.  For SH, SI and SJ.                     
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
!      CHARACTER*2 ALPHA                                                  
!      DIMENSION ALPHA(3)                                                    
!      DIMENSION A(6),C(6)                                                
!      DATA ALPHA/'SH','SI','SJ'/                                            
!
      CHARACTER(2) :: ALPHA(3)= (/'SH','SI','SJ'/)
!
      INTEGER :: I, J, K, KE, L, LCT, N
!
      REAL :: A(6), C(6)
      REAL :: D, D1, D3, VX
      REAL :: H1P1, H2I, H2P1, H3P1 
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 300                                              
      WRITE(NQ,5)(ITR(I),I=65,67)                                        
    5 FORMAT(1X,'ZTRTB called',3I3)                                      
  300 ER=ERMF                                                               
      MARK=MARK2                                                         
      N=MARK                                                                
      KE=NNE                                                             
!-----Collect I molecule positions for double standoff                   
      DO 309 I=1,NMOD                                                       
      DO 308 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=W(1,J)-W(1,I)                                                      
      D3=W(3,J)-W(3,I)                                                      
  301 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 301                                              
  302 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 302                                              
      IF(D3.LT.-CN(K))GO TO 308                                             
  303 D1=D1+H1M                                                             
      IF(D1.LE.CN(K))GO TO 303                                              
  304 D1=D1-H1M                                                             
      IF(D1.GT.CN(K))GO TO 304                                              
      IF(D1.LT.-CN(K))GO TO 302                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 304                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)-W(2,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 304                                                             
  308 CONTINUE                                                              
  309 CONTINUE                                                              
      NARK=N                                                             
!-----Go to first designated ZTRTB space group                           
      IF(ITR(65)-2)320,380,390                                           
!     Code SH - P1(1/2,1/2,1/2)                                          
  320 LCT=1                                                              
      H1P1=H1N                                                           
      H3P1=H3N                                                           
!-----Calculate standoff of grid of P1 molecules in axis-2 direction     
  330 N=NARK                                                             
      DO 339 I=1,NMOD                                                       
      DO 338 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=H1P1-W(1,J)-W(1,I)                                                 
      D3=H3P1+W(3,J)-W(3,I)                                                 
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
      CALL MINHI(H2P1,0.,100.,H2M,KE)                                       
      IF(KCT)904,340,906                                                 
  340 H2I=2.*H2P1                                                        
      VX=H1N*H2I*H3M                                                     
      IF(V3-VX.LT..01)GO TO 360                                             
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 360                                             
      HT1(1)=H1M                                                            
      HT1(2)=0.                                                          
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                         
      HT1(5)=0.                                                             
      HT1(6)=0.                                                             
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                  
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=26                                                             
      MLOW=LCT+64                                                        
  360 IF(IPR.LT.2)GO TO 370                                              
      A(1)=H1M                                                           
      A(2)=H2I                                                           
      A(3)=H3M                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=0.                                                            
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      A(6)=H3P1                                                          
      A(5)=H2P1                                                          
      A(4)=H1P1                                                          
!-----Then to fractional coordinates                                     
      D=SQRT(H1P1**2+H2P1**2+H3P1**2)                                    
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=90.                                                           
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,365)ALPHA(LCT),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)      

  365 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3/1X,'  COSINES' &     
     &,3F7.4,'  ANGLES',3F7.2)                                           
      WRITE(NQ,366)H1P1,H2P1,H3P1,(A(I),I=4,6),D                         
  366 FORMAT(1X,'P3 AT',3F7.3,';  ',3F7.4,';  D =',F6.2)          
  370 IF(LCT-2)375,385,1000                                              
!-----Are Codes SI or SJ to be calculated?                               
  375 IF(ITR(66)-2)380,390,1000                                          
!     Code SI - P1(1/2,0,1/2)                                            
  380 LCT=2                                                              
      H1P1=H1N                                                           
      H3P1=0.                                                            
      GO TO 330                                                          
!     Is Code SJ to be calculated?                                       
  385 IF(ITR(67).GE.2)GO TO 1000                                         
!     Code SJ - P1(0,0,1/2)                                              
  390 LCT=3                                                              
      H1P1=0.                                                            
      H3P1=0.                                                            
      GO TO 330                                                          
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTD,5X,3(I9,I3),I9)                  
      GO TO 999                                                             
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)                                                      
  905 FORMAT('0PARAMETER TOO SMALL')                                     
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT('0PARAMETER TOO LARGE')                                     
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)                                                         
  921 FORMAT(16H0NO INTERACTIONS)                                           
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                             
      END SUBROUTINE ZSRTD  

