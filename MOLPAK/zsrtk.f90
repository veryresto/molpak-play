      SUBROUTINE ZSRTK     
!                                                                                                                                                
!-----This routine finds structures in space group Pbcn when the         
!      molecule contains a center of symmetry.  The coordination sphere  
!      contains 2 I molecules on axis-3. For SV.                                
!
      USE molpakCommonMod
!
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
!      CHARACTER*2 ALPHA,SYM                                              
!      DIMENSION SM(6),SYM(3)                                                
!      DATA ALPHA/'SV'/                                                      
!      DATA SYM/'A1','A2','A3'/                                           
!
      CHARACTER(2) :: ALPHA = 'SV' 
      CHARACTER(2) :: SYM(3) = (/'A1','A2','A3'/)
!
      INTEGER :: I, J, K, L, N
!
      REAL :: D, D1, D2, D3
      REAL :: H1I, H1IH, H2I, H2IH
      REAL :: SM(6), VX
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)                                                        
    5 FORMAT(1X,'ZSRTK called')                                          
!-----TRIPLE AXIS CASE                                                      
!     Calculate the offset distance of a line of A3 screw axis molecules 
!     along axis-2.                                                      
   10 N=MARK2                                                            
      ER=ERMT                                                            
      MARK=N                                                             
      DO 19 I=MARK1,NEND1,NSTP1                                          
      K=-IT(I)                                                           
      IF(K.LE.0)GO TO 19                                                 
      L=IT(I+4)                                                          
      D=TI(I+1)+TI(I+5)**2                                               
      IF(D.GT.CN(L))GO TO 19                                             
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+6)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
   19 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      CALL MINHI(H2IH,0.,100.,H2N,NNE)                                   
      IF(KCT)904,20,906                                                  
   20 H2I=2.*H2IH                                                        
!-----Calculate the offset along axis-1 of a grid of A1 and A2 screw     
!     axis molecules.                                                    
      N=MARK                                                             
      ER=ERMF                                                            
      DO 39 I=1,NMOD                                                     
      DO 38 J=1,NMOD                                                     
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D3=H3N-W(3,J)-W(3,I)                                               
!-----Set D1 and D2 to collect distances for A1 molecules                
      D1=1.                                                              
      D2=H2IH-W(2,J)-W(2,I)                                              
   21 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 21                                            
   22 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 22                                            
      IF(D3.LT.-CN(K))GO TO 29                                           
   23 D2=D2+H2I                                                          
      IF(D2.LE.CN(K))GO TO 23                                            
   24 D2=D2-H2I                                                          
      IF(D2.GT.CN(K))GO TO 24                                            
      IF(D2.LT.-CN(K))GO TO 22                                           
      D=D2**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 24                                             
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=D1*W(1,J)-W(1,I)                                           
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 24                                                           
!-----Have distances for A2 molecules been collected?                    
   29 IF(D1.LT.0.)GO TO 38                                               
!     Set D1, D2, and D3 for A2 molecules.                               
      D1=-1.                                                             
      D2=W(2,J)-W(2,I)                                                   
      D3=D3-H3N                                                          
      GO TO 21                                                            
   38 CONTINUE                                                           
   39 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      CALL MINHI(H1IH,0.,100.,H1N,NNE)                                   
      IF(KCT)904,40,906                                                  
   40 H1I=2.*H1IH                                                        
      VX=H1IH*H2IH*H3M                                                   
      IF(V3.LT.VX)GO TO 400                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA                                                    
      IF(V3.GE.HT1(8))GO TO 400                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                            
      HT1(5)=0.                                                             
      HT1(6)=0.                                                             
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA                                                         
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=36                                                               
      MLOW=83                                                               
  400 IF(IPR.LT.2)GO TO 1000                                             
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,401)ALPHA,A1,A2,A3,VX,H1I,H2I,H3M                         
  401 FORMAT(1X,A2,3F7.1,'  V =',F8.2,'  AXES',3F8.3)                    
      DO 819 I=1,3                                                       
      IF(I-2)805,806,807                                                 
  805 SM(1)=H1IH                                                         
      SM(2)=H2IH                                                         
      SM(3)=H3N                                                          
      GO TO 812                                                          
  806 SM(1)=H1IH                                                         
      SM(2)=0.                                                           
      SM(3)=0.                                                           
      GO TO 812                                                          
  807 SM(1)=0.                                                           
      SM(2)=H2IH                                                         
      SM(3)=H3N                                                          
  812 SM(4)=SM(1)/H1I                                                    
      SM(5)=SM(2)/H2I                                                    
      SM(6)=SM(3)/H3M                                                    
      D=SQRT(SM(1)**2+SM(2)**2+SM(3)**2)                                 
      WRITE(NQ,815)SYM(I),SM,D                                           
  815 FORMAT(1X,A2,' AT',3X,3F8.3,';  ',3F8.4,'  D =',F6.2)              
  819 CONTINUE                                                           
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK3,NSTP3,MARK4,NSTP4,MARK,NSTP,MAXT                   
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTG,5X,3(I9,I3),I9)                  
      GO TO 999                                                             
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)                                                      
  905 FORMAT(20H0PARAMETER TOO SMALL)                                    
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT(1H0,'PARAMETER TOO LARGE')                                  
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)KCT,N,MARK3,NARK,MARK,NEND2                              
  921 FORMAT(16H0NO INTERACTIONS,I3,5I5)                                    
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END  SUBROUTINE ZSRTK  

