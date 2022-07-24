      SUBROUTINE ZSRTJ                                                      
!
!     For ST and SU.
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
!      DIMENSION SM(6)                                                    
!      CHARACTER*2 ALPHA(2)                                               
!      CHARACTER*3 SYM(7)                                                 
!      DATA ALPHA/'ST','SU'/                                              
!      DATA SYM/'T23','T13','T12','D00','D23','D13','D12'/                
!
      CHARACTER(2) :: ALPHA(2) = (/'ST','SU'/)
      CHARACTER(3) :: SYM(7) = (/'T23','T13','T12','D00','D23', &
     &                           'D13','D12'/)
!
      INTEGER :: I, ICT, J, K, L, LCT, N, NERR
!
      REAL :: D, D1, D2, D3
      REAL :: H1I, H1IH, H1IQ
      REAL :: H2I, H2IH, H2IQ
      REAL :: H3I, H3IH, H3IQ
      REAL :: SM(6)
      REAL :: VX, VH
      REAL :: X1I, X1IH, X2IQ, X3IH
      REAL :: Y2IQ
!
      NARK=0                                                             
      VH=1000000.                                                        
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)ITR(81),ITR(82)                                         
    5 FORMAT(1X,'ZSRTJ called',2I3)                                      
   10 IF(ITR(81)-2)100,200,1000                                          
!-----Calculate area-1,3 with I molecules in contact along axis-3        
  100 LCT=1                                                              
      H3I=H3M                                                            
      H3IH=H3N                                                           
      H3IQ=H3M/4.                                                        
!-----Determine H1IH and H1I by placing a row of T13 molecules parallel  
!     to axis-3 - H2T13=0 - H1I=2*H1T13                                  
!-----Collect I molecules for double standoff                            
      MARK=MARK2                                                         
      N=MARK                                                             
      DO 105 I=MARK1,NEND1,NSTP1                                         
      K=IT(I)                                                            
      IF(K.LE.0)GO TO 105                                                
      D2=TI(I+3)                                                         
      IF(ABS(D2).GT.CN(K))GO TO 105                                      
      D=TI(I+1)+D2**2                                                    
      L=IT(I+4)                                                          
      IF(D.GT.CN(L))GO TO 105                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 891                                               
  105 CONTINUE                                                              
      NARK=N                                                             
      DO 109 I=MARK1,NEND1,NSTP1                                         
      K=-IT(I)                                                           
      IF(K.LE.0)GO TO 109                                                
      D2=TI(I+3)                                                         
      IF(ABS(D2).GT.CN(K))GO TO 109                                      
      D=TI(I+1)+D2**2                                                    
      L=IT(I+4)                                                          
      IF(D.GT.CN(L))GO TO 109                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 892                                               
  109 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERMT                                                               
      NEND=N-1                                                              
      H1IH=H2N                                                              
      CALL MINHI(H1IH,0.,100.,H1M,NNE)                                      
      IF(KCT)904,110,906                                                 
  110 H1I=2.*H1IH                                                        
      H1IQ=H1IH/2.                                                       
      IF(IPR.LT.3)GO TO 120                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,111)H1IH,H1I                                              
  111 FORMAT(1X,'H1IH =',F7.3,',  H1I =',F7.3)                           
  120 IF(ITR(82).NE.1)GO TO 300                                          
!-----Calculate area-1,3 with molecules in contact along axis-1          
  200 NARK=0                                                             
      MARK=MARK2                                                         
      N=MARK                                                             
      DO 205 I=MARK1,NEND1,NSTP1                                         
      IF(IT(I+7).NE.0)GO TO 205                                          
      K=IT(I)                                                            
      IF(K.LE.0)GO TO 205                                                
      D2=TI(I+3)                                                         
      IF(ABS(D2).GT.CN(K))GO TO 205                                      
      D=TI(I+1)+D2**2                                                    
      L=IT(I+4)                                                          
      IF(D.GT.CN(L))GO TO 205                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 893                                               
  205 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      X1I=H1M                                                            
      CALL MINHI(X1I,0.,100.,H1M,NNE)                                    
      IF(KCT)904,206,906                                                 
!-----Determine H3IH and H3I by placing a row of T13 molecules parallel  
!     to axis-1 - H2T13=0 - H3I=2*H3T13                                  
!-----Collect I molecules for double standoff                            
  206 X1IH=X1I/2.                                                        
      N=MARK                                                             
      DO 210 I=1,NMOD                                                    
      DO 209 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=W(2,J)-W(2,I)                                                   
      IF(D2.GT.ABS(CN(K)))GO TO 209                                      
      D1=W(1,J)-W(1,I)                                                   
  207 D1=D1+X1I                                                          
      IF(D1.LE.CN(K))GO TO 207                                           
  208 D1=D1-X1I                                                          
      IF(D1.LT.-CN(K))GO TO 209                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L))GO TO 208                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(3,J)-W(3,I)                                              
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 894                                               
      GO TO 208                                                          
  209 CONTINUE                                                             
  210 CONTINUE                                                           
      NARK=N                                                             
      DO 219 I=1,NMOD                                                    
      DO 218 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=W(2,J)-W(2,I)                                                   
      IF(D2.GT.ABS(CN(K)))GO TO 218                                      
      D1=X1IH+W(1,J)-W(1,I)                                              
  215 D1=D1+X1I                                                          
      IF(D1.LE.CN(K))GO TO 215                                           
  216 D1=D1-X1I                                                          
      IF(D1.LT.-CN(K))GO TO 218                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L))GO TO 216                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(3,J)-W(3,I)                                              
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 895                                               
      GO TO 216                                                          
  218 CONTINUE                                                             
  219 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERMT                                                            
      X3IH=H3N                                                           
      CALL MINHI(X3IH,0.,100.,H3M,NNE)                                   
      IF(KCT)904,221,906                                                 
  221 IF(IPR.LT.3)GO TO 230                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,225)X1I,X3IH                                              
  225 FORMAT(' X1I =',F7.3,', X3IH =',F7.3)                              
  230 IF(ITR(81).NE.1)GO TO 240                                          
      IF(H1I*H3IH.LE.X1I*X3IH)GO TO 300                                  
  240 LCT=2                                                              
      H1I=X1I                                                            
      H1IH=X1I/2.                                                        
      H1IQ=X1I/4.                                                        
      H3IH=X3IH                                                          
      H3I=2.*X3IH                                                        
      H3IQ=X3IH/2.                                                       
!-----Determine H2IQ and H2I by placing a plane of D00 and D13           
!     molecules above and D23 and D12 molecules below.                   
!-----Collect T23 and T12 molecules for axis-2 double standoff.          
  300 N=MARK                                                             
      DO 309 I=1,NMOD                                                    
      DO 308 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
!-----Set flag for T23 molecules                                         
      ICT=1                                                              
      D1=W(1,J)-W(1,I)                                                   
      D3=H3IH+W(3,J)-W(3,I)                                              
  301 D3=D3+H3I                                                          
      IF(D3.LT.CN(K))GO TO 301                                           
  302 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 305                                          
  303 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 303                                           
  304 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 302                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 304                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 896                                            
      GO TO 304                                                          
!-----Have T12 molecules been added?                                     
  305 IF(ICT.LT.0)GO TO 308                                              
      ICT=-1                                                             
      D1=D1+H1IH                                                         
      D3=D3-H3IH                                                         
      GO TO 301                                                          
  308 CONTINUE                                                           
  309 CONTINUE                                                           
      NARK=N                                                             
      DO 329 I=1,NMOD                                                    
      DO 328 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
!-----Set flag for D00 molecules                                         
      ICT=1                                                              
      D1=H1IQ-W(1,J)-W(1,I)                                              
      D2=W(2,J)-W(2,I)                                                   
      D3=H3IQ+W(3,J)-W(3,I)                                              
  311 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 311                                           
  312 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 320                                          
  313 D3=D3+H3I                                                          
      IF(D3.LT.CN(K))GO TO 313                                           
  314 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 312                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 314                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=D2                                                         
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 897                                            
      GO TO 314                                                          
!-----Have D13 molecule positions been added?                            
  320 IF(ICT.LT.0)GO TO 328                                              
      D1=D1+H1IH                                                         
      D3=D3+H3IH                                                         
      ICT=-1                                                             
      GO TO 311                                                          
  328 CONTINUE                                                           
  329 CONTINUE                                                           
      IF(N.LE.NARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERMF                                                            
      X2IQ=H2N/2.                                                        
      CALL MINHI(X2IQ,0.,100.,H2N,NNE)                                   
      IF(KCT)904,330,906                                                 
!-----Approach of D23 and D12 molecules from below.                      
  330 N=NARK                                                             
      DO 349 I=1,NMOD                                                    
      DO 348 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
!-----Set flag for D12 molecules                                         
      ICT=1                                                              
      D1=H1IH+H1IQ-W(1,J)-W(1,I)                                             
      D2=W(2,J)-W(2,I)                                                   
      D3=H3IQ+W(3,J)-W(3,I)                                              
  331 D1=D1+H1I                                                          
      IF(D1.LT.CN(K))GO TO 331                                           
  332 D1=D1-H1I                                                          
      IF(D1.LT.-CN(K))GO TO 340                                          
  333 D3=D3+H3I                                                          
      IF(D3.LT.CN(K))GO TO 333                                           
  334 D3=D3-H3I                                                          
      IF(D3.LT.-CN(K))GO TO 332                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 334                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=D2                                                         
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 898                                            
      GO TO 334                                                          
!-----Have D23 molecule positions been added?                            
  340 IF(ICT.LT.0)GO TO 348                                              
      ICT=-1                                                             
      D1=D1-H1IH                                                         
      D3=D3+H3IH                                                         
      GO TO 331                                                          
  348 CONTINUE                                                           
  349 CONTINUE                                                           
      IF(N.LE.NARK)GO TO 920                                             
      NEND=N-1                                                           
      Y2IQ=-H2N/2.                                                       
      CALL MINHI(Y2IQ,-X2IQ,-100.,-H2N,NNE)                                  
      IF(KCT)351,350,906                                                 
  350 H2IQ=-Y2IQ                                                         
      GO TO 352                                                          
  351 H2IQ=X2IQ                                                          
  352 H2IH=2.*H2IQ                                                       
      H2I=4.*H2IQ                                                        
      IF(IPR.LT.3)GO TO 400                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,355)X2IQ,Y2IQ                                             
  355 FORMAT(' X2IQ =',F7.3,', Y2IQ =',F7.3)                             
  400 VX=H1I*H2I*H3I/8.                                                     
      IF(VX.GE.V3)GO TO 500                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 500                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
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
      NLOW=35                                                               
      MLOW=81                                                            
  500 IF(IPR.LT.2)GO TO 1000                                             
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,525)ALPHA(LCT),A1,A2,A3,V3,H1I,H2I,H3I                    
  525 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3)                    
  526 FORMAT(1X,A3,' AT',3F7.3,'  D =',F7.3)                             
      D=SQRT(H2IH**2+H3IH**2)                                            
      WRITE(NQ,526)SYM(1),0.,H2IH,H3IH,D                                 
      D=SQRT(H1IH**2+H3IH**2)                                            
      WRITE(NQ,526)SYM(2),H1IH,0.,H3IH,D                                 
      D=SQRT(H1IH**2+H2IH**2)                                            
      WRITE(NQ,526)SYM(3),H1IH,H2IH,0.,D                                 
      D=SQRT(H1IQ**2+H2IQ**2+H3IQ**2)                                    
      WRITE(NQ,526)SYM(4),H1IQ,H2IQ,H3IQ,D                               
      GO TO 1000                                                         
  891 NERR=1                                                             
      GO TO 902                                                          
  892 NERR=2                                                             
      GO TO 902                                                          
  893 NERR=3                                                             
      GO TO 902                                                          
  894 NERR=4                                                             
      GO TO 902                                                          
  895 NERR=5                                                             
      GO TO 902                                                          
  896 NERR=6                                                             
      GO TO 902                                                          
  897 NERR=7                                                             
      GO TO 902                                                          
  898 NERR=8                                                             
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)NERR,MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT           
  903 FORMAT(26H0STORAGE EXCEEDED BY ZSRTJ,2X,I2,3(I9,I3),I9)               
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
      WRITE(NQ,921)KCT,MARK1,MARK2,NEND2,MARK                               
  921 FORMAT(16H0NO INTERACTIONS,I3,4I9)                                    
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZSRTJ  
