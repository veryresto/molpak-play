      SUBROUTINE ZORTA    
!                                                                        
!     For AP and AQ.                                                   
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
!      COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000), &
!     &      CODE(5000)
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
!      DATA ALPHA/'AP','AQ'/                                                 
!      DATA SYM/'A1','A2','A3'/  
!
      CHARACTER(2) :: ALPHA(2) = (/'AP','AQ'/)
      CHARACTER(2) :: SYM(3) = (/'A1','A2','A3'/)
!
      INTEGER :: I, ICT, IE, J, JCT, JE, K, KE, L, LCT
      INTEGER :: MCT, N, NZS
!
      REAL :: D, D1, D2, D3, DX2A1, DY1A2, DY3A1
      REAL :: H1A1, H1A2, H1A3, H1I,  HMIN, HL 
      REAL :: H2A1, H2A1R, H2A1S, H2A2, H2A3, H2I 
      REAL :: H3A1, H3A2, H3A3   
      REAL :: SM(6), VX, VY 
      REAL :: X1I, X1A1, X1A2, X2A1, X2A2, X3A1, X3A2
      REAL :: Y1A1, Y1A2, Y2A2, Y2A1, Y3A1, Z2A2
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)ITR(17),ITR(18)                                         
    5 FORMAT(1X,'ZORTA called',2I3)                                      
!-----TRIPLE AXIS CASE                                                      
!     Determine range of possible values of H2A1 from minimum positive   
!     and negative values with H1IH set at H1N                           
   10 KE=NE                                                                 
!-----Collect short interatomic distances along axis-2 and axis-3 for    
!     screw axis molecule along axis-1 (A1 molecule)                     
      N=MARK2                                                               
      NSTP2=5                                                            
      DO 89 I=1,NMOD                                                        
      DO 88 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D1=H1N+W(1,J)-W(1,I)                                                  
      IF(ABS(D1).GT.CN(K))GO TO 88                                          
      IT(N)=K                                                               
      TI(N+1)=D1**2                                                         
      TI(N+2)=-W(2,J)-W(2,I)                                                
      TI(N+3)=-W(3,J)-W(3,I)                                                
      IT(N+4)=L                                                             
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
   88 CONTINUE                                                              
   89 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
      ER=ERMT                                                               
!-----Set up for positive offset calculation                             
      HL=H1M                                                                
      X2A1=100.                                                             
      Y2A1=.75*H2M                                                          
      JCT=0                                                              
   92 DY3A1=H3M/NV                                                          
      Y3A1=-H3N                                                             
      MCT=0                                                                 
      JE=1                                                                  
   93 N=MARK                                                                
      DO 96 I=MARK2,NEND2,NSTP2                                             
      K=IT(I)                                                               
      D3=Y3A1+TI(I+3)                                                       
   94 D3=D3+H3M                                                             
      IF(D3.LE.CN(K))GO TO 94                                               
   95 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 95                                               
      IF(D3.LT.-CN(K))GO TO 96                                              
      D=TI(I+1)+D3**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 95                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 95                                                              
   96 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(Y2A1,0.,X2A1,HL,KE)                                        
      IF(KCT)904,151,152                                                 
  151 X3A1=Y3A1                                                          
      X2A1=Y2A1                                                          
      IF(MCT)156,153,156                                                 
  152 IF(MCT)156,153,155                                                 
  153 IF(JE.GE.NV)GO TO 154                                              
      JE=JE+1                                                            
      Y3A1=Y3A1+DY3A1                                                    
      GO TO 93                                                           
  154 MCT=1                                                              
      JE=2.*JE                                                           
      DY3A1=H3M/JE                                                       
      Y3A1=X3A1+DY3A1                                                    
      GO TO 93                                                           
  155 MCT=-1                                                             
      Y3A1=X3A1-DY3A1                                                    
      GO TO 93                                                           
  156 IF(JE.LT.KE)GO TO 154                                              
      IF(JE.GE.NNE)GO TO 158                                             
      KE=NNE                                                             
      GO TO 154                                                          
  158 CONTINUE                                                           
!-----Has negative offset been calculated?                               
      IF(JCT.NE.0)GO TO 160                                                 
      H2A1R=X2A1                                                            
!     Set up for negative offset calculation                             
      HL=-H1M                                                               
      X2A1=-100.                                                            
      Y2A1=-H2A1R                                                           
      JCT=1                                                              
      GO TO 92                                                              
!-----Set start of H2A1 range, H2A1S, and length of range, H2A1R         
  160 H2A1S=X2A1                                                            
      H2A1R=H2A1R-H2A1S                                                  
      IF(IPR.LT.3)GO TO 162                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,161)H2A1S,H2A1R                                              
  161 FORMAT(1X,7HH2A1S =,F8.3,9H, H2A1R =,F8.3)                            
  162 IF(ITR(17)-2)200,201,1000                                             
!-----One twofold axis - space group P2(1)2(1)2 - Code AP                
  200 LCT=1                                                                 
!     Set atom type multiplier to pick up A3 twofold axis distances      
      NZS=1                                                                 
      H3A3=0.                                                               
      GO TO 202                                                             
!-----Three screw axes - space group P2(1)2(1)2(1) - Code AQ             
  201 LCT=2                                                                 
!     Set atom type multiplier to pick up A3 screw axis distances        
      NZS=-1                                                                
      H3A3=H3N                                                              
!-----Set dummy value of cell volume per molecule                        
  202 VX=1000000.       ! 03-15-06                                                            
      KE=NE                                                                 
      ICT=0                                                                 
      IE=1                                                                  
      DX2A1=H2A1R/NV                                                        
      X2A1=H2A1S                                                            
!-----TOP OF A1 MOLECULE PLACEMENT LOOP                                     
  205 N=MARK2                                                               
      DO 209 I=1,NMOD                                                       
      DO 208 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=X2A1-W(2,J)-W(2,I)                                                 
      IF(ABS(D2).GT.CN(K))GO TO 208                                         
      IT(N)=K                                                               
      TI(N+1)=D2**2                                                         
      TI(N+2)=W(1,J)-W(1,I)                                                 
      TI(N+3)=-W(3,J)-W(3,I)                                                
      IT(N+4)=L                                                             
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
  208 CONTINUE                                                              
  209 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 260                                               
      NEND2=N-1                                                             
      MARK=N                                                                
      ER=ERMT                                                               
!-----Set up for calculation of A1 offset along axis-1                   
      X1A1=100.                                                             
      Y1A1=.75*H1M                                                          
      HMIN=H1N                                                           
      HL=H1M                                                             
      DY3A1=H3M/NV                                                          
      Y3A1=-H3N                                                             
      JE=1                                                                  
      MCT=0                                                                 
  216 N=MARK                                                                
      NARK=0                                                             
      DO 229 I=MARK2,NEND2,NSTP2                                            
      K=IT(I)                                                               
      D3=Y3A1+TI(I+3)                                                       
  221 D3=D3+H3M                                                             
      IF(D3.LT.CN(K))GO TO 221                                              
  222 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 222                                              
      IF(D3.LT.-CN(K))GO TO 229                                             
      D=TI(I+1)+D3**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 222                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 222                                                             
  229 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 260                                                
      NEND=N-1                                                              
      CALL MINHI(Y1A1,HMIN,X1A1,HL,KE)                                      
      IF(KCT)260,251,252                                                 
  251 X3A1=Y3A1                                                          
      X1A1=Y1A1                                                          
      IF(MCT)256,253,256                                                 
  252 IF(MCT)256,253,255                                                 
  253 IF(JE.GE.NV)GO TO 254                                              
      JE=JE+1                                                            
      Y3A1=Y3A1+DY3A1                                                    
      GO TO 216                                                          
  254 MCT=1                                                              
      JE=2.*JE                                                           
      DY3A1=H3M/JE                                                       
      Y3A1=X3A1+DY3A1                                                    
      GO TO 216                                                          
  255 MCT=-1                                                             
      Y3A1=X3A1-DY3A1                                                    
      GO TO 216                                                          
  256 IF(JE.LT.KE)GO TO 254                                              
      GO TO 265                                                          
  260 X1A1=H1N                                                              
      X3A1=Y3A1                                                             
!-----The vector from the origin to the A2 molecule is equal to the sum  
!     of those to the A1 and A3 molecules                                
  265 X3A2=X3A1+H3A3                                                     
!     The screw displacement is half of the cell length                  
      X1I=2.*X1A1                                                           
      IF(IPR.LT.3)GO TO 300                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,270)X1A1,X2A1,X3A1                                        
  270 FORMAT(1X,'X1A1 =',F8.3,',  X2A1 =',F8.3,'  X3A1 =',F8.3)          
!-----Find short A2 interatomic distances along axis-3                   
  300 N=MARK2                                                               
      DO 304 I=1,NMOD                                                       
      DO 303 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D3=X3A2-W(3,J)-W(3,I)                                                 
  301 D3=D3+H3M                                                             
      IF(D3.LE.CN(K))GO TO 301                                              
  302 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 302                                              
      IF(D3.LT.-CN(K))GO TO 303                                             
      IT(N)=K                                                               
      TI(N+1)=D3**2                                                         
      TI(N+2)=-W(1,J)-W(1,I)                                                
      TI(N+3)=W(2,J)-W(2,I)                                                 
      IT(N+4)=L                                                             
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 302                                                             
  303 CONTINUE                                                              
  304 CONTINUE                                                              
      NEND2=N-1                                                             
      MARK3=N                                                               
!-----Store A1 positions for double-standoff calculation by MINHI        
      MARK=MARK3                                                         
      N=MARK                                                             
      DO 308 I=1,NMOD                                                    
      DO 307 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D1=X1A1+W(1,J)-W(1,I)                                              
      IF(D1.GT.ABS(CN(K)))GO TO 307                                      
      D3=X3A1-W(3,J)-W(3,I)                                              
  305 D3=D3+H3M                                                          
      IF(D3.LT.CN(K))GO TO 305                                           
  306 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 306                                           
      IF(D3.LT.-CN(K))GO TO 307                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 306                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=X2A1-W(2,J)-W(2,I)                                         
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 306                                                          
  307 CONTINUE                                                           
  308 CONTINUE                                                           
      NARK=N                                                             
!-----TOP OF A2 MOLECULE PLACEMENT LOOP                                     
      DY1A2=X1I/NV                                                          
      Y1A2=-X1I/2.                                                          
      X2A2=100.                                                             
      Y2A2=1.5*H2N                                                          
      HL=H2M                                                                
      HMIN=H2N                                                           
      ER=ERMF                                                               
      JE=1                                                                  
      MCT=0                                                                 
  310 N=NARK                                                                
!-----PICK UP PREVIOUSLY STORED A3 MOLECULE DISTANCES, IF ANY            
      IF(NEND1.LE.MARK1)GO TO 314                                           
      DO 313 I=MARK1,NEND1,NSTP1                                            
!-----Pick up either twofold or screw axis distances                     
      K=NZS*IT(I)                                                           
      IF(K.LE.0)GO TO 313                                                   
!     The vector from an A2 to an A3 molecule is equal to that from the  
!     origin to an A1 molecule (X1A1, X2A1, X3A1)                        
      D1=Y1A2+X1A1+TI(I+5)                                                  
  311 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 311                                              
  312 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 312                                              
      IF(D1.LT.-CN(K))GO TO 313                                             
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 312                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+6)+X2A1                                                  
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 312                                                             
  313 CONTINUE                                                              
!-----ADD A2 MOLECULE DISTANCES, IF ANY                                     
  314 IF(NEND2.LE.MARK2)GO TO 318                                           
      DO 317 I=MARK2,NEND2,NSTP2                                            
      K=IT(I)                                                               
      D1=TI(I+2)+Y1A2                                                       
  315 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 315                                              
  316 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 316                                              
      IF(D1.LT.-CN(K))GO TO 317                                             
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 316                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 316                                                             
  317 CONTINUE                                                              
  318 IF(N.LE.MARK)GO TO 321                                                
      NEND=N-1                                                              
      CALL MINHI(Y2A2,HMIN,X2A2,HL,KE)                                      
      IF(KCT)320,325,342                                                    
!-----The minimum possible standoff is found in the positive axis-2      
!      direction - try negative direction                                
  320 Z2A2=-Y2A2                                                         
      CALL MINHI(Z2A2,-Y2A2,-X2A2,-HL,KE)                                
      IF(KCT)321,340,342                                                 
  321 X2A2=H2N                                                           
      X1A2=Y1A2                                                             
      GO TO 357                                                             
!-----Try negative axis-2 direction because of A1 and A3 molecules       
  325 Z2A2=-Y2A2                                                         
      CALL MINHI(Z2A2,-Y2A2,-X2A2,-HL,KE)                                
      IF(KCT)341,340,342                                                 
!-----Negative standoff is larger - change Y2A2                          
  340 Y2A2=-Z2A2                                                         
!-----This is the smallest value of Y1A2 so far                          
  341 X1A2=Y1A2                                                          
      X2A2=Y2A2                                                          
      IF(MCT)346,343,346                                                 
  342 IF(MCT)346,343,345                                                 
  343 IF(JE.GE.NV)GO TO 344                                              
      JE=JE+1                                                            
      Y1A2=Y1A2+DY1A2                                                    
      GO TO 310                                                          
  344 MCT=1                                                              
      JE=2.*JE                                                           
      DY1A2=X1I/JE                                                       
      Y1A2=X1A2+DY1A2                                                    
      GO TO 310                                                          
  345 MCT=-1                                                             
      Y1A2=X1A2-DY1A2                                                    
      GO TO 310                                                          
  346 IF(JE.LT.KE)GO TO 344                                              
  357 VY=X1A1*X2A2*H3M                                                      
      IF(IPR.LT.3)GO TO 359                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,358)IE,X1A1,X2A2,VY                            
  358 FORMAT(1X,'IE =',I5,', X1A1 =',F8.3,', X2A2 ='     &
     &,F8.3,', VY =',F8.2)                                               
  359 IF(VX-VY.GT..01)GO TO 361                                             
      IF(ICT)362,364,360                                                    
  360 ICT=-1                                                                
      X2A1=H2A1-DX2A1                                                       
      GO TO 205                                                             
  361 H1A1=X1A1                                                          
      H2A1=X2A1                                                          
      H3A1=X3A1                                                          
      H1A2=X1A2                                                          
      H2A2=X2A2                                                          
      VX=VY                                                                 
      IF(ICT)362,364,362                                                 
  362 IF(IE.GE.KE)GO TO 365                                                 
  363 IE=2*IE                                                               
      DX2A1=H2A1R/IE                                                        
      ICT=1                                                                 
      X2A1=H2A1+DX2A1                                                       
      GO TO 205                                                             
  364 IF(IE.GE.NV)GO TO 363                                                 
      IE=IE+1                                                               
      X2A1=X2A1+DX2A1                                                       
      GO TO 205                                                             
  365 IF(KE.GE.NNE)GO TO 366                                                
      KE=NNE                                                                
      GO TO 363                                                             
  366 H1I=2.*H1A1                                                        
      H2I=2.*H2A2                                                        
      H3A2=H3A1+H3A3                                                     
      H1A3=H1A1+H1A2                                                     
      H2A3=H2A1+H2A2                                                     
      IF(V3-VX.LT..1)GO TO 400                                              
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 400                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                            
      HT1(5)=H1A2                                                           
      HT1(6)=H2A3                                                           
      HT1(7)=H3A1                                                           
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=12                                                               
      MLOW=LCT+16                                                           
  400 IF(IPR.LT.2)GO TO 405                                              
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,401)ALPHA(LCT),A1,A2,A3,VX,H1I,H2I,H3M                    
  401 FORMAT(1X,A2,3F7.1,'  V =',F8.2,'  AXES',3F8.3)                    
      DO 819 I=1,3                                                       
      IF(I-2)805,806,807                                                 
  805 SM(1)=H1A1                                                         
      SM(2)=H2A1                                                         
      SM(3)=H3A1                                                         
      SM(4)=SM(1)                                                        
      SM(5)=SM(2)-HT1(6)                                                 
      SM(6)=SM(3)-HT1(7)                                                 
      GO TO 812                                                          
  806 SM(1)=H1A2                                                         
      SM(2)=H2A2                                                         
      SM(3)=H3A2                                                         
      SM(4)=SM(1)-HT1(5)                                                 
      SM(5)=SM(2)                                                        
      SM(6)=SM(3)-HT1(7)                                                 
      GO TO 812                                                          
  807 SM(1)=H1A3                                                         
      SM(2)=H2A3                                                         
      SM(3)=H3A3                                                         
      SM(4)=SM(1)-HT1(5)                                                 
      SM(5)=SM(2)-HT1(6)                                                 
      SM(6)=SM(3)                                                        
  812 SM(4)=SM(4)/H1I                                                    
      SM(5)=SM(5)/H2I                                                    
      SM(6)=SM(6)/H3M                                                    
      IF(SM(4).GT..5)SM(4)=SM(4)-1.                                      
      IF(SM(4).LE.-.5)SM(4)=SM(4)+1.                                     
      IF(SM(5).GT..5)SM(5)=SM(5)-1.                                      
      IF(SM(5).LE.-.5)SM(5)=SM(5)+1.                                     
      IF(SM(6).GT..5)SM(6)=SM(6)-1.                                      
      IF(SM(6).LE.-.5)SM(6)=SM(6)+1.                                     
      D=SQRT(SM(1)**2+SM(2)**2+SM(3)**2)                                 
      WRITE(NQ,815)SYM(I),SM,D                                           
  815 FORMAT(1X,A2,' AT',3X,3F8.3,';  ',3F8.4,'  D =',F6.2)              
  819 CONTINUE                                                           
      SM(1)=HT1(5)/2.                                                    
      SM(2)=HT1(6)/2.                                                    
      SM(3)=HT1(7)/2.                                                    
      SM(4)=SM(1)/H1I                                                    
      SM(5)=SM(2)/H2I                                                    
      SM(6)=SM(3)/H3M                                                    
      WRITE(NQ,820)SM                                                    
  820 FORMAT(1X,'ORIGIN',2X,3F8.3,';  ',3F8.4)                           
  405 IF(LCT.NE.1)GO TO 1000                                             
      IF(ITR(18)-1)1000,201,1000                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK3,NSTP3,MARK4,NSTP4,MARK,NSTP,MAXT         
  903 FORMAT(26H0STORAGE EXCEEDED BY ZORTA,5X,3(I9,I3),I9)      
      GO TO 999                                                             
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)                                                      
  905 FORMAT(20H0PARAMETER TOO SMALL)                                    
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)KCT,MCT,N,MARK3,NARK,MARK,NEND2                          
  921 FORMAT(16H0NO INTERACTIONS,2I3,5I5)                                   
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZORTA  
