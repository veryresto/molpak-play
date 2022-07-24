      SUBROUTINE ZNRTA     
!                                                                        
!     For DA, DB, DC, DD, DE                                                   
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
      CHARACTER(2) :: ALPHA(5) = (/'DA','DB','DC','DD','DE'/)
      CHARACTER(3) :: SYM(7) =(/'P3 ','A3 ','C  ','IC ','P3C', &
     & 'A3C','CC '/) 
!      DIMENSION ALPHA(5),A(3),C(6),SM(6),ORIG(3),SYM(7)                  
!      DATA ALPHA/'DA','DB','DC','DD','DE'/                               
!      DATA SYM/'P3 ','A3 ','C  ','IC ','P3C','A3C','CC '/                
!-----This subroutine treats C-centered monoclinic space groups          
! 
      INTEGER :: I, ICT, IE, J, JCT, JE, K, KE, L, LE, LCT
      INTEGER :: M, MCT, N, NCT, NAREA, NERR, NSRCH, NVT 
!
      REAL :: A(3), AR, AREA
      REAL :: C(6), D, D1, D2, D3, D1AJ, D2AJ, D2IJ
      REAL :: DX1IC, DX1P3, DX3IC, DX3P3, DY1A3, DY1C, DY3C 
      REAL :: EGI, EGY
      REAL :: H1A3, H1C, H1I, H1IC, H1ICF, H1ICFF, H1ICR, H1ICS,  H1P3 
      REAL :: H1A3N, H1A3P, H2A3, H2A3N, H2A3P, H2C, H21I, H2I, H2P3, HL
      REAL :: H3A3, H3C, H3I, H3IC, H3ICF, H3ICR, H3P3, HH1P3, HH3P3 
      REAL :: ORIG(3), SINE3, SM(6)
      REAL :: VX, VY
      REAL :: X1A3, X1C, X1CT, X1I, X1IC, X1P3 
      REAL :: X2A3, X2A3P, X2C, X2I, X2CT, X21I, X2P3
      REAL :: X3C, X3IC, X3I, X3P3
      REAL :: XX1C, XX3C
      REAL :: Y1A3, Y1A3P, Y1A3T, Y1AMC, Y1C, Y1TAMC 
      REAL :: Y2A3, Y2A3P, Y2C, Y3C
      REAL :: Z2A3

      NARK=0                                                             
      NVT=2*NV                                                           
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=44,48)                                        
    5 FORMAT(1X,'ZNRTA called',5I3)                                      
!     Find minimum value of H1IC - molecules in contact along axis-1     
   10 N=MARK2                                                               
      NSTP2=4                                                            
      DO 15 I=1,NMOD                                                        
      DO 14 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      D2=W(2,J)-W(2,I)                                                      
      IF(ABS(D2).GT.CN(K))GO TO 14                                          
      IT(N)=K                                                            
      TI(N+1)=D2**2                                                      
      TI(N+2)=W(1,J)-W(1,I)                                              
      TI(N+3)=W(3,J)-W(3,I)                                              
      N=N+NSTP2                                                          
      IF(N.GT.LIMIT)GO TO 902                                            
   14 CONTINUE                                                           
   15 CONTINUE                                                           
      IF(N.LE.MARK2)GO TO 930                                            
      NEND2=N-1                                                          
      MARK3=N                                                            
      MARK=N                                                             
      DO 19 I=MARK2,NEND2,NSTP2                                          
      K=IT(I)                                                            
      IF(ABS(TI(I+3)).GT.CN(K))GO TO 19                                  
      D=TI(I+1)+TI(I+3)**2                                               
      IF(D.GT.CN(K+10))GO TO 19                                          
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
   19 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 930                                             
      NEND=N-1                                                           
      X1I=H1M                                                            
      CALL MINHI(X1I,0.,100.,H1M,NNE)                                    
      IF(KCT)911,20,921                                                  
   20 H1ICS=.5*X1I                                                       
!-----Find value of X1IC with molecular contact along axis-3             
   40 MARK=MARK3                                                         
      N=MARK                                                             
      DO 44 I=MARK2,NEND2,NSTP2                                          
      K=IT(I)                                                            
      D3=H3N+TI(I+3)                                                     
   41 D3=D3+H3M                                                          
      IF(D3.LE.CN(K))GO TO 41                                            
   42 D3=D3-H3M                                                          
      IF(D3.GT.CN(K))GO TO 42                                            
      IF(D3.LT.-CN(K))GO TO 44                                           
      D=TI(I+1)+D3**2                                                    
      IF(D.GT.CN(K+10))GO TO 42                                          
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 42                                                           
   44 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 930                                             
      NEND=N-1                                                           
      H1ICF=H1M                                                          
      CALL MINHI(H1ICF,H1ICS,100.,H1M,NNE)                               
      IF(KCT)50,45,922                                                   
   45 NAREA=0                                                            
      H1ICR=H1ICF-H1ICS                                                  
      IF(IPR.LT.3)GO TO 100                                              
      WRITE(NQ,46)H1ICS,H1ICF                                            
   46 FORMAT(1X,'H1ICS =',F7.3,', H1ICF =',F7.3)                         
      GO TO 100                                                          
!-----There is no range of possible X1IC values                          
   50 NAREA=-1                                                           
      H1IC=H1ICS                                                         
      H1I=2.*H1IC                                                        
      H3IC=H3N                                                           
      H3I=H3M                                                            
      IF(IPR.LT.3)GO TO 100                                              
      WRITE(NQ,55)H1IC,H1I                                               
   55 FORMAT(' H1IC =',F8.3,'   H1I =',F8.3)                             
  100 LCT=0                                                              
!     NAREA=0                                                            
!-----Transfer to next operation based on space group/sub-group          
  101 LCT=LCT+1                                                          
      GO TO(105,110,115,120,125,1000),LCT                                
  105 GO TO(106,111,116,121,126,1000),ITR(44)                            
!-----Space group Cc - Code DA                                           
  106 LCT=1                                                              
      GO TO 160                                                          
  110 GO TO(111,116,121,126,1000),ITR(45)                                
!-----Space group C2 - Code DB                                           
  111 LCT=2                                                              
      IF(NAREA)300,160,300                                               
  115 GO TO(116,121,126,1000),ITR(46)                                    
!-----Space group C2/c - Subgroup 1 - Code DC                            
  116 LCT=3                                                              
      X3I=H3M                                                            
      X3IC=H3N                                                           
      X1IC=H1ICF                                                         
      X1I=2.*X1IC                                                        
      KE=NNE                                                             
      GO TO 450                                                          
  120 GO TO(121,126,1000),ITR(47)                                        
!-----Space group C2/c - Subgroup 2 - Code DD                            
  121 LCT=4                                                              
      IF(NAREA)425,160,425                                               
  125 IF(ITR(48).GT.1)GO TO 1000                                         
!-----Space group C2/c - Subgroup 3 - Code DE                            
  126 LCT=5                                                              
      GO TO 400                                                          
!-----Find value of H1IC, between H1ICS and H1ICF, which gives min. AREA 
  160 NAREA=1                                                            
      N=MARK3                                                            
      MARK=N                                                             
      KE=NE                                                                 
      JE=0                                                                  
      MCT=0                                                              
      X3I=H3M                                                               
      X3IC=H3N                                                           
      X1IC=H1ICS                                                            
      X1I=2.*X1IC                                                        
      DX1IC=H1ICR/NV                                                     
      ER=ERMT                                                               
!-----Set dummy starting face half area                                  
      AREA=10000.                                                           
  170 N=MARK                                                             
!-----Collect cell corner distances for double standoff                  
      DO 173 I=MARK2,NEND2,NSTP2                                         
      K=IT(I)                                                            
      L=K+10                                                             
      D1=TI(I+2)                                                         
  171 D1=D1+X1I                                                          
      IF(D1.LE.CN(K))GO TO 171                                           
  172 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 172                                           
      IF(D1.LT.-CN(K))GO TO 173                                          
      D=TI(I+1)+D1**2                                                    
      IF(D.GT.CN(L))GO TO 172                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+3)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 172                                                          
  173 CONTINUE                                                           
      NARK=N                                                             
!-----Add face center I molecules                                        
      DO 176 I=MARK2,NEND2,NSTP2                                         
      K=IT(I)                                                            
      L=K+10                                                             
      D1=X1IC+TI(I+2)                                                       
  174 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 174                                              
  175 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 175                                              
      IF(D1.LT.-CN(K))GO TO 176                                             
      D=TI(I+1)+D1**2                                                       
      IF(D.GT.CN(L))GO TO 175                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 175                                                             
  176 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 930                                                
      NEND=N-1                                                              
      CALL MINHI(X3IC,0.,100.,H3M,KE)                                       
      AR=X1I*X3IC                                                        
      IF(IPR.LT.3)GO TO 178                                              
      WRITE(NQ,177)X1IC,X3IC,AR                                          
  177 FORMAT(1X,'X1IC =',F7.3,', X3IC =',F7.3,', AR =',F8.3)             
!-----Is this the smallest face half area found so far?                  
  178 IF(AR.GE.AREA)GO TO 191                                            
      DX3IC=X3IC/NNE                                                     
  180 X3I=2.*X3IC                                                        
!-----Check center molecules in surrounding cells                        
      EGY=0.                                                             
      DO 185 I=MARK2,NEND2,NSTP2                                         
      K=IT(I)                                                            
      L=K+10                                                             
      D1=X1IC+TI(I+2)                                                    
      D3=X3IC+TI(I+3)                                                    
  181 D3=D3+X3I                                                          
      IF(D3.LE.CN(K))GO TO 181                                           
  182 D3=D3-X3I                                                          
      IF(D3.GT.CN(K))GO TO 182                                           
      IF(D3.LT.-CN(K))GO TO 185                                          
  183 D1=D1+X1I                                                          
      IF(D1.LE.CN(K))GO TO 183                                           
  184 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 184                                           
      IF(D1.LT.-CN(K))GO TO 182                                          
      D=TI(I+1)+D1**2+D3**2                                              
      IF(D.GT.CN(L))GO TO 184                                            
      M=IFIX((CN(L)-D)/CN(K+11))+1                                      
      IF(M.LE.9)EGI=CN(K+M)                               
      IF(M.GT.9)EGI=CN(K+9)+(CN(K+9)-CN(K+8))*(M-9)         
      EGY=EGY+EGI                                          
      IF(EGY.GT.ERMT)GO TO 186                                           
      GO TO 184                                                          
  185 CONTINUE                                                           
      GO TO 190                                                          
  186 X3IC=X3IC+DX3IC                                                    
      AR=X3IC*X1I                                                        
      IF(AR.GE.AREA)GO TO 191                                            
  190 AREA=AR                                                            
      H1IC=X1IC                                                          
      H3IC=X3IC                                                          
      H1I=2.*H1IC                                                        
      H3I=2.*H3IC                                                        
      IF(MCT)193,192,193                                                 
  191 IF(MCT)193,192,195                                                 
  192 IF(JE.EQ.0)H3ICF=X3IC                                              
      IF(LCT.EQ.4)GO TO 425                                              
      IF(JE.GE.NV)GO TO 194                                              
      JE=JE+1                                                            
      X1IC=X1IC+DX1IC                                                    
      X1I=2.*X1IC                                                        
      GO TO 170                                                          
  193 IF(JE.GE.KE)GO TO 196                                              
  194 JE=2*JE                                                            
      DX1IC=H1ICR/JE                                                     
      X1IC=H1IC+DX1IC                                                    
      IF(X1IC.GT.H1ICF)GO TO 195                                         
      X1I=2.*X1IC                                                        
      MCT=1                                                              
      GO TO 170                                                          
  195 X1IC=H1IC-DX1IC                                                    
      IF(X1IC.LT.H1ICS)GO TO 193                                         
      X1I=2.*X1IC                                                        
      MCT=-1                                                             
      GO TO 170                                                          
  196 IF(JE.GE.NNE)GO TO 197                                             
      KE=NNE                                                             
      GO TO 194                                                          
  197 IF(IPR.LT.3)GO TO 200                                                 
      WRITE(NQ,198)A3,H1I,H3I                                               
  198 FORMAT(1X,4HA3 =,F5.1,', H1I =',F7.3,', H3I =',F7.3)                  
  200 IF(LCT.GE.2)GO TO 300                                              
!-----P3 molecules                                                       
      H2P3=100.                                                             
      X2P3=H2N                                                              
      HL=H2M                                                                
      ER=ERMF                                                               
      JE=1                                                                  
      LE=1                                                                  
      KE=NE                                                                 
      DX1P3=H1I/NV                                                          
      DX3P3=H3I/NV                                                           
      X1P3=-H1IC                                                            
      X3P3=-H3IC                                                            
      MCT=1                                                                 
      NCT=1                                                              
  205 N=MARK3                                                               
      X21I=2.*X1P3                                                       
!-----Collect I molecules for double standoff                            
      DO 219 I=1,NMOD                                                    
      DO 218 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      JCT=0                                                              
      D1=X21I+W(1,J)-W(1,I)                                              
      D3=W(3,J)-W(3,I)                                                   
  211 D1=D1+H1I                                                          
      IF(D1.LE.CN(K))GO TO 211                                           
  212 D1=D1-H1I                                                          
      IF(D1.GT.CN(K))GO TO 212                                           
      IF(D1.LT.-CN(K))GO TO 215                                          
  213 D3=D3+H3I                                                          
      IF(D3.LE.CN(K))GO TO 213                                           
  214 D3=D3-H3I                                                          
      IF(D3.GT.CN(K))GO TO 214                                           
      IF(D3.LT.-CN(K))GO TO 212                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 214                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 214                                                          
  215 IF(JCT.NE.0)GO TO 218                                              
      JCT=1                                                              
      D1=D1+H1IC                                                         
      D3=D3+H3IC                                                         
      GO TO 211                                                          
  218 CONTINUE                                                           
  219 CONTINUE                                                           
      NARK=N                                                             
      DO 229 I=1,NMOD                                                       
      DO 228 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                             
      JCT=0                                                              
      D1=X1P3+W(1,J)-W(1,I)                                                 
      D3=X3P3-W(3,J)-W(3,I)                                                 
  221 D1=D1+H1I                                                              
      IF(D1.LE.CN(K))GO TO 221                                              
  222 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 222                                              
      IF(D1.LT.-CN(K))GO TO 225                                             
  223 D3=D3+H3I                                                          
      IF(D3.LE.CN(K))GO TO 223                                           
  224 D3=D3-H3I                                                          
      IF(D3.GT.CN(K))GO TO 224                                           
      IF(D3.LT.-CN(K))GO TO 222                                          
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 224                                            
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=W(2,J)-W(2,I)                                                 
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 224                                                             
  225 IF(JCT.NE.0)GO TO 228                                              
      JCT=1                                                              
      D1=D1+H1IC                                                         
      D3=D3+H3IC                                                         
      GO TO 221                                                          
  228 CONTINUE                                                              
  229 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 930                                                
      NEND=N-1                                                              
      CALL MINHI(X2P3,0.,H2P3,HL,KE)                                       
      IF(KCT)913,242,241                                                 
  241 IF(MCT.LE.1)GO TO 244                                              
      J=MCT-NCT                                                          
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 257                                 
      GO TO 243                                                          
  242 H1P3=X1P3                                                          
      H3P3=X3P3                                                          
      H2P3=X2P3                                                          
      NCT=MCT                                                            
  243 GO TO(244,249,250,251,252,253,254,255,248),MCT                     
  244 IF(LE.GE.NVT)GO TO 245                                             
      LE=LE+1                                                            
      X3P3=X3P3+DX3P3                                                    
      GO TO 205                                                          
  245 IF(JE.GE.NV)GO TO 246                                              
      JE=JE+1                                                            
      X1P3=X1P3+DX1P3                                                    
      LE=1                                                               
      X3P3=0.                                                            
      GO TO 205                                                          
  246 NCT=6                                                              
  247 JE=2.*JE                                                           
      DX1P3=H1I/JE                                                       
      DX3P3=H3I/JE                                                       
      HH1P3=H1P3                                                         
      HH3P3=H3P3                                                         
      NCT=NCT+1                                                          
      IF(NCT.GT.9)NCT=NCT-8                                              
      GO TO(923,250,251,252,253,254,255,248,249),NCT                     
  248 MCT=2                                                              
      X1P3=HH1P3+DX1P3                                                   
      X3P3=HH3P3                                                         
      GO TO 205                                                          
  249 MCT=3                                                              
      X1P3=HH1P3+DX1P3                                                   
      X3P3=HH3P3+DX3P3                                                   
      GO TO 205                                                          
  250 MCT=4                                                              
      X1P3=HH1P3                                                         
      X3P3=HH3P3+DX3P3                                                   
      GO TO 205                                                          
  251 MCT=5                                                              
      X1P3=HH1P3-DX1P3                                                   
      X3P3=HH3P3+DX3P3                                                   
      GO TO 205                                                          
  252 MCT=6                                                              
      X1P3=HH1P3-DX1P3                                                   
      X3P3=HH3P3                                                         
      GO TO 205                                                          
  253 MCT=7                                                              
      X1P3=HH1P3-DX1P3                                                   
      X3P3=HH3P3-DX3P3                                                   
      GO TO 205                                                          
  254 MCT=8                                                              
      X1P3=HH1P3                                                         
      X3P3=HH3P3-DX3P3                                                   
      GO TO 205                                                          
  255 MCT=9                                                              
      X1P3=HH1P3+DX1P3                                                   
      X3P3=HH3P3-DX3P3                                                   
      GO TO 205                                                          
  257 IF(JE.LT.KE)GO TO 247                                              
      IF(JE.GE.NNE)GO TO 258                                             
      KE=NNE                                                             
      GO TO 247                                                          
  258 CONTINUE                                                           
      VX=H1IC*H2P3*H3I                                                      
      H2I=2.*H2P3                                                           
      H21I=2.*H1P3                                                       
      IF(IPR.LT.2)GO TO 270                                                 
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
      SINE3=SQRT(1.-C(3)**2)                                             
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      ORIG(1)=0.                                                         
      ORIG(2)=0.                                                         
      ORIG(3)=H3P3                                                       
      WRITE(NQ,260)ALPHA(LCT),A1,A2,A3,VX,A,C                            
  260 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3  &                  
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
  261 FORMAT(1X,A3,' AT',3F7.3,5X,3F7.4,', D =',F7.2)                    
      DO 265 I=1,3                                                       
      GO TO(262,263,264),I                                               
  262 SM(3)=H3P3                                                         
      SM(2)=H2P3/SINE3                                                   
      SM(1)=H1P3-SM(2)*C(3)                                              
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=SM(2)/A(2)                                                   
      SM(6)=(SM(3)-ORIG(3))/A(3)                                         
      D=SQRT(H1P3**2+H2P3**2+H3P3**2)                                    
      WRITE(NQ,261)SYM(1),SM,D                                           
      GO TO 265                                                          
  263 SM(3)=H3IC                                                         
      SM(2)=0.                                                           
      SM(1)=H1IC                                                         
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=SM(2)/A(2)                                                   
      SM(6)=SM(3)/A(3)                                                   
      D=SQRT(H1IC**2+H3IC**2)                                            
      WRITE(NQ,261)SYM(4),SM,D                                           
      GO TO 265                                                          
  264 D1=H1P3+H1IC                                                       
      D2=H2P3                                                            
      SM(3)=H3P3+H3IC                                                    
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=SM(2)/A(2)                                                   
      SM(6)=(SM(3)-ORIG(3))/A(3)                                         
      D=SQRT(D1**2+D2**2+SM(3)**2)                                       
      WRITE(NQ,261)SYM(5),SM,D                                           
  265 CONTINUE                                                           
  270 IF(VX.GE.V3)GO TO 101                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 101                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=0.                                                             
      HT1(6)=0.                                                             
      HT1(7)=H3P3                                                           
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=16                                                               
      MLOW=44                                                            
      GO TO 101                                                             
!-----Space group C2 - Code DB - Locate A3 molecules                        
  300 N=MARK2                                                            
      MARK=N                                                                
      NARK=0                                                             
      JCT=0                                                                 
      X2A3=100.                                                             
      Y2A3=1.5*H2N                                                          
      HL=H2M                                                                
      ER=ERMF                                                               
      GO TO 304                                                             
  301 X2A3=-100.                                                            
      Y2A3=-1.5*H2N                                                         
      HL=-H2M                                                               
      D1AJ=-H1A3P
	D2AJ=-H2A3P
      GO TO 304
  302 X2A3=100.
      Y2A3=Y2A3P
	HL=H2M
	D1AJ=-H1A3N
	D2AJ=-H2A3N
  304 JE=1                                                                  
      KE=NE                                                                 
      DY1A3=H1I/NV                                                          
      Y1A3=-H1IC                                                            
      MCT=0                                                                 
  305 N=MARK                                                                
      NARK=0                                                             
!-----Collect second layer of A3 molecules on second and third passes    
      IF(JCT.EQ.0)GO TO 316                                              
      Y1A3T=2.*Y1A3      
      DO 312 I=1,NMOD                                                    
      DO 311 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      ICT=0                                                                 
      D1=Y1A3T+D1AJ-W(1,J)-W(1,I)                                           
      D2=D2AJ-W(2,J)-W(2,I)
      D3=W(3,J)-W(3,I)                                                   
  306 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 306                                              
  307 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 307                                              
      IF(D3.LT.-CN(K))GO TO 310                                             
  308 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 308                                              
  309 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 309                                              
      IF(D1.LT.-CN(K))GO TO 307                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 309                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                         
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 309                                                             
  310 IF(ICT.NE.0)GO TO 311                                                 
      ICT=1                                                                 
      D1=D1+H1IC                                                            
      D3=D3+H3IC                                                            
      GO TO 306                                                             
  311 CONTINUE                                                              
  312 CONTINUE                                                           
      NARK=N                                                             
  316 DO 331 I=1,NMOD                                                    
      DO 330 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      ICT=1                                                                 
      D1=Y1A3-W(1,J)-W(1,I)                                                 
      D2=-W(2,J)-W(2,I)                                                  
      D3=W(3,J)-W(3,I)                                                   
  317 D3=D3+H3I                                                             
      IF(D3.LE.CN(K))GO TO 317                                              
  318 D3=D3-H3I                                                             
      IF(D3.GT.CN(K))GO TO 318                                              
      IF(D3.LT.-CN(K))GO TO 327                                             
  325 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 325                                              
  326 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 326                                              
      IF(D1.LT.-CN(K))GO TO 318                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 326                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                         
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 326                                                             
  327 GO TO(328,329,328,330),ICT                                            
!-----Add face center molecules			
  328 ICT=ICT+1                                                             
      D1=D1+H1IC                                                            
      D3=D3+H3IC                                                            
      GO TO 317                                                             
!-----Add I molecules on second and third passes                         
  329 IF(JCT.EQ.0)GO TO 330                                              
      ICT=3                                                              
      D1=Y1A3+D1AJ+W(1,J)-W(1,I)                                         
      D2=D2AJ+W(2,J)-W(2,I)                                              
      D3=D3-H3IC                                                         
      GO TO 317                                                          
  330 CONTINUE                                                              
  331 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 930                                                
      NEND=N-1                                                              
      CALL MINHI(Y2A3,0.,X2A3,HL,KE)                                        
      IF(KCT)914,351,352                                                 
  351 X1A3=Y1A3                                                          
      X2A3=Y2A3                                                          
      IF(MCT)356,353,356                                                 
  352 IF(MCT)356,353,355                                                 
  353 IF(JE.GE.NV)GO TO 354                                              
      JE=JE+1                                                            
      Y1A3=Y1A3+DY1A3                                                    
      GO TO 305                                                          
  354 MCT=1                                                              
      JE=2*JE                                                            
      DY1A3=H1I/JE                                                       
      Y1A3=X1A3+DY1A3                                                    
      GO TO 305                                                          
  355 MCT=-1                                                             
      Y1A3=X1A3-DY1A3                                                    
      GO TO 305                                                          
  356 IF(JE.LT.KE)GO TO 354                                              
      IF(JE.GE.NNE)GO TO 358                                             
      KE=NNE                                                             
      GO TO 354                                                          
  358 CONTINUE                                                           
      IF(JCT)361,359,360                                                    
  359 H1A3P=X1A3                                                            
      H2A3P=X2A3                                                            
      JCT=1                                                                 
      IF(NCTR.EQ.0)GO TO 301                                             
      H1A3N=-H1A3P                                                       
      H2A3N=-H2A3P                                                       
      JCT=-1
      GO TO 302                                                             
  360 H1A3N=X1A3                                                         
      H2A3N=X2A3                                                         
	JCT=-1
	GO TO 302
  361 H1A3P=X1A3
      H2A3P=X2A3
	H21I=H1A3P-H1A3N
	H2I=H2A3P-H2A3N
  380 VX=H1IC*H2I*H3IC                                                      
      IF(IPR.LT.2)GO TO 390                                                 
      H3A3=0.                                                               
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
      SINE3=SQRT(1.-C(3)**2)                                             
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      ORIG(3)=0.                                                         
      ORIG(2)=H2A3P/SINE3                                                
      ORIG(1)=H1A3P-ORIG(2)*C(3)                                         
      WRITE(NQ,260)ALPHA(LCT),A1,A2,A3,VX,A,C                            
      DO 385 I=1,3                                                       
      GO TO(382,383,384),I                                               
  382 SM(3)=0.                                                           
      SM(2)=H2A3P/SINE3                                                  
      SM(1)=H1A3P-SM(2)*C(3)                                             
      SM(4)=(SM(1)-ORIG(1))/A(1)                                         
      SM(5)=(SM(2)-ORIG(2))/A(2)                                         
      SM(6)=0.                                                           
      D=SQRT(H1A3P**2+H2A3P**2)                                          
      WRITE(NQ,261)SYM(2),SM,D                                           
      GO TO 385                                                          
  383 SM(3)=H3IC                                                         
      SM(2)=0.                                                           
      SM(1)=H1IC                                                         
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=SM(2)/A(2)                                                   
      SM(6)=SM(3)/A(3)                                                   
      D=SQRT(H1IC**2+H3IC**2)                                            
      WRITE(NQ,261)SYM(4),SM,D                                           
      GO TO 385                                                          
  384 D1=H1A3P+H1IC                                                      
      D2=H2A3P                                                           
      SM(3)=H3IC                                                         
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=(SM(1)-ORIG(1))/A(1)                                         
      SM(5)=(SM(2)-ORIG(2))/A(2)                                         
      SM(6)=SM(3)/A(3)                                                   
      D=SQRT(D1**2+D2**2+SM(3)**2)                                       
      WRITE(NQ,261)SYM(6),SM,D                                           
  385 CONTINUE                                                           
  390 IF(VX.GE.V3)GO TO 101                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 101                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=H1A3P                                                          
      HT1(6)=H2A3P                                                          
      HT1(7)=0.                                                             
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=17                                                               
      MLOW=45                                                            
      GO TO 101                                                             
!-----Space group C2/c - Subgroup 3                                      
  400 VX=1000000.                                                        
      X3I=H3M                                                            
      X3IC=H3N                                                           
!     Search for best value of X1IC between its minimum, H1ICF, when     
!      molecules are in contact along axis-3 and twice this value.       
      H1ICFF=2.*H1ICF                                                    
      X1IC=H1ICF                                                         
      X1I=2.*X1IC                                                        
      DX1IC=H1ICF/NV                                                     
      KE=NE                                                              
      IE=1                                                               
      ICT=0                                                              
      GO TO 450                                                          
!-----Space group C2/c - Subgroup 2                                      
  425 VX=1000000.                                                        
!     Search for best value of X3IC between its minimum, H3N, when       
!      molecules are in contact along axis-3 and its maximum, H3ICF,     
!      when molecules are in contact along axis-1.                       
      H3ICR=H3ICF-H3N                                                    
      X3IC=H3N                                                           
      X3I=H3M                                                            
      DX3IC=H3ICR/NV                                                     
      X1IC=H1N                                                           
      KE=NE                                                              
      IE=1                                                               
      ICT=0                                                              
!-----Find value of X1IC for this value of X3IC                          
!     Collect corners for double standoff                                
  440 MARK=MARK3                                                         
      N=MARK                                                             
      DO 443 I=MARK2,NEND2,NSTP2                                         
      K=IT(I)                                                            
      L=K+10                                                             
      D3=TI(I+3)                                                         
  441 D3=D3+X3I                                                          
      IF(D3.LT.CN(K))GO TO 441                                           
  442 D3=D3-X3I                                                          
      IF(D3.GT.CN(K))GO TO 442                                           
      IF(D3.LT.-CN(K))GO TO 443                                          
      D=TI(I+1)+D3**2                                                    
      IF(D.GT.CN(L))GO TO 442                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 442                                                          
  443 CONTINUE                                                           
      NARK=N                                                             
!-----Add face center molecule distances                                 
      DO 446 I=MARK2,NEND2,NSTP2                                         
      K=IT(I)                                                            
      L=K+10                                                             
      D3=TI(I+3)+X3IC                                                    
  444 D3=D3+X3I                                                          
      IF(D3.LT.CN(K))GO TO 444                                           
  445 D3=D3-X3I                                                          
      IF(D3.GT.CN(K))GO TO 445                                           
      IF(D3.LT.-CN(K))GO TO 446                                          
      D=TI(I+1)+D3**2                                                    
      IF(D.GT.CN(L))GO TO 445                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+2)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 445                                                          
  446 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 930                                             
      NEND=N-1                                                           
      ER=ERMT                                                            
      CALL MINHI(X1IC,H1ICS,H1ICF,H1M,KE)                                
      IF(KCT)447,449,447                                                 
  447 IF(IPR.LT.3)GO TO 449                                              
      WRITE(NQ,448)X1IC,H1ICS,H1ICF                                      
  448 FORMAT(1X,'X1IC set at limit',3F8.3)                               
  449 X1I=2.*X1IC                                                        
  450 ER=ERMF                                                               
      JE=1                                                                  
      LE=1                                                                  
      DY1C=X1IC/NV                                                          
      DY3C=X3IC/NV                                                          
      Y1C=-X1IC                                                             
      Y3C=0.                                                                
      X2C=100.                                                           
      Y2C=H2N                                                            
      MCT=1                                                                 
      NCT=1                                                              
      MARK=MARK3                                                         
  451 N=MARK                                                                
      NARK=0                                                             
      DO 460 I=1,NMOD                                                       
      DO 459 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      JCT=0                                                              
      D1=Y1C-W(1,J)-W(1,I)                                                  
      D2=-W(2,J)-W(2,I)                                                     
      D3=Y3C-W(3,J)-W(3,I)                                                  
  454 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 454                                              
  455 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 455                                              
      IF(D1.LT.-CN(K))GO TO 458                                             
  456 D3=D3+X3I                                                             
      IF(D3.LE.CN(K))GO TO 456                                              
  457 D3=D3-X3I                                                             
      IF(D3.GT.CN(K))GO TO 457                                              
      IF(D3.LT.-CN(K))GO TO 455                                             
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 457                                            
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                         
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 457                                                             
  458 IF(JCT.NE.0)GO TO 459                                              
      JCT=1                                                              
      D1=D1+X1IC                                                         
      D3=D3+X3IC                                                         
      GO TO 454                                                          
  459 CONTINUE                                                              
  460 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 550                                                
      NEND=N-1                                                              
      CALL MINHI(Y2C,0.,X2C,H2M,KE)                                      
      IF(IPR.LT.3)GO TO 2460                                             
      WRITE(NQ,1460)Y1C,Y2C,Y3C                                          
 1460 FORMAT(1X,'Y1C =',F7.3,', Y2C =',F7.3,', Y3C =',F7.3)              
 2460 IF(KCT)552,462,461                                                 
  461 IF(MCT.LE.1)GO TO 464                                              
      J=MCT-NCT                                                          
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 476                                 
      GO TO 463                                                          
  462 X1C=Y1C                                                            
      X3C=Y3C                                                            
      X2C=Y2C                                                            
      NCT=MCT                                                            
  463 GO TO(464,469,470,471,472,473,474,475,468),MCT                     
  464 IF(LE.GE.NVT)GO TO 465                                             
      LE=LE+1                                                            
      Y1C=Y1C+DY1C                                                       
      GO TO 451                                                          
  465 IF(JE.GE.NV)GO TO 466                                              
      JE=JE+1                                                            
      Y3C=Y3C+DY3C                                                       
      LE=1                                                               
      Y1C=-X1IC                                                          
      GO TO 451                                                          
  466 NCT=6                                                              
  467 JE=2.*JE                                                           
      DY1C=X1IC/JE                                                       
      DY3C=X3IC/JE                                                       
      XX1C=X1C                                                           
      XX3C=X3C                                                           
      NCT=NCT+1                                                          
      IF(NCT.GT.9)NCT=NCT-8                                              
      GO TO(924,470,471,472,473,474,475,468,469),NCT                     
  468 MCT=2                                                              
      Y1C=XX1C+DY1C                                                      
      Y3C=XX3C                                                           
      GO TO 451                                                          
  469 MCT=3                                                              
      Y1C=XX1C+DY1C                                                      
      Y3C=XX3C+DY3C                                                      
      GO TO 451                                                          
  470 MCT=4                                                              
      Y1C=XX1C                                                           
      Y3C=XX3C+DY3C                                                      
      GO TO 451                                                          
  471 MCT=5                                                              
      Y1C=XX1C-DY1C                                                      
      Y3C=XX3C+DY3C                                                      
      GO TO 451                                                          
  472 MCT=6                                                              
      Y1C=XX1C-DY1C                                                      
      Y3C=XX3C                                                           
      GO TO 451                                                          
  473 MCT=7                                                              
      Y1C=XX1C-DY1C                                                      
      Y3C=XX3C-DY3C                                                      
      GO TO 451                                                          
  474 MCT=8                                                              
      Y1C=XX1C                                                           
      Y3C=XX3C-DY3C                                                      
      GO TO 451                                                          
  475 MCT=9                                                              
      Y1C=XX1C+DY1C                                                      
      Y3C=XX3C-DY3C                                                      
      GO TO 451                                                          
  476 IF(JE.LT.KE)GO TO 467                                                
      NSRCH=0                                                            
  477 IF(IPR.LT.3)GO TO 479                                              
      WRITE(NQ,478)X1C,X2C,X3C                                           
  478 FORMAT(1X,'X1C =',F7.3,', X2C =',F7.3,', X3C =',F7.3)              
  479 X1CT=2.*X1C                                                        
      X2CT=2.*X2C                                                        
      X2A3=-100.                                                            
      X2A3P=100.                                                         
      MCT=0                                                                 
      DY1A3=X1I/NV                                                          
      JE=1                                                                  
      Y1A3=-X1IC                                                            
      Y2A3=-1.5*H2N                                                         
      NARK=0                                                             
  480 N=MARK                                                                
      Y1AMC=Y1A3-X1C                                                     
      Y1TAMC=Y1AMC+Y1A3                                                  
!-----Collect C molecule distances on negative axis-2 side for double    
!     standoff.                                                          
      DO 487 I=1,NMOD                                                       
      DO 486 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      JCT=0                                                                 
      D1=Y1TAMC-W(1,J)-W(1,I)                                               
      D2=-W(2,J)-W(2,I)-X2C                                              
      D3=-W(3,J)-W(3,I)+X3C                                             
  481 D3=D3+X3I                                                             
      IF(D3.LE.CN(K))GO TO 481                                              
  482 D3=D3-X3I                                                             
      IF(D3.GT.CN(K))GO TO 482                                              
      IF(D3.LT.-CN(K))GO TO 485                                             
  483 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 483                                              
  484 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 484                                              
      IF(D1.LT.-CN(K))GO TO 482                                             
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 484                                            
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 484                                                             
!-----Has face center offset molecule been added?                        
  485 IF(JCT.NE.0)GO TO 486                                              
      JCT=1                                                              
      D1=D1+X1IC                                                         
      D3=D3+X3IC                                                         
      GO TO 481                                                          
  486 CONTINUE                                                           
  487 CONTINUE                                                           
      NARK=N                                                             
!-----Collect A and P molecules for negative standoff                    
      DO 500 I=1,NMOD                                                       
      DO 499 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      JCT=1                                                                 
      D1=Y1A3-W(1,J)-W(1,I)                                              
      D2=-W(2,J)-W(2,I)                                                  
      D3=W(3,J)-W(3,I)                                                      
  491 D3=D3+X3I                                                             
      IF(D3.LE.CN(K))GO TO 491                                              
  492 D3=D3-X3I                                                             
      IF(D3.GT.CN(K))GO TO 492                                              
      IF(D3.LT.-CN(K))GO TO 495                                             
  493 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 493                                              
  494 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 494                                              
      IF(D1.LT.-CN(K))GO TO 492                                             
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 494                                            
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 494                                                             
  495 GO TO(496,497,498,499),JCT                                         
!-----Add face center A3 molecule distances                              
  496 JCT=2                                                              
      D1=D1+X1IC                                                         
      D3=D3+X3IC                                                         
      GO TO 491                                                          
!-----Add P3 molecule distances                                          
  497 JCT=3                                                              
      D1=Y1AMC+W(1,J)-W(1,I)                                             
      D2=W(2,J)-W(2,I)-X2C                                               
      D3=X3C-W(3,J)-W(3,I)                                               
      GO TO 491                                                          
!-----Add face center P3 molecule distances                              
  498 JCT=4                                                                 
      D1=D1+X1IC                                                            
      D3=D3+X3IC                                                            
      GO TO 491                                                             
  499 CONTINUE                                                              
  500 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 554                                                
      NEND=N-1                                                              
      CALL MINHI(Y2A3,0.,X2A3,-H2M,KE)                                      
      IF(KCT)556,501,542                                                 
  501 NARK=0                                                             
!-----Check A3 on positive axis-2 side                                   
      Y1A3P=X1CT-Y1A3                                                    
      Y2A3P=X2CT-Y2A3                                                    
      N=MARK                                                             
      DO 508 I=1,NMOD                                                       
      DO 507 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      JCT=0                                                                 
      D1=Y1A3P-W(1,J)-W(1,I)                                                
      D2=-W(2,J)-W(2,I)                                                  
      D3=W(3,J)-W(3,I)                                                      
  502 D3=D3+X3I                                                             
      IF(D3.LE.CN(K))GO TO 502                                              
  503 D3=D3-X3I                                                             
      IF(D3.GT.CN(K))GO TO 503                                              
      IF(D3.LT.-CN(K))GO TO 506                                             
  504 D1=D1+X1I                                                             
      IF(D1.LE.CN(K))GO TO 504                                              
  505 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 505                                              
      IF(D1.LT.-CN(K))GO TO 503                                             
      D=D1**2+D3**2                                                      
      IF(D.GT.CN(L))GO TO 505                                            
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 505                                                             
  506 IF(JCT.NE.0)GO TO 507                                              
!-----Add face center A3 molecule distances                              
      JCT=1                                                              
      D1=D1+X1IC                                                         
      D3=D3+X3IC                                                         
      GO TO 502                                                          
  507 CONTINUE                                                              
  508 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 510                                                
      NEND=N-1                                                              
      Z2A3=Y2A3P                                                         
      CALL MINHI(Z2A3,Y2A3P,X2A3P,H2M,KE)                                
      IF(KCT)510,509,542                                                 
  509 Y2A3=-Z2A3                                                         
  510 CONTINUE                                                           
  541 X1A3=Y1A3                                                          
      X2A3=Y2A3                                                          
      X2A3P=X2CT-X2A3                                                    
      IF(MCT)546,543,546                                                 
  542 IF(MCT)546,543,545                                                 
  543 IF(JE.GE.NV)GO TO 544                                              
      JE=JE+1                                                            
      Y1A3=Y1A3+DY1A3                                                    
      GO TO 480                                                          
  544 MCT=1                                                              
      JE=2.*JE                                                           
      DY1A3=X1I/JE                                                       
      Y1A3=X1A3+DY1A3                                                    
      GO TO 480                                                          
  545 MCT=-1                                                             
      Y1A3=X1A3-DY1A3                                                    
      GO TO 480                                                          
  546 IF(JE.LT.KE)GO TO 544                                              
      X2P3=X2C-X2A3                                                      
      X2I=2.*X2P3                                                        
      VY=X1IC*X2P3*X3IC                                                     
      IF(IPR.LT.4)GO TO 549                                                 
      WRITE(NQ,548)X1I,X2I,X3I,VY                                        
  548 FORMAT(' X1I =',F7.3,', X2I =',F7.3,', X3I =',F7.3,', VY =',F8.2)  
  549 IF(LCT-4)640,600,610                                               
!-----A reduced limit for the X1IC scan of Code DE may be required       
  550 IF(LCT.LT.4)GO TO 570                                              
      H1ICFF=X1IC                                                        
      IF(IPR.LT.3)GO TO 560                                              
      WRITE(NQ,551)H1ICFF                                                
  551 FORMAT(1X,'H1ICFF =',F7.3,', Set by no interactions for C')        
      GO TO 560                                                          
  552 IF(LCT.LT.4)GO TO 570                                              
      H1ICFF=X1IC                                                        
      IF(IPR.LT.3)GO TO 560                                              
      WRITE(NQ,553)H1ICFF                                                
  553 FORMAT(1X,'H1ICFF =',F7.3,', Set by H2C becoming zero')            
      GO TO 560                                                          
  554 IF(LCT.LT.4)GO TO 930                                              
      H1ICFF=X1IC                                                        
      IF(IPR.LT.3)GO TO 560                                              
      WRITE(NQ,555)H1ICFF                                                
  555 FORMAT(1X,'H1ICFF =',F7.3,', Set by no interactions for A3')       
      GO TO 560                                                          
  556 IF(LCT.LT.4)GO TO 915                                              
      H1ICFF=X1IC                                                        
      IF(IPR.LT.3)GO TO 560                                              
      WRITE(NQ,557)H1ICFF                                                
  557 FORMAT(1X,'H1ICFF =',F7.3,', Set by H2A3 becoming zero')           
  560 X1C=Y1C                                                            
      X2C=0.                                                             
      X3C=Y3C                                                            
      NSRCH=1                                                            
      IF(ICT.EQ.0)IE=NV                                                  
      GO TO 477                                                          
!-----Code DA, DB, or DC - A value of zero was found for Y2C             
  570 X1C=Y1C                                                            
      X2C=0.                                                             
      X3C=Y3C                                                            
      GO TO 477                                                          
!-----Is this the smallest volume found by Code DD so far?               
  600 IF(VY.GE.VX)GO TO 601                                              
      VX=VY                                                              
      H1I=X1I                                                            
      H2I=X2I                                                            
      H3I=X3I                                                            
      H1IC=X1IC                                                          
      H3IC=X3IC                                                          
      H1C=X1C                                                            
      H2C=X2C                                                            
      H3C=X3C                                                            
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
      IF(NSRCH.GT.0)GO TO 650                                            
      IF(ICT)603,602,603                                                 
  601 IF(NSRCH.GT.0)GO TO 650                                            
      IF(ICT)603,602,605                                                 
  602 IF(IE.GE.NV)GO TO 604                                              
      IE=IE+1                                                            
      X3IC=X3IC+DX3IC                                                    
      X3I=2.*X3IC                                                        
      GO TO 440                                                          
  603 IF(IE.GE.KE)GO TO 606                                              
  604 IE=2*IE                                                            
      DX3IC=H3ICR/IE                                                     
      X3IC=H3IC+DX3IC                                                    
      X3I=2.*X3IC                                                        
      ICT=1                                                              
      GO TO 440                                                          
  605 X3IC=H3IC-DX3IC                                                    
      IF(X3IC.LT.H3N)GO TO 603                                           
      X3I=2.*X3IC                                                        
      ICT=-1                                                             
      GO TO 440                                                          
  606 IF(IE.GE.NNE)GO TO 650                                             
      KE=NNE                                                             
      GO TO 604                                                          
!-----Is this the smallest volume found by Code DE so far?               
  610 IF(VY.GE.VX)GO TO 611                                              
      VX=VY                                                              
      H1I=X1I                                                            
      H2I=X2I                                                            
      H3I=X3I                                                            
      H1IC=X1IC                                                          
      H3IC=X3IC                                                          
      H1C=X1C                                                            
      H2C=X2C                                                            
      H3C=X3C                                                            
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
      IF(NSRCH.GT.0)GO TO 650                                            
      IF(ICT)613,612,613                                                 
  611 IF(NSRCH.GT.0)GO TO 650                                            
      IF(ICT)613,612,615                                                 
  612 IF(IE.GE.NV)GO TO 614                                              
      IE=IE+1                                                            
      X1IC=X1IC+DX1IC                                                    
      X1I=2.*X1IC                                                        
      GO TO 450                                                          
  613 IF(IE.GE.KE)GO TO 616                                              
  614 IE=2*IE                                                            
      DX1IC=H1ICF/IE                                                     
      IF(X1IC.GT.H1ICFF)GO TO 613                                        
      X1IC=H1IC+DX1IC                                                    
      X1I=2.*X1IC                                                        
      ICT=1                                                              
      GO TO 450                                                          
  615 X1IC=H1IC-DX1IC                                                    
      IF(X1IC.LT.H1ICF)GO TO 613                                         
      X1I=2.*X1IC                                                        
      ICT=-1                                                             
      GO TO 450                                                          
  616 IF(IE.GE.NNE)GO TO 650                                             
      KE=NNE                                                             
      GO TO 614                                                          
  640 VX=VY                                                              
      H1I=X1I                                                            
      H2I=X2I                                                            
      H3I=X3I                                                            
      H1IC=X1IC                                                          
      H3IC=X3IC                                                          
      H1C=X1C                                                            
      H2C=X2C                                                            
      H3C=X3C                                                            
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
  650 H1P3=H1C-H1A3                                                      
      IF(H1P3.GT.H1IC)H1P3=H1P3-H1I                                      
      IF(H1P3.LT.-H1IC)H1P3=H1P3+H1I                                     
      H2P3=H2C-H2A3                                                      
      H3P3=H3C                                                           
      H21I=2.*H1P3                                                       
      IF(IPR.LT.2)GO TO 670                                                 
      H3A3=0.                                                               
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
      SINE3=SQRT(1.-C(3)**2)                                             
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
      ORIG(3)=H3C                                                        
      ORIG(2)=H2C/SINE3                                                  
      ORIG(1)=H1C-ORIG(2)*C(3)                                           
      WRITE(NQ,260)ALPHA(LCT),A1,A2,A3,VX,A,C                           
      DO 660 I=1,7                                                       
      GO TO(652,653,654,655,656,657,658),I                               
  652 D1=H1P3                                                            
      D2=H2P3                                                            
      SM(3)=H3P3                                                         
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=SM(2)/A(2)                                                   
      SM(6)=(SM(3)-ORIG(3))/A(3)                                         
      GO TO 659                                                          
  653 D1=H1A3                                                            
      D2=H2A3                                                            
      SM(3)=0.                                                           
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=(SM(1)-ORIG(1))/A(1)                                         
      SM(5)=(SM(2)-ORIG(2))/A(2)                                         
      SM(6)=0.                                                           
      GO TO 659                                                          
  654 D1=H1C                                                             
      D2=H2C                                                             
      SM(3)=H3C                                                          
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=(SM(1)-ORIG(1))/A(1)                                         
      SM(5)=(SM(2)-ORIG(2))/A(2)                                         
      SM(6)=(SM(3)-ORIG(3))/A(3)                                         
      GO TO 659                                                          
  655 D1=H1IC                                                            
      D2=0.                                                              
      SM(3)=H3IC                                                         
      SM(2)=0.                                                           
      SM(1)=D1                                                           
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=0.                                                           
      SM(6)=SM(3)/A(3)                                                   
      GO TO 659                                                          
  656 D1=H1IC+H1P3                                                       
      D2=H2P3                                                            
      SM(3)=H3IC+H3P3                                                    
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=SM(1)/A(1)                                                   
      SM(5)=SM(2)/A(2)                                                   
      SM(6)=(SM(3)-ORIG(3))/A(3)                                         
      GO TO 659                                                          
  657 D1=H1IC+H1A3                                                       
      D2=H2A3                                                            
      SM(3)=H3IC                                                         
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=(SM(1)-ORIG(1))/A(1)                                         
      SM(5)=(SM(2)-ORIG(2))/A(2)                                         
      SM(6)=SM(3)/A(3)                                                   
      GO TO 659                                                          
  658 D1=H1IC+H1C                                                        
      D2=H2C                                                             
      SM(3)=H3IC+H3C                                                     
      SM(2)=D2/SINE3                                                     
      SM(1)=D1-SM(2)*C(3)                                                
      SM(4)=(SM(1)-ORIG(1))/A(1)                                         
      SM(5)=(SM(2)-ORIG(2))/A(2)                                         
      SM(6)=(SM(3)-ORIG(3))/A(3)                                         
  659 D=SQRT(D1**2+D2**2+SM(3)**2)                                       
      WRITE(NQ,261)SYM(I),SM,D                                           
  660 CONTINUE                                                           
  670 IF(VX.GE.V3)GO TO 101                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 101                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=H1C                                                            
      HT1(6)=H2C                                                            
      HT1(7)=H3C                                                         
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(LCT)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=18                                                               
      MLOW=LCT+43                                                        
      GO TO 101                                                             
  902 WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK3,NSTP3,MARK,NSTP,LIMIT   
      WRITE(17,903)MARK1,NSTP1,MARK2,NSTP2,MARK3,NSTP3,MARK,NSTP,LIMIT   
  903 FORMAT(26H0STORAGE EXCEEDED BY ZNRTA,5X,4(I9,I3),I9)                  
      GO TO 999                                                             
  911 NERR=1                                                             
      GO TO 918                                                          
  913 NERR=3                                                             
      GO TO 918                                                          
  914 NERR=4                                                             
      GO TO 918                                                          
  915 NERR=5                                                             
  918 WRITE(NQ,919)NERR,A1,A2,A3                                         
      WRITE(17,919)NERR,A1,A2,A3                                         
  919 FORMAT(' PARAMETER TOO SMALL - ZNRTA Location',I2,3F7.1)           
      GO TO 999                                                          
  921 NERR=1                                                             
      GO TO 928                                                          
  922 NERR=2                                                             
      GO TO 928                                                          
  923 NERR=3                                                             
      GO TO 928                                                          
  924 NERR=4                                                             
  928 WRITE(NQ,929)NERR,A1,A2,A3                                         
      WRITE(17,929)NERR,A1,A2,A3                                         
  929 FORMAT(' PARAMETER TOO LARGE - ZNRTA Location',I2,3F7.1)           
      GO TO 999                                                          
  930 WRITE(NQ,931)KCT,MARK1,MARK2,NEND2,MARK                               
      WRITE(17,931)KCT,MARK1,MARK2,NEND2,MARK                               
  931 FORMAT(' NO INTERACTIONS in ZNRTA',I3,4I9)                         
  999 continue
! 999 KILL=1              ! 3-9-2010 temp                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZNRTA         
