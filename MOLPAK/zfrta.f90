      SUBROUTINE ZFRTA   
!                                                                        
!     For AI, AJ, AK, Al
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
!      COMMON ICOUNT,MCOUNT,IORDER(500)                                 
!      COMMON ICOUNT,MCOUNT,IORDER(5000)                                 
!      COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500) 
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
      CHARACTER(2) :: ALPHA(4) =(/'AI','AJ','AK','AL'/)
      CHARACTER(2) :: SYM(3) =  (/'C ','A3','P3'/)      
!      DIMENSION ALPHA(4),A(12),C(6),SM(3)                                   
!      DATA ALPHA/'AI','AJ','AK','AL'/                                       
!      DATA SYM/'C ','A3','P3'/                                           
!
      INTEGER :: I, ICT, J, JCT, JE, K, KE, L, LCT, LCTM, LE
      INTEGER :: MCT, N, NCT, NERR
!
      REAL :: A(12), C(6) 
      REAL :: D, D1, D2, D3, DX3P3, DY1A3, DZ1A3 
      REAL :: H1A3, H1C, H1I, H1IH, H1P3 
      REAL :: H2A3, H2C, H2I, H21I  
      REAL :: H3A3, H3C, H3P3, HL
      REAL :: SM(3), SINE3
      REAL :: VX, VY, V3S
      REAL :: X1A3, X1I, X1IH, X1P3
      REAL :: X2I, X21I, X2A3, X3P3 
      REAL :: Y1A3, Y2A3
      REAL :: Z1A3, Z1A3T, Z2A3, Z3CT
!
!-----Set starting value of V3                                           
      V3S=V3                                                             
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)(ITR(I),I=10,13)                                        
    5 FORMAT(1X,'ZFRTA called',4I3)                                      
   10 GO TO(100,100,300,301,401),ITR(10)                                    
!-----IDENTITIES IN CONTACT ALONG B AXIS - AXIS-3                           
!                                                                           
!     CALC. DUAL CONTACT P DISTANCE ALONG AXIS-1                            
  100 KE=NNE                                                                
      N=MARK2                                                               
      NSTP2=6                                                            
      DO 109 I=1,NMOD                                                       
      DO 108 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=W(2,J)-W(2,I)                                                      
      IF(ABS(D2).GT.CN(K))GO TO 108                                         
      IT(N)=K                                                               
      TI(N+1)=D2**2                                                         
      TI(N+2)=W(1,J)-W(1,I)                                                 
      TI(N+3)=-W(3,J)-W(3,I)                                                
      IT(N+4)=L                                                             
      TI(N+5)=W(3,J)-W(3,I)                                              
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
  108 CONTINUE                                                              
  109 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
      H1P3=100.                                                             
      X1P3=H1N                                                           
      ER=ERMT                                                               
      KE=NE                                                              
      JE=1                                                                  
      DX3P3=H3M/NV                                                          
      X3P3=-H3N                                                             
      MCT=0                                                                 
  116 N=MARK                                                                
!-----Collect I molecule distances for double standoff                   
      DO 119 I=MARK2,NEND2,NSTP2
	K=IT(I)
	D3=TI(I+5)
  117 D3=D3+H3M
      IF(D3.LE.CN(K))GO TO 117
  118 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 118
	IF(D3.LT.-CN(K))GO TO 119
	D=TI(I+1)+D3**2
	L=IT(I+4)
	IF(D.GT.CN(L))GO TO 118
	IT(N)=K
	TI(N+1)=D
	TI(N+2)=TI(I+2)
	N=N+NSTP
      IF(N.GT.LIMIT)GO TO 902
      GO TO 118
  119 CONTINUE
      NARK=N
  120 N=NARK
      DO 129 I=MARK2,NEND2,NSTP2                                            
      K=IT(I)                                                               
      D3=X3P3+TI(I+3)                                                       
  121 D3=D3+H3M                                                             
      IF(D3.LE.CN(K))GO TO 121                                              
  122 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 122                                              
      IF(D3.LT.-CN(K))GO TO 129                                             
      D=TI(I+1)+D3**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 122                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 122                                                             
  129 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(X1P3,0.,H1P3,H1M,KE)                                       
      IF(KCT)1904,151,152                                               
  151 H3P3=X3P3                                                         
      H1P3=X1P3                                                         
      IF(MCT)156,153,156                                                
  152 IF(MCT)156,153,155                                                
  153 IF(JE.GE.NV)GO TO 154                                             
      JE=JE+1                                                           
      X3P3=X3P3+DX3P3                                                   
      GO TO 120                                                         
  154 MCT=1                                                             
      JE=2*JE                                                           
      DX3P3=H3M/JE                                                      
      X3P3=H3P3+DX3P3                                                   
      GO TO 120                                                         
  155 MCT=-1                                                            
      X3P3=H3P3-DX3P3                                                   
      IF(X3P3.LT.-H3N)X3P3=X3P3+H3M                                      
      GO TO 120                                                         
  156 IF(JE.LT.KE)GO TO 154                                             
      IF(JE.GE.NNE)GO TO 158                                            
      KE=NNE                                                            
      GO TO 154                                                         
  158 CONTINUE                                                          
      NARK=0
      H1I=2.*H1P3                                                        
      H1IH=H1P3                                                          
      IF(IPR.LT.4)GO TO 161                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,160)H1P3,H3P3                                                
  160 FORMAT(15X,6HH1P3 =,F7.3,8H, H3P3 =,F7.3)                             
  161 CONTINUE
!-----CALC. AXIS-2 SPACING OF PLANE OF AXIS AND CENTER MOLECULES            
  200 MARK=MARK2                                                            
      ER=ERMT                                                            
      IF(ITR(10)-2)201,202,300                                              
!-----Screw axis - space group P21/c - Code AI                           
  201 H3A3=H3N                                                              
      LCT=1                                                                 
      GO TO 210                                                             
!-----Twofold axis - space group P2/c - Code AJ                          
  202 H3A3=0.                                                               
      LCT=2                                                                 
  210 ER=ERMF                                                               
!-----Set up to calculate standoff in positive axis-2 direction          
      JCT=3                                                              
      Y2A3=H2N                                                              
      HL=H2M                                                                
      X2A3=100.                                                             
  211 JE=1                                                                  
      DY1A3=H1I/NV                                                          
      Y1A3=-H1P3                                                            
      MCT=0                                                                 
      KE=NE                                                                 
  212 N=MARK                                                                
      NARK=0
      IF(JCT.EQ.3)GO TO 220
!-----Collect distances for second layer of A3 and C molecules           
      Z1A3T=2.*Y1A3-H1A3
      Z3CT=H3A3+H3P3
      DO 219 I=1,NMOD
      DO 218 J=1,NMOD
	K=IA(I)+IAA(J)
	L=K+10
	ICT=0
	D1=Z1A3T-W(1,J)-W(1,I)
	D2=-H2A3-W(2,J)-W(2,I)
	D3=H3A3+W(3,J)-W(3,I)
  213 D1=D1+H1I
      IF(D1.LT.CN(K))GO TO 213
  214 D1=D1-H1I
      IF(D1.GT.CN(K))GO TO 214
	IF(D1.LT.-CN(K))GO TO 217
  215 D3=D3+H3M
      IF(D3.LT.CN(K))GO TO 215
  216 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 216
      IF(D3.LT.-CN(K))GO TO 214
	D=D1**2+D3**2
	IF(D.GT.CN(L))GO TO 216
	IT(N)=K
	TI(N+1)=D
  	TI(N+2)=D2
 	N=N+NSTP
	IF(N.GT.LIMIT)GO TO 902
	GO TO 216
!-----Have C molecule distances been calculated                          
  217 IF(ICT.GT.0)GO TO 218
      ICT=1
	D1=D1+H1P3
	D3=Z3CT-W(3,J)-W(3,I)
	GO TO 213
  218 CONTINUE
  219 CONTINUE
      NARK=N
  220 DO 230 I=1,NMOD                                                       
      DO 229 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=-W(2,J)-W(2,I)                                                     
      ICT=JCT                                                               
      D1=Y1A3-W(1,J)-W(1,I)                                                 
      D3=H3A3+W(3,J)-W(3,I)                                                 
  221 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 221                                              
  222 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 222                                              
      IF(D1.LT.-CN(K))GO TO 225                                             
  223 D3=D3+H3M                                                             
      IF(D3.LE.CN(K))GO TO 223                                              
  224 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 224                                              
      IF(D3.LT.-CN(K))GO TO 222                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 224                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 224                                                             
!-----Have OTHER molecule contacts been calculated?                      
  225 GO TO(226,227,228,229),ICT                                         
!-----PLANE CONTACT                                                         
  226 D1=H1P3+Y1A3-H1A3+W(1,J)-W(1,I)                                       
      D2=-H2A3+W(2,J)-W(2,I)                                             
      D3=H3P3-W(3,J)-W(3,I)                                              
      ICT=2                                                              
      GO TO 221                                                             
!-----IDENTITY CONTACT                                                      
  227 D1=D1-H1P3                                                            
      D3=W(3,J)-W(3,I)                                                      
      ICT=3                                                              
      GO TO 221                                                             
!-----CENTER CONTACT                                                        
  228 D1=Y1A3+H1P3-W(1,J)-W(1,I)                                            
      D2=-W(2,J)-W(2,I)                                                     
      D3=H3A3+H3P3-W(3,J)-W(3,I)                                            
      ICT=4                                                              
      GO TO 221                                                             
  229 CONTINUE                                                              
  230 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(Y2A3,0.,X2A3,HL,KE)                                        
      IF(KCT)2904,251,252                                               
  251 X1A3=Y1A3                                                         
      X2A3=Y2A3                                                         
      IF(MCT)256,253,256                                                
  252 IF(MCT)256,253,255                                                
  253 IF(JE.GE.NV)GO TO 254                                             
      JE=JE+1                                                           
      Y1A3=Y1A3+DY1A3                                                   
      GO TO 212                                                         
  254 MCT=1                                                             
      JE=2.*JE                                                          
      DY1A3=H1I/JE                                                      
      Y1A3=X1A3+DY1A3                                                   
      GO TO 212                                                         
  255 MCT=-1                                                            
      Y1A3=X1A3-DY1A3                                                   
      GO TO 212                                                         
  256 IF(JE.LT.KE)GO TO 254                                             
      IF(JE.GE.NNE)GO TO 258                                            
      KE=NNE                                                            
      GO TO 254                                                         
  258 CONTINUE                                                          
!-----Has negative standoff been determined                              
      IF(JCT.EQ.1)GO TO 265
      H1A3=X1A3
	H2A3=X2A3
	JCT=1
	Y2A3=-H2N
	X2A3=-100.
	HL=-H2M
      GO TO 211
  265 H21I=H1A3-X1A3
      H2I=H2A3-X2A3
	VX=H1IH*H2I*H3N
      IF(VX.GE.V3)GO TO 297                                                 
!-----Check upper layers on positive axis-2 side                         
      Z1A3T=2.*H1A3-X1A3
      Z3CT=H3A3+H3P3
	N=MARK
      DO 279 I=1,NMOD
      DO 278 J=1,NMOD
	K=IA(I)+IAA(J)
	L=K+10
	ICT=0
	D1=Z1A3T-W(1,J)-W(1,I)
	D2=-X2A3-W(2,J)-W(2,I)
	D3=H3A3+W(3,J)-W(3,I)
  271 D1=D1+H1I
      IF(D1.LT.CN(K))GO TO 271
  272 D1=D1-H1I
      IF(D1.GT.CN(K))GO TO 272
	IF(D1.LT.-CN(K))GO TO 275
  273 D3=D3+H3M
      IF(D3.LT.CN(K))GO TO 273
  274 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 274
      IF(D3.LT.-CN(K))GO TO 272
	D=D1**2+D3**2
	IF(D.GT.CN(L))GO TO 274
	IT(N)=K
	TI(N+1)=D
	TI(N+2)=D2
	N=N+NSTP
	IF(N.GT.LIMIT)GO TO 902
	GO TO 274
!-----Have C molecule distances been calculated                          
  275 IF(ICT.GT.0)GO TO 278
      ICT=1
	D1=D1+H1P3
	D3=Z3CT-W(3,J)-W(3,I)
	GO TO 271
  278 CONTINUE
  279 CONTINUE
      NARK=N
      DO 290 I=1,NMOD                                                       
      DO 289 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      ICT=1                                                                
!-----AXIS CONTACT                                                          
      D1=H1A3-W(1,J)-W(1,I)                                                 
      D2=-W(2,J)-W(2,I)                                                  
      D3=H3A3+W(3,J)-W(3,I)                                                 
  281 D1=D1+H1I                                                             
      IF(D1.LT.CN(K))GO TO 281                                              
  282 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 282                                              
      IF(D1.LT.-CN(K))GO TO 285                                             
  283 D3=D3+H3M                                                             
      IF(D3.LT.CN(K))GO TO 283                                              
  284 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 284                                              
      IF(D3.LT.-CN(K))GO TO 282                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 284                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 284                                                             
!-----Have OTHER molecule contacts been calculated?                      
  285 GO TO(286,287,288,289),ICT                                         
!-----PLANE CONTACT                                                         
  286 D1=H1P3+H1A3-X1A3+W(1,J)-W(1,I)                                       
      ICT=2                                                              
      D2=-X2A3+W(2,J)-W(2,I)                                             
      GO TO 281                                                             
!-----IDENTITY CONTACT                                                      
  287 D1=D1-H1P3                                                            
      D3=W(3,J)-W(3,I)                                                      
      ICT=3                                                              
      GO TO 281                                                             
!-----CENTER CONTACT                                                        
  288 D1=H1A3+H1P3-W(1,J)-W(1,I)                                            
      D2=-W(2,J)-W(2,I)                                                     
      D3=H3A3+H3P3-W(3,J)-W(3,I)                                            
      ICT=4                                                              
      GO TO 281                                                             
  289 CONTINUE                                                              
  290 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      Y2A3=H2A3
      CALL MINHI(Y2A3,H2A3,100.,H2M,KE)                                     
      NARK=0
      IF(KCT)292,291,906
  291 H2A3=Y2A3
      H2I=H2A3-X2A3
	VX=H1IH*H2I*H3N
	IF(VX.GE.V3)GO TO 297
  292 V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 297                                             
      H1C=H1A3+H1P3                                                      
      H2C=H2A3                                                           
      H3C=H3A3+H3P3                                                      
!     IF(H1C.GT.H1IH)H1C=H1C-H1I                                         
      IF(H3C.GT.H3N)H3C=H3C-H3M                                          
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                            
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
      NLOW=8                                                                
      MLOW=LCT+9                                                         
  297 IF(IPR.LT.2)GO TO 298                                              
      LCTM=LCT                                                           
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3M                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      SINE3=SQRT(1.-C(3)**2)                                             
      A(6)=H3C                                                           
      A(5)=H2C/SINE3                                                     
      A(4)=H1C-A(5)*C(3)                                                 
      A(9)=H3P3                                                          
      A(8)=0.                                                            
      A(7)=H1P3                                                          
      A(12)=H3A3                                                         
      A(11)=H2A3/SINE3                                                   
      A(10)=H1A3-A(11)*C(3)                                              
!-----Then to fractional coordinates                                     
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      A(7)=A(7)/A(1)                                                     
      A(8)=A(8)/A(2)                                                     
      A(9)=A(9)/A(3)                                                     
      A(10)=A(10)/A(1)                                                   
      A(11)=A(11)/A(2)                                                   
      A(12)=A(12)/A(3)                                                   
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
  298 IF(LCT.NE.1)GO TO 299                                                 
      GO TO(202,300,301,401)ITR(11)                                         
  299 IF(ITR(12)-2)300,301,401                                           
!-----Screw axis - space group P21/c - Code AK                           
  300 H3A3=H3N                                                           
      LCT=3                                                              
      GO TO 302                                                          
!-----Twofold axis - space group P2/c - Code AL                          
  301 LCT=4                                                              
      H3A3=0.                                                            
  302 KE=NE                                                                 
      N=MARK2                                                               
      NSTP2=5                                                            
!-----Collect close interatomic distances along axis-2                   
      DO 304 I=1,NMOD                                                       
      DO 303 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=W(2,J)-W(2,I)                                                      
      IF(ABS(D2).GT.CN(K))GO TO 303                                         
      IT(N)=K                                                               
      TI(N+1)=D2**2                                                         
      TI(N+2)=W(1,J)-W(1,I)                                                 
      TI(N+3)=-W(3,J)-W(3,I)                                                
      IT(N+4)=L                                                             
      N=N+NSTP2                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
  303 CONTINUE                                                              
  304 CONTINUE                                                              
      IF(N.LE.MARK2)GO TO 920                                               
      NEND2=N-1                                                             
      MARK=N                                                                
      X1P3=H1N                                                           
      DX3P3=H3M/NV                                                          
      X3P3=-H3N                                                             
!-----Top of P molecule axis-3 offset placement loop                     
! 305 ER=ERM                                                             
  305 NCT=0                                                                 
      VX=1000000.                                                           
      LE=1                                                                  
!-----CALC. SINGLE CONTACT P DISTANCE                                       
  306 N=MARK                                                                
      ER=ERM                                                             
      DO 309 I=MARK2,NEND2,NSTP2                                            
      K=IT(I)                                                               
      D3=X3P3+TI(I+3)                                                       
  307 D3=D3+H3M                                                             
      IF(D3.LE.CN(K))GO TO 307                                              
  308 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 308                                              
      IF(D3.LT.-CN(K))GO TO 309                                             
      D=TI(I+1)+D3**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 308                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 308                                                             
  309 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(X1P3,H1N,100.,H1M,KE)                                      
      IF(IPR.LT.4)GO TO 310                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,160)X1P3,X3P3                                                
!-----Calculate offset of A and C molecules in positive axis-2           
!     direction then A, C, P, and I in negative direction                
  310 X1I=2.*X1P3                                                           
      X1IH=X1P3                                                          
!-----Set up for calculation of standoff in positive direction           
      Z2A3=H2N                                                              
      Y2A3=100.                                                             
      Y1A3=H1A3                                                          
      HL=H2M                                                             
      ER=ERMF                                                               
      JCT=3                                                                 
  312 MCT=0                                                                 
!-----Skip scan of Z1A3 range in last stages of refinement               
      IF(LE.LE.NE)GO TO 315                                              
      JE=LE/2                                                            
      DZ1A3=X1I/JE                                                       
      Z1A3=Y1A3                                                          
      GO TO 316                                                          
  315 JE=1                                                                  
      DZ1A3=X1I/NV                                                          
      Z1A3=-X1IH                                                            
  316 N=MARK                                                                
      NARK=0
      IF(JCT.EQ.3)GO TO 340
!-----Collect distances for second layer of A3 and C molecules           
      Z1A3T=2.*Z1A3-X1A3
      Z3CT=H3A3+X3P3
!     WRITE(NQ,317)Z1A3T,Z3CT
! 317 FORMAT(' Z1A3T =',F8.3,'  Z3CT =',F8.3)
      DO 329 I=1,NMOD
      DO 328 J=1,NMOD
      K=IA(I)+IAA(J)
	  L=K+10
	  ICT=0
	  D1=Z1A3T-W(1,J)-W(1,I)
	  D2=-X2A3-W(2,J)-W(2,I)
	  D3=H3A3+W(3,J)-W(3,I)
  321 D1=D1+X1I
      IF(D1.LT.CN(K))GO TO 321
  322 D1=D1-X1I
      IF(D1.GT.CN(K))GO TO 322
	  IF(D1.LT.-CN(K))GO TO 325
  323 D3=D3+H3M
      IF(D3.LT.CN(K))GO TO 323
  324 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 324
      IF(D3.LT.-CN(K))GO TO 322
	  D=D1**2+D3**2
 	  IF(D.GT.CN(L))GO TO 324
	  IT(N)=K
	  TI(N+1)=D
	  TI(N+2)=D2
	  N=N+NSTP
	  IF(N.GT.LIMIT)GO TO 902
	  GO TO 324
!-----Have C molecule distances been calculated                          
  325 IF(ICT.GT.0)GO TO 328
      ICT=1
	  D1=D1+X1P3
	  D3=Z3CT-W(3,J)-W(3,I)
	  GO TO 321
  328 CONTINUE
  329 CONTINUE
!     WRITE(NQ,330)MARK,N
! 330 FORMAT(' MARK =',I10,'   N =',I10)
      NARK=N
  340 DO 350 I=1,NMOD                                                       
      DO 349 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      ICT=JCT                                                              
!-----AXIS CONTACT                                                          
      D1=Z1A3-W(1,J)-W(1,I)                                                 
      D2=-W(2,J)-W(2,I)                                                  
      D3=H3A3+W(3,J)-W(3,I)                                                 
  341 D1=D1+X1I                                                             
      IF(D1.LT.CN(K))GO TO 341                                              
  342 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 342                                              
      IF(D1.LT.-CN(K))GO TO 345                                             
  343 D3=D3+H3M                                                             
      IF(D3.LT.CN(K))GO TO 343                                              
  344 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 344                                              
      IF(D3.LT.-CN(K))GO TO 342                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 344                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 344                                                             
!-----Have OTHER molecule contacts been calculated?                      
  345 GO TO(346,347,348,349),ICT                                         
!-----PLANE CONTACT                                                         
  346 D1=X1P3+Z1A3-X1A3+W(1,J)-W(1,I)                                       
      D2=-X2A3+W(2,J)-W(2,I)                                             
      D3=X3P3-W(3,J)-W(3,I)                                              
      ICT=2                                                              
      GO TO 341                                                             
!-----IDENTITY CONTACT                                                      
  347 D1=D1-X1P3                                                            
      D3=W(3,J)-W(3,I)                                                      
      ICT=3                                                              
      GO TO 341                                                             
!-----CENTER CONTACT                                                        
  348 D1=Z1A3+X1P3-W(1,J)-W(1,I)                                            
      D2=-W(2,J)-W(2,I)                                                     
      D3=H3A3+X3P3-W(3,J)-W(3,I)                                            
      ICT=4                                                              
      GO TO 341                                                             
  349 CONTINUE                                                              
  350 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(Z2A3,0.,Y2A3,HL,KE)                                        
      IF(KCT)358,351,352
  351 Y1A3=Z1A3                                                         
      Y2A3=Z2A3                                                         
      IF(MCT)356,353,356                                                
  352 IF(MCT)356,353,355                                                
  353 IF(JE.GE.NV)GO TO 354                                             
      JE=JE+1                                                           
      Z1A3=Z1A3+DZ1A3                                                   
      GO TO 316                                                         
  354 MCT=1                                                             
      JE=2*JE                                                           
      DZ1A3=X1I/JE                                                      
      Z1A3=Y1A3+DZ1A3                                                   
      GO TO 316                                                         
  355 MCT=-1                                                            
      Z1A3=Y1A3-DZ1A3                                                   
      IF(Z1A3.LT.-X1IH)Z1A3=Z1A3+X1I                                     
      GO TO 316                                                         
  356 IF(JE.LT.KE)GO TO 354                                             
      GO TO 359                                                          
  358 CONTINUE                                                          
      Y1A3=Z1A3                                                         
      Y2A3=Z2A3                                                         
!1359 WRITE(NQ,1358)Z1A3,Z2A3
!1358 FORMAT(' Z1A3 =',F8.3,'  Z2A3 =',F8.3)
!-----Has negative standoff been calculated?                             
  359 IF(JCT.EQ.1)GO TO 362                                              
  360 X1A3=Y1A3                                                          
      X2A3=Y2A3                                                          
!     Set up for negative standoff calculation                           
      JCT=1                                                              
      Z2A3=-H2N                                                          
      HL=-H2M                                                               
      Y2A3=-100.                                                            
      Y1A3=H1A3-H21I                                                     
      GO TO 312                                                             
  362 X21I=X1A3-Y1A3                                                     
      X2I=X2A3-Y2A3                                                         
      NARK=0
      VY=X1IH*X2I*H3N                                                       
      IF(IPR.LT.4)GO TO 365                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,364)A3,ALPHA(LCT),X1I,X2I,H3M,VY,X1A3,X2A3,H3A3,X1P3,X3P3 &
     &,X21I                                                                 
  364 FORMAT(1X,F7.1,2X,A2,3F7.3,F7.2,6F7.3)                             
  365 IF(VY.GE.VX)GO TO 392                                                 
!-----Check upper layers on positive axis-2 side                         
      Z1A3T=2.*X1A3-Y1A3
      Z3CT=H3A3+X3P3
	  N=MARK
      DO 379 I=1,NMOD
      DO 378 J=1,NMOD
	  K=IA(I)+IAA(J)
	  L=K+10
	  ICT=0
	  D1=Z1A3T-W(1,J)-W(1,I)
	  D2=-Y2A3-W(2,J)-W(2,I)
	  D3=H3A3+W(3,J)-W(3,I)
  371 D1=D1+X1I
      IF(D1.LT.CN(K))GO TO 371
  372 D1=D1-X1I
      IF(D1.GT.CN(K))GO TO 372
	  IF(D1.LT.-CN(K))GO TO 375
  373 D3=D3+H3M
      IF(D3.LT.CN(K))GO TO 373
  374 D3=D3-H3M
      IF(D3.GT.CN(K))GO TO 374
      IF(D3.LT.-CN(K))GO TO 372
	  D=D1**2+D3**2
	  IF(D.GT.CN(L))GO TO 374
	  IT(N)=K
	  TI(N+1)=D
	  TI(N+2)=D2
	  N=N+NSTP
	  IF(N.GT.LIMIT)GO TO 902
	  GO TO 374
!-----Have C molecule distances been calculated                          
  375 IF(ICT.GT.0)GO TO 378
      ICT=1
 	  D1=D1+X1P3
	  D3=Z3CT-W(3,J)-W(3,I)
	  GO TO 371
  378 CONTINUE
  379 CONTINUE
      NARK=N
	  DO 390 I=1,NMOD
	  DO 389 J=1,NMOD
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      ICT=1                                                                
!-----AXIS CONTACT                                                          
      D1=X1A3-W(1,J)-W(1,I)                                                 
      D2=-W(2,J)-W(2,I)                                                  
      D3=H3A3+W(3,J)-W(3,I)                                                 
  381 D1=D1+X1I                                                             
      IF(D1.LT.CN(K))GO TO 381                                              
  382 D1=D1-X1I                                                             
      IF(D1.GT.CN(K))GO TO 382                                              
      IF(D1.LT.-CN(K))GO TO 385                                             
  383 D3=D3+H3M                                                             
      IF(D3.LT.CN(K))GO TO 383                                              
  384 D3=D3-H3M                                                             
      IF(D3.GT.CN(K))GO TO 384                                              
      IF(D3.LT.-CN(K))GO TO 382                                             
      D=D1**2+D3**2                                                         
      IF(D.GT.CN(L))GO TO 384                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D2                                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 384                                                             
!-----Have OTHER molecule contacts been calculated?                      
  385 GO TO(386,387,388,389),ICT                                         
!-----PLANE CONTACT                                                         
  386 D1=X1P3+X1A3-Y1A3+W(1,J)-W(1,I)                                       
      D3=X3P3-W(3,J)-W(3,I)                                              
      ICT=2                                                              
      D2=-Y2A3+W(2,J)-W(2,I)                                             
      GO TO 381                                                             
!-----IDENTITY CONTACT                                                      
  387 D1=D1-X1P3                                                            
      D3=W(3,J)-W(3,I)                                                      
      ICT=3                                                              
      GO TO 381                                                             
!-----CENTER CONTACT                                                        
  388 D1=X1A3+X1P3-W(1,J)-W(1,I)                                            
      D2=-W(2,J)-W(2,I)                                                     
      D3=H3A3+X3P3-W(3,J)-W(3,I)                                            
      ICT=4                                                              
      GO TO 381                                                             
  389 CONTINUE                                                              
  390 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      Z2A3=X2A3
      CALL MINHI(Z2A3,X2A3,100.,H2M,KE)                                     
      NARK=0
      IF(KCT)391,1390,906
 1390 CONTINUE
!	WRITE(NQ,2390)X2A3,Z2A3
!2390 FORMAT(' X2A3 =',F8.3,'   Z2A3 =',F8.3)
      X2A3=Z2A3
      X2I=X2A3-Y2A3
      VY=X1IH*X2I*H3N
	  IF(VY.GE.VX)GO TO 392
  391 VX=VY                                                                 
      H1I=X1I                                                            
      H1IH=X1IH                                                          
      H21I=X21I                                                          
      H2I=X2I                                                            
      H1P3=X1P3                                                          
      H3P3=X3P3                                                          
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
      IF(NCT)395,393,395                                                 
  392 IF(NCT)395,393,394                                                 
  393 IF(LE.GE.NV)GO TO 395                                                 
      LE=LE+1                                                               
      X3P3=X3P3+DX3P3                                                       
      GO TO 306                                                             
! 394 X3P3=H3P3-X3P3                                                        
  394 X3P3=H3P3-DX3P3                                                    
      NCT=-1                                                                
      GO TO 306                                                             
  395 IF(LE.GE.KE)GO TO 397                                                 
  396 NCT=1                                                                 
      LE=2*LE                                                               
      DX3P3=H3M/LE                                                          
      X3P3=H3P3+DX3P3                                                       
      GO TO 306                                                             
  397 IF(KE.GE.NNE)GO TO 398                                                
      KE=NNE                                                                
      GO TO 396                                                             
  398 IF(VX.GE.V3)GO TO 400                                                 
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(LCT)                                               
      IF(V3.GE.HT1(8))GO TO 399                                             
      H1C=H1A3+H1P3                                                      
      H2C=H2A3                                                           
      H3C=H3A3+H3P3                                                      
!     IF(H1C.GT.H1IH)H1C=H1C-H1I                                         
      IF(H3C.GT.H3N)H3C=H3C-H3M                                          
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3M                                                            
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
      NLOW=9                                                                
      MLOW=LCT+9                                                         
  399 IF(IPR.LT.2)GO TO 400                                              
      LCTM=LCT                                                           
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3M                                                           
!     Cosines of cell angles                                             
      C(1)=0.                                                            
      C(2)=0.                                                            
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates  
!-----Transform first to Angstrom coordinates                            
      SINE3=SQRT(1.-C(3)**2)                                             
      A(6)=H3C                                                           
      A(5)=H2C/SINE3                                                     
      A(4)=H1C-A(5)*C(3)                                                 
      A(9)=H3P3                                                          
      A(8)=0.                                                            
      A(7)=H1P3                                                          
      A(12)=H3A3                                                         
      A(11)=H2A3/SINE3                                                   
      A(10)=H1A3-A(11)*C(3)                                              
!-----Then to fractional coordinates                                     
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      A(7)=A(7)/A(1)                                                     
      A(8)=A(8)/A(2)                                                     
      A(9)=A(9)/A(3)                                                     
      A(10)=A(10)/A(1)                                                   
      A(11)=A(11)/A(2)                                                   
      A(12)=A(12)/A(3)                                                   
      C(4)=90.                                                           
      C(5)=90.                                                           
      C(6)=57.296*ACOS(C(3))                                             
  400 IF(LCT.NE.3)GO TO 401                                                 
      IF(ITR(13).EQ.1)GO TO 301                                             
  401 IF(IPR.LT.2)GO TO 1000                                             
      IF(V3.GE.V3S)GO TO 1000                                            
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,405)ALPHA(LCTM),A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)     
  405 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3   &                  
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
  406 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)
      SM(1)=A(4)-A(4)
      SM(2)=A(5)-A(5)
      SM(3)=A(6)-A(6)
      WRITE(NQ,406)SYM(1),(A(I),I=4,6),SM
      SM(1)=A(10)-A(4)
      SM(2)=A(11)-A(5)
      SM(3)=A(12)
      WRITE(NQ,406)SYM(2),(A(I),I=10,12),SM
      SM(1)=A(7)
      SM(2)=A(8)
      SM(3)=A(9)-A(6)
      WRITE(NQ,406)SYM(3),(A(I),I=7,9),SM
      GO TO 1000
  902 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZFRTA,5X,3(I9,I3),I9)                  
      GO TO 999                                                             
 1904 NERR=151
      GO TO 904
 2904 NERR=251
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)NERR                                                  
  905 FORMAT(' PARAMETER TOO SMALL NEAR',I4)                             
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT('0PARAMETER TOO LARGE')                                     
      GO TO 999                                                          
  920 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,921)KCT,MCT,NCT,MARK1,MARK2,MARK,NEND,NARK,N                 
  921 FORMAT(16H0NO INTERACTIONS,3I3,6I6)                                   
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZFRTA                              
