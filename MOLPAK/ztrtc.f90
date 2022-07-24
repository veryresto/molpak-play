      SUBROUTINE ZTRTC  
!                                                        
!-----This subroutine finds structures in the P1bar space group with I   
!     molecules in contact along all three triclinic axes.  First, the   
!     axis-1 I separation is determined by Procedure 1.  Then an I line  
!     offset in the axis-2 direction is determined (with one I molecule  
!     arbitrarily placed on axis-2) by a shortened version of Procedure  
!     2.  That is, the offset of the I molecule along the axis-1         
!     direction is determined by iteration, not by Procedure 2.  The     
!     rest of the "iteration loop" consists of approach of a layer of C  
!     molecules from both the positive and negative axis-3 directions by 
!     Procedure 3.  A vector between a negative C and a positive C       
!     places the remainder of the I molecules.  If one or more I         
!     molecules are too close, both the I and C molecules are backed off.
!     Iteration continues to find the I offset from axis-2 which produces
!     the minimum cell volume. For CA. 
!                                                                        
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
!-----Linear storage array and markers                                   
!     COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000),
!     &       CODE(5000)
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
!
      CHARACTER(2) ::  ALPHA ='CA' 
! 
      INTEGER :: I, ICT, IE, J, JCT, JE, K, KE, L, LCT, LE
      INTEGER :: N, NARC, NCT, NERR, MCT
!     
      REAL, DIMENSION (6) ::  A, C
      REAL :: D, D1, D2, D3, DX21I, DZ1C, DZ2C
      REAL :: E, F, FXX, FXY
      REAL :: H1C, H1I, H1IH, H1IT
      REAL :: H21I, H2C, H2I, H3C 
      REAL :: H31I, H32I, H3I, HL
      REAL :: SINE2, SINE3
      REAL :: X1C, X2C, X2I, X2IH, X21I, X21IT
      REAL :: X3C, X31I, X32I, X3I
      REAL :: Y1C, Y2C, Y3C, YY1C, YY2C
      REAL :: VX, VY, WM
      REAL :: Z1C, Z2C, Z3C, Z1CF, Z2CF
!

      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)ITR(9)                                                  
    5 FORMAT(1X,'ZTRTC called',I3)                                       
!-----Set dummy value of the minimum cell volume                         
   10 VX=1000000.                                                        
!-----Calculate axis-1 length                                            
      MARK=MARK2                                                            
      N=MARK                                                             
      DO 20 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 20                                                 
      IF(IT(I+7).NE.0)GO TO 20                                           
      D2=TI(I+3)                                                            
      IF(ABS(D2).GT.CN(K))GO TO 20                                          
      D=TI(I+1)+D2**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 20                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+2)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
   20 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                               
      ER=ERM                                                                
      H1I=H1M                                                                
      CALL MINHI(H1I,0.,100.,H1M,NNE)                                       
      IF(KCT)1904,21,906                                                 
   21 H1IH=H1I/2.                                                        
      H1IT=2.*H1I                                                        
      IF(IPR.LT.3)GO TO 23                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,22)H1I                                                    
   22 FORMAT(1X,'H1I =',F6.3)                                            
!-----Top of axis iteration loop                                         
   23 X2I=H2M                                                            
      KE=NE                                                                 
      IE=1                                                                  
      DX21I=H1I/NV                                                          
      X21I=-H1IH                                                            
      LCT=0                                                                 
   25 ER=ERMT                                                            
      N=MARK                                                                
!-----Collect second line for double standoff                            
      X21IT=2.*X21I                                                       
      DO 30 I=MARK1,NEND1,NSTP1	
	K=IT(I)					    
	IF(K.LE.0)GO TO 30		    
	IF(IT(I+7).NE.0)GO TO 30	    
	D1=TI(I+2)+X21IT		    
   28 D1=D1+H1I				    
      IF(D1.LE.CN(K))GO TO 28		    
   29 D1=D1-H1I					    
      IF(D1.GT.CN(K))GO TO 29	    
	IF(D1.LT.-CN(K))GO TO 30	    
	D=TI(I+1)+D1**2			    
      L=IT(I+4)				    
	IF(D.GT.CN(L))GO TO 29		    
	IT(N)=K				    
	TI(N+1)=D				    
	TI(N+2)=TI(I+3)			    
	N=N+NSTP				    
	IF(N.GT.LIMIT)GO TO 902			    
	GO TO 29			    
   30 CONTINUE			
      NARK=N					    
      DO 33 I=MARK1,NEND1,NSTP1                                             
      K=IT(I)                                                               
      IF(K.LE.0)GO TO 33                                                 
      IF(IT(I+7).NE.0)GO TO 33                                           
      D1=TI(I+2)+X21I                                                       
   31 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 31                                               
   32 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 32                                               
      IF(D1.LT.-CN(K))GO TO 33                                              
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 32                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 32                                                              
   33 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(X2I,0.,100.,H2M,KE)                                        
      NARK=0					    
      IF(KCT)2904,35,906                                                 
!-----Calculate closest approach of C molecules along axis-3             
!     Setup for approach from positive axis-3 direction                  
   35 Z3C=H3M                                                               
      Y3C=100.                                                              
      HL=H3M                                                             
      ER=ERMF                                                               
      JCT=0                                                                 
      NARK=0                                                              
      X2IH=X2I/2.                                                         
	N=MARK                                                              
      IF(IPR.LT.3)GO TO 40                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,36)X2I                                                    
   36 FORMAT(1X,'X2I =',F7.3)                                            
   40 CONTINUE	
      JE=1                                                                  
      LE=1                                                                  
      DZ1C=H1I/NV                                                           
      DZ2C=X2I/NV                                                           
      Z2C=-X2IH                                                             
      Z1C=-H1IH                                                             
      NCT=1                                                                 
      MCT=1                                                                 
   42 N=MARK                                                                
      NARK=0					    
      IF(JCT.LE.0)GO TO 242	    
!-----Collect distances for second layer of C molecules                  
!     WRITE(NQ,41)		    
!  41 FORMAT(' Start negative standoff')
      Z1CF=2.*Z1C-X1C			    
	Z2CF=2.*Z2C-X2C			    
      DO 109 I=1,NMOD			    
	DO 108 J=1,NMOD			    
	K=IA(I)+IAA(J)			    
	L=K+10					    
	D1=Z1CF-W(1,J)-W(1,I)	    
	D2=Z2CF-W(2,J)-W(2,I)	    
  101 D2=D2+X2I				    
      D1=D1+X21I					    
	IF(D2.LT.CN(K))GO TO 101 	    
  102 D2=D2-X2I					     
      D1=D1-X21I					    
	IF(D2.GT.CN(K))GO TO 102    	
  	IF(D2.LT.-CN(K))GO TO 108	
  103 D1=D1+H1I				    
      IF(D1.LT.CN(K))GO TO 103	
  104 D1=D1-H1I				    
      IF(D1.GT.CN(K))GO TO 104	
	IF(D1.LT.-CN(K))GO TO 102  
      D=D1**2+D2**2			    
	IF(D.GT.CN(L))GO TO 104		     
	IT(N)=K						    
	TI(N+1)=D				    
	TI(N+2)=-X3C-W(3,J)-W(3,I)   
	N=N+NSTP					    
	IF(N.GE.LIMIT)GO TO 902		    
	GO TO 104				    
  108 CONTINUE				    
  109 CONTINUE	 			    
      NARK=N					    
  242 DO 250 I=1,NMOD                                                       
      DO 249 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      ICT=JCT                                                             
      D1=Z1C-W(1,J)-W(1,I)                                                  
      D2=Z2C-W(2,J)-W(2,I)                                                  
      D3=0.                                                               
      WM=-1.                                                              
  243 D2=D2+X2I                                                             
      D1=D1+X21I                                                            
      IF(D2.LE.CN(K))GO TO 243                                              
  244 D2=D2-X2I                                                             
      D1=D1-X21I                                                            
      IF(D2.GT.CN(K))GO TO 244                                              
      IF(D2.LT.-CN(K))GO TO 247                                             
  245 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 245                                              
  246 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 246                                              
      IF(D1.LT.-CN(K))GO TO 244                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 246                                               
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=D3+WM*W(3,J)-W(3,I)                                            
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 246                                                             
  247 IF(ICT.LE.0)GO TO 249                                               
!-----Add I molecule distances on second standoff only                   
      D1=Z1C-X1C+W(1,J)-W(1,I)                                            
	D2=Z2C-X2C+W(2,J)-W(2,I)                                            
	D3=-X3C                                                             
	ICT=0                                                               
      WM=1.                                                               
	GO TO 243                                                           
  249 CONTINUE                                                              
  250 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                              
      CALL MINHI(Z3C,0.,Y3C,HL,KE)                                          
      NARK=0					    
!     IF(KCT)3904,352,351                                                
      IF(KCT)368,352,351		    
  351 IF(MCT.LE.1)GO TO 354                                              
      J=MCT-NCT                                                          
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 367                                 
      GO TO 353                                                          
  352 Y1C=Z1C                                                            
      Y2C=Z2C                                                            
      Y3C=Z3C                                                            
      NCT=MCT                                                            
  353 GO TO(354,359,360,361,362,363,364,365,358),MCT                     
  354 IF(LE.GE.NV)GO TO 355                                              
      LE=LE+1                                                            
      Z2C=Z2C+DZ2C                                                       
      GO TO 42                                                           
  355 IF(JE.GE.NV)GO TO 356                                              
      JE=JE+1                                                            
      Z1C=Z1C+DZ1C                                                       
      LE=1                                                               
      Z2C=-X2IH                                                          
      GO TO 42                                                           
  356 NCT=6                                                              
  357 JE=2.*JE                                                           
      DZ1C=H1I/JE                                                        
      DZ2C=X2I/JE                                                        
      YY1C=Y1C                                                           
      YY2C=Y2C                                                           
      NCT=NCT+1                                                          
      IF(NCT.GT.9)NCT=NCT-8                                              
      GO TO(906,360,361,362,363,364,365,358,359),NCT                     
  358 MCT=2                                                              
      Z1C=YY1C+DZ1C                                                      
      Z2C=YY2C                                                           
      GO TO 42                                                           
  359 MCT=3                                                              
      Z1C=YY1C+DZ1C                                                      
      Z2C=YY2C+DZ2C                                                      
      GO TO 42                                                           
  360 MCT=4                                                              
      Z1C=YY1C                                                           
      Z2C=YY2C+DZ2C                                                      
      GO TO 42                                                           
  361 MCT=5                                                              
      Z1C=YY1C-DZ1C                                                      
      IF(Z1C.LT.-H1IH)Z1C=Z1C+H1I                                        
      Z2C=YY2C+DZ2C                                                      
      GO TO 42                                                           
  362 MCT=6                                                              
      Z1C=YY1C-DZ1C                                                      
      IF(Z1C.LT.-H1IH)Z1C=Z1C+H1I                                        
      Z2C=YY2C                                                           
      GO TO 42                                                           
  363 MCT=7                                                              
      Z1C=YY1C-DZ1C                                                      
      IF(Z1C.LT.-H1IH)Z1C=Z1C+H1I                                        
      Z2C=YY2C-DZ2C                                                      
      IF(Z2C.LT.-X2IH)Z2C=Z2C+X2I                                        
      GO TO 42                                                           
  364 MCT=8                                                              
      Z1C=YY1C                                                           
      Z2C=YY2C-DZ2C                                                      
      IF(Z2C.LT.-X2IH)Z2C=Z2C+X2I                                        
      GO TO 42                                                           
  365 MCT=9                                                              
      Z1C=YY1C+DZ1C                                                      
      Z2C=YY2C-DZ2C                                                      
      IF(Z2C.LT.-X2IH)Z2C=Z2C+X2I                                        
      GO TO 42                                                           
  367 IF(JE.LT.KE)GO TO 357                                              
      GO TO 400                                                           
  368 CONTINUE                                                           
!-----Stop if zero offset found after first positive determination       
!     IF(JCT.GT.0)GO TO 3904                                              
      Y1C=Z1C                                                             
	Y2C=Z2C                                                             
	Y3C=Z3C                                                             
!-----Have both directions been fully calculated?                        
  400 IF(JCT-1)410,420,475                                                
  410 X3C=Y3C                                                               
      X2C=Y2C                                                            
      X1C=Y1C                                                            
!-----Setup for approach from negative axis-3 direction                  
      Y3C=-100.                                                             
      Z3C=-H3M                                                              
      HL=-H3M                                                               
      JCT=1                                                                 
      GO TO 40                                                              
!-----Setup for second approach from positive axis-3 direction           
  420 X1C=Y1C					    
      X2C=Y2C					    
	X3C=Y3C					    
	Y3C=100.				    
	Z3C=H3M					    
	HL=H3M						    
	JCT=2					    
	GO TO 40				    
!-----Both directions along axis-3 have been completed                   
  475 IF(IPR.LT.3)GO TO 477                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,476)X1C,X2C,X3C,Y1C,Y2C,Y3C                               
  476 FORMAT(1X,'X and YC =',3F8.3,3X,3F8.3)                             
  477 X3I=Y3C-X3C                                                            
      X32I=Y2C-X2C                                                       
      X31I=Y1C-X1C                                                       
  491 VY=H1I*X2I*X3I/2.                                                     
      IF(IPR.LT.3)GO TO 494                                                 
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,492)X21I,X2I,X3I,D                                           
  492 FORMAT(7H X21I =,F7.3,7H, X2I =,F7.3,7H, X3I =,F7.3,5H, D =,F8.2)     
!-----Is this the smallest cell volume per molecule found so far?        
  494 IF(VY.GT.VX)GO TO 495                                              
!-----Yes, this is the smallest                                          
      VX=VY                                                              
      H2I=X2I                                                            
      H21I=X21I                                                          
      H3I=X3I                                                            
      H31I=X31I                                                          
      H32I=X32I                                                          
      H1C=X1C                                                            
      H2C=X2C                                                            
      H3C=X3C                                                            
      IF(LCT)498,496,498                                                 
!-----No, this is not the smallest                                       
  495 IF(LCT)498,496,497                                                 
  496 IF(IE.GE.NV)GO TO 498                                                 
      IE=IE+1                                                               
      X21I=X21I+DX21I                                                       
      GO TO 25                                                              
  497 LCT=-1                                                             
      X21I=H21I-DX21I                                                    
      IF(X21I.LT.-H1IH)X21I=X21I+H1I                           
      GO TO 25                                                              
  498 IF(IE.GE.KE)GO TO 580                                                 
  499 IE=2*IE                                                               
      DX21I=H1I/IE                                                          
      X21I=H21I+DX21I                                                       
      LCT=1                                                                 
      GO TO 25                                                              
  580 IF(KE.GE.NNE)GO TO 590                                             
!-----Shift to maximum accuracy                                          
      KE=NNE                                                             
      GO TO 499                                                          
!-----Is this the smallest volume found at this molecular rotation?      
  590 IF(VX.GE.V3)GO TO 600                                              
      V3=VX                                                              
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA                                                    
!-----Is this the smallest volume found in this SEEK operation?          
      IF(V3.GE.HT1(8))GO TO 600                                             
      HT1(1)=H1I                                                            
      HT1(2)=H21I                                                           
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=H1C                                                            
      HT1(6)=H2C                                                            
      HT1(7)=H3C                                                            
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA                                                         
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=H31I                                                          
      HT1(14)=H32I                                                          
      NLOW=7                                                                
      MLOW=9                                                             
  600 IF(IPR.LT.2)GO TO 1000                                             
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=SQRT(H3I**2+H31I**2+H32I**2)                                  
!     Cosines of cell angles                                             
      C(1)=(H21I*H31I+H2I*H32I)/(A(2)*A(3))                              
      C(2)=H31I/A(3)                                                     
      C(3)=H21I/A(2)                                                     
!-----Transform molecular positions from orthogonal to cell coordinates  
      SINE2=SQRT(1.-C(2)**2)                                             
      SINE3=SQRT(1.-C(3)**2)                                             
      E=(C(2)*C(3)-C(1))/(SINE2*SINE3)                                      
      F=SQRT(1.-E**2)                                                    
      FXX=SINE2*F                                                           
      FXY=SINE2*E/SINE3                                                     
!-----Transform first to Angstrom coordinates                            
      A(6)=H3C/FXX                                                       
      A(5)=H2C/SINE3+A(6)*FXY                                            
      A(4)=H1C-A(6)*C(2)-A(5)*C(3)                                       
!-----Then to fractional coordinates                                     
      A(4)=A(4)/A(1)                                                     
      A(5)=A(5)/A(2)                                                     
      A(6)=A(6)/A(3)                                                     
      C(4)=57.296*ACOS(C(1))                                             
      C(5)=57.296*ACOS(C(2))                                             
      C(6)=57.296*ACOS(C(3))                                             
      D=SQRT(H1C**2+H2C**2+H3C**2)                                       
      CALL PAGE(2,2,0)                                                   
      WRITE(NQ,605)ALPHA,A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)           
  605 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3   &                  
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
      WRITE(NQ,606)H1C,H2C,H3C,(A(I),I=4,6),D                            
  606 FORMAT(1X,'C AT',3F7.3,';  ',3F7.4,';  D =',F6.2)                  
      GO TO 1000                                                         
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZTRTC,5X,3(I9,I3),I9)                  
      GO TO 999                                                          
 1904 NERR=21					    
      GO TO 904						    
 2904 NERR=35					    
  904 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,905)NERR,JCT,H1I,X21I,X2I                                 
  905 FORMAT(' PARAMETER TOO SMALL NEAR',I4,I2,3F8.3)                    
      GO TO 999                                                          
  906 CALL PAGE(2,2,0)                                                   
      WRITE(NQ,907)                                                      
  907 FORMAT(20H0PARAMETER TOO LARGE)                                    
      GO TO 999                                                          
  920 WRITE(NQ,921)                                                      
  921 FORMAT(21H0NO ATOM INTERACTIONS)                                   
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZTRTC           

