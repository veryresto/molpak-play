      SUBROUTINE ZFRTB   
!
!-----This subroutine finds structures in space group P21/c with both I  
!     and P molecules in contact along the axis-1 direction.             
!     The AM coordination sphere routine (like all others) assumes that  
!     the structure collapses along the three orthogonal axes until      
!     stopped by intermolecular replusion.  Axis-1, 2, and 3 are treated 
!     separately in order.  Identities, if present, are used first;      
!     therefore, axis-1 (monoclinic a) is determined by I contact.  Axis-
!     2 is determined by a row of P molecules along the axis-1 direction.
!     Their offset along the axis-3 direction (monoclinic b) is unknown; 
!     therefore it is arbitrarily set and (later) determined by          
!     iteration.  Monoclinic á for each iterative step is determined by  
!     sliding the P row along axis-1 to the point which gives a minimum  
!     axis-2 half-length (Procedure 2).  The axis-3 half-length for each 
!     iteration is determined by lowering a configuration of A and !     
!     molecules.  The shape of the configuration is specified by axis-1  
!     and the position of the P molecules.  The vector between an A and  
!     C molecule is equal to the vector between the central I and a P    
!     molecule. For AM.                                                         
!
!
      USE molpakCommonMod
!
      IMPLICIT NONE
!
!
!     CHARACTER*80 BUF                                                       
!      CHARACTER*60 HEAD                                                  
!      CHARACTER*4 AN,CDIJ,CHT1                                           
!      CHARACTER*2 AT,CODE                                                 
!-----General program communication and transfer variables               
!     COMMON BUF,KILL,NPAGE,LINE,NEW,JTR(30),ITR(90)                     
!-----Energy constants                                                       
!      COMMON NSEP,CN(1200),ERM,FERM,ER,ERMT,ERMF,NV,NE,NNE                   
!-----Compound input data                                                
!     COMMON NMOD,NCTR,IA(200),IAA(200),X(3,200),AT(200),AN(200)           
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
      CHARACTER(2) :: ALPHA ='AM'
      CHARACTER(2) :: SYM(3) = (/'C ','A3','P3'/)
!      DATA ALPHA/'AM'/                                                      
!      DATA SYM/'C ','A3','P3'/                                           
!
      INTEGER :: I, ICT, IE, J, JCT, JE, K, KE, L, LE
      INTEGER :: N, NCT, MCT
!
      REAL :: A(12),C(6)
      REAL :: D, D1, D2, D3, DX3P3
      REAL :: DY1A3, DY2A3, DY1P3
      REAL :: H1C, H2C, H3C, H1I, H2I, H3I, H21I, H1IH, H1IT 
      REAL :: H1A3, H2A3, H3A3, H1P3, H2P3, H3P3, H3P3R, H3P3S
      REAL :: SM(3), SINE3
      REAL :: VX, VY
      REAL :: X1A3, X1P3, X2A3, X2I, X21I, X2IH,  X2P3
      REAL :: X3A3, X3P3, XX1A3, XX2A3
      REAL :: Y1A3, Y1P3, Y2A3, Y2P3, Y3A3  
      REAL :: Z3A3
!
      NARK=0                                                             
      IF(IPR.LT.2)GO TO 10                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)ITR(14)                                                 
    5 FORMAT(1X,'ZFRTB called',I3)                                       
!-----Set dummy value of unit cell volume                                
   10 VX=1000000.                                                        
!     Set to intermediate accuracy                                       
      KE=NE                                                              
!-----Separate I contacts with both short axis-3 and axis-2 distances    
      NSTP2=5                                                            
      N=MARK2                                                            
      DO 15 I=MARK1,NEND1,NSTP1                                          
      K=IT(I)                                                            
      IF(K.LE.0)GO TO 15                                                 
      IF(IT(I+7).NE.0)GO TO 15                                           
      IT(N)=K                                                               
      TI(N+1)=TI(I+1)                                                       
      TI(N+2)=TI(I+2)                                                       
      TI(N+3)=TI(I+3)                                                    
      IT(N+4)=IT(I+4)                                                    
      N=N+NSTP2                                                          
      IF(N.GT.LIMIT)GO TO 902                                            
   15 CONTINUE                                                           
      IF(N.LE.MARK2)GO TO 920                                            
      NEND2=N-1                                                          
      MARK3=N                                                            
!-----Calculate axis-1 length                                            
      MARK=N                                                             
      DO 20 I=MARK2,NEND2,NSTP2                                             
      K=IT(I)                                                               
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
!-----Make sure some interactions were found                             
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                               
      ER=ERM                                                                
      H1I=H1M                                                            
   21 CALL MINHI(H1I,0.,100.,H1M,NNE)                                       
      IF(KCT)904,22,906                                                  
!-----Assume that the range of possible P3 offsets is from -H3N to H3N   
   22 H1IT=2.*H1I                                                        
      H1IH=H1I/2.                                                        
      IF(IPR.LT.3)GO TO 25                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,23)H1I                                                    
   23 FORMAT(1X,'H1I =',F8.3)                                            
   25 H3P3S=-H3N                                                         
      H3P3R=H3M                                                          
!-----INITIALIZE PLANE OFFSET =LOOP=                                        
   30 NSTP3=NSTP2                                                        
      JCT=0                                                                 
      IE=1                                                                  
      DX3P3=H3P3R/NV                                                        
      X3P3=H3P3S                                                            
!-----TOP OF PLANE OFFSET =LOOP=                                            
!-----FIND POSSIBLE PLANE CONTACTS AT THIS OFFSET                           
   31 N=MARK3                                                               
      DO 33 I=1,NMOD                                                        
      DO 32 J=1,NMOD                                                        
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D3=X3P3-W(3,J)-W(3,I)                                                 
      IF(ABS(D3).GT.CN(K))GO TO 32                                          
      IT(N)=K                                                               
      TI(N+1)=D3**2                                                         
      TI(N+2)=W(1,J)-W(1,I)                                                 
      TI(N+3)=W(2,J)-W(2,I)                                                 
      IT(N+4)=L                                                             
      N=N+NSTP3                                                             
      IF(N.GT.LIMIT)GO TO 902                                               
   32 CONTINUE                                                              
   33 CONTINUE                                                              
      IF(N.LE.MARK3)GO TO 920                                               
      NEND3=N-1                                                             
!-----INITIALIZE SEARCH FOR MIN. PLANE-1,2 AREA USING IDENT AND PLANE       
!        CONTACTS                                                           
      MARK=N                                                                
      ER=ERMT                                                               
      X2P3=100.                                                             
      Y2P3=H2N                                                           
      DY1P3=H1I/NV                                                       
      MCT=0                                                                 
!-----Skip scan of range of possible Y1P3 values if in final stages of   
!     refinement                                                         
      IF(IE.LE.NE)GO TO 35                                               
      Y1P3=H1P3                                                          
      JE=IE/2                                                            
      GO TO 40                                                           
!-----Initial stages, make scan                                          
   35 Y1P3=-H1IH                                                         
      JE=1                                                                  
!-----TOP OF BETA FINDING =LOOP=                                            
   40 N=MARK                                                                
!     Pick up both I and P molecule contacts - start at MARK2            
      DO 45 I=MARK2,NEND3,NSTP2                                             
      K=IT(I)                                                               
      D1=TI(I+2)+Y1P3                                                       
      IF(I-MARK3)41,42,43                                                   
!-----IDENT CONTACT - twice as far along axis-1                          
   41 D1=D1+Y1P3                                                            
      GO TO 43                                                              
!-----SET NARK AS SIGNAL TO MINHI OF SHIFT TO PLANE CONTACTS                
   42 NARK=N                                                                
   43 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 43                                               
   44 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 44                                               
      IF(D1.LT.-CN(K))GO TO 45                                              
      D=TI(I+1)+D1**2                                                       
      L=IT(I+4)                                                             
      IF(D.GT.CN(L))GO TO 44                                                
      IT(N)=K                                                               
      TI(N+1)=D                                                             
      TI(N+2)=TI(I+3)                                                       
      N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 44                                                              
   45 CONTINUE                                                              
      IF(N.LE.MARK)GO TO 920                                                
      NEND=N-1                                                              
      CALL MINHI(Y2P3,0.,X2P3,H2M,KE)                                        
      IF(KCT)904,51,52                                                   
   51 X1P3=Y1P3                                                          
      X2P3=Y2P3                                                          
      IF(MCT)56,53,56                                                    
   52 IF(MCT)56,53,55                                                    
   53 IF(JE.GE.NV)GO TO 54                                               
      JE=JE+1                                                            
      Y1P3=Y1P3+DY1P3                                                    
      GO TO 40                                                           
   54 MCT=1                                                              
      JE=2.*JE                                                           
      DY1P3=H1I/JE                                                       
      Y1P3=X1P3+DY1P3                                                    
      GO TO 40                                                           
   55 MCT=-1                                                             
      Y1P3=X1P3-DY1P3                                                    
      IF(Y1P3.LT.-H1IH)Y1P3=Y1P3+H1I                                     
      GO TO 40                                                           
   56 IF(JE.LT.KE)GO TO 54                                               
   58 CONTINUE                                                           
      IF(IPR.LT.3)GO TO 60                                                  
      CALL PAGE(1,1,0)                                                      
      WRITE(NQ,59)X1P3,X2P3,X3P3                                         
   59 FORMAT(1X,'X1P3 =',F8.3,', X2P3 =',F8.3,', X3P3 =',F8.3)              
   60 MARK=MARK3                                                            
!-----Find axis-3 length by placing a grid of A and C molecules          
!        INITIALIZE                                                         
      NARK=0                                                                
      X2I=2.*X2P3                                                           
      X2IH=X2P3                                                          
      X21I=2.*X1P3                                                          
!-----FOUR TIMES ER BECAUSE INCLUDED PLANE CONTACTS ARE RARE             
      ER=ERMF                                                            
      H1C=0.                                                                
      H2C=0.                                                                
      X3A3=100.                                                             
      Y3A3=H3N                                                              
      MCT=1                                                              
      NCT=1                                                              
      DY1A3=H1I/NV                                                          
      DY2A3=X2I/NV                                                          
!-----Skip scan of Y1A3, Y2A3 range in final stages of refinement        
      IF(IE.LE.KE)GO TO 61                                               
      Y1A3=H1A3                                                          
      Y2A3=H2A3                                                          
      JE=IE/2                                                            
      LE=IE/2                                                            
      GO TO 62                                                           
!-----Initial stages of refinement                                       
   61 Y1A3=-H1IH                                                            
      Y2A3=-X2IH                                                            
      JE=1                                                                  
      LE=1                                                                  
!-----TOP OF AXIS-CENTER PLACING =LOOP=                                     
!-----FIND PLANE CONTACTS - DO NOT CHANGE WITH AXIS PLACEMENT            
   62 N=MARK                                                             
      DO 69 I=1,NMOD                                                     
      DO 68 J=1,NMOD                                                     
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D2=X2P3+W(2,J)-W(2,I)                                              
      IF(ABS(D2).GT.CN(K)) GO TO 68                                      
      D1=X1P3+W(1,J)-W(1,I)                                              
   63 D1=D1+H1I                                                          
      IF(D1.LE.CN(K)) GO TO 63                                           
   64 D1=D1-H1I                                                          
      IF(D1.GT.CN(K)) GO TO 64                                           
      IF(D1.LT.-CN(K)) GO TO 68                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L)) GO TO 64                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=X3P3-W(3,J)-W(3,I)                                         
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 64                                                           
   68 CONTINUE                                                           
   69 CONTINUE                                                           
!-----Set signal for MINHI to switch from P contacts for which the       
!     standoff is twice that for A contacts.                             
      NARK=N                                                             
  116 N=NARK                                                             
      DO 129 I=1,NMOD                                                       
      DO 128 J=1,NMOD                                                       
      K=IA(I)+IAA(J)                                                        
      L=K+10                                                                
      D2=Y2A3-W(2,J)-W(2,I)                                                 
      D1=Y1A3-W(1,J)-W(1,I)                                                 
      ICT=0                                                                 
  121 D2=D2+X2I                                                             
      D1=D1+X21I                                                            
      IF(D2.LE.CN(K))GO TO 121                                              
  122 D2=D2-X2I                                                             
      D1=D1-X21I                                                            
      IF(D2.GT.CN(K))GO TO 122                                              
      IF(D2.LT.-CN(K))GO TO 127                                             
  123 D1=D1+H1I                                                             
      IF(D1.LE.CN(K))GO TO 123                                              
  124 D1=D1-H1I                                                             
      IF(D1.GT.CN(K))GO TO 124                                              
      IF(D1.LT.-CN(K))GO TO 122                                             
      D=D1**2+D2**2                                                         
      IF(D.GT.CN(L))GO TO 124                                               
      IT(N)=K                                                            
      TI(N+1)=D                                                             
!-----If ICT is not zero, this is a C contact                            
      IF(ICT.EQ.0) GO TO 125                                             
      TI(N+2)=X3P3-W(3,J)-W(3,I)                                            
      GO TO 126                                                             
  125 TI(N+2)=W(3,J)-W(3,I)                                                 
  126 N=N+NSTP                                                              
      IF(N.GT.LIMIT)GO TO 902                                               
      GO TO 124                                                             
!-----Have the C molecules been collected? - non-zero means yes.         
  127 IF(ICT.NE.0)GO TO 128                                                 
!-----Switch to collecting C molecules                                   
      ICT=1                                                                 
!-----The vector between an A and a C molecule is equal to that from the 
!     origin to a P molecule.                                            
      D2=D2+X2P3                                                            
      D1=D1+X1P3                                                            
      GO TO 121                                                             
  128 CONTINUE                                                              
  129 CONTINUE                                                              
      NEND=N-1                                                              
      CALL MINHI(Y3A3,H3N,X3A3,H3M,KE)                                      
      IF(IPR.LT.4)GO TO 134                                              
      CALL PAGE(1,1,0)                                                   

      WRITE(NQ,131)KCT,Y3A3,H3N,X3A3,H3M,KE                              
  131 FORMAT(1X,'KCT =',I3,4F8.3,', KE =',I5)                            
  134 IF(KCT)135,139,141                                                    
  135 IF(IPR.LT.3)GO TO 158                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,136)MARK,NARK,NEND,KE,Y3A3,X3A3                           
  136 FORMAT(1X,'MARKS =',3I6,', KE =',I5,', Y3A3 =',2F8.3)              
      GO TO 158                                                          
!-----Try the negative axis-3 direction because of C and P molecules     
  139 Z3A3=-Y3A3                                                         
      CALL MINHI(Z3A3,-Y3A3,-100.,-H3M,KE)                               
      IF(KCT)142,140,906                                                 
!-----Adjust if the negative standoff is larger.                         
  140 Y3A3=-Z3A3                                                         
      IF(Y3A3.LT.X3A3)GO TO 142                                          
  141 IF(MCT.LE.1)GO TO 144                                              
      J=MCT-NCT                                                          
      IF((J.EQ.1).OR.(J.EQ.-7))GO TO 157                                 
      GO TO 143                                                          
  142 X1A3=Y1A3                                                          
      X2A3=Y2A3                                                          
      X3A3=Y3A3                                                          
      NCT=MCT                                                            
  143 GO TO(144,149,150,151,152,153,154,155,148),MCT                     
  144 IF(LE.GE.NV)GO TO 145                                              
      LE=LE+1                                                            
      Y2A3=Y2A3+DY2A3                                                    
      GO TO 116                                                          
  145 IF(JE.GE.NV)GO TO 146                                              
      JE=JE+1                                                            
      Y1A3=Y1A3+DY1A3                                                    
      LE=1                                                               
      Y2A3=-X2IH                                                         
      GO TO 116                                                          
  146 NCT=6                                                              
  147 JE=2.*JE                                                           
      DY1A3=H1I/JE                                                       
      DY2A3=X2I/JE                                                       
      XX1A3=X1A3                                                         
      XX2A3=X2A3                                                         
      NCT=NCT+1                                                          
      IF(NCT.GT.9)NCT=NCT-8                                              
      GO TO(906,150,151,152,153,154,155,148,149),NCT                     
  148 MCT=2                                                              
      Y1A3=XX1A3+DY1A3                                                   
      Y2A3=XX2A3                                                         
      GO TO 116                                                          
  149 MCT=3                                                              
      Y1A3=XX1A3+DY1A3                                                   
      Y2A3=XX2A3+DY2A3                                                   
      GO TO 116                                                          
  150 MCT=4                                                              
      Y1A3=XX1A3                                                         
      Y2A3=XX2A3+DY2A3                                                   
      GO TO 116                                                          
  151 MCT=5                                                              
      Y1A3=XX1A3-DY1A3                                                   
      IF(Y1A3.LT.-H1IH)Y1A3=Y1A3+H1I                                     
      Y2A3=XX2A3+DY2A3                                                   
      GO TO 116                                                          
  152 MCT=6                                                              
      Y1A3=XX1A3-DY1A3                                                   
      IF(Y1A3.LT.-H1IH)Y1A3=Y1A3+H1I                                     
      Y2A3=XX2A3                                                         
      GO TO 116                                                          
  153 MCT=7                                                              
      Y1A3=XX1A3-DY1A3                                                   
      IF(Y1A3.LT.-H1IH)Y1A3=Y1A3+H1I                                     
      Y2A3=XX2A3-DY2A3                                                   
      IF(Y2A3.LT.-X2IH)Y2A3=Y2A3+X2I                                     
      GO TO 116                                                          
  154 MCT=8                                                              
      Y1A3=XX1A3                                                         
      Y2A3=XX2A3-DY2A3                                                   
      IF(Y2A3.LT.-X2IH)Y2A3=Y2A3+X2I                                     
      GO TO 116                                                          
  155 MCT=9                                                              
      Y1A3=XX1A3+DY1A3                                                   
      Y2A3=XX2A3-DY2A3                                                   
      IF(Y2A3.LT.-X2IH)Y2A3=Y2A3+X2I                                     
      GO TO 116                                                          
  157 IF(JE.LT.KE)GO TO 147                                              
      GO TO 160                                                          
!-----H3N is the smallest possible standoff, check negative value        
  158 Z3A3=-H3N                                                          
      CALL MINHI(Z3A3,-H3N,-100.,-H3M,KE)                                
      IF(KCT)159,140,906                                                 
!     Smallest value confirmed, terminate search                         
  159 X3A3=H3N                                                           
      X2A3=Y2A3                                                          
      X1A3=Y1A3                                                          
  160 VY=H1I*X2P3*X3A3                                                      
      IF(IPR.LT.3)GO TO 162                                                 
      WRITE(NQ,161)X1A3,X2A3,X3A3,VY                                        
  161 FORMAT(1X,'X1A3 =',F8.3,', X2A3 =',F8.3,', X3A3 =',F8.3,', VY =' &  
     &,F8.2)                                                             
  162 IF(VX-VY.LT..01)GO TO 165                                             
      VX=VY                                                                 
      H1P3=X1P3                                                          
      H2P3=X2P3                                                          
      H3P3=X3P3                                                          
      H1A3=X1A3                                                          
      H2A3=X2A3                                                          
      H3A3=X3A3                                                          
  163 IF(IE.GE.NV)GO TO 167                                                 
      IE=IE+1                                                               
      X3P3=X3P3+DX3P3                                                    
      GO TO 31                                                           
  164 X3P3=H3P3-DX3P3                                                    
      JCT=-1                                                                
      GO TO 31                                                              
  165 IF(JCT)166,163,164                                                    
  166 IF(IE.GE.KE)GO TO 168                                                 
  167 JCT=1                                                                 
      IE=2*IE                                                               
      DX3P3=H3P3R/IE                                                        
      X3P3=H3P3+DX3P3                                                    
      GO TO 31                                                              
  168 IF(KE.GE.NNE)GO TO 169                                                
!-----Shift to maximum accuracy                                          
      KE=NNE                                                                
      GO TO 167                                                          
  169 VX=H1I*H2P3*H3A3                                                      
!-----Prepare lowest volume found output                                 
      H21I=2.*H1P3                                                       
      H2I=2.*H2P3                                                        
      H3I=2.*H3A3                                                        
      H1C=H1A3+H1P3                                                      
      H2C=H2A3+H2P3                                                      
      H3C=H3A3+H3P3                                                      
      IF(H2C.LE.H2I/2.)GO TO 170                                         
      H2C=H2C-H2I                                                        
      H1C=H1C-H21I                                                       
  170 IF(H2C.GE.-H2I/2.)GO TO 171                                        
      H2C=H2C+H2I                                                        
      H1C=H1C+H21I                                                       
  171 IF(H1C.GT.H1I/2.)H1C=H1C-H1I                                       
      IF(H1C.LT.-H1I/2.)H1C=H1C+H1I                                      
      IF(V3-VX.LT..01)GO TO 174                                             
      V3=VX                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA                                                    
      IF(V3.GE.HT1(8))GO TO 174                                             
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
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=10                                                               
      MLOW=14                                                            
  174 IF(IPR.LT.2)GO TO 1000                                             
      A(1)=H1I                                                           
      A(2)=SQRT(H21I**2+H2I**2)                                          
      A(3)=H3I                                                           
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
      A(8)=H2P3/SINE3                                                    
      A(7)=H1P3-A(8)*C(3)                                                
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
      CALL PAGE(4,4,0)                                                   
      WRITE(NQ,175)ALPHA,A1,A2,A3,VX,(A(I),I=1,3),(C(I),I=1,6)           
  175 FORMAT(1X,A2,3F7.1,'  V =',F7.2,'  AXES',3F7.3   &
     &/1X,'  COSINES',3F7.4,'  ANGLES',3F7.2)                            
  176 FORMAT(1X,A2,' AT',3F7.4,5X,3F7.4)
      SM(1)=A(4)-A(4)
      SM(2)=A(5)-A(5)
      SM(3)=A(6)-A(6)
      WRITE(NQ,176)SYM(1),(A(I),I=4,6),SM
      SM(1)=A(10)-A(4)
      SM(2)=A(11)-A(5)
      SM(3)=A(12)
      WRITE(NQ,176)SYM(2),(A(I),I=10,12),SM
      SM(1)=A(7)
      SM(2)=A(8)
      SM(3)=A(9)-A(6)
      WRITE(NQ,176)SYM(3),(A(I),I=7,9),SM
      GO TO 1000                                                            
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZFRTB,5X,3(I9,I3),I9)                  
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
      WRITE(NQ,921)KCT,MCT,N,MARK3,MARK,NEND3                               
  921 FORMAT(16H0NO INTERACTIONS,2I3,4I5)                                   
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZFRTB            

