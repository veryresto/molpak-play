      SUBROUTINE ZERTB            
!
      USE molpakCommonMod
      
      IMPLICIT NONE      
!                                                          
!      CHARACTER*80 BUF                                                      
!      CHARACTER*60 HEAD                                                 
!      CHARACTER*4 AN,CDIJ,CHT1                                          
!      CHARACTER*2 AT,CODE                                                
!----General program communication and transfer variables              
!     COMMON BUF,KILL,NPAGE,LINE,NEW,JTR(30),ITR(90)                      
!----Energy constants                                                      
!      COMMON NSEP,CN(1200),ERM,FERM,ER,ERMT,ERMF,NV,NE,NNE                  
!----Compound input data                                               
!      COMMON NMOD,NCTR,IA(200),IAA(200),X(3,200),AT(200),AN(200)          
!----Compound rotation and structure parameters                        
!      COMMON CS2(37),SN2(37),CS3(19),SN3(19),V(3,100),W(3,100)          
!      COMMON NLOW,MLOW,NP2,KCT                                          
!      COMMON H1M,H2M,H3M,H1N,H2N,H3N,A1,A2,A3,V3                        
!----Output flags and storage                                          
!      COMMON NQ,NR,NS,NWM,NSGT,IPR,HEAD,DIJ(500),HT1(14)                
!      DIMENSION IJD(500),CDIJ(500),CHT1(14)                             
!      EQUIVALENCE (DIJ(1),IJD(1)),(DIJ(1),CDIJ(1)),(HT1(1),CHT1(1))     
!      COMMON ICOUNT,MCOUNT,IORDER(500)                                  
!      COMMON ICOUNT,MCOUNT,IORDER(5000)                                  
!      COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500)  
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
!
      CHARACTER (3) :: ALPHA(3) = (/'CD', 'CE', 'CF'/)
      CHARACTER (2) :: SYM(7)=(/'P1','P2','P3','A1','A2','A3','C '/)
!
      INTEGER :: I, ICT, IE, J, JE, JCL, JCT 
      INTEGER :: K, KE, KTR, L, LE, LOC, LTR 
      INTEGER :: N, NCT, NNI, NSK1, NSK2, NSK3, NY2P2
!
      REAL :: D, D1, D2, D3, DX1P1, DY2P2, DZ3P3
      REAL :: H1C, H1I, H1P1, H1P3, H1P1F, H1P1S, H1P1R 
      REAL :: H2C, H2I, H2P1, H2P2, H2P2R, H2P2S  
      REAL :: H3C, H3I, H3P2, H3P3, H3P3R, H3P3S
      REAL :: SM(6), VH, VX, VY, VZ
      REAL :: X1I, X1P1, X1P3, X2P1 
      REAl :: X2I, X2P2, X3I, X3P2, X3P3
      REAL :: Y1I, Y3I, Y1P3, Y2P2, Y3P2, Y3P3 
      REAL :: Z1I, Z1IH, Z1P3, Z1P3N, Z2I, Z2IH, Z3P3   
!
      NARK=0                                                            
      NNI=0                                                             
      VH=1000000.                                                       
      KTR=ITR(41)                                                       
      IF(IPR.LT.2)GO TO 10                                              
      WRITE(NQ,5)(ITR(I),I=41,43)                                       
    5 FORMAT(1X,'ZERTB called',3I3)                                     
   10 IF(KTR-3)11,910,1000                                              
!-----Calculate minimum I molecule separation along axis-1               
   11 MARK=MARK2                                                        
      N=MARK                                                            
      DO 22 I=MARK1,NEND1,NSTP1                                            
      K=IT(I)                                                              
      IF(K.LE.0)GO TO 22                                                
      IF(IT(I+7).NE.0)GO TO 22                                          
      D2=TI(I+3)                                                           
      IF(ABS(D2).GT.CN(K))GO TO 22                                         
      D=TI(I+1)+D2**2                                                      
      L=IT(I+4)                                                            
      IF(D.GT.CN(L))GO TO 22                                               
      IT(N)=K                                                              
      TI(N+1)=D                                                            
      TI(N+2)=TI(I+2)                                                      
      N=N+NSTP                                                             
      IF(N.GT.LIMIT)GO TO 902                                              
   22 CONTINUE                                                             
      IF(N.LE.MARK)GO TO 912                                            
      ER=ERM                                                               
      NEND=N-1                                                             
      Z1I=H1M                                                              
      CALL MINHI(Z1I,0.,100.,H1M,NNE)                                      
      IF(KCT)1903,25,906                                                
   25 Z1IH=.5*Z1I                                                       
      IF(IPR.LT.3)GO TO 30                                              
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,26)Z1I                                                   
   26 FORMAT(1X,'Z1I =',F7.3)                                           
!-----Calculate minimum I molecule separation along axis-2               
   30 N=MARK                                                            
      DO 32 I=MARK1,NEND1,NSTP1                                            
      K=IT(I)                                                              
      IF(K.LE.0)GO TO 32                                                
      IF(IT(I+7).NE.0)GO TO 32                                          
      D1=TI(I+2)                                                           
      IF(ABS(D1).GT.CN(K))GO TO 32                                         
      D=TI(I+1)+D1**2                                                      
      L=IT(I+4)                                                            
      IF(D.GT.CN(L))GO TO 32                                               
      IT(N)=K                                                              
      TI(N+1)=D                                                            
      TI(N+2)=TI(I+3)                                                      
      N=N+NSTP                                                             
      IF(N.GT.LIMIT)GO TO 902                                              
   32 CONTINUE                                                             
      IF(N.LE.MARK)GO TO 913                                            
      ER=ERM                                                               
      NEND=N-1                                                             
      Z2I=H2M                                                              
      CALL MINHI(Z2I,0.,100.,H2M,NNE)                                     
      IF(KCT)2903,40,906                                                
   40 Z2IH=.5*Z2I                                                       
      IF(IPR.LT.3)GO TO 43                                              
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,41)Z2I                                                   
   41 FORMAT(1X,'Z2I =',F7.3)                                           
!-----Calculate range of possible values of H1P1 by finding the          
!     positive and negative axis-1 contacts with the axis-2 offset set   
!     at its minimum value, Z2IH                                         
   43 N=MARK                                                            
      DO 45 I=MARK1,NEND1,NSTP1                                         
      K=IT(I)                                                           
      IF(K.LE.0)GO TO 45                                                
      IF(IT(I+7).NE.0)GO TO 45                                          
      D2=Z2IH+TI(I+3)                                                   
      IF(ABS(D2).GT.CN(K))GO TO 45                                      
      D=TI(I+1)+D2**2                                                   
      L=IT(I+4)                                                         
      IF(D.GT.CN(L))GO TO 45                                            
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=TI(I+5)                                                   
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
   45 CONTINUE                                                          
      IF(N.LE.MARK)GO TO 914                                            
      NEND=N-1                                                          
      ER=ERM                                                            
      H1P1F=H1N                                                         
      CALL MINHI(H1P1F,0.,100.,H1M,NNE)                                 
      IF(KCT)46,46,906                                                  
   46 H1P1S=-H1N                                                        
      CALL MINHI(H1P1S,0.,-100.,-H1M,NNE)                               
      IF(KCT)47,47,906                                                  
   47 H1P1R=H1P1F-H1P1S                                                 
      IF(IPR.LT.3)GO TO 50                                              
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,48)H1P1S,H1P1F                                           
   48 FORMAT(1X,'H1P1S =',F7.3,',  H1P1F =',F7.3)                       
!-----Calculate range of possible values of H3P3 by finding the          
!     positive and negative axis-3 contacts with the axis-1 offset set   
!     at its minimum value, Z1IH                                         
   50 N=MARK                                                            
      DO 55 I=1,NMOD                                                    
      DO 54 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                    
      D1=Z1IH+W(1,J)-W(1,I)                                             
      IF(ABS(D1).GT.CN(K))GO TO 54                                      
      D2=W(2,J)-W(2,I)                                                  
      IF(ABS(D2).GT.CN(K))GO TO 54                                      
      D=D1**2+D2**2                                                     
      IF(D.GT.CN(K+10))GO TO 54                                         
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=-W(3,J)-W(3,I)                                            
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
   54 CONTINUE                                                          
   55 CONTINUE                                                          
      IF(N.LE.MARK)GO TO 915                                            
      NEND=N-1                                                          
      ER=ERM                                                            
      H3P3R=H3N                                                         
      CALL MINHI(H3P3R,0.,100.,H3M,NNE)                                 
      IF(KCT)56,56,906                                                  
   56 H3P3S=-H3N                                                        
      CALL MINHI(H3P3S,0.,-100.,-H3M,NNE)                               
      IF(KCT)57,57,906                                                  
   57 H3P3R=H3P3R-H3P3S                                                 
      IF(IPR.LT.3)GO TO 60                                              
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,59)H3P3S,H3P3R                                           
   59 FORMAT(1X,'H3P3S =',F7.3,',  H3P3R =',F7.3)                       
!-----Calculate range of possible values of H2P2 by finding the          
!     positive and negative axis-2 contacts with the axis-3 and axis-1   
!     offsets set at their minimum values, H3N and Z1IH.                 
   60 N=MARK                                                            
      NSK1=0
	NSK2=0
	NSK3=0
      DO 65 I=1,NMOD                                                    
      DO 64 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                    
      D3=H3N+W(3,J)-W(3,I)                                              
      IF(ABS(D3).GT.CN(K))GO TO 64                                      
      D1=Z1IH+W(1,J)-W(1,I)                                             
      IF(ABS(D1).GT.CN(K))GO TO 64                                      
      D=D1**1+D3**2                                                     
      IF(D.GT.CN(K+10))GO TO 64                                         
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=-W(2,J)-W(2,I)                                            
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
   64 CONTINUE                                                          
   65 CONTINUE                                                          
      IF(N.LE.MARK)GO TO 70                                             
      NEND=N-1                                                          
      ER=ERM                                                            
      NY2P2=0                                                            
      H2P2R=H2N                                                         
      CALL MINHI(H2P2R,0.,100.,H2M,NNE)                                 
      IF(KCT)66,66,906                                                  
   66 H2P2S=-H2N                                                        
      CALL MINHI(H2P2S,0.,-100.,-H2M,NNE)                               
      IF(KCT)77,77,906                                                  
   70 H2P2S=0.                                                           
      H2P2R=0.                                                           
      NY2P2=1                                                            
   77 H2P2R=H2P2R-H2P2S                                                 
      IF(IPR.LT.3)GO TO 100                                             
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,78)H2P2S,H2P2R                                           
   78 FORMAT(1X,'H2P2S =',F7.3,',  H2P2R =',F7.3)                       
!-----Determination of H1P1 and H2P1                                     
  100 ICT=0                                                             
      KE=NNE                                                            
!-----Do value of X1P1 where I molecules are in conact at the far end    
!     of its possible range.                                              
      X1P1=H1P1F                                                        
      X2P1=Z2IH                                                         
      X2I=Z2I                                                           
      LTR=2                                                             
      IE=0                                                              
      NARK=0                                                            
      GO TO 110                                                         
!-----Set up X1P1 placement iteration loop - Code CD                     
  103 X1P1=H1P1S                                                        
      DX1P1=H1P1R/NV                                                    
      X2P1=Z2IH                                                         
      X2I=Z2I                                                           
      IE=1                                                              
      NARK=0                                                            
!-----I molecules are in contact at the start of X1P1 range.             
      GO TO 110                                                         
!-----Top of X1P1 placement iteration loop - Code CD                     
  105 MARK=MARK2                                                        
      NARK=0                                                            
      N=MARK                                                            
      DO 106 I=MARK1,NEND1,NSTP1                                        
      K=IT(I)                                                           
      IF(K.LE.0)GO TO 106                                               
      IF(IT(I+7).NE.0)GO TO 106                                         
      D1=X1P1+TI(I+5)                                                   
      IF(ABS(D1).GT.CN(K))GO TO 106                                     
      D=TI(I+1)+D1**2                                                   
      IF(D.GT.CN(K+10))GO TO 106                                        
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=TI(I+3)                                                   
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
  106 CONTINUE                                                          
      IF(N.LE.MARK)GO TO 108                                            
      NEND=N-1                                                          
      ER=ERM                                                            
      CALL MINHI(X2P1,Z2IH,100.,H2M,KE)                                 
      IF(KCT)109,109,906                                                
  108 X2P1=Z2IH                                                         
  109 X2I=2.*X2P1                                                       
  110 IF(IPR.LT.4)GO TO 115                                             
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,111)X1P1,X2P1                                            
  111 FORMAT(1X,'X1P1 =',F7.3,',  X2P1 =',F7.3)                         
!-----Set up Y2P2 placement iteration loop                               
  115 IF(IE.LE.NV)GO TO 119                                             
      IF(KE.LE.NE)GO TO 119                                             
      Y3P2=H3P2                                                         
      LE=NE                                                             
      DY2P2=H2P2R/LE                                                    
      Y2P2=H2P2+DY2P2                                                   
      JCT=1                                                             
      GO TO 120                                                         
  119 Y3P2=H3N                                                          
      Y2P2=H2P2S                                                        
      DY2P2=H2P2R/NV                                                    
      JCT=0                                                             
      LE=1                                                              
  120 VX=1000000.           ! 03-15-06                                                      
!-----Top of Y2P2 placement iteration loop                               
  121 N=MARK                                                            
      NARK=0                                                            
      DO 127 I=1,NMOD                                                   
      DO 126 J=1,NMOD                                                   
      K=IA(I)+IAA(J)                                                    
      L=K+10                                                            
!-----First enter P2 3 molecules                                         
      JCL=1                                                             
      D1=W(1,J)-W(1,I)                                                  
      D2=Y2P2-W(2,J)-W(2,I)                                             
  122 IF(ABS(D1).GT.CN(K))GO TO 125                                     
  123 D2=D2+X2I                                                         
      IF(D2.LT.CN(K))GO TO 123                                          
  124 D2=D2-X2I                                                         
      IF(D2.GT.CN(K))GO TO 124                                          
      IF(D2.LT.-CN(K))GO TO 125                                         
      D=D1**2+D2**2                                                     
      IF(D.GT.CN(L))GO TO 124                                           
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=W(3,J)-W(3,I)                                             
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
      GO TO 124                                                         
!-----Have S3 molecules been added?                                      
  125 IF(JCL.LT.0)GO TO 126                                             
      JCL=-1                                                            
      D1=X1P1-W(1,J)-W(1,I)                                             
      D2=D2+X2P1                                                        
      GO TO 122                                                         
  126 CONTINUE                                                          
  127 CONTINUE                                                          
      IF(N.LE.MARK)GO TO 128                                            
      NEND=N-1                                                          
      ER=ERM                                                            
      CALL MINHI(Y3P2,H3N,100.,H3M,KE)                                  
      IF(KCT)129,129,906                                                
  128 Y3P2=H3N                                                          
  129 Y3I=2.*Y3P2                                                       
      IF(IPR.LT.4)GO TO 131                                             
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,130)Y2P2,Y3P2                                            
  130 FORMAT(1X,'Y2P2 =',F7.3,',  Y3P2 =',F7.3)                         
!-----Collect P1 2, P2 3, S3, and I molecules for double standoff        
  131 N=MARK                                                            
      DO 141 I=1,NMOD                                                   
      DO 140 J=1,NMOD                                                   
      K=IA(I)+IAA(J)                                                    
      L=K+10                                                            
      D1=X1P1-W(1,J)-W(1,I)                                             
      D2=X2P1+W(2,J)-W(2,I)                                             
      D3=W(3,J)-W(3,I)                                                  
!-----Set flag for P1 molecule contacts                                 
      JCL=1                                                             
  132 D3=D3+Y3I                                                         
      IF(D3.LT.CN(K))GO TO 132                                          
  133 D3=D3-Y3I                                                         
      IF(D3.GT.CN(K))GO TO 133                                          
      IF(D3.LT.-CN(K))GO TO 136                                         
  134 D2=D2+X2I                                                         
      IF(D2.LT.CN(K))GO TO 134                                          
  135 D2=D2-X2I                                                         
      IF(D2.GT.CN(K))GO TO 135                                          
      IF(D2.LT.-CN(K))GO TO 133                                         
      D=D2**2+D3**2                                                     
      IF(D.GT.CN(L))GO TO 135                                           
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=D1                                                        
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
      GO TO 135                                                         
  136 GO TO(137,138,139,140),JCL                                        
!-----S3 molecule contacts                                               
  137 JCL=2                                                             
      D2=X2P1+Y2P2-W(2,J)-W(2,I)                                        
      D3=D3+Y3P2                                                        
      GO TO 132                                                         
!-----P2 3 molecule contacts                                             
  138 JCL=3                                                             
      D1=W(1,J)-W(1,I)                                                  
      D2=D2-X2P1                                                        
      GO TO 132                                                         
!-----I molecule contacts                                                
  139 JCL=4                                                             
      D2=W(2,J)-W(2,I)                                                  
      D3=D3-Y3P2                                                        
      GO TO 132                                                         
  140 CONTINUE                                                          
  141 CONTINUE                                                          
      NARK=N                                                            
!-----Set up Z3P3 placement iteration loop                               
      IF(IE.LE.NV)GO TO 149                                             
      IF(KE.LE.NE)GO TO 149                                             
      Z1P3=H1P3                                                         
      JE=NE                                                             
      DZ3P3=H3P3R/JE                                                    
      Z3P3=H3P3+DZ3P3                                                   
      NCT=1                                                             
      GO TO 150                                                         
  149 Z1P3=Z1IH                                                         
      Z3P3=H3P3S                                                        
      DZ3P3=H3P3R/NV                                                    
      NCT=0                                                             
      JE=1                                                              
  150 VY=1000000.          ! 03-15-06                                                        
      ER=ERMF                                                           
!-----Top of Z3P3 placement iteration loop                               
  151 N=NARK                                                            
      DO 165 I=1,NMOD                                                   
      DO 164 J=1,NMOD                                                   
      K=IA(I)+IAA(J)                                                    
      L=K+10                                                            
!-----Set flag for P3 n contacts                                         
      JCL=1                                                             
      D1=W(1,J)-W(1,I)                                                  
      D2=X2P1+W(2,J)-W(2,I)                                             
      D3=Z3P3-W(3,J)-W(3,I)                                             
  152 D2=D2+X2I                                                         
      IF(D2.LT.CN(K))GO TO 152                                          
  153 D2=D2-X2I                                                         
      IF(D2.GT.CN(K))GO TO 153                                          
      IF(D2.LT.-CN(K))GO TO 159                                         
  154 D3=D3+Y3I                                                         
      IF(D3.LT.CN(K))GO TO 154                                          
  155 D3=D3-Y3I                                                         
      IF(D3.GT.CN(K))GO TO 155                                          
      IF(D3.LT.-CN(K))GO TO 153                                         
      D=D2**2+D3**2                                                     
      IF(D.GT.CN(L))GO TO 155                                           
      IT(N)=K                                                           
      TI(N+1)=D                                                         
      TI(N+2)=D1                                                        
      N=N+NSTP                                                          
      IF(N.GT.LIMIT)GO TO 902                                           
      GO TO 155                                                         
!-----Have S1, A2, and C contacts been calculated?                       
  159 GO TO(160,161,162,164),JCL                                        
!-----S1 molecule contacts                                               
  160 JCL=2                                                             
      D2=X2P1+Y2P2-W(2,J)-W(2,I)                                        
      D3=D3+Y3P2                                                        
      GO TO 152                                                         
!-----C molecule contacts                                                
  161 JCL=3                                                             
      D1=X1P1-W(1,J)-W(1,I)                                             
      D2=D2-X2P1                                                        
      GO TO 152                                                         
!-----A2 molecule contacts                                               
  162 JCL=4                                                             
      D2=W(2,J)-W(2,I)                                                  
      D3=D3-Y3P2                                                        
      GO TO 152                                                         
  164 CONTINUE                                                          
  165 CONTINUE                                                          
      IF(N.LE.NARK)GO TO 169                                            
      NEND=N-1                                                          
      CALL MINHI(Z1P3,Z1IH,100.,Z1IH,KE)                                
      IF(KCT)168,168,906                                                
  168 Z1P3N=-Z1P3                                                       
      CALL MINHI(Z1P3N,-Z1P3,-100.,-Z1IH,KE)                            
!-----If negative standoff is greater than positive, reset               
      IF(KCT)171,170,906                                                
  169 Z1P3=Z1IH                                                         
      GO TO 171                                                         
  170 Z1P3=-Z1P3N                                                       
  171 Z1I=2.*Z1P3                                                       
!-----Bottom of P3 n molecule placement iteration loop                   
      VZ=Z1P3*X2P1*Y3P2                                                 
      IF(IPR.LT.4)GO TO 205                                             
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,204)Z3P3,Z1P3,VZ                                         
  204 FORMAT(1X,'Z3P3 =',F7.3,', Z1P3 =',F7.3,', VZ =',F9.2)            
  205 IF(VZ.GE.VY)GO TO 210                                             
      VY=VZ                                                             
      Y3P3=Z3P3                                                         
      Y1P3=Z1P3                                                         
      Y1I=Z1I                                                           
      IF(NCT)212,211,212                                                
  210 IF(NCT)212,211,214                                                
  211 IF(JE.GE.NV)GO TO 212                                             
  213 JE=JE+1                                                           
      Z3P3=Z3P3+DZ3P3                                                   
      GO TO 151                                                         
  212 IF(JE.GE.KE)GO TO 219                                             
      NCT=1                                                             
      JE=2*JE                                                           
      DZ3P3=H3P3R/JE                                                    
      Z3P3=Y3P3+DZ3P3                                                   
      GO TO 151                                                         
  214 NCT=-1                                                            
      Z3P3=Y3P3-DZ3P3                                                   
      IF(Z3P3-H3P3S)212,151,151                                         
!-----Bottom of P2 3 molecule placement iteration loop                   
  219 IF(VY.GE.VX)GO TO 230                                             
      VX=VY                                                             
      X2P2=Y2P2                                                         
      X3P2=Y3P2                                                         
      X3I=Y3I                                                           
      X3P3=Y3P3                                                         
      X1P3=Y1P3                                                         
      X1I=Y1I                                                           
      IF(NY2P2.EQ.1)GO TO 235                                            
      IF(JCT)232,231,232                                                
  230 IF(JCT)232,231,234                                                
  231 IF(LE.GE.NV)GO TO 232                                             
      LE=LE+1                                                           
      Y2P2=Y2P2+DY2P2                                                   
      GO TO 121                                                         
  232 IF(LE.GE.KE)GO TO 235                                             
      JCT=1                                                             
      LE=2*LE                                                           
      DY2P2=H2P2R/LE                                                    
      Y2P2=X2P2+DY2P2                                                   
      GO TO 121                                                         
  234 JCT=-1                                                            
      Y2P2=X2P2-DY2P2                                                   
      IF(Y2P2-H2P2S)232,121,121                                         
  235 IF(IPR.LT.3)GO TO 239                                             
      CALL PAGE(1,1,0)                                                  
      WRITE(NQ,236)X1P1,X2P2,X3P3,X1P3,X2P1,X3P2,VX                     
  236 FORMAT(' ',6F7.3,', VX =',F9.2)                                   
!-----Bottom of P1 2 molecule placement iteration loop                   
  239 IF(VX.GE.VH)GO TO 280                                             
      IF(IE.GT.1)LTR=1                                                  
      VH=VX                                                             
      H2P2=X2P2                                                         
      H3P2=X3P2                                                         
      H1P3=X1P3                                                         
      H1P1=X1P1                                                         
      H2P1=X2P1                                                         
      H2P2=X2P2                                                         
      H3P3=X3P3                                                         
      H1I=X1I                                                           
      H2I=X2I                                                           
      H3I=X3I                                                           
      IF(ICT)287,281,287                                                
  280 IF(ICT)287,281,289                                                
!-----Have the ends of the H1P1 range been calculated?                   
  281 IF(IE-1)103,282,286                                               
!     If Code CE was specified, only the ends are calculated.            
  282 IF(KTR-2)285,300,1000                                             
  285 KE=NE                                                             
  286 IF(IE.GE.NV)GO TO 287                                             
      IE=IE+1                                                           
      X1P1=X1P1+DX1P1                                                   
      GO TO 105                                                         
  287 IF(IE.GE.KE)GO TO 290                                             
  288 ICT=1                                                             
      IE=2*IE                                                           
      DX1P1=H1P1R/IE                                                    
      X1P1=H1P1+DX1P1                                                   
      IF(X1P1-H1P1F)105,105,289                                         
  289 ICT=-1                                                            
      X1P1=H1P1-DX1P1                                                   
      IF(X1P1-H1P1S)287,105,105                                         
  290 IF(IE.GE.NNE)GO TO 300                                            
      KE=NNE                                                            
      GO TO 288                                                         
  300 H1C=H1P1+H1P3                                                     
      H2C=H2P2                                                          
      H3C=H3P3+H3P2                                                     
      IF(H1C.GT.H1I/2.)H1C=H1C-H1I                                      
      IF(H1C.LT.-H1I/2.)H1C=H1C+H1I                                     
      IF(H2C.GT.H2I/2.)H2C=H2C-H2I                                      
      IF(H2C.LT.-H2I/2.)H2C=H2C+H2I                                     
      IF(H3C.GT.H3I/2.)H3C=H3C-H3I                                      
      IF(H3C.LT.-H3I/2.)H3C=H3C+H3I                                     
      IF(IPR.LT.2)GO TO 320                                             
      CALL PAGE(8,8,0)                                                  
      WRITE(NQ,301)ALPHA(LTR),A1,A2,A3,VH,H1I,H2I,H3I                   
  301 FORMAT(1X,A2,3F7.1,'  V =',F8.2,'  AXES',3F8.3)                   
      DO 316 I=1,7                                                      
      GO TO (305,306,307,308,309,310,311),I                             
  305 SM(1)=H1P1                                                        
      SM(2)=H2P1                                                        
      SM(3)=0.                                                          
      GO TO 312                                                         
  306 SM(1)=0.                                                          
      SM(2)=H2P2                                                        
      SM(3)=H3P2                                                        
      GO TO 312                                                         
  307 SM(1)=H1P3                                                        
      SM(2)=H2P1                                                        
      SM(3)=H3P3                                                        
      GO TO 312                                                         
  308 SM(1)=H1P3                                                        
      SM(2)=H2P1+H2P2                                                   
      SM(3)=H3P2+H3P3                                                   
      GO TO 312                                                         
  309 SM(1)=H1P1+H1P3                                                   
      SM(2)=0.                                                          
      SM(3)=H3P3                                                        
      GO TO 312                                                         
  310 SM(1)=H1P1                                                        
      SM(2)=H2P1+H2P2                                                   
      SM(3)=H3P2                                                        
      GO TO 312                                                         
  311 SM(1)=H1C                                                         
      SM(2)=H2C                                                         
      SM(3)=H3C                                                         
  312 SM(4)=SM(1)/H1I                                                   
      SM(5)=SM(2)/H2I                                                   
      SM(6)=SM(3)/H3I                                                   
      IF(SM(4).GT..5)SM(4)=SM(4)-1.                                     
      IF(SM(4).LE.-.5)SM(4)=SM(4)+1.                                    
      IF(SM(5).GT..5)SM(5)=SM(5)-1.                                     
      IF(SM(5).LE.-.5)SM(5)=SM(5)+1.                                    
      IF(SM(6).GT..5)SM(6)=SM(6)-1.                                     
      IF(SM(6).LE.-.5)SM(6)=SM(6)+1.                                    
  314 D=SQRT(SM(1)**2+SM(2)**2+SM(3)**2)                                
      WRITE(NQ,315)SYM(I),SM,D                                          
  315 FORMAT(1X,A2,' AT',3X,3F8.3,';  ',3F8.4,'  D =',F6.2)             
  316 CONTINUE                                                          
      SM(1)=H1C/2.                                                      
      SM(2)=H2C/2.                                                      
      SM(3)=H3C/2.                                                      
      SM(4)=SM(1)/H1I                                                   
      SM(5)=SM(2)/H2I                                                   
      SM(6)=SM(3)/H3I                                                   
      WRITE(NQ,319)SM                                                   
  319 FORMAT(1X,'ORIGIN',2X,3F8.3,';  ',3F8.4)                          
  320 IF(V3.LE.VH)GO TO 1000                                               
      V3=VH                                                                
      DIJ(NP2)=V3                                                          
      CDIJ(NP2+20)=ALPHA(LTR)                                              
      IF(V3.GE.HT1(8))GO TO 1000                                           
      HT1(1)=H1I                                                           
      HT1(2)=0.                                                            
      HT1(3)=H2I                                                           
      HT1(4)=H3I                                                           
      HT1(5)=H1C                                                           
      HT1(6)=H2C                                                           
      HT1(7)=H3C                                                           
      HT1(8)=V3                                                            
      CHT1(9)=ALPHA(LTR)                                                   
      HT1(10)=A1                                                           
      HT1(11)=A2                                                           
      HT1(12)=A3                                                           
      HT1(13)=0.                                                           
      HT1(14)=0.                                                           
      NLOW=32                                                              
      MLOW=LTR+40                                                          
      GO TO 1000                                                        
  902 CALL PAGE(2,2,0)                                                     
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                 
  903 FORMAT(26H0STORAGE EXCEEDED BY ZERTA,5X,3(I9,I3),I9)                 
      GO TO 999                                                            
 1903 LOC=25
      GO TO 904
 2903 LOC=40
  904 CALL PAGE(2,2,0)                                                  
      WRITE(NQ,905)LOC                                                  
  905 FORMAT(' PARAMETER TOO SMALL NEAR',I3)                            
      GO TO 999                                                         
  906 CALL PAGE(2,2,0)                                                  
      WRITE(NQ,907)                                                     
  907 FORMAT(20H0PARAMETER TOO LARGE)                                   
      GO TO 999                                                         
  910 CALL PAGE(2,2,0)                                                  
      WRITE(NQ,911)                                                     
  911 FORMAT(1H0,'Code CF, Pnma, programmed only with PLAN line (Z=4)') 
      GO TO 999                                                         
  912 NNI=1                                                             
      GO TO 920                                                         
  913 NNI=2                                                             
      GO TO 920                                                         
  914 NNI=3                                                             
      GO TO 920                                                         
  915 NNI=4                                                             
  920 CALL PAGE(2,2,0)                                                     
      WRITE(NQ,921)NNI,MARK1,MARK2,MARK,N,NSK1,NSK2,NSK3,H3N,Z1IH          
  921 FORMAT(16H0NO INTERACTIONS,I3,4I6,3I5,2F6.2)                         
      GO TO 999                                                            
  999 KILL=1                                                               
 1000 RETURN                                                               
      END SUBROUTINE ZERTB
