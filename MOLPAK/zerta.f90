      SUBROUTINE ZERTA     
!            
!     For CB and CC.
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
      CHARACTER (2) :: ALPHA(2) = (/'CB', 'CC'/)
      CHARACTER (2) :: SYM(7) = (/'P1','P2','P3','A1','A2','A3','C '/)
      
      INTEGER :: I, ICT, IE, J, JCL, JCT, JE,  K, KE, KTR
      INTEGER :: L, LE, N, MCT, NCT
!
      REAL :: D, D1, D2, D3, DX3P3, DY1P1, DZ1P1, DZ2P2
      REAL :: H1C, H1I, H1P1, H1P1R, H1P1S, H1P3
      REAL :: H2C, H2I, H2P1, H2P2, H2P2R, H2P2S
      REAL :: H3C, H3I, H3P2, H3P3, H3P3R, H3P3S
      REAL :: SM(6), VH, VX, VY, VZ
      REAL :: X1I, X1IH, X1P1, X1P3, X2P1, X2P2, X2I
      REAl :: X3I, X3P2, X3P3
      REAL :: Y1P1, Y2P1, Y2I, Y3I, Y2P2, Y3P2 
      REAL :: Z1I, Z1IH, Z1P1, Z2I, Z2IH, Z2P1
      REAL :: Z2P2, Z3I, Z3P2, Z3P2N 
!
      NARK=0                                                             
      VH=1000000.                                                        
      KTR=ITR(39)                                                        
      IF(IPR.LT.2)GO TO 10                                               
      WRITE(NQ,5)(ITR(I),I=39,40)                                        
    5 FORMAT(1X,'ZERTA called',2I3)                                      
   10 IF(KTR.GT.2)GO TO 1000                                             
!----Calculate minimum I molecule separation along axis-1              
      MARK=MARK2                                                         
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
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      Z1I=H1M                                                               
      CALL MINHI(Z1I,0.,100.,H1M,NNE)                                       
      IF(KCT)904,25,906                                                  
   25 Z1IH=.5*Z1I                                                        
      IF(IPR.LT.3)GO TO 30                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,26)Z1I                                                    
   26 FORMAT(1X,'Z1I =',F7.3)                                            
!----Calculate minimum I molecule separation along axis-2              
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
      IF(N.LE.MARK)GO TO 920                                             
      ER=ERM                                                                
      NEND=N-1                                                              
      Z2I=H2M                                                               
      CALL MINHI(Z2I,0.,100.,H2M,NNE)                                      
      IF(KCT)904,40,906                                                  
   40 Z2IH=.5*Z2I                                                        
      IF(IPR.LT.3)GO TO 43                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,41)Z2I                                                    
   41 FORMAT(1X,'Z2I =',F7.3)                                            
!----Calculate range of possible values of H1P1 by finding the         
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
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      H1P1R=H1N                                                          
      CALL MINHI(H1P1R,0.,100.,H1M,NNE)                                  
      IF(KCT)904,46,906                                                  
   46 H1P1S=-H1N                                                         
      CALL MINHI(H1P1S,0.,-100.,-H1M,NNE)                                
      IF(KCT)904,47,906                                                  
   47 H1P1R=H1P1R-H1P1S                                                  
      IF(IPR.LT.3)GO TO 50                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,49)H1P1S,H1P1R                                            
   49 FORMAT(1X,'H1P1S =',F7.3,',  H1P1R =',F7.3)                        
!----Calculate range of possible values of H3P3 by finding the         
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
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      H3P3R=H3N                                                          
      CALL MINHI(H3P3R,0.,100.,H3M,NNE)                                  
      IF(KCT)904,56,906                                                  
   56 H3P3S=-H3N                                                         
      CALL MINHI(H3P3S,0.,-100.,-H3M,NNE)                                
      IF(KCT)904,57,906                                                  
   57 H3P3R=H3P3R-H3P3S                                                  
      IF(IPR.LT.3)GO TO 60                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,58)H3P3S,H3P3R                                            
   58 FORMAT(1X,'H3P3S =',F7.3,',  H3P3R =',F7.3)                        
!----Calculate range of possible values of H2P2 by finding the         
!     positive and negative axis-2 contacts with the axis-3 offset set  
!     at its minimum value, H3N                                         
   60 N=MARK                                                             
      DO 65 I=1,NMOD                                                     
      DO 64 J=1,NMOD                                                     
      K=IA(I)+IAA(J)                                                     
      D3=H3N+W(3,J)-W(3,I)                                               
      IF(ABS(D3).GT.CN(K))GO TO 64                                       
      D1=W(1,J)-W(1,I)                                                   
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
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      H2P2R=H2N                                                          
      CALL MINHI(H2P2R,0.,100.,H2M,NNE)                                  
      IF(KCT)904,66,906                                                  
   66 H2P2S=-H2N                                                         
      CALL MINHI(H2P2S,0.,-100.,-H2M,NNE)                                
      IF(KCT)904,67,906                                                  
   67 H2P2R=H2P2R-H2P2S                                                  
      IF(IPR.LT.3)GO TO 70                                               
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,68)H2P2S,H2P2R                                            
   68 FORMAT(1X,'H2P2S =',F7.3,',  H2P2R =',F7.3)                        
   70 IF(KTR-2)100,500,1000                                              
!----Set up P3 placement iteration loop - Code CB                      
  100 X3P3=H3P3S                                                         
      DX3P3=H3P3R/NV                                                     
      ICT=0                                                              
      IE=1                                                               
      X1P3=Z1IH                                                          
      KE=NE                                                              
!----Top of P3 placement iteration loop                                
  105 MARK=MARK2                                                         
      NARK=0                                                             
      N=MARK                                                             
      DO 109 I=1,NMOD                                                    
      DO 108 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      D3=X3P3-W(3,J)-W(3,I)                                              
      D2=W(2,J)-W(2,I)                                                   
      IF(ABS(D3).GT.CN(K))GO TO 108                                      
      IF(ABS(D2).GT.CN(K))GO TO 108                                      
      D=D3**2+D2**2                                                      
      IF(D.GT.CN(K+10))GO TO 108                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(1,J)-W(1,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
  108 CONTINUE                                                           
  109 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(X1P3,Z1IH,100.,H1M,KE)                                  
      X1I=2.*X1P3                                                        
      IF(IPR.LT.3)GO TO 115                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,110)X3P3,X1P3                                             
  110 FORMAT(1X,'X3P3 =',F7.3,',  X1P3 =',F7.3)                          
!----Set up P1 placement iteration loop                                
  115 IF(KE.LE.NE)GO TO 119                                              
      Y2P1=H2P1                                                          
      LE=NE                                                              
      DY1P1=H1P1R/LE                                                     
      Y1P1=H1P1+DY1P1                                                    
      JCT=1                                                              
      GO TO 120                                                          
  119 Y2P1=Z2IH                                                          
      Y1P1=H1P1S                                                         
      DY1P1=H1P1R/NV                                                     
      JCT=0                                                              
      LE=1                                                               
  120 VX=1000000.        ! 03-15-06                                                         
!----Top of P1 placement iteration loop                                
  121 N=MARK                                                             
      NARK=0                                                             
!     Find P1 molecule contacts                                         
      DO 125 I=MARK1,NEND1,NSTP1                                         
      K=IT(I)                                                            
      IF(K.LE.0)GO TO 125                                                
      IF(IT(I+7).NE.0)GO TO 125                                          
      D1=Y1P1+TI(I+5)                                                    
  122 D1=D1+X1I                                                          
      IF(D1.LT.CN(K))GO TO 122                                           
  123 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 123                                           
      IF(D1.LT.-CN(K))GO TO 125                                          
      D=TI(I+1)+D1**2                                                    
      IF(D.GT.CN(K+10))GO TO 123                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+3)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 123                                                          
  125 CONTINUE                                                           
!----Add A2 molecule contacts                                          
      DO 129 I=1,NMOD                                                    
      DO 128 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      D3=X3P3-W(3,J)-W(3,I)                                              
      IF(ABS(D3).GT.CN(K))GO TO 128                                      
      D1=X1P3+Y1P1-W(1,J)-W(1,I)                                         
  126 D1=D1+X1I                                                          
      IF(D1.LT.CN(K))GO TO 126                                           
  127 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 127                                           
      IF(D1.LT.-CN(K))GO TO 128                                          
      D=D3**2+D1**2                                                      
      IF(D.GT.CN(K+10))GO TO 127                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 127                                                          
  128 CONTINUE                                                           
  129 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      ER=ERM                                                             
      CALL MINHI(Y2P1,Z2IH,100.,Z2I,KE)                                  
      Y2I=2.*Y2P1                                                        
      IF(IPR.LT.3)GO TO 131                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,130)Y1P1,Y2P1                                             
  130 FORMAT(1X,'Y1P1 =',F7.3,',  Y2P1 =',F7.3)                          
!----Collect P1, P3, A2, and I molecules for double standoff           
  131 N=MARK                                                             
      DO 141 I=1,NMOD                                                    
      DO 140 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      L=K+10                                                             
      D1=X1P3+W(1,J)-W(1,I)                                              
      D2=W(2,J)-W(2,I)                                                   
      D3=X3P3-W(3,J)-W(3,I)                                              
!----Set flag for P3 molecule contacts                                  
      JCL=1                                                              
  132 D1=D1+X1I                                                          
      IF(D1.LT.CN(K))GO TO 132                                           
  133 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 133                                           
      IF(D1.LT.-CN(K))GO TO 136                                          
  134 D2=D2+Y2I                                                          
      IF(D2.LT.CN(K))GO TO 134                                           
  135 D2=D2-Y2I                                                          
      IF(D2.GT.CN(K))GO TO 135                                           
      IF(D2.LT.-CN(K))GO TO 133                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(L))GO TO 135                                            
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=D3                                                         
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 135                                                          
  136 GO TO(137,138,139,140),JCL                                         
!----A2 molecule contacts                                              
  137 JCL=2                                                              
      D1=X1P3+Y1P1-W(1,J)-W(1,I)                                         
      D2=D2+Y2P1                                                         
      GO TO 132                                                          
!----P1 molecule contacts                                              
  138 JCL=3                                                              
      D1=D1-X1P3                                                         
      D3=W(3,J)-W(3,I)                                                   
      GO TO 132                                                          
!----I molecule contacts                                               
  139 JCL=4                                                              
      D1=W(1,J)-W(1,I)                                                   
      D2=D2-Y2P1                                                         
      GO TO 132                                                          
  140 CONTINUE                                                           
  141 CONTINUE                                                           
      NARK=N                                                             
!----Set up P2 placement iteration loop                                
      IF(KTR.EQ.2)GO TO 149                                              
      IF(KE.LE.NE)GO TO 149                                              
      Z3P2=H3P2                                                          
      JE=NE                                                              
      DZ2P2=H2P2R/JE                                                     
      Z2P2=H2P2+DZ2P2                                                    
      NCT=1                                                              
      GO TO 150                                                          
  149 Z3P2=H3N                                                           
      Z2P2=H2P2S                                                         
      DZ2P2=H2P2R/NV                                                     
      NCT=0                                                              
      JE=1                                                               
  150 VY=1000000.       ! 03-15-06                                                     
      ER=ERMF                                                            
!----Top of P2 placement iteration loop                                
  151 N=NARK                                                             
      DO 169 I=1,NMOD                                                    
      DO 168 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
!----Set flag for P2 contacts                                          
      JCL=1                                                              
      D2=Z2P2-W(2,J)-W(2,I)                                              
      D1=W(1,J)-W(1,I)                                                   
      D3=W(3,J)-W(3,I)                                                   
  152 D2=D2+Y2I                                                          
      IF(D2.LT.CN(K))GO TO 152                                           
  153 D2=D2-Y2I                                                          
      IF(D2.GT.CN(K))GO TO 153                                           
      IF(D2.LT.-CN(K))GO TO 159                                          
  154 D1=D1+X1I                                                          
      IF(D1.LT.CN(K))GO TO 154                                           
  155 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 155                                           
      IF(D1.LT.-CN(K))GO TO 153                                          
      D=D1**2+D2**2                                                      
      IF(D.GT.CN(K+10))GO TO 155                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=D3                                                         
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 155                                                          
!----Have A1, A3, and C contacts been calculated                       
  159 GO TO(160,161,162,168),JCL                                         
!----A3 molecule contacts (P2 + P1)                                    
  160 JCL=2                                                              
      D2=D2+Y2P1                                                         
      D1=Y1P1-W(1,J)-W(1,I)                                              
      GO TO 152                                                          
!----C molecule contacts (A3 + P3)                                     
  161 JCL=3                                                              
      D2=Y2P1+Z2P2-W(2,J)-W(2,I)                                         
      D1=D1+X1P3                                                         
      D3=X3P3-W(3,J)-W(3,I)                                              
      GO TO 152                                                          
!----A1 molecule contacts (P2 + P3)                                    
  162 JCL=4                                                              
      D2=D2-Y2P1                                                         
      D1=X1P3+W(1,J)-W(1,I)                                              
      GO TO 152                                                          
  168 CONTINUE                                                           
  169 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      CALL MINHI(Z3P2,H3N,100.,H3M,KE)                                   
      Z3P2N=-Z3P2                                                        
      CALL MINHI(Z3P2N,-Z3P2,-100.,-H3M,KE)                              
!----If negative standoff is greater than positive, reset              
      IF(KCT)171,170,906                                                 
  170 Z3P2=-Z3P2N                                                        
  171 Z3I=2.*Z3P2                                                        
!----Bottom of P2 molecule placement iteration loop                    
      VZ=X1P3*Y2P1*Z3P2                                                  
      IF(IPR.LT.3)GO TO 200                                              
      WRITE(NQ,195)Z2P2,Z3P2,VZ                                          
  195 FORMAT('Z2P2 =',F8.3,', Z3P2 =',F8.3,', VZ =',F10.3)               
  200 IF(VZ.GE.VY)GO TO 210                                              
      VY=VZ                                                              
      Y2P2=Z2P2                                                          
      Y3P2=Z3P2                                                          
      Y3I=Z3I                                                            
      IF(NCT)212,211,212                                                 
  210 IF(NCT)212,211,214                                                 
  211 IF(JE.GE.NV)GO TO 212                                              
  213 JE=JE+1                                                            
      Z2P2=Z2P2+DZ2P2                                                    
      GO TO 151                                                          
  212 IF(JE.GE.KE)GO TO 219                                              
      NCT=1                                                              
      JE=2*JE                                                            
      DZ2P2=H2P2R/JE                                                     
      Z2P2=Y2P2+DZ2P2                                                    
      GO TO 151                                                          
  214 NCT=-1                                                             
      Z2P2=Y2P2-DZ2P2                                                    
      IF(Z2P2-H2P2S)212,151,151                                          
!----Bottom of P1 molecule placement loop                              
  219 IF(VY.GE.VX)GO TO 230                                              
      VX=VY                                                              
      X1P1=Y1P1                                                          
      X2P1=Y2P1                                                          
      X2I=Y2I                                                            
      X2P2=Y2P2                                                          
      X3P2=Y3P2                                                          
      X3I=Y3I                                                            
      IF(JCT)232,231,232                                                 
  230 IF(JCT)232,231,234                                                 
  231 IF(KTR.EQ.2)GO TO 239                                              
      IF(LE.GE.NV)GO TO 232                                              
  233 LE=LE+1                                                            
      Y1P1=Y1P1+DY1P1                                                    
      GO TO 121                                                          
  232 IF(LE.GE.KE)GO TO 239                                              
      JCT=1                                                              
      LE=2*LE                                                            
      DY1P1=H1P1R/LE                                                     
      Y1P1=X1P1+DY1P1                                                    
      GO TO 121                                                          
  234 JCT=-1                                                             
      Y1P1=X1P1-DY1P1                                                    
      IF(Y1P1-H1P1S)232,121,121                                          
!----Bottom of P3 molecule placement loop                              
  239 IF(VX.GE.VH)GO TO 285                                              
      VH=VX                                                              
      H1P3=X1P3                                                          
      H3P2=X3P2                                                          
      H2P1=X2P1                                                          
      H1P1=X1P1                                                          
      H2P2=X2P2                                                          
      H3P3=X3P3                                                          
      H1I=X1I                                                            
      H2I=X2I                                                            
      H3I=X3I                                                            
!----Is this the special case with I molecules in contact?             
  280 IF(KTR.EQ.2)GO TO 501                                              
      IF(ICT)287,286,287                                                 
  285 IF(KTR.EQ.2)GO TO 501                                              
      IF(ICT)287,286,289                                                 
  286 IF(IE.GE.NV)GO TO 287                                              
      IE=IE+1                                                            
      X3P3=X3P3+DX3P3                                                    
      GO TO 105                                                          
  287 IF(IE.GE.KE)GO TO 290                                              
  288 ICT=1                                                              
      IE=2*IE                                                            
      DX3P3=H3P3R/IE                                                     
      X3P3=H3P3+DX3P3                                                    
      GO TO 105                                                          
  289 ICT=-1                                                             
      X3P3=H3P3-DX3P3                                                    
      IF(X3P3-H3P3S)287,105,105                                          
  290 IF(IE.GE.NNE)GO TO 300                                             
      KE=NNE                                                             
      GO TO 288                                                          
  300 H1C=H1P1+H1P3                                                      
      H2C=H2P1+H2P2                                                      
      H3C=H3P2+H3P3                                                      
      IF(H1C.GT.H1I/2.)H1C=H1C-H1I                                       
      IF(H1C.LT.-H1I/2.)H1C=H1C+H1I                                      
      IF(H2C.GT.H2I/2.)H2C=H2C-H2I                                       
      IF(H2C.LT.-H2I/2.)H2C=H2C+H2I                                      
      IF(H3C.GT.H3I/2.)H3C=H3C-H3I                                       
      IF(H3C.LT.-H3I/2.)H3C=H3C+H3I                                      
      IF(IPR.LT.2)GO TO 320                                              
      CALL PAGE(8,8,0)                                                   
      WRITE(NQ,301)ALPHA(KTR),A1,A2,A3,VH,H1I,H2I,H3I                    
  301 FORMAT(1X,A2,3F7.1,'  V =',F8.2,'  AXES',3F8.3)                    
      DO 319 I=1,7                                                       
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
      SM(2)=0.                                                           
      SM(3)=H3P3                                                         
      GO TO 312                                                          
  308 SM(1)=H1P3                                                         
      SM(2)=H2P2                                                         
      SM(3)=H3P2+H3P3                                                    
      GO TO 312                                                          
  309 SM(1)=H1P1+H1P3                                                    
      SM(2)=H2P1                                                         
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
  319 CONTINUE                                                           
      SM(1)=H1C/2.                                                       
      SM(2)=H2C/2.                                                       
      SM(3)=H3C/2.                                                       
      SM(4)=SM(1)/H1I                                                    
      SM(5)=SM(2)/H2I                                                    
      SM(6)=SM(3)/H3I                                                    
      WRITE(NQ,820)SM                                                    
  820 FORMAT(1X,'ORIGIN',2X,3F8.3,';  ',3F8.4)                           
  320 IF(V3.LE.VH)GO TO 400                                                 
      V3=VH                                                                 
      DIJ(NP2)=V3                                                           
      CDIJ(NP2+20)=ALPHA(KTR)                                               
      IF(V3.GE.HT1(8))GO TO 400                                             
      HT1(1)=H1I                                                            
      HT1(2)=0.                                                             
      HT1(3)=H2I                                                            
      HT1(4)=H3I                                                            
      HT1(5)=H1C                                                            
      HT1(6)=H2C                                                            
      HT1(7)=H3C                                                            
      HT1(8)=V3                                                             
      CHT1(9)=ALPHA(KTR)                                                    
      HT1(10)=A1                                                            
      HT1(11)=A2                                                            
      HT1(12)=A3                                                            
      HT1(13)=0.                                                            
      HT1(14)=0.                                                            
      NLOW=15                                                               
      MLOW=KTR+38                                                           
  400 IF(KTR.GE.2)GO TO 1000                                         
      IF(ITR(40).NE.1)GO TO 1000               
      KTR=2                                              
!----Set X1I and X1P3 for I molecules in contact along axis-1          
  500 X1I=Z1I                                                            
      X1IH=Z1IH                                                          
      X1P3=Z1IH                                                          
!----Set flag for negative axis-3 placement                            
      ICT=0                                                              
      VX=1000000.        ! 03-15-06                                                      
      X3P3=H3P3S                                                         
      GO TO 505                                                          
!----Has positive axis-3 placement been calculated?                    
  501 IF(ICT.NE.0)GO TO 300                                              
      ICT=1                                                              
      X3P3=H3P3S+H3P3R                                                   
!----Calculate axis-2 standoff, Y2P1, of a line of P1 molecules - CC   
  505 DZ1P1=X1I/NV                                                       
      Z1P1=-X1IH                                                         
      Y2P1=100.                                                          
      Z2P1=Z2IH                                                          
      JE=1                                                               
      KE=NNE                                                             
      MCT=0                                                              
      ER=ERMT                                                            
      NARK=0                                                             
  515 N=MARK                                                             
      DO 518 I=MARK1,NEND1,NSTP1                                         
      K=IT(I)                                                            
      IF(K.LE.0)GO TO 518                                                
      IF(IT(I+7).NE.0)GO TO 518                                          
      D1=Z1P1+TI(I+5)                                                    
  516 D1=D1+X1I                                                          
      IF(D1.LT.CN(K))GO TO 516                                           
  517 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 517                                           
      IF(D1.LT.-CN(K))GO TO 518                                          
      D=TI(I+1)+D1**2                                                    
      IF(D.GT.CN(K+10))GO TO 517                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=TI(I+3)                                                    
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 517                                                          
  518 CONTINUE                                                           
!----Add A2 molecule contacts                                          
      DO 529 I=1,NMOD                                                    
      DO 528 J=1,NMOD                                                    
      K=IA(I)+IAA(J)                                                     
      D3=X3P3-W(3,J)-W(3,I)                                              
      IF(ABS(D3).GT.CN(K))GO TO 528                                      
      D1=X1P3+Z1P1-W(1,J)-W(1,I)                                         
  526 D1=D1+X1I                                                          
      IF(D1.LT.CN(K))GO TO 526                                           
  527 D1=D1-X1I                                                          
      IF(D1.GT.CN(K))GO TO 527                                           
      IF(D1.LT.-CN(K))GO TO 528                                          
      D=D3**2+D1**2                                                      
      IF(D.GT.CN(K+10))GO TO 527                                         
      IT(N)=K                                                            
      TI(N+1)=D                                                          
      TI(N+2)=W(2,J)-W(2,I)                                              
      N=N+NSTP                                                           
      IF(N.GT.LIMIT)GO TO 902                                            
      GO TO 527                                                          
  528 CONTINUE                                                           
  529 CONTINUE                                                           
      IF(N.LE.MARK)GO TO 920                                             
      NEND=N-1                                                           
      CALL MINHI(Z2P1,Z2IH,Y2P1,Z2I,KE)                                  
  530 IF(KCT)531,531,532                                                
  531 Y1P1=Z1P1                                                         
      Y2P1=Z2P1                                                         
      IF(MCT)536,533,536                                                
  532 IF(MCT)536,533,535                                                
  533 IF(JE.GE.NV)GO TO 534                                             
      JE=JE+1                                                           
      Z1P1=Z1P1+DZ1P1                                                   
      GO TO 515                                                         
  534 MCT=1                                                             
      JE=2*JE                                                           
      DZ1P1=X1I/JE                                                      
      Z1P1=Y1P1+DZ1P1                                                   
      GO TO 515                                                         
  535 MCT=-1                                                            
      Z1P1=Y1P1-DZ1P1                                                   
      IF(Z1P1+X1IH)536,536,515                                           
  536 IF(JE.LT.KE)GO TO 534                                             
      IF(JE.GE.NNE)GO TO 538                                            
      KE=NNE                                                            
      GO TO 534                                                         
  538 Y2I=2.*Y2P1                                                        
      JCT=0                                                              
      IF(IPR.LT.3)GO TO 131                                              
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,539)Y1P1,Y2P1,X3P3                                        
  539 FORMAT(1X,'Y1P1 =',F8.3,', Y2P1 =',F8.3,', X3P3 =',F8.3)           
      GO TO 131                                                          
  902 CALL PAGE(2,2,0)                                                      
      WRITE(NQ,903)MARK1,NSTP1,MARK2,NSTP2,MARK,NSTP,LIMIT                  
  903 FORMAT(26H0STORAGE EXCEEDED BY ZERTA,5X,3(I9,I3),I9)                  
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
      GO TO 999                                                             
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZERTA                                                              
