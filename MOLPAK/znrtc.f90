      SUBROUTINE ZNRTC 
!
!     For DF and DG.  Space Group Fdd2, programmed only for   
!      molecules with an internal two-fold axis.                                                                                                    
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
!      CHARACTER*2 ALPHA(2),SYM(7)                                        
!      DATA ALPHA/'DF','DG'/                                              
!      DATA SYM/'P1','P2','P3','A1','A2','A3','C '/ 
!                      
      CHARACTER(2):: ALPHA(2) = (/'DF','DG'/)
      CHARACTER(2):: SYM(7) = (/'P1','P2','P3','A1','A2','A3','C '/)
!
      IF(ITR(49).NE.1)GO TO 10                                           
      CALL PAGE(1,1,0)                                                   
      WRITE (NQ,5) ALPHA(1)                                                
    5 FORMAT(1X,'Code ',A2,', Space Group Fdd2, programmed only for', &   
     &' molecules with an internal two-fold axis')                       
   10 IF(ITR(50).NE.1)GO TO 999                                          
      CALL PAGE(1,1,0)                                                   
      WRITE(NQ,5)ALPHA(2)                                                
  999 KILL=1                                                                
 1000 RETURN                                                                
      END SUBROUTINE ZNRTC       
