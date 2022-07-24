      SUBROUTINE PAGE (NEWP, MLIN, MHED)   
                                                    
!                                                                        
!     This subroutine keeps track of the number of lines which have      
!     been printed on the current page and starts a new page with the    
!     designated heading when it is full (60 lines).                     
!         NEWP - A new page is called if fewer than this number of lines 
!                 are available on the current page - should be equal to 
!                 or greater than MLIN.                                  
!         MLIN - Number of lines which are to be printed                 
!         MHED - Number of lines associated with the subheading which    
!                 is to be printed on a new page (usually = 0).          
!               NOTE: PAGE should not be called before the subheading    
!                      is printed by the calling program.                
!                                                                        
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
!!     COMMON ANGLE1(500),ANGLE2(500),ANGLE3(500),VOLUME(500),CODE(500)   
!      COMMON ANGLE1(5000),ANGLE2(5000),ANGLE3(5000),VOLUME(5000),
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
      INTEGER :: NEWP, MLIN, MHED

!-----Tentitively set new-page flag to NO new page was called            
      NEW=0                                                                 
!     LINE is the number of lines still available on the current page.   
!     Is it sufficient to prevent a new page from being called?          
      IF(LINE.GE.NEWP)GO TO 15                                              
!     Increment the page number to be printed.                           
      NPAGE=NPAGE+1                                                         
!     Set the new-page flag - used by the calling program to determine   
!     whether a subheading should be printed.                            
      NEW=1                                                                 
!     Call a new page and print the heading.                             
      WRITE(NQ,10)HEAD,NPAGE                                                
   10 FORMAT(1H1,A60,30X,12HMOLPAK  PAGE,I3)                                
!     Calculate the lines which will be available on the new page after  
!     the subheading is printed.                                         
      LINE=59-MHED                                                          
!     Decrement the number of lines available on the current page.   
   15 LINE=LINE-MLIN                                                        
      RETURN                                                                
      END SUBROUTINE PAGE                                                               
