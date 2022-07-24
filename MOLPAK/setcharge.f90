      SUBROUTINE SETCHARGE              ! last change 04-10-2010 
!
      USE molpakCommonMod
      
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: N_PC=89          ! 04-30-2010 

!      CHARACTER*80 BUF                        
!      CHARACTER*60 HEAD                                                  
!      CHARACTER*4 AN,CDIJ,CHT1                                           
!      CHARACTER*2 AT,CODE                                                 
!-----General program communication and transfer variables               
!     Only 14 JTR and 65 ITR locations are currently in use.             
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

      LOGICAL :: IS_FILE_READ = .false. 
 
      CHARACTER(40) :: PCOMMENT(N_PC)
      
      INTEGER :: I, J, K, NUMBER
      
      REAL  ::  CHARGE(200)                        
      REAL, DIMENSION(N_PC)   :: AMASS, BR 
      REAL, DIMENSION(2,N_PC) :: G92   
   
      DATA AMASS/1.0079, 1.0079,  1.0079, 12.0110, 12.0110, 15.9994,&
     &         14.0067, 14.0067, 14.0067, 14.0067, 15.9994, 79.904, &
     &         18.998 , 15.9994, 15.9994, 14.0067, 15.9994, 14.0067,&
     &         14.0067, 14.0067, 14.0067, 15.9994, 15.9994, 14.0067,&
     &         15.9994, 14.0067, 14.0067, 14.0067, 14.0067, 14.0067,&
     &         14.0067, 14.0067, 14.0067, 12.0110, 10.8110, 32.066, &
     &         15.9994, 14.0067,  1.0079, 14.0067, 15.9994, 14.0067,&
     &          1.0079, 14.0067, 14.0067, 12.0110,  1.0079, 14.0067,&
     &         15.9994, 15.9994, 14.0067, 15.9994,  18.998, 14.0067,&
     &         15.9994, 14.0067, 14.0067, 15.9994,  1.0079, 126.9045,&
     &          1.0079, 14.0067, 15.9994, 1.0079,  15.9994, 14.0067, & ! 12/31/03
     &         12.0110, 1.0079,  14.0067, 4*14.0067, 1.0079,2*14.0067, &
     &         32.0660, 14.0067, 12.0110, 2*15.9994, 14.0067,14.0067, &
     &         15.9994, 14.0067, 1.0079 , 35.4527, 14.0067, 1.0079/    ! 04-30-2010
!
!----Should WMIN potential params be read from POTENTIAL_FILE
      IF (IS_FILE_READ) GO TO 10
         OPEN (UNIT=55, FILE=POTENTIAL_FILE, STATUS='OLD', ERR=1000)
         READ (55,*)        ! skip first title line of potential file
         DO I=1,N_PC   
            READ (55,100) NUMBER, (G92(J,I),J=1,2), BR(I), PCOMMENT(I)
100            FORMAT (I3,3F15.5,A)
            IF (I .EQ. NUMBER) CYCLE
            WRITE (*,105) POTENTIAL_FILE, I, NUMBER
105            FORMAT ('**Inconsistency in potential file ',A,'...'/&
     &                 '  parameter number expected =',I3,', number found =',I3)
            STOP
         ENDDO
         CLOSE (UNIT = 55)
         IS_FILE_READ = .true.
         GO TO 10
!----using DMACRYS, 
1000  OPEN (UNIT=56, FILE='pote.dat', STATUS='OLD', ERR=1002)
      DO I=1,N_PC
       G92(1,I) = 0.0
       G92(2,I) = 0.0
       BR(I)    = 0.0
      END DO
      CLOSE (UNIT = 56)
      GO TO 10 
1002  CLOSE (UNIT = 56)
      OPEN (UNIT=56, FILE='will01.pots', STATUS='OLD', ERR=1001)
      DO I=1,N_PC
       G92(1,I) = 0.0
       G92(2,I) = 0.0
       BR(I)    = 0.0
      END DO
      CLOSE (UNIT = 56)
      GO TO 10
!----Can't find potential file, required for WMIN 
1001  WRITE (*,1005) POTENTIAL_FILE
1005  FORMAT ('**Cannot locate WMIN potential parameter file ',A)
      STOP
10    DO I=1,NMOD                                                      
         K = IDATOM(I)                                                       
         WRITE (NR,30) AT(I), IDATOM(I), G92CHARGE(I), (G92(J,K),J=1,2),&    
     &                 BR(IDATOM(I)), AMASS(IDATOM(I))                      
30          FORMAT (A2,I4,3X,5F9.4)                                             
      ENDDO                                                           
      RETURN                                                             
      END SUBROUTINE SETCHARGE
