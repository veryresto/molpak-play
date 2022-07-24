 Program Main
 use nbsLatticeMod, only : ai,alphi,beti,bi,ci,gami,ick011,ick012,ick013,&
                         & ick014,ick021,ick022,ick023,ick024,ick025,&
                         & ick026,ick027,ick028,ick029,ierr1,ierr2,ifcn,&
                         & iunita,iunitb,iunitc,iunitd,ipch,iprob,iprrss,&
                         & ntd,nunk,radian,tol,tolv,v1,v2,v3,voli,u1,u2,&
                         & u3,w1,w2,w3

 implicit none

 character(1) :: asrch(5)
 character(5) :: afcn
 character(5) :: ALM='LM   ', ARSS='RSS  ', ATRANS='TRANS', AINV='INV  '
 character(1) :: IHI='I', IHO='O'

 integer :: i,ictfcn,ierrin,ihead1,ihead2,ihead3,ihead4,isrchi,isrcho
 real :: temrec

 interface
   subroutine ckpt02(iflag2)
   implicit none
    integer :: iflag2
   end subroutine ckpt02
   subroutine determ
   end subroutine determ
   subroutine head0
   end subroutine head0
   subroutine head2
   end subroutine head2
   subroutine head3
   end subroutine head3
   subroutine head4
   end subroutine head4
   subroutine invers(iprinv)
   implicit none
   integer :: iprinv
   end subroutine invers
   subroutine lmpre
   end subroutine lmpre
   subroutine lmseto
   end subroutine lmseto
   subroutine lmsrch
   end subroutine lmsrch
   subroutine outpt1(iflag1)
   implicit none
   integer :: iflag1
   end subroutine outpt1
   subroutine reduce
   end subroutine reduce
   subroutine trans(iprtra)
   implicit none
   integer :: iprtra
   end subroutine trans
   subroutine volume(ax,bx,cx,alphx,betx,gamx,volx)
   implicit none
   real :: ax,bx,cx,alphx,betx,gamx,volx
   end subroutine volume
 
 end interface

!      CHARACTER ASRCH(5),IHI,IHO
!      CHARACTER*5 AFCN,AINV,ALM,ARSS,ATRANS
!      COMMON /CELLI/ AI,BI,CI,ALPHI,BETI,GAMI,VOLI
!      COMMON /CONST1/ RADIAN
!      COMMON /ERR1/ IERR1
!      COMMON /ERR2/ IERR2
!      COMMON /FCN1/ IFCN
!      COMMON /MATR1/ U1,V1,W1,U2,V2,W2,U3,V3,W3
!      COMMON /PRINT1/ IPRRSS
!      COMMON /PROB1/ NTD,IPROB
!      COMMON /PROB2/ IPCH,NUNK
!      COMMON /TOLER1/ TOL,TOLV
!      COMMON /UNIT1/ IUNITA
!      COMMON /UNIT2/ IUNITB
!      COMMON /UNIT3/ IUNITC
!      COMMON /UNIT4/ IUNITD
!**
!      COMMON /CK01/ ICK011,ICK012,ICK013,ICK014
!     COMMON /CK02/ ICK021,ICK022,ICK023,ICK024,ICK025,ICK026,ICK027,
!     $              ICK028,ICK029




!     --- OPENS, CLOSES :
!         IUNITC (= 7)  (ASSOCIATES FILENAMES: NBSI, NBSO)
!         IUNITD (=10)  (ASSOCIATED FILENAME:  NBS10)

!     --- INITIALIZE PROGRAM CONSTANTS
      RADIAN = 180.0/3.14159265

!     --- INITIALIZE VARIABLES
      ICTFCN = 0
      IERRIN = 0
      IHEAD1 = 0
      IHEAD2 = 0
      IHEAD3 = 0
      IHEAD4 = 0
      IPRRSS = 0
      ISRCHI = 0
      ISRCHO = 0
!**
!     --- FOR CHECKING ... READ OR INITIALIZE TEST PRINT/RUN CONTROL
!         VARIABLES
!     READ(IUNITA,1000) ICK011,ICK012,ICK013,ICK014,
!    $    ICK021,ICK022,ICK023,ICK024,ICK025,ICK026,ICK027,ICK028,ICK029
!**
      ICK011 = 0
      ICK012 = 0
      ICK013 = 0
      ICK014 = 0
      ICK021 = 0
      ICK022 = 0
      ICK023 = 0
      ICK024 = 0
      ICK025 = 0
      ICK026 = 0
      ICK027 = 0
      ICK028 = 0
      ICK029 = 0
!**
!     --- WHEN NOT EXECUTING SPECIAL CHECK RUN, GO TO SECTION
!         OF CODE FOR NORMAL EXECUTION OF PROGRAM
      IF(ICK029.NE.1) GO TO 10
!**
!        --- THIS CHECK CODE ALLOWS FOR THE REDUCTION OF A LARGE
!            NUMBER OF CELLS PER RUN (UP TO 99999).  THE REDUCED
!            CELL, TRANSFORMATION MATRIX FROM THE INITIAL CELL TO
!            THE REDUCED CELL, AND THE PROBLEM NUMBER WILL BE
!            WRITTEN TO IUNITD.  IUNITD MUST BE PRE-ASSIGNED TO
!            THE COMPUTER RUN. THE PRINT OUT IS LIMITED TO THE
!            PROBLEM NUMBER FOLLOWED BY ANY WARNING AND/OR
!            ERROR MESSAGES.

!        --- INITIALIZE VARIABLES
         IFCN = 2
         IPCH = 0
         IPRRSS = 1
         NUNK = 0

!        --- READ NUMBER OF PROBLEMS
         READ(IUNITA,980) NTD

!        --- WRITE SPECIAL HEADING, NUMBER OF PROBLEMS
         CALL HEAD0
         CALL CKPT02(12)

!        --- CALCULATE REDUCED CELLS
         CALL REDUCE
         GO TO 800
   10 CONTINUE

!     --- WRITE GENERAL HEADING FOR PROGRAM
      CALL HEAD0

!     --- STARTING POINT FOR EACH PROGRAM FUNCTION TO BE PROCESSED
   50 CONTINUE
      ICTFCN = ICTFCN + 1

!     --- READ PROGRAM CONTROL CARD
!         (TYPE OF PROGRAM FUNCTION, NUMBER OF PROBLEMS)
      READ(IUNITA,1100,END=800) AFCN, NTD

!     --- ASSIGN A VARIABLE (IFCN) TO INDICATE THE TYPE OF
!         PROGRAM FUNCTION
      IFCN = 0
      IF(AFCN.EQ.ALM)    IFCN = 1
      IF(AFCN.EQ.ARSS)   IFCN = 2
      IF(AFCN.EQ.ATRANS) IFCN = 3
      IF(AFCN.EQ.AINV)   IFCN = 4

!     --- CHECK FOR VALID PROGRAM FUNCTION (IF NOT VALID,
!         WRITE ERROR MESSAGE AND STOP EXECUTION OF PROGRAM)
      IF(IFCN.LT.1.OR.IFCN.GT.4) GO TO 700

!     --- GO TO SECTION OF CODE FOR SPECIFIED FUNCTION
      GO TO (100,200,300,400) IFCN

  100 CONTINUE


!     *** LM FUNCTION ***

!     --- IDENTIFY UNKNOWN CELLS BASED ON LATTICE
!         MATCHING WITH THE NBS CRYSTAL DATA FILE

!     --- CHECK THAT THIS IS THE FIRST LM FUNCTION TO BE PROCESSED
      IF(IHEAD1.EQ.0) GO TO 120

!        --- ONLY ONE LM OPERATION PER RUN IS ALLOWED ... SKIP EXPECTED
!            INPUT RECORDS AND PROCEED TO NEXT PROGRAM FUNCTION
         WRITE(IUNITB,1200)
         DO 110 I = 1,NTD+1
            READ(IUNITA,1300) TEMREC
  110    CONTINUE
         GO TO 50
  120 CONTINUE

!     --- READ LM PARAMETER CARD
      READ(IUNITA,1400) TOL, TOLV, ASRCH, IPRRSS

!     --- ASSIGN VARIABLES TO INDICATE WHICH 1-LINE DATA BASES ARE
!         TO BE SEARCHED (DEFAULT = INORGANIC)
      DO 130 I = 1,5
         IF(ASRCH(I).EQ.IHI) ISRCHI = 1
         IF(ASRCH(I).EQ.IHO) ISRCHO = 1
  130 CONTINUE
      IF(ISRCHI.EQ.0.AND.ISRCHO.EQ.0) ISRCHI = 1

!     --- ALWAYS WRITE IUNITD, CONTROL WRITE OF RSS HEADING
      IPCH = 1
      IF(ICTFCN.GT.1) WRITE(IUNITB,5000)
      IF(IPRRSS.EQ.0) CALL HEAD2

!     --- CHECK NUMBER OF PROBLEMS INPUT ... IF NECESSARY, RESET TO 1
      IF(NTD.GE.1.AND.NTD.LE.20) GO TO 140
         IERRIN = NTD
         WRITE(IUNITB,1500)
         NTD = 1
  140 CONTINUE

!     --- WRITE NUMBER OF INDEPENDENT PROBLEMS TO STUDY
      CALL OUTPT1(1)

!     --- OUTPUT FOR LATTICE MATCHING OF REDUCED CELLS
      CALL OUTPT1(2)

!     --- INITIALIZE VARIABLE
      NUNK = 0

!     --- CALCULATE REDUCED AND DERIVATIVE CELLS
      OPEN(UNIT=IUNITD,FILE='NBS10',STATUS='NEW')
      CALL REDUCE
      CLOSE(IUNITD)

!     --- IDENTIFICATION OF REDUCED CELL(S) BASED ON LATTICE MATCHING
!         (FIRST PREPARE INPUT CELLS)
      CALL LMPRE

!     --- SEARCH INORGANIC 1-LINE DATA BASE IF INDICATED
      IF(ISRCHI.EQ.0) GO TO 150

!        --- WRITE HEADING, THEN SEARCH DATA BASE
         WRITE(IUNITB,1600)
         OPEN(UNIT=IUNITC,FILE='NBSI',STATUS='OLD',ACCESS='DIRECT',&
     &        FORM='FORMATTED',RECL=132)
         CALL LMSRCH
         CLOSE(IUNITC)
  150 CONTINUE

!     --- SEARCH ORGANIC 1-LINE DATA BASE IF INDICATED
      IF(ISRCHO.EQ.0) GO TO 160

!        --- WRITE HEADING
         WRITE(IUNITB,1700)

!        --- RESET VARIABLES, POINTERS FOR ORGANIC SEARCH
         CALL LMSETO

!        --- ORGANIC SEARCH
         OPEN(UNIT=IUNITC,FILE='NBSO',STATUS='OLD',ACCESS='DIRECT',&
     &        FORM='FORMATTED',RECL=132)
         CALL LMSRCH
         CLOSE(IUNITC)
  160 CONTINUE

!     --- IF INPUT ERROR, SKIP EXPECTED INPUT RECORDS AND PROCEED TO
!         NEXT PROGRAM FUNCTION
      IF(IERRIN.EQ.0) GO TO 180
         DO 170 I = 1,IERRIN-1
            READ(IUNITA,1300) TEMREC
  170    CONTINUE
         IERRIN = 0
  180 CONTINUE

      IHEAD1 = IHEAD1 + 1
      IPRRSS = 0
      GO TO 50

  200 CONTINUE


!     *** RSS FUNCTION ***

!     --- NEVER WRITE IUNITD, CONTROL WRITE OF RSS HEADING
      IPCH = 0
      IF(ICTFCN.GT.1) WRITE(IUNITB,5000)
      IF(IHEAD2.EQ.0) CALL HEAD2

!     --- CHECK NUMBER OF PROBLEMS INPUT ... IF NECESSARY, RESET TO 1
      IF(NTD.GE.1.AND.NTD.LE.20) GO TO 210
         IERRIN = NTD
         WRITE(IUNITB,1500)
         NTD = 1
  210 CONTINUE

!     --- WRITE NUMBER OF INDEPENDENT PROBLEMS TO STUDY
      CALL OUTPT1(1)

!     --- CALCULATE REDUCED AND DERIVATIVE CELLS
      CALL REDUCE

!     --- IF INPUT ERROR, SKIP EXPECTED INPUT RECORDS AND PROCEED TO
!         NEXT PROGRAM FUNCTION
      IF(IERRIN.EQ.0) GO TO 230
         DO 220 I = 1,IERRIN-1
            READ(IUNITA,1300) TEMREC
  220    CONTINUE
         IERRIN = 0
  230 CONTINUE

      IHEAD2 = IHEAD2 + 1
      GO TO 50

  300 CONTINUE


!     *** TRANS FUNCTION ***

!     --- CONTROL WRITE OF TRANS HEADINGS
      IF(ICTFCN.GT.1) WRITE(IUNITB,5000)
      IF(IHEAD3.EQ.0) CALL HEAD3
      CALL OUTPT1(14)

!     --- INITIALIZE VARIABLES
      IPCH  = 0
      IPROB = 0

!     --- TRANSFORM EACH INPUT CELL BY SPECIFIED MATRIX
      DO 320 I = 1,NTD
         IPROB = I

!        --- WRITE SUB-HEADING
         CALL OUTPT1(15)

!        --- READ TRANSFORMATION MATRIX AND CELL
         READ(IUNITA,1800) U1, V1, W1, U2, V2, W2, U3, V3, W3
         READ(IUNITA,1900) AI, BI, CI, ALPHI, BETI, GAMI

!        --- CALCULATE VOLUME OF THE INPUT CELL
!            (AND CHECK FOR VALID INPUT DATA)
         CALL VOLUME(AI,BI,CI,ALPHI,BETI,GAMI,VOLI)

!        --- CHECK THAT INPUT CELL IS VALID
!            (INDICATED BY ERROR FLAG SET IN *VOLUME*) ...
!            IF NOT, WRITE ERROR MESSAGE AND GO TO NEXT PROBLEM
         IF(IERR1.EQ.0) GO TO 310
            WRITE(IUNITB,2000)
            GO TO 320
  310    CONTINUE

!        --- TRANSFORM CELL AND WRITE
         CALL TRANS(1)
  320 CONTINUE
      IHEAD3 = IHEAD3 + 1
      GO TO 50

  400 CONTINUE


!     *** INV FUNCTION ***

!     --- CONTROL WRITE OF INV HEADINGS
      IF(ICTFCN.GT.1) WRITE(IUNITB,5000)
      IF(IHEAD4.EQ.0) CALL HEAD4
      CALL OUTPT1(16)

!      --- INITIALIZE VARIABLES
      IPROB = 0

!     --- CALCULATE INVERSE FOR EACH 3X3 MATRIX
      DO 420 I = 1,NTD
         IPROB = I
         READ(IUNITA,1800) U1, V1, W1, U2, V2, W2, U3, V3, W3
         CALL OUTPT1(17)
         CALL OUTPT1(11)

!        --- CALCULATE DETERMINANT AND CHECK FOR VALID INPUT MATRIX
         CALL DETERM
         IF(IERR2.EQ.0) GO TO 410

!           --- INVALID MATRIX (DETERMINANT IS ZERO, >= 100, OR
!               <= -100), WRITE ERROR MESSAGE AND GO TO NEXT PROBLEM
            WRITE(IUNITB,2100)
            GO TO 420
  410    CONTINUE
         CALL INVERS(1)
  420 CONTINUE
      IHEAD4 = IHEAD4 + 1
      GO TO 50
  700 CONTINUE

!     --- ERROR IN SPECIFIED PROGRAM FUNCTION
      WRITE(IUNITB,4000)
  800 CONTINUE
      STOP
  980 FORMAT(5X,I5)
!**
!1000 FORMAT(4I1,6X,9I1,1X)
 1100 FORMAT(A5,3X,I2)
 1200 FORMAT(1H1////1X,'INPUT ERROR ... Only one Lattice Matching operat&
     &ion per run is allowed.')
 1300 FORMAT(A80)
 1400 FORMAT(2F10.2,5X,5A1,9X,I1)
 1500 FORMAT(///1X,'INPUT ERROR ... Check columns 9-10 on Program Contro&
     &l Line.'/1X,16X,'Number of independent problems set to 1.')
 1600 FORMAT(1H1,41X,'INORGANIC')
 1700 FORMAT(1H1,42X,'ORGANIC')
 1800 FORMAT(9F8.2)
 1900 FORMAT(10X,6F10.5)
 2000 FORMAT(/1X,'INPUT ERROR ... Input cell has illegal cell parameter(&
     &s) and/or cell volume.'/)
 2100 FORMAT(/1X,'INPUT ERROR ... Invalid matrix, check determinant.'/)
 4000 FORMAT(1H1////1X,'INPUT ERROR ... Invalid program function was spe&
     &cified.')
 5000 FORMAT(1H1)
 end Program Main
 subroutine center

! **********************************************************************
!  Purpose: Lattice type of input cell is identified and then a call to
!           subroutine outpt1 results in the writing of the elements of
!           the matrix that transforms the input cell to the primitive 
!           cell to file unit number iunitb
!  
!  Algorithm:
!
!  Parameters:
!
!  Variables: icntr specifies input cell lattice type
!             icntr=1 denotes primitive lattice
!             icntr=2 denotes a-centered lattice
!             icntr=3 denotes b-centered lattice
!             icntr=4 denotes c-centered lattice
!             icntr=5 denotes f-centered lattice
!             icntr=6 denotes i-centered lattice
!             icntr=7 denotes rhombohdral cell set on hexagonal axes
!             
!             u1, v1, w1, u2, v2, w2, u3, v3, w2 are the elements of
!             the transformation matrix
! 
! **********************************************************************
!
 use nbsLatticeMod, only:  icntr1
 use nbsLatticeMod, only:  u1,v1,w1,u2,v2,w2,u3,v3,w3

 implicit none

 character( 6 ), parameter :: procedureName = 'center'

 select case (icntr1)

   case(1)      !     --- primitive lattice
      u1 = 1.0
      v1 = 0.0
      w1 = 0.0
      u2 = 0.0
      v2 = 1.0
      w2 = 0.0
      u3 = 0.0
      v3 = 0.0
      w3 = 1.0
      call outpt1(10) 
      return

   case(2)      !     --- a-centered lattice
      u1 =  0.0
      v1 =  0.5
      w1 = -0.5
      u2 =  0.0
      v2 =  0.5
      w2 =  0.5
      u3 =  1.0
      v3 =  0.0
      w3 =  0.0
      call outpt1(10) 
      return

   case(3)      !     --- b-centered lattice
      u1 = -0.5
      v1 =  0.0
      w1 =  0.5
      u2 =  0.5
      v2 =  0.0
      w2 =  0.5
      u3 =  0.0
      v3 =  1.0
      w3 =  0.0
      call outpt1(10) 
      return

   case(4)      !     --- c-centered lattice
      u1 =  0.5
      v1 = -0.5
      w1 =  0.0
      u2 =  0.5
      v2 =  0.5
      w2 =  0.0
      u3 =  0.0
      v3 =  0.0
      w3 =  1.0
      call outpt1(10) 
      return

   case(5)      !     --- f-centered lattice
      u1 = 0.5
      v1 = 0.5
      w1 = 0.0
      u2 = 0.0
      v2 = 0.5
      w2 = 0.5
      u3 = 0.5
      v3 = 0.0
      w3 = 0.5
      call outpt1(10) 
      return

   case(6)      !     --- i-centered lattice
      u1 =  0.5
      v1 =  0.5
      w1 = -0.5
      u2 = -0.5
      v2 =  0.5
      w2 =  0.5
      u3 =  0.5
      v3 = -0.5
      w3 =  0.5
      call outpt1(10) 
      return

   case(7)      !     --- rhombohedral cell set on hexagonal axes
      u1 =  1.0/3.0
      v1 = -1.0/3.0
      w1 = -1.0/3.0
      u2 = -2.0/3.0
      v2 = -1.0/3.0
      w2 = -1.0/3.0
      u3 =  1.0/3.0
      v3 =  2.0/3.0
      w3 = -1.0/3.0
      call outpt1(10) 

 end select

 return
 end subroutine center
      SUBROUTINE CKPT02(IFLAG2)

 use nbsLatticeMod, only : iprob,itype,iunitb,ntd,s11,s22,s33,s23,s13,&
                         & s12,u1,u2,u3,v1,v2,v3,w1,w2,w3

 implicit none

 integer, intent(inout) :: iflag2
 integer :: icond

!     --- IF THE CALL ARGUMENT (IFLAG2) IS GREATER THAN OR EQUAL TO
!         31, IT MUST BE DECODED TO OBTAIN THE NUMBER OF THE SPECIAL
!         CONDITION AND THE 'CALL ARGUMENT' FOR THE OUTPUT
      ICOND = 0
      IF(IFLAG2.GE.31) ICOND = IFLAG2 - 30
      IF(IFLAG2.GE.31) IFLAG2 = 10

!     --- GO TO APPROPRIATE SECTION OF OUTPUT
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130) IFLAG2
   10 CONTINUE

!     *** CALLED FROM 'MNCOND'

!     --- WRITE CELL MATRIX AND TRANSFORMATION MATRIX AT BEGINNING OF
!         SUBROUTINE 'MNCOND'
      WRITE(IUNITB,1100) IFLAG2, S11, S22, S33, S23, S13, S12,&
     &                   U1, V1, W1, U2, V2, W2, U3, V3, W3
      GO TO 900
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE

!     *** CALLED FROM 'MNCOND'

!     --- WRITE PRIOR TO CALLS TO 'SHORTV'
!         IFLAG2 = 2    ...   PRIOR TO FIRST CALLS
!         IFLAG2 = 3    ...   PRIOR TO SECOND CALLS
!         IFLAG2 = 4    ...   PRIOR TO THIRD CALLS
      WRITE(IUNITB,1100) IFLAG2
      GO TO 900
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE

!     *** CALLED FROM 'MNCOND'

!     --- WRITE CELL MATRIX AND TRANSFORMATION MATRIX AT
!         INTERMEDIATE STAGES
!         IFLAG2 = 5 ... AFTER CELL EDGES A AND B HAVE BEEN INTERCHANGED
!         IFLAG2 = 6 ... AFTER CELL EDGES B AND C HAVE BEEN INTERCHANGED
!         IFLAG2 = 7 ... AFTER CALCULATING THE BODY DIAGONAL OF THE CELL
      WRITE(IUNITB,1100) IFLAG2, S11, S22, S33, S23, S13, S12,&
     &                   U1, V1, W1, U2, V2, W2, U3, V3, W3
      GO TO 900
   80 CONTINUE

!     *** CALLED FROM 'SHORTV'

!     --- WRITE CELL MATRIX AND TRANSFORMATION MATRIX AT END
!         OF SUBROUTINE 'SHORTV'
      WRITE(IUNITB,1200) S11, S22, S33, S23, S13, S12,&
     &                   U1, V1, W1, U2, V2, W2, U3, V3, W3
      GO TO 900
   90 CONTINUE

!     *** CALLED FROM 'NORMAL'

!     --- WRITE CELL TYPE, CELL MATRIX, TRANSFORMATION MATRIX
      WRITE(IUNITB,1300) ITYPE, S11, S22, S33, S23, S13, S12,&
     &                   U1, V1, W1, U2, V2, W2, U3, V3, W3
      GO TO 900
  100 CONTINUE

!     *** CALLED FROM 'SPCON2'

!     --- WRITE SPECIAL CONDITION NUMBER, CELL TYPE, CELL MATRIX,
!         TRANSFORMATION MATRIX
      WRITE(IUNITB,1400) ICOND, ITYPE, S11, S22, S33, S23, S13, S12,&
     &                   U1, V1, W1, U2, V2, W2, U3, V3, W3
      GO TO 900
  110 CONTINUE

!     *** CALLED FROM 'SPCON2'

!     --- WRITE CELL MATRIX AND TRANSFORMATION MATRIX AFTER THE
!         TRANSFORMATION TO SATISFY THE SPECIAL CONDITION HAS BEEN
!         APPLIED
      WRITE(IUNITB,1500) S11, S22, S33, S23, S13, S12,&
     &                   U1, V1, W1, U2, V2, W2, U3, V3, W3
      GO TO 900
  120 CONTINUE

!     *** CALLED FROM 'MAIN' FOR SPECIAL CHECK RUN ONLY

!     --- WRITE HEADING FOR CHECK RUN, NUMBER OF PROBLEMS
      WRITE(IUNITB,2000)
      WRITE(IUNITB,2100) NTD
      GO TO 900
  130 CONTINUE

!     *** CALLED FROM 'OUTPT1' FOR CHECK RUN ONLY

!     --- WRITE PROBLEM NUMBER
      WRITE(IUNITB,2200) IPROB
  900 CONTINUE
      RETURN
 1100 FORMAT(1X,'MNCOND',I2,1X,6F12.6,3X,9F5.2)
 1200 FORMAT(1X,'SHORTV  ',1X,6F12.6,3X,9F5.2)
 1300 FORMAT(1X,'NORMAL',1X,I1,1X,6F12.6,3X,9F5.2)
 1400 FORMAT(1X,'SPCON2',2I1,1X,6F12.6,3X,9F5.2)
 1500 FORMAT(1X,'SPCON2  ',1X,6F12.6,3X,9F5.2)
 2000 FORMAT(///1X,'This mode of operation allows for the reduction of a&
     & large number'     /1X,'of cells (up to 99999) in a single run.  I&
     &n this mode of operation,'     /1X,'the reduced cell, the transfor&
     &mation matrix from the initial cell to'     /1X,'the reduced cell,&
     & and the problem number will be written on unit 10'     /1X,'(Form&
     &at = 6F10.6,2X,9F7.2,2X,I5).  Unit 10 must be pre-assigned to the'&
     &     /1X,'computer run.  The printed output is limited to the prob&
     &lem number'     /1X,'followed by any warning and/or error messages&
     &.')
 2100 FORMAT(//1X,'Number of cells to be reduced = ',I5//////)
 2200 FORMAT(1X,I5)
 end subroutine ckpt02
      SUBROUTINE DERIV
 use nbsLatticeMod, only : ideriv,idet,ifinis,ilatno,irss,iqmatf, &
                         & isq11,isq12,isq13,isq22,isq23,isq33,&
                         & ui1,ui2,ui3,vi1,vi2,vi3,wi1,wi2,wi3,&
                         & u1,u2,u3,v1,v2,v3,w1,w2,w3
 implicit none

 integer :: idpre
 real :: det

 interface 

  subroutine invers(iprinv)
    implicit none
    integer :: iprinv
  end subroutine invers

 end interface

!     --- INITIALIZE VARIABLES FOR EACH NEW PROBLEM
      IF(IDERIV.NE.1) GO TO 100
         IDET = 0
         IFINIS = 0
         ILATNO = 0
  100 CONTINUE

!     --- INCREMENT DERIVATIVE LATTICE NUMBER, AND COUNT THE NUMBER
!         OF UPPER TRIANGULAR MATRICES SELECTED (IDERIV = IQMATF WHEN
!         LAST MATRIX IS SELECTED)
      ILATNO = ILATNO + 1
      IDERIV = IDERIV + 1

!     --- SELECT AN UPPER TRIANGULAR MATRIX
      U1 = ISQ11(IDERIV-1)
      V1 = ISQ12(IDERIV-1)
      W1 = ISQ13(IDERIV-1)
      U2 = 0.0
      V2 = ISQ22(IDERIV-1)
      W2 = ISQ23(IDERIV-1)
      U3 = 0.0
      V3 = 0.0
      W3 = ISQ33(IDERIV-1)

!     --- CALCULATE DETERMINANT FOR UPPER TRIANGULAR MATRIX
!         (SAVE PREVIOUS VALUE FOR LATER COMPARISON)
      IDPRE = IDET
      DET  = U1*V2*W3 + V1*W2*U3 + W1*U2*V3 - U3*V2*W1 - V3*W2*U1 -&
     &       W3*U2*V1
      IDET = NINT(DET)

!     --- DO NOT CALCULATE THE TRANSPOSE OF THE INVERSE FOR THE UPPER
!         TRIANGULAR MATRIX WHEN SUPERCELLS ARE BEING CALCULATED
!         (WHEN BOTH SUPERCELL AND SUBCELL CALCULATIONS ARE TO BE
!         DONE, SUPERCELLS ARE CALCULATED FIRST)
      IF(IRSS.EQ.1.OR.IRSS.EQ.3) GO TO 300

!        --- CALCULATE THE TRANSPOSE OF THE INVERSE FOR AN UPPER
!            TRIANGULAR MATRIX (USED FOR CALCULATING SUBCELLS)
         CALL INVERS(0)
         U1 = UI1
         V1 = UI2
         W1 = UI3
         U2 = VI1
         V2 = VI2
         W2 = VI3
         U3 = WI1
         V3 = WI2
         W3 = WI3
  300 CONTINUE

!     --- CHECK FOR A CHANGE IN THE VALUE OF THE DETERMINANT OF
!         THE UPPER TRIANGULAR MATRIX ... ENABLES APPROPRIATE
!         SUB-HEADING TO BE WRITTEN
      IF(IDET.EQ.IDPRE) GO TO 400

!        --- CHANGE IN DELTA ENCOUNTERED ... WRITE SUB-HEADING FOR
!            DERIVATIVE LATTICES
         ILATNO = 1
         IF(IRSS.NE.2) CALL OUTPT1(7)
         IF(IRSS.EQ.2) CALL OUTPT1(8)
  400 CONTINUE

!     --- WRITE DERIVATIVE LATTICE NUMBER AND TRANSFORMATION MATRIX
!         TO DERIVATIVE CELL
      CALL OUTPT1(9)
      CALL OUTPT1(10)

!     --- ONLY WHEN THE LAST UPPER TRIANGULAR MATRIX HAS BEEN SELECTED,
!         INITIALIZE VARIABLES FOR THE NEXT PROBLEM, OR FOR THE SUBCELL
!         CALCULATIONS WHEN BOTH SUPERCELL AND SUBCELL CALCULATIONS ARE
!         ARE TO BE DONE
      IF(IDERIV-1.LT.IQMATF) GO TO 700

!        --- LAST UPPER TRIANGULAR MATRIX, CHECK IF FINISHED PROBLEM
         IF(IRSS.NE.3) IFINIS = 1
         IF(IRSS.NE.3) GO TO 700

!           --- FINISHED SUPERCELLS, INITIALIZE VARIABLES FOR SUBCELL
!               CALCULATION ON NEXT CALL TO SUBROUTINE 'DERIV'
            IRSS = 2
            IDET = 0
            IDERIV = 1
            ILATNO = 0
  700 CONTINUE
      RETURN
      END
      subroutine determ

 use nbsLatticeMod, only : ierr2,ifcn,iunitb,u1,u2,u3,&
                         & v1,v2,v3,w1,w2,w3
 implicit none

 real :: det

!     --- INITIALIZE VARIABLES
      DET = 0
      IERR2 = 0

!     --- CALCULATE DETERMINANT
      DET = U1*V2*W3 + V1*W2*U3 + W1*U2*V3 - U3*V2*W1 - V3*W2*U1 -&
     &      W3*U2*V1

!     --- CHECK VALUE OF DETERMINANT FOR ERROR OR WARNING
!         SITUATIONS (DIFFERENT ACTIONS ARE REQUIRED FOR
!         DIFFERENT PROGRAM FUNCTIONS)
      IF(.NOT.(IFCN.EQ.1.OR.IFCN.EQ.2)) GO TO 100

!        --- LM,RSS PROGRAM FUNCTIONS ARE BEING EXECUTED

!        --- FOR ZERO OR NEGATIVE DETERMINANT, STOP PROGRAM
!            EXECUTION (SHOULD NOT OCCUR)
         IF(DET.LE.0.01) STOP&
     &      '*DETERM* ERROR ... Determinant is zero or negative.'

!        --- FOR DETERMINANTS GREATER THAN 9, STOP PROGRAM
!            EXECUTION (SHOULD NOT OCCUR)
         IF(DET.GE.9.01) STOP&
     &      '*DETERM* ERROR ... Determinant is greater than nine.'
         GO TO 400
  100 CONTINUE
      IF(.NOT.(IFCN.EQ.3.OR.IFCN.EQ.4)) GO TO 200

!        --- TRANS,INV PROGRAM FUNCTIONS ARE BEING EXECUTED

!        --- FOR NEGATIVE DETERMINANT, WRITE WARNING AND CONTINUE
         IF(DET.LT.-0.01) WRITE(IUNITB,1000)

!        --- FOR ZERO DETERMINANT, SET ERROR FLAG
!            (DO NOT HAVE A VALID MATRIX)
         IF(DET.GE.-0.01.AND.DET.LE.0.01) IERR2 = 1

!        --- FOR DETERMINANTS >= 100 OR <= -100, SET ERROR FLAG
!            (DO NOT HAVE A VALID MATRIX)
         IF(ABS(DET).GE.100.0) IERR2 = 1
         GO TO 400
  200 CONTINUE

!     --- STOP PROGRAM EXECUTION SINCE A VALID PROGRAM FUNCTION
!         HAS NOT BEEN ENCOUNTERED (SHOULD NOT OCCUR)
      STOP '*DETERM* ERROR ... Invalid progam function.'
  400 CONTINUE
      RETURN
 1000 FORMAT(/1X,'*DETERM* WARNING ... Determinant of matrix is &
     &negative.'/)
 end subroutine determ
 subroutine dot (x1,x2,x3,y1,y2,y3,xi,yi,zi,ar,br,gr,dotxy)

 implicit none

 real, intent(in) :: ar,br,gr,xi,x1,x2,x3,yi,y1,y2,y3,zi
 real, intent(out) :: dotxy
 

!     --- vector x = (x1a + x2b + x3c)
!         vector y = (y1a + y2b + y3c) , where (a,b,c) are (xi,yi,zi)

 dotxy = x1*y1*(xi**2) + x2*y2*(yi**2) + x3*y3*(zi**2) +&
       &(x3*y1+x1*y3)*zi*xi*cos(br) + (x3*y2+x2*y3)*zi*yi*cos(ar) +&
       &(x1*y2+x2*y1)*xi*yi*cos(gr)

 return
 end subroutine dot
 subroutine head0

 use nbsLatticemod, only: iunitb

 implicit none

!     --- write heading for output

      write(iunitb,1100)
      write(iunitb,1110)
      write(iunitb,1120)
      write(iunitb,1130)
      write(iunitb,1140)
      write(iunitb,1150)
      write(iunitb,1160)
      write(iunitb,2000)
      return
 1100 FORMAT(1H1,////T45,'***  NBS*LATTICE  ***'/)
 1110 FORMAT(1X,T35,'A PROGRAM TO ANALYZE LATTICE RELATIONSHIPS')
 1120 FORMAT(1X,T45,'Version of Summer, 1985'/)
 1130 FORMAT(1X,T44,'NBS Crystal Data Center')
 1140 FORMAT(1X,T42,'National Bureau of Standards')
 1150 FORMAT(1X,T42,'Reactor Radiation Division')
 1160 FORMAT(1X,T44,'Gaithersburg, MD  20899')
 2000 FORMAT(////)
 end subroutine head0
 subroutine head2
 use nbsLatticeMod, only:  iunitb

 implicit none

!     --- WRITE HEADING FOR OUTPUT
      WRITE(IUNITB,1100)
      WRITE(IUNITB,1110)
      WRITE(IUNITB,1120)
      WRITE(IUNITB,1130)
      WRITE(IUNITB,1140)
      WRITE(IUNITB,1150)
      WRITE(IUNITB,1160)
      WRITE(IUNITB,1170)
      WRITE(IUNITB,1180)
      WRITE(IUNITB,1190)
      WRITE(IUNITB,1200)
      WRITE(IUNITB,1210)
      WRITE(IUNITB,1220)
      WRITE(IUNITB,1230)
      WRITE(IUNITB,1240)
      WRITE(IUNITB,1250)
      WRITE(IUNITB,1260)
      WRITE(IUNITB,1270)
      WRITE(IUNITB,2000)
      RETURN
 1100 FORMAT(2X,25X,'** REDUCTION AND DERIVATIVE LATTICE **'/)
 1110 FORMAT(2X,'These calculations fall into two categories:'/)
 1120 FORMAT(2X,' I. Reduction of an input cell.'/)
 1130 FORMAT(2X,'    CELL 1 = Input cell.  This cell may be primitive or&
     & centered (A,B,C,I,F,RR,RH).')
 1140 FORMAT(2X,'    CELL 2 = Reduced primitive cell of the lattice.')
 1150 FORMAT(2X,'    T 1    = A matrix that transforms CELL 1 to a primi&
     &tive cell of the lattice.')
 1160 FORMAT(2X,'    T 2    = A matrix that transforms CELL 1 to CELL 2.&
     &'//)
 1170 FORMAT(2X,'II. Calculation and reduction of a series of derivative&
     & supercells and/or subcells.')
 1180 FORMAT(2X,'    These derivative cells are calculated from the redu&
     &ced cell of the lattice')
 1190 FORMAT(2X,'    (i.e. to carry out the Type II calculation, the pro&
     &gram first carries out'/2X,'    the Type I calculation).'/)
 1200 FORMAT(2X,'    CELL 1 = Reduced primitive cell (i.e. CELL 2 from P&
     &art I).')
 1210 FORMAT(2X,'    CELL 2 = Reduced supercell or subcell.')
 1220 FORMAT(2X,'    T 1    = A matrix that transforms CELL 1 to a super&
     &cell or subcell of the lattice.')
 1230 FORMAT(2X,'    T 2    = A matrix that transforms CELL 1 to CELL 2.&
     &'/)
 1240 FORMAT(2X,'For CELL 1 and CELL 2, the output parameters given are:&
     & a, b, c, alpha, beta, gamma')
 1250 FORMAT(2X,'    and volume.  Cell edges are in angstroms and angles&
     & in degrees.'/)
 1260 FORMAT(2X,'The reduced cell matrix is of the form:           a.a&
     & b.b   c.c')
 1270 FORMAT(1X,51X,'b.c   a.c   a.b     .')
 2000 FORMAT(1H1)
 end subroutine head2
      subroutine head3

 use nbsLatticeMod, only: iunitb

 implicit none

!     --- WRITE HEADING FOR OUTPUT
      WRITE(IUNITB,1100)
      WRITE(IUNITB,1110)
      WRITE(IUNITB,1120)
      WRITE(IUNITB,1130)
      WRITE(IUNITB,1140)
      WRITE(IUNITB,1150)
      WRITE(IUNITB,1160)
      WRITE(IUNITB,1170)
      WRITE(IUNITB,1180)
      WRITE(IUNITB,1190)
      WRITE(IUNITB,2000)
      RETURN
 1100 FORMAT(19X,'** CELL TRANSFORMATION **'/)
 1110 FORMAT(5X,'    CELL 1   =   Input cell.')
 1120 FORMAT(5X,'    CELL 2   =   Transformed cell.')
 1130 FORMAT(5X,'    T 2      =   Input transformation matrix.')
 1140 FORMAT(5X,'    T 2 INV  =   Inverse matrix for T 2.'/)
 1150 FORMAT(5X,'For CELL 1 and CELL 2, the output parameters given are:&
     &')
 1160 FORMAT(5X,'    a, b, c, alpha, beta, gamma and volume.  Cell edges&
     &')
 1170 FORMAT(5X,'    are in angstroms and angles in degrees.'/)
 1180 FORMAT(5X,'The cell matrix is of the form:      a.a   b.b   c.c')
 1190 FORMAT(34X,'        b.c   a.c   a.b     .')
 2000 FORMAT(//)
 end subroutine head3
      SUBROUTINE HEAD4
 use nbsLatticeMod, only:  iunitb

!     --- WRITE HEADING FOR OUTPUT
      WRITE(IUNITB,1100)
      WRITE(IUNITB,1110)
      WRITE(IUNITB,1120)
      WRITE(IUNITB,2000)
      RETURN
 1100 FORMAT(21X,'** MATRIX INVERSION **'/)
 1110 FORMAT(7X,'    T 2      =   Input transformation matrix.')
 1120 FORMAT(7X,'    T 2 INV  =   Inverse matrix for T 2.')
 2000 FORMAT(//)
 end subroutine head4
      SUBROUTINE INVERS(IPRINV)
 use nbsLatticeMod, only : ui1,ui2,ui3,u1,u2,u3,vi1,vi2,vi3, &
                         & v1,v2,v3,wi1,wi2,wi3,w1,w2,w3
 implicit none

 integer, intent(in) :: iprinv
 real :: det

 interface
   subroutine outpt1(iflag)
   implicit none
   integer, intent(in) :: iflag
   end subroutine outpt1
 end interface

!     --- CALCULATE DETERMINANT
   DET = U1*V2*W3 + V1*W2*U3 + W1*U2*V3 - U3*V2*W1 - V3*W2*U1 -&
     &      W3*U2*V1

!     --- CALCULATE INVERSE
      UI1 =  (V2*W3 - V3*W2)/DET
      VI1 = -(V1*W3 - V3*W1)/DET
      WI1 =  (V1*W2 - V2*W1)/DET
      UI2 = -(U2*W3 - U3*W2)/DET
      VI2 =  (U1*W3 - U3*W1)/DET
      WI2 = -(U1*W2 - U2*W1)/DET
      UI3 =  (U2*V3 - U3*V2)/DET
      VI3 = -(U1*V3 - U3*V1)/DET
      WI3 =  (U1*V2 - U2*V1)/DET

!     --- PRINT OUT FOR INVERSE MATRIX
      IF(IPRINV.EQ.1) CALL OUTPT1(12)
      RETURN
 end subroutine invers
 subroutine lmpre

 use nbslatticemod, only : ad,alphad,bd,betad,cd,gammad,iamax,iamin, &
                         & ibmax,ibmin,ick014,icmax,icmin,ie,iseqd, &
                         & iunitb,iunitd,ivmax,ivmin,nunk,&
                         & tol,tolv,vd

 implicit none
 
 integer :: ia0(900),ialph0(900),ib0(900),ibeta0(900),ic0(900),&
          & igamm0(900),iseq0(900)

 integer :: i,k0,idelet,irset1,irset2

 real :: a0(900),alpha0(900),b0(900),beta0(900),c0(900),&
       & gamma0(900), v0(900)
 real :: tolv1
 
 interface
   subroutine sort1(a1,a2,a3,a4,a5,a6,a7,n)
   implicit none
   integer,intent(in) :: n
   integer :: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:)
   end subroutine sort1

 end interface

!      dimension A0(900), B0(900), C0(900), ALPHA0(900), BETA0(900)
!      DIMENSION GAMMA0(900), V0(900)
!      DIMENSION IA0(900), IB0(900), IC0(900), IALPH0(900), IBETA0(900)
!      DIMENSION IGAMM0(900), ISEQ0(900)

!      COMMON /LMCELL/ AD(900),BD(900),CD(900),ALPHAD(900),BETAD(900),
!     $                GAMMAD(900),VD(900),ISEQD(900)
!      COMMON /MAX1/ IAMAX(900),IBMAX(900),ICMAX(900),IVMAX(900)
!      COMMON /MIN1/ IAMIN(900),IBMIN(900),ICMIN(900),IVMIN(900)
!      COMMON /PROB2/ IPCH,NUNK
!      COMMON /PROB3/ IE
!      COMMON /TOLER1/ TOL,TOLV
!      COMMON /UNIT2/ IUNITB
!      COMMON /UNIT4/ IUNITD
!**
!     COMMON /CK01/ ICK011,ICK012,ICK013,ICK014


!     A0,B0,C0,ALPHA0,BETA0,GAMMA0,V0 = CELL PARAMETERS AND VOLUME FOR
!          UNKNOWN ENTRIES AS READ FROM INPUT UNIT IUNITD
!     IA0,IB0,IC0,IALPH0,IBETA0,IGAMM0 = ROUNDED OFF CELL PARAMETERS
!          FOR UNKNOWN ENTRIES (ARRAYS ARE SORTED AFTER ROUNDING
!          OFF AND ADDING SEQUENCE NUMBERS)
!     NUNK = NUMBER OF UNKNOWN ENTRIES READ
!     TOL = TOLERANCES FOR A, B, C (ANGSTROMS)
!     TOLV = TOLERANCE FOR VOLUME (PERCENTAGE OF ANGSTROMS**3)
!     ISEQ0,ISEQD = SEQUENCE NUMBER OF UNKNOWN ENTRY IN INPUT FILE (IUNI
!     IDELET = NUMBER OF DUPLICATE ENTRIES ELIMINATED
!     IE = NUMBER OF ENTRIES TO BE CHECKED AGAINST THE 1-LINE DATA BASE
!          (DUPLICATES HAVE BEEN ELIMINATED)
!     AD,BD,CD,ALPHAD,BETAD,GAMMAD,VD = CELL PARAMETERS AND VOLUMES FOR
!           UNKNOWN ENTRIES (DUPLICATES HAVE BEEN ELIMINATED)
!     IAMIN,IBMIN,ICMIN,IVMIN = MINIMUM VALUE OF A,B,C,V FOR A MATCH
!     IAMAX,IBMAX,ICMAX,IVMAX = MAXIMUM VALUE OF A,B,C,V FOR A MATCH


!     --- INITIALIZE VARIABLES
      IE = 1
      K0 = 0
      IDELET = 0
      IRSET1 = 0
      IRSET2 = 0

!     --- IF NECESSARY, RESET NUMBER OF UNKNOWNS, TOLERANCES TO
!         MORE REASONABLE VALUES (ERROR/WARNING MESSAGES ARE
!         WRITTEN WITH SUMMARY DATA)
      IF(NUNK.LT.1.OR.NUNK.GT.900) IRSET1 = 1
      IF(NUNK.GT.900) NUNK = 900
      IF(TOL.LT.0.02.OR.TOLV.LT.5.0) IRSET2 = 1
      IF(TOL.LT.0.02) TOL  = 0.02
      IF(TOLV.LT.5.0) TOLV = 5.0

!     --- READ INPUT UNKNOWNS, ROUND OFF NUMBERS AND ADD SEQUENCE NUMBER
!         (READS IUNITD)
      OPEN(UNIT=IUNITD,FILE='NBS10',STATUS='OLD')
      DO 100 I = 1,NUNK
         READ(IUNITD,1100) A0(I), B0(I), C0(I), ALPHA0(I), BETA0(I),&
     &                GAMMA0(I), V0(I)
         IA0(I)    = NINT(A0(I)*1000)
         IB0(I)    = NINT(B0(I)*1000)
         IC0(I)    = NINT(C0(I)*1000)
         IALPH0(I) = NINT(ALPHA0(I)*100)
         IBETA0(I) = NINT(BETA0(I)*100)
         IGAMM0(I) = NINT(GAMMA0(I)*100)
         ISEQ0(I)  = I
  100 CONTINUE
!**
      IF(ICK014.NE.1) GO TO 160

!        --- FOR CHECKING, WRITE INPUT UNKNOWNS,
!            FOLLOWED BY WRITE OF INPUT UNKNOWNS AFTER ROUNDING
!            OFF AND SEQUENCE NUMBERS ADDED (BEFORE SORTING)
         WRITE(IUNITB,910)
         DO 120 I = 1,NUNK
            WRITE(IUNITB,920) A0(I), B0(I), C0(I), ALPHA0(I), BETA0(I),&
     &                        GAMMA0(I), V0(I)
  120    CONTINUE
         WRITE(IUNITB,910)
         DO 140 I = 1,NUNK
            WRITE(IUNITB,940) IA0(I), IB0(I), IC0(I), IALPH0(I),&
     &                        IBETA0(I), IGAMM0(I), ISEQ0(I)
  140    CONTINUE
  160 CONTINUE

!     --- SORT ROUNDED OFF CELLS ON A,B,C,ALPHA,BETA,GAMMA,SEQUENCE NUMB
      CALL SORT1(IA0,IB0,IC0,IALPH0,IBETA0,IGAMM0,ISEQ0,NUNK)
!      CALL SORT1( IA0(1:nunk),IB0(1:nunk),IC0(1:nunk),IALPH0(1:nunk),&
!                & IBETA0(1:nunk),IGAMM0(1:nunk),ISEQ0(1:nunk),NUNK)
!**
      IF(ICK014.NE.1) GO TO 200

!        --- FOR CHECKING, WRITE INPUT UNKNOWNS AFTER ROUNDING
!            OFF, ADDING SEQUENCE NUMBERS, SORTING
         WRITE(IUNITB,910)
         DO 180 I = 1,NUNK
            WRITE(IUNITB,940) IA0(I), IB0(I), IC0(I), IALPH0(I),&
     &                        IBETA0(I), IGAMM0(I), ISEQ0(I)
  180    CONTINUE
  200 CONTINUE

!     --- WRITE FIRST ENTRY AND ASSIGN VALUES TO ARRAY TO BE
!         CHECKED AGAINST THE 1-LINE DATA BASE (FIRST ENTRY IS
!         NEVER CONSIDERED A DUPLICATE ENTRY)
      K0 = ISEQ0(1)
      WRITE(IUNITB,1300)
      WRITE(IUNITB,1400)
      WRITE(IUNITB,1500) IE, ISEQ0(1), A0(K0), B0(K0), C0(K0),&
     &                   ALPHA0(K0), BETA0(K0), GAMMA0(K0), V0(K0)
      AD(IE)     = A0(K0)
      BD(IE)     = B0(K0)
      CD(IE)     = C0(K0)
      ALPHAD(IE) = ALPHA0(K0)
      BETAD(IE)  = BETA0(K0)
      GAMMAD(IE) = GAMMA0(K0)
      VD(IE)     = V0(K0)
      ISEQD(IE)  = ISEQ0(1)

!     --- BUILD INPUT ARRAY OF CELLS WHICH IS TO BE CHECKED AGAINST
!         THE 1-LINE DATA BASE FOR A MATCH
      DO 500 I = 2,NUNK

!        --- CHECK WHETHER A,B,C,ALPHA,BETA,GAMMA (ROUNDED OFF) ARE EQUA
!            FOR TWO ADJACENT ENTRIES IN SORTED ARRAY OF UNKNOWNS
         IF(IC0(I).NE.IC0(I-1)) GO TO 300
         IF(IB0(I).NE.IB0(I-1).OR.IA0(I).NE.IA0(I-1)) GO TO 300
         IF(IALPH0(I).NE.IALPH0(I-1).OR.IBETA0(I).NE.IBETA0(I-1).OR.&
     &      IGAMM0(I).NE.IGAMM0(I-1)) GO TO 300

!           --- HAVE A DUPLICATE ENTRY ... DO NOT INCLUDE IN ARRAY
!               TO BE CHECKED AGAINST THE 1-LINE DATA BASE
            K0 = ISEQ0(I)
            WRITE(IUNITB,1600) ISEQ0(I), A0(K0), B0(K0), C0(K0),&
     &                         ALPHA0(K0), BETA0(K0), GAMMA0(K0), V0(K0)
            IDELET = IDELET + 1
            GO TO 400
  300    CONTINUE

!        --- NOT A DUPLICATE ENTRY, ASSIGN VALUES TO ARRAYS TO BE
!            CHECKED AGAINST THE 1-LINE DATA BASE
         IE = IE + 1
         K0 = ISEQ0(I)
         WRITE(IUNITB,1500) IE, ISEQ0(I), A0(K0), B0(K0), C0(K0),&
     &                      ALPHA0(K0), BETA0(K0), GAMMA0(K0), V0(K0)
         AD(IE)     = A0(K0)
         BD(IE)     = B0(K0)
         CD(IE)     = C0(K0)
         ALPHAD(IE) = ALPHA0(K0)
         BETAD(IE)  = BETA0(K0)
         GAMMAD(IE) = GAMMA0(K0)
         VD(IE)     = V0(K0)
         ISEQD(IE)  = ISEQ0(I)
  400 CONTINUE
  500 CONTINUE

!     --- WRITE SUMMARY DATA FOR PREPARATION OF INPUT FILE OF UNKNOWNS
!         (AND, IF INDICATED, ANY ERROR/WARNING MESSAGES FOR THE
!         INPUT TOLERANCES, NUMBER OF CELLS TO BE CHECKED
!         AGAINST THE 1-LINE DATA BASE)
      WRITE(IUNITB,1700)
      IF(IRSET1.EQ.1) WRITE(IUNITB,1720)
      IF(IRSET2.EQ.1) WRITE(IUNITB,1740)
      WRITE(IUNITB,1800) NUNK, IDELET, IE, TOL, TOLV
!**
      IF(ICK014.NE.1) GO TO 540

!        --- FOR CHECKING, WRITE UNKNOWNS TO BE CHECKED AGAINST
!            THE 1-LINE DATA BASE FOR A MATCH (AFTER DUPLICATES
!            HAVE BEEN ELIMINATED)
         WRITE(IUNITB,910)
         DO 520 I = 1,IE
            WRITE(IUNITB,920) AD(I), BD(I), CD(I), ALPHAD(I), BETAD(I),&
     &                        GAMMAD(I), VD(I), ISEQD(I)
  520    CONTINUE
  540 CONTINUE

!     --- CALCULATE MINIMUM AND MAXIMUM VALUES FOR A,B,C,V
      DO 600 I = 1,IE
         TOLV1 = 0.01*TOLV*VD(I)
         IAMIN(I) = NINT((AD(I) - TOL)*100)
         IBMIN(I) = NINT((BD(I) - TOL)*100)
         ICMIN(I) = NINT((CD(I) - TOL)*100)
         IVMIN(I) = NINT(VD(I) - TOLV1)
         IAMAX(I) = NINT((AD(I) + TOL)*100)
         IBMAX(I) = NINT((BD(I) + TOL)*100)
         ICMAX(I) = NINT((CD(I) + TOL)*100)
         IVMAX(I) = NINT(VD(I) + TOLV1)
  600 CONTINUE

!     --- DELETE FILE
      CLOSE(UNIT=IUNITD,STATUS='DELETE')
  800 CONTINUE
      RETURN
  910 FORMAT(1H1)
  920 FORMAT(1X,6F10.5,F10.2,I10)
  940 FORMAT(1X,7I10)
 1100 FORMAT(6F10.5,F10.2,I10)
 1300 FORMAT(1H1,35X,'** LATTICE MATCHING **'//1X,20X,'Identification of&
     & unknown lattices using data from'/1X,20X,'the National Bureau of&
     &Standards Crystal Data File'/////)
 1400 FORMAT(1X,37X,'***  UNKNOWNS  ***'&
     &/1X,91(1H-)//,1X,'  Search',T11,' Original',T27,'a',T37,'b'&
     &,T47,'c',T55,'Alpha',T66,'Beta',T75,'Gamma',T87,'Vol'/1X,' Sequenc&
     &e',T11,' Sequence'/1X,91(1H-)/)
 1500 FORMAT(1X,1X,I5,3X,I5,6X,6F10.5,F10.2)
 1600 FORMAT(1X,' Delete ',5X,I5,2X,6F10.5,F10.2)
 1700 FORMAT(/1X,91(1H-)/)
 1720 FORMAT(1H1,'*LMPRE* ERROR ... Check number of unknowns generated (&
     &input to unit 10).'//)
 1740 FORMAT(////1X,'*LMPRE* WARNING ... Input tolerances have been rese&
     &t to minimum allowable values.  It is still'/1X,20X,&
     &'likely that the reset tolerances are unacceptable since the known&
     & data in'/1X,20X,&
     &'the Data Base may not be as accurate and/or as precise as the exp&
     &erimental'/1X,20X,&
     &'data input to the program.  Consult the manual for recommended va&
     &lues.'//)
 1800 FORMAT(////1X,'Original number of  unknowns  =',T35,I5/1X,&
     &'Number of duplicate unknowns  =',T35,I5/1X,&
     &'Number of unknowns for search =',T35,I5//1X,&
     &'Tolerance for cell edges  =',1X,F10.2/1X,&
     &'Tolerance for cell volume =',1X,F10.2)

 end subroutine lmpre
      SUBROUTINE LMSETO

 use nbsLatticeMod, only : idata1,ipoint,istart,nogr,point

!      INTEGER POINT,IPOINT
!      COMMON /DATA1/ IDATA1,IDATA2,IDATA3
!      COMMON /POINT1/ POINT(62),IPOINT(62)
!     COMMON /SRCH1/ ISTART(900),NOGR(900)

      IDATA1 = 62

      IPOINT(1)  =     1
      IPOINT(2)  =   497
      IPOINT(3)  =  1003
      IPOINT(4)  =  1499
      IPOINT(5)  =  1996
      IPOINT(6)  =  2505
      IPOINT(7)  =  3006
      IPOINT(8)  =  3499
      IPOINT(9)  =  4006
      IPOINT(10) =  4484
      IPOINT(11) =  5008
      IPOINT(12) =  5502
      IPOINT(13) =  6001
      IPOINT(14) =  6506
      IPOINT(15) =  6977
      IPOINT(16) =  7482
      IPOINT(17) =  7995
      IPOINT(18) =  8504
      IPOINT(19) =  8999
      IPOINT(20) =  9523
      IPOINT(21) =  9989
      IPOINT(22) = 10515
      IPOINT(23) = 10987
      IPOINT(24) = 11505
      IPOINT(25) = 12001
      IPOINT(26) = 12478
      IPOINT(27) = 13011
      IPOINT(28) = 13476
      IPOINT(29) = 14000
      IPOINT(30) = 14517
      IPOINT(31) = 15017
      IPOINT(32) = 15510
      IPOINT(33) = 15991
      IPOINT(34) = 16513
      IPOINT(35) = 16986
      IPOINT(36) = 17496
      IPOINT(37) = 18006
      IPOINT(38) = 18516
      IPOINT(39) = 18980
      IPOINT(40) = 19523
      IPOINT(41) = 19990
      IPOINT(42) = 20505
      IPOINT(43) = 20997
      IPOINT(44) = 21499
      IPOINT(45) = 21996
      IPOINT(46) = 22506
      IPOINT(47) = 23015
      IPOINT(48) = 23518
      IPOINT(49) = 23995
      IPOINT(50) = 24495
      IPOINT(51) = 24999
      IPOINT(52) = 25503
      IPOINT(53) = 26008
      IPOINT(54) = 26515
      IPOINT(55) = 27011
      IPOINT(56) = 27501
      IPOINT(57) = 28002
      IPOINT(58) = 28511
      IPOINT(59) = 29002
      IPOINT(60) = 29501
      IPOINT(61) = 30000
      IPOINT(62) = 30232

      POINT(1)  =  160
      POINT(2)  =  388
      POINT(3)  =  419
      POINT(4)  =  474
      POINT(5)  =  497
      POINT(6)  =  516
      POINT(7)  =  536
      POINT(8)  =  556
      POINT(9)  =  573
      POINT(10) =  586
      POINT(11) =  599
      POINT(12) =  611
      POINT(13) =  623
      POINT(14) =  635
      POINT(15) =  646
      POINT(16) =  658
      POINT(17) =  670
      POINT(18) =  681
      POINT(19) =  691
      POINT(20) =  702
      POINT(21) =  711
      POINT(22) =  720
      POINT(23) =  728
      POINT(24) =  737
      POINT(25) =  745
      POINT(26) =  753
      POINT(27) =  763
      POINT(28) =  771
      POINT(29) =  780
      POINT(30) =  789
      POINT(31) =  799
      POINT(32) =  809
      POINT(33) =  819
      POINT(34) =  830
      POINT(35) =  840
      POINT(36) =  851
      POINT(37) =  863
      POINT(38) =  876
      POINT(39) =  887
      POINT(40) =  899
      POINT(41) =  911
      POINT(42) =  924
      POINT(43) =  936
      POINT(44) =  949
      POINT(45) =  962
      POINT(46) =  974
      POINT(47) =  988
      POINT(48) = 1003
      POINT(49) = 1018
      POINT(50) = 1034
      POINT(51) = 1054
      POINT(52) = 1076
      POINT(53) = 1101
      POINT(54) = 1129
      POINT(55) = 1157
      POINT(56) = 1192
      POINT(57) = 1234
      POINT(58) = 1278
      POINT(59) = 1349
      POINT(60) = 1452
      POINT(61) = 1679
      POINT(62) = 3595

      DO 100 I = 1,900
         NOGR(I) = 0
         ISTART(I) = 0
  100 CONTINUE

      RETURN
      END subroutine lmseto
 subroutine lmsrch

 use nbsLatticeMod, only : ad,alphad,bd,betad,cd,gammad,iamax,iamin,&
                         & ibmax,ibmin,ick011,ick012,ick013,ick014,&
                         & icmax,icmin,idata1,&
                         & idata2,idata3,ie,ipoint,iseqd,istart,&
                         & iunitb,iunitc,ivmax,ivmin,nogr,point,vd

 implicit none

 character(132)::   reslts
 integer :: krecf(2000), kunkno(2000), kseqno(2000)
 integer :: i,ia1l,ib1l,ic1l,idscan,ig,imatch,index1,index2,irecf,&
          & irecno,iseqno,iset1,iset2,iskip,iv1l,iunkno,ix,j,k,k1
 real :: a1l,b1l,c1l

 interface
  subroutine sort2(a1,a2,a3,n)
  implicit none
  integer :: a1(:), a2(:), a3(:), n
  end subroutine sort2
 end interface

!      INTEGER POINT,IPOINT
!      COMMON /DATA1/ IDATA1,IDATA2,IDATA3
!      COMMON /LMCELL/ AD(900),BD(900),CD(900),ALPHAD(900),BETAD(900),
!     $                GAMMAD(900),VD(900),ISEQD(900)
!      COMMON /MAX1/ IAMAX(900),IBMAX(900),ICMAX(900),IVMAX(900)
!      COMMON /MIN1/ IAMIN(900),IBMIN(900),ICMIN(900),IVMIN(900)
!      COMMON /POINT1/ POINT(62),IPOINT(62)
!      COMMON /PROB3/ IE
!      COMMON /SRCH1/ ISTART(900),NOGR(900)
!      COMMON /UNIT2/ IUNITB
!      COMMON /UNIT3/ IUNITC
!**
!      COMMON /CK01/ ICK011,ICK012,ICK013,ICK014


!     ISEQD = SEQUENCE NUMBER OF UNKNOWN ENTRY IN INPUT FILE (IUNITD)
!     IE = NUMBER OF ENTRIES TO BE CHECKED AGAINST THE 1-LINE
!          DATA BASE FOR A MATCH
!          (AFTER SEQUENCED, SORTED, DUPLICATE CELLS ELIMINATED)
!     AD,BD,CD,ALPHAD,BETAD,GAMMAD,VD = CELL PARAMETERS AND VOLUMES FOR
!           UNKNOWN ENTRIES (AFTER DUPLICATES HAVE BEEN ELIMINATED)
!     A1L,B1L,C1L,IV1L = CELL EDGES AND VOLUME READ FROM 1-LINE DATA BAS
!     IAMIN,IBMIN,ICMIN,IVMIN = MINIMUM VALUE OF A,B,C,V FOR A MATCH
!     IAMAX,IBMAX,ICMAX,IVMAX = MAXIMUM VALUE OF A,B,C,V FOR A MATCH
!     POINT(I) = VALUE OF A AT RECORD NUMBER IPOINT(I)
!     IPOINT(I) = SEQUENCE NUMBER IN 1-LINE FILE
!     IRECNO = FOR EACH GROUP OF UNKNOWN ENTRIES, RECORD NUMBER TO START
!              SEARCHING FOR A MATCH
!     IMATCH = TOTAL NUMBER OF MATCHES FOUND FOR ENTIRE RUN
!     ISTART = FOR EACH UNKNOWN ENTRY, ENTRY NUMBER IN 1-LINE DATA BASE
!              TO START SEARCHING FOR A MATCH
!     ISKIP = FOR EACH GROUP OF ENTRIES, NUMBER OF TIMES GROUPS OF X ENT
!             ARE ELIMINATED FROM FINAL SCAN  (X=IDATA2)
!     IDSCAN = FOR EACH GROUP OF ENTRIES, NUMBER OF ENTRIES CHECKED WHEN
!              LOOKING FOR A MATCH IN THE 1-LINE FILE
!     IG = INPUT ENTRIES ARE DIVIDED INTO IG NUMBER OF GROUPS
!     NOGR = FOR EACH GROUP, NUMBER OF ENTRIES IN THE GROUP
!     INDEX1 = STARTING UNKNOWN NUMBER TO DEFINE A GROUP
!     INDEX2 = ENDING UNKNOWN NUMBER TO DEFINE A GROUP


!     --- INITIALIZE VARIABLES
      IMATCH = 0

!     --- IF INPUT DATA IS NOT WITHIN RANGE OF A VALUES IN THE DATA BASE
!         NO MATCHES ARE POSSIBLE ... DO NOT READ FILE
      IF(IAMIN(1).GT.POINT(IDATA1).OR.IAMAX(IE).LT.POINT(1))&
     &   WRITE(IUNITB,1000)
      IF(IAMIN(1).GT.POINT(IDATA1).OR.IAMAX(IE).LT.POINT(1)) GO TO 700

!     --- CALCULATE APPROXIMATE STARTING VALUES FOR EACH UNKNOWN
      DO 100 J = 1,IE
         K1 = 0
         DO 40 I = 1,IDATA1
            IF(IAMIN(J).LT.POINT(I)) K1 = I
            IF(IAMIN(J).LT.POINT(I)) GO TO 80
   40    CONTINUE
   80    CONTINUE
         IF(K1.EQ.1) ISTART(J) = 1
         IF(IAMIN(J).GE.POINT(IDATA1)) ISTART(J) = IPOINT(IDATA1-1)
         IF(K1.GT.1) ISTART(J) = IPOINT(K1-1)
  100 CONTINUE
!**
!     --- FOR CHECKING, WRITE INTERMEDIATE VARIABLES
      IF(ICK012.EQ.1) WRITE(IUNITB,740)
      IF(ICK012.EQ.1) WRITE(IUNITB,741) (I,IAMIN(I),IAMAX(I),ISTART(I),&
     &   IVMIN(I),IVMAX(I),I=1,IE)

!     --- DEFINE GROUPS OF UNKNOWNS AND COUNT NUMBER OF GROUPS FORMED
      IG = 1
      NOGR(1) = 1
      DO 200 J = 1,IE-1
         IF(IAMAX(J).LT.IAMIN(J+1)) GO TO 180
            NOGR(IG) = NOGR(IG) + 1
            GO TO 200
  180    CONTINUE
         IG = IG + 1
         NOGR(IG) = 1
  200 CONTINUE
!**
!     --- FOR CHECKING, WRITE INTERMEDIATE VARIABLES
      IF(ICK011.EQ.1) WRITE(IUNITB,745) IG
      IF(ICK011.EQ.1) WRITE(IUNITB,746) (I,NOGR(I),I=1,IG)

!     --- START PROCEDURE TO FIND MATCHES WITH DATA BASE ... EACH GROUP
!         OF UNKNOWN ENTRIES IS TREATED SEPARATELY
      INDEX2 = 0
      DO 500 J = 1,IG
         ISKIP  = 0
         IDSCAN = 0
         IX = INDEX2
         INDEX1 = IX + 1
         INDEX2 = IX + NOGR(J)
         IRECNO = ISTART(INDEX1)
!**
!        --- FOR CHECKING, WRITE INTERMEDIATE VARIABLES
         IF(ICK011.EQ.1) WRITE(IUNITB,755) J, INDEX1, INDEX2

!        --- PRE-SCAN OF 1-LINE DATA BASE TO FIND STARTING RECORD NUMBER
!            CLOSE TO MINIMUM VALUE OF A FOR EACH GROUP ... READS EVERY
!            X-TH VALUE OF A IN 1-LINE DATA BASE (X=IDATA2)
         DO 300 I = 1,IDATA3
            READ(IUNITC,1100,REC=IRECNO) A1L
!**
!           --- FOR CHECKING, WRITE INTERMEDIATE VARIABLES
            IF(ICK013.EQ.1) WRITE(IUNITB,760) I, ISKIP, IRECNO, A1L
            IA1L = NINT(A1L*100)
            IF(IAMIN(INDEX1).LE.IA1L) GO TO 340
               IRECNO = IRECNO + IDATA2
               IF(IRECNO.GE.IPOINT(IDATA1)) IRECNO = IPOINT(IDATA1)
               IF(IRECNO.GE.IPOINT(IDATA1)) GO TO 340
                  ISKIP = ISKIP + 1
  300    CONTINUE
  340    CONTINUE
         IRECNO = IRECNO - IDATA2
         IF(IRECNO.LT.1) IRECNO = 1
!**
!        --- FOR CHECKING, WRITE INTERMEDIATE VARIABLES
         IF(ICK012.EQ.1) WRITE(IUNITB,765) ISKIP, IRECNO

!        --- SEARCH SELECTED REGION OF 1-LINE DATA BASE FOR A MATCH
!            OF A,B,C,V WITHIN SPECIFIED TOLERANCES
  380    CONTINUE
         READ(IUNITC,1100,REC=IRECNO) A1L, B1L, C1L, IV1L
         IA1L = NINT(A1L*100)
         IB1L = NINT(B1L*100)
         IC1L = NINT(C1L*100)

!        --- CHECK THAT VALUE OF A FROM 1-LINE DATA BASE IS LESS THAN
!            MAXIMUM VALUE OF A IN GROUP ... IF NOT, STOP SEARCHING FOR
!            A MATCH
         IF(IA1L.GT.IAMAX(INDEX2)) GO TO 440
            IDSCAN = IDSCAN + 1

!           --- CHECK FOR A MATCH OF A,B,C,V FOR EACH UNKNOWN ENTRY IN G
            DO 400 K = INDEX1,INDEX2
               IF(IC1L.LT.ICMIN(K).OR.IC1L.GT.ICMAX(K).OR.&
     &            IB1L.LT.IBMIN(K).OR.IB1L.GT.IBMAX(K).OR.&
     &            IV1L.LT.IVMIN(K).OR.IV1L.GT.IVMAX(K).OR.&
     &            IA1L.LT.IAMIN(K).OR.IA1L.GT.IAMAX(K)) GO TO 400

!                 --- HAVE A MATCH, ASSIGN ARRAYS SPECIFYING
!                     RECORD NUMBER (OF MATCH) IN THE 1-LINE
!                     DATA BASE, UNKNOWN NUMBER (SEARCH SEQUENCE),
!                     AND ORIGINAL SEQUENCE NUMBER OF UNKNOWN
                  IMATCH = IMATCH + 1
                  IF(IMATCH.GT.2000) GO TO 510
                     KRECF(IMATCH)  = IRECNO
                     KUNKNO(IMATCH) = K
                     KSEQNO(IMATCH) = ISEQD(K)
  400       CONTINUE
            IRECNO = IRECNO + 1
            IF(IRECNO.GT.IPOINT(IDATA1)) GO TO 440
               GO TO 380
  440    CONTINUE
!**
!        --- FOR CHECKING, WRITE INTERMEDIATE VARIABLES
         IF(ICK011.EQ.1) WRITE(IUNITB,770) IDSCAN, IMATCH
  500 CONTINUE
  510 CONTINUE

!     --- CHECK NUMBER OF MATCHES FOUND
      IF(IMATCH.GT.0.AND.IMATCH.LE.2000) GO TO 530
         IF(IMATCH.LE.2000) GO TO 520

!           --- TOO MANY MATCHES GENERATED,
!               WRITE ERROR MESSAGE AND STOP PROGRAM EXECUTION
            WRITE(IUNITB,1200)
            STOP
  520    CONTINUE

!        --- NO MATCHES FOUND, WRITE MESSAGE AND CONTINUE
         WRITE(IUNITB,1300)
         GO TO 700
  530 CONTINUE
!**
      IF(ICK014.NE.1) GO TO 550

!        --- FOR CHECKING, WRITE ARRAYS WITH INDICES TO
!            MATCHES AND UNKNOWNS (BEFORE SORTING)
         WRITE(IUNITB,720)
         DO 540 I = 1,IMATCH
            WRITE(IUNITB,780) KRECF(I), KUNKNO(I), KSEQNO(I)
  540    CONTINUE
  550 CONTINUE

!     --- SORT ARRAYS WITH INDICES TO MATCHES AND UNKNOWNS
!         (FIRST SORT ON ORIGINAL SEQUENCE NUMBER, SECOND SORT
!         ON RECORD NUMBER OF MATCH FROM 1-LINE DATA BASE, THIRD
!         SORT ON SEARCH SEQUENCE (UNKNOWN NUMBER FOR SEARCH))

!     --- THIRD SORT IS NOT REQUIRED
      CALL SORT2(KSEQNO,KRECF,KUNKNO,IMATCH)
!**
      IF(ICK014.NE.1) GO TO 570

!        --- FOR CHECKING, WRITE ARRAYS WITH INDICES TO MATCHES
!            AND UNKNOWNS (AFTER SORTING)
         WRITE(IUNITB,720)
         DO 560 I = 1,IMATCH
             WRITE(IUNITB,780) KRECF(I), KUNKNO(I), KSEQNO(I)
  560    CONTINUE
  570 CONTINUE

!     --- WRITE RESULTS AND SUMMARY DATA
      WRITE(IUNITB,1400)
      ISET1 = 0
      ISET2 = 1
      DO 600 I = 1,IMATCH
         IRECF  = KRECF(I)
         IUNKNO = KUNKNO(I)
         ISEQNO = KSEQNO(I)

         READ(IUNITC,1500,REC=IRECF) RESLTS
         IF(ISEQNO.GT.ISET2) ISET1 = 0
         IF(ISET1.EQ.0) WRITE(IUNITB,1600) IUNKNO, ISEQD(IUNKNO),&
     &      AD(IUNKNO), BD(IUNKNO), CD(IUNKNO), ALPHAD(IUNKNO),&
     &      BETAD(IUNKNO), GAMMAD(IUNKNO), VD(IUNKNO)
         WRITE(IUNITB,1700) RESLTS
         ISET1 = 1
         ISET2 = ISEQNO
  600 CONTINUE
      WRITE(IUNITB,1800)
  700 CONTINUE
      RETURN
  720 FORMAT(1H1)
  740 FORMAT(/1X,'UNK NUMBER',10X,'IAMIN',10X,'IAMAX',9X,'ISTART' ,11X,&
     &'IVMIN',10X,'IVMAX'/)
  741 FORMAT(1X,I10,5X,I10,5X,I10,5X,I10,5X,I10,5X,I10)
  745 FORMAT(/1X,'745   NUMBER OF GROUPS = ',I10//1X,'GROUP NUMBER',5X,'&
     &NUMBER IN GROUP')
  746 FORMAT(1X,I10,9X,I10)
  755 FORMAT(/1X,'755   Group number = ',I10,10X,'INDEX1 = ',I10,10X,'IN&
     &DEX2 = ',I10)
  760 FORMAT(/1X,'760   I = ',I10,10X,'ISKIP = ',I10,10X,'IRECNO = ',I10&
     &,10X,'A1L = ',F6.2)
  765 FORMAT(/1X,'765   ISKIP = ',I10,10X,'IRECNO = ',I10)
  770 FORMAT(/1X,'Number of entries searched = ',I7,10X,'Number of match&
     &es = ',I7)
  780 FORMAT(1X,3I10)
 1000 FORMAT(////1X,'** Input Data Not Within Range of Data Base ... No&
     &Matches Possible **')
 1100 FORMAT(3F6.2,18X,I6)
 1200 FORMAT(/////1X,'*LMSRCH* ERROR ... More than 2000 matches have bee&
     &n found.'/1X,19X,'Input parameters and/or program limits must be c&
     &hanged.')
 1300 FORMAT(////1X,'** NONE of Unknown Entries Matched Those in the Dat&
     &a Base **')
 1400 FORMAT(1X,37X,'***  RESULTS  ***'/1X,&
     &91(1H-)//,4X,'  Search',T15,' Original',T28,'a',T38,'b'&
     &,T48,'c',T56,'Alpha',T67,'Beta',T76,'Gamma',T88,'Vol'/4X,' Sequenc&
     &e',T15,' Sequence'/1X,91(1H-)/)
 1500 FORMAT(A132)
 1600 FORMAT(/1X,5X,I5,I10,3F10.3,3F10.2,F10.1)
 1700 FORMAT(1X,A132)
 1800 FORMAT(/1X,91(1H-)/)
 end subroutine lmsrch
      SUBROUTINE MNCOND

 use nbsLatticeMod, only : ick021,ick022,ick023,ick024,ick025,ick026,&
                         & ick027,ick028,ick029,itype,iunitb, &
                         & s11,s22,s33,s23,s13,s12,u

 implicit none

 integer :: ict1,ict2,ict3
 real :: a12,a13,a23,q1,q2,q3,tx1,tx2
 real :: pc1=0.00006

 interface
   subroutine ckpt02(iflag2)
    implicit none
    integer :: iflag2
   end subroutine ckpt02

   subroutine multip
   end subroutine multip

   subroutine normal
   end subroutine normal

   subroutine set
    end subroutine set

   subroutine shortv(sxx,sxy,elem1,elem5,elem9,elemx)
   implicit none
   real :: elem1,elem5,elem9,elemx,sxx,sxy
   end subroutine shortv

   subroutine trans(iprtra)
    implicit none
    integer :: iprtra
   end subroutine trans

 end interface

!     --- INITIALIZE VARIABLES
      ICT1 = 0
      ICT2 = 0
      ICT3 = 0


!      --- START CALCULATIONS TO SATISFY MAIN CONDITIONS
   50 CONTINUE
!**
!     --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE VARIABLES
      IF(ICK021.EQ.1) CALL CKPT02(1)

      A23 = ABS(S23)
      A13 = ABS(S13)

      A12 = ABS(S12)

!     --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY WHEN
!         THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET FOR THE
!         SPECIFIC COMPUTER
      ICT1 = ICT1 + 1
      IF(ICT1.GT.21) WRITE(IUNITB,1000)
      IF(ICT1.GT.21) GO TO 300

!     --- CONTINUE WHEN HAVE THE SHORTEST VECTORS IN THE BC PLANE;
!         OTHERWISE, CALCULATE
      IF(.NOT.((2.0*A23-S22-PC1).GT.0.0.OR.(2.0*A23-S33-PC1).GT.0.0))&
     &   GO TO 100
!**
!        --- FOR CHECKING, WRITE EXECUTION POINT
         IF(ICK022.EQ.1) CALL CKPT02(2)

!        --- CALCULATE SHORTEST VECTORS IN BC PLANE,
!            UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY
!            THE RESULTING MATRIX TO THE INPUT CELL
         CALL SHORTV(S22,S23,U(1),U(5),U(9),U(8))
         CALL MULTIP
         CALL TRANS(0)
         CALL SHORTV(S33,S23,U(1),U(5),U(9),U(6))
         CALL MULTIP
         CALL TRANS(0)
         GO TO 50
  100 CONTINUE

!     --- CONTINUE WHEN HAVE THE SHORTEST VECTORS IN THE AC PLANE;
!         OTHERWISE, CALCULATE
      IF(.NOT.((2.0*A13-S11-PC1).GT.0.0.OR.(2.0*A13-S33-PC1).GT.0.0))&
     &   GO TO 200
!**
!        --- FOR CHECKING, WRITE EXECUTION POINT
         IF(ICK022.EQ.1) CALL CKPT02(3)

!        --- CALCULATE SHORTEST VECTORS IN AC PLANE,
!            UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY
!            THE RESULTING MATRIX TO THE INPUT CELL
         CALL SHORTV(S11,S13,U(1),U(5),U(9),U(7))
         CALL MULTIP
         CALL TRANS(0)
         CALL SHORTV(S33,S13,U(1),U(5),U(9),U(3))
         CALL MULTIP
         CALL TRANS(0)
         GO TO 50
  200 CONTINUE

!     --- CONTINUE WHEN HAVE THE SHORTEST VECTORS IN THE AB PLANE;
!         OTHERWISE, CALCULATE
      IF(.NOT.((2.0*A12-S11-PC1).GT.0.0.OR.(2.0*A12-S22-PC1).GT.0.0))&
     &   GO TO 300
!**
!        --- FOR CHECKING, WRITE EXECUTION POINT
         IF(ICK022.EQ.1) CALL CKPT02(4)

!        --- CALCULATE SHORTEST VECTORS IN AB PLANE,
!            UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY
!            THE RESULTING MATRIX TO THE INPUT CELL
         CALL SHORTV(S11,S12,U(1),U(5),U(9),U(4))
         CALL MULTIP
         CALL TRANS(0)
         CALL SHORTV(S22,S12,U(1),U(5),U(9),U(2))
         CALL MULTIP
         CALL TRANS(0)
         GO TO 50
  300 CONTINUE

  400 CONTINUE

!     --- ARRANGE SO A IS <= TO B AND B IS <= TO C

!     --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY WHEN
!         THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET FOR THE
!         SPECIFIC COMPUTER
      ICT2 = ICT2 + 1
      IF(ICT2.GT.11) WRITE(IUNITB,1100)
      IF(ICT2.GT.11) GO TO 500

!     --- CHECK THAT A IS <= TO B
      IF(S11.LE.(S22+PC1)) GO TO 420

!        --- INTERCHANGE VECTORS A AND B
         CALL SET
         U(2) =  1.0
         U(4) =  1.0
         U(9) = -1.0
         CALL MULTIP
         CALL TRANS(0)
!**
!        --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE
!            VARIABLES
         IF(ICK023.EQ.1) CALL CKPT02(5)

  420 CONTINUE

!     --- CHECK THAT B IS <= TO C
      IF(S22.LE.(S33+PC1)) GO TO 440

!        --- INTERCHANGE VECTORS B AND C
         CALL SET
         U(1) =  1.0
         U(6) =  1.0
         U(8) = -1.0
         CALL MULTIP
         CALL TRANS(0)
!**
!        --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE
!            VARIABLES
         IF(ICK023.EQ.1) CALL CKPT02(6)
  440 CONTINUE

!     --- CHECK THAT A IS <= TO B AND B IS <= TO C
      IF(.NOT.(S11.LE.(S22+PC1).AND.S22.LE.(S33+PC1))) GO TO 400
  500 CONTINUE

!     --- PUT CELL MATRIX IN NORMAL REPRESENTATION
      CALL NORMAL
      IF(ITYPE.EQ.1) GO TO 700

!        --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY
!            WHEN THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET
!            FOR THE SPECIFIC COMPUTER
         ICT3 = ICT3 + 1
         IF(ICT3.GT.5) WRITE(IUNITB,1200)
         IF(ICT3.GT.5) GO TO  700

!           --- TYPE 2 CELL MATRIX ... CHECK PART 2 OF MAIN CONDITIONS
            TX1 = S11 + S22 + PC1
            TX2 = 2.0*(ABS(S23) + ABS(S13) + ABS(S12))
            IF(TX1.GE.TX2) GO TO 700
               CALL SET
               U(1) = 1.0
               U(5) = 1.0
               U(7) = 1.0
               U(8) = 1.0
               U(9) = 1.0
               CALL MULTIP
               CALL TRANS(0)
!**
!              --- FOR CHECKING, WRITE EXECUTION POINT AND
!                  INTERMEDIATE VARIABLES
               IF(ICK024.EQ.1) CALL CKPT02(7)

               ICT2 = 0
               GO TO 400
  700 CONTINUE

!     --- MAIN CONDITIONS SHOULD BE SATISFIED ...
!         RE-CHECK, WRITE WARNING IF NOT SATISFIED AND CONTINUE
      Q1 = 2.0*ABS(S23) - 2.0*PC1
      Q2 = 2.0*ABS(S13) - 2.0*PC1
      Q3 = 2.0*ABS(S12) - 2.0*PC1
      IF(.NOT.(S11.LE.(S22+2.0*PC1).AND.S22.LE.(S33+2.0*PC1).AND.Q1.LE.&
     &   S22.AND.Q2.LE.S11.AND.Q3.LE.S11)) WRITE(IUNITB,1300)
      IF(ITYPE.EQ.1) GO TO 800

!        --- TYPE 2 CELL MATRIX (---), TEST EXTRA MAIN CONDITION
         TX1 = S11 + S22 + 2.0*PC1
         TX2 = 2.0*(ABS(S23) + ABS(S13) + ABS(S12))
         IF(TX1.LT.TX2) WRITE(IUNITB,1400)
  800 CONTINUE
      RETURN
 1000 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a&
     &ll conditions for reduction.'/1X,21X,'Program error when determini&
     &ing the shortest vectors in the BC, AC, AB planes.'&
     &/1X,21X,'The program constant (PC1) may be incorrectly set for thi&
     &s computer.'/)
 1100 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a&
     &ll conditions for reduction.'/1X,21X,'Program error when ordering&
     &the cell edges so that a<=b<=c.'&
     &/1X,21X,'The program constant (PC1) may be incorrectly set for thi&
     &s computer.'/)
 1200 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a&
     &ll conditions for reduction.'/1X,21X,'Program error when checking&
     &Part II of Main Conditions for a type II cell.'&
     &/1X,21X,'The program constant (PC1) may be incorrectly set for thi&
     &s computer.'/)
 1300 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a&
     &ll conditions for reduction.'/1X,21X,'Program error in the final c&
     &heck of Part I of the Main Conditions.'&
     &/1X,21X,'The program constant (PC1) may be incorrectly set for thi&
     &s computer.'/)
 1400 FORMAT(/1X,'*MNCOND* WARNING ... Cell may not be reduced - Check a&
     &ll conditions for reduction.'/1X,21X,'Program error in the final c&
     &heck of Part II of the Main Conditions.'&
     &/1X,21X,'The program constant (PC1) may be incorrectly set for thi&
     &s computer.'/)
      END
      subroutine multip

 use nbsLatticeMod, only : u,u1,u2,u3,v1,v2,v3,t,&
                         & w1,w2,w3
 implicit none

!     --- (PRODUCT MATRIX) = (U-MATRIX) (T-MATRIX)
      U1 = T(1)*U(1) + T(4)*U(2) + T(7)*U(3)
      V1 = T(2)*U(1) + T(5)*U(2) + T(8)*U(3)
      W1 = T(3)*U(1) + T(6)*U(2) + T(9)*U(3)
      U2 = T(1)*U(4) + T(4)*U(5) + T(7)*U(6)
      V2 = T(2)*U(4) + T(5)*U(5) + T(8)*U(6)
      W2 = T(3)*U(4) + T(6)*U(5) + T(9)*U(6)
      U3 = T(1)*U(7) + T(4)*U(8) + T(7)*U(9)
      V3 = T(2)*U(7) + T(5)*U(8) + T(8)*U(9)
      W3 = T(3)*U(7) + T(6)*U(8) + T(9)*U(9)
      RETURN
 end subroutine multip
      SUBROUTINE NORMAL

 use nbsLatticeMod, only : cosa,cosb,cosg,ick021,ick022,ick023,&
                         & ick024,ick025,ick026,ick027, &
                         & ick028,ick029,itype,u,var90

 implicit none

 integer :: iprod,irev12,irev13,irev23,isgn12,isgn13,isgn23
 real :: acosa,acosb,acosg

 interface
  subroutine ckpt02(iflag2)
  implicit none
  integer :: iflag2
  end subroutine ckpt02

  subroutine multip
  end subroutine multip

  subroutine set
  end subroutine set

  subroutine trans(iprtra)
  implicit none
  integer :: iprtra
  end subroutine trans

 end interface


!     THE CELL IS SET AS A TYPE 1 OR TYPE 2 CELL ...
!     TYPE 1 CELL = (+++) ALL ANGLES LESS THAN 90 DEGREES
!                         ( COS(ANGLES) ARE ALL > THAN ZERO )
!     TYPE 2 CELL = (---) ALL ANGLES EQUAL TO OR GREATER THAN 90 DEGREES
!                         ( COS(ANGLES) ARE ALL <= TO ZERO )

!     IF NECESSARY, THIS SUBROUTINE CHANGES THE SIGNS ON THE UNSYMMETRIC
!     SCALARS (S23, S13, S12) SO THAT THEY ARE ALL POSITIVE, ALL
!     NEGATIVE, OR ZERO AND NEGATIVE.  ALSO, THE TOTAL TRANSFORMATION
!     MATRIX IS UPDATED.


!     --- INITIALIZE VARIABLES
      ITYPE = 0

      ACOSA = ABS(COSA)
      ACOSB = ABS(COSB)
      ACOSG = ABS(COSG)

!     --- VARIABLES ISGN-- INDICATE WHETHER THE UNSYMMETRICAL SCALARS
!         (S23, S13, S12) ARE ZERO, POSITIVE OR NEGATIVE (THIS CODE
!         USES THE COS(ANGLES) AS AN INDICATION OF THE SIGNS)
      ISGN23 = -1
      ISGN13 = -1
      ISGN12 = -1
      IF(COSA.GT.VAR90)  ISGN23 = 1
      IF(COSB.GT.VAR90)  ISGN13 = 1
      IF(COSG.GT.VAR90)  ISGN12 = 1
      IF(ACOSA.LE.VAR90) ISGN23 = 0
      IF(ACOSB.LE.VAR90) ISGN13 = 0
      IF(ACOSG.LE.VAR90) ISGN12 = 0

!     --- INITIALIZE MATRIX ARRAY U TO IDENTITY MATRIX
      CALL SET
      U(1) = 1.0
      U(5) = 1.0
      U(9) = 1.0

!     --- TYPE 1 CELL IF PRODUCT OF SIGNS OF UNSYMMETRICAL SCALARS IS +1
!         TYPE 2 CELL IF PRODUCT IS -1 OR 0
      IPROD = ISGN23*ISGN13*ISGN12

!     --- PUT THE CELL MATRIX IN NORMAL REPRESENTATION
      IF(IPROD.EQ.0.OR.IPROD.EQ.-1) GO TO 100

!        --- TYPE 1 CELL ... POSITIVE REDUCED FORM
!            SET UNSYMMETRICAL SCALARS TO BE POSITIVE (ALL ANGLES
!            LESS THAN 90 DEGREES)
         ITYPE = 1

!        --- UPDATE TRANSFORMATION MATRIX
         IF(ISGN23.EQ.-1) U(1) = -1.0
         IF(ISGN13.EQ.-1) U(5) = -1.0
         IF(ISGN12.EQ.-1) U(9) = -1.0
         GO TO 200
  100 CONTINUE

!     --- TYPE 2 CELL ... NEGATIVE REDUCED FORM
!         SET UNSYMMETRICAL SCALARS TO BE LESS THAN OR EQUAL TO ZERO
!         (ALL ANGLES >= 90 DEGREES)
      ITYPE = 2

!     --- UPDATE TRANSFORMATION MATRIX
      IF(ISGN23.EQ.1) U(1) = -1.0
      IF(ISGN13.EQ.1) U(5) = -1.0
      IF(ISGN12.EQ.1) U(9) = -1.0

!     --- HAVE ONE OR MORE 90 DEGREE ANGLES IF IPROD = 0
      IF(IPROD.NE.0) GO TO 200

!        --- CHECK WHETHER OR NOT MUST REVERSE DIRECTION OF ONE
!            ADDITIONAL VECTOR TO MAINTAIN A RIGHT-HANDED SYSTEM
         IREV23 = ISGN23
         IREV13 = ISGN13
         IREV12 = ISGN12
         IF(ISGN23.EQ.0) IREV23 = -1
         IF(ISGN13.EQ.0) IREV13 = -1
         IF(ISGN12.EQ.0) IREV12 = -1
         IPROD = IREV23*IREV13*IREV12
         IF(IPROD.EQ.-1) GO TO 200

!           --- MUST REVERSE DIRECTION OF ONE VECTOR TO MAINTAIN A
!               RIGHT-HANDED SYSTEM
            IF(ISGN23.EQ.0) U(1) = -1.0
            IF(ISGN23.EQ.0) GO TO 200
            IF(ISGN13.EQ.0) U(5) = -1.0
            IF(ISGN13.EQ.0) GO TO 200
            IF(ISGN12.EQ.0) U(9) = -1.0
  200 CONTINUE

!     --- UPDATE THE TOTAL TRANSFORMATION MATRIX AND APPLY THE
!         RESULTING MATRIX TO THE INPUT CELL.  TO UPDATE THE TOTAL
!         TRANSFORMATION MATRIX, THE PROGRAM PREMULTIPLIES THE
!         CURRENT MATRIX BY AN UPDATE MATRIX U.  IN THIS ROUTINE, U
!         MAY BE (-1 0 0 / 0 -1 0 / 0 0 1), (-1 0 0 / 0 1 0 / 0 0 -1)
!         OR (1 0 0 / 0 -1 0 / 0 0 -1).
      CALL MULTIP
      CALL TRANS(0)
!**
!     --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE VARIABLES
      IF(ICK025.EQ.1) CALL CKPT02(9)

!     --- STOP PROGRAM EXECUTION IF PROGRAM FAILS TO TYPE
!         CELL MATRIX (SHOULD NOT OCCUR)
      IF(ITYPE.EQ.0) STOP&
     &  '*NORMAL* ERROR ... Failed to type cell matrix.'
      RETURN
end subroutine normal
      subroutine outpt1(iflag1)
 use nbsLatticeMod, only : af,ai,alphf,alphi,bf,bi,betf,beti,cf,ci,&
                         & gamf,gami,ick021,ick022,ick023,ick024,&
                         & ick025,ick026,ick027,ick028,ick029,&
                         & icntr1,idet,ilatno,iprob,iprrss,irss,&
                         & iunitb,ntd,s11,s12,s13,s22,s23,s33,&
                         & ui1,ui2,ui3,u1,u2,u3,vi1,vi2,vi3,&
                         & v1,v2,v3,volf,voli,wi1,wi2,wi3,w1,&
                         & w2,w3

 implicit none
 integer, intent(in) :: iflag1

 interface
   subroutine ckpt02(iflag2)
   implicit none
   integer :: iflag2
   end subroutine ckpt02
 end interface

!**
!     --- WHEN NOT EXECUTING SPECIAL CHECK RUN, GO TO SECTION
!         OF CODE FOR NORMAL EXECUTION OF PROGRAM
      IF(ICK029.NE.1) GO TO 5

!        --- WRITE PROBLEM NUMBER IF INDICATED
         IF(IFLAG1.EQ.3) CALL CKPT02(13)
         GO TO 900
    5 CONTINUE

!     --- RETURN IF NO REDUCTION AND DERIVATIVE LATTICE OUTPUT
      IF(IPRRSS.EQ.1.AND.IFLAG1.GT.2) GO TO 900

!     --- GO TO APPROPRIATE SECTION OF OUTPUT (INDICATED BY IFLAG1)
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170)&
     & IFLAG1
   10 CONTINUE

!     *** CALLED FROM 'MAIN'

!     --- WRITE NUMBER OF INDEPENDENT PROBLEMS TO STUDY
      WRITE(IUNITB,1100)
      WRITE(IUNITB,1150) NTD
      GO TO 900
   20 CONTINUE

!     *** CALLED FROM 'MAIN'

!     --- WRITE FOR LATTICE MATCHING OF REDUCED CELLS
      WRITE(IUNITB,1200)
      GO TO 900
   30 CONTINUE

!     *** CALLED FROM 'REDUCE'

!     --- WRITE PROBLEM NUMBER
      IF(IPROB.NE.1) WRITE(IUNITB,1300)
      WRITE(IUNITB,1350) IPROB
      GO TO 900
   40 CONTINUE

!     *** CALLED FROM 'READ'

!     --- HEADINGS PRODUCED BY PROBLEM SPECIFICATION AND CELL CARD
      IF(IRSS.EQ.0) WRITE(IUNITB,1400)
      IF(IRSS.EQ.1) WRITE(IUNITB,1500)
      IF(IRSS.EQ.2) WRITE(IUNITB,1600)
      IF(IRSS.EQ.3) WRITE(IUNITB,1700)
      GO TO 900
   50 CONTINUE

!     *** CALLED FROM 'READ'

!     --- HEADING FOR CENTERING OF INPUT CELL
      IF(ICNTR1.EQ.0.OR.ICNTR1.EQ.1) WRITE(IUNITB,1800)
      IF(ICNTR1.EQ.2) WRITE(IUNITB,1900)
      IF(ICNTR1.EQ.3) WRITE(IUNITB,2000)
      IF(ICNTR1.EQ.4) WRITE(IUNITB,2100)
      IF(ICNTR1.EQ.5) WRITE(IUNITB,2200)
      IF(ICNTR1.EQ.6) WRITE(IUNITB,2300)
      IF(ICNTR1.EQ.7) WRITE(IUNITB,2400)
      GO TO 900
   60 CONTINUE
   70 CONTINUE

!     *** CALLED FROM 'DERIV'

!     --- HEADING FOR SUPERLATTICES
      WRITE(IUNITB,2500) IDET
      GO TO 900
   80 CONTINUE

!     *** CALLED FROM 'DERIV'

!     --- HEADING FOR SUBLATTICES
      WRITE(IUNITB,2600) IDET
      GO TO 900
   90 CONTINUE

!     *** CALLED FROM 'DERIV'

!     --- WRITE LATTICE NUMBER IF SUPERCELL OR SUBCELL
      IF(IRSS.EQ.1.OR.IRSS.EQ.3) WRITE(IUNITB,2700) ILATNO
      IF(IRSS.EQ.2) WRITE(IUNITB,2800) ILATNO
      GO TO 900
  100 CONTINUE

!     *** CALLED FROM 'CENTER' OR 'DERIV'

!     --- WRITE TRANSFORMATION MATRIX
      WRITE(IUNITB,2900) U1, V1, W1, U2, V2, W2, U3, V3, W3
      WRITE(11,2900)     U1, V1, W1, U2, V2, W2, U3, V3, W3    ! DU ADD
      GO TO 900
  110 CONTINUE

!     *** CALLED FROM 'TRANS'

!     --- WRITE TRANSFORMATION MATRIX
      WRITE(IUNITB,3000) U1, V1, W1, U2, V2, W2, U3, V3, W3
      WRITE(11,3000)     U1, V1, W1, U2, V2, W2, U3, V3, W3    ! DU ADD
      GO TO 900
  120 CONTINUE

!     *** CALLED FROM 'INVERS'

!     --- WRITE INVERSE MATRIX
      WRITE(IUNITB,3100) UI1, VI1, WI1, UI2, VI2, WI2,&
     &      UI3, VI3, WI3
      WRITE(11,3100)     UI1, VI1, WI1, UI2, VI2, WI2,&         ! DU ADD
     &      UI3, VI3, WI3
      GO TO 900
  130 CONTINUE

!     *** CALLED FROM 'TRANS'

!     --- HEADING FOR CELL MATRIX
      WRITE(IUNITB,3200)
      WRITE(11,3200)                                           ! DU ADD

!     --- WRITE INITIAL CELL, FINAL CELL AND CELL MATRIX
      WRITE(IUNITB,3300) AI, BI, CI, ALPHI, BETI, GAMI, VOLI,&
     &      S11, S22, S33
      WRITE(11,3300)     AI, BI, CI, ALPHI, BETI, GAMI, VOLI,&  ! DU ADD
     &      S11, S22, S33
      WRITE(IUNITB,3400) AF, BF, CF, ALPHF, BETF, GAMF, VOLF,&
     &      S23, S13, S12
      WRITE(11,3400)     AF, BF, CF, ALPHF, BETF, GAMF, VOLF,&  ! DU ADD
     &      S23, S13, S12
      GO TO 900
  140 CONTINUE

!     *** CALLED FROM 'MAIN'

!     --- HEADING FOR CELL TRANSFORMATION
      WRITE(IUNITB,3500)
      GO TO 900
  150 CONTINUE

!     *** CALLED FROM 'MAIN'

!     --- SUB-HEADING FOR CELL TRANSFORMATION
      WRITE(IUNITB,3600) IPROB
      GO TO 900
  160 CONTINUE

!     *** CALLED FROM 'MAIN'

!     --- HEADING FOR MATRIX INVERSION
      WRITE(IUNITB,3700)
      GO TO 900
  170 CONTINUE

!     *** CALLED FROM 'MAIN'

!     --- SUB-HEADING FOR MATRIX INVERSION
      WRITE(IUNITB,3800) IPROB

  900 CONTINUE
      RETURN
 1100 FORMAT(////1X,10X,'REDUCTION AND DERIVATIVE LATTICE'///)
 1150 FORMAT(2X,'Number of independent problems to study =',I3//)
 1200 FORMAT(2X,'Lattice Matching of reduced cell(s) with the'/2X,&
     &'          NBS  Crystal  Data  File'//)
 1300 FORMAT(1H1)
 1350 FORMAT(////10X,I3,'.  ')
 1400 FORMAT(1H+,16X,'REDUCTION'/)
 1500 FORMAT(1H+,16X,'SUPERLATTICES'/)
 1600 FORMAT(1H+,16X,'SUBLATTICES'/)
 1700 FORMAT(1H+,16X,'SUPERLATTICES for a given delta followed by SUBLAT&
     &TICES for 1/delta'/)
 1800 FORMAT(//33X,'** Initial Cell is Primitive **'/)
 1900 FORMAT(//33X,'** Initial Cell is A-Centered **'/)
 2000 FORMAT(//33X,'** Initial Cell is B-Centered **'/)
 2100 FORMAT(//33X,'** Initial Cell is C-Centered **'/)
 2200 FORMAT(//33X,'** Initial Cell is F-Centered **'/)
 2300 FORMAT(//33X,'** Initial Cell is I-Centered **'/)
 2400 FORMAT(//32X,'** Initial Cell is RH-Centered **'/)
 2500 FORMAT(///2X,'Superlattices for Delta = ',I2)
 2600 FORMAT(///2X,'Sublattices for Delta = 1/',I2)
 2700 FORMAT(//36X,'***** Supercell ',I3,' *****'/)
 2800 FORMAT(//37X,'***** Subcell ',I3,' *****'/)
 2900 FORMAT(1X,'T 1=',4X,3F7.2,'/',3F7.2,'/',3F7.2)
 3000 FORMAT(1X,'T 2=',4X,3F7.2,'/',3F7.2,'/',3F7.2)
 3100 FORMAT(1X,'T 2 INV=',3F7.2,'/',3F7.2,'/',3F7.2)
 3200 FORMAT(1H+,108X,'** Cell Matrix **'/)
 3300 FORMAT(1X,'CELL  1=',3F10.4,1X,3F10.3,14X,'V1=',F10.2,5X,3F10.3)
 3400 FORMAT(1X,'CELL  2=',3F10.4,1X,3F10.3,14X,'V2=',F10.2,5X,3F10.3)
 3500 FORMAT(////1X,10X,'CELL TRANSFORMATION'/)
 3600 FORMAT(//33X,'** Cell Transformation ',I2,' **'/)
 3700 FORMAT(////1X,10X,'MATRIX INVERSION'/)
 3800 FORMAT(//33X,'** Matrix Inversion ',I2,' **'/)
 end subroutine outpt1
 subroutine qmatri

 use nbsLatticeMod, only : idel1,idel2,iqmatf,isq11,isq12,&
                             & isq13,isq22,isq23,isq33

 implicit none

 integer :: i,idel,ii,j,jj,k,kk,iprod,iq11,iq12,iq13,iq21,iq22,&
          & iq23,iq31,iq32,iq33

!     --- IDEL IS THE VALUE OF THE DETERMINANT AND IMATT IS THE TOTAL
!         NON-EQUIVALENT MATRICES CONSISTENT WITH A GIVEN IDEL
!     --- CALCULATE ALL THE UNIQUE MATRICES (IQMATF) CONSISTENT WITH
!         WITH A GIVEN DELTA (IDEL) AND STORE THEM IN THE ARRAYS.
!         EACH UNIQUE MATRIX IS IN THE FORM :
!               Q11     Q12     Q13
!               0       Q22     Q23
!               0       0       Q33

!     --- A TOTAL OF 519 UNIQUE MATRICES WILL BE GENERATED FOR
!         DELTAS 2 THROUGH 9 (I.E. 7+13+35+31+91+57+155+130 = 519).

      DO 300  IDEL = IDEL1,IDEL2
         DO 200 I = 1,IDEL
         DO 200 J = 1,IDEL
         DO 200 K = 1,IDEL
            IPROD = I*J*K
            IF(IPROD.NE.IDEL) GO TO 200

!              --- SET UP DIAGONAL ELEMENTS
               IQ11 = I
               IQ22 = J
               IQ33 = K

!              --- SET BOTTOM TRIANGULAR ELEMENTS EQUAL TO ZERO
               IQ21 = 0
               IQ31 = 0
               IQ32 = 0

!              --- SET TOP TRIANGULAR ELEMENTS EQUAL TO ZERO
               IQ12 = 0
               IQ13 = 0
               IQ23 = 0

!              --- GENERATE ALL MATRICES CONSISTENT WITH GIVEN DIAGONAL
               DO 100 II = 1,IDEL
               DO 100 JJ = 1,IDEL
               DO 100 KK = 1,IDEL
                  IQ12 = II - 1
                  IQ13 = JJ - 1
                  IQ23 = KK - 1

!                 --- THE VALUE OF EACH NON-DIAGONAL ELEMENT IS RESTRICT
!                     TO BE LESS THAN THE DIAGONAL ELEMENT IN THE SAME
!                     COLUMN.  IN DOING SO, A UNIQUE UPPER TRIANGULAR
!                     MATRIX IS GENERATED FOR A GIVEN DELTA.
                  IF(.NOT.(IQ12.LT.IQ22.AND.IQ13.LT.IQ33.AND.IQ23.LT.&
     &               IQ33)) GO TO 100

!                    --- MATRIX IS ACCEPTABLE
                     IQMATF = IQMATF + 1

!                    --- STOP PROGRAM EXECUTION IF TOO MANY UPPER
!                        TRIANGULAR MATRICES HAVE BEEN GENERATED
!                        (SHOULD NOT OCCUR)
                     IF(IQMATF.GT.519) STOP  '*QMATRI* ... Error in the&
     &generation of upper triangular matrices.'

                     ISQ11(IQMATF) = IQ11
                     ISQ12(IQMATF) = IQ12
                     ISQ13(IQMATF) = IQ13
                     ISQ22(IQMATF) = IQ22
                     ISQ23(IQMATF) = IQ23
                     ISQ33(IQMATF) = IQ33
  100          CONTINUE
  200    CONTINUE
  300 CONTINUE
      RETURN
end subroutine qmatri
      subroutine read

 use nbsLatticeMod, only : ai,alphi,beti,bi,ci,gami,icntr1,idel1,&
                      & idel2,ierr1,irss,iunita,iunitb,voli

!      CHARACTER AHP,AHA,AHB,AHC,AHF,AHI,AHH,AHR,ACNTR1,ACNTR2

!      COMMON /CELLI/ AI,BI,CI,ALPHI,BETI,GAMI,VOLI
!      COMMON /CNTR1/ ICNTR1
!      COMMON /DELTA1/ IRSS
!      COMMON /DELTA2/ IDEL1,IDEL2
!      COMMON /ERR1/ IERR1
!      COMMON /UNIT1/ IUNITA
!      COMMON /UNIT2/ IUNITB

 implicit none

 character(1) acntr1, acntr2
 character(1) :: AHP='P', AHA='A', AHB='B', AHC='C', AHF='F', AHI='I'
 character(1) :: AHR='R', AHH='H'

 interface

  subroutine center
  end subroutine center

  subroutine outpt1(iflag1)
  implicit none
  integer :: iflag1
  end subroutine outpt1

  subroutine volume(ax,bx,cx,alphx,betx,gamx,volx)
  implicit none
  real :: ax,bx,cx,alphx,betx,gamx,volx
  end subroutine volume
 end interface

!     --- READ PROBLEM SPECIFICATIONS AND CELL PARAMETERS
!         (ON RSS CONTROL LINE)
      READ(IUNITA,1000) IRSS, IDEL1, IDEL2, ACNTR1, ACNTR2,&
     &             AI, BI, CI, ALPHI, BETI, GAMI

!     --- IF NOT REDUCTION ONLY, CHECK AND, IF NECESSARY, RESET
!         SPECIFICATIONS FOR DERIVATIVE LATTICE CALCULATIONS
      IF(IRSS.EQ.0) GO TO 100
         IF(IRSS.GT.3.OR.IRSS.LT.0) WRITE(IUNITB,1100)
         IF(IRSS.GT.3.OR.IRSS.LT.0) IRSS = 0
         IF(IDEL1.GT.IDEL2.OR.IDEL1.LT.2.OR.IDEL2.GT.9)&
     &      WRITE(IUNITB,1200)
         IF(IDEL1.GT.IDEL2.OR.IDEL1.LT.2.OR.IDEL2.GT.9) IRSS = 0
  100 CONTINUE

!     --- WRITE HEADINGS FOR RSS PROBLEM SPECIFICATIONS
      CALL OUTPT1(4)

!     --- ASSIGN VARIABLE (ICNTR1) TO INDICATE CENTERING OF INPUT CELL
      ICNTR1 = 0
      IF(ACNTR1.EQ.AHP)  ICNTR1 = 1
      IF(ACNTR1.EQ.AHA)  ICNTR1 = 2
      IF(ACNTR1.EQ.AHB)  ICNTR1 = 3
      IF(ACNTR1.EQ.AHC)  ICNTR1 = 4
      IF(ACNTR1.EQ.AHF)  ICNTR1 = 5
      IF(ACNTR1.EQ.AHI)  ICNTR1 = 6
      IF(ACNTR1.EQ.AHR.AND.ACNTR2.EQ.AHH)   ICNTR1 = 7
      IF(ACNTR1.EQ.AHR.AND.ACNTR2.EQ.AHR)   ICNTR1 = 1
      IF(ICNTR1.EQ.0) WRITE(IUNITB,1300)
      IF(ICNTR1.EQ.0) ICNTR1 = 1

!     --- WRITE HEADINGS FOR CENTERING OF INPUT CELL
      CALL OUTPT1(5)

!     --- FIND MATRIX TO TRANSFORM THE INPUT CELL TO
!         A PRIMITIVE CELL OF THE LATTICE
      CALL CENTER

!     --- CALCULATE THE VOLUME OF THE INPUT CELL
!         (AND CHECK FOR VALID INPUT DATA)
      CALL VOLUME(AI,BI,CI,ALPHI,BETI,GAMI,VOLI)

!     --- IF INPUT CELL HAS ILLEGAL CELL PARAMETER(S) AND/OR CELL VOLUME
!         WRITE ERROR MESSAGE (ERROR FLAG IERR1 IS SET IN *VOLUME*)
      IF(IERR1.EQ.1) WRITE(IUNITB,1400)
      RETURN
 1000 FORMAT(I1,2X,2I1,3X,2A1,6F10.2)
 1100 FORMAT(/1X,'*READ* ERROR ... Check RSS Control Line, column 1.'&
     &/1X,17X,'Only cell reduction calculations will be carried out.'//)
 1200 FORMAT(/1X,'*READ* ERROR ... Check RSS Control Line, columns 4-5.'&
     &/1X,17X,'Only cell reduction calculations will be carried out.'//)
 1300 FORMAT(/1X, '*READ* WARNING ... Input centering symbol not legal (&
     &i.e. P,A,B,C,I,F,RR,RH).',&
     &/20X,'Program assumes that the input cell is primitive.')
 1400 FORMAT(/1X,'*READ* ERROR ... Input cell has illegal cell parameter&
     &(s) and/or cell volume.'/)
 end subroutine read
 subroutine reduce
 
 use nbsLatticeMod, only : af,ai,alphf,alphi,betf,beti,&
                         & bf,bi,cf,ci,gamf,gami,ideriv,&
                         & ierr1,ifinis,irss,iprob,iqmatf,&
                         & ntd,var90,volf,voli

!      COMMON /CELLI/ AI,BI,CI,ALPHI,BETI,GAMI,VOLI
!      COMMON /CELLF/ AF,BF,CF,ALPHF,BETF,GAMF,VOLF
!      COMMON /DELTA1/ IRSS
!      COMMON /DELTA4/ IDERIV,IFINIS,IQMATF
!      COMMON /ERR1/ IERR1
!      COMMON /PROB1/ NTD,IPROB
!      COMMON /VAR1/ VAR90

 implicit none

 interface

  subroutine deriv
  end subroutine deriv

  subroutine mncond
  end subroutine mncond

  subroutine normal
  end subroutine normal

  subroutine outpt1(iflag1)
  implicit none
  integer :: iflag1
  end subroutine outpt1

  subroutine qmatri
  end subroutine qmatri

  subroutine read
  end subroutine read
  
  subroutine spcond
  end subroutine spcond

 
  subroutine trans(iprtra)
  implicit none
  integer :: iprtra
  end subroutine trans
  
  
 end interface

!     --- INITIALIZE VARIABLES
      IPROB = 1

  100 CONTINUE

!     --- BEGINNING OF LOOP TO PROCESS EACH INPUT CELL

!     --- WRITE PROBLEM NUMBER
      CALL OUTPT1(3)

!     --- INITIALIZE VARIABLES
      IDERIV = 0
      IQMATF = 0
      IFINIS = 0

!     --- READ AND INTERPRET RSS CONTROL LINE
!         (RSS PROBLEM SPECIFICATIONS, CELL CENTERING, CELL PARAMETERS)
      CALL READ

!     --- GO TO NEXT PROBLEM IF INPUT CELL HAS ILLEGAL CELL PARAMETER(S)
!         AND/OR CELL VOLUME (ERROR FLAG IERR1 IS SET IN *VOLUME*)
      IF(IERR1.EQ.1) GO TO 400

!        --- GENERATE UPPER TRIANGULAR MATRICES IF DERIVATIVE LATTICES
!            ARE TO BE CALCULATED
         IF(IRSS.NE.0) CALL QMATRI

  200    CONTINUE

!        --- STARTING POINT FOR REDUCTION CALCULATIONS

!        --- THIS STATEMENT DOES NOT APPLY TO REDUCTION-ONLY
!            CALCULATIONS ... GO TO NEXT THE PROBLEM IF ALL THE
!            DERIVATIVE LATTICES HAVE BEEN CALCULATED AND REDUCED
!            (IFINIS IS SET TO 1 IN DERIV AFTER THE LAST UPPER
!            TRIANGULAR MATRIX HAS BEEN SELECTED)
         IF(IFINIS.EQ.1) GO TO 400

!           --- FOR EACH INPUT CELL, REDUCE THE CELL AND, IF APPROPRIATE
!               REDUCE ITS ASSOCIATED DERIVATIVE CELLS
            IF(IDERIV.GT.0) CALL DERIV
            CALL TRANS(0)

!           --- VAR90 DEFINES RANGE OF ANGLES TO BE CONSIDERED EQUAL
!               TO 90 DEGREES IN SUBROUTINE 'NORMAL'
            VAR90 = 0.000001
            CALL MNCOND
            VAR90 = 0.000349
            CALL NORMAL
            CALL SPCOND
            CALL TRANS(1)

!           --- IF DERIVATIVE LATTICES ARE NOT BEING CALCULATED,
!               GO TO NEXT PROBLEM
            IF(IRSS.EQ.0) GO TO 400
               IF(IDERIV.NE.0) GO TO 200

!              --- RESET INITIAL CELL TO BE THE REDUCED (PRIMITIVE) CELL
!                  OF THE LATTICE
               AI = AF
               BI = BF
               CI = CF
               ALPHI = ALPHF
               BETI  = BETF
               GAMI  = GAMF
               VOLI  = VOLF

!              --- FIRST UPPER TRIANGULAR MATRIX WILL BE SELECTED
!                  TO CALCULATE DERIVATIVE LATTICES
               IDERIV = 1
               GO TO 200

  400 CONTINUE

!     --- FINISHED RSS CALCULATIONS FOR A SINGLE PROBLEM,
!         PROCEED TO NEXT PROBLEM IF REQUIRED
      IPROB = IPROB + 1
      IF(IPROB.LE.NTD) GO TO 100
      RETURN
       END
      subroutine set
 use nbsLatticeMod, only : t,u,u1,u2,u3,v1,v2,v3,w1,w2,w3

 implicit none
 integer :: i

!     --- SAVE THE CURRENT VALUES OF THE ELEMENTS FOR THE TOTAL
!         TRANSFORMATION MATRIX IN THE ARRAY T AND SET THE MATRIX
!         ELEMENTS IN THE UPDATE MATRIX U TO ZERO
      T(1) = U1
      T(2) = V1
      T(3) = W1
      T(4) = U2
      T(5) = V2
      T(6) = W2
      T(7) = U3
      T(8) = V3
      T(9) = W3
      DO 100 I = 1,9
         U(I) = 0.0
  100 CONTINUE
      RETURN
      END subroutine set
 subroutine shortv(sxx,sxy,elem1,elem5,elem9,elemx)
 use nbsLatticemod, only : ick021,ick022,ick023,ick024,ick025,ick026,ick027,&
                     & ick028,ick029
 implicit none

 real, intent(in) :: sxx,sxy
 real, intent(out) :: elem1,elem5,elem9,elemx 
 real :: axy,rxy,xn,xn1
 interface
 subroutine ckpt02(iflag2)
 implicit none
 integer :: iflag2
 end subroutine ckpt02
 end interface 
 
!**
!     --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE VARIABLES
      IF(ICK026.EQ.1) CALL CKPT02(8)

!     --- INITIALIZE VARIABLES
      AXY = ABS(SXY)
      XN = 0.0
      IF(SXY.LE.0.0) XN =  1.0
      IF(SXY.GT.0.0) XN = -1.0
      XN1 = 0.0

!     --- DETERMINE THE NUMBER OF VECTOR ADDITIONS OR SUBTRACTIONS
!         REQUIRED TO OBTAIN THE SHORTEST Y' RELATIVE TO X.  IF THE
!         INITIAL ANGLE IS < 90 DEGREES, VECTOR X IS SUBTRACTED XN1
!         TIMES, IF THE INITIAL ANGLE IS >= 90 DEGREES, VECTOR X IS
!         ADDED XN1 TIMES.
  100 CONTINUE
      RXY = ABS(SXY + XN*SXX)
      IF((RXY-AXY).GT.0.0) GO TO 160
         IF(SXY.GE.0.0) GO TO 120
            XN  = XN  + 1.0
            XN1 = XN1 + 1.0
            GO TO 140
  120    CONTINUE
         XN  = XN  - 1.0
         XN1 = XN1 - 1.0
  140    CONTINUE
         AXY = RXY
         GO TO 100
  160 CONTINUE

!     --- SET THE MATRIX ELEMENTS FOR THE MATRIX USED TO UPDATE
!         THE TOTAL TRANSFORMATION MATRIX (U)
      CALL SET
      ELEM1 = 1.0
      ELEM5 = 1.0
      ELEM9 = 1.0
      ELEMX = XN1

      RETURN
      END
      subroutine sort1(a1,a2,a3,a4,a5,a6,a7,n)
 implicit none

 integer :: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:),a7(:)
 integer :: i,ii,iim,im,index(900),j,k,m,n,q,temp1,temp2,&
          & temp3,temp4,temp5,temp6,temp7

!     A1, ... A7 = ARRAYS TO BE SORTED
!     N = LENGTH OF ARRAYS A1, ... A7
!     INDEX = ARRAY SPECIFYING ORDER OF INPUT ARRAYS AFTER SORT
!     Q = FLAG TO DETERMINE WHETHER OR NOT TO SORT ARRAY
!         IN THIS SUBROUTINE
!        Q = 0 ... SUBROUTINE SORTS ARRAYS (ALSO CALCULATES INDEX
!                  SPECIFYING HOW ARRAY WAS SORTED)
!        Q.NE.0 .. SUBROUTINE DOES NOT SORT ARRAYS BUT CALCULATES
!                  INDEX SPECIFYING HOW TO SORT

!     *** NOTE THAT THE VARIABLES INDEX AND Q MAY BE
!         INCLUDED IN THE CALL ARGUMENT




!     --- THIS SUBROUTINE RETURNS SORTED ARRAYS
      Q = 0

!     --- INITIALIZE VARIABLES
      DO 10 I = 1,N
         INDEX(I) = I
   10 CONTINUE

      M = N

!     --- BEGINNING OF CODE TO DETERMINE ORDER OF INPUT ARRAYS
!         AFTER SORTING
   20 CONTINUE

!     --- MUST HAVE AT LEAST 2 ITEMS TO SORT
      IF((M-2).LT.0) GO TO 60
         M = 2*(M/4) + 1
         K = N - M
         J = 1

   30    CONTINUE
         I = J

   40    CONTINUE
         IM = I + M
         II = INDEX(I)
         IIM = INDEX(IM)
         IF(A1(II).LT.A1(IIM)) GO TO 50
         IF(A1(II).GT.A1(IIM)) GO TO 45
         IF(A2(II).LT.A2(IIM)) GO TO 50
         IF(A2(II).GT.A2(IIM)) GO TO 45
         IF(A3(II).LT.A3(IIM)) GO TO 50
         IF(A3(II).GT.A3(IIM)) GO TO 45
         IF(A4(II).LT.A4(IIM)) GO TO 50
         IF(A4(II).GT.A4(IIM)) GO TO 45
         IF(A5(II).LT.A5(IIM)) GO TO 50
         IF(A5(II).GT.A5(IIM)) GO TO 45
         IF(A6(II).LT.A6(IIM)) GO TO 50
         IF(A6(II).GT.A6(IIM)) GO TO 45
         IF(A7(II).LE.A7(IIM)) GO TO 50
         IF(A7(II).GT.A7(IIM)) GO TO 45
   45    CONTINUE
            INDEX(I) = IIM
            INDEX(IM) = II
            I = I - M
            IF(I.GT.0) GO TO 40

   50    CONTINUE
         J = J + 1
         IF((K-J).LT.0) GO TO 20
            GO TO 30

   60 CONTINUE

!     --- BEGINNING OF CODE TO SORT THE ARRAYS (IF INDICATED)
      IF(Q.NE.0) GO TO 200

!     --- INDEX(A)=B MEANS THAT THE 'B'TH ENTRY BELONGS IN POSITION 'A'
      DO 80 I = 1,N
         IF(INDEX(I).EQ.0) GO TO 80
         IF((INDEX(I) - I).EQ.0) GO TO 80
            J = I
            TEMP1 = A1(I)
            TEMP2 = A2(I)
            TEMP3 = A3(I)
            TEMP4 = A4(I)
            TEMP5 = A5(I)
            TEMP6 = A6(I)
            TEMP7 = A7(I)

   70       CONTINUE
            K = INDEX(J)
            INDEX(J) = 0
            A1(J) = A1(K)
            A2(J) = A2(K)
            A3(J) = A3(K)
            A4(J) = A4(K)
            A5(J) = A5(K)
            A6(J) = A6(K)
            A7(J) = A7(K)
            J = K
            IF((INDEX(J) - I).NE.0) GO TO 70
               INDEX(J) = 0
               A1(J) = TEMP1
               A2(J) = TEMP2
               A3(J) = TEMP3
               A4(J) = TEMP4
               A5(J) = TEMP5
               A6(J) = TEMP6
               A7(J) = TEMP7

   80 CONTINUE

  200 CONTINUE
      RETURN
      END subroutine sort1
      SUBROUTINE SORT2(A1,A2,A3,N)
 implicit none
 integer :: a1(:),a2(:),a3(:)
 integer :: i,ii,iim,im,j,k,index(2000),m,n,q,temp1,temp2,temp3



!     THIS SUBROUTINE IS THE SAME AS SORT1 EXCEPT ONLY THREE
!       ARRAYS ARE SORTED


!     A1, A2, A3 = ARRAYS TO BE SORTED
!     N = LENGTH OF ARRAYS A1, A2, A3
!     INDEX = ARRAY SPECIFYING ORDER OF INPUT ARRAYS AFTER SORT
!     Q = FLAG TO DETERMINE WHETHER OR NOT TO SORT ARRAY
!         IN THIS SUBROUTINE
!        Q = 0 ... SUBROUTINE SORTS ARRAYS (ALSO CALCULATES INDEX
!                  SPECIFYING HOW ARRAY WAS SORTED)
!        Q.NE.0 .. SUBROUTINE DOES NOT SORT ARRAYS BUT CALCULATES
!                  INDEX SPECIFYING HOW TO SORT

!     *** NOTE THAT THE VARIABLES INDEX AND Q MAY BE
!         INCLUDED IN THE CALL ARGUMENT




!     --- THIS SUBROUTINE RETURNS SORTED ARRAYS
      Q = 0

!     --- INITIALIZE VARIABLES
      DO 10 I = 1,N
         INDEX(I) = I
   10 CONTINUE

      M = N

!     --- BEGINNING OF CODE TO DETERMINE ORDER OF INPUT ARRAYS
!         AFTER SORTING
   20 CONTINUE

!     --- MUST HAVE AT LEAST 2 ITEMS TO SORT
      IF((M-2).LT.0) GO TO 60
         M = 2*(M/4) + 1
         K = N - M
         J = 1

   30    CONTINUE
         I = J

   40    CONTINUE
         IM = I + M
         II = INDEX(I)
         IIM = INDEX(IM)
         IF(A1(II).LT.A1(IIM)) GO TO 50
         IF(A1(II).GT.A1(IIM)) GO TO 45
         IF(A2(II).LT.A2(IIM)) GO TO 50
         IF(A2(II).GT.A2(IIM)) GO TO 45
         IF(A3(II).LE.A3(IIM)) GO TO 50
         IF(A3(II).GT.A3(IIM)) GO TO 45
   45    CONTINUE
            INDEX(I) = IIM
            INDEX(IM) = II
            I = I - M
            IF(I.GT.0) GO TO 40

   50    CONTINUE
         J = J + 1
         IF((K-J).LT.0) GO TO 20
            GO TO 30

   60 CONTINUE

!     --- BEGINNING OF CODE TO SORT THE ARRAYS (IF INDICATED)
      IF(Q.NE.0) GO TO 200

!     --- INDEX(A)=B MEANS THAT THE 'B'TH ENTRY BELONGS IN POSITION 'A'
      DO 80 I = 1,N
         IF(INDEX(I).EQ.0) GO TO 80
         IF((INDEX(I) - I).EQ.0) GO TO 80
            J = I
            TEMP1 = A1(I)
            TEMP2 = A2(I)
            TEMP3 = A3(I)

   70       CONTINUE
            K = INDEX(J)
            INDEX(J) = 0
            A1(J) = A1(K)
            A2(J) = A2(K)
            A3(J) = A3(K)
            J = K
            IF((INDEX(J) - I).NE.0) GO TO 70
               INDEX(J) = 0
               A1(J) = TEMP1
               A2(J) = TEMP2
               A3(J) = TEMP3
   80 CONTINUE

  200 CONTINUE
      RETURN
 end subroutine sort2
  subroutine spcon2(icond,exp1,exp2,exp3,ispfix,ele1,ele2,ele3,ele4,&
     &                  ele5,ele6)

 use nbsLatticeMod, only:  ick021,ick022,ick023,ick024,ick025,ick026,ick027,&
                    & ick028,ick029

 implicit none

 integer :: icond,ispfix
 real :: pc2=0.0005
 real :: ele1,ele2,ele3,ele4,ele5,ele6,exp1,exp2,exp3

 interface
  subroutine ckpt02(iflag2)
  implicit none
  integer :: iflag2
  end subroutine ckpt02
  subroutine multip
  end subroutine multip
  subroutine normal
  end subroutine normal
  subroutine trans(iprtra)
  implicit none
  integer :: iprtra
  end subroutine trans
 
 end interface
!**
!      COMMON /CK02/ ICK021,ICK022,ICK023,ICK024,ICK025,ICK026,ICK027,
!     $              ICK028,ICK029
!      DATA PC2/0.0005/


!     --- WHEN A SPECIAL RELATIONSHIP IN THE CELL MATRIX OCCURS,
!         CHECK, AND IF NECESSARY SATISFY, THE SPECIAL CONDITION
      IF(EXP1.GT.PC2) GO TO 100
!**
!        --- FOR CHECKING, WRITE EXECUTION POINT, SPECIAL CONDITION
!            NUMBER, INTERMEDIATE VARIABLES
         IF(ICK027.EQ.1) CALL CKPT02(30+ICOND)

         IF(EXP2.GE.ABS(EXP3)) GO TO 100

!           --- FAILED ... NOW SATISFY SPECIAL CONDITION
            ELE1 = -1.0
            ELE2 = -1.0
            ELE3 = -1.0
            ELE4 =  1.0
            ELE5 =  1.0
            ELE6 =  1.0

!           --- UPDATE THE TOTAL TRANSFORMATION MATRIX AND
!               APPLY THE RESULTING MATRIX TO THE INPUT CELL
            CALL MULTIP
            CALL TRANS(0)
!**
!           --- FOR CHECKING, WRITE EXECUTION POINT AND INTERMEDIATE
!               VARIABLES
            IF(ICK028.EQ.1) CALL CKPT02(11)

            CALL NORMAL
            ISPFIX = 1
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SPCOND
 use nbsLatticeMod, only : cosb,cosg,itype,iunitb,s11,s12,s13,&
                         & s22,s23,s33,u,var90


!      COMMON /COSANG/ COSA,COSB,COSG
!      COMMON /DOTP/ S11,S22,S33,S23,S13,S12
!      COMMON /MATR2/ U(9),T(9)
!      COMMON /TYPE/ ITYPE
!      COMMON /UNIT2/ IUNITB
!      COMMON /VAR1/ VAR90
!      DATA PC1/0.00006/

 implicit none
 
 integer :: icond,icount,ispfix
 real :: dum1,dum2,dum3,ts1,ts2,ts3,ts4,ts5,ts6,ts7,tem1,tem2
  real ::   pc1=0.00006

 interface
  subroutine set
  end subroutine set
  subroutine spcon2( icond,exp1,exp2,exp3,ispfix,ele1,ele2,ele3,ele4,&
                   & ele5,ele6)
  implicit none
 integer :: icond,ispfix
 real :: ele1,ele2,ele3,ele4,ele5,ele6,exp1,exp2,exp3
  end subroutine spcon2
 end interface

!     --- INITIALIZE VARIABLES
      ICOUNT = 0
      DUM1 = 0.0
      DUM2 = 0.0
      DUM3 = 0.0

  100 CONTINUE

!     --- STARTING POINT TO CHECK, AND IF NECESSARY SATISFY, EACH
!         SPECIAL CONDITION

!     --- INITIALIZE VARIABLES
      CALL SET
      ISPFIX = 0

!     --- PROTECTIVE CODE ... THIS CODE SHOULD BE ACTIVE ONLY
!         WHEN THE PROGRAM CONSTANT (PC1) IS INCORRECTLY SET
!         FOR THE SPECIFIC COMPUTER
      ICOUNT = ICOUNT + 1
      IF(ICOUNT.LE.20) GO TO 200
         WRITE(IUNITB,1000)
         GO TO  600
  200 CONTINUE

!     --- CALCULATE TEMPORARY VARIABLES USED TO CHECK SPECIAL CONDITIONS
      TS1  = ABS(S11 - S22)/S11
      TS2  = ABS(S22 - S33)/S22
      TS3  = ABS(ABS(2.0*S23) - S22)/S22
      TS4  = ABS(ABS(2.0*S13) - S11)/S11
      TS5  = ABS(ABS(2.0*S12) - S11)/S11
      TEM1 = S11 + S22
      TEM2 = 2.0*(ABS(S23) + ABS(S13) + ABS(S12))
      TS6  = ABS(TEM1 - TEM2)/TEM1
      TS7  = 2.0*ABS(S13) + ABS(S12) + PC1


!     --- SPECIAL CONDITION (A) ...
!         TYPE 1 CELL - IF A.A = B.B THEN  B.C  <=  A.C
!         TYPE 2 CELL - IF A.A = B.B THEN /B.C/ <= /A.C/
      ICOND = 1
      CALL SPCON2(ICOND,TS1,ABS(S13)+PC1,S23,ISPFIX,U(2),U(4),U(9),&
     &            DUM1,DUM2,DUM3)
      IF(ISPFIX.EQ.1) GO TO 100

!     --- SPECIAL CONDITION (B) ...
!         TYPE 1 CELL - IF B.B = C.C THEN  A.C  <=  A.B
!         TYPE 2 CELL - IF B.B = C.C THEN /A.C/ <= /A.B/
      ICOND = 2
      CALL SPCON2(ICOND,TS2,ABS(S12)+PC1,S13,ISPFIX,U(1),U(6),U(8),&
     &            DUM1,DUM2,DUM3)
      IF(ISPFIX.EQ.1) GO TO 100

!     --- SPECIAL CONDITION (C) ...
!         TYPE 1 CELL - IF  B.C  = 1/2B.B THEN A.B <= 2A.C
!         TYPE 2 CELL - IF /B.C/ = 1/2B.B THEN A.B = 0
      ICOND = 3
      IF(ITYPE.EQ.1) CALL SPCON2(ICOND,TS3,2.0*S13+PC1,S12,ISPFIX,&
     &                           U(1),U(5),U(8),U(9),DUM1,DUM2)
      IF(ITYPE.EQ.2) CALL SPCON2(ICOND,TS3,VAR90,COSG,ISPFIX,&
     &                           U(5),U(8),U(9),U(1),DUM1,DUM2)
      IF(ISPFIX.EQ.1) GO TO 100

!     --- SPECIAL CONDITION (D) ...
!         TYPE 1 CELL - IF  A.C  = 1/2A.A THEN A.B <= 2B.C
!         TYPE 2 CELL - IF /A.C/ = 1/2A.A THEN A.B = 0
      ICOND = 4
      IF(ITYPE.EQ.1) CALL SPCON2(ICOND,TS4,2.0*S23+PC1,S12,ISPFIX,&
     &                           U(1),U(5),U(7),U(9),DUM1,DUM2)
      IF(ITYPE.EQ.2) CALL SPCON2(ICOND,TS4,VAR90,COSG,ISPFIX,&
     &                           U(1),U(7),U(9),U(5),DUM1,DUM2)
      IF(ISPFIX.EQ.1) GO TO 100

!     --- SPECIAL CONDITION (E) ...
!         TYPE 1 CELL - IF  A.B  = 1/2A.A THEN A.C <= 2B.C
!         TYPE 2 CELL - IF /A.B/ = 1/2A.A THEN A.C = 0
      ICOND = 5
      IF(ITYPE.EQ.1) CALL SPCON2(ICOND,TS5,2.0*S23+PC1,S13,ISPFIX,&
     &                           U(1),U(4),U(9),U(5),DUM1,DUM2)
      IF(ITYPE.EQ.2) CALL SPCON2(ICOND,TS5,VAR90,COSB,ISPFIX,&
     &                           U(1),U(4),U(5),U(9),DUM1,DUM2)
      IF(ISPFIX.EQ.1) GO TO 100

!     --- SPECIAL CONDITIONS (A)-(E) HAVE BEEN SATISFIED.
!         REDUCTION IS COMPLETE FOR A TYPE 1 CELL.
      IF(ITYPE.EQ.1) GO TO 600

!        --- SPECIAL CONDITION (F) ... APPLIES ONLY FOR TYPE 2 CELL
!            IF (/B.C/ + /A.C/ + /A.B/) = 1/2(A.A + B.B) ,
!            THEN A.A <= (2/A.C/ + A.B)
         ICOND = 6
         CALL SPCON2(ICOND,TS6,TS7,S11,ISPFIX,&
     &              U(1),U(5),DUM1,U(7),U(8),U(9))
         IF(ISPFIX.EQ.1) GO TO 100


  600 CONTINUE
      RETURN
 1000 FORMAT(/1X,'*SPCOND* WARNING ... Cell may not be reduced - Check a&
     &ll conditions for reduction.'/1X,21X,'Program error when satisfyin&
     &g Special Conditions.'&
     &/1X,21X,'The program constant (PC1) may be incorrectly set for thi&
     &s computer.'/)
      END
 subroutine trans(iprtra)

 use nbsLatticeMod, only : af,ai,alphf,alphi,betf,beti,bf,bi,&
                          & cf,ci,cosa,cosb,cosg,gamf,gami,ick021,&
                          & ick022,ick023,ick024,ick025,ick026,ick027,&
                          & ick028,ick029,ierr2,ipch,iunitb,iunitd,&
                          & nunk,radian,s11,s22,s33,s23,s13,s12,&
                          & u1,u2,u3,v1,v2,v3,volf,voli,w1,w2,w3
 implicit none

 integer :: iprtra
 real :: ar, br, gr

 interface
  subroutine determ
  end subroutine determ
  subroutine dot(x1,x2,x3,y1,y2,y3,xi,yi,&
                & zi,ar,br,gr,dotxy)
 implicit none
 real, intent(in) :: ar,br,gr,xi,x1,x2,x3,yi,y1,y2,y3,zi
 real, intent(out) :: dotxy
  end subroutine dot
  subroutine invers(iprinv)
  implicit none
  integer :: iprinv
  end subroutine invers
  subroutine outpt1(iflag1)
  implicit none
  integer :: iflag1
  end subroutine outpt1
 end interface


!      COMMON /CELLI/ AI,BI,CI,ALPHI,BETI,GAMI,VOLI
!      COMMON /CELLF/ AF,BF,CF,ALPHF,BETF,GAMF,VOLF
!      COMMON /CONST1/ RADIAN
!      COMMON /COSANG/ COSA,COSB,COSG
!      COMMON /DOTP/ S11,S22,S33,S23,S13,S12
!      COMMON /ERR2/ IERR2
!      COMMON /MATR1/ U1,V1,W1,U2,V2,W2,U3,V3,W3
!      COMMON /PROB2/ IPCH,NUNK
!      COMMON /UNIT2/ IUNITB
!      COMMON /UNIT4/ IUNITD
!**
!      COMMON /CK02/ ICK021,ICK022,ICK023,ICK024,ICK025,ICK026,ICK027,
!     $              ICK028,ICK029



!     --- CALCULATE DETERMINANT AND CHECK FOR A VALID
!         TRANSFORMATION MATRIX. FOR THE LM,RSS PROGRAM FUNCTIONS, THE
!         PROGRAM EXECUTION IS STOPPED IN *DETERM* IF THE DETERMINANT
!         IS NOT GREATER THAN ZERO (SHOULD NOT OCCUR)
      CALL DETERM
      IF(IERR2.EQ.0) GO TO 100

!        --- INVALID MATRIX WHEN EXECUTING TRANS PROGRAM FUNCTION
!            (DETERMINANT OF MATRIX IS ZERO, >= 100, OR <= -100),
!            WRITE TRANSFORMATION MATRIX, ERROR MESSAGE,
!            AND GO TO NEXT PROBLEM
         CALL OUTPT1(11)
         WRITE(IUNITB,1000)
         GO TO 400
  100 CONTINUE

      AR = ALPHI/RADIAN
      BR = BETI/RADIAN
      GR = GAMI/RADIAN

!     --- CALCULATE A
      CALL DOT (U1,V1,W1,U1,V1,W1,AI,BI,CI,AR,BR,GR,S11)
      AF = SQRT(S11)

!     --- CALCULATE B
      CALL DOT (U2,V2,W2,U2,V2,W2,AI,BI,CI,AR,BR,GR,S22)
      BF = SQRT(S22)

!     --- CALCULATE C
      CALL DOT (U3,V3,W3,U3,V3,W3,AI,BI,CI,AR,BR,GR,S33)
      CF = SQRT(S33)

!     --- CALCULATE ANGLE BETWEEN AF AND BF (GAMMA)
      CALL DOT (U1,V1,W1,U2,V2,W2,AI,BI,CI,AR,BR,GR,S12)
      COSG = S12/(AF*BF)
      GAMF = (ACOS(COSG))*RADIAN

!     --- CALCULATE ANGLE BETWEEN AF AND CF (BETA)
      CALL DOT (U1,V1,W1,U3,V3,W3,AI,BI,CI,AR,BR,GR,S13)
      COSB  = S13/(AF*CF)
      BETF = (ACOS(COSB))*RADIAN

!     --- CALCULATE ANGLE BETWEEN BF AND CF (ALPHA)
      CALL DOT (U2,V2,W2,U3,V3,W3,AI,BI,CI,AR,BR,GR,S23)
      COSA = S23/(BF*CF)
      ALPHF = (ACOS(COSA))*RADIAN

!     --- CALCULATE VOLUME
      VOLF = AF*BF*CF*SQRT(1.0 - COSG**2 - COSB**2 - COSA**2 +&
     &       2.0*COSA*COSB*COSG)

      IF(IPRTRA.EQ.0) GO TO 400

!        --- WRITE OUTPUT
         CALL OUTPT1(11)
         CALL INVERS(1)
         CALL OUTPT1(13)

!        --- OPTIONAL WRITE OF REDUCED CELLS ON IUNITD
!            (ONLY FOR LATTICE MATCHING)
         IF(IPCH.EQ.1) WRITE(IUNITD,1100) AF, BF, CF, ALPHF,&
     &      BETF, GAMF, VOLF
         IF(IPCH.EQ.1) NUNK = NUNK + 1
!**
!        --- ONLY FOR SPECIAL CHECK RUN, WRITE REDUCED CELL,
!            TRANSFORMATION MATRIX FROM INITIAL CELL TO REDUCED
!            CELL, PROBLEM SEQUENCE NUMBER ON IUNITD
         IF(ICK029.NE.1) GO TO 400
            NUNK = NUNK + 1
            WRITE(IUNITD,1200) AF,BF,CF,ALPHF,BETF,GAMF,&
     &            U1,V1,W1,U2,V2,W2,U3,V3,W3, NUNK
  400 CONTINUE
      RETURN
 1000 FORMAT(/1X,'*TRANS* ERROR ... Invalid matrix, check determinant.'/&
     &)
 1100 FORMAT(6F10.5,F10.2)
!**
 1200 FORMAT(6F10.5,2X,9F7.2,2X,I5)
       END
      SUBROUTINE VOLUME(AX,BX,CX,ALPHX,BETX,GAMX,VOLX)

 use nbsLatticeMod, only: radian, ierr1, iunitb

 implicit none

 real :: alphx,ax,betx,bx,cos1,cos2,cos3,cx,gamx,tem,volx
 

!      COMMON /CONST1/ RADIAN
!      COMMON /ERR1/ IERR1
!      COMMON /UNIT2/ IUNITB



!     --- INITIALIZE VARIABLES
      VOLX = 0.0
      IERR1 = 0

!     --- WHEN ONE OR MORE CELL EDGES ARE <= 0.01,
!         SET ERROR FLAG AND DO NOT CALCULATE VOLUME
      IF(AX.GT.0.01.AND.BX.GT.0.01.AND.CX.GT.0.01) GO TO 100
         IERR1 = 1
         GO TO 400
  100 CONTINUE

!     --- WHEN ONE OR MORE CELL ANGLES ARE <= 0.1 OR >=180.0,
!         SET ERROR FLAG AND DO NOT CALCULATE VOLUME
      IF(ALPHX.GT.0.1.AND.ALPHX.LT.180.0.AND.BETX.GT.0.1.AND.BETX.LT.&
     &   180.0.AND.GAMX.GT.0.1.AND.GAMX.LT.180.0) GO TO 200
         IERR1 = 1
         GO TO 400
  200 CONTINUE

!      --- CALCULATE THE COSINES OF THE CELL ANGLES
      COS1 = COS(ALPHX/RADIAN)
      COS2 = COS(BETX/RADIAN)
      COS3 = COS(GAMX/RADIAN)

!     --- IF THE INTERMEDIATE EXPRESSION (TEM) IS GREATER THAN ZERO,
!         CALCULATE THE CELL VOLUME ... OTHERWISE, SET ERROR FLAG
!         SINCE DO NOT HAVE A LEGAL CELL
      TEM  = (1.0 - COS1**2 - COS2**2 - COS3**2 + 2.0*COS1*COS2*COS3)
      IF(TEM.GT.0.0) GO TO 300
         IERR1 = 1
         GO TO 400
  300 CONTINUE
      VOLX = AX*BX*CX*SQRT(TEM)

!     --- WRITE WARNING MESSAGE IF CALCULATED CELL VOLUME IS
!         EITHER VERY SMALL OR VERY LARGE
      IF(VOLX.LT.8.0.OR.VOLX.GT.50000.0) WRITE(IUNITB,1000)
  400 CONTINUE
      RETURN
 1000 FORMAT(/1X,'*VOLUME* WARNING ... Unusual cell volume has been calc&
     &ulated.'/)
 end subroutine volume
