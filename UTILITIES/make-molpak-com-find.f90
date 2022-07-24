      PROGRAM MAKE_MOLPAK_COM_FIND              ! 1/21/03
!
!----Read the 1st MOLPAK input file and prepare a file with the
!     FIND instruction for a fine grid search. 
!
       IMPLICIT NONE
!
       INTEGER, PARAMETER :: NSET = 60000      ! ~39**3
       CHARACTER (4) :: J(20)
       CHARACTER (1) :: JUNK, FNAME(20) 
       CHARACTER (20) :: FNAME20, FN2
       CHARACTER (4) :: FNAME7
       CHARACTER (6) :: FN1, INPUT
       REAL :: EANGLE(NSET,3), VOL(NSET), ZNTANG(3), ANGLE(NSET,3) 
       INTEGER :: IMINA(NSET,3), IMAXA(NSET,3), NMOL, I, N, IMOL
       INTEGER :: NSOL, II, NS
       REAL :: SANG(3), STEPS(3)  
       LOGICAL INDEX
       DATA INPUT /'input.'/
       DATA FNAME /20*' '/
       EQUIVALENCE  (FNAME20,FNAME(1),FN1,FN2), (FNAME(7),FNAME7) 
!
       WRITE (16,566)
       WRITE (*,566)
566    FORMAT (' Start making individual MOLPAK input files ')
       DO I=1,3
           READ (15,*)   ! discard 1st 3 lines of ....data file
       ENDDO
!----If this a a rerun of MOLPAK+WMIN, NSOL is the solution number
       READ (15,155) NMOL, (ZNTANG(I),I=1,3), NSOL 
155       FORMAT (I5,3F5.1,I5)    
       DO I=1,3
           SANG(I) = ZNTANG(I)  
       END DO
       IF (NSOL .GT. 0) WRITE (16,53) NSOL
53        FORMAT (' Elaborate soln #',I4)
       write (*,53) nsol                                ! <-----$$$$$$$$$$
       WRITE (16,55) ((ZNTANG(I)/2),I=1,3) ! fine step = initial/2
       WRITE (*,55) ((ZNTANG(I)/2),I=1,3) ! fine step = initial/2   <----$$$$$$
55        FORMAT(' Initial fine step increments for the three angles =',  &
     &           2(F5.1,','),',',F5.1)     
!----Read 3 angles from volume.find file created after run second MOLPAK run
      INDEX = .false.
      OPEN (UNIT=7, FILE='volume.find', STATUS='OLD', ERR=900 )
      INDEX = .true.
      N = 0
30    N = N + 1                                       ! read all angles
      READ (7,35,END=100) (ANGLE(N,I),I=1,3), VOL(N)   ! and volumes
35       FORMAT(9X,3F8.1,F10.3)
      GO TO 30
 900  CLOSE (UNIT=7)
       OPEN (UNIT = 8, FILE = 'volume.min', STATUS = 'OLD')
!----Read the minimum volume file produced by 1st MOLPAK run
       READ (8,*)   ! discard header line 
       N = 0
!----Read 2nd line and pick up the 1st MOLPAK angle steps.  
       READ (8,9) (STEPS(I),I=1,3)
9         FORMAT (9X,3F8.1) 
10     N = N + 1
       READ (8,11,END=100) NS, (EANGLE(N,I),I=1,3), VOL(N)
11         FORMAT (1X,I8,3F8.1,F10.3)
       DO I=1,3
           ANGLE(N,I) = EANGLE(N,I)    
       END DO
       GO TO 10
100    N = N - 1
       CLOSE (UNIT = 8)
!----Which is smaller...N or NMOL
       IF (N .LT. NMOL) NMOL = N
       WRITE (16,577) NMOL
577    FORMAT (' Number of MOLPAK runs =',I4)                
!----Now read the 1st molpak input file (unit # 10) and modify it for
!     the FIND minimum volume point search.
       DO IMOL=1,NMOL          ! cycle thru NMOL solutions
           FN1 = INPUT
          IF(IMOL .GE. 1000) THEN           ! create
              WRITE (FNAME7,'(I4)') IMOL    ! numerous
          ELSE
           IF(IMOL .GE. 100) THEN           ! create
              WRITE (FNAME7,67) IMOL        ! numerous
67            FORMAT (I3)                   ! MOLPAK 
           ELSE                             ! input
              IF(IMOL .GE. 10) THEN         ! files
                  WRITE (FNAME7,65) IMOL    ! named
65                FORMAT(I2)                ! input.nnn
              ELSE
                  WRITE (FNAME(7),66) IMOL
66                FORMAT(I1)
              ENDIF
           ENDIF
          END IF
           OPEN (UNIT=26, FILE=FN2, STATUS='NEW') ! open input.nnn
           OPEN (UNIT=10, FILE='fort.10', STATUS='OLD')   ! 10 contains
           READ (10,1) J      ! various misc information; header line       
           WRITE (26,1) J
           WRITE (26,3)
3             FORMAT ('LIST 1')      ! input instructions and FIND vol sequence
20         READ (10,1,END=999) J   ! read lines
1             FORMAT (20A4)
           IF (J(1) .EQ. 'VOLS') GO TO 20   ! eliminate the VOLS line
           IF (J(1) .EQ. 'SEEK') THEN       ! SEEK line...replace with FIND
               IF (INDEX) THEN    ! true...have volume.find
                   IF(NMOL.EQ.1.AND.NSOL.NE.0) THEN   ! rerun M-W with the NSOL solution
                       WRITE (26,45) (ANGLE(NSOL,I),I=1,3), (STEPS(II),II=1,3) 
 45                    FORMAT('FIND 0',6F6.1)    
                   ELSE
                       WRITE (26,45) (ANGLE(IMOL,I),I=1,3), (STEPS(II),II=1,3)  
                   END IF
               ELSE
                  IF(NMOL.EQ.1.AND.NSOL.NE.0) THEN   ! rerun M-W with the NSOL solution
                      WRITE (26,5)  (ANGLE(NSOL,I),I=1,3), (STEPS(II),II=1,3)      ! 6-12-02
 5                    FORMAT ('FIND 2',6F6.1)   
                  ELSE
                      WRITE (26,5)  (ANGLE(IMOL,I),I=1,3), (STEPS(II),II=1,3)      ! 6-12-02
                  END IF
               END IF
           ELSE                             
               WRITE (26,1) J
           ENDIF
       GO TO 20
999    CLOSE (UNIT=26, STATUS='KEEP')
!      N = 0
       CLOSE (UNIT=10)   ! close misc info file
       ENDDO          ! done with NMOL solutions       
       WRITE (*,510) NMOL
510        FORMAT (' Finished producing ',I3,' MOLPAK input files',  &
     &             ' for small angle step refinement')
       STOP
       END PROGRAM MAKE_MOLPAK_COM_FIND
