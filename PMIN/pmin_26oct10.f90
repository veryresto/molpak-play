!---Program for crystal lattice energy minimization by optimization of the 
!    unit cell parameters and positions and orientations of the  
!    molecules...by Rosenbrock step search and least-squares.
!---Conformational refinement concomitant with the cell refinment...the
!    intrmolecular energies (PM3 or b3lyl/631g*) are obtained by external 
!    calls to G03. 
!    
!---Multiple rigid bodies introduced 
!    Rosenbrock step search (RSS) introduced , primarily by Nicolae Albu
!    Conformational optimization for minimum E (rotation of bond(s)) introduced 
!    Major developments by Satya M. Prasad.
!    Major tweaking, testing and general trouble shooting by Zuyue Du and Herman Ammon.
!
!    *******  The latest version made on OCT. 2010.  **********
!
      PROGRAM PMIN   
!
      USE PMIN_MODULE , ONLY:NTX_min , NTY_min , NTZ_min , NTX_max ,    &
     &    NTY_max , NTZ_max , NMOl , NATm , IDX , ISYstem , N , ILS ,   &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , STHl , SYM , Q , A1 , B1 , C1 , CK , WT_mol ,   &
     &    DMAx , DDMax , GNI , CE12 , VE12 , K1 , N11 , IHKl , IRIj ,   &
     &    IRSs_call , I_Cross , A_Cross , B_Cross , C_Cross , ANAme1 ,  &
     &    ANAme2 ,  PI , DIFfpc , PRDifpc , TWOPI , ispl, ispl_1,       &  ! from S.P Feb.5, 09
     &    ANG_TO_RADIAN , NCYc_rss, cell_input, f_factor_rss, vab_name, &  ! 1-8-09 
     &    bond_type, dd_cross                                              ! 5-4-10 DU
      USE RSS1_MODULE
      USE RSS2_MODULE
      USE RSS3_MODULE
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!
! Local variables
!
      REAL*8 :: a , alpha , anij , are = 0.00001 , b , beta , c , c11 , &
     &          c22 , c33 , cmpr = 0.1 , del_egymol_rss , dens_diff ,   &
     &          dens_diff_per , dens_ini , dens_last , dgan , ec ,      &
     &          egymol , egymol_ori , elim = 0.0 , emol_1st ,           &
     &          emol_diff , er , ev , e_diff , e_diff_ls_test = 1.0D-9 ,&
     &          e_nmol , gamma , param_fract = 1.0 , sah , sin1 , sin2 ,&
     &          sin3 , sins , slim = 0.000001 , step1 , step2 , step3 , &
     &          step4 , sumderiv_1st , sumx , sumx_new , sumy ,         &
     &          sumy_new , sumz , sumz_new , sum_deriv1, sum_deriv_previous2, & ! 7/2/09
     &          sum_deriv_previous = 0.0 , vol , vol_diff ,             &
     &          vol_diff_per , vol_ini , w01 , wt_mol_asym , zz,        &
     &          ec_ini,ev_ini,er_ini,e_ini,e_nmol_ini,ts,rt,f_factor,   &             
     &          sumxx, sumyy, sumzz, xc1, yc1, zc1, dxy, dyz            ! from S.P Feb. 5,09
      REAL*8 , DIMENSION(5) :: egy_ini               
      REAL*8 , DIMENSION(NMOLD,MATM) :: ag , bg , cg , qg , wg , xg ,   &
     &                                  yg , zg , IWCodeg    ! 1-14-08 DU
      INTEGER :: allow_change , cycle2 , i , ib , icell , iend , ii ,   &
     &           ij , ik , im , imol , ipp , irosen , irun , istart ,   &
     &           istop , iv1 , j , ji , jj , jkk , k , k2 , k55 , ki ,  &
     &           kj , kjk , natom , ncycle , ne , np , nset , nsym2 ,   &
     &           nv , nv1 , nvr , nvu , nvv, coef_table
      REAL*8 , DIMENSION(NVD*NVD) :: an , btr , evec
      CHARACTER(6) , DIMENSION(10,MATM) :: atomg
      REAL*8 , DIMENSION(MATM) :: atw , fx , fy , fz , pc , xo , x_cor ,&
     &                            x_prev , yo , y_cor , y_prev , zo ,   &
     &                            z_cor , z_prev,                       &
     &                            A1_tab, B1_tab, C1_tab   ! 1-5-09 DU
      REAL*8 , DIMENSION(6) :: cell_1st , cell_cor , cell_diff ,        &
     &                         cell_diff_per , cell_ini , cell_last ,   &
     &                         cell_prev
      CHARACTER(6) , DIMENSION(6) :: cell_name , p_name
      REAL*8 , DIMENSION(NVD,6) :: cell_new , cell_new1
      REAL*8 , DIMENSION(NVD) :: corr , corr_sum , deriv1 , dkt , dyc , &
     &                           eval , p , pcmax , pcmin , pd , pdj ,  &
     &                           psv , rt_sum , san , tq , vn , wn , wp
      REAL*8 :: DCOS , DSQRT
      REAL*8 , DIMENSION(NVD,NVD) :: deriv2 , wpp
      REAL*8 , DIMENSION(4) :: ecvr_1st , ecvr_diff , ecvr_last ,       &
     &                         frac_change
      REAL*8 , DIMENSION(500) :: emol_prev , e_prev
      INTEGER , DIMENSION(MATM) :: igroup , katom , mcity
      INTEGER , DIMENSION(20,4) :: jatom
      INTEGER , DIMENSION(4) :: k11
      INTEGER , DIMENSION(NMOLD,MATM) :: k1atom
      CHARACTER(80) :: line
      INTEGER , DIMENSION(100) :: listatm
      LOGICAL :: l_allow_change = .TRUE.
      INTEGER , DIMENSION(20) :: nsym1
      REAL :: REAL
      REAL*8 , DIMENSION(NMOLD) :: xc , yc , zc
      REAL*8 , DIMENSION(NVD,MATM) :: xnew , xnew1 , ynew , ynew1 ,     &
     &                                znew , znew1
!
      character(1) :: flaggy
!
!     DOUBLE PRECISION :: sumderiv_1st,dens_last,dens_diff,vol_last,vol_diff
!
!----The following are defined in data statements above...
!
!     are = Acceptable ratio of eigenvalues.  0.0001 may be satisfactory.  For ill-conditioned problems    ! added
!           this value is used to prevent blow-up by holding poorly determined eigenvectors constant.      ! 7/9/08
!     slim = Scaling limit (eg 0.000001).  Diagonal elements smaller than slim are treated as 0 when       ! HLA
!            normal equations are scaled
!     cmpr = Constant for Marquardt's compromise.  Added to diagonal elements of the scaled matrix of
!            normal equations (otherwise are unity).  Has the effect of damping an ill-conditioned
!            problem by added cmpr to each eigenvalue before its reciprocal is taken.  Reduces the
!            contribution of each eigenvalue, especially those with small eigenvalues.  Non-0 values
!            up to 0.1 may provide suitable damping.
!
      DATA cell_name/'a     ' , 'b     ' , 'c     ' ,& 
                     'cos_1 ' , 'cos_2 ' , 'cos_3 '/  
      DATA p_name/'RX    ' , 'RY    ' , 'RZ    ' , &
                  'TX    ' , 'TY    ' , 'TZ    '/
!
      OPEN (UNIT=7,FILE='PMIN.inp',STATUS='OLD')
      OPEN (UNIT=12,FILE='PMIN.summary_last',STATUS='UNKNOWN')
      OPEN (UNIT=13,FILE='PMIN.save',STATUS='UNKNOWN')
      OPEN (UNIT=23,FILE='pmin_LSQ.save',STATUS='UNKNOWN')
!
      del_egymol_rss = 0.0                           ! 9-30-08 DU
      DEL_egymol = 0.0                               ! 9-30-08 DU
!----Line # 1 : title line
!
      READ (7,99055) TITle
      WRITE (12,1281) TITle
1281     format (' PMIN.inp line # 1, title : ',a)
      WRITE (13,99056) TITle
      WRITE (61,99056) TITle                         ! added, 7/4/08 HLA
!
!----Line # 2 : principal control parameters (iswitch, irosen, iuser_E, isystem, nmol, sthl,
!                                             ck, dmax, ddmax, cross_term, allow_change,
!                                             coef_table)
!
!       iswitch = 1, usual wmin.input format
!       iswitch = 2, if PMIN.inp has everything and wmin.input file is not used
!
!       irosen = 0, if Rosenbrock step search is to be bypassed.
!       irosen = 1, if only Rosenbrock step search is to be done
!       irosen = 2, if Rosenbrock step search is to be followed by l.s. refinement (typical)
!
!       iuser_E = 0, E's calculated with subroutine POT_E
!               = 1, 2, 3, E's calculated with subroutine USER_E_1 (see below)
!
!       iuser_E = 0, Ec, Ev from Ewald summation in subroutine POT_E (same as wmin function)
!               = 1, Ec, Ev from Ewald summation in subroutine USER_E_1
!               = 2, Ec, Ev by direct summation distance = ddmax in USER_E_1
!               = 3, Ec from Ewald summation, Ev by direct summation to ddmax in USER_E_1
!
!          IROT = 0, if no bond to be rotated, bent or turned
!               = 1, if bond(s) to be rotated, bent or turned in RSS & LSQ with intra E 
!                    calculated with PM3 or b3lyp_631g*
!
!       ISYSTEM = 0 for triclinic,monoclinic and orthorhombic
!               = 1 for tetragonal and hexagonal
!               = 2 for cubic
!               = 3 for trigonal
!
!          nmol = # of molecules/independent fragments (usually 1)
!
!          sthl = max sin theta/lambda (usually 0.5)
!
!            ck = Ewald constant (noramlly 0.25, must be > 0)
!
!          dmax = maximum radius for intermolecular interactions around atoms of asymmetric unit
!                 (normally 6 - 10 Angs)
!
!         ddmax = radius for interactions for direct summation if iuser_E = 2 or 3 (50 Angs might work);
!                 all atoms in surrounding molecules whose centers are .le. ddmax of center of parent
!                 molecule are included
!
!       i_cross = # of potential coefficient lines (line # 5 and onward); if 2 then lines # 5 and 6 would
!                 contain potential data for two atomi...atomj intermolecular interactions
!
!  allow_change = (0/1) for (do not/do) adjust l.s. param changes to reflect changes in sum of 1st 
!                 derivatives; if 1, may speed up convergence; probably best to leave at 0
!
!    coef_table = 0/1 for (no/yes) read A, B, C potential coefs from file = coefficient_table
!
!----Line # 2: principal calculation control parameters
!
      READ (7,*) ISWitch , irosen , IUSer_e , IROt , ISYstem , NMOl ,   &
                 STHl , CK , DMAx , DDMax , I_Cross , allow_change, &
                 coef_table
      WRITE (12,99001) ISWitch , irosen , IUSer_e , IROt , ISYstem ,    &
                       NMOl , STHl , CK , DMAx , DDMax , I_Cross ,      &
                       allow_change, coef_table
99001    FORMAT (' PMIN.inp line # 2: iswitch,irosen,iuser_E,irot,isystem,nmol,sthl,ck,dmax,ddmax,',&
                 'i_cross,allow_change,coef_table:'/'                   '6I3,2F6.2,2F5.1,3I3)       
!
!----Lines # 3 & 4 : parameters for LS refinement
!
      READ (7,*) ncycle , NCYc_rss , step1 , step2 , step3 , step4 ,    &
                 frac_change                                                 
      WRITE (12,99002) ncycle , NCYc_rss , step1 , step2 , step3 ,      &
                       step4 , frac_change                                     
99002    FORMAT (' PMIN.inp line # 3: LS cycles, RSS cycles, step1...4, frac_change1...4:'/ &
                 '                    '2I5,4F10.5,4F7.2)
      WRITE (12,99003) are , slim , cmpr      ! defined in data statements
99003    FORMAT (' are, slim, cmpr:',3F10.6/)
      IMOde = IUSer_e
      NROtbond = 0
      IF ( allow_change==0 ) l_allow_change = .FALSE.      ! 0 = false = make l.s. param changes
      DO i = 1 , 20                                        ! 1 = true  = do not make
         INDx(i) = 0                                       !             l. s. param changes
         DO j = 1 , 100
            IATom(i,j) = 0
         ENDDO
      ENDDO
      nvr = 0
      ispl = 0                                             ! from S.P Feb.5,09
      do i = 1,5                                           ! from S.P Feb.5,09
        ispl_1(i) = 0                                      ! from S.P Feb.5,09
      enddo                                                ! from S.P Feb.5,09
!
!----Information on bonds to be modified (conformational refinement)
!      nrotbond must be .le 20
!      E_weight = weight for intramolecular energy term in total E sum
!      NE_type = 1/2/3/4 for intramol E calcd with PM3 from g03/ b3lyp/631g* from g03 / 
!                PM3 from MOPAC2009 code/ PM6 from MOPAC2009 code                          ! 4-28-09 DU
!
Rotbnds:IF ( IROt/=0 ) THEN
           READ (7,*) NROtbond , E_Weight , NE_type  
           IF ( E_Weight<0.00001 ) E_Weight = 1.0                                          ! 9-2-08 DU
           IF ( NE_type==1 ) WRITE (61,1288) NROtbond , E_Weight       
1288             format ('Number of bonds to be modified',i3/  &
                         '  PM3 intramolecular energy evaluation'/ &
                         '  Intramolecular E weight=',f6.2) 
           if (ne_type==2)  WRITE (61,1289) NROtbond , E_Weight                   
1289           format ('Number of bonds to be modified',i3/  &
                       '  b3lyp-631gs*  intramolecular energy evaluation'/ &
                       '  Intramolecular E weight=',f6.2) 
           if (ne_type==3)  WRITE (61,1287) NROtbond , E_Weight                            ! 4-28-09 DU
1287           format ('Number of bonds to be modified',i3/  &                             ! 4-28-09 DU
                       '  PM3 intramolecular energy evaluation from MOPAC code'/ &         ! 4-28-09 DU
                       '  Intramolecular E weight=',f6.2)                                  ! 4-28-09 DU
           if (ne_type==4)  WRITE (61,1286) NROtbond , E_Weight                            ! 4-28-09 DU
1286           format ('Number of bonds to be modified',i3/  &                             ! 4-28-09 DU
                       '  PM6 intramolecular energy evaluation from MOPAC code'/ &         ! 4-28-09 DU
                       '  Intramolecular E weight=',f6.2)                                  ! 4-28-09 DU
           nvr = NROtbond
           DO i = 1 , NROtbond
              READ (7,*) INDx(i) , (IATom(i,ik),ik=1,4)
              WRITE (61,"('index :',4I4)") (IATom(i,ik),ik=1,4)                            ! 3-5-08 Du
           ENDDO
        ENDIF Rotbnds
!----Read information if cross-terms are used
!--- The normal energy coefficients A, B, C are rplaced by special coefficients            ! 10-28-09 
!    A_Cross, B_Cross, C_Cross                                                             ! 10-28-09
!    for O-H... O and N-H... O Hydrogen bonds. These coefficients are read as input lines. ! 10-28-09
!    I_Cross is the no. of H-bonds with special coefficients. ANAme1 and ANAme2 refer to   ! 10-28-09
!    the donor and acceptor atoms respectively. This treatment of H-bonds can be extended  ! 10-18-09
!    to any particular interaction for which the special coefficients A_Coss, B_Cross,     ! 10-31-09
!    C_Cross are read as input.                                                            ! 10-31-09
!
Crosscoef:IF ( I_Cross/=0 ) THEN                                                           ! 5-1-08 DU
             WRITE (61,99004)                                                              ! 5-1-08 DU
99004           FORMAT (' Special A, B and C cross-term energy coefficients...'/)           ! 5-1-08 DU
             write (12,99004)                                                              ! 6/29/09
             DO i = 1 , I_Cross                                                            ! 5-1-08 DU
                READ (7,99055) line                                                        ! 5-1-08 DU
                READ (line,99005) ANAme1(i) , ANAme2(i)                                    ! 5-1-08 DU
99005              FORMAT (2A6)                                                            ! 5-1-08 DU
                READ (line(13:80),*) A_Cross(i) , B_Cross(i) , C_Cross(i), DD_Cross(i)     ! 5-4-10 DU
                WRITE (20,99057) ANAme1(i) , ANAme2(i) , A_Cross(i) ,       &              ! 5-1-08 DU
                                 B_Cross(i) , C_Cross(i), DD_Cross(i)                      ! 5-4-10 DU
                WRITE (61,99057) ANAme1(i) , ANAme2(i) , A_Cross(i) ,       &              ! 5-1-08 DU
                                 B_Cross(i) , C_Cross(i), DD_Cross(i)                      ! 5-4-10 DU 
                WRITE (12,99057) ANAme1(i) , ANAme2(i) , A_Cross(i) ,       &              ! 6/29/09 
                                 B_Cross(i) , C_Cross(i), DD_Cross(i)                      ! 5-4-10 DU
             ENDDO                                                                         ! 5-1-08 DU
          ENDIF Crosscoef 
!----Finished with cross-term input
      nv = 6 + 6*NMOl
      DO i = 1 , MATM
         ATWt(i) = 0.0
         atw(i) = 0.0
         igroup(i) = 1
      ENDDO
      DO i = 1 , nv
         deriv1(i) = 0.0
         DO j = 1 , nv
            deriv2(i,j) = 0.0
         ENDDO
      ENDDO
!
!----Replace the A,B,C coeffecients from a table
!
C_table:IF (coef_table .ne. 0) THEN  ! read pot coefs from coefficient_table 
           open (unit=65, file = 'coefficient_table', status='old')
           read (65,*)      ! first line is header, skip it
           do j=1,matm
              read (65,8811,end=8822) i, A1_tab(i), B1_tab(i), C1_tab(i)      ! 1-5-09 DU
8811             format (i3,3f15.5)
           enddo
8822       write (12,8833) j-1
8833          format (i4,' potential coefficients read from coefficient_table')
           close (unit = 65)
        ENDIF C_table
!----End of section to read pot coefs from table
!
!----Section of read cell params, coords, etc
!
      go to (8787, 8788, 8788), iswitch
!
!----iswitch = 2, read all info from PMIN.inp...wmin.input NOT used
!     
8788     READ (7,99959) REFcode, SYMcd     ! refcode
99959    FORMAT(A8,2X,A2)                  ! 1-29-09 DU
         READ (7,*) AH                     ! cell params
         DO j = 1 , 6                      ! 1-5-09 DU save initial cell
            cell_ini(j) = AH(j)
         ENDDO
         READ (7,*) natom, NSYm, wt_mol    ! # atoms, # of sym ops, molec wt
         N = natom                         ! 1-5-09 DU
         DO j = 1 , NSYm                   ! read nsym sym op lines
            READ (7,*) (SYM(j,i),i=1,12)
         ENDDO
         zz = NSYm
!----Read the atom parameters and (maybe) potential coefficients (std wmin type)
         DO i = 1,natom   ! natom = # atoms 
            if (coef_table .eq. 0) then
               READ (7,99006) ATOm(i) , IWCode(i) , X(i) , Y(i) , Z(i), &  ! everything on
                           Q(i) , A1(i) , B1(i) , C1(i) , ATWt(i)          !  this line
99006             FORMAT (5x,a6,i5,8F10.5)
               WRITE(61,99006) ATOm(i) , IWCode(i) , X(i) , Y(i) , Z(i), & ! 1-5-09 DU
                           Q(i) , A1(i) , B1(i) , C1(i) , ATWt(i)          ! 1-5-09 DU 
            else
              if (iswitch == 2) then                                       ! 1-8-09 DU
               read  (7,8793) ATOm(i) , IWCode(i) , X(i) , Y(i) , Z(i), &  ! no coefs on
                             Q(i), ATWt(i)                                 !  this line
               write (61,8793) ATOm(i) , IWCode(i) , X(i) , Y(i) , Z(i), & ! 1-5-09 DU 
                             Q(i), ATWt(i)
8793              format (5x,a6,i5,5f10.5)                                 ! no coefs on the line
              end if   
               A1(i) = A1_tab(IWCode(i))                                   ! 1-5-09 DU
               B1(i) = B1_tab(IWCode(i))                                   ! 1-5-09 DU
               C1(i) = C1_tab(IWCode(i))                                   ! 1-5-09 DU
            endif
         ENDDO
!----Read refinement codes (0/1 for no/yes); to refine 6 cell params + orientation &
!     translation params for 1 molecule would be 111111111111.  For orthorhombic
!     unit cell and 2 independent molecules, 111000111111111111.
         READ (7,99060) (IDX(i),i=1,6+nmol*6) ! format = 36I1
         read (7,99011) (igroup(i),i=1,natom) ! indicates to which group(s) atoms belong; 
         go to 8899
!
8787   OPEN (UNIT=9,FILE='wmin.input',STATUS='OLD')  ! iswitch = 1...usaal wmin.input file  
!
         DO i = 1 , 4          ! skip 4 lines at top
            READ (9,*)
         ENDDO
         GOTO 50
 50      READ (9,99059) SYMcd , REFcode             ! 6/11/07
         WRITE(99,"(A8)") REFcode                   ! 1-5-09 DU
         READ (9,99008) natom , NSYm
99008    FORMAT (3x,2I3)
         N = natom
         READ (9,*)
         READ (9,*) WT_mol
         DO i = 1 , 2
            READ (9,*)
         ENDDO
         READ (9,*) AH            ! a, b, c, alpha, beta, gamma
         WRITE(99,"(6F9.4)") AH                     ! 1-5-09 DU
         WRITE(99,"(2I3,F8.3)") natom, NSYm, WT_mol ! 1-5-09 DU
         DO j = 1 , 6                                            ! 9-16-08 DU save initial cell
            cell_ini(j) = AH(j)
         ENDDO                                                   ! 9-16-08 DU
         READ (9,*)
         DO j = 1 , NSYm
            READ (9,*) (SYM(j,i),i=1,12)
            WRITE(99,"(3(10X,F5.2,3F4.0))") (SYM(j,i),i=1,12)    ! 1-5-09 DU
         ENDDO
         zz = NSYm
!
!----Read the atom parameters and wmin coefficients
!
         DO i = 1 , natom
            if (coef_table .eq. 0) then
               READ (9,99009) ATOm(i) , IWCode(i) , Q(i) , A1(i) ,    &
                              B1(i) , C1(i) , ATWt(i)
99009          FORMAT (a6,i3,f9.6,4F9.4,i3)
            else
               read (9,99099) ATOm(i) , IWCode(i) , Q(i), atwt(i)
99099          FORMAT (a6,i3,f9.6,27X,F9.4)                                 ! 1-5-09 DU
               A1(i) = A1_tab(IWCode(i))                                    ! 1-5-09 DU
               B1(i) = B1_tab(IWCode(i))                                    ! 1-5-09 DU
               C1(i) = C1_tab(IWCode(i))                                    ! 1-5-09 DU
            endif
            atw(i) = ATWt(i)
         enddo
         DO i=1,natom    ! fractl coords from latter part of wmin.input file
            ATWt(i) = atw(i)
            READ (9,99010) ATOm(i) , X(i) , Y(i) , Z(i)
99010       FORMAT (a6,21x,3F9.5)
         ENDDO
!
         DO i = 1 , NMOl
            READ (9,*)
         ENDDO
         READ (9,99011) (igroup(i),i=1,natom)            ! 1-8-09 DU
99011       FORMAT (24I3)                                ! + nmol extra atoms (XTRA)
         DO i = 1 , NMOl + 1
            READ (9,*)
         ENDDO
         READ (9,99060) (IDX(i),i=1,nv)
!
8899     IF ( IROt/=0 ) THEN
            DO i = 1 , NROtbond
               DO ik = 1 , 4
                  k11(ik) = IATom(i,ik)
               ENDDO
               CALL ATOM_LINKS(N,AH,X,Y,Z,k11,k55,listatm)
               K5(i) = k55
               DO ik = 1 , k55
                  LINkatom(i,ik) = listatm(ik)
               ENDDO
               WRITE (61,*) ' i,k11,k55,linkatom ' , i , k11 , K5(i) ,  &
     &                      (LINkatom(i,ik),ik=1,k55)       
               WRITE (61,"(20X,I4,4A4,I4,20A4)") i ,         &
     &                (ATOm(k11(ik)),ik=1,4) , K5(i) ,       &
     &                (ATOm(LINkatom(i,ik)),ik=1,k55)        
            ENDDO
         ENDIF
!
         IF ( NMOl>1 ) THEN            ! nmol = 1, only one fragment
!
!         set atoms to proper molecules
!
            DO i = 1 , NMOl            ! more than one fragment
               NATm(i) = 0
            ENDDO
            DO i = 1 , natom
               NATm(igroup(i)) = NATm(igroup(i)) + 1
               k1atom(igroup(i),NATm(igroup(i))) = i
               xg(igroup(i),NATm(igroup(i))) = X(i)
               yg(igroup(i),NATm(igroup(i))) = Y(i)
               zg(igroup(i),NATm(igroup(i))) = Z(i)
               qg(igroup(i),NATm(igroup(i))) = Q(i)
               ag(igroup(i),NATm(igroup(i))) = A1(i)
               bg(igroup(i),NATm(igroup(i))) = B1(i)
               cg(igroup(i),NATm(igroup(i))) = C1(i)
               wg(igroup(i),NATm(igroup(i))) = ATWt(i)
               IWCodeg(igroup(i),NATm(igroup(i))) = IWCode(i)     ! 1-14-09 DU
               atomg(igroup(i),NATm(igroup(i))) = ATOm(i)
            ENDDO
!
            istart = 0
            DO imol = 1 , NMOl
               DO i = 1 , NATm(imol)
                  katom(istart+i) = k1atom(imol,i)
                  X(istart+i) = xg(imol,i)
                  Y(istart+i) = yg(imol,i)
                  Z(istart+i) = zg(imol,i)
                  Q(istart+i) = qg(imol,i)
                  A1(istart+i) = ag(imol,i)
                  B1(istart+i) = bg(imol,i)
                  C1(istart+i) = cg(imol,i)
                  ATWt(istart+i) = wg(imol,i)
                  ATOm(istart+i) = atomg(imol,i)
                  IWCode(istart+i) = IWCodeg(imol,i)             ! 1-14-09 DU  
               ENDDO
               istart = istart + NATm(imol)
            ENDDO
!
            IF ( IROt/=0 ) THEN
               DO i = 1 , NROtbond
                  DO ik = 1 , 4
                     DO j = 1 , natom
                        IF ( IATom(i,ik)==katom(j) ) THEN
                           jatom(i,ik) = j
                           IATom(i,ik) = jatom(i,ik)
                           EXIT
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
!
            IF ( IROt/=0 ) THEN
               DO i = 1 , NROtbond
                  DO ik = 1 , 4
                     k11(ik) = IATom(i,ik)
                  ENDDO
                  CALL ATOM_LINKS(N,AH,X,Y,Z,k11,k55,listatm)
                  K5(i) = k55
                  DO ik = 1 , k55
                     LINkatom(i,ik) = listatm(ik)
                  ENDDO
               ENDDO
            ENDIF
         ELSE
            NATm(1) = natom
         ENDIF
!
         istop = 0
         IF ( istop==1 ) GOTO 99999
         N = natom
!      write(8,516)idx                                           ! removed 7/11/08 HLA
!
!------ END LOOP TO READ ATOM COORDINATES
!
      WT_mol = 0.0
      DO i = 1 , natom
         WT_mol = WT_mol + ATWt(i)
      ENDDO
      N = natom
!
      IF ( ISYstem==2 ) THEN
         CALL UNIQUE_SYM(SYM,X,Y,Z,nsym1)
         nsym2 = NSYm
         NSYm = 0
         DO j = 1 , nsym2
            IF ( nsym1(j)==1 ) THEN
               NSYm = NSYm + 1
               DO i = 1 , 12
                  SYM(NSYm,i) = SYM(j,i)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      CALL MULTIPLICITY(SYM,ATWt,X,Y,Z,mcity,wt_mol_asym)
      WT_mol = wt_mol_asym
!
!---- END LOOP TO READ ATOM COORDINATES
!
      sum_deriv1 = 0.0
!
      IRSs_call = 0                          ! 4-30-08, HLA, added
      ICYcle = 0
      DO i = 1 , nv      !modified
         corr_sum(i) = 0.0
      ENDDO
      IF ( IROt>0 ) THEN                      ! 9-23-08 from S.P
         CALL ORTHO_COD(N,AH,X,Y,Z,xo,yo,zo)  ! 9-23-08 from S.P
         IF ( NE_type==1 ) THEN               ! 9-30-08 DU
            CALL ENERGY_PM3(N,xo,yo,zo,egymol)
         ELSE IF ( NE_type==2 ) THEN
            CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)
         ELSE 
           IF ( NE_type==3 .or. ne_type==4 ) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)   ! PM3 or PM6 energy from MOPAC code 4-28-09 DU
         ENDIF
         egymol_ori = egymol                  ! 9-23-08 from S.P
         EGYmol0 = egymol                     ! 9-23-08 from S.P
         WRITE (13,7124) egymol_ori
7124        format (4X,'E of initial conformation =',E16.10)
      ENDIF
!
 100  ICYcle = ICYcle + 1
      irun = 0
!
      a = AH(1)
      b = AH(2)
      c = AH(3)
!
      IF ( AH(4)>10.0 ) THEN      ! angles or cosines?
         AH(4) = COSD(AH(4))      ! changed to COSD from COS   6/11/07
         AH(5) = COSD(AH(5))      ! ah(4-6) are cosines        6/11/07
         AH(6) = COSD(AH(6))      ! 6/11/07
      ENDIF
      c11 = AH(4)                 ! c11, c22, c33 are cosines
      c22 = AH(5)
      c33 = AH(6)
      alpha = ACOSD(c11)          ! alpha, beta, gamma
      beta = ACOSD(c22)           ! are angles in degrees
      gamma = ACOSD(c33)
      CELl1(1) = a                ! lengths
      CELl1(2) = b
      CELl1(3) = c
      CELl1(4) = c11              ! cell1 are cosines
      CELl1(5) = c22
      CELl1(6) = c33
      IF ( ICYcle==1 ) THEN
         DO ik = 1 , 6
            cell_1st(ik) = CELl1(ik)
         ENDDO
      ENDIF
      CELl(1) = a
      CELl(2) = b
      CELl(3) = c
      CELl(4) = alpha             ! cell(4-6) are angles
      CELl(5) = beta              ! in degrees
      CELl(6) = gamma
!----Parameters to be adjusted are p(i)'s
      DO i = 1 , 6
         p(i) = CELl1(i)
      ENDDO
      DO i = 7 , nv               ! modified
         p(i) = 0.0
      ENDDO
!----Initial cell params, etc
      IF ( ICYcle==1 ) THEN                  
          WRITE (12,99015) CELl1 , CELl       
99015        FORMAT (/' Initial cell parameters, lengths, cosines and angles...' & 
                    //3f10.4,3F12.6/3f10.4,3F12.6/)     
      do i=1,6
         cell_input(i) = cell1(i)     ! cell_input = inital values of cell params
      enddo                           ! lengths & cosines
!----Initial list of atoms, etc
      WRITE (12,99012)
99012 FORMAT (                                                          &
     &' Atoms, fracl coords, charges, potential coefs, at wts, fragment &
     &#, Z''...'//                                                      &
     &'  #...ID...code......x.........y.........z........q....',        &
     &'....A........B......C.....atwt..mol#'''/)
      istart = 1
      DO imol = 1 , NMOl
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            WRITE (12,99013) i , ATOm(i) , IWCode(i) , X(i) , Y(i) ,    &
     &                       Z(i) , Q(i) , A1(i) , B1(i) , C1(i) ,      &
     &                       ATWt(i) , imol
99013       FORMAT (i3,2x,a6,i3,3F10.5,3F9.4,f6.2,f9.4,2I4)
            WRITE (99,99006) ATOm(i) , IWCode(i) , X(i) , Y(i) ,        &    ! 1-5-08 DU
     &                       Z(i) , Q(i) , A1(i) , B1(i) , C1(i) ,      &
     &                       ATWt(i) 
         ENDDO
         istart = istart + NATm(imol)
      ENDDO
      WRITE(99,99060)(IDX(i),i=1,nv)                  ! 1-5-09 DU  
      WRITE(99,99011)(igroup(i),i=1,natom)            ! 1-5-09 DU
!----Calc centers of independent nmol fragments!
      istart = 1
      DO imol = 1 , NMOl      
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = X(i) + sumx
            sumy = Y(i) + sumy
            sumz = Z(i) + sumz
         ENDDO
         xc(imol) = REAL(sumx/NATm(imol))
         yc(imol) = REAL(sumy/NATm(imol))
         zc(imol) = REAL(sumz/NATm(imol))
         WRITE (12,99014) xc(imol) , yc(imol) , zc(imol) , imol
99014       FORMAT ('     XTRA',5x,3F10.5,43x,i3)
         dxy = abs(abs(xc(imol))-abs(yc(imol)))             ! from S.P. Feb.5,09
         dyz = abs(abs(yc(imol))-abs(zc(imol)))             ! from S.P. Feb.5,09
         if(abs(dxy) < 0.0001 .and. abs(dyz) < 0.0001)then  ! from S.P. Feb.5,09
         ispl_1(imol) = 1                                   ! from S.P. Feb.5,09
         ispl = 1                                           ! from S.P. Feb.5,09
         endif                                              ! from S.P. Feb.5,09
         istart = istart + NATm(imol)
      ENDDO
      ENDIF       
!
      sumxx = 0.0                                           ! from S.P. Feb.5,09
      sumyy = 0.0                                           ! from S.P. Feb.5,09
      sumzz = 0.0                                           ! from S.P. Feb.5,09
      DO i = 1 , N
         X_Ori(i) = X(i)
         Y_Ori(i) = Y(i)
         Z_Ori(i) = Z(i)
         sumxx = sumxx + x(i)                               ! from S.P. Feb.5,09
         sumyy = sumyy + y(i)                               ! from S.P. Feb.5,09
         sumzz = sumzz + z(i)                               ! from S.P. Feb.5,09
      ENDDO
         xc1 = sumxx/n                                      ! from S.P. Feb.5,09
         yc1 = sumyy/n                                      ! from S.P. Feb.5,09
         zc1 = sumzz/n                                      ! from S.P. Feb.5,09
       write(8,6756)icycle,xc1,yc1,zc1                      ! from S.P. Feb.5,09
6756   format(/,' icycle, centriod of molecule',i5,3f12.6,/)! from S.P. Feb.5,09
!
      sah = 0.5*(alpha+beta+gamma)
      sins = SIN(ANG_TO_RADIAN*sah)                         ! 5/13/08, HLA, altered
      sin1 = SIN(ANG_TO_RADIAN*(sah-alpha))                 ! 5/13/08, HLA, altered
      sin2 = SIN(ANG_TO_RADIAN*(sah-beta))                  ! 5/13/08, HLA, altered
      sin3 = SIN(ANG_TO_RADIAN*(sah-gamma))                 ! 5/13/08, HLA, altered
      vol = 2.0*AH(1)*AH(2)*AH(3)*DSQRT(sins*sin1*sin2*sin3)
      DENsity = 1.6605*zz*WT_mol/vol
!
      IF ( ICYcle==1 .AND. IRSs_call==0 ) THEN
         DENs_1st = DENsity
         VOL_1st = vol
         dens_ini = DENsity                                 ! 9-16-08 DU
         vol_ini = vol                                      ! 9-16-08 DU
         write(61,4237) icycle,zz,wt_mol,vol,density                       
 4237    format(/,5x,'icycle, z_volue, wt_mol;  initial vol & density ',i5,2f10.2,1X,2f10.5,/) 
      ENDIF
!
      IF ( ICYcle<21 .OR. MOD(ICYcle,10)==1 ) NENtry = 1
      ILS = 0
      IF ( IUSer_e>0 ) THEN
         CALL USER_E_1(CELl,X,Y,Z,ec,ev,er,e_nmol,E)        ! 9-23-08 from S.P
      ELSE
         CALL POT_E(CELl,X,Y,Z,ec,ev,er,e_nmol,E)           ! 9-23-08 from S.P
      ENDIF
!
     IF ( ICYcle==1 .AND. IRSs_call==0 ) THEN               ! 12-16-08 
      ec_ini = ec                                           ! 12-16-08 DU
      ev_ini = ev                                           ! 12-16-08 DU
      er_ini = er                                           ! 12-16-08 DU
      e_ini  = e                                            ! 12-16-08 DU
      egy_ini(1) = ec_ini                                   ! 12-16-08 DU
      egy_ini(2) = ev_ini                                   ! 12-16-08 DU
      egy_ini(3) = er_ini                                   ! 12-16-08 DU
      egy_ini(4) = e_ini                                    ! 12-16-08 DU
      e_nmol_ini = e_nmol                                   ! 12-16-08 DU
      write(8,"(5F16.6)") ec_ini,ev_ini,er_ini,e_ini,ec_ini ! temp
      write(8,"(2F16.4)") vol_ini, dens_ini                 ! temp
     END IF                                                 ! 12-16-08 DU
      w01 = E
      WRITE (8,*) ' E without egymol ' , E                  ! 9-23-08 from S.P
!
      IF ( IROt==2 ) THEN                                   ! 9-23-08 DU
         CALL ORTHO_COD(N,CELl,X,Y,Z,xo,yo,zo)              ! 9-23-08 from S.P
         IF ( NE_type==1 ) THEN                             ! 9-30-08 DU
            CALL ENERGY_PM3(N,xo,yo,zo,egymol)              ! 9-23-08 from S.P
         ELSE IF ( NE_type==2 ) THEN
           CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)             ! 9-30-08 DU
         ELSE 
           IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)            ! 4-28-09 DU
         ENDIF                                                               ! 9-30-08 DU
         DEL_egymol = egymol - egymol_ori                                    ! 9-23-08 from S.P
         w01 = E + DEL_egymol                                                ! 9-23-08 from S.P
      ENDIF                                                                  ! 9-23-08 from S.P
!
      WRITE (8,99018) ICYcle , NENtry , N11 , E , DEL_egymol , w01           ! 9-23-08 from S.P
99018 FORMAT ('  icycle,nentry,n11,E,del_egymol,W01 ',3I6,3F16.8)            ! 9-23-08 from S.P
!
!       write(8,1212) icycle,nentry,n11,W01                                  ! removed 7/11/08 HLA
!1212    format('  icycle,nentry,n11,W01 ',3I6,F16.8)                        ! removed 7/11/08 HLA
!
      IF ( irosen==0 .OR. ICYcle>1 ) THEN                                    ! 9-30-08 DU
         IF ( IROt==1 ) IROt = 2                                             ! to refine conformation only in LSQ 9-30-08
         GOTO 200
      ENDIF
!
      IF ( ICYcle==1 .AND. irosen==1 ) THEN
         IF ( IROt/=0 ) IROt = 1                                             ! 10-10-08 DU
         CALL ROSENBROCK(fx,fy,fz)                                           ! Rosenbrock step search
         GOTO 99999
      ENDIF
!
      IF ( ICYcle==1 .AND. irosen==2 ) THEN
         IF ( IROt/=0 ) IROt = 1                                             ! 10-10-08 DU
         CALL ROSENBROCK(fx,fy,fz)
         IF ( IROt==1 ) THEN
            IROt = 2                                       ! to continue refinement conformation in LSQ after RSS 9-30-08 DU
            del_egymol_rss = EGYmol0 - egymol_ori          ! 9-30-08 DU new conf. diff E after RSS
            egymol_ori = EGYmol0                           ! 9-30-08 DU new conf. E after RSS as initial conf. E for LS
         ENDIF
      ENDIF
 
      NENtry = 1
      WRITE (13,99019)
99019 FORMAT (/'****Start least squares refinement...')
!
      a = CELl(1)
      b = CELl(2)
      c = CELl(3)
      alpha = CELl(4)
      beta = CELl(5)
      gamma = CELl(6)
!
      sah = 0.5*(alpha+beta+gamma)
      sins = SIN(ANG_TO_RADIAN*sah)                                    ! 5/13/08, HLA, added
      sin1 = SIN(ANG_TO_RADIAN*(sah-alpha))                            ! 5/13/08, HLA, added
      sin2 = SIN(ANG_TO_RADIAN*(sah-beta))                             ! 5/13/08, HLA, added
      sin3 = SIN(ANG_TO_RADIAN*(sah-gamma))                            ! 5/13/08, HLA, added
      vol = 2.0*AH(1)*AH(2)*AH(3)*DSQRT(sins*sin1*sin2*sin3)
      DENsity = 1.6605*zz*WT_mol/vol
!
!     IF ( ICYcle==1 .AND. IRSs_call==0 ) THEN
      IF ( ICYcle==1) THEN                                             ! 1-22-09 DU
         DENs_1st = DENsity
         VOL_1st = vol
!       write(8,4237) icycle -1,zz,wt_mol,vol,density                  ! removed 7/11/08 HLA
      ENDIF
!        write(8,5223)cell                                             ! removed 7/11/08 HLA
!5223     format(5x,'cell after RSS ',6f14.8)                          ! removed 7/11/08 HLA
!
      DO i = 1 , 3
         cell_1st(i) = CELl(i)
      ENDDO
      DO i = 4 , 6
         cell_1st(i) = DCOS(CELl(i)*ANG_TO_RADIAN)                              ! 5/12/08, HLA, modification
!        write(8,*)' After RSS cell(i), cell_1st(4,5,6) ',cell(i),cell_1st(i)   ! removed 7/11/08 HLA
      ENDDO
      DO i = 1 , N
!        write(8,513)atom(i),x(i),y(i),z(i)                                     ! removed 7/11/08 HLA
         X_Ori(i) = X(i)
         Y_Ori(i) = Y(i)
         Z_Ori(i) = Z(i)
      ENDDO
!
      istart = 1
      DO imol = 1 , NMOl
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = X(i) + sumx
            sumy = Y(i) + sumy
            sumz = Z(i) + sumz
         ENDDO
! Centriod of atoms for atoms of the group imol xc,yc,zc
         xc(imol) = REAL(sumx/NATm(imol))
         yc(imol) = REAL(sumy/NATm(imol))
         zc(imol) = REAL(sumz/NATm(imol))
!
!        write(8,8163)imol,xc(imol),yc(imol),zc(imol)               ! removed 7/11/08 HLA
!8163     format('XTR',i1,5x,18x,3f12.8)                            ! removed 7/11/08 HLA
         istart = istart + NATm(imol)
      ENDDO
!
      WRITE (13,*) ' '
      WRITE (13,7231) icycle-1
7231     format ('  After RSS, least squares cycle',i3,'...')
      WRITE (13,*) ' '
!
      IF ( IROt>0 ) THEN                                            ! 9-23-08 from S.P
         CALL ORTHO_COD(N,CELl,X,Y,Z,xo,yo,zo)                      ! 9-23-08 from S.P
         IF ( NE_type==1 ) THEN                                     ! 9-30-08 DU
            CALL ENERGY_PM3(N,xo,yo,zo,egymol)                      ! 9-23-08 from S.P
         ELSE IF ( NE_type==2 ) THEN
            CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)                    ! 9-30-08 DU
         ELSE
            IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)  ! 4-28-09 DU
         ENDIF                                                      ! 9-30-08 DU
!         egymol0 = egymol                                          ! 9-23-08 from S.P
!         del_egymol = 0.0                                          ! 9-23-08 from S.P
         DEL_egymol = egymol - egymol_ori                           ! 9-23-08 from S.P
         WRITE (8,*) 'del_egymol after RSS egymol,egymol_ori ' , &
     &               egymol , egymol_ori , DEL_egymol               ! 9-23-08 from S.P
      ENDIF                                                         ! 9-23-08 from S.P
!
      IF ( IMOde>0 ) THEN
         CALL USER_E_1(CELl,X,Y,Z,ec,ev,er,e_nmol,E)                ! 9-23-08 from S.P
      ELSE
         CALL POT_E(CELl,X,Y,Z,ec,ev,er,e_nmol,E)                   ! 9-23-08 from S.P
      ENDIF
!
      w01 = E
      WRITE (13,99020) ICYcle - 1 , NENtry , N11 , ec , ev , er , &
     &                 e_nmol , w01                                 ! 9-23-08 from S.P
99020 FORMAT (2x,'Icycle, nentry, n11, Ec, Ev, Er, E after RSS...'/ &
              I5,i3,i5,5F13.8)
!
 200  ILS = 1
      nv = 6*NMOl + 6                                               ! exclude nvr
!
      IF ( IROt==2 ) nv = 6*NMOl + 6 + nvr                          ! 9-23-08 from S.P
!
!        write(8,4332)icycle,nentry,W01                             ! removed 7/11/08 HLA
!4332     format(/5x,'CYCLE,nentry,W01 ',2i5,f20.12,/)              ! removed 7/11/08 HLA
!
!        set parameters limits
!
      DO i = 1 , 3                                                  ! unit cell lengths; if length < 10 Angs
         IF ( p(i)>10. ) THEN
            pcmin(i) = p(i) - 1.0                                   ! make limit just 10%
            pcmax(i) = p(i) + 1.0
         ELSE
            pcmin(i) = 0.9*p(i)
            pcmax(i) = 1.1*p(i)
         ENDIF
      ENDDO
!
      DO i = 4 , 6                                                 ! unit cell cosine limits
         pcmin(i) = p(i) - 0.1                                     ! changed from 0.2 to 0.1, 7/4/08, HLA
         pcmax(i) = p(i) + 0.1               
         IF ( pcmin(i)<-0.8 ) pcmin(i) = -0.82
         IF ( pcmax(i)>0.8 ) pcmax(i) = 0.82
      ENDDO
!
      DO imol = 1 , NMOl                                          ! model orientation cosines
         DO i = 6*imol + 1 , 6*imol + 3
            pcmin(i) = -0.1                                       ! changed from .2 to .1, 7/4/08, HLA
            pcmax(i) = +0.1                    
         ENDDO
         DO i = 6*imol + 4 , 6*imol + 6
            pcmin(i) = -0.25                                      ! reduced translation from 0.5 to 0.25, 7/4/08, HLA
            pcmax(i) = 0.25                                             
         ENDDO
      ENDDO
      IF ( IROt>0 ) THEN                                          ! 9-23-08 from S.P
!        nvr = nrotbond
         DO j = nv - nvr + 1 , nv
            pcmin(j) = -PI/2.0   ! min & max rotation angles in radian for rotbond
            pcmax(j) = PI/2.0
            IDX(j) = 1                                            ! 9-23-08 from S.P
         ENDDO
      ENDIF
!
      IF ( ICYcle==1 ) THEN
         ecvr_1st(1) = ec
         ecvr_1st(2) = ev
         ecvr_1st(3) = er
         ecvr_1st(4) = E
         emol_1st = 0.0                                           ! 9-23-08 from S.P
         IF ( irosen==2 .AND. IROt==1 ) THEN                      ! 9-23-08 from S.P
            emol_1st = DEL_egymol                                 ! 9-23-08 from S.P
            ecvr_1st(4) = ecvr_1st(4) + DEL_egymol                ! 9-23-08 from S.P
         ENDIF                                                    ! 9-23-08 from S.P
      ENDIF
      e_prev(ICYcle) = w01
      emol_prev(ICYcle) = DEL_egymol                              ! 9-23-08 from S.P
      IF ( ICYcle==1 ) THEN
!        E_diff = 0.0
         emol_diff = 0.0                                          ! 9-23-08 from S.P
         e_diff = 0.111111111                                     ! test, 6/30/08
      ELSE
         e_diff = e_prev(ICYcle) - e_prev(ICYcle-1)
         emol_diff = emol_prev(ICYcle) - emol_prev(ICYcle-1)      ! 9-23-08 from S.P
         cycle2 = ICYcle
      ENDIF
!
      IF ( IRSs_call==0 .AND. ICYcle==1 ) WRITE (12,99021)        ! 5-25-08 DU
99021 FORMAT (/' Start least squares refinement...')
!     if(E_prev(icycle) > 0.0 ) then
      IF ( e_prev(ICYcle)>0.0 .OR. (e_diff>1.0 .AND. irosen==0) ) THEN      ! 9-30-08 DU
         ICYcle = ICYcle - 1
         NENtry = 1
         ILS = 0
         DO i = 1 , 6
            CELl(i) = cell_prev(i)
            AH(i) = cell_prev(i)
         ENDDO
         DO i = 1 , N
            X(i) = x_prev(i)
            Y(i) = y_prev(i)
            Z(i) = z_prev(i)
         ENDDO
         ncycle = ICYcle
         IRSs_call = IRSs_call + 1
         IF ( IRSs_call>=2 ) e_diff_ls_test = 1.0D-2                       ! 1-20-09 DU
         IF ( e_prev(ICYcle)>0.0 ) WRITE (62,99022) ICYcle , &
     &        e_prev(ICYcle)                                               ! 9-30-08 DU
99022    FORMAT (' ----> Warning...E positive in LS; Cycle # & E : ',I3,&
     &           F16.6)                                                    ! 9-30-08 DU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF ( e_diff>1.0 .AND. irosen==0 ) THEN         ! 9-30-08 DU
!         write (62,8008) E_diff, cycle2, E_prev(cycle2)! 9-30-08 DU
!         print 8008, E_diff, cycle2, E_prev(cycle2)    ! 9-30-08 DU
            irosen = 2                                  ! 9-30-08 DU  reset "irosen" toRestart PMIN with RSS plus LS
            ICYcle = 0                                  ! 9-30-08 DU
            DO i = 1 , nv                               ! 10-10-08 DU
               corr(i) = 0.0                            ! 10-10-08 DU
               corr_sum(i) = 0.0                        ! 10-10-08 DU
            ENDDO                                       ! 10-10-08 DU
            IRSs_call = 0                               ! 9-3008 DU
            IF ( IROt>0 ) THEN                          ! 9-30-08 DU
               CALL ORTHO_COD(N,AH,X,Y,Z,xo,yo,zo)      ! 9-30-08 DU
               IF ( NE_type==1 ) THEN                   ! 9-30-08 DU
                  CALL ENERGY_PM3(N,xo,yo,zo,egymol)    ! 9-30-08 DU
               ELSE IF ( NE_type==2 ) THEN
                  CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)  ! 9-30-08 DU
               ELSE
                  IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)   ! 4-28-09
               ENDIF                                    ! 9-30-08 DU
               egymol_ori = egymol                      ! 9-30-08 DU
               EGYmol0 = egymol                         ! 9-30-08 DU
!          print "('new conf. E ',E16.10)", egymol      ! 9-30-08 DU
               WRITE (62,99061) e_diff , cycle2 , e_prev(cycle2) , &
     &                          egymol
                                                        ! 9-30-08 DU
               WRITE (13,99061) e_diff , cycle2 , e_prev(cycle2) , &
     &                          egymol
                                                        ! 9-30-08 DU
               PRINT 99061 , e_diff , cycle2 , e_prev(cycle2) , egymol
                                                        ! 9-30-08 DU
            ENDIF                                       ! 9-30-08 DU
         ENDIF                                          ! 9-30-08 DU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         GOTO 100
      ENDIF
!     if( E_diff > 1.0 ) then                           ! 4-30-08, HLA, changed from 0.0 to 1.0
!     IF ( e_diff>0.0 ) then
!      print *, '(failure in Least sqr, Job terminated)'
!      GO TO 99999                                      ! END job   10-26-10
!     END IF                                         
      IF ( e_diff>1.0 .AND. irosen==2 ) THEN            ! 9-30-08 DU
         NENtry = 1
         ILS = 0
         DO i = 1 , nv                                  ! 10-10-08 DU
            corr(i) = 0.0                               ! 10-10-08 DU
            corr_sum(i) = 0.0                           ! 10-10-08 DU
         ENDDO
         IF ( IROt/=0 ) IROt = 1                        ! 10-7-08 DU for test
!          WRITE(8,8001)                                ! removed 7/11/08 HLA
!8001          FORMAT('**E_difference > 1 kcal/mol...continue new RSS and LS')   ! removed 7/11/08 HLA
         WRITE (61,99062)                                                        ! 5-25-08 DU
         WRITE (12,99023) e_diff , icycle-1 , e_prev(cycle2)                     ! 1-22-09
         WRITE (61,99023) e_diff , icycle-1 , e_prev(cycle2)
         WRITE (62,99023) e_diff , icycle-1 , e_prev(cycle2)
99023    FORMAT (/'--> Warning...E increase > 1 kcal/mol,', &
     &           '...E diff =',e11.4,', cycle',i3,', current E =',f12.7,&
     &           '...repeat RSS calcns')
         WRITE (61,99062)                                                        ! 5-25-08 DU
         IRSs_call = IRSs_call + 1                                               ! 7-7-08 DU
         IF ( IRSs_call>=2 ) e_diff_ls_test = 1.0D-2                             ! 1-20-09 DU
         IF ( IRSs_call>=5 ) then
         DO i = 1 , 6
            CELl(i) = cell_prev(i)
            AH(i) = cell_prev(i)
         ENDDO
         DO i = 1 , N
            X(i) = x_prev(i)
            Y(i) = y_prev(i)
            Z(i) = z_prev(i)
         ENDDO
         ncycle = ICYcle
         GO TO 101                                                               ! 2-3-09
         END IF
         ICYcle =0
         GOTO 100
      ENDIF
!
!       print "('ICLCLE, E_diff', I6,D16.8)", icycle-1, E_diff                   ! 10-7-08 DU temp
101   DO i = 1 , 6
         cell_prev(i) = CELl(i)
      ENDDO
      DO i = 1 , N
         x_prev(i) = X(i)
         y_prev(i) = Y(i)
         z_prev(i) = Z(i)
      ENDDO
!
!     write(8,421)icycle-1, nentry,a,b,c,alpha,beta,gamma,vol,density, &         ! removed 7/11/08 HLA
!                 Ec,Ev,Er,E_nmol,E,E_diff,sum_deriv1                            ! 9-23-08 from S.P
      IF ( IROt==0 ) THEN                                                        ! 9-30-08 DU
         WRITE (18,99064) ICYcle - 1 , NENtry , a , b , c , alpha ,     &
     &                    beta , gamma , vol , DENsity , ec , ev , er , &
     &                    e_nmol , E , e_diff , sum_deriv1                       ! 9-23-08 from S.P
         WRITE (13,99064) ICYcle - 1 , NENtry , a , b , c , alpha ,     &
     &                    beta , gamma , vol , DENsity , ec , ev , er , &
     &                    e_nmol , E , e_diff , sum_deriv1                       ! 9-23-08 from S.P
      ELSE
         WRITE (18,99065) ICYcle - 1 , NENtry , a , b , c , alpha ,     &
     &                    beta , gamma , vol , DENsity , ec , ev , er , &
     &                    DEL_egymol , e_nmol , E , e_diff , sum_deriv1          ! 9-30-08 DU
!        WRITE(18,423)  Ec,Ev,Er,Ec+Ev+Er,del_egymol,E,E_diff                    ! 9-30-08 DU
         WRITE (13,99065) ICYcle - 1 , NENtry , a , b , c , alpha ,     &
     &                    beta , gamma , vol , DENsity , ec , ev , er , &
     &                    DEL_egymol , e_nmol , E , e_diff , sum_deriv1
      ENDIF                                                                      ! 9-30-08 DU
99024 FORMAT (/8X,'WC',12X,'WV',12X,'WR',12X,'W ',8X,'Wconf',12X,'WT',  &
     &        8x,'E_diff'/f12.6,5F14.6,D16.8)                                    ! 9-30-08 DU
      DO i=1,nv
         corr_sum(i) = corr_sum(i) + corr(i)
      ENDDO
!      write(13,*) '    corrections to parameters (corr & corr_sum) applied to  this cycle'
!      write(13,4222)(corr(i)*frac_change(1),i=1,nv)
      WRITE (13,6777) icycle-1           
6777     format ('  Corrections to parameters (corr) applied to cycle',i3,' are...')
      WRITE (13,99066) (corr(i),i=1,nv)                                          ! 9-23-08 from S.P
      WRITE (13,6778) icycle-1   
6778     format ('  Corrections to parameters (corr_sum) applied to cycle',i3,' are...')
      WRITE (13,99066) (corr_sum(i),i=1,nv)                                      ! 9-30-08 DU
      IF ( IROt/=0 ) WRITE (13,"(4X,'New conformation in LS ',3F16.8)") &
     &                      (corr_sum(i)*57.29578,i=nv-nvr+1,nv)                 ! 10-10-08 DU
!
!      Write the final results of LS refinement                                  ! 9-23-08 from S.P
!
      IF ( ICYcle>3 ) THEN                                                       ! do at least 4 cycles of l.s, 6/29/08
! skip comparison if icycle = 1 6/11/07
! abs(E_diff) < 1.0D-9)) THEN   6/27/08...changed from -8 to -9
! test on # cycle & E_diff      6/11/07
         IF ( .NOT.(ICYcle==1) .AND.                                    &
     &        (ICYcle>ncycle .OR. ABS(e_diff)<e_diff_ls_test) ) THEN
                                                                                 ! 7/1/08...changed
            WRITE (13,99025) ICYcle - 1
99025       FORMAT (//,5x,                                              &
     &            'Cell parameters, E, VOL and density after LS  cycle '&
     &            ,i5)
            WRITE (13,99026) REFcode , a , b , c , alpha , beta ,       &
     &                       gamma , w01 , vol , DENsity
99026       FORMAT (/,5x,a8,3F10.4,3F10.5,f12.6,2F10.3,/)
            DO ik = 1 , 6
               cell_last(ik) = CELl1(ik)
               cell_diff(ik) = cell_last(ik) - cell_1st(ik)
            ENDDO
            if (IROSEN == 0) then                                                ! 12-22-08 DU
            write(12,99127)ICYcle - 1                                            ! 12-22-08 DU
            else                                                                 ! 12-22-08 DU
            write(12,99126)ICYcle - 1                                            ! 12-22-08 DU
            end if                                                               ! 12-22-08 DU
99126       FORMAT (/' Parameter changes and values after LS cycles =',i3,  &    ! 12-22-08 DU
     &        '...'//5X,'#',2X,'Parameter',13X,'after RSS',9X,'change',9X,'value',11X,'diff %'/)     ! 1-8-09 DU
99127       FORMAT (/' Parameter changes and values after LS cycles =',i3,  &                        ! 12-22-08 DU
     &        '...'//5X,'#',2X,'Parameter',13X,' initial ',9X,'change',9X,'value',11X,'diff %'/)     ! 1-8-09 DU
           DO ik = 1, 6                                                                              ! 12-22-08 DU
           if (IDX(ik) /= 0) then                                                                    ! 1-8-09   DU
           if (ik < 4) then                                                                          ! 12-22-08 DU
           WRITE (12,"(I6,1X,A15,1X,4F16.6)") IK, vab_name(ik),cell_1st(IK),&                        ! 1-8-09 DU
     &     cell_diff(ik), cell_last(IK),(cell_diff(ik)*100/cell_1st(IK))                             ! 1-8-09 DU
           else
           WRITE (12,"(I6,1X,A15,1X,4F16.6)") IK, vab_name(ik),ACOSD(cell_1st(IK)),&                 ! 1-8-09 DU
     &     (ACOSD(cell_last(IK)) - ACOSD(cell_1st(IK))), ACOSD(cell_last(IK)), &                     ! 1-8-09 DU
     &     ((ACOSD(cell_last(IK)) - ACOSD(cell_1st(IK)))*100/ACOSD(cell_1st(IK)))                    ! 1-8-09 DU
           end if
           end if                                                                                    ! 1-8-09
           END DO                                                                                    ! 12-22-08 DU
           DO ik = 7, nv                                                                             ! 12-22-08 DU
           if (IDX(ik) /= 0) then                                                                    ! 1-8-09   DU
            if ( IROt/=0 .AND. ik>(6+6*nmol) ) THEN                                                  ! 12-25-08 DU
            ii = INDX(ik-(6+6*nmol)) + 1                                                             ! 1-22--9 DU
           WRITE (12,"(I6,1X,A10,I2,', deg',F15.6,3F16.6)") &                                        ! 1-22-09 DU
     &     IK, bond_type(ii),ik-(6+6*nmol), diffpc(IK)*57.29578, corr_sum(IK)*57.29578, &            ! 1-9-09 DU
     &     (diffpc(IK)+corr_sum(IK))*57.29578, &
     &     (corr_sum(IK)*57.29578)*100/(diffpc(IK)*57.29578)                                         ! 1-8-09 DU
            else                                                                                     ! 1-8-09 DU
             if ( vab_name(ik)(5:9) .EQ. 'rotn,') then                                               ! 1-8-09 DU
             WRITE (12,"(I6,1X,A15,1X,4F16.6)") IK, vab_name(ik),diffpc(IK)*57.29578, &              ! 1-8-09 DU
     &       corr_sum(IK)*57.29578, (diffpc(IK)*57.29578+corr_sum(IK)*57.29578), &                   ! 1-8-09 DU
     &       ((corr_sum(IK)*57.29578)/(diffpc(IK)*57.29578))*100                                          
             else                                                                                    ! 1-8-09 DU
             WRITE (12,"(I6,1X,A15,1X,4F16.6)") IK, vab_name(ik),diffpc(IK), corr_sum(IK), &         ! 1-8-09 DU
     &       diffpc(IK)+corr_sum(IK),(corr_sum(IK)/diffpc(IK))*100                                   ! 1-8-09 DU
             end if                                                                                  ! 1-8-09 DU
            end if                                                                                   ! 1-8-09 DU
           end if                                                                                    ! 1-8-09   DU
           END DO                                                                                    ! 12-22-08 DU
            WRITE (12,99067) ICYcle - 1
            WRITE (13,99067) ICYcle - 1
            WRITE (12,99068) cell_last
            WRITE (13,99068) cell_last
            WRITE (23,"(6f9.5)") cell_last                                                           ! 3-31-08 DU
!        sumx5 = 0                                                                                   ! deleted, 11/15/08,HLA
!        sumy5 = 0                                                                                   ! deleted, 11/15/08,HLA
!        sumz5 = 0                                                                                   ! deleted, 11/15/08,HLA
            DO i = 1 , N
               WRITE (12,99069) ATOm(i) , i , X(i) , Y(i) , Z(i)
               WRITE (13,99069) ATOm(i) , i , X(i) , Y(i) , Z(i)
               WRITE (23,99069) ATOm(i) , i , X(i) , Y(i) , Z(i)
            ENDDO
!
            istart = 1
            DO imol = 1 , NMOl
               sumx = 0.0
               sumy = 0.0
               sumz = 0.0
               iend = istart + NATm(imol) - 1
               DO i = istart , iend
                  sumx = X(i) + sumx
                  sumy = Y(i) + sumy
                  sumz = Z(i) + sumz
               ENDDO
! Centriod of atoms for atoms of the group imol xc,yc,zc
               xc(imol) = REAL(sumx/NATm(imol))
               yc(imol) = REAL(sumy/NATm(imol))
               zc(imol) = REAL(sumz/NATm(imol))
!
               WRITE (12,99027) imol , xc(imol) , yc(imol) , zc(imol)
99027          FORMAT ('XTR',I1,5x,18x,3F9.5)
               WRITE (23,"('XTRA',5x,18x,3f9.5)") xc(imol) , yc(imol) , &
     &                zc(imol)                                                           ! 3-31-08 DU
               istart = istart + NATm(imol)
            ENDDO
!
!           WRITE (12,99028)                                                             ! 12-16-08
!99028      FORMAT (/,51x,'Total',6x,'Last cycle ',/50X,'after RSS',/)                   
!           WRITE (12,99070) (i,cell_name(i),corr_sum(i),corr(i),i=1,6)
!           WRITE (13,99070) (i,cell_name(i),corr_sum(i),corr(i),i=1,6)
!
            IF ( IROt==2 ) THEN                                                          ! 9-23-08 from S.P
               DO i = nv - nvr + 1 , nv                                                  ! 9-23-08 from S.P
                  corr_sum(i) = corr_sum(i)*180.0/PI                                     ! 9-23-08 from S.P
                  corr(i) = corr(i)*180.0/PI                                             ! 9-23-08 from S.P
               ENDDO                                                                     ! 9-23-08 from S.P
            ENDIF                                                                        ! 9-23-08 from S.P
!
            WRITE (12,99073) sumderiv_1st , sum_deriv1
            WRITE (13,99073) sumderiv_1st , sum_deriv1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                          9-16-08 DU
            WRITE (13,*)                                                &
     &   '                     initial         final           diffs'
            WRITE (13,*)                                                &
     &'                     exppt.           LS           A & degree    &
     &      %'
            DO i = 1 , 6
               IF ( i>3 ) THEN
                  cell_diff(i) = ACOSD(cell_last(i))                    &
     &                           - ACOSD(cell_ini(i))
                  cell_diff_per(i) = cell_diff(i)*100/ACOSD(cell_ini(i))
               ELSE
                  cell_diff(i) = cell_last(i) - cell_ini(i)
                  cell_diff_per(i) = cell_diff(i)*100/cell_ini(i)
               ENDIF
            ENDDO
            WRITE (13,99074) (cell_ini(i),cell_last(i),cell_diff(i),    &
     &                       cell_diff_per(i),i=1,6)
            if (IROSEN == 0) then
            WRITE (12,1277)
            else
            WRITE (12,1266)
            end if
1266           format (/28x,'initial        after RSS        final     difference',5X,'%')
1277           format (/28x,'initial         initial          final    difference',5X,'%')
            WRITE (12,1278)
1278           format ('                           values         values           values    Angs & degs')
            WRITE (12,*)
            WRITE (12,1253) (vab_name(i), cell_ini(i),cell_1st(i),cell_last(i),  &                         ! 1-22-09 DU
                             cell_diff(i),cell_diff_per(i),i=1,3)                                          ! 1-9-09 DU
            WRITE (12,1253) (vab_name(i), ACOSD(cell_ini(i)), ACOSD(cell_1st(i)),ACOSD(cell_last(i)), &    ! 1-22-09 DU
                             cell_diff(i),cell_diff_per(i),i=4,6)                                          ! 1-9-09 DU
1253            FORMAT (5x,a15,3F16.5,2X,2f10.5)
            DO imol = 1 , NMOl
               DO j = 1 , 6
                  rt_sum(j+6*imol) = DIFfpc(j+6*imol)                   &
     &                               + corr_sum(j+6*imol)
               ENDDO
               WRITE (12,99072) imol
               WRITE (13,99072) imol
               WRITE (12,99076)
               WRITE (13,99076)
               DO i = 1 , 6
                  IF ( i<4 ) THEN
                     WRITE (12,99075) p_name(i) , DIFfpc(i+6*imol)*57.29578 ,    &  ! 1-22-09
     &                                corr_sum(i+6*imol)*57.29578 ,              &  ! 1-22-09
     &                                rt_sum(i+6*imol)*57.29578                     ! 1-22-09
                     WRITE (13,99075) p_name(i) , DIFfpc(i+6*imol)*57.29578 ,    &
     &                                corr_sum(i+6*imol)*57.29578 ,              &
     &                                rt_sum(i+6*imol)*57.29578                   
                  ELSE
                     WRITE (12,99075) p_name(i) , DIFfpc(i+6*imol) ,    &
     &                                corr_sum(i+6*imol) ,              &
     &                                rt_sum(i+6*imol)                   
                     WRITE (13,99075) p_name(i) , DIFfpc(i+6*imol) ,    &
     &                                corr_sum(i+6*imol) ,              &
     &                                rt_sum(i+6*imol)                   
                  ENDIF
               ENDDO
            ENDDO
            rt = SQRT(((rt_sum(7)*57.29578)**2)+((rt_sum(8)*57.29578)**2) &         ! 12-16-08 DU
     &               +((rt_sum(9)*57.29578)**2))
            ts = SQRT((rt_sum(10)**2)+(rt_sum(11)**2)+(rt_sum(12)**2))              ! 12-16-08 DU
            f_factor = (0.5*rt)**2 + (10*ts)**2 + cell_diff_per(1)**2 + &           ! 12-16-08 DU
     &                 cell_diff_per(2)**2 + cell_diff_per(3)**2   + &              ! 12-16-08 DU
     &                 cell_diff(4)**2 + cell_diff(5)**2 + cell_diff(6)**2          ! 12-16-08 DU
            vol_diff_per = (vol-vol_ini)*100/vol_ini
            dens_diff_per = (DENsity-dens_ini)*100/dens_ini
!           WRITE (12,99029) vol_ini , vol , vol - vol_ini ,            &           ! 12-16-08
!    &                       vol_diff_per , dens_ini , DENsity ,        &
!    &                       DENsity - dens_ini , dens_diff_per
99029       FORMAT (/,5x,'vol     ',2F16.8,2F16.6/5x,'D       ',2F16.5, &
     &              2F16.5,/)
!                                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            WRITE (13,*)                                                &
     &   '                   initial           final         difference'
            WRITE (13,*) '                  after RSS '
            WRITE (13,99077) (cell_1st(i),cell_last(i),cell_diff(i),i=1,6)
            vol_diff = vol - VOL_1st
            dens_last = DENsity
            dens_diff = dens_last - DENs_1st
            ecvr_last(1) = ec
            ecvr_last(2) = ev
            ecvr_last(3) = er
            ecvr_last(4) = E
            IF ( IROt>0 ) ecvr_last(4) = E + DEL_egymol                   ! 9-23-08 from S.P
            DO i = 1 , 4
               ecvr_diff(i) = ecvr_last(i) - ecvr_1st(i)
            ENDDO
            WRITE (12,*) ' '
            WRITE (13,*) ' '
!
            WRITE (12,99330)
99330       FORMAT(19X,'initial',10X,'after RSS',8X,'final')              ! 12-16-08 DU
            WRITE (12,99030) (egy_ini(i),ecvr_1st(i),ecvr_last(i),i=1,3)  ! 12-16-08 DU
99030       FORMAT (/5x,'Ec    ',F16.6,2X,2F16.6/ &                       ! 12-16-08 DU
     &               5x,'Ev    ',F16.6,2X,2F16.6/ &                       ! 12-16-08 DU
     &               5x,'Er    ',F16.6,2X,2F16.6)                         ! 12-16-08 DU
            WRITE (12,99033) egy_ini(4),ecvr_1st(4),E,           &        ! 12-16-08 DU
     &                       vol_ini   ,VOL_1st    ,vol         ,&        ! 12-16-08 DU   
     &                       Dens_ini  ,DENs_1st   ,dens_last             ! 12-16-08 DU
99033       FORMAT (5x,'E_latt',F16.6,2X,2F16.6//  &                      ! 12-16-08 DU
     &              5x,'vol   ',F16.6,2X,2F16.6/   &                      ! 12-16-08 DU
     &              5x,'D     ',F16.6,2X,2F16.6/)                         ! 12-16-08 DU                       
            WRITE(12,99333) f_factor_rss,f_factor                         ! 12-16-08 DU
99333       FORMAT(5x,'F factor',16X,2F16.6)                              ! 12-16-08 DU 
!
            OPEN (UNIT=21,POSITION='append')
                                          ! unit # 21 = SUMRY file                 ! 6/11/07
            WRITE (21,99034) VOL_1st , vol , DENs_1st , dens_last ,     &
     &                       ecvr_last(4)                                          ! 6/11/07
99034       FORMAT (4x,'PMIN initial & final V, D, E =',2X,2F8.2,2F8.3, &
     &              F9.2,' MODE=  1')                                              ! 6/11/07
            CLOSE (21,STATUS='keep')                                               ! 6/11/07
!
!     if (irot >= 1) ECVR_last(4) = ECVR_last(4) + del_egymol                      ! 3-31-08 DU
!     if (irot == 1) ECVR_last(4) = ECVR_last(4) + del_egymol                      ! 9-23-08 from S.P
!     print 1252, icycle-1, dens_last, ECVR_last(4), E_diff, del_egymol            ! 3-31-08 DU
!1252 format ('PMIN Least Sqar refn: cycles =',i3, &                               ! 6/11/07
!                '; density =',f6.3,'; Ecalc =',f11.6, &
!                '; final E diff =',E13.6,' (E_conf =', F10.5,' )')
!
            IF ( IROt>0 ) THEN                                                     ! 9-30-08 DU
               PRINT 99035 , ICYcle - 1 , dens_last , E , DEL_egymol ,  &
     &               ecvr_last(4) , e_diff                                         ! 9-30-08 DU
99035          FORMAT ('PMIN Least Sqar refn: cycles =',i3,             &
     &                 '; density =',f6.3,'; Lat E =',f11.6,'; Con E =',&
     &                 f11.6,'; Tol E =',f11.6,'; final E diff =',E13.6)           ! 9-30-08 DU
            ELSE
               PRINT 99036 , ICYcle - 1 , dens_last , ecvr_last(4) ,    &
     &               e_diff                                                        ! 9-30-08 DU
99036          FORMAT ('PMIN Least sqrs refn: cycles =',i3,             &
     &                 '; density =',f6.3,'; Lat E =',f11.6,            &
     &                 '; final E diff =',E13.6)
            ENDIF                                                                  ! 9-30-08 DU
            WRITE (62,99078) ICYcle - 1 , sum_deriv1 ,                  &
     &                       sum_deriv_previous , E , param_fract ,     &
     &                       e_diff                                                ! 9-30-08 DU
            IF ( IMOde>0 ) PRINT 99037 , (ecvr_last(i),i=1,4)                      ! 6-23-08 DU
99037       FORMAT ('EC, EV, ER, ET: ',10X,4F11.6)                                 ! 6-23-08 DU
!                                                                     
!----Write fort.31 for nbslattice calcns                              6/11/07
!     SYMCD is the coordn geom. code from line # 5 of wmin.input      6/11/07
!                                                                     
            DO i = 4 , 6                                            
               cell_last(i) = ACOSD(cell_last(i))
            ENDDO                                                   
            OPEN (UNIT=31,FORM='formatted')                         
            IF ( SYMcd=='DA' .OR. SYMcd=='DB' .OR. SYMcd=='DC' .OR.     &
     &           SYMcd=='DD' .OR. SYMcd=='DE' ) THEN                
               WRITE (31,99038) cell_last
!  <--- final cell parameters                                        
99038          FORMAT ('RSS',6X,'1'/8X,'C',1X,6F10.4)            
            ELSEIF ( SYMcd=='DF' .OR. SYMcd=='DG' ) THEN
               WRITE (31,99039) cell_last                     
99039          FORMAT ('RSS',6X,'1'/8X,'F',1X,6F10.4)        
            ELSE                                           
               WRITE (31,99040) cell_last                
99040          FORMAT ('RSS',6X,'1'/8X,'P',1X,6F10.4)      
            ENDIF
            CLOSE (UNIT=31,STATUS='keep')                 
            CLOSE (UNIT=31,STATUS='keep')               
!
            GOTO 99999    ! END JOB
         ENDIF
      ENDIF
!
      DO i = 1 , nv
         DEL_param(i) = step1
         IF ( IROt==2 .AND. i>nv-nvr ) DEL_param(i) = 0.001         ! 9-23-08 from S.P
      ENDDO
      DO i = 4 , 6
         DEL_param(i) = step2
      ENDDO
      DO i = 7 , 9
         DEL_param(i) = step3
      ENDDO
      DO i = 10 , 12
         DEL_param(i) = step4
      ENDDO
!
! if(irot == 2) write(3,"('icycle nset irun     del_param      egymol     egymol_ori     del_egymol    Total')")     
      DO irun = 1 , 2
         DO i = 1 , nv                                                       ! modified
            DEL_param(i) = -DEL_param(i)
         ENDDO
!       for irun = 2, del_param is set positive
         DO i = 1 , N
            X(i) = X_Ori(i)
            Y(i) = Y_Ori(i)
            Z(i) = Z_Ori(i)
         ENDDO
         CELl(1) = a
         CELl(2) = b
         CELl(3) = c
         CELl(4) = alpha
         CELl(5) = beta
         CELl(6) = gamma
!
!        write(8,*) 'icycle,irun,step1,del_param,cell1(1) ', &               ! removed 7/11/08 HLA
!                   icycle,irun,step1,del_param(1),cell1(1)                  ! removed 7/11/08 HLA
!
         CALL NEW_XYZ(CELl1,X,Y,Z,cell_new1,xnew1,ynew1,znew1)
!
!       write(8,*) ' Just after call to NEW_XYZ'                             ! removed 7/11/08 HLA
!        write(8,*)'icycle,irun,step1,del_param,cell1(1),cell_new1(1,1) ',&  ! removed 7/11/08 HLA
!                   icycle,irun,step1,del_param(1),cell1(1),cell_new1(1,1)   ! removed 7/11/08 HLA
!        write(8,*)'icycle,irun,x(1),x_ori(1),xnew1(1,1) ',icycle,irun, &    ! removed 7/11/08 HLA
!                  x(1),x_ori(1),xnew1(1,1)                                  ! removed 7/11/08 HLA
!
!        GRAND LOOP STARTS HERE
!
         DO nset = 1 , nv            ! modified
            IF ( IDX(nset)==0 ) THEN
               IF ( irun==1 ) wn(nset) = w01
               IF ( irun==2 ) wp(nset) = w01
               CYCLE
                      ! if pc invariant, skip calculation
            ENDIF
            DO i = 1 , 6
               CELl(i) = cell_new1(nset,i)
            ENDDO
!           write(8,*)' nset = ',nset
            sumx_new = 0.0
            sumy_new = 0.0
            sumz_new = 0.0
            DO i = 1 , N
               X(i) = xnew1(nset,i)
               Y(i) = ynew1(nset,i)
               Z(i) = znew1(nset,i)
               sumx_new = sumx_new + X(i)
               sumy_new = sumy_new + Y(i)
               sumz_new = sumz_new + Z(i)
            ENDDO
!           xc_new = sumx_new/FLOAT(N)                    ! removed, 11/16/08
!           yc_new = sumy_new/FLOAT(N)                    ! removed, 11/16/08
!           zc_new = sumz_new/FLOAT(N)                    ! removed, 11/16/08
!
            NENtry = NENtry + 1
!
            IF ( IMOde>0 ) THEN
!        CALL USER_E_1(cell,x,y,z,Ec,Ev,Er,E_mol,E)
               CALL USER_E_1(CELl,X,Y,Z,ec,ev,er,e_nmol,E)
            ELSE
               CALL POT_E(CELl,X,Y,Z,ec,ev,er,e_nmol,E)   ! 9-23-08 from S.P
            ENDIF
!
!     if(icycle == 1 .and. nset <= 6)then
!        write(8,5821)icycle,irun,nentry,nset,cell,W01,E,W01-E                      ! removed 7/11/08 HLA
!5821     format(/,'icycle,irun,nentry,nset,cell,W01,E,del_E',4i5,6f12.8,3f16.10,/) ! removed 7/11/08 HLA
!        do i= 1,n                                                                  ! removed 7/11/08 HLA
!           write(8,5823)i,x_ori(i),y_ori(i),z_ori(i),x(i),y(i),z(i), &             ! removed 7/11/08 HLA
!                        x_ori(i)-x(i),y_ori(i)-y(i),z_ori(i)-z(i)                  ! removed 7/11/08 HLA
!5823           format('i,xyz_ori,xyz,diff xyz ',i5,9f13.9)                         ! removed 7/11/08 HLA
!        enddo                                                                      ! removed 7/11/08 HLA
!        write(8,582)Ec,Ev,Er,E,ECVR_1st,ECVR_1st(1)-Ec,ECVR_1st(2)-Ev,ECVR_1st(3)-Er,ECVR_1st(4)-E   ! removed 7/11/08 HLA
!582        format(/,5x, 'Energy = ',8f14.8,/)                                                        ! removed 7/11/08 HLA
!      endif
!
            IF ( IROt==2 .AND. nset>nv-nvr ) THEN
               CALL ORTHO_COD(N,CELl,X,Y,Z,xo,yo,zo)
               IF ( NE_type==1 ) THEN                                               ! 9-30-08 DU
                  CALL ENERGY_PM3(N,xo,yo,zo,egymol)                                ! 9-23-08 from S.P
               ELSE IF ( NE_type==2 ) THEN
                  CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)                              ! 9-30-08 DU
               ELSE
                  IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)            ! 4-28-09 DU
               ENDIF                                                                ! 9-30-08 DU
               DEL_egymol = egymol - egymol_ori
               E = E + DEL_egymol
               WRITE (3,99041) ICYcle , nset , irun , DEL_param(nset) , &
     &                         egymol , egymol_ori , DEL_egymol , E
!        Print 8436,icycle,nset,irun,del_param(nset),egymol,egymol_ori,del_egymol,E
99041          FORMAT (                                                 &
     &     ' icycle,nset,irun,del_param,egymol,egymol_ori,del_egymol,E '&
     &     ,3I5,5F14.8)
            ENDIF
!
!        Save energy value
            IF ( irun==1 ) wn(nset) = E
            IF ( irun==2 ) wp(nset) = E
         ENDDO
      ENDDO
!
      DO ii = 1 , nv                              ! modified
         IF ( IDX(ii)/=0 ) THEN                   ! 9-23-08 from S.P
            DO icell = 1 , 6
               CELl1(icell) = cell_new1(ii,icell)
            ENDDO
            DO i = 1 , N
               X(i) = xnew1(ii,i)
               Y(i) = ynew1(ii,i)
               Z(i) = znew1(ii,i)
            ENDDO
            CALL NEW_XYZ(CELl1,X,Y,Z,cell_new,xnew,ynew,znew)
!
            DO jj = 1 , nv
               IF ( IDX(jj)==0 ) THEN
                  wpp(ii,jj) = wp(ii)
                  CYCLE   ! if pc invariant, skip calculation
               ENDIF
               DO icell = 1 , 6
                  CELl(icell) = cell_new(jj,icell)
               ENDDO
               DO i = 1 , N
                  X(i) = xnew(jj,i)
                  Y(i) = ynew(jj,i)
                  Z(i) = znew(jj,i)
               ENDDO
               NENtry = NENtry + 1
!
               IF ( IMOde>0 ) THEN
                  CALL USER_E_1(CELl,X,Y,Z,ec,ev,er,e_nmol,E)                  ! 9-23-08 from S.P
               ELSE
                  CALL POT_E(CELl,X,Y,Z,ec,ev,er,e_nmol,E)                     ! 9-23-08 from S.P
               ENDIF
!
               IF ( IROt==2 .AND. (ii>nv-nvr .OR. jj>nv-nvr) ) THEN
                  CALL ORTHO_COD(N,CELl,X,Y,Z,xo,yo,zo)                        ! 9-23-08 from S.P
                  IF ( NE_type==1 ) THEN                                       ! 9-30-08 DU
                     CALL ENERGY_PM3(N,xo,yo,zo,egymol)                        ! 9-23-08 from S.P
                  ELSE IF ( NE_type==2 ) THEN
                     CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)                      ! 9-30-08 DU
                  ELSE
                     IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)    ! 4-28-09 DU
                  ENDIF                                                        ! 9-30-08 DU
!                 del_egymol = egymol - egymol0                                ! 9-23-08 from S.P
                  DEL_egymol = egymol - egymol_ori                             ! 9-23-08 from S.P
                  E = E + DEL_egymol                                           ! 9-23-08 from S.P
               ENDIF                                                           ! 9-23-08 from S.P
               wpp(ii,jj) = E
            ENDDO
         ENDIF
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Now calculate the 1st and 2nd derivatives     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      i = 0
      j = 0
!     sum_dias = 0.0                                                          ! deleted, 11/21/08, HLA
      sum_deriv1 = 0.0
      nset = nv                                                               ! modified
      WRITE (18,"('ICYCLE = ', I6)") ICYcle                                   ! 5-1-08 DU
      WRITE (18,99042)                                                        ! 9-23-08 DU
99042 FORMAT ('icycle',3X,'i',5X,'W01',12X,'WP',13X,'WN',13X,'W01-WP',  &
     &        9X,'W01-WN',9X,'WP-WN',10X,'deriv1')                            ! 9-23-08 DU
      DO iv1 = 1 , nv                                                         ! modified
         i = i + 1
         IF ( IDX(iv1)==0 ) THEN
            deriv1(i) = 0.0
            CYCLE      ! CHECK later
         ENDIF
         deriv1(i) = (wp(i)-wn(i))/(2*DEL_param(i))
         IF ( ICYcle<=5 ) THEN                                                ! 5-1-08 DU
            WRITE (18,99043) ICYcle , i , w01 , wp(i) , wn(i) ,         &
     &                       w01 - wp(i) , w01 - wn(i) , wp(i) - wn(i) ,&
     &                       deriv1(i)
99043       FORMAT (2I5,7F15.10)                                              ! 9-23-08 from S.P
         ENDIF
         sum_deriv1 = sum_deriv1 + ABS(deriv1(i))
      ENDDO
      IF ( ICYcle==1 ) THEN
         sumderiv_1st = sum_deriv1
         sum_deriv_previous = sum_deriv1                                      ! 6/2/08, HLA, added
      ENDIF
!        write(8,677)icycle,sum_deriv1                                        ! removed 7/11/08 HLA
!        write(8,674)(deriv1(i),i=1,nv)                                       ! removed 7/11/08 HLA
      WRITE (13,99044) ICYcle , sum_deriv1                                    ! modified
99044 FORMAT (/'Sum of 1st derivatives for cycle',i3,' is',f12.6,'....'/) 
      WRITE (13,99045) (deriv1(i),i=1,nv)                                     ! modified
99045 FORMAT (4x,6F15.10)                                                        ! 5-1-08 DU
99046 FORMAT (12F11.4)                                                        ! newline
!
!----Start least squares
      IF ( ICYcle==1 ) then
         WRITE (62,99047)
99047       FORMAT ('Start least squares...')
         write (62,1256) E
1256        format (' Start with geometry from previous RSS optimization, RSS E =',f13.8)
         go to 1258
      endif
      WRITE (62,99078) ICYcle - 1 , sum_deriv1 , sum_deriv_previous ,   &
                       E , param_fract , e_diff
      sum_deriv_previous2 = sum_deriv_previous                                ! 7/2/09
1258  sum_deriv_previous = sum_deriv1                                         ! 6/2/08, HLA, added
      IF ( ICYcle>=10 ) e_diff_ls_test = 1.0D-6                               ! 7/1/08, set e-diff test to d-6
!----If have done 99 cycles of l.s. w/o convergence, go back to RSS           ! 6/30/08, HLA
!     IF (ICYcle .le. 98 .and. icycle .gt. 1) THEN                            ! 7/2/09 
      IF (ICYcle .le. 98) THEN                                                ! 7/2/09 
     write (62,1818) icycle, sum_deriv1, sum_deriv_previous, param_fract    ! temp   7/2/09
1818    format ('at 1557..',i4,2f14.5,f5.1)                                 ! temp   7/2/09
         IF (l_allow_change) THEN  ! if true then allow param_fract to change   7/2/09
            IF ( (sum_deriv1 >= 1.0) .AND.                              &
     &           (sum_deriv1 >= 1.1*sum_deriv_previous2) ) THEN    ! changed from 2.0 to 1.1  7/2/09
               param_fract = 0.5                                              ! reduce param changes to 0.5    6/4/08, HLA
               WRITE (62,99048) sum_deriv1 , sum_deriv_previous ,       &
     &                          param_fract                                   ! 6/3/08, HLA, added
99048          FORMAT (' ...current 1st deriv sum (',f10.5,             &
     &                 ') > 2*previous sum (',f10.5,')...',             &
     &                 'return to RSS'/                                 &
     &                 '    l.s. parameter change multiplier now',f4.1)       ! 6/03/08, HLA, added
               DO i = 1 , 4                                                   ! 6/4/08, HLA, added
                  frac_change(i) = param_fract                                ! 6/3/08, HLA, added
               ENDDO                                                          ! 6/3/08, HLA, added
               IRSs_call = IRSs_call + 1                                      ! 7-7-08 DU
               IF ( IRSs_call>=2 ) e_diff_ls_test = 1.0D-2                    ! 1-20-09 DU
               ICYcle = 0                                                     ! 6/3/08, HLA, added
               GOTO 100                                                       ! 6/3/08, HLA, added
            ENDIF                                                             ! 6/3/08, HLA, added
            sum_deriv_previous = sum_deriv1                                   ! 6/3/08, HLA, added
!        write(18,*)' i,j,w01,wp(i),wp(j),wpp(i,j),deriv2(i,j)'               ! 9-23-08 from S.P
            WRITE (18,99049)                                                  ! 9-23-08 DU
99049       FORMAT (4X,'i',4X,'j',5X,'w01',14X,'wp(i)',12X,'wp(j)',10X, &
     &              'wpp(i,j)',15X,'deriv2(i,j)',1X,'del_param(i)',1X,  &
     &              'del_param(j)')                                           ! 9-23-08 DU
         ENDIF
         DO i = 1 , nset
            DO j = 1 , nset
               IF ( IDX(i)==0 .OR. IDX(j)==0 ) THEN
                  deriv2(i,j) = 0.0
                  CYCLE
               ENDIF
               IF ( i==j ) THEN
                  deriv2(i,i) = (wp(i)+wn(i)-2*w01)                     &
     &                          /(DEL_param(i)*DEL_param(i))
                  IF ( ICYcle<=5 ) THEN
!                    write(18,7721)i,j,w01,wn(i),wp(i),deriv2(i,i),del_param(i)
99050                FORMAT (4x,                                        &
     &                   'i,j,w01,wn(i),wp(i),deriv2(i,i),del_param(i):'&
     &                   ,2I5,3F17.12,F20.8,f10.6)                                      ! 9-23-08 DU
                  ENDIF
!              ELSE                                                                     ! 9-23-08 S.P
!               if(i == 1 .and. j == 2)then                                             ! 9-23-08 S.P
!               diff_W = WPP(i,j) + W01 - WP(i) - WP(j)                                 ! 9-23-08 S.P
!               write(18,7743)i,j,diff_W                                                ! 9-23-08 S.P
!               end if                                                                  ! 9-23-08 S.P
!7743             format(' i,j,diff_W ',2i5,e18.12)                                     ! 9-23-08 S.P
               ENDIF
!
               deriv2(i,j) = (wpp(i,j)+w01-wp(i)-wp(j))                 &
     &                       /(DEL_param(i)*DEL_param(j))
!              s_deriv2    =  WPP(i,j)+W01-WP(i)-WP(j)                                  ! 5-1-08 DU TEMP
!
               IF ( ICYcle<5 ) THEN                                                     ! 9-23-08 from S.P
                  WRITE (18,99051) i , j , w01 , wp(i) , wp(j) ,        &
     &                             wpp(i,j) , deriv2(i,j) , DEL_param(i)&
     &                             , DEL_param(j)                                       ! 9-23-08 from S.P
99051             FORMAT (2I5,4F17.12,F20.8,2F10.6)                                     ! 9-23-08 from S.P
!                 write(18,7724) s_deriv2,del_param(i),del_param(j)                     ! 9-23-08 from S.P
!7724             format('WPP(i,j)+W01-WP(i)-WP(j), del_param(i) & del_param(j) = ', &  ! 9-23-08 from S.P
!                 F20.16,2F10.6)                                                        ! 9-23-08 from S.P
               ENDIF
            ENDDO
         ENDDO
!
!     if(icycle <= 5) then                                ! removed 7/11/08 HLA
!        write(8,*) '   1st & 2nd derivatives are:'       ! removed 7/11/08 HLA
!        write(8,*) ' '                                   ! removed 7/11/08 HLA
!        write(8,674) (deriv1(i),i=1,nset)                ! removed 7/11/08 HLA
!        write(8,*) ' '                                   ! removed 7/11/08 HLA
!        do i=1,nset                                      ! removed 7/11/08 HLA
!           write(8,672) (deriv2(i,j),j=1,nset)           ! removed 7/11/08 HLA
!        enddo                                            ! removed 7/11/08 HLA
!     endif                                               ! removed 7/11/08 HLA
!
!     New insertions to be compatible with WMIN
!
         im = 0
         k2 = 0
         DO i = 1 , nv
            dyc(i) = -deriv1(i)
            DO j = 1 , nv
               im = im + 1
               btr(im) = deriv2(i,j)
            ENDDO
         ENDDO
         CALL MCPY(btr,an,nv,nv,0)
         CALL MCPY(dyc,vn,nv,1,0)
         nvv = nv*nv
         im = 1
         DO i = 1 , nv
            IF ( slim>=0.0 ) THEN
               dgan = an(im)
               IF ( dgan<=slim ) THEN
                  san(i) = 0.0
               ELSE
                  san(i) = 1.0/(DSQRT(ABS(dgan)))
               ENDIF
            ELSE
               san(i) = 1.0
            ENDIF
!         write(8,7112)icycle, i,im,an(im),dgan,slim,san(i)                ! removed 7/11/08 HLA
!7112         format(5x,'icycle,i,im,an(im),dgan,slim,san(i) ',3i5,4f12.6) ! removed 7/11/08 HLA
            im = im + nv + 1
         ENDDO
!
!    Scale matrix and vector
!
         im = 0
         DO i = 1 , nv
            vn(i) = vn(i)*san(i)
            DO j = 1 , nv
               im = im + 1
               an(im) = an(im)*san(i)*san(j)
            ENDDO
         ENDDO
!   Add constant (unity) to diagonal elements for Marquardt's compromise
         im = 1
         DO i = 1 , nv
            an(im) = an(im) + cmpr !cmpr = 0.1
            im = im + nv + 1
         ENDDO
! Determine Eigen values and Eigen vectors of matrix
         nvv = nv*nv
         DO i = 1 , nvv
            evec(i) = an(i)
         ENDDO
         IF ( nv==1 ) THEN
            eval(1) = an(1)
            GOTO 300
         ENDIF
!
! Temporay insertion of evec.dat from WMIN
!
99052    FORMAT (5x,f14.8)
!
         DO ik = 1 , nv
            IF ( IDX(ik)==0 ) THEN
               eval(ik) = 0.0
               evec((ik-1)*nv+ik) = 0.0
            ENDIF
         ENDDO
         IF ( ICYcle<=2 ) THEN
!        write(8,*) ' Before HOUSEH j,eval(j),evec(k)',icycle          ! removed 7/11/08 HLA
            k2 = 0
            DO j = 1 , nv
               K1 = k2 + 1
               k2 = k2 + nv
!     WRITE(8,200)J,EVAL(J),(EVEC(K),K=K1,K2)                          ! removed 7/11/08 HLA
! 200 FORMAT(1X,I3,E11.3,9F11.6/(15X,9F11.6))                          ! removed 7/11/08 HLA
            ENDDO
         ENDIF
         nv1 = nv
         CALL HOUSEH(evec,nv,nv,eval,tq)
         IF ( ICYcle<=2 ) THEN
!        write(8,*) ' After HOUSEH icycle j,eval(j),evec(k) ',icycle   ! removed 7/11/08 HLA
            k2 = 0
            DO j = 1 , nv
               K1 = k2 + 1
               k2 = k2 + nv
!     WRITE(8,200)J,EVAL(J),(EVEC(K),K=K1,K2)                          ! removed 7/11/08 HLA
            ENDDO
         ENDIF
!
!        New lines      04Jan08
!
         DO ik = nv , 1 , -1
            IF ( IDX(ik)==1 ) THEN
               kjk = ik
               EXIT
            ENDIF
         ENDDO
         IF ( kjk<nv ) THEN
            DO jkk = kjk + 1 , nv
               eval(jkk) = 0.0
!        write(8,*)kjk,jkk,eval(jkk)                                  ! removed 7/11/08 HLA
            ENDDO
         ENDIF
         nv = nv1
         CALL VALVEC(eval,tq,evec,ne,nv,nv)
         IF ( ICYcle<=2 ) THEN
!        write(8,*)' after VALVEC icycle, ne ',icycle,ne              ! removed 7/11/08 HLA
            k2 = 0
            DO j = 1 , nv
               K1 = k2 + 1
               k2 = k2 + nv
!     WRITE(8,200)J,EVAL(J),(EVEC(K),K=K1,K2)                         ! removed 7/11/08 HLA
            ENDDO
         ENDIF
      ELSE
         WRITE (62,99053) ICYcle , E , e_diff
99053    FORMAT (' ~~~l.s. cycle =',i3,', E =',f12.7,', E_diff =',e12.6,&
     &           '...return to RSS')
         WRITE (61,99054) ICYcle , E , e_diff
99054    FORMAT (/'~~~l.s. cycle =',i3,', E =',f12.7,', E_diff =',e12.6,&
     &           '...repeat RSS')
         param_fract = 1.0
         DO i = 1 , 4
            frac_change(i) = param_fract
         ENDDO
         IRSs_call = IRSs_call + 1                                    ! 7-7-08 DU
         IF ( IRSs_call>=2 ) e_diff_ls_test = 1.0D-2                  ! 1-20-09 DU
         ICYcle = 0
         GOTO 100
      ENDIF
!      if(ne .ne. 0) write(8,*) '*** VALVEC fails ***'                     ! removed 7/11/08 HLA
 300  nvu = nv
!            write(8,*)'icycle,eval(1),are,elim ',icycle,eval(1),are,elim  ! removed 7/11/08 HLA
!            write(8,710)eval                                              ! removed 7/11/08 HLA
!710             format(12f11.6)                                           ! removed 7/11/08 HLA
      elim = eval(1)*are
!
      DO k = 1 , nv
         IF ( eval(k)<0.0 ) eval(k) = 0.1
         IF ( eval(k)>elim ) eval(k) = 1.0/eval(k)
      ENDDO
      DO
!
!       write(8,710) (eval(i),i=1,nv)                                     ! removed 7/11/08 HLA
         ij = 0
         DO j = 1 , nv
            ji = j
            DO i = 1 , nv
               IF ( i>=j ) THEN
                  ki = i
                  kj = j
                  anij = 0.0
                  DO k = 1 , nvu
                     anij = anij + eval(k)*evec(ki)*evec(kj)
                     ki = ki + nv
                     kj = kj + nv
                  ENDDO
               ELSE
                  anij = an(ji)
               ENDIF
               ij = ij + 1
               an(ij) = anij
               ji = ji + nv
            ENDDO
         ENDDO
!
!       calculate parameter changes and their errors
!
!       write(8,*)'   After inversion icycle ',icycle         ! removed 7/11/08 HLA
!       write(8,*) ' '                                        ! removed 7/11/08 HLA
!       write(8,710)(vn(i),i=1,nv)                            ! removed 7/11/08 HLA
!       write(8,*) ' '                                        ! removed 7/11/08 HLA
!       write(8,710)(an(i),i=1,nvv)                           ! removed 7/11/08 HLA
!       write(8,*) ' '                                        ! removed 7/11/08 HLA
         im = 0
         DO j = 1 , nv
            pdj = 0.0
            DO i = 1 , nv
               im = im + 1
               pdj = pdj + vn(i)*an(im)
!           IF ( i.EQ.j ) THEN                                ! deleted, 11/21/08, HLA
!              diag = an(im)                                  ! deleted, 11/21/08, HLA
!              sigma = 0.01                                   ! deleted, 11/21/08, HLA
!              err(j) = DSQRT(diag)*san(j)*sigma              ! deleted, 11/21/08, HLA
!           ENDIF                                             ! deleted, 11/21/08, HLA
            ENDDO
            pd(j) = pdj(j)*san(j)
         ENDDO
!
!        Adjust parameters
!
!        write(8,235)icycle                                   ! removed 7/11/08 HLA
!235      format(/' Parameters after cycle',i3/)              ! removed 7/11/08 HLA
!     id = 0                                                  ! removed, 11/16/08, HLA
         ib = 0
         DO i = 1,nv                                        ! modified
            dkt(i) = frac_change(1)
         ENDDO
         DO i = 1 , nset
            psv(i) = p(i)
            IF ( IDX(i)==0 ) THEN
               corr(i) = 0.0
            ELSE
               corr(i) = pd(i)                                ! Change later
            ENDIF
            p(i) = p(i) + corr(i)*dkt(i)
         ENDDO
         ib = 0
!     id = 0               ! deleted, 11/20/08, HLA
         ik = 0
         DO i = 1 , nv
!          write(8,1478)i,p(i),pcmin(i),pcmax(i),ib                               ! removed 7/11/08 HLA
!1478       format(/,'parameter exceeds limit i,p,pcmin,pcmax,ib ',i5,3f14.8,i5)  ! removed 7/11/08 HLA
            IF ( p(i)<pcmin(i) .OR. p(i)>pcmax(i) ) ib = 1
         ENDDO
!       write(8,*) ' '                                                            ! removed 7/11/08 HLA
!         write(8,*)' Parameters and Corrections for Cycle,ib ',icycle,ib         ! removed 7/11/08 HLA
         IF ( ib==1 ) THEN
            DO ik = nvu , 1 , -1
               IF ( IDX(ik)/=0 ) THEN
                  nvu = ik - 1
                         ! No. of variable parameters reduced by 1
                  EXIT
               ENDIF
            ENDDO
!             icycle = icycle - 1
            DO i = 1 , nv
               p(i) = psv(i)
            ENDDO
!          write(8,*)' No. of variable parameters reduced to ',nvu               ! removed 7/11/08 HLA
            IF ( nvu==0 ) THEN
               DO i = 1 , nv
                  corr(i) = 0.0
               ENDDO
               EXIT
            ENDIF
            CYCLE
         ENDIF
         EXIT
      ENDDO
!
      IF ( ISYstem==1 ) corr(2) = corr(1)
      IF ( ISYstem==2 ) THEN
         corr(2) = corr(1)
         corr(3) = corr(1)
         if(ispl ==1)then                 ! from S.P Feb.5,09
          do imol = 1,nmol                ! from S.P Feb.5,09
           if(ispl_1(imol) == 1)then      ! from S.P Feb.5,09
           j = 6*imol + 4                 ! from S.P Feb.5,09
         corr(j+1) = corr(j)              ! from S.P Feb.5,09
         corr(j+2) = corr(j)              ! from S.P Feb.5,09
           endif                          ! from S.P Feb.5,09
         enddo                            ! from S.P Feb.5,09
         endif                            ! from S.P Feb.5,09
      ENDIF
      IF ( ISYstem==3 ) THEN
         corr(2) = corr(1)
         corr(3) = corr(1)
         corr(5) = corr(4)
         corr(6) = corr(4)
      ENDIF
!
!     alp1 = ACOSD(p(4))                                                ! deleted, 11/15/08,HLA
!     bet1 = ACOSD(p(5))                                                ! deleted, 11/15/08,HLA
!     gam1 = ACOSD(p(6))                                                ! deleted, 11/15/08,HLA
!
!         write(8,*)' Parameters and Corrections for Cycle ',icycle     ! removed 7/11/08 HLA
!
      DO i = 1 , N
         X(i) = X_Ori(i)
         Y(i) = Y_Ori(i)
         Z(i) = Z_Ori(i)
      ENDDO
!
      CALL CORR_XYZ(corr,frac_change,CELl1,X,Y,Z,cell_cor,x_cor,y_cor,  &
     &              z_cor)
!
!       write(8,522)icycle,cell_cor                                             ! removed 7/11/08 HLA
!522        format(/5x,'corrected cell parameters after cycle  ',i5/5x,6f12.6)  ! removed 7/11/08 HLA
!      if(icycle == 1) then                                                     ! removed 7/11/08 HLA
!         write(8,*)  'coordinates before and after CORR_XYZ'                   ! removed 7/11/08 HLA
!         do i=1,n                                                              ! removed 7/11/08 HLA
!            write(8,192) atom(i),i,x(i),y(i),z(i),x_cor(i),y_cor(i),z_cor(i)   ! removed 7/11/08 HLA
!192          format(a4,1x,i4,18x,6f9.5)                                        ! removed 7/11/08 HLA
!         enddo                                                                 ! removed 7/11/08 HLA
!      endif                                                                    ! removed 7/11/08 HLA
!
      DO i = 1 , 6
         AH(i) = cell_cor(i)
      ENDDO
!
      DO i = 1 , N
         X(i) = x_cor(i)
         Y(i) = y_cor(i)
         Z(i) = z_cor(i)
      ENDDO
      write(8,522)icycle,CELl1,cell_cor                                              ! 1-22-09 temp
522   format(/'corrected cell parameters after cycle ##### ',i5/5x,6f10.5,2X,6f10.5) ! 1-22-09 temp
      write(8,192) atom(1),x(1),y(1),z(1),x_cor(1),y_cor(1),z_cor(1)                 ! 1-22-09 temp
192  format(a4,1x,'   1',18x,6f9.5)                                                  ! 1-22-09 temp
      GOTO 100
99055 FORMAT (a)
99056 FORMAT (a/)
99057 FORMAT (1x,a6,'...',a6,1x,2f14.6,f8.3)                                         ! 10/6/09 
99058 FORMAT (a6)
!      write (8,519) symcd, refcode                                                  ! removed 7/11/08 HLA
!519    format(/4x,a2,1a8,/)                                                         ! removed 7/11/08 HLA
99059 FORMAT (4x,a2,1A8)                                                             ! 6/11/07
99060 FORMAT (36I1)
                                                       !                             ! 9-30-08 DU
99061 FORMAT ('PMIN Least squares ..warning E increase in LS of',e11.4, &
     &        ' > 1 kcal/mol, ','LS cycle =',i2,' and E =',f11.5/       &
     &'Restart PMIN with RSS followed by LS calculations...E of starting&
     & conf. ',F16.8)                                                                ! 9-30-08 DU
99062 FORMAT ('|',71('-'),'|')                                                       ! 5-25-08 DU
99063 FORMAT (' ----> Warning...E increase in LS > 1 kcal/mol',         &
     &        '...E diff =',e11.4,'...cycle',i3,' E =',                 &
     &        f15.7/'       Restart RSS and LS calculations...')
99064 FORMAT (/'  Icycle, nentry, unit cell params, volume, density...'/ &
              4x,2I3,6F11.4,F10.3,f7.3/                                  &
              '  Ec, Ev, Er, E_nmol, E_total, E_diff, sum_deriv1...'/   &
              4x,5F12.6,d16.8,f9.5/)
!422    format(5x,'icycle, nentry, Cell, Vol, Density ',2i5,6f11.6,2f10.5/ &
!             5x,'Ec, Ev, Er, Econf, E_nmol, E_total, E_diff, sum_deriv1 ',/6f12.8,d16.8,f10.5/)    ! 9-30-08 DU
99065 FORMAT (5x,'icycle, nentry, Cell, Vol, Density ',2I5,6F11.6,      &
     &        2F10.5/8X,'EC',14X,'EV',14X,'ER',14X,'Econf',11X,'E_nmol',&
     &        10X,'E_total',5X,'E_diff',11X,'sum_deriv1'/6F16.8,d16.8,  &
     &        f10.5/)                                                                               ! 9-30-08 DU
99066 FORMAT (5x,12F10.6)
99067 FORMAT (/' After least squares cycle #',i3,','  &
     &        ' new cell params & coordinates in ',                   &
     &        'wmin format follow...'/5x,'(file name = pmin_LSQ.save)'/)
99068 FORMAT (6F9.5)
!514         format(a6,i3,18x,3f9.5)
99069 FORMAT (a6,i3,18x,6F9.5)                   ! 9-23-08 from S.P
99070 FORMAT (5x,'Change in parameter No.',i5,3x,a6,f16.6,E16.6)
!9071 FORMAT (5x,'Change in parameter No.',i5,3x,'bond-rotation ',i2,   &
!    &        f16.6,E16.6,' in degrees ')                                                          ! 9-23-08 from S.P
99072 FORMAT (/3x,'For molecule',i2,'...'/)
99073 FORMAT (/' After least squares, initial and final sums of ', &
               '1st_derivatives are',f8.4,','f8.4)
99074 FORMAT (/,5x,'a       ',4F16.5/5x,'b       ',4F16.5/5x,'c       ',&
     &        4F16.5/5x,'cos_1   ',4F16.6/5x,'cos_2   ',4F16.6/5x,      &
     &        'cos_3   ',4F16.6)
99075 FORMAT (5x,a6,8X,'0.000000',2x,4F16.6)                                                       ! 12-16-08 DU
99076 FORMAT (19X,'initial',11X,'from RSS',8X,'from LS',7X,'sum (Angs/degs)'/)                     ! 12-16-08 DU
99077 FORMAT (/5x,'a     ',2F16.5,2X,f16.5/5x,'b     ',2F16.5,2X,      &
     &        f16.5/5x,'c     ',2F16.5,2X,f16.5/5x,'cos_1 ',2F16.6,2X,  &
     &        f16.5/5x,'cos_2 ',2F16.6,2X,f16.5/5x,'cos_3 ',2F16.6,2X,  &
     &        f16.5/)
99078 FORMAT (' Cycle =',i4,', sum 1st derivs =',f11.6,                 &
     &        ', previous sum =',f11.6,', E =',F13.8,                   &
     &        ', l.s. shift multr =',f4.1,', E_diff =',e12.6)
99999 END PROGRAM PMIN                          ! 4/9/08
!
!end----------------------program pmin-------------------------------end
!
      SUBROUTINE RECIP(X,Y)
!
!----GET RECIPROCAL CELL PARAMETERS FROM CRYSTAL LATTICE AND VICE VERSA
!    INPUT IS X OUTPUT IS Y
!
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 , DIMENSION(6) :: X , Y
      INTENT (IN) X
      INTENT (OUT) Y
!
! Local variables
!
      REAL*8 :: arg , rv , sina , sinb , sinc , v , x4 , x5 , x6
      REAL*8 :: DSQRT
!
      x4 = X(4)**2
      x5 = X(5)**2
      x6 = X(6)**2
      arg = 1.0 - x4 - x5 - x6 + 2.0*X(4)*X(5)*X(6)
      v = X(1)*X(2)*X(3)*DSQRT(arg)
      rv = 1.0/v
      sina = DSQRT(1.0-x4)
      sinb = DSQRT(1.0-x5)
      sinc = DSQRT(1.0-x6)
      Y(1) = X(2)*X(3)*sina*rv
      Y(2) = X(3)*X(1)*sinb*rv
      Y(3) = X(1)*X(2)*sinc*rv
      Y(4) = (X(5)*X(6)-X(4))/(sinb*sinc)
      Y(5) = (X(6)*X(4)-X(5))/(sinc*sina)
      Y(6) = (X(4)*X(5)-X(6))/(sina*sinb)
      END SUBROUTINE RECIP
!
!end---------------------subroutine recip-----------------------------end
!
! Program to get new fractional-coordinates rotated and translated
! orthogonal coordinates with  x' along a, z' along c* and y' normal to x'z'.
! if the input coordinates are in A and orthogonal, put a,b,c =1.0 & cell angles as 90.0
!
      SUBROUTINE NEW_XYZ(Cell1,X,Y,Z,Cell_new,Xnew,Ynew,Znew)
!
      USE PMIN_MODULE , ONLY:NMOl , NATm , IDX , ISYstem , N , ILS ,    &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , PI, ispl, ispl_1                                  ! from S.P Feb 5,09
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!*--NEW_XYZ2017
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      REAL*8 , DIMENSION(6) :: Cell1
      REAL*8 , DIMENSION(NVD,6) :: Cell_new
      REAL*8 , DIMENSION(MATM) :: X , Y , Z
      REAL*8 , DIMENSION(NVD,MATM) :: Xnew , Ynew , Znew
      INTENT (OUT) Xnew , Ynew , Znew
      INTENT (INOUT) Cell_new
!
! Local variables
!
      REAL*8 :: a , b , c , c1 , c11 , c2 , c22 , c3 , c33 , c51 , cc1 ,&
     &          cc2 , cc3 , dxc_all , dyc_all , dzc_all , egymol , phi ,&
     &          psi , r1 , r2 , s1 , s11 , s2 , s3 , s33 , ss1 , ss3 ,  &
     &          sumx , sumx_all , sumy , sumy_all , sumz , sumz_all ,   &
     &          theta , xc1 , xc1_all , xc2 , xc_all , yc1 , yc1_all ,  &
     &          yc2 , yc_all , zc1 , zc1_all , zc2 , zc_all
      REAL*8 , DIMENSION(6) :: cell2
      REAL*8 :: DCOS , DSIN , DSQRT
      REAL*8 , DIMENSION(NMOLD) :: dxc , dyc , dzc , xc , xco , yc ,    &
     &                             yco , zc , zco
      REAL :: FLOAT , REAL
      INTEGER :: i , iend , imol , istart , ivr , j , nv , nvr
      REAL*8 , DIMENSION(MATM) :: x2 , x3 , x4 , x6 , x7 , xo , y2 ,    &
     &                            y3 , y4 , y6 , y7 , yo , z2 , z3 ,    &
     &                            z4 , z6 , z7 , zo
      REAL*8 , DIMENSION(3,3) :: xtm
      REAL*8 , DIMENSION(MATM,3) :: xyz
!
!*** End of declarations rewritten by SPAG
!
!      integer :: i,j,imol,nv,nvt,istart,iend            ! modified
!
      nvr = NROtbond                                     ! 9-23-08 from S.P
      a = Cell1(1)
      b = Cell1(2)
      c = Cell1(3)
      c1 = Cell1(4)
      c2 = Cell1(5)
      c3 = Cell1(6)
!     alpha = ACOSD(c1)            ! deleted, 11/18/08, HLA  ! 5/13/08, HLA, changed to acosd
!     beta = ACOSD(c2)             ! deleted, 11/18/08, HLA  ! 5/13/08, HLA, changed to acosd
!     gamma = ACOSD(c3)            ! deleted, 11/18/08, HLA  ! 5/13/08, HLA, changed to acosd
      s1 = DSQRT(1.0-c1*c1)
      s2 = DSQRT(1.0-c2*c2)
      s3 = DSQRT(1.0-c3*c3)
!
!      treat the nv variables separately                 ! modified
!
      istart = 1
!
      nv = 6 + 6*NMOl
                                                         ! 9-23-08 from S.P
                                                         ! 9-23-08 from S.P
      IF ( IROt==2 ) nv = 6 + 6*NMOl + nvr               ! 9-23-08 from S.P
!
      DO imol = 1 , NMOl
!
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = X(i) + sumx
            sumy = Y(i) + sumy
            sumz = Z(i) + sumz
         ENDDO
! Centriod of atoms for atoms of the group imol xc,yc,zc
         xc(imol) = REAL(sumx/NATm(imol))
         yc(imol) = REAL(sumy/NATm(imol))
         zc(imol) = REAL(sumz/NATm(imol))
         IF ( NENtry==1 ) THEN
!      write(8,*)'imol,natm,sumx,sumy,sumz ',imol,natm(imol),sumx,sumy,sumz              ! removed 7/11/08 HLA
!      write(8,311)icycle,imol,xc(imol),yc(imol),zc(imol)                                ! removed 7/11/08 HLA
!311     format(/,5x,'icycle,imol and Centriod of the atoms of this group',2i5,3f10.5,/) ! removed 7/11/08 HLA
         ENDIF
!
! transform the centriod of atoms to origin
!
         xc2 = a*xc(imol)
         yc2 = b*yc(imol)
         zc2 = c*zc(imol)
         xco(imol) = xc2 + yc2*c3 + zc2*c2
         yco(imol) = yc2*s3 + zc2*(c1-c2*c3)/s3
         zco(imol) = zc2*DSQRT(s1*s1-((c2-c1*c3)/s3)**2)
!
         DO i = istart , iend
            x2(i) = a*(X(i)-xc(imol))
            y2(i) = b*(Y(i)-yc(imol))
            z2(i) = c*(Z(i)-zc(imol))
         ENDDO
         DO i = istart , iend
            x3(i) = x2(i) + y2(i)*c3 + z2(i)*c2
            y3(i) = y2(i)*s3 + z2(i)*(c1-c2*c3)/s3
            z3(i) = z2(i)*DSQRT(s1*s1-((c2-c1*c3)/s3)**2)
         ENDDO
!
         istart = istart + NATm(imol)
      ENDDO
!
      DO i = 1 , nv                                      ! modified
         DO j = 1 , 6
            Cell_new(i,j) = Cell1(j)
         ENDDO
         DO j = 1 , N
            Xnew(i,j) = X(j)
            Ynew(i,j) = Y(j)
            Znew(i,j) = Z(j)
         ENDDO
      ENDDO
!
      sumx_all = 0.0
      sumy_all = 0.0
      sumz_all = 0.0
      DO i = 1 , N
         sumx_all = X(i) + sumx_all
         sumy_all = Y(i) + sumy_all
         sumz_all = Z(i) + sumz_all
      ENDDO
!
      xc_all = sumx_all/FLOAT(N)
      yc_all = sumy_all/FLOAT(N)
      zc_all = sumz_all/FLOAT(N)
!
Labc: DO i = 1 , 3                                                       ! loop for varying a, b, c
!
         IF ( IDX(i)/=0 ) THEN
            IF ( i==1 ) THEN
               Cell_new(1,1) = a + DEL_param(1)
               xc1_all = xc_all*a/(a+DEL_param(1))
               dxc_all = xc1_all - xc_all
               DO j = 1 , N
                  Xnew(i,j) = X(j)*(a/(a+DEL_param(1))) - dxc_all
               ENDDO
!
               IF ( ISYstem==1 .OR. ISYstem==2 ) THEN                    !  For tetragonal,hexagonal, and Cubic systems, b = a
                  Cell_new(1,2) = Cell_new(1,1)
                  yc1_all = yc_all*a/(a+DEL_param(1))
                  dyc_all = yc1_all - yc_all
                  DO j = 1 , N
                     Ynew(i,j) = Y(j)*(a/(a+DEL_param(1))) - dyc_all
                  ENDDO
               ENDIF
!
               IF ( ISYstem==2 ) THEN                                    ! For Cubic system, c = a
                  Cell_new(1,3) = Cell_new(1,1)
                  zc1_all = zc_all*a/(a+DEL_param(1))
                  dzc_all = zc1_all - zc_all
                  DO j = 1 , N
                     Znew(i,j) = Z(j)*(a/(a+DEL_param(1))) - dzc_all
                  ENDDO
               ENDIF
!
               IF ( ISYstem==3 ) THEN                                    ! For trigonal system, a = b = c
                  CALL NEW_XYZ_TRIGONAL(i,N,Cell_new(1,1),c1,           &
     &                                  DEL_param(i),ATOm,X,Y,Z,cell2,  &
     &                                  xyz)
                  DO j = 1 , N
                     Xnew(i,j) = xyz(j,1)
                     Ynew(i,j) = xyz(j,2)
                     Znew(i,j) = xyz(j,3)
                  ENDDO
               ENDIF
            ENDIF
!
            IF ( i==2 ) THEN
               Cell_new(2,2) = b + DEL_param(2)
               DO imol = 1 , NMOl
                  yc1 = yc(imol)*b/(b+DEL_param(2))
                  dyc(imol) = yc1 - yc(imol)
               ENDDO
               yc1_all = (yc_all*b)/(b+DEL_param(2))
               dyc_all = yc1_all - yc_all
               istart = 1
               DO imol = 1 , NMOl
                  iend = istart + NATm(imol) - 1
                  DO j = istart , iend
                     Ynew(i,j) = Y(j)*(b/(b+DEL_param(2))) - dyc_all          !  - dyc(imol)
                  ENDDO
                  istart = istart + NATm(imol)
               ENDDO
            ENDIF
            IF ( i==3 ) THEN
               Cell_new(3,3) = c + DEL_param(3)
               DO imol = 1 , NMOl
                  zc1 = zc(imol)*c/(c+DEL_param(3))
                  dzc(imol) = zc1 - zc(imol)
               ENDDO
               zc1_all = (zc_all*c)/(c+DEL_param(3))
               dzc_all = zc1_all - zc_all
               istart = 1
               DO imol = 1 , NMOl
                  iend = istart + NATm(imol) - 1
                  DO j = istart , iend
                     Znew(i,j) = Z(j)*(c/(c+DEL_param(3))) - dzc_all          ! - dzc(imol)
                  ENDDO
                  istart = istart + NATm(imol)
               ENDDO
            ENDIF
         ENDIF
!
      ENDDO Labc
!
!     get new coordinates after increasing cosine of cell angles by del_param
!
      Cell_new(4,4) = Cell1(4)
      Cell_new(5,5) = Cell1(5)
      Cell_new(6,6) = Cell1(6)
!
Labg: DO i =4,6                                                              ! loop for varying alpha, beta, gamma
!
         cc1 = c1
         cc2 = c2
         cc3 = c3
         ss1 = s1
!        ss2 = s2                                                            ! deleted, 11/18/08, HLA
         ss3 = s3
         IF ( IDX(i)/=0 ) THEN
!
            IF ( i==4 ) THEN
               IF ( ISYstem/=3 ) THEN
                  Cell_new(4,4) = Cell1(4) + DEL_param(4)
                  cc1 = Cell_new(4,4)
                  ss1 = DSQRT(1.0-cc1*cc1)
                  r1 = (cc1-cc2*cc3)/ss3
                  r2 = SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
!
                  istart = 1
                  DO imol = 1 , NMOl
                     zc1 = zco(imol)/r2
                     yc1 = (yco(imol)-zc1*r1)/ss3
                     xc1 = xco(imol) - yc1*cc3 - zc1*cc2
                     dxc(imol) = xc1/a - xc(imol)
                     dyc(imol) = yc1/b - yc(imol)
                     dzc(imol) = zc1/c - zc(imol)
                     iend = istart + NATm(imol) - 1
                     DO j = istart , iend
                        z6(j) = z3(j)/r2
                        y6(j) = (y3(j)-z6(j)*r1)/ss3
                        x6(j) = x3(j) - y6(j)*cc3 - z6(j)*cc2
                        Xnew(4,j) = (x6(j)+xc1)/a - dxc(imol)
                        Ynew(4,j) = (y6(j)+yc1)/b - dyc(imol)
                        Znew(4,j) = (z6(j)+zc1)/c - dzc(imol)
                     ENDDO
                     istart = istart + NATm(imol)
                  ENDDO
               ENDIF
!
!
               IF ( ISYstem==3 ) THEN  !   for trigonal system only
                  CALL NEW_XYZ_TRIGONAL(i,N,Cell_new(1,1),c1,           &
     &                                  DEL_param(i),ATOm,X,Y,Z,cell2,  &
     &                                  xyz)
                  DO j = 1 , 6
                     Cell_new(i,j) = cell2(j)
                  ENDDO
                  DO j = 1 , N
                     Xnew(i,j) = xyz(j,1)
                     Ynew(i,j) = xyz(j,2)
                     Znew(i,j) = xyz(j,3)
                  ENDDO
               ENDIF
            ENDIF
!
            IF ( i==5 ) THEN
               Cell_new(5,5) = Cell1(5) + DEL_param(5)
               cc2 = Cell_new(5,5)
!           ss2 = DSQRT(1.0-cc1*cc1)                                         ! deleted, 11/18/08, HLA
               r1 = (cc1-cc2*cc3)/ss3
               r2 = SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
!
               istart = 1
               DO imol = 1 , NMOl
                  zc1 = zco(imol)/r2
                  yc1 = (yco(imol)-zc1*r1)/ss3
                  xc1 = xco(imol) - yc1*cc3 - zc1*cc2
                  dxc(imol) = xc1/a - xc(imol)
                  dyc(imol) = yc1/b - yc(imol)
                  dzc(imol) = zc1/c - zc(imol)
                  iend = istart + NATm(imol) - 1
                  DO j = istart , iend
                     z6(j) = z3(j)/r2
                     y6(j) = (y3(j)-z6(j)*r1)/ss3
                     x6(j) = x3(j) - y6(j)*cc3 - z6(j)*cc2
                     Xnew(5,j) = (x6(j)+xc1)/a - dxc(imol)
                     Ynew(5,j) = (y6(j)+yc1)/b - dyc(imol)
                     Znew(5,j) = (z6(j)+zc1)/c - dzc(imol)
                  ENDDO
                  istart = istart + NATm(imol)
               ENDDO
!
            ENDIF
!
            IF ( i==6 ) THEN
               Cell_new(6,6) = Cell1(6) + DEL_param(6)
               cc3 = Cell_new(6,6)
               ss3 = DSQRT(1.0-cc3*cc3)
               r1 = (cc1-cc2*cc3)/ss3
               r2 = SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
!
               istart = 1
               DO imol = 1 , NMOl
                  zc1 = zco(imol)/r2
                  yc1 = (yco(imol)-zc1*r1)/ss3
                  xc1 = xco(imol) - yc1*cc3 - zc1*cc2
                  dxc(imol) = xc1/a - xc(imol)
                  dyc(imol) = yc1/b - yc(imol)
                  dzc(imol) = zc1/c - zc(imol)
                  iend = istart + NATm(imol) - 1
                  DO j = istart , iend
                     z6(j) = z3(j)/r2
                     y6(j) = (y3(j)-z6(j)*r1)/ss3
                     x6(j) = x3(j) - y6(j)*cc3 - z6(j)*cc2
                     Xnew(6,j) = (x6(j)+xc1)/a - dxc(imol)
                     Ynew(6,j) = (y6(j)+yc1)/b - dyc(imol)
                     Znew(6,j) = (z6(j)+zc1)/c - dzc(imol)
                  ENDDO
                  istart = istart + NATm(imol)
               ENDDO
!
            ENDIF
         ENDIF
!
      ENDDO Labg
!
      s11 = s1
!     s22 = s2                                                               ! deleted, 11/18/08, HLA
      s33 = s3
      c11 = c1
      c22 = c2
      c33 = c3
      r1 = (c11-c22*c33)/s33
      r2 = DSQRT(s11*s11-((c22-c11*c33)/s33)**2)
!
!   r1 and r2 are used later for converting rotated orthogonal coordinates to Angstrom coordinates
!   with respect to the cell edges.
!
      istart = 1
!
      DO imol = 1 , NMOl
!
!
! Now Calculate the orthognal coordinates x2,y2,z2 rotated by Euler angles phi,theta,psi
! Take the coordinates referred to centriod as origin
!
         DO j = 6*imol + 1 , 6*imol + 3
!
            phi = 0.0
            theta = 0.0
            psi = 0.0
            IF ( IDX(j)/=0 ) THEN
                                  ! Do not generate coordinates for variables kept constant
               IF ( j==6*imol+1 ) phi = DEL_param(6*imol+1)
               IF ( j==6*imol+2 ) theta = DEL_param(6*imol+2)
               IF ( j==6*imol+3 ) psi = DEL_param(6*imol+3)
!
!      the angles phi, theta, psi are in radians
!
               s1 = DSIN(phi)
               s2 = DSIN(theta)
               s3 = DSIN(psi)
               c51 = DCOS(phi)
               c2 = DCOS(theta)
               c3 = DCOS(psi)
!
               xtm(1,1) = c3*c2
               xtm(1,2) = -s3*c2
               xtm(1,3) = s2
               xtm(2,1) = c3*s2*s1 + s3*c51
               xtm(2,2) = c3*c51 - s3*s2*s1
               xtm(2,3) = -c2*s1
               xtm(3,1) = -c3*s2*c51 + s3*s1
               xtm(3,2) = c3*s1 + s3*s2*c51
               xtm(3,3) = c2*c51
!
!
!       The matrix below gives the coordinates rotated  by angles phi, theta, psi
!       which are inremental angles and only one is varied at a time.
!
               iend = istart + NATm(imol) - 1
               DO i = istart , iend
                  x4(i) = x3(i)*xtm(1,1) + y3(i)*xtm(2,1) + z3(i)       &
     &                    *xtm(3,1)
                  y4(i) = x3(i)*xtm(1,2) + y3(i)*xtm(2,2) + z3(i)       &
     &                    *xtm(3,2)
                  z4(i) = x3(i)*xtm(1,3) + y3(i)*xtm(2,3) + z3(i)       &
     &                    *xtm(3,3)
               ENDDO
!
!  convert rotated orthogonal coordinates to crystallographic cell in Angs
!
               iend = istart + NATm(imol) - 1
               DO i = istart , iend
                  z6(i) = z4(i)/r2
                  y6(i) = (y4(i)-z6(i)*r1)/s33
                  x6(i) = x4(i) - y6(i)*c33 - z6(i)*c22
!
!     convert to fractional  and move the molecule centre to the original position
!
                  Xnew(j,i) = x6(i)/a + xc(imol)
                  Ynew(j,i) = y6(i)/b + yc(imol)
                  Znew(j,i) = z6(i)/c + zc(imol)
               ENDDO
            ENDIF
!
         ENDDO
!
!      Now get new-coordinates for translation of the rigid body by del_param in Angs
         if(isystem == 2 .and. ispl == 1)then                                    ! from S.P Feb 5,09
            if(ispl_1(imol) == 1)then                                            ! from S.P Feb 5,09
             j = 6*imol + 4                                                      ! from S.P Feb 5,09
             if(idx(j) == 0) go to 5033                                          ! from S.P Feb 5,09
             iend = istart + natm(imol)-1                                        ! from S.P Feb 5,09
           do i = istart,iend                                                    ! from S.P Feb 5,09
             xnew(j,i) = x(i) + del_param(j)*xc(imol)/abs(xc(imol))/sqrt(3.0)    ! from S.P Feb 5,09
             ynew(j,i) = y(i) + del_param(j)*yc(imol)/abs(yc(imol))/sqrt(3.0)    ! from S.P Feb 5,09
             znew(j,i) = z(i) + del_param(j)*zc(imol)/abs(zc(imol))/sqrt(3.0)    ! from S.P Feb 5,09
           enddo                                                                 ! from S.P Feb 5,09
             endif                                                               ! from S.P Feb 5,09
5033      continue                                                               ! from S.P Feb 5,09
          RETURN                                                                 ! from S.P Feb 5,09
         endif                                                                   ! from S.P Feb 5,09
!
         DO j = 6*imol + 4 , 6*imol + 6
            IF ( IDX(j)/=0 ) THEN
               IF ( j==6*imol+4 ) THEN
                  iend = istart + NATm(imol) - 1
                  DO i = istart , iend
                     Xnew(j,i) = X(i) + DEL_param(6*imol+4)/a
                     Ynew(j,i) = Y(i)
                     Znew(j,i) = Z(i)
                  ENDDO
               ENDIF
               IF ( j==6*imol+5 ) THEN
                  iend = istart + NATm(imol) - 1
                  DO i = istart , iend
                     Xnew(j,i) = X(i)
                     Ynew(j,i) = Y(i) + DEL_param(6*imol+5)/b
                     Znew(j,i) = Z(i)
                  ENDDO
               ENDIF
               IF ( j==6*imol+6 ) THEN
                  iend = istart + NATm(imol) - 1
                  DO i = istart , iend
                     Xnew(j,i) = X(i)
                     Ynew(j,i) = Y(i)
                     Znew(j,i) = Z(i) + DEL_param(6*imol+6)/c
                  ENDDO
               ENDIF
            ENDIF
!
         ENDDO
!
         istart = istart + NATm(imol)
      ENDDO
!
      IF ( IROt==2 ) THEN                                             ! 9-23-08 from S.P
         DO j = nv - nvr + 1 , nv                                     ! 9-23-08 from S.P
            IF ( IDX(j)/=0 ) THEN                                     ! 9-23-08 from S.P
               ivr = j - (nv-nvr)                                     ! 9-23-08 from S.P
               DTHeta = DEL_param(j)                                  ! 9-23-08 from S.P
!             Print *,'j,del_param,dtheta ',j,del_param(j),dtheta     ! 9-23-08 from S.P
               CALL BEND_BOND_NEW(ivr,N,Cell1,X,Y,Z,xo,yo,zo,egymol)  ! 9-23-08 from S.P
               CALL FRACTIONAL_COD(N,Cell1,xo,yo,zo,x7,y7,z7)         ! 9-23-08 from S.P
!             del_egymol = egymol - egymol0                           ! 9-23-08 from S.P
!             write(3,8436)ils,ik,dtheta,egymol,egymol0,del_egymol    ! 9-23-08 from S.P
!             Print 8436,ivr,dtheta,egymol                            ! 9-23-08 from S.P
!             write(8,8436)ivr,dtheta,egymol                          ! 9-23-08 from S.P
99001          FORMAT (' ivr,dtheta,egymol ',i5,f9.5,3F10.6,f14.6)    ! 9-23-08 from S.P
               DO i = 1 , N                                           ! 9-23-08 from S.P
                  Xnew(j,i) = x7(i)                                   ! 9-23-08 from S.P
                  Ynew(j,i) = y7(i)                                   ! 9-23-08 from S.P
                  Znew(j,i) = z7(i)                                   ! 9-23-08 from S.P
               ENDDO                                                  ! 9-23-08 from S.P
            ENDIF
         ENDDO                                                        ! 9-23-08 from S.P
      ENDIF                                                           ! 9-23-08 from S.P
!
!
      END SUBROUTINE NEW_XYZ
!
!end-----------------------subroutine new_xyz-----------------------------end
!
      SUBROUTINE CORR_XYZ(Corr,Frac_change,Cell1,X,Y,Z,Cell_new,X_cor,  &
     &                    Y_cor,Z_cor)
!
      USE PMIN_MODULE , ONLY:NMOl , NATm , IDX , ISYstem , N , ILS ,    &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , PI , ANG_TO_RADIAN          ! 9-23-08 from S.P
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 , DIMENSION(6) :: Cell1 , Cell_new
      REAL*8 , DIMENSION(NVD) :: Corr
      REAL*8 , DIMENSION(4) :: Frac_change
      REAL*8 , DIMENSION(MATM) :: X , X_cor , Y , Y_cor , Z , Z_cor
      INTENT (IN) Corr , Frac_change
      INTENT (INOUT) Cell_new , X_cor , Y_cor , Z_cor
!
! Local variables
!
      REAL*8 :: a , b , c , c2 , c3 , c61 , cc1 , cc2 , cc3 , egymol ,  &
     &          phi , psi , r1 , r2 , s1 , s2 , s3 , ss1 , ss3 , sum1 , &
     &          sum2 , sum3 , sumx , sumy , sumz , theta
      INTEGER :: i , iend , ik , imol , istart , ivr , j , nv , nvr
      REAL*8 , DIMENSION(MATM) :: x3 , x4 , x6 , xo , xx2 , y3 , y4 ,   &
     &                            y6 , yo , yy2 , z3 , z4 , z6 , zo ,   &
     &                            zz2
      REAL*8 , DIMENSION(NMOl) :: xcc , ycc , zcc
      REAL*8 , DIMENSION(3,3) :: xtm
!
      nv = 6 + 6*NMOl                          ! 9-23-08 from S.P
      IF ( IROt==2 ) THEN                      ! For LS only      9-23-08 from S.P
         nvr = NROtbond                        ! 9-23-08 from S.P
         nv = 6 + 6*NMOl + nvr                 ! 9-23-08 from S.P
      ENDIF                                    ! 9-23-08 from S.P
!        write (fmt2(36:37),'(i2)') nv         ! removed 7/11/08 HLA
!        write (8,fmt2) nv,(corr(i),i=1,nv)    ! removed 7/11/08 HLA
      a = Cell1(1)
      b = Cell1(2)
      c = Cell1(3)
!     alpha = ACOS(Cell1(4)*ANG_TO_RADIAN)     ! deleted, 11/18/08, ! 5/13/08, HLA, altered
!     beta = ACOS(Cell1(5)*ANG_TO_RADIAN)      ! deleted, 11/18/08, ! 5/13/08, HLA, altered
!     gamma = ACOS(Cell1(6)*ANG_TO_RADIAN)     ! deleted, 11/18/08, ! 5/13/08, HLA, altered
      cc1 = Cell1(4)
      cc2 = Cell1(5)
      cc3 = Cell1(6)
      ss1 = SQRT(1.0-cc1*cc1)
!     ss2 = SQRT(1.0-cc2*cc2)                  ! deleted, 11/18/08, HLA
      ss3 = SQRT(1.0-cc3*cc3)
!
!    First apply the corrections for bond rotation and bending                  ! 9-23-08 from S.P
!
      IF ( IROt==2 ) THEN                                                       ! 9-23-08 from S.P
         DO ik = 1 , nv                                                         ! 9-23-08 from S.P
            IF ( ik>nv-nvr ) THEN                                               ! 9-23-08 from S.P
               ivr = ik - (nv-nvr)                                              ! 9-23-08 from S.P
               DTHeta = Corr(ik)                                                ! 9-23-08 from S.P
               CALL BEND_BOND_NEW(ivr,N,Cell1,X,Y,Z,xo,yo,zo,egymol)            ! 9-23-08 from S.P
               CALL FRACTIONAL_COD(N,Cell1,xo,yo,zo,X,Y,Z)                      ! 9-23-08 from S.P
!              del_egymol = egymol - egymol0                                    ! 9-23-08 from S.P
!              write(3,8436)ik,dtheta,egymol,egymol0,del_egymol                 ! 9-23-08 from S.P
99001          FORMAT (' i,dtheta,egymol,egymol0,del_egymol ',i5,f9.5,  &
     &                 3F14.8)                                                  ! 9-23-08 from S.P
            ENDIF                                                               ! 9-23-08 from S.P
         ENDDO                                                                  ! 9-23-08 from S.P
      ENDIF                                                                     ! 9-23-08 from S.P
!
! Now Calculate the orthognal coordinates xx2,yy2,zz2 rotated by Euler angles phi,theta,psi
! Take the coordinates referred to centriod as origin
!
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
!
      istart = 1
      DO imol = 1 , NMOl
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = X(i) + sumx
            sumy = Y(i) + sumy
            sumz = Z(i) + sumz
         ENDDO
         xcc(imol) = sumx/NATm(imol)
         ycc(imol) = sumy/NATm(imol)
         zcc(imol) = sumz/NATm(imol)
!        xco(imol) = xcc(imol)*a                                                      ! deleted, 11/18/08, HLA
!        yco(imol) = ycc(imol)*b                                                      ! deleted, 11/18/08, HLA
!        zco(imol) = zcc(imol)*c                                                      ! deleted, 11/18/08, HLA
         IF ( NENtry==1 ) THEN
!           write(8,338)imol,natm(imol),sumx,sumy,sumz,xcc(imol),ycc(imol),zcc(imol)  ! removed 7/11/08 HLA
!338            format(/5x,'imol,natm(imol),sumx,sumy,sumz,xcc,ycc,zcc ',2i5,6f12.8/) ! removed 7/11/08 HLA
         ENDIF
         DO i = istart , iend
            xx2(i) = a*(X(i)-xcc(imol))
            yy2(i) = b*(Y(i)-ycc(imol))
            zz2(i) = c*(Z(i)-zcc(imol))
         ENDDO
         DO i = istart , iend
            x3(i) = xx2(i) + yy2(i)*cc3 + zz2(i)*cc2
            y3(i) = yy2(i)*ss3 + zz2(i)*(cc1-cc2*cc3)/ss3
            z3(i) = zz2(i)*SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
         ENDDO
!        xco3 = xco(imol) + yco(imol)*cc3 + zco(imol)*cc2                             ! deleted, 11/18/08, HLA
!        yco3 = yco(imol)*ss3 + zco(imol)*(cc1-cc2*cc3)/ss3                           ! deleted, 11/18/08, HLA
!        zco3 = zco(imol)*SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)                        ! deleted, 11/18/08, HLA
!
         phi = 0.0
         theta = 0.0
         psi = 0.0
         phi = Corr(6*imol+1)*Frac_change(3)
         theta = Corr(6*imol+2)*Frac_change(3)
         psi = Corr(6*imol+3)*Frac_change(3)
!
!----These angles are in radians
         s1 = SIN(phi)
         s2 = SIN(theta)
         s3 = SIN(psi)
         c61 = COS(phi)
         c2 = COS(theta)
         c3 = COS(psi)
!
         xtm(1,1) = c3*c2
         xtm(1,2) = -s3*c2
         xtm(1,3) = s2
         xtm(2,1) = c3*s2*s1 + s3*c61
         xtm(2,2) = c3*c61 - s3*s2*s1
         xtm(2,3) = -c2*s1
         xtm(3,1) = -c3*s2*c61 + s3*s1
         xtm(3,2) = c3*s1 + s3*s2*c61
         xtm(3,3) = c2*c61
!
!----Rotate by phi, theta, psi and translate to xcc+corr(11) etc
         DO j = istart , iend
            x4(j) = x3(j)*xtm(1,1) + y3(j)*xtm(2,1) + z3(j)*xtm(3,1)
            y4(j) = x3(j)*xtm(1,2) + y3(j)*xtm(2,2) + z3(j)*xtm(3,2)
            z4(j) = x3(j)*xtm(1,3) + y3(j)*xtm(2,3) + z3(j)*xtm(3,3)
         ENDDO
!
         IF ( ILS==0 ) THEN
            DO i = istart , iend
               x4(i) = x4(i) + Corr(6*imol+4)*Frac_change(4)
               y4(i) = y4(i) + Corr(6*imol+5)*Frac_change(4)
               z4(i) = z4(i) + Corr(6*imol+6)*Frac_change(4)
            ENDDO
         ENDIF
         istart = istart + NATm(imol)
      ENDDO
!
!----Convert rotated orthogonal coordinates to crystallographic cell in Angs
      Cell_new(1) = a + Corr(1)*Frac_change(1)
      Cell_new(2) = b + Corr(2)*Frac_change(1)
      Cell_new(3) = c + Corr(3)*Frac_change(1)
      Cell_new(4) = Cell1(4) + Corr(4)*Frac_change(2)
      Cell_new(5) = Cell1(5) + Corr(5)*Frac_change(2)
      Cell_new(6) = Cell1(6) + Corr(6)*Frac_change(2)
!
      cc1 = Cell_new(4)
      cc2 = Cell_new(5)
      cc3 = Cell_new(6)
      ss1 = SQRT(1.0-cc1*cc1)
!     ss2 = SQRT(1.0-cc2*cc2)              ! deleted, 11/18/08, HLA
      ss3 = SQRT(1.0-cc3*cc3)
!
      r1 = (cc1-cc2*cc3)/ss3
      r2 = SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
!
      istart = 1
      DO imol = 1 , NMOl
         iend = istart + NATm(imol) - 1
!
         DO j = istart , iend
            z6(j) = z4(j)/r2
            y6(j) = (y4(j)-z6(j)*r1)/ss3
            x6(j) = x4(j) - y6(j)*cc3 - z6(j)*cc2
         ENDDO
!
         DO j = istart , iend
            IF ( ILS==1 ) THEN
               x6(j) = x6(j) + Corr(6*imol+4)*Frac_change(4)
               y6(j) = y6(j) + Corr(6*imol+5)*Frac_change(4)
               z6(j) = z6(j) + Corr(6*imol+6)*Frac_change(4)
            ENDIF
            X_cor(j) = x6(j)/Cell_new(1) + xcc(imol)
            Y_cor(j) = y6(j)/Cell_new(2) + ycc(imol)
            Z_cor(j) = z6(j)/Cell_new(3) + zcc(imol)
            sum1 = sum1 + X_cor(j)
            sum2 = sum2 + Y_cor(j)
            sum3 = sum3 + Z_cor(j)
         ENDDO
         istart = istart + NATm(imol)
      ENDDO
!
!     extra_x = sum1/N                     ! deleted, 11/18/08, HLA
!     extra_y = sum2/N                     ! deleted, 11/18/08, HLA
!     extra_z = sum3/N                     ! deleted, 11/18/08, HLA
!
      END SUBROUTINE CORR_XYZ
!
!end-------------------------corr_xyz-----------------------------end
!
      SUBROUTINE CORR_XYZ_ROT(Corr,Frac_change,Cell1,X,Y,Z,Cell_new,    &
     &                        X_cor,Y_cor,Z_cor)
      USE PMIN_MODULE , ONLY:NMOl , NATm , IDX , ISYstem , N , ILS ,    &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , PI , ANG_TO_RADIAN                                               ! 8-19-08
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 , DIMENSION(6) :: Cell1 , Cell_new
      REAL*8 , DIMENSION(NVD) :: Corr
      REAL*8 , DIMENSION(4) :: Frac_change
      REAL*8 , DIMENSION(MATM) :: X , X_cor , Y , Y_cor , Z , Z_cor
      INTENT (IN) Cell1 , Corr , Frac_change , X , Y , Z
      INTENT (INOUT) Cell_new , X_cor , Y_cor , Z_cor
!
! Local variables
!
      REAL*8 :: a , b , c , c2 , c3 , c61 , cc1 , cc2 , cc3 , egymol ,  &
     &          phi , psi , r1 , r2 , s1 , s2 , s3 , ss1 , ss3 , sum1 , &
     &          sum2 , sum3 , sumx , sumy , sumz , sx , sy , sz , theta
      INTEGER :: i , iend , ik , imol , istart , ivr , j , nv , nvr
      REAL*8 , DIMENSION(MATM) :: x3 , x4 , x6 , xo , xx2 , x_temp ,    &
     &                            y3 , y4 , y6 , yo , yy2 , y_temp ,    &
     &                            z3 , z4 , z6 , zo , zz2 , z_temp
      REAL*8 , DIMENSION(NMOl) :: xcc , ycc , zcc
      REAL*8 , DIMENSION(3,3) :: xtm
!
      nvr = NROtbond                           ! 8-19-08
      nv = 6 + 6*NMOl + nvr                    ! 8-19-08
      a = Cell1(1)
      b = Cell1(2)
      c = Cell1(3)
      cc1 = Cell1(4)
      cc2 = Cell1(5)
      cc3 = Cell1(6)
      ss1 = SQRT(1.0-cc1*cc1)
!     ss2 = SQRT(1.0-cc2*cc2)                  ! deleted, 11/18/08, HLA
      ss3 = SQRT(1.0-cc3*cc3)
!
! Now Calculate the orthognal coordinates xx2,yy2,zz2 rotated by Euler angles phi,theta,psi
! Take the coordinates referred to centriod as origin
!
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
!
      istart = 1
      DO imol = 1 , NMOl
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = X(i) + sumx
            sumy = Y(i) + sumy
            sumz = Z(i) + sumz
         ENDDO
         xcc(imol) = sumx/NATm(imol)
         ycc(imol) = sumy/NATm(imol)
         zcc(imol) = sumz/NATm(imol)
         DO i = istart , iend
            xx2(i) = a*(X(i)-xcc(imol))
            yy2(i) = b*(Y(i)-ycc(imol))
            zz2(i) = c*(Z(i)-zcc(imol))
         ENDDO
         DO i = istart , iend
            x3(i) = xx2(i) + yy2(i)*cc3 + zz2(i)*cc2
            y3(i) = yy2(i)*ss3 + zz2(i)*(cc1-cc2*cc3)/ss3
            z3(i) = zz2(i)*SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
         ENDDO
!
         phi = 0.0
         theta = 0.0
         psi = 0.0
         phi = Corr(6*imol+1)*Frac_change(3)
         theta = Corr(6*imol+2)*Frac_change(3)
         psi = Corr(6*imol+3)*Frac_change(3)
!
!----These angles are in radians
         s1 = SIN(phi)
         s2 = SIN(theta)
         s3 = SIN(psi)
         c61 = COS(phi)
         c2 = COS(theta)
         c3 = COS(psi)
!
         xtm(1,1) = c3*c2
         xtm(1,2) = -s3*c2
         xtm(1,3) = s2
         xtm(2,1) = c3*s2*s1 + s3*c61
         xtm(2,2) = c3*c61 - s3*s2*s1
         xtm(2,3) = -c2*s1
         xtm(3,1) = -c3*s2*c61 + s3*s1
         xtm(3,2) = c3*s1 + s3*s2*c61
         xtm(3,3) = c2*c61
!
!----Rotate by phi, theta, psi and translate to xcc+corr(11) etc
         DO j = istart , iend
            x4(j) = x3(j)*xtm(1,1) + y3(j)*xtm(2,1) + z3(j)*xtm(3,1)
            y4(j) = x3(j)*xtm(1,2) + y3(j)*xtm(2,2) + z3(j)*xtm(3,2)
            z4(j) = x3(j)*xtm(1,3) + y3(j)*xtm(2,3) + z3(j)*xtm(3,3)
         ENDDO
!
         IF ( ILS==0 ) THEN
            DO i = istart , iend
               x4(i) = x4(i) + Corr(6*imol+4)*Frac_change(4)
               y4(i) = y4(i) + Corr(6*imol+5)*Frac_change(4)
               z4(i) = z4(i) + Corr(6*imol+6)*Frac_change(4)
            ENDDO
         ENDIF
         istart = istart + NATm(imol)
      ENDDO
!
!----Convert rotated orthogonal coordinates to crystallographic cell in Angs
      Cell_new(1) = a + Corr(1)*Frac_change(1)
      Cell_new(2) = b + Corr(2)*Frac_change(1)
      Cell_new(3) = c + Corr(3)*Frac_change(1)
      Cell_new(4) = Cell1(4) + Corr(4)*Frac_change(2)
      Cell_new(5) = Cell1(5) + Corr(5)*Frac_change(2)
      Cell_new(6) = Cell1(6) + Corr(6)*Frac_change(2)
!
      cc1 = Cell_new(4)
      cc2 = Cell_new(5)
      cc3 = Cell_new(6)
      ss1 = SQRT(1.0-cc1*cc1)
!     ss2 = SQRT(1.0-cc2*cc2)                  ! deleted, 11/18/09, HLA
      ss3 = SQRT(1.0-cc3*cc3)
!
      r1 = (cc1-cc2*cc3)/ss3
      r2 = SQRT(ss1*ss1-((cc2-cc1*cc3)/ss3)**2)
!
      istart = 1
      DO imol = 1 , NMOl
         iend = istart + NATm(imol) - 1
!
         DO j = istart , iend
            z6(j) = z4(j)/r2
            y6(j) = (y4(j)-z6(j)*r1)/ss3
            x6(j) = x4(j) - y6(j)*cc3 - z6(j)*cc2
         ENDDO
!
         DO j = istart , iend
            IF ( ILS==1 ) THEN
               x6(j) = x6(j) + Corr(6*imol+4)*Frac_change(4)
               y6(j) = y6(j) + Corr(6*imol+5)*Frac_change(4)
               z6(j) = z6(j) + Corr(6*imol+6)*Frac_change(4)
            ENDIF
            X_cor(j) = x6(j)/Cell_new(1) + xcc(imol)
            Y_cor(j) = y6(j)/Cell_new(2) + ycc(imol)
            Z_cor(j) = z6(j)/Cell_new(3) + zcc(imol)
            sum1 = sum1 + X_cor(j)
            sum2 = sum2 + Y_cor(j)
            sum3 = sum3 + Z_cor(j)
         ENDDO
         istart = istart + NATm(imol)
      ENDDO
!
! ----------------------------- added following lines on 8-19-08 DU ----------------------------
      IF ( IROt/=0 ) THEN
         DO ik = nv - nvr + 1 , nv
            ivr = ik - (nv-nvr)
            DTHeta = Corr(ik)
            DO i = 1 , N
               x_temp(i) = X_cor(i)
               y_temp(i) = Y_cor(i)
               z_temp(i) = Z_cor(i)
            ENDDO
            CALL BEND_BOND_NEW(ivr,N,Cell_new,x_temp,y_temp,z_temp,xo,  &
     &                         yo,zo,egymol)
            CALL FRACTIONAL_COD(N,Cell_new,xo,yo,zo,X_cor,Y_cor,Z_cor)
            IF ( I_Nv==1 .OR. I_Nv>(nv-nvr) ) THEN                           ! altered, 11/16/08, HLA
               DEL_egymol = egymol - EGYmol0                                 ! altered, 11/16/08, HLA
                                                                      !      ! altered, 11/16/08, HLA
               WRITE (61,99001) ivr , Corr(ik) , Corr(ik)*57.29578 ,    &
     &                          EGYmol0 , egymol , DEL_egymol                ! altered, 11/16/08, HLA
99001          FORMAT ('intra-E for each bond:',I4,F16.6,F9.3,3F16.6)        ! altered, 11/16/08, HLA
            ENDIF                                                            ! altered, 11/16/08, HLA
         ENDDO
!
         sx = 0.0
         sy = 0.0
         sz = 0.0
99002    FORMAT (6F9.5)
         DO ik = 1 , N
            sx = sx + X_cor(ik)
            sy = sy + Y_cor(ik)
            sz = sz + Z_cor(ik)
         ENDDO
         sx = sx/N
         sy = sy/N
         sz = sz/N
         ATOm(N+1) = 'XTRA  '
      ELSE
         DEL_egymol = 0.0
      ENDIF
99003 FORMAT (a6,i3,18x,3F9.6)
      END SUBROUTINE CORR_XYZ_ROT
!
!end-----------------------subroutine corr_xyz_rot----------------------------end
!
! The subroutine POT_E below calculates the total lattice energy for atoms within one asymmetric unit
! If the asymmetric unit contains one molecule, the calculated lattice energy is per molecule.
! The total lattice energy has three parts given by 6-exp formula. All intermolecular interactions within
! a radius of dmax (6 to 10 Angstrom) are considered around the atoms of the asymmetric unit. The atoms
! surrounding the molecule are obtained by the space group symmetry.
! The Coulomb and van der Waals' terms are calculated using the reciprocal space summation.
!
      SUBROUTINE POT_E(Cell,X,Y,Z,Ec,Ev,Er,E_nmol,Pe)
      USE PMIN_MODULE , ONLY:NTX_min , NTY_min , NTZ_min , NTX_max ,    &
     &    NTY_max , NTZ_max , NMOl , NATm , IDX , ISYstem , N , ILS ,   &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , STHl , SYM , Q , A1 , B1 , C1 , CK , WT_mol ,   &
     &    DMAx , DDMax , GNI , CE12 , VE12 , K1 , N11 , IHKl , IRIj ,   &
     &    KK , I_Cross , A_Cross , B_Cross , C_Cross , ANAme1 , ANAme2 ,&
     &    A12 , B12 , C12 , PI , TWOPI , ANG_TO_RADIAN, dd_cross            ! 5-4-10 DU 
      USE RSS1_MODULE , ONLY:REFcode, Icall_pot_e                           ! 2-3-09 
      USE RSSPM3_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 :: Ec , Er , Ev , E_nmol , Pe
      REAL*8 , DIMENSION(6) :: Cell
      REAL*8 , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN) Cell , X , Y , Z
      INTENT (INOUT) Ec , Er , Ev , E_nmol , Pe
!
! Local variables
!
      REAL*8 :: a , alpha , asqrt , asquare , b , beta , bsq , c , c11 ,&
     &          c22 , c33 , cb1 , cb2 , ce1 , ce11 , ce2 , ce22 , ce3 , &
     &          ck11 , ck22 , ck33 , ck44 , coulomb_en , d , distmax ,  &
     &          distmin , dx , dy , dz , e11 , e12 , e2 , e22 , e3 ,    &
     &          e33 , ec_nmol , er_nmol , ev_nmol , gamma , gi , r1 ,   &
     &          r2 , r3 , r4 , r5 , r6 , repulse_en , rij , sah , sin1 ,&
     &          sin2 , sin3 , sins , sqrtpi , sum123 , sum_en , v11 ,   &
     &          v12 , vander_en , ve , ve1 , ve11 , ve2 , ve22 , ve3 ,  &
     &          ve4 , vol , x11 , x2 , x22 , x4 , xni , y11 , y2 , y22 ,&
     &          y4 , z11 , z2 , z22 , z4
      REAL*8 , DIMENSION(6) :: ah , rh
      REAL*8 :: DERFC
      REAL*8 :: DSQRT
      REAL*8 , DIMENSION(1000) :: fc , fcb , fv , fvb , qh
      INTEGER :: i , i1 , iend1 , iend2 , ig , ih , ihmax , ihmin , ik ,&
     &           ik3 , ik4 , il , imol1 , imol2 , istart1 , istart2 ,   &
     &           j , jh , jk , jl , k , k11 , k111 , k2 , k22 , k3 ,    &
     &           k33 , kmax , kmin , ktx , kty , ktz , lmax , lmin ,    &
     &           n12 , natom , ne1 , nrij , ntotal , ntx , ntx2 , nty , &
     &           nty2 , ntz , ntz2
      INTEGER :: INT
      REAL*8 :: zz
!
!     WRITE(20,*)                                                           ! 2-26-09
      e22 = 0.0                                                             ! 4-18-08
      distmin = 999.0
      distmax = -999.0
      short_dd =0.0                                                         ! 2-26-09
      DO i = 1 , 1000
         qh(i) = 0.0
      ENDDO
      zz = NSYm
      natom = N
      DO i = 1 , 6
         ah(i) = Cell(i)
      ENDDO
!
      IF ( ah(4)>10.0 ) THEN      ! angles are degrees or cosines?
         ah(4) = COS(ah(4)*ANG_TO_RADIAN)                                   ! 5/13/08, HLA,added
         ah(5) = COS(ah(5)*ANG_TO_RADIAN)                                   ! 5/13/08, HLA,added
         ah(6) = COS(ah(6)*ANG_TO_RADIAN)                                   ! 5/13/08, HLA,added
      ENDIF
      c11 = ah(4)                 ! c11-c33 are cosines
      c22 = ah(5)
      c33 = ah(6)
      a = ah(1)
      b = ah(2)
      c = ah(3)
      alpha = ACOSD(ah(4))        ! alpha, beta,
      beta = ACOSD(ah(5))         ! gamma are degrees
      gamma = ACOSD(ah(6))
      sah = 0.5*(alpha+beta+gamma)
      sins = SIN(ANG_TO_RADIAN*sah)                                         ! 5/13/08, HLA, altered
      sin1 = SIN(ANG_TO_RADIAN*(sah-alpha))                                 ! 5/13/08, HLA, altered
      sin2 = SIN(ANG_TO_RADIAN*(sah-beta))                                  ! 5/13/08, HLA, altered
      sin3 = SIN(ANG_TO_RADIAN*(sah-gamma))                                 ! 5/13/08, HLA, altered
      vol = 2.0*ah(1)*ah(2)*ah(3)*DSQRT(sins*sin1*sin2*sin3)
!
      CALL RECIP(ah,rh)
!
      r1 = rh(1)**2
      r2 = rh(2)**2
      r3 = rh(3)**2
      r4 = 2.0*rh(1)*rh(2)*rh(6)
      r5 = 2.0*rh(1)*rh(3)*rh(5)
      r6 = 2.0*rh(3)*rh(2)*rh(4)
!
IF_NENtry_LE_1:      IF ( NENtry<=1 ) THEN
!----Excludes calculation of hkl's & gni except on 1st entry
!
         ihmax = INT(2.0*STHl/rh(1))                ! h, k, l maxes for reciprocal lattice
         kmax = INT(2.0*STHl/rh(2))                 ! summations
         lmax = INT(2.0*STHl/rh(3))
         ihmin = -ihmax
         kmin = -kmax
         lmin = 0
         gi = 2.0
         IF ( ABS(rh(4))<=0.001 .AND. ABS(rh(6))<=0.001 ) THEN
            kmin = 0
            gi = 2.0*gi
         ENDIF
         IF (ABS(rh(5))<=0.000001 ) THEN            ! 2-15-09 DU      
            ihmin = 0
            gi = 2.0*gi
         ENDIF 
!
!----Determine set of hkl values for given value of sthl
!
         K1 = 0
         DO ih = ihmin , ihmax
            DO ik = kmin , kmax
               DO il = lmin , lmax
                  xni = 1.0
                  IF ( ih/=0 .OR. ik/=0 .OR. il/=0 ) THEN
                     IF ( ih==0 .AND. ihmin==0 ) xni = 0.5*xni
                     IF ( ik==0 .AND. kmin==0 ) xni = 0.5*xni
                     IF ( il==0 ) xni = 0.5*xni
                     d = SQRT((r1*ih**2)+(r2*ik**2)+(r3*il**2) &
     &                   +(r5*ih*il)+(r4*ih*ik)+(r6*ik*il))
                     IF ( d<=STHl ) THEN
                        K1 = K1 + 1
                        IHKl(K1,1) = ih
                        IHKl(K1,2) = ik
                        IHKl(K1,3) = il
                        qh(K1) = d
                        GNI(K1) = gi*xni
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      END IF IF_NENtry_LE_1 
      if ( Icall_pot_e ==  1) then
      write(61,"('STHl,rh : ',/7F12.6)") STHl, rh
      write(61,"('ihmin , ihmax, kmin , kmax, lmin , lmax, K1 :'/6I6,1X,I6)")  &
    &             ihmin , ihmax, kmin , kmax, lmin , lmax, K1
      end if

!
DO_LOOP_k:      DO k = 1 , K1
         jh = IHKl(k,1)
         jk = IHKl(k,2)
         jl = IHKl(k,3)
         IF ( NENtry>1 ) THEN
            d = r1*jh**2 + r2*jk**2 + r3*jl**2 + r5*jh*jl + r4*jh*jk +  &
     &          r6*jk*jl
            qh(k) = SQRT(d)
         ENDIF
         fc(k) = 0.0
         fcb(k) = 0.0
         fv(k) = 0.0
         fvb(k) = 0.0
         DO ig = 1 , NSYm
            DO j = 1 , N                                                  ! loop over all atoms
               x2 = SYM(ig,1) + X(j)*SYM(ig,2) + Y(j)*SYM(ig,3) + Z(j)  &
     &              *SYM(ig,4)
               y2 = SYM(ig,5) + X(j)*SYM(ig,6) + Y(j)*SYM(ig,7) + Z(j)  &
     &              *SYM(ig,8)
               z2 = SYM(ig,9) + X(j)*SYM(ig,10) + Y(j)*SYM(ig,11) + Z(j)&
     &              *SYM(ig,12)
               fc(k) = fc(k) + Q(j)*COS(TWOPI*(jh*x2+jk*y2+jl*z2))
               fv(k) = fv(k) + A1(j)*COS(TWOPI*(jh*x2+jk*y2+jl*z2))
               fcb(k) = fcb(k) + Q(j)*SIN(TWOPI*(jh*x2+jk*y2+jl*z2))
               fvb(k) = fvb(k) + A1(j)*SIN(TWOPI*(jh*x2+jk*y2+jl*z2))
            ENDDO
         ENDDO
         fc(k) = SQRT(fc(k)**2+fcb(k)**2)
         fv(k) = SQRT(fv(k)**2+fvb(k)**2)
      ENDDO DO_LOOP_k
!
!----Calculation of structure factors finished
!
!----Calculate the reciprocal sum for Coulomb and Van der Waals energies
!
      ce2 = 0.0
      ve2 = 0.0
      sqrtpi = DSQRT(PI)
      zz = NSYm
      ck11 = 1.0/(TWOPI*vol*zz)                     ! 5/13/08, HLA added
      ck22 = PI**4.5/(3*vol*zz)
      ck33 = (PI**3.0)*(CK**6)/12.0
      ck44 = (PI**3.0)*(CK**3)*zz/(6.0*vol)
!
      DO i = 1 , K1
         bsq = PI*qh(i)*qh(i)/(CK*CK)
         cb1 = SQRT(bsq)   ! Check if taking sqare root for getting constant B is valid or not
         cb2 = sqrtpi*DERFC(cb1) + (1.0/(2.0*bsq*cb1)-1.0/cb1)*EXP(-bsq)
         ce22 = 332.17*fc(i)*fc(i)*GNI(i)*(EXP(-bsq))/(qh(i)*qh(i))
         ce22 = ck11*ce22
         ce2 = ce2 + ce22
         ve22 = GNI(i)*fv(i)*fv(i)*(qh(i)**3)*cb2
         ve22 = -ck22*ve22
         ve2 = ve2 + ve22
      ENDDO
!
!----No special position assumed presently
      ce3 = 0.0
      ve3 = 0.0
      ve4 = 0.0
      DO i = 1 , N
         ce3 = ce3 + Q(i)*Q(i)
         ve3 = ve3 + A1(i)*A1(i)      ! No special position of atom assumed n(i) = 1.0
         ve4 = ve4 + A1(i)
      ENDDO
      ce3 = -ce3*CK*332.17
      ve3 = ck33*ve3
      ve4 = -ck44*ve4*ve4
!
      e2 = 0.0
      e3 = 0.0
      n12 = 0
      ve11 = 0.0
      ce11 = 0.0
      ntotal = 0
      sum123 = 0.0
IF_NENtry_GT_1_1:      IF ( NENtry>1 ) THEN
!
         nrij = 0
!
         CE12 = 0.0
         VE12 = 0.0
         ce11 = 0.0
         ve11 = 0.0
         e3 = 0.0
         n12 = 0
!
IF_NENtry_GT_1_2:      IF ( NENtry>1 ) THEN
!
            CE12 = 0.0
            VE12 = 0.0
            ce11 = 0.0
            ve11 = 0.0
            e3 = 0.0
!
!----Use n11 interactions only
DO_loop_ik: DO ik = 1 , N11
               i = IRIj(ik,1)
               i1 = IRIj(ik,2)
               j = IRIj(ik,3)
               ktx = IRIj(ik,4)
               kty = IRIj(ik,5)
               ktz = IRIj(ik,6)
               x11 = X(i)*a
               y11 = Y(i)*b
               z11 = Z(i)*c
               x4 = a*(X(i1)*SYM(j,2)+Y(i1)*SYM(j,3)+Z(i1)*SYM(j,4)     &
     &              +SYM(j,1)+ktx)
               dx = x11 - x4
               y4 = b*(X(i1)*SYM(j,6)+Y(i1)*SYM(j,7)+Z(i1)*SYM(j,8)     &
     &              +SYM(j,5)+kty)
               dy = y11 - y4
               z4 = c*(X(i1)*SYM(j,10)+Y(i1)*SYM(j,11)+Z(i1)*SYM(j,12)  &
     &              +SYM(j,9)+ktz)
               dz = z11 - z4
!
               rij = DSQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+2*dy*dz*c11+   &
     &               2*dx*dz*c22)
!
!      Now calculate the potential energy
!
               asquare = PI*CK*CK*rij*rij
               asqrt = DSQRT(asquare)
!
               IF ( j==1 .AND. ktx==0 .AND. kty==0 .AND. ktz==0 ) THEN
                 e12 = 332.17*Q(i)*Q(i1)*(DERFC(asqrt)-1.0)/rij
                 IF ( I_Cross/=0 ) THEN                                                    ! 10-16-09 DU
                  DO KK =1,I_Cross                                                         ! 10-16-09 DU
                     IF ((ATOm(i)==ANAme1(KK) .AND. ATOm(i1)==ANAme2(KK).AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)==ANAme1(KK).AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1)==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK))) THEN       ! 5-4-10 DU    
                        if (ICYcle==1 .AND. NCY<=2 ) &                                     ! 10-16-09&
     &                  write(20, '(A6,A6,3F9.4,1X,F6.4)') ATOm(i),ATOm(i1),A1(i),A1(i1),A12,rij           !10-16-09 DU
                        A12 = A_Cross(KK)                                                  ! 10-16-09 DU
                        B12 = B_Cross(KK)                                                  ! 10-16-09 DU
                        C12 = C_Cross(KK)                                                  ! 10-16-09 DU
                        GO TO 99029                                                        ! found cross term, exit DO Loop I_croos
                     ELSE                                                                  ! 10-16-09 DU
                        A12 = A1(i)*A1(i1)                                                 ! 10-16-09 DU
                        B12 = B1(i)*B1(i1)                                                 ! 10-16-09 DU
                        C12 = C1(i) + C1(i1)                                               ! 10-16-09 DU
                     ENDIF                                                                 ! 10-16-09 DU
                  ENDDO                                                                    ! 10-16-09 DU
                 ELSE                                                                      ! 10-16-09 DU
                  A12 = A1(i)*A1(i1)                                                       ! 10-16-09 DU
                  B12 = B1(i)*B1(i1)                                                       ! 10-16-09 DU
                  C12 = C1(i) + C1(i1)                                                     ! 10-16-09 DU
                 ENDIF
99029             v12 = -A12                                            &
     &                  *((1.0+asquare+0.5*asquare*asquare)*EXP         &
     &                  (-asquare)-1.0)/(rij**6)
                  IF ( i==i1 ) THEN
                     e12 = 0.5*e12
                     v12 = 0.5*v12
                  ENDIF
                  CE12 = CE12 + e12
                  VE12 = VE12 + v12
                  CYCLE
               ENDIF
               e11 = 332.17*Q(i)*Q(i1)*DERFC(asqrt)/rij
!----Cross-term coefficients ??
               IF ( I_Cross/=0 ) THEN                                                      ! 5-1-08 DU
                  DO KK =1,I_Cross                                                         ! 5-1-08 DU
!                  IF ( rij<3.0 .AND. IMOde==0 .AND. ICYcle==1 .AND.       &
!    &              NCY<=2 ) THEN                                                          
!                   WRITE (20,'(4A6,1X)') ATOm(i), ATOm(i1), ANAme1(KK), ANAme2(KK)          
!                  END IF                                                                  ! 10-13-09
                     IF ((ATOm(i)==ANAme1(KK) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)==ANAme1(KK) .AND. Rij <= dd_cross(KK)) .OR.  & 
     &                   (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1)==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK))) THEN       ! 5-4-10 DU     
                        A12 = A_Cross(KK)                                                  ! 5-1-08 DU
                        B12 = B_Cross(KK)                                                  ! 5-1-08 DU
                        C12 = C_Cross(KK)                                                  ! 5-1-08 DU
!                       write(20,"(6H######,F6.4, 2F16.8,6H######)") Rij, A12, B12                          ! 5-4-10
                        GO TO 99009                                                        ! found cross term, exit DO Loop I_croos 
                     ELSE                                                                  ! 5-1-08 DU
                        A12 = A1(i)*A1(i1)                                                 ! 5-1-08 DU
                        B12 = B1(i)*B1(i1)                                                 ! 5-1-08 DU
                        C12 = C1(i) + C1(i1)                                               ! 5-1-08 DU
                     ENDIF                                                                 ! 5-1-08 DU
                  ENDDO                                                                    ! 5-1-08 DU
               ELSE                                                                        ! 5-1-08 DU
                  A12 = A1(i)*A1(i1)                                                       ! 5-1-08 DU
                  B12 = B1(i)*B1(i1)                                                       ! 5-1-08 DU
                  C12 = C1(i) + C1(i1)                                                     ! 5-1-08 DU
               ENDIF                                                                       ! 5-1-08 DU
99009          v11 = A12*((1.0+asquare+0.5*asquare*asquare)            &                   ! 10-13-09
     &               *EXP(-asquare))/(rij**6)                                              ! 5-1-08 DU
               e33 = B12*EXP(-C12*rij)                                                     ! 5-1-08 DU
!
               IF ( rij<3.0 .AND. IMOde==0 .AND. ICYcle==1 .AND.       &
     &              NCY<=2 ) THEN                                                          ! 5-1-08 DU
                  WRITE (20,99001) NENtry , ATOm(i) , ATOm(i1) , rij , &
     &                             A1(i) , A1(i1) , B1(i) , B1(i1) ,   &
     &                             A12 , B12 , e11 , e22 , e33,I_Cross                     ! 10-13-09 DU
99001             FORMAT (I6,1X,A6,'---',A6,2X,f6.4,4F9.4,1X,2F16.8,   &
     &                    3F12.6,' #3',I3)                                                 ! 10-13-09 DU
               ENDIF                                                                       ! 5-1-08 DU
               IF ( i==i1 ) THEN
                  e11 = 0.5*e11
                  v11 = 0.5*v11
                  e22 = 0.5*e22
                  e33 = 0.5*e33
               ENDIF
               ce11 = ce11 + e11
               ve11 = ve11 - v11
               e3 = e33 + e3
            ENDDO DO_loop_ik
         ELSE                        ! IF NENtry_GT_1.2                                    ! 10-16-09
         GO TO 2211                  ! 10-31-09 SMP
!
!----Loop over all the atoms
            DO i = 1 , natom
               x11 = X(i)*a
               y11 = Y(i)*b
               z11 = Z(i)*c
!----Over all the second atoms
               DO i1 = i , natom
                                ! Check for index i
!----Over all the symmetries
                  DO j = 1 , NSYm
!----Over all the translations
                     DO ktx = NTX_min , NTX_max
                        x4 = a*(X(i1)*SYM(j,2)+Y(i1)*SYM(j,3)+Z(i1)     &
     &                       *SYM(j,4)+SYM(j,1)+ktx)
                        dx = x11 - x4
                        DO kty = NTY_min , NTY_max
                           y4 = b*(X(i1)*SYM(j,6)+Y(i1)*SYM(j,7)+Z(i1)  &
     &                          *SYM(j,8)+SYM(j,5)+kty)
                           dy = y11 - y4
                           DO ktz = NTZ_min , NTZ_max
                              z4 = c*(X(i1)*SYM(j,10)+Y(i1)*SYM(j,11)   &
     &                             +Z(i1)*SYM(j,12)+SYM(j,9)+ktz)
                              dz = z11 - z4
                              rij = DSQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+&
     &                              2*dy*dz*c11+2*dx*dz*c22)
                              IF ( rij>0.00001 ) THEN
                                 nrij = nrij + 1
!
!---Calculate the potential energy
!
                                 asquare = PI*CK*CK*rij*rij
                                 asqrt = DSQRT(asquare)
                                 IF ( j==1 .AND. ktx==0 .AND.           &
     &                                kty==0 .AND. ktz==0 ) THEN
                                    e12 = 332.17*Q(i)*Q(i1)             &
     &                                 *(DERFC(asqrt)-1.0)/rij
                                 IF ( I_Cross/=0 ) THEN                                           ! 10-16-09 DU
                                  DO KK =1,I_Cross                                                ! 10-16-09 DU
                                   IF ( (ATOm(i)==ANAme1(KK) .AND. ATOm(i1)==ANAme2(KK).AND. Rij <= dd_cross(KK)) .OR.  &
     &                                 (ATOm(i)==ANAme2(KK)  .AND. ATOm(i1)==ANAme1(KK).AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1)==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN       ! 5-4-10 DU
                                        if (ICYcle==1 .AND. NCY<=2 ) &
     &                                  write(20, '(A6,A6,3F9.4,1X,F6.4)') ATOm(i),ATOm(i1),A1(i),A1(i1),A12,rij           !10-16-09 DU
                                    A12 = A_Cross(KK)                                             ! 10-16-09 DU
                                    B12 = B_Cross(KK)                                             ! 10-16-09 DU
                                    C12 = C_Cross(KK)                                             ! 10-16-09 DU
                                    GO TO 99039                                                   ! found cross term, exit DO Loop I_croos
                                   ELSE                                                           ! 10-16-09 DU
                                    A12 = A1(i)*A1(i1)                                            ! 10-16-09 DU
                                    B12 = B1(i)*B1(i1)                                            ! 10-16-09 DU
                                    C12 = C1(i) + C1(i1)                                          ! 10-16-09 DU
                                   ENDIF                                                          ! 10-16-09 DU
                                  ENDDO                                                           ! 10-16-09 DU
                                 ELSE                                                             ! 10-16-09 DU
                                    A12 = A1(i)*A1(i1)                                            ! 10-16-09 DU
                                    B12 = B1(i)*B1(i1)                                            ! 10-16-09 DU
                                    C12 = C1(i) + C1(i1)                                          ! 10-16-09 DU
                                 ENDIF
99039                               v12 = -A12                          &
     &                                 *((1.0+asquare+0.5*asquare*      &
     &                                 asquare)*EXP(-asquare)-1.0)      &
     &                                 /(rij**6)
                                    IF ( i==i1 ) THEN
                                       e12 = 0.5*e12
                                       v12 = 0.5*v12
                                    ENDIF
                                    CE12 = CE12 + e12
                                    VE12 = VE12 + v12
                                    CYCLE
                                 ENDIF
!
                                 IF ( rij<=DMAx ) THEN
!
                                    e11 = 332.17*Q(i)*Q(i1)*DERFC(asqrt)/rij
!----Use cross-term coefficients ??
                                    IF ( I_Cross/=0 ) THEN                                 ! 5-1-08 DU
                                       DO KK = 1 , I_Cross                                 ! 5-1-08 DU
                                         IF ( (ATOm(i)==ANAme1(KK) &
     &                                      .AND. ATOm(i1)         &
     &                                      ==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.     &
     &                                      (ATOm(i)==ANAme2(KK)   &
     &                                      .AND. ATOm(i1)         &
     &                                      ==ANAme1(KK) .AND. Rij <= dd_cross(KK)) .OR.     &
     &                   (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1)==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN       ! 5-4-10 DU
                                         A12 = A_Cross(KK)                                 ! 5-1-08 DU
                                         B12 = B_Cross(KK)                                 ! 5-1-08 DU
                                         C12 = C_Cross(KK)                                 ! 5-1-08 DU
                                         GO TO  99922                                      ! found cross term, exit DO Loop I_croos
                                         ELSE                                              ! 5-1-08 DU
                                         A12 = A1(i)*A1(i1)                                ! 5-1-08 DU
                                         B12 = B1(i)*B1(i1)                                ! 5-1-08 DU
                                         C12 = C1(i) + C1(i1)                              ! 5-1-08 DU
                                         ENDIF                                             ! 5-1-08 DU
                                       ENDDO                                               ! 5-1-08 DU
                                    ELSE                                                   ! 5-1-08 DU
                                       A12 = A1(i)*A1(i1)                                  ! 5-1-08 DU
                                       B12 = B1(i)*B1(i1)                                  ! 5-1-08 DU
                                       C12 = C1(i) + C1(i1)                                ! 5-1-08 DU
                                    ENDIF                                                  ! 5-1-08 DU
99922                               v11 = A12*                           &                 ! 10-13-09
     &                                 ((1.0+asquare+0.5*asquare*asquare &
     &                                 )*EXP(-asquare))/(rij**6)                           ! 5-1-08 DU
                                    e33 = B12*EXP(-C12*rij)                                ! 5-1-08 DU
!
                                    IF ( rij<3.0 .AND. IMOde==0 .AND.    &
     &                                 ICYcle==1 .AND. NCY<=2 ) THEN                       ! 5-1-08 DU
                                       WRITE (20,99002) NENtry , ATOm(i) &
     &                                    , ATOm(i1) , rij , A1(i) ,     &
     &                                    A1(i1) , B1(i) , B1(i1) ,      &
     &                                    A12 , B12 , e11 , e22 , e33
99002                                  FORMAT (I6,1X,A6,'---',A6,2X,     &
     &                                    f6.4,4F9.4,1X,2F16.8,3F12.6,   &
     &                                    ' #2')                                           ! 5-1-08 DU
                                    ENDIF                                                  ! 5-1-08 DU
!
                                    IF ( i==i1 ) THEN
                                       e11 = 0.5*e11
                                       v11 = 0.5*v11
                                       e22 = 0.5*e22
                                       e33 = 0.5*e33
                                    ENDIF
                                    ce11 = ce11 + e11
                                    ve11 = ve11 - v11
                                    e3 = e33 + e3
                                 ENDIF
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
!
            ENDDO
         ENDIF   IF_NENtry_GT_1_2
2211     continue     ! 10-31-09 SMP                 
      ELSE            ! IF NENtry_GT_1_1   10-16-09 DU    
!
!----For translations along x, y, z
         N11 = 0
         VE12 = 0.0
         CE12 = 0.0
         ntx = INT(DMAx/a) + 2
         nty = INT(DMAx/b) + 2
         ntz = INT(DMAx/c) + 2
         NTX_min = 0
         NTX_max = 0
         NTY_min = 0
         NTY_max = 0
         NTZ_min = 0
         NTZ_max = 0
         ntx2 = 2*ntx + 1
         nty2 = 2*nty + 1
         ntz2 = 2*ntz + 1
!
         IF (IMOde==0 .AND.ICYcle==1 .AND. NCY<=2 ) write(20,99033)                       ! 2-3-09 
99033    format('   Names        R      Q(i)    Q(i1)   A(i)    A(i1)  B(i)     B(i1)', &
     &          '   Qi*Qi1 DERFC(asqrt) Ai*Ai1    Bi*Bi1   Ci+Ci1  CNB     VDW     REP'/)
!----Over all of the atoms
loop_i:  DO i = 1 , natom   ! over all of the atoms
            x11 = X(i)*a
            y11 = Y(i)*b
            z11 = Z(i)*c
!----Over all the second atoms
loop_i1:  DO i1 = i , natom  ! Check for index i
!----Over all the symmetries
loop_j:    DO j = 1 , NSYm
!----Over all the translations
loop_k111:  DO k111 = 1 , ntx2
                     x4 = a*(X(i1)*SYM(j,2)+Y(i1)*SYM(j,3)+Z(i1)        &
     &                    *SYM(j,4)+SYM(j,1)+(ntx+1-k111))
                     dx = x11 - x4
loop_k2:             DO k2 = 1 , nty2
                        y4 = b*(X(i1)*SYM(j,6)+Y(i1)*SYM(j,7)+Z(i1)     &
     &                       *SYM(j,8)+SYM(j,5)+(nty+1-k2))
                        dy = y11 - y4
loop_k3:                DO k3 = 1 , ntz2
                           z4 = c*(X(i1)*SYM(j,10)+Y(i1)*SYM(j,11)+Z(i1)&
     &                          *SYM(j,12)+SYM(j,9)+(ntz+1-k3))
                           dz = z11 - z4
!
                           rij = DSQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+   &
     &                           2*dy*dz*c11+2*dx*dz*c22)
                           ntotal = ntotal + 1
                           IF ( rij>distmax ) distmax = rij
                           IF ( rij>0.000001 ) THEN
!
                              k11 = 5 + (ntx+1-k111)
                              k22 = 5 + (nty+1-k2)
                              k33 = 5 + (ntz+1-k3)
 
                              asquare = PI*CK*CK*rij*rij
                              asqrt = DSQRT(asquare)
!
                              IF ( j==1 .AND. k11==5 .AND. k22==5 .AND. &
     &                             k33==5 ) THEN
                                 e12 = 332.17*Q(i)*Q(i1)                &
     &                                 *(DERFC(asqrt)-1.0)/rij
                               IF ( I_Cross/=0 ) THEN                                                    ! 10-16-09 DU
                                DO KK =1,I_Cross                                                         ! 10-16-09 DU
                                IF ( (ATOm(i)==ANAme1(KK) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR.  &
     &                         (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)==ANAme1(KK) .AND. Rij <= dd_cross(KK))      .OR.  &
     &                   (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK))  .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1)==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN       ! 5-4-10 DU
                                if (ICYcle==1 .AND. NCY<=2 ) &
     &                          write(20, '(A6,A6,3F9.4,1X,F6.4)') ATOm(i),ATOm(i1),A1(i),A1(i1),A12,rij           !10-16-09 DU
                                A12 = A_Cross(KK)                                                  ! 10-16-09 DU
                                B12 = B_Cross(KK)                                                  ! 10-16-09 DU
                                C12 = C_Cross(KK)                                                  ! 10-16-09 DU
                                GO TO 99049                                                        ! found cross term, exit DO Loop I_croos
                                ELSE                                                               ! 10-16-09 DU
                                A12 = A1(i)*A1(i1)                                                 ! 10-16-09 DU
                                B12 = B1(i)*B1(i1)                                                 ! 10-16-09 DU
                                C12 = C1(i) + C1(i1)                                               ! 10-16-09 DU
                                ENDIF                                                              ! 10-16-09 DU
                                ENDDO                                                              ! 10-16-09 DU
                               ELSE                                                                ! 10-16-09 DU
                                A12 = A1(i)*A1(i1)                                                 ! 10-16-09 DU
                                B12 = B1(i)*B1(i1)                                                 ! 10-16-09 DU
                                C12 = C1(i) + C1(i1)                                               ! 10-16-09 DU
                               ENDIF
99049                           v12 = -A12                              &
     &                                 *((1.0+asquare+0.5*asquare*      &
     &                                 asquare)*EXP(-asquare)-1.0)      &
     &                                 /(rij**6)
                                 IF ( i==i1 ) THEN
                                    e12 = 0.5*e12
                                    v12 = 0.5*v12
                                 ENDIF
                                 CE12 = CE12 + e12
                                 VE12 = VE12 + v12
                                 n12 = n12 + 1
!
                                 N11 = N11 + 1
                                 IRIj(N11,1) = i
                                 IRIj(N11,2) = i1
                                 IRIj(N11,3) = j
                                 IRIj(N11,4) = ntx + 1 - k111
                                 IRIj(N11,5) = nty + 1 - k2
                                 IRIj(N11,6) = ntz + 1 - k3
                                 CYCLE
                              ENDIF
!
                              IF ( rij>0.000001 .AND. rij<=DMAx ) THEN
                                 IF ( NENtry==1 ) THEN
                                    ktx = ntx + 1 - k111
                                    kty = nty + 1 - k2
                                    ktz = ntz + 1 - k3
                                    N11 = N11 + 1
                                    IRIj(N11,1) = i
                                    IRIj(N11,2) = i1
                                    IRIj(N11,3) = j
                                    IRIj(N11,4) = ntx + 1 - k111
                                    IRIj(N11,5) = nty + 1 - k2
                                    IRIj(N11,6) = ntz + 1 - k3
                                    IF ( NTX_min>ktx ) NTX_min = ktx
                                    IF ( NTX_max<ktx ) NTX_max = ktx
                                    IF ( NTY_min>kty ) NTY_min = kty
                                    IF ( NTY_max<kty ) NTY_max = kty
                                    IF ( NTZ_min>ktz ) NTZ_min = ktz
                                    IF ( NTZ_max<ktz ) NTZ_max = ktz
                                 ENDIF
!----Minimum contact distance is distmin
                                 IF ( rij<10.0 .AND. rij<distmin ) distmin = rij
!
!----Calculate the potential energy
                                 e22 = 0.0                                ! 4/29/08, HLA, e22 initialized
                                 e11 = 332.17*Q(i)*Q(i1)*DERFC(asqrt)/rij
!
! cross-term coefs?? 
                                 IF (I_Cross/=0) THEN
cross_do:                            DO KK = 1,I_Cross                   ! 5-1-08 DU
                                        IF ( (ATOm(i)==ANAme1(KK) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK)) .OR. &
                                             (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)==ANAme1(KK) .AND. Rij <= dd_cross(KK)) .OR. &
     &                   (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1)==ANAme2(KK) .AND. Rij <= dd_cross(KK))  .OR.  &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1)==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN       ! 5-4-10 DU
                                           A12 = A_Cross(KK)                ! 5-1-08 DU
                                           B12 = B_Cross(KK)                ! 5-1-08 DU
                                           C12 = C_Cross(KK)                ! 5-1-08 DU
                                          go to 5684                        ! found cross term, exit DO Loop I_croos
                                       ELSE                                 ! 5-1-08 DU
                                          A12 = A1(i)*A1(i1)  ! no cross terms for i and i1
                                          B12 = B1(i)*B1(i1)                ! 5-1-08 DU
                                          C12 = C1(i) + C1(i1)              ! 5-1-08 DU
                                       ENDIF                                ! 5-1-08 DU
                                    ENDDO cross_do                          ! 5-1-08 DU
                                 ELSE                                       ! 5-1-08 DU
                                    A12 = A1(i)*A1(i1) ! no cross-terms for anything      
                                    B12 = B1(i)*B1(i1)                    ! 5-1-08 DU
                                    C12 = C1(i) + C1(i1)                  ! 5-1-08 DU
                                 ENDIF                                    ! 5-1-08 DU
5684                             v11 = A12*((1.0+asquare+0.5*asquare*asquare) & ! HLA...7/26/09
     &                                 *EXP(-asquare))/(rij**6)           ! 5-1-08 DU
                                 e22 = -A12/(rij**6)                      ! 5-1-08 DU
                                 e33 = B12*EXP(-C12*rij)                  ! 5-1-08 DU
!----End of calculation of three terms in potential
                                 IF ( rij<3.0 .AND. IMOde==0 .AND.      & ! 2-3-09 DU
     &                                ICYcle==1 .AND. NCY<=2 ) THEN                      
                                    WRITE (20,99003) NENtry ,           &
     &                                 ATOm(i) , ATOm(i1) ,   &
     &                                 rij , Q(i),Q(i1),A1(i) , A1(i1) , B1(i) ,                &     
     &                                 B1(i1) , Q(i)*Q(i1),DERFC(asqrt),A12 , B12 , C12 , e11 , &
     &                                 e22 , e33                                         
99003                               FORMAT (I2,1X,2A5,1X,f6.4,4F8.4,2F9.4,                      &
     &                                 2f8.4,2F12.4,f5.2,3F8.4,' #1')                     
                                 ENDIF                                    ! 2-3-09 DU
                                    IF ( rij<2.4 .AND. IMOde==0 .AND. ICYcle==1 .AND. ncy_rss_end > 0 ) THEN    ! 2-26-09 DU
!                                   IF ( rij<1.2 .AND. IMOde==0 .AND. ICYcle==1 .AND. ncy_rss_end > 0 ) THEN
                                     if ( ATOm(i)(1:1) == 'H' .AND. ATOm(i1)(1:1) == 'H' .AND. rij > 1.3) go to 99335  ! 5-4-10 
                                      WRITE (22,99003) ncy ,                                    &
     &                                 ATOm(i) , ATOm(i1) ,                                     &
     &                                 rij , Q(i),Q(il),A1(i) , A1(i1) , B1(i) ,                &
     &                                 B1(i1) , Q(i)*Q(il),DERFC(asqrt),A12 , B12 , C12 , e11 , &
     &                                 e22 , e33
                                       short_dd = rij                                                 ! 2-26-09 DU
                                       print 99333, ATOm(i) , ATOm(i1) ,rij, ncy                      ! 2-26-09
                                       write(61,99333) ATOm(i) , ATOm(i1) ,rij, ncy                   ! 2-26-09
99333                                  FORMAT('** The structure-> with short intermolecular distance: ',A5,'...',A5,' = 'f6.3, &
!    &                                 ' Angs, RSS cycle ',I2,' job terminated**')
     &                                 ' Angs, RSS cycle ',I2)                                        ! 5-4-10
99335                                  continue                                                       ! 5-4-10
                                    END IF                                                            ! 2-26-09 DU
                                 IF ( i==i1 ) THEN
! same atoms??
                                    e11 = 0.5*e11
                                    v11 = 0.5*v11
                                    e22 = 0.5*e22
                                    e33 = 0.5*e33
                                 ENDIF
                                 ce11 = ce11 + e11
                                 ve11 = ve11 - v11
                                 e2 = e22 + e2
                                 e3 = e33 + e3
                                 sum123 = sum123 + e11 + e22 + e33
                              ENDIF
                           ENDIF
                        ENDDO loop_k3
                     ENDDO loop_k2
                  ENDDO loop_k111
               ENDDO loop_j
            ENDDO loop_i1
         ENDDO loop_i
      ENDIF  IF_NENtry_GT_1_1
!
!
      ne1 = N11 + n12
      ce1 = ce11 + CE12
      coulomb_en = ce1 + ce2 + ce3
      ve1 = ve11 + VE12
!
      ve = ve1 + ve2 + ve3 + ve4
      vander_en = ve
      repulse_en = e3
      sum_en = coulomb_en + vander_en + repulse_en
      Ec = coulomb_en
      Ev = vander_en
      Er = repulse_en
!
      Pe = sum_en
!
      IF ( NENtry==1 .AND. ILS==1 ) THEN
         WRITE (13,99004) ICYcle - 1 , NENtry , K1 , n12 , N11 , ne1 ,  &
     &                    ntotal
99004    FORMAT ('From POT_E...icycle,nentry,refls & interactions',     &
     &           ' - intra, inter & total =',7I6)
         WRITE (8,99005) ce11 , CE12 , ce1 , ce2 , ce3 , coulomb_en
99005    FORMAT (5x,'CE11,CE12,CE1, CE2, CE3, total = ',6F14.8)
         WRITE (8,99006) ve11 , VE12 , ve1 , ve2 , ve3 , ve4 , ve
99006    FORMAT (5x,'VE11,VE12,VE1,VE2,VE3,VE4,VE = ',7F14.8)
      ENDIF
!                                                                                  ! newline
!----Additional energy contribution from nmol molecules in the assymmetric unit    ! newline
      ec_nmol = 0.0
      ev_nmol = 0.0
      er_nmol = 0.0                                                                ! newline
      E_nmol = 0.0                                                                 ! newline
!
      IF ( NMOl>1 ) THEN
!
         istart1 = 1
         DO imol1 = 1 , NMOl - 1
            iend1 = istart1 + NATm(imol1) - 1
            DO ik3 = istart1 , iend1
               istart2 = 0
               DO ik = 1 , imol1
                  istart2 = istart2 + NATm(ik)
               ENDDO
               istart2 = istart2 + 1
               DO imol2 = imol1 + 1 , NMOl
                  iend2 = istart2 + NATm(imol2) - 1
                  DO ik4 = istart2 , iend2
                     x11 = a*X(ik3)                                                ! newline
                     y11 = b*Y(ik3)                                                ! newline
                     z11 = c*Z(ik3)                                                ! newline
                     x22 = a*X(ik4)                                                ! newline
                     y22 = b*Y(ik4)                                                ! newline
                     z22 = c*Z(ik4)                                                ! newline
                     dx = x11 - x22                                                ! newline
                     dy = y11 - y22                                                ! newline
                     dz = z11 - z22                                                ! newline
                     rij = SQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+          &
     &                     2*dy*dz*c11+2*dx*dz*c22)                                ! newline
                     ec_nmol = ec_nmol + 332.17*Q(ik3)*Q(ik4)/rij
                     ev_nmol = ev_nmol - A1(ik3)*A1(ik4)/(rij**6)
                     er_nmol = er_nmol + B1(ik3)*B1(ik4)                &
     &                         *EXP(-(C1(ik3)+C1(ik4))*rij)                        ! newline
                  ENDDO
                  istart2 = istart2 + NATm(imol2)
               ENDDO
            ENDDO
            istart1 = istart1 + NATm(imol1)
         ENDDO
      ENDIF
!
      E_nmol = ec_nmol + ev_nmol + er_nmol
      Ec = Ec + ec_nmol
      Ev = Ev + ev_nmol
      Er = Er + er_nmol
      Pe = Pe + E_nmol
!
      END SUBROUTINE POT_E
!
!end------------------------subroutine pot_e-----------------------------end
!
      SUBROUTINE MCPY(A,R,N,M,Ms)                                       !  MCPY   5
!                                                                       !  MCPY  30
!        PURPOSE                                                        !  MCPY  35
!           COPY ENTIRE MATRIX                                          !  MCPY  40
!                                                                       !  MCPY  45
!        USAGE                                                          !  MCPY  50
!           CALL MCPY (A,R,N,M,MS)                                      !  MCPY  55
!                                                                       !  MCPY  60
!        DESCRIPTION OF PARAMETERS                                      !  MCPY  65
!           A - NAME OF INPUT MATRIX                                    !  MCPY  70
!           R - NAME OF OUTPUT MATRIX                                   !  MCPY  75
!           N - NUMBER OF ROWS IN A OR R                                !  MCPY  80
!           M - NUMBER OF COLUMNS IN A OR R                             !  MCPY  85
!           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) !  MCPY  90
!                  0 - GENERAL                                          !  MCPY  95
!                  1 - SYMMETRIC                                        !  MCPY 100
!                  2 - DIAGONAL                                         !  MCPY 100
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  !  MCPY 130
!           LOK (FORMERLY LOC)                                          !  MCPY 135
!        METHOD                                                         !  MCPY 145
!           EACH ELEMENT OF MATRIX A IS MOVED TO THE CORRESPONDING      !  MCPY 150
!           ELEMENT OF MATRIX R                                         !  MCPY 155
!                                                                       !  MCPY 160
      USE F77KINDS                        
      IMPLICIT NONE
!*--MCPY3721
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: M , Ms , N
      REAL*8 , DIMENSION(1) :: A , R
      INTENT (IN) A
      INTENT (OUT) R
!
! Local variables
!
      INTEGER :: i , it
!
!*** End of declarations rewritten by SPAG
!
!
!                                                                       !  MCPY 190
      CALL LOK(N,M,it,N,M,Ms)                                           !  MCPY 195
!                                                                       !  MCPY 200
!        COPY MATRIX                                                    !  MCPY 205
      DO i = 1 , it                                                     !  MCPY 215
         R(i) = A(i)                                                    !  MCPY 220
      ENDDO
!
      END SUBROUTINE MCPY                                               !  MCPY 230
!
!end-------------------subroutine mcpy--------------------------------end
!
      SUBROUTINE VALVEC(D,E,A,Ic,N,Id)
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Ic , Id , N
      REAL*8 , DIMENSION(Id,1) :: A
      REAL*8 , DIMENSION(1) :: D , E
      INTENT (IN) Id , N
      INTENT (OUT) Ic
      INTENT (INOUT) A , D , E
!
! Local variables
!
      REAL*8 :: b , c , f , g , h , p , q , r , s , tol
      REAL*8 :: DABS , DSQRT
      INTEGER :: i , ii , j , k , l , m
!
!     GET EIGENVALUES AND EIGENVECTORS OF TRIDIAGONAL MATRIX FROM HOUSEH
! PROGRAM AUTHOR G. W. WESTLEY
!
      Ic = 0
      tol = 16.D0**(-14)
      DO i = 2 , N
         E(i-1) = E(i)
      ENDDO
      E(N) = 0.0D0
      f = 0.0D0
      b = 0.0D0
      DO l = 1 , N
         j = 0
         h = tol*(DABS(D(l))+DABS(E(l)))
         IF ( b<h ) b = h
         DO m = l , N
            IF ( DABS(E(m))<=b ) EXIT
         ENDDO
         IF ( m/=l ) THEN
            DO WHILE ( j<30 )
               j = j + 1
               p = (D(l+1)-D(l))/(2.0D0*E(l))
               r = DSQRT(p*p+1.0D0)
               q = p + r
               IF ( p<0.0D0 ) q = q - 2.0D0*r
               h = D(l) - E(l)/q
               DO i = l , N
                  D(i) = D(i) - h
               ENDDO
               f = f + h
               p = D(m)
               c = 1.0D0
               s = 0.0D0
               i = m
               DO
                  i = i - 1
                  IF ( i<l ) THEN
                     E(l) = s*p
                     D(l) = c*p
                     IF ( DABS(E(l))<=b ) GOTO 50
                     EXIT
                  ELSE
                     g = c*E(i)
                     h = c*p
                     IF ( DABS(p)>=DABS(E(i)) ) THEN
                        c = E(i)/p
                        r = DSQRT(c*c+1.0D0)
                        E(i+1) = s*r*p
                        s = c/r
                        c = 1.0D0/r
                     ELSE
                        c = p/E(i)
                        r = DSQRT(c*c+1.0D0)
                        E(i+1) = s*E(i)*r
                        s = 1.0D0/r
                        c = c/r
                     ENDIF
                     p = c*D(i) - s*g
                     D(i+1) = h + s*(c*g+s*D(i))
                     DO k = 1 , N
                        h = A(k,i+1)
                        A(k,i+1) = s*A(k,i) + c*h
                        A(k,i) = c*A(k,i) - s*h
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            GOTO 100
         ENDIF
 50      D(l) = D(l) + f
      ENDDO
      DO i = 1 , N
         k = i
         p = D(i)
         IF ( i/=N ) THEN
            ii = i + 1
            DO j = ii , N
               IF ( D(j)>p ) THEN
                  k = j
                  p = D(j)
               ENDIF
            ENDDO
         ENDIF
         IF ( k/=i ) THEN
            D(k) = D(i)
            D(i) = p
            DO j = 1 , N
               p = A(j,i)
               A(j,i) = A(j,k)
               A(j,k) = p
            ENDDO
         ENDIF
      ENDDO
      RETURN
 100  Ic = 1
!
!
      END SUBROUTINE VALVEC
!
!end--------------------subroutine valvec---------------------------end
!
      SUBROUTINE LOK(I,J,Ir,N,M,Ms)                                     !  LOK    5
!
!        SUBROUTINE LOK (FORMERLY LOC)                                  !  LOK   25
!                                                                       !  LOK   30
!        PURPOSE                                                        !  LOK   35
!           COMPUTE A VECTOR SUBSCRIPT FOR AN ELEMENT IN A MATRIX OF    !  LOK   40
!           SPECIFIED STORAGE MODE                                      !  LOK   45
!                                                                       !  LOK   50
!        USAGE                                                          !  LOK   55
!                                                                       !  LOK   65
!        DESCRIPTION OF PARAMETERS                                      !  LOK   70
!           I   - ROW NUMBER OF ELEMENT                                 !  LOK   75
!           J   - COLUMN NUMBER  OF ELEMENT                             !  LOK   80
!           IR  - RESULTANT VECTOR SUBSCRIPT                            !  LOK   85
!           N   - NUMBER OF ROWS IN MATRIX                              !  LOK   90
!           M   - NUMBER OF COLUMNS IN MATRIX                           !  LOK   95
!           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX           !  LOK  100
!                  0 - GENERAL                                          !  LOK  105
!                  1 - SYMMETRIC                                        !  LOK  110
!                  2 - DIAGONAL                                         !  LOK  115
!                                                                       !  LOK  120
!        REMARKS                                                        !  LOK  125
!           NONE                                                        !  LOK  130
!                                                                       !  LOK  135
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  !  LOK  140
!           NONE                                                        !  LOK  145
!                                                                       !  LOK  150
!        METHOD                                                         !  LOK  155
!           MS=0   SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N*M ELEMENTS !  LOK  160
!                  IN STORAGE (GENERAL MATRIX)                          !  LOK  165
!           MS=1   SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N*(N+1)/2 IN !  LOK  170
!                  STORAGE (UPPER TRIANGLE OF SYMMETRIC MATRIX). IF     !  LOK  175
!                  ELEMENT IS IN LOWER TRIANGULAR PORTION, SUBSCRIPT IS !  LOK  180
!                  CORRESPONDING ELEMENT IN UPPER TRIANGLE.             !  LOK  185
!           MS=2   SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N ELEMENTS   !  LOK  190
!                  IN STORAGE (DIAGONAL ELEMENTS OF DIAGONAL MATRIX).   !  LOK  195
!                  IF ELEMENT IS NOT ON DIAGONAL (AND THEREFORE NOT IN  !  LOK  200
!                  STORAGE), IR IS SET TO ZERO.                         !  LOK  205
!                                                                       !  LOK  225
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: I , Ir , J , M , Ms , N
      INTENT (IN) I , J , Ms , N
      INTENT (OUT) Ir
!
! Local variables
!
      INTEGER :: irx , ix , jx
!
      ix = I                                                            !  LOK  230
      jx = J                                                            !  LOK  235
      IF ( Ms<1 ) THEN                                                  !  LOK  240
         irx = N*(jx-1) + ix                                            !  LOK  245
      ELSEIF ( Ms==1 ) THEN
         IF ( ix<jx ) THEN                                              !  LOK  255
            irx = ix + (jx*jx-jx)/2                                     !  LOK  260
         ELSE
            irx = jx + (ix*ix-ix)/2                                     !  LOK  270
         ENDIF
      ELSE
         irx = 0                                                        !  LOK  280
                                                                        !  LOK  285
         IF ( ix==jx ) irx = ix                                         !  LOK  290
      ENDIF
      Ir = irx                                                          !  LOK  295
      END SUBROUTINE LOK                                                !  LOK  305
!
!end-----------------------subroutine lok----------------------------end
!
      SUBROUTINE HOUSEH(A,Id,N,D,E)
!     HOUSEHOLDER REDUCTION OF MATRIX TO TRIDIAGONAL FORM BEFORE VALVEC
! PROGRAM AUTHOR G. W. WESTLEY
!
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Id , N
      REAL*8 , DIMENSION(Id,1) :: A
      REAL*8 , DIMENSION(1) :: D , E
      INTENT (IN) Id , N
      INTENT (INOUT) A , D , E
!
! Local variables
!
      REAL*8 :: DSQRT
      REAL*8 :: f , g , h , hh , tol
      INTEGER :: i , j , jp1 , k , l , n1 , nv
!
      tol = 16.D0**(-51.)
      i = N
      nv = N
      n1 = 0                                    ! 4/29/08, HLA, initialized
 100  l = i - 2
      f = A(i,i-1)
      g = 0.0D0
      IF ( l/=0 ) THEN
         DO k = 1 , l
            g = g + A(i,k)*A(i,k)
         ENDDO
      ENDIF
      h = g + f*f
      IF ( g<=tol ) THEN
         E(i) = f
         h = 0.0D0
      ELSE
         l = l + 1
         g = -DSQRT(h)
         IF ( f<0.0 ) g = -g
         E(i) = g
         h = h - f*g
         A(i,i-1) = f - g
         f = 0.0D0
         DO j = 1 , l
            A(j,i) = A(i,j)/h
            g = 0.0D0
            DO k = 1 , j
               g = g + A(j,k)*A(i,k)
            ENDDO
            IF ( j/=l ) THEN
               jp1 = j + 1
               DO k = jp1 , l
                  g = g + A(k,j)*A(i,k)
               ENDDO
            ENDIF
            E(j) = g/h
            f = f + g*A(j,i)
         ENDDO
         hh = f/(h+h)
         DO j = 1 , l
            f = A(i,j)
            E(j) = E(j) - hh*f
            g = E(j)
            DO k = 1 , j
               A(j,k) = A(j,k) - f*E(k) - g*A(i,k)
            ENDDO
         ENDDO
      ENDIF
      D(i) = h
      i = i - 1
      n1 = n1 + 1
      IF ( i>1 ) GOTO 100
      D(1) = 0.0D0
      E(1) = 0.0D0
      DO i = 1 , nv
         l = i - 1
         IF ( ABS(D(i))>1.0D-20 ) THEN          ! 4/28/08, HLA, line used in place of previous one
            DO j = 1 , l
               g = 0.0D0
               DO k = 1 , l
                  g = g + A(i,k)*A(k,j)
               ENDDO
               DO k = 1 , l
                  A(k,j) = A(k,j) - g*A(k,i)
               ENDDO
            ENDDO
         ENDIF
         D(i) = A(i,i)
         A(i,i) = 1.0D0
         IF ( l/=0 ) THEN
            DO j = 1 , l
               A(i,j) = 0.0D0
               A(j,i) = 0.0D0
            ENDDO
         ENDIF
      ENDDO
!
      END SUBROUTINE HOUSEH
!
!end--------------------subroutine househ-------------------------end
!
!----------------------------------------------------------------------------------------
!      Rosenbrock Step Search Subroutine                                                !
!      N. M. Albu 06.13.2005  Modified by S. M. Prasad 10. 26. 2006 & Z.Du 2008         !
!----------------------------------------------------------------------------------------
!
      SUBROUTINE ROSENBROCK(Fx,Fy,Fz)
      USE PMIN_MODULE , ONLY:NTX_min , NTY_min , NTZ_min , NTX_max ,    &
     &    NTY_max , NTZ_max , NMOl , NATm , IDX , ISYstem , N , ILS ,   &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , STHl , SYM , Q , A1 , B1 , C1 , CK , WT_mol ,   &
     &    DMAx , DDMax , GNI , CE12 , VE12 , K1 , N11 , IHKl , IRIj ,   &
     &    IRSs_call , PI , DIFfpc , PRDifpc , NCYc_rss, f_factor_rss,   &
     &    vab_name, bond_type, ispl, ispl_1                             ! from S.P Feb 5,09    
      USE RSS1_MODULE
      USE RSS2_MODULE
      USE RSS3_MODULE
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE RSSPM3_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 , DIMENSION(MATM) :: Fx , Fy , Fz
      INTENT (INOUT) Fx , Fy , Fz
!
! Local variables
!
      REAL*8 :: anglemax , anglemin , coef , dens_per , ec , egymol ,   &
     &          er , ev , e_mol , f_factor , iw , nmc , ro , rt , sm ,  &
     &          sumx , sumy , sumz , sum_epv , sum_x , sum_y , sum_z ,  &
     &          sx , sy , sz , ts , vbsq , volume , vol_per , w ,       &
     &          wd_lim , w_1 , w_2 , w_cor , w_diff , w_previous ,      &
     &          w_save , w_test , w_vt , xcr , xtra_x , xtra_y ,        &
     &          xtra_z , ycr , zcr
      REAL*8 , DIMENSION(6) :: buff , cell_new , cell_ori , save_cell
      CHARACTER(26) :: cline1
      REAL*8 , DIMENSION(4) :: cm
      REAL*8 , DIMENSION(NVD) :: corr , epv , epv_555 , epv_save , fpc ,&
     &                           ipc , new_pcitm , pcitm , pcmax ,      &
     &                           pcmin , pc_ori , pqc , pqc_1st ,       &
     &                           psave , psave_rc , sum_corr , vbest
      LOGICAL :: done = .FALSE.
      REAL*8 :: DSQRT
      REAL*8 , DIMENSION(0:999) :: econf_cyc , epack_cyc , wcyc
      INTEGER :: fa , i , iend , ik , imol , in_cyc , istart , istop ,  &
     &           ivr , j , jj , kill , nka , nrc , ns , nv , nvr , su,  &
     &           kmol, ii
      REAL*8 , DIMENSION(MATM) :: ix , iy , iz , save_x , save_y ,      &
     &                            save_z , xo , xsave , xx , xx_nw ,    &
     &                            yo , ysave , yy , yy_nw , zo , zsave ,&
     &                            zz , zz_nw
      CHARACTER(80) :: line
      CHARACTER(27) :: line1
      REAL :: REAL
      REAL*8 , DIMENSION(NMOl) :: xc , yc , zc
!
      WRITE (61,99001)       
99001 FORMAT (/'Start Rosenbrock step refinement...'/)
      egymol = 0.0D0                                                        ! 4/29/08, HLA, initialization
      EGYmol0 = 0.0D0                                                       ! 4/29/08, HLA, initialization
      nvr = NROtbond
      nv = 6 + 6*NMOl + nvr
      WRITE (61,99002) nv , NMOl , nvr     !  3-5-08 DU
99002 FORMAT ('Initial numbers... # variables (nv) =',i3,               &
     &        ', # molecules (nmol) =',i2,                              &
     &        ', # bonds to be rotated (nvr) =',i2)
!-------------------------------------------------------------------------------
!   Input  variables: AH, x, y, z, cell, cell1, nentry, wt_mol, nsym           ]
!   Output variables: AH, x, y, z, cell, cell1, nentry, ro (density)           ]
!   Rosenbrock Step Search variables:                                          ]
!        - global  ones : ipc, fpc, ix, iy, iz, fx, fy, fz, W, Wcyc            ]
!        - changing ones: pcitm, xx, yy, zz, xx_nw, yy_nw, zz_nw               ]
!        - report  ones : ipc, fpc, ix, iy, iz, fx, fy, fz, difxyz*, W         ]
!        - buffer  ones : buff, buff_new, diffpc, corr, W_cor, W_vt            ]
!        - calcul. ones : su, fa, cm, pqc, vbest, vbsq                         ]
!   *'difxyz' should be understand like 'difx', 'dify', 'difz'.                ]
!                                                                              ]
!   We have: - 'i'      stand for initial                                      ]
!            - 'f'      stand for final                                        ]
!            - 'W'      stand for potential energy                             ]
!            - 'cor'    stand for corrected/correction                         ]
!            - 'itm'    stand for intermediar                                  ]
!            - 'pc'     stand for parameter cell                               ]
!            - 'su'     stand for success                                      ]
!            - 'fa'     stand for failure                                      ]
!            - 'vbest'  stand for vector best direction                        ]
!            - 'vbsq'   stand for square of vector (distance in 12 dimmensions)]
!            - 'qpc'    stand for quantity parameter change                    ]
!            - 'cm'     stand for class multiplier                             ]
!            - 'rc'     stand for Rosenbrock stage counter                     ]
!            - 'ncy'    stand for Rosenbrock cycle counter                     ]
!            - 'difxyz' stand for difference between initial and final xyz     ]
!-------------------------------------------------------------------------------
!      Step 1: Initializing variables                                          ]
!-------------------------------------------------------------------------------
      wd_lim = 1.0D0       ! 4-13-08
      w_diff = 1.0D-6      ! 6-30-08 DU
      kill = 0             ! 2-6-08 DU
      NENtry = 1           ! 2-6-08 DU
      ncy_rss_end = 0      ! 2-26-09 DU
      sum_x = 0.0
      sum_y = 0.0
      sum_z = 0.0
      xtra_x = 0.0
      xtra_y = 0.0
      xtra_z = 0.0         ! all this for calculating the XTRA coordinates
!----- set initial value
      DO i = 1 , nv
         ipc(i) = 0.0
         fpc(i) = 0.0
         pcitm(i) = 0.0
         sum_corr(i) = 0.0
!        if (IRSS_call == 0) DIFfpc(i) = 0.0        ! 1-20-09 DU       
         DIFfpc(i) = 0.0                            ! 1-20-09 DU
         vbest(i) = 0.0
         IF ( i>=4 .AND. i<=6 ) THEN
            pqc(i) = 0.100000D-04
         ELSE
            pqc(i) = 0.100000D-03
         ENDIF
         IF ( IROt>=1 ) THEN
            DO ik = nv - nvr + 1 , nv
!              pqc(ik) = 0.01    ! initial angle in radians for rot_bond
               pqc(ik) = 0.1     ! initial angle in radians for rot_bond
               IDX(ik) = 1       ! refinement code for rot_bond
            ENDDO
         ENDIF
      ENDDO
      DO i = 1 , N
         xx(i) = 0.0
         yy(i) = 0.0
         zz(i) = 0.0
 
         xx_nw(i) = 0.0
         yy_nw(i) = 0.0
         zz_nw(i) = 0.0
 
!        xx_s1(i) = 0.0          ! deleted, 11/20/08, HLA
!        yy_s1(i) = 0.0          ! deleted, 11/20/08, HLA
!        zz_s1(i) = 0.0          ! deleted, 11/20/08, HLA
!        xx_s2(i) = 0.0          ! deleted, 11/20/08, HLA
!        yy_s2(i) = 0.0          ! deleted, 11/20/08, HLA
!        zz_s2(i) = 0.0          ! deleted, 11/20/08, HLA
 
         ix(i) = 0.0
         iy(i) = 0.0
         iz(i) = 0.0
         Fx(i) = 0.0
         Fy(i) = 0.0
         Fz(i) = 0.0
!        difx(i) = 0.0           ! deleted, 11/20/08, HLA
!        dify(i) = 0.0           ! deleted, 11/20/08, HLA
!        difz(i) = 0.0           ! deleted, 11/20/08, HLA
      ENDDO
      DO i = 1 , 6
         buff(i) = 0.0
!        buff_new(i) = 0.0       ! deleted, 11/20/08, HLA
!        buff_s1(i) = 0.0        ! deleted, 11/20/08, HLA
      ENDDO
      DO i = 1 , 4
         cm(i) = 1.0
      ENDDO
      DO i = 1 , nv
         DEL_param(i) = pqc(i)
      ENDDO
!
      DO i = 0 , 999
         wcyc(i) = 0.0
         econf_cyc(i) = 0.0      ! 8-28-08 DU
         epack_cyc(i) = 0.0      ! 8-28-08 DU
      ENDDO
 
      su = 0
      fa = 0
      RC = 1
      nrc = 5
      in_cyc = 0                 ! 7-16-08 DU
      NCY = 1
!     flag = 0                   ! deleted, 11/20/08, HLA
      vbsq = 0.0
      w = 0.0
      w_cor = 0.0
      w_vt = 0.0
      iw = 0.0
!     fw = 0.0                   ! deleted, 11/20/08, HLA
!     zero = 0.0                 ! deleted, 11/20/08, HLA
      nmc = NSYm                 ! nmc is DP
      ro = 0.0
      coef = 0.0
 
!-------------------------------------------------------------------------------
!      Step 2: Writing in 'nrss.output' and 'nrss.sumary'                      ]
!-------------------------------------------------------------------------------
 
      DO i = 1 , 6
         ipc(i) = AH(i)
         fpc(i) = AH(i)
         pcitm(i) = ipc(i)
         cell_ori(i) = AH(i)
      ENDDO
      DO i = 1 , nv
         pc_ori(i) = pcitm(i)
      ENDDO
      DO i = 1 , N
         xx(i) = X(i)
         yy(i) = Y(i)
         zz(i) = Z(i)
         ix(i) = X(i)
         iy(i) = Y(i)
         iz(i) = Z(i)
         X_Ori(i) = X(i)
         Y_Ori(i) = Y(i)
         Z_Ori(i) = Z(i)
      ENDDO
 
!-------------------------------------------------------------------------------
!      Step 3: Calculation                                                     ]
!-------------------------------------------------------------------------------
 
!----Initialize min and max values for cell parameters
      DO i = 1 , 3                                           ! unit cell lengths, +/- 1 Angs or 10%
         IF ( pcitm(1)<10.0 ) THEN
            pcmin(i) = pcitm(i) - 0.1*pcitm(i)
            pcmax(i) = pcitm(i) + 0.1*pcitm(i)
         ELSE
            pcmin(i) = pcitm(i) - 1.0
            pcmax(i) = pcitm(i) + 1.0
         ENDIF
      ENDDO
      anglemin = -0.866025
      anglemax = +0.866025
      DO i = 4 , 6                                           ! unit cell cosines
         pcmin(i) = pcitm(i) - 0.1                           ! changed from .2 to .1, 7/4/08, HLA
         pcmax(i) = pcitm(i) + 0.1                           !                 "
         IF ( pcmin(i)<=anglemin ) pcmin(i) = anglemin       ! 2-6-08 DU
         IF ( pcmax(i)>=anglemax ) pcmax(i) = anglemax       ! 2-6-08 DU
      ENDDO
!
      DO imol = 1 , NMOl
         DO i = 6*imol + 1 , 6*imol + 3                      ! model rotation
            pcmin(i) = -0.1                                  ! changed from .2 to .1, 7/4/08 HLA
            pcmax(i) = 0.1                                   !                "
         ENDDO
         DO i = 6*imol + 4 , 6*imol + 6
            pcmin(i) = -0.25                                 ! changed from .5 to .25, 7/4/08 HLA
            pcmax(i) = 0.25                                  !                 "
         ENDDO
      ENDDO
      IF ( IROt>=1 ) THEN
         DO j = nv - nvr + 1 , nv
! min & max rotation angles in radian for rotbond
            pcmin(j) = -PI/2.0
            pcmax(j) = PI/2.0
         ENDDO
      ENDIF
!-------------------------------------------------------------------------------
 
      DO i = 1 , 6
         buff(i) = pcitm(i)
      ENDDO
!
      Icall_pot_e = 1                                                               ! 2-3-09 DU
      IF ( IUSer_e==1 ) THEN
         CALL USER_E_1(cell_ori,X_Ori,Y_Ori,Z_Ori,ec,ev,er,e_mol,w)
      ELSE
         CALL POT_E(cell_ori,X_Ori,Y_Ori,Z_Ori,ec,ev,er,e_mol,w)
      ENDIF
      WRITE (61,99003)                                                              ! 8-19-08 DU
99003 FORMAT (/'Initial lattice energies before RSS...')         
      WRITE (61,99004) ec , ev , er , w, k1, n11                                    ! 2-3-09 DU
99004    format ('   coul (WC) ',2x,'   vdW (WV) ',2x,'  repl (WR) ', &
                 2x,'  total (W)',6x,'K1',5x,'N11'/f12.6,3F14.6,2I8)     
      Icall_pot_e = 0                                                               ! 2-3-09 DU
!
      w_vt = w
      iw = w
!     fw = w                                                                        ! deleted, 11/20/08, HLA
      wcyc(1) = w
      NENtry = 1
!      if(irot == 0)nentry = 1
 
!-------------------------------------------------------------------------------
!----Start big loop for RSS calculations
!-------------------------------------------------------------------------------
!
!----Save the initial parameters
      DO i = 1 , nv
         IF ( i<=6 ) THEN
            cell_ori(i) = pcitm(i)
            pc_ori(i) = pcitm(i)
         ELSE
!           pcyc(i) = 0.0                                                          ! deleted, 11/20/08, HLA
            pc_ori(i) = 0.0
         ENDIF
         epv(i) = 0.0
         psave(i) = pcitm(i)
         epv_save(i) = epv(i)
         psave_rc(i) = pcitm(i)
      ENDDO
!
      DO i = 1 , N
         X_Ori(i) = xx(i)
         Y_Ori(i) = yy(i)
         Z_Ori(i) = zz(i)
      ENDDO
!
!----Calculate intramolecular E for the initial conformation
      IF ( IROt/=0 ) THEN
         CALL ORTHO_COD(N,cell_ori,X_Ori,Y_Ori,Z_Ori,xo,yo,zo)
         IF ( NE_type==1 ) THEN                                                   ! 9-30-08 DU
            CALL ENERGY_PM3(N,xo,yo,zo,egymol)                                    ! PM3 energy
         ELSE IF ( NE_type==2 ) THEN
            CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)                                  ! b3lyp/631g** energy 9-30-08 DU
         ELSE
            IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)                ! 4-28-09 DU
         ENDIF
         EGYmol0 = egymol
         DEL_egymol = EGYmol0 - EGYmol0                                           ! 3-5-08 DU
         WRITE (61,"('Initial intra E before RSS (egymol0): ', F20.10)")&
     &          EGYmol0                                                           ! 8-19-08 DU
      ENDIF
!
!        do i=1,nv
!           IF (i.GE.4.AND.i.LE.6) THEN
!              pqc(i) = 0.100000D-04
!           ELSE
!              pqc(i) = 0.100000D-03
!           END IF
!        enddo
!
!     if(irot >= 1) then
!        do ik = nv-nvr+1,nv
!        pqc(ik) = 0.001  ! starting angle twist in radian for rot_bond
!        pqc(ik) = 0.2    ! 8-13-08 DU
!        enddo
!      endif
!
      DO i = 1 , nv
         pqc_1st(i) = pqc(i)
      ENDDO
      DO
!
!     do i = 1,nv
!      epv(i) = 0.0
!     enddo
!-----------------------------------------------------------------!
!----Start of a NEW CYCLE
!-----------------------------------------------------------------!
!
!      if (ncy > 100) W_diff = 1.0D-8                    ! 4-8-08 DU
!      if (ncy > 100) W_diff = 1.0D-2                    ! 6-30-08 DU
         IF ( NCY>NCYc_rss ) w_diff = 1.0D-1             ! 10-31-08 copied from S.P
         WRITE (61,99005) NCY
99005       FORMAT (/'Cycle #',I3,'...')                 ! 12/3/08
         CALL PC_LIMITS(ncy,NMOl,IROt,nvr,nv,pcitm,pcmin,pcmax)
         NCY = NCY + 1
         NENtry = 1
!
         IF ( ISYstem/=0 ) THEN
            epv(2) = epv(1)
            psave(2) = psave(1)
            IF ( ISYstem>1 ) THEN
               epv(3) = epv(1)
               psave(3) = psave(1)
            ENDIF
!
       IF ( ISYstem==2 .and. ispl==1) THEN          ! from S.P Feb.5,09
          do imol = 1,nmol                          ! from S.P Feb.5,09
           if(ispl_1(imol)==1)then                  ! from S.P Feb.5,09
               j = 6*imol + 4                       ! from S.P Feb.5,09
               epv(j+1) = epv(j)                    ! from S.P Feb.5,09
               psave(j+1) = psave(j)                ! from S.P Feb.5,09
               epv(j+2) = epv(j)                    ! from S.P Feb.5,09
               psave(j+2) = psave(j)                ! from S.P Feb.5,09
           endif                                    ! from S.P Feb.5,09
          enddo                                     ! from S.P Feb.5,09
        ENDIF                                       ! from S.P Feb.5,09
!
            IF ( ISYstem==3 ) THEN
               epv(5) = epv(4)
               epv(6) = epv(4)
               psave(5) = psave(4)
               psave(6) = psave(4)
            ENDIF
         ENDIF
!
!     WRITE (fmt4(13:14),'(i2)') nv                                         ! deleted, 11/20/08, HLA
!     WRITE (fmt4(32:33),'(i2)') nv                                         ! deleted, 11/20/08, HLA
!
         WRITE (61,99006)                                                   ! 8-19-08 DU
99006    FORMAT ('Parameters, values, min & max ranges...'/ &
     &           '  #       value     minimum     maximum',10X,'epv')
         DO i = 1 , nv                                                      ! 8-19-08 DU
            WRITE (61,99007) i , pcitm(i) , pcmin(i) , pcmax(i) , epv(i)
99007       FORMAT (i3,3F12.6,D16.6)
         ENDDO
 
!        write(61,"(10X,'epv'/)")                                           ! 8-19-08 DU
!       DO ik = 1,nv
!        write(61,"(1X,D16.6)") epv(ik)                                     ! 8-19-08 DU
!       END DO
         IF ( IROt==0 ) THEN                                                ! 8-19-08 DU
            CALL CORR_XYZ(epv,cm,cell_ori,X_Ori,Y_Ori,Z_Ori,cell_new,   &
     &                    xx_nw,yy_nw,zz_nw)
         ELSE
            I_Nv = 1                                                        ! 8-19-08 DU
            CALL CORR_XYZ_ROT(epv,cm,cell_ori,X_Ori,Y_Ori,Z_Ori,        &
     &                        cell_new,xx_nw,yy_nw,zz_nw)                   ! 8-19-08 DU
         ENDIF
!
         IF ( IMOde>0 ) THEN
            CALL USER_E_1(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol,w)
         ELSE
            CALL POT_E(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol,w)
         ENDIF
         IF ( IROt/=0 ) w = w + DEL_egymol
!
!        IF ( IRSs_call >=5 .AND. (w>0.0 .OR. ABS(w_cor-w)>5.0) ) RETURN                      
!        IF ( IRSs_call >=5 ) RETURN                                                            ! 2-3-09
         IF ( IRSs_call>0 ) then                                                                ! 1-20-09
          WRITE(13,"('recall RSS IRSS_call: ', I3, 2F16.6)") IRSs_call, w, w_cor                ! 1-20-09
         END IF                                                                                 ! 1-20-09
!
         WRITE (61,99008) NCY - 1
99008       FORMAT (/'Initial energies cycle #',I3,'...')
!       if((ncy-1) > 1)  WRITE(61,1133) Ec,Ev,Er,Ec+Ev+Er,del_egymol,W, abs(w_cor-w), nentry    ! 6/30/08
         WRITE (61,99009) ec , ev , er , ec + ev + er , DEL_egymol , w ,&
     &                    NENtry , I_Nv                                                         ! 8-19-08
99009       format ('   coul (WC)      vdW (WV)     repl (WR)     lattE (W)', &
                    '        conf E    total (WT)     nentry'/ &      
                    f12.6,5F14.6,2I6)           
         IF ( (NCY-1)==1 ) THEN
            econf_cyc(NCY-1) = 0.0D-12                                                          ! 8-28-08 DU
         ELSE
            econf_cyc(NCY-1) = DEL_egymol                                                       ! 8-28-08 DU
         ENDIF
         epack_cyc(NCY-1) = ec + ev + er                                                        ! 8-28-08 DU
!
         WRITE (28,99010) NCY - 1 , ATOm(1) , X_Ori(1) , Y_Ori(1) ,     &
     &                    Z_Ori(1) , xx_nw(1) , yy_nw(1) , zz_nw(1) , w
99010    FORMAT (/,' ncy,atom1,xyz_ori,xyz_nw,W',I2,2x,a6,6F12.8,f14.8, &
     &           ' <-------------',/)
!
         sx = 0.0
         sy = 0.0
         sz = 0.0
         DO ik = 1 , 6                                  ! 2-6-08
            save_cell(ik) = pcitm(ik)                   ! 2-6-08
         ENDDO                                          ! 2-6-08
         DO ik = 1 , N
            save_x(ik) = xx_nw(ik)                      ! 2-6-08
            save_y(ik) = yy_nw(ik)                      ! 2-6-08
            save_z(ik) = zz_nw(ik)                      ! 2-6-08
            sx = sx + xx_nw(ik)
            sy = sy + yy_nw(ik)
            sz = sz + zz_nw(ik)
         ENDDO
         xcr = sx/N
         ycr = sy/N
         zcr = sz/N
         ATOm(N+1) = 'XTRA  '
!       in_cyc = 0
!
!----------------------- Start of loop 555 for stage search ---------------------------------
         DO RC = 1 , nrc
            NENtry = 2                                  ! 6-30-08
            I_Nv = 0
            WRITE (61,99011) NCY - 1 , RC  ! 3-5-08
99011       FORMAT (/'Start of Rosenbrock cycle & stage =',2I3/)
            vbsq = 0.0
!8045       format (' ncy, rc, in_cyc, nentry, n11, calcd energy:',5i5,F12.6)
!
!        in_cyc = 0                                     ! 7-16-08 DU
!767     in_cyc = in_cyc + 1                            ! 7-16-08 DU
!
            IF ( in_cyc==1 ) THEN                       ! 7-16-08 DU
               DO j = 1 , nv
                  IF ( j>=4 .AND. j<=6 ) THEN
                     pqc(j) = 0.100000D-06
                  ELSE
                     pqc(j) = 0.100000D-05
                  ENDIF
               ENDDO
            ENDIF                                       ! 7-16-08 DU
            WRITE (61,99012)
99012       FORMAT (3X,'I    delta param   +',5X,'epv',6X,              &     ! altered, 7/4/08 HLA
     &              '= new param chng')                 
            DO i = 1 , nv
!          if (irot >= 1 .and. i > nv-nvr) then         ! 8-19-08 DU
!             del_param(i) = pqc_1st(i)
!          else
!             del_param(i) = pqc(i)
!          end if
               DEL_param(i) = pqc(i)
               epv_555(i) = epv(i)                      ! 8-19-08 DU
               WRITE (61,99013) i , DEL_param(i) , epv(i) , DEL_param(i)&
     &                          + epv(i)
99013          FORMAT (I4,3D16.6)
            ENDDO
            WRITE (61,99014)
99014       FORMAT (                                                    &
     &' ** del_param is from best del_param of previously stage or cycle&
     &'/' ** epv is from best epv of previously vector search'/         &
     &' ** new param (new epv) is used to next satge or cycle')
            IF ( ISYstem/=0 ) THEN
               epv(2) = epv(1)
               pcitm(2) = pcitm(1)
               IF ( ISYstem>1 ) THEN
                  epv(3) = epv(1)
                  pcitm(3) = pcitm(1)
               ENDIF
!
       IF ( ISYstem==2 .and. ispl==1) THEN              ! from S.P Feb 5,09
          do imol = 1,nmol                              ! from S.P Feb 5,09
           if(ispl_1(imol)==1)then                      ! from S.P Feb 5,09
               j = 6*imol + 4                           ! from S.P Feb 5,09
               epv(j+1) = epv(j)                        ! from S.P Feb 5,09
               pcitm(j+1) = pcitm(j)                    ! from S.P Feb 5,09
               epv(j+2) = epv(j)                        ! from S.P Feb 5,09
               pcitm(j+2) = pcitm(j)                    ! from S.P Feb 5,09
           endif                                        ! from S.P Feb 5,09
          enddo                                         ! from S.P Feb 5,09
        ENDIF                                           ! from S.P Feb 5,09
!
               IF ( ISYstem==3 ) THEN
                  epv(5) = epv(4)
                  epv(6) = epv(4)
                  pcitm(5) = pcitm(4)
                  pcitm(6) = pcitm(4)
               ENDIF
            ENDIF
!
            IF ( RC<=1 ) THEN
!          endif
               IF ( IMOde>0 ) THEN
                  CALL USER_E_1(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,    &
     &                          e_mol,w_1)
               ELSE
                  CALL POT_E(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol, &
     &                       w_1)
               ENDIF
               WRITE (61,99015) RC
99015          FORMAT (/'Starting energies for new stage #',I3,'...')
               IF ( IROt/=0 ) w_1 = w_1 + DEL_egymol
               WRITE (61,99009) ec , ev , er , ec + ev + er ,           &
     &                          DEL_egymol , w_1 , NENtry
               w_test = w_1
!        WRITE (fmt5(12:13),'(i2)') nv     ! deleted, 11/20/08, HLA
               w_save = w_test
               w = w_test
!        WRITE(61,8884) W_save, W
99017          FORMAT (/'Total energies, W_save & W = ',2F11.6)
!        Wum_deriv1, param_fract               ! 6/4/08, HLA, added
99018          FORMAT (' ...sum 1st deriv =',f11.7,                     &
     &                 /'    l.s. parameter change multiplier reset to',&
     &                 f4.1)
!        wrITE(61,8891) cell_new
99019          FORMAT (/'Unit cell:',6F11.6)                               ! 3-5-08 DU
            ENDIF
!--------- write out new cell and coords -------------------------------------
            sx = 0.0
            sy = 0.0
            sz = 0.0
            DO ik = 1 , N
               sx = sx + xx_nw(ik)
               sy = sy + yy_nw(ik)
               sz = sz + zz_nw(ik)
            ENDDO
            xcr = sx/N
            ycr = sy/N
            zcr = sz/N
            ATOm(N+1) = 'XTRA  '
!------------------------------------------------------------------------------
            DO ik = 1 , N
               xx(ik) = xx_nw(ik)
               yy(ik) = yy_nw(ik)
               zz(ik) = zz_nw(ik)
            ENDDO
!
!--------- write out new cell and coords -------------------------------------
            WRITE (61,99020)
99020       FORMAT (/'Unit cell parameters and coordinates...'/)
            sx = 0.0
            sy = 0.0
            sz = 0.0
            WRITE (61,99056) (cell_new(ik),ik=1,6)                         ! 2-6-08 DU
            DO ik = 1 , N
               WRITE (61,99057) ATOm(ik) , ik , xx(ik) , yy(ik) , zz(ik)   ! 2-6-08 DU
               sx = sx + xx(ik)
               sy = sy + yy(ik)
               sz = sz + zz(ik)
            ENDDO
            xcr = sx/N
            ycr = sy/N
            zcr = sz/N
            ATOm(N+1) = 'XTRA  '
            WRITE (61,99057) ATOm(N+1) , N - N , xcr , ycr , zcr           ! 2-6-08 DU
            WRITE (61,99021) w_1 - DEL_egymol , DEL_egymol , w             ! 8-19-08 DU
99021       FORMAT (/'Search for minimum energy...starting lattice, intramolecular and total energies =', &
                    3F12.6/)
            WRITE (61,99022) (pc_ori(i),i=1,6)                             ! 7-16-08 DU
99022       FORMAT ('Original cell parameters:',3f9.4,3F10.5/)      
            WRITE (61,99023)                                               ! 2-6-08
99023       FORMAT (                                                    &
     &'   #  delta param + started epv -/+ change(epv)   new        tota&
     &l E',5X,'E_conf    ',1X,'Ec',9X,'Ev',9X,'Er',8X,'k1',3X,'n11',    &
     &' N SU FA NV')
 
!--------------------- loop 666 for each parameter ---------------------------------------
            DO i = 1 , nv
               DO ik = 1 , N
                  xsave(ik) = xx(ik)
                  ysave(ik) = yy(ik)
                  zsave(ik) = zz(ik)
!              IF ( IROt.GE.1 .AND. i.GT.nv-nvr ) THEN    ! deleted, 11/20/08, HLA
!                 x_rot_sv(ik) = xx(ik)                   ! deleted, 11/20/08, HLA
!                 y_rot_sv(ik) = yy(ik)                   ! deleted, 11/20/08, HLA
!                 z_rot_sv(ik) = zz(ik)                   ! deleted, 11/20/08, HLA
!              ENDIF                                      ! deleted, 11/20/08, HLA
               ENDDO
               su = 0
               fa = 0
!          del_param(i) = pqc(i)                          ! 8-19-08 check?
               IF ( IDX(i)/=0 ) THEN
!----------------- loop 777 for minmun E search --------------------------------
                  DO jj = 1 , 10   ! 3-5-08 DU
                     epv(i) = epv(i) + DEL_param(i)
                     pcitm(i) = pc_ori(i) + epv(i)
                     IF ( (pcitm(i)<=pcmin(i)) .OR. (pcitm(i)>=pcmax(i))&
     &                    ) THEN
                        IF ( pcitm(i)<pcmin(i) ) WRITE (61,99024) i ,   &
     &                       DEL_param(i) , pcitm(i) , pcmin(i)
99024                   FORMAT ('**',I4,'...param change (',D12.5,      &
     &                          ') -> new param (',f9.5,')',            &
     &                          '...param < allowed param minimum (',   &
     &                          F9.5,')')
                        IF ( pcitm(i)>pcmax(i) ) WRITE (61,99025) i ,   &
     &                       DEL_param(i) , pcitm(i) , pcmax(i)
99025                   FORMAT ('**',I4,'...param change (',D12.5,      &
     &                          ') -> new param (',f9.5,')',            &
     &                          '...param > allowed param maximum (',   &
     &                          F9.5,')')
                        pcitm(i) = pcitm(i) - DEL_param(i)
                        epv(i) = epv(i) - DEL_param(i)
                        IF ( i>=4 .AND. i<=6 ) THEN
                           DEL_param(i) = DEL_param(i)/3.0
                           EXIT
                        ENDIF
                        IF ( i<=3 .OR. i>6 ) DEL_param(i) = pqc_1st(i)
                        pcitm(i) = pc_ori(i) + epv(i)
                        EXIT
                     ENDIF
!--------------------------------------------------------------------------------------------------------------------------
!          if(irot >= 1 .and. i > nv-nvr) then            ! 3-5-08 DU
!             ivr = i-(nv-nvr)
!             dtheta = epv(i)
!             if (ivr == 1 .AND. nvr == 1) then
!                do j=1,n
!                   x_rot(j) = xx(j)                      ! 3-25-08 DU
!                   y_rot(j) = yy(j)                      ! 3-25-08 DU
!                   z_rot(j) = zz(j)                      ! 3-25-08 DU
!                end do
!             else
!                do j=1,n
!                    x_rot(j) = x_rot_sv(j)               ! 3-28-08 DU
!                    y_rot(j) = y_rot_sv(j)               ! 3-28-08 DU
!                    z_rot(j) = z_rot_sv(j)               ! 3-28-08 DU
!                end do
!             end if
!             CALL BEND_BOND_NEW(ivr,n,cell_new,x_rot,y_rot,z_rot,xo,yo,zo,egymol)
!             CALL FRACTIONAL_COD(n,cell_new,xo,yo,zo,xx_nw,yy_nw,zz_nw)
!             del_egymol = egymol - egymol0
!             go to 7721
!          endif
!
                     IF ( ISYstem/=0 ) THEN
                        epv(2) = epv(1)
                        IF ( ISYstem>1 ) epv(3) = epv(1)
!
                        IF ( ISYstem==2 .and. ispl==1) THEN       ! from S.P Feb 5,09
                         do imol = 1,nmol                         ! from S.P Feb 5,09
                          if(ispl_1(imol)==1)then                 ! from S.P Feb 5,09
                           j = 6*imol + 4                         ! from S.P Feb 5,09
                           epv(j+1) = epv(j)                      ! from S.P Feb 5,09
                           epv(j+2) = epv(j)                      ! from S.P Feb 5,09
                          endif                                   ! from S.P Feb 5,09
                         enddo                                    ! from S.P Feb 5,09
                        ENDIF                                     ! from S.P Feb 5,09
!
                        IF ( ISYstem==3 ) THEN
                           epv(5) = epv(4)
                           epv(6) = epv(4)
                        ENDIF
                     ENDIF
                     IF ( IROt==0 ) THEN                                    ! 8-19-08 DU
                        CALL CORR_XYZ(epv,cm,cell_ori,X_Ori,Y_Ori,Z_Ori,&
     &                                cell_new,xx_nw,yy_nw,zz_nw)
                     ELSE
                        I_Nv = i                                            ! 8-19-08 DU
                        CALL CORR_XYZ_ROT(epv,cm,cell_ori,X_Ori,Y_Ori,  &
     &                     Z_Ori,cell_new,xx_nw,yy_nw,zz_nw)                ! 8-19-08 DU
                     ENDIF
!
                     NENtry = NENtry + 1                                    ! 2july08
                     w_cor = 0.0
                     IF ( IMOde>0 ) THEN
                        CALL USER_E_1(cell_new,xx_nw,yy_nw,zz_nw,ec,ev, &
     &                                er,e_mol,w_cor)
                     ELSE
                        CALL POT_E(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er, &
     &                             e_mol,w_cor)
                     ENDIF
                                           ! 3-5-08 DU
                     IF ( IROt>=1 ) w_cor = DEL_egymol + w_cor
                                           ! egymol = 0.0 if irot = 0
                     IF ( jj==1 ) THEN                                      ! 8-19-08 DU
                        IF ( IROt>=1 .AND. i>nv-nvr ) THEN
!           WRITE(61,"(I6,F13.6,D16.6,6F13.6,2I6,I5,1X,2I3)") &
                           WRITE (61,                                   &
     &                            "(I4,3D14.6,6F11.6,2I6,2I2,1X,I2,I3)")&
     &                            i , DEL_param(i) , epv_555(i) , epv(i)&
     &                            , pcitm(i)*57.29578 , w_cor ,         &
     &                            DEL_egymol , ec , ev , er , K1 , N11 ,&
     &                            RC , su , fa , I_Nv                       ! 8-19-08
                        ELSE
!           WRITE(61,"(I6,D13.6,D16.6,6F13.6,2I6,I5,1X,2I3)") &
                           WRITE (61,                                   &
     &                            "(I4,3D14.6,6F11.6,2I6,2I2,1X,I2,I3)")&
     &                            i , DEL_param(i) , epv_555(i) , epv(i)&
     &                            , pcitm(i) , w_cor , DEL_egymol , ec ,&
     &                            ev , er , K1 , N11 , RC , su , fa ,   &
     &                            I_Nv                                      ! 8-19-08
                        ENDIF
                     ELSEIF ( IROt>=1 .AND. i>nv-nvr ) THEN
!           WRITE(61,"(I6,F13.6,D16.6,6F13.6,2I6,I5,1X,2I3)") &
!      &       i,pc_ori(i),epv(i),pcitm(i)*57.29578,W_cor,del_egymol,Ec,Ev,Er,k1,n11,rc,SU,FA
                        WRITE (61,                                      &
     &                    "(I4,D14.6,14X,D14.6,6F11.6,2I6,2I2,1X,I2,I3)"&
     &                    ) i , DEL_param(i) , epv(i) , pcitm(i)        &
     &                      *57.29578 , w_cor , DEL_egymol , ec , ev ,  &
     &                      er , K1 , N11 , RC , su , fa , I_Nv             ! 8-19-08
                     ELSE
!           WRITE(61,"(I6,D13.6,D16.6,6F13.6,2I6,I5,1X,2I3)") &
!      &       i,pc_ori(i),epv(i),pcitm(i),W_cor,del_egymol,Ec,Ev,Er,k1,n11,rc,SU,FA
                        WRITE (61,                                      &
     &                    "(I4,D14.6,14X,D14.6,6F11.6,2I6,2I2,1X,I2,I3)"&
     &                    ) i , DEL_param(i) , epv(i) , pcitm(i) ,      &
     &                      w_cor , DEL_egymol , ec , ev , er , K1 ,    &
     &                      N11 , RC , su , fa , I_Nv                       ! 8-19-08
                     ENDIF
!
                     IF ( (w_cor-w_save)>0.0 .AND. su==1 )              &
     &                    WRITE (61,99026) i , DEL_param(i) , pcitm(i) ,&
     &                    w_cor , w_save
99026                FORMAT ('**',I2,'...param change (',D12.5,         &
     &                       ') -> new param (',f9.5,')',               &
     &                       '...current E (',f11.6,') > previous E (', &
     &                       f11.6,')')                                    ! 7/4/08 HLA
! -----------------------------------------------------------------------------------------------
!            IF (W_cor - W_save)7135,7135,7136
!7135        if ((W_cor - W_save) .ge.  1.0 ) then     ! decreases by more than 1 (-9 to -11),
                     IF ( w_cor>w_save ) THEN
                        fa = 1
                        epv(i) = epv(i) - DEL_param(i)
                        pcitm(i) = pc_ori(i) + epv(i)
                        DO ik = 1 , N
                           xx(ik) = xsave(ik)
                           yy(ik) = ysave(ik)
                           zz(ik) = zsave(ik)
                        ENDDO
                        w = w_save
                        IF ( ABS(w_cor-w)<=w_diff ) WRITE (61,99027) i ,&
     &                       w_cor , w , ABS(w_cor-w) , w_diff ,        &
     &                       DEL_egymol                                        ! altered 7/10/08 HLA
99027                   FORMAT ('**',I4,"  current & previous E's =",   &
     &                          2D16.8,', E diff =',d15.8,' <',d11.4,   &
     &                          ', conf E =',f12.5)                            ! altered 7/10/08 HLA
                        IF ( ABS(w_cor-w)<=w_diff ) EXIT                       ! 2-6-08 DU
                        IF ( su==1 ) THEN
                           DEL_param(i) = DEL_param(i)/3.0
                           EXIT
                        ENDIF
                        IF ( jj==10 ) THEN                                     ! 3-5-08 DU
                           DEL_param(i) = pqc_1st(i)
                           EXIT
                        ENDIF
!               if(irot >= 1 .and. i >  nv-nvr)then                            ! 8-19-08 DU
!                if (mod(jj,2) == 0) then
!                 del_param(i) = -0.5*del_param(i)
!                 GO TO 8135
!                else
!                 del_param(i) = -del_param(i)
!                 GO TO 8135
!                end if
!               else                                                           ! 3-5-08 DU
!                del_param(i) = -0.5*del_param(i)
!                GOTO 8135
!               end if
                        DEL_param(i) = -0.5*DEL_param(i)                       ! 8-19-08 DU
                     ELSE
!                       IF ( (w_cor<0.0) .AND. ((w_cor-w_save)>=1.0) )  &
                        IF ( (w_cor<0.0) .AND. (ABS(w_cor-w_save)>=1.0) .AND. NCY > 3)  &   ! 5-4-10 DU-
     &                       THEN                                              ! if W_cor < 0 & E decreases
                                                                               ! by more than 1 ( eg, -9 to -11),
                           WRITE (61,99028) i , w_cor , w_save                 ! use previous param value
99028                      FORMAT ('** ',I3,'...current E (',f11.6,     &
     &                             ') < previous E (',f11.6,')',        &
     &                      '...big jump in E, use previous param value'&
     &                      )
                           epv(i) = epv(i) - DEL_param(i)                      ! 7-16-08 DU
                           pcitm(i) = pc_ori(i) + epv(i)                       ! 7-16-08 DU
!               DO J = 1,nv
!                  IF (J.GE.4.AND.J.LE.6) THEN
!                     del_param(J) = 0.100000D-04
!                  ELSE
!                     del_param(J) = 0.100000D-03
!                  END IF
!               END DO
                           DO ik = 1 , N                                       ! 7-16-08 DU
                              xx(ik) = xsave(ik)                               ! 7-16-08 DU
                              yy(ik) = ysave(ik)                               ! 7-16-08 DU
                              zz(ik) = zsave(ik)                               ! 7-16-08 DU
                           ENDDO                                               ! 7-16-08 DU
                           w = w_save                                          ! 7-16-08 DU
                           in_cyc = 1                                          ! 7-16-08 DU
                           EXIT                                                ! 7-16-08 DU
                        ENDIF                                                  ! 7-16-08 DU
!7135              su = 1
!                  su = 1
                        su = su + 1                                            ! 7-16-08 DU
                        DO j = 1 , N
                           xx(j) = xx_nw(j)
                           yy(j) = yy_nw(j)
                           zz(j) = zz_nw(j)
                        ENDDO
                        w = w_cor
                        psave(i) = pcitm(i)
                        DO ik = 1 , N
                           xsave(ik) = xx(ik)
                           ysave(ik) = yy(ik)
                           zsave(ik) = zz(ik)
                        ENDDO
                        w_save = w
                        IF ( fa/=0 ) EXIT
                        DEL_param(i) = 3*DEL_param(i)
                        IF ( in_cyc==1 ) THEN
                           IF ( su>=5 ) EXIT
                        ENDIF                                                  ! 7-16-08 DU
                     ENDIF
                  ENDDO
               ENDIF
 
!        if (i ==  nv-nvr) then
!           write(61,"('cell_new ',6F16.6)") (pcitm(ik),ik=1,6)
!           write(61,"('cell_new ',6F16.6)") (psave(ik),ik=1,6)
!           write(61,"('cell_new ',6F16.6)") (cell_new(ik),ik=1,6)
!        end if
            ENDDO
            DO j = 1 , N
               Fx(j) = xx(j)
               Fy(j) = yy(j)
               Fz(j) = zz(j)
            ENDDO
            WRITE (61,                                                  &
     &"('   #         VBEST(I)   =    P(I)      -     PBASE(I)         e&
     &pv(I)')")                                                                ! 8-19-08 DU
            DO i = 1 , nv
               vbest(i) = pcitm(i) - psave_rc(i)
               WRITE (61,"(I4,1X,D16.6,2F16.6,D16.6)") i , vbest(i) ,   &
     &                pcitm(i) , psave_rc(i) , epv(i)                          ! 6-23-08 DU
               vbsq = vbsq + (vbest(i))**2
               fpc(i) = pcitm(i)
               psave_rc(i) = pcitm(i)
            ENDDO
            IF ( IROt==0 ) THEN                                                ! 8-19-08
               CALL CORR_XYZ(epv,cm,cell_ori,X_Ori,Y_Ori,Z_Ori,cell_new,&
     &                       xx_nw,yy_nw,zz_nw)
            ELSE
               I_Nv = 1                                                        ! 8-19-08
               CALL CORR_XYZ_ROT(epv,cm,cell_ori,X_Ori,Y_Ori,Z_Ori,     &
     &                           cell_new,xx_nw,yy_nw,zz_nw)                   ! 8-19-08
            ENDIF
            w_cor = 0.0
            IF ( IMOde>0 ) THEN
               CALL USER_E_1(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol, &
     &                       w_cor)
            ELSE
               CALL POT_E(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol,    &
     &                    w_cor)
            ENDIF
!--------- write out new cell and coords 8-19-08 DU -----------------------------
            WRITE (61,99029) RC
99029       FORMAT (/'After stage...cell parameters and coordinates...',&
     &              I4/)
            sx = 0.0
            sy = 0.0
            sz = 0.0
            WRITE (61,99056) cell_new
                                   ! cell parameters
            DO ik = 1 , N
               WRITE (61,99057) ATOm(ik) , ik , xx_nw(ik) , yy_nw(ik) , &
     &                          zz_nw(ik)                               
               sx = sx + xx_nw(ik)
               sy = sy + yy_nw(ik)
               sz = sz + zz_nw(ik)
            ENDDO
            xcr = sx/N
            ycr = sy/N
            zcr = sz/N
            ATOm(N+1) = 'XTRA  '
            WRITE (61,99057) ATOm(N+1) , N - N , xcr , ycr , zcr
!-------------------------------------------------------------------------------
            DO i = 1 , 6
               buff(i) = pcitm(i)
            ENDDO
!-------------- the following lines removed on 8-19-08 ----------------------------
!       del_egymol = 0.0
!       if(irot > 0)then
!          CALL ORTHO_COD(n,buff,xx,yy,zz,xo,yo,zo)
!          CALL ENERGY_pm3(n,xo,yo,zo,egymol)
!          del_egymol = egymol - egymol0
!       endif
!       if(imode > 0) then
!          CALL USER_E_1(buff,xx,yy,zz,Ec,Ev,Er,E_mol,W_cor)
!       else
!          CALL POT_E(buff,xx,yy,zz,Ec,Ev,Er,E_mol,W_cor)
!       endif
            w = DEL_egymol + w_cor                    ! egymol = 0.0 if irot = 0
!--------- write out new cell and coords 3-5-08 DU -----------------------------
!      write (61,8381)
99030       FORMAT (/                                                   &
     &       'Before vector search...cell parameters and coordinates...'&
     &       /)
            sx = 0.0
            sy = 0.0
            sz = 0.0
!        WRITE(61,5646) buff    ! cell parameters
            DO ik = 1 , N
!          WRITE(61,81321)atom(ik),ik,xx(ik),yy(ik),zz(ik)     
               sx = sx + xx(ik)
               sy = sy + yy(ik)
               sz = sz + zz(ik)
            ENDDO
            xcr = sx/N
            ycr = sy/N
            zcr = sz/N
            ATOm(N+1) = 'XTRA  '
!       write(61,81321)atom(n+1),n-n,xcr,ycr,zcr
!-------------------------------------------------------------------------------
! - here is starting the Vector Step Search (inside RSS calculation)           ]
!-------------------------------------------------------------------------------
!
!        IF (vbsq .le. 1.0D-15) go to 8065                                ! 4-23-08 for testing DU
            WRITE (61,99031) w_cor , DEL_egymol , w
99031       FORMAT (/                                                   &
     &'Start vector search, inter-, intra, Total energy before changes =&
     &',3F11.6,' -->'/)
!        WRITE (fmt6(30:31),'(i2)') nv     ! deleted, 11/20/08, HLA
!        WRITE (fmt6(57:58),'(i2)') nv     ! deleted, 11/20/08, HLA
!        WRITE (fmt7(14:15),'(i2)') nv     ! deleted, 11/20/08, HLA
!        WRITE (fmt7(34:35),'(i2)') nv     ! deleted, 11/20/08, HLA
!        WRITE (fmt7(54:55),'(i2)') nv     ! deleted, 11/20/08, HLA
!        WRITE (fmt7(74:75),'(i2)') nv     ! deleted, 11/20/08, HLA
            DO i = 1 , nv
               pqc(i) = DEL_param(i)
            ENDDO
            NENtry = 1
            sm = 0.25
            su = 0
            fa = 0
!        kid = 0                           ! deleted, 11/20/08, HLA
            w_vt = 0.0
!        ism = 0                           ! deleted, 11/20/08, HLA
            DO jj = 1 , 10          
               DO i = 1 , nv
                  epv_save(i) = epv(i)
                  psave(i) = pcitm(i)
               ENDDO
               DO i = 1 , N
                  xsave(i) = xx(i)
                  ysave(i) = yy(i)
                  zsave(i) = zz(i)
               ENDDO
               w_save = w
               DO i = 1 , nv
                  corr(i) = epv(i)
!                 IF ( i>=4 .AND. i<=6 ) THEN
!                    epv(i) = 0.1*sm*vbest(i) + epv(i)
!                 ELSE
!                    epv(i) = sm*vbest(i) + epv(i)
!                 ENDIF
                  epv(i) = sm*vbest(i) + epv(i)                            !  2-3-09
                  IF ( IDX(i)==0 ) epv(i) = 0.0
               ENDDO
!
               DO i = 1 , nv
                  IF ( (pcitm(i)<=pcmin(i)) .OR. (pcitm(i)>=pcmax(i)) ) &
     &                 THEN
!                    epv(i) = corr(i)
!                    pcitm(i) = pc_ori(i) + epv(i)
                     epv(i) = epv_save(i)                                  ! 2-3-09
                     pcitm(i) = psave(i)                                   ! 2-3-09
                     vbest(i) = 0.0                                        ! 2-3-09
                     if (pcitm(i) <= pcmin(i)) then                        ! 2-3-09
                     WRITE(61,"('i, pcitm(i) <= pcmin(i) :',I4,3f16.6)") & ! 2-3-09
     &                           i,pcitm(i),pcmin(i),epv(i)                ! 2-3-09
!                    pcitm(i) = pc_ori(i)                                  ! 5-4-10
!                    psave(i) = pc_ori(i)                                  ! 5-4-10
!                    epv(i)   = pc_ori(i)                                  ! 5-4-10
                     else                                                  ! 2-3-09
                     WRITE(61,"('i, pcitm(i) >= pcmax(i) :',I4,3f16.6)") & ! 2-3-09
     &                           i,pcitm(i),pcmax(i),epv(i)                ! 2-3-09
!                    pcitm(i) = pc_ori(i)                                  ! 5-4-10
!                    psave(i) = pc_ori(i)                                  ! 5-4-10
!                    epv(i)   = pc_ori(i)                                  ! 5-4-10
                     end if
                  ENDIF
               ENDDO
               DO i = 1 , 6
                  buff(i) = pcitm(i)
               ENDDO
               IF ( ISYstem/=0 ) THEN
                  buff(2) = buff(1)
                  epv(2) = epv(1)
                  IF ( ISYstem>1 ) THEN
                     buff(3) = buff(1)
                     epv(3) = epv(1)
                  ENDIF
!
       IF ( ISYstem==2 .and. ispl==1) THEN                                 ! from S.P Feb 5,09
          do imol = 1,nmol                                                 ! from S.P Feb 5,09
           if(ispl_1(imol)==1)then                                         ! from S.P Feb 5,09
               j = 6*imol + 4                                              ! from S.P Feb 5,09
               epv(j+1) = epv(j)                                           ! from S.P Feb 5,09
               buff(j+1) = buff(j)                                         ! from S.P Feb 5,09
               epv(j+2) = epv(j)                                           ! from S.P Feb 5,09
               buff(j+2) = buff(j)                                         ! from S.P Feb 5,09
           endif                                                           ! from S.P Feb 5,09
          enddo                                                            ! from S.P Feb 5,09
        ENDIF                                                              ! from S.P Feb 5,09
!
                  IF ( ISYstem==3 ) THEN
                     buff(5) = buff(4)
                     buff(6) = buff(4)
                     epv(5) = epv(4)
                     epv(6) = epv(4)
                  ENDIF
               ENDIF
!           WRITE (fmt8(14:15),'(i2)') nv                                  ! deleted, 11/20/08, HLA
!           WRITE (fmt8(34:35),'(i2)') nv                                  ! deleted, 11/20/08, HLA
!           WRITE (fmt8(54:55),'(i2)') nv                                  ! deleted, 11/20/08, HLA
               WRITE (61,99032)                                            ! 8-19-08 DU
!           WRITE(61,2140)
99032          FORMAT (6X,'i',5X,                                       &
     &'new value(i) =  vbest(I)  *   sm   +  current(i)= new epv(i)  +  &
     &original(i)')
               DO i = 1 , nv
!              new_pcitm = psave(i)+vbest(i)*sm                            ! 3-25-08 DU
               new_pcitm(i) = psave(i) + vbest(i)*sm                       ! 8-19-08 DU
                  WRITE (61,99033) i , new_pcitm(i) , vbest(i) , sm ,   &
     &                             psave(i) , epv(i) , pc_ori(i)
99033             FORMAT (I7,f14.6,f15.6,f11.6,f12.6,D14.6,f12.6)
               ENDDO
 
               CALL CORR_XYZ(epv,cm,cell_ori,X_Ori,Y_Ori,Z_Ori,cell_new,&
     &                       xx_nw,yy_nw,zz_nw)
               DO ik = 1 , N
               ENDDO
               egymol = EGYmol0
               IF ( IROt>=1 ) THEN
                  DO ik = nv - nvr + 1 , nv
                     I_Nv = ik                                             ! 8-28-08 DU
                     ivr = ik - (nv-nvr)
                     DTHeta = epv(ik)
!             write(61,"('cell_ori :',6F16.6)") cell_ori                   ! 3-25-08 DU
!             write(61,"('cell_new :',6F16.6)") cell_new                   ! 3-25-08 DU
                     DO i = 1 , N
                        X(i) = xx_nw(i)
                        Y(i) = yy_nw(i)
                        Z(i) = zz_nw(i)
                     ENDDO
                     CALL BEND_BOND_NEW(ivr,N,cell_new,X,Y,Z,xo,yo,zo,  &
     &                                  egymol)
                     CALL FRACTIONAL_COD(N,cell_new,xo,yo,zo,xx_nw,     &
     &                  yy_nw,zz_nw)
                     DEL_egymol = egymol - EGYmol0
                     WRITE (61,                                         &
     &               "('intra-E for each bond : ',I4,F16.6,F9.3,3F16.6)"&
     &               ) ivr , epv(ik) , epv(ik)*57.29578 , EGYmol0 ,     &
     &                 egymol , DEL_egymol
                  ENDDO
                  DEL_egymol = egymol - EGYmol0
!--------- write out new cell and coords after bonds rotated 3-25-08 DU -
                  sx = 0.0
                  sy = 0.0
                  sz = 0.0
                  WRITE (61,99056) cell_new
                  DO ik = 1 , N
                     WRITE (61,99057) ATOm(ik) , ik , xx_nw(ik) ,       &
     &                                yy_nw(ik) , zz_nw(ik)
                     sx = sx + xx_nw(ik)
                     sy = sy + yy_nw(ik)
                     sz = sz + zz_nw(ik)
                  ENDDO
                  xcr = sx/N
                  ycr = sy/N
                  zcr = sz/N
                  ATOm(N+1) = 'XTRA  '
                  WRITE (61,99057) ATOm(N+1) , N - N , xcr , ycr , zcr
!-----------------------------------------------------------------------
               ELSE
                  DEL_egymol = 0.0
               ENDIF
               w_cor = 0.0
               IF ( IMOde>0 ) THEN
                  CALL USER_E_1(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,    &
     &                          e_mol,w_cor)
               ELSE
                  CALL POT_E(cell_new,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol, &
     &                       w_cor)
               ENDIF
               w_vt = DEL_egymol + w_cor
                                     ! egymol = 0.0 if irot = 0
               WRITE (61,99034) w_cor , DEL_egymol , w_cor + DEL_egymol
99034          FORMAT (/'Inter, intra & total energies =',3F11.6)
               WRITE (61,99035)
99035          FORMAT (17X,'sm',11X,'current E',5X,'   new E  ',11X,    &
     &                 'Ec',13X,'Ev',15X,'Er',11X,'inter E',5X,'nentry')
               WRITE (61,99036) sm , w , w_vt , ec , ev , er , w_cor ,  &
     &                          NENtry
99036          FORMAT (5X,'**',D16.4,6F16.6,I6/120('*'))
               IF ( w_vt<=w ) THEN
                  IF ( ABS(w_vt-w)<=w_diff ) EXIT                          ! 2-6-08
! ---------------------------------------------------------------------
                  IF ( IROt/=0 ) wd_lim = 5.0D0                            ! 4-13-08 DU
                  IF ( ABS(w_vt-w)>wd_lim .OR. iw>0.0 ) THEN
                     DO i = 1 , nv
                        pcitm(i) = psave(i)
                        epv(i) = epv_save(i)
                     ENDDO
                     DO i = 1 , N
                        xx(i) = xsave(i)
                        yy(i) = ysave(i)
                        zz(i) = zsave(i)
                     ENDDO
                     w = w_save
                     EXIT
                  ENDIF
!----------------------------------------------------------------------
!7158        do i = 1,6
                  DO i = 1 , 6
                     pcitm(i) = cell_new(i)
                  ENDDO
                  DO i = 1 , nv                                           ! 8-28-08  DU
                     pcitm(i) = new_pcitm(i)                              ! 8-28-08
                  ENDDO                                                   ! 8-28-08
                  DO i = 1 , N
                     xx(i) = xx_nw(i)
                     yy(i) = yy_nw(i)
                     zz(i) = zz_nw(i)
                     xsave(i) = xx(i)
                     ysave(i) = yy(i)
                     zsave(i) = zz(i)
                  ENDDO
                  w = w_vt
                  IF ( fa==1 ) EXIT
                  sm = 3.0*sm
                  su = 1
               ELSE
                  DO i = 1 , nv
                     pcitm(i) = psave(i)
                     epv(i) = epv_save(i)
                  ENDDO
                  DO i = 1 , N
                     xx(i) = xsave(i)
                     yy(i) = ysave(i)
                     zz(i) = zsave(i)
                  ENDDO
!           del_w = w_vt - w       ! deleted, 11/20/08, HLA
                  w = w_save
                  IF ( ABS(w_vt-w)<=w_diff ) EXIT
                  DO i = 1 , nv
                     pcitm(i) = psave(i)
                  ENDDO
                  DO i = 1 , N
                     xx(i) = xsave(i)
                     yy(i) = ysave(i)
                     zz(i) = zsave(i)
                  ENDDO
                  w = w_save
                  sm = -0.5*sm
                  fa = 1
                  IF ( su>=1 ) EXIT
               ENDIF
 
!
            ENDDO
!-------------------------------------------------------------------------------
! - here is ending the Vector Step Search (inside RSS calculation)             ]
!-------------------------------------------------------------------------------
 
            DO i = 1 , nv
               psave(i) = pcitm(i)
               epv_save(i) = epv(i)
               IF ( i<=6 ) cell_new(i) = pcitm(i)
            ENDDO
            IF ( IROt>=1 ) THEN
               CALL ORTHO_COD(N,psave,xx,yy,zz,xo,yo,zo)
               IF ( NE_type==1 ) THEN                                        ! 9-30-08 DU
                  CALL ENERGY_PM3(N,xo,yo,zo,egymol)                         ! this gets PM3 energy
               ELSE IF ( NE_type==2 ) THEN
                  CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)                       ! this gets b3lyp/631g** energy 9-30-08 DU
               ELSE
                  IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)     ! 4-28-09 DU
               ENDIF
            ENDIF
            w_cor = 0.0
            IF ( IMOde>0 ) THEN
               CALL USER_E_1(cell_new,xx,yy,zz,ec,ev,er,e_mol,w_1)
            ELSE
               CALL POT_E(cell_new,xx,yy,zz,ec,ev,er,e_mol,w_1)
            ENDIF
!
            DEL_egymol = egymol - EGYmol0
            w_2 = DEL_egymol + w_1      ! egymol = 0.0 if irot = 0
            w_save = w_2
            w = w_2
            DO i = 1 , nv
               epv(i) = epv_save(i)
               pcitm(i) = psave(i)
               DEL_param(i) = pqc(i)
            ENDDO
            DO i = 1 , N
               xx_nw(i) = xx(i)
               yy_nw(i) = yy(i)
               zz_nw(i) = zz(i)
            ENDDO
            WRITE (61,"('Vector search complete, go to next stage...')")
            WRITE (61,99037) ec , ev , er , w_1 , DEL_egymol , w_2
99037       FORMAT ('     **   Ec, Ev, Er, intra E, inter E, total E =',&
     &              6F13.6/103('*'))
            WRITE (61,99038)
99038       FORMAT (4X,'delta param',10X,'epv',14X,'pcitm')
            DO i = 1 , nv
               WRITE (61,99039) DEL_param(i) , epv(i) , pcitm(i)
99039          FORMAT (2D16.6,F16.6)
            ENDDO
!        if (rc ==  1 ) goto 8888                                            ! 8-19-08 DU
!        GO TO               8888                                            ! 8-19-08 DU
!           if (irot >= 1) then
!              do i = 1,n
!                 xr(i) = xx_nw(i)
!                 yr(i) = yy_nw(i)
!                 zr(i) = zz_nw(i)
!              enddo
!              do i = nv-nvr+1,nv
!                 ivr = i-(nv-nvr)
!                 dtheta = epv(i)
!                 CALL BEND_BOND_NEW(ivr,n,cell_new,xr,yr,zr,xo,yo,zo,egymol)
!                 CALL  FRACTIONAL_COD(n,cell_new,xo,yo,zo,xx_nw,yy_nw,zz_nw)
!                 WRITE(61,3435) i,epv(i),epv(i)*57.29578,egymol-egymol0
99040       FORMAT ('Angle =',I4,F9.6,F9.3,F16.6)
!                 do ik = 1,n
!                    xr(ik) = xx_nw(ik)
!                    yr(ik) = yy_nw(ik)
!                    zr(ik) = zz_nw(ik)
!                 enddo
!              enddo
!              del_egymol = egymol - egymol0
!             WRITE (61,3440) egymol, egymol0, del_egymol
99041       FORMAT ('PM3 intramolecular E, initial, difference:',3F16.6)
!
         ENDDO
         wcyc(NCY) = w
         WRITE (28,99042) NCY - 1 , wcyc(NCY) , (NCY-2) , wcyc(NCY-1) , &
     &                    (wcyc(NCY)-wcyc(NCY-1)) , vbsq
99042    FORMAT ('Wcyc(',I3,') =',F15.10,2x,'Wcyc(',I3,') =',F15.10,2x, &
     &           'Diff =',F15.10,2x,' VBSQ ',e16.8/)
99043    FORMAT ('Wcyc(',I3,') =',F15.10,2x,'Wcyc(',I3,') =',F15.10,2x, &
     &           'Diff =',F15.10,2x,' VBSQ ',e16.8,f12.6)
!
!      IF (abs(Wcyc(ncy) - Wcyc(ncy-1)) .le. 1.0D-12) then                   ! 4-8-08 DU
         IF ( ABS(wcyc(NCY)-wcyc(NCY-1))<=w_diff ) EXIT                      ! 6-30-08 DU
         IF ( vbsq<=1.0D-15 .AND. IROt==0 ) THEN                             ! 8-13-08 DU
            WRITE (61,99058) vbsq
            EXIT
         ENDIF
         IF ( IROt/=0 ) THEN                                                 ! 8-13-08 DU
            IF ( vbsq<=1.0D-15 .AND.                                    &
     &           (ABS(wcyc(NCY)-wcyc(NCY-1))<=w_diff) ) THEN                 ! 8-13-08 DU
               WRITE (61,99058) vbsq                                         ! 8-13-08 DU
               EXIT                                                          ! 8-13-08 DU
            ENDIF                                                            ! 8-13-08 DU
         ENDIF                                                               ! 8-13-08 DU
 
!      IF ( in_cyc == 1 ) GO TO 8063
         IF ( (ABS(wcyc(NCY)-wcyc(NCY-1))>5.0 .AND. (wcyc(NCY)>wcyc(NCY-&
     &        1)) .AND. NCY>3) .OR.                                     &
     &        (ABS(wcyc(NCY)-wcyc(NCY-1))>5.0 .AND. (wcyc(NCY)<wcyc(NCY-&    ! 3-11-09 DU
     &        1)) .AND. NCY>3) .OR.                                     &    ! 3-11-09 DU             
     &        (wcyc(NCY)>0.0 .AND. (wcyc(NCY)>wcyc(NCY-1))) ) THEN           ! 7-21-08 DU
            NCY = NCY - 1
            DO i = 1 , 6
               cell_new(i) = save_cell(i)
               AH(i) = save_cell(i)
               fpc(i) = save_cell(i)
            ENDDO
            WRITE (61,"(6F9.6)") (save_cell(i),i=1,6)     
            DO i = 1 , N
               xx(i) = save_x(i)
               yy(i) = save_y(i)
               zz(i) = save_z(i)
               Fx(i) = save_x(i)
               Fy(i) = save_y(i)
               Fz(i) = save_z(i)
               WRITE (61,"(3F12.5)") save_x(i) , save_y(i) , save_z(i)
            ENDDO
            kill = 1
            IF ( IMOde>0 ) THEN
               CALL USER_E_1(save_cell,save_x,save_y,save_z,ec,ev,er,   &
     &                       e_mol,w_cor)                                          ! 2-15-08
            ELSE
               CALL POT_E(save_cell,save_x,save_y,save_z,ec,ev,er,e_mol,&
     &                    w_cor)                                                   ! 2-15-08
               WRITE (61,"(I6,5F16.6)") NENtry , ec , ev , er , e_mol , &
     &                                  w_cor
            ENDIF
            EXIT
         ENDIF
         w_previous = w                                                            ! 6/11/07
         w = w_vt
         DO i = 1 , nv
            IF ( i<=6 ) buff(i) = pcitm(i)
         ENDDO
         IF ( IMOde>0 ) THEN
            CALL USER_E_1(buff,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol,w_cor)
         ELSE
            CALL POT_E(buff,xx_nw,yy_nw,zz_nw,ec,ev,er,e_mol,w_cor)
         ENDIF
         sum_epv = 0.0
         DO i = 1 , nv
            sum_epv = sum_epv + epv(i)
            IF ( i>6 ) pcitm(i) = epv(i)                                          ! 8-19-08 DU 
         ENDDO
!      WRITE (61,8041) ncy, w_cor
99044    FORMAT ('Go to next cycle...cycle =',I3,', E =',F13.6)
      ENDDO
!
!-------------------------------------------------------------------------------
! - here is ending the big loop for RSS calculation                            ]
!-------------------------------------------------------------------------------
!
      DO i = 1 , nv
         DEL_param(i) = 0.0
      ENDDO
 
!     DO i = 1 , N                     ! deleted, 11/20/08, HLA
!        difx(i) = Fx(i) - ix(i)       ! deleted, 11/20/08, HLA
!        dify(i) = Fy(i) - iy(i)       ! deleted, 11/20/08, HLA
!        difz(i) = Fz(i) - iz(i)       ! deleted, 11/20/08, HLA
!     ENDDO                            ! deleted, 11/20/08, HLA
!
      DO i = 1 , nv
         IF ( i>(nv-nvr) ) fpc(i) = sum_corr(i)
         DIFfpc(i) = fpc(i) - ipc(i)          
         IF ( ipc(i)<=0.000001 ) THEN
            PRDifpc(i) = DIFfpc(i)
         ELSE
            PRDifpc(i) = (DIFfpc(i)/ipc(i))*100
         ENDIF
      ENDDO
      DO i = 1 , nv
         DIFfpc(i) = fpc(i) - ipc(i)          
         IF ( ipc(i)<=0.000001 ) THEN
            PRDifpc(i) = DIFfpc(i)
         ELSE
            PRDifpc(i) = (DIFfpc(i)/ipc(i))*100
         ENDIF
      ENDDO
!
      istart = 1
      DO imol = 1 , NMOl
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = ix(i) + sumx
            sumy = iy(i) + sumy
            sumz = iz(i) + sumz
         ENDDO
! Centriod of atoms for atoms of the group imol xc,yc,zc
         xc(imol) = REAL(sumx/NATm(imol))
         yc(imol) = REAL(sumy/NATm(imol))
         zc(imol) = REAL(sumz/NATm(imol))
         WRITE (13,99045) imol , xc(imol) , yc(imol) , zc(imol)
99045    FORMAT (' Molecular center(s)...'/'  XTR',i1,3F10.5)
!        WRITE (12,99059) imol , xc(imol) , yc(imol) , zc(imol)           ! 12-16-08
         istart = istart + NATm(imol)
      ENDDO
!
! Final Energy Calculation after RSS
      IF ( kill==0 ) NENtry = 1
      IF ( IROt>=1 ) THEN
         CALL ORTHO_COD(N,psave,xx,yy,zz,xo,yo,zo)
         IF ( NE_type==1 ) THEN                                           ! 9-30-08 DU
            CALL ENERGY_PM3(N,xo,yo,zo,egymol)                            ! this gets PM3 energy
         ELSE IF ( NE_type==2 ) THEN
            CALL ENERGY_B3LYP(N,xo,yo,zo,egymol)                          ! this gets b3lyp/631g** energy 9-30-08 DU
         ELSE
            IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)        ! 4-28-09 DU
         ENDIF
      ENDIF
      w_cor = 0.0
      WRITE (61,"('final cell ',6F9.5)") (cell_new(i),i=1,6)              ! 2-15-08
      DO i = 1 , N
         WRITE (61,"(3F12.6)") xx(i) , yy(i) , zz(i)
      ENDDO
      IF ( IMOde>0 ) THEN
         CALL USER_E_1(cell_new,xx,yy,zz,ec,ev,er,e_mol,w_1)
      ELSE
         ncy_rss_end = ncy                                                ! 2-26-09
         CALL POT_E(cell_new,xx,yy,zz,ec,ev,er,e_mol,w_1)                 ! 2-26-09
!        if (short_dd > 0) then                                           ! short intermolecular dis at RSS final cycle
!          stop                                                           ! 2-26-09  5-4-10 
!        end if                                                           ! 2-26-09  5-4-10
         ncy_rss_end = 0                                                  ! 2-26-09  5-4-10
      ENDIF
      DEL_egymol = egymol - EGYmol0
      w = DEL_egymol + w_1       ! egymol = 0.0 if irot = 0
      WRITE (61,"(6F16.8)") ec , ev , er , DEL_egymol , w_1 , w           ! 2-15-08
!
!       del_egymol = egymol-egymol0
!       W = del_egymol  + W_1    ! egymol = 0.0 if irot = 0
      IF ( IROt>=1 ) THEN                                                 ! 3-31-08 DU
         DO i = nv - nvr + 1 , nv
            WRITE (61,99046) epv(i) , epv(i)*57.29578                     ! 3-31-08 DU
            DIFfpc(i) = epv(i)       
99046       FORMAT ('Angles : ',F9.6,2X,F9.3)
         ENDDO
      ENDIF
      WRITE (61,"('intra_E & inter_E & total E : ',32F16.8)")           &
     &       DEL_egymol , w_1 , w                                         ! 3-31-08 DU
      WRITE (61,99060) NCY - 1                                            ! 9-16-08 DU
      WRITE (61,99061)                                                    ! 9-16-08 DU
      WRITE (13,99060) NCY - 1
      WRITE (13,99061)
      WRITE (12,99060) NCY - 1                                            ! 9-30-08 DU
      WRITE (12,99061)                                                    ! 9-30-08 DU
      kmol = 0                                                            ! 12-22-08 DU
      rt = 0.0
      ts = 0.0
      DO i = 1 , nv
         IF ( IDX(i)/=0 ) THEN
!        write(8,818) i,ipc(i),diffpc(i),fpc(i),prdifpc(i),epv(i)         ! removed 7/11/08 HLA
            IF ( i>3 .AND. i<7 ) PRDifpc(i)                             &
     &           = ((ACOSD(fpc(i))-ACOSD(ipc(i)))/ACOSD(ipc(i)))*100                                                
            IF ( IROt/=0 .AND. i>(6+6*nmol) ) THEN                        ! 12-22-08 DU
               ii = INDX(i-(6+6*nmol)) + 1                                ! 1-22--9 DU 
               WRITE (13,99163) i,bond_type(ii),i-(6+6*nmol),epv(i),epv(i), 57.29578*epv(i)                     ! 1-22-09 DU
               WRITE (61,99163) i,bond_type(ii),i-(6+6*nmol),epv(i)*57.29578,epv(i)*57.29578, 57.29578*epv(i)   ! 1-22-09 DU
               WRITE (12,99163) i,bond_type(ii),i-(6+6*nmol),epv(i)*57.29578,epv(i)*57.29578, 57.29578*epv(i)   ! 1-22-09 DU
99163          FORMAT(3X,I3,1X,A10,I2,', deg',4X,'0.000000',F14.6,F14.6,F10.4)                                  ! 1-8-09  DU
            ELSE
             IF ( vab_name(i)(5:9) .EQ. 'rotn,' ) THEN
               WRITE (13,99079) i , vab_name(i), ipc(i) , DIFfpc(i)*57.29578 , fpc(i)*57.29578 , &      ! 1-8-08 DU               
     &                          PRDifpc(i)*57.29578                                                     ! 1-8-08 DU
               WRITE (61,99079) i , vab_name(i), ipc(i) , DIFfpc(i)*57.29578 , fpc(i)*57.29578 , &      ! 1-8-08 DU
     &                          PRDifpc(i)*57.29578                                                     ! 1-8-09 DU
               WRITE (12,99079) i , vab_name(i), ipc(i) , DIFfpc(i)*57.29578 , fpc(i)*57.29578 , &      ! 1-8-09 DU
     &                          PRDifpc(i)*57.29578                                                     ! 12-22-08 DU
             ELSE
               if ( i>3 .AND. i<7 ) then
               WRITE (13,99079) i , vab_name(i), ipc(i) , DIFfpc(i) , fpc(i) ,       &                  ! 1-8-09 DU
     &                          PRDifpc(i)                                                              ! 12-22-08 DU
               WRITE (61,99079) i , vab_name(i), ACOSD(ipc(i)), (ACOSD(fpc(i))-ACOSD(ipc(i))), ACOSD(fpc(i)),&  ! 1-8-09 DU 
     &                          PRDifpc(i)                                                                          
               WRITE (12,99079) i , vab_name(i), ACOSD(ipc(i)), (ACOSD(fpc(i))-ACOSD(ipc(i))), ACOSD(fpc(i)),&  ! 1-8-09 DU
     &                          PRDifpc(i)
               else
               WRITE (13,99079) i , vab_name(i), ipc(i) , DIFfpc(i) , fpc(i) ,       &                  ! 1-8-09 DU
     &                          PRDifpc(i)                                                              ! 12-22-08 DU
               WRITE (61,99079) i , vab_name(i), ipc(i) , DIFfpc(i) , fpc(i) ,       &                  ! 1-8-09 DU
     &                          PRDifpc(i)
               WRITE (12,99079) i , vab_name(i), ipc(i) , DIFfpc(i) , fpc(i) ,       &                  ! 1-8-09 DU
     &                          PRDifpc(i)
               end if
!     ts = ts + SQRT((PRDifpc(10+6*kmol)**2)+(PRDifpc(11+6*kmol)**2)+(PRDifpc(12+6*kmol)**2)) ! translation A
!    &                          PRDifpc(i) , epv(i)                          
99079             FORMAT (i6,1X,A15,2x,F12.6,2x,F12.6,2x,F12.6,2x,F8.4)
             ENDIF
            END IF
         ENDIF
      ENDDO
      rt = SQRT(((PRDifpc(7)*57.29578)**2)+((PRDifpc(8)*57.29578)**2)   &
     &     +((PRDifpc(9)*57.29578)**2))                                             ! rotation degree
      ts = SQRT((PRDifpc(10)**2)+(PRDifpc(11)**2)+(PRDifpc(12)**2))                 ! translation Angs
      f_factor = (0.5*rt)**2 + (10*ts) + PRDifpc(1)**2 + PRDifpc(2)     &
     &           **2 + PRDifpc(3)**2 + (ACOSD(fpc(4))-ACOSD(ipc(4)))    &
     &           **2 + (ACOSD(fpc(5))-ACOSD(ipc(5)))                    &
     &           **2 + (ACOSD(fpc(6))-ACOSD(ipc(6)))**2
      WRITE (13,"( 'F factor = ', F13.6)") f_factor                                 ! 9-16-08 DU
      f_factor_rss = f_factor                                                       ! 12-16-08 DU
      IF ( IROt>0 ) THEN
         WRITE (12,99065) iw , w , DEL_egymol , NCY
         WRITE (13,99065) iw , w , DEL_egymol , NCY
      ELSE
         WRITE (13,99064) NCY - 1 , iw , w , wcyc(NCY) - wcyc(NCY-1)
         WRITE (61,99064) NCY - 1 , iw , w , wcyc(NCY) - wcyc(NCY-1)                ! 9-16-08 DU
      ENDIF
!     WRITE (61,"(/12X,'#',9X,'E_conf',10X,'E_pack',10X,'E_toatl'/)")
      WRITE (61,"(/12X,'#',9X,'E_conf',10X,'E_toatl',10X,'E diff'/)")               ! 12-16-08
      DO ik = 1 , NCY                                                               ! 12-16-08
!        WRITE (61,"('CYCLE : ',I5,1X ,3F16.8)") ik - 1 , econf_cyc(ik) &
!    &          , epack_cyc(ik) , wcyc(ik)
         if (ik == 1) wcyc(ik-1) = wcyc(ik)
         if (ik == ncy) econf_cyc(ncy) = econf_cyc(ncy-1)                           ! 1-22-09
         WRITE (61,"('CYCLE : ',I5,1X ,2F16.8, 2X,E16.8)") ik - 1 , econf_cyc(ik) & ! 12-16-08
     &          , wcyc(ik), wcyc(ik)-wcyc(ik-1)
99047    FORMAT (/,5x,'cycle, E',i5,f16.8)
      ENDDO
!
!-----------------------------------------------------------------------------------------
!        - calculating density                                                           ]
!-----------------------------------------------------------------------------------------
!                                                                                        
      coef = DSQRT(1-(fpc(4))**2-(fpc(5))**2-(fpc(6))**2+2*fpc(4)*fpc(5)&
     &       *fpc(6))
      volume = fpc(1)*fpc(2)*fpc(3)*coef
      ro = (nmc*WT_mol)/(0.6023*volume)
      dens_per = ((ro-DENs_1st)/DENs_1st)*100                                       ! 9-16-08 DU
      vol_per = ((volume-VOL_1st)/VOL_1st)*100                                      ! 9-16-08 DU
      IF ( done ) WRITE (61,99068) VOL_1st , volume , vol_per                       ! 9-16-08 DU
      IF ( done ) WRITE (12,99068) VOL_1st , volume , vol_per                       ! 9-16-08 DU
      IF ( done ) WRITE (61,99066) DENs_1st , ro , dens_per                         ! 9-16-08 DU
      IF ( done ) WRITE (12,99066) DENs_1st , ro , dens_per
      WRITE (13,99067) DENs_1st , ro , dens_per                                     ! 9-16-08 DU
      WRITE (61,99067) DENs_1st , ro , dens_per                                     ! 9-16-08 DU
      WRITE (61,99068) VOL_1st , volume , vol_per                                   ! 9-16-08 DU
      WRITE (13,99068) VOL_1st , volume , vol_per                                   ! 9-16-08 DU
 
!
      OPEN (UNIT=21,POSITION='append')                                              ! for SUMRY file 6/11/07
      IF ( IRSs_call==0 ) WRITE (21,99048) VOL_1st , volume , DENs_1st ,&
     &                           ro , w                                             ! 4-30-08, HLA
99048 FORMAT (4x,'PMIN initial & final V, D, E =',2x,2F8.2,2F8.3,F9.2,  &
     &        ' MODE=  3')                                                          ! 6/11/07
      CLOSE (21,STATUS='keep')                                                      ! 6/11/07
!     if (irss_call == 0) print 1250, ncy-1,ro, w, (W_previous - w),vbsq            ! 4-30-08, HLA
!     print 1250, ncy-1,ro, w, (W_previous - w),vbsq                                ! 4-30-08, HLA
!1250 format ('PMIN Rosenbrock refn: cycles =',i3, &                                ! 6/11/07
!                '; density =',f6.3,'; Ecalc =',f11.6,'; final E diff =',E13.6,', vbsq =', D13.6)
      IF ( IROt>0 ) THEN                                                            ! 9-30-08 DU
         PRINT 99049 , NCY - 1 , ro , w_1 , DEL_egymol , w ,            &
     &         (w_previous-w)                                                       ! 9-30-08 DU
                                                                                    ! 9-30-08 DU
!        write(12, 99049) NCY - 1 , ro , w_1 , DEL_egymol , w ,         &
!    &         (w_previous-w)                                                       ! 9-30-08 DU
!        write(13, 99049) NCY - 1 , ro , w_1 , DEL_egymol , w ,         &
!    &         (w_previous-w)                                                       ! 9-30-08 DU
99049    FORMAT ('PMIN Rosenbrock refn: cycles =',i3,'; density =',f6.3,&
     &           '; Lat E =',f11.6,'; Con E =',f11.6,'; Tol E =',f11.6, &
     &           '; final E diff =',E13.6)                                          ! 9-30-08 DU
      ELSE
         PRINT 99050 , NCY - 1 , ro , w , (w_previous-w)                            ! 9-30-08 DU
!        write(12,99050) NCY - 1 , ro , w , (w_previous-w)                          ! 9-30-08 DU
!        write(13,99050) NCY - 1 , ro , w , (w_previous-w)                          ! 9-30-08 DU
99050    FORMAT ('PMIN Rosenbrock refn: cycles =',i3,'; density =',f6.3,&
     &           '; Lat E =',f11.6,'; final E diff =',E13.6)                        ! 9-30-08 DU
      ENDIF                                                                         ! 9-30-08 DU
!
!-----------------------------------------------------------------------------------------
      DENsity = ro
!     DO i = 1 , N                 ! deleted, 11/20/08, HLA
!        difx(i) = Fx(i) - ix(i)   ! deleted, 11/20/08, HLA
!        dify(i) = Fy(i) - iy(i)   ! deleted, 11/20/08, HLA
!        difz(i) = Fz(i) - iz(i)   ! deleted, 11/20/08, HLA
!     ENDDO                        ! deleted, 11/20/08, HLA
      E = w
      DO i = 1 , 6
         AH(i) = fpc(i)
         CELl1(i) = fpc(i)
      ENDDO
      DO i = 1 , 6
         IF ( i>=4 .AND. i<=6 ) THEN
            CELl(i) = ACOSD(CELl1(i))                                               ! 5/12/08, HLA, using acosd
         ELSE
            CELl(i) = CELl1(i)
         ENDIF
      ENDDO
      DO i = 1 , N
         X(i) = Fx(i)
         Y(i) = Fy(i)
         Z(i) = Fz(i)
      ENDDO
!
      DO i = 1 , N
         X_Ori(i) = Fx(i)
         Y_Ori(i) = Fy(i)
         Z_Ori(i) = Fz(i)
      ENDDO
!-----------------------------------------------------------------------------------------
!        - calculating the XTRA coordinates                                              ]
!-----------------------------------------------------------------------------------------
!
      istart = 1
      DO imol = 1 , NMOl
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         iend = istart + NATm(imol) - 1
         DO i = istart , iend
            sumx = Fx(i) + sumx
            sumy = Fy(i) + sumy
            sumz = Fz(i) + sumz
         ENDDO
! Centriod of atoms for atoms of the group imol xc,yc,zc
         xc(imol) = REAL(sumx/NATm(imol))
         yc(imol) = REAL(sumy/NATm(imol))
         zc(imol) = REAL(sumz/NATm(imol))
         istart = istart + NATm(imol)
      ENDDO
      DO i = 1 , N
         sum_x = sum_x + ix(i)
         sum_y = sum_y + iy(i)
         sum_z = sum_z + iz(i)
      ENDDO
      xtra_x = sum_x/N
      xtra_y = sum_y/N
      xtra_z = sum_z/N
      sum_x = 0.0
      sum_y = 0.0
      sum_z = 0.0
      DO i = 1 , N
         sum_x = sum_x + Fx(i)
         sum_y = sum_y + Fy(i)
         sum_z = sum_z + Fz(i)
      ENDDO
      xtra_x = sum_x/N
      xtra_y = sum_y/N
      xtra_z = sum_z/N
!-----------------------------------------------------------------------------------------
      WRITE (12,99051)
99051 FORMAT (/                                                         &
     &     ' New cell params & coordinates from RSS in wmin file format'&
     &     ,' follow...'/5x,'(file name = pmin_RSS.save)'/)
      WRITE (12,99069) AH
      WRITE (13,*)                                                      &
     &    '    The NEW   cell parameters and coordinates after RSS are:'
      WRITE (13,99069) AH
      DO i = 1 , N
99052    FORMAT (5x,a6,i5,10F10.5)
         WRITE (12,99070) ATOm(i) , i , Fx(i) , Fy(i) , Fz(i)
         WRITE (13,99070) ATOm(i) , i , Fx(i) , Fy(i) , Fz(i)
      ENDDO
!
      istop = 0
      IF ( istop/=1 ) THEN
!
!-----------------------------------------------------------------------------------------
!        - writing the WMIN type input file after RSS
!        CLOSE (UNIT=9)
!        OPEN (UNIT=9,FILE='wmin.input',STATUS='OLD')
!        OPEN (UNIT=99,FILE='wmin.input_RSS',STATUS='UNKNOWN')
         OPEN (UNIT=22,FILE='pmin_RSS.save',STATUS='UNKNOWN')                  
 
!        DO i = 1 , 5
!           READ (9,99071) line
!           WRITE (99,99071) line
!        ENDDO
!        READ (9,99071) line
!        WRITE (99,99071) line
!        READ (line,'(3x,2I3)') nka , ns
!        DO i = 1 , 4
!           READ (9,99071) line
!           WRITE (99,99071) line
!        ENDDO
!        READ (9,99071) line
!        READ (line,'(54x,A26)') cline1
!        WRITE (99,99053) AH , cline1
!99053   FORMAT (3F9.4,3F9.5,A26)
!        READ (9,99071) line
!        WRITE (99,99071) line
!        DO i = 1 , ns
!           READ (9,99071) line
!           WRITE (99,99071) line
!        ENDDO
!        DO i = 1 , nka
!           READ (9,99071) line
!           WRITE (99,99071) line
!        ENDDO
!        DO i = 1 , nka
!           READ (9,99072) line1
!           WRITE (99,99054) line1 , Fx(i) , Fy(i) , Fz(i) , cline1
!99054      FORMAT (A27,3F9.5,A26)
!        ENDDO
!
         DO imol = 1 , NMOl
            WRITE (12,99059) imol , xc(imol) , yc(imol) , zc(imol)
            WRITE (13,99059) imol , xc(imol) , yc(imol) , zc(imol)
!           WRITE (99,99059) imol , xc(imol) , yc(imol) , zc(imol)
         ENDDO
!
!        READ (9,99072) line1
!99055   FORMAT (A27,3F9.5,A26)
!        DO i = 1 , 6
!           READ (9,99071) line
!           WRITE (99,99071) line
!        ENDDO
         WRITE (22,"(6F9.5)") AH
         DO i = 1 , N
            WRITE (22,"(A6,I3,18X,3F9.5)") ATOm(i) , i , Fx(i) , Fy(i) ,&
     &             Fz(i)
         ENDDO
         DO imol = 1 , NMOl
            WRITE (22,"('XTRA    0',18X,3F9.5)") xc(imol) , yc(imol) , zc(imol)
         ENDDO
!        WRITE (22,"('XTRA    0',18X,3F9.5)") xtra_x , xtra_y , xtra_z
      ENDIF
!
      IF ( IROt/=0 ) THEN                                             ! 9-30-08 DU
         WRITE (62,"('conf. E before RSS: ',D16.8)") EGYmol0          ! 9-30-08 DU
         EGYmol0 = egymol                                             ! 9-30-08 DU new conf. E after RSS as initial conf. E for
         WRITE (62,"('conf. E and Angles after RSS: ',D16.8,10F8.3/)") &   ! 4-28-09
     &          EGYmol0 , (epv(i)*57.29578,i=nv-nvr+1,nv)
      ENDIF                                                           ! 9-30-08 DU
      CLOSE (UNIT=22)
!
99056 FORMAT (6F9.5)
99057 FORMAT (a6,i3,18x,3F9.5)
99058 FORMAT (/'VBSQ effectively NULL --> STOP...vbsq =',D15.8)
99059 FORMAT ('XTR',i1,5x,18x,3F9.5)
!       write(8,8217) ncy-1                                                            ! removed 7/11/08 HLA
99060 FORMAT (/' Initial and new parameters after RSS, # cycles =',i3,  &
     &        '...'/)
!       write(8,8215)                                                                  ! removed 7/11/08 HLA
99061 FORMAT (5X,                                                       &
     &      '#  Parameter            initial       change        new P      diff %'/)  ! 1-8-09  DU
99062 FORMAT (i6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F8.4,2x,f10.4,2X,f10.6)
99063 FORMAT (i6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F8.4,2x,f10.4,2X,f10.6,  &
     &        '(',f10.6,')')                                                           ! 9-30-08 DU
99064 FORMAT (/' RSS refinement',i3,' cycles,',                         &
     &        ' initial, final and E differences =',2(f12.6,','),e13.6)
99065 FORMAT (/,'Initial Energy before RSS =',F13.8,/,                  &
     &          'Final   Energy after  RSS =',F13.8,2x,                 &
     &        '  del_egymol =  ',f13.8,2x,'--  calculated in ncy =',I3, &
     &        ' cycles')
99066 FORMAT (/' Initial and new densities diff % =',f8.4,',',f8.4,',', &
     &        f8.4)                                                                    ! 9-16-08 DU
99067 FORMAT (/'Initial, final and density differences =',2(f10.4,','),f9.4,'% g/cc') 
99068 FORMAT ( 'Initial, final and volume differences  =',2(f10.4,','),f9.4,'% Angs**3')
99069 FORMAT (6F9.5)
99070 FORMAT (a6,i3,18x,3F9.5)
99071 FORMAT (A80)
99072 FORMAT (A27)
!
      END SUBROUTINE ROSENBROCK
!
!end--------------------subroutine rosenbrock-----------------------------end
!
      SUBROUTINE NEW_XYZ_TRIGONAL(Icode,N,A,C1,Dp,Atom,X,Y,Z,Cell,Xyz)
!
      USE PMIN_MODULE , ONLY:MATM , PI                         
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 :: A , C1 , Dp
      INTEGER :: Icode , N
      CHARACTER(6) , DIMENSION(MATM) :: Atom
      REAL*8 , DIMENSION(6) :: Cell
      REAL*8 , DIMENSION(MATM) :: X , Y , Z
      REAL*8 , DIMENSION(MATM,3) :: Xyz
      INTENT (IN) Atom , C1 , Dp , Icode , N , X , Y , Z
      INTENT (INOUT) A , Cell , Xyz
!
! Local variables
!
      REAL*8 :: b , c , c2 , c3 , delx1 , delx2 , delx3 , dely1 ,       &
     &          dely2 , dely3 , delz1 , delz2 , delz3 , del_d ,         &
     &          del_param , del_r , dist , dl , dl1 , dm , dm1 , dn ,   &
     &          dn1 , dx , dy , dz , r , r1 , r12 , r2 , s1 , s3 , xc2 ,&
     &          xco , yc2 , yco , zc2 , zco
      REAL*8 :: DSQRT
      INTEGER :: i , i1 , ij , ik , itrigonal , j , k
      REAL*8 , DIMENSION(3) :: x1 , y1 , z1
!
!
!       open (unit = 7,file = 'new_xyz.inp', status = 'old')
!       open (unit = 8,file = 'trigonal.out', status = 'unknown')
!       read(7,*) (x1(i),y1(i),z1(i),i=1,3)
!       write(8,2) (x1(i),y1(i),z1(i),i=1,3)
99001 FORMAT (5x,'Cell defining coordinates are :',/5x,3F10.5/5x,       &
     &        3F10.5/5x,3F10.5,/)
!      go to 10
!
      itrigonal = 1
!
!       ortho coordinates for cell vectors for trigonal case only
!
      del_param = Dp
      IF ( Icode==1 ) A = A + Dp
      b = A
      c = A
      c2 = C1
      c3 = C1
      s1 = SQRT(1-C1*C1)
!     s2 = s1           ! deleted, 11/18/08, HLA
      s3 = s1
!           do i = 1,n
!       x(i) = x5(i) + y5(i)*c3+z5(i)*c2
!       y(i) = y5(i)*s3 + z5(i)*(c1 - c2*c3)/s3
!       z(i) = z5(i)*dsqrt(s1*s1 -((c2-c1*c3)/s3)**2)
!           enddo
 
      x1(1) = A         ! for [100]
      y1(1) = 0.0
      z1(1) = 0.0
!
      x1(2) = b*c3      ! for [010]
      y1(2) = b*s3
      z1(2) = 0.0
!
      x1(3) = c*c2      ! for [001]
      y1(3) = c*(C1-c2*c3)/s3
      z1(3) = c*DSQRT(s1*s1-((c2-C1*c3)/s3)**2)
      IF ( Icode/=1 ) THEN
!
         xc2 = 0.5*A     ! for 1/2,1/2,1/2
         yc2 = 0.5*b
         zc2 = 0.5*c
         xco = xc2 + yc2*c3 + zc2*c2
         yco = yc2*s3 + zc2*(C1-c2*c3)/s3
         zco = zc2*DSQRT(s1*s1-((c2-C1*c3)/s3)**2)
!
         r1 = A*SQRT(2*(1-C1))
         r2 = A*SQRT(2*(1-C1-del_param))
         del_r = r2 - r1
!
         delx1 = xco - x1(1)
         dely1 = yco - y1(1)
         delz1 = zco - z1(1)
         dist = SQRT(delx1*delx1+dely1*dely1+delz1*delz1)
         del_d = dist*(del_r/r1)
         x1(1) = x1(1) + del_d*delx1/dist
         y1(1) = y1(1) + del_d*dely1/dist
         z1(1) = z1(1) + del_d*delz1/dist
!
         delx2 = xco - x1(2)
         dely2 = yco - y1(2)
         delz2 = zco - z1(2)
         x1(2) = x1(2) + del_d*delx2/dist
         y1(2) = y1(2) + del_d*dely2/dist
         z1(2) = z1(2) + del_d*delz2/dist
!
         delx3 = xco - x1(3)
         dely3 = yco - y1(3)
         delz3 = zco - z1(3)
         x1(3) = x1(3) + del_d*delx3/dist
         y1(3) = y1(3) + del_d*dely3/dist
         z1(3) = z1(3) + del_d*delz3/dist
      ENDIF
!
!
      DO k = 1 , 3
         i = k
         j = k + 1
         ij = k + 2
         IF ( j==4 ) j = 1
         IF ( ij>3 ) ij = MOD(ij,3)
!      write(8,4)k,i,j,ij
99002    FORMAT (4I6)
!
!      coefficients (direction ratios) defining plane through (0,0,0),(x1,y1,z1),(x2,y2.z2)
!      Equation of the plane is
!      dL*x + dm*y + dn*z = 0 where
!
         dl = y1(i)*z1(j) - y1(j)*z1(i)
         dm = z1(i)*x1(j) - z1(j)*x1(i)
         dn = x1(i)*y1(j) - x1(j)*y1(i)
         r = SQRT(x1(ij)*x1(ij)+y1(ij)*y1(ij)+z1(ij)*z1(ij))
         Cell(ij) = r
!     direction cosines of the third axis
         dl1 = x1(ij)/r
         dm1 = y1(ij)/r
         dn1 = z1(ij)/r
!
         DO i1 = 1 , N
            Xyz(i1,ij) = (dl*X(i1)+dm*Y(i1)+dn*Z(i1))                   &
     &                   /(dl*dl1+dm*dm1+dn*dn1)
            IF ( itrigonal==1 ) THEN
               Xyz(i1,ij) = Xyz(i1,ij)/A
            ELSE
               Xyz(i1,ij) = Xyz(i1,ij)/r
            ENDIF
         ENDDO
      ENDDO
!      write(8,*) ' '
!
! Get cell angles
!
      DO i = 1 , 3
         j = i + 1
         IF ( j==4 ) j = 1
         ij = MOD(i+4,3) + 4
         dx = x1(i) - x1(j)
         dy = y1(i) - y1(j)
         dz = z1(i) - z1(j)
         r12 = SQRT(dx*dx+dy*dy+dz*dz)
         Cell(ij) = ACOS((-r12*r12+Cell(i)*Cell(i)+Cell(j)*Cell(j))     &
     &              /(2*Cell(i)*Cell(j)))
         Cell(ij) = Cell(ij)*180.0/PI
         if (ij > 3) Cell(ij) = cosd(Cell(ij))                          
!        write(8,8)i,j,ij,r12,cell(i),cell(j),cell(ij)
 8       format(5x,'i,j,ij,r12,cell(i),cell(j),cell(ij)',3i5,4f10.4,/)
      ENDDO
      IF ( itrigonal==1 ) THEN
         Cell(1) = A
         Cell(2) = A
         Cell(3) = A
      ENDIF
      WRITE (8,99003) (Cell(i),i=1,6)
99003 FORMAT (5x,' Cell from trigonal subroutine: ',6F16.10,/)
      DO i = 1 , N
!       write(8,1)atom(i),x5(i),y5(i),z5(i),x(i),y(i),z(i),xyz(i,1),xyz(i,2),xyz(i,3)
         WRITE (18,99004) Atom(i) , X(i) , Y(i) , Z(i) ,                &
     &                    (Xyz(i,ik),ik=1,3)
!
!10    n = n + 1
!      read(7,1,end=100) atom(n),x(n),y(n),z(n)
!      write(8,1) atom(n),x(n),y(n),z(n)
99004    FORMAT (5x,a6,3F9.5,6F15.10)
      ENDDO
      RETURN
      END SUBROUTINE NEW_XYZ_TRIGONAL
!*==unique_sym.spg  processed by SPAG 6.56Rc at 11:22 on 24 Nov 2008
!
!end-------------------------new_xyz_trigonal---------------------------------end
!
      SUBROUTINE UNIQUE_SYM(Sym,X,Y,Z,Nsym1)
!
      USE PMIN_MODULE , ONLY:NMOl , NATm , IDX , ISYstem , N , ILS ,    &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , MATM , ATOm                    
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER , DIMENSION(20) :: Nsym1
      REAL*8 , DIMENSION(20,12) :: Sym
      REAL*8 , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN) Sym , X , Y , Z
      INTENT (INOUT) Nsym1
!
! Local variables
!
      REAL*8 :: dx , dy , dz , rij , x4 , x5 , y4 , y5 , z4 , z5
      INTEGER :: i , i1 , j , j1
      INTEGER , DIMENSION(NSYm,NSYm) :: ksym
!
!
      WRITE (8,99001) ATOm
99001 FORMAT (100(1x,a6))
      WRITE (8,99002) NATm , ISYstem , N , NSYm
99002 FORMAT (/,'natm,isystem,n,nsym ',20I5,/)
      DO i = 1 , NSYm
         Nsym1(i) = 1
      ENDDO
!
      DO j = 1 , NSYm
         DO j1 = 1 , NSYm
            ksym(j,j1) = 0
         ENDDO
      ENDDO
      DO i = 1 , N
!      Print 212,atom(i),x(i),y(i),z(i)       
99003    FORMAT (a6,3x,3F9.5)
      ENDDO
! Loop over all the atoms
!
      DO j = 1 , NSYm
         DO i = 1 , N
            x4 = X(i)*Sym(j,2) + Y(i)*Sym(j,3) + Z(i)*Sym(j,4)          &
     &           + Sym(j,1)
            y4 = X(i)*Sym(j,6) + Y(i)*Sym(j,7) + Z(i)*Sym(j,8)          &
     &           + Sym(j,5)
            z4 = X(i)*Sym(j,10) + Y(i)*Sym(j,11) + Z(i)*Sym(j,12)       &
     &           + Sym(j,9)
!
            DO j1 = 1 , NSYm
               IF ( j1/=j ) THEN
                  DO i1 = 1 , N
                     IF ( i/=i1 ) THEN
                        x5 = X(i1)*Sym(j1,2) + Y(i1)*Sym(j1,3) + Z(i1)  &
     &                       *Sym(j1,4) + Sym(j,1)
                        dx = x5 - x4
                        y5 = X(i1)*Sym(j1,6) + Y(i1)*Sym(j1,7) + Z(i1)  &
     &                       *Sym(j1,8) + Sym(j,5)
                        dy = y5 - y4
                        z5 = X(i1)*Sym(j1,10) + Y(i1)*Sym(j1,11) + Z(i1)&
     &                       *Sym(j1,12) + Sym(j,9)
                        dz = z5 - z4
                        rij = SQRT(dx*dx+dy*dy+dz*dz)
!
!     Print 210, atom(i),atom(i1),j,j1,rij
!       if(i /= i1 .and. j /= j1 .and. rij == 0.0 )then
                        IF ( rij<=0.000001 ) THEN
!        write(8,210) atom(i),atom(i1),j,j1,rij
!        Print 210,atom(i),atom(i1),j,j1,rij
                           ksym(j,j1) = ksym(j,j1) + 1
99004                      FORMAT (5x,a6,3x,a6,2I5,f9.5)
!         k1 = j
!        nsym1(j) = 1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      DO j = 1 , NSYm
         DO j1 = 1 , NSYm
            IF ( ksym(j,j1)/=0 ) THEN
!      Print 211,j,j1,ksym(j,j1)               
               WRITE (8,99005) j , j1 , ksym(j,j1)
99005          FORMAT (/5x,'j,j1,ksym(j,j1) ',3I6)
               IF ( j1>j ) Nsym1(j1) = 0
               IF ( j>j1 ) Nsym1(j) = 0
            ENDIF
         ENDDO
      ENDDO
      WRITE (8,99006) (Nsym1(i),i=1,NSYm)
99006 FORMAT (20I1)
      END SUBROUTINE UNIQUE_SYM
!
!end-------------------------unique_sym----------------------------------end
!
      SUBROUTINE MULTIPLICITY(Sym,Atwt,X,Y,Z,Mcity,Wt_mol_asym)
!         
      USE PMIN_MODULE , ONLY:NMOl , NATm , IDX , ISYstem , N , ILS ,    &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , MATM , ATOm                    
      USE F77KINDS                        
      IMPLICIT NONE
!*--MULTIPLICITY6074
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      REAL*8 :: Wt_mol_asym
      REAL*8 , DIMENSION(MATM) :: Atwt , X , Y , Z
      INTEGER , DIMENSION(MATM) :: Mcity
      REAL*8 , DIMENSION(20,12) :: Sym
      INTENT (IN) Atwt , Sym , X , Y , Z
      INTENT (INOUT) Mcity , Wt_mol_asym
!
! Local variables
!
      REAL*8 :: dx , dy , dz , rij , x4 , y4 , z4
      REAL :: FLOAT
      INTEGER :: i , j
!
!*** End of declarations rewritten by SPAG
!
!
!       Print *,'MULTIPLICITY entered',n,nsym
!
      Wt_mol_asym = 0.0
      DO j = 1 , N
         Mcity(j) = 0
      ENDDO
!
!       do i = 1,n
!       Print 212,atom(i),x(i),y(i),z(i)
99001 FORMAT (a6,3x,3F9.5)
!        enddo
! Loop over all the atoms
!
      DO i = 1 , N
         DO j = 1 , NSYm
            x4 = X(i)*Sym(j,2) + Y(i)*Sym(j,3) + Z(i)*Sym(j,4)          &
     &           + Sym(j,1)
            y4 = X(i)*Sym(j,6) + Y(i)*Sym(j,7) + Z(i)*Sym(j,8)          &
     &           + Sym(j,5)
            z4 = X(i)*Sym(j,10) + Y(i)*Sym(j,11) + Z(i)*Sym(j,12)       &
     &           + Sym(j,9)
!
            dx = X(i) - x4
            dy = Y(i) - y4
            dz = Z(i) - z4
!
            rij = SQRT(dx*dx+dy*dy+dz*dz)
!
            IF ( rij<=0.000001 ) Mcity(i) = Mcity(i) + 1
!        write(8,210) i,atom(i),j,rij,mcity(i)
!        Print 210,i,atom(i),j,rij,mcity(i)
99002       FORMAT (5x,'i,atom,j,rij,mcity  ',i5,3x,a6,3x,i5,f9.5,i5)
         ENDDO
      ENDDO
      DO i = 1 , N
         Wt_mol_asym = Wt_mol_asym + Atwt(i)/(FLOAT(Mcity(i)))
!         write(8,213)atom(i),x(i),y(i),z(i),atwt(i),mcity(i)
99003    FORMAT ('From MULTIPLICITY ',a6,3x,4F9.5,i5)
      ENDDO
!        write(8,211)wt_mol_asym                                          ! removed 7/11/08 HLA
!211      format(/,5x,'Molecular weight of the asymmetric unit ',f9.4,/)  ! removed 7/11/08 HLA
      END SUBROUTINE MULTIPLICITY
!*==bend_bond_new.spg  processed by SPAG 6.56Rc at 11:22 on 24 Nov 2008
 
!end------------------------multiplicity--------------------------------end
!
!END MULTIPLICITY
!
      SUBROUTINE BEND_BOND_NEW(Ivr,N,Cell,X,Y,Z,Xo,Yo,Zo,Egymol)
!
!  Program to rotate a bond RS , bend and turn  a bond RS out of the plane PRQ.
!  P(x1,y1,z1), Q(x2,y2,z2), R(x3,y3,z3), S(x4,y4,z4)
!  Orthogonal coordinates assumed.  08. 31. 04
!
      USE PMIN_MODULE , ONLY:MATM , ATOm , PI , NMOl                                      
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE RSSPM3_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!*--BEND_BOND_NEW6157
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      REAL*8 :: Egymol
      INTEGER :: Ivr , N
      REAL*8 , DIMENSION(6) :: Cell
      REAL*8 , DIMENSION(MATM) :: X , Xo , Y , Yo , Z , Zo
      INTENT (IN) Cell , Ivr , X , Y , Z
      INTENT (INOUT) Xo , Yo , Zo
!
! Local variables
!
      REAL*8 :: alpha , beta , c1 , c2 , c3 , dist_in , dist_out , dx , &
     &          dx1 , dy , dy1 , dz , dz1 , gamma , gamma1 , gamma2 ,   &
     &          k_dist , s1 , s3 , x1 , x2 , x21 , x22 , x3 , x31 ,     &
     &          x32 , x33 , x34 , x340 , x341 , x342 , x4 , x5 , y1 ,   &
     &          y2 , y21 , y3 , y31 , y32 , y33 , y34 , y340 , y341 ,   &
     &          y342 , y4 , y5 , z1 , z2 , z21 , z22 , z23 , z3 , z31 , &
     &          z32 , z33 , z34 , z340 , z341 , z342 , z4 , z5
      REAL*8 :: DSQRT
      INTEGER :: i , ik , jk , k , k2 , kdx , nobond , nv , nvr
      INTEGER , DIMENSION(4) :: k1
      INTEGER , DIMENSION(100) :: linkatm
      REAL*8 , DIMENSION(MATM) :: x1o , y1o , z1o
      REAL*8 :: dXAB,dYAB,dZAB, dXBC,dYBC,dZBC, DISAB_in, DISBC_in, CA_in, CA_out, & ! 4-28-09 DU temp
     &          dX1AB,dY1AB,dZ1AB, dX1BC,dY1BC,dZ1BC, DISAB_out, DISBC_out, &        ! 4-28-09 DU temp
     &          org_angle, new_angle                                                 ! 4-28-09 DU temp
      INTEGER :: IDIS_ab                                                             ! 4-28-09 DU temp
!
!*** End of declarations rewritten by SPAG
!
!
!        write(8,*)'From BEND_BOND ivr,irot,nrotbond',ivr,irot,nrotbond
!        write(8,454) dtheta
99001 FORMAT (/,5x,'From BEND_BOND_NEW dtheta ',f10.6,/)
                               ! 8-28-08 DU
!        write(8,*)'k5,linkatom',k5,linkatom
!        write(8,*)' indx,iatom',indx,iatom
!
!        open(unit = 7, file = 'bend_bond_new.inp', status = 'old')
!        open(unit = 8, file = 'bend_bond_new.out', status = 'unknown')
!        open(unit = 10, file = 'bend_bond_new.xyz', status = 'old')
!
!         index = 0 to rotate about bond 1-2
!               = 1 to bend the bond 3 about an axis normal to bond 3 and in the plane of bonds 1-2.
!               = 2 to bend the bond 3 about an axis normal to the plane of bond 1-2
!               = 3 to bend the bond 3 about an axis normal to bond 3 and nearly normal to plane of bonds 1-2.
!
!        read(7,*)nrotbond,dtheta
!        read(7,*) index,(k1(i),i=1,4)    ! No. of 4 defining atoms, bond 3-4 to be rotated, bent or turned
!
! Now read the atom coordinates from unit 10
!
!         pi = 3.141592653589793D0
!         read(10,*)cell
      nvr = NROtbond                                                ! 8-28-08 DU
      nv = 6 + 6*NMOl + nvr                                         ! 8-28-08 DU
!        if (i_nv > nv-nvr) write(8,2231)cell                       ! 9-2-08  DU
99002 FORMAT (5x,'cell  ',6F9.5)
      c1 = Cell(4)
      c2 = Cell(5)
      c3 = Cell(6)
      s1 = SQRT(1-c1*c1)
!     s2 = SQRT(1-c2*c2)                                            ! deleted, 11/16/08, HLA
      s3 = SQRT(1-c3*c3)
!     r1 = (c1-c2*c3)/s3                                            ! deleted, 11/16/08, HLA
!     r2 = SQRT(s1*s1-((c2-c1*c3)/s3)**2)                           ! deleted, 11/16/08, HLA
!
!        n= 0
!10      n = n + 1
!        read(10,1,end=100) atom(n),xi(n),yi(n),zi(n)
!        write(8,1)atom(n),xi(n),yi(n),zi(n)
!           xout(n) = xi(n)
!           yout(n) = yi(n)
!           zout(n) = zi(n)
99003 FORMAT (2x,a6,3F12.6)
!           go to 10
!100     n = n-1
!
      DO i = 1 , N
!        xout(i) = X(i)                                             ! deleted, 11/16/08, HLA
!        yout(i) = Y(i)                                             ! deleted, 11/16/08, HLA
!        zout(i) = Z(i)                                             ! deleted, 11/16/08, HLA
      ENDDO
!
!        n1 = 1
!        n2 = 1
!        angle = 0.0
!
!      CALL ENERGY_pm3(n)
!        write(8,501)egymol
99004 FORMAT (/,'Intramol energy from ENERGY_pm3',f20.15,/)
!
!        get  orthogonal coordinates
      DO i = 1 , N
         Xo(i) = Cell(1)*X(i) + Cell(2)*Y(i)*c3 + Cell(3)*Z(i)*c2
         Yo(i) = Cell(2)*Y(i)*s3 + Cell(3)*Z(i)*(c1-c2*c3)/s3
         Zo(i) = Cell(3)*Z(i)*DSQRT(s1*s1-((c2-c1*c3)/s3)**2)
         x1o(i) = Xo(i)
         y1o(i) = Yo(i)
         z1o(i) = Zo(i)
 
!       write(8,5112)atom(i),xi(i),yi(i),zi(i),x(i),y(i),z(i)
99005    FORMAT (a6,9F10.6)
      ENDDO
!
!        do 700 i = 1,nrotbond
!
      DO i = 1 , 4
         k1(i) = IATom(Ivr,i)
      ENDDO
      k2 = K5(Ivr)
      DO i = 1 , k2
         linkatm(i) = LINkatom(Ivr,i)
      ENDDO
!        write(8,*)'From BEND_BOND  k1,k2,linkatm',k1,k2,(linkatm(i),i=1,k2)
!
      gamma1 = DTHeta
!       write(8,121)gamma1
99006 FORMAT (/,5x,' dtheta ',f9.5,/)
!        gamma = (gamma1*pi)/180.0  ! for angle in degrees
      gamma = gamma1                ! for angle in radian
!
      kdx = INDx(Ivr)
      IF ( kdx==0 ) THEN
         x1 = Xo(k1(3))
         y1 = Yo(k1(3))
         z1 = Zo(k1(3))
         x2 = Xo(k1(4))
         y2 = Yo(k1(4))
         z2 = Zo(k1(4))
         GOTO 100
      ENDIF
!
      x1 = Xo(k1(1))
      y1 = Yo(k1(1))
      z1 = Zo(k1(1))
      x2 = Xo(k1(2))
      y2 = Yo(k1(2))
      z2 = Zo(k1(2))
      x3 = Xo(k1(3))
      y3 = Yo(k1(3))
      z3 = Zo(k1(3))
!
      x4 = Xo(k1(4))
      y4 = Yo(k1(4))
      z4 = Zo(k1(4))
!
!        write(8,*) ' Before Call to Cross_Product '
!        write(8,21)x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
99007 FORMAT (10x,'INPUT COORDINATES ARE ',/,4(10x,3F9.5/))
!
      CALL CROSS_PRODUCT(kdx,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5, &
     &                   z5)
!
      x1 = x5
      y1 = y5
      z1 = z5
      x2 = x3
      y2 = y3
      z2 = z3
!       write(8,51)x1,y1,z1
99008 FORMAT (/5x,'Coordinates of the atom positon normal to the bond ',&
     &        3F9.5,/)
!        endif
!
!
 100  DO k = 1 , k2
         x3 = Xo(linkatm(k))
         y3 = Yo(linkatm(k))
         z3 = Zo(linkatm(k))
!
!        write(8,2)k,x1,y1,z1,x2,y2,z2,x3,y3,z3,gamma
99009    FORMAT (10x,'INPUT COORDINATES  for bending ',/,i5,            &
     &           3(10x,3F9.5/),10x,'Rotation angle = ',2x,f8.3/)
         x21 = x2 - x1
         y21 = y2 - y1
         z21 = z2 - z1
         x31 = x3 - x1
         y31 = y3 - y1
         z31 = z3 - z1
!        IF ( ABS(z21)<=0.000001 ) THEN                             ! 4-28-09
         IF ( z21 == 0.0 ) THEN                                     ! 4-28-09
            alpha = 0.5*3.141592653589793D0
!           if ((ncy-1) == 1 .and. rc ==1 .and.i_nv > 12) &         ! 4-28-09
!    &      write(63,"('alpha & z21 : ',I4,F9.3, F20.10)") i_nv, alpha*57.29578, z21
            GOTO 150
         ENDIF
         alpha = ATAN(y21/z21)
 150     x22 = x21
!        y22 = y21*COS(alpha) - z21*SIN(alpha)                      ! deleted, 11/20/08, HLA
         z22 = y21*SIN(alpha) + z21*COS(alpha)
         x32 = x31
         y32 = y31*COS(alpha) - z31*SIN(alpha)
         z32 = y31*SIN(alpha) + z31*COS(alpha)
!       write(8,3) x22,y22,z22,x32,y32,z32
99010    FORMAT (/,'after rotation by alpha about x-axis ',/,           &
     &           2(10x,3F9.5/))
!
!        IF ( ABS(z22)<=0.000001 ) THEN                            ! 4-28-09
         IF ( z22 == 0.0 ) THEN                                    ! 4-28-09
            beta = 0.5*3.141592653589793D0
!           if ((ncy-1) == 1 .and. rc ==1 .and.i_nv > 12) &         ! 4-28-09
!    &      write(63,"('beta & z22 : ',I4,F9.3, F20.10)") i_nv, beta*57.29578, z22
            GOTO 200
         ENDIF
         beta = -ATAN(x22/z22)
!        x23 = z22*SIN(beta) + x22*COS(beta)
!        y23 = y22
 200     z23 = z22*COS(beta) - x22*SIN(beta)
         x33 = z32*SIN(beta) + x32*COS(beta)
         y33 = y32
         z33 = z32*COS(beta) - x32*SIN(beta)
!       write(8,331) x23,y23,z23,x33,y33,z33
99011    FORMAT (/,'after rotation by beta about y-axis ',/,            &
     &           2(10x,3F9.5/))
!
!
         gamma2 = gamma
         IF ( z23<0.0 ) gamma2 = -gamma
!
         x34 = x33*COS(gamma2) - y33*SIN(gamma2)
         y34 = x33*SIN(gamma2) + y33*COS(gamma2)
         z34 = z33
!       write(8,332) x34,y34,z34
99012    FORMAT (/,'after rotation by gamma  about z-axis ',/,          &
     &           2(10x,3F9.5/))
!
!       now rotate back by -beta about y and then by -alpha about x 
!
!        x232 = -z23*SIN(beta) + x23*COS(beta)                      ! deleted, 11/20/08, HLA
!        y232 = y23                                                 ! deleted, 11/20/08, HLA
!        z232 = z23*COS(beta) + x23*SIN(beta)                       ! deleted, 11/20/08, HLA
!        x332 = -z33*SIN(beta) + x33*COS(beta)                      ! deleted, 11/20/08, HLA
!        y332 = y33                                                 ! deleted, 11/20/08, HLA
!        z332 = z33*COS(beta) + x33*SIN(beta)                       ! deleted, 11/20/08, HLA
         x342 = -z34*SIN(beta) + x34*COS(beta)
         y342 = y34
         z342 = z34*COS(beta) + x34*SIN(beta)
!
!        write(8,333)x232,y232,z232,x332,y332,z332,x342,y342,z342
99013    FORMAT (/,'after back rotation by beta about y-axis ',/,       &
     &           3(10x,3F9.5/))
!        x231 = x232                                                ! deleted, 11/18/08, HLA
!        y231 = y232*COS(alpha) + z232*SIN(alpha)                   ! deleted, 11/18/08, HLA
!        z231 = -y232*SIN(alpha) + z232*COS(alpha)                  ! deleted, 11/18/08, HLA
!        x331 = x332                                                ! deleted, 11/18/08, HLA
!        y331 = y332*COS(alpha) + z332*SIN(alpha)                   ! deleted, 11/18/08, HLA
!        z331 = -y332*SIN(alpha) + z332*COS(alpha)                  ! deleted, 11/18/08, HLA
         x341 = x342
         y341 = y342*COS(alpha) + z342*SIN(alpha)
         z341 = -y342*SIN(alpha) + z342*COS(alpha)
!
!        x230 = x231 + x1                                           ! deleted, 11/16/08, HLA
!        y230 = y231 + y1                                           ! deleted, 11/16/08, HLA
!        z230 = z231 + z1                                           ! deleted, 11/16/08, HLA
!
!        x330 = x331 + x1                                           ! deleted, 11/16/08, HLA
!        y330 = y331 + y1                                           ! deleted, 11/16/08, HLA
!        z330 = z331 + z1                                           ! deleted, 11/16/08, HLA
!
         x340 = x341 + x1
         y340 = y341 + y1
         z340 = z341 + z1
!
!      if ((ncy-1) == 1 .and. rc ==1) &
!    &   write(8,1) atom(linkatm(k)),x340,y340,z340  !   Final coordinates after rotation or bending  4-28-09 DU
         Xo(linkatm(k)) = x340
         Yo(linkatm(k)) = y340
         Zo(linkatm(k)) = z340
1        format('Final coordinates after rotation or bending :',3F9.5)
!
      ENDDO
!
      nobond = 1 
!     if ((ncy-1) == 1 .and. rc ==1) nobond = 0             ! 4-28-09 DU temp
      IF ( nobond /= 1 ) THEN
!        WRITE (63,99014)
99014    FORMAT (/,5x,'Bondlengths before and after ',/)
         k_dist = 0
         DO ik = 1 , N - 1
            DO jk = ik + 1 , N
               dx1 = x1o(ik) - x1o(jk)
               dy1 = y1o(ik) - y1o(jk)
               dz1 = z1o(ik) - z1o(jk)
               dist_in = SQRT(dx1*dx1+dy1*dy1+dz1*dz1)
               dx = Xo(ik) - Xo(jk)
               dy = Yo(ik) - Yo(jk)
               dz = Zo(ik) - Zo(jk)
               dist_out = SQRT(dx*dx+dy*dy+dz*dz)
!        calculate interbond angles  A-B-C                ! 4-28-09 DU
               if (ik == k1(1) .and. jk == k1(3)) then    ! 4-28-09 DU distance A-B
               dx1ab = dx1
               dy1ab = dy1
               dz1ab = dz1
               disab_in = dist_in
               dxab  = dx
               dyab  = dy
               dzab  = dz
               disab_out  = dist_out
               IDIS_ab = 1
               else if (ik == k1(3) .and. jk == k1(4)) then ! 4-28-09 distance B-C
               dx1bc = dx1
               dy1bc = dy1
               dz1bc = dz1
               disbc_in = dist_in
               dxbc  = dx
               dybc  = dy
               dzbc  = dz
               disbc_out  = dist_out
                if (IDIS_ab == 1 .and. i_nv > 12) then
                 CA_in  = -(dx1ab*dx1bc + dy1ab*dy1bc + dz1ab*dz1bc)/(disbc_in*disab_in)
!                write(63,"('CA_in',9F9.3)") dx1ab,dx1bc,dy1ab,dy1bc,dz1ab,dz1bc,disbc_in,disab_in,CA_in
                 org_angle = ATAN(SQRT(1-CA_in**2)/CA_in) * 57.29578
                 if ( org_angle < 0.0 ) org_angle = 180 +  org_angle
                 write(63, "('angle A-B-C before rot. or bend :',3I4, 3F9.3)")  k1(1),K1(3),K1(4), &
     &                disab_in, disbc_in, org_angle
                 CA_out = -(dxab*dxbc + dyab*dybc + dzab*dzbc)/(disbc_out*disab_out)
                 new_angle = ATAN(SQRT(1-CA_out**2)/CA_out) * 57.29578
                 if ( new_angle < 0.0 ) new_angle = 180 +  new_angle
                 write(63, "('angle A-B-C after  rot. or bend :',3I4, 3F9.3)")  k1(1),K1(3),K1(4), & 
     &                disab_out, disbc_out, new_angle
                end if
               end if
!
               IF ( dist_out<=1.7 ) THEN
                  IF ( ABS(dist_in-dist_out)>0.00001 ) THEN
                     k_dist = 1
!                    WRITE (63,99015) ATOm(ik) , ATOm(jk) , dist_in ,    &
!    &                               dist_out , k_dist
99015                FORMAT (/,5x,a6,2x,a6,2x,2F10.6,i5,/)
                  ENDIF
!                 WRITE (63,99016) ATOm(ik) , ATOm(jk) , dist_in ,       &
!    &                            dist_out
!
99016             FORMAT (5x,a6,2x,a6,3x,2F10.5)
               ENDIF
            ENDDO
         ENDDO
         if (IDIS_ab == 1 .and. i_nv > 12) WRITE(63,"(100('*'))")
      ENDIF
!
!
      IF ( I_Nv>=(nv-nvr+1) .OR. I_Nv==1 ) THEN                     ! 8-28-08  DU
!        write (61,"('nv & nvr & i_nv = ', 3I4)") nv, nvr, i_nv     ! 8-19-08  DU  temp
         IF ( NE_type==1 ) THEN                                     ! 9-30-08 DU
            CALL ENERGY_PM3(N,Xo,Yo,Zo,Egymol)                      ! this gets PM3 energy
         ELSE IF ( NE_type==2 ) THEN
            CALL ENERGY_B3LYP(N,Xo,Yo,Zo,Egymol)  ! this gets b3lyp/631g** energy 9-30-08 DU
         ELSE
            IF ( NE_type==3 .or. ne_type==4) CALL ENERGY_MOPAC(N,xo,yo,zo,egymol)  ! 4-28-09 DU 
         ENDIF
      ENDIF                                                         ! 8-28-08  DU
!
!         Convert to fractional coordinates
!
!       write(8,123)
99017 FORMAT (/,'input and output coordinates',/)
!      do  j = 1,n
!      zout(j) = zo(j)/r2
!      yout(j) = (yo(j)-zout(j)*r1)/s3
!      xout(j) = xo(j)-yout(j)*c3-zout(j)*c2
!      write(8,5112)atom(j),x(j),y(j),z(j),xout(j),yout(j),zout(j)
!      enddo
!
!      CALL ENERGY_pm3(n)
!       write(8,501)egymol
!
!         stop
      END SUBROUTINE BEND_BOND_NEW
!
!end--------------------------bend_bond_new----------------------------------end
!
      SUBROUTINE CROSS_PRODUCT(Kdx,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4, &
     &                         X5,Y5,Z5)
!
      USE PMIN_MODULE , ONLY:MATM , ATOm
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Kdx
      REAL*8 :: X1 , X2 , X3 , X4 , X5 , Y1 , Y2 , Y3 , Y4 , Y5 , Z1 ,  &
     &          Z2 , Z3 , Z4 , Z5
      INTENT (IN) Kdx , X1 , X2 , X3 , X4 , Y1 , Y2 , Y3 , Y4 , Z1 ,    &
     &            Z2 , Z3 , Z4
      INTENT (OUT) X5 , Y5 , Z5
!
! Local variables
!
      REAL*8 :: cccx , cccy , cccz , ccx , ccy , ccz , cx , cy , cz ,   &
     &          x31 , x32 , x34 , y31 , y32 , y34 , z31 , z32 , z34
      INTEGER :: index1
!
      index1 = Kdx
!
!              Vector CA
!
      x31 = X1 - X3
      y31 = Y1 - Y3
      z31 = Z1 - Z3
!
!              Vector CB
!
      x32 = X2 - X3
      y32 = Y2 - Y3
      z32 = Z2 - Z3
!
!              Vector CD
!
      x34 = X4 - X3
      y34 = Y4 - Y3
      z34 = Z4 - Z3
!
!              Vector CA x CB
!
      cx = (y31*z32-y32*z31)
      cy = (z31*x32-z32*x31)
      cz = (x31*y32-x32*y31)
      IF ( index1==2 ) THEN
         X5 = cx + X3
         Y5 = cy + Y3
         Z5 = cz + Z3
         GOTO 99999
      ENDIF
!
!
!             Vector CE = CD x (CA x CB)
!
      ccx = (y34*cz-cy*z34)
      ccy = (z34*cx-cz*x34)
      ccz = (x34*cy-cx*y34)
!
      IF ( index1==1 ) THEN
         X5 = ccx + X3
         Y5 = ccy + Y3
         Z5 = ccz + Z3
         GOTO 99999
      ENDIF
!
!             Vector   CD x CE
!
      cccx = (-y34*ccx+ccy*z34)
      cccy = (-z34*ccx+ccz*x34)
      cccz = (-x34*ccy+ccx*y34)
!
      IF ( index1==3 ) THEN
         X5 = cccx + X3
         Y5 = cccy + Y3
         Z5 = cccz + Z3
      ENDIF
!
99999 END SUBROUTINE CROSS_PRODUCT
!
!end-----------------------cross_product---------------------------end
!
      SUBROUTINE ATOM_LINKS(N,Cell,X,Y,Z,K11,K2,Linkatom)
!
      USE PMIN_MODULE , ONLY:MATM , ATOm
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: K2 , N
      REAL*8 , DIMENSION(6) :: Cell
      INTEGER , DIMENSION(4) :: K11
      INTEGER , DIMENSION(MATM) :: Linkatom
      REAL*8 , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN) Cell , K11 , N , X , Y , Z
      INTENT (OUT) Linkatom
      INTENT (INOUT) K2
!
! Local variables
!
      REAL*8 :: a , b , c , c1 , c2 , c3 , dist , dx , dy , dz
      INTEGER :: i , ik , j , k , kfind , nfind
      INTEGER , DIMENSION(MATM,MATM) :: idist
      INTEGER , DIMENSION(MATM) :: idx , idx1 , idx2 , list
      INTEGER , DIMENSION(4) :: k1
!
      DO i = 1 , 4
         k1(i) = K11(i)
      ENDDO
!
! Now read the atom coordinates from unit 10
!
 
!        read(7,*)a,b,c,c1,c2,c3
      a = Cell(1)
      b = Cell(2)
      c = Cell(3)
      c1 = Cell(4)
      c2 = Cell(5)
      c3 = Cell(5)
!        write(8,120)a,b,c,c1,c2,c3
!       write(8,120)cell
99001 FORMAT (/,5x,'From ATOM_LINKS cell ',6F9.5,/)
!        read(7,*)index,(k1(i),i=1,4)
!        write(8,121) index,(k1(i),i=1,4)
99002 FORMAT (5x,'index,k1(i) ',6I5,/)
!        read(10,*)natom
!        n= 0
!10      n = n + 1
!        read(10,1,end=100) atom(n),x(n),y(n),z(n)
!        write(8,1)atom(n),x(n),y(n),z(n)
99003 FORMAT (2x,a6,3F12.6,5A5)
99004 FORMAT ('ATOM',1x,a2,a4,3F10.6,i5,2F10.6)
!           go to 10
!100     n = n-1
!
!         Calculate bondlengths
!
      DO i = 1 , N
         DO j = 1 , N
            idist(i,j) = 0
         ENDDO
      ENDDO
!
      K2 = 0
      DO i = 1 , N
!              linkatom(i) = 0
         list(i) = 0
      ENDDO
!
      DO i = 1 , N
         DO j = 1 , N
            dx = a*(X(i)-X(j))
            dy = b*(Y(i)-Y(j))
            dz = c*(Z(i)-Z(j))
            dist = SQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c3+2*dy*dz*c1+        &
     &             2*dx*dz*c2)
            IF ( ATOm(i)(1:1)/='H' .OR. ATOm(j)(1:1)/='H' ) THEN          ! 11-24-08 DU
               IF ( dist<1.6 .AND. dist>0.6 ) THEN
                  idist(i,j) = 1
!                 write(8,4)atom(i),atom(j),dist
99005             FORMAT (5x,a6,2x,a6,3x,f10.4)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!          do 40 i = 1,n
!          write(8,41) (idist(i,j),j=1,n)
99006 FORMAT (40I3)
!40        continue
!
      idx1(1) = k1(4)
      kfind = 1
      DO
         nfind = 0
         DO k = 1 , kfind
            idx(k) = idx1(k)
            DO i = 1 , N
               DO j = 1 , N
!
                  IF ( i/=k1(3) .AND. j/=k1(3) ) THEN                     ! Exclude the parent link
                     IF ( idist(i,j)/=0 ) THEN
                        IF ( idist(i,j)/=99 ) THEN
                           IF ( idist(idx(k),j)==1 ) THEN
                              idist(idx(k),j) = 99
!
                              nfind = nfind + 1
                              idx2(nfind) = j
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         IF ( nfind>0 ) THEN
            kfind = nfind
            DO ik = 1 , nfind
               idx1(ik) = idx2(ik)
            ENDDO
            CYCLE
         ENDIF
!
         DO i = 1 , N
            DO j = 1 , N
               IF ( idist(i,j)==99 ) THEN
!            write(8,21)i,j,idist(i,j)
!             k2 = k2 + 1
                  list(j) = 1
99007             FORMAT (/,5x,'i,j,idist(i,j) ',3I5)
               ENDIF
            ENDDO
         ENDDO
         list(k1(4)) = 1
         DO i = 1 , N
            IF ( list(i)==1 ) THEN
               K2 = K2 + 1
               Linkatom(K2) = i
            ENDIF
         ENDDO
!          write(8,23) atom(k1(4)),k1(4),(linkatom(i),i=1,k2)
99008    FORMAT (/,5x,'list of atom numbers linked to terminal atom ',  &
     &           a6,2x,'no.',i5,/,5x,40I3)
         EXIT
      ENDDO
!           Print *, k2                 
!           STOP
      END SUBROUTINE ATOM_LINKS
!
!end-----------------------atom_links----------------------------------end
!
      SUBROUTINE ENERGY_PM3(N,Xo,Yo,Zo,Egymol)
!
!-----Program for calculation of Intramolecular energy, ENERGY_pm3
!
      USE PMIN_MODULE , ONLY:MATM , ATOm                                                  
      USE RSSPM3_MODULE
      USE RSS5_MODULE , ONLY:I_Nv , E_Weight
      USE BEND1_MODULE , ONLY:EGYmol0
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 :: Egymol
      INTEGER :: N
      REAL*8 , DIMENSION(MATM) :: Xo , Yo , Zo
      INTENT (IN) N , Xo , Yo , Zo
      INTENT (INOUT) Egymol
!
! Local variables
!
      REAL*8 :: angle
      CHARACTER(68) :: e_line
      INTEGER :: i , k
      INTEGER :: INDEX, NCONT
!
!      n = NMOD
       NCONT = 0                                                         ! 1-22-09
 
!---- PM3 energy calcd with G03, make necessary input file on unit # 60   
!
      OPEN (UNIT=60,FILE='fort.60',STATUS='unknown')
      angle = 0.0
      WRITE (60,99001) 
99001 FORMAT ('RUNGAUSS'/'#p pM3'/'         '/                          &
     &        'PM3 intermolecular energy',                              &
     &       /'       '/'   0   1')
      DO i = 1 , N
         WRITE (60,99002) ATOm(i) , Xo(i) , Yo(i) , Zo(i)
99002    FORMAT (a1,4x,3F15.9)
      ENDDO
      WRITE (60,99003)
99003 FORMAT ('      ')
      CLOSE (UNIT=60)
      Egymol = 0.0
!
!   Call system ('cat fort.60 >> fort.88')
!
      CALL SYSTEM('g03 < fort.60 >& fort.26')   
101   CLOSE (UNIT=60,STATUS='delete')
      OPEN (UNIT=26,STATUS='old')
      DO
         READ (26,'(A)',END=100) e_line                                  ! find intramolecular PM3 energy on G03 log file
         k = INDEX(e_line,' Energy=')
         IF ( k/=0 ) THEN
            READ (e_line(9:26),99004) Egymol                             ! PM3 E in Hartrees
99004       FORMAT (F18.12)
            Egymol = Egymol*627.5095                                     ! convert to kcal
            Egymol = Egymol*E_Weight                                     ! 8-13-08
!      write(8,315)egymol, egymol-egymol0, ncy-1, rc, i_nv
99005       FORMAT (/,5x,'intramolecular energy from pm3 ',2F15.6,3I6)
            CLOSE (UNIT=26,STATUS='delete')
            RETURN
         ENDIF
      ENDDO
 100  NCONT = NCONT + 1                                                  ! 1-22-09 DU
      WRITE (*,99006) NCY , RC, NCONT                                    ! 1-22-09 DU 
99006 FORMAT ('**cannot find PM3 energy on fort.26,ncy,rc,ncont ; RE_RUN g03** '&
     &        ,3I5)
      if (NCONT <= 2) then
      CALL SYSTEM('rm -f fort.26 ')                                      ! 1-22-09 DU
      CALL SYSTEM('g03 < fort.60 >& fort.26')                            ! 1-22-09 DU
      GO TO 101                                                          ! re-run g03 1-22-09 DU 
      else
      return
      end if                                                             ! 1-22-09 DU 
!     STOP
!
!     RETURN
!
      END SUBROUTINE ENERGY_PM3
!
!end---------------------energy_pm3----------------------------------end
!
      SUBROUTINE ENERGY_MOPAC(N,Xo,Yo,Zo,Egymol)           ! 4-28-09
!
!-----Program for calculation of Intramolecular energy, ENERGY_pm3
!
      USE PMIN_MODULE , ONLY:MATM , ATOm                                                  
      USE RSSPM3_MODULE
      USE RSS5_MODULE , ONLY:I_Nv , E_Weight, NE_type
      USE BEND1_MODULE , ONLY:EGYmol0
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 :: Egymol
      INTEGER :: N
      REAL*8 , DIMENSION(MATM) :: Xo , Yo , Zo
      INTENT (IN) N , Xo , Yo , Zo
      INTENT (INOUT) Egymol
!
! Local variables
!
      REAL*8 :: angle
      CHARACTER(68) :: e_line
      INTEGER :: i , k
      INTEGER :: INDEX, NCONT
!
!      n = NMOD
       NCONT = 0                                                         ! 1-22-09
 
!---- PM3 energy calcd with G03, make necessary input file on unit # 60   
!
      OPEN (UNIT=60,FILE='conf_E.mop',STATUS='unknown')
      angle = 0.0
      if ( NE_type == 3) then
      WRITE (60,99001) 
      else if ( NE_type == 4) then
      WRITE (60,99089)
      end if
99089 FORMAT ('PM6'/'PM6 intermolecular energy'/                        &
     &        'All coordinates are Cartesian')
99001 FORMAT ('PM3'/'PM3 intermolecular energy'/                        &
     &        'All coordinates are Cartesian')     
      DO i = 1 , N
         WRITE (60,99002) ATOm(i) , Xo(i) , Yo(i) , Zo(i)
99002    FORMAT (a1,4x,3(F15.9,1X,'0'))
      ENDDO
      WRITE (60,99003)
99003 FORMAT ('      ')
      CLOSE (UNIT=60)
      Egymol = 0.0
!
!
      CALL SYSTEM('/export/software/MOPAC2009.exe conf_E')              ! 11-4-2010
!     CALL SYSTEM('/opt/mopac/MOPAC2009.exe conf_E')   
!     CALL SYSTEM('./mopac.com')
!     CALL SYSTEM('wait')
 101  OPEN (UNIT=26,FILE='conf_E.arc',STATUS='unknown')
!101  OPEN (UNIT=26,STATUS='OLD')
      DO
         READ (26,'(A)',END=100) e_line                                  ! find intramolecular PM3 energy on MOPAC *.arc file
         k = INDEX(e_line,'          HEAT OF FORMATION')
         IF ( k/=0 ) THEN
            READ (e_line(36:52),99004) Egymol                            ! PM3 E Kcal
99004       FORMAT (F17.5)
            Egymol = Egymol*E_Weight                                     ! 8-13-08
!           write(8,315)egymol, egymol-egymol0, ncy-1, rc, i_nv
99005       FORMAT (/,5x,'intramolecular energy from pm3 ',2F15.6,3I6)
            CLOSE (UNIT=26,STATUS='delete')
!           CLOSE (UNIT=60,STATUS='delete')
            CALL SYSTEM('rm -f *.out')
            RETURN
         ENDIF
      ENDDO
 100  NCONT = NCONT + 1                                                  ! 1-22-09 DU
      WRITE (*,99006) NCY , RC, NCONT                                    ! 1-22-09 DU 
99006 FORMAT ('**cannot find PM3 energy on fort.26,ncy,rc,ncont ; RE_RUN MOPAC '&
     &        ,3I5)
      if (NCONT <= 2) then
      CALL SYSTEM('rm -f *.arc *.out ')                                  ! 4-28-09 DU
      CALL SYSTEM('/opt/mopac/MOPAC2009.exe conf_E')                     ! 4-28-09 DU
      GO TO 101                                                          ! re-run MOPAC 4-28-09 DU 
      else
      return
      end if                                                             ! 1-22-09 DU 
!     STOP
!
      END SUBROUTINE ENERGY_MOPAC
!
!end---------------------energy_mopac----------------------------------end
!
      SUBROUTINE ENERGY_B3LYP(N,Xo,Yo,Zo,Egymol)
!
!----Program for evaluation of intramolecular b3lyp/631g* energy
!
      USE PMIN_MODULE , ONLY:MATM , ATOm                                                  
      USE RSSPM3_MODULE
      USE RSS5_MODULE , ONLY:I_Nv , E_Weight
      USE BEND1_MODULE , ONLY:EGYmol0
      USE F77KINDS                        
      IMPLICIT NONE
!
      REAL*8 , PARAMETER :: E_CONVERT = 627.5095
!
! Dummy arguments
!
      REAL*8 :: Egymol
      INTEGER :: N
      REAL*8 , DIMENSION(MATM) :: Xo , Yo , Zo
      INTENT (IN) N , Xo , Yo , Zo
      INTENT (INOUT) Egymol
!
! Local variables
!
      REAL*8 :: angle
      REAL*8 :: egymol_81
      CHARACTER(68) :: e_line
      INTEGER :: i , k
      INTEGER :: INDEX
!
!     CHARACTER(37) :: E_PM3
!
!      n = NMOD
 
!---- b3lyp/631g* energy calcd with G03, make necessary input file on unit # 60
!
101   OPEN (UNIT=60,FILE='fort.60',STATUS='unknown')
      angle = 0.0
      WRITE (60,99001) 
99001 FORMAT ('RUNGAUSS'/'#p b3lyp/6-31G*'/'         '/                 &
     &        ' intramolecular energy w/b3lyp-6-31G*  ',                &
     &       /'       '/'   0   1')
      DO i = 1 , N
         WRITE (60,99002) ATOm(i) , Xo(i) , Yo(i) , Zo(i)
99002    FORMAT (a1,4x,3F15.9)
      ENDDO
      WRITE (60,99003)
99003 FORMAT ('      ')
      CLOSE (UNIT=60)
      Egymol = 0.0
!
!   Call system ('cat fort.60 >> fort.88')
!   for linux chem-69 or chem-65 or 91a-454
!
102   CALL SYSTEM('g03 < fort.60 >& fort.26')   
      CLOSE (UNIT=60,STATUS='delete')
      OPEN (UNIT=26,STATUS='old')
      DO
         READ (26,'(A)',END=100) e_line                                  ! find intramolecular b3lyp/631g* energy on G03 log file
         k = INDEX(e_line,' SCF Done:')
         IF ( k/=0 ) THEN
            READ (e_line(27:42),99004) egymol_81                         ! g03 E in Hartrees
99004       FORMAT (F16.9)
            Egymol = egymol_81*E_CONVERT
            Egymol = Egymol*E_Weight                                     ! 8-13-08
            WRITE (8,99005) Egymol , Egymol - EGYmol0 , NCY - 1 , RC ,  &
     &                      I_Nv
99005       FORMAT (/,5x,'intramolecular energy from pm3 ',2F20.8,3I6)
            CLOSE (UNIT=26,STATUS='delete')
            RETURN
         ENDIF
      ENDDO
 100  WRITE (*,99006) NCY , RC , INV , IJJ
99006 FORMAT (                                                          &
     &'**Quit, cannot find b3lyp/631g* energy on fort.26,ncy,rc,inv,ijj'&
     &,4I5)
      CALL SYSTEM('rm -f fort.26')                                       ! 1-22-09 DU
      GO TO 102                                                          ! re-run g03 1-22-09 DU
!     STOP
!
!     RETURN
!
      END SUBROUTINE ENERGY_B3LYP
!
!end---------------------energy_b3lyp----------------------------------end
!
      SUBROUTINE ORTHO_COD(N,Cell,X,Y,Z,Xo,Yo,Zo)
!
!----x,y,z = fractional coordinates
!    xo,yo,zo = orthoigoinal coordinates in Angs
!    cell = a,b,c, cosines of angles
!
      USE PMIN_MODULE , ONLY:MATM , ATOm , PI
      USE F77KINDS                        
      IMPLICIT NONE
!*--ORTHO_COD6940
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: N
      REAL*8 , DIMENSION(6) :: Cell
      REAL*8 , DIMENSION(MATM) :: X , Xo , Y , Yo , Z , Zo
      INTENT (IN) Cell , N , X , Y , Z
      INTENT (OUT) Xo , Yo , Zo
!
! Local variables
!
      REAL*8 :: c1 , c2 , c3 , s1 , s3
      REAL*8 :: DSQRT
      INTEGER :: i
      REAL*8 , DIMENSION(MATM) :: xa , ya , za
!
!*** End of declarations rewritten by SPAG
!
!
!        c1 = cell(4)
!        c2 = cell(5)
!        c3 = cell(6)
!        s1 = sqrt(1-c1*c1)
!        s2 = sqrt(1-c2*c2)
!        s3 = sqrt(1-c3*c3)
!        r1 = (c1 -  c2*c3)/s3
!        r2 = sqrt(s1*s1 -((c2-c1*c3)/s3)**2)
!     pi = 3.141592653589793D0                                    ! deleted, 11/21/08, HLA
      c1 = Cell(4)                                                ! 9-23-08
      c2 = Cell(5)                                                ! 9-23-08
      c3 = Cell(6)                                                ! 9-23-08
      IF ( c1>10.0 ) THEN                                         ! 9-23-08
         c1 = COS(c1*PI/180.0)                                    ! 9-23-08
         c2 = COS(c2*PI/180.0)                                    ! 9-23-08
         c3 = COS(c3*PI/180.0)                                    ! 9-23-08
      ENDIF
!
      s1 = SQRT(1-c1*c1)                                          ! 9-23-08
!     s2 = SQRT(1-c2*c2)                                          ! deleted, 11/18/08, HLA 
      s3 = SQRT(1-c3*c3)                                          ! 9-23-08
!     r1 = (c1-c2*c3)/s3                                          ! deleted, 11/18/08, HLA
!     r2 = SQRT(s1*s1-((c2-c1*c3)/s3)**2)                         ! deleted, 11/18/08, HLA 
!
      DO i = 1 , N
         xa(i) = Cell(1)*X(i)
         ya(i) = Cell(2)*Y(i)
         za(i) = Cell(3)*Z(i)
      ENDDO
!        get  orthogonal coordinates
      DO i = 1 , N
         Xo(i) = xa(i) + ya(i)*c3 + za(i)*c2
         Yo(i) = ya(i)*s3 + za(i)*(c1-c2*c3)/s3
         Zo(i) = za(i)*DSQRT(s1*s1-((c2-c1*c3)/s3)**2)
      ENDDO
      END SUBROUTINE ORTHO_COD
!
!end---------------------ortho_cod------------------------------end
!
      SUBROUTINE FRACTIONAL_COD(N,Cell,Xo,Yo,Zo,X,Y,Z)
!
!----x,y,z = fractional coordinates
!    xo,yo,zo = orthogonal coordinates in Angs
!    cell = a,b,c, cosines of angles
!
      USE PMIN_MODULE , ONLY:MATM , ATOm
      USE F77KINDS                        
      IMPLICIT NONE
!*--FRACTIONAL_COD7014
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: N
      REAL*8 , DIMENSION(6) :: Cell
      REAL*8 , DIMENSION(MATM) :: X , Xo , Y , Yo , Z , Zo
      INTENT (IN) Cell , N , Xo , Yo , Zo
      INTENT (INOUT) X , Y , Z
!
! Local variables
!
      REAL*8 :: c1 , c2 , c3 , r1 , r2 , s1 , s3
      INTEGER :: i , j
!
!*** End of declarations rewritten by SPAG
!
!
      c1 = Cell(4)
      c2 = Cell(5)
      c3 = Cell(6)
      c3 = Cell(6)
      s1 = SQRT(1-c1*c1)
!     s2 = SQRT(1-c2*c2)                                          ! deleted, 11/18/08, HLA
      s3 = SQRT(1-c3*c3)
      r1 = (c1-c2*c3)/s3
      r2 = SQRT(s1*s1-((c2-c1*c3)/s3)**2)
!
!        get  fractional coordinates
      DO j = 1 , N
         Z(j) = Zo(j)/r2
         Y(j) = (Yo(j)-Z(j)*r1)/s3
         X(j) = Xo(j) - Y(j)*c3 - Z(j)*c2
      ENDDO
      DO i = 1 , N
         X(i) = X(i)/Cell(1)
         Y(i) = Y(i)/Cell(2)
         Z(i) = Z(i)/Cell(3)
      ENDDO
      END SUBROUTINE FRACTIONAL_COD
!
!end--------------------fractional_cod---------------------------end
!
      SUBROUTINE PC_LIMITS(ncy,Nmol,Irot,Nvr,Nv,Pcitm,Pcmin,Pcmax)
!
      USE PMIN_MODULE , ONLY:NMOLD , NVD , PI                               
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: Irot , Nmol , Nv , Nvr, ncy
      REAL*8 , DIMENSION(NVD) :: Pcitm , Pcmax , Pcmin
      INTENT (IN) Irot , Nmol , Nv , Nvr , Pcitm
      INTENT (INOUT) Pcmax , Pcmin
!
! Local variables
!
      REAL*8 :: anglemax , anglemin
      INTEGER :: i , imol , j
!
      DO i = 1 , 3                                    ! unit cell lengths
         IF ( Pcitm(i)<10.0 ) THEN
            Pcmin(i) = Pcitm(i) - 0.1*Pcitm(i)        ! changed from 0.25 to 0.1 7/10/08 HLA
            Pcmax(i) = Pcitm(i) + 0.1*Pcitm(i)        ! changed from 0.25 to 0.1 7/10/08 HLA
         ELSE
            Pcmin(i) = Pcitm(i) - 1.0
            Pcmax(i) = Pcitm(i) + 1.0
         ENDIF
      ENDDO
      anglemin = -0.866025
      anglemax = +0.866025
      DO i = 4 , 6                                    ! unit cell cosines
         Pcmin(i) = Pcitm(i) - 0.1                    ! changed from 0.2 to 0.1, 7/4/08 HLA
         Pcmax(i) = Pcitm(i) + 0.1                                  
         IF ( Pcmin(i)<=anglemin ) Pcmin(i) = anglemin
         IF ( Pcmax(i)>=anglemax ) Pcmax(i) = anglemax 
      ENDDO
!
      DO imol = 1 , Nmol                              ! model rotation
         DO i = 6*imol + 1 , 6*imol + 3
            Pcmin(i) = Pcitm(i) - 0.1                 ! changed from 0.2 to 0.2, 7/4/08 HLA
            Pcmax(i) = Pcitm(i) + 0.1                             
         ENDDO
         DO i = 6*imol + 4 , 6*imol + 6               ! model translation
            Pcmin(i) = Pcitm(i) - 0.25                ! reduced translation from 0.5 to 0.25, 7/4/08, HLA
            Pcmax(i) = Pcitm(i) + 0.25             
         ENDDO
      ENDDO
      IF ( Irot>=1 ) THEN
         DO j = Nv - Nvr + 1 , Nv
            Pcmin(j) = -PI/2.0                        ! min & max rotation angles in radian for rotbond
            Pcmax(j) = PI/2.0
         ENDDO
      ENDIF
!      write(61,333)
!333       format ('Parameters, values, min & max ranges'/ &
!                 '  #       value     minimum     maximum')
!      DO i=1,nv                                         
!         write(61,334) i, pcitm(i), pcmin(i), pcmax(i)
!334          format (i3,3f12.5)
!      END DO
      END SUBROUTINE PC_LIMITS
!
!end---------------------subroutine pc_limits----------------------------end
!
! The subroutine USER_E_1 below calculates the total lattice energy for atoms within one asymmetric unit
! If the asymmetric unit contains one molecule, the calculated lattice energy is per molecule.
! The total lattice energy has three parts given by 6-exp formula. All intermolecular interactions within
! a radius of dmax (6 to 10 Angstrom) are considered around the atoms of the asymmetric unit. The atoms
! surrounding the molecule are obtained by the space group symmetry.
! The Coulomb and van der Waals' terms are calculated using the reciprocal space summation if IEC = 0 and
! IEV = 0. If IEC = 1 and IEV = 1, they are calculated by summation in direct space.
!
      SUBROUTINE USER_E_1(Cell,X,Y,Z,Ec,Ev,Er,E_nmol,Pe)
!
      USE PMIN_MODULE , ONLY:NTX_min , NTY_min , NTZ_min , NTX_max ,    &
     &    NTY_max , NTZ_max , NMOl , NATm , IDX , ISYstem , N , ILS ,   &
     &    NSYm , ICYcle , NENtry , IMOde , NMOLD , NVD , DEL_param ,    &
     &    MATM , ATOm , STHl , SYM , Q , A1 , B1 , C1 , CK , WT_mol ,   &
     &    DMAx , DDMax , GNI , CE12 , VE12 , K1 , N11 , IHKl , IRIj ,   &
     &    KK , I_Cross , A_Cross , B_Cross , C_Cross , ANAme1 , ANAme2 ,&
     &    A12 , B12 , C12 , PI , TWOPI , ANG_TO_RADIAN, dd_cross        ! 5-4-10 DU                      
      USE RSSPM3_MODULE
      USE F77KINDS                        
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL*8 :: Ec , Er , Ev , E_nmol , Pe
      REAL*8 , DIMENSION(6) :: Cell
      REAL*8 , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN) Cell , X , Y , Z
      INTENT (INOUT) Ec , Er , Ev , E_nmol , Pe
!
! Local variables
!
      REAL*8 :: a , alpha , asqrt , asquare , b , beta , bsq , c , c11 ,&
     &          c22 , c33 , cb1 , cb2 , ce1 , ce11 , ce2 , ce22 , ce3 , &
     &          ck11 , ck22 , ck33 , ck44 , coulomb_en , d , dc ,       &
     &          distmax = -999.0 , distmin = 999.0 , dx , dy , dz , e1 ,&
     &          e11 , e12 , e123 , e123_all , e1_all , e2 , e22 ,       &
     &          e2_all , e3 , e33 , e3_all , ec_nmol , energy , eng1 ,  &
     &          eng2 , er_nmol , ev_nmol , gamma , gi , r1 , r2 , r3 ,  &
     &          r4 , r5 , r6 , repulse_en , rij , sah , sin1 , sin2 ,   &
     &          sin3 , sins , sqrtpi , sum123 , sumx , sumy , sumz ,    &
     &          sum_en , v11 , v12 , vander_en , ve , ve1 , ve11 , ve2 ,&
     &          ve22 , ve3 , ve4 , vol
      REAL*8 , DIMENSION(6) :: ah , rh
      REAL*8 :: DERFC
      REAL*8 :: DSQRT
      REAL*8 , DIMENSION(1000) :: fc , fcb , fv , fvb , qh
      INTEGER :: i , i1 , iend1 , iend2 , ig , ih , ihmax , ihmin , ik ,&
     &           ik3 , ik4 , il , imol1 , imol2 , index1 , istart1 ,    &
     &           istart2 , j , jh , jk , jl , k , k11 , k111 , k2 ,     &
     &           k22 , k3 , k33 , kk1 , kmax , kmin , ktx , kty , ktz , &
     &           lmax , lmin , n12 , natom , ne1 , ne1_all , ntotal ,   &
     &           ntx , ntx2 , nty , nty2 , ntz , ntz2
      INTEGER :: INT
      REAL*8 :: x11 , x2 , x22 , x4 , xc , xni , y11 , y2 , y22 , y4 ,  &
     &          yc , z11 , z2 , z22 , z4 , zc , zz
!
      DO i = 1 , 1000                                                                 
         qh(i) = 0.0
      ENDDO
!
      zz = NSYm
      natom = N
!
      DO i = 1 , 6
         ah(i) = Cell(i)
      ENDDO
!
      IF ( ah(4)>10.0 ) THEN
         ah(4) = COS(ah(4)*ANG_TO_RADIAN)                       ! 5/13/08, HLA, altered
         ah(5) = COS(ah(5)*ANG_TO_RADIAN)                       ! 5/13/08, HLA, altered
         ah(6) = COS(ah(6)*ANG_TO_RADIAN)                       ! 5/13/08, HLA, altered
99001    FORMAT (7X,A4,A2,3F8.3,3F9.5)
      ENDIF
      c11 = ah(4)
      c22 = ah(5)
      c33 = ah(6)
      a = ah(1)
      b = ah(2)
      c = ah(3)
      alpha = ACOSD(ah(4))                                      ! 5/14/08, HLA, changed to acosd
      beta = ACOSD(ah(5))                                                  
      gamma = ACOSD(ah(6))                                     
!         write(8,520)(ah(i),i=1,3),alpha,beta,gamma
!
      sah = 0.5*(alpha+beta+gamma)
      sins = SIN(ANG_TO_RADIAN*sah)                             ! 5/14/08, HLA, altered
      sin1 = SIN(ANG_TO_RADIAN*(sah-alpha))                     ! 5/14/08, HLA, altered
      sin2 = SIN(ANG_TO_RADIAN*(sah-beta))                      ! 5/14/08, HLA, altered
      sin3 = SIN(ANG_TO_RADIAN*(sah-gamma))                     ! 5/14/08, HLA, altered
      vol = 2.0*ah(1)*ah(2)*ah(3)*DSQRT(sins*sin1*sin2*sin3)
!        write(8,5126) vol
99002 FORMAT (/,' vol = ',f10.4)
!
!     density = 1.6605*zz*WT_mol/vol                            ! deleted, 11/20/08, HLA
!  101 FORMAT (' CRYSTAL AND RECIP LATTICE PARAMETERS')
!      write(8,517)AH,vol,density
99003 FORMAT (/,5x,' Cell, Volume & density = '/8F10.4,/)
!
      CALL RECIP(ah,rh)
!       WRITE (8,107) RH
99004 FORMAT (/,'reciprocal cell ',8F10.5,/)
!
      IF ( IMOde/=2 ) THEN
!
         r1 = rh(1)**2
         r2 = rh(2)**2
         r3 = rh(3)**2
         r4 = 2.0*rh(1)*rh(2)*rh(6)
         r5 = 2.0*rh(1)*rh(3)*rh(5)
         r6 = 2.0*rh(3)*rh(2)*rh(4)
!
!       print *, ' nentry ',nentry
         IF ( NENtry<=1 ) THEN   ! Excludes calculation of hkl's & gni except on 1st entry
!
            ihmax = INT(2.0*STHl/rh(1))
            kmax = INT(2.0*STHl/rh(2))
            lmax = INT(2.0*STHl/rh(3))
!       write(8,108)sthl,rh(1),rh(2),rh(3),ihmax,kmax,kmin
99005       FORMAT (5x,'sthl,rh1,2,3,ihmax,kmax,kmin ',4F8.4,3I5,/)
            ihmin = -ihmax
            kmin = -kmax
            lmin = 0
            gi = 2.0
            IF ( ABS(rh(4))<=0.001 .AND. ABS(rh(6))<=0.001 ) THEN
               kmin = 0
               gi = 2.0*gi
            ENDIF
            IF ( ABS(rh(5))<=0.000001 ) THEN                    ! 2-15-09 DU
               ihmin = 0
               gi = 2.0*gi
            ENDIF
!      WRITE (8,103) IHMAX, KMAX, LMAX,IHMIN,KMIN,LMIN, STHL, gi
99006       FORMAT (/,' H,K,L MAX & MIN VALUES ',6I5,                   &
     &              ' FOR SINTH/L 'F6.3,' gi = ',f8.2/)
!
!  MAKE A set of reflections or hkl values for given value of sthl
!
!     sumfc = 0.0                                               ! deleted, 11/20/08, HLA
            K1 = 0
            DO ih = ihmin , ihmax
               DO ik = kmin , kmax
                  DO il = lmin , lmax
                     xni = 1.0
                     IF ( ih/=0 .OR. ik/=0 .OR. il/=0 ) THEN
                        IF ( ih==0 .AND. ihmin==0 ) xni = 0.5*xni
                        IF ( ik==0 .AND. kmin==0 ) xni = 0.5*xni
                        IF ( il==0 ) xni = 0.5*xni
                        d = SQRT((r1*ih**2)+(r2*ik**2)+(r3*il**2)       &
     &                      +(r5*ih*il)+(r4*ih*ik)+(r6*ik*il))
!       write(8,531)ih,ik,il,D
99007                   FORMAT (3I5,f9.3)
!        if(D <= STHL) qh1(k1) = D
                        IF ( d<=STHl ) THEN
                           K1 = K1 + 1
                           IHKl(K1,1) = ih
                           IHKl(K1,2) = ik
                           IHKl(K1,3) = il
                           qh(K1) = d
                           GNI(K1) = gi*xni
!        write(8,541)k1,ih,ik,il,d,qh(k1)   ! comment out later
99008                      FORMAT (4I6,2F10.4)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
!
!
!  PICK UP EACH REFLECTION ONE BY ONE
         DO k = 1 , K1
            jh = IHKl(k,1)
            jk = IHKl(k,2)
            jl = IHKl(k,3)
            IF ( NENtry>1 ) THEN
               d = r1*jh**2 + r2*jk**2 + r3*jl**2 + r5*jh*jl +          &
     &             r4*jh*jk + r6*jk*jl
               qh(k) = SQRT(d)
            ENDIF
!
            fc(k) = 0.0
            fcb(k) = 0.0
            fv(k) = 0.0
            fvb(k) = 0.0
!
            DO ig = 1 , NSYm
!
!    LOOP OVER ALL ATOMS
!
               DO j = 1 , N
!
!      X2 = X(J)*SYM(IG,1) + SYM(IG,2)
!      Y2 = Y(J)*SYM(IG,3) + SYM(IG,4)
!      Z2 = Z(J)*SYM(IG,5) + SYM(IG,6)
                  x2 = SYM(ig,1) + X(j)*SYM(ig,2) + Y(j)*SYM(ig,3)      &
     &                 + Z(j)*SYM(ig,4)
                  y2 = SYM(ig,5) + X(j)*SYM(ig,6) + Y(j)*SYM(ig,7)      &
     &                 + Z(j)*SYM(ig,8)
                  z2 = SYM(ig,9) + X(j)*SYM(ig,10) + Y(j)*SYM(ig,11)    &
     &                 + Z(j)*SYM(ig,12)
                  fc(k) = fc(k) + Q(j)*COS(TWOPI*(jh*x2+jk*y2+jl*z2))
                  fv(k) = fv(k) + A1(j)*COS(TWOPI*(jh*x2+jk*y2+jl*z2))
!         IF(ICENT == 1) THEN
                  fcb(k) = fcb(k) + Q(j)*SIN(TWOPI*(jh*x2+jk*y2+jl*z2))
                  fvb(k) = fvb(k) + A1(j)*SIN(TWOPI*(jh*x2+jk*y2+jl*z2))
!         ENDIF
               ENDDO
            ENDDO
!
!         IF(ICENT == 0) THEN
!            FC(K) = 2*FC(K)
!         ELSE
!      if(icent == 1) THEN
            fc(k) = SQRT(fc(k)**2+fcb(k)**2)
            fv(k) = SQRT(fv(k)**2+fvb(k)**2)
!      ENDIF
!      SUMFC = SUMFC + ABS(FC(K))
!        if(abs(FC(k)) > 0.001) then
!      WRITE(8,345)K,jh(k),jk(k),jl(k),qh(k),FC(K),FV(K)
99009       FORMAT (5X,4I5,3F10.4)
!        endif
         ENDDO
!
!---- LOOP TO CALCULATE ST FACTOR ENDS
!
!     Now calculate the reciprocal sum for Coulomb and Van der Waals Energy
!
!      write(8,*) 'i,hkl,qh(i),bsq,fc(i),gni(i),clmb,coulomb'
         ce2 = 0.0
         ve2 = 0.0
         sqrtpi = DSQRT(PI)
         zz = NSYm
         ck11 = 1.0/(2.0*PI*vol*zz)
         ck22 = PI**4.5/(3*vol*zz)
         ck33 = (PI**3.0)*(CK**6)/12.0
         ck44 = (PI**3.0)*(CK**3)*zz/(6.0*vol)
!      write(8,711) zz,vol,ck11,ck22,ck33,ck44
99010    FORMAT (/,5x,'zz,vol,ck11,ck22,ck33,ck44 = ',2F10.4,4F15.10,/)
!
         DO i = 1 , K1
            bsq = PI*qh(i)*qh(i)/(CK*CK)
            cb1 = SQRT(bsq)
                        ! Check if taking sqare root for getting constant B is valid or not
            cb2 = sqrtpi*DERFC(cb1) + (1.0/(2*bsq*cb1)-1.0/cb1)         &
     &            *EXP(-bsq)
            ce22 = 332.17*fc(i)*fc(i)*GNI(i)*(EXP(-bsq))/(qh(i)*qh(i))
            ce22 = ck11*ce22
            ce2 = ce2 + ce22
            ve22 = GNI(i)*fv(i)*fv(i)*(qh(i)**3)*cb2
            ve22 = -ck22*ve22
            ve2 = ve2 + ve22
         ENDDO
!      No special position assumed presently
         ce3 = 0.0
         ve3 = 0.0
         ve4 = 0.0
         DO i = 1 , N
            ce3 = ce3 + Q(i)*Q(i)
            ve3 = ve3 + A1(i)*A1(i) ! No special position of atom assumed n(i) = 1.0
            ve4 = ve4 + A1(i)
         ENDDO
         ce3 = -ce3*CK*332.17
         ve3 = ck33*ve3
         ve4 = -ck44*ve4*ve4
!
! Now calculate the 1st part (Direct space sum) of Coulomb energy
! and van der Waals attraction and repulsion terms
!
         eng2 = 0.0
         e1 = 0.0
         e2 = 0.0
         e3 = 0.0
         e123 = 0.0
!     sum_e = 0.0                                                        ! deleted, 11/20/08, HLA
         n12 = 0
         ve11 = 0.0
         ce11 = 0.0
         ntotal = 0
         sum123 = 0.0
!&&&       if(nentry > 1) go to 3333
         IF ( NENtry>1 ) THEN
!        if(nentry > 0) go to 333                                        ! always go, by pass nentry > 1
!
!
!      For nentry > 1
!
            CE12 = 0.0
            VE12 = 0.0
            ce11 = 0.0
            ve11 = 0.0
            e3 = 0.0
!
!  Use n11 interactions only
!
            DO ik = 1 , N11
               i = IRIj(ik,1)
               i1 = IRIj(ik,2)
               j = IRIj(ik,3)
               ktx = IRIj(ik,4)
               kty = IRIj(ik,5)
               ktz = IRIj(ik,6)
               x11 = X(i)*a
               y11 = Y(i)*b
               z11 = Z(i)*c
               x4 = a*(X(i1)*SYM(j,2)+Y(i1)*SYM(j,3)+Z(i1)*SYM(j,4)     &
     &              +SYM(j,1)+ktx)
               dx = x11 - x4
               y4 = b*(X(i1)*SYM(j,6)+Y(i1)*SYM(j,7)+Z(i1)*SYM(j,8)     &
     &              +SYM(j,5)+kty)
               dy = y11 - y4
               z4 = c*(X(i1)*SYM(j,10)+Y(i1)*SYM(j,11)+Z(i1)*SYM(j,12)  &
     &              +SYM(j,9)+ktz)
               dz = z11 - z4
!
               rij = DSQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+2*dy*dz*c11+   &
     &               2*dx*dz*c22)
!        if(rij == 0.0) go to 220
!
!      Now calculate the potential energy
!
               asquare = PI*CK*CK*rij*rij
               asqrt = DSQRT(asquare)
!
               IF ( j==1 .AND. ktx==0 .AND. kty==0 .AND. ktz==0 ) THEN
!       Print *,'ils,icycle,nentry,ik,n11 intra',ils,icycle,nentry,ik,n11
                  e12 = 332.17*Q(i)*Q(i1)*(DERFC(asqrt)-1.0)/rij
                  v12 = -A1(i)*A1(i1)                                   &
     &                  *((1.0+asquare+0.5*asquare*asquare)*EXP         &
     &                  (-asquare)-1.0)/(rij**6)
                  IF ( i==i1 ) THEN
                     e12 = 0.5*e12
                     v12 = 0.5*v12
                  ENDIF
                  CE12 = CE12 + e12
                  VE12 = VE12 + v12
!        n12 = n12 + 1
!
                  IF ( NENtry<=5 ) THEN
                     WRITE (28,99011) NENtry , n12 , N11 , i , i1 ,     &
     &                                rij , e12 , v12 , CE12 , VE12
!        write(28,768)nentry,n12,n11,i,i1,rij,e12,v12,ce12,ve12
99011                FORMAT ('nentry,n12,n11,i1,rij,e12,v12,ce12,ve12', &
     &                       5I5,f9.5,4F16.10)
                  ENDIF
                  CYCLE
               ENDIF
!
!          if(rij > dmax)go to 2225
!
               e11 = 332.17*Q(i)*Q(i1)*DERFC(asqrt)/rij
! -----------------------5-1-08 using Cross-term PC --------------------------------------
               IF ( I_Cross/=0 ) THEN                                    ! 5-1-08 DU
                  DO KK = 1 , I_Cross                                    ! 5-1-08 DU
                     IF ( (ATOm(i)==ANAme1(KK) .AND. ATOm(i1) &
     &                    ==ANAme2(KK).AND. Rij <= dd_cross(KK)) .OR. &
     &                    (ATOm(i)==ANAme2(KK) .AND. ATOm(i1) &
     &                    ==ANAme1(KK) .AND. Rij <= dd_cross(KK)).OR. &
     &                    (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1) &
     &                    ==ANAme2(KK).AND. Rij <= dd_cross(KK)) .OR. &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1) &
     &                    ==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN    ! 5-4-10 DU
                        A12 = A_Cross(KK)                                ! 5-1-08 DU
                        B12 = B_Cross(KK)                                ! 5-1-08 DU
                        C12 = C_Cross(KK)                                ! 5-1-08 DU
                        GO TO 99922                                      ! found cross term, exit DO Loop I_croos 
                     ELSE                                                ! 5-1-08 DU
                        A12 = A1(i)*A1(i1)                               ! 5-1-08 DU
                        B12 = B1(i)*B1(i1)                               ! 5-1-08 DU
                        C12 = C1(i) + C1(i1)                             ! 5-1-08 DU
                     ENDIF                                               ! 5-1-08 DU
                  ENDDO                                                  ! 5-1-08 DU
               ELSE                                                      ! 5-1-08 DU
                  A12 = A1(i)*A1(i1)                                     ! 5-1-08 DU
                  B12 = B1(i)*B1(i1)                                     ! 5-1-08 DU
                  C12 = C1(i) + C1(i1)                                   ! 5-1-08 DU
               ENDIF                                                     ! 5-1-08 DU
99922          v11 = A12*((1.0+asquare+0.5*asquare*asquare)             &
     &               *EXP(-asquare))/(rij**6)                            ! 5-1-08 DU
               e33 = B12*EXP(-C12*rij)                                   ! 5-1-08 DU
!        v11 = a1(i)*a1(i1)*((1.0 + asquare + 0.5*asquare*asquare)*exp(-asquare))/(rij**6) ! 5-1-08 DU
!        e22 = -A1(i)*A1(i1)/(rij**6)
!        e33 = B1(i)*B1(i1)*exp(-(C1(i)+C1(i1))*rij)                     ! 5-1-08 DU
!------------------------------------------------------------------------! 5-1-08 DU
!##################### test for cross-term 5-1-08 #######################! 5-1-08 DU
               IF ( rij<3.0 .AND. IMOde==0 .AND. ICYcle==1 .AND.        &
     &              NCY<=2 ) THEN                                        ! 5-1-08 DU
!         write(20, 2202) nentry,atom(i),atom(i1),rij,A1(i),A1(i1),B1(i),B1(i1),B12,C1(i),C1(i1),e11,e22,e33     ! 5-1-08 DU
                  WRITE (20,99012) NENtry , ATOm(i) , ATOm(i1) , rij ,  &
     &                             A1(i) , A1(i1) , B1(i) , B1(i1) ,    &
     &                             A12 , B12 , e11 , e22 , e33           ! 5-1-08 DU
99012             FORMAT (I6,1X,A6,'---',A6,2X,f6.4,4F9.4,1X,2F16.8,    &
     &                    3F12.6,' #2')                                  ! 5-1-08 DU
               ENDIF                                                     ! 5-1-08 DU
!########################################################################! 5-1-08 DU
               IF ( i==i1 ) THEN
                  e11 = 0.5*e11
                  v11 = 0.5*v11
                  e22 = 0.5*e22
                  e33 = 0.5*e33
               ENDIF
               ce11 = ce11 + e11
               ve11 = ve11 - v11
!         e2 = e22 + e2
               e3 = e33 + e3
!         if(nentry == 113 .or. nentry == 114) then
!        write(28,2209)icycle,nentry,ik,atom(i),atom(i1),rij,e33,e3
!          if(ik == 77)then
!        write(28,6567)icycle,nentry,atom(i),atom(i1),x(i),y(i),z(i),x(i1),y(i1),z(i1),rij
!6567    format(/,'icycle,nentry,atoms,xyz(i),xyz(i1),rij',2i5,2x,a6,2x,a6,6f10.6,f10.6,/)
!          endif
 
!         endif
!        sum123 = sum123 + e11 + e22 + e33
!        if( nentry == 353 .or. nentry == 354 .or. nentry == 362) then
!         if(nentry <= 4 .or. ils == 1)then
!        write(8,2201)nentry,ik,nrij,i,i1,j,rij,e11,CE11
!         Print 2201,ils,nentry,ik,nrij,i,i1,j,rij,e33,e3
99013          FORMAT (5x,'nentry,ik,nrij,i,i1,j,rij,e33,e3 ',7I5,      &
     &                 3F14.8)
!         endif
            ENDDO
!
!         if(icycle > 5) then
!      write(28,2208)icycle,nentry,n12,n11, CE12,VE12,CE11,VE11,e3
99014       FORMAT (/,5x,'icycle,nentry,n12,n11,CE12,VE12,CE11,VE11,e3 '&
     &              ,4I6,5F16.10,/)
         ELSE
!
!      for  translations along x, y, z
!
            N11 = 0
            VE12 = 0.0
            CE12 = 0.0
            ntx = INT(DMAx/a) + 2
            nty = INT(DMAx/b) + 2
            ntz = INT(DMAx/c) + 2
            NTX_min = 0
            NTX_max = 0
            NTY_min = 0
            NTY_max = 0
            NTZ_min = 0
            NTZ_max = 0
            ntx2 = 2*ntx + 1
            nty2 = 2*nty + 1
            ntz2 = 2*ntz + 1
!       write(8,67)ntx,nty,ntz
99015       FORMAT (/,5x,                                               &
     &           'translations along x,y,z both +ve and -ve directions '&
     &           ,3I5,/)
!
! Loop over all the atoms
!
            DO i = 1 , natom
               x11 = X(i)*a
               y11 = Y(i)*b
               z11 = Z(i)*c
! Loop over all the second atoms
               DO i1 = i , natom
! Loop over all the symmetries
                  DO j = 1 , NSYm
! Loop over all the translations
                     DO k111 = 1 , ntx2
                        x4 = a*(X(i1)*SYM(j,2)+Y(i1)*SYM(j,3)+Z(i1)     &
     &                       *SYM(j,4)+SYM(j,1)+(ntx+1-k111))
                        dx = x11 - x4
                        DO k2 = 1 , nty2
                           y4 = b*(X(i1)*SYM(j,6)+Y(i1)*SYM(j,7)+Z(i1)  &
     &                          *SYM(j,8)+SYM(j,5)+(nty+1-k2))
                           dy = y11 - y4
                           DO k3 = 1 , ntz2
                              z4 = c*(X(i1)*SYM(j,10)+Y(i1)*SYM(j,11)   &
     &                             +Z(i1)*SYM(j,12)+SYM(j,9)+(ntz+1-k3))
                              dz = z11 - z4
!
                              rij = DSQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+&
     &                              2*dy*dz*c11+2*dx*dz*c22)
                              ntotal = ntotal + 1
                              IF ( rij>distmax ) distmax = rij
!
                              IF ( rij>0.000001 ) THEN
!
                                 k11 = 5 + (ntx+1-k111)
                                 k22 = 5 + (nty+1-k2)
                                 k33 = 5 + (ntz+1-k3)
 
                                 asquare = PI*CK*CK*rij*rij
                                 asqrt = DSQRT(asquare)
!
                                 IF ( j==1 .AND. k11==5 .AND.           &
     &                                k22==5 .AND. k33==5 ) THEN
!            if(rij < 2.0)then                                           ! removed 7/11/08 HLA
!          write(8,166)i,i1,atom(i),atom(i1),rij                         ! removed 7/11/08 HLA
!            endif                                                       ! removed 7/11/08 HLA
!166       format(2i5,2x,a6,2x,a6,2x,f10.5)                              ! removed 7/11/08 HLA
                                    e12 = 332.17*Q(i)*Q(i1)             &
     &                                 *(DERFC(asqrt)-1.0)/rij
                                    v12 = -A1(i)*A1(i1)                 &
     &                                 *((1.0+asquare+0.5*asquare*      &
     &                                 asquare)*EXP(-asquare)-1.0)      &
     &                                 /(rij**6)
                                    IF ( i==i1 ) THEN
                                       e12 = 0.5*e12
                                       v12 = 0.5*v12
                                    ENDIF
                                    CE12 = CE12 + e12
                                    VE12 = VE12 + v12
                                    n12 = n12 + 1
!
                                    N11 = N11 + 1
                                    IRIj(N11,1) = i
                                    IRIj(N11,2) = i1
                                    IRIj(N11,3) = j
                                    IRIj(N11,4) = ntx + 1 - k111
                                    IRIj(N11,5) = nty + 1 - k2
                                    IRIj(N11,6) = ntz + 1 - k3
!         if(n11 == 77)then
!       write(28,6567)icycle,nentry,atom(i),atom(i1),x(i),y(i),z(i),x(i1),y(i1),z(i1),rij
99016                               FORMAT (/,                          &
     &                          'icycle,nentry,atoms,xyz(i),xyz(i1),rij'&
     &                          ,2I5,2x,a6,2x,a6,6F10.6,f10.6,/)
!         endif
!
                                    CYCLE
                                 ENDIF
!
                                 IF ( rij>0.000001 .AND. rij<=DMAx )    &
     &                                THEN
                                    IF ( NENtry==1 ) THEN
                                       ktx = ntx + 1 - k111
                                       kty = nty + 1 - k2
                                       ktz = ntz + 1 - k3
                                       N11 = N11 + 1
                                       IRIj(N11,1) = i
                                       IRIj(N11,2) = i1
                                       IRIj(N11,3) = j
                                       IRIj(N11,4) = ntx + 1 - k111
                                       IRIj(N11,5) = nty + 1 - k2
                                       IRIj(N11,6) = ntz + 1 - k3
                                       IF ( NTX_min>ktx ) NTX_min = ktx
                                       IF ( NTX_max<ktx ) NTX_max = ktx
                                       IF ( NTY_min>kty ) NTY_min = kty
                                       IF ( NTY_max<kty ) NTY_max = kty
                                       IF ( NTZ_min>ktz ) NTZ_min = ktz
                                       IF ( NTZ_max<ktz ) NTZ_max = ktz
                                    ENDIF
!      if(n11 < 10) then
!        write(8,727)n11,(irij(n11,ij),ij=1,6)
!727     format(5x,'n11,irij ',8i5)
!      endif
!
! Pick up the minimum contact distance distmin
!                          k1min = k11                                   ! deleted, 11/20/08, HLA
!                          k2min = k22                                   ! deleted, 11/20/08, HLA
!                          k3min = k33                                   ! deleted, 11/20/08, HLA
!                          jmin = j                                      ! deleted, 11/20/08, HLA
!                          atom1m = ATOm(i)                              ! deleted, 11/20/08, HLA
!                          atom2m = ATOm(i1)                             ! deleted, 11/20/08, HLA
                                    IF ( rij<10.0 .AND. rij<distmin )   &
     &                                 distmin = rij
!
!      Now calculate the potential energy
!
                                    e11 = 332.17*Q(i)*Q(i1)*DERFC(asqrt)&
     &                                 /rij
! -----------------------5-1-08 using Cross-term PC --------------------------------------
                                    IF ( I_Cross/=0 ) THEN               ! 5-1-08 DU
                                       DO KK = 1 , I_Cross               ! 5-1-08 DU
                                         IF ( (ATOm(i)==ANAme1(KK)      &
     &                                      .AND. ATOm(i1)              &
     &                                      ==ANAme2(KK).AND. Rij <= dd_cross(KK) ) .OR.          &
     &                                      (ATOm(i)==ANAme2(KK)        &
     &                                      .AND. ATOm(i1)              &
     &                                      ==ANAme1(KK).AND.Rij <= dd_cross(KK))   .OR.          &
     &                    (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1) &
     &                    ==ANAme2(KK).AND. Rij <= dd_cross(KK)) .OR. &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1) &
     &                    ==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN    ! 5-4-10 DU
                                         A12 = A_Cross(KK)               ! 5-1-08 DU
                                         B12 = B_Cross(KK)               ! 5-1-08 DU
                                         C12 = C_Cross(KK)               ! 5-1-08 DU
                                         GO TO 5684                      ! found cross term, exit DO Loop I_croos 
                                         ELSE                            ! 5-1-08 DU
                                         A12 = A1(i)*A1(i1)              ! 5-1-08 DU
                                         B12 = B1(i)*B1(i1)              ! 5-1-08 DU
                                         C12 = C1(i) + C1(i1)            ! 5-1-08 DU
                                         ENDIF                           ! 5-1-08 DU
                                       ENDDO                             ! 5-1-08 DU
                                    ELSE                                 ! 5-1-08 DU
                                       A12 = A1(i)*A1(i1)                ! 5-1-08 DU
                                       B12 = B1(i)*B1(i1)                ! 5-1-08 DU
                                       C12 = C1(i) + C1(i1)              ! 5-1-08 DU
                                    ENDIF                                ! 5-1-08 DU
5684                                v11 = A12*                          &
     &                                 ((1.0+asquare+0.5*asquare*asquare&
     &                                 )*EXP(-asquare))/(rij**6)         ! 5-1-08 DU
                                    e22 = -A12/(rij**6)                  ! 5-1-08 DU
                                    e33 = B12*EXP(-C12*rij)              ! 5-1-08 DU
!        v11 = a1(i)*a1(i1)*((1.0 + asquare + 0.5*asquare*asquare)*exp(-asquare))/(rij**6) ! 5-1-08 DU
!        e22 = -A1(i)*A1(i1)/(rij**6)                                    ! 5-1-08 DU
!        e33 = B1(i)*B1(i1)*exp(-(C1(i)+C1(i1))*rij)                     ! 5-1-08 DU
!##################### test for cross-term 5-1-08 #######################! 5-1-08 DU
                                    IF ( rij<5.0 .AND. IMOde==0 .AND.   &
     &                                 ICYcle==1 .AND. NCY<=2 ) THEN     ! 5-1-08
                                       WRITE (20,99017) NENtry ,        &
     &                                    ATOm(i) , ATOm(i1) ,&
     &                                    rij , A1(i) , A1(i1) , B1(i) ,&
     &                                    B1(i1) , A12 , B12 , C12 ,    &
     &                                    e11 , e22 , e33                 ! 5-1-08
!         write(20, 2202) nentry,atom(i),atom(i1),rij,A1(i),A1(i1),B1(i),B1(i1),A12,B12,e11,e22,e33  ! 5-1-08
!2202     format(I6,1X,A6,'---',A6,2X,f6.4,4f9.4,1X,2F16.8,3f12.6,' #1')  ! 5-1-08
99017                                  FORMAT (I6,1X,2A6,1X,2X,f6.4,    &
     &                                    4F9.4,1X,2F16.8,f5.2,3F12.6,  &
     &                                    ' #1')                          ! 5-1-08
                                    ENDIF                                 ! 5-1-08
!#########################################################################! 5-1-08
! ---------------------------------------------------------------------------------------
                                    IF ( i==i1 ) THEN
                                       e11 = 0.5*e11
                                       v11 = 0.5*v11
                                       e22 = 0.5*e22
                                       e33 = 0.5*e33
                                    ENDIF
                                    ce11 = ce11 + e11
                                    ve11 = ve11 - v11
                                    e2 = e22 + e2
                                    e3 = e33 + e3
                                    sum123 = sum123 + e11 + e22 + e33
!       write(28,2209)icycle,nentry,n11,atom(i),atom(i1),rij,e33,e3
99018                               FORMAT (5x,                         &
     &                             'icycle,nentry,n11,atoms,rij,e33,e3 '&
     &                             ,3I5,2x,a6,2x,a6,f10.5,2F14.8)
                                 ENDIF
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
!
!
!          if(nentry <= 5) then
!       write(28,2208)icycle,nentry,n12,n11, CE12,VE12,CE11,VE11,e3
!       write(28,*) 'ntx_min,ntx_max,nty_min,nty_max,ntz_min,ntz_max',ntx_min,ntx_max,nty_min,nty_max,ntz_min,ntz_max
!          endif
!
            ENDDO
         ENDIF
!         endif
!
!        Now print the result
!
         ne1 = N11 + n12
!        if(nentry == 1 .and. icycle == 1) then
!         if(nentry == 353 .or. nentry == 354 .or. nentry == 362)then
!        print 817, k1,n12, n11,ne1,ntotal
!        write(8,817) k1,n12, n11,ne1
!        write(11,817) k1,n12, n11,ne1
!        write(12,817) k1,n12, n11,ne1
!        write(13,817) icycle,nentry,k1,n12, n11,ne1,ntotal
!        write(8,817) icycle,nentry,k1,n12,n11,ne1,ntotal
!817     format(' From POT_E icycle,nentry,reflections &  interactions - intra, inter & total ', 7i6)
!        endif
         ce1 = ce11 + CE12
         coulomb_en = ce1 + ce2 + ce3
         ve1 = ve11 + VE12
!
         ve = ve1 + ve2 + ve3 + ve4
!     vander_en  = e2
         vander_en = ve
         repulse_en = e3
         sum_en = coulomb_en + vander_en + repulse_en
         Ec = coulomb_en
!      Ev = vander_en
         Ev = ve
!      Er = repulse_en
         Er = e3
!
         Pe = sum_en
!
!      write(8,817) icycle-1,nentry,k1,n12, n11,ne1,ntotal                      ! removed 7/11/08 HLA
!      write(8,417) CE11,CE12,CE1,CE2,CE3, coulomb_en                           ! removed 7/11/08 HLA
!      write(8,419)VE11,VE12,VE1,VE2,VE3,VE4,VE                                 ! removed 7/11/08 HLA
!      write(8,817) icycle-1,nentry,k1,n12,n11,ne1,ntotal                       ! removed 7/11/08 HLA
!817    format('  From LATT_E icycle,nentry,refls & interactions - intra, inter & total ', 7i6)  ! removed 7/11/08 HLA
!417    format(5x,'CE11,CE12,CE1, CE2, CE3, total = ',6f14.8)                   ! removed 7/11/08 HLA
!419    format(5x,'VE11,VE12,VE1,VE2,VE3,VE4,VE = ',7f14.8)  ! removed 7/11/08 HLA
!      write(8,412)dmax,Ec,Ev,Er,PE                                             ! removed 7/11/08 HLA
!412    format(/,5x,'dmax,Ec,EV,Er,PE from Ewald smmation for Ec and Ev ',5f12.6,/)              ! removed 7/11/08 HLA
!  
!---- Additional energy contribution from nmol molecules in the asymmetric unit  
         IF ( NMOl>1 ) THEN
            ec_nmol = 0.0
            ev_nmol = 0.0
            er_nmol = 0.0                                          
            E_nmol = 0.0                                
!
            istart1 = 1
            DO imol1 = 1 , NMOl - 1
               iend1 = istart1 + NATm(imol1) - 1
               DO ik3 = istart1 , iend1
                  istart2 = 0
                  DO ik = 1 , imol1
                     istart2 = istart2 + NATm(ik)
                  ENDDO
                  istart2 = istart2 + 1
                  DO imol2 = imol1 + 1 , NMOl
                     iend2 = istart2 + NATm(imol2) - 1
                     DO ik4 = istart2 , iend2
                        x11 = a*X(ik3)                   
                        y11 = b*Y(ik3)                   
                        z11 = c*Z(ik3)                   
                        x22 = a*X(ik4)                   
                        y22 = b*Y(ik4)                   
                        z22 = c*Z(ik4)                   
                        dx = x11 - x22                      
                        dy = y11 - y22                  
                        dz = z11 - z22                 
                        rij = SQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+       &
     &                        2*dy*dz*c11+2*dx*dz*c22)                                      
                        ec_nmol = ec_nmol + 332.17*Q(ik3)*Q(ik4)/rij
                        ev_nmol = ev_nmol - A1(ik3)*A1(ik4)/(rij**6)
                        er_nmol = er_nmol + B1(ik3)*B1(ik4)             &
     &                            *EXP(-(C1(ik3)+C1(ik4))*rij)                  
99019                   FORMAT (4x,a6,2x,a6,f9.4,4F12.6)
                     ENDDO
                     istart2 = istart2 + NATm(imol2)
                  ENDDO
               ENDDO
               istart1 = istart1 + NATm(imol1)
            ENDDO
         ENDIF
 
         E_nmol = ec_nmol + ev_nmol + er_nmol
         Ec = Ec + ec_nmol
         Ev = Ev + ev_nmol
         Er = Er + er_nmol
         Pe = Pe + E_nmol
      ENDIF
!
!     write(8,5367)dmax,n11                                              ! removed 7/11/08 HLA
!      Print 5367,dmax,n11
!5367  format(/,' Ec,Ev from Ewald summation, Er from direct summation' &! removed 7/11/08 HLA
!             ' upto dmax, No. interactions = ',f8.2,' A',i10)           ! removed 7/11/08 HLA
!     write(8,5368)nmol,Ec,Ev,Er,E_nmol,PE                               ! removed 7/11/08 HLA
!      Print 5368,nmol,Ec,Ev,Er,E_nmol,PE
5368    format(5x,'nmol,Ec,Ev,Er,E_nmol,E_total ',i5,2x,5f12.6,/)        ! removed 7/11/08 HLA
!        endif
!
!
      IF ( IMOde/=1 ) THEN
!
!-------------------------------------------------------------------------------
!    For direct summatio in real space
!-------------------------------------------------------------------------------
!
!      write(8,*) ' '                                                    ! removed 7/11/08 HLA
!       Print *,'             --------> From direct summation<---------- '
!       write(8,*)'           --------> From direct summation<---------- '       
!      write(8,*) ' '                                                    ! removed 7/11/08 HLA
         ntx = INT(DDMax/a) + 4
         nty = INT(DDMax/b) + 4
         ntz = INT(DDMax/c) + 4
         ntx2 = 2*ntx - 1
         nty2 = 2*nty - 1
         ntz2 = 2*ntz - 1
         sumx = 0.0
         sumy = 0.0
         sumz = 0.0
         DO i = 1 , N
            sumx = sumx + X(i)
            sumy = sumy + Y(i)
            sumz = sumz + Z(i)
         ENDDO
         xc = sumx/N
         yc = sumy/N
         zc = sumz/N
!    write(8,105)ddmax,ntx,nty,ntz,xc,yc,zc
99020    FORMAT (/'ddmax,ntx,nty,ntz,centriod xc,yc,zc ',f9.3,3I5,      &
     &           3F10.6,/)
!
! Loop over all the atoms
!
!     write(8,661)
99021    FORMAT (/,5x,                                                  &
     &'atom      e       e_sum       e1       e2        e3       e_sum  &
     & No. interactions',/)
!        n2 = 2*ntrans - 1
         e1_all = 0.0
         e2_all = 0.0
         e3_all = 0.0
         e123_all = 0.0
         ne1_all = 0
         energy = 0.0
!     eg2 = 0.0                                                          ! deleted, 11/20/08, HLA
         DO i = 1 , N
            x11 = X(i)*a
            y11 = Y(i)*b
            z11 = Z(i)*c
            eng2 = 0.0
            e1 = 0.0
            e2 = 0.0
            e3 = 0.0
            e123 = 0.0
            ne1 = 0
! Loop over all the second atoms
!        do 20 i1 = 1,n
! Loop over all the symmetries
            DO j = 1 , NSYm
!loop over all the translations
!       do 20 k1 = 1,19
!       do 20 k1 = 1,n2
               DO kk1 = 1 , ntx2
                  DO k2 = 1 , nty2
                     DO k3 = 1 , ntz2
!
!         if(i == 1)then
                        x4 = a*(xc*SYM(j,2)+yc*SYM(j,3)+zc*SYM(j,4)     &
     &                       +SYM(j,1)+(ntx-kk1))
                        dx = a*xc - x4
                        y4 = b*(xc*SYM(j,6)+yc*SYM(j,7)+zc*SYM(j,8)     &
     &                       +SYM(j,5)+(nty-k2))
                        dy = b*yc - y4
                        z4 = c*(xc*SYM(j,10)+yc*SYM(j,11)+zc*SYM(j,12)  &
     &                       +SYM(j,9)+(ntz-k3))
                        dz = c*zc - z4
                        dc = SQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+        &
     &                       2*dy*dz*c11+2*dx*dz*c22)
                        IF ( dc>0.000001 .AND. dc<=DDMax ) THEN
!         endif
!
!-----Loop over all the second atoms
                           DO i1 = 1 , N
                              x4 = a*(X(i1)*SYM(j,2)+Y(i1)*SYM(j,3)     &
     &                             +Z(i1)*SYM(j,4)+SYM(j,1)+(ntx-kk1))
                              dx = x11 - x4
                              y4 = b*(X(i1)*SYM(j,6)+Y(i1)*SYM(j,7)     &
     &                             +Z(i1)*SYM(j,8)+SYM(j,5)+(nty-k2))
                              dy = y11 - y4
                              z4 = c*(X(i1)*SYM(j,10)+Y(i1)*SYM(j,11)   &
     &                             +Z(i1)*SYM(j,12)+SYM(j,9)+(ntz-k3))
                              dz = z11 - z4
!
                              d = SQRT(dx*dx+dy*dy+dz*dz+2*dx*dy*c33+   &
     &                            2*dy*dz*c11+2*dx*dz*c22)
                              IF ( d>distmax ) distmax = d
!
                              IF ( d>0.000001 ) THEN
!
                                 k11 = 5 + (ntx-kk1)
                                 k22 = 5 + (nty-k2)
                                 k33 = 5 + (ntz-k3)
                                 index1 = j*k11*k22*k33
                                 IF ( index1/=125 ) THEN
                                       ! Exclude intramolecular distances
! Pick up the minimum contact distance distmin
!                          k1min = k11                                   ! deleted, 11/20/08, HLA
!                          k2min = k22                                   ! deleted, 11/20/08, HLA
!                          k3min = k33                                   ! deleted, 11/20/08, HLA
!                          jmin = j                                      ! deleted, 11/20/08, HLA
!                          atom1m = ATOm(i)                              ! deleted, 11/20/08, HLA
!                          atom2m = ATOm(i1)                             ! deleted, 11/20/08, HLA
                                    IF ( d<10.0 .AND. d<distmin )       &
     &                                 distmin = d
!
!      Now calculate the potential energy
 
!
                                    e11 = 0.5*332.17*Q(i)*Q(i1)/d
                                    e1 = e11 + e1                        ! 5/9/08, HLA
! -----------------------5-1-08 using Cross-term PC --------------------------------------
                                    IF ( I_Cross/=0 ) THEN               ! 5-1-08 DU
                                       DO KK = 1 , I_Cross               ! 5-1-08 DU
                                         IF ( (ATOm(i)==ANAme1(KK)      &
     &                                      .AND. ATOm(i1)              &
     &                                      ==ANAme2(KK).AND. d <= dd_cross(KK)) .OR.          &
     &                                      (ATOm(i)==ANAme2(KK)        &
     &                                      .AND. ATOm(i1)              &
     &                                      ==ANAme1(KK).AND. d <= dd_cross(KK)) .OR.          &  
     &                    (ATOm(i)(1:1)==ANAme1(KK)(1:1) .AND. ATOm(i1) &
     &                    ==ANAme2(KK).AND. Rij <= dd_cross(KK)) .OR. &
     &                   (ATOm(i)==ANAme2(KK) .AND. ATOm(i1)(1:1) &
     &                    ==ANAme1(KK)(1:1) .AND. Rij <= dd_cross(KK)) ) THEN    ! 5-4-10 DU
                                         A12 = A_Cross(KK)               ! 5-1-08 DU
                                         B12 = B_Cross(KK)               ! 5-1-08 DU
                                         C12 = C_Cross(KK)               ! 5-1-08 DU
                                         GO TO 99921                     ! found cross term, exit DO Loop I_croos 
                                         ELSE                            ! 5-1-08 DU
                                         A12 = A1(i)*A1(i1)              ! 5-1-08 DU
                                         B12 = B1(i)*B1(i1)              ! 5-1-08 DU
                                         C12 = C1(i) + C1(i1)            ! 5-1-08 DU
                                         ENDIF                           ! 5-1-08 DU
                                       ENDDO                             ! 5-1-08 DU
                                    ELSE                                 ! 5-1-08 DU
                                       A12 = A1(i)*A1(i1)                ! 5-1-08 DU
                                       B12 = B1(i)*B1(i1)                ! 5-1-08 DU
                                       C12 = C1(i) + C1(i1)              ! 5-1-08 DU
                                    ENDIF                                ! 5-1-08 DU
99921                               e22 = -0.5*A12/(d**6)                ! 10-13-09
                                    e2 = -0.5*A12/(d**6) + e2            ! 5-1-08 DU
                                    e33 = 0.5*B12*EXP(-C12*d)            ! 5-1-08 DU
                                    e3 = 0.5*B12*EXP(-C12*d) + e3        ! 5-1-08 DU
                                    ne1 = ne1 + 1
                                    eng1 = 0.5*(332.17*Q(i)*Q(i1)       &
     &                                 /d-A12/(d**6)+B12*EXP(-C12*d))    ! 5-1-08 DU
!        e22 = -0.5*A1(i)*A1(i1)/(d**6)                                  ! 5-1-08 DU
!        e2 = -0.5*A1(i)*A1(i1)/(d**6) + e2                              ! 5-1-08 DU
!        e33 = 0.5*B1(i)*B1(i1)*exp(-(C1(i)+C1(i1))*d)                   ! 5-1-08 DU
!        e3 =  0.5*B1(i)*B1(i1)*exp(-(C1(i)+C1(i1))*d) + e3              ! 5-1-08 DU
!        ne1 = ne1 + 1                                                   ! 5-1-08 DU
!        eng1 = 0.5*(332.17*q(i)*q(i1)/d -A1(i)*A1(i1)/(d**6) + B1(i)*B1(i1)*exp(-(C1(i)+C1(i1))*d)) ! 5-1-08 DU
! -----------------------------------------------------------------------! 5-1-08 DU
!##################### test for cross-term 5-1-08 #######################! 5-1-08 DU
                                    IF ( rij<3.0 .AND. IMOde==0 .AND.   &
     &                                 ICYcle==1 .AND. NCY<=2 ) THEN     ! 5-1-08 DU
                                       WRITE (20,99022) NENtry , ATOm(i)&
     &                                    , ATOm(i1) , rij , A1(i) ,    &
     &                                    A1(i1) , B1(i) , B1(i1) ,     &
     &                                    A12 , B12 , e11 , e22 , e33    ! 5-1-08 DU
99022                                  FORMAT (I6,1X,A6,'---',A6,2X,    &
     &                                    f6.4,4F9.4,1X,2F16.8,3F12.6,  &
     &                                    ' #3')                         ! 5-1-08 DU
!         write(20, 2202) nentry,atom(i),atom(i1),rij,A1(i),A1(i1),B1(i),B1(i1),B12,C1(i),C1(i1),e11,e22,e33 ! 5-1-08 DU
                                    ENDIF                                ! 5-1-08 DU
!########################################################################! 5-1-08 DU
                                    energy = energy + eng1
                                    eng2 = eng2 + eng1
                                 ENDIF
                              ENDIF
!
!
!           if(i == 1) then
!           write(8,1111)ne1,i,i1,atom(i),atom(i1),q(i),q(i1),d,e11,e1
!1111       format(5x,3i5,a4,2x,a4,5f10.4)
!           endif
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            e123 = e1 + e2 + e3
            e123_all = e123_all + e123
            e1_all = e1_all + e1
            e2_all = e2_all + e2
            e3_all = e3_all + e3
            ne1_all = ne1_all + ne1
!        write(8,811) atom(i),eng2,energy, e1, e2, e3, e123, ne1,ne1_all,e1_all,e2_all,e3_all,e123_all
!
99023       FORMAT (5x,a4,6F10.4,2I10,3x,4F10.4)
         ENDDO
!       write(8,815) atom1m,atom2m,jmin,k1min,k2min,k3min,distmin
99024    FORMAT (/,5x,'Minimum intermol contact : ',2A4,4I1,f9.3,/)
!        write(8,816) distmax
99025    FORMAT (/,5x,' Maximum calculated distance = ',f10.4,/)
!
         e1_all = e1_all + ec_nmol
         e2_all = e2_all + ev_nmol
         e3_all = e3_all + er_nmol
         e123_all = e123_all + E_nmol
!
!
      write (8,99026) ddmax, ne1_all
      Print 99026, ddmax, ne1_all                 ! altered, 7/13/09
99026    FORMAT (' Ec, Ev, Er, E_nmol, PE from direct summation to',f6.2,' Angs, # terms =', & ! altered, 7/13/09
                 i12)
      write(8,5368) nmol, e1_all, e2_all, e3_all, E_nmol, e123_all
      Print 5368, nmol,e1_all,e2_all,e3_all,E_nmol,e123_all
99027    FORMAT (5x,'nmol,imode,Ec,Ev,Er,E_nmol,PE ',2I5,5F12.6,/)
!          write(8,*)' imode,Ec ',imode,Ec
         IF ( IMOde==2 ) THEN
            Ec = e1_all
            Ev = e2_all
            Er = e3_all
            Pe = e123_all + E_nmol
!         Print 5365, nmol,imode,Ec,Ev,Er,E_nmol,PE
          write(8,99027) nmol,imode,Ec,Ev,Er,E_nmol,PE                   ! from S.P Feb 5,09
            GOTO 99999
         ENDIF
!
!        Ewald summation for Ec, but for others direct summation upto ddmax
!
!     e_total = Ec + e2_all + e3_all + E_nmol                            ! deleted, 11/20/08, HLA
!
         IF ( IMOde==3 ) THEN
            Ec = Ec
            Ev = e2_all
            Er = e3_all
            Pe = Ec + e2_all + e3_all + E_nmol
         ENDIF
!
!     write(8,5369)dmax,ddmax,n11,ne1_all                                ! removed 7/11/08 HLA
!      Print 5369,dmax,ddmax,n11,ne1_all
99028    FORMAT (/                                                      &
     &' Ec from Ewald summation + Ev, Er and E_nmol from direct summatio&
     &n',' up to ',2F8.2,' A',2I12)
!        write(8,5370) nmol,Ec,e2_all,e3_all,E_nmol,E_total              ! removed 7/11/08 HLA
!         Print 5370, nmol,Ec,e2_all,e3_all,E_nmol,PE
99029    FORMAT (5x,'nmol,Ec,Ev,Er,E_nmol,PE ',i5,5F12.6,/)
      ENDIF
!        all done
!
99999 END SUBROUTINE USER_E_1
