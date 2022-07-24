!*==s_recip.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_RECIP
      IMPLICIT NONE
!*--S_RECIP4
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE RECIP(X,Y)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , DIMENSION(6) :: X , Y
      INTENT (IN)X
      INTENT (OUT)Y
      END SUBROUTINE RECIP
      END INTERFACE
      END MODULE S_RECIP
!*==s_new_xyz.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_NEW_XYZ
      IMPLICIT NONE
!*--S_NEW_XYZ23
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE NEW_XYZ(CELL1,X,Y,Z,CELL_NEW,XNEW,YNEW,ZNEW)
      USE PMIN_MODULE , ONLY:NMOL , NATM , IDX , ISYSTEM , N , ILS ,    &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , PI              ! 9-23-08 from S.P
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , DIMENSION(6) :: CELL1
      REAL(R8KIND) , DIMENSION(NVD,6) :: CELL_NEW
      REAL(R8KIND) , DIMENSION(MATM) :: X , Y , Z
      REAL(R8KIND) , DIMENSION(NVD,MATM) :: XNEW , YNEW , ZNEW
      INTENT (OUT)XNEW , YNEW , ZNEW
      INTENT (INOUT)CELL_NEW
      END SUBROUTINE NEW_XYZ
      END INTERFACE
      END MODULE S_NEW_XYZ
!*==s_corr_xyz.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_CORR_XYZ
      IMPLICIT NONE
!*--S_CORR_XYZ50
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CORR_XYZ(CORR,FRAC_CHANGE,CELL1,X,Y,Z,CELL_NEW,X_COR,  &
     &                    Y_COR,Z_COR)
      USE PMIN_MODULE , ONLY:NMOL , NATM , IDX , ISYSTEM , N , ILS ,    &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , PI , ANG_TO_RADIAN          ! 9-23-08 from S.P
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , DIMENSION(6) :: CELL1 , CELL_NEW
      REAL(R8KIND) , DIMENSION(NVD) :: CORR
      REAL(R8KIND) , DIMENSION(4) :: FRAC_CHANGE
      REAL(R8KIND) , DIMENSION(MATM) :: X , X_COR , Y , Y_COR , Z ,     &
     &                                  Z_COR
      INTENT (IN)CORR , FRAC_CHANGE
      INTENT (INOUT)CELL_NEW , X_COR , Y_COR , Z_COR
      END SUBROUTINE CORR_XYZ
      END INTERFACE
      END MODULE S_CORR_XYZ
!*==s_corr_xyz_rot.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_CORR_XYZ_ROT
      IMPLICIT NONE
!*--S_CORR_XYZ_ROT79
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CORR_XYZ_ROT(CORR,FRAC_CHANGE,CELL1,X,Y,Z,CELL_NEW,    &
     &                        X_COR,Y_COR,Z_COR)
      USE PMIN_MODULE , ONLY:NMOL , NATM , IDX , ISYSTEM , N , ILS ,    &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , PI , ANG_TO_RADIAN                                               ! 8-19-08
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , DIMENSION(6) :: CELL1 , CELL_NEW
      REAL(R8KIND) , DIMENSION(NVD) :: CORR
      REAL(R8KIND) , DIMENSION(4) :: FRAC_CHANGE
      REAL(R8KIND) , DIMENSION(MATM) :: X , X_COR , Y , Y_COR , Z ,     &
     &                                  Z_COR
      INTENT (IN)CELL1 , CORR , FRAC_CHANGE , X , Y , Z
      INTENT (INOUT)CELL_NEW , X_COR , Y_COR , Z_COR
      END SUBROUTINE CORR_XYZ_ROT
      END INTERFACE
      END MODULE S_CORR_XYZ_ROT
!*==s_pot_e.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_POT_E
      IMPLICIT NONE
!*--S_POT_E109
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE POT_E(CELL,X,Y,Z,EC,EV,ER,E_NMOL,PE)
      USE PMIN_MODULE , ONLY:NTX_MIN , NTY_MIN , NTZ_MIN , NTX_MAX ,    &
     &    NTY_MAX , NTZ_MAX , NMOL , NATM , IDX , ISYSTEM , N , ILS ,   &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , STHL , SYM , Q , A1 , B1 , C1 , CK , WT_MOL ,   &
     &    DMAX , DDMAX , GNI , CE12 , VE12 , K1 , N11 , IHKL , IRIJ ,   &
     &    KK , I_CROSS , A_CROSS , B_CROSS , C_CROSS , ANAME1 , ANAME2 ,&
     &    ATOM_CODE , A12 , B12 , C12 , PI , TWOPI , ANG_TO_RADIAN                      ! 5-13-08
      USE RSS1_MODULE , ONLY:REFCODE
      USE RSSPM3_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) :: EC , ER , EV , E_NMOL , PE
      REAL(R8KIND) , DIMENSION(6) :: CELL
      REAL(R8KIND) , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN)CELL , X , Y , Z
      INTENT (INOUT)EC , ER , EV , E_NMOL , PE
      END SUBROUTINE POT_E
      END INTERFACE
      END MODULE S_POT_E
!*==s_mcpy.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_MCPY
      IMPLICIT NONE
!*--S_MCPY139
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE MCPY(A,R,N,M,MS)                                       !  MCPY   5
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: M , MS , N
      REAL(R8KIND) , DIMENSION(1) :: A , R
      INTENT (IN)A
      INTENT (OUT)R
      END SUBROUTINE MCPY
      END INTERFACE
      END MODULE S_MCPY
!*==s_valvec.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_VALVEC
      IMPLICIT NONE
!*--S_VALVEC159
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE VALVEC(D,E,A,IC,N,ID)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: IC , ID , N
      REAL(R8KIND) , DIMENSION(ID,1) :: A
      REAL(R8KIND) , DIMENSION(1) :: D , E
      INTENT (IN)ID , N
      INTENT (OUT)IC
      INTENT (INOUT)A , D , E
      END SUBROUTINE VALVEC
      END INTERFACE
      END MODULE S_VALVEC
!*==s_lok.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_LOK
      IMPLICIT NONE
!*--S_LOK181
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE LOK(I,J,IR,N,M,MS)                                     !  LOK    5
      IMPLICIT NONE
      INTEGER :: I , IR , J , M , MS , N
      INTENT (IN)I , J , MS , N
      INTENT (OUT)IR
      END SUBROUTINE LOK
      END INTERFACE
      END MODULE S_LOK
!*==s_househ.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_HOUSEH
      IMPLICIT NONE
!*--S_HOUSEH199
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE HOUSEH(A,ID,N,D,E)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: ID , N
      REAL(R8KIND) , DIMENSION(ID,1) :: A
      REAL(R8KIND) , DIMENSION(1) :: D , E
      INTENT (IN)ID , N
      INTENT (INOUT)A , D , E
      END SUBROUTINE HOUSEH
      END INTERFACE
      END MODULE S_HOUSEH
!*==s_rosenbrock.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_ROSENBROCK
      IMPLICIT NONE
!*--S_ROSENBROCK220
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ROSENBROCK(FX,FY,FZ)
      USE PMIN_MODULE , ONLY:NTX_MIN , NTY_MIN , NTZ_MIN , NTX_MAX ,    &
     &    NTY_MAX , NTZ_MAX , NMOL , NATM , IDX , ISYSTEM , N , ILS ,   &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , STHL , SYM , Q , A1 , B1 , C1 , CK , WT_MOL ,   &
     &    DMAX , DDMAX , GNI , CE12 , VE12 , K1 , N11 , IHKL , IRIJ ,   &
     &    IRSS_CALL , PI , DIFFPC , PRDIFPC , NCYC_RSS                                   ! 4-23-08
      USE RSS1_MODULE
      USE RSS2_MODULE
      USE RSS3_MODULE
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE RSSPM3_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , DIMENSION(MATM) :: FX , FY , FZ
      INTENT (INOUT)FX , FY , FZ
      END SUBROUTINE ROSENBROCK
      END INTERFACE
      END MODULE S_ROSENBROCK
!*==s_new_xyz_rss.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_NEW_XYZ_RSS
      IMPLICIT NONE
!*--S_NEW_XYZ_RSS251
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE NEW_XYZ_RSS(IDX_CODE,CELL1,X,Y,Z,CELL_NEW1,XNEW1,YNEW1,&
     &                       ZNEW1)
      USE PMIN_MODULE , ONLY:NMOL , NATM , IDX , ISYSTEM , N , ILS ,    &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , PI                                                               ! 5/13/08
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: IDX_CODE
      REAL(R8KIND) , DIMENSION(6) :: CELL1 , CELL_NEW1
      REAL(R8KIND) , DIMENSION(MATM) :: X , XNEW1 , Y , YNEW1 , Z ,     &
     &                                  ZNEW1
      INTENT (IN)CELL1 , IDX_CODE
      INTENT (OUT)XNEW1 , YNEW1 , ZNEW1
      INTENT (INOUT)CELL_NEW1
      END SUBROUTINE NEW_XYZ_RSS
      END INTERFACE
      END MODULE S_NEW_XYZ_RSS
!*==s_new_xyz_trigonal.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_NEW_XYZ_TRIGONAL
      IMPLICIT NONE
!*--S_NEW_XYZ_TRIGONAL278
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE NEW_XYZ_TRIGONAL(ICODE,N,A,C1,DP,ATOM,X,Y,Z,CELL,XYZ)
      USE PMIN_MODULE , ONLY:MATM , PI                         ! 5/13/08
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) :: A , C1 , DP
      INTEGER :: ICODE , N
      CHARACTER(6) , DIMENSION(MATM) :: ATOM
      REAL(R8KIND) , DIMENSION(6) :: CELL
      REAL(R8KIND) , DIMENSION(MATM) :: X , Y , Z
      REAL(R8KIND) , DIMENSION(MATM,3) :: XYZ
      INTENT (IN)ATOM , C1 , DP , ICODE , N , X , Y , Z
      INTENT (INOUT)A , CELL , XYZ
      END SUBROUTINE NEW_XYZ_TRIGONAL
      END INTERFACE
      END MODULE S_NEW_XYZ_TRIGONAL
!*==s_unique_sym.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_UNIQUE_SYM
      IMPLICIT NONE
!*--S_UNIQUE_SYM303
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE UNIQUE_SYM(SYM,X,Y,Z,NSYM1)
      USE PMIN_MODULE , ONLY:NMOL , NATM , IDX , ISYSTEM , N , ILS ,    &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , MATM , ATOM                    ! 4-18-08
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , DIMENSION(20) :: NSYM1
      REAL(R8KIND) , DIMENSION(20,12) :: SYM
      REAL(R8KIND) , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN)SYM , X , Y , Z
      INTENT (INOUT)NSYM1
      END SUBROUTINE UNIQUE_SYM
      END INTERFACE
      END MODULE S_UNIQUE_SYM
!*==s_multiplicity.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_MULTIPLICITY
      IMPLICIT NONE
!*--S_MULTIPLICITY326
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE MULTIPLICITY(SYM,ATWT,X,Y,Z,MCITY,WT_MOL_ASYM)
      USE PMIN_MODULE , ONLY:NMOL , NATM , IDX , ISYSTEM , N , ILS ,    &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , MATM , ATOM                    ! 4-18-08
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) :: WT_MOL_ASYM
      REAL(R8KIND) , DIMENSION(MATM) :: ATWT , X , Y , Z
      INTEGER , DIMENSION(MATM) :: MCITY
      REAL(R8KIND) , DIMENSION(20,12) :: SYM
      INTENT (IN)ATWT , SYM , X , Y , Z
      INTENT (INOUT)MCITY , WT_MOL_ASYM
      END SUBROUTINE MULTIPLICITY
      END INTERFACE
      END MODULE S_MULTIPLICITY
!*==s_bend_bond_new.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_BEND_BOND_NEW
      IMPLICIT NONE
!*--S_BEND_BOND_NEW350
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE BEND_BOND_NEW(IVR,N,CELL,X,Y,Z,XO,YO,ZO,EGYMOL)
      USE PMIN_MODULE , ONLY:MATM , ATOM , PI , NMOL                                      ! 4-23-08
      USE RSS4_MODULE
      USE RSS5_MODULE
      USE BEND1_MODULE
      USE RSSPM3_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) :: EGYMOL
      INTEGER :: IVR , N
      REAL(R8KIND) , DIMENSION(6) :: CELL
      REAL(R8KIND) , DIMENSION(MATM) :: X , XO , Y , YO , Z , ZO
      INTENT (IN)CELL , IVR , X , Y , Z
      INTENT (INOUT)XO , YO , ZO
      END SUBROUTINE BEND_BOND_NEW
      END INTERFACE
      END MODULE S_BEND_BOND_NEW
!*==s_cross_product.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_CROSS_PRODUCT
      IMPLICIT NONE
!*--S_CROSS_PRODUCT377
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CROSS_PRODUCT(KDX,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4, &
     &                         X5,Y5,Z5)
      USE PMIN_MODULE , ONLY:MATM , ATOM
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: KDX
      REAL(R8KIND) :: X1 , X2 , X3 , X4 , X5 , Y1 , Y2 , Y3 , Y4 , Y5 , &
     &                Z1 , Z2 , Z3 , Z4 , Z5
      INTENT (IN)KDX , X1 , X2 , X3 , X4 , Y1 , Y2 , Y3 , Y4 , Z1 , Z2 ,&
     &        Z3 , Z4
      INTENT (OUT)X5 , Y5 , Z5
      END SUBROUTINE CROSS_PRODUCT
      END INTERFACE
      END MODULE S_CROSS_PRODUCT
!*==s_atom_links.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_ATOM_LINKS
      IMPLICIT NONE
!*--S_ATOM_LINKS401
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ATOM_LINKS(N,CELL,X,Y,Z,K11,K2,LINKATOM)
      USE PMIN_MODULE , ONLY:MATM , ATOM
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: K2 , N
      REAL(R8KIND) , DIMENSION(6) :: CELL
      INTEGER , DIMENSION(4) :: K11
      INTEGER , DIMENSION(MATM) :: LINKATOM
      REAL(R8KIND) , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN)CELL , K11 , N , X , Y , Z
      INTENT (OUT)LINKATOM
      INTENT (INOUT)K2
      END SUBROUTINE ATOM_LINKS
      END INTERFACE
      END MODULE S_ATOM_LINKS
!*==s_energy_pm3.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_ENERGY_PM3
      IMPLICIT NONE
!*--S_ENERGY_PM3426
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ENERGY_PM3(N,XO,YO,ZO,EGYMOL)
      USE PMIN_MODULE , ONLY:MATM , ATOM                                                  ! 9-2-08
      USE RSSPM3_MODULE
      USE RSS5_MODULE , ONLY:I_NV , E_WEIGHT
      USE BEND1_MODULE , ONLY:EGYMOL0
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) :: EGYMOL
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(MATM) :: XO , YO , ZO
      INTENT (IN)N , XO , YO , ZO
      INTENT (INOUT)EGYMOL
      END SUBROUTINE ENERGY_PM3
      END INTERFACE
      END MODULE S_ENERGY_PM3
!*==s_energy_b3lyp.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_ENERGY_B3LYP
      IMPLICIT NONE
!*--S_ENERGY_B3LYP451
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ENERGY_B3LYP(N,XO,YO,ZO,EGYMOL)
      USE PMIN_MODULE , ONLY:MATM , ATOM                                                  ! 9-2-08
      USE RSSPM3_MODULE
      USE RSS5_MODULE , ONLY:I_NV , E_WEIGHT
      USE BEND1_MODULE , ONLY:EGYMOL0
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER :: E_CONVERT = 627.5095
      REAL(R8KIND) :: EGYMOL
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(MATM) :: XO , YO , ZO
      INTENT (IN)N , XO , YO , ZO
      INTENT (INOUT)EGYMOL
      END SUBROUTINE ENERGY_B3LYP
      END INTERFACE
      END MODULE S_ENERGY_B3LYP
!*==s_ortho_cod.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_ORTHO_COD
      IMPLICIT NONE
!*--S_ORTHO_COD480
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ORTHO_COD(N,CELL,X,Y,Z,XO,YO,ZO)
      USE PMIN_MODULE , ONLY:MATM , ATOM
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(6) :: CELL
      REAL(R8KIND) , DIMENSION(MATM) :: X , XO , Y , YO , Z , ZO
      INTENT (IN)CELL , N , X , Y , Z
      INTENT (OUT)XO , YO , ZO
      END SUBROUTINE ORTHO_COD
      END INTERFACE
      END MODULE S_ORTHO_COD
!*==s_fractional_cod.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_FRACTIONAL_COD
      IMPLICIT NONE
!*--S_FRACTIONAL_COD502
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE FRACTIONAL_COD(N,CELL,XO,YO,ZO,X,Y,Z)
      USE PMIN_MODULE , ONLY:MATM , ATOM
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(6) :: CELL
      REAL(R8KIND) , DIMENSION(MATM) :: X , XO , Y , YO , Z , ZO
      INTENT (IN)CELL , N , XO , YO , ZO
      INTENT (INOUT)X , Y , Z
      END SUBROUTINE FRACTIONAL_COD
      END INTERFACE
      END MODULE S_FRACTIONAL_COD
!*==s_pc_limits.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_PC_LIMITS
      IMPLICIT NONE
!*--S_PC_LIMITS524
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE PC_LIMITS(NMOL,IROT,NVR,NV,PCITM,PCMIN,PCMAX)
      USE PMIN_MODULE , ONLY:NMOLD , NVD , PI                               ! 5/13/08
      USE F77KINDS
      IMPLICIT NONE
      INTEGER :: IROT , NMOL , NV , NVR
      REAL(R8KIND) , DIMENSION(NVD) :: PCITM , PCMAX , PCMIN
      INTENT (IN)IROT , NMOL , NV , NVR , PCITM
      INTENT (INOUT)PCMAX , PCMIN
      END SUBROUTINE PC_LIMITS
      END INTERFACE
      END MODULE S_PC_LIMITS
!*==s_user_e_1.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
      MODULE S_USER_E_1
      IMPLICIT NONE
!*--S_USER_E_1545
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE USER_E_1(CELL,X,Y,Z,EC,EV,ER,E_NMOL,PE)
      USE PMIN_MODULE , ONLY:NTX_MIN , NTY_MIN , NTZ_MIN , NTX_MAX ,    &
     &    NTY_MAX , NTZ_MAX , NMOL , NATM , IDX , ISYSTEM , N , ILS ,   &
     &    NSYM , ICYCLE , NENTRY , IMODE , NMOLD , NVD , DEL_PARAM ,    &
     &    MATM , ATOM , STHL , SYM , Q , A1 , B1 , C1 , CK , WT_MOL ,   &
     &    DMAX , DDMAX , GNI , CE12 , VE12 , K1 , N11 , IHKL , IRIJ ,   &
     &    KK , I_CROSS , A_CROSS , B_CROSS , C_CROSS , ANAME1 , ANAME2 ,&
     &    ATOM_CODE , A12 , B12 , C12 , PI , TWOPI , ANG_TO_RADIAN                      ! 4-23-08
      USE RSSPM3_MODULE
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) :: EC , ER , EV , E_NMOL , PE
      REAL(R8KIND) , DIMENSION(6) :: CELL
      REAL(R8KIND) , DIMENSION(MATM) :: X , Y , Z
      INTENT (IN)CELL , X , Y , Z
      INTENT (INOUT)EC , ER , EV , E_NMOL , PE
      END SUBROUTINE USER_E_1
      END INTERFACE
      END MODULE S_USER_E_1
