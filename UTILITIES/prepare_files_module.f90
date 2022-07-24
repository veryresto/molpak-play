      MODULE mod_preppot
!This module is generated for preppot.f90
!
      IMPLICIT NONE
!
!  COMMON /MOPAC_STUFF/
      CHARACTER (72) :: MTITLE
      character (1)  :: char_atom_order=' '
      INTEGER        :: MOPAC_TYPE,I_CH_TYPE
      integer        :: atom_order=0
!                
!  COMMON /CHEM3D_STUFF/ 
      CHARACTER (11) :: DATNM2 
      CHARACTER (19) :: OUTFIL, CHEM3D_NAME
      INTEGER        :: NCC
!
!  COMMON /COM1/ 
      LOGICAL :: CONNECT_FLAG
      INTEGER :: ICON(6,200), NCON(200)
!
!  COMMON /ADJUST2/
      CHARACTER (LEN=5) :: LABEL(500)
! RLABEL(500) removed, not in use
      INTEGER :: ISPECIAL, NATOMS, N_HY, N_OTHER, N_OX
      INTEGER :: NEW, N_2, N_CENTRAL
      REAL    :: CELL(6), C1(3)    
      REAL    :: OR_XYZ(3,500), R(3,3)
      REAL    :: XYZ(3,500), XYZR(3,3), XYZT(3), XMAT2(3,3)           
! 
!from subroutine ADJUST
! COMMON/T_MATRIX/ T
      REAL    :: T(3,3)
! COMMON /CONNECT/ N_NEAR, N_ATOMS, J1, J2
      INTEGER :: N_NEAR(100), N_ATOMS(4,100), J1(3), J2(3)
!
!from subroutine MAKE_CHEM3D
!      COMMON /CONNECT_INFO/ ICON, NCON, AT2
      CHARACTER (LEN=1)  :: AT2(500)
!    
      END MODULE mod_preppot
