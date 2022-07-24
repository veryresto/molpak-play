! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 module transformCommonMod
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------------------------------
!  Common blocks that are used to construct variables and parameters
!  for this module.

!     COMMON /COM1/ NATOMS, CELL(6), LABEL(200), O_XYZ(3,200),
!     x              N_NEAR(200), N_ATOMS(6,200)
!      COMMON /COM2/XYZ 
!   COMMON /COM3/ CONNECT_FLAG, ICHARGE_FLAG, G92_CHARGE,
!     X              XMNDO_CHARGE
!-----------------------------------------------------------------------

 implicit none
 integer, parameter :: single=kind(0.0e0)

 character(5) :: label(200)

 logical :: connect_flag
 
 integer :: icharge_flag
 integer :: natoms
 integer :: n_near(200)
 integer :: n_atoms(6,200)

 real (single) :: cell(6)=(/1.0, 1.0, 1.0, 90.0, 90.0, 90.0/)
 real (single) :: g92_charge
 real (single) :: o_xyz(3,200)
 real (single) :: xmndo_charge
 real (single) :: xyz(3,200)


 end module transformCommonMod

