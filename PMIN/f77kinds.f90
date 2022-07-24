      MODULE F77KINDS
!
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER :: I1KIND = 1 , I2KIND = 2 , I4KIND = 4 ,     &
     &                       L1KIND = 1 , L4KIND = 4 , R4KIND = 4 ,     &
     &                       R8KIND = 8 , DPKIND = 8 , CX8KIND = 4 ,    &
     &                       CX16KIND = 8
!
!
! This module was generated automatically by MKKIND.F90
!    (supplied with the plusFORT toolkit from Polyhedron Software).
! It maps common non-standard Fortran 77 types to Fortran 90 kinds.
! Because different Fortran 90 compilers use different kind
! numbers, MKKIND.F90 should be compiled and run with the
! compiler you intend to use.
!
                        ! ! INTEGER*1
                        ! ! INTEGER*2
                        ! ! INTEGER*4
                        ! ! LOGICAL*1
                        ! ! LOGICAL*4
                        ! ! REAL*4
                        ! ! REAL*8
                        ! ! DOUBLE PRECISION
                        ! ! COMPLEX*8
      END MODULE F77KINDS ! COMPLEX*16
