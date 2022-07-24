      MODULE PMIN_MODULE
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , PARAMETER :: MATM = 200 , NMOLD = 4 ,                   &
     &                       NVD = 6 + 6*NMOLD + 20
      REAL(R8KIND) , PARAMETER :: PI = 3.141592653589793D0 ,            &
     &                            TWOPI = 2.0*PI ,                      &
     &                            ANG_TO_RADIAN = PI/180.0
!
! Local variables
!
      REAL(R8KIND) , DIMENSION(MATM) :: a1 , a_cross , b1 , b_cross ,   &
     &                                  c1 , c_cross , q, dd_cross      ! 5-4-10
      REAL(R8KIND) :: a12 , b12 , c12 , ce12 , ck , ddmax , dmax ,      &
     &                sthl , ve12 , wt_mol, f_factor_rss
      CHARACTER(6) , DIMENSION(MATM) :: aname1 , aname2 , atom ,        &
     &                                  atom_code
      CHARACTER(15) , DIMENSION(42) ::  vab_name
      CHARACTER(10) , DIMENSION(3)  ::  bond_type
      REAL(R8KIND) , DIMENSION(NVD) :: del_param , diffpc , prdifpc
      REAL(R8KIND) , DIMENSION(10000) :: gni
      INTEGER :: icycle , ils , imode , irss_call , isystem , i_cross , &
     &           k1 , kk , n , n11 , ncyc_rss , nentry , nmol , nsym ,  &
     &           ntx_max , ntx_min , nty_max , nty_min , ntz_max ,      &
     &           ntz_min, ispl, ispl_1(5)
      INTEGER , DIMENSION(NVD) :: idx
      real(r8kind), dimension(6) :: cell_input
      REAL(R8KIND) , DIMENSION(10000,3) :: ihkl
      REAL(R8KIND) , DIMENSION(50000,6) :: irij
      INTEGER , DIMENSION(NMOLD) :: natm
      REAL(R8KIND) , DIMENSION(20,12) :: sym
      DATA vab_name/'a length, Angs ',  'b length, Angs ',  'c length, Angs ',  &
     &              'alpha,    deg  ',  'beta,     deg  ',  'gamma,    deg  ',  &
     &            3*'RB1 rotn,  deg ',3*'RB1 trans, Angs',3*'RB2 rotn,  deg ',3*'RB2 trans, Angs', &
     &            3*'RB3 rotn,  deg ',3*'RB3 trans, Angs',3*'RB4 rotn,  deg ',3*'RB4 trans, Angs', &
     &            3*'RB5 rotn,  deg ',3*'RB5 trans, Angs',3*'RB6 rotn,  deg ',3*'RB6 trans, Angs'/
      DATA bond_type/'bond rota ','bond bend ','angle bend'/
      END MODULE PMIN_MODULE                                                               ! 9-16-08 DU
