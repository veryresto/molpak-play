!*==rss5_module.spg  processed by SPAG 6.56Rc at 11:14 on 24 Nov 2008
 
 
      MODULE RSS5_MODULE
!
      USE F77KINDS                        
      IMPLICIT NONE
!*--RSS5_MODULE8
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
      REAL*8 :: e_weight
      INTEGER , DIMENSION(20,100) :: iatom , linkatom
      INTEGER , DIMENSION(20) :: indx , k5
      INTEGER :: irot , i_nv , ne_type , nrotbond
!
!*** End of declarations rewritten by SPAG
!
                                                              !   ! (20,4) changed to
                                            !                     ! (20,100)
                                                                  ! 9-30-08 DU
!
      END MODULE RSS5_MODULE                                      ! 9-2-08 DU
