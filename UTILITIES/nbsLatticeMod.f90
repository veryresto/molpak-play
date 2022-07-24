!  Programming notes:  1/16/03.  Today I compared the variable names
!  contained within all named common blocks among all subroutines, to
!  ensure that variable names were the same throughout the program, and
!  thus could be used with confidence with "use commonMod, only : "
!  statements
! -----------------------------------------------------------------------

 module nbsLatticeMod

 implicit none


!      INTEGER POINT,IPOINT
!      common /data1/ idata1,idata2,idata3
!      common /point1/ point(62),ipoint(62)
!      common /srch1/ istart(900),nogr(900)
!      common /unit1/ iunita
!      common /unit2/ iunitb
!      common /unit3/ iunitc
!      common /unit4/ iunitd
!      data idata1/55/, idata2/25/, idata3/24/
!      data iunita/15/, iunitb/16/, iunitc/7/, iunitd/10/
!      data istart/900*0/, nogr/900*0/
!From center.f
!      common /cntr1/ icntr1
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!From ckpt02.f90
!      common /dotp/ s11,s22,s33,s23,s13,s12
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /prob1/ ntd,iprob
!      common /type/ itype
!      common /unit2/ iunitb
!From deriv.f90
!      common /delta1/ irss
!      common /delta3/ idet,ilatno
!      common /delta4/ ideriv,ifinis,iqmatf
!      common /inver1/ ui1,vi1,wi1,ui2,vi2,wi2,ui3,vi3,wi3
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /qmat/ isq11(519),isq12(519),isq13(519),
!     $              isq22(519),isq23(519),isq33(519)
!From determ.f90
!      common /err2/ ierr2
!      common /fcn1/ ifcn
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /unit2/ iunitb
!From head0.f90
!     common /unit2/ iunitb
!From head2.f90
!      common /unit2/ iunitb
!From head3.f90
!      common /unit2/ iunitb
!From head4.f90
!      common /unit2/ iunitb
!From invers.f90
!      common /inver1/ ui1,vi1,wi1,ui2,vi2,wi2,ui3,vi3,wi3
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!Fromlmpre.f90
!      common /lmcell/ ad(900),bd(900),cd(900),alphad(900),betad(900),
!     $                gammad(900),vd(900),iseqd(900)
!      common /max1/ iamax(900),ibmax(900),icmax(900),ivmax(900)
!      common /min1/ iamin(900),ibmin(900),icmin(900),ivmin(900)
!      common /prob2/ ipch,nunk
!      common /prob3/ ie
!      common /toler1/ tol,tolv
!      common /unit2/ iunitb
!      common /unit4/ iunitd
!      common /ck01/ ick011,ick012,ick013,ick014
!From lmseto.f90
!      common /data1/ idata1,idata2,idata3
!      common /point1/ point(62),ipoint(62)
!      common /srch1/ istart(900),nogr(900)
!From lmsrch.f90
!      common /data1/ idata1,idata2,idata3
!      common /lmcell/ ad(900),bd(900),cd(900),alphad(900),betad(900),
!     $                gammad(900),vd(900),iseqd(900)
!      common /max1/ iamax(900),ibmax(900),icmax(900),ivmax(900)
!      common /min1/ iamin(900),ibmin(900),icmin(900),ivmin(900)
!      common /point1/ point(62),ipoint(62)
!      common /prob3/ ie
!      common /srch1/ istart(900),nogr(900)
!      common /unit2/ iunitb
!      common /unit3/ iunitc
!      common /ck01/ ick011,ick012,ick013,ick014
!from mncond.f90
!      common /dotp/ s11,s22,s33,s23,s13,s12
!      common /matr2/ u(9),t(9)
!      common /type/ itype
!      common /unit2/ iunitb
!      common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029
!from multip.f90
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /matr2/ u(9),t(9)
!from normal.f90
!      common /cosang/ cosa,cosb,cosg
!      common /matr2/ u(9),t(9)
!      common /type/ itype
!      common /var1/ var90
!      common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029
!From outpt1.f90
!      common /celli/ ai,bi,ci,alphi,beti,gami,voli
!      common /cellf/ af,bf,cf,alphf,betf,gamf,volf
!      common /cntr1/ icntr1
!      common /delta1/ irss
!      common /delta3/ idet,ilatno
!      common /dotp/ s11,s22,s33,s23,s13,s12
!      common /inver1/ ui1,vi1,wi1,ui2,vi2,wi2,ui3,vi3,wi3
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /print1/ iprrss
!      common /prob1/ ntd,iprob
!      common /unit2/ iunitb
!      common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029
!fromqmatri.f90
!      common /delta2/ idel1,idel2
!      common /delta4/ ideriv,ifinis,iqmatf
!      common /qmat/ isq11(519),isq12(519),isq13(519),
!     $              isq22(519),isq23(519),isq33(519)
!from read.f90
!      common /celli/ ai,bi,ci,alphi,beti,gami,voli
!      common /cntr1/ icntr1
!      common /delta1/ irss
!      common /delta2/ idel1,idel2
!      common /err1/ ierr1
!      common /unit1/ iunita
!      common /unit2/ iunitb
!from reduce.f90 
!      common /celli/ ai,bi,ci,alphi,beti,gami,voli
!      common /cellf/ af,bf,cf,alphf,betf,gamf,volf
!      common /delta1/ irss
!      common /delta4/ ideriv,ifinis,iqmatf
!      common /err1/ ierr1
!      common /prob1/ ntd,iprob
!      common /var1/ var90
!from set.f90
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /matr2/ u(9),t(9)
!from shortv.f90
!      common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029
!from spcon2.f90
!      common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029
!from spcond.f90
!      common /cosang/ cosa,cosb,cosg
!      common /dotp/ s11,s22,s33,s23,s13,s12
!      common /matr2/ u(9),t(9)
!      common /type/ itype
!      common /unit2/ iunitb
!      common /var1/ var90
!from trans.f90
!      common /celli/ ai,bi,ci,alphi,beti,gami,voli
!      common /cellf/ af,bf,cf,alphf,betf,gamf,volf
!      common /const1/ radian
!      common /cosang/ cosa,cosb,cosg
!      common /dotp/ s11,s22,s33,s23,s13,s12
!      common /err2/ ierr2
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /prob2/ ipch,nunk
!      common /unit2/ iunitb
!      common /unit4/ iunitd
!      common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029
!from volume.f90
!      common /const1/ radian
!      common /err1/ ierr1
!      common /unit2/ iunitb
!from main program
!      common /celli/ ai,bi,ci,alphi,beti,gami,voli
!      common /const1/ radian
!      common /err1/ ierr1
!      common /err2/ ierr2
!      common /fcn1/ ifcn
!      common /matr1/ u1,v1,w1,u2,v2,w2,u3,v3,w3
!      common /print1/ iprrss
!      common /prob1/ ntd,iprob
!      common /prob2/ ipch,nunk
!      common /toler1/ tol,tolv
!      common /unit1/ iunita
!      common /unit2/ iunitb
!      common /unit3/ iunitc
!      common /unit4/ iunitd
!**
!      common /ck01/ ick011,ick012,ick013,ick014
 !     common /ck02/ ick021,ick022,ick023,ick024,ick025,ick026,ick027,
!     $              ick028,ick029


 integer :: i
 integer :: iamax(900),ibmax(900),icmax(900),ivmax(900)
 integer :: iamin(900),ibmin(900),icmin(900),ivmin(900)
 integer :: ick011,ick012,ick013,ick014
 integer :: ick021,ick022,ick023,ick024,ick025,ick026,ick027, &
          & ick028,ick029
 integer :: icntr1
 integer :: idata1=55
 integer :: idata2=25
 integer :: idata3=24
 integer :: idel1, idel2
 integer :: ideriv
 integer :: idet
 integer :: ie
 integer :: ierr1
 integer :: ierr2
 integer :: ifcn
 integer :: ifinis
 integer :: ilatno
 integer :: iqmatf
 integer :: ipch
 integer :: ipoint(62)=&
     &(/1,503,990,1501,2007,2498,3010,3527,3987,4473,4976,&
     &5479,5983,6469,7000,7476,8007,8499,9026,9508,10015,10475,11021,&
     &11518,12021,12486,13001,13487,13992,14483,14977,15471,16012,16506,&
     &17013,17493,17993,18506,18976,19500,20011,20492,21007,21491,21994,&
     &22505,22989,23493,23999,24507,25001,25498,26000,26500,26988,0,0,0,&
     &0,0,0,0/)
 integer :: iprob
 integer :: iprrss
 integer :: irss
 integer :: iseqd(900)
 integer :: isq11(519),isq12(519),isq13(519),&
          & isq22(519),isq23(519),isq33(519)
 integer :: istart(900)=(/ (0,i=1,900) /)
 integer :: itype
 integer :: iunita=15
 integer :: iunitb=16
 integer :: iunitc=7
 integer :: iunitd=10
 integer :: nogr(900)=(/ (0, i=1,900) /)
 integer :: ntd
 integer :: nunk
 integer :: point(62)=&
     &(/207,292,308,328,345,359,370,379,384,390,396,401,407,&
     &413,422,431,444,458,472,482,494,502,511,518,525,531,539,548,558,&
     &568,577,588,599,612,626,640,655,673,688,702,715,727,740,753,773,&
     &796,819,850,891,942,1007,1055,1125,1272,3840,0,0,0,0,0,0,0/)

 real :: ai,bi,ci,alphi,beti,gami,voli
 real :: af,bf,cf,alphf,betf,gamf,volf
 real :: ad(900),bd(900),cd(900),alphad(900),betad(900),gammad(900),vd(900)
 real :: cosa, cosb, cosg
 real :: radian
 real :: s11,s22,s33,s23,s13,s12
 real :: tol, tolv
 real :: u1,v1,w1,u2,v2,w2,u3,v3,w3
 real :: ui1,vi1,wi1,ui2,vi2,wi2,ui3,vi3,wi3
 real :: u(9), t(9) 
 real :: var90

 end module nbsLatticeMod
