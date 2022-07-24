!
!  Program to check if two molecules are enantiomorphs 
!  May 5, 2005
!
!  Orthogonal coordinates assumed.
!
!       implicit real*8 (a-h,o-z)
        SUBROUTINE HAND(NMOD,AA,CC,XYZ,W,AT,AN,LH)
        logical :: pair
        character AT(200)*2, AN(200)*4
        character*6  atom1(100),atom2(100)
        dimension x1(100),y1(100),z1(100),x2(100),y2(100),z2(100), x11(100),y11(100),z11(100)
        dimension x12(100),y12(100),z12(100), x13(100),y13(100),z13(100), x14(100),y14(100),z14(100)
        dimension x22(100),y22(100),z22(100)
        DIMENSION AA(3),CC(3),XYZ(3,200),W(3,200),V(3,200),RV(3,3)
        common atom1,atom2
!        open(unit = 7, file = 'enantio.inp', status = 'old')
!       open(unit = 88, file = 'enantio.out', status = 'unknown')
!       open(unit = 87, file = 'molpak2.xyz', status = 'unknown')
!       open(unit = 10, file = 'input.cod1', status = 'old')
!       open(unit = 11, file = 'input.cod2', status = 'old')
!
! Now read the atom coordinates from unit 10
!
         LH = 0
         n = NMOD
         DO I = 1, NMOD
          x1(I) = XYZ(1,I)
          y1(I) = XYZ(2,I)
          Z1(I) = XYZ(3,I)
!         write(88,1) AT(I),AN(I),x1(I),y1(I),Z1(I)
         END DO 
 1       FORMAT(12x,a2,a4,2x,3f10.5)
!10      n = n + 1
!       read(10,1,end = 100)atom1(n),x1(n),y1(n),z1(n)
!       write(8,1)atom1(n),x1(n),y1(n),z1(n)
!       format(12x,a6,2x,3f10.5)
!          go to 10
!100     n = n-1
!       write(88,2)n
!       covert W to v  (transform from fractional coords. to Cartesian) 
!       WRITE(88,'(6F9.5)') AA,CC
         DO I = 1,3
          IF(CC(I).LT. 1.0) CC(I) = ACOSD(CC(I))
         END DO
!        WRITE(88,'(6F9.5)') AA,CC
        CALL MATRXT(AA,CC,RV)
      DO K=1,NMOD
      V(1,K)=RV(1,1)*W(1,K)+RV(1,2)*W(2,K)+RV(1,3)*W(3,K)
      V(2,K)=RV(2,1)*W(1,K)+RV(2,2)*W(2,K)+RV(2,3)*W(3,K)
      V(3,K)=RV(3,1)*W(1,K)+RV(3,2)*W(2,K)+RV(3,3)*W(3,K)
      END DO
2       format(/,5x,'No. of atoms ',i5,/)
          DO I = 1, NMOD
          x2(I) = V(1,I)
          y2(I) = V(2,I)
          Z2(I) = V(3,I)
!         write(88,3) AT(I),AN(I),(w(j,I),J=1,3),x2(I),y2(I),Z2(I)
!         write(87,5) AT(I),AN(I),x2(I),y2(I),Z2(I)
         END DO
3       format(12x,a2,a4,6f10.5)
5       format(5x,a2,a4,3f10.6)
!         do i = 1,n
!       read(11,3)  atom2(i),x2(i),y2(i),z2(i)
!       the following three lines is for testing 
!            x2(i) = -x2(i)
!            y2(i) = -y2(i)   
!            z2(i) = -z2(i)
!       write(8,1) atom2(i),x2(i),y2(i),z2(i)
!3       format(12x,a6,66x,3f10.5)
!         enddo
!        32  H_9        0.64549  -3.25263  -0.02195     1.20050  -3.09113  -0.02195    -1.20050  -3.09113   0.02195
!
200      continue
!
         CALL ALIGN(n,x1,y1,z1,x11,y11,z11)
         CALL ALIGN(n,x2,y2,z2,x22,y22,z22)
          do i = 1,n
!        write(88,4) atom1(i),x11(i),y11(i),z11(i),x22(i),y22(i),z22(i)
!4        format(5x,a6,2x,6f10.5)
!        write(88,4) AT(i),AN(i),x11(i),y11(i),z11(i),x22(i),y22(i),z22(i)
 4        format(5x,a2,a4,2x,6f10.5)
          enddo
            do i = 1,n
          sum_y1 = sum_y1 + y11(i)
          sum_y2 = sum_y2 + y22(i)
            enddo
          sum_y1y2 =  sum_y1 + sum_y2
          ratio = sum_y1y2/(abs(sum_y1) + abs(sum_y2))
!         write(88,111)sum_y1,sum_y2,sum_y1y2,ratio
!         print 111,sum_y1,sum_y2,sum_y1y2,ratio
!111      format(/,'sum_y1,sum_y2,sum_y1y2,ratio ',4f10.5,/)
!
          if(abs(ratio) < 0.1) then
          pair = .FALSE.
          LH = 1
!         write(88,*) '   The two molecules are enantiomorphs'
          print 111,sum_y1,sum_y2
 111      format('sum_y1, sum_y2 =',2f10.5,': molecules are enantiomorphs')
          else
          LH = 0
          pair = .TRUE.
!         write(88,*) '   The two molecules are not enantiomorph of each other'
          endif
!
!         write(88,*) '   Are the two molecules of the same hand:  ',pair
        RETURN
        END
!
        SUBROUTINE ALIGN(n,x1,y1,z1,x2,y2,z2)
!       implicit real*8 (a-h,o-z)
!
      dimension  x1(100),y1(100),z1(100), x2(100),y2(100),z2(100), x13(100),y13(100),z13(100)
      dimension  x11(100),y11(100),z11(100), x12(100),y12(100),z12(100)
      common atom1(100),atom2(100)
      pi = 3.141592653
!
      x11(2) = x1(2) - x1(1)
      y11(2) = y1(2) - y1(1)
      z11(2) = z1(2) - z1(1)
         if(z11(2) == 0.0) then
             alpha =0.5*pi
           go to 44
           endif
           alpha = atan(y11(2)/z11(2))
44     continue
!
!     Rotate the entire molecule by alpha about x-axis with atom 1 at center
!
       do i = 1,n
      x11(i) = x1(i) - x1(1)
      y11(i) = y1(i) - y1(1)
      z11(i) = z1(i) - z1(1)
      x12(i) = x11(i)
      y12(i) = y11(i)*cos(alpha) - z11(i)*sin(alpha)
      z12(i) = y11(i)*sin(alpha) + z11(i)*cos(alpha)
!     write(8,3) x22,y22,z22,x32,y32,z32
3    format(/,'after rotation by alpha about x-axis ',/,2(10x,3f9.5/))
     enddo
!
     if(z12(2) == 0.0) then
       beta =0.5*pi
       go to 45
     endif
       beta = -atan(x12(2)/z12(2))
45   continue
!
!     Rotate the entire molecule by beta about y-axis with atom 1 at center
!
       z13(2) =z12(2)*cos(beta) - x12(2)*sin(beta)
        if(z13(2) < 0.0) beta = pi + beta
       do i = 1,n
      x13(i) = z12(i)*sin(beta) + x12(i)*cos(beta)
      y13(i) = y12(i)
      z13(i) = z12(i)*cos(beta) - x12(i)*sin(beta)
       enddo
!
!     write(8,331) x23,y23,z23,x33,y33,z33
331    format(/,'after rotation by beta about y-axis ',/,2(10x,3f9.5/))
!
     if(x13(3) == 0.0) then
       gamma =0.5*pi
       go to 46
     endif
       gamma = -atan(y13(3)/x13(3))
46     continue
!
!     Rotate the entire molecule by gamma about z-axis to bring atom 3 in xz-plane
!
        x2(3) =x13(3)*cos(gamma) - y13(2)*sin(gamma)
        if(x2(3) < 0.0) gamma = pi + gamma
!
       do i = 1,n
!      x14(i) = x13(i)*cos(gamma) - y13(i)*sin(gamma)
!      y14(i) = x13(i)*sin(gamma) + y13(i)*cos(gamma)
!      z14(i) = z13(i)
       x2(i) = x13(i)*cos(gamma) - y13(i)*sin(gamma)
       y2(i) = x13(i)*sin(gamma) + y13(i)*cos(gamma)
       z2(i) = z13(i)
       enddo
!
!    Now write the result
!
      alpha1 = 180.0*alpha/pi
      beta1  = 180.0*beta/pi
      gamma1 = 180.0*gamma/pi
    write(88,15)alpha1,beta1,gamma1
15   format(/,5x,'rotation angles alpha,beta,gamma ',3f10.4,/)
!      do i= 1,n
!      write(8,12)atom1(i),x1(i),y1(i),z1(i),x12(i),y12(i),z12(i),x13(i),y13(i),z13(i),x2(i),y2(i),z2(i)
!12    format(2x,a6,2x,4(1x,3f9.5))
!      enddo
      RETURN
       end

      SUBROUTINE MATRXT (AA,CC,T)
      DIMENSION T(3,3),A(6),AA(3),CC(3)
      do n = 1,3
       a(n) = aa(n)
       a(n+3) = cc (n)
      end do
!      print '(6F9.5)',a

      T(1,1)=1.0*A(1)
      T(1,2)=COSD(A(6))*A(2)
      T(1,3)=COSD(A(5))*A(3)

      T(2,1)=0.0
      T(2,2)=SIND(A(6))*A(2)
      T(2,3)=A(3)*(COSD(A(4))-COSD(A(5))*COSD(A(6)))/SIND(A(6))

      W=SQRT(1-COSD(A(4))**2-COSD(A(5))**2-COSD(A(6))**2+ &
     & 2*COSD(A(4))*COSD(A(5))*COSD(A(6)))
      T(3,1)=0.0
      T(3,2)=0.0
      T(3,3)=A(3)*W/SIND(A(6))

      WRITE(88,100)((T(I,J),J=1,3),I=1,3)
  100 FORMAT(3X,3F8.5)
      RETURN
      END

