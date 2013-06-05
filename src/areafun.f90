subroutine areafun(A,m,out)

implicit none
integer :: m,i
double precision :: A(m,12),out(m),z(3),AC(3),BD(3)

do i=1,m
   AC = A(i,1:3)-A(i,4:6)
   BD = A(i,7:9)-A(i,10:12)
   call crosspr(AC,BD,z)
out(i) = 0.5*sqrt(dot_product(z,z))
end do

contains
subroutine crosspr (x,y,z)
   real*8 :: x(3),y(3),z(3),lz
        
   z(1) = x(2)*y(3)-x(3)*y(2)
   z(2) = x(3)*y(1)-x(1)*y(3)
   z(3) = x(1)*y(2)-x(2)*y(1)
   lz = sqrt(dot_product(z,z))
   
    end subroutine crosspr
end subroutine areafun
