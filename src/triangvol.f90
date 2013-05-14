subroutine trianvol(vb1,vb2,it,dimvb,dimit,V)
  implicit none
  integer :: i,j,dimvb,dimit,it(3,dimit),subvols(3,4)
  real*8 :: vb1(4,dimvb),vb2(4,dimvb),allvol(4,6), V,tmp,Vtmp

  subvols(1,:) = (/1,2,3,4/)
  subvols(2,:) = (/1,2,4,5/)
  subvols(3,:) = (/2,4,5,6/)
  V = 0.0d0
  do i = 1,dimit
     Vtmp = 0.0d0
     allvol(:,1:3) = vb1(:,it(:,i))
     allvol(:,(/5,6,4/)) = vb2(:,it(:,i))
     
     
     do j = 1,3
        call detfour(transpose(allvol(:,subvols(j,:))),tmp)
        Vtmp = Vtmp + abs(tmp)
     end do
     V = V + Vtmp/6

  end do
end subroutine trianvol

subroutine detfour(A,det)
 real*8:: A(4,4),det
DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
             A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
             A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
             A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
             A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))


end subroutine detfour
