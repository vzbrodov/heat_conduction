program Qz
implicit none
integer i, m
real, allocatable :: A(:,:)
open(unit=1,file='Qs')

m=10

allocate(A(m,m))
write(1,*) '#',m
A=0
!A(50,2)=-1000
!A(50,99)=-1000
!A(2,50)=1000
!A(99,50)=1000
!A(50,45)=-40
A(50,2)=2
A(2,50)=-2
A(40,60)=-3
A(40,100)=3
!(m/2-10,m/2+10)=-1
!(m/2+10,m/2-10)=1
!A(10,10)=-2


do i=1,m
    write(1,*) A(1:m,i)
enddo
end program Qz
