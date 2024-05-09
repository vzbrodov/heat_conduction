program heat_conduction  
    use :: mpi
implicit none

real, dimension(:,:), allocatable :: qh2, u_old, u_new, u_matrix, ql
integer,dimension(MPI_STATUS_SIZE) :: STAT
integer, parameter :: IN = 1, OUT = 2
integer, parameter :: MAX_ITER = 100, OUTPUT_FREQ = 10
integer :: i, j, k, n, ierr, my_rank, num_proc,u,d,l,r
character(len=15) :: filename

call MPI_Init(ierr)

STAT(MPI_SOURCE)=MPI_ANY_SOURCE
STAT(MPI_TAG)=MPI_ANY_TAG
STAT(MPI_ERROR)=ierr

call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierr)

if (num_proc .eq. 0) then
    open(IN, file='./sources.dat')  
    read(IN,'(2x,i10)') n
    allocate(qh2(n,n))
    read(IN,*) qh2
    close(IN)
end if
call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

allocate(u_old(n/2+2,n/2+2))
allocate(u_new(n/2+2,n/2+2))
allocate(ql(n/2,n/2))
allocate(u_matrix(n,n))


if (num_proc .eq. 0) ql=qh2(1:n/2,1:n/2)
if (num_proc .eq. 1) ql=qh2(1:n/2,n/2+1:n)
if (num_proc .eq. 2) ql=qh2(n/2+1:n,1:n/2)
if (num_proc .eq. 3) ql=qh2(n/2+1:n,n/2+1:n)
u_old=0
u_new=0
u_matrix=0

do k = 1, MAX_ITER

        if (num_proc .eq. 2) then
             call MPI_Send(u_old(2,2:n/2+1), n/2, MPI_REAL, 0, 11, MPI_COMM_WORLD, ierr)
             call MPI_Send(u_old(2:n/2+1,n/2+1), n/2, MPI_REAL, 3, 88, MPI_COMM_WORLD, ierr)  
        endif
        if (num_proc .eq. 3) then
            call MPI_Send(u_old(n/2+1,2:n/2+1), n/2, MPI_REAL, 1, 33, MPI_COMM_WORLD, ierr)
            call MPI_Send(u_old(2:n/2+1,2), n/2, MPI_REAL, 2, 66, MPI_COMM_WORLD, ierr)
        endif
         if (num_proc .eq. 1) then
            call MPI_Send(u_old(2:n/2+1,1), n/2, MPI_REAL, 0, 22, MPI_COMM_WORLD, ierr)
            call MPI_Send(u_old(n/2+1,2:n/2+1), n/2, MPI_REAL, 3, 77, MPI_COMM_WORLD, ierr) 
        endif 
        if (num_proc .eq. 0) then
            call MPI_Send(u_old(2:n/2+1,n/2+1), n/2, MPI_REAL, 1, 44, MPI_COMM_WORLD, ierr)
            call MPI_Send(u_old(n/2+1,2:n/2+1), n/2, MPI_REAL, 2, 55, MPI_COMM_WORLD, ierr)
        endif 

        if (num_proc .eq. 0) then
           call MPI_Recv(u_old(2:n/2+1,n/2+2), n/2, MPI_REAL, 1, 22, MPI_COMM_WORLD, STAT,ierr)
           call MPI_Recv(u_old(n/2+2,2:n/2+1), n/2, MPI_REAL, 2, 11, MPI_COMM_WORLD, STAT,ierr)
        endif 

        if (num_proc .eq. 1) then
           call MPI_Recv(u_old(2:n/2+1,1), n/2, MPI_REAL, 0, 44, MPI_COMM_WORLD, STAT,ierr)
            call MPI_Recv(u_old(n/2+2,2:n/2+1), n/2, MPI_REAL, 3, 33, MPI_COMM_WORLD, STAT,ierr)
        endif
        if (num_proc .eq. 2) then
            call MPI_Recv(u_old(1,2:n/2+1), n/2, MPI_REAL, 0, 55, MPI_COMM_WORLD, STAT,ierr)
           call MPI_Recv(u_old(2:n/2+1,n/2+2), n/2, MPI_REAL, 3, 66, MPI_COMM_WORLD, STAT,ierr)
        endif
        if (num_proc .eq. 3) then
            call MPI_Recv(u_old(1,2:n/2+1), n/2, MPI_REAL, 1, 77, MPI_COMM_WORLD, STAT,ierr)
            call MPI_Recv(u_old(2:n/2+1,1), n/2, MPI_REAL, 2, 88, MPI_COMM_WORLD, STAT,ierr)
        endif
         

        do i=2,n/2+1
            write(*,*)'111',k
            do j=2, n/2+1
                write(*,*)'222',k
                u_new(i,j)=1.0/4.0*(u_old(i-1,j)+u_old(i,j-1)+u_old(i+1,j)+u_old(i,j+1)+ql(i-1,j-1))  
            enddo
        enddo

        u_old=u_new

        if (num_proc .ne. 0) then
        call MPI_Send(u_old(2:n/2+1,2:n/2+1), n*n/4, MPI_REAL, 0, 100, MPI_COMM_WORLD, ierr)
    endif
        
        if (num_proc .eq. 0) then
            if (mod(k, OUTPUT_FREQ) == 0) then
                u_matrix(1:n/2,1:n/2) = u_old(2:n/2+1,2:n/2+1)        
                    call MPI_Recv(u_matrix(n/2+1:n,1:n/2), n*n/4, MPI_REAL, 1, 100, MPI_COMM_WORLD, STAT,ierr)
                    call MPI_Recv(u_matrix(1:n/2,n/2+1:n), n*n/4, MPI_REAL, 2, 100, MPI_COMM_WORLD, STAT,ierr)
                    call MPI_Recv(u_matrix(n/2+1:n,n/2+1:n), n*n/4, MPI_REAL, 3 , 100, MPI_COMM_WORLD, STAT,ierr)

                write(filename,'(A,I3.3,A)') "DATA/data", k, ".dat"
                open(OUT, file=trim(filename))

                do i = 1, n
                    write(OUT,'((4X,F16.8))') u_matrix(i,:)
                end do
                close(OUT)
            endif
        endif

end do

call MPI_Finalize(ierr)

end program heat_conduction
