program heat_conduction  
    use :: mpi
implicit none

real, dimension(:,:), allocatable :: qh2, u_old, u_new, u_matrix,ql 
integer,dimension(MPI_STATUS_SIZE) :: status
integer, parameter :: IN = 1, OUT = 2
integer, parameter :: MAX_ITER = 20000, OUTPUT_FREQ = 1000               !Формат принимаемого файла: в первой строке # и размер матрицы источников
integer :: n, i, j, k, ierr, my_rank                                    ! в следующих строках сама матрица
character(len=15) :: filename
                                                                        !OUTPUT_FREQ - частота вывода результата в файл
call MPI_Init(ierr)

call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

if (my_rank .eq. 0) then
    open(IN, file='./sources.dat') 
    read(IN, '(2x, i10)')n 
    allocate(qh2(n,n))
    read(IN,*) qh2
    close(IN)
endif

call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                                    ! нумерация рангов на сетке
  if (my_rank .eq. 0 ) write(*,*)"N=",n               !     | 0 | 1 |
                                                      !     ---------
                                                      !     | 3 | 2 |

allocate(u_old(n/2+2,n/2+2))                            
allocate(u_new(n/2+2,n/2+2))                            
allocate(u_matrix(n,n))                             
allocate(ql(n/2,n/2))
if (my_rank .eq. 0) then
    ql=qh2(1:n/2,1:n/2)
    call MPI_Send(qh2(1:n/2,n/2+1:n), n*n/4, MPI_REAL, 1, 910, MPI_COMM_WORLD, ierr)
    call MPI_Send(qh2(n/2+1:n,n/2+1:n), n*n/4, MPI_REAL, 2, 910, MPI_COMM_WORLD, ierr)
    call MPI_Send(qh2(2:n/2+1,n/2+1), n*n/4, MPI_REAL, 3, 910, MPI_COMM_WORLD, ierr)
endif
if (my_rank .ne. 0) call MPI_Recv(ql, n*n/4, MPI_REAL, 0, 910, MPI_COMM_WORLD, status,ierr) 
u_old=0
u_new=0
u_matrix=0
 write(*,*)"Программа работает..." 
do k = 1, MAX_ITER

    if (my_rank .eq. 0) then

        call MPI_Send(u_old(2:n/2+1,n/2+1), n/2, MPI_REAL, 1, 111, MPI_COMM_WORLD, ierr)
        call MPI_Send(u_old(n/2+1,2:n/2+1), n/2, MPI_REAL, 3, 222, MPI_COMM_WORLD, ierr)

        call MPI_Recv(u_old(2:n/2+1,n/2+2), n/2, MPI_REAL, 1, 333, MPI_COMM_WORLD, status,ierr) 
        call MPI_Recv(u_old(n/2+2,2:n/2+1), n/2, MPI_REAL, 3, 888, MPI_COMM_WORLD, status,ierr)

    
        if (mod(k, OUTPUT_FREQ) == 0) then
            u_matrix(1:n/2,1:n/2) = u_old(2:n/2+1,2:n/2+1) 
        endif
    endif



if (my_rank .eq. 1) then

    call MPI_Send(u_old(2:n/2+1,2), n/2, MPI_REAL, 0, 333, MPI_COMM_WORLD, ierr)
    call MPI_Send(u_old(n/2+1,2:n/2+1), n/2, MPI_REAL, 2, 444, MPI_COMM_WORLD, ierr)

    call MPI_Recv(u_old(2:n/2+1,1), n/2, MPI_REAL, 0, 111, MPI_COMM_WORLD, status,ierr)
    call MPI_Recv(u_old(n/2+2,2:n/2+1), n/2, MPI_REAL, 2, 666, MPI_COMM_WORLD, status,ierr)

endif

if (my_rank .eq. 2) then
    
    call MPI_Send(u_old(2:n/2+1,2), n/2, MPI_REAL, 3, 555, MPI_COMM_WORLD, ierr)
    call MPI_Send(u_old(2,2:n/2+1), n/2, MPI_REAL, 1, 666, MPI_COMM_WORLD, ierr)

    call MPI_Recv( u_old(1,2:n/2+1), n/2, MPI_REAL, 1, 444, MPI_COMM_WORLD, status,ierr)
    call MPI_Recv(u_old(2:n/2+1,1), n/2, MPI_REAL, 3, 777, MPI_COMM_WORLD, status,ierr)

endif



if (my_rank .eq. 3) then
    
    call MPI_Send(u_old(2:n/2+1,n/2+1), n/2, MPI_REAL, 2, 777, MPI_COMM_WORLD, ierr)
    call MPI_Send(u_old(2,2:n/2+1), n/2, MPI_REAL, 0, 888, MPI_COMM_WORLD, ierr)

    call MPI_Recv(u_old(1,2:n/2+1), n/2, MPI_REAL, 0, 222, MPI_COMM_WORLD, status,ierr)
    call MPI_Recv(u_old(2:n/2+1,n/2+2), n/2, MPI_REAL, 2, 555, MPI_COMM_WORLD, status,ierr)

endif


    do i=2,n/2+1
            do j=2, n/2+1
                u_new(i,j)=(u_old(i-1,j)+u_old(i,j-1)+u_old(i+1,j)+u_old(i,j+1)+ql(i-1,j-1))/4.0   
            enddo
        enddo
        
        u_old=u_new 
    if (my_rank .ne. 0) then
        if (mod(k, OUTPUT_FREQ) == 0) then
            call MPI_Send(u_old(2:n/2+1,2:n/2+1), n*n/4, MPI_REAL, 0, 100, MPI_COMM_WORLD, ierr)
        endif
    endif

        
    if (my_rank .eq. 0) then
        if (mod(k, OUTPUT_FREQ) == 0) then

            call MPI_Recv(u_matrix(1:n/2,n/2+1:n), n*n/4, MPI_REAL, 1, 100, MPI_COMM_WORLD, status,ierr)
            call MPI_Recv(u_matrix(n/2+1:n,n/2+1:n), n*n/4, MPI_REAL, 2 , 100, MPI_COMM_WORLD, status,ierr)
            call MPI_Recv(u_matrix(n/2+1:n,1:n/2), n*n/4, MPI_REAL, 3, 100, MPI_COMM_WORLD, status,ierr)

            write(filename,'(A,I5.5,A)')"data",k, ".dat"
            open(OUT, file=trim(filename))

            do i = 1, n
                write(OUT,*) u_matrix(i,:)
            end do
            close(OUT)
        endif
    endif

end do
    if (my_rank .eq. 0 ) write(*,*)"End." 
call MPI_Finalize(ierr)

end program heat_conduction
