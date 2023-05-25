program main
    implicit none

    character(len=64) :: filename, logfile
    real:: start, finish

    ! for input matrix
    integer :: nmax, nzmax, n, nnz, ierr, iounit, i
    real*8, dimension (:), allocatable :: a, rhs
    integer, dimension (:), allocatable :: ia, ja

    ! for ILU(k)
    real*8, dimension (:), allocatable :: alu, w
    integer, dimension (:), allocatable :: jlu, ju, levs, jw
    integer :: lfil, iwk

    ! for PGMRES
    real*8, dimension (:), allocatable :: sol, vv
    real*8 :: eps
    integer :: maxits, im, iout

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read matrix
    nmax = 264751
    nzmax = 3541545
    allocate(ia(nmax + 1))
    allocate(ja(nzmax), a(nzmax))

    call get_command_argument(1, filename)
    iounit = 1
    open(iounit, file=filename)

    call readsk(nmax, nzmax, n, nnz, a, ja, ia, iounit, ierr)

    print *, "Matrix size = ", n
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ILU(K)
    iwk = nnz * 100
    allocate(w(n), ju(n))
    allocate(jw(3 * n))
    allocate(alu(iwk), jlu(iwk), levs(iwk))

    lfil = 0
    call cpu_time(start)
    call iluk(n, a, ja, ia, lfil, alu, jlu, ju, levs, iwk, w, jw, ierr)
    call cpu_time(finish)

    print *, "Time for iluk = ", finish - start
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PGMRES
    allocate(rhs(n))
    do i = 1, n
        rhs(i) = sin(real(i))
    end do

    im = 12 ! amount of Krylov subspaces
    eps = 1e-7 ! residual decrease
    maxits = n * n 
    iout = 2

    call get_command_argument(2, logfile)
    open(iout, file=logfile)

    allocate(sol(n)) ! solution
    allocate(vv(n * (im + 1))) ! work space

    call cpu_time(start)
    call pgmres(n, im, rhs, sol, vv, eps, maxits, iout, a, ja, ia, alu, jlu, ju, ierr)
    call cpu_time(finish)

    print *, "Time for pgmres = ", finish - start
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(ia)
    deallocate(ja, a)
    deallocate(w, ju)
    deallocate(jw)
    deallocate(alu, jlu, levs)
    deallocate(rhs)
    deallocate(sol)
    deallocate(vv)

end program main
