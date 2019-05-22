!********************************************************************************
!   BACKERR_WELLCOND: Compare the backward error of FPML and FPML-CP.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 28 April 2019
!********************************************************************************
!   
!********************************************************************************
program backerr_wellcond
    use fpml_comp
    use mproutines
    implicit none
    ! testing variables
    integer                                     :: deg, it, j, startDegree, endDegree, itnum, clock_rate, clock_start, clock_stop
    real(kind=dp), allocatable                  :: results(:,:)
    ! FPML variables
    integer, parameter                          :: nitmax = 30
    logical, allocatable                        :: conv(:)
    real(kind=dp), allocatable                  :: berr(:), cond(:)   
    complex(kind=dp), allocatable               :: p(:), roots(:)
    ! external function
    double precision                            :: dznrm2
    external                                    :: dznrm2
    
    call mpinit
    
    ! Test: random coefficients
    startDegree=50
    endDegree=550
    itnum=128
    call init_random_seed()
    open(unit=1,file="data_files/backerr_wellcond_rand.dat")
    write(1,'(A)') 'Degree, FPML-time, FPML-err, FPML-CP-time, FPML-CP-err'
    allocate(results(itnum,4))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ','
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
        do it=1,itnum
            ! complex polynomial with random coefficients
            call cmplx_rand(deg+1,p)
            ! FPML, No Polish
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond, conv, nitmax)
            call system_clock(count=clock_stop)
            results(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2) = error(p, roots, deg)
            ! FPML, Compensated Laguerre Polish
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call MLPolish(p, deg, roots)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,4) = error(p, roots, deg)
        end do
        deallocate(p, roots, berr, cond, conv)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,3))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,4))/itnum
        ! update deg
        deg = deg + 10
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
    
    ! Test: natural number coefficients
    startDegree=50
    endDegree=550
    itnum=128
    open(unit=1,file="data_files/backerr_wellcond_nat.dat")
    write(1,'(A)') 'Degree, FPML-time, FPML-err, FPML-CP-time, FPML-CP-err'
    allocate(results(itnum,6))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
        do it=1,itnum
            ! polynomial with natural number coefficients
            p = (/ (j, j=1,deg+1)/)
            ! FPML, No Polish
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call main(p, deg, roots, berr, cond, conv, nitmax)
            call system_clock(count=clock_stop)
            results(it,1)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,2) = error(p, roots, deg)
            ! FPML, Compensated Laguerre Polish
            call system_clock(count_rate=clock_rate)
            call system_clock(count=clock_start)
            call MLPolish(p, deg, roots)
            call system_clock(count=clock_stop)
            results(it,3)=dble(clock_stop-clock_start)/dble(clock_rate)
            results(it,4) = error(p, roots, deg)
        end do
        deallocate(p, roots, berr, cond, conv)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,3))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,4))/itnum
        ! update deg
        deg = deg + 10
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
    
contains
    !************************************************
    !                   init_random_seed            *
    !************************************************
    ! Initiate random seed using system_clock. This
    ! seed is then available for the random number
    ! generator in random_number for the life of
    ! the program.
    !************************************************
    subroutine init_random_seed()
        implicit none
        ! local variables
        integer                             :: i, n, clock
        integer, dimension(:), allocatable  :: seed
        ! intrinsic subroutines
        intrinsic                           :: random_seed, system_clock
        
        ! main
        call random_seed(size = n)
        allocate(seed(n))
        
        call system_clock(count = clock)
        seed = clock + 37 * (/ (i - 1, i = 1,n) /)
        call random_seed(put = seed)
        
        deallocate(seed)
    end subroutine init_random_seed
    !************************************************
    !                   cmplx_rand                  *
    !************************************************
    ! Creates array of random complex numbers of
    ! size n. Each complex number has real and
    ! imaginary parts uniformly distributed on (-1,1).
    !************************************************
    subroutine cmplx_rand(n,x)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(out)   :: x(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: r1, r2
        
        ! main
        do k=1,n
            call random_number(r1)
            call random_number(r2)
            x(k) = cmplx(-1+2*r1,-1+2*r2,kind=dp)
        end do
    end subroutine cmplx_rand
    !************************************************
    !                   error                       *
    !************************************************
    ! Computes the backward error of the polynomial
    ! with respect to the 1-norm.
    !************************************************
    function error(p, roots, deg) result(err)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: p(:), roots(:)
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: err, rho
        complex(kind=dp)                :: alpha(deg+1)
        
        alpha(deg+1) = p(deg+1)
        call RootCoeffCplx2(deg, roots, alpha)
        rho = maxval((/ (abs(p(k)), k=1,deg+1)/))
        err = maxval((/ (abs(alpha(k)-p(k))/rho, k=1,deg+1)/))
        return
    end function error
end program backerr_wellcond