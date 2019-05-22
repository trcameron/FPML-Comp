!********************************************************************************
!   BACKERR_ILLCOND: Compare the backward error of FPML, FPML-P, and FPML-CP.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 28 April 2019
!********************************************************************************
!   
!********************************************************************************
program backerr_illcond
    use fpml_comp
    use mproutines
    implicit none
    ! testing variables
    integer                                     :: deg, it, j, startDegree, endDegree, itnum
    real(kind=dp), allocatable                  :: results(:,:)
    ! FPML variables
    integer, parameter                          :: nitmax = 30
    logical, allocatable                        :: conv(:)
    real(kind=dp), allocatable                  :: berr(:), cond(:)   
    complex(kind=dp), allocatable               :: p(:), roots(:)
    
    call mpinit
    
    ! Test: random roots in unit disk
    startDegree=10
    endDegree=150
    itnum=128
    call init_random_seed()
    open(unit=1,file="data_files/backerr_illcond_unit.dat")
    write(1,'(A)') 'Degree, FPML_err, FPML-P_err, FPML-CP_err'
    allocate(results(itnum,3))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
        do it=1,itnum
            ! polynomial with random roots in unit disk
            call cmplx_rand_unit_roots(deg,roots)
            call RootCoeffCplx(deg,roots,p)
            ! FPML, No Polish
            call main(p, deg, roots, berr, cond, conv, nitmax)
            results(it,1) = error(p, roots, deg)
            ! FPML, Simple Polish
            call NPolish(p, deg, roots)
            results(it,2) = error(p, roots, deg)
            ! FPML, Compensated Laguerre Polish
            call MLPolish(p, deg, roots)
            results(it,3) = error(p, roots, deg)
        end do
        deallocate(p, roots, berr, cond, conv)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,3))/itnum
        ! update deg
        deg = deg + 10
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
    
    ! Test: truncated exponential
    startDegree=10
    endDegree=100
    itnum=128
    open(unit=1,file="data_files/backerr_illcond_trunc.dat")
    write(1,'(A)') 'Degree, FPML_err, FPML-P_err, FPML-CP_err'
    allocate(results(itnum,3))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
        do it=1,itnum
            ! truncated exponential
            p(1) = 1.0_dp
            p(2) = 1.0_dp
            do j=3,deg+1
                p(j) = p(j-1)/(j-1)
            end do
            ! FPML, No Polish
            call main(p, deg, roots, berr, cond, conv, nitmax)
            results(it,1) = error(p, roots, deg)
            ! FPML, Simple Polish
            call NPolish(p, deg, roots)
            results(it,2) = error(p, roots, deg)
            ! FPML, Compensated Laguerre Polish
            call MLPolish(p, deg, roots)
            results(it,3) = error(p, roots, deg)
        end do
        deallocate(p, roots, berr, cond, conv)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,3))/itnum
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
    !                   cmplx_rand_unit_roots       *
    !************************************************
    ! Creates array of random complex numbers of
    ! size n that lies in the unit disk.
    !************************************************
    subroutine cmplx_rand_unit_roots(n,x)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(out)   :: x(:)
        ! local variables
        integer                         :: k
        real(kind=dp), parameter        :: pi2 = 6.2831853071795865_dp
        real(kind=dp)                   :: r, t
        
        ! main
        do k=1,n
            call random_number(r)
            call random_number(t)
            x(k) = r*cmplx(cos(t*pi2),sin(t*pi2),kind=dp)
        end do
    end subroutine cmplx_rand_unit_roots
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
end program backerr_illcond