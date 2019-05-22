!********************************************************************************
!   FORWERR_SIMPLE: Compare the forward error of FPML, FPML-CP, and C02AFF.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 28 April 2019
!********************************************************************************
!   
!********************************************************************************
program forwerr_simple
    use fpml_comp
    use mproutines
    implicit none
    ! testing variables
    integer                                     :: deg, it, j, startDegree, endDegree, itnum
    real(kind=dp), allocatable                  :: results(:,:)
    complex(kind=dp), allocatable               :: exact_roots(:)
    ! FPML variables
    integer, parameter                          :: nitmax = 30
    logical, allocatable                        :: conv(:)
    real(kind=dp), allocatable                  :: berr(:), cond(:)   
    complex(kind=dp), allocatable               :: p(:), roots(:)
    ! NAG variables
    logical, parameter                          :: scal = .false.
    integer                                     :: ifail
    real(kind=dp), allocatable                  :: a(:,:), w(:), z(:,:)
    complex(kind=dp), allocatable               :: zeros(:)
    
    call mpinit
    
    ! Test1: Chebyshev Polynomial
    startDegree = 5
    endDegree = 80
    itnum = 10
    open(unit=1,file="data_files/forwerr_cheby.dat")
    write(1,'(A)') 'Degree, FPML_err, FPML-CP_err, C02AFF_err'
    allocate(results(itnum,3))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg), zeros(deg))
        ! polynomial
        call make_poly(1, deg, exact_roots, p)
        do j=0,deg
            a(1,j) = real(p(deg+1-j))
            a(2,j) = aimag(p(deg+1-j))
        end do
        do it=1,itnum
            ! FPML, No Polish
            call main(p, deg, roots, berr, cond, conv, nitmax)
            results(it,1) = max_rel_err(roots, exact_roots, deg)
            ! FPML, Compensated Laguerre Polish
            call MLPolish(p, deg, roots)
            results(it,2) = max_rel_err(roots, exact_roots, deg)
            ! C02AFF
            ifail = 0
            call c02aff(a,deg,scal,z,w,ifail)
            zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
            results(it,3) = max_rel_err(zeros, exact_roots, deg)
        end do
        deallocate(roots, berr, cond, conv, a, w, z, zeros)
        deallocate(exact_roots, p)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,3))/itnum
        ! update deg
        deg = deg + 1
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
    
    ! Test2: Prescribed Roots of Varying Scale - 3
    startDegree = 5
    endDegree = 20
    itnum = 10
    open(unit=1,file="data_files/forwerr_presc.dat")
    write(1,'(A)') 'Degree, FPML_err, FPML-CP_err, C02AFF_err'
    allocate(results(itnum,3))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg), zeros(deg))
        ! polynomial
        call make_poly(2, deg, exact_roots, p)
        do j=0,deg
            a(1,j) = real(p(deg+1-j))
            a(2,j) = aimag(p(deg+1-j))
        end do
        do it=1,itnum
            ! FPML, No Polish
            call main(p, deg, roots, berr, cond, conv, nitmax)
            results(it,1) = max_rel_err(roots, exact_roots, deg)
            ! FPML, Compensated Laguerre Polish
            call MLPolish(p, deg, roots)
            results(it,2) = max_rel_err(roots, exact_roots, deg)
            ! C02AFF
            ifail = 0
            call c02aff(a,deg,scal,z,w,ifail)
            zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
            results(it,3) = max_rel_err(zeros, exact_roots, deg)
        end do
        deallocate(roots, berr, cond, conv, a, w, z, zeros)
        deallocate(exact_roots, p)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,3))/itnum
        ! update deg
        deg = deg + 1
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
    
    ! Test3: Small Imaginary Part
    startDegree = 5
    endDegree = 25
    itnum = 10
    open(unit=1,file="data_files/forwerr_smallimag.dat")
    write(1,'(A)') 'Degree, FPML_err, FPML-CP_err, C02AFF_err'
    allocate(results(itnum,3))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg), zeros(deg))
        ! polynomial
        call make_poly(3, deg, exact_roots, p)
        do j=0,deg
            a(1,j) = real(p(deg+1-j))
            a(2,j) = aimag(p(deg+1-j))
        end do
        do it=1,itnum
            ! FPML, No Polish
            call main(p, deg, roots, berr, cond, conv, nitmax)
            results(it,1) = max_rel_err(roots, exact_roots, deg)
            ! FPML, Compensated Laguerre Polish
            call MLPolish(p, deg, roots)
            results(it,2) = max_rel_err(roots, exact_roots, deg)
            ! C02AFF
            ifail = 0
            call c02aff(a,deg,scal,z,w,ifail)
            zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
            results(it,3) = max_rel_err(zeros, exact_roots, deg)
        end do
        deallocate(roots, berr, cond, conv, a, w, z, zeros)
        deallocate(exact_roots, p)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,3))/itnum
        ! update deg
        deg = deg + 1
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
    
    ! Test4: Wilkinson Polynomial
    startDegree = 5
    endDegree = 20
    itnum = 10
    open(unit=1,file="data_files/forwerr_wilk.dat")
    write(1,'(A)') 'Degree, FPML_err, FPML-CP_err, C02AFF_err'
    allocate(results(itnum,3))
    deg = startDegree
    do while(deg<=endDegree)
        write(1,'(I10,A)', advance='no') deg, ', '
        allocate(roots(deg), berr(deg), cond(deg), conv(deg))
        allocate(a(2,0:deg),w(4*(deg+1)), z(2,deg), zeros(deg))
        ! polynomial
        call make_poly(4, deg, exact_roots, p)
        do j=0,deg
            a(1,j) = real(p(deg+1-j))
            a(2,j) = aimag(p(deg+1-j))
        end do
        do it=1,itnum
            ! FPML, No Polish
            call main(p, deg, roots, berr, cond, conv, nitmax)
            results(it,1) = max_rel_err(roots, exact_roots, deg)
            ! FPML, Compensated Laguerre Polish
            call MLPolish(p, deg, roots)
            results(it,2) = max_rel_err(roots, exact_roots, deg)
            ! C02AFF
            ifail = 0
            call c02aff(a,deg,scal,z,w,ifail)
            zeros = (/ (cmplx(z(1,j),z(2,j),kind=dp), j=1,deg)/)
            results(it,3) = max_rel_err(zeros, exact_roots, deg)
        end do
        deallocate(roots, berr, cond, conv, a, w, z, zeros)
        deallocate(exact_roots, p)
        ! write results to file
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,1))/itnum, ','
        write(1,'(ES15.2,A)', advance='no') sum(results(1:itnum,2))/itnum, ','
        write(1,'(ES15.2)') sum(results(1:itnum,3))/itnum
        ! update deg
        deg = deg + 1
    end do
    ! deallocate results and close file
    deallocate(results)
    close(1)
contains
    !************************************************
    !                       make_poly               *
    !************************************************
    ! Make special polynomial corresponding to 
    ! poly_num and deg > 1.
    !************************************************
    subroutine make_poly(poly_num, deg, exact_roots, p)
        implicit none
        ! argument variables
        integer, intent(in)                         :: poly_num, deg
        complex(kind=dp), allocatable, intent(out)  :: exact_roots(:), p(:)
        ! local variables
        integer                                     :: j
        real(kind=dp)                               :: tmp, a(deg+1), b(deg+1), c(deg+1)
        real(kind=dp), parameter                    :: pi = 3.1415926535897932385E0_dp
        
        select case (poly_num)
            case(1)
            ! Chebyshev Polynomial
            allocate(exact_roots(deg), p(deg+1))
            do j=1,deg
                tmp = cos(pi*(2.0_dp*j-1.0_dp)/(2.0_dp*deg))
                if(abs(tmp)<mu) tmp = 0.0_dp
                exact_roots(j) = cmplx(tmp,0,kind=dp)
            end do
            a = 0.0_dp; a(1) = 1.0_dp
            b = 0.0_dp; b(2) = 1.0_dp
            do j=2,deg
                c = 0.0_dp
                c(2:deg+1) = b(1:deg)
                p = cmplx(2.0_dp*c - a,0,kind=dp)
                a = b
                b = real(p,kind=dp)
            end do
            case(2)
            ! Prescribed Roots of Varying scale
            allocate(exact_roots(deg), p(deg+1))
            if(mod(deg,2)==0) then
                exact_roots = (/ (cmplx(2.0_dp**(-deg/2+j)-3.0_dp,0,kind=dp), j=0,deg-1)/)
            else
                exact_roots = (/ (cmplx(2.0_dp**(-(deg-1)/2+j)-3.0_dp,0,kind=dp), j=0,deg-1)/)
            end if
            call RootCoeffCplx(deg,exact_roots,p)
            case(3)
            ! Small Imaginary Part
            allocate(exact_roots(deg), p(deg+1))
            do j=1,deg
                if(mod(j,2)==0) then
                    exact_roots(j) = cmplx(j,-10*mu,kind=dp)
                else
                    exact_roots(j) = cmplx(j,10*mu,kind=dp)
                end if
            end do
            call RootCoeffCplx(deg,exact_roots,p)
            case(4)
            ! Wilkinson Polynomial
            allocate(exact_roots(deg), p(deg+1))
            exact_roots = (/ (cmplx(j,0,kind=dp), j=1,deg)/)
            call RootCoeffCplx(deg,exact_roots,p)
        end select
    end subroutine make_poly
    !************************************************
    !                       max_rel_err             *
    !************************************************
    ! Compute maximum relative error between roots 
    ! and exact roots.
    !************************************************
    function max_rel_err(roots, exact_roots, deg) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: exact_roots(:), roots(:)
        ! local variables
        logical                         :: m(deg)
        integer                         :: i, j, k
        real(kind=dp)                   :: delta, err, res, relerr, rmax
        ! intrinsic procedures
        intrinsic                       :: abs, max, huge
        
        ! initialize
        m = .True.
        rmax = huge(1.0_dp)
        res = mu
        ! main loop
        do i=1,deg
            ! find corresponding exact root (yet to be matched)
            err = rmax
            do j=1,deg
                if(m(j)) then
                    delta = abs(exact_roots(j) - roots(i))
                    if(delta < err) then
                        err = delta
                        k = j
                    end if
                end if
            end do
            ! mark corresponding root as matched
            m(k) = .False.
            ! calculate relative error on this root and update max
            ! zero roots give unhelpful relative errors, so revert to absolute error
            if(abs(roots(i))<mu .or. abs(exact_roots(k))<mu) then
                relerr = err
            else
                relerr = err/abs(exact_roots(k))
            end if
            res = max(relerr, res)
        end do
        res = max(res,mu)
    end function max_rel_err
end program forwerr_simple