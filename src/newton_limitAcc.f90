!********************************************************************************
!   NEWTON_TEST: Test Compensated Newton Routine implemented as it would be in root solver
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 22 March 2019
!********************************************************************************
program newton_limitAcc
    use eft
    use mproutines
    implicit none
    ! parameters
    integer, parameter              :: itnum = 5
    ! polynomial variables
    integer                         :: deg
    real(kind=dp)                   :: x
    complex(kind=dp)                :: z
    real(kind=dp), allocatable      :: poly(:)
    complex(kind=dp), allocatable   :: cpoly(:)
    ! testing variables
    real(kind=dp)                   :: a, b, g
    real(kind=dp), allocatable      :: bound(:), error(:), compbound(:), comperror(:), cond(:), exact(:)
    complex(kind=dp), allocatable   :: cexact(:)
    
    call mpinit
    
    ! random seed
    call init_random_seed()
    ! allocate real test variables
    allocate(bound(40),error(40),compbound(40),comperror(40),cond(40),exact(40))
    ! run real test
    do deg=1,40
        ! allocate polynomial variable
        allocate(poly(deg+1))
        ! polynomial
        call createpoly(poly, deg)
        ! exact root
        exact(deg) = createroot(poly,deg)
        ! condition number and error bounds
        g = (2*deg*mu/(1-2*deg*mu))
        cond(deg) = compcond(poly,exact(deg),deg)
        bound(deg) = cond(deg)*g
        compbound(deg) = mu + cond(deg)*g**2
        if(bound(deg)>1) bound(deg) = 1
        if(compbound(deg)>1) compbound(deg) = 1
        ! standard Newton's method
        call random_number(a)
        x = exact(deg) + max(compbound(deg),1E-4_dp)*a/abs(a)
        call SNewton(poly,x,deg,itnum)
        error(deg) = abs(x - exact(deg))/abs(exact(deg))
        if (error(deg)<mu) then
            error(deg) = mu
        else if (error(deg)>1) then
            error(deg) = 1
        end if
        ! compensated Newton's method
        x = exact(deg) + max(compbound(deg),1E-4_dp)*a/abs(a)
        call CNewton(poly,x,deg,itnum)
        comperror(deg) = abs(x - exact(deg))/abs(exact(deg))
        if (comperror(deg)<mu) then
            comperror(deg) = mu
        else if (comperror(deg)>1) then
            comperror(deg) = 1
        end if
        ! deallocate polynomial variable
        deallocate(poly)
    end do
    ! write real test results
    open(unit=1,file="data_files/newton_limitAcc_real.dat")
    write(1,'(A)') 'cond, stan_err_bound, stan_err, comp_err_bound, comp_err'
    do deg=1,40
        write(1,'(ES15.2)', advance='no') cond(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') bound(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') error(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') compbound(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') comperror(deg)
    end do
    ! close file
    close(1)
    ! deallocate real test variables
    deallocate(bound,error,compbound,comperror,cond,exact)
    
    ! allocate complex test variables
    allocate(bound(40),error(40),compbound(40),comperror(40),cond(40),cexact(40))
    ! run complex test
    do deg=1,40
        ! allocate polynomial variable
        allocate(cpoly(deg+1))
        ! polynomial
        call createpolycplx(cpoly, deg)
        ! exact root
        cexact(deg) = createrootcplx(cpoly, deg)
        ! condition number and error bounds
        g = 2*(4*deg+2)*mu*sqrt(2.0_dp)
        g = g/((1-2*mu)-g)
        cond(deg) = compcondcplx(cpoly,cexact(deg),deg)
        bound(deg) = cond(deg)*g
        compbound(deg) = mu + cond(deg)*g**2
        if(bound(deg)>1) bound(deg) = 1
        if(compbound(deg)>1) compbound(deg) = 1
        ! standard Newton's method
        call random_number(a)
        call random_number(b)
        z = cexact(deg) + max(compbound(deg),1E-4_dp)*cmplx(a,b,kind=dp)/sqrt(a*a+b*b)
        call SNewtonCplx(cpoly,z,deg,itnum)
        error(deg) = abs(z - cexact(deg))/abs(cexact(deg))
        if (error(deg)<mu) then
            error(deg) = mu
        else if (error(deg)>1) then
            error(deg) = 1
        end if
        ! compensated Newton's method
        z = cexact(deg) + max(compbound(deg),1E-4_dp)*cmplx(a,b,kind=dp)/sqrt(a*a+b*b)
        call CNewtonCplx(cpoly,z,deg,itnum)
        comperror(deg) = abs(z - cexact(deg))/abs(cexact(deg))
        if (comperror(deg)<mu) then
            comperror(deg) = mu
        else if (comperror(deg)>1) then
            comperror(deg) = 1
        end if
        ! deallocate polynomial variable
        deallocate(cpoly)
    end do
    ! write complex test results
    open(unit=1,file="data_files/newton_limitAcc_complex.dat")
    write(1,'(A)') 'cond, stan_err_bound, stan_err, comp_err_bound, comp_err'
    do deg=1,40
        write(1,'(ES15.2)', advance='no') cond(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') bound(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') error(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') compbound(deg)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') comperror(deg)
    end do
    ! close file
    close(1)
    ! deallocate complex test variables
    deallocate(bound,error,compbound,comperror,cond,cexact)
    
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
    !********************************************************
    !                   Standard Newton                     *
    !********************************************************
    ! Computes a root approximation of a polynomial using
    ! standard Newton's method in dp-precision. The result
    ! is stored in x. 
    !********************************************************
    subroutine SNewton(poly,x,deg,itnum)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, itnum
        real(kind=dp), intent(in)       :: poly(:)
        real(kind=dp), intent(inout)    :: x
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: alpha, beta
        
        ! Newton's method
        do k=1,itnum
            ! compute alpha and beta
            alpha = poly(deg+1)
            beta = 0
            do j=deg,1,-1
                beta = x*beta + alpha
                alpha = x*alpha + poly(j)
            end do
            ! check residual
            if(abs(alpha)<tiny(1.0_dp)) return
            ! update x
            x = x - (alpha/beta)
        end do
    end subroutine SNewton
    !********************************************************
    !                   Standard Newton Cplx                *
    !********************************************************
    ! Computes a root approximation of a complex polynomial 
    ! using standard Newton's method in dp-precision. The 
    ! result is stored in x. 
    !********************************************************
    subroutine SNewtonCplx(poly,x,deg,itnum)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, itnum
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: x
        ! local variables
        integer                         :: j, k
        complex(kind=dp)                :: alpha, beta
        
        ! Newton's method
        do k=1,itnum
            ! compute alpha and beta
            alpha = poly(deg+1)
            beta = 0
            do j=deg,1,-1
                beta = x*beta + alpha
                alpha = x*alpha + poly(j)
            end do
            ! update x
            x = x - (alpha/beta)
        end do
    end subroutine SNewtonCplx
    !********************************************************
    !              Compensated Newton                       *
    !********************************************************
    ! Computes a root approximation of a polynomial using
    ! Newton's method with compensated Horner routine.
    ! We start with initial estiamte x and update for
    ! itmax iterations. The result is stored in x. 
    !********************************************************
    subroutine CNewton(poly,x,deg,itnum)
        ! argument variables
        integer, intent(in)             :: deg, itnum
        real(kind=dp), intent(in)       :: poly(:)
        real(kind=dp), intent(inout)    :: x
        ! local variables
        integer                         :: i, k
        real(kind=dp)                   :: errBound, g, r
        real(kind=dp)                   :: a, compPoly, compDer, errPoly, errDer, z
        type(REFT)                      :: prod, sum
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! allocate memory for compPoly and compDer polynomials
        do i=1,itnum
            z = x
            r = abs(z)
            compDer = 0.0_dp
            errPoly = 0.0_dp
            errDer = 0.0_dp
            errBound = 0.0_dp
            if(r>1) then
                z = 1/z
                r = 1/r
                ! EFT Horner applied to reversal polynomial
                compPoly = poly(1)
                do k=2,deg+1
                    ! product and sum for Der
                    prod = TwoProduct(compDer,z)
                    sum = TwoSum(prod%x,compPoly)
                    ! update compDer
                    compDer = sum%x
                    ! update errDer
                    errDer = z*errDer + errPoly + (prod%y + sum%y)
                    ! product and sum for Poly
                    prod = TwoProduct(compPoly,z)
                    sum = TwoSum(prod%x,poly(k))
                    ! update CompPoly
                    compPoly = sum%x
                    ! update errPoly
                    errPoly = z*errPoly + (prod%y + sum%y)
                    ! update errBound
                    errBound = r*errBound + (abs(prod%y) + abs(sum%y))
                end do
                ! add error back into result
                compDer = compDer + errDer
                compPoly = compPoly + errPoly
                ! compute Newton correction term
                a = 1/(z*(deg-z*compDer/compPoly))
            else
                ! EFT Horner applied to polynomial
                compPoly = poly(deg+1)
                do k=deg,1,-1
                    ! product and sum for Der
                    prod = TwoProduct(compDer,z)
                    sum = TwoSum(prod%x,compPoly)
                    ! update compDer
                    compDer = sum%x
                    ! update errDer
                    errDer = z*errDer + errPoly + (prod%y + sum%y)
                    ! product and sum for Poly
                    prod = TwoProduct(compPoly,z)
                    sum = TwoSum(prod%x,poly(k))
                    ! update CompPoly
                    compPoly = sum%x
                    ! update errPoly
                    errPoly = z*errPoly + (prod%y + sum%y)
                    ! update errBound
                    errBound = r*errBound + (abs(prod%y) + abs(sum%y))
                end do
                ! add error back into result
                compDer = compDer + errDer
                compPoly = compPoly + errPoly
                ! compute Newton correction term
                a = compPoly/compDer
            end if
            ! errBound
            g = (4*deg+2)*mu/(1-(4*deg+2)*mu)
            errBound = mu*abs(compPoly) + (g*errBound + 2*mu**2*abs(compPoly))
            ! update x
            if(errBound<abs(compPoly) .and. abs(a)>mu*abs(x)) then
                x = x - a
            else
                return
            end if
        end do
    end subroutine CNewton
    !********************************************************
    !              Compensated Newton Cplx                  *
    !********************************************************
    ! Computes a root approximation of a complex polynomial 
    ! using Newton's method with compensated Horner routine.
    ! We start with initial estiamte x and update for
    ! itmax iterations. The result is stored in x. 
    !********************************************************
    subroutine CNewtonCplx(poly,x,deg,itnum)
        ! argument variables
        integer, intent(in)             :: deg, itnum
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: x
        ! local variables
        integer                         :: i, k
        real(kind=dp)                   :: errBound, g, r
        complex(kind=dp)                :: a, compPoly, compDer, errPoly, errDer, z
        type(CEFTSum)                   :: sum
        type(CEFTProd)                  :: prod
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! allocate memory for compPoly and compDer polynomials
        do i=1,itnum
            z = x
            r = abs(z)
            compDer = cmplx(0,0,kind=dp)
            errPoly = cmplx(0,0,kind=dp)
            errDer = cmplx(0,0,kind=dp)
            errBound = 0.0_dp
            if(r>1) then
                z = 1/z
                r = 1/r
                ! EFT Horner applied to reversal polynomial
                compPoly = poly(1)
                do k=2,deg+1
                    ! product and sum for Der
                    prod = TwoProductCplx(compDer,z)
                    sum = TwoSumCplx(prod%p,compPoly)
                    ! update compDer
                    compDer = sum%x
                    ! update errDer
                    errDer = z*errDer + errPoly + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
                    ! product and sum for Poly
                    prod = TwoProductCplx(compPoly,z)
                    sum = TwoSumCplx(prod%p,poly(k))
                    ! update CompPoly
                    compPoly = sum%x
                    ! update errPoly
                    errPoly = z*errPoly + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
                    ! update errBound
                    errBound = r*errBound + FaithSum(abs(prod%e),abs(prod%f),abs(prod%g),abs(sum%y))
                end do
                ! add error back into result
                compDer = compDer + errDer
                compPoly = compPoly + errPoly
                ! compute Newton correction term
                a = 1/(z*(deg-z*compDer/compPoly))
            else
                ! EFT Horner applied to polynomial
                compPoly = poly(deg+1)
                do k=deg,1,-1
                    ! product and sum for Der
                    prod = TwoProductCplx(compDer,z)
                    sum = TwoSumCplx(prod%p,compPoly)
                    ! update compDer
                    compDer = sum%x
                    ! update errDer
                    errDer = z*errDer + errPoly + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
                    ! product and sum for Poly
                    prod = TwoProductCplx(compPoly,z)
                    sum = TwoSumCplx(prod%p,poly(k))
                    ! update compPoly
                    compPoly = sum%x
                    ! update errPoly
                    errPoly = z*errPoly + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
                    ! update errBound
                    errBound = r*errBound + FaithSum(abs(prod%e),abs(prod%f),abs(prod%g),abs(sum%y))
                end do
                ! add error back into result
                compDer = compDer + errDer
                compPoly = compPoly + errPoly
                ! compute Newton correction term
                a = compPoly/compDer
            end if
            ! errBound
            g = 2*(4*deg+2)*mu*sqrt(2.0_dp)
            g = g/((1-2*mu)-g)
            errBound = mu*abs(compPoly) + (g*errBound + 2*mu**2*abs(compPoly))
            ! update x
            if(errBound<abs(compPoly) .and. abs(a)>mu*abs(x)) then
                x = x - a
            else
                return
            end if
        end do
    end subroutine CNewtonCplx
end program newton_limitAcc