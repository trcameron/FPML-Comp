!********************************************************************************
!   LAG_TEST: Test Compensated Newton Routine vs Compensated Laguerre Routine
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 22 March 2019
!********************************************************************************
program lag_limitAcc
    use eft
    use mproutines
    implicit none
    ! parameters
    integer, parameter              :: itnum = 5
    ! polynomial variables
    integer                         :: deg
    complex(kind=dp)                :: z
    complex(kind=dp), allocatable   :: poly(:)
    ! testing variables
    real(kind=dp)                   :: a, b, g
    real(kind=dp), allocatable      :: bound(:), error(:), compbound(:), comperror(:), cond(:)
    complex(kind=dp), allocatable   :: exact(:)
    
    call mpinit
    
    ! random seed
    call init_random_seed()
    ! allocate test variables
    allocate(bound(40),error(40),compbound(40),comperror(40),cond(40),exact(40))
    ! run test
    do deg=1,40
        ! allocate polynomial variable
        allocate(poly(deg+1))
        ! polynomial
        call createpolycplx(poly, deg)
        ! exact root
        exact(deg) = createrootcplx(poly, deg)
        ! condition number and comp error bound
        g = 2*(4*deg+2)*mu*sqrt(2.0_dp)
        g = g/((1-2*mu)-g)
        cond(deg) = compcondcplx(poly,exact(deg),deg)
        bound(deg) = cond(deg)*g
        compbound(deg) = mu + cond(deg)*g**2
        if(bound(deg)>1) bound(deg) = 1
        if(compbound(deg)>1) compbound(deg) = 1
        ! standard Laguerre's method
        call random_number(a)
        call random_number(b)
        z = exact(deg) + max(compbound(deg),1E-4_dp)*cmplx(a,b,kind=dp)/sqrt(a*a+b*b)
        call SLaguerreCplx(poly,z,deg,itnum)
        error(deg) = abs(z - exact(deg))/abs(exact(deg))
        if (error(deg)<mu) then
            error(deg) = mu
        else if (error(deg)>1) then
            error(deg) = 1
        end if
        ! compensated Laguerre's method
        z = exact(deg) + max(compbound(deg),1E-4_dp)*cmplx(a,b,kind=dp)/sqrt(a*a+b*b)
        call CLaguerreCplx(poly,z,deg,itnum)
        comperror(deg) = abs(z - exact(deg))/abs(exact(deg))
        if (comperror(deg)<mu) then
            comperror(deg) = mu
        else if (comperror(deg)>1) then
            comperror(deg) = 1
        end if
        ! deallocate polynomial variable
        deallocate(poly)
    end do
    ! write test results
    open(unit=1,file="data_files/lag_limitAcc.dat")
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
    ! deallocate test variables
    deallocate(bound,error,compbound,comperror,cond,exact)
    
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
    !                   Standard Laguerre Cplx              *
    !********************************************************
    ! Computes a root approximation of a complex polynomial 
    ! using standard Laguerre method in dp-precision. The 
    ! result is stored in x. 
    !********************************************************
    subroutine SLaguerreCplx(poly,x,deg,itnum)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, itnum
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: x
        ! local variables
        integer                         :: j, k
        complex(kind=dp)                :: a, b, c
        
        ! Laguerre's method
        do k=1,itnum
            ! compute a, b, c
            a = poly(deg+1)
            b = 0
            c = 0
            do j=deg,1,-1
                c = x*c + b
                b = x*b + a
                a = x*a + poly(j)
            end do
            ! check residual
            if(abs(a)<tiny(1.0_dp)) return
            ! compute Laguerre correction
            a = b/a
            b = a**2 - 2*(c/a)
            c = sqrt((deg-1)*(deg*b-a**2))
            b = a + c
            a = a - c
            if(abs(a)>abs(b)) then
                c = deg/a
            else
                c = deg/b 
            end if
            ! update x
            x = x - c
        end do
    end subroutine SLaguerreCplx
    !********************************************************
    !              Compensated Laguerre Cplx                *
    !********************************************************
    ! Computes a root approximation of a complex polynomial 
    ! using Newton's method with compensated Horner routine.
    ! We start with initial estiamte x and update for
    ! itmax iterations. The result is stored in x. 
    !********************************************************
    subroutine CLaguerreCplx(poly,x,deg,itnum)
        ! argument variables
        integer, intent(in)             :: deg, itnum
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: x
        ! local variables
        integer                         :: i, k
        real(kind=dp)                   :: errBound, g, r
        complex(kind=dp)                :: a, b, compPoly, compDer, compDer2, errPoly, errDer, errDer2, z
        type(CEFTSum)                   :: sum
        type(CEFTProd)                  :: prod
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! allocate memory for compPoly and compDer polynomials
        do i=1,itnum
            z = x
            r = abs(z)
            compDer2 = cmplx(0,0,kind=dp)
            compDer = cmplx(0,0,kind=dp)
            errPoly = cmplx(0,0,kind=dp)
            errDer = cmplx(0,0,kind=dp)
            errDer2 = cmplx(0,0,kind=dp)
            errBound = 0.0_dp
            if(r>1) then
                z = 1/z
                r = 1/r
                ! EFT Horner applied to reversal polynomial
                compPoly = poly(1)
                do k=2,deg+1
                    ! product and sum for Der2
                    prod = TwoProductCplx(compDer2,z)
                    sum = TwoSumCplx(prod%p,compDer)
                    ! update compDer2
                    compDer2 = sum%x
                    ! update errDer2
                    errDer2 = z*errDer2 + errDer + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
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
                compDer2 = compDer2 + errDer2
                compDer = compDer + errDer
                compPoly = compPoly + errPoly
                ! compute Laguerre correction terms
                a = compDer/compPoly
                b = 2*(compDer2/compPoly)
                b = z**2*(deg-2*z*a+z**2*(a**2-b))
                a = z*(deg-z*a)
            else
                ! EFT Horner applied to polynomial
                compPoly = poly(deg+1)
                do k=deg,1,-1
                    ! product and sum for Der2
                    prod = TwoProductCplx(compDer2,z)
                    sum = TwoSumCplx(prod%p,compDer)
                    ! update compDer2
                    compDer2 = sum%x
                    ! update errDer2
                    errDer2 = z*errDer2 + errDer + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
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
                compDer2 = compDer2 + errDer2
                compDer = compDer + errDer
                compPoly = compPoly + errPoly
                ! compute Laguerre correction terms
                a = compDer/compPoly
                b = a**2 - 2*(compDer2/compPoly)
            end if
            ! errBound
            g = 2*(4*deg+2)*mu*sqrt(2.0_dp)
            g = g/((1-2*mu)-g)
            errBound = mu*abs(compPoly) + (g*errBound + 2*mu**2*abs(compPoly))
            ! update x
            if(errBound<abs(compPoly)) then
                z = sqrt((deg-1)*(deg*b-a**2))
                b = a + z
                a = a - z
                if(abs(a)>abs(b)) then
                    z = deg/a
                else
                    z = deg/b 
                end if
                if(abs(z)>mu*abs(x)) then
                    x = x - z
                else
                    return
                end if
            else
                return
            end if
        end do
    end subroutine CLaguerreCplx
end program lag_limitAcc