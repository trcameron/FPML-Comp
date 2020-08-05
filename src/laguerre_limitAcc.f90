!********************************************************************************
!   Laguerre_LIMITACC: Test limiting accuracy of compensated Laguerre method
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 22 March 2019
!********************************************************************************
program laguerre_limitAcc
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
    real(kind=dp)                   :: bound, error, compBound, compError, cond, exact, g
    complex(kind=dp)                :: cexact
    
    ! random seed
    call init_random_seed()
    ! open real test file
    open(unit=1,file="data_files/laguerre_limitAcc_real.dat")
    write(1,'(A)') 'cond, stan_err_bound, stan_err, comp_err_bound, comp_err'
    ! run real test
    do deg=1,40
        ! allocate polynomial variable
        allocate(poly(deg+1))
        ! polynomial
        call limAccPoly(poly,deg)
        ! exact root
        exact = limAccRoot(poly,deg)
        ! condition number and error bounds
        g = 2*deg*mu/(1-2*deg*mu)
        cond = compCond(poly,exact,deg)
        bound = cond*g
        compBound = mu + cond*g**2
        if(bound>1.0_dp) bound = 1.0_dp
        if(compBound>1.0_dp) compBound = 1.0_dp
        ! standard Laguerre method
        x = exact + exact*compBound
        call SLaguerre(poly,x,deg)
        error = abs(x - exact)/abs(exact)
        if(error<mu) then
            error = mu
        else if(error>1.0_dp) then
            error = 1.0_dp
        end if    
        ! compensated Laguerre method
        x = exact + exact*compBound
        call CLaguerre(poly,x,deg)
        compError = abs(x - exact)/abs(exact)
        if(compError<mu) then
            compError = mu
        else if(compError>1.0_dp) then
            compError = 1.0_dp
        end if
        ! deallocate polynomial
        deallocate(poly)
        ! write to file
        write(1,'(ES15.2)', advance='no') cond
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') bound
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') error
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') compBound
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') compError
    end do
    ! close real test file
    close(1)
    
    ! open complex test file
    open(unit=1,file="data_files/laguerre_limitAcc_cmplx.dat")
    write(1,'(A)') 'cond, stan_err_bound, stan_err, comp_err_bound, comp_err'
    ! run complex text
    do deg=1,40
        ! allocate polynomial variable
        allocate(cpoly(deg+1))
        ! polynomial
        call limAccPolyCplx(cpoly,deg)
        ! exact root
        cexact = limAccRootCplx(cpoly,deg)
        ! condition number and error bounds
        g = 2*mu/(1-2*mu)
        g = 2*deg*sqrt(2.0_dp)*g/(1-2*deg*sqrt(2.0_dp)*g)
        cond = compCondCplx(cpoly,cexact,deg)
        bound = cond*g
        compBound = mu + cond*g**2
        if(bound>1.0_dp) bound = 1.0_dp
        if(compBound>1.0_dp) compBound = 1.0_dp
        ! standard Laguerre method
        z = cexact + cexact*compBound
        call SLaguerreCplx(cpoly,z,deg)
        error = abs(z - cexact)/abs(cexact)
        if(error<mu) then
            error = mu
        else if(error>1.0_dp) then
            error = 1.0_dp
        end if
        ! compensated Laguerre method
        z = cexact + cexact*compBound
        call CLaguerreCplx(cpoly,z,deg)
        compError = abs(z - cexact)/abs(cexact)
        if(compError<mu) then
            compError = mu
        else if(compError>1.0_dp) then
            compError = 1.0_dp
        end if
        ! deallocate polynomial
        deallocate(cpoly)
        ! write to file
        write(1,'(ES15.2)', advance='no') cond
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') bound
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') error
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') compBound
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') compError
    end do
    ! close complex test file
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
    !********************************************************
    !                   Standard Laguerre                   *
    !********************************************************
    ! Computes a root approximation of a polynomial using
    ! standard Laguerre method in dp-precision. The result
    ! is stored in x. 
    !********************************************************
    subroutine SLaguerre(poly,x,deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(in)       :: poly(:)
        real(kind=dp), intent(inout)    :: x
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: alpha, beta, gamma
        
        ! Laguerre's method
        do k=1,itnum
            ! compute alpha and beta
            alpha = poly(deg+1)
            beta = 0.0_dp
            gamma = 0.0_dp
            do j=deg,1,-1
                gamma = x*gamma + beta
                beta = x*beta + alpha
                alpha = x*alpha + poly(j)
            end do
            ! check residual
            if(abs(alpha)<small) return
            ! compute Laguerre correction terms
            alpha = beta/alpha
            beta = alpha**2 - 2*(gamma/alpha)
            gamma = sqrt((deg-1)*(deg*beta-alpha**2))
            beta = alpha + gamma
            alpha = alpha - gamma
            ! update x
            if(abs(alpha)>abs(beta)) then
                x = x - (deg/alpha)
            else
                x = x - (deg/beta)
            end if
        end do
    end subroutine SLaguerre
    !********************************************************
    !                   Standard Laguerre Cplx              *
    !********************************************************
    ! Complex version of Standard Laguerre.
    !********************************************************
    subroutine SLaguerreCplx(poly,x,deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: x
        ! local variables
        integer                         :: j, k
        complex(kind=dp)                :: alpha, beta, gamma
        
        ! Laguerre's method
        do k=1,itnum
            ! compute alpha and beta
            alpha = poly(deg+1)
            beta = cmplx(0.0_dp,0.0_dp,kind=dp)
            gamma = cmplx(0.0_dp,0.0_dp,kind=dp)
            do j=deg,1,-1
                gamma = x*gamma + beta
                beta = x*beta + alpha
                alpha = x*alpha + poly(j)
            end do
            ! check residual
            if(abs(alpha)<small) return
            ! compute Laguerre correction terms
            alpha = beta/alpha
            beta = alpha**2 - 2*(gamma/alpha)
            gamma = sqrt((deg-1)*(deg*beta-alpha**2))
            beta = alpha + gamma
            alpha = alpha - gamma
            ! update x
            if(abs(alpha)>abs(beta)) then
                x = x - (deg/alpha)
            else
                x = x - (deg/beta)
            end if
        end do
    end subroutine SLaguerreCplx
    !********************************************************
    !                   Compensated Laguerre                *
    !********************************************************
    ! Computes a root approximation of a polynomial using
    ! compensated Laguerre method in dp-precision, i.e.,
    ! the polynomial evaluation is done using a compensated
    ! Horner's method. The result is stored in x. 
    !********************************************************
    subroutine CLaguerre(poly,x,deg)
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(in)       :: poly(:)
        real(kind=dp), intent(inout)    :: x
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: alpha, alphaErr, beta, betaErr, gamma, gammaErr
        type(REFT)                      :: prod, sum
        
        ! Laguerre method
        do k=1,itnum
            ! compute alpha and beta
            alpha = poly(deg+1)
            alphaErr = 0.0_dp
            beta = 0.0_dp
            betaErr = 0.0_dp
            gamma = 0.0_dp
            gammaErr = 0.0_dp
            do j=deg,1,-1
                ! product and sum for gamma
                prod = TwoProduct(gamma,x)
                sum = TwoSum(prod%x,beta)
                ! update gamma and gammaErr
                gamma = sum%x
                gammaErr = x*gammaErr + betaErr + (prod%y + sum%y)
                ! product and sum for beta
                prod = TwoProduct(beta,x)
                sum = TwoSum(prod%x,alpha)
                ! update beta and betaErr
                beta = sum%x
                betaErr = x*betaErr + alphaErr + (prod%y + sum%y)
                ! product and sum for alpha
                prod = TwoProduct(alpha,x)
                sum = TwoSum(prod%x,poly(j))
                ! update alpha and alphaErr
                alpha = sum%x
                alphaErr = x*alphaErr + (prod%y + sum%y)
            end do
            ! add error back into result
            gamma = gamma + gammaErr
            beta = beta + betaErr
            alpha = alpha + alphaErr
            ! check residual
            if(abs(alpha)<small) return
            ! compute Laguerre correction terms
            alpha = beta/alpha
            beta = alpha**2 - 2*(gamma/alpha)
            gamma = sqrt((deg-1)*(deg*beta-alpha**2))
            beta = alpha + gamma
            alpha = alpha - gamma
            ! update x
            if(abs(alpha)>abs(beta)) then
                x = x - (deg/alpha)
            else
                x = x - (deg/beta)
            end if
        end do
    end subroutine CLaguerre
    !********************************************************
    !                   Compensated Laguerre Cplx           *
    !********************************************************
    ! Complex version of Compensated Laguerre.
    !********************************************************
    subroutine CLaguerreCplx(poly,x,deg)
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: x
        ! local variables
        integer                         :: j, k
        complex(kind=dp)                :: alpha, alphaErr, beta, betaErr, gamma, gammaErr
        type(CEFTSum)                   :: sum
        type(CEFTProd)                  :: prod
        
        ! Laguerre method
        do k=1,itnum
            ! compute alpha and beta
            alpha = poly(deg+1)
            alphaErr = cmplx(0.0_dp,0.0_dp,kind=dp)
            beta = cmplx(0.0_dp,0.0_dp,kind=dp)
            betaErr = cmplx(0.0_dp,0.0_dp,kind=dp)
            gamma = cmplx(0.0_dp,0.0_dp,kind=dp)
            gammaErr = cmplx(0.0_dp,0.0_dp,kind=dp)
            do j=deg,1,-1
                ! product and sum for gamma
                prod = TwoProductCplx(gamma,x)
                sum = TwoSumCplx(prod%p,beta)
                ! update gamma and gammaErr
                gamma = sum%x
                gammaErr = x*gammaErr + betaErr + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
                ! product and sum for beta
                prod = TwoProductCplx(beta,x)
                sum = TwoSumCplx(prod%p,alpha)
                ! update beta and betaErr
                beta = sum%x
                betaErr = x*betaErr + alphaErr + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
                ! product and sum for alpha
                prod = TwoProductCplx(alpha,x)
                sum = TwoSumCplx(prod%p,poly(j))
                ! update alpha and alphaErr
                alpha = sum%x
                alphaErr = x*alphaErr + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            end do
            ! add error back into result
            gamma = gamma + gammaErr
            beta = beta + betaErr
            alpha = alpha + alphaErr
            ! check residual
            if(abs(alpha)<small) return
            ! compute Laguerre correction terms
            alpha = beta/alpha
            beta = alpha**2 - 2*(gamma/alpha)
            gamma = sqrt((deg-1)*(deg*beta-alpha**2))
            beta = alpha + gamma
            alpha = alpha - gamma
            ! update x
            if(abs(alpha)>abs(beta)) then
                x = x - (deg/alpha)
            else
                x = x - (deg/beta)
            end if
        end do
    end subroutine CLaguerreCplx
end program laguerre_limitAcc