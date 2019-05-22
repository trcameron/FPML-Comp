!********************************************************************************
!   HORNER_TEST: Test Compensated Horner Routines
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 27 November 2018
!********************************************************************************
program horner_aprioriErr
    use eft
    use mproutines
    implicit none
    ! polynomial variables
    integer                         :: deg, k
    real(kind=dp)                   :: x
    complex(kind=dp)                :: z
    real(kind=dp), allocatable      :: poly(:), roots(:)
    complex(kind=dp), allocatable   :: cpoly(:), croots(:)
    ! testing variables
    real(kind=dp)                   :: g
    real(kind=dp), allocatable      :: bound(:), error(:), compbound(:), comperror(:), cond(:), exact(:)
    complex(kind=dp), allocatable   :: cexact(:)
    
    call mpinit
    
    ! allocate real test variables
    allocate(bound(40),error(40),compbound(40),comperror(40),cond(40),exact(40))
    ! run real test
    x = 1.333_dp
    do deg=3,42
        ! allocate polynomial variables
        allocate(poly(deg+1),roots(deg))
        ! roots (all ones)
        roots = (/ (1.0_dp, k=1,deg)/)
        ! polynomial
        call RootCoeff(deg,roots,poly)
        ! multi-precision real Horner
        exact(deg-2) = MPHorner(poly,x,deg)
        ! standard real Horner
        error(deg-2) = abs(SHorner(poly,x,deg) - exact(deg-2))/abs(exact(deg-2))
        if (error(deg-2)<mu) then
            error(deg-2) = mu
        else if (error(deg-2)>1) then
            error(deg-2) = 1
        end if
        ! compensated Horner
        comperror(deg-2) = abs(CHorner(poly,x,deg) - exact(deg-2))/abs(exact(deg-2))
        if (comperror(deg-2)<mu) then
            comperror(deg-2) = mu
        else if (comperror(deg-2)>1) then
            comperror(deg-2) = 1
        end if
        ! error bounds
        g = (2*deg*mu/(1-2*deg*mu))
        cond(deg-2) = MPHorner(abs(poly),abs(x),deg)/abs(exact(deg-2))
        bound(deg-2) = cond(deg-2)*g
        compbound(deg-2) = mu + cond(deg-2)*g**2
        if(bound(deg-2)>1) bound(deg-2) = 1
        if(compbound(deg-2)>1) compbound(deg-2) = 1
        ! deallocate polynomial variables
        deallocate(poly,roots)
    end do
    ! write real test results
    open(unit=1,file="data_files/horner_aprioriErr_real.dat")
    write(1,'(A)') 'cond, stan_err_bound, stan_err, comp_err_bound, comp_err'
    do k=1,40
        write(1,'(ES15.2)', advance='no') cond(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') bound(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') error(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') compbound(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') comperror(k)
    end do
    ! close file
    close(1)
    ! deallocate real testing variables
    deallocate(bound,error,compbound,comperror,cond,exact)
    
    ! allocate complex test variables
    allocate(bound(40),error(40),compbound(40),comperror(40),cond(40),cexact(40))
    ! run complex test
    z = cmplx(1.333,1.333,kind=dp)
    do deg=3,42
        ! allocate polynomial variables
        allocate(cpoly(deg+1),croots(deg))
        ! roots (all 1+i)
        croots = (/ (cmplx(1,1,kind=dp), k=1,deg)/)
        ! polynomial
        call RootCoeffCplx(deg,croots,cpoly)
        ! multi-precision complex Horner
        cexact(deg-2) = MPHornerCplx(cpoly,z,deg)
        ! standard complex Horner
        error(deg-2) = abs(SHornerCplx(cpoly,z,deg) - cexact(deg-2))/abs(cexact(deg-2))
        if (error(deg-2)<mu) then
            error(deg-2) = mu
        else if (error(deg-2)>1) then
            error(deg-2) = 1
        end if
        ! compensated complex Horner
        comperror(deg-2) = abs(CHornerCplx(cpoly,z,deg) - cexact(deg-2))/abs(cexact(deg-2))
        if (comperror(deg-2)<mu) then
            comperror(deg-2) = mu
        else if (comperror(deg-2)>1) then
            comperror(deg-2) = 1
        end if
        ! error bounds
        g = 2*(4*deg+2)*mu*sqrt(2.0_dp)
        g = g/((1-2*mu)-g)
        cond(deg-2) = MPHorner(abs(cpoly),abs(z),deg)/abs(cexact(deg-2))
        bound(deg-2) = cond(deg-2)*g
        compbound(deg-2) = mu + cond(deg-2)*g**2
        if(bound(deg-2)>1) bound(deg-2) = 1
        if(compbound(deg-2)>1) compbound(deg-2) = 1
        ! deallocate polynomial variables
        deallocate(cpoly,croots)
    end do
    ! write complex test results
    open(unit=1,file="data_files/horner_aprioriErr_complex.dat")
    write(1,'(A)') 'cond, stan_err_bound, stan_err, comp_err_bound, comp_err'
    do k=1,40
        write(1,'(ES15.2)', advance='no') cond(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') bound(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') error(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)', advance='no') compbound(k)
        write(1,'(A)', advance='no') ','
        write(1,'(ES15.2)') comperror(k)
    end do
    ! close file
    close(1)
    ! deallocate complex test variables
    deallocate(bound,error,compbound,comperror,cond,cexact)

contains

    !********************************************************
    !                   Standard Horner                     *
    !********************************************************
    ! Evaluate polynomial at floating-point number in dp
    ! precision. The result is returned in comp which has 
    ! type real(kind=dp). 
    !********************************************************
    function SHorner(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        real(kind=dp)                   :: poly(:), x
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: comp
        
        ! Horner's method
        comp = poly(deg+1)
        do k=deg,1,-1
            comp = x*comp + poly(k)
        end do
        return
    end function SHorner
    !********************************************************
    !                   Standard Complex Horner             *
    !********************************************************
    ! Evaluate complex polynomial at complex floating-point 
    ! number in dp-precision. The result is returned in comp
    ! which has type complex(kind=dp). 
    !********************************************************
    function SHornerCplx(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                 :: deg
        complex(kind=dp)        :: poly(:), x
        ! local variables
        integer                 :: k
        complex(kind=dp)        :: comp
        
        ! Horner's method
        comp = poly(deg+1)
        do k=deg,1,-1
            comp = x*comp + poly(k)
        end do
        return
    end function SHornerCplx
    !****************************************************
    !                       HornerSum                   *
    !****************************************************
    ! Computes the evaluation of the sum of two
    ! polynomials of degree (deg-1), at a floating-point 
    ! number. Result is returned in comp which has type 
    ! real(kind=dp).
    !****************************************************
    function HSum(p,q,x,deg) result(comp)
        implicit none
        ! argument variables
        integer             :: deg
        real(kind=dp)       :: p(:), q(:), x
        ! local variables
        integer             :: k
        real(kind=dp)       :: comp
        
        ! Horner's method
        comp = p(deg) + q(deg)
        do k=deg-1,1,-1
            comp = x*comp + (p(k) + q(k))
        end do
        return
    end function HSum
    !****************************************************
    !              Compensated Horner                   *
    !****************************************************
    ! Computes the evaluation of a polynomial at
    ! a floating point number using the compensated
    ! Horner routine. The result is returned in comp
    ! which has type real(kind=dp).
    !****************************************************
    function CHorner(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer             :: deg
        real(kind=dp)       :: poly(:), x
        ! local variables
        real(kind=dp)       :: comp
        type(REFTHorner)    :: eftcomp
            
        ! compute EFTHorner
        eftcomp = EFTHorner(poly,x,deg)
        ! compute error using HSum and add back to result
        comp = eftcomp%h + HSum(eftcomp%p,eftcomp%q,x,deg)
        ! deallocate error polynomials in eftcomp
        deallocate(eftcomp%p,eftcomp%q)
        return
    end function CHorner
    !****************************************************
    !                       HornerSumCplx               *
    !****************************************************
    ! Computes the evaluation of the sum of four
    ! complex polynomials of degree deg-1, at a complex 
    ! floating-point number. Result is returned in comp 
    ! which has type complex(kind=dp).
    !****************************************************
    function HSumCplx(p,q,r,s,x,deg) result(comp)
        implicit none
        ! argument variables
        integer             :: deg
        complex(kind=dp)    :: p(:), q(:), r(:), s(:), x
        ! local variables
        integer             :: k
        complex(kind=dp)    :: comp
        
        ! Horner's method
        comp = FaithSumCplx(p(deg),q(deg),r(deg),s(deg))
        do k=deg-1,1,-1
            comp = x*comp + FaithSumCplx(p(k),q(k),r(k),s(k))
        end do
        return
    end function HSumCplx
    !****************************************************
    !              Compensated Horner Cplx              *
    !****************************************************
    ! Computes the evaluation of a complex polynomial at
    ! a complex floating point number using the 
    ! compensated Horner routine. The result is returned 
    ! in comp which has type real(complex=dp).
    !****************************************************
    function CHornerCplx(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer             :: deg
        complex(kind=dp)    :: poly(:), x
        ! local variables
        complex(kind=dp)    :: comp
        type(CEFTHorner)    :: eftcomp
            
        ! compute EFTHornerCplx
        eftcomp = EFTHornerCplx(poly,x,deg)
        ! compute error using HSumCplx and add back to result
        comp = eftcomp%h + HSumCplx(eftcomp%p,eftcomp%q,eftcomp%r,eftcomp%s,x,deg)
        ! deallocate error polynomials in eftcomp
        deallocate(eftcomp%p,eftcomp%q,eftcomp%r,eftcomp%s)
        return
    end function CHornerCplx
end program horner_aprioriErr