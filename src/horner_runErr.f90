!********************************************************************************
!   HORNER_RUNERR: Test compensated Horner running error bound
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 29 March 2020
!********************************************************************************
program horner_runErr
    use eft 
    use mproutines
    implicit none
    ! paramters
    integer, parameter              :: deg = 5, n = 2000
    real(kind=dp), parameter        :: dx = 1E-5_dp
    ! polynomial variables
    real(kind=dp), allocatable      :: poly(:), roots(:)
    complex(kind=dp), allocatable   :: cpoly(:), croots(:)
    ! testing variables
    integer                         :: k
    real(kind=dp)                   :: abound, comp, error, exact, g, rbound, x
    complex(kind=dp)                :: ccomp, cexact, z
    
    ! allocate polynomial variables
    allocate(poly(deg+1), roots(deg))
    roots = (/ (1.0_dp, k=1,deg)/)
	poly(deg+1) = 1.0_dp
    call rootCoeff(deg,roots,poly)
    ! open file
    open(unit=1,file="data_files/horner_runErr_real.dat")
    write(1,'(A)') 'x, a priori error bound, running error bound, forward error'
    ! run real test
    x = 0.99_dp
    do k=1,n+1
        ! multi-precision real Horner
        exact = MPHorner(poly,x,deg)
        ! compensated Horner and running error bound
        call CHorner(poly,x,deg,comp,rbound)
        ! a priori error bound
        g = 2*deg*mu/(1-2*deg*mu)
        abound = mu*abs(exact) + MPHorner(abs(poly),abs(x),deg)*g**2
        ! forward error
        error = CHornerErr(poly,x,deg,comp)
        ! write real test results
        write(1,'(ES15.5,A)', advance='no') x, ', '
        write(1,'(ES15.5,A)', advance='no') abound, ', '
        write(1,'(ES15.5,A)', advance='no') rbound, ', '
        write(1,'(ES15.5)') error
        ! update x
        x = x + dx
    end do
    ! deallocate polynomial variables
    deallocate(poly, roots)
    ! close file
    close(1)
    
    ! allocate polynomial variables
    allocate(cpoly(deg+1), croots(deg))
    croots = (/ (cmplx(1.0_dp,1.0_dp,kind=dp), k=1,deg)/)
	cpoly(deg+1) = cmplx(1.0_dp,0.0_dp,kind=dp)
    call rootCoeffCplx(deg,croots,cpoly)
    ! open file
    open(unit=1,file="data_files/horner_runErr_cmplx.dat")
    write(1,'(A)') 'x, a priori error bound, running error bound, forward error'
    ! run complex test
    z = cmplx(0.99_dp,1.0_dp,kind=dp)
    do k=1,n+1
        ! multi-precision complex Horner
        cexact = MPHornerCplx(cpoly,z,deg)
        ! compensated Horner and running error bound
        call CHornerCplx(cpoly,z,deg,ccomp,rbound)
        ! a priori error bound
        g = 2*mu/(1-2*mu)
        g = 2*deg*sqrt(2.0_dp)*g/(1-2*deg*sqrt(2.0_dp)*g)
        abound = mu*abs(cexact) + MPHorner(abs(cpoly),abs(z),deg)*g**2
        ! forward error
        error = CHornerErrCplx(cpoly,z,deg,ccomp)
        ! write real test results
        write(1,'(ES15.5,A)', advance='no') real(z), ', '
        write(1,'(ES15.5,A)', advance='no') abound, ', '
        write(1,'(ES15.5,A)', advance='no') rbound, ', '
        write(1,'(ES15.5)') error
        ! update z
        z = z + cmplx(dx,0,kind=dp)
    end do
    ! deallocate polynomial variables
    deallocate(cpoly, croots)
    ! close file
    close(1)
contains
    
    !****************************************************
    !                       Horner Sum                  *
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
    !              Compensated Horner                   *
    !****************************************************
    ! Computes the evaluation of a polynomial at
    ! a floating point number using the compensated
    ! Horner routine. The result is stored in comp and the
    ! running relative forward error is stored in error.
    !****************************************************
    subroutine CHorner(poly,x,deg,comp,error)
        implicit none
        ! argument variables
        integer, intent(in)         :: deg
        real(kind=dp), intent(in)   :: poly(:), x
        real(kind=dp), intent(out)  :: comp, error
        ! local variables
        real(kind=dp)               :: g
        type(REFTHorner)            :: eft
            
        ! compute EFTHorner
        eft = EFTHorner(poly,x,deg)
        ! compute error using HSum and add back to result
        comp = eft%h + HSum(eft%p,eft%q,x,deg)
        !  compute running forward error
        g = (4*deg+2)*mu/(1-(4*deg+2)*mu)
        error = mu*abs(comp) + (2*mu**2*abs(comp) + g*HSum(abs(eft%p),abs(eft%q),abs(x),deg))
        ! deallocate error polynomials in eft
        deallocate(eft%p,eft%q)
    end subroutine CHorner
    !****************************************************
    !              Compensated Horner Cplx              *
    !****************************************************
    ! Computes the evaluation of a polynomial at
    ! a floating point number using the compensated
    ! Horner routine. The result is stored in comp and the
    ! running relative forward error is stored in error.
    !****************************************************
    subroutine CHornerCplx(poly,x,deg,comp,error)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: poly(:), x
        real(kind=dp), intent(out)      :: error
        complex(kind=dp), intent(out)   :: comp
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: g
        type(CEFTHorner)                :: eft
            
        ! compute EFTHornerCplx
        eft = EFTHornerCplx(poly,x,deg)
        ! compute error using HSumCplx and add back to result
        comp = eft%h + HSumCplx(eft%p,eft%q,eft%r,eft%s,x,deg)
        !  compute running forward error
        g = 2*mu/(1-2*mu)
        g = (4*deg+2)*sqrt(2.0_dp)*g/(1-(4*deg+2)*sqrt(2.0_dp)*g)
        error = FaithSum(abs(eft%p(deg)),abs(eft%q(deg)),abs(eft%r(deg)),abs(eft%s(deg)))
        do k=deg-1,1,-1
            error = abs(x)*error + FaithSum(abs(eft%p(k)),abs(eft%q(k)),abs(eft%r(k)),abs(eft%s(k)))
        end do
        error = mu*abs(comp) + (g*error + 2*mu**2*abs(comp))
        ! deallocate error polynomials in eft
        deallocate(eft%p,eft%q,eft%r,eft%s)
    end subroutine CHornerCplx
    !****************************************************
    !              Compensated Horner Error             *
    !****************************************************
    ! Compute error in compensated Horner evaluation vs.
    ! MP Horner evaluation. 
    !****************************************************
    function CHornerErr(poly,x,deg,comp) result(err)
        implicit none
        ! argument variables
        integer                 :: deg
        real(kind=dp)           :: poly(:), x, comp
        ! local variables
        integer                 :: k
        real(kind=dp)           :: err
        type(mp_real)           :: mp_comp, mp_err, mp_horner, mp_poly(deg+1), mp_x
            
        ! convert to mp-type
        mp_comp = mpreald(comp)
        do k=1,deg+1
            mp_poly(k) = mpreald(poly(k))
        end do
        mp_x = mpreald(x)
        ! Horner's method
        mp_horner = mp_poly(deg+1)
        do k=deg,1,-1
            mp_horner = mp_x*mp_horner + mp_poly(k)
        end do
        ! error
        mp_err = abs(mp_comp - mp_horner)
        ! convert to dp-precision
        err = mpreal(mp_err)
    end function CHornerErr
    !****************************************************
    !              Compensated Horner Error Cplx        *
    !****************************************************
    ! Compute error in compensated Horner evaluation vs.
    ! MP Horner evaluation. 
    !****************************************************
    function CHornerErrCplx(poly,x,deg,comp) result(err)
        implicit none
        ! argument variables
        integer                 :: deg
        complex(kind=dp)        :: poly(:), x, comp
        ! local variables
        integer                 :: k
        real(kind=dp)           :: err
        type(mp_real)           :: mp_err
        type(mp_complex)        :: mp_comp, mp_horner, mp_poly(deg+1), mp_x
            
        ! convert to mp-type
        mp_comp = mpcmplxdc(comp)
        do k=1,deg+1
            mp_poly(k) = mpcmplxdc(poly(k))
        end do
        mp_x = mpcmplxdc(x)
        ! Horner's method
        mp_horner = mp_poly(deg+1)
        do k=deg,1,-1
            mp_horner = mp_x*mp_horner + mp_poly(k)
        end do
        ! error
        mp_err = abs(mp_comp - mp_horner)
        ! convert to dp-precision
        err = mpreal(mp_err)
    end function CHornerErrCplx
end program horner_runErr