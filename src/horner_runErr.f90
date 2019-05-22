!********************************************************************************
!   HORNER_DEB_TEST: Test Compensated Horner Dynamic Error Bound
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 21 March 2019
!********************************************************************************
program horner_runErr
    use eft
    use mproutines
    implicit none
    ! paramters
    integer, parameter              :: deg = 5, n = 2000        ! degree and number of trials
    real(kind=dp), parameter        :: dx = 1E-5_dp             ! change in x after each trial
    ! polynomial variables
    integer                         :: k
    real(kind=dp)                   :: x
    complex(kind=dp)                :: z
    real(kind=dp), allocatable      :: poly(:), roots(:)
    complex(kind=dp), allocatable   :: cpoly(:), croots(:)
    ! testing variables
    real(kind=dp)                   :: abound, comp, error, exact, g, rbound
    complex(kind=dp)                :: ccomp, cexact
    
    call mpinit
    
    ! allocate polynomial variables
    allocate(poly(deg+1), roots(deg))
    roots = (/ (1.0_dp, k=1,deg)/)
    call RootCoeff(deg,roots,poly)
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
        error = MPHornerErr(poly,x,deg,comp)
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
    croots = (/ (cmplx(1,1,kind=dp), k=1,deg)/)
    call RootCoeffCplx(deg,croots,cpoly)
    ! open file
    open(unit=1,file="data_files/horner_runErr_complex.dat")
    write(1,'(A)') 'x, a priori error bound, running error bound, forward error'
    ! run complex test
    z = cmplx(0.99_dp,1.0_dp,kind=dp)
    do k=1,n+1
        ! multi-precision complex Horner
        cexact = MPHornerCplx(cpoly,z,deg)
        ! compensated Horner and running error bound
        call CHornerCplx(cpoly,z,deg,ccomp,rbound)
        ! a priori error bound
        g = 4*deg*mu*sqrt(2.0_dp)
        g = g/((1-2*mu)-g)
        abound = mu*abs(cexact) + MPHorner(abs(cpoly),abs(z),deg)*g**2
        ! forward error
        error = MPHornerCplxErr(cpoly,z,deg,ccomp)
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
    !                       MPHornerErr                 *
    !****************************************************
    function MPHornerErr(poly,x,deg,comp) result(err)
        implicit none
        ! argument variables
        integer                         :: deg
        real(kind=dp)                   :: poly(:), x, comp
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: err
        type(mp_real)                   :: mp_comp, mp_err, mp_horn, mp_poly(deg+1), mp_x
            
        ! convert to mp-type
        mp_comp = comp
        do k=1,deg+1
            mp_poly(k) = poly(k)
        end do
        mp_x = x
        ! Horner's method
        mp_horn = mp_poly(deg+1)
        do k=deg,1,-1
            mp_horn = mp_x*mp_horn + mp_poly(k)
        end do
        ! error
        mp_err = abs(mp_comp - mp_horn)
        ! convert to dp-precision
        err = mp_err
        return 
    end function MPHornerErr
    !****************************************************
    !                       MPHornerErr                 *
    !****************************************************
    function MPHornerCplxErr(poly,x,deg,comp) result(err)
        implicit none
        ! argument variables
        integer                         :: deg
        complex(kind=dp)                :: poly(:), x, comp
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: err
        type(mp_real)                   :: mp_err
        type(mp_complex)                :: mp_comp, mp_horn, mp_poly(deg+1), mp_x
            
        ! convert to mp-type
        mp_comp = comp
        do k=1,deg+1
            mp_poly(k) = poly(k)
        end do
        mp_x = x
        ! Horner's method
        mp_horn = mp_poly(deg+1)
        do k=deg,1,-1
            mp_horn = mp_x*mp_horn + mp_poly(k)
        end do
        ! error
        mp_err = abs(mp_comp - mp_horn)
        ! convert to dp-precision
        err = mp_err
        return 
    end function MPHornerCplxErr
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
        type(REFTHorner)            :: eftcomp
            
        ! compute EFTHorner
        eftcomp = EFTHorner(poly,x,deg)
        ! compute error using HSum and add back to result
        comp = eftcomp%h + HSum(eftcomp%p,eftcomp%q,x,deg)
        !  compute running forward error
        g = (4*deg+2)*mu/(1-(4*deg+2)*mu)
        error = mu*abs(comp) + (2*mu**2*abs(comp) + g*HSum(abs(eftcomp%p),abs(eftcomp%q),abs(x),deg))
        ! deallocate error polynomials in eftcomp
        deallocate(eftcomp%p,eftcomp%q)
    end subroutine CHorner
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
        type(CEFTHorner)                :: eftcomp
            
        ! compute EFTHornerCplx
        eftcomp = EFTHornerCplx(poly,x,deg)
        ! compute error using HSum and add back to result
        comp = eftcomp%h + HSumCplx(eftcomp%p,eftcomp%q,eftcomp%r,eftcomp%s,x,deg)
        !  compute running forward error
        g = 2*(4*deg+2)*mu*sqrt(2.0_dp)
        g = g/((1-2*mu)-g)
        error = FaithSum(abs(eftcomp%p(deg)),abs(eftcomp%q(deg)),abs(eftcomp%r(deg)),abs(eftcomp%s(deg)))
        do k=deg-1,1,-1
            error = abs(x)*error + FaithSum(abs(eftcomp%p(k)),abs(eftcomp%q(k)),abs(eftcomp%r(k)),abs(eftcomp%s(k)))
        end do
        error = mu*abs(comp) + (g*error + 2*mu**2*abs(comp))
        ! deallocate error polynomials in eftcomp
        deallocate(eftcomp%p,eftcomp%q,eftcomp%r,eftcomp%s)
    end subroutine CHornerCplx
end program horner_runErr