!********************************************************************************
!   MPROUTINES: Multi-Precision Routines
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 29 March 2020
!********************************************************************************
module mproutines
    use eft
    use mpmodule
    implicit none
    
contains
    !********************************************************
    !                  Compute Condition                    *
    !********************************************************
    ! Given the polynomial, use multi-precision Horner method
    ! to compute the condition number of poly at x. The
    ! result is returned in comp which has type real(kind=dp). 
    !********************************************************
    function compCond(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        real(kind=dp)                   :: poly(:), x
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: comp
        type(mp_real)                   :: mp_comp, mp_poly(deg+1), mp_apoly(deg+1), alpha, beta, rho
            
        ! convert to mp_type
        do k=1,deg+1
            mp_poly(k) = mpreald(poly(k))
            mp_apoly(k) = abs(mp_poly(k))
        end do
        ! compute rho
        mp_comp = mpreald(abs(x))
        rho = mp_apoly(deg+1)
        do k=deg,1,-1
            rho = mp_comp*rho + mp_apoly(k)
        end do
        ! compute alpha and beta
        mp_comp = mpreald(x)
        alpha = mp_poly(deg+1)
        beta = mpreald(0.0_dp)
        do k=deg,1,-1
            beta = mp_comp*beta + alpha
            alpha = mp_comp*alpha + mp_poly(k)
        end do
        ! compute cond
        if(beta==mpreald(0.0_dp)) then
            mp_comp = mpreald(big)
        else
            mp_comp = rho/abs(mp_comp*beta)
        end if
        ! convert to dp-precision
        comp = mp_comp
        return
    end function compCond
    !********************************************************
    !                  Compute Condition   Cplx             *
    !********************************************************
    ! Complex version of Compute Condition.
    !********************************************************
    function compCondCplx(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        complex(kind=dp)                :: poly(:), x
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: comp
        type(mp_real)                   :: mp_comp, mp_apoly(deg+1), rho
        type(mp_complex)                :: cmp_comp, mp_poly(deg+1), alpha, beta
            
        ! convert to mp_type
        do k=1,deg+1
            mp_poly(k) = mpcmplxdc(poly(k))
            mp_apoly(k) = abs(mp_poly(k))
        end do
        ! compute rho
        mp_comp = mpreald(abs(x))
        rho = mp_apoly(deg+1)
        do k=deg,1,-1
            rho = mp_comp*rho + mp_apoly(k)
        end do
        ! compute alpha and beta
        cmp_comp = mpcmplxdc(x)
        alpha = mp_poly(deg+1)
        beta = mpcmplxdc(cmplx(0.0_dp,0.0_dp,kind=dp))
        do k=deg,1,-1
            beta = cmp_comp*beta + alpha
            alpha = cmp_comp*alpha + mp_poly(k)
        end do
        ! compute cond
        if(beta==mpcmplxdc(cmplx(0.0_dp,0.0_dp,kind=dp))) then
            mp_comp = mpreald(big)
        else
            mp_comp = rho/abs(cmp_comp*beta)
        end if
        ! convert to dp-precision
        comp = mp_comp
        return
    end function compCondCplx
    !********************************************************
    !                   limAccPoly                          *
    !********************************************************
    ! Creates test polynomial (x-1)^deg-10^(-8) for the 
    ! limiting accuracy testing. The result is stored in poly. 
    !********************************************************
    subroutine limAccPoly(poly,deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(inout)    :: poly(:)
        ! local variables
        type(mp_real)                   :: r(deg), c(deg), alpha
        integer                         :: i, j
        
        ! store mp-type arrays
        do i=1,deg
            r(i) = mpreald(1.0_dp)
            c(i) = mpreald(0.0_dp)
        end do
        ! create poly
        do i=1,deg
            alpha = -r(i)
            do j=i,1,-1
                if(j==1) then
                    c(j) = c(j) + alpha
                else
                    c(j) = c(j) + alpha*c(j-1)
                end if
            end do
        end do
        ! shift c(deg) by 10E-8
        c(deg) = c(deg) - mpreald(1.0E-8_dp)
        ! convert to dp-precision
        do i=1,deg
            poly(i) = c(deg+1-i)
        end do
        poly(deg+1) = 1.0_dp
    end subroutine limAccPoly
    !********************************************************
    !                   limAccPolyCplx                      *
    !********************************************************
    ! Complex version of limAccPoly.
    !********************************************************
    subroutine limAccPolyCplx(poly,deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(inout) :: poly(:)
        ! local variables
        type(mp_complex)                :: r(deg), c(deg), alpha
        integer                         :: i, j
        
        ! store mp-type arrays
        do i=1,deg
            r(i) = mpcmplxdc(cmplx(1.0_dp,1.0_dp,kind=dp))
            c(i) = mpcmplxdc(cmplx(0.0_dp,0.0_dp,kind=dp))
        end do
        ! create poly
        do i=1,deg
            alpha = -r(i)
            do j=i,1,-1
                if(j==1) then
                    c(j) = c(j) + alpha
                else
                    c(j) = c(j) + alpha*c(j-1)
                end if
            end do
        end do
        ! shift c(deg) by 10E-8
        c(deg) = c(deg) - mpcmplxdc(cmplx(1.0E-8_dp,0.0_dp,kind=dp))
        ! convert to dp-precision
        do i=1,deg
            poly(i) = c(deg+1-i)
        end do
        poly(deg+1) = cmplx(1.0_dp,0.0_dp,kind=dp)
    end subroutine limAccPolyCplx
    !********************************************************
    !                   limAccRoot                          *
    !********************************************************
    ! Given coeeficients of (x-1)^deg-10^(-8) stored in 
    ! dp-precision, use multi-precision Newton's method to 
    ! compute exact root. The result is returned in comp which
    ! has type real(kind=dp).
    !********************************************************
    function limAccRoot(poly,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        real(kind=dp)                   :: poly(:)
        ! local variables
        integer                         :: j, k
        real(kind=dp)                   :: comp
        type(mp_real)                   :: mp_comp, mp_poly(deg+1), alpha, beta
            
        ! convert to mp-type
        do k=1,deg+1
            mp_poly(k) = mpreald(poly(k))
        end do
        ! root approximation
        mp_comp = mpreald(1.0_dp + 10.0_dp**(-8.0_dp/real(deg,kind=dp)))
        ! Newton's method
        do k=1,5
            ! MPHorner applied to mp_mpoly
            alpha = mp_poly(deg+1)
            beta = mpreald(0.0_dp)
            do j=deg,1,-1
                beta = mp_comp*beta + alpha
                alpha = mp_comp*alpha + mp_poly(j)
            end do
            ! update root approximation
            mp_comp = mp_comp - (alpha/beta)
        end do
        ! convert to dp-precision
        comp = mp_comp
        return
    end function limAccRoot
    !********************************************************
    !                   limAccRootCplx                      *
    !********************************************************
    ! Complex version of limAccRoot.
    !********************************************************
    function limAccRootCplx(poly,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        complex(kind=dp)                :: poly(:)
        ! local variables
        integer                         :: j, k
        complex(kind=dp)                :: comp
        type(mp_complex)                :: mp_comp, mp_poly(deg+1), alpha, beta
            
        ! convert to mp-type
        do k=1,deg+1
            mp_poly(k) = mpcmplxdc(poly(k))
        end do
        ! root approximation
        mp_comp = mpcmplxdc(cmplx(1.0_dp + 10.0_dp**(-8.0_dp/real(deg,kind=dp)),1.0_dp,kind=dp))
        ! Newton's method
        do k=1,5
            ! MPHorner applied to mp_mpoly
            alpha = mp_poly(deg+1)
            beta = mpcmplxdc(cmplx(0.0_dp,0.0_dp,kind=dp))
            do j=deg,1,-1
                beta = mp_comp*beta + alpha
                alpha = mp_comp*alpha + mp_poly(j)
            end do
            ! update root approximation
            mp_comp = mp_comp - (alpha/beta)
        end do
        ! convert to dp-precision
        comp = mp_comp
        return
    end function limAccRootCplx
    !********************************************************
    !                   Multi-precision Horner              *
    !********************************************************
    ! Evaluate polynomial at floating-point number in multi-
    ! precision using MPFun and Horner's method. The result
    ! is returned in comp which has type real(kind=dp). 
    !********************************************************
    function MPHorner(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        real(kind=dp)                   :: poly(:), x
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: comp
        type(mp_real)                   :: mp_comp, mp_poly(deg+1), mp_x
        
        ! convert to mp-type
        do k=1,deg+1
            mp_poly(k) = mpreald(poly(k))
        end do
        mp_x = mpreald(x)
        ! Horner's method
        mp_comp = mp_poly(deg+1)
        do k=deg,1,-1
            mp_comp = mp_x*mp_comp + mp_poly(k)
        end do
        ! convert to dp-precision
        comp = mp_comp
        return
    end function MPHorner
    !********************************************************
    !                   Multi-precision Horner Cplx         *
    !********************************************************
    ! Complex version of MPHorner.
    !********************************************************
    function MPHornerCplx(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        complex(kind=dp)                :: poly(:), x
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: comp
        type(mp_complex)                :: mp_comp, mp_poly(deg+1), mp_x
        
        ! convert to mp-type
        do k=1,deg+1
            mp_poly(k) = mpcmplxdc(poly(k))
        end do
        mp_x = mpcmplxdc(x)
        ! Horner's method
        mp_comp = mp_poly(deg+1)
        do k=deg,1,-1
            mp_comp = mp_x*mp_comp + mp_poly(k)
        end do
        ! convert to dp-precision
        comp = mp_comp
        return
    end function MPHornerCplx    
    !********************************************************
    !                   Roots to Coeffs                     *
    !********************************************************
    ! Uses mpfun to compute the coefficients of a polynomial
    ! whose roots are given in roots. The result is stored
    ! in coeffs.
    !********************************************************
    subroutine rootCoeff(deg,roots,coeffs)
        implicit none
        ! argument variables
        integer, intent(in)           :: deg
        real(kind=dp), intent(in)     :: roots(:)
        real(kind=dp), intent(out)    :: coeffs(:)
        ! local variables
        type(mp_real)                 :: r(deg), c(deg+1), alpha
        integer                       :: i, j
        
        ! convert to mp-type
        do i=1,deg
            r(i) = mpreald(roots(i))
            c(i) = mpreald(0.0_dp)
        end do
		! expand polynomial
        do i=1,deg
            alpha = -r(i)
            do j=i,1,-1
                if(j==1) then
                    c(j) = c(j) + alpha
                else
                    c(j) = c(j) + alpha*c(j-1)
                end if
            end do
        end do
        ! convert to dp-precision
        do i=1,deg
            coeffs(i) = c(deg+1-i)
        end do
        ! multiply by leading coefficient
        do i=1,deg
            coeffs(i) = coeffs(deg+1)*coeffs(i)
        end do
    end subroutine rootCoeff
    !********************************************************
    !                   Roots to Coeffs Cplx                *
    !********************************************************
    ! Complex version of rootCoeff.
    !********************************************************
    subroutine rootCoeffCplx(deg,roots,coeffs)
        implicit none
        ! argument variables
        integer, intent(in)           :: deg
        complex(kind=dp), intent(in)  :: roots(:)
        complex(kind=dp), intent(out) :: coeffs(:)
        ! local variables
        type(mp_complex)              :: r(deg), c(deg+1), alpha
        integer                       :: i, j
        
        ! convert to mp-type
        do i=1,deg
            r(i) = mpcmplxdc(roots(i))
            c(i) = mpcmplxdc(cmplx(0.0_dp,0.0_dp,kind=dp))
        end do
		! expand polynomial
        do i=1,deg
            alpha = -r(i)
            do j=i,1,-1
                if(j==1) then
                    c(j) = c(j) + alpha
                else
                    c(j) = c(j) + alpha*c(j-1)
                end if
            end do
        end do
        ! convert to dp-precision
        do i=1,deg
            coeffs(i) = c(deg+1-i)
        end do
        ! multiply by leading coefficient
        do i=1,deg
            coeffs(i) = coeffs(deg+1)*coeffs(i)
        end do
    end subroutine rootCoeffCplx
end module mproutines