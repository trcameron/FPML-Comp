!********************************************************************************
!   MPROUTINES: Multi-Precision Routines
!   Authors: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 21 March 2019
!********************************************************************************
module mproutines
    use eft
    use mpmodule
    implicit none
    
contains
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
            mp_poly(k) = poly(k)
        end do
        mp_x = x
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
    !                   Multi-precision Complex Horner      *
    !********************************************************
    ! Evaluate complex polynomial at complex floating-point 
    ! number in multi-precision using MPFun and Horner's method.
    ! The result is returned in comp which has type complex(kind=dp). 
    !********************************************************
    function MPHornerCplx(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                     :: deg
        complex(kind=dp)            :: poly(:), x
        ! local variables
        integer                     :: k
        complex(kind=dp)            :: comp
        type(mp_complex)            :: mp_comp, mp_poly(deg+1), mp_x
            
        ! covert to mp-type
        do k=1,deg+1
            mp_poly(k) = poly(k)
        end do
        mp_x = x
        ! Horner's method
        mp_comp = mp_poly(deg+1)
        do k=deg,1,-1
            mp_comp = mp_x*mp_comp + mp_poly(k)
        end do
        ! covert to dp-precision
        comp = mp_comp
        return
    end function MPHornerCplx
    !********************************************************
    !                   Roots to Coeffs                     *
    !********************************************************
    ! Uses mpfun to computes the coefficients of a polynomial
    ! whose roots are given in roots. The result is stored
    ! in coeffs. 
    !********************************************************
    subroutine RootCoeff(deg,roots,coeffs)
      implicit none
      ! argument variables
      integer, intent(in)           :: deg
      real(kind=dp), intent(in)     :: roots(:)
      real(kind=dp), intent(out)    :: coeffs(:)
      ! local variables
      type(mp_real)                 :: rr(deg), c(deg+1), cx(deg+1), cz(deg+1)
      integer                       :: ii, jj
      
      ! convert to mp-type
      do ii=1,deg
          rr(ii) = roots(ii)
          c(ii) = 0.0_dp
      end do
      c(deg+1) = 1.0_dp
      
      do ii=1,deg
          ! shift coefficients left
          cx(1:deg) = c(2:deg+1)
          cx(deg+1) = 0.0_dp
          ! mulltiply coefficients by new root
          do jj=1,deg+1
              cz(jj) = c(jj)*rr(ii)
          end do
          ! adjust c as difference between cx and cz
          do jj=1,deg+1
              c(jj) = cx(jj) - cz(jj)
          end do
      end do
      
      ! convert to dp-precision
      do ii=1,deg+1
          coeffs(ii) = c(deg+2-ii)
      end do
    end subroutine RootCoeff
    !********************************************************
    !                   Roots to Coeffs Cplx                *
    !********************************************************
    ! Uses mpfun to compute the coefficients of a complex
    ! polynomial whose roots are given in roots. The result
    ! is stored in coeffs (monic). 
    !********************************************************
    subroutine RootCoeffCplx(deg,roots,coeffs)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(out)   :: coeffs(:)
        ! local variables
        type(mp_complex)                :: rr(deg), c(deg+1), cx(deg+1), cz(deg+1)
        integer                         :: ii, jj
        
        ! convert to mp-type
        do ii=1,deg
            rr(ii) = roots(ii)
            c(ii) = cmplx(0,0,kind=dp)
        end do
        c(deg+1) = cmplx(1,0,kind=dp)
        
        do ii=1,deg
            ! shift coefficients left
            cx(1:deg) = c(2:deg+1)
            cx(deg+1) = cmplx(0,0,kind=dp)
            ! mulltiply coefficients by new root
            do jj=1,deg+1
                cz(jj) = c(jj)*rr(ii)
            end do
            ! adjust c as difference between cx and cz
            do jj=1,deg+1
                c(jj) = cx(jj) - cz(jj)
            end do
        end do
        
        ! convert to dp-precision
        do ii=1,deg+1
            coeffs(ii) = c(deg+2-ii)
        end do
    end subroutine RootCoeffCplx
    !********************************************************
    !                   Roots to Coeffs Cplx 2              *
    !********************************************************
    ! Uses mpfun to computes the coefficients of a complex
    ! polynomial whose roots are given in roots. The result
    ! is stored in coeffs (non-monic).
    !********************************************************
    subroutine RootCoeffCplx2(deg,roots,coeffs)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: roots(:)
        complex(kind=dp), intent(inout) :: coeffs(:)
        ! local variables
        type(mp_complex)                :: rr(deg), c(deg+1), cx(deg+1), cz(deg+1), alpha
        integer                         :: ii, jj
        
        ! convert to mp-type
        do ii=1,deg
            rr(ii) = roots(ii)
            c(ii) = cmplx(0,0,kind=dp)
        end do
        c(deg+1) = cmplx(1,0,kind=dp)
        alpha = coeffs(deg+1)
        
        do ii=1,deg
            ! shift coefficients left
            cx(1:deg) = c(2:deg+1)
            cx(deg+1) = cmplx(0,0,kind=dp)
            ! mulltiply coefficients by new root
            do jj=1,deg+1
                cz(jj) = c(jj)*rr(ii)
            end do
            ! adjust c as difference between cx and cz
            do jj=1,deg+1
                c(jj) = cx(jj) - cz(jj)
            end do
        end do
        ! multiply by leading coefficient
        do ii=1,deg+1
            c(ii) = alpha*c(ii)
        end do
        
        ! convert to dp-precision
        do ii=1,deg+1
            coeffs(ii) = c(deg+2-ii)
        end do
    end subroutine RootCoeffCplx2
    !********************************************************
    !                   Create Test Poly                    *
    !********************************************************
    ! Creates test polynomial (x-1)^deg-10^(-8). The result 
    ! is stored in poly. 
    !********************************************************
    subroutine createpoly(poly,deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        real(kind=dp), intent(inout)    :: poly(:)
        ! local variables
        type(mp_real)                   :: rr(deg), c(deg), alpha
        integer                         :: ii, jj
        
        ! store mp-type arrays
        do ii=1,deg
            rr(ii) = 1.0_dp
            c(ii) = 0.0_dp
        end do
        ! create poly
        do ii=1,deg
            alpha = -rr(ii)
            do jj=ii,1,-1
                if(jj==1) then
                    c(jj) = c(jj) + alpha*1.0_dp
                else
                    c(jj) = c(jj) + alpha*c(jj-1)
                end if
            end do
        end do
        ! shift c(deg) by 10^-8
        alpha = 1.0E-8_dp
        c(deg) = c(deg) - alpha
        ! convert to dp-precision
        do ii=1,deg
            poly(ii) = c(deg+1-ii)
        end do
        poly(deg+1) = 1.0_dp
    end subroutine createpoly
    !********************************************************
    !                   Create Test Poly cplx               *
    !********************************************************
    ! Creates test polynomial (x-(1+i))^deg-10^(-8). 
    ! The result is stored in poly. 
    !********************************************************
    subroutine createpolycplx(poly,deg)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(inout) :: poly(:)
        ! local variables
        type(mp_complex)                :: rr(deg), c(deg), alpha
        integer                         :: ii, jj
        
        ! store mp-type arrays
        do ii=1,deg
            rr(ii) = cmplx(1,1,kind=dp)
            c(ii) = cmplx(0,0,kind=dp)
        end do
        ! create poly
        do ii=1,deg
            alpha = -rr(ii)
            do jj=ii,1,-1
                if(jj==1) then
                    c(jj) = c(jj) + alpha*cmplx(1,0,kind=dp)
                else
                    c(jj) = c(jj) + alpha*c(jj-1)
                end if
            end do
        end do
        ! shift c(deg) by 10^-8
        alpha = cmplx(1.0E-8,0,kind=dp)
        c(deg) = c(deg) - alpha
        ! convert to dp-precision
        do ii=1,deg
            poly(ii) = c(deg+1-ii)
        end do
        poly(deg+1) = cmplx(1,0,kind=dp)
    end subroutine createpolycplx
    !********************************************************
    !                   Create Root                         *
    !********************************************************
    ! Given the polynomial, use multi-precision Newton method 
    ! to compute an exact root. The result is returned in comp 
    ! which has type real(kind=dp). 
    !********************************************************
    function createroot(poly,deg) result(comp)
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
            mp_poly(k) = poly(k)
        end do
        ! root approximation
        mp_comp = 1.0_dp + 10**(-8.0_dp/real(deg,kind=dp))
        ! Newton's method
        do k=1,10
            ! MPHorner applied to mp_poly
            alpha = mp_poly(deg+1)
            beta = 0.0_dp
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
    end function createroot
    !********************************************************
    !                   Create Root Cplx                    *
    !********************************************************
    ! Given the  complex polynomial, use multi-precision 
    ! Newton method to compute an exact complex root. The
    ! result is returned in comp which has type complex(kind=dp). 
    !********************************************************
    function createrootcplx(poly,deg) result(comp)
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
            mp_poly(k) = poly(k)
        end do
        ! root approximation
        mp_comp = cmplx(1.0_dp + 10**(-8.0_dp/real(deg,kind=dp)),1,kind=dp)
        ! Newton's method
        do k=1,10
            ! MPHorner applied to mp_poly
            alpha = mp_poly(deg+1)
            beta = cmplx(0,0,kind=dp)
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
    end function createrootcplx
    !********************************************************
    !                  Compute Condition                    *
    !********************************************************
    ! Given the polynomial, use multi-precision Horner method
    ! to compute the condition number of poly at x. The
    ! result is returned in comp which has type real(kind=dp). 
    !********************************************************
    function compcond(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        real(kind=dp)                   :: poly(:), x
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: comp
        type(mp_real)                   :: mp_comp, mp_poly(deg+1), mp_apoly(deg+1), alpha, beta, rho
            
        ! convert to mp-type
        do k=1,deg+1
            mp_poly(k) = poly(k)
            mp_apoly(k) = abs(poly(k))
        end do
        ! compute rho
        mp_comp = abs(x)
        rho = mp_apoly(deg+1)
        do k=deg,1,-1
            rho = mp_comp*rho + mp_apoly(k)
        end do
        ! compute alpha and beta
        mp_comp = x
        alpha = mp_poly(deg+1)
        beta = 0.0_dp
        do k=deg,1,-1
            beta = mp_comp*beta + alpha
            alpha = mp_comp*alpha + mp_poly(k)
        end do
        ! compute cond
        if(beta==0) then
            mp_comp = huge(1.0_dp)
        else
            mp_comp = rho/abs(mp_comp*beta)
        end if
        ! convert to dp-precision
        comp = mp_comp
        return
    end function compcond
    !********************************************************
    !                  Compute Condition Cplx               *
    !********************************************************
    ! Given the complex polynomial, use multi-precision 
    ! Horner method to compute the condition number of poly at x.
    ! The result is returned in comp which has type real(kind=dp).
    !********************************************************
    function compcondcplx(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer                         :: deg
        complex(kind=dp)                :: poly(:), x
        ! local variables
        integer                         :: k
        real(kind=dp)                   :: comp
        type(mp_real)                   :: mp_comp, mp_poly(deg+1), rho
        type(mp_complex)                :: mp_compcplx, mp_polycplx(deg+1), alpha, beta
            
        ! convert to mp-type
        do k=1,deg+1
            mp_polycplx(k) = poly(k)
            mp_poly(k) = abs(poly(k))
        end do
        ! compute rho
        mp_comp = abs(x)
        rho = mp_poly(deg+1)
        do k=deg,1,-1
            rho = mp_comp*rho + mp_poly(k)
        end do
        ! compute alpha and beta
        mp_compcplx = x
        alpha = mp_polycplx(deg+1)
        beta = cmplx(0,0,kind=dp)
        do k=deg,1,-1
            beta = mp_compcplx*beta + alpha
            alpha = mp_compcplx*alpha + mp_polycplx(k)
        end do
        ! compute cond
        if(beta==0) then
            mp_comp = huge(1.0_dp)
        else
            mp_comp = rho/abs(mp_compcplx*beta)
        end if
        ! convert to dp-precision
        comp = mp_comp
        return
    end function compcondcplx
end module mproutines