!********************************************************************************
!   FPML_COMP: Compensated Polishing Technique for FPML based on Modified Laguerre's Method
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 31 March 2020
!********************************************************************************
module fpml_comp
	use eft
    use fpml, only : check_nan_inf, dp, main
    implicit none
    
contains
    !************************************************
    !           main comp                           *
    !************************************************
    subroutine main_comp(poly, deg, roots, itmax)
        ! argument variables
        integer, intent(in)             :: deg, itmax
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        real(kind=dp)                   :: alpha(deg+1), r
		complex(kind=dp)                :: b, c, z
        ! fpml variables
        integer                         :: conv(deg)
        real(kind=dp)                   :: berr(deg), cond(deg)
        
        ! call fpml main
        call main(poly, deg, roots, berr, cond, conv, itmax)
        ! polish root approximations
        alpha = abs(poly)
        conv = (/ (0, i=1,deg)/)
        nz = 0
        do i=1,(itmax/2)
            do j=1,deg
                if(conv(j)==0) then
                    z = roots(j)
                    r = abs(z)
                    if(r > 1.0_dp) then
						call rcheck_lag(poly, deg, b, c, z, r, conv(j))
                    else
						call check_lag(poly, deg, b, c, z, r, conv(j))
                    end if
                    if(conv(j)==0) then
                        call modify_lag(deg, b, c, z, j, roots, conv(j))
                        roots(j) = z - c
                    else
                        nz = nz + 1
                        if(nz==deg) return
                    end if
                end if
            end do
        end do
    end subroutine main_comp
    !************************************************
    !                       rcheck_lag              *
    !************************************************
    subroutine rcheck_lag(poly, deg, b, c, z, r, conv)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        integer, intent(out)            :: conv
        real(kind=dp), intent(in)       :: r
        complex(kind=dp), intent(in)    :: poly(:), z
        complex(kind=dp), intent(out)   :: b, c
        ! local variables
        integer                         :: k
		real(kind=dp)					:: g, rr, err
        complex(kind=dp)                :: a, a_err, b_err, c_err, zz
        type(CEFTSum)                   :: sum
        type(CEFTProd)                  :: prod
        ! intrinsic functions
        intrinsic                       :: abs
		
        ! evaluate polynomial and derivatives
        zz = 1/z
        rr = 1/r
        a = poly(1)
        a_err = cmplx(0.0_dp,0.0_dp,kind=dp)
        b = cmplx(0.0_dp,0.0_dp,kind=dp)
        b_err = cmplx(0.0_dp,0.0_dp,kind=dp)
        c = cmplx(0.0_dp,0.0_dp,kind=dp)
        c_err = cmplx(0.0_dp,0.0_dp,kind=dp)
		err = 0.0_dp
        do k=2,deg+1
            ! product and sum for c
            prod = TwoProductCplx(c,zz)
            sum = TwoSumCplx(prod%p,b)
            ! update c and c_err
            c = sum%x
            c_err = zz*c_err + b_err + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            ! product and sum for b
            prod = TwoProductCplx(b,zz)
            sum = TwoSumCplx(prod%p,a)
            ! update a and a_berr
            b = sum%x
            b_err = zz*b_err + a_err + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            ! product and sum for a
            prod = TwoProductCplx(a,zz)
            sum = TwoSumCplx(prod%p,poly(k))
            ! update a and a_err
            a = sum%x
            a_err = zz*a_err + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            ! update error bound
            err = rr*err + FaithSum(abs(prod%e),abs(prod%f),abs(prod%g),abs(sum%y))
        end do
        ! add error back into result
        c = c + c_err
        b = b + b_err
        a = a + a_err
        ! berr
        g = (4*deg+2)*sqrt(2.0_dp)*eps/(1-eps)
        g = g/(1-g)
		err = mu*abs(a) + (g*err + 2*mu**2*abs(a))
		! laguerre correction terms/ convergence
		if(abs(a)>err) then
            b = b/a
            c = 2*(c/a)
            c = zz**2*(deg-2*zz*b+zz**2*(b**2-c))
            b = zz*(deg-zz*b)
            if(check_nan_inf(b) .or. check_nan_inf(c)) conv = -1
		else
			conv = 1
		end if
    end subroutine rcheck_lag
    !************************************************
    !                       check_lag              	*
    !************************************************
    subroutine check_lag(poly, deg, b, c, z, r, conv)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        integer, intent(out)            :: conv
        real(kind=dp), intent(in)       :: r
        complex(kind=dp), intent(in)    :: poly(:), z
        complex(kind=dp), intent(out)   :: b, c
        ! local variables
        integer                         :: k
		real(kind=dp)					:: g, err
        complex(kind=dp)                :: a, a_err, b_err, c_err
        type(CEFTSum)                   :: sum
        type(CEFTProd)                  :: prod
        ! intrinsic functions
        intrinsic                       :: abs
		
        ! evaluate polynomial and derivatives
        a = poly(deg+1)
        a_err = cmplx(0.0_dp,0.0_dp,kind=dp)
        b = cmplx(0.0_dp,0.0_dp,kind=dp)
        b_err = cmplx(0.0_dp,0.0_dp,kind=dp)
        c = cmplx(0.0_dp,0.0_dp,kind=dp)
        c_err = cmplx(0.0_dp,0.0_dp,kind=dp)
		err = 0.0_dp
        do k=deg,1,-1
            ! product and sum for c
            prod = TwoProductCplx(c,z)
            sum = TwoSumCplx(prod%p,b)
            ! update c and c_err
            c = sum%x
            c_err = z*c_err + b_err + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            ! product and sum for b
            prod = TwoProductCplx(b,z)
            sum = TwoSumCplx(prod%p,a)
            ! update b and b_err
            b = sum%x
            b_err = z*b_err + a_err + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            ! product and sum for a
            prod = TwoProductCplx(a,z)
            sum = TwoSumCplx(prod%p,poly(k))
            ! update a and a_err
            a = sum%x
            a_err = z*a_err + FaithSumCplx(prod%e,prod%f,prod%g,sum%y)
            ! update berr
            err = r*err + FaithSum(abs(prod%e),abs(prod%f),abs(prod%g),abs(sum%y))
        end do
        ! add error back into result
        c = c + c_err
        b = b + b_err
        a = a + a_err
        ! berr
        g = (4*deg+2)*sqrt(2.0_dp)*eps/(1-eps)
        g = g/(1-g)
		err = mu*abs(a) + (g*err + 2*mu**2*abs(a))
		! laguerre correction terms/ convergence
		if(abs(a)>err) then
            b = b/a
            c = b**2 - 2*(c/a)
            if(check_nan_inf(b) .or. check_nan_inf(c)) conv = -1
		else
			conv = 1
		end if
    end subroutine check_lag
    !************************************************
    !                       modify_lag              *
    !************************************************
    subroutine modify_lag(deg, b, c, z, j, roots, conv)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, j
		integer, intent(out)            :: conv
        complex(kind=dp), intent(in)    :: roots(:), z
        complex(kind=dp), intent(inout) :: b, c
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: t
        ! intrinsic functions
        intrinsic                       :: abs, sqrt

        ! Aberth correction terms
        do k=1,j-1
            t = 1/(z - roots(k))
            b = b - t
            c = c - t**2
        end do
        do k=j+1,deg
            t = 1/(z - roots(k))
            b = b - t
            c = c - t**2
        end do
        ! Laguerre correction/ convergence
        t = sqrt((deg-1)*(deg*c-b**2))
        c = b + t
        b = b - t
        if(abs(b)>abs(c)) then
            c = deg/b
        else
            c = deg/c
        end if
		if(abs(c) <= mu*abs(z)) conv = 1
    end subroutine modify_lag
end module fpml_comp