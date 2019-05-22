!********************************************************************************
!   FPML_COMP: Compensated Polishing Technique for FPML based on Newton's method,
!   The Ehrlich-Aberth method, and the modified Laguerre's method.
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 24 April 2019
!********************************************************************************
module fpml_comp
    use fpml
    implicit none
    
contains
    !************************************************
    !                       NPolish                 *
    !************************************************
    ! Single iteration of Newton's method in the 
    ! original polynomial. 
    !************************************************
    subroutine NPolish(poly, deg, roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        logical                         :: conv(deg)
        integer                         :: j, k
        integer, parameter              :: itmax = 30
        real(kind=dp)                   :: berr(deg), cond(deg)
        complex(kind=dp)                :: a, b, z
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! call fpml main
        call main(poly, deg, roots, berr, cond, conv, itmax)
        ! polish root approximations from fpml
        do j=1,deg
            z = roots(j)
            b = 0
            if(abs(z)>1) then
                z = 1/z
                a = poly(1)
                do k=2,deg+1
                    b = z*b + a
                    a = z*a + poly(k)
                end do
                a = a/b
                a = a/(z*(deg*a - z))
            else
                a = poly(deg+1)
                do k=deg,1,-1
                    b = z*b + a
                    a = z*a + poly(k)
                end do
                a = a/b
            end if
            roots(j) = roots(j) - a
        end do
    end subroutine NPolish
    !************************************************
    !           MLPolish                            *
    !************************************************
    ! Iterative refinement of root approximations
    ! using compensated Horner's method and the
    ! modified Laguerre method to avoid unecessary
    ! multiple convergence to same root.
    !************************************************
    subroutine MLPolish(poly, deg, roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        complex(kind=dp), intent(in)    :: poly(:)
        complex(kind=dp), intent(inout) :: roots(:)
        ! local variables
        logical                         :: conv(deg)
        integer                         :: i, j, k, nz
        integer, parameter              :: itnum = 17
        real(kind=dp)                   :: errBound, g, r
        complex(kind=dp)                :: a, b, compPoly, compDer, compDer2, errPoly, errDer, errDer2, z
        type(CEFTSum)                   :: sum
        type(CEFTProd)                  :: prod
        ! fpml variables
        integer, parameter              :: itmax = 30
        real(kind=dp)                   :: berr(deg), cond(deg)
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! call fpml main
        call main(poly, deg, roots, berr, cond, conv, itmax)
        ! polish root approximations from fpml
        conv = .false.
        nz = 0
        do i=1,itnum
            do j=1,deg
                if(.not.conv(j)) then
                    z = roots(j)
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
                    ! update roots(j) and conv(j)
                    if(2*errBound<abs(compPoly)) then
                        ! compute modified Laguerre correction terms
                        do k=1,j-1
                            ! compute sum term
                            z = 1/(roots(j) - roots(k))
                            a = a - z
                            b = b - z**2
                        end do
                        do k=j+1,deg
                            ! compute sum term
                            z = 1/(roots(j) - roots(k))
                            a = a - z
                            b = b - z**2
                        end do
                        ! update roots(j)
                        z = sqrt((deg-1)*(deg*b-a**2))
                        b = a + z
                        a = a - z
                        if(abs(a)>abs(b)) then
                            z = deg/a
                        else
                            z = deg/b 
                        end if
                        if(abs(z)>2*mu*abs(roots(j))) then
                            roots(j) = roots(j) - z
                        else
                            ! update conv(j) and nz
                            conv(j) = .true.
                            nz = nz + 1
                        end if
                    else
                        ! update conv(j) and nz
                        conv(j) = .true.
                        nz = nz + 1
                    end if
                    if(nz==deg) return
                end if
            end do
        end do
    end subroutine MLPolish
end module fpml_comp