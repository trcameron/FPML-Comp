!********************************************************************************
!   FPML_QUAD: Fourth order Parallelizable Modification of Laguerre's method in Quadruple Precision
!   Author: Thomas R. Cameron, Davidson College
!   Last Modified: 1 November 2018
!********************************************************************************
module fpml_quad
    use eft, only : dp
    implicit none
    integer, parameter                  :: qp = selected_real_kind(2*precision(1.0_dp))
    real(kind=qp), parameter            :: eps = epsilon(1.0_qp), big = huge(1.0_qp), small = tiny(1.0_qp)
    
contains
    !************************************************
    !                       main quad               *
    !************************************************
    ! Computes the roots of a polynomial of degree
    ! deg whose coefficients are stored in p. The
    ! root approximations are stored in roots, the
    ! backward error in each approximation in berr,
    ! and the condition number of each root
    ! approximation is stored in cond.
    !************************************************
    subroutine main_quad(poly, deg, roots, berr, cond, conv, itmax)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, itmax
        integer, intent(out)            :: conv(:)
        real(kind=qp), intent(out)      :: berr(:), cond(:)
        complex(kind=qp), intent(in)    :: poly(:)
        complex(kind=qp), intent(out)   :: roots(:)
        ! local variables
        integer                         :: i, j, nz
        real(kind=qp)                   :: r
        real(kind=qp), dimension(deg+1) :: alpha
        complex(kind=qp)                :: b, c, z
        ! intrinsic functions
        intrinsic                       :: abs, sqrt

        ! precheck
        alpha = abs(poly)
        if(alpha(deg+1)<small) then
            write(*,'(A)') 'Warning: leading coefficient too small.'
            return
        elseif(deg==1) then
            roots(1) = -poly(1)/poly(2)
            conv = 1
            berr = 0.0_qp
            cond(1) = (alpha(1) + alpha(2)*abs(roots(1)))/(abs(roots(1))*alpha(2))
            return
        elseif(deg==2) then
            b = -poly(2)/(2*poly(3))
            c = sqrt(poly(2)**2-4*poly(3)*poly(1))/(2*poly(3))
            roots(1) = b - c
            roots(2) = b + c
            conv = 1
            berr = 0.0_qp
            cond(1) = (alpha(1)+alpha(2)*abs(roots(1))+alpha(3)*abs(roots(1))**2)/(abs(roots(1))* &
                        abs(poly(2)+2.0_qp*poly(3)*roots(1)))
            cond(2) = (alpha(1)+alpha(2)*abs(roots(2))+alpha(3)*abs(roots(2))**2)/(abs(roots(2))* &
                        abs(poly(2)+2.0_qp*poly(3)*roots(2)))
            return
        end if
        ! initial estimates
        conv = (/ (0, i=1,deg)/)
        nz = 0
        call estimates(alpha, deg, roots, conv, nz)
        ! main loop
        alpha = (/ (alpha(i)*(3.8_qp*(i-1)+1),i=1,deg+1)/)
        do i=1,itmax
            do j=1,deg
                if(conv(j)==0) then
                    z = roots(j)
                    r = abs(z)
                    if(r > 1.0_qp) then
                        call rcheck_lag(poly, alpha, deg, b, c, z, r, conv(j), berr(j), cond(j))
                    else
                        call check_lag(poly, alpha, deg, b, c, z, r, conv(j), berr(j), cond(j))
                    end if
                    if(conv(j)==0) then
                        call modify_lag(deg, b, c, z, j, roots)
                        roots(j) = roots(j) - c
                    else
                        nz = nz + 1
                        if(nz==deg) go to 10
                    end if
                end if
            end do
        end do
        ! final check
        10 continue
        if(minval(conv)==1) then
            return
        else
            ! display warrning
            write(*,'(A)') 'Some root approximations did not converge or experienced overflow/underflow.'
            ! compute backward error and condition number for roots that did not converge;
            ! note that this may produce overflow/underflow.
            do j=1,deg
                if(conv(j) .ne. 1) then
                    z = roots(j)
                    r = abs(z)
                    if(r>1.0_qp) then
                        z = 1/z
                        r = 1/r
                        c = 0
                        b = poly(1)
                        berr(j) = alpha(1)
                        do i=2,deg+1
                            c = z*c + b
                            b = z*b + poly(i)
                            berr(j) = r*berr(j) + alpha(i)
                        end do
                        cond(j) = berr(j)/abs(deg*b-z*c)
                        berr(j) = abs(b)/berr(j)
                    else
                        c = 0
                        b = poly(deg+1)
                        berr(j) = alpha(deg+1)
                        do i=deg,1,-1
                            c = z*c + b
                            b = z*b + poly(i)
                            berr(j) = r*berr(j) + alpha(i)
                        end do
                        cond(j) = berr(j)/(r*abs(c))
                        berr(j) = abs(b)/berr(j)
                    end if
                end if
            end do
        end if
    end subroutine main_quad
    !************************************************
    !                       rcheck_lag              *
    !************************************************
    ! Computes backward error of root approximation
    ! with moduli greater than 1.
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the Laguerre correction terms
    ! are computed and stored in variables b and c.
    !************************************************
    subroutine rcheck_lag(p, alpha, deg, b, c, z, r, conv, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        integer, intent(out)            :: conv
        real(kind=qp), intent(in)       :: alpha(:), r
        real(kind=qp), intent(out)      :: berr, cond
        complex(kind=qp), intent(in)    :: p(:), z
        complex(kind=qp), intent(out)   :: b, c
        ! local variables
        integer                         :: k
        real(kind=qp)                   :: rr
        complex(kind=qp)                :: a, zz
        ! intrinsic functions
        intrinsic                       :: abs

        ! evaluate polynomial and derivatives
        zz = 1/z
        rr = 1/r
        a = p(1)
        b = 0
        c = 0
        berr = alpha(1)
        do k=2,deg+1
            c = zz*c + b
            b = zz*b + a
            a = zz*a + p(k)
            berr = rr*berr + alpha(k)
        end do
        ! laguerre correction/ backward error and condition
        if(abs(a)>eps*berr) then
            b = b/a
            c = 2*(c/a)
            c = zz**2*(deg-2*zz*b+zz**2*(b**2-c))
            b = zz*(deg-zz*b)
            if(check_nan_inf(b) .or. check_nan_inf(c)) conv = -1
        else
            cond = berr/abs(deg*a-zz*b)
            berr = abs(a)/berr
            conv = 1
        end if
    end subroutine rcheck_lag
    !************************************************
    !                       check_lag               *
    !************************************************
    ! Computes backward error of root approximation
    ! with moduli less than or equal to 1.
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the Laguerre correction terms
    ! Gj and Hj are computed and stored in variables
    ! b and c, respectively.
    !************************************************
    subroutine check_lag(p, alpha, deg, b, c, z, r, conv, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        integer, intent(out)            :: conv
        real(kind=qp), intent(in)       :: alpha(:), r
        real(kind=qp), intent(out)      :: berr, cond
        complex(kind=qp), intent(in)    :: p(:), z
        complex(kind=qp), intent(out)   :: b, c
        ! local variables
        integer                         :: k
        complex(kind=qp)                :: a
        ! intrinsic functions
        intrinsic                       :: abs

        ! evaluate polynomial and derivatives
        a = p(deg+1)
        b = 0
        c = 0
        berr = alpha(deg+1)
        do k=deg,1,-1
            c = z*c + b
            b = z*b + a
            a = z*a + p(k)
            berr = r*berr + alpha(k)
        end do
        ! laguerre correction/ backward error and condition
        if(abs(a)>eps*berr) then
            b = b/a
            c = b**2 - 2*(c/a)
            if(check_nan_inf(b) .or. check_nan_inf(c)) conv = -1
        else
            cond = berr/(r*abs(b))
            berr = abs(a)/berr
            conv = 1
        end if
    end subroutine check_lag
    !************************************************
    !                       modify_lag              *
    !************************************************
    ! Computes modified Laguerre correction term of
    ! the jth rooot approximation.
    ! The coefficients of the polynomial of degree
    ! deg are stored in p, all root approximations
    ! are stored in roots. The values b, and c come
    ! from rcheck_lag or check_lag, c will be used
    ! to return the correction term.
    !************************************************
    subroutine modify_lag(deg, b, c, z, j, roots)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg, j
        complex(kind=qp), intent(in)    :: roots(:), z
        complex(kind=qp), intent(inout) :: b, c
        ! local variables
        integer                         :: k
        complex(kind=qp)                :: t
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
        ! Laguerre correction
        t = sqrt((deg-1)*(deg*c-b**2))
        c = b + t
        b = b - t
        if(abs(b)>abs(c)) then
            c = deg/b
        else
            c = deg/c
        end if
    end subroutine modify_lag
    !************************************************
    !                       estimates               *
    !************************************************
    ! Computes initial estimates for the roots of an
    ! univariate polynomial of degree deg, whose
    ! coefficients moduli are stored in alpha. The
    ! estimates are returned in the array roots.
    ! The computation is performed as follows: First
    ! the set (i,log(alpha(i))) is formed and the
    ! upper envelope of the convex hull of this set
    ! is computed, its indices are returned in the
    ! array h (in descending order). For i=c-1,1,-1
    ! there are h(i) - h(i+1) zeros placed on a
    ! circle of radius alpha(h(i+1))/alpha(h(i))
    ! raised to the 1/(h(i)-h(i+1)) power.
    !************************************************
    subroutine estimates(alpha, deg, roots, conv, nz)
        implicit none
        ! argument variables
        integer, intent(in)             :: deg
        integer, intent(inout)          :: conv(:), nz
        real(kind=qp), intent(in)       :: alpha(:)
        complex(kind=qp), intent(inout) :: roots(:)
        ! local variables
        integer                         :: c, i, j, k, nzeros
        real(kind=qp)                   :: a1, a2, ang, r, th
        integer, dimension(deg+1)       :: h
        real(kind=qp), dimension(deg+1) :: a
        real(kind=qp), parameter        :: pi2 = 6.2831853071795865_qp, sigma = 0.7_qp
        ! intrinsic functions
        intrinsic                       :: log, cos, sin, cmplx

        ! Log of absolute value of coefficients
        do i=1,deg+1
            if(alpha(i)>0) then
                a(i) = log(alpha(i))
            else
                a(i) = -1E+30_qp
            end if
        end do
        call conv_hull(deg+1, a, h, c)
        k=0
        th=pi2/deg
        ! Initial Estiamtes
        do i=c-1,1,-1
            nzeros = h(i)-h(i+1)
            a1 = alpha(h(i+1))**(1.0_qp/nzeros)
            a2 = alpha(h(i))**(1.0_qp/nzeros)
            if(a1 .le. a2*small) then
                ! r is too small
                r = 0.0_qp
                nz = nz + nzeros
                conv(k+1:k+nzeros) = -1
                roots(k+1:k+nzeros) = cmplx(0,0,kind=qp)
            else if(a1 .ge. a2*big) then
                ! r is too big
                r = big
                nz = nz+nzeros
                conv(k+1:k+nzeros) = -1
                ang = pi2/nzeros
                do j=1,nzeros
                    roots(k+j) = r*cmplx(cos(ang*j+th*h(i)+sigma),sin(ang*j+th*h(i)+sigma),kind=qp)
                end do
            else
                ! r is just right
                r = a1/a2
                ang = pi2/nzeros
                do j=1,nzeros
                    roots(k+j) = r*cmplx(cos(ang*j+th*h(i)+sigma),sin(ang*j+th*h(i)+sigma),kind=qp)
                end do
            end if
            k = k+nzeros
        end do
    end subroutine estimates
    !************************************************
    !                       conv_hull               *
    !************************************************
    ! Computex upper envelope of the convex hull of
    ! the points in the array a, which has size n.
    ! The number of vertices in the hull is equal to
    ! c, and they are returned in the first c entries
    ! of the array h.
    ! The computation follows Andrew's monotone chain
    ! algorithm: Each consecutive three pairs are
    ! tested via cross to determine if they form
    ! a clockwise angle, if so that current point
    ! is rejected from the returned set.
    !************************************************
    subroutine conv_hull(n, a, h, c)
        implicit none
        ! argument variables
        integer, intent(in)         :: n
        integer, intent(inout)      :: c
        integer, intent(inout)      :: h(:)
        real(kind=qp), intent(in)   :: a(:)
        ! local variables
        integer                     :: i

        ! covex hull 
        c=0
        do i=n,1,-1
            do while(c>=2 .and. cross(h, a, c, i)<eps)
                c = c - 1
            end do
            c = c + 1
            h(c) = i
        end do
    end subroutine conv_hull
    !************************************************
    !                       cross                   *
    !************************************************
    ! Returns 2D cross product of OA and OB vectors,
    ! where
    ! O=(h(c-1),a(h(c-1))),
    ! A=(h(c),a(h(c))),
    ! B=(i,a(i)).
    ! If det>0, then OAB makes counter-clockwise turn.
    !************************************************
    function cross(h, a, c, i) result(det)
        implicit none
        ! argument variables
        integer, intent(in)         :: c, i
        integer, intent(in)         :: h(:)
        real(kind=qp), intent(in)   :: a(:)
        ! local variables
        real(kind=qp)               :: det

        ! determinant
        det = (a(i)-a(h(c-1)))*(h(c)-h(c-1)) - (a(h(c))-a(h(c-1)))*(i-h(c-1))
        return
    end function cross
    !************************************************
    !                       check_nan_inf           *
    !************************************************
    ! Check if real or imaginary part of complex
    ! number a is either NaN or Inf.
    !************************************************
    function check_nan_inf(a) result(res)
        implicit none
        ! argument variables
        complex(kind=qp)            :: a
        ! local variables
        logical                     :: res
        real(kind=qp)               :: re_a, im_a
        ! intrinsic functions
        intrinsic                   :: abs

        ! check for nan and inf
        re_a = real(a,kind=qp)
        im_a = aimag(a)
        res = isnan(re_a) .or. isnan(im_a) .or. (abs(re_a)>big) .or. (abs(im_a)>big)
        return
    end function check_nan_inf
end module fpml_quad