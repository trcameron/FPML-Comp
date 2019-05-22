!********************************************************************************
!   EFT: Error-Free Transformations and Compensated Arithmetic Routines
!   Author: Thomas R. Cameron and Aidan O'Neill
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 21 March 2019
!********************************************************************************
module eft
    implicit none
    integer, parameter                          :: dp = kind(1.d0), factor = 2**27 + 1
    real(kind=dp), parameter                    :: mu = 2.0_dp**(-53), eps = 2.0_dp**(-52)
    !********************************************************
    !   REFT: Real Error Free Transformation type.
    !********************************************************
    !   Made up of two real(kind=dp) variables x and y
    !   called primary and secondary, respectively. 
    !   For TwoSum and TwoProd, x will denote floating-
    !   point result and y will denote error in result.
    !   For Split, x and y are the two parts of the
    !   floating-point number. 
    !********************************************************                 
    type REFT
        real(kind=dp)                           :: x        ! primary
        real(kind=dp)                           :: y        ! secondary
    end type REFT
    !********************************************************
    !   CEFTSum: Complex Error Free Transformation Sum.
    !********************************************************
    !   Made up of two complex(kind=dp) variables x and y
    !   called primary and secondary, respectively. 
    !   For TwoSumCplx, x denotes the floating-point 
    !   result and y denotes the error.
    !********************************************************
    type CEFTSum
        complex(kind=dp)                        :: x        ! primary
        complex(kind=dp)                        :: y        ! secondary
    end type CEFTSum
    !********************************************************
    !   CEFTProd: Complex Error Free Transformation Product.
    !********************************************************
    !   Made up of four complex(kind=dp) variables p, e, f,
    !   and g, which we call primary, secondary, tertiary,
    !   and quaternary, respectively. For TwoProductCplx,
    !   p denotes the floating-point result, and e, f, and g
    !   denote errors associated with the sum and product
    !   of real and imaginary parts.                        
    !********************************************************
    type CEFTProd
        complex(kind=dp)                        :: p        ! primary
        complex(kind=dp)                        :: e        ! secondary
        complex(kind=dp)                        :: f        ! tertiary
        complex(kind=dp)                        :: g        ! quaternary
    end type CEFTProd
    !********************************************************
    !   REFTHorner: Real Error Free Transformation Horner.
    !********************************************************
    !   Made up of a real(kind=dp) variable h to store result
    !   of standard Horner method and two real(kind=dp) 
    !   allocatable variables to store polynomial coefficients
    !   p and q, where pi is the error in the product term
    !   and qi is the error in the sum term on the ith
    !   iteration of Horner's method.                       
    !********************************************************
    type REFTHorner
        real(kind=dp)                           :: h        ! result
        real(kind=dp), allocatable              :: p(:)     ! error in product
        real(kind=dp), allocatable              :: q(:)     ! error in sum
    end type REFTHorner
    !********************************************************
    !   CEFTHorner: Complex Error Free Transformation Horner.
    !********************************************************
    !   Made up of a complex(kind=dp) variable h to store 
    !   result of standard Horner method and four allocatable
    !   complex(kind=dp) variables to store polynomial 
    !   coefficients p, q, r, and s, where pi, qi, ri is the
    !   error in the product term and si is the error in the
    !   sum term on the ith iteration of Horner's method.   
    !********************************************************
    type CEFTHorner
        complex(kind=dp)                        :: h        ! result
        complex(kind=dp), allocatable           :: p(:)     ! error in product
        complex(kind=dp), allocatable           :: q(:)     ! error in product
        complex(kind=dp), allocatable           :: r(:)     ! error in product
        complex(kind=dp), allocatable           :: s(:)     ! error in sum
    end type CEFTHorner

contains

    !********************************************************
    !                       TwoSum                          *
    !********************************************************
    ! Computes the sum of two floating-point numbers
    ! a and b. The result and error is returned in
    ! comp which has type REFT. 
    !********************************************************
    function TwoSum(a,b) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)       :: a, b
        ! local variables
        real(kind=dp)       :: z
        type(REFT)          :: comp

        ! compute floating-point sum
        comp%x = a + b
        z = comp%x - a
        ! compute error in floating-point sum
        comp%y = (a - (comp%x - z)) + (b - z)
        return
    end function TwoSum
    !********************************************************
    !                       TwoSumCplx                      *
    !********************************************************
    ! Computes the sum of two complex floating-point 
    ! numbers a and b. The result and error is 
    ! returned in comp which has type CEFTSum. 
    !********************************************************
    function TwoSumCplx(a,b) result(comp)
        implicit none
        ! argument variables
        complex(kind=dp)        :: a, b
        ! local variables
        type(REFT)              :: rcomp, icomp
        type(CEFTSum)           :: comp
            
        ! compute sum of real and imaginary parts
        rcomp = TwoSum(real(a),real(b))
        icomp = TwoSum(aimag(a),aimag(b))
        ! store primary and secondary complex numbers in comp
        comp%x = cmplx(rcomp%x,icomp%x,kind=dp)
        comp%y = cmplx(rcomp%y,icomp%y,kind=dp)
        return
    end function TwoSumCplx
    !********************************************************
    !                       Split                           *
    !********************************************************
    ! Computes the splitting of a floating-point 
    ! number a. Both parts have at most 26 nonzero bits 
    ! and are returned in comp which has type REFT.
    !********************************************************
    function Split(a) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)       :: a
        ! local variables   
        real(kind=dp)       :: c
        type(REFT)          :: comp
        
        ! factor a
        c = factor * a
        ! split into two 26-bit numbers
        comp%x = c - (c - a)
        comp%y = a - comp%x
        return
    end function Split
    !********************************************************
    !                       TwoProduct                      *
    !********************************************************
    ! Computes the product of two floating-point
    ! numbers a and b. The result and error is 
    ! returned in comp which has type REFT. 
    !********************************************************
    function TwoProduct(a,b) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)       :: a, b
        ! local variables
        type(REFT)          :: comp, sa, sb
        
        ! compute floating-point product
        comp%x = a * b
        ! split a and b
        sa = Split(a)
        sb = Split(b)
        ! compute error in floating-point product
        comp%y = sa%y * sb%y - (((comp%x - sa%x * sb%x) - sa%y * sb%x) - sa%x * sb%y)
        return
    end function TwoProduct
    !********************************************************
    !                       TwoProductCplx                  *
    !********************************************************
    ! Computes the product of two complex floating-point 
    ! numbers a and b. The result and error is returned 
    ! in comp which has type CEFTP. 
    !********************************************************
    function TwoProductCplx(a,b) result(comp)
        implicit none
        ! argument variables
        complex(kind=dp)    :: a, b
        ! local variables
        real(kind=dp)       :: z1, z2, z3, z4
        type(REFT)          :: sra, sia, srb, sib, sumr, sumi
        type(CEFTProd)      :: comp
            
        ! perform splitting of real(a), imag(a), real(b), imag(b)
        sra = Split(real(a)); sia = Split(aimag(a))
        srb = Split(real(b)); sib = Split(aimag(b))
        ! compute product of real and imaginary parts
        z1 = real(a) * real(b); z2 = aimag(a) * aimag(b)
        z3 = real(a) * aimag(b); z4 = aimag(a) * real(b)
        ! error in z1 and z3
        comp%e = cmplx(sra%y * srb%y - (((z1 - sra%x * srb%x) - sra%y * srb%x) - sra%x * srb%y),&
                    sra%y * sib%y - (((z3 - sra%x * sib%x) - sra%y * sib%x) - sra%x * sib%y),kind=dp)
        ! error in z2 and z4
        comp%f = cmplx(-(sia%y * sib%y - (((z2 - sia%x * sib%x) - sia%y * sib%x) - sia%x * sib%y)),&
                    sia%y * srb%y - (((z4 - sia%x * srb%x) - sia%y * srb%x) - sia%x * srb%y),kind=dp)
        ! compute sum z1-z2
        sumr = TwoSum(z1,-z2)
        ! compute sum z3+z4
        sumi = TwoSum(z3,z4)
        ! store floating-point result
        comp%p = cmplx(sumr%x,sumi%x,kind=dp)
        ! error in sums
        comp%g = cmplx(sumr%y,sumi%y,kind=dp)
        return
    end function TwoProductCplx
    !********************************************************
    !                       EFT Horner                      *
    !********************************************************
    ! Compute the evaluation of a polynomial at a
    ! floating-point number. The result and error
    ! is returned in comp which has type REFTHorner. 
    !********************************************************
    function EFTHorner(poly,x,deg) result(comp)
        implicit none
        ! argument variables
        integer             :: deg
        real(kind=dp)       :: poly(:), x
        ! local variables
        integer             :: k
        type(REFT)          :: prod, sum
        type(REFTHorner)    :: comp
            
        ! allocate memory for comp polynomials
        allocate(comp%p(deg),comp%q(deg))
        ! Horner's method
        comp%h = poly(deg+1)
        do k=deg,1,-1
            ! product and sum
            prod = TwoProduct(comp%h,x)
            sum = TwoSum(prod%x,poly(k))
            ! update comp
            comp%h = sum%x
            comp%p(k) = prod%y
            comp%q(k) = sum%y
        end do
        return
    end function EFTHorner
    !********************************************************
    !                       EFT Horner Cplx                 *
    !********************************************************
    ! Compute the evaluation of a polynomial at a complex
    ! floating-point number. The result and error is returned
    ! in comp which has type CEFTHorner.        
    !********************************************************
    function EFTHornerCplx(poly,x,deg)  result(comp)
        implicit none
        ! argument variables
        integer             :: deg
        complex(kind=dp)    :: poly(:), x
        ! local variables
        integer             :: k
        type(CEFTSum)       :: sum
        type(CEFTProd)      :: prod
        type(CEFTHorner)    :: comp
            
        ! allocate memory for comp polynomials
        allocate(comp%p(deg),comp%q(deg),comp%r(deg),comp%s(deg))
        ! Horner's method
        comp%h = poly(deg+1)
        do k=deg,1,-1
            ! product and sum
            prod = TwoProductCplx(comp%h,x)
            sum = TwoSumCplx(prod%p,poly(k))
            ! update comp
            comp%h = sum%x
            comp%p(k) = prod%e
            comp%q(k) = prod%f
            comp%r(k) = prod%g
            comp%s(k) = sum%y
        end do
        return
    end function EFTHornerCplx
    !****************************************************
    !                       NextPowerTwo                *
    !****************************************************
    ! Computation of 2^ceil(log2(abs(p))), for p \neq 0.          
    !****************************************************
    function NextPowerTwo(p) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)               :: p
        ! local variables
        real(kind=dp)               :: comp, q
        
        q = p/mu
        comp = abs((q+p)-q)
        if(comp==0.0_dp) then
            comp = abs(p)
        end if
        return
    end function NextPowerTwo
    !****************************************************
    !                       FaithSum                    *
    !****************************************************
    ! Computes the faithful sum of the real floating-
    ! point numbers a, b, c, and d. Result is returned
    ! in comp which has type real(kind=dp).          
    !****************************************************
    function FaithSum(a,b,c,d) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)                   :: a, b, c, d
        ! local variables
        integer                         :: n
        real(kind=dp)                   :: fac, m, Ms, phi, sigma
        real(kind=dp)                   :: comp, t, tau, tau1, tau2
        real(kind=dp), allocatable      :: p(:), q(:)
        real(kind=dp), parameter        :: realmin = tiny(1.0_dp)
            
        ! initialize arrays
        n = 4
        allocate(p(n),q(n))
        p(1)=a; p(2)=b; p(3)=c; p(4)=d
        ! check size and max
        10 m = maxval(abs(p))
        if(n==0 .or. m==0.0_dp) then
            comp = 0.0_dp
            deallocate(p,q)
            return
        end if
        ! initialize variables
        Ms = NextPowerTwo((n+2)*1.0_dp)
        sigma = Ms*NextPowerTwo(m)
        phi = mu*Ms
        fac = eps*Ms*Ms
        ! main loop
        t = 0
        do while(.True.)
            q = (sigma + p) - sigma
            tau = sum(q)
            p = p - q
            tau1 = t + tau
            ! check new approximation tau1
            if(abs(tau1)>=fac*sigma .or. sigma<=realmin) then
                tau2 = tau - (tau1 - t)
                comp = tau1 + (tau2 + sum(p))
                deallocate(p,q)
                return
            end if
            t = tau1
            ! accelerate zero sum case
            if(t==0.0_dp) then
                p = pack(p, p /= 0.0_dp)
                n = size(p)
                deallocate(q)
                allocate(q(n))
                go to 10
            end if
            ! new extraction unit
            sigma = phi*sigma
        end do
    end function FaithSum
    !****************************************************
    !                       FaithSumCplx                *
    !****************************************************
    ! Computes the faithful sum of the complex floating-
    ! point numbers a, b, c, and d. Result is returned
    ! in comp which has type complex(kind=dp).          
    !****************************************************
    function FaithSumCplx(a,b,c,d) result(comp)
        implicit none
        ! argument variables
        complex(kind=dp)                :: a, b, c, d
        ! local variables
        complex(kind=dp)                :: comp
        
        comp = cmplx(FaithSum(real(a),real(b),real(c),real(d)),&
                    FaithSum(aimag(a),aimag(b),aimag(c),aimag(d)),kind=dp)
    end function FaithSumCplx
end module eft