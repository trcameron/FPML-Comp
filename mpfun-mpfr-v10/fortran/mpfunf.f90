!*****************************************************************************

!  MPFUN-MPFR: An MPFR-based Fortran arbitrary precision computation package
!  Precision level declaration module (module MPFUNF)

!  Version date:  19 Apr 2017

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2017 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    This version differs from the MPFUN-Fort version by the same author in that
!    it is based on MPFR, which presently is the fastest available low-level
!    package for high-precision floating-point computation.  Thus most user
!    applications typically run 3X faster.  In addition, the developers of the
!    MPFR package have taken considerable pains to ensure that the many functions
!    return correctly rounded (to the last bit) results for each input.  At the
!    Fortran user level, application codes written for MPFUN-Fort may be compiled
!    and executed with MPFUN-MPFR --i.e., MPFUN-MPFR is "plug compatible" with
!    MPFUN-Fort.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNF):
!    This module defines the default precision level in digits (mpipl) and
!    the equivalent precision level in words (mpwds), which is calculated
!    below according to the formula:
!       mpwds = int (mpipl / mpdpw + 2)
!    where mpdpw = 19.26591972249... is an approximation to log10(2^64)
!    (mpdpw is set in module MPFUNA). This precision level is the working
!    precision level for all operations that use module MPFUNG, unless one
!    dynamically varies the medium precision level, in which case this is the
!    maximum working precision level.

module mpfunf
use mpfuna

implicit none
integer, public:: mpipl

!  *** Set the default precision level (in digits) here.

parameter (mpipl = 1200)

!----------------------------------------------------------------------------

!  *** Do not change the following code (in normal usage).

integer, public:: mpwds, mpwdsbt, mpwds6
parameter (mpwds = int (mpipl / mpdpw + 2.d0), mpwdsbt = mpnbt * mpwds, &
  mpwds6 = mpwds + 6)

end module mpfunf
