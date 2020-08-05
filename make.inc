#############################################
#      CompRoots Make Include File          #
#############################################

# fortran compiler and flags
FC := gfortran
FFLAGS := -O2

# latex compiler and flags
LC := pdflatex

# mpfun objects
MPOBJS := mpfun-mpfr-v10/fortran/*.o