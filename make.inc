#############################################
#      CompRoots Make Include File          #
#############################################

# fortran compiler and flags
FC := gfortran
FFLAGS := -O2

# latex compiler and flags
LC := pdflatex

# mpfun objects
MPOBJS := mpfun/mpfun90.o mpfun/mpmod90.o mpfun/mpmodm90.o mpfun/mpmodx90.o

# nag info
NAGDIR := /Users/thcameron/Documents/Research/c02aff
NAGINC := $(NAGDIR)/include
NAGOBJ := $(NAGDIR)/objects
FCOMPILE := gfortran -O2 -I $(NAGINC) $(NAGOBJ)/*.o