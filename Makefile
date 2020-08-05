#####################################################
#      CompRoots Make File Fortran Root directory	#
#####################################################
include make.inc

backerr_illcond:
	$(FC) $(FFLAGS) -o backerr_illcond src/eft.f90 src/fpml.f90 src/fpml_comp.f90 src/fpml_quad.f90 src/mproutines.f90 src/backerr_illcond.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
backerr_wellcond:
	$(FC) $(FFLAGS) -o backerr_wellcond src/eft.f90 src/fpml.f90 src/fpml_comp.f90 src/fpml_quad.f90 src/mproutines.f90 src/backerr_wellcond.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
clean:
	@rm -f *.mod
	
compile:
	@$(MAKE) all -C TeX
	
forwerr_multi:
	$(FC) $(FFLAGS) -o forwerr_multi src/eft.f90 src/fpml.f90 src/fpml_comp.f90 src/fpml_quad.f90 src/mproutines.f90 src/forwerr_multi.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
forwerr_simple:
	$(FC) $(FFLAGS) -o forwerr_simple src/eft.f90 src/fpml.f90 src/fpml_comp.f90 src/fpml_quad.f90 src/mproutines.f90 src/forwerr_simple.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib

horner_runErr:
	$(FC) $(FFLAGS) -o horner_runErr src/eft.f90 src/mproutines.f90 src/horner_runErr.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
install: backerr_illcond backerr_wellcond forwerr_multi forwerr_simple horner_runErr laguerre_limitAcc newton_limitAcc special_poly
	
laguerre_limitAcc:
	$(FC) $(FFLAGS) -o laguerre_limitAcc src/eft.f90 src/mproutines.f90 src/laguerre_limitAcc.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
newton_limitAcc:
	$(FC) $(FFLAGS) -o newton_limitAcc src/eft.f90 src/mproutines.f90 src/newton_limitAcc.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
run:
	@echo "Backward Error Ill Conditioned ..."
	@./backerr_illcond
	@echo "Backward Error Well Conditioned ..."
	@./backerr_wellcond
	@echo "Forward Error Multi Roots ..."
	@./forwerr_multi
	@echo "Forward Error Simple Roots ..."
	@./forwerr_simple
	@echo "Horner Running Error ..."
	@./horner_runErr
	@echo "Laguerre Limiting Accuracy ..."
	@./laguerre_limitAcc
	@echo "Newton Limiting Accuracy ..."
	@./newton_limitAcc
	@echo "Special Polynomials"
	@./special_poly
	
special_poly:
	$(FC) $(FFLAGS) -o special_poly src/eft.f90 src/fpml.f90 src/fpml_comp.f90 src/fpml_quad.f90 src/mproutines.f90 src/special_poly.f90 $(MPOBJS) -I mpfun-mpfr-v10/fortran -lmpfr -L/usr/local/lib
	
uninstall: clean
	@rm -f backerr_illcond
	@rm -f backerr_wellcond
	@rm -f forwerr_multi
	@rm -f forwerr_simple
	@rm -f horner_runErr
	@rm -f laguerre_limitAcc
	@rm -f newton_limitAcc
	@rm -f special_poly