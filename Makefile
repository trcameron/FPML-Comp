#####################################################
#      CompRoots Make File Fortran Root directory	#
#####################################################
include make.inc
	
clean:
	@rm -f *.mod

compile:
	@$(MAKE) all -C TeX
	
backerr_illcond: mpfunsrc
	$(FC) $(FFLAGS) -o backerr_illcond src/eft.f90 src/mproutines.f90 src/fpml.f90 src/fpml_comp.f90 src/backerr_illcond.f90 $(MPOBJS) -I mpfun
	
backerr_wellcond: mpfunsrc
	$(FC) $(FFLAGS) -o backerr_wellcond src/eft.f90 src/mproutines.f90 src/fpml.f90 src/fpml_comp.f90 src/backerr_wellcond.f90 $(MPOBJS) -I mpfun

#forwerr_multi: mpfunsrc nag
#	$(FCOMPILE) $(FFLAGS) -o forwerr_multi src/eft.f90 src/mproutines.f90 src/fpml.f90 src/fpml_comp.f90 src/forwerr_multi.f90 $(MPOBJS) -I mpfun
	
#forwerr_simple: mpfunsrc nag
#	$(FCOMPILE) $(FFLAGS) -o forwerr_simple src/eft.f90 src/mproutines.f90 src/fpml.f90 src/fpml_comp.f90 src/forwerr_simple.f90 $(MPOBJS) -I mpfun

horner_aprioriErr: mpfunsrc
	$(FC) $(FFLAGS) -o horner_aprioriErr src/eft.f90 src/mproutines.f90 src/horner_aprioriErr.f90 $(MPOBJS) -I mpfun
	
horner_runErr: mpfunsrc
	$(FC) $(FFLAGS) -o horner_runErr src/eft.f90 src/mproutines.f90 src/horner_runErr.f90 $(MPOBJS) -I mpfun
	
install: backerr_illcond backerr_wellcond horner_aprioriErr horner_runErr lag_limitAcc newton_limitAcc
	
lag_limitAcc: mpfunsrc
	$(FC) $(FFLAGS) -o lag_limitAcc src/eft.f90 src/mproutines.f90 src/lag_limitAcc.f90 $(MPOBJS) -I mpfun

#nag:
#	@$(MAKE) install -C $(NAGDIR)

newton_limitAcc: mpfunsrc
	$(FC) $(FFLAGS) -o newton_limitAcc src/eft.f90 src/mproutines.f90 src/newton_limitAcc.f90 $(MPOBJS) -I mpfun
	
#special_poly: mpfunsrc nag
#	$(FCOMPILE) $(FFLAGS) -o special_poly src/eft.f90 src/mproutines.f90 src/fpml.f90 src/fpml_comp.f90 src/special_poly.f90 $(MPOBJS) -I mpfun
	
mpfunsrc:
	@$(MAKE) all -C mpfun
	
run:
	@echo "Backward Error Ill-Conditioned Test ..."
	@./backerr_illcond
	@echo "Backward Error Well-Conditioned Test ..."
	@./backerr_wellcond
	#@echo "Forward Error Multi Test ..."
	#@./forwerr_multi
	#@echo "Forward Error Simple Test ..."
	#@./forwerr_simple
	@echo "Horner A Priori Error Bound Test ..."
	@./horner_aprioriErr
	@echo "Horner Running Error Bound Test ..."
	@./horner_runErr
	@echo "Laguerre Limiting Accuracy Test ..."
	@./lag_limitAcc
	@echo "Newton Limiting Accuracy Test ..."
	@./newton_limitAcc
	#@echo "Special Polynomials Test ..."
	#@./special_poly
	
uninstall: clean
	@$(MAKE) clean -C mpfun
	#@$(MAKE) clean -C $(NAGDIR)
	@rm -f backerr_illcond
	@rm -f backerr_wellcond
	#@rm -f forwerr_multi
	#@rm -f forwerr_simple
	@rm -f horner_aprioriErr
	@rm -f horner_runErr
	@rm -f lag_limitAcc
	@rm -f lag_time
	@rm -f newton_limitAcc
	#@rm -f special_poly