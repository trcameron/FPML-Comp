#####################################################
#      CompRoots Make File Fortran TeX directory.   #
#####################################################
include ../make.inc

all: backerr_illcond backerr_wellcond horner_runErr laguerre_limitAcc newton_limitAcc clean
	
backerr_illcond:
	$(LC) backerr_illcond.tex

backerr_wellcond:
	$(LC) backerr_wellcond_err.tex && $(LC) backerr_wellcond_time.tex
	
clean: copy
	@rm -f *.gz
	@rm -r *.log
	@rm -r *.aux
	@rm -r *.pdf

copy:
	@cp -r *.pdf ../figures/
	
forwerr_simple:
	$(LC) forwerr_simple.tex

horner_runErr:
	$(LC) horner_runErr_real.tex && $(LC) horner_runErr_cmplx.tex
	
laguerre_limitAcc:
	$(LC) laguerre_limitAcc_real.tex && $(LC) laguerre_limitAcc_cmplx.tex
	
limitAcc:
	$(LC) limitAcc.tex

newton_limitAcc:
	$(LC) newton_limitAcc_real.tex && $(LC) newton_limitAcc_cmplx.tex