#####################################################
#      CompRoots Make File Fortran TeX directory.   #
#####################################################
include ../make.inc

all: backErr-illCond backErr-wellCond forwErr \
	horner_aprioriErr horner_runErr limitAcc \
	newton_limitAcc special_poly clean

backErr-illCond:
	$(LC) backErr-illCond.tex
	
backErr-wellCond:
	$(LC) backErr-wellCond-err.tex && $(LC) backErr-wellCond-time.tex
	
forwErr:
	$(LC) forwErr-simple.tex
	
clean: copy
	@rm -f *.gz
	@rm -r *.log
	@rm -r *.aux
	@rm -r *.pdf

copy:
	@cp -r *.pdf ../figures/
	
horner_aprioriErr:
	$(LC) horner_aprioriErr_complex.tex && $(LC) horner_aprioriErr_real.tex

horner_runErr:
	$(LC) horner_runErr_complex.tex && $(LC) horner_runErr_real.tex
	
limitAcc:
	$(LC) limitAcc.tex

newton_limitAcc:
	$(LC) newton_limitAcc_real.tex
	
special_poly:
	$(LC) special_poly.tex