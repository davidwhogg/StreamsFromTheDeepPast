RM=/bin/rm -vf
CP=/bin/cp -vf

all: lsr_1_lsr.ps lsr_1big_lsr.ps lsr_1fake_lsr.ps lsr_bestngauss_lsr.ps
	$(CP) *.ps *.tex ../../tex/lsr/

lsr_1_lsr.ps: lsr.pro
	echo lsr | idl

lsr_1big_lsr.ps: lsr.pro
	echo lsr,/bighalo | idl

lsr_1fake_lsr.ps: lsr.pro
	echo lsr,/fakedata | idl

lsr_bestngauss_lsr.ps: lsr.pro
	echo lsr,/bestngauss | idl

clean:
	$(RM) lsr*.ps lsr*.tex

spotless: clean
