all: rv.pdf

.secondary: rv.ps

%.pdf: %.ps
	ps2pdf $<

%.ps: %.dvi
	dvips -Ppdf $< -o

%.dvi: %.tex
	latex $<
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "
