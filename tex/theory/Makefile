.SUFFIXES: .tex .dvi .ps .pdf .mp .eps .png
figpath= ../../pro/theory

all: stream_perturbation.pdf stream_direction.pdf

figs: $(figpath)/*.eps
#	cd $(figpath) ; make
	cp $(figpath)/*.eps .
	touch figs

%.ps: %.dvi
	dvips -t letter $< -o

%.dvi: %.tex figs
	latex $<
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "

%.pdf: %.dvi
	dvipdf $<
