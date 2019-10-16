RM=/bin/rm -fv

all: paper1.ps.gz

paper1.ps.gz: paper1.dvi
	dvips paper1 -o
	gzip -fv paper1.ps

paper1.dvi: paper1.tex paper1.bbl
	latex paper1
	latex paper1
	latex paper1

clean:
	$(RM) *.aux *.dvi *.log *.blg *.bbl

spotless: clean
	$(RM) paper1.ps.gz
