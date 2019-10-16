RM= /bin/rm -vf

targets= bovy_abstract.pdf bovy_abstract2.pdf

all: $(targets)

.SECONDARY: $(subst .pdf,.dvi,$(targets))

%.pdf: %.ps
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true $<
	$(RM) $<

%.gz: %
	gzip -f $<

%.ps: %.dvi
	dvips -t letter $< -o

%.dvi: %.tex Makefile
	latex $<
	latex $<
#	bibtex $(subst .tex,,$<)
	latex $<
	latex $<

.PHONY: clean spotless

clean:
	$(RM) *.aux *.bbl *.blg *.dvi *.log

spotless: clean
	$(RM) $(targets)
