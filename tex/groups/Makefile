RM=/bin/rm -vf
.SUFFIXES: .dvi .pdf .eps
.SECONDARY: clusterhist_ngc1901.eps clusterhist_sirius.eps clusterhist_pleiades.eps clusterhist_hyades.eps hyadesClusterHist.eps model_selection.eps fit_ssp_fixpij_fixalpha_ngc1901.eps fit_ssp_fixpij_fixalpha_sirius.eps fit_ssp_fixpij_fixalpha_pleiades.eps fit_ssp_fixpij_fixalpha_hyades.eps fit_ssp_fixpij_fixalpha_hercules.eps fit_ssp_fixpij_pleiades.eps fit_ssp_fixpij_hyades.eps fit_ssp_fixpij_hercules.eps fit_ssp_fixpij_pleiades_alpha.eps fit_ssp_fixpij_hyades_alpha.eps fit_ssp_fixpij_hercules_alpha.eps fit_ssp_fixpij_ngc1901.eps fit_ssp_fixpij_ngc1901_alpha.eps fit_ssp_fixpij_sirius.eps fit_ssp_fixpij_sirius_alpha.eps predict_plx_foreground_sirius.eps veldist.eps

all: ms.pdf

%.pdf: %.dvi
	dvips -t letter $< -o
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true $(subst .dvi,.ps,$<)

%.gz: %
	gzip -f $<

%.dvi: %.tex clusterhist_ngc1901.eps clusterhist_sirius.eps clusterhist_pleiades.eps clusterhist_hyades.eps hyadesClusterHist.eps model_selection.eps fit_ssp_fixpij_fixalpha_ngc1901.eps fit_ssp_fixpij_fixalpha_sirius.eps fit_ssp_fixpij_fixalpha_pleiades.eps fit_ssp_fixpij_fixalpha_hyades.eps fit_ssp_fixpij_fixalpha_hercules.eps fit_ssp_fixpij_pleiades.eps fit_ssp_fixpij_hyades.eps fit_ssp_fixpij_hercules.eps fit_ssp_fixpij_pleiades_alpha.eps fit_ssp_fixpij_hyades_alpha.eps fit_ssp_fixpij_hercules_alpha.eps fit_ssp_fixpij_ngc1901.eps fit_ssp_fixpij_ngc1901_alpha.eps fit_ssp_fixpij_sirius.eps fit_ssp_fixpij_sirius_alpha.eps predict_plx_foreground_sirius.eps veldist.eps
	latex $<
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "

model_selection.eps: model_selection.ps
	cp $< tmp.ps
	echo "1,\$$s/%%BoundingBox: 36 126 576 666/%%BoundingBox: 88 394 525 616/g" > edcmd
	echo "w" >> edcmd
	ed $< < edcmd
	cp $< $@
	cp tmp.ps $< 
	rm edcmd tmp.ps

veldist.eps: veldist.ps
	cp $< tmp.ps
	echo "1,\$$s/%%BoundingBox: 27 195 585 596/%%BoundingBox: 27 360 311 561/g" > edcmd
	echo "w" >> edcmd
	ed $< < edcmd
	cp $< $@
	cp tmp.ps $< 
	rm edcmd tmp.ps

predict_plx_foreground_sirius.eps: predict_plx_foreground_sirius.ps
	cp $< tmp.ps
	echo "1,\$$s/%%BoundingBox: 36 126 576 666/%%BoundingBox: 55 318 526 628/g" > edcmd
	echo "w" >> edcmd
	ed $< < edcmd
	cp $< $@
	cp tmp.ps $< 
	rm edcmd tmp.ps

fit_ssp_fixpij_%_alpha.eps: fit_ssp_fixpij_%_alpha.ps
	cp $< tmp.ps
	echo "1,\$$s/%%BoundingBox: 36 126 576 666/%%BoundingBox: 73 172 484 596/g" > edcmd
	echo "w" >> edcmd
	ed $< < edcmd
	cp $< $@
	cp tmp.ps $< 
	rm edcmd tmp.ps

fit_ssp_fixpij_%.eps: fit_ssp_fixpij_%.ps
	cp $< tmp.ps
	echo "1,\$$s/%%BoundingBox: 36 126 576 666/%%BoundingBox: 56 170 481 600/g" > edcmd
	echo "w" >> edcmd
	ed $< < edcmd
	cp $< $@
	cp tmp.ps $< 
	rm edcmd tmp.ps

%.eps: %.ps
	cp $< tmp.ps
	echo "1,\$$s/%%BoundingBox: 36 126 576 666/%%BoundingBox: 46 340 563 621/g" > edcmd
	echo "w" >> edcmd
	ed $< < edcmd
	cp $< $@
	cp tmp.ps $< 
	rm edcmd tmp.ps

%.dvi: %.tex

.PHONY: clean spotless

clean:
	$(RM) *.dvi *.aux

spotless: clean
	$(RM) *.pdf *.eps