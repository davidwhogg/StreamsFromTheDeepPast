RM= /bin/rm -vf

targets= fitting_gaussians.pdf ms.pdf shadowms.pdf ms2.pdf

figures= hipparcos2_properties.ps hipparcos2_err_properties.ps hipvship2.ps hipvstyc.ps veldensXY.ps veldensXZ.ps veldensYZ.ps loglike.ps crossval.ps aic.ps partition.ps mdl.ps mml.ps bayes.ps gcs_cross.ps predict_gcs_rand.ps predict_gcs_good.ps predict_gcs_bad.ps hist_gcs_like.ps predict_gcs_high.ps predict_gcs_low.ps hist_gcs_ent.ps checkquants.ps radecpatches.ps tile_onepred.ps info_hip_high.ps info_hip_low.ps hist_hip_ent.ps info_entropies_low.ps info_entropies_high.ps 

all: $(targets)

.SECONDARY: $(subst .pdf,.dvi,$(targets))

%.pdf: %.ps
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true $<
	$(RM) $<

%.gz: %
	gzip -f $<

%.ps: %.dvi
	dvips -t letter $< -o

%.dvi: %.tex $(figures) Makefile EM.bib bovy_streams.bib MML.bib
	latex $<
	latex $<
	bibtex $(subst .tex,,$<)
	latex $<
	latex $<
	latex $<

$(figures):
	touch $(figures)

ms.tex: paperI.tex paperI_figures.tex parameters.txt high_info_hip.txt low_info_hip.txt
	cat paperI.tex parameters.txt low_info_hip.txt high_info_hip.txt paperI_figures.tex > ms.tex

submit_paper1: ms.tex
	mkdir $@
	cp ms.tex $@
	cp hipparcos2_properties.ps $@/f1.eps
	cp veldensXY.ps $@/f2.eps
	cp veldensXZ.ps $@/f3.eps
	cp veldensYZ.ps $@/f4.eps
	cp loglike.ps $@/f5.eps
	cp predict_gcs_rand.ps $@/f6.eps
	cp tile_onepred.ps $@/f7.eps
	cp crossval.ps $@/f8a.eps
	cp aic.ps $@/f8b.eps
	cp mdl.ps $@/f8c.eps
	cp mml.ps $@/f8d.eps
	cp bayes.ps $@/f8e.eps
	cp gcs_cross.ps $@/f8f.eps
	cp checkquants.ps $@/f9.eps
	cp predict_gcs_good.ps $@/f10.eps
	cp predict_gcs_bad.ps $@/f11.eps
	cp hist_gcs_like.ps $@/f12.eps
	cp predict_gcs_low.ps $@/f13.eps
	cp predict_gcs_high.ps $@/f14.eps
	cp hist_gcs_ent.ps $@/f15.eps
	cp radecpatches.ps $@/f16.eps
	cp info_hip_low.ps $@/f17.eps
	cp info_hip_high.ps $@/f18.eps
	cp hist_hip_ent.ps $@/f19.eps
	cp annotated_veldist.ps $@/f20.eps
	cp annotated_veldist2.ps $@/f21.eps

.PHONY: clean spotless

clean:
	$(RM) *.aux *.bbl *.blg *.dvi *.log *.out

spotless: clean
	$(RM) $(targets)
