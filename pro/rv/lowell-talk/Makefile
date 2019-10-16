K= 10
w= 4
basedir= '../ALMSv2_reconvergev2/'

SRCDIR= ~/streams/pro/rv/

all: veldist_3_4.jpg veldist_5_4.jpg veldist_10_4.jpg veldist_15_1.jpg veldist_20_4.jpg illustrate_photplx.jpg

veldist_3_4.ps: $(SRCDIR)plot_projected_gaussians.pro $(SRCDIR)annotated_veldist.pro
	echo "annotated_veldist, K=3, w=4, filename='$@', basedir=$(basedir)" | idl

veldist_5_4.ps: $(SRCDIR)plot_projected_gaussians.pro $(SRCDIR)annotated_veldist.pro
	echo "annotated_veldist, K=5, w=4, filename='$@', basedir=$(basedir)" | idl

veldist_10_4.ps: $(SRCDIR)plot_projected_gaussians.pro $(SRCDIR)annotated_veldist.pro
	echo "annotated_veldist, K=10, w=4, filename='$@', basedir=$(basedir)" | idl

veldist_15_1.ps: $(SRCDIR)plot_projected_gaussians.pro $(SRCDIR)annotated_veldist.pro
	echo "annotated_veldist, K=15, w=1, filename='$@', basedir=$(basedir)" | idl

veldist_20_4.ps: $(SRCDIR)plot_projected_gaussians.pro $(SRCDIR)annotated_veldist.pro
	echo "annotated_veldist, K=20, w=4, filename='$@', basedir=$(basedir)" | idl

illustrate_photplx.ps: illustrate_photplx.pro
	echo "illustrate_photplx,plotfilename='$@'" | idl

%.jpg: %.ps
	convert $< $@
