IDL=idl
VELDIST='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
CUT=0.1
LOO_NWS=101
LOO_LOG_WRANGE=[-2.,-.5]
LOO_WRANGE=[.01,.05]
LOO_XSIZE=4.375
LOO_YSIZE=4.375
BACKGROUND_KERNEL_WIDTH=0.05
BACKGROUND_KERNEL='tricube'
BACK_LOGLIKE_FILENAME='back_loglike_all.sav'
SSPSPREAD=0.1
FIXEDALPHA_FITSSP_GRID=50
FITSSP_GRID=50
FITSSP_AGRID=50
FOREGROUND_SIRIUS_AGE=8.62
FOREGROUND_SIRIUS_Z= 0.023
FOREGROUND_SIRIUS_ALPHA=0.9
FOREGROUND_SIRIUS_FIXEDAGE=8.54
FOREGROUND_SIRIUS_FIXEDZ= 0.026
GCS_COLORCUT=1.
LOO_NWS_GCS=101
LOO_LOG_WRANGE_GCS=[-2.,-.5]
GCS_ZCUT= -0.132429
GCS_KERNEL='tricube'
GCS_KERNEL_WIDTH_FULL=0.050699069
GCS_KERNEL_WIDTH_HIGHZ=0.048977880
GCS_KERNEL_WIDTH_LOWZ=0.052480744
SELLWOOD_HARDCUT=0.5
allcmddependon=group_cmd.pro hip2-aumer-giants.sav
paperfigs=veldist.ps cmd.ps cmd_ngc1901.ps cmd_pleiades.ps cmd_hyades.ps cmd_arcturus.ps\
	 cmd_hercules.ps cmd_sirius.ps clusterhist_ngc1901.ps\
	 clusterhist_sirius.ps clusterhist_pleiades.ps\
	 clusterhist_hyades.ps hyadesClusterCMD.ps hyadesClusterHist.ps\
	 model_selection.ps loo_back_log.ps predict_plx_background_random.ps\
	 back_quants.ps\
	 back_loglike_all.ps back_loglike_ngc1901.ps back_loglike_sirius.ps\
	 back_loglike_pleiades.ps back_loglike_hyades.ps\
	 back_loglike_hercules.ps fit_ssp_fixpij_fixalpha_ngc1901.ps\
	 ssp_fixpij_fixalpha_sirius.ps ssp_fixpij_fixalpha_pleiades.ps\
	 ssp_fixpij_fixalpha_hyades.ps ssp_fixpij_fixalpha_hercules.ps\
	 predict_plx_foreground_sirius.ps fit_ssp_fixpij_sirius.ps\
	 fit_ssp_fixpij_ngc1901.ps fit_ssp_fixpij_pleiades.ps\
	 fit_ssp_fixpij_hyades.ps fit_ssp_fixpij_hercules.ps\
	 gcs_cmd.ps gcs_zdist.ps\
	 loo_gcs_all.ps 
papertables=model_selection.txt model_selection_highlowz.txt\
	 avg_metal.txt
extrafigs=pijshist_ngc1901.ps pijshist_sirius.ps pijshist_pleiades.ps\
	 pijshist_hyades.ps pijshist_hercules.ps pijshist_arcturus.ps\
	 extra_clusterhist_ngc1901.ps extra_clusterhist_sirius.ps\
	 extra_clusterhist_pleiades.ps extra_clusterhist_hyades.ps\
	 back_quants_gaussian.ps gcs_zdist_3gauss.ps loo_epanechnikov_log_gcs.ps\
	 loo_gaussian_log_gcs.ps loo_epanechnikov_log_gcs_highz.ps\
	 loo_epanechnikov_log_gcs_lowz.ps loo_gaussian_log_gcs_highz.ps\
	 loo_gaussian_log_gcs_lowz.ps sellwood_gcs_raw.ps\
	 hyadesClusterPij.ps sellwood_gcs_mg_color.ps\
	 sellwood_gcs_mg_gray.ps

all: $(paperfigs) $(extrafigs)
extra: $(extrafigs)
paper: $(paperfigs)
	cp $(paperfigs) ../../tex/groups
.PRECIOUS: gcs_loglike_ngc1901_fixalpha.sav gcs_loglike_sirius_fixalpha.sav gcs_loglike_pleiades_fixalpha.sav gcs_loglike_hyades_fixalpha.sav gcs_loglike_hercules_fixalpha.sav 
.PHONY: model_selection.txt model_selection_highlowz.txt avg_metal.txt
#PAPER
veldist.ps: plot_veldist.pro
	echo "plot_veldist, plotfilename='$@'" | $(IDL)
cmd.ps: ../lsr2/select_sample.pro plot_mscuts.pro
	echo "hip=select_sample(tycho_savefile='../lsr2/tycho.sav',/savefile,/giants,/cmd,filename='hip2-aumer-giants.sav',/overplotmscuts)" | $(IDL)
cmd_ngc1901.ps: $(allcmddependon)
	echo "group_cmd, 0, savefilename=$(VELDIST),cut=$(CUT)" | $(IDL)
cmd_pleiades.ps: $(allcmddependon)
	echo "group_cmd, 3, savefilename=$(VELDIST),cut=$(CUT)" | $(IDL)
cmd_hyades.ps: $(allcmddependon)
	echo "group_cmd, 6, savefilename=$(VELDIST),cut=$(CUT)" | $(IDL)
cmd_arcturus.ps: $(allcmddependon)
	echo "group_cmd, 2, savefilename=$(VELDIST),cut=$(CUT)" | $(IDL)
cmd_hercules.ps: $(allcmddependon)
	echo "group_cmd, 7, savefilename=$(VELDIST),cut=$(CUT)" | $(IDL)
cmd_sirius.ps: $(allcmddependon)
	echo "group_cmd, 1, savefilename=$(VELDIST),cut=$(CUT)" | $(IDL)
clusterhist_ngc1901.ps: iscluster.pro
	echo "iscluster, 0" | $(IDL)
clusterhist_sirius.ps: iscluster.pro
	echo "iscluster, 1" | $(IDL)
clusterhist_pleiades.ps: iscluster.pro
	echo "iscluster, 3" | $(IDL)
clusterhist_hyades.ps: iscluster.pro
	echo "iscluster, 6" | $(IDL)
hyadesClusterCMD.ps: plot_hyades_cmd.pro
	echo "plot_hyades_cmd,savefilename='hyadesMembers.sav',plotfilename='$@'" | $(IDL)
hyadesClusterHist.ps: hyades_cluster_hist.pro
	echo "hyades_cluster_hist, savefilename='hyadesMembers.sav',plotfilename='$@',/colorcut" | $(IDL)
model_selection.ps: illustrate_model_selection.pro
	echo "illustrate_model_selection, plotfilename='$@'" | $(IDL)
loo_back_log.ps: plot_loo_kernel_all.pro Makefile
	echo 'plot_loo_kernel_all, plotfilename="$@",xsize=3.675' | $(IDL)
predict_plx_background_random.ps: predictplx.pro tile_predictplx.pro kernel_phot_plx.pro Makefile
	echo "tile_predictplx, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), mplot=[0,3,3], plotfilename='$@', /random,seed=-1L" | $(IDL)
back_quants.ps: check_quant_back_plx.pro kernel_phot_plx.pro Makefile
	echo "check_quant_back_plx, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), plotfilename='$@', savefilename='$(subst .ps,.sav,$@)',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE)" | $(IDL)
back_loglike_all.ps: hist_loglike_back.pro kernel_phot_plx.pro
	echo "hist_loglike_back, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), plotfilename='$@', savefilename='$(subst .ps,.sav,$@)',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE), group=20" | $(IDL)
back_loglike_ngc1901.ps: hist_loglike_back.pro kernel_phot_plx.pro
	echo "hist_loglike_back, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), savefilename=$(BACK_LOGLIKE_FILENAME),xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE), group=0" | $(IDL)
back_loglike_sirius.ps: hist_loglike_back.pro kernel_phot_plx.pro
	echo "hist_loglike_back, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), savefilename=$(BACK_LOGLIKE_FILENAME),xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE), group=1" | $(IDL)
back_loglike_pleiades.ps: hist_loglike_back.pro kernel_phot_plx.pro
	echo "hist_loglike_back, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), savefilename=$(BACK_LOGLIKE_FILENAME),xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE), group=3" | $(IDL)
back_loglike_hyades.ps: hist_loglike_back.pro kernel_phot_plx.pro
	echo "hist_loglike_back, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), savefilename=$(BACK_LOGLIKE_FILENAME),xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE), group=6" | $(IDL)
back_loglike_hercules.ps: hist_loglike_back.pro kernel_phot_plx.pro
	echo "hist_loglike_back, width=$(BACKGROUND_KERNEL_WIDTH), kernel=$(BACKGROUND_KERNEL), savefilename=$(BACK_LOGLIKE_FILENAME),xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE), group=7" | $(IDL)
fit_ssp_fixpij_fixalpha_ngc1901.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro
	echo "fit_mg_fixpij, 0, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),/fixalpha,grid=$(FIXEDALPHA_FITSSP_GRID)" | $(IDL)
fit_ssp_fixpij_fixalpha_sirius.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro
	echo "fit_mg_fixpij, 1, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),/fixalpha,grid=$(FIXEDALPHA_FITSSP_GRID)" | $(IDL)
fit_ssp_fixpij_fixalpha_pleiades.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro
	echo "fit_mg_fixpij, 3, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),/fixalpha,grid=$(FIXEDALPHA_FITSSP_GRID)" | $(IDL)
fit_ssp_fixpij_fixalpha_hyades.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro
	echo "fit_mg_fixpij, 6, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),/fixalpha,grid=$(FIXEDALPHA_FITSSP_GRID)" | $(IDL)
fit_ssp_fixpij_fixalpha_hercules.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro
	echo "fit_mg_fixpij, 7, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),/fixalpha,grid=$(FIXEDALPHA_FITSSP_GRID)" | $(IDL)
fit_ssp_fixpij_ngc1901.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro Makefile
	echo "fit_mg_fixpij, 0, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',alphaplotfilename='$(subst .ps,_alpha.ps,$@)',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),/dontlogalpha" | $(IDL)
fit_ssp_fixpij_sirius.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro Makefile
	echo "fit_mg_fixpij, 1, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',alphaplotfilename='$(subst .ps,_alpha.ps,$@)',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),/dontlogalpha" | $(IDL)
fit_ssp_fixpij_pleiades.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro Makefile
	echo "fit_mg_fixpij, 3, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',alphaplotfilename='$(subst .ps,_alpha.ps,$@)',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),/dontlogalpha" | $(IDL)
fit_ssp_fixpij_hyades.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro Makefile
	echo "fit_mg_fixpij, 6, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',alphaplotfilename='$(subst .ps,_alpha.ps,$@)',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),/dontlogalpha" | $(IDL)
fit_ssp_fixpij_hercules.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro Makefile
	echo "fit_mg_fixpij, 7, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',alphaplotfilename='$(subst .ps,_alpha.ps,$@)',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),/dontlogalpha" | $(IDL)
predict_plx_foreground_sirius.ps: predictplx.pro predictplx_iso.pro tile_predictplx.pro
	echo "tile_predictplx_iso, 1, $(FOREGROUND_SIRIUS_ALPHA), $(FOREGROUND_SIRIUS_AGE), $(FOREGROUND_SIRIUS_Z), $(FOREGROUND_SIRIUS_FIXEDAGE), $(FOREGROUND_SIRIUS_FIXEDZ), sspspread=$(SSPSPREAD),mplot=[0,3,2],indices=[5942,6274],plotfilename='$@',width=$(BACKGROUND_KERNEL_WIDTH),kernel=$(BACKGROUND_KERNEL),/usegcs" | $(IDL)
model_selection.txt: tabulate_model_selection.pro calc_gcs_prob.pro
	echo "tabulate_model_selection, kernel=$(BACKGROUND_KERNEL),width=$(BACKGROUND_KERNEL_WIDTH),colorcut=$(GCS_COLORCUT),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),sspspread=$(SSPSPREAD)" | $(IDL)
gcs_cmd.ps: plot_gcs_cmd.pro
	echo "plot_gcs_cmd, plotfilename='$@'" | $(IDL)
gcs_zdist.ps: plot_gcs_zdist.pro
	echo "plot_gcs_zdist, plotfilename='$@',/reload" | $(IDL)
loo_gcs_all.ps: plot_loo_kernel_gcs_all.pro Makefile
	echo 'plot_loo_kernel_gcs_all, plotfilename="$@",xsize=3.675' | $(IDL)
model_selection_highlowz.txt: tabulate_model_selection_highlowz.pro
	echo "tabulate_model_selection_highlowz, kernel=$(GCS_KERNEL),full_width=$(GCS_KERNEL_WIDTH_FULL),highz_width=$(GCS_KERNEL_WIDTH_HIGHZ),lowz_width=$(GCS_KERNEL_WIDTH_LOWZ),zcut=$(GCS_ZCUT)" | $(IDL)
sellwood_gcs_mg.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=5.5,ysize=5.,plotfilename='$@',output=3,hardcut=$(SELLWOOD_HARDCUT)" | $(IDL)
avg_metal.txt: avg_metal_groups.pro
	echo "avg_metal_groups" | $(IDL)
#Prep
gcs_loglike_ngc1901_fixalpha.sav: tabulate_model_selection.pro calc_gcs_prob.pro
	echo "tabulate_model_selection, kernel=$(BACKGROUND_KERNEL),width=$(BACKGROUND_KERNEL_WIDTH),colorcut=$(GCS_COLORCUT),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),sspspread=$(SSPSPREAD),/justonegroup,group=0" | $(IDL)
gcs_loglike_sirius_fixalpha.sav: tabulate_model_selection.pro calc_gcs_prob.pro
	echo "tabulate_model_selection, kernel=$(BACKGROUND_KERNEL),width=$(BACKGROUND_KERNEL_WIDTH),colorcut=$(GCS_COLORCUT),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),sspspread=$(SSPSPREAD),/justonegroup,group=1" | $(IDL)
gcs_loglike_pleiades_fixalpha.sav: tabulate_model_selection.pro calc_gcs_prob.pro
	echo "tabulate_model_selection, kernel=$(BACKGROUND_KERNEL),width=$(BACKGROUND_KERNEL_WIDTH),colorcut=$(GCS_COLORCUT),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),sspspread=$(SSPSPREAD),/justonegroup,group=3" | $(IDL)
gcs_loglike_hyades_fixalpha.sav: tabulate_model_selection.pro calc_gcs_prob.pro
	echo "tabulate_model_selection, kernel=$(BACKGROUND_KERNEL),width=$(BACKGROUND_KERNEL_WIDTH),colorcut=$(GCS_COLORCUT),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),sspspread=$(SSPSPREAD),/justonegroup,group=6" | $(IDL)
gcs_loglike_hercules_fixalpha.sav: tabulate_model_selection.pro calc_gcs_prob.pro
	echo "tabulate_model_selection, kernel=$(BACKGROUND_KERNEL),width=$(BACKGROUND_KERNEL_WIDTH),colorcut=$(GCS_COLORCUT),grid=$(FITSSP_GRID),agrid=$(FITSSP_AGRID),sspspread=$(SSPSPREAD),/justonegroup,group=7" | $(IDL)
#EXTRA
pijshist_ngc1901.ps: hist_pijs.pro
	echo "hist_pijs, 0, savefilename=$(VELDIST),plotfilename='$@'" | $(IDL)
pijshist_sirius.ps: hist_pijs.pro
	echo "hist_pijs, 1, savefilename=$(VELDIST),plotfilename='$@'" | $(IDL)
pijshist_pleiades.ps: hist_pijs.pro
	echo "hist_pijs, 3, savefilename=$(VELDIST),plotfilename='$@'" | $(IDL)
pijshist_hyades.ps: hist_pijs.pro
	echo "hist_pijs, 6, savefilename=$(VELDIST),plotfilename='$@'" | $(IDL)
pijshist_hercules.ps: hist_pijs.pro
	echo "hist_pijs, 7, savefilename=$(VELDIST),plotfilename='$@'" | $(IDL)
pijshist_arcturus.ps: hist_pijs.pro
	echo "hist_pijs, 2, savefilename=$(VELDIST),plotfilename='$@'" | $(IDL)
extra_clusterhist_ngc1901.ps: iscluster.pro
	echo "iscluster, 0,cut=0.1, bvcut=0.9,plotfilename='$@'" | $(IDL)
extra_clusterhist_sirius.ps: iscluster.pro
	echo "iscluster, 1,cut=0.1, bvcut=0.9,plotfilename='$@'" | $(IDL)
extra_clusterhist_pleiades.ps: iscluster.pro
	echo "iscluster, 3,cut=0.1, bvcut=0.9,plotfilename='$@'" | $(IDL)
extra_clusterhist_hyades.ps: iscluster.pro
	echo "iscluster, 6,cut=0.1, bvcut=0.9,plotfilename='$@'" | $(IDL)
back_quants_gaussian.ps: check_quant_back_plx.pro kernel_phot_plx.pro Makefile
	echo "check_quant_back_plx, width=0.025, kernel='gaussian', plotfilename='$@', savefilename='$(subst .ps,.sav,$@)',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE)" | $(IDL)
gcs_zdist_3gauss.ps: plot_gcs_zdist.pro
	echo "plot_gcs_zdist, plotfilename='$@',ngauss=3,/reload" | $(IDL)
sellwood_gcs_raw.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=5,ysize=5.,plotfilename='$@',output=1" | $(IDL)
hyadesClusterPij.ps: hyades_high_prob.pro
	echo "hyades_high_prob,savefilename='hyadesMembers_full.sav',plotfilename='$@'" | $(IDL)
sellwood_gcs_mg_gray.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=5,ysize=5.,plotfilename='$@',output=2,hardcut=$(SELLWOOD_HARDCUT)" | $(IDL)
#DEPRECATED
loo_tricube_log.ps: plot_loo_kernel.pro loo_kernel.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel, wrange=$(LOO_LOG_WRANGE),nws=$(LOO_NWS),/rangelog,/log,kernel='tricube',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_epanechnikov_log.ps: plot_loo_kernel.pro loo_kernel.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel, wrange=$(LOO_LOG_WRANGE),nws=$(LOO_NWS),/rangelog,/log,kernel='epanechnikov',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_gaussian_log.ps: plot_loo_kernel.pro loo_kernel.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel, wrange=$(LOO_LOG_WRANGE),nws=$(LOO_NWS),/rangelog,/log,kernel='gaussian',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_tricube.ps: plot_loo_kernel.pro loo_kernel.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel, wrange=$(LOO_WRANGE),nws=$(LOO_NWS),/log,kernel='tricube',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_epanechnikov.ps: plot_loo_kernel.pro loo_kernel.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel, wrange=$(LOO_WRANGE),nws=$(LOO_NWS),/log,kernel='epanechnikov',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_gaussian.ps: plot_loo_kernel.pro loo_kernel.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel, wrange=$(LOO_WRANGE),nws=$(LOO_NWS),/log,kernel='gaussian',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_tricube_log_gcs.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='tricube',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_tricube_log_gcs_highz.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='tricube',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)',/highz,zcut=$(GCS_ZCUT)" | $(IDL)
loo_tricube_log_gcs_lowz.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='tricube',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)',/lowz,zcut=$(GCS_ZCUT)" | $(IDL)
loo_epanechnikov_log_gcs.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='epanechnikov',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_gaussian_log_gcs.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='gaussian',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
loo_epanechnikov_log_gcs_highz.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='epanechnikov',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)',/highz,zcut=$(GCS_ZCUT)" | $(IDL)
loo_epanechnikov_log_gcs_lowz.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='epanechnikov',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)',/lowz,zcut=$(GCS_ZCUT)" | $(IDL)
loo_gaussian_log_gcs_highz.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='gaussian',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)',/highz,zcut=$(GCS_ZCUT)" | $(IDL)
loo_gaussian_log_gcs_lowz.ps: plot_loo_kernel_gcs.pro loo_kernel_gcs.pro kernel_phot_plx.pro Makefile
	echo "plot_loo_kernel_gcs, wrange=$(LOO_LOG_WRANGE_GCS),nws=$(LOO_NWS_GCS),/rangelog,/log,kernel='gaussian',plotfilename='$@',xsize=$(LOO_XSIZE),ysize=$(LOO_YSIZE),savefilename='$(subst .ps,.sav,$@)',/lowz,zcut=$(GCS_ZCUT)" | $(IDL)
#TEST
fit_ssp_fixpij_fixalpha_sirius_test.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro
	echo "fit_mg_fixpij, 1, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),/fixalpha,grid=5" | $(IDL)
fit_ssp_fixpij_sirius_test.ps: calc_logpost_fixedpij.pro fit_mg_fixpij.pro phot_plx.pro read_isochrone.pro get_isoname.pro Makefile
	echo "fit_mg_fixpij, 1, savefilename='$(subst .ps,.sav,$@)',plotfilename='$@',alphaplotfilename='$(subst .ps,_alpha.ps,$@)',back_loglike_filename=$(BACK_LOGLIKE_FILENAME),sspspread=$(SSPSPREAD),grid=5,agrid=4,/dontlogalpha" | $(IDL)
