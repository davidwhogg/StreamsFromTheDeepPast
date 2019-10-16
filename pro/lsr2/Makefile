IDL=idl

paper: cmd.ps lsr_dropblue_groupcut0.5.ps lsr_dropblue_groupcut1.1.ps SMvsPMG.ps

extra: lsr_dropblue_groupcut0.9.ps lsr_dropblue_groupcut1.0.ps lsr_s500_groupcut0.5.ps SMvsPMG_S500.ps SMvsPMG_NGC1901.ps SMvsPMG_S500_NGC1901.ps SMvsPMG_ex1000.ps SMvsPMG_ex1000_NGC1901.ps gcs_vphi_dist.ps gcs_vphi_dist_cumul.ps

#Plots for paper
lsr_dropblue_groupcut0.5.ps: lsr.pro
	echo "lsr,groupcut=0.5,prefix='lsr',/dropblue" | $(IDL)
lsr_dropblue_groupcut1.1.ps: lsr.pro
	echo "lsr,groupcut=1.1,prefix='lsr',/dropblue" | $(IDL)
cmd.ps: select_sample.pro
	echo "hip=select_sample(tycho_savefile='tycho.sav',/savefile,/cmd)" | $(IDL)
SMvsPMG.ps: solarmotion_vs_groupcut.pro lsr.pro
	echo "solarmotion_vs_groupcut, groupcut=[.2,1,.1],/dropblue, plotfilename='$@', savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
#Extra
lsr_dropblue_groupcut0.9.ps: lsr.pro
	echo "lsr,groupcut=0.9,prefix='lsr',/dropblue" | $(IDL)
lsr_dropblue_groupcut1.0.ps: lsr.pro
	echo "lsr,groupcut=1.0,prefix='lsr',/dropblue" | $(IDL)
lsr_s500_groupcut0.5.ps: lsr.pro
	echo "lsr,groupcut=0.5,prefix='lsr',/s500" | $(IDL)
SMvsPMG_S500.ps: solarmotion_vs_groupcut.pro
	echo "solarmotion_vs_groupcut, groupcut=[.2,1,.1],/S500, plotfilename='$@', savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
SMvsPMG_NGC1901.ps: solarmotion_vs_groupcut.pro
	echo "solarmotion_vs_groupcut, groupcut=[.2,1,.1],/dropblue,/remove_ngc1901, plotfilename='$@', savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
SMvsPMG_S500_NGC1901.ps: solarmotion_vs_groupcut.pro
	echo "solarmotion_vs_groupcut, groupcut=[.2,1,.1],/S500,/remove_ngc1901, plotfilename='$@', savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
SMvsPMG_ex1000.ps: solarmotion_vs_groupcut.pro
	echo "solarmotion_vs_groupcut, groupcut=[.2,1,.1],/ex1000, plotfilename='$@', savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
SMvsPMG_ex1000_NGC1901.ps: solarmotion_vs_groupcut.pro
	echo "solarmotion_vs_groupcut, groupcut=[.2,1,.1],/ex1000,/remove_ngc1901, plotfilename='$@', savefilename='$(subst .ps,.sav,$@)'" | $(IDL)
gcs_vphi_dist.ps: gcs_vphi_dist.pro
	echo "gcs_vphi_dist, plotfilename='$@'" | $(IDL)
gcs_vphi_dist_cumul.ps: gcs_vphi_dist.pro
	echo "gcs_vphi_dist, plotfilename='$@',/cum" | $(IDL)
#Samples
hip2-aumer.sav: select_sample.pro
	echo "hip=select_sample(tycho_savefile='tycho.sav',/savefile,/cmd)" | $(IDL)
hip2-aumer_formal.sav: select_sample.pro
	echo "hip=select_sample(tycho_savefile='tycho.sav',/savefile,filename='hip2-aumer_formal.sav',/formal)" | $(IDL)
BovySample.dat: list_sample.pro
	echo "list_sample, filename='$@'" | $(IDL)

