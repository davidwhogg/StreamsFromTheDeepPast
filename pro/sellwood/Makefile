IDL=idl
VELDIST='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
SELLWOOD_HARDCUT=0.5
all: sellwood_gcs_mg_zoom.ps sellwood_gcs_zoom.ps sellwood_gcs_mg_zoom_6panel.ps\
	 sellwood_gcs_hyades.ps
sellwood_gcs_mg_zoom.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=5.5,ysize=5.,plotfilename='$@',output=2,hardcut=$(SELLWOOD_HARDCUT),/colour" | $(IDL)
sellwood_gcs_zoom.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=5.5,ysize=5.,plotfilename='$@',output=1,hardcut=$(SELLWOOD_HARDCUT)" | $(IDL)
sellwood_gcs_mg_zoom_6panel.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=5.5,ysize=5.,plotfilename='$@',output=3,hardcut=$(SELLWOOD_HARDCUT)" | $(IDL)
sellwood_gcs_hyades.ps: sellwood_model.pro Makefile
	echo "sellwood_model, sample='GCS',disk='mestel',xsize=3.375,ysize=5.,plotfilename='$@',/justhyades,hardcut=$(SELLWOOD_HARDCUT)" | $(IDL)
