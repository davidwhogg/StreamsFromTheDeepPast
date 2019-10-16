basedir= '../ALMSv2_reconvergev2/'
ks= [3,5,10,13,15,20]
ws= [1.,4.,9.]
modelks_high= 20
modelks_low= 4

gcs_cross.ps: plot_kwgcs.pro
	echo "plot_kwgcs, filename='gcs_cross.sav',plotfilename='gcs_cross.ps',/hor,basedir=$(basedir),ws=$(ws), ks_low=$(modelks_low), ks_high=$(modelks_high)" | idl

theory_veldist.ps: annotated_veldist.pro
	echo "annotated_veldist, K=1, w=1, sample=6,filename='$@'" | idl
