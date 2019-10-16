;+
;   NAME:
;      iscluster
;   PURPOSE:
;      compare plxs from isochrone and observation to see whether the
;      hypothesis that all the stars are part of an evaporating
;      cluster is correct. Weights stars in the histogram with the
;      probability that they are part of the cluster
;   CALLING SEQUENCE:
;   INPUT:
;      group - which group?
;      cut   - member prob has to be higher than this to be considered
;      bvcut - bv has to be small than this to be considered
;      savefilename - filename that holds the velocity distribution
;   OUTPUT:
;   REVISION HISTORY:
;      2009-09-11 - Written - Bovy (NYU)
;-
PRO ISCLUSTER, group, savefilename=savefilename, cut=cut, bvcut=bvcut, plotfilename=plotfilename
IF ~keyword_set(cut) THEN cut=0.
IF ~keyword_set(bvcut) THEN bvcut= 10.
IF ~keyword_set(savefilename) THEN $
  savefilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
basefilename='clusterhist_'
IF group EQ 0 THEN BEGIN
    name='NGC 1901'
    IF ~keyword_set(plotfilename) THEN plotfilename=basefilename+'ngc1901'
ENDIF ELSE IF group EQ 1 THEN BEGIN
    name='Sirius'
    IF ~keyword_set(plotfilename) THEN plotfilename=basefilename+'sirius'
ENDIF ELSE IF group EQ 3 THEN BEGIN
    name='Pleiades'
    IF ~keyword_set(plotfilename) THEN plotfilename=basefilename+'pleiades'
ENDIF ELSE IF group EQ 6 THEN BEGIN
    name='Hyades'
    IF ~keyword_set(plotfilename) THEN plotfilename=basefilename+'hyades'
ENDIF
IF ~keyword_set(plotfilename) THEN plotfilename+='.ps'
;;Restore data
restore, filename='../lsr2/hip2-aumer.sav'
nhip= n_elements(hip.hip)
ydata= transpose([[hip.vl], [hip.vb]])
ycovar= hip.vlvbc
projection= hip.nsm
;;Restore solution
restore, filename=savefilename
;;First compute all the probabilities
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
;;Sort the stars
IF group EQ 3 THEN BEGIN
    logpost[group,*]= alog(exp(logpost[group,*])+exp(logpost[group+1,*]))
ENDIF
;;Select stars that we are going to consider
indx= where(hip.bvcolor LT bvcut AND exp(logpost[group,*]) GT cut)
n_in_sample= n_elements(indx)
splog, 'n_in_sample = '+strtrim(string(n_in_sample),2)
;;Now compute the photometric plxs
IF group EQ 3 THEN BEGIN
    datafilename='pleiades_isochrone.dat'
ENDIF ELSE IF group EQ 6 THEN BEGIN
    datafilename='hyades_isochrone.dat'
ENDIF ELSE IF group EQ 1 THEN BEGIN
    datafilename='sirius_isochrone.dat'
ENDIF ELSE IF group EQ 0 THEN BEGIN
    datafilename='ngc1901_isochrone.dat'
ENDIF
isochrone= read_isochrone(datafilename)
iso= dblarr(2,40)
iso[0,*]= isochrone.h09[0:39]
iso[1,*]= isochrone.h10[0:39]
plxs= dblarr(n_in_sample)
FOR ii=0L, n_in_sample-1 DO Begin
    BV= hip[indx[ii]].bvcolor
    V= hip[indx[ii]].vmag
;    IF BV GE 1.5 OR BV LE .2 THEN CONTINUE
    plxs[ii]= phot_plx(BV,V,iso)
ENDFOR
;;Now histogram the observed and the model plxs, weighted by the
;;probability of being part of the moving group
plotthese= where(plxs NE -1)
weights= exp(logpost[group,indx[plotthese]])
k_print, filename=plotfilename
hogg_plothist, hip[indx[plotthese]].plx,weight=weights, xrange=[0,29],linestyle=2,$
  xtitle=textoidl('\pi_{obs}, \pi_{iso} [mas]'), position=[.1,.5,.5,.9]
hogg_plothist, plxs[plotthese],weight=weights,/overplot
legend, [name], /right,box=0
;;Then plot the normalized difference
hogg_plothist, (hip[indx[plotthese]].plx-plxs[plotthese])/hip[indx[plotthese]].e_plx,$
  weight=weights, xrange=[-19,10], $
  xtitle=textoidl('(\pi_{obs}-\pi_{iso})\sigma^{-1}_\pi'), $
  position=[.5,.5,.9,.9],/noerase, ytickname=REPLICATE(' ',30), hist=hist
rightytitle, position=[.5,.5,.9,.9], yrange=[-0.1,1.1]*max(hist), xticks=1,$
  xtickname=REPLICATE(' ',30), xrange=[-19,10]
;;this yrange comes from hogg_plothist

k_end_print


END
