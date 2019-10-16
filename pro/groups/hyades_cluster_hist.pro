;+
;   NAME:
;      hyades_cluster_hist
;   PURPOSE:
;      histogram the isochrone plxs and the observed plxs for the
;      Hyades open cluster
;   CALLING SEQUENCE:
;   INPUT:
;      savefilename - filename for Hyades members data
;      plotfilename - name for plot
;   KEYWORDS:
;      colorcut     - cut everything out that has B-V > 1
;   OUTPUT:
;   REVISION HISTORY:
;      2009-09-14 - Written - Bovy (NYU)
;-
PRO HYADES_CLUSTER_HIST, savefilename=savefilename, $
                         plotfilename=plotfilename, $
                         colorcut=colorcut
;;Restore Hyades members
IF ~file_test(savefilename) THEN BEGIN
    hyades= read_hyades(savefilename=savefilename,/single,/tenpc)
ENDIF ELSE BEGIN
    restore, savefilename
ENDELSE
;;select (possible) MS members
IF keyword_set(colorcut) THEN usethese= where(hyades.vmag-5.*alog10(100./hyades.plx) GE 0.5 $
                                              AND hyades.bvcolor LT 1.) $
  ELSE usethese= where(hyades.vmag-5.*alog10(100./hyades.plx) GE 0.5)
;;Now get the photometric plxs
nhyades= n_elements(usethese)
splog, 'Using '+strtrim(string(nhyades),2)+' Hyades members'
plxs= dblarr(nhyades)
absmags= dblarr(nhyades)
;;Read the isochrone
datafilename='hyades_isochrone2.dat'
isochrone= read_isochrone(datafilename)
iso= dblarr(2,40)
iso[0,*]= isochrone.h09[0:39]
iso[1,*]= isochrone.h10[0:39]
FOR ii=0L, nhyades-1 DO Begin
    BV= hyades[usethese[ii]].bvcolor
    V= hyades[usethese[ii]].vmag
    plxs[ii]= phot_plx(BV,V,iso,absmag=absmag)
    absmags[ii]= absmag
ENDFOR
;;Calculate sigma_M
varM= variance(absmags-hyades[usethese].vmag+5.*alog10(100./hyades[usethese].plx))
extravarPi= varM*(alog(10)/5.)^2.
splog, 'extra uncertainty', sqrt(varM), sqrt(extravarPi)
;;Now histogram the observed and the model plxs, weighted by the
;;probability of being part of the moving group
plotthese= where(plxs NE -1)
npix=19
k_print, filename=plotfilename
hogg_plothist, hyades[usethese[plotthese]].plx, xrange=[0,49],$
  linestyle=2, npix=npix,$
  xtitle=textoidl('\pi_{obs}, \pi_{iso} [mas]'), position=[.1,.5,.5,.9]
hogg_plothist, plxs[plotthese],/overplot
legend, ['Hyades cluster'], /right,box=0
;;Then plot the normalized difference
;;Overplot what happens when you add in the extra variance
weight=dblarr(nhyades)+1/3.
hogg_plothist, (hyades[usethese[plotthese]].plx-plxs[plotthese])/hyades[usethese[plotthese]].e_plx,$
  xtitle=textoidl('(\pi_{obs}-\pi_{iso})\sigma^{-1}_\pi'), npix=npix,$
  position=[.5,.5,.9,.9],/noerase, ytickname=REPLICATE(' ',30), hist=hist, xrange=[-9.,9.]
hogg_plothist, (hyades[usethese[plotthese]].plx-plxs[plotthese])/(sqrt(hyades[usethese[plotthese]].e_plx^2.+extravarPi*(hyades[usethese[plotthese]].plx)^2.)),linestyle=2, npix=floor(1.5*npix), weight=weight, /overplot
rightytitle, position=[.5,.5,.9,.9], yrange=[-0.1,1.1]*max(hist), xticks=1,$
  xtickname=REPLICATE(' ',30), xrange=[-9,9]
;;this yrange comes from hogg_plothist

indx=where((hyades[usethese[plotthese]].plx-plxs[plotthese])/(sqrt(hyades[usethese[plotthese]].e_plx^2.+extravarPi*(hyades[usethese[plotthese]].plx)^2.)) GT -3.)
splog, 'mean and variance', mean((hyades[usethese[plotthese[indx]]].plx-plxs[plotthese[indx]])/(sqrt(hyades[usethese[plotthese[indx]]].e_plx^2.+extravarPi*(hyades[usethese[plotthese[indx]]].plx)^2.))), mean(((hyades[usethese[plotthese[indx]]].plx-plxs[plotthese[indx]])/(sqrt(hyades[usethese[plotthese[indx]]].e_plx^2.+extravarPi*(hyades[usethese[plotthese[indx]]].plx)^2.)))^2)


k_end_print

END
