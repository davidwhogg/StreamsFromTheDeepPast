;+
;   NAME:
;      hist_pijs
;   PURPOSE:
;      make a histogram of the distribution of pijs for group j
;   CALLING SEQUENCE:
;   INPUT:
;      group - which moving group
;      savefilename - filename with the velocity field
;      datafilename - filename with the data
;      plotfilename - filename for plot
;   OUTPUT:
;      histogram
;   REVISION HISTORY:
;      2009-09-15 - Written - Bovy (NYU)
;-
PRO HIST_PIJS, group, savefilename=savefilename, $
               datafilename=datafilename, plotfilename=plotfilename
IF group EQ 0 THEN BEGIN
    name='NGC 1901'
ENDIF ELSE IF group EQ 1 THEN BEGIN
    name='Sirius'
ENDIF ELSE IF group EQ 2 THEN BEGIN
    name='Arcturus'
ENDIF ELSE IF group EQ 3 THEN BEGIN
    name='Pleiades'
ENDIF ELSE IF group EQ 5 THEN BEGIN
    name='Background'
ENDIF ELSE IF group EQ 6 THEN BEGIN
    name='Hyades'
ENDIF ELSE IF group EQ 7 THEN BEGIN
    name='Hercules'
ENDIF ELSE IF group EQ 8 THEN BEGIN
    name='Background'
ENDIF ELSE IF group EQ 9 THEN BEGIN
    name='Background'
ENDIF
;;Restore data
IF ~keyword_set(savefilename) THEN $
  savefilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
;;Restore data
IF ~keyword_set(datafilename) THEN datafilename='../lsr2/hip2-aumer.sav'
restore, filename=datafilename
nhip= n_elements(hip.hip)
ydata= transpose([[hip.vl], [hip.vb]])
ycovar= hip.vlvbc
projection= hip.nsm
;;Restore solution
restore, filename=savefilename
;;First compute all the probabilities
splog, 'Assigning group probabilities...'
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
IF group EQ 3 THEN BEGIN;;Pleiades is a special case
    logpost[group,*]= alog(exp(logpost[group,*])+exp(logpost[group+1,*]))
ENDIF
pijs= exp(logpost[group,*])
;;Compute cumulative sums of pijs
splog, 'Calculating cumulative sum of pijs...'
nsamples=1001
xs=dindgen(nsamples)/(nsamples-1)
cumpijs= dblarr(nsamples)
sortedpijs= pijs[sort(pijs)]
cumsum= 0.
indx=0.
FOR ii=0L, nsamples-1 DO BEGIN
    WHILE sortedpijs[indx] LT xs[ii] DO BEGIN
        indx+= 1
        cumsum+= sortedpijs[indx]
        IF indx EQ n_elements(hip.hip)-1 THEN BREAK
    ENDWHILE
    cumpijs[ii]= cumsum
    IF indx EQ n_elements(hip.hip)-1 THEN BEGIN
        cumpijs[ii+1:nsamples-1]= cumsum
        BREAK
    ENDIF
ENDFOR
splog, 'Plotting...'
;;First plot the histogram
charsize=1
hogg_plothist, pijs, xrange=[0,.99], hist=hist,/dontplot
weight= dblarr(n_elements(hip.hip))+1./n_elements(hip.hip)/n_elements(hist)
k_print, filename=plotfilename;, xsize=3.375, ysize=3.5
hogg_plothist, pijs,xrange=[0,.99],xtitle=textoidl('p_{i}'), $
  ytitle='Fraction', weight=weight, charsize=charsize, position=[.1,.5,.5,.9]
legend, [name],box=0.,/right
;;Then plot the cumulative sum next to it
yrange=[0,1]
djs_plot, xs, cumpijs/cumsum, xtitle='x', xrange=[0,1], yrange=yrange, $
  position=[.5,.5,.9,.9],/noerase, ytickname=REPLICATE(' ',30), $
  charsize=charsize
rightytitle, position=[.5,.5,.9,.9], yrange=yrange, xticks=1,$
  xtickname=REPLICATE(' ',30), xrange=[0,1], $
  ytitle=textoidl('\Sigma_{p_i < x} p_i / \Sigma_i p_i')
k_end_print
END
