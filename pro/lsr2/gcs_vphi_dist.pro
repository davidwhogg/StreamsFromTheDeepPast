;+
;   NAME:
;      gcs_vphi_dist
;   PURPOSE:
;      plot the distribution of GCS vphi velocities (like in Binney,
;      2009), leaving out moving group members
;   INPUT:
;      plotfilename - filename for plot
;      groupcuts    - cuts for the plots [.4,.7]
;      Vsun         - phi (V) component of the Solar motion
;      circv        - assumed circular velocity
;   KEYWORDS:
;      cum - plot cumulative distribution
;   OUTPUT:
;      plot
;   REVISION HISTORY:
;      2009-10-09 - Written - Bovy (NYU)
;-
PRO GCS_VPHI_DIST, plotfilename=plotfilename, groupcut1=groupcuts, Vsun=Vsun, $
                   cum=cum

IF ~keyword_set(groupcuts) THEN groupcuts=[0.5,0.2]
IF~keyword_set(Vsun) THEN Vsun= 4.4
IF~keyword_set(circv) THEN circv= 220.

xrange=[0.,319]
nsteps= 101
xcumul= dindgen(nsteps)/double(nsteps)*(xrange[1]-xrange[0])+xrange[0]
cumul= dblarr(nsteps)

restore, filename='fitv10_bestV_V4.0_sample6.sav'
ngauss= n_elements(amp)
;;Get the GCS data
gcsfilename='/global/data/rv/gcs/gcs.sav'
;;Restore the GCS sample
IF file_test(gcsfilename) THEN BEGIN
    splog, 'Restoring the GCS savefile: '+gcsfilename
    restore, gcsfilename
ENDIF ELSE BEGIN
    splog, 'GCS savefile does not exist in '+gcsfilename
    splog, 'GCS savefile must exist'
    splog, 'Returning...'
    RETURN
ENDELSE
ngcs= n_elements(gcs.hip)
;;Calculate all of the V velocities
v=dblarr(ngcs)
FOR ii=0L, ngcs-1 DO BEGIN
    v[ii]= (invert(gcs[ii].sm)##[gcs[ii].vr,gcs[ii].vl,gcs[ii].vb])[1]+Vsun+circv
ENDFOR
;;First histogram: all GCS
k_print, filename=plotfilename
charsize=.6
charthick=2
IF keyword_set(cum) THEN BEGIN
    FOR ii=0L, nsteps-1 DO BEGIN
        cumul[ii]= n_elements(where(v LE xcumul[ii]))/double(ngcs)
        if cumul[ii] EQ -1 THEN cumul[ii]= 0.
    ENDFOR
    djs_plot, xcumul, cumul, xrange=xrange, yrange=[0.,1.],$
    position=[0.1,0.5,0.5,0.9], charsize=charsize, $
      charthick=charthick, xtickname=strarr(30)+' '
    allcumul= cumul
    charlegend=.8
    legend, ['No moving-group'], /left, box=0., charsize=charlegend, position=[20.,.95]
    legend, ['stars excluded'], /left, box=0., charsize=charlegend, position=[20.,.9]
ENDIF ELSE BEGIN
    npix= 320/5+1
    hogg_plothist, v, xrange=xrange, /totalweight, npix=npix, /dontplot, hist=hist
    yrange= [0.,1.1]*max(hist)
    hogg_plothist, v, xrange=xrange, /totalweight, npix=npix,$
      position=[0.1,0.5,0.5,0.9], charsize=charsize, $
      charthick=charthick, xtickname=strarr(30)+' ', yrange=yrange
    charlegend=.8
    legend, ['No moving-group'], /left, box=0., charsize=charlegend, position=[20.,1025]
    legend, ['stars excluded'], /left, box=0., charsize=charlegend, position=[20.,975]
ENDELSE


;;Now leave out groupcuts[0]
ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
ycovar= gcs.vrvlvbc
projection= gcs.sm
;;First construct the sample, determine the data points which aren't
;;part of a moving group (based on a given probability)
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
;;Make the sample cut
;groups=[0,1,2,3,4,6,7];;include NGCS1901
groups=[1,2,3,4,6,7]
ngroups= n_elements(groups)
sample_indx= lonarr(ngcs)
nsample_indx= 0L
groupcut= groupcuts[0]
FOR ii=0L, ngcs-1 DO BEGIN
    ingroup= 0
    FOR jj=0L, ngroups-1 DO BEGIN
        IF keyword_set(relcut) THEN BEGIN
            IF logpost[groups[jj],ii] NE max(logpost[*,ii]) THEN CONTINUE
            sortindx= sort(logpost[*,ii])
            IF exp(logpost[groups[jj],ii]) GE $
              relcut*exp(logpost[sortindx[ngauss-2],ii]) THEN BEGIN
                ingroup= 1
                BREAK
            ENDIF
        ENDIF ELSE IF exp(logpost[groups[jj],ii]) GE groupcut THEN BEGIN
            ingroup= 1
            BREAK
        ENDIF
    ENDFOR
    IF ingroup EQ 0 THEN BEGIN
        sample_indx[nsample_indx]= ii
        nsample_indx+= 1
    ENDIF
ENDFOR
sample_indx= sample_indx[0:nsample_indx-1]
frac_excluded= (ngcs-nsample_indx)/double(ngcs)
IF ~keyword_set(silent) THEN splog, strtrim(string(nsample_indx),2)+" stars out of " +strtrim(string(ngcs),2)+" not considered to be in a moving group"
;;Now plot these
IF keyword_set(cum) THEN BEGIN
    FOR ii=0L, nsteps-1 DO BEGIN
        cumul[ii]= n_elements(where(v[sample_indx] LE xcumul[ii]))/double(n_elements(sample_indx))
        if cumul[ii] EQ -1 THEN cumul[ii]= 0.
    ENDFOR
    djs_plot, xcumul, cumul, xrange=xrange, yrange=[0.,1.],$
    position=[0.5,0.5,0.9,0.9], charsize=charsize, $
      charthick=charthick, xtickname=strarr(30)+' ',/noerase, ytickname=strarr(30)+' '
    djs_oplot, xcumul, allcumul, color=djs_icolor('gray')
    allcumul= cumul
    legend, ['Moving-group stars'], /left, box=0., charsize=charlegend, position=[20.,.95]
    legend, [textoidl('excluded at p_{MG} = ')+strtrim(string(groupcuts[0],'(F3.1)'),2)], /left, box=0., charsize=charlegend, position=[20.,.9]
ENDIF ELSE BEGIN
    weights= dblarr(n_elements(sample_indx))+double(ngcs)/n_elements(sample_indx)
    hogg_plothist, v[sample_indx], xrange=xrange, /totalweight, npix=npix, $
      position=[0.5,0.5,0.9,0.9], /NOERASE, charsize=charsize, charthick=charthick, $
      xtickname=strarr(30)+' ', ytickname=strarr(30)+' ', yrange=yrange, weight=weights
    legend, ['Moving-group stars'], /left, box=0., charsize=charlegend, position=[20.,1025]
    legend, [textoidl('excluded at p_{MG} = ')+strtrim(string(groupcuts[0],'(F3.1)'),2)], /left, box=0., charsize=charlegend, position=[20.,975]
ENDELSE


groups=[0,1,2,3,4,6,7];;include NGCS1901
;groups=[1,2,3,4,6,7]
ngroups= n_elements(groups)
sample_indx= lonarr(ngcs)
nsample_indx= 0L
groupcut= groupcuts[0]
FOR ii=0L, ngcs-1 DO BEGIN
    ingroup= 0
    FOR jj=0L, ngroups-1 DO BEGIN
        IF keyword_set(relcut) THEN BEGIN
            IF logpost[groups[jj],ii] NE max(logpost[*,ii]) THEN CONTINUE
            sortindx= sort(logpost[*,ii])
            IF exp(logpost[groups[jj],ii]) GE $
              relcut*exp(logpost[sortindx[ngauss-2],ii]) THEN BEGIN
                ingroup= 1
                BREAK
            ENDIF
        ENDIF ELSE IF exp(logpost[groups[jj],ii]) GE groupcut THEN BEGIN
            ingroup= 1
            BREAK
        ENDIF
    ENDFOR
    IF ingroup EQ 0 THEN BEGIN
        sample_indx[nsample_indx]= ii
        nsample_indx+= 1
    ENDIF
ENDFOR
sample_indx= sample_indx[0:nsample_indx-1]
frac_excluded= (ngcs-nsample_indx)/double(ngcs)
IF ~keyword_set(silent) THEN splog, strtrim(string(nsample_indx),2)+" stars out of " +strtrim(string(ngcs),2)+" not considered to be in a moving group"
;;Now plot these
IF keyword_set(cum) THEN BEGIN
    FOR ii=0L, nsteps-1 DO BEGIN
        cumul[ii]= n_elements(where(v[sample_indx] LE xcumul[ii]))/double(n_elements(sample_indx))
        if cumul[ii] EQ -1 THEN cumul[ii]= 0.
    ENDFOR
    djs_plot, xcumul, cumul, xrange=xrange, yrange=[0.,1.],$
    position=[0.1,0.1,0.5,0.5], charsize=charsize, $
      charthick=charthick, /noerase
    djs_oplot, xcumul, allcumul, color=djs_icolor('gray')
    allcumul= cumul
    legend, ['Moving-group stars'], /left, box=0., charsize=charlegend, position=[20.,.95]
    legend, ['(including NGC 1901)'], /left, box=0., charsize=charlegend, position=[20.,.9]
    legend, [textoidl('excluded at p_{MG} = ')+strtrim(string(groupcuts[0],'(F3.1)'),2)], /left, box=0., charsize=charlegend, position=[20.,.85]
ENDIF ELSE BEGIN
    weights= dblarr(n_elements(sample_indx))+double(ngcs)/n_elements(sample_indx)
    hogg_plothist, v[sample_indx], xrange=xrange, /totalweight, npix=npix, $
      position=[0.1,0.1,0.5,0.5], /NOERASE, charsize=charsize, charthick=charthick, $
      yrange=yrange, weight=weights
    legend, ['Moving-group stars'], /left, box=0., charsize=charlegend, position=[20.,1025]
    legend, ['(including NGC 1901)'], /left, box=0., charsize=charlegend, position=[20.,975]
    legend, [textoidl('excluded at p_{MG} = ')+strtrim(string(groupcuts[0],'(F3.1)'),2)], /left, box=0., charsize=charlegend, position=[20.,925]
ENDELSE




;groups=[0,1,2,3,4,6,7];;include NGCS1901
groups=[1,2,3,4,6,7]
ngroups= n_elements(groups)
sample_indx= lonarr(ngcs)
nsample_indx= 0L
groupcut= groupcuts[1]
FOR ii=0L, ngcs-1 DO BEGIN
    ingroup= 0
    FOR jj=0L, ngroups-1 DO BEGIN
        IF keyword_set(relcut) THEN BEGIN
            IF logpost[groups[jj],ii] NE max(logpost[*,ii]) THEN CONTINUE
            sortindx= sort(logpost[*,ii])
            IF exp(logpost[groups[jj],ii]) GE $
              relcut*exp(logpost[sortindx[ngauss-2],ii]) THEN BEGIN
                ingroup= 1
                BREAK
            ENDIF
        ENDIF ELSE IF exp(logpost[groups[jj],ii]) GE groupcut THEN BEGIN
            ingroup= 1
            BREAK
        ENDIF
    ENDFOR
    IF ingroup EQ 0 THEN BEGIN
        sample_indx[nsample_indx]= ii
        nsample_indx+= 1
    ENDIF
ENDFOR
sample_indx= sample_indx[0:nsample_indx-1]
frac_excluded= (ngcs-nsample_indx)/double(ngcs)
IF ~keyword_set(silent) THEN splog, strtrim(string(nsample_indx),2)+" stars out of " +strtrim(string(ngcs),2)+" not considered to be in a moving group"
;;Now plot these
IF keyword_set(cum) THEN BEGIN
    FOR ii=0L, nsteps-1 DO BEGIN
        cumul[ii]= n_elements(where(v[sample_indx] LE xcumul[ii]))/double(n_elements(sample_indx))
        if cumul[ii] EQ -1 THEN cumul[ii]= 0.
    ENDFOR
    djs_plot, xcumul, cumul, xrange=[xrange[0]+0.0001,xrange[1]], yrange=[0.,1.],$
    position=[0.5,0.1,0.9,0.5], ytickname=strarr(30)+' ', charsize=charsize, $
      charthick=charthick,/noerase
    djs_oplot, xcumul, allcumul, color=djs_icolor('gray')
    allcumul= cumul
    legend, ['Moving-group stars'], /left, box=0., charsize=charlegend, position=[20.,.95]
    legend, [textoidl('excluded at p_{MG} = ')+strtrim(string(groupcuts[1],'(F3.1)'),2)], /left, box=0., charsize=charlegend, position=[20.,.9]
ENDIF ELSE BEGIN
    weights= dblarr(n_elements(sample_indx))+double(ngcs)/n_elements(sample_indx)
    hogg_plothist, v[sample_indx], xrange=[xrange[0]+0.0001,xrange[1]], /totalweight, npix=npix, $
      position=[0.5,0.1,0.9,0.5], /NOERASE, charsize=charsize, charthick=charthick, $
      yrange=yrange, ytickname=strarr(30)+' ', weight=weights
    legend, ['Moving-group stars'], /left, box=0., charsize=charlegend, position=[20.,1025]
    legend, [textoidl('excluded at p_{MG} = ')+strtrim(string(groupcuts[1],'(F3.1)'),2)], /left, box=0., charsize=charlegend, position=[20.,975]
ENDELSE


k_end_print
END
