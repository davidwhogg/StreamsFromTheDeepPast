;+
;   NAME:
;      predictplx_iso
;   PURPOSE:
;      plot the prediction of the photometric parallax from the
;      foreground model
;   INPUT:
;      id - HIP number
;      alpha - alpha parameter of the foreground model
;      age - age parameter of the foreground model (log Age)
;      z - metallicity parameter of the foreground model
;      datafilename - filename with the data
;      nplx - number of plxs to evaluate the probability at
;      xrange - range of plx
;      yrange - range of p(plx)
;      confidence - confidence interval (default .95) on the
;                   plot
;      isodir - directory that holds the isochrones
;      sspspread - spread to add to the SSP prediction
;      width - width parameter of the kernel
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;   KEYWORDS:
;      hiplabel - put Hipparcos identifier as a label
;      byindx - if set, the id is the indx in the hipparcos struct
;      noxlabels - no x labels
;      noylabels - no y labels
;      dontplot - don't plot
;      gzip - if set, gunzip and gzip the isochrones
;      usegcs - if set, use the GCS catalog and not the Hipparcos catalog
;   OUTPUT:
;      plot
;      quantile at which the true parallax is found if asked for
;   REVISION HISTORY:
;      2009-11-10 - Written - Bovy (NYU)
;-
PRO PREDICTPLX_ISO, id, group, alpha, age, z, $
                    datafilename=datafilename, xrange=xrange, yrange=yrange, $
                    nplx=nplx, confidence=confidence, hiplabel=hiplabel,$
                    isodir=isodir, sspspread=sspspread,$
                    _EXTRA=KeywordsForPlot, byindx=byindx, $
                    noxlabels=noxlabels, noylabels=noylabels, quantile=quantile, $
                    dontplot=dontplot, kernel=kernel, width=width, gzip=gzip, $
                    solutionfilename=solutionfilename, usegcs=usegcs
IF ~keyword_set(xrange) THEN xrange=[.05,60.]
IF ~keyword_set(nplx) THEN nplx= 1001
IF ~keyword_set(confidence) THEN confidence= .95
IF ~keyword_set(sspspread) THEN sspspread= 0.
IF ~keyword_set(isodir) THEN isodir='isochrones/'
IF ~keyword_set(width) THEN width= 0.05
IF ~keyword_set(kernel) THEN kernel='tricube'
;;Restore the Hipparcos data
thisdatafilename='../lsr2/hip2-aumer.sav'
restore, filename=thisdatafilename
nhip= n_elements(hip.hip)
data= dblarr(nhip,2)
data[*,0]= hip.bvcolor
data[*,1]= hip.vmag-5.*alog10(100./hip.plx)

IF ~keyword_set(datafilename) THEN BEGIN
    IF keyword_set(usegcs) THEN BEGIN
        thisdatafilename='gcs.sav'
    ENDIF ELSE BEGIN
        thisdatafilename='../lsr2/hip2-aumer.sav'
    ENDELSE
ENDIF ELSE BEGIN
    thisdatafilename=datafilename
ENDELSE
restore, filename=thisdatafilename
IF keyword_set(usegcs) THEN hip= gcs
nhip= n_elements(hip.hip)
IF ~keyword_set(solutionfilename) THEN $
  solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
IF keyword_set(byindx) THEN thisid= hip[id].hip ELSE thisid= id
hipindx= (where(hip.hip EQ thisid))[0]
IF keyword_set(usegcs) THEN BEGIN
    ydata= transpose([[hip[hipindx].vr]])
    ycovar= hip[hipindx].vrc
    ycovar= reform(ycovar,1,1,1)
    projection= hip[hipindx].sm[*,0]
    projection= reform(projection,3,1,1)
ENDIF ELSE BEGIN
    ydata= transpose([[hip[hipindx].vl], [hip[hipindx].vb]])
    ycovar= [hip[hipindx].vlvbc]
    projection= [hip[hipindx].nsm]
ENDELSE
;;Restore solution
restore, filename=solutionfilename
;;First compute all the probabilities
assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
IF group EQ 3 THEN BEGIN
    grouplogpost[group,*]= alog(exp(grouplogpost[group,*])+exp(grouplogpost[group+1,*]))
ENDIF

;;Read the isochrones
isoname= get_isoname(z,metal=zout)
IF keyword_set(gzip) THEN BEGIN
    cmd= 'gunzip '+isodir+isoname+'.gz'
    spawn, cmd
ENDIF
iso= read_isochrone(isodir+isoname)
IF keyword_set(gzip) THEN BEGIN
    cmd= 'gzip '+isodir+isoname
    spawn, cmd
ENDIF
ages= iso.h01[UNIQ(iso.h01, SORT(iso.h01))]
ii=0L
WHILE ii LT n_elements(ages)-1 DO BEGIN
    IF ages[ii] GE age THEN BREAK
    ii+= 1
ENDWHILE
IF ii GE n_elements(ages)-1 THEN BEGIN
    isoage= ages[ii]
ENDIF ELSE IF ii EQ 0 THEN BEGIN
    isoage= ages[0]
ENDIF ELSE IF abs(age-ages[ii]) GT abs(age-ages[ii-1]) THEN BEGIN
    isoage= ages[ii-1]
ENDIF ELSE BEGIN
    isoage= ages[ii]
ENDELSE
thisiso_indx= where(iso.h01 EQ isoage)
nthisiso= n_elements(thisiso_indx)
thisiso= dblarr(2,nthisiso)
thisiso[0,*]= iso.h09[thisiso_indx]
thisiso[1,*]= iso.h10[thisiso_indx]

;;Calculate quantile
IF arg_present(quantile) THEN BEGIN
    ;;Sample the distribution finer, on a wider grid
    quant_grid= 10*nplx+1
    quant_range=[0.01,200.]
    quant_step= double(quant_range[1]-quant_range[0])/(10*nplx)
    quant_plxs= dindgen(quant_grid)*quant_step + quant_range[0]
    quant_dist= dblarr(quant_grid)
    FOR kk=0L, quant_grid-1 DO BEGIN
        back_loglike= kernel_phot_plx(quant_plxs[kk],$
                                        hip[hipindx].bvcolor,$
                                        hip[hipindx].vmag,$
                                        hip[hipindx].e_plx,data,$
                                        width=width,kernel=kernel)
        quant_dist[kk]= exp(iso_phot_plx(quant_plxs[kk],$
                                     hip[hipindx].bvcolor,$
                                     hip[hipindx].vmag,$
                                     hip[hipindx].e_plx,thisiso,$
                                     alpha=alpha,back_loglike=back_loglike, $
                                     sspspread=sspspread)*(exp(grouplogpost[group,0])))*exp(back_loglike*(1-exp(grouplogpost[group,0])))
    ENDFOR
    ;;Normalize to one
    quant_dist=quant_dist/(total(quant_dist,/NaN)*quant_step)
    ;;Find where vr is
    IF quant_range[1] LT hip[hipindx].plx THEN BEGIN
        quantile= 1.
    ENDIF ELSE IF quant_range[0] GE hip[hipindx].plx THEN BEGIN
        quantile= 0.
    ENDIF ELSE BEGIN
        indx= (hip[hipindx].plx-quant_range[0]) / double(quant_step)
        quantile= total(quant_dist[0:indx],/NaN)*quant_step
    ENDELSE
ENDIF

IF keyword_set(dontplot) THEN RETURN

step= (xrange[1]-xrange[0])/(nplx-1.)
plxs= dindgen(nplx)*step+xrange[0]
probs= dblarr(nplx)
FOR kk=0L, nplx-1 DO BEGIN
    back_loglike= alog(kernel_phot_plx(plxs[kk],$
                                  hip[hipindx].bvcolor,$
                                  hip[hipindx].vmag,$
                                  hip[hipindx].e_plx,data,$
                                  width=width,kernel=kernel))
    probs[kk]= exp(iso_phot_plx(plxs[kk],$
                                     hip[hipindx].bvcolor,$
                                     hip[hipindx].vmag,$
                                     hip[hipindx].e_plx,thisiso,$
                                     alpha=alpha,back_loglike=back_loglike, $
                                     sspspread=sspspread)*(exp(grouplogpost[group,0])))*exp(back_loglike*(1-exp(grouplogpost[group,0])))
ENDFOR
probs= probs/(total(probs)*step)

IF ~keyword_set(yrange) THEN yrange=[0.,1.1*max(probs)]
IF keyword_set(noxlabels) THEN BEGIN
    IF keyword_set(noylabels) THEN BEGIN
        xtickname=REPLICATE(' ',30)
        ytickname=REPLICATE(' ',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, $
          xtickname=xtickname, ytickname=ytickname, _EXTRA=KeywordsForPlot
    ENDIF ELSE BEGIN
        xtickname=REPLICATE(' ',30)
        ytickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, $
          ytitle='p(\pi)', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDELSE
ENDIF ELSE BEGIN
    IF keyword_set(noylabels) THEN BEGIN
        ytickname=REPLICATE(' ',30)
        xtickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, $
          xtitle='\pi [mas]', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDIF ELSE BEGIN
        xtickname=REPLICATE('',30)
        ytickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, xtitle='\pi [mas]', $
          ytitle='p(\pi)', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDELSE
ENDELSE
djs_oplot, plxs, probs
;;Calculate confidence intervals
calc_conf_int, probs, confidence, cll=cll, cul=cul, /multimodal
cll= xrange[0]+step*cll
cul= xrange[0]+step*cul
oplotbarx, [cll,cul], color=djs_icolor('gray'), $
  xtickname=xtickname, ytickname=ytickname
oplotbarx, hip[hipindx].plx, thick=2, $
xtickname=xtickname, ytickname=ytickname
;;Put labels on
IF keyword_set(hiplabel) THEN BEGIN
    idposx= xrange[1] - (xrange[1]-xrange[0])*(.4)
    idposy= yrange[0] + (yrange[1]-yrange[0])*0.95
    identify= 'HIP'+strtrim(string(hip[hipindx].hip),2)
    legend, [identify], pos=[idposx,idposy], box=0, charsize=charsize
ENDIF
END

