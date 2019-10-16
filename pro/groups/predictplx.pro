;+
;   NAME:
;      predictplx
;   PURPOSE:
;      plot the prediction of the photometric parallax from the
;      background model
;   INPUT:
;      id - HIP number
;      width - width parameter of the kernel
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      datafilename - filename with the data
;      nplx - number of plxs to evaluate the probability at
;      xrange - range of plx
;      yrange - range of p(plx)
;      confidence - confidence interval (default .95) on the
;                   plot
;   KEYWORDS:
;      loo - leave the star out of the predictive set      
;      hiplabel - put Hipparcos identifier as a label
;      byindx - if set, the id is the indx in the hipparcos struct
;      noxlabels - no x labels
;      noylabels - no y labels
;      dontplot - don't plot
;      labelposition - if set, shift the label to the left
;      usegcs - if set, the star is in the GCS catalog
;   OUTPUT:
;      plot
;      quantile at which the true parallax is found if asked for
;   REVISION HISTORY:
;      2009-11-05 - Written - Bovy (NYU)
;-
PRO PREDICTPLX, id, width=width, loo=loo, kernel=kernel, $
                datafilename=datafilename, xrange=xrange, yrange=yrange, $
                nplx=nplx, confidence=confidence, hiplabel=hiplabel,$
                _EXTRA=KeywordsForPlot, byindx=byindx, $
                noxlabels=noxlabels, noylabels=noylabels, quantile=quantile, $
                dontplot=dontplot, labelposition=labelposition, charsize=charsize, $
                usegcs=usegcs
IF ~keyword_set(charsize) THEN charsize=1
IF ~keyword_set(xrange) THEN xrange=[.05,60.]
IF ~keyword_set(nplx) THEN nplx= 1001
IF ~keyword_set(confidence) THEN confidence= .95
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN thisdatafilename='../lsr2/hip2-aumer.sav' ELSE thisdatafilename=datafilename
restore, filename=thisdatafilename
nhip= n_elements(hip.hip)
IF keyword_set(byindx) THEN thisid= hip[id].hip ELSE thisid= id
IF keyword_set(loo) THEN BEGIN
    indx= where(hip.hip NE thisid)
ENDIF ELSE BEGIN
    indx= lindgen(nhip)
ENDELSE
nindx= n_elements(indx)
data= dblarr(nindx,2)
data[*,0]= hip[indx].bvcolor
data[*,1]= hip[indx].vmag-5.*alog10(100./hip[indx].plx)
hipindx= (where(hip.hip EQ thisid))[0]

IF keyword_set(usegcs) THEN BEGIN
    thisdatafilename='gcs.sav'
    restore, filename=thisdatafilename
    hip=gcs
    IF keyword_set(byindx) THEN thisid= hip[id].hip ELSE thisid= id
    hipindx= (where(hip.hip EQ thisid))[0]
ENDIF

;;Calculate quantile
IF arg_present(quantile) THEN BEGIN
    ;;Sample the distribution finer, on a wider grid
    quant_grid= 10*nplx+1
    quant_range=[0.01,200.]
    quant_step= double(quant_range[1]-quant_range[0])/(10*nplx)
    quant_plxs= dindgen(quant_grid)*quant_step + quant_range[0]
    quant_dist= dblarr(quant_grid)
    FOR kk=0L, quant_grid-1 DO BEGIN
        quant_dist[kk]= kernel_phot_plx(quant_plxs[kk],$
                                        hip[hipindx].bvcolor,$
                                        hip[hipindx].vmag,$
                                        hip[hipindx].e_plx,data,$
                                        width=width,kernel=kernel)
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
probs= kernel_phot_plx(plxs,hip[hipindx].bvcolor,hip[hipindx].vmag,$
                      hip[hipindx].e_plx,data,width=width,kernel=kernel)
probs= probs/(total(probs)*step)

IF ~keyword_set(yrange) THEN yrange=[0.,1.1*max(probs)]
IF keyword_set(noxlabels) THEN BEGIN
    IF keyword_set(noylabels) THEN BEGIN
        xtickname=REPLICATE(' ',30)
        ytickname=REPLICATE(' ',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, charsize=charsize,$
          xtickname=xtickname, ytickname=ytickname, _EXTRA=KeywordsForPlot
    ENDIF ELSE BEGIN
        xtickname=REPLICATE(' ',30)
        ytickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, charsize=charsize,$
          ytitle='p(\pi)', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDELSE
ENDIF ELSE BEGIN
    IF keyword_set(noylabels) THEN BEGIN
        ytickname=REPLICATE(' ',30)
        xtickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, charsize=charsize,$
          xtitle='\pi [mas]', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDIF ELSE BEGIN
        xtickname=REPLICATE('',30)
        ytickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, xtitle='\pi [mas]', $
          ytitle='p(\pi)', xtickname=xtickname, ytickname=ytickname, charsize=charsize,$
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
    IF ~keyword_set(labelposition) THEN BEGIN
        idposx= xrange[1] - (xrange[1]-xrange[0])*(0.4)
        idposy= yrange[0] + (yrange[1]-yrange[0])*0.95
        pos=[idposx,idposy]
    ENDIF ELSE BEGIN
        idposx= xrange[1] - (xrange[1]-xrange[0])*(0.5)
        idposy= yrange[0] + (yrange[1]-yrange[0])*0.95
        pos=[idposx,idposy]
    ENDELSE
    identify= 'HIP'+strtrim(string(hip[hipindx].hip),2)
    legend, [identify], pos=pos, box=0, charsize=charsize
ENDIF
END
