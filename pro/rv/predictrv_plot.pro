;+
;   NAME:
;      predictrv_plot
;
;   PURPOSE:
;      plot the prediction for the radial velocity of a star
;
;   CALLING SEQUENCE:
;      predictrv_plot, xmean, xcovar, xamp, ll, bb,
;      confidence=confidence, measuredrv= measuredrv, measuredvar=
;      measurederr, color=color, rrange= rrange, sample= sample,
;      confintwidth=confintwidth, /calconly, /one_interval
;
;   INPUTS:
;      xmean, xcovar, xamp - parameters of the underlying distribution
;      ll, bb              - galactic coordinates (in deg)
;
;   OPTIONAL INPUTS:
;      confidence          - confidence interval limit (default .95)
;      measuredrv          - measured value of the radial velocity
;      measurederr         - measurement error
;      rrange              - r-range to use
;      sample              - number of samples to take
;
;   OUTPUT:
;      confintwidth        - confidence interval width
;
;   KEYWORDS:
;      color               - plot in color? default no.
;      calconly            - only calculate (the confidence interval)
;      one_interval        - specify that you
;                            want a single confidence interval
;
;   REVISION HISTORY:
;      2008-09-02 - Written Bovy
;-
PRO CALC_CONF_INT, ds, confidence, cul=cul, cll=cll, confintwidth=confintwidth, $
                   one_interval=one_interval
;;Search for highest likelihood interval that holds confidence% of the
;;probability (~plot_projected_gaussians contouring algo)
IF ~keyword_set(one_interval) THEN BEGIN
    cumindx= reverse(sort(ds))
    grid= n_elements(ds)
    cumds= dblarr(grid)
    cumds[cumindx]= total(ds[cumindx],/cumulative)
;    cumds= cumds*step
    cumds= cumds/max(cumds)
    confindx= where(cumds LE confidence)
    nconf= n_elements(confindx)
    cul= bytarr(nconf)
    cll= bytarr(nconf)
    IF (cumds[0] LE confidence) THEN BEGIN
        print, cumds[0]
        splog, 'ERROR: interval is not wide enough, increase size and run again'
        return
    ENDIF
    cll[0]= 1B
    cul[nconf-1]= 1B
    up= 1B
    FOR ii= 1, nconf-1 DO BEGIN
        IF confindx[ii] NE (confindx[ii-1]+1) THEN BEGIN
            cul[ii-1]= 1B
            cll[ii]= 1B
        ENDIF
    ENDFOR
    cul= confindx[where(cul EQ 1B)]
    cll= confindx[where(CLL EQ 1B)]
    IF arg_present(confintwidth) THEN confintwidth= nconf
ENDIF ELSE BEGIN
;;old interval method
;;Search for minimal interval that holds 95% probability
;;For every point, look ahead to where 95% is reached
    cds= total(ds,/cumulative,/double)/total(ds,/double)
    candidates= where(cds LE 1D0-confidence)
    ncand= n_elements(candidates)
    pairs= lonarr(3,ncand)
    FOR ii=0, ncand-1 DO BEGIN
        candidates2= where(cds GT cds[ii])
        ncand2= n_elements(candidates2)
        tmpcds= cds[candidates2]- cds[candidates[ii]]
        pairs[0,ii]= ii
        pairs[1,ii]= min(where(tmpcds GE confidence))+ii
        pairs[2,ii]= pairs[1,ii]-pairs[0,ii]
    ENDFOR
    mindist= min(pairs[2,*],minii)
    cll= pairs[0,minii]
    cul= pairs[1,minii]
;;Return the interval?
    IF arg_present(confintwidth) THEN confintwidth= cul-cll
ENDELSE
END
PRO PREDICTRV_PLOT, xmean, xcovar, xamp, ll, bb, confidence=confidence, $
                    measuredrv= measuredrv, measurederr= measurederr, $
                    color=color, rrange=rrange, sample=sample, $
                    confintwidth=confintwidth, calconly= calconly, $
                    one_interval= one_interval

ON_ERROR, 2

IF ~keyword_set(confidence) THEN confidence = .95D
IF ~keyword_set(rrange) THEN rrange= [-200D,200D]
IF ~keyword_set(sample) THEN sample= 1000L


;Construct 1D means and covariances
llrad= ll/180.0*!DPI
bbrad= bb/180.0*!DPI
rr= [[cos(llrad)*cos(bbrad)],[sin(llrad)*cos(bbrad)],[sin(bbrad)]]
nmeans= n_elements(xamp)
rmean= dblarr(1,nmeans)
rcovar= dblarr(1,nmeans)
FOR ii= 0L, nmeans-1 DO BEGIN
    rmean[0,ii]=rr#xmean[*,ii]
    rcovar[0,ii]= rr#xcovar[*,*,ii]#transpose(rr)
ENDFOR

;;sample distribution
grid= sample+1
step= double(rrange[1]-rrange[0])/sample
rs= dindgen(grid)*step + rrange[0]
ds= dblarr(grid)
FOR kk=0L, grid-1 DO ds[kk]= oned_sum_gaussians(rs[kk],rmean, rcovar, xamp)

calc_conf_int, ds, confidence, cll=cll, cul=cul, one_interval=one_interval, $
  confintwidth=confintwidth

cll= rrange[0]+step*cll
cul= rrange[0]+step*cul
IF arg_present(confintwidth) THEN confintwidth= confintwidth*step

IF keyword_set(calconly) THEN RETURN
    

yrange= [0,max(ds)]

;;If measurederr is given, convolve the distribution with the errors
;;and compute confidence limits again
;;Don't do this if the measurederror is smaller than the sampling
IF keyword_set(measurederr) AND keyword_set(measuredrv) THEN BEGIN
    IF measurederr LT step THEN BEGIN
        splog, 'Warning: step is larger than measurederr'
        splog, 'Skipping convolution'
;;plot
        djs_plot, rs, ds, xtitle='v_r [km s^{-1}]', ytitle= 'P(v_r)', yrange=yrange
    ENDIF ELSE BEGIN
        kernel= dindgen(sample+1)*step + rrange[0]
        kernel= 1/sqrt(2.*!DPI)/measurederr*exp(-kernel^2/(2.*measurederr^2))
        convds= convol(ds,kernel,/edge_truncate)*step
        IF max(convds) GT yrange[1] THEN yrange[1]= max(convds)
        IF keyword_set(color) THEN djs_plot, rs, convds, $
          color=djs_icolor('cyan green'), xtitle='v_r [km s^{-1}]', $
          ytitle= 'P(v_r)', yrange=yrange $
        ELSE djs_plot, rs, convds, $
          color=djs_icolor('gray'), xtitle='v_r [km s^{-1}]', $
          ytitle= 'P(v_r)', yrange=yrange
        ;;also calculate confidence limits for the conv-distribution
        calc_conf_int, convds, confidence, cll=cllconv, cul=culconv, one_interval=one_interval
        culconv= rrange[0]+step*culconv
        cllconv= rrange[0]+step*cllconv
;        cds= total(convds,/cumulative,/double)/total(convds,/double)
;        candidates= where(cds LE 1D0-confidence)
;        ncand= n_elements(candidates)
;        pairs= lonarr(3,ncand)
;        FOR ii=0, ncand-1 DO BEGIN
;            candidates2= where(cds GT cds[ii])
;            ncand2= n_elements(candidates2)
;            tmpcds= cds[candidates2]- cds[candidates[ii]]
;            pairs[0,ii]= ii
;            pairs[1,ii]= min(where(tmpcds GE confidence))+ii
;            pairs[2,ii]= pairs[1,ii]-pairs[0,ii]
;        ENDFOR
;        mindist= min(pairs[2,*],minii)
;        cllconv= rrange[0]+step*pairs[0,minii]
;        culconv= rrange[0]+step*pairs[1,minii]
        IF keyword_set(color) THEN oplotbarx, [cllconv,culconv], $
          color=djs_icolor('cyan green') $
        ELSE oplotbarx, [cllconv,culconv], color=djs_icolor('gray')
    ENDELSE
ENDIF ELSE BEGIN
    djs_plot, rs, ds, xtitle='v_r [km s^{-1}]', ytitle= 'P(v_r)', yrange=yrange
ENDELSE


;;Plot
djs_oplot, rs, ds, xtitle='v_r [km s^{-1}]', ytitle= 'P(v_r)', yrange=yrange
IF keyword_set(color) THEN oplotbarx, [cll,cul], color=djs_icolor('red') $
ELSE oplotbarx, [cll,cul]

;;Plot measured value
IF keyword_set(measuredrv) THEN IF keyword_set(color) THEN oplotbarx, measuredrv, color=djs_icolor('blue') $
ELSE oplotbarx, measuredrv, linestyle=3

END
