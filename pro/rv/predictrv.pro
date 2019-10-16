;+
;   NAME:
;      predictrv
;   PURPOSE:
;      predict the radial velocity of a starr given the velocity
;      distribution, the locus of the star and (optionally) its
;      tangential velocity (+errors); (optionally) plot the prediction
;   CALLING SEQUENCE:
;      predictrv, xamp, xmean, xcovar, 
;   INPUT:
;      xamp, xmean, xcovar - parameters of the vlocity distribution in
;                            (U,V,W)
;      projection          - projection from (U,V,W) to (vr,vl,vb)
;   OPTIONAL INPUT:
;      vr                  - measured radial velocity
;      vlvb=[vl,vb]        - measured tangential velocity
;      vlvbc               - error covariance for tangential
;                            components
;      cvr                 - radial velocity squared error
;      confidence          - confidence interval (default .95) on the
;                            plot
;      rrange              - radial velocity range to use (default
;                            [-150,150])
;      yrange              - speaks for itself
;      sample              - number of samples to use for plotting
;      identity            - identifier string to use as a label
;      seed                - seed for random number generator
;      maxvr               - vr with max probability
;      cll                 - 95% lower limit
;      cul                 - 95% upper limit      
;   KEYWORDS:
;      noxlabels           - no x labels
;      noylabels           - no y labels
;      putloglikelabels    - put the loglike as a label in the plot
;      putentropylabels    - put the entropies
;                            of the distributions in the plot
;      noerr               - plot the predictions not convolved
;                            with the errors
;   OUTPUT:
;      loglike             - logarithm of the probability of the
;                            radial velocity given the model and the
;                            tangential velocities (or just the model
;                            if the tan. velocities are not given);
;                            Not returned if the radial velocity is
;                            not given
;      altloglike          - logarithm of the proability of the
;                            rad. vel. given the model
;      entropy             - entropy of the conditional distribution
;      altentropy          - entropy of the marginal distribution
;      quantile            - quantile of the conditional dist. @ which
;                            the true velocity is
;   KEYWORDS:
;      plot                - send a plot of the prediction to the
;                            plotting device (?)
;      color               - plot in color
;   REVISION HISTORY:
;      2008-10-26 - Written Bovy (NYU)
;-
;PRO CALC_CONF_INT, ds, confidence, cul=cul, cll=cll, confintwidth=confintwidth, $
;                   one_interval=one_interval
;;;Search for highest likelihood interval that holds confidence% of the
;;;probability (~plot_projected_gaussians contouring algo)
;IF ~keyword_set(one_interval) THEN BEGIN
;    cumindx= reverse(sort(ds))
;    grid= n_elements(ds)
;    cumds= dblarr(grid)
;    cumds[cumindx]= total(ds[cumindx],/cumulative)
;;    cumds= cumds*step
;    cumds= cumds/max(cumds)
;    confindx= where(cumds LE confidence)
;    nconf= n_elements(confindx)
;    cul= bytarr(nconf)
;    cll= bytarr(nconf)
;    IF (cumds[0] LE confidence) THEN BEGIN
;        print, cumds[0]
;        splog, 'ERROR: interval is not wide enough, increase size and run again'
;        return
;    ENDIF
;    cll[0]= 1B
;    cul[nconf-1]= 1B
;    up= 1B
;    FOR ii= 1, nconf-1 DO BEGIN
;        IF confindx[ii] NE (confindx[ii-1]+1) THEN BEGIN
;            cul[ii-1]= 1B
;            cll[ii]= 1B
;        ENDIF
;    ENDFOR
;    cul= confindx[where(cul EQ 1B)]
;    cll= confindx[where(CLL EQ 1B)]
;    IF arg_present(confintwidth) THEN confintwidth= nconf
;ENDIF ELSE BEGIN
;;;old interval method
;;;Search for minimal interval that holds 95% probability
;;;For every point, look ahead to where 95% is reached
;    cds= total(ds,/cumulative,/double)/total(ds,/double)
;    candidates= where(cds LE 1D0-confidence)
;    ncand= n_elements(candidates)
;    pairs= lonarr(3,ncand)
;    FOR ii=0, ncand-1 DO BEGIN
;        candidates2= where(cds GT cds[ii])
;        ncand2= n_elements(candidates2)
;        tmpcds= cds[candidates2]- cds[candidates[ii]]
;        pairs[0,ii]= ii
;        pairs[1,ii]= min(where(tmpcds GE confidence))+ii
;        pairs[2,ii]= pairs[1,ii]-pairs[0,ii]
;    ENDFOR
;    mindist= min(pairs[2,*],minii)
;    cll= pairs[0,minii]
;    cul= pairs[1,minii]
;;;Return the interval?
;    IF arg_present(confintwidth) THEN confintwidth= cul-cll
;ENDELSE
;END
PRO PREDICTRV, xamp, xmean, xcovar, projection, vr=vr, vlvb=vlvb, $
               cvlvb=cvlvb,cvr=cvr, $
               confidence=confidence, rrange=rrange, sample=sample, $
               loglike=loglike, altloglike=altloglike, plot=plot, $
               color=color, entropy=entropy, altentropy=altentropy,$
               yrange=yrange, noxlabels=noxlabels, noylabels=noylabels, $
               identify=identify, putloglikelabels=putloglikelabels, $
               seed=seed, putentropylabels=putentropylabels, $
               quantile=quantile, maxvr=maxvr, cul=cul, cll=cll, $
               _EXTRA=KeywordsForPlot, noerr=noerr


;;Defaults
IF ~keyword_set(confidence) THEN confidence = .95D
IF ~keyword_set(rrange) THEN rrange= [-150D,150D]
IF ~keyword_set(sample) THEN sample= 1000L
IF ~keyword_set(seed) THEN seed=1L

IF ~keyword_set(cvr) THEN cvr=0D

;;keyword abbreviations in IDL are annoying
vlvbc= cvlvb
vrc= cvr

;;Perhaps do the following for lists of inputs?

;;Number of Gaussians used
K= n_elements(xamp)
d= 3

twopiterm=-0.5*alog(2.*!DPI)

;;First rotate from (U,V,W) to (vr,vl,vb)
rmean= dblarr(d,K)
rcovar= dblarr(d,d,K)
FOR kk=0L, K-1 DO BEGIN
    rmean[*,kk]= projection##xmean[*,kk]
    rcovar[*,*,kk]= (projection##xcovar[*,*,kk])##transpose(projection)
    rcovar[*,*,kk]= 0.5D * (rcovar[*,*,kk] + transpose(rcovar[*,*,kk]))
ENDFOR

;;Then find the following (1) marginalized without vr errors (2)
;;marginalized with errors (3) conditioned on vl,vb without vr errors
;;(4) conditioned on vl,vb with vr errors

;;(1)+(2)
rmean_marg= dblarr(1,K)
rcovar_marg_noerr= dblarr(1,K)
rcovar_marg= dblarr(1,K)
loglike_tmp= dblarr(K)
FOR kk=0L,K-1 DO BEGIN
    rmean_marg[0,kk] = rmean[0,kk]
    rcovar_marg_noerr[0,kk]= rcovar[0,0,kk]
    rcovar_marg[0,kk]= rcovar[0,0,kk]
    IF keyword_set(cvr) THEN rcovar_marg[0,kk]+= vrc
    IF keyword_set(vr) THEN BEGIN
        delta= vr - rmean_marg[0,kk]
        IF keyword_set(cvr) THEN loglike_tmp[kk]= alog(xamp[kk]) + twopiterm $
          - 0.5 * alog(rcovar_marg[0,kk]) -0.5 * delta^2 / rcovar_marg[0,kk] $ 
        ELSE loglike_tmp[kk]= alog(xamp[kk]) + twopiterm - 0.5 * alog(rcovar_marg_noerr[0,kk]) $
          -0.5 * delta^2 / rcovar_marg_noerr[0,kk]
    ENDIF
ENDFOR
IF keyword_set(vr) THEN loglike= logsum(loglike_tmp,/double)

;;(3)+(4)
IF keyword_set(vlvb) THEN BEGIN
    IF keyword_set(vr) THEN altloglike= loglike
    ramp_cond= dblarr(K);first will hold log like alpha_j N(vlvb|mlb,Tlb)
    rmean_cond= dblarr(1,K)
    rcovar_cond_noerr= dblarr(1,K)
    rcovar_cond= dblarr(1,K)
    FOR kk=0L, K-1 DO BEGIN
        mlb= rmean[1:2,kk]
        Tlb= rcovar[1:2,1:2,kk]
        IF keyword_set(vlvbc) THEN Tlb += vlvbc
        Vrlb= rcovar[1:2,0,kk]
        delta= (vlvb-mlb)
        rmean_cond[0,kk]= rmean[0,kk] + (Vrlb##special_invert(Tlb))##delta
        rcovar_cond_noerr[0,kk]= rcovar[0,0,kk] - (Vrlb##special_invert(Tlb))##transpose(Vrlb)
        rcovar_cond[0,kk]= rcovar_cond_noerr[0,kk]
        IF keyword_set(cvr) THEN rcovar_cond[0,kk]+= vrc
        ramp_cond[kk]= alog(xamp[kk]) - 0.5D * alog(determ(Tlb,/double,/check)) - $
          0.5 * transpose(delta)##special_invert(Tlb)##delta
        IF keyword_set(vr) THEN BEGIN
            deltavr= vr - rmean_cond[0,kk]
            IF keyword_set(cvr) THEN loglike_tmp[kk]= ramp_cond[kk] + twopiterm $
              - 0.5 * alog(rcovar_cond[0,kk]) -0.5 * deltavr^2/rcovar_cond[0,kk] $
            ELSE loglike_tmp[kk]= alog(ramp_cond[kk]) + twopiterm $
              - 0.5 * alog(rcovar_cond_noerr[0,kk]) -0.5 * deltavr^2/rcovar_cond_noerr[0,kk]
        ENDIF
    ENDFOR
    normalize= logsum(ramp_cond,/double)
    ramp_cond -= normalize
    ramp_cond = exp(ramp_cond)
    IF keyword_set(vr) THEN loglike= logsum(loglike_tmp,/double) - normalize
ENDIF

;;Calculate quantile
IF arg_present(quantile) THEN BEGIN
    IF ~keyword_set(vr) THEN BEGIN
        splog, "Can't calculate the quantile if the true velocity isn't given"
        quantile= -1
    ENDIF ELSE BEGIN
        ;;Sample the distribution
        quant_grid= 10*sample+1
        quant_step= double(rrange[1]-rrange[0])/(10*sample)
        quant_rs= dindgen(quant_grid)*quant_step + rrange[0]
        quant_dist_cond= dblarr(quant_grid)
        FOR kk=0L, quant_grid-1 DO quant_dist_cond[kk]= oned_sum_gaussians(quant_rs[kk],rmean_cond,rcovar_cond,ramp_cond)
        ;;Normalize to one
        quant_dist_cond=quant_dist_cond/(total(quant_dist_cond,/NaN)*quant_step)
        ;;Find where vr is
        IF rrange[1] LT vr THEN BEGIN
            quantile= 1.
        ENDIF ELSE IF rrange[0] GE vr THEN BEGIN
            quantile= 0.
        ENDIF ELSE BEGIN
            indx= (vr-rrange[0]) / double(quant_step)
            quantile= total(quant_dist_cond[0:indx],/NaN)*quant_step
        ENDELSE
    ENDELSE
ENDIF

;;Calculate entropies
IF (arg_present(entropy) OR keyword_set(plot)) THEN BEGIN
    calc_entropy_sum_gaussians, ramp_cond, rmean_cond, rcovar_cond, $
      nboot=100, nsamplings=1000, entropy=entropy, sigma_H=sigma_H, seed=seed
    IF sigma_H/entropy GT 0.001 THEN BEGIN
        splog, 'Warning: entropy sampling has variation GT 0.001'
        splog, 'Trying with 10x more samples...'
        calc_entropy_sum_gaussians, ramp_cond, rmean_cond, rcovar_cond, $
          nboot=100, nsamplings=10000, entropy=entropy, sigma_H=sigma_H, $
          seed=seed
        IF sigma_H/entropy GT 0.001 THEN $
          splog, 'Warning: entropy sampling *still* has variation GT 0.001'
    ENDIF
ENDIF

IF (arg_present(altentropy) OR keyword_set(plot)) THEN BEGIN
    calc_entropy_sum_gaussians, xamp, rmean_marg, rcovar_marg, nboot=100, $
      nsamplings=1000, entropy=altentropy, sigma_H=sigma_H, seed=seed
    IF sigma_H/altentropy GT 0.001 THEN BEGIN
        splog, 'Warning: entropy sampling has variation GT 0.001'
        splog, 'Trying with 10x more samples...'
        calc_entropy_sum_gaussians, xamp, rmean_marg, rcovar_marg, nboot=100, $
          nsamplings=10000, entropy=altentropy, sigma_H=sigma_H, seed=seed
        IF sigma_H/entropy GT 0.001 THEN $
          splog, 'Warning: entropy sampling *still* has variation GT 0.001'
    ENDIF
ENDIF

;;If plot is not set return
IF (~keyword_set(plot) AND ~arg_present(maxvr) AND ~arg_present(cul) AND ~arg_present(cll) )THEN RETURN

;;Plotting:

;;sample distribution
grid= sample+1
step= double(rrange[1]-rrange[0])/sample
rs= dindgen(grid)*step + rrange[0]
ds_marg= dblarr(grid)
ds_cond= dblarr(grid)
ds_marg_noerr= dblarr(grid)
ds_cond_noerr= dblarr(grid)
FOR kk=0L, grid-1 DO ds_marg[kk]= oned_sum_gaussians(rs[kk],rmean_marg,rcovar_marg,xamp)
FOR kk=0L, grid-1 DO ds_cond[kk]= oned_sum_gaussians(rs[kk],rmean_cond,rcovar_cond,ramp_cond)
FOR kk=0L, grid-1 DO ds_marg_noerr[kk]= oned_sum_gaussians(rs[kk],rmean_marg,rcovar_marg_noerr,xamp)
FOR kk=0L, grid-1 DO ds_cond_noerr[kk]= oned_sum_gaussians(rs[kk],rmean_cond,rcovar_cond_noerr,ramp_cond)


IF arg_present(maxvr) THEN BEGIN
    ;;Calculate where the max of the probability occurs
    max_indx= (where(ds_cond_noerr EQ max(ds_cond_noerr,/NaN)))[0]
    maxvr= rs[max_indx]
ENDIF

;;Calculate confidence intervals
IF ((arg_present(cll) OR arg_present(cul)) AND ~keyword_set(plot)) THEN BEGIN
    calc_conf_int, ds_cond, confidence, cll=cll, cul=cul
    cll= rrange[0]+step*cll
    cul= rrange[0]+step*cul
ENDIF ELSE BEGIN
    calc_conf_int, ds_cond, confidence, cll=cll, cul=cul, /multimodal
    cll= rrange[0]+step*cll
    cul= rrange[0]+step*cul
ENDELSE

IF ~keyword_set(plot) THEN RETURN

IF ~keyword_set(yrange) THEN yrange= [0,max([ds_marg,ds_cond,ds_marg_noerr,ds_cond_noerr])]

;;Then plot everything

;;Make axes
IF keyword_set(noxlabels) THEN BEGIN
    IF keyword_set(noylabels) THEN BEGIN
        xtickname=REPLICATE(' ',30)
        ytickname=REPLICATE(' ',30)
        djs_plot, [0.], [0.], xrange=rrange, yrange=yrange, $
          xtickname=xtickname, ytickname=ytickname, _EXTRA=KeywordsForPlot
    ENDIF ELSE BEGIN
        xtickname=REPLICATE(' ',30)
        ytickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=rrange, yrange=yrange, $
          ytitle= 'P(v_r)', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDELSE
ENDIF ELSE BEGIN
    IF keyword_set(noylabels) THEN BEGIN
        ytickname=REPLICATE(' ',30)
        xtickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=rrange, yrange=yrange, $
          xtitle='v_r [km s^{-1}]', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDIF ELSE BEGIN
        xtickname=REPLICATE('',30)
        ytickname=REPLICATE('',30)
        djs_plot, [0.], [0.], xrange=rrange, yrange=yrange, xtitle='v_r [km s^{-1}]', $
          ytitle= 'P(v_r)', xtickname=xtickname, ytickname=ytickname, $
          _EXTRA=KeywordsForPlot
    ENDELSE
ENDELSE

;And plot!
IF keyword_set(noerr) THEN BEGIN
    IF ~keyword_set(color) THEN djs_oplot, rs, ds_cond_noerr, $
      linestyle=2, xtickname=xtickname, ytickname=ytickname $
    ELSE djs_oplot, rs, ds_cond_noerr, color=djs_icolor('magenta'), $
      xtickname=xtickname, ytickname=ytickname
ENDIF
IF ~keyword_set(color) THEN djs_oplot, rs, ds_cond, xtickname=xtickname, $
  ytickname=ytickname $
  ELSE djs_oplot, rs, ds_cond, color=djs_icolor('black'), $
  xtickname=xtickname, ytickname=ytickname
IF keyword_set(noerr) THEN BEGIN
    IF ~keyword_set(color) THEN djs_oplot, rs, ds_marg_noerr, $
      linestyle=3, xtickname=xtickname, ytickname=ytickname $
    ELSE djs_oplot, rs, ds_marg_noerr, color=djs_icolor('yellow'),$
      xtickname=xtickname, ytickname=ytickname
ENDIF
IF ~keyword_set(color) THEN djs_oplot, rs, ds_marg, $
  color=djs_icolor('gray'), xtickname=xtickname, ytickname=ytickname $
ELSE djs_oplot, rs, ds_marg, color=djs_icolor('green'), $
  xtickname=xtickname, ytickname=ytickname


;;Plot confidence intervals
IF keyword_set(color) THEN oplotbarx, [cll,cul], color=djs_icolor('red'), xtickname=xtickname, ytickname=ytickname $
ELSE oplotbarx, [cll,cul], color=djs_icolor('gray'), xtickname=xtickname, ytickname=ytickname

;;confidence intervals for marginal distribution (DEPRECATED)
;calc_conf_int, ds_marg, confidence, cll=cll, cul=cul,/multimodal
;cll= rrange[0]+step*cll
;cul= rrange[0]+step*cul
;IF keyword_set(color) THEN oplotbarx, [cll,cul], color=djs_icolor('orange'), xtickname=xtickname, ytickname=ytickname $
;ELSE oplotbarx, [cll,cul], linestyle=5, xtickname=xtickname, ytickname=ytickname


;;Plot measured value
IF keyword_set(vr) THEN BEGIN
    IF keyword_set(color) THEN oplotbarx, vr, color=djs_icolor('blue'), thick=2, xtickname=xtickname, ytickname=ytickname $
    ELSE oplotbarx, vr, thick=2, xtickname=xtickname, ytickname=ytickname
ENDIF


;;Put labels on
idposx= rrange[0] + (rrange[1]-rrange[0])*(-0.02)
idposy= yrange[0] + (yrange[1]-yrange[0])*0.95

IF keyword_set(identify) THEN legend, [identify], $
  pos=[idposx,idposy], box=0, charsize=.8


;;Put entropies or log likelihoods as labels
IF keyword_set(putentropylabels) THEN BEGIN
    IF entropy GE 10 THEN $
      lstring= textoidl('H_1 = ')+strtrim(string(entropy,format='(f6.2)'),2) $
    ELSE $
      lstring= textoidl('H_1 = ')+strtrim(string(entropy,format='(f6.3)'),2)
    IF altentropy GE 10 THEN $
      lstring2= textoidl('H_2 = ')+strtrim(string(altentropy,format='(f6.2)'),2) $
    ELSE $
      lstring2= textoidl('H_2 = ')+strtrim(string(altentropy,format='(f6.3)'),2)
    
    IF entropy GT 99.999 THEN lstring= textoidl('H_1 > ')+strtrim(string(99.999,format='(f6.3)'),2)
    IF altentropy GT 99.999 THEN lstring2= textoidl('H_2 > ')+strtrim(string(99.999,format='(f6.3)'),2)
ENDIF ELSE IF keyword_set(putloglikelabels) THEN BEGIN
    IF loglike GE 10 THEN $
      lstring= textoidl('L_1 = ')+strtrim(string(loglike,format='(f6.2)'),2) $
    ELSE $
      lstring= textoidl('L_1 = ')+strtrim(string(loglike,format='(f6.3)'),2)
    IF altloglike GE 10 THEN $
      lstring2= textoidl('L_2 = ')+strtrim(string(altloglike,format='(f6.2)'),2) $
    ELSE $
      lstring2= textoidl('L_2 = ')+strtrim(string(altloglike,format='(f6.3)'),2)
    
    IF loglike LT -9.999 THEN lstring= textoidl('L_1 < ')+strtrim(string(-9.999,format='(f6.3)'),2)
    IF altloglike LT -9.999 THEN lstring2= textoidl('L_2 < ')+strtrim(string(-9.999,format='(f6.3)'),2)
ENDIF

lposx= rrange[0] + (rrange[1]-rrange[0])*0.55
lposy= yrange[0] + (yrange[1]-yrange[0])*0.95

l2posx= rrange[0] + (rrange[1]-rrange[0])*0.55
l2posy= yrange[0] + (yrange[1]-yrange[0])*0.85

IF keyword_set(putloglikelabels) OR keyword_set(putentropylabels) THEN BEGIN
    legend, [lstring], pos=[lposx,lposy], box=0, charsize=.8
    legend, [lstring2], pos=[l2posx,l2posy], box=0, charsize=.8
ENDIF



END
