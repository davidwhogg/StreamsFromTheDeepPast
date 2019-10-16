;+
;   NAME:
;      plot_radecpatch
;   PURPOSE:
;      plot the probability distribution for the radial velocity in a
;      path of the sky and compare it to the observed distribution
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      plotfilename - filename to plot to
;      basedir - basedirectory for density estimate (ends in /)
;      mplot  - multiple plots array (e.g [0,1,2])
;      seed   - seed for random
;      ra_span - box around central ra of this size
;      dec_span - box around central dec of this size
;      xsize    - xsize for k_print
;      ysize    - ysize for k_print
;      nsamples - number of samples to use to plot the dist.
;      rrange   - v_r range to plot
;      yrange   - yrange to plot
;   KEYWORDS:
;      /QA     - Print QA info
;      meandist- instead of calculating the expected distribution at
;                the center of the box, calculate it at the mean ra
;                and dec
;      median  - instead of calculating the expected distribution at
;                the center of the box, calculate it at the median ra
;                and dec
;   OUTPUT:
;      plot with the distribution(s)
;   REVISION HISTORY:
;      2009-01-16 - Written Bovy (NYU)
;-
PRO PLOT_RADECPATCH, K=K, w=w, sample=sample, plotfilename=plotfilename, $
                     basedir=basedir, QA=QA, mplot=mplot, seed=seed, $
                     ra_span=ra_span, dec_span=dec_span, xsize=xsize, $
                     ysize=ysize, nsamples= nsamples, rrange=rrange, $
                     yrange=yrange, meandist=meandist, median=median

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(seed) THEN seed= 4L
IF ~keyword_set(ra_span) THEN ra_span=40;;Box around the patch center of this size
IF ~keyword_set(dec_span) THEN dec_span= 20
IF ~keyword_set(xsize) THEN xsize=3.375;;Actual size of ApJ column!
IF ~keyword_set(ysize) THEN ysize= 3.375;;Careful with this
IF ~keyword_set(nsamples) THEN nsamples= 1000L
IF ~keyword_set(rrange) THEN rrange= [-85.,85.]
IF ~keyword_set(yrange) THEN yrange=[0,0.044]
IF ~keyword_set(mplot) THEN mplot= [0,1,2]
d= 3

;;# of patches
npatch= mplot[1] * mplot[2] > 2
IF npatch EQ 2 THEN mplot= [0,1,2]

;;ra and decs of the patches, NEP, SEP, and random ones
patches= dblarr(2,npatch)
patches[*,0]= [270.,66.5];;North ecliptic pole in deg
;patches[*,1]= [90.,-66.5];;South ecliptic pole in deg
FOR ii= 1, npatch-1 DO patches[*,ii]= [randomu(seed)*(360.-ra_span)$
                                       +ra_span/2.,$
                                       randomu(seed)*(90.-dec_span)+$
                                       dec_span/2.]

;;Create the string to restore the density distribution
densityfile=basedir+'fitv'
densityfile+= strtrim(string(K),2)
densityfile+='_bestV_V'
densityfile+=strtrim(string(w,format='(f4.1)'),2)
densityfile+= '_sample'
densityfile+= strtrim(string(sample),2)
densityfile+= '.sav'
IF file_test(densityfile) THEN BEGIN
    IF keyword_set(QA) THEN splog, 'reading '+densityfile
    restore, densityfile
ENDIF ELSE BEGIN
    splog, 'Density estimate not found in file '+densityfile
    splog, 'Returning...'
    RETURN
ENDELSE

;;Restore gcs
genevafilename= '/global/data/rv/gcs/gcs.sav'
IF file_test(genevafilename) THEN BEGIN
    restore, genevafilename
ENDIF ELSE BEGIN
    splog, 'Error: GCS-savefile does not exist in ', $
      genevafilename
    splog, 'GCS-savefile must exist, returning...'
    RETURN
ENDELSE
ngcs= n_elements(gcs.hip)

IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=ysize

!p.multi=mplot
!p.charsize=.5
!p.charthick=2

FOR ii=0L, npatch-1 DO BEGIN
    ;;Find stars in this patch
    indx_patch1= where((gcs.ra GE patches[0,ii]-ra_span/2. AND $
                       gcs.ra LE patches[0,ii]+ra_span/2. AND $
                       gcs.dec GE patches[1,ii]-dec_span/2. AND $
                       gcs.dec LE patches[1,ii]+dec_span/2.))
    indx_patch2= where(((180.-gcs.ra) GE patches[0,ii]-ra_span/2. AND $
                       (180.-gcs.ra) LE patches[0,ii]+ra_span/2. AND $
                       -gcs.dec GE patches[1,ii]-dec_span/2. AND $
                       -gcs.dec LE patches[1,ii]+dec_span/2.) OR $
                     ((180.+gcs.ra) GE patches[0,ii]-ra_span/2. AND $
                       (180.+gcs.ra) LE patches[0,ii]+ra_span/2. AND $
                       -gcs.dec GE patches[1,ii]-dec_span/2. AND $
                       -gcs.dec LE patches[1,ii]+dec_span/2.))
    n_in_patch= n_elements(indx_patch1)+n_elements(indx_patch2)
    ;;First find the rotation from (U,V,W) to (vr,vl,vb) and rotate
    ;;the velocity distribution accordingly
    IF keyword_set(QA) THEN splog, "mean ra = "+strtrim(string((moment(gcs[indx_patch1].ra))[0],format='(F5.1)'),2)
    IF keyword_set(QA) THEN splog, "mean dec = "+strtrim(string((moment(gcs[indx_patch1].dec))[0],format='(F5.1)'),2)
    IF keyword_set(QA) THEN splog, "skewness ra = "+strtrim(string((moment(gcs[indx_patch1].ra))[2],format='(F5.1)'),2)
    IF keyword_set(QA) THEN splog, "skewness dec = "+strtrim(string((moment(gcs[indx_patch1].dec))[2],format='(F5.1)'),2)
    IF keyword_set(meandist) THEN BEGIN
        radec_to_lb, (moment(gcs[indx_patch1].ra))[0],(moment(gcs[indx_patch1].dec))[0],0,0,ll,bb,pmll,pmbb,1991.25 ;;Hipparcos epoch
    ENDIF ELSE IF keyword_set(median) THEN BEGIN
        radec_to_lb, (median(gcs[indx_patch1].ra)),(median(gcs[indx_patch1].dec)),0,0,ll,bb,pmll,pmbb,1991.25 ;;Hipparcos epoch
    ENDIF ELSE BEGIN
        radec_to_lb, patches[0,ii],patches[1,ii],0,0,ll,bb,pmll,pmbb,1991.25 ;;Hipparcos epoch
    ENDELSE
    ll= ll*!DPI/1.80D2
    bb= bb*!DPI/1.80D2
    projection= dblarr(d,d)
    projection[0,0]=  cos(bb)*cos(ll)
    projection[1,0]= -sin(ll)
    projection[2,0]= -sin(bb)*cos(ll)
    projection[0,1]=  cos(bb)*sin(ll)
    projection[1,1]=  cos(ll)
    projection[2,1]= -sin(bb)*sin(ll)
    projection[0,2]=  sin(bb)
    projection[1,2]=  0
    projection[2,2]=  cos(bb)    
    projection= transpose(projection)
    
    rmean= dblarr(d,K)
    rcovar= dblarr(d,d,K)
    FOR kk=0L, K-1 DO BEGIN
        rmean[*,kk]= projection##mean[*,kk]
        rcovar[*,*,kk]= (projection##covar[*,*,kk])##transpose(projection)
        rcovar[*,*,kk]= 0.5D * (rcovar[*,*,kk] + transpose(rcovar[*,*,kk]))
    ENDFOR
    rmean_marg= dblarr(1,K)
    rcovar_marg_noerr= dblarr(1,K)
    FOR kk=0L,K-1 DO BEGIN
        rmean_marg[0,kk] = rmean[0,kk]
        rcovar_marg_noerr[0,kk]= rcovar[0,0,kk]
    ENDFOR
    grid= nsamples+1
    step= double(rrange[1]-rrange[0])/nsamples
    rs= dindgen(grid)*step + rrange[0]
    ds_marg= dblarr(grid)
    FOR kk=0L, grid-1 DO ds_marg[kk]= $
      oned_sum_gaussians(rs[kk],rmean_marg,rcovar_marg_noerr,amp)

    ;;Calculate histogram, well, not really
    npix= 25;ceil(0.3*sqrt(n_in_patch)) > 10
    binwidth= (rrange[1]-rrange[0])/double(npix+1);;Conforms to the way hogg_plothist works
    norm= total(ds_marg,/NaN)*step
    weight = make_array(n_in_patch,/float,Value=norm/float(n_in_patch))

    ;;Plot
    xtitle= 'v_r [km s^{-1}]'
    ytitle= 'P(v_r!10#!X\alpha,\delta)'
    hist_this= dblarr(n_in_patch)
    hist_this[0:n_elements(indx_patch1)-1]= gcs[indx_patch1].vr
    hist_this[n_elements(indx_patch1):n_in_patch-1]= -gcs[indx_patch2].vr
    IF (ii MOD mplot[1]) EQ 0 THEN BEGIN
        IF ii/mplot[1] EQ mplot[2]-1 THEN BEGIN
            djs_plot, rs, ds_marg, xrange=rrange, yrange=yrange, $
              xtitle=xtitle, ytitle=ytitle
            hogg_plothist, hist_this, npix=npix, /overplot, $
              weight=weight
        ENDIF ELSE BEGIN
            djs_plot, rs, ds_marg, xrange=rrange, yrange=yrange, $
              xtickname=REPLICATE(' ',30), ytitle=ytitle
            hogg_plothist, hist_this, npix=npix, /overplot, $
              weight=weight
        ENDELSE
    ENDIF ELSE BEGIN
        IF ii/mplot[1] EQ mplot[2]-1 THEN BEGIN
            djs_plot, rs, ds_marg, xrange=rrange, yrange=yrange, $
              ytickname=REPLICATE(' ',30), xtitle=xtitle
            hogg_plothist, hist_this, npix=npix, /overplot, $
              weight=weight
        ENDIF ELSE BEGIN
            djs_plot, rs, ds_marg, xrange=rrange, yrange=yrange, $
              ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
            hogg_plothist, hist_this, npix=npix, /overplot, $
              weight=weight
        ENDELSE
    ENDELSE
    ;;put labels
    lstring= textoidl('\alpha = ')+strtrim(string(patches[0,ii],$
                                                  format='(F5.1)'))+$
      textoidl('^\circ')
    lstring2= textoidl('\delta = ')+strtrim(string(patches[1,ii],$
                                                  format='(F5.1)'))+$
      textoidl('^\circ')
    lstring3= strtrim(string(n_in_patch),2)+' stars'

    lposx= rrange[0] + (rrange[1]-rrange[0])*0.38
    lposy= yrange[0] + (yrange[1]-yrange[0])*1.05;0.95
    
    l2posx= rrange[0] + (rrange[1]-rrange[0])*0.38
    l2posy= yrange[0] + (yrange[1]-yrange[0])*.95;0.85

    idposx= rrange[0] + (rrange[1]-rrange[0])*(-0.12)
    idposy= yrange[0] + (yrange[1]-yrange[0])*1.05
    
    legend, [lstring], pos=[lposx,lposy], box=0, charsize=.6
    legend, [lstring2], pos=[l2posx,l2posy], box=0, charsize=.6
    legend, [lstring3], pos=[idposx,idposy], box=0, charsize=.6

    IF keyword_set(QA) THEN BEGIN
        splog, "This patch has "+strtrim(string(n_in_patch),2)+" members"
    ENDIF


ENDFOR


IF keyword_set(plotfilename) THEN k_end_print

END
