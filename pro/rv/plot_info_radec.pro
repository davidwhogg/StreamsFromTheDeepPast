;+
;   NAME:
;      plot_info_radec
;   PURPOSE:
;      Calculate the most informative direction on the celestial
;      sphere, and plot the predicted distribution + observed
;      distribution from GCS for this direction
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      filename - filename to save the entropies
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
;      patches  - number of patches to calculate the entropy in
;   KEYWORDS:
;      QA      - print QA info
;   OUTPUT:
;   REVISION HISTORY:
;      2009-01-20 - Written Bovy (NYU)
;-
PRO PLOT_INFO_RADEC, K=K, w=w, sample=sample, filename=filename, $
                     plotfilename=plotfilename, basedir=basedir, $
                     mplot=mplot, seed=seed, ra_span=ra_span, $
                     dec_span=dec_span, xsize=xsize, ysize=ysize, $
                     nsamples=nsamples, rrange=rrange, yrange=yrange, $
                     patches=patches, QA=QA, highentropy=highentropy

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(ra_span) THEN ra_span=40;;Box around the patch center of this size
IF ~keyword_set(dec_span) THEN dec_span= 20
IF ~keyword_set(xsize) THEN xsize=3.375;;Actual size of ApJ column!
IF ~keyword_set(ysize) THEN ysize= 3.375;;Careful with this
IF ~keyword_set(nsamples) THEN nsamples= 1000L
IF ~keyword_set(rrange) THEN rrange= [-85.,85.]
IF ~keyword_set(yrange) THEN yrange=[0,0.044]
IF ~keyword_set(mplot) THEN mplot= [0,2,2]
IF ~keyword_set(patches) THEN patches= 16200L
d= 3

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


IF file_test(filename) THEN BEGIN
    IF keyword_set(QA) THEN splog, "Restoring "+filename
    restore, filename
ENDIF ELSE BEGIN
    ra_bin= sqrt(360.*180/patches)
    nras= long(360/ra_bin)
    ras= dindgen(nras)*ra_bin
    dec_bin= sqrt(360.*180/patches)
    ndecs= long(180/dec_bin)
    decs= dindgen(ndecs)*dec_bin-90.+dec_bin/2

    altentropies= dblarr(nras,ndecs)
    altentropiesoned= dblarr(nras*ndecs)

    ;;Now calculate the entropies
    FOR ii= 0L, nras-1 DO FOR jj= 0L, ndecs-1 DO BEGIN
        radec_to_lb, ras[ii],decs[jj],0,0,ll,bb,pmll,pmbb,1991.25 ;;Hipparcos epoch
        
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
        ;;Actual calculation of the entropy
        calc_entropy_sum_gaussians, amp, rmean_marg, rcovar_marg_noerr, nboot=100, $
          nsamplings=1000, entropy=altentropy, sigma_H=sigma_H, seed=seed
        IF sigma_H/altentropy GT 0.001 THEN BEGIN
            splog, 'Warning: entropy sampling has variation GT 0.001'
            splog, 'Trying with 10x more samples...'
            calc_entropy_sum_gaussians, amp, rmean_marg, rcovar_marg_noerr, nboot=100, $
              nsamplings=10000, entropy=altentropy, sigma_H=sigma_H, seed=seed
            IF sigma_H/altentropy GT 0.001 THEN $
              splog, 'Warning: entropy sampling *still* has variation GT 0.001'
        ENDIF

        altentropies[ii,jj]= altentropy
        altentropiesoned[ii*ndecs+jj]= altentropy
    ENDFOR
    IF keyword_set(filename) THEN BEGIN
        IF keyword_set(QA) THEN splog, "Saving entropies in "+filename
        save, filename=filename, ras, decs, nras, ndecs, ra_bin, $
          dec_bin, altentropies, altentropiesoned
    ENDIF
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

;;Then plot the most informative one + the observed distribution
npatch= mplot[1]*mplot[2]

;;Most informative one
sort_indx= sort(altentropiesoned)
sorted_altentropiesoned= altentropiesoned[sort_indx]
plotted_ras= dblarr(npatch)-1.
plotted_decs= dblarr(npatch)-91.
;max_info_indx= where(altentropiesoned EQ min(altentropiesoned))
;ra= ras[max_info_indx/ndecs]
;dec= decs[max_info_indx MOD ndecs]

    
IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=ysize


!p.multi=mplot
!p.charsize=.5
!p.charthick=2


FOR ii= 0L, npatch-1 DO BEGIN
    ;;Find the stars in gcs in the relevant patch
    range_in_range= bytarr(nras*ndecs)+1B
    FOR pp= 0L, npatch-1 DO BEGIN
        IF plotted_ras[pp] EQ -1. THEN CONTINUE
        FOR rr= 0L, nras-1 DO FOR dd= 0L, ndecs-1 DO BEGIN
            IF decs[dd] LT 0. THEN BEGIN
                range_in_range[rr*ndecs+dd]= 0
                CONTINUE
            ENDIF
            range_in_range[rr*ndecs+dd]*= (1B-bovy_in_range_periodic_nd([(ras[rr]-ra_span/2.+360.) MOD 360.,((decs[dd]-dec_span/2.+270.) MOD 180.)-90.],$
                                                                        [[plotted_ras[pp]-ra_span/2.,plotted_ras[pp]+ra_span/2.],$
                                                                         [plotted_decs[pp]-dec_span/2.,plotted_decs[pp]+dec_span/2.]],$
                                                                        [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([ras[rr]+ra_span/2. MOD 360.,((decs[dd]-dec_span/2.+270.) MOD 180.)-90.],$
                                            [[plotted_ras[pp]-ra_span/2.,plotted_ras[pp]+ra_span/2.],$
                                             [plotted_decs[pp]-dec_span/2.,plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([(ras[rr]-ra_span/2.+360.) MOD 360.,((decs[dd]+dec_span/2.+90.) MOD 180.)-90. ],$
                                            [[plotted_ras[pp]-ra_span/2.,plotted_ras[pp]+ra_span/2.],$
                                             [plotted_decs[pp]-dec_span/2.,plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([ras[rr]+ra_span/2. MOD 360.,((decs[dd]+dec_span/2.+90.) MOD 180.)-90.],$
                                            [[plotted_ras[pp]-ra_span/2.,plotted_ras[pp]+ra_span/2.],$
                                             [plotted_decs[pp]-dec_span/2.,plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([(ras[rr]-ra_span/2.+360.) MOD 360.,((decs[dd]-dec_span/2.+270.) MOD 180.)-90.],$
                                            [[180.+plotted_ras[pp]-ra_span/2.,180.+plotted_ras[pp]+ra_span/2.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([ras[rr]+ra_span/2. MOD 360.,((decs[dd]-dec_span/2.+270.) MOD 180.)-90.],$
                                            [[180.+plotted_ras[pp]-ra_span/2,180.+plotted_ras[pp]+ra_span/2.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([(ras[rr]-ra_span/2.+360.) MOD 360.,((decs[dd]+dec_span/2.+90.) MOD 180.)-90. ],$
                                            [[180.+plotted_ras[pp]-ra_span/2.,180.+plotted_ras[pp]+ra_span/2.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([ras[rr]+ra_span/2. MOD 360.,((decs[dd]+dec_span/2.+90.) MOD 180.)-90.],$
                                            [[180.+plotted_ras[pp]-ra_span/2.,180.+plotted_ras[pp]+ra_span/2.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([(ras[rr]-ra_span/2.+360.) MOD 360.,((decs[dd]-dec_span/2.+270.) MOD 180.)-90.],$
                                            [[plotted_ras[pp]-ra_span/2.-180.,plotted_ras[pp]-180.+ra_span/2.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([ras[rr]+ra_span/2. MOD 360.,((decs[dd]-dec_span/2.+270.) MOD 180.)-90.],$
                                            [[plotted_ras[pp]-ra_span/2.-180.,plotted_ras[pp]-180.+ra_span/2.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([(ras[rr]-ra_span/2.+360.) MOD 360.,((decs[dd]+dec_span/2.+90.) MOD 180.)-90. ],$
                                            [[plotted_ras[pp]-ra_span/2.-180.,plotted_ras[pp]+ra_span/2.-180.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))*$
              (1B-bovy_in_range_periodic_nd([ras[rr]+ra_span/2. MOD 360.,((decs[dd]+dec_span/2.+90.) MOD 180.)-90.],$
                                            [[plotted_ras[pp]-ra_span/2.-180.,plotted_ras[pp]+ra_span/2.-180.],$
                                             [-plotted_decs[pp]-dec_span/2.,-plotted_decs[pp]+dec_span/2.]],$
                                            [[0.,360.],[-90.,90.]]))

        ENDFOR
    ENDFOR
    FOR rr=0L, nras*ndecs-1 DO range_in_range[rr]= 1B-range_in_range[rr]
    ;;Sort this range
    range_in_range= range_in_range[sort_indx]
    IF keyword_set(QA) THEN splog, strtrim(string(n_elements(where(range_in_range EQ 0B))),2)+" possible ra and decs in this iteration"
    IF keyword_set(highentropy) THEN $
      possible_stars= (where(range_in_range EQ 0B))[n_elements(where(range_in_range EQ 0B))-1] ELSE $
      possible_stars= (where(range_in_range EQ 0B))[0]
    plot_this_entropy= sorted_altentropiesoned[possible_stars]
    IF keyword_set(QA) THEN splog, "Plotting entropy "+strtrim(string(plot_this_entropy,format='(F5.3)'),2)
    ra= ras[where((altentropiesoned EQ plot_this_entropy))/ndecs]
    dec= decs[where((altentropiesoned EQ plot_this_entropy)) MOD ndecs]
    plotted_ras[ii]= ra
    plotted_decs[ii]= dec

    ra_in_range= dblarr(ngcs)
    dec_in_range= dblarr(ngcs)
    FOR gg=0L, ngcs-1 DO BEGIN
        ra_in_range[gg]= bovy_in_range_periodic(gcs[gg].ra,[ra-ra_span/2.,ra+ra_span/2.],[0.,360.])
        dec_in_range[gg]= bovy_in_range_periodic(gcs[gg].dec,[dec-dec_span/2.,dec+dec_span/2.],[-90.,90.])
    ENDFOR
    indx_patch= where(ra_in_range EQ 1 AND dec_in_range EQ 1)
    ra_in_range2= dblarr(ngcs)
    dec_in_range2= dblarr(ngcs)
    FOR gg=0L, ngcs-1 DO BEGIN
        ra_in_range2[gg]= bovy_in_range_periodic(180.-gcs[gg].ra,[ra-ra_span/2.,ra+ra_span/2.],[0.,360.])+$
          bovy_in_range_periodic(180.+gcs[gg].ra,[ra-ra_span/2.,ra+ra_span/2.],[0.,360.])
        dec_in_range2[gg]= bovy_in_range_periodic(-gcs[gg].dec,[dec-dec_span/2.,dec+dec_span/2.],[-90.,90.])
    ENDFOR
    indx_patch2= where(ra_in_range2 EQ 1 AND dec_in_range2 EQ 1)
    n_in_patch= n_elements(indx_patch)+n_elements(indx_patch2)

    ;;Find the distribution
    radec_to_lb, ra,dec,$
      0,0,ll,bb,pmll,pmbb,1991.25 ;;Hipparcos epoch
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
    npix= 25                    ;ceil(0.3*sqrt(n_in_patch)) > 10
    binwidth= (rrange[1]-rrange[0])/double(npix+1) ;;Conforms to the way hogg_plothist works
    norm= total(ds_marg,/NaN)*step
    weight = make_array(n_in_patch,/float,Value=norm/float(n_in_patch))
    
;;Plot
    xtitle= 'v_r [km s^{-1}]'
    ytitle= 'P(v_r!10#!X\alpha,\delta)'
    hist_this= dblarr(n_in_patch)
    hist_this[0:n_elements(indx_patch)-1]= gcs[indx_patch].vr
    hist_this[n_elements(indx_patch):n_in_patch-1]= -gcs[indx_patch2].vr
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
    lstring= textoidl('\alpha = ')+strtrim(string(ra,$
                                                  format='(F5.1)'))+$
      textoidl('^\circ')
    lstring2= textoidl('\delta = ')+strtrim(string(dec,$
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
