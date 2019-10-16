;+
;   NAME:
;      plot_kwMML
;   PURPOSE:
;      plot the MML surface in (k,w) FOR THE ALMS sample,
;      new reduction!
;   CALLING SEQUENCE:
;   INPUT:
;      filename       - save filename
;      plotfilename   - save the plot here
;      ks_low         - smallest K to include
;      ks_high        - highest K to include
;      ks             - explicitly set all the ks
;      ws             - explicitly set all the ws
;      basedir        - dir that holds the fitv savefiles
;   KEYWORDS:
;      horizontal     - if set, draw a horizontal colorbar
;      upper          - if set, use the upper bound on the lattice const
;   OUTPUT:
;   REVISION HISTORY:
;      2009-01-14 - Written Bovy (NYU)
;-
PRO PLOT_KWMML, filename=filename, plotfilename=plotfilename, $
                    horizontal=horizontal, upper=upper, basedir=basedir, $
                     ks_low=ks_low, ks_high=ks_high, $
                     ks=ks, ws=ws

sample= 6;;ALMS
charsize=.6
d= 3

IF ~keyword_set(basedir) THEN basedir= '../ALMSv2_reconverge/'

IF keyword_set(ks_low) THEN BEGIN
    IF keyword_set(ks_high) THEN BEGIN
        nks= ks_high-ks_low+1
        ks= lindgen(nks)+ks_low
    ENDIF ELSE BEGIN
        ks_high= 15
        nks= ks_high-ks_low+1
        ks= lindgen(nks)+ks_low
    ENDELSE
ENDIF ELSE IF keyword_set(ks_high) THEN BEGIN
    ks_low= 1
    nks= ks_high-ks_low+1
    ks= lindgen(nks)+ks_low
ENDIF ELSE BEGIN
    nks= 14
    ks= lindgen(nks)+1L
ENDELSE

IF ~keyword_set(ws) THEN ws= [1.,4.,9.]
nws= n_elements(ws)

mltwo= dblarr(nks,nws)
mlone= dblarr(nks*nws)

;;Calculate the ML
IF keyword_set(filename) AND file_test(filename) THEN BEGIN
    splog, "Restoring savefile: "+filename
    restore, filename
ENDIF ELSE BEGIN
    ;;Restore data
    dehnenfilename= '/global/data/hipparcos/hipnewredux/hip2-dehnen.sav'
    splog, 'QA: Using the new Hipparcos reduction'
    IF file_test(dehnenfilename) THEN BEGIN
        restore, dehnenfilename
    ENDIF ELSE BEGIN
        splog, 'Data file not found in '+dehnenfilename
        splog, 'Data file must exist!'
        splog, 'Returning...'
        RETURN
    ENDELSE
    nhip= n_elements(hip)
    splog, 'QA: Data sample has '+strtrim(string(nhip),2)+' members'
    minamp= 1D0/nhip

    ;;Run EM with k=1 to find the sigma_pops
    ydata= transpose([[dblarr(nhip)],[hip.vl], [hip.vb]])
    ycovar_tmp= hip.vlvbc
    projection= hip.sm
    ;;Add radial velocities
    ycovar= dblarr(3,3,nhip)
    FOR gg= 0, nhip-1 DO BEGIN
        ;ycovar[*,*,gg]=[[10000.^2,0.,0.],[0.,ycovar_tmp[0,0,gg],ycovar_tmp[0,1,gg]],[0.,ycovar_tmp[1,0,gg],ycovar_tmp[1,1,gg]]]
        ycovar[0,0,gg]=10000.^2
    ENDFOR
    ;;set-up initial conditions
    xamp= 1.
    xmean= dblarr(d)
    xcovar= dblarr(d,d)
    FOR dd= 0L, d-1 DO xcovar[dd,dd]= 100.^2
    ;;Setting up projected_gauss_mixtures
    logfile= ''
    splitnmerge= 0
    fixmean= bytarr(1)
    fixcovar= bytarr(1)
    fixamp= bytarr(1)+1
    tol=1D-6
    constraint= 0.
    ;;Run projected_gauss_mixtures
    projected_gauss_mixtures_c, 1, ydata, ycovar, projection, $
      xamp, xmean, xcovar, fixmean=fixmean, fixcovar=fixcovar, $
      tol=tol, avgloglikedata=avgloglikedata, w=constraint, logfile=logfile, $
      splitnmerge=splitnmerge,/quiet
    ;;Now extract the sigma_pops (squared): diagonalize xcovar 
    sigmapops= EIGENQL(xcovar,/double)
    ;;Now load the right covariance matrices
    FOR gg= 0, nhip-1 DO BEGIN
        ycovar[*,*,gg]=[[10000.^2,0.,0.],[0.,ycovar_tmp[0,0,gg],ycovar_tmp[0,1,gg]],[0.,ycovar_tmp[1,0,gg],ycovar_tmp[1,1,gg]]]
    ENDFOR
    
    ;;Now run through the different models to find the ML for each
    FOR kk=0L, nks-1 DO FOR ww=0L, nws-1 DO Begin
        ;;Restore the right file
        savefilename= basedir+'fitv'+strtrim(string(ks[kk]),2)+$
          '_bestV_V'+strtrim(string(ws[ww],format='(f4.1)'),2)+$
          '_sample'+strtrim(string(sample),2)+$
          '.sav'
        restore, savefilename
        ;Compute total log likelihood (this is needed in the ML)
        avgloglikedata= 0.
        projected_gauss_mixtures_C, ks[kk], ydata, ycovar, $
          projection, amp, mean, covar, /likeonly, $
          avgloglikedata=avgloglikedata, /quiet
        ;;Calculate ML
        ML= 0D0
        FOR dd= 0L, d-1 DO ML+= alog(2D0*sigmapops[dd])
        ML*= ks[kk]
        ML-= alog(GAMMA(ks[kk]))
        Np= double(ks[kk])*(1D0+d+d*(d+1)/2)-1;;Number of parameters
        lattice= lattice_const(Np,/zador,upper=upper);;Uses the Zador lower bound for now
        ML+= Np/2.*alog(lattice)
        ;;Make sure the amps are finite
        amp= amp > minamp
        FOR kk2= 0L, ks[kk]-1 DO BEGIN
            ;;Diagonalize the covariance of this Gaussian
            sigmas= EIGENQL(covar[*,*,kk2],/double);;Actually sigma-squared
            FOR dd= 0L, d-1 DO ML+= alog(sqrt(2.)*amp[kk2]*nhip/sigmas[dd])
        ENDFOR
        ML+= 0.5D0*alog(nhip)
        FOR kk2=0L, ks[kk]-1 DO ML-= 0.5D0*alog(amp[kk2])
        ML-= avgloglikedata*nhip
        ML+= Np/double(2)
        mltwo[kk,ww]= ML
        mlone[ww+kk*nws]= ML
        ;splog, ML, mltwo[kk,ww], mlone[ww+kk*nws]
    ENDFOR
    splog, mlone
    IF keyword_set(filename) THEN BEGIN
        save, mltwo, ks,ws,mlone, filename=filename
        splog, "Saved MML surface in "+filename
    ENDIF
ENDELSE


loglike= mltwo
ml= mlone
;;Print the min
min_indx= where(ml EQ min(ml))
splog, "Min= "+strtrim(string(ml(min_indx[0])),2)+$
  " @: K= "+strtrim(string(ks[min_indx[0]/nws]),2)+$
  ", w= "+strtrim(string(ws[min_indx[0] MOD nws],format='(F4.2)'),2)

;;Plot
charthick= 2
ytickinterval=1.
xsize=2.3;one third of ApJ width
IF keyword_set(horizontal) THEN $
  tvposition= [.1,.1,.9,.8] $
  ELSE tvposition= [.1,.1,.8,.9]

;;Make axes
;;We assume that the ks are equally spaced
startk= 0
kspacing= (ks[startk+1]-ks[startk])/2.
xrange= [ks[startk]-kspacing,ks[n_elements(ks)-1]+kspacing]
;;we assume that the ws are equally spaced in sqrt(w)
wspacing= (sqrt(ws[1])-sqrt(ws[0]))/2.
yrange= [sqrt(ws[0])-wspacing,sqrt(ws[n_elements(ws)-1])+wspacing]

IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=xsize

djs_plot, [0], [1], xrange=xrange,yrange=yrange, charsize=charsize, $
  position=tvposition, charthick=charthick, ytickinterval=ytickinterval


tvimage= loglike[startk:n_elements(ks)-1,*]
;tvimage= dblarr(n_elements(ks),n_elements(ws))
;tvimage[0,0]=10
tvrange= minmax(tvimage)
ncolors= 190
tvimage= (256-ncolors)+((floor(ncolors*(tvimage-tvrange[0])/ $
                         (tvrange[1]-tvrange[0])) > 0) $
                  < (ncolors-1))
tv, tvimage, !X.RANGE[0]+kspacing+2*kspacing*(ks[startk]-1),$
  !Y.RANGE[0]+wspacing,/data, $
  xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
  ysize=(!Y.CRANGE[1]-!Y.CRANGE[0])

;;Remake axes, and put axis-labels
!P.MULTI[0]= !P.MULTI[0]+1
djs_plot, [ks[min_indx[0]/nws]],[sqrt(ws[min_indx[0] MOD nws])],$
  xrange=xrange,yrange=yrange, charthick=charthick,$
  xtitle='K', ytitle=TEXTOIDL('!9r!xw [km s^{-1}]'), charsize=charsize, $
  position=tvposition, ytickinterval=ytickinterval, $
  psym=7, color=djs_icolor('white')

;;Plot colorbar
nsamples= 256
ticklen=1
nticks= 3
IF keyword_set(horizontal) THEN BEGIN
    cb= dblarr(nsamples,2)
    cb[*,0]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cb[*,1]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cbposition= [!X.WINDOW[0],!Y.WINDOW[1]+.02,!X.WINDOW[1],!Y.WINDOW[1]+.05]
    djs_plot, [0.], [0.], yrange=[0,1], xrange=tvrange, xticks=nticks, $
      xtickname=REPLICATE(' ',30),charthick=charthick,$
      position=cbposition, ytickname=REPLICATE(' ',30), xthick=2,$
      yticklen=.0001, /NOERASE, xticklen=ticklen, xstyle=5, charsize=charsize
ENDIF ELSE BEGIN
    cb= dblarr(2,nsamples)
    cb[0,*]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cb[1,*]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cbposition= [!X.WINDOW[1]+.02,!Y.WINDOW[0],!X.WINDOW[1]+.05,!Y.WINDOW[1]]
    djs_plot, [0.], [0.], xrange=[0,1], yrange=tvrange, yticks=nticks, $
      ytickname=REPLICATE(' ',30),charthick=charthick,$
      position=cbposition, xtickname=REPLICATE(' ',30), $
      xticklen=.0001, /NOERASE, yticklen=.02, ystyle=5, charsize=charsize
ENDELSE
cbimage= (256-ncolors)+((floor(ncolors*(cb-tvrange[0])/ $
                         (tvrange[1]-tvrange[0])) > 0) $
                  < (ncolors-1))
tv, cbimage, cbposition[0],cbposition[1],/normal, $
  xsize=(cbposition[2]-cbposition[0]), $
  ysize=(cbposition[3]-cbposition[1])
;!X.RANGE[0],!Y.RANGE[0],/data, $
;      xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
;      ysize=(!Y.CRANGE[1]-!Y.CRANGE[0])
!P.MULTI[0]= !P.MULTI[0]+1
IF keyword_set(horizontal) THEN BEGIN
    djs_plot, [0],[1],yrange=[0,1], xrange=tvrange, charthick=charthick,$
      ytickname=REPLICATE(' ',30),  charsize=charsize, $
      position=cbposition, yticklen=.0001, /NOERASE, xstyle=5
    axis, xaxis=0, xtickformat='(A1)', charsize=.45*charsize, xticklen=ticklen, $
      xticks=nticks, charthick=.6*charthick
    axis, xaxis=1, charsize=.45*charsize, xtickformat='(E10.4)', $
      xticks=nticks, charthick=.6*charthick
    djs_plot, [0.],[0.], /NOERASE, title='MML', $
      ystyle=4, xstyle=4, charsize=charsize, position=[0.1,0.1,.9,.9], $
      charthick=charthick
ENDIF ELSE BEGIN
    djs_plot, [0],[1],xrange=[0,1], yrange=tvrange, charthick=charthick,$
      xtickname=REPLICATE(' ',30),  charsize=charsize, $
      position=cbposition, xticklen=.0001, /NOERASE, ystyle=5
    axis, yaxis=0, ytickformat='(A1)', charsize=.45*charsize, yticklen=ticklen, $
      yticks= nticks, charthick=.6*charthick
    axis, yaxis=1, charsize=.45*charsize, ytickformat='(E10.4)', $
      yticks= nticks, charthick=.6*charthick
    djs_plot, [0.],[0.], /NOERASE, title='MML', $
      ystyle=4, xstyle=4, charsize=charsize, position=[0.1,0.1,.9,.9],$
      charthick=charthick
ENDELSE

IF keyword_set(plotfilename) THEN k_end_print


END
