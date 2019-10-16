;+
;   NAME:
;      plot_kwcrossval
;   PURPOSE:
;      plot the cross validation surface in (k,w) FOR THE ALMS sample,
;      new reduction!
;   CALLING SEQUENCE:
;   INPUT:
;      plotfilename   - save the plot here
;      horizontal     - if set, draw a horizontal colorbar
;      ks_low         - smallest K to include
;      ks_high        - highest K to include
;      ks             - explicitly set all the ks
;      ws             - explicitly set all the ws
;      basedir        - dir that holds the fitv savefiles
;   OUTPUT:
;   REVISION HISTORY:
;      2009-01-12 - Written Bovy
;-
PRO PLOT_KWCROSSVAL, plotfilename=plotfilename, $
                    horizontal=horizontal, basedir=basedir, $
                     ks_low=ks_low, ks_high=ks_high, $
                     ks=ks, ws=ws

sample= 6;;ALMS
charsize=.6

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

;;Read the data
cmd= '../make_crossval_file.sh '+strtrim(string(ks[0]),2)+' '+strtrim(string(ks[nks-1]),2)+' '+basedir
spawn, cmd
readcol, 'crossval.dat', crossval, format='D'
cmd= '\rm crossval.dat'
spawn, cmd

;;Put the data in the right format
loglike= dblarr(nks,nws)
FOR kk=0L, nks-1 DO FOR ww= 0L, nws-1 DO BEGIN
    loglike[kk,ww]= crossval[ww+kk*nws]
ENDFOR


;;Print the max
max_indx= where(crossval EQ max(crossval))
splog, "Max= "+strtrim(string(crossval(max_indx[0])),2)+" @: K= "+strtrim(string(ks[max_indx[0]/nws]),2)+", w= "+strtrim(string(ws[max_indx[0] MOD nws],format='(F4.2)'),2)

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
tvimage= (256-ncolors)+(ncolors-1-((floor(ncolors*(tvimage-tvrange[0])/ $
                         (tvrange[1]-tvrange[0])) > 0) $
                  < (ncolors-1)))
tv, tvimage, !X.RANGE[0]+kspacing+2*kspacing*(ks[startk]-1),$
  !Y.RANGE[0]+wspacing,/data, $
  xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
  ysize=(!Y.CRANGE[1]-!Y.CRANGE[0])

;;Remake axes, and put axis-labels
!P.MULTI[0]= !P.MULTI[0]+1
djs_plot, [ks[max_indx[0]/nws]],[sqrt(ws[max_indx[0] MOD nws])],$
  xrange=xrange,yrange=yrange, charthick=charthick,$
  xtitle='K', ytitle=TEXTOIDL('!9r!xw [km s^{-1}]'), charsize=charsize, $
  position=tvposition, ytickinterval=ytickinterval, $
  psym=7, color=djs_icolor('white')

;;Plot colorbar
nsamples= 256
ticklen=1
IF keyword_set(horizontal) THEN BEGIN
    cb= dblarr(nsamples,2)
    cb[*,0]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cb[*,1]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cbposition= [!X.WINDOW[0],!Y.WINDOW[1]+.02,!X.WINDOW[1],!Y.WINDOW[1]+.05]
    djs_plot, [0.], [0.], yrange=[0,1], xrange=tvrange, $
      xtickname=REPLICATE(' ',30),charthick=charthick,$
      position=cbposition, ytickname=REPLICATE(' ',30), xthick=2,$
      yticklen=.0001, /NOERASE, xticklen=ticklen, xstyle=5, charsize=charsize
ENDIF ELSE BEGIN
    cb= dblarr(2,nsamples)
    cb[0,*]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cb[1,*]= dindgen(nsamples)*(tvrange[1]-tvrange[0])/(nsamples-1)+tvrange[0]
    cbposition= [!X.WINDOW[1]+.02,!Y.WINDOW[0],!X.WINDOW[1]+.05,!Y.WINDOW[1]]
    djs_plot, [0.], [0.], xrange=[0,1], yrange=tvrange, $
      ytickname=REPLICATE(' ',30),charthick=charthick,$
      position=cbposition, xtickname=REPLICATE(' ',30), $
      xticklen=.0001, /NOERASE, yticklen=.02, ystyle=5, charsize=charsize
ENDELSE
cbimage= (256-ncolors)+ncolors-1-((floor(ncolors*(cb-tvrange[0])/ $
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
    djs_plot, [0],[1],yrange=[0,1], xrange=tvrange, $
      ytickname=REPLICATE(' ',30),  charsize=charsize, charthick=charthick,$
      position=cbposition, yticklen=.0001, /NOERASE, xstyle=5
    axis, xaxis=0, xtickformat='(A1)', charsize=0.35*charsize, $
      xticklen=ticklen, charthick=0.6*charthick
    axis, xaxis=1, $
      charsize=.35*charsize, xtickformat='(E10.3)', charthick=.6*charthick
    djs_plot, [0.],[0.], /NOERASE, title='ln(cross-validation likelihood)', $
      ystyle=4, xstyle=4, charsize=charsize, position=[0.1,0.1,.9,.9], $
      charthick=charthick
ENDIF ELSE BEGIN
    djs_plot, [0],[1],xrange=[0,1], yrange=tvrange, charthick=charthick,$
      xtickname=REPLICATE(' ',30),  charsize=charsize, $
      position=cbposition, xticklen=.0001, /NOERASE, ystyle=5
    axis, yaxis=0, ytickformat='(A1)', charsize=0.4*charsize, $
      yticklen=ticklen, charthick=.6*charthick
    axis, yaxis=1, $
      charsize=.4*charsize, ytickformat='(E10.3)', charthick=0.6*charthick
    djs_plot, [0.],[0.], /NOERASE, title='ln(cross-validation likelihood)', $
      ystyle=4, xstyle=4, charsize=charsize, position=[0.1,0.1,.9,.9],$
      charthick=charthick
ENDELSE

IF keyword_set(plotfilename) THEN k_end_print


END
