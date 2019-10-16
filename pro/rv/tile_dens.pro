;+
;   NAME:
;      tile_dens.pro
;   PURPOSE:
;      plot the velocity distribution as a function of k, w
;   CALLING SEQUENCE:
;   INPUT:
;      sample     - sample to use (default ALMS=6, see fitv for
;                   different samples)
;      xrange     - xrange for plot
;      yrange     - yrange for plot
;      zrange     - zrange for plot
;      xvsz       - plot v_x vs v_z (default v_x vs v_y)
;      yvsz       - plot v_y vs v_z
;      grid       - grid size
;      quiet      - be quiet!
;      basedir    - directory that holds the savefiles
;      ks         - ks to plot
;      ws         - ws to plot
;   OUTPUT:
;   REVISION HISTORY:
;      2009-01-12 - Written Bovy (NYU)
;-
PRO TILE_DENS, sample=sample, xrange=xrange, yrange=yrange,$
               zrange=zrange, filename=filename, xvsz=xvsz,yvsz=yvsz, $
               grid=grid, quiet=quiet, basedir=basedir, ks=ks, ws=ws

d=3

IF ~keyword_set(sample) THEN sample= 6
IF ~keyword_set(xrange) THEN xrange=[-130,120]
IF ~keyword_set(yrange) THEN yrange=[-120,60]
IF ~keyword_set(zrange) THEN zrange=[-70,70]
IF ~keyword_set(ks) THEN ks= [3,5,7,10,13,15];ks= [3,5,10,15,20,25]
nks= n_elements(ks)
IF ~keyword_set(ws) THEN ws= [1.,4.,9.]
nws= n_elements(ws)

;Projection matrices
;projection[*,*,[x,y,z]] = projection perpendicular to x-direction. etc.
projection= dblarr(d,d,d)
FOR ii=0L, d-1 DO FOR jj=0L, d-1 DO BEGIN
    IF ii NE jj THEN projection[jj,jj,ii]= double(3-jj)
ENDFOR
IF keyword_set(xvsz) THEN BEGIN
    projection= projection[*,*,1]
    xlab= TEXTOIDL('v_x [km s^{-1}]')
    ylab= TEXTOIDL('v_z [km s^{-1}]')
    xrange=xrange
    yrange=zrange
    suffix='XZ.ps'
ENDIF ELSE IF keyword_set(yvsz) THEN BEGIN
    projection= projection[*,*,0]
    xlab= TEXTOIDL('v_y [km s^{-1}]')
    ylab= TEXTOIDL('v_z [km s^{-1}]')
    xrange=yrange
    yrange=zrange
    suffix='YZ.ps'
ENDIF ELSE BEGIN
    projection= projection[*,*,2]
    xlab= TEXTOIDL('v_x [km s^{-1}]')
    ylab= TEXTOIDL('v_y [km s^{-1}]')
    xrange=xrange
    yrange=yrange
    suffix='XY.ps'
ENDELSE


IF ~keyword_set(filename) THEN filename='veldens'+suffix


;;For legend later
xspan= max(xrange)-min(xrange)
yspan= max(yrange)-min(yrange)


;Set up contours
cntrlevels= dblarr(10)
ncntr= n_elements(cntrlevels)
FOR ii= 0L, ncntr-1 DO cntrlevels[ncntr-1-ii]=10D^(-ii)
nsub=4
ndiv=6
denslevels= dblarr(ndiv*nsub)
ndens= n_elements(denslevels)
FOR ii= 0L, ndiv-1 DO FOR jj=0L, nsub-1 DO denslevels[nsub*ii+jj]= 10D^(-ndiv+ii+1)*double(jj+1)/double(nsub)
denslevels= denslevels[UNIQ(denslevels, SORT(denslevels))]


cntrlevels=[.02D,.06D,.12D,.21D,.33D,.5D,.68D,.8D,.9D,.95D,.99D,.999D]




IF keyword_set(filename) THEN k_print, filename=filename

xmargin= .1
ymargin= .1
ywidth= .13
xwidth= ywidth*(xrange[1]-xrange[0])/(yrange[1]-yrange[0])

charsize= .4





FOR kk= 0L, nks-1 DO FOR ww= 0L, nws-1 DO BEGIN
    ;;Restore the right file
    savefilename= basedir+'fitv'+strtrim(string(ks[kk]),2)+$
      '_bestV_V'+strtrim(string(ws[ww],format='(f4.1)'),2)+$
      '_sample'+strtrim(string(sample),2)+$
      '.sav'
    restore, savefilename
    IF ww EQ 0 THEN BEGIN
        ytickname=REPLICATE('',30)
        ylabel= ylab
    ENDIF ELSE BEGIN
        ytickname=REPLICATE(' ',30)
        ylabel= " "
    ENDELSE
    IF kk EQ nks-1 THEN BEGIN
        xtickname=REPLICATE('',30) 
        xlabel=xlab
    ENDIF ELSE BEGIN
        xtickname=REPLICATE(' ',30)
        xlabel= " "
    ENDELSE
    plot_projected_gaussians, mean, covar, projection, xrange=xrange, $
      yrange=yrange, amp=amp, xlabel=xlabel, ylabel=ylabel, filename=fileXY, $
      cntrlevels=cntrlevels, denslevels=denslevels, quiet=quiet,grid=grid, $
      position=[xmargin+ww*xwidth,ymargin+(nks-kk-1)*ywidth,$
      xmargin+(ww+1)*xwidth,ymargin+(nks-kk)*ywidth], /NOERASE, $
      xtickname=xtickname, ytickname=ytickname, charsize=charsize
    
    IF kk EQ 0 THEN BEGIN
        legend, ["w="+strtrim(string(ws[ww],format='(f4.1)'),2)+textoidl(' km s^{-1}')],$
          pos=[xmargin+ww*xwidth+xwidth/2-.07,ymargin+(nks-kk)*ywidth+.03], $
          box=0, charsize=1.5*charsize,/normal
    ENDIF
    IF ww EQ nws-1 THEN BEGIN
        legend, ["K="+strtrim(string(ks[kk]),2)], $
          pos=[xmargin+(ww+1)*xwidth-.01,ymargin+(nks-kk)*ywidth-ywidth/2.+.02], $
          box=0, charsize=1.5*charsize,/normal
    ENDIF
ENDFOR

splog, "Success!"


IF keyword_set(filename) THEN k_end_print

END
