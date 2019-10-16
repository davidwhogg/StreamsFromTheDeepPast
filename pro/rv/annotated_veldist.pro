;+
;   NAME:
;      annotated_veldist
;   PURPOSE:
;      plot the three projections of the final velocity distribution
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      filename - filename for plot
;      basedir - basedirectory for density estimate (ends in /)
;      xrange
;      yrange
;      zrange
;   KEYWORDS:
;      QA    - print QA info
;   OUTPUT:
;   REVISION HISTORY:
;      2009-04-08 - Written - Bovy (NYU)
;-
PRO ANNOTATED_VELDIST, K=K, w=w, sample=sample, filename=filename, $
                       basedir=basedir, QA=QA, xrange=xrange, yrange=yrange, $
                       zrange=zrange

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(xrange) THEN xrange= [-130,120]
IF ~keyword_set(yrange) THEN yrange= [-120,60]
IF ~keyword_set(zrange) THEN zrange= [-70,70]
lsr= [-10.1,-4.0,-6.7]

;;Restore the velocity distribution
;;Restore density estimate
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


;;Prepare for plots
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
;Set up labels
xlab= TEXTOIDL('v_x [km s^{-1}]')
ylab= TEXTOIDL('v_y [km s^{-1}]')
zlab= TEXTOIDL('v_z [km s^{-1}]')
;Define the projection matrices
;projection[*,*,[x,y,z]] = projection perpendicular to x-direction. etc.
d=3
projection= dblarr(d,d,d)
FOR ii=0L, d-1 DO FOR jj=0L, d-1 DO BEGIN
    IF ii NE jj THEN projection[jj,jj,ii]= double(3-jj)
ENDFOR


k_print, filename=filename
charsize=.6
;;lsr symbol
xarr = [-1,0,1,-1]
yarr = [-1,1,-1,-1]


;;First plot X-Y
position=[0.1,.324,0.9,0.9]
plot_projected_gaussians, mean, covar, projection[*,*,2], xrange=xrange, $
  yrange=yrange, amp=amp, xlabel=xlab, ylabel=ylab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=400, charsize=charsize, $
  position=position, xstyle=5, /keepbangX
;axis, 0.1,.324,xaxis=0, xtickformat='(A1)'
axis,xaxis=1, xtitle=xlab, /NOERASE, charsize=charsize
;;Add the lsr label
usersym, xarr, yarr,/fill, color=djs_icolor('white')
djs_oplot, [lsr[0]],[lsr[1]], psym=8
usersym, xarr, yarr
djs_oplot, [lsr[0]],[lsr[1]], psym=8

;;Then plot X-Z
position= [0.1,.1,0.5,.324]
plot_projected_gaussians, mean, covar, projection[*,*,1], xrange=xrange, $
  yrange=zrange, amp=amp, xlabel=xlab, ylabel=zlab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=grid, charsize=charsize, $
  position=position, /NOERASE,/keepbangX
;;Add the lsr label
usersym, xarr, yarr,/fill, color=djs_icolor('white')
djs_oplot, [lsr[0]],[lsr[2]], psym=8
usersym, xarr, yarr
djs_oplot, [lsr[0]],[lsr[2]], psym=8

;;Finally plot Y-Z
position= [0.5,.1,0.9,.324]
plot_projected_gaussians, mean, covar, projection[*,*,0], xrange=xrange, $
  yrange=zrange, amp=amp, xlabel=ylab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=grid, charsize=charsize, $
  ytickname=REPLICATE(' ',30), $
  position=position, /NOERASE,/keepbangX
;;Add the lsr label
usersym, xarr, yarr,/fill, color=djs_icolor('white')
djs_oplot, [lsr[1]],[lsr[2]], psym=8
usersym, xarr, yarr
djs_oplot, [lsr[1]],[lsr[2]], psym=8


k_end_print


    


END
