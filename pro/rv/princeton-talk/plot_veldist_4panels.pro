;+
;   NAME:
;      plot_veldist_4panels
;   PURPOSE:
;      plot 3 projections of the velocity distribution + 4th which has
;      the actual Gaussians
;   CALLING SEQUENCE:
;   INPUT:
;      k   - number of Gaussians
;      w   - regularization parameter
;      xrange
;      yrange
;      zrange
;   KEYWORDS:
;      plotdata - overplot a sampling of the data as gray lines (it's
;                 a mess)
;   OUTPUT:
;   REVISION HISTORY:
;      2008-12-10 - Written Bovy (NYU)
;      2008-12-16 - Thickness of ellipses propto their amplitudes
;-
PRO PLOT_VELDIST_4PANELS, K, w, sample=sample, xrange=xrange, $
                          yrange=yrange, zrange=zrange, filename=filename, $
                          plotdata=plotdata

IF ~keyword_set(sample) THEN sample= 6
IF ~keyword_set(xrange) THEN xrange=[-130,120]
IF ~keyword_set(yrange) THEN yrange=[-120,60]
IF ~keyword_set(zrange) THEN zrange=[-70,70]
IF ~keyword_set(filename) THEN filename='4veldist'

savefilename= '../ALMSv2_reconvergev2/'+'fitv'+strtrim(string(K),2)+$
  '_bestV_V'+strtrim(string(w,format='(f4.1)'),2)+$
  '_sample'+strtrim(string(sample),2)+$
  '.sav'

restore, savefilename

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

k_print, filename=filename+'1.ps'
!p.multi=[2,2,1]

charsize=.6
;;First plot X-Y
plot_projected_gaussians, mean, covar, projection[*,*,2], xrange=xrange, $
  yrange=yrange, amp=amp, xlabel=xlab, ylabel=ylab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=400, charsize=charsize,/keepbangX

;;Overplot the measurements
;;First get the measurements
restorefilename='/global/data/hipparcos/hipnewredux/hip2-dehnen.sav'
restore, filename=restorefilename
ndata= n_elements(hip.hip)
UVW= dblarr(ndata,3)
radialdir=dblarr(ndata,3)
boundaries= dblarr(ndata,3,2) ;;For each data point, calculate where the projected velocity reaches the boundary for each coordinate
FOR ii=0L, 0 DO BEGIN
    l= hip[ii].l
    b= hip[ii].b
    ;;Create the rotation matrices
    rot1= dblarr(3,3)
    rot2= dblarr(3,3)
    rot1[0:2,0]= [cos(l),-sin(l),0]
    rot1[0:2,1]= [sin(l),cos(l),0]
    rot1[0:2,2]= [0,0,1]
    rot2[0:2,0]= [cos(b),0,-sin(b)]
    rot2[0:2,1]= [0,1,0]
    rot2[0:2,2]= [sin(b),0,cos(b)]
    ;;Rotate the velocities
    vrlb= [[0],[hip[ii].vl],[hip[ii].vb]]
    UVW[ii,0:2]= (rot1##rot2)##vrlb
    radialdir[ii,0:2]= [cos(l)*cos(b),sin(l)*cos(b),sin(b)]
    boundaries[ii,0,0]= (xrange[0]-UVW[ii,0])/radialdir[ii,0]
    boundaries[ii,0,1]= (xrange[1]-UVW[ii,0])/radialdir[ii,0]
    boundaries[ii,1,0]= (yrange[0]-UVW[ii,1])/radialdir[ii,1]
    boundaries[ii,1,1]= (yrange[1]-UVW[ii,1])/radialdir[ii,1]
    boundaries[ii,2,0]= (zrange[0]-UVW[ii,2])/radialdir[ii,2]
    boundaries[ii,2,1]= (zrange[1]-UVW[ii,2])/radialdir[ii,2]
    djs_oplot,xrange,[radialdir[ii,1]/radialdir[ii,0]*(xrange[0]-UVW[ii,0])+UVW[ii,1],$
                      radialdir[ii,1]/radialdir[ii,0]*(xrange[1]-UVW[ii,0])+UVW[ii,1]], $
      psym=-3, thick=3.
ENDFOR

;;Plot X-Y with ellipses
djs_plot, [0.],[0.], xrange=xrange, yrange=yrange, /isotropic, $
  xtitle=xlab, charsize=charsize, ytickname=REPLICATE(' ',30)
minthick=0.01
maxthick=4.
minalogamp= max(alog(amp))-2.
thicks= (alog(amp)-minalogamp)/(max(alog(amp)-minalogamp))*(maxthick-minthick)+minthick > minthick
FOR kk=0L, K-1 DO hogg_oplot_covar, mean[0,kk], mean[1,kk], covar[0:1,0:1,kk], thick=thicks[kk]
;FOR kk=0L, K-1 DO hogg_oplot_covar, mean[0,kk], mean[1,kk], covar[0:1,0:1,kk], thick=alog(amp[kk])-min(alog(amp))-4

k_end_print

k_print, filename=filename+'2.ps'
!p.multi=[2,2,1]

;;Then plot X-Z
plot_projected_gaussians, mean, covar, projection[*,*,1], xrange=xrange, $
  yrange=zrange, amp=amp, xlabel=xlab, ylabel=zlab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=grid, charsize=charsize

;;Finally plot Y-Z
plot_projected_gaussians, mean, covar, projection[*,*,0], xrange=xrange, $
  yrange=zrange, amp=amp, xlabel=ylab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=grid, charsize=charsize, $
  ytickname=REPLICATE(' ',30)


k_end_print

;;Done!
END
