;+
;   NAME:
;      plot_projected_gaussians_wrap
;
;   PURPOSE:
;      plots 2D projections of 3D gaussians
;
;   CALLING SEQUENCE:
;      plot_projected_gaussians_wrap, mean, covar, xrange=xrange,
;      yrange=yrange, zrange=zrange [, amp=amp, quiet=quiet
;      basefilename=basefilename, label=label]
;
;   INPUTS:
;      mean       : means of the gaussians (can be a vector, i.e. 
;                  dblarr(3,ngauss)
;      covar      : covariances of the gaussians (as a vector:
;                  dblarr(3,3,ngauss)
;      [xyz]range : x-range, etc.
;
;   OPTIONAL INPUTS:
;      amp        : relative amplitudes of the gaussians (should add up
;                   to 1; default: 1D/ngauss for all)
;      basefilename : Sets the base-filename. XY/YZ/XZ is added to
;                     this to make the different plots
;      label       : label to be put in the top-left corner
;      grid        : number of pixels in each direction to use
;
;   KEYWORDS:
;      quiet           : be (relatively) quiet.
;      log             : plot projected log density, not linear (which
;                        is the default)
;
;   OUTPUTS:
;      ps-files with the plots
;
;   REVISION HISTORY:
;      06-03-08  - Written Bovy
;-
;
PRO PLOT_PROJECTED_GAUSSIANS_WRAP, mean, covar, xrange=xrange, $
                                   yrange=yrange, zrange=zrange,amp=amp, $
                                   basefilename=basefilename, label=label, $
                                   quiet=quiet, grid=grid, log=log

;Set dimension 3D
d=3

;Check inputs
ngauss= n_elements(mean[0,*])
IF n_elements(covar[0,0,*]) NE ngauss THEN BEGIN
    splog, 'Error: Number of covariances has to equal the number of means'
    RETURN
ENDIF
IF NOT keyword_set(amp) THEN amp= dblarr(ngauss)+1D/double(ngauss)
IF n_elements(amp) NE ngauss THEN BEGIN
    splog, 'Error: number of given amplitudes should equal the number of means'
    RETURN
ENDIF


;Define the projection matrices
;projection[*,*,[x,y,z]] = projection perpendicular to x-direction. etc.
projection= dblarr(d,d,d)
FOR ii=0L, d-1 DO FOR jj=0L, d-1 DO BEGIN
    IF ii NE jj THEN projection[jj,jj,ii]= double(3-jj)
ENDFOR

;Set up filenames
fileXY= basefilename+'XY.ps'
fileXZ= basefilename+'XZ.ps'
fileYZ= basefilename+'YZ.ps'

;Set up labels
xlab= TEXTOIDL('v_x [km s^{-1}]')
ylab= TEXTOIDL('v_y [km s^{-1}]')
zlab= TEXTOIDL('v_z [km s^{-1}]')

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

;axis stuffz
xspan= max(xrange)-min(xrange)
yspan= max(yrange)-min(yrange)
zspan= max(zrange)-min(zrange)

;Make plots
k_print, filename=fileXY

plot_projected_gaussians, mean, covar, projection[*,*,2], xrange=xrange, $
  yrange=yrange, amp=amp, xlabel=xlab, ylabel=ylab, filename=fileXY, $
  cntrlevels=cntrlevels, denslevels=denslevels, quiet=quiet,grid=grid, log=log

;IF keyword_set(label) THEN legend, [label],/left_legend, box=0
;IF keyword_set(label) THEN legend, [label],pos=[0,1], box=0
IF keyword_set(label) THEN IF xspan GE yspan THEN legend, [label],pos=[0,yspan/double(xspan)-0.03], box=0 $
  ELSE legend, [label],pos=[xspan/double(yspan)-0.03,1], box=0

k_end_print

;STOP

k_print, filename=fileXZ

plot_projected_gaussians, mean, covar, projection[*,*,1], xrange=xrange, $
  yrange=zrange, amp=amp, xlabel=xlab, ylabel=zlab, filename=fileXZ, $
  cntrlevels=cntrlevels, denslevels=denslevels, quiet=quiet, grid=grid, log=log

IF keyword_set(label) THEN IF xspan GE zspan THEN legend, [label],pos=[0,zspan/double(xspan)-0.03], box=0 $
  ELSE legend, [label],pos=[xspan/double(zspan)-0.03,1], box=0

k_end_print


k_print, filename=fileYZ

plot_projected_gaussians, mean, covar, projection[*,*,0], xrange=yrange, $
  yrange=zrange, amp=amp, xlabel=ylab, ylabel=zlab, filename=fileYZ, $
  cntrlevels=cntrlevels, denslevels=denslevels, quiet=quiet,grid=grid, log=log

IF keyword_set(label) THEN IF yspan GE zspan THEN legend, [label],pos=[0,zspan/double(yspan)-0.03], box=0 $
  ELSE legend, [label],pos=[yspan/double(zspan)-0.03,1], box=0

k_end_print



END
