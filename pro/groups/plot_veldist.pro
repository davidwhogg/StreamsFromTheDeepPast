;+
;   NAME:
;      plot_veldist
;   PURPOSE:
;      plot the velocity distribution in vx, vy, annotated
;   CALLING SEQUENCE:
;   INPUT:
;      plotfilename - filename for plot
;   OUTPUT:
;   REVISION HISTORY:
;      2009-11-28 - Written - Bovy (NYU)
;-
PRO PLOT_VELDIST, plotfilename=plotfilename
solutionfilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'

IF ~keyword_set(xrange) THEN xrange= [-130,120]
IF ~keyword_set(yrange) THEN yrange= [-120,60]
lsr= [-10.1,-4.0,-6.7]

restore, solutionfilename

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
cntrlevels=[.02D,.06D,.12D,.21D,.33D,.5D,.68D,.8D,.9D,.95D];,.99D,.999D]
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

xsize=7.75
k_print, filename=plotfilename, xsize=xsize, ysize=xsize/250.*180.
charsize=.6
;;lsr symbol
xarr = [-1,0,1,-1]
yarr = [-1,1,-1,-1]

;;First plot X-Y
position=[0.1,.5,0.656,0.9];;position of top panel
plot_projected_gaussians, mean, covar, projection[*,*,2], xrange=xrange, $
  yrange=yrange, amp=amp, ylabel=ylab, cntrlevels=cntrlevels, $
  denslevels=denslevels, quiet=quiet,grid=400, charsize=charsize, $
  position=position, xlabel=xlab,/keepbangX
;;Add the lsr label
usersym, xarr, yarr,/fill, color=djs_icolor('white')
djs_oplot, [lsr[0]],[lsr[1]], psym=8
usersym, xarr, yarr
djs_oplot, [lsr[0]],[lsr[1]], psym=8
;;Label them
labelcharsize= .8
;;Sirius/Uma
plots, [9.,55], [3.,30.]
legend, ['Sirius/UMa'], pos=[45.,45.],/data, box=0, charsize=labelcharsize
;;NGC 1901
plots, [-23.,-70.], [-11.,20.]
legend, ['NGC 1901'], pos=[-130.,35.],/data, box=0, charsize=labelcharsize
;;Coma Berenices
;plots, [-10.,-70], [-5.,40.]
;legend, ['Coma Berenices'], pos=[-130.,55.],/data, box=0, charsize=labelcharsize
;;Hyades
plots, [-40.,-80], [-19.,-10.]
legend, ['Hyades'], pos=[-130.,5.],/data, box=0, charsize=labelcharsize
;;Pleiades
plots, [-12.,55], [-21.,-50.]
legend, ['Pleiades'], pos=[45.,-45.],/data, box=0, charsize=labelcharsize
;;Hercules
plots, [-20.,-70], [-35.,-75.]
legend, ['Hercules'], pos=[-130.,-67.],/data, box=0, charsize=labelcharsize
;;Arcturus
plots, [2.,55.], [-104.,-90.]
legend, ['Arcturus'], pos=[45.,-75.],/data, box=0, charsize=labelcharsize


k_end_print


    


END
