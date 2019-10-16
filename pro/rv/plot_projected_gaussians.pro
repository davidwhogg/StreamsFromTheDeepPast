;+
;   NAME:
;      plot_projected_gaussians
;
;   PURPOSE:
;      plots 1D/2D projections of N-D gaussians
;
;   CALLING SEQUENCE:
;      plot_projected_gaussians, mean, covar, projection, $
;      xrange=xrange, yrange=yrange, cntrlevels=cntrlevels, $ 
;      denslevels=denslevels [, amp=amp, grid=grid, $
;      xlabel=xlabel, ylabel=ylabel,/quiet]
;
;   INPUTS:
;      mean       : means of the gaussians (can be a vector, i.e. 
;                   dblarr(N,ngauss)
;      covar      : covariances of the gaussians (as a vector:
;                   dblarr(N,N,ngauss)
;      projection : projection matrix (NxN); the plot is made with the
;                   axes in the directions of the eigenvectors
;                   corresponding to the non-zero eigenvalues
;                   (i.e. eigenvalue eq 1). This leaves some ambiguity
;                   as to what the eigenvectors are. The routine
;                   allows to give any matrix with only 2 non-zero
;                   eigenvalues. Use the larger eigenvalue for the eigenvector
;                   which will be the x-direction, the smaller
;                   eigenvalue (>0) for the eigenvector which will be
;                   the y-direction.
;      [xy]range  : x-range, etc.
;      cntrlevels : levels at which contours should be drawn
;      denslevels : levels which should be filled in (grayscale)
;
;   OPTIONAL INPUTS:
;      amp        : relative amplitudes of the gaussians (should add up
;                   to 1; default: 1D/ngauss for all)
;      grid       : grid on which to calculate values (grid x grid);
;                   default 100
;      [xy]label  : label for x-axis resp. y-axis
;      quiet      : be quiet
;      + all other keywords allowed by CONTOUR
;      keepbangX  : don't revert to previous graphics settings
;
;   OUTPUTS:
;      for 2D: a contour-plot of the Gaussians
;      for 1D: a regular plot of the projected sum of the Gaussians
;
;   REVISION HISTORY:
;      06-03-08  - Written Bovy
;-
;
PRO PLOT_PROJECTED_GAUSSIANS, mean, covar, projection, xrange=xrange, $
                              yrange=yrange, amp=amp, grid= grid, $
                              xlabel=xlabel, ylabel=ylabel, $
                              cntrlevels=cntrlevels, denslevels=denslevels, $
                              _EXTRA= KeywordsForPlot, quiet=quiet, $
                              keepbangX=keepbangX

IF ~keyword_set(grid) THEN grid= 100L

;Set the higher dimension
d=n_elements(projection[0,*])

;Number of Gaussians in the sum
ngauss= n_elements(mean[0,*])

;Check inputs
IF n_elements(covar[0,0,*]) NE ngauss THEN BEGIN
    splog, 'Error: Number of covariances has to equal the number of means'
    RETURN
ENDIF
IF NOT keyword_set(amp) THEN amp= dblarr(ngauss)+1D/double(ngauss)
IF n_elements(amp) NE ngauss THEN BEGIN
    splog, 'Error: number of given amplitudes should equal the number of means'
    RETURN
ENDIF
IF (abs(1-TOTAL(amp)) GT 1D-7 AND NOT (keyword_set(quiet))) THEN BEGIN
    splog, 'Warning: relative amplitudes dont add up to 1'
    splog, 'They add up to '+strtrim(string(total(amp)),2)
    stop
ENDIF

;Determine the projection, i.e is it 1D or 2D; what directions are
;involved, directions are set by the eigenvectors corresponding to the
;non-zero eigenvalues.
eigenva= EIGENQL(float(projection),/DOUBLE,EIGENVECTORS=eigenve)

;project the means and covariances
np= n_elements(where(eigenva NE 0)); dimension of projected Gaussian
pmean= dblarr(np,ngauss)
pcovar= dblarr(np,np,ngauss)
FOR ii=0L, np-1 DO FOR gg=0L, ngauss-1 DO BEGIN
    pmean[ii,gg]= eigenve[*,ii]##transpose(mean[*,gg])
    FOR jj=ii, np-1 DO BEGIN
        pcovar[ii,jj,gg]= eigenve[*,ii]##covar[*,*,gg]##TRANSPOSE(eigenve[*,jj])
        pcovar[jj,ii,gg]= pcovar[ii,jj,gg]
    ENDFOR
ENDFOR


;Plot either 1D or 2D
IF np EQ 2 THEN BEGIN
    IF NOT keyword_set(quiet) THEN splog, '2D projection'
    delta_x= (xrange[1]-xrange[0])/double(grid)
    delta_y= (yrange[1]-yrange[0])/double(grid)
    x=xrange[0]+delta_x*(dindgen(grid)+0.5D)
    y=yrange[0]+delta_y*(dindgen(grid)+0.5D)
    z=dblarr(grid,grid)
    FOR ii=0L, grid-1 DO FOR jj= 0L, grid-1 DO z[ii,jj]= $
      twod_sum_gaussians([x[ii],y[jj]],pmean,pcovar,amp)
    LOADCT, 0, BOTTOM=0, SILENT=QUIET

;;Save old settings
    bangX= !X
    bangY= !Y
    !X.RANGE= xrange
    !Y.RANGE= yrange
    nticks= 5
    xinterval= hogg_interval(!X.RANGE,nticks=nticks)
    yinterval= hogg_interval(!Y.RANGE,nticks=nticks)
;;Make axes
    IF ~keyword_set(xtickname) THEN $
      djs_plot, [0],[1],xtickinterval=xinterval,ytickinterval=yinterval, $
      /isotropic, _EXTRA= KeywordsForPlot $
      ELSE $
      djs_plot, [0],[1],xrange=xrange,yrange=yrange, $
      /isotropic, _EXTRA= KeywordsForPlot
      


;;Make image
    tvimage= dblarr(grid,grid)
    IF keyword_set(log) THEN BEGIN
        notzero= where(z NE 0D, zcount)
        iszero= where(z EQ 0D,nzcount)
        IF (zcount NE 0) THEN BEGIN
            tvimage[notzero]= alog(z[notzero])
            IF (nzcount NE 0) THEN tvimage[iszero]= alog(min(z[notzero]))
        ENDIF
    ENDIF ELSE BEGIN
        tvimage= z
    ENDELSE
    IF keyword_set(log) THEN BEGIN
        IF ~keyword_set(quiet) THEN splog, 'log-range = ',minmax(tvimage)
        tvrange= [max(tvimage)-10D,1.0*max(tvimage)]
        toosmall= where(tvimage LE tvrange[0])
        IF (toosmall[0] NE -1) THEN tvimage[toosmall]= tvrange[0]
    ENDIF ELSE BEGIN
        tvrange= minmax(tvimage)
    ENDELSE
    tvimage= 255-((floor(256*(tvimage-tvrange[0])/ $
                         (tvrange[1]-tvrange[0])) > 0) $
                  < 255)
    tv, tvimage, !X.RANGE[0],!Y.RANGE[0],/data, $
      xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
      ysize=(!Y.CRANGE[1]-!Y.CRANGE[0])

;;Remake axes, and put axis-labels
    !P.MULTI[0]= !P.MULTI[0]+1
    IF ~keyword_set(xtickname) THEN $
      djs_plot, [0],[1],xtickinterval=xinterval,ytickinterval=yinterval, $
      xtitle=xlabel, ytitle=ylabel, /isotropic, _EXTRA= KeywordsForPlot $
      ELSE $
      djs_plot, [0],[1],xrange=xrange,yrange=yrange, $
      xtitle=xlabel, ytitle=ylabel, /isotropic, _EXTRA= KeywordsForPlot
      

;;Compute contours
    cumindex= reverse(sort(z))
    cumimage= dblarr(grid,grid)
    cumimage[cumindex]= total(z[cumindex],/cumulative)
    cumimage= cumimage*abs(delta_x*delta_y)
    cumimage= cumimage/max(cumimage)

;;split up contours in > .5 and <= .5
    splitindx= 0
    FOR ii=0L, n_elements(cntrlevels)-1 DO $
      IF (cntrlevels[splitindx++] GE 0.5D) THEN BREAK
    splitindx= splitindx-2
    
    cntrs1= dblarr(splitindx+1)
    cntrs2= dblarr(n_elements(cntrlevels)-splitindx-1)
    cntrs1= cntrlevels[0:splitindx]
    cntrs2= cntrlevels[splitindx+1:*]
;;Plot contours
    contour, cumimage,x,y,levels=cntrs1, $
      thick=3,/overplot, color=djs_icolor('white')
    contour, cumimage,x,y,levels=cntrs2, $
      thick=3,/overplot


    IF ~keyword_set(keepbangX) THEN BEGIN 
        !X= bangX
        !Y= bangY
    ENDIF
ENDIF ELSE BEGIN $
  IF np EQ 1 THEN BEGIN
    splog, '1D projection'
ENDIF 
ENDELSE

END
