;+
;   NAME:
;      plot_convergence
;   PURPOSE:
;      plot the convergence of the EM split and merge etc algorithm
;   CALLING SEQUENCE:
;   INPUT:
;      k      - number of Gaussians
;      w      - regularization parameter
;      plotfilename
;   KEYWORDS:
;   OUTPUT:
;   REVISION HISTORY
;      2008-12-10 - Written Bovy (NYU)
;      2008-12-16 - Plot total log likelihood
;-
PRO PLOT_CONVERGENCE, k, w, sample=sample, plotfilename=plotfilename

IF ~keyword_set(sample) THEN sample=6

IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=4, ysize=4

filename= 'paper0/fitv'+$
  strtrim(string(k),2)+$
  '_bestV_V'+strtrim(string(w,format='(f4.1)'),2)+$
  '_sample'+$
  strtrim(string(sample),2)+$
  '_loglike.log'

readcol, filename, loglike,format='D'

;;Also restore the original data file to find out how many objects
;;this represents (if the sample is not '6' then we need to do
;;something more complicated here)
IF sample EQ 6 THEN BEGIN
    dehnenfilename= '/global/data/hipparcos/hip-dehnen.sav'
    restore, dehnenfilename
    ndata= n_elements(hip.hip)
ENDIF ELSE BEGIN
    ndata= 1
ENDELSE

print, n_elements(loglike)
print, minmax(loglike*ndata)


;;Go through the likelihood to find the locations of the split and
;;merge steps and where they reach the previous level
maxima= dblarr(n_elements(loglike))
maxima_ind= lonarr(n_elements(loglike))
previous_maxima_ind= lonarr(n_elements(loglike))
jj= 0L
FOR ii=1L, n_elements(loglike)-1 DO BEGIN
    IF loglike[ii] LT loglike[ii-1] THEN BEGIN
        maxima_ind[jj]= ii-1
        maxima[jj]= loglike[ii-1]
    ENDIF ELSE IF (loglike[ii] GE maxima[jj] AND loglike[ii-1] LT maxima[jj]) THEN BEGIN
        previous_maxima_ind[jj]= ii
        jj+= 1
    ENDIF
ENDFOR

xs= lindgen(n_elements(loglike))+1L
small= (loglike[n_elements(loglike)-1]*ndata-maxima[0]*ndata)*.05
djs_plot, [0.],[0.], ytitle='ln(likelihood)', yticks=2,$
  xtickformat="(I)", ytickformat="(G11.5)",yminor=4,$
  xtitle='No. of iterations', charsize=.5, $
  yrange=[maxima[0]*ndata-small,loglike[n_elements(loglike)-1]*ndata+small], $
  xrange=[0,max(xs)],thick=4
;;yrange=[-2.255D5,-2.24D5]
breaknow= 0
;;First plot the verticals
FOR ii= 0L, n_elements(loglike)-1 DO BEGIN
    IF maxima[ii+1] EQ 0 THEN BEGIN
        maxima_ind[ii+1]= n_elements(loglike)-1
        breaknow= 1
    ENDIF
    oplotbarx, xs[maxima_ind[ii]], color=djs_icolor('light gray')
    IF breaknow EQ 1 THEN break
ENDFOR
;;Then plot the function, gray part
djs_oplot, xs[0:maxima_ind[0]],ndata*loglike[0:maxima_ind[0]]
breaknow= 0
FOR ii= 0L, n_elements(loglike)-1 DO BEGIN
    IF maxima[ii+1] EQ 0 THEN BEGIN
        maxima_ind[ii+1]= n_elements(loglike)-1
        breaknow= 1
    ENDIF
    djs_oplot, [xs[maxima_ind[ii]],xs[previous_maxima_ind[ii]]], $
      [ndata*loglike[maxima_ind[ii]],ndata*loglike[maxima_ind[ii]]], color=djs_icolor('light gray'),psym=-3
    IF breaknow EQ 1 THEN break
ENDFOR
breaknow= 0
FOR ii= 0L, n_elements(loglike)-1 DO BEGIN
    IF maxima[ii+1] EQ 0 THEN BEGIN
        maxima_ind[ii+1]= n_elements(loglike)-1
        breaknow= 1
    ENDIF
    djs_oplot, xs[previous_maxima_ind[ii]:maxima_ind[ii+1]],$
      ndata*loglike[previous_maxima_ind[ii]:maxima_ind[ii+1]], $
      thick=4
    IF breaknow EQ 1 THEN break
ENDFOR

IF keyword_set(plotfilename) THEN k_end_print

END
