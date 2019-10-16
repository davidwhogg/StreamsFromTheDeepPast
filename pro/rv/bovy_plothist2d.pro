;+
;   NAME:
;      bovy_plothist2D
;   PURPOSE:
;      make a 2D histogram plot, with contours
;   CALLING SEQUENCE:
;      bovy_plothist2d,...
;   INPUT:
;      xdata  - x data
;      ydata  - y data
;      nbins  - [2] array of bin numbers (default, [16,16]
;      nticks - desired number of ticks on the x axis
;   KEYWORDS:
;      log    - plot log of the histogram
;   OUTPUT:
;   REVISION HISTORY:
;      2008-11-26 - Written Bovy (NYU)
;-
PRO BOVY_PLOTHIST2D, xdata, ydata, nbins=nbins, $
                     _EXTRA= KeywordsForPlot, nticks=nticks, $
                     xrange=xrange, yrange=yrange, log=log

IF ~keyword_set(nbins) THEN nbins= [16,16]
IF ~keyword_set(nticks) THEN nticks= 5
bangX= !X
bangY= !Y
IF ~keyword_set(xrange) THEN BEGIN
    !X.RANGE= minmax(xdata)
    xrange= minmax(xdata)
ENDIF ELSE BEGIN
    !X.RANGE= xrange
ENDELSE
IF ~keyword_set(yrange) THEN BEGIN
    !Y.RANGE= minmax(ydata) 
    yrange= minmax(ydata) 
ENDIF ELSE BEGIN
    !Y.RANGE= yrange
ENDELSE


;;Call hogg_histogram
hist= hogg_histogram([transpose(xdata),transpose(ydata)],[transpose([xrange[0],yrange[0]]), transpose([xrange[1],yrange[1]])],nbins)
xs= dindgen(nbins[0])/double(nbins[0])*(xrange[1]-xrange[0])+(xrange[1]-xrange[0])/double(2*nbins[0])+xrange[0]
ys= dindgen(nbins[1])/double(nbins[1])*(yrange[1]-yrange[0])+(yrange[1]-yrange[0])/double(2*nbins[1])+yrange[0]

IF keyword_set(log) THEN BEGIN
    nonzerohist= where(hist NE 0.)
    zerohist= where(hist EQ 0.)
    hist[nonzerohist]= alog(hist[nonzerohist])
    IF zerohist[0] NE -1 THEN hist[zerohist]= min(hist[nonzerohist])
ENDIF

tvrange= minmax(hist)
tvimage= 255-((floor(256*(hist-tvrange[0])/ $
                         (tvrange[1]-tvrange[0])) > 0) $
                  < 255)
;;Make axes
IF ~keyword_set(xtickinterval) THEN $
  xinterval= hogg_interval(xrange,nticks=nticks) $
  ELSE xinterval= xtickinterval
IF ~keyword_set(ytickinterval) THEN $
  yinterval= hogg_interval(yrange,nticks=nticks) $
  ELSE yinterval= ytickinterval

;;Make axes
djs_plot, [0],[1],xtickinterval=xinterval,ytickinterval=yinterval, $
  _EXTRA= KeywordsForPlot

tv, tvimage, !X.RANGE[0],!Y.RANGE[0],/data, $
      xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
      ysize=(!Y.CRANGE[1]-!Y.CRANGE[0])

;;Remake axes, and put axis-labels
djs_plot, [0],[1],xtickinterval=xinterval,ytickinterval=yinterval, $
  _EXTRA= KeywordsForPlot, /noerase

;;Contour
contour, hist,xs,ys,_EXTRA=KeywordsForPlot, /noerase


;;Revert to old settings
!X= bangX
!Y= bangY

END
