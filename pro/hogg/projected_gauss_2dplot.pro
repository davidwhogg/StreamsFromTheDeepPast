;+
; NAME:
;   projected_gauss_2dplot
; PURPOSE:
;   two-d projections of n-d gaussians
; INPUTS:
; OUTPUTS:
; BUGS:
;   this should be linked to projected_gauss_step
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro projected_gauss_2dplot, ydata, ycovar, projection, amp, xmean, xcovar, $
                            filename=filename,nogreyscale=nogreyscale

if(NOT keyword_set(filename)) then $
  filename='qa_project_gauss_2dplot.ps'

ndimy=n_elements(ycovar)/n_elements(ydata)
ndata=n_elements(ydata)/ndimy
ndimx=n_elements(projection)/ndata/ndimy

xdata=dblarr(ndimx,ndata)
for ii=0L, ndata-1L do $
  xdata[*,ii]=projection[*,*,ii]#ydata[*,ii]
weight=dblarr(ndata)+1.

ex_max_plot,weight,xdata,amp,xmean,xcovar,filename,nsig=2, $
  label=['U','V','W'],/nodata,range=[[-100.,100.],[-100.,100.], $
                                     [-40.,40.]],nogreyscale=nogreyscale

end
