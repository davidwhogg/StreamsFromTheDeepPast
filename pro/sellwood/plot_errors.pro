;+
;
;-
PRO PLOT_ERRORS

restore, '../groups/phot_gcs.sav'
ngcs= n_elements(gcs.hip)
UVW=dblarr(ngcs,3)
UVWc= dblarr(ngcs,3,3)
for ii=0L, ngcs-1 do begin
    deproject= special_invert(gcs[ii].sm)
    loslb= double([gcs[ii].vr,gcs[ii].vl,gcs[ii].vb])
    UVW[ii,*]= deproject##loslb
    UVWc[ii,*,*]= deproject##(gcs[ii].vrvlvbc##transpose(deproject))
endfor

yrange=[0,1]

k_print, filename='Uerr.ps'
xrange=[0,10]
hogg_plothist, sqrt(UVWc[*,0,0]),/dontplot,xvec=xvec,hist=hist,xrange=xrange
cumul= total(hist,/cumulative)
cumul=cumul/cumul[n_elements(cumul)-1]
djs_plot,xvec,cumul,xrange=xrange,yrange=yrange,$
  xtitle='\sigma_U [km s^{-1}]',ytitle='fraction of stars with uncertainty < \sigma_U'
k_end_print

k_print, filename='Verr.ps'
xrange=[0,10]
hogg_plothist, sqrt(UVWc[*,1,1]),/dontplot,xvec=xvec,hist=hist,xrange=xrange
cumul= total(hist,/cumulative)
cumul=cumul/cumul[n_elements(cumul)-1]
djs_plot,xvec,cumul,xrange=xrange,yrange=yrange,$
  xtitle='\sigma_V [km s^{-1}]',ytitle='fraction of stars with uncertainty < \sigma_V'
k_end_print

k_print, filename='Ufracerr.ps'
xrange=[0,1]
hogg_plothist, sqrt(UVWc[*,0,0])/UVW[*,0],/dontplot,xvec=xvec,hist=hist,xrange=xrange
cumul= total(hist,/cumulative)
cumul=cumul/cumul[n_elements(cumul)-1]
djs_plot,xvec,cumul,xrange=xrange,yrange=yrange,$
  xtitle='\sigma_U/U',ytitle='fraction of stars with frac. uncertainty < \sigma_U/U'
k_end_print


k_print, filename='Vfracerr.ps'
xrange=[0,1]
hogg_plothist, sqrt(UVWc[*,1,1])/UVW[*,1],/dontplot,xvec=xvec,hist=hist,xrange=xrange
cumul= total(hist,/cumulative)
cumul=cumul/cumul[n_elements(cumul)-1]
djs_plot,xvec,cumul,xrange=xrange,yrange=yrange,$
  xtitle='\sigma_V/U',ytitle='fraction of stars with frac. uncertainty < \sigma_V/V'
k_end_print


k_print, filename='pmrafracerr.ps'
xrange=[0,1]
hogg_plothist, gcs.e_pmra/gcs.pmra,/dontplot,xvec=xvec,hist=hist,xrange=xrange
cumul= total(hist,/cumulative)
cumul=cumul/cumul[n_elements(cumul)-1]
djs_plot,xvec,cumul,xrange=xrange,yrange=yrange,$
  xtitle='\sigma_\mu/\mu_\alpha',ytitle='fraction of stars with frac. uncertainty < \sigma_\mu/\mu_\alpha'
k_end_print


k_print, filename='plxfracerr.ps'
xrange=[0,1]
hogg_plothist, gcs.e_plx/gcs.plx,/dontplot,xvec=xvec,hist=hist,xrange=xrange
cumul= total(hist,/cumulative)
cumul=cumul/cumul[n_elements(cumul)-1]
djs_plot,xvec,cumul,xrange=xrange,yrange=yrange,$
  xtitle='\sigma_\pi/\pi',ytitle='fraction of stars with frac. uncertainty < \sigma_\pi/\pi'
k_end_print




end
