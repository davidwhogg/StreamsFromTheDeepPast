pro plot_qso_pm
qsopm= mrdfits('../../data/qsopm*.fits',1)
seed= -1
nqso= n_elements(qsopm)
qsopm.pmra= qsopm.pmra+randomu(seed,nqso)
qsopm.pmdec= qsopm.pmdec+randomu(seed,nqso)
set_plot, 'ps'
device, filename='qso_pm.ps', $
  xsize=6.5,ysize=6.5,xoffset=1.0,yoffset=2.0, $
  /inches,/color
xrange= [-1,1]*50.0
hogg_plothist, qsopm.pmra,npix=100,xrange=xrange,xtitle='pm RA'
hogg_plothist, qsopm.pmdec,npix=100,xrange=xrange,xtitle='pm Dec'
plot, qsopm.pmra,qsopm.pmdec,psym=1,symsize=0.2,thick=4, $
  xrange=xrange,xtitle='pm RA', $
  yrange=xrange,ytitle='pm Dec'
device,/close
return
end
