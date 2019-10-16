function az_spherical_powerlaw,nstars=nstars
;+
;Adi Zolotov NYU
;April 2008
;combines az_make_radii(makes radii according to large scale nfw
;powerlaw)with az_make_spherical (unitvectors)
;
;-
nstars=nstars
radii=az_make_radii(nstars=nstars)
sphere=az_make_spherical_vectors(nstars=nstars)
distr=dblarr(3,nstars)
distr[0,*]=radii*sphere[*,0]
distr[1,*]=radii*sphere[*,1]
distr[2,*]=radii*sphere[*,2]

xrange=[-15,15]
yrange=[-15,15]
xsize=10
ysize=10
set_plot,'ps'
device,file='spher_power.ps';,xsize=xsize,ysize=ysize
erase & multiplot, [3,1]
plot,distr[0,*],distr[1,*],psym=3,xrange=xrange,yrange=yrange,symsize=4,/isotropic
multiplot
plot,distr[0,*],distr[1,*],psym=3,xrange=xrange,yrange=yrange,symsize=4,/isotropic
multiplot
plot,distr[0,*],distr[1,*],psym=3,xrange=xrange,yrange=yrange,symsize=4,/isotropic
multiplot,/reset
device,/close


return,distr
end
