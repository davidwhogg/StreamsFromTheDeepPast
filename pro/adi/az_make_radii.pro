function az_make_radii, nstars=nstars
;+
;Adi Zolotov NYU
;April 2008
;
;-
nstars=nstars
rmin=1.
rmax=40.
seed=-7
z=abs(randomu(seed,nstars))
rr=dblarr(1,nstars)
rr[*]=rmin+exp(alog(rmax/rmin)*z)

set_plot,'ps'
device,file='az_radii.ps',/inches
hogg_plothist,rr,hist=hist,xvec=xvec
plot,alog(xvec),alog(hist),psym=10
device,/close

return,rr
end

