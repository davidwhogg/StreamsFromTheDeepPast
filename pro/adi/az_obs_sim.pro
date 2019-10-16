pro az_obs_sim
;Adi Zolotov
;Jan 21,08
;observe simulation like Carollo et al
;
;

mw1=mrdfits('/global/data/scr/adi2/mw1.fits',1)
main=where((mw1.haloid eq 1));and (mw1.d gt 5) and (mw1.d lt 10))
mw1=mw1[main]
nelem=n_elements(mw1)
dist=sqrt(((mw1.x-8.)*(mw1.x-8.))+(mw1.y*mw1.y)+(mw1.z*mw1.z))

unitv=[[mw1.x-8.],[mw1.y],[mw1.z]]/(replicate(1d0,3)##dist)
vrad=total([[mw1.vx],[mw1.vy],[mw1.vz]]*unitv,2)

dec=asin(mw1.z/dist)
ra=atan((mw1.y/dist),((mw1.x-8.)/dist))
dec=dec*180/!DPI
ra=ra*180/!DPI
for ii=0L,nelem-1 do begin
if (ra[ii] lt 0) then (ra[ii]=ra[ii]+360)
endfor

glactc,ra,dec,2000,gl,gb,1,/degree
vgsr=vrad+10*cos(gl)*cos(gb)+7.2*sin(gb)+225.2*sin(gl)*cos(gb)
gl=(gl*!dpi)/180.
gb=(gb*!dpi)/180.
rr=sqrt((8.-(dist*cos(gb)*cos(gl)))^2+(dist*dist*sin(gb)*sin(gb))+ $
  (dist*dist*cos(gb)*cos(gb)*sin(gl)*sin(gl)))

platelist,plist=str_sdss
spherematch,ra,dec,str_sdss.ra,str_sdss.dec,1.0,mall,m12,dis11,maxmatch=0
rafoot=ra[mall]

decfoot=dec[mall]
sortinx=sort(rafoot)
rafoot_sort=rafoot(sortinx)
decfoot_sort=decfoot(sortinx)
uniq_inx=uniq(rafoot_sort)
rauniq=rafoot_sort(uniq_inx)
decuniq=decfoot_sort(uniq_inx)


begplot,'carollo.ps'
plot,rauniq,decuniq,psym=4
plot,rr[mall],mw1[mall].z,psym=4,xrange=[0,20],yrange=[-15,15]
endplot
stop
end
