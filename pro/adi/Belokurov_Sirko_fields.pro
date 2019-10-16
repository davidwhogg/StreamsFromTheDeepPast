pro Belokurov_Sirko_fields

file=mrdfits('/global/data/scr/adi/spectro.fits',1)
nelem=n_elements(file)
str1=create_struct('g',0.,'u',0.,'i',0.,'r',0.,'ra',0.D,'dec',0.D,'velocity',0.D)
mags_struct=replicate(str1,nelem)
mags_struct.g=22.5-(2.5*alog10(file.psfflux[1]))-file.extinction[1]
mags_struct.i=22.5-(2.5*alog10(file.psfflux[3]))-file.extinction[3]
mags_struct.u=22.5-(2.5*alog10(file.psfflux[0]))-file.extinction[0]
mags_struct.r=22.5-(2.5*alog10(file.psfflux[2]))-file.extinction[2]
mags_struct.ra=file.ra
mags_struct.dec=file.dec
mags_struct.velocity=file.z*(3d5)
gmr=mags_struct.g-mags_struct.r
umg=mags_struct.u-mags_struct.g
rmi=mags_struct.r-mags_struct.i
gmi=mags_struct.g-mags_struct.i



cuts=where((gmr gt (-0.4)) and (gmr lt 0.0) and (umg gt 0.8) $
           and (umg lt 1.35),ncut)
mags_struct=mags_struct[cuts]
print, ncut
racen14=150.0
deccen14=18.4
if (ncut gt 1) then begin
cutA= where ((abs(mags_struct.ra-150.0) lt (3.0)) and $
             (abs(mags_struct.dec-18.4) lt 3.0),ncutA)
;cutA= where ((abs(mags_struct.ra-150.0) lt (6.0*cos(18.4*!DPI/180.0))) and $
             ;(abs(mags_struct.dec-18.4) lt 3.0),ncutA)
; (mags_struct.ra lt 150) and (mags_struct.dec gt 18.4) $
          ; and (mags_struct.dec lt 19),ncutA)
cutB= where ((abs(mags_struct.ra-145.0) lt (3.0*cos(19.0*!DPI/180.0))) and $
             (abs(mags_struct.dec-19.0) lt 3.0),ncutB)
cutC= where ((abs(mags_struct.ra-140.0) lt (3.0*cos(19.4*!DPI/180.0))) and $
             (abs(mags_struct.dec-19.4) lt 3.0),ncutC)
cutD= where ((abs(mags_struct.ra-130.0) lt (3.0*cos(19.5*!DPI/180.0))) and $
             (abs(mags_struct.dec-19.5) lt 3.0),ncutD)







if (ncutB ge 1) then begin
print, ncutB
boxcutsB=mags_struct[cutB]
stB=create_struct('gmr',0.,'umg',0.,'rmi',0.,'ra',0.D,'dec',0.,'g',0.,'velocity',0.D, 'gmi',0.,'i',0.)

boxB=replicate(stB,ncutB)
boxB.gmr=boxcutsB.g-boxcutsB.r
boxB.umg=boxcutsB.u-boxcutsB.g
boxB.rmi=boxcutsB.r-boxcutsB.i
boxB.gmi=boxcutsB.g-boxcutsB.i
boxB.ra=boxcutsB.ra
boxB.velocity=boxcutsB.velocity
boxB.dec=boxcutsB.dec
boxB.g=boxcutsB.g
boxB.i=boxcutsB.i
endif

if (ncutA ge 1) then begin
print, ncutA


boxcuts=mags_struct[cutA]
st=create_struct('gmr',0.,'umg',0.,'rmi',0.,'gmi',0.,'i',0.,'ra',0.D,'dec',0.,'g',0.,'velocity',0.D)

box=replicate(st,ncutA)
box.gmr=boxcuts.g-boxcuts.r
box.umg=boxcuts.u-boxcuts.g
box.rmi=boxcuts.r-boxcuts.i
box.gmi=boxcuts.g-boxcuts.i
box.ra=boxcuts.ra
box.velocity=boxcuts.velocity
box.dec=boxcuts.dec
box.g=boxcuts.g
box.i=boxcuts.i
endif

if (ncutC ge 1) then begin
print, ncutC


boxcutsC=mags_struct[cutC]
stC=create_struct('gmr',0.,'umg',0.,'gmi',0.,'i',0.,'rmi',0.,'ra',0.D,'dec',0.,'g',0.,'velocity',0.D)

boxC=replicate(stC,ncutC)
boxC.gmr=boxcutsC.g-boxcutsC.r
boxC.umg=boxcutsC.u-boxcutsC.g
boxC.rmi=boxcutsC.r-boxcutsC.i
boxC.gmi=boxcutsC.g-boxcutsC.i
boxC.ra=boxcutsC.ra
boxC.velocity=boxcutsC.velocity
boxC.dec=boxcutsC.dec
boxC.g=boxcutsC.g
boxC.i=boxcutsC.i
endif

if (ncutD ge 1) then begin
print, ncutD


boxcutsD=mags_struct[cutD]
stD=create_struct('gmr',0.,'umg',0.,'gmi',0.,'i',0.,'rmi',0.,'ra',0.D,'dec',0.,'g',0.,'velocity',0.D)

boxD=replicate(stD,ncutD)
boxD.gmr=boxcutsD.g-boxcutsD.r
boxD.umg=boxcutsD.u-boxcutsD.g
boxD.rmi=boxcutsD.r-boxcutsD.i
boxD.gmi=boxcutsD.g-boxcutsD.i
boxD.ra=boxcutsD.ra
boxD.velocity=boxcutsD.velocity
boxD.dec=boxcutsD.dec
boxD.g=boxcutsD.g
boxD.i=boxcutsD.i
endif



endif


;begplot, 'cmd_bel.ps'
plotsym,0,0.7,/fill
;!P.multi=[0,2,2,0,0] 
;multiplot
;plot,box.gmi, box.i,psym=8,xtitle='g-i',xrange=[-0.2,1.2],yrange=[22,14]
;multiplot
;plot,boxB.gmi,boxB.i, psym=8,xrange=[-0.2,1.2],yrange=[22,14];,xyout='field A15'
;multiplot
;plot,boxC.gmi,boxC.i, psym=8,ytitle='imag',xrange=[-0.2,1.2],yrange=[22,14]
;multiplot
;plot,boxD.gmi,boxD.i, psym=8,xrange=[-0.2,1.2],yrange=[22,14]
;endplot

;plot, box.velocity,box.g, xtitle='A14',psym=8
;plot, boxB.velocity,boxB.g, xtitle='A15',psym=8
;plot, boxC.velocity,boxC.g, xtitle='A16',psym=8
;plot, boxD.velocity,boxD.g, xtitle='A17',psym=8
;plot,box.gmr, box.g,psym=3
totall=total(nelem)
totsirko=total(ncutA)


x1=[-300,300]
y1=[17.2,17.2]


x2=[-300,300]
y2=[19.2,19.2]

y3=[17.15,17.15]
y4=[17.00,17.0]
y5=[16.85,16.85]
y6=[19.15,19.15]
y7=[19.0,19.0]
y8=[18.85,18.85]
begplot, 'Belokurov_Sirko1.ps'

plot,box.velocity,box.g,xtitle='velocity (km/s)',ytitle='g',title='Field A14',psym=8
endplot
begplot, 'Belokurov_Sirko2.ps'
plot,boxB.velocity,boxB.g,xtitle='velocity (km/s)',ytitle='g',title='Field A15',psym=8
endplot
;begplot, 'Belokurov_Sirko3.ps'
;plot,boxC.velocity,boxC.g,xtitle='velocity (km/s)',ytitle='g',title='Field A16',psym=8
;endplot
begplot, 'Belokurov_Sirko4.ps'
plot,boxD.velocity,boxD.g,xtitle='velocity (km/s)',ytitle='g',title='Field A17',psym=8
endplot

begplot, 'allBfields.ps'

yrange=[22,16]
xrange=[-300,300]

!P.multi=[0,2,2] 
multiplot
plot,box.velocity, box.g,psym=8,title='Belokurov Fields With Sirko''s color cut',yrange=[22,14],xrange=[-300,300]
xyouts, -270, 21.5, 'A14'
plots,x1,y1,linestyle=0
plots,x2,y2,linestyle=2
legend,'BHB',linestyle=0
multiplot
;yrange=[22,14]
;xrange=[-300,300]

plot,boxB.velocity,boxB.g, psym=8,yrange=[22,14],xrange=[-300,300]
xyouts,-270,21.5,'A15'
plots,x1,y3,linestyle=0
plots,x2,y6,linestyle=2

legend,'BS',linestyle=2,/right
multiplot
plot,boxC.velocity,boxC.g, psym=8,ytitle='gmag',xtitle='velocity (km/s)',yrange=[22,14],xrange=[-300,300]
xyouts, -270,21.5,'A16'
plots,x1,y4,linestyle=0
plots,x2,y7,linestyle=2
multiplot
plot,boxD.velocity,boxD.g, psym=8,yrange=[22,14],xrange=[-300,300]
xyouts, -270,21.5,'A18'
plots,x1,y5,linestyle=0
plots,x2,y8,linestyle=2
multiplot,/default

endplot


end

