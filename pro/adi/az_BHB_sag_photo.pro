pro az_BHB_sag_photo
;+
;Adi Zolotov
;April 4, 2007
;Use photometric data (with color cuts for BHBs) to look for sagittarius
;
;
;-

BHB=mrdfits('/global/data/scr/adi/BHB_photo.fits',1)

readcol,"/global/data/scr/adi/prolate.dat",lbd,bta,l,b,xg,yg,zg,xs $
  ,ys,zs,u,v,w,d,vgsr,pcol


g_bhb=5*(alog10(d*10^3))-5+0.7
g_bs=5*(alog10(d*10^3))-5+2.7

;r_bhb=((BHB.g)-0.7+5.)/5.
;rad_bhb=10^(r_bhb)
;dist_bhb=rad_bhb/(1d3)
glactc,BHB.ra,BHB.dec,2000,gl,gb,1,/degree
gl=gl
gb=gb

nelem=n_elements(BHB)
r=dblarr(1,nelem)
rad=dblarr(1,nelem)
x=dblarr(1,nelem)
y=dblarr(1,nelem)
z=dblarr(1,nelem)
dist=dblarr(1,nelem)
lambda=dblarr(1,nelem)
lam=dblarr(1,nelem)
beta=dblarr(1,nelem)
correct=dblarr(1,nelem)


deg2ra=double(!DPI/180)

phi=double((180+3.75)*deg2ra)
theta=double((90-13.46)*deg2ra)
psi=double((180+14.111534)*deg2ra)


a11=double(cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi))
a12=double(cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi))
a13=double(sin(psi)*sin(theta))
a21=double(-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi))
a22=double(-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi))
a23=double(cos(psi)*sin(theta))
a31=double(sin(theta)*sin(phi))
a32=double(-sin(theta)*cos(phi))
a33=double(cos(theta))
                                                   

r=((BHB.g)-0.7+5.)/5.
rad=10^(r)

r_bs=((BHB.g)-2.7+5.)/5
rad_bs=10^(r_bs)

x=rad*cos((gl*!DPI)/180)*cos((gb*!DPI)/180)
y=rad*sin((gl*!DPI)/180)*cos((gb*!DPI)/180)
z=rad*sin((gb*!DPI)/180)

xs=a11*x+a12*y+a13*z
ys=a21*x+a22*y+a23*z
zs=a31*x+a32*y+a33*z

dist=sqrt((xs*xs)+(ys*ys)+(zs*zs))

lam=atan(ys,xs)
lambda=lam/deg2ra

for ii=0L,nelem-1 do begin
if (lambda[ii] gt 360)then (lambda[ii]=lambda[ii]-360)
;if ((lambda[ii] gt 0) and (lambda[ii] lt 180)) then (lambda[ii]=lambda[ii]+180)
if (lambda[ii] lt 0)then (lambda[ii]=lambda[ii]+360)
endfor
beta=(acos(zs/sqrt(xs*xs+ys*ys+zs*zs)))/deg2ra
beta=beta-90


begplot,"lambda_photo_set.ps"
for ii=0,24 do begin
d_sag=where((BHB.g gt (14.5+ii*0.25)) and (BHB.g lt (15+ii*0.25)))
a=strn(14.5+ii*0.25,length=5)
b=strn(15+ii*0.25,length=5)
plot,lambda[d_sag],beta[d_sag],psym=4,xrange=[340,60],yrange=[-40,40],symsize=0.5,xtitle="lambda",ytitle="beta",title=a+'< g < '+b
endfor
endplot

begplot,"radec_photo_set.ps"
for ii=0,24 do begin
d_sag=where((BHB.g gt (14.5+ii*0.25)) and (BHB.g lt (15+ii*0.25)))
a=strn(14.5+ii*0.25,length=5)
b=strn(15+ii*0.25,length=5)
plot,BHB[d_sag].ra,BHB[d_sag].dec,psym=4,xrange=[300,100],yrange=[0,60],symsize=0.5,xtitle="RA",ytitle="DEC",title=a+'< g < '+b
endfor
endplot


betacut=where((beta gt -10) and (beta lt 10))
sag=where((BHB.g gt 18.0) and (BHB.g lt 20.0) $
          and ( beta gt -10.) and (beta lt 10.),nsagg)

;begplot,"birfurcation.ps"
;!P.Multi=[0,2,2]
;multiplot
;plot,BHB[d_sag1].ra,BHB[d_sag1].dec,psym=4,xrange=[220,120],yrange=[0,50],symsize=0.4;,color=200
;legend,'19.2-19.7'
;multiplot
;plot,BHB[d_sag2].ra,BHB[d_sag2].dec,psym=4,symsize=0.4,xrange=[220,120],yrange=[0,50];,color=300
;legend,'18.7-19.2'
;multiplot
;plot,BHB[d_sag3].ra,BHB[d_sag3].dec,psym=4,symsize=0.4,xrange=[220,120],yrange=[0,50];,color=300;,color=400
;legend,'18.2-18.7'
;multiplot
;plot,BHB[d_sag4].ra,BHB[d_sag4].dec,psym=4,symsize=0.4,xrange=[220,120],yrange=[0,50];,color=300;,color=500
;legend,'17.7-18.2'
;multiplot,/default
;oplot,BHB[mod_d_bhb5].ra,BHB[mod_d_bhb5].dec,psym=6;,symsize=0.4,color=600
;oplot,BHB[mod_d_bhb6].ra,BHB[mod_d_bhb6].dec,psym=7;,symsize=0.4,color=700
;endplot


;begplot,'forgroup.ps'
;!P.Multi=[0,1,2]
;multiplot

;plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[35,65],ytitle='d(kpc)'
;legend,"BHB"
;oplot,lambda[sag],rad[sag]/1d3,psym=4,symsize=0.4
;multiplot

;plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,35],ytitle='d(kpc)',xtitle="prolate model"
;legend,'BS stars'
;oplot,lambda[sag],rad_bs[sag]/1d3,psym=6,symsize=0.4

;multiplot,/default
;endplot

;begplot,'canusee.ps'
;!P.Multi=[0,1,2]
;multiplot
;plot,BHB[sag].ra,BHB[sag].dec,xrange=[220,120],yrange=[0,60],psym=4

;endplot


end
