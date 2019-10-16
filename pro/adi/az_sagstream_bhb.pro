pro az_sagstream_bhb
;+
;Adi Zolotov
;
;takes in spectro file with stars that make Sirko's BHB cuts (this
;should be changed to the strict cuts)
;converts (l,b) to (lambda,beta) defined for sag dwarf
;converts helio velocity to gsr velocity
;
;-

;read in sgr models 
readcol,"/global/data/scr/adi/prolate.dat",lbd,bta,l,b,xg,yg,zg,xs $
  ,ys,zs,u,v,w,d,vgsr,pcol
readcol,"/global/data/scr/adi/oblate.dat",lbd2,bta2,l2,b2,xg2,yg2,zg2,xs2 $
  ,ys2,zs2,u2,v2,w2,d2,vgsr2,pcol2

; for the models:
g_bhb=5*(alog10(d*10^3))-5+0.7
g_bs=5*(alog10(d*10^3))-5+2.7
mod_d_bhb1=where((g_bhb gt 18.4) and (g_bhb lt 19.5))
mod_d_bs1=where((g_bs gt 18.4) and (g_bs lt 19.5))
mod_d_bhb2=where((g_bhb gt 17.9) and (g_bhb lt 19.0))
mod_d_bs2=where((g_bs gt 17.9) and (g_bs lt 19.))
mod_d_bhb3=where((g_bhb gt 17.4) and (g_bhb lt 18.5))
mod_d_bs3=where((g_bs gt 17.4) and (g_bs lt 18.5))
mod_d_bhb4=where((g_bhb gt 16.9) and (g_bhb lt 18.))
mod_d_bs4=where((g_bs gt 16.9) and (g_bs lt 18.))
mod_d_bhb5=where((g_bhb gt 16.4) and (g_bhb lt 17.5))
mod_d_bs5=where((g_bs gt 16.4) and (g_bs lt 17.5))
mod_d_bhb6=where((g_bhb gt 15.9) and (g_bhb lt 17.))
mod_d_bs6=where((g_bs gt 15.9) and (g_bs lt 17.))

BHB=mrdfits('/global/data/scr/adi/BHBspectro.fits',1)


sdist1=where((BHB.g gt 18.8) and (BHB.g lt 19.2),nsdist)
sdist2=where((BHB.g gt 18.1) and (BHB.g lt 18.5),nsdist2)
print,nsdist
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

x=rad*cos((BHB.gl*!DPI)/180)*cos((BHB.gb*!DPI)/180)
y=rad*sin((BHB.gl*!DPI)/180)*cos((BHB.gb*!DPI)/180)
z=rad*sin((BHB.gb*!DPI)/180)

xs=a11*x+a12*y+a13*z
ys=a21*x+a22*y+a23*z
zs=a31*x+a32*y+a33*z

dist=sqrt((xs*xs)+(ys*ys)+(zs*zs))

lam=atan(ys,xs)
lambda=lam/deg2ra
for ii=0,nelem-1 do begin
if (lambda[ii] gt 360)then (lambda[ii]=lambda[ii]-360)
;if ((lambda[ii] gt 0) and (lambda[ii] lt 180)) then (lambda[ii]=lambda[ii]+180)
if (lambda[ii] lt 0)then (lambda[ii]=lambda[ii]+360)
endfor
beta=(acos(zs/sqrt(xs*xs+ys*ys+zs*zs)))/deg2ra
beta=beta-90

betacut=where((beta gt -10.) and (beta lt 10.))
sag=where((BHB.g gt 18.7) and (BHB.g lt 19.2) $
          and ( beta gt -9.) and (beta lt 9.),nsagg)
offsag=where((BHB.g gt 18.9) and (BHB.g lt 19.2) $
             and (BHB.ra gt 140) and (BHB.ra lt 200)$
             and (BHB.dec gt -3.0) and (BHB.dec lt 10.0),noff)
sag1=where((BHB.g gt 18.2) and (BHB.g lt 18.7) $
           and ( beta gt -9.) and (beta lt 9.))  
sag2=where((BHB.g gt 17.7) and (BHB.g lt 18.2) $
          and ( beta gt -9.) and (beta lt 9.))       
sag3=where((BHB.g gt 17.2) and (BHB.g lt 17.7) $
          and ( beta gt -9.) and (beta lt 9.))
sag4=where((BHB.g gt 16.7) and (BHB.g lt 17.2) $
           and ( beta gt -9.) and (beta lt 9.))
sag5=where((BHB.g gt 16.2) and (BHB.g lt 16.7) $
           and ( beta gt -9.) and (beta lt 9.))
;******polar coordinates & plot


;radial=sqrt((BHB.g*BHB.g)+((!dpi*lambda/180)*(!dpi*lambda/180)))
;theta=atan(BHB.g,lambda)
vel1=where((BHB.velocity gt -200) and (BHB.velocity lt -100))
vel2=where((BHB.velocity gt -100) and (BHB.velocity lt 00))
vel3=where((BHB.velocity gt 0) and (BHB.velocity lt 100))
vel4=where((BHB.velocity gt 100) and (BHB.velocity lt 200))
vel5=where((BHB.velocity gt 0) and (BHB.velocity lt 50))
vel6=where((BHB.velocity gt 50) and (BHB.velocity lt 100))
vel7=where((BHB.velocity gt 100) and (BHB.velocity lt 150))
vel8=where((BHB.velocity gt 150) and (BHB.velocity lt 200))

;polar coordinates definition
theta=dblarr(1,nelem)
lam_rad=(!dpi*lambda/180)
xx=(Bhb.g-14) *cos(lam_rad)
yy=(bhb.g-14)*sin(lam_rad)

begplot,"polar.ps"
!x.ticks=15
!x.tickname=['14','15','16','17','18','19','20','21','14','15','16','17','18','19','20','21']
hogg_scatterplot,XX,yy,xrange=[-6,6],yrange=[-6,6],/nocontours,grid=grid,darkest=70,satfrac=0.05,xnpix=120,ynpix=120,xsty=4,ysty=4
axis,0,0,xax=0,/data
axis,0,0,yax=0,/data
plot,(BHB[vel1].g-14),(!dpi*lambda[vel1]/180),/polar,psym=4,symsize=0.2,xrange=[-6,6],yrange=[-6,6]
;legend,['-200<vel<-100','-100<vel<0','0<vel<100','100<vel<200','leading arm'],psym=[4,4,4,4,8],color=[djs_icolor('black'),djs_icolor('red'),djs_icolor('blue'),djs_icolor('green'),djs_icolor('yellow')]

oplot,BHB[vel2].g-14,(!dpi*lambda[vel2]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('red')
oplot,BHB[vel3].g-14,(!dpi*lambda[vel3]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('blue')
;oplot,BHB[vel4].g-14,(!dpi*lambda[vel4]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('green')
oplot,[5],[!dpi*200/180],psym=8,/polar,thick=4.0,color=djs_icolor('red')
oplot,[3.2],[!dpi*200/180],psym=8,/polar,thick=4.0,color=djs_icolor('pink')
oplot,[5],[!dpi*280/180],psym=8,/polar,thick=4.0,color=djs_icolor('yellow')
oplot,[5.25],[!dpi*100/180],psym=8,/polar,thick=4.0,color=djs_icolor('yellow')
oplot,[3.3],[!dpi*190/180],psym=8,/polar,thick=4.0,color=djs_icolor('yellow')
endplot
stop
mvel1=where((vgsr gt -200) and (vgsr lt -100))
mvel2=where((vgsr gt -100) and (vgsr lt 00))
mvel3=where((vgsr gt 0) and (vgsr lt 100))
mvel4=where((vgsr gt 100) and (vgsr lt 200))

begplot,"polar_models.ps",/color
plot,g_bhb[mvel1]-14,(!dpi*lbd[mvel1]/180),/polar,psym=4,symsize=0.2
oplot,g_bhb[mvel2]-14,(!dpi*lbd[mvel2]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('red')
oplot,g_bhb[mvel3]-14,(!dpi*lbd[mvel3]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('blue')
oplot,g_bhb[mvel4]-14,(!dpi*lbd[mvel4]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('green')
endplot
;stop

begplot,"vel_lambda_set.ps"
for ii=0,24 do begin
d_sag=where((BHB.g gt (14.5+ii*0.25)) and (BHB.g lt (15+ii*0.25)) and (beta lt 15) and (beta gt -15))
a=strn(14.5+ii*0.25,length=5)
b=strn(15+ii*0.25,length=5)
plot,lambda[d_sag],BHB[d_sag].velocity,psym=4,xrange=[340,60],yrange=[-300,300],symsize=0.5,xtitle="lambda",ytitle="radial velocity",title=a+'< g < '+b+' beta cut +-15'
endfor
endplot



begplot,"vel_ra_set.ps"
for ii=0,24 do begin
d_sag=where((BHB.g gt (14.5+ii*0.25)) and (BHB.g lt (15+ii*0.25)) and (beta lt 15) and (beta gt -15))
a=strn(14.5+ii*0.25,length=5)
b=strn(15+ii*0.25,length=5)
plot,BHB[d_sag].ra,BHB[d_sag].velocity,psym=4,xrange=[340,10],yrange=[-300,300],symsize=0.5,xtitle="RA",ytitle="radial velocity",title=a+'< g < '+b+' beta cut +-15'
endfor
endplot

rabeta=where((beta gt -10.) and (beta lt 10.)and (BHB.ra gt 130) and (BHB.ra lt 180))

begplot,'g_v_ra.ps'
plot,BHB[rabeta].velocity,BHB[rabeta].g,psym=4,xrange=[-400,400],yrange=[15,21.5]
endplot

on=BHB[sag]
on1=BHB[sag1]
on2=BHB[sag2]
on3=BHB[sag3]
off=BHB[offsag]

correct=(sin(BHB.gl*deg2ra)*cos(BHB.gb*deg2ra)*225)+(7*sin(BHB.gb*deg2ra))-(10*cos(BHB.gl*deg2ra)*cos(BHB.gb*deg2ra))

gsr=BHB.velocity+correct

gsr1=where((gsr gt -200) and (gsr lt -100))
gsr2=where((gsr gt -100) and (gsr lt 00))
gsr3=where((gsr gt 0) and (gsr lt 100))
gsr4=where((gsr gt 100) and (gsr lt 200))


begplot,"polar_gsr.ps",/color
;plot,(BHB[betacut].g-14),(!dpi*lambda[betacut]/180),/polar,psym=4,symsize=0.4,xrange=[-6,6]
;axis,0,0,xax=0
;axis,0,0,yax=0
plot,(BHB[gsr1].g-14),(!dpi*lambda[gsr]/180),/polar,psym=4,symsize=0.2,xrange=[-6,6],yrange=[-6,6]
legend,['-200<vel<-100','-100<vel<0','0<vel<100','100<vel<200'],psym=[4,4,4,4],color=[djs_icolor('black'),djs_icolor('red'),djs_icolor('blue'),djs_icolor('green')]

oplot,BHB[gsr2].g-14,(!dpi*lambda[gsr2]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('red')
oplot,BHB[gsr3].g-14,(!dpi*lambda[gsr3]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('blue')
oplot,BHB[gsr4].g-14,(!dpi*lambda[gsr4]/180),/polar,psym=4,symsize=0.2,color=djs_icolor('green')
endplot

begplot,"g_lambda_set.ps"
for ii=0,38 do begin
    vel=where((BHB.velocity gt (-300.+15.*ii)) and (BHB.velocity lt (-270.+15.*ii)) and (beta gt -15) and (beta lt 15))
    c=strn(-300+15*ii,length=4)
    e=strn(-270+15*ii,length=4)
    plot,lambda[vel],BHB[vel].g,psym=4,xrange=[340,60],yrange=[14,20.5],$
  title=c+'< vel <'+e+'beta cut +-15',xtitle="lambda",ytitle="g mag"
endfor
endplot


begplot,"g_RA_set.ps"
for ii=0,38 do begin
    vel=where((BHB.velocity gt (-300.+15.*ii)) and (BHB.velocity lt (-270.+15.*ii)) and (beta gt -15) and (beta lt 15))
    c=strn(-300+15*ii,length=4)
    e=strn(-270+15*ii,length=4)
    plot,BHB[vel].ra,BHB[vel].g,psym=4,xrange=[340,10],yrange=[14,20.5],$
  title=c+'< vel <'+e+'beta cut +-15',xtitle="lambda",ytitle="g mag"
endfor
endplot



;stop
begplot,"sagstream.ps"
plot,lambda[sdist1],beta[sdist1],psym=4,xrange=[0,360],yrange=[-90,90]
;oplot,lambda[sag],beta[sag],psym=6,color=200
plot,BHB[sdist1].gl,BHB[sdist1].gb,psym=4,xrange=[0,360],yrange=[-90,90]
;oplot,BHB[sag].gl,BHB[sag].gb,psym=6,color=200
plot,on.gl,on.gb,psym=4
plot,lambda[sag1],rad[sag1]/(10^3),psym=4,yrange=[30,55]
plot,on.gl,on.gb,psym=4
plot,lambda[betacut],BHB[betacut].g,psym=1,symsize=0.5,xrange=[360,0],yrange=[15,20]
endplot
;stop

begplot,"g_sag1.ps"
!P.Multi=[0,2,2]
multiplot
plot,lbd[mod_d_bhb1],vgsr[mod_d_bhb1],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],ytitle="Vgsr (km/s)"
oplot,lambda[sag],gsr[sag],symsize=0.4,psym=4
legend,"BHB"
multiplot
plot,lbd[mod_d_bs1],vgsr[mod_d_bs1],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag],gsr[sag],symsize=0.4,psym=4
legend,"BMS"
multiplot
plot,lbd[mod_d_bhb1],g_bhb[mod_d_bhb1],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20],xtitle="lambda",ytitle="g mag"
oplot,lambda[sag],BHB[sag].g,psym=4,symsize=0.4
multiplot
plot,lbd[mod_d_bs1],g_bs[mod_d_bs1],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20]
oplot,lambda[sag],BHB[sag].g,psym=4,symsize=0.4
multiplot,/default
endplot


begplot,"g_sag2.ps"
!P.Multi=[0,2,2]
multiplot
plot,lbd[mod_d_bhb2],vgsr[mod_d_bhb2],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],ytitle="Vgsr"
oplot,lambda[sag1],gsr[sag1],symsize=0.4,psym=4
legend,"BHB"
multiplot
plot,lbd[mod_d_bs2],vgsr[mod_d_bs2],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag1],gsr[sag1],symsize=0.4,psym=4
legend,"BMS"
multiplot
plot,lbd[mod_d_bhb2],g_bhb[mod_d_bhb2],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20],xtitle="lambda",ytitle="g mag"
oplot,lambda[sag1],BHB[sag1].g,psym=4,symsize=0.4
multiplot
plot,lbd[mod_d_bs2],g_bs[mod_d_bs2],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20]
oplot,lambda[sag1],BHB[sag1].g,psym=4,symsize=0.4
multiplot,/default
endplot




begplot,"g_sag3.ps"
!P.Multi=[0,2,2]
multiplot
plot,lbd[mod_d_bhb3],vgsr[mod_d_bhb3],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],ytitle="Vgsr"
oplot,lambda[sag2],gsr[sag2],symsize=0.4,psym=4
legend,"BHB"
multiplot
plot,lbd[mod_d_bs3],vgsr[mod_d_bs3],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag2],gsr[sag2],symsize=0.4,psym=4
legend,"BMS"
multiplot
plot,lbd[mod_d_bhb3],g_bhb[mod_d_bhb3],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20],ytitle="g mag",xtitle="lambda"
oplot,lambda[sag2],BHB[sag2].g,psym=4,symsize=0.4
multiplot
plot,lbd[mod_d_bs3],g_bs[mod_d_bs3],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20]
oplot,lambda[sag2],BHB[sag2].g,psym=4,symsize=0.4
multiplot,/default
endplot


begplot,"g_sag4.ps"
!P.Multi=[0,2,2]
multiplot
plot,lbd[mod_d_bhb4],vgsr[mod_d_bhb4],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],ytitle="Vgsr"
oplot,lambda[sag3],gsr[sag3],symsize=0.4,psym=4
legend,"BHB"
multiplot
plot,lbd[mod_d_bs4],vgsr[mod_d_bs4],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag3],gsr[sag3],symsize=0.4,psym=4
legend,"BMS"
multiplot
plot,lbd[mod_d_bhb4],g_bhb[mod_d_bhb4],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20],ytitle="g mag",xtitle="lambda"
oplot,lambda[sag3],BHB[sag3].g,psym=4,symsize=0.4
multiplot
plot,lbd[mod_d_bs4],g_bs[mod_d_bs4],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20]
oplot,lambda[sag3],BHB[sag3].g,psym=4,symsize=0.4
multiplot,/default
endplot



begplot,"g_sag5.ps"
!P.Multi=[0,2,2]
multiplot
plot,lbd[mod_d_bhb5],vgsr[mod_d_bhb5],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],ytitle="Vgsr"
oplot,lambda[sag4],gsr[sag4],symsize=0.4,psym=4
legend,"BHB"
multiplot
plot,lbd[mod_d_bs5],vgsr[mod_d_bs5],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag4],gsr[sag4],symsize=0.4,psym=4
legend,"BMS"
multiplot
plot,lbd[mod_d_bhb5],g_bhb[mod_d_bhb5],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20],ytitle="g mag",xtitle="lambda"
oplot,lambda[sag4],BHB[sag4].g,psym=4,symsize=0.4
multiplot
plot,lbd[mod_d_bs5],g_bs[mod_d_bs5],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20]
oplot,lambda[sag4],BHB[sag4].g,psym=4,symsize=0.4
multiplot,/default
endplot


begplot,"g_sag6.ps"
!P.Multi=[0,2,2]
multiplot
plot,lbd[mod_d_bhb6],vgsr[mod_d_bhb6],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],ytitle="Vgsr"
oplot,lambda[sag5],gsr[sag5],symsize=0.4,psym=4
legend,"BHB"
multiplot
plot,lbd[mod_d_bs6],vgsr[mod_d_bs6],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag5],gsr[sag5],symsize=0.4,psym=4
legend,"BMS"
multiplot
plot,lbd[mod_d_bhb6],g_bhb[mod_d_bhb6],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20],ytitle="g mag",xtitle="lambda"
oplot,lambda[sag5],BHB[sag5].g,psym=4,symsize=0.4
multiplot
plot,lbd[mod_d_bs6],g_bs[mod_d_bs6],psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[14,20]
oplot,lambda[sag5],BHB[sag5].g,psym=4,symsize=0.4
multiplot,/default
endplot



begplot,'forgroup.ps'
!P.Multi=[0,2,3]
multiplot
plot,lbd,vgsr,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],$
  ytitle="vgsr"
oplot,lambda[sag],gsr[sag],symsize=0.4,psym=4
legend,"d~45kpc"
multiplot
plot,lbd,vgsr,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag2],gsr[sag2],symsize=0.4,psym=4
legend,'d~33kpc'
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)'
legend,"BHB"
oplot,lambda[sag],rad[sag]/1d3,psym=4,symsize=0.4
multiplot

plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)'
oplot,lambda[sag2],rad[sag2]/1d3,psym=4,symsize=0.4
legend,"BHB"
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)',xtitle="prolate model"
legend,'BS stars'
oplot,lambda[sag],rad_bs[sag]/1d3,psym=6,symsize=0.4
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],xtitle="lambda",ytitle='d(kpc)'
legend,'BS stars'
oplot,lambda[sag2],rad_bs[sag2]/1d3,psym=6,symsize=0.4

multiplot,/default
endplot

begplot,"forgroup_oblate.ps"

!P.Multi=[0,2,3]
multiplot
plot,lbd2,vgsr2,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],$
  ytitle="vgsr"
oplot,lambda[sag],gsr[sag],symsize=0.4,psym=4
legend,"d~45kpc"
multiplot
plot,lbd2,vgsr2,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300]
oplot,lambda[sag2],gsr[sag2],symsize=0.4,psym=4
legend,'d~33kpc'
multiplot
plot,lbd2,d2,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)'
legend,"BHB"
oplot,lambda[sag],rad[sag]/1d3,psym=4,symsize=0.4
multiplot

plot,lbd2,d2,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)'
oplot,lambda[sag2],rad[sag2]/1d3,psym=4,symsize=0.4
legend,"BHB"
multiplot
plot,lbd2,d2,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)',xtitle="oblate model"
legend,'BS stars'
oplot,lambda[sag],rad_bs[sag]/1d3,psym=6,symsize=0.4
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],xtitle="lambda",ytitle='d(kpc)'
legend,'BS stars'
oplot,lambda[sag2],rad_bs[sag2]/1d3,psym=6,symsize=0.4

multiplot,/default
endplot

begplot,'sag_models.ps'
!P.Multi=[0,0,3]
multiplot
plot,lbd,vgsr,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],$
  ytitle="vgsr"
oplot,lambda[betacut],gsr[betacut],symsize=0.4,psym=4
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)'
legend,"BHB"
oplot,lambda[betacut],rad[betacut]/1d3,psym=4,symsize=0.4
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],xtitle="lambda",ytitle='d(kpc)'
legend,'BS stars'
oplot,lambda[betacut],rad_bs[betacut]/1d3,psym=6,symsize=0.4

multiplot,/default
endplot

begplot,"models_sag.ps"

!P.Multi=[0,0,3]
multiplot
plot,lbd,vgsr,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[-300,300],$
  ytitle="vgsr"
oplot,lambda[sag],gsr[sag],symsize=0.4,psym=4
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],ytitle='d(kpc)'
legend,"BHB"
oplot,lambda[sag],rad[sag]/1d3,psym=4,symsize=0.4
multiplot
plot,lbd,d,psym=1,symsize=0.5,color=200,xrange=[360,0],yrange=[0,100],xtitle="lambda",ytitle='d(kpc)'
legend,'BS stars'
oplot,lambda[sag],rad_bs[sag]/1d3,psym=6,symsize=0.4

multiplot,/default
endplot

begplot,'sag_dist_comp.ps'
!P.multi=[0,2,2]
multiplot
plot,lambda[sag],gsr[sag],psym=4,xrange=[90,310],yrange=[-400,400],ytitle='velocity(gsr)'

legend,'39-50kpc'
multiplot
plot,lambda[sag1],gsr[sag1],psym=4,xrange=[90,310],yrange=[-600,600]
legend,'17-26kpc'
multiplot
plot,lambda[sag2],gsr[sag2],psym=4,xrange=[90,310],yrange=[-400,400]
legend,'30-36kpc'
multiplot
plot,lambda[sag3],gsr[sag3],psym=4,xrange=[90,310],yrange=[-400,400],xtitle='lambda'
legend,'36-44kpc'
multiplot,/default
endplot

begplot,'sag_loc_comp.ps'
!P.multi=[0,0,2]
multiplot
plot,lambda[sag],beta[sag],psym=4,xrange=[00,360],yrange=[-40,20],ytitle="Beta"
multiplot
plot,lambda[offsag],beta[offsag],psym=4,xrange=[00,360],yrange=[-40,20],xtitle='lambda'
multiplot,/default
endplot

begplot,'sag_off_vel.ps'
!P.multi=[0,0,2]
multiplot
plot,lambda[sag],gsr[sag],psym=4,xrange=[0,320],yrange=[-400,400],ytitle="gsr velocity"
multiplot
plot,lambda[offsag],gsr[offsag],psym=4,xrange=[0,320],yrange=[-400,400],xtitle='lambda'
multiplot,/default
endplot

begplot,'sag_paper.ps'
plot,lambda[sag],beta[sag],psym=4,xrange=[00,360],yrange=[-20,20],ytitle="Beta",xtitle="lambda"
endplot

begplot,"sag_lb.ps"
plot,BHB.gl,BHB.gb,psym=4,xrange=[0,360],yrange=[-90,90],xtitle="gl",ytitle="gb"
endplot

;begplot,"birfurcation.ps"
;!P.Multi=[0,2,2]
;multiplot
;plot,BHB[d_sag1].ra,BHB[d_sag1].dec,psym=4,xrange=[360,0],symsize=0.4;,color=200
;legend,'19.2-19.7'
;multiplot
;plot,BHB[d_sag2].ra,BHB[d_sag2].dec,psym=4,symsize=0.4;,color=300
;legend,'18.7-19.2'
;multiplot
;plot,BHB[d_sag3].ra,BHB[d_sag3].dec,psym=4,symsize=0.4;,color=400
;legend,'18.2-18.7'
;multiplot
;plot,BHB[d_sag4].ra,BHB[d_sag4].dec,psym=4,symsize=0.4;,color=500
;legend,'17.7-18.2'
;multiplot,/default
;oplot,BHB[mod_d_bhb5].ra,BHB[mod_d_bhb5].dec,psym=6;,symsize=0.4,color=600
;oplot,BHB[mod_d_bhb6].ra,BHB[mod_d_bhb6].dec,psym=7;,symsize=0.4,color=700
;endplot

end

