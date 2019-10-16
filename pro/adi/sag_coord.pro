pro sag_coord
;+
;Adi Zolotov
;
;takes in spectro file with stars that make Sirko's BHB cuts (this
;should be changed to the strict cuts)
;converts (l,b) to (lambda,beta) defined for sag dwarf
;converts helio velocity to gsr velocity
;
;-



BHB=mrdfits('/global/data/scr/adi/BHBspectro.fits',1)

sag=where((BHB.g gt 18.9) and (BHB.g lt 19.2) $
          and (BHB.ra gt 200) and (BHB.ra lt 230)$
          and (BHB.dec gt -3.0) and (BHB.dec lt 5.0),nsag)
offsag=where((BHB.g gt 18.9) and (BHB.g lt 19.2) $
          and (BHB.ra gt 140) and (BHB.ra lt 200)$
          and (BHB.dec gt -3.0) and (BHB.dec lt 10.0),noff)

sag1=where((BHB.g gt 16.9) and (BHB.g lt 17.7) $
          and (BHB.ra gt 200) and (BHB.ra lt 230)$
          and (BHB.dec gt -3.0) and (BHB.dec lt 5.0),nsag1)

sag2=where((BHB.g gt 18.1) and (BHB.g lt 18.5) $
          and (BHB.ra gt 200) and (BHB.ra lt 230)$
          and (BHB.dec gt -3.0) and (BHB.dec lt 5.0),nsag2)

sag3=where((BHB.g gt 18.5) and (BHB.g lt 18.9) $
          and (BHB.ra gt 200) and (BHB.ra lt 230)$
          and (BHB.dec gt -3.0) and (BHB.dec lt 5.0),nsag3)

on=BHB[sag]
on1=BHB[sag1]
on2=BHB[sag2]
on3=BHB[sag3]
off=BHB[offsag]


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

x=rad*cos((BHB.gl*!DPI)/180)*cos((BHB.gb*!DPI)/180)
y=rad*sin((BHB.gl*!DPI)/180)*cos((BHB.gb*!DPI)/180)
z=rad*sin((BHB.gb*!DPI)/180)

xs=a11*x+a12*y+a13*z
ys=a21*x+a22*y+a23*z
zs=a31*x+a32*y+a33*z

dist=sqrt((xs*xs)+(ys*ys)+(zs*zs))

lam=atan(ys/xs)
lambda=lam/deg2ra
for ii=0,nelem-1 do begin
if (lambda[ii] gt 360)then (lambda[ii]=lambda[ii]-360)
if ((lambda[ii] gt 0) and (lambda[ii] lt 180)) then (lambda[ii]=lambda[ii]+180)
if (lambda[ii] lt 0)then (lambda[ii]=lambda[ii]+360)
endfor
beta=(asin(zs/sqrt(xs*xs+ys*ys+zs*zs)))/deg2ra


correct=(sin(BHB.gl*deg2ra)*cos(BHB.gb*deg2ra)*225)+(7*sin(BHB.gb*deg2ra))-(10*cos(BHB.gl*deg2ra)*cos(BHB.gb*deg2ra))

gsr=BHB.velocity+correct


begplot,'sag_dist_comp.ps'
!P.multi=[0,2,2]
multiplot
plot,lambda[sag],gsr[sag],psym=4,xrange=[210,310],yrange=[-400,400],ytitle='velocity(gsr)'
legend,'44-50kpc'
multiplot
plot,lambda[sag1],gsr[sag1],psym=4,xrange=[210,310],yrange=[-600,600]
legend,'17-26kpc'
multiplot
plot,lambda[sag2],gsr[sag2],psym=4,xrange=[210,310],yrange=[-400,400]
legend,'30-36kpc'
multiplot
plot,lambda[sag3],gsr[sag3],psym=4,xrange=[210,310],yrange=[-400,400],xtitle='lambda'
legend,'36-44kpc'
multiplot,/default
endplot

begplot,'sag_loc_comp.ps'
!P.multi=[0,0,2]
multiplot
plot,lambda[sag],beta[sag],psym=4,xrange=[200,320],yrange=[-40,20],ytitle="Beta"
multiplot
plot,lambda[offsag],beta[offsag],psym=4,xrange=[200,320],yrange=[-40,20],xtitle='lambda'
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
end

