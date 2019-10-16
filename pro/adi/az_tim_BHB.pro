pro az_tim_BHB
;+
;Adi Zolotov
;July 10,2007
;
;Read in Beers' BHB catalogue for DR6 and see whats what
;does sag coordinates conversion
;-


readcol,"/global/data/scr/adi/timBHB.dat",id,Teff,eTeff,Logg,elogg,fe,efe,vv,dis,edist,ra,dec
nelem=2490
seed=7
glactc,ra,dec,2000,gl,gb,1,/degree

dist=dis*(1d-3)
deg2ra=!dpi/180
readcol,"/global/data/scr/adi/prolate.dat",lbd,bta,l,b,xg,yg,zg,xs $
  ,ys,zs,u,v,w,d,vgsr,pcol




correct=(sin(gl*deg2ra)*cos(gb*deg2ra)*225)+(7*sin(gb*deg2ra))-(10*cos(gl*deg2ra)*cos(gb*deg2ra))

vel=vv+correct


;*****************************sag coordinates

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
                        

rad=dis*(1d-3)
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
;if ((lambda[ii] gt 0) and (lambda[ii] lt 180)) then
;(lambda[ii]=lambda[ii]+180)
if (lambda[ii] lt 0)then (lambda[ii]=lambda[ii]+360)
endfor
beta=(acos(zs/sqrt(xs*xs+ys*ys+zs*zs)))/deg2ra
beta=beta-90

tvlct,255,255,0,200
begplot,'tim_sag.ps'
plot,lbd,d,psym=4,color=200,xtitle='lambda',ytitle='distance (kpc)',symsize=0.5,xrange=[360,0],yrange=[0,100]
oplot,lambda,dist,psym=4
plot,lbd,vgsr,psym=4,color=200,xtitle='lambda',ytitle='Vgsr',symsize=0.5,xrange=[360,0],yrange=[-300,300]
oplot,lambda,vel,psym=4
endplot

set_plot,'ps'
xszie=7.5 & ysize=7.5
device, file='timbhb.ps',/inches,xsize=xsize,ysize=ysize
plot,ra,dec,psym=4
plot,dist,vel,psym=4
device,/close

xx=dist*cos(ra*!dpi/180)
yy=dist*sin(ra*!dpi/180)

;*******************velocity bins*************************
indx=where((vel gt -200) and (vel lt -100))
indx2=where((vel gt -100) and (vel lt 0))
indx3=where((vel gt 0) and (vel lt 100))
indx4=where((vel gt 100) and (vel lt 200))
color=replicate('black',nelem)
color[indx]=djs_icolor('black')
color[indx2]='red'
color[indx3]='blue'
color[indx4]='green'
sindex=shuffle_indx(nelem,seed=seed)
plotsym,0,0.3,/fill

set_plot,'ps'
device,file='bhb_beers.ps',/color,bits=8
djs_plot,ra[sindex],dist[sindex],psym=8,color=color[sindex],xsty=5,ysty=5,xrange=[0,360],yrange=[0,100],xtitle='ra',ytitle='distance (kpc)',title='Beers BHB'
axis,0,0,xax=0
axis,0,0,yax=0
legend,['200<v<-100','-100<v<0','0<v<100','100<v<200'],psym=[8,8,8,8],color=[djs_icolor('black'),djs_icolor('red'),djs_icolor('blue'),djs_icolor('green')],symsize=[1.5,1.5,1.5,1.5]
device,/close
end
