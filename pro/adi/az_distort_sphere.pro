function az_distort_sphere, nstars=nstars,lambda1=lambda1,lambda2=lambda2,lambda3=lambda3

;+
;Adi Zolotov NYU
;April 23,2008
;This code takes in az_spherical_powerlaw, which produces nstars
;distributed according to a nfw powerlaw (at large distances)
;and distorts the spherical distibution by altering the axis ratios (lambdas)
;and by changing the principal axes orientation (change xhat,yhat,zhat
;by hand for now)
;-

nstars=nstars
;sphere is 3x nstars matrix
sphere=az_spherical_powerlaw(nstars=nstars)
;begin to distort here
lambda1=lambda1
lambda2=lambda2
lambda3=lambda3

xhat=[0.7,0.7,0]
yhat=[0,0.866,0.5]
zhat=[0,0.5,.866]
;xhat=[1,0,0]
;yhat=[0,1,0]
;zhat=[0,0,1]


x1=lambda1*((sphere[0,*]*xhat[0])+(sphere[1,*]*xhat[1])+(sphere[2,*]*xhat[2]))
x2=lambda2*((sphere[0,*]*yhat[0])+(sphere[1,*]*yhat[1])+(sphere[2,*]*yhat[2]))
x3=lambda3*((sphere[0,*]*zhat[0])+(sphere[1,*]*zhat[1])+(sphere[2,*]*zhat[2]))
dis_sphere=dblarr(3,nstars)
dis_sphere[0,*]=x1[*]
dis_sphere[1,*]=x2[*]
dis_sphere[2,*]=x3[*]

;make some plots to check everything
xrange=[-50,50]
yrange=[-50,50]

set_plot,'ps'
device,file='distort_sphere2.ps'
erase & multiplot, [3,1]
plot,dis_sphere[0,*],dis_sphere[1,*],psym=3,symsize=4,/isotropic,xrange=xrange,yrange=yrange
multiplot
plot,dis_sphere[0,*],dis_sphere[2,*],psym=3,symsize=4,/isotropic,xrange=xrange,yrange=yrange
multiplot
plot,dis_sphere[1,*],dis_sphere[2,*],psym=3,symsize=4,/isotropic,xrange=xrange,yrange=yrange
multiplot,/reset
device,/close

return,dis_sphere

end
