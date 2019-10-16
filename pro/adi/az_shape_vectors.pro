pro az_shape_vectors
;
;
;
;
;
;

MW1=mrdfits('/global/data/scr/adi/fakedata3.fits',1)
nelem=n_elements(mw1)

;define initial variables

xhat=dblarr(1,3)
xhat[0,0]=1
xhat[0,1]=0
xhat[0,2]=0

yhat=dblarr(1,3)
yhat[0,0]=0
yhat[0,1]=1
yhat[0,2]=0

zhat=dblarr(1,3)
zhat[0,0]=0
zhat[0,1]=0
zhat[0,2]=1



newx=(mw1.x*xhat[0,0])+(mw1.y*xhat[0,1])$
  +(mw1.z*xhat[0,2])
newy=(mw1.x*yhat[0,0])+(mw1.y*yhat[0,1])$
  +(mw1.z*yhat[0,2]);
newz=(mw1.x*zhat[0,0])+(mw1.y*zhat[0,1])$
  +(mw1.z*zhat[0,2])            

;calculate initial radius of each particle in sphere
q=1
s=1
r_n2=(newx^2)+((newy^2)/q)+((newz^2)/s)

set_plot,'ps'
device,file='fakedata_1d2_iter.ps',/inches,xsize=7,ysize=7,/color

nbins=5
rmax=150

rr=dblarr(1,nbins)
rr_in=dblarr(1,nbins)
qq=dblarr(1,nbins)
ss=dblarr(1,nbins)
;lambda1=dblarr(1,nbins)
;lambda2=dblarr(1,nbins)
;lambda3=dblarr(1,nbins)


for jj=1,nbins-1 do begin

lambda1=1
lambda2=1
lambda3=1 
difference=1
rr[jj]=rmax*0.18*jj
rr_in[jj]=rr[jj]-40.

while (difference gt 1d-3) DO BEGIN

;is this particle within rmax squared?
insphere=where((r_n2 lt rr[jj]) and (r_n2 gt rr_in[jj]))
;find variance of particles in ellipsoid
var=az_variance(xx=mw1[insphere].x,yy=mw1[insphere].y,zz=mw1[insphere].z,$
                q=(lambda2/lambda1),s=(lambda3/lambda1))
;find eigenvalues and vectors
eigen=eigenql(var,eigenvectors=evector,/double)

lambda1_new=eigen[0]
lambda2_new=eigen[1]
lambda3_new=eigen[2]

xhat=evector[0,*]
yhat=evector[1,*]
zhat=evector[2,*]

;new rotated coordinate system
newx=(mw1.x*xhat[0,0])+(mw1.y*xhat[0,1])$
  +(mw1.z*xhat[0,2])
newy=(mw1.x*yhat[0,0])+(mw1.y*yhat[0,1])$
  +(mw1.z*yhat[0,2]);
newz=(mw1.x*zhat[0,0])+(mw1.y*zhat[0,1])$
  +(mw1.z*zhat[0,2])            

;set variables
q=lambda2_new/lambda1_new
s=lambda3_new/lambda1_new

;recalculate radius to use for cut
r_n2=(newx^2)+((newy^2)/q)+((newz^2)/s)

difference=abs(lambda1_new-lambda1)

print,(lambda2_new/lambda1_new)
lambda1=lambda1_new
lambda2=lambda2_new
lambda3=lambda3_new


xrange=[-10,10]
yrange=[-10,10]
bin=strn(jj,length=2)

erase & multiplot, [3,3]
plot,mw1[insphere].x,mw1[insphere].y,xrange=xrange,yrange=yrange,psym=3,title='bin'+bin
multiplot
plot,mw1[insphere].x,mw1[insphere].z,xrange=xrange,yrange=yrange,psym=3
multiplot
plot,mw1[insphere].y,mw1[insphere].z,xrange=xrange,yrange=yrange,psym=3
multiplot
plot,newx[insphere],newy[insphere],xrange=xrange,yrange=yrange,psym=3, $
  xtitle='x',ytitle='y'
oplot,xhat,[0,0,0],color=djs_icolor('red'),thick=2
oplot,[0,0,0],yhat,color=djs_icolor('red'),thick=2
multiplot
plot,newx[insphere],newz[insphere],xrange=xrange,yrange=yrange,psym=3, $
  xtitle='x'
oplot,xhat,[0,0,0],color=djs_icolor('red')
oplot,[0,0,0],zhat,color=djs_icolor('red')
multiplot
plot,newy[insphere],newz[insphere],xrange=xrange,yrange=yrange,psym=3, $
  xtitle='y'
oplot,yhat,[0,0,0],color=djs_icolor('red')
oplot,[0,0,0],zhat,color=djs_icolor('red')
multiplot,/reset

endwhile
qq[jj]=sqrt(lambda2/lambda1)
ss[jj]=sqrt(lambda3/lambda1)

;plot q and s

endfor
device,/close

print,sqrt(lambda2/lambda1)

set_plot,'ps'
device,file='fakedata_1d2.ps',/inches,xsize=7,ysize=7,/color
;erase & multiplot,[1,0]
plot,[0],[0],/nodata,xrange=[0,1],yrange=[0,1],xtitle='r/rmax',ytitle='q  (b/a) & s(c/a)'

line=dblarr(200)
line[*]=0.6666667
xline=[0,200]
oplot,xline,line

line2=dblarr(200)
line2[*]=0.333333
xline2=[0,200]
oplot,xline2,line2
plots,rr/rmax,qq,psym=4
plots,rr/rmax,ss,psym=4
device,/close
stop

end
