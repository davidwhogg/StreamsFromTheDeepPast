pro inertia_tensor 
;+
;Adi Zolotov
;July 31,2007
;Using simulations, measure the velocity ellipsoid of the halo  using the variance tensor.
;
;
;-

MW1=mrdfits('/global/data/scr/adi2/mw1.fits',1)

;mw_halo=where((mw1.haloid eq 0) or (mw1.haloid eq 1))  ; only use stars indentified as halo members
mw_halo=where(mw1.comp eq 2)  ; only use stars indentified as halo members
mw1=mw1[mw_halo]

;start at the full virial radius and define smaller concentric spheres
rvir=268
;for ii=0L,nelem-1 do begin
mom=dblarr(3,3)
mom1=dblarr(3,3)


;initial conditions
aa=1
bb=1
cc=1
q=bb/aa
s=cc/aa
radius=sqrt((mw1.x)^2+(mw1.y/q)^2+(mw1.z/s)^2)
shell=where((radius lt 0.3*rvir),nrad)
position=dblarr(nrad,3)
tot_sph=dblarr(nrad,3,3)
var=dblarr(3,3)
rtot=total(sqrt((mw1[shell].x)^2+(mw1[shell].y/q)^2+(mw1[shell].z/s)^2))

position[*,0]=mw1[shell].x
position[*,1]=mw1[shell].y
position[*,2]=mw1[shell].z


pos_t=transpose(position)
for ii=0L,nrad-1 do begin
    tot_sph[ii,*,*]=position[ii,*]##pos_t[*,ii]
endfor


var[0,0]=total(tot_sph[*,0,0])/rtot
var[0,1]=total(tot_sph[*,0,1])/rtot
var[0,2]=total(tot_sph[*,0,2])/rtot
var[1,0]=total(tot_sph[*,1,0])/rtot
var[1,1]=total(tot_sph[*,1,1])/rtot
var[1,2]=total(tot_sph[*,1,2])/rtot
var[2,0]=total(tot_sph[*,2,0])/rtot
var[2,1]=total(tot_sph[*,2,1])/rtot
var[2,2]=total(tot_sph[*,2,2])/rtot

eigen=eigenql(var,eigenvectors=evector)
aa=sqrt(eigen[0])
bb=sqrt(eigen[1])
cc=sqrt(eigen[2])

xprime0=position[*,0]*evector[0,0]+position[*,1]*evector[0,1]+position[*,2]*evector[0,2]
;xprime[*,0]=[xprime0,xprime1,xprime2]
;for jj=0L,nrad do begin

;xprime[jj]=az_dotp(position[jj,*],evector[0,*])
;endfor
;yprime=az_dotp(position[*,1],evector[1,*])
;zprime=az_dotp(position[*,2],evector[2,*])


stop

;*************************IGNORE everything after this **********************


mom[0,0]=total((mw1[shell].y)^2+(mw1[shell].z)^2)/r
mom[1,1]=total((mw1[shell].x)^2+(mw1[shell].z)^2)/r
mom[2,2]=total((mw1[shell].y)^2+(mw1[shell].x)^2)/r
mom[0,1]=-total((mw1[shell].y)*(mw1[shell].x))/r
mom[1,0]=-total((mw1[shell].y)*(mw1[shell].x))/r
mom[0,2]=-total((mw1[shell].x)*(mw1[shell].z))/r
mom[2,0]=-total((mw1[shell].x)*(mw1[shell].z))/r
mom[1,2]=-total((mw1[shell].y)*(mw1[shell].z))/r
mom[2,1]=-total((mw1[shell].y)*(mw1[shell].z))/r

eigen=eigenql(mom,eigenvector=evector)
if ((eigen[0] gt eigen[1]) and (eigen[0] gt eigen[2])) then a=sqrt(eigen[0])
if ((eigen[1] gt eigen[2]) and (eigen[1] lt eigen[0])) then b=sqrt(eigen[1])
if ((eigen[2] lt eigen[0]) and (eigen[2] lt eigen[1])) then c=sqrt(eigen[2])
q=b/a
s=c/a

r1=total(sqrt((mw1[shell].x)^2+(mw1[shell].y/q)^2+(mw1[shell].z/s)^2))
mom1[0,0]=total((mw1[shell].y)^2+(mw1[shell].z)^2)/r1
mom1[1,1]=total((mw1[shell].x)^2+(mw1[shell].z)^2)/r1
mom1[2,2]=total((mw1[shell].y)^2+(mw1[shell].x)^2)/r1
mom1[0,1]=total((mw1[shell].y)*(mw1[shell].x))/r1
mom1[1,0]=total((mw1[shell].y)*(mw1[shell].x))/r1
mom1[0,2]=total((mw1[shell].x)*(mw1[shell].z))/r1
mom1[2,0]=total((mw1[shell].x)*(mw1[shell].z))/r1
mom1[1,2]=total((mw1[shell].y)*(mw1[shell].z))/r1
mom1[2,1]=total((mw1[shell].y)*(mw1[shell].z))/r1
eigen1=eigenql(mom1,eigenvector=evector1)
if ((eigen1[0] gt eigen1[1]) and (eigen1[0] gt eigen1[2])) then a1=sqrt(eigen1[0])
if ((eigen1[1] gt eigen1[2]) and (eigen1[1] lt eigen1[0])) then b1=sqrt(eigen1[1])
if ((eigen1[2] lt eigen1[0]) and (eigen1[2] lt eigen1[1])) then c1=sqrt(eigen1[2])
q1=b1/a1
s1=c1/a1
stop

q=b/a
s=c/a








rvir=200       ;virial radius of 200 kpc
xsize=7.5
ysize=7.5

; define the mean of the velocities
plotsym,0,0.6,/fill
set_plot,'ps'
device, file='vel_var_stars.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
plot,[0],[0],/nodata,xrange=[0,350],yrange=[0,1]
for jj=2,2 do begin
    cond=where(mw1.d lt jj*0.1*200)

   mw=mw1[cond]
    nelem=n_elements(MW)
    vel=dblarr(nelem,3)
    mean=dblarr(1,3)
    
    vel[*,0]=mw.vx
    vel[*,1]=mw.vy
    vel[*,2]=mw.vz
    mean[0,0]=total(vel[*,0])/nelem
    mean[0,1]=total(vel[*,1])/nelem
    mean[0,2]=total(vel[*,2])/nelem
    
    
    
; define a variance tensor
    var=dblarr(nelem,3)
    tot=dblarr(nelem,3,3)
    variance=dblarr(3,3)
    
    
    var[*,0]=vel[*,0]-mean[0,0]
    var[*,1]=vel[*,1]-mean[0,1]
    var[*,2]=vel[*,2]-mean[0,2]
    
    var_t=transpose(var)
    for ii=0L,nelem-1 do begin
        tot[ii,*,*]=var[ii,*]##var_t[*,ii]
    endfor
    
    variance[0,0]=total(tot[*,0,0])/nelem
    variance[0,1]=total(tot[*,0,1])/nelem
    variance[0,2]=total(tot[*,0,2])/nelem
    variance[1,0]=total(tot[*,1,0])/nelem
    variance[1,1]=total(tot[*,1,1])/nelem
    variance[1,2]=total(tot[*,1,2])/nelem
    variance[2,0]=total(tot[*,2,0])/nelem
    variance[2,1]=total(tot[*,2,1])/nelem
    variance[2,2]=total(tot[*,2,2])/nelem
    
    sigmax=sqrt(variance[0,0])
    sigmay=sqrt(variance[1,1])
    sigmaz=sqrt(variance[2,2])
    radius=jj*0.1*200

    plots,radius,sigmax/sigmay,psym=8,color=djs_icolor('red')
    plots,radius,sigmax/sigmaz,psym=8,color=djs_icolor('black')

endfor

device,/close
;**************************spatial variance tensor*************





set_plot,'ps'
device, file='shape_var_stellar.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
plot,[0],[0],/nodata,xrange=[0,350],yrange=[0,1],xtitle='radius (kpc)',title='halo stars only'

for kk=2,15 do begin

    cond2=where(mw1.d lt kk*0.1*200)
    mw2=mw1[cond2]
    nelem2=n_elements(mw2)
    rho=sqrt(mw1.x*mw1.x+mw1.y*mw1.y+mw1.z*mw1.z)
    phi=acos(mw1.z/rho)
    theta=atan(mw1.y/mw1.x)
    
    
    s_mean=dblarr(1,3)
    
    s_mean[0,0]=total(mw2.x)/nelem2
    s_mean[0,1]=total(mw2.y)/nelem2
    s_mean[0,2]=total(mw2.z)/nelem2
    
; define a variance tensor
    var_sph=dblarr(nelem2,3)
    tot_sph=dblarr(nelem2,3,3)
    variance_sph=dblarr(3,3)
    
    var_sph[*,0]=mw2.x-s_mean[0,0]
    var_sph[*,1]=mw2.y-s_mean[0,1]
    var_sph[*,2]=mw2.z-s_mean[0,2]
    
    var_t_sph=transpose(var_sph)
    for ii=0L,nelem2-1 do begin
        tot_sph[ii,*,*]=var_sph[ii,*]##var_t_sph[*,ii]
        
        
    endfor
    variance_sph[0,0]=total(tot_sph[*,0,0])/nelem2
    variance_sph[0,1]=total(tot_sph[*,0,1])/nelem2
    variance_sph[0,2]=total(tot_sph[*,0,2])/nelem2
    variance_sph[1,0]=total(tot_sph[*,1,0])/nelem2
    variance_sph[1,1]=total(tot_sph[*,1,1])/nelem2
    variance_sph[1,2]=total(tot_sph[*,1,2])/nelem2
    variance_sph[2,0]=total(tot_sph[*,2,0])/nelem2
    variance_sph[2,1]=total(tot_sph[*,2,1])/nelem2
    variance_sph[2,2]=total(tot_sph[*,2,2])/nelem2
    
    eigen=eigenql(variance_sph)
    aa=sqrt(eigen[0])   
    bb=sqrt(eigen[1])
    cc=sqrt(eigen[2])
    radii=kk*0.1*200
    tri=(aa*aa-bb*bb)/(aa*aa-cc*cc) ;triaxial parameter
    plots,radii,bb/aa,psym=8,color=djs_icolor('red')
    plots,radii,cc/aa,psym=8,color=djs_icolor('blue')
    plots,radii,tri,psym=8,color=djs_icolor('black')
    plots,radii,cc/bb,psym=8,color=djs_icolor('green')


endfor
legend,['T','q=b/a','s=c/a','p=c/b'],psym=[8,8,8,8],color=[djs_icolor('black'),djs_icolor('red'),djs_icolor('blue'),djs_icolor('green')],/bottom
device,/close
stop
end
