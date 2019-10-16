pro vel_ellip
;+
;Adi Zolotov
;July 31,2007
;Using simulations, measure the velocity ellipsoid of the halo  using the variance tensor.
;
;
;-

MW1=mrdfits('/global/data/scr/adi/mw1.fits',1)

mw_halo=where(mw1.comp eq 2)  ; only use stars indentified as halo members
mw1=mw1[mw_halo]


rvir=200       ;virial radius of 200 kpc
xsize=7.5
ysize=7.5

; define the mean of the velocities

set_plot,'ps'
device, file='vel_var.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
plot,[0],[0],/nodata,xrange=[0,200],yrange=[0,1]
for jj=4,8 do begin
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
    radius=jj*0.1*200

;plots,radius,sigmax/sigmay,psym=4
endfor

device,/close
;**************************spatial variance tensor*************
set_plot,'ps'
device, file='shape_var.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
plot,[0],[0],/nodata,xrange=[0,200],yrange=[0,1]

for kk=2,8 do begin

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
    plots,radii,bb/aa,psym=4
endfor
device,/close

end
