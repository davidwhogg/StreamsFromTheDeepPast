function az_clumps,deltav=deltav,radius=radius
;+
;Adi Zolotov
;June 14,2007
;Find clumps in distance (mag,ra,dec) and velocity in SDSS
;spectroscopic data of BHB stars (selected with Sirko's color cuts)
;      
;INPUTS: delta v input is in km/s
;        radius input in kpc  
;
;returns a structure with all parent and child info that matched
;deltav and radius criteria and outputs fits file with the structure saved.
;-

file=mrdfits('/global/data/scr/adi2/BHB_spectro.fits',1)
magcut=where(file.g gt 15.5) ;cut out saturated data at g~15.5
file=file[magcut]
nelem=n_elements(file)
kk=0
;********************define distance to stars using M=0.8, and x,y,z**********
rad=10^((file.g+4.3)/5.)
rad=rad*1d-3  ;distance to stars in kpc

x=rad*cos((file.ra*!DPI)/180)*cos((file.dec*!DPI)/180)
y=rad*sin((file.ra*!DPI)/180)*cos((file.dec*!DPI)/180)
z=rad*sin((file.dec*!DPI)/180)
vec=[x,y,z]

;***************define all arrays*********************************
cen=dblarr(nelem)
xcen=dblarr(nelem)
ycen=dblarr(nelem)
zcen=dblarr(nelem)
rcen=dblarr(nelem)
cenvel=dblarr(nelem)
cvec=dblarr(3,nelem)
veldiff=dblarr(3,nelem)
sphere=dblarr(nelem)
str1=create_struct('parent',0.,'child',0.,'par_g',0.,'par_u',0.,'par_i',0.,'par_r',0.,'par_ra',0.D,'par_dec',0.D,'par_velocity',0.D,'ch_g',0.,'ch_u',0.,'ch_i',0.,'ch_r',0.,'ch_ra',0.D,'ch_dec',0.D,'ch_velocity',0.D,'veldiff',0.,'raddiff',0.)
clump=replicate(str1,100000)

;*************** pick star to be parent************
for ii=0L,nelem-1 do begin
 cen[ii]=ii
 rcen[ii]=10^((file[ii].g+4.3)/5.)
 rcen[ii]=rcen[ii]*(1d-3)
 xcen[ii]=rcen[ii]*cos((file[ii].ra*!DPI)/180)*cos((file[ii].dec*!DPI)/180)
 ycen[ii]=rcen[ii]*sin((file[ii].ra*!DPI)/180)*cos((file[ii].dec*!DPI)/180)
 zcen[ii]=rcen[ii]*sin((file[ii].dec*!DPI)/180)
 cvec[ii]=[xcen[ii],ycen[ii],zcen[ii]]
 cenvel[ii]=file[ii].velocity

;***********************velocity pairs******************************
 for jj=0L,nelem-1 do begin
     if (jj ne ii) then begin
         veldiff[jj]=file[jj].velocity-cenvel[ii]
         if (abs(veldiff[jj]) le deltav) then begin
             
;**************** define radius/sphere around parent star to be searched******
             
             sphere[jj]=(x[jj]-xcen[ii])^2+(y[jj]-ycen[ii])^2+(z[jj]-zcen[ii])^2-(radius)^2.
             if sphere[jj] le 0 then begin
                 clump[ii+kk].parent=ii
                 clump[ii+kk].child=jj
                 clump[ii+kk].par_ra=file[ii].ra
                 clump[ii+kk].par_dec=file[ii].dec
                 clump[ii+kk].par_r=file[ii].r
                 clump[ii+kk].par_u=file[ii].u
                 clump[ii+kk].par_i=file[ii].i
                 clump[ii+kk].par_g=file[ii].g
                 clump[ii+kk].par_velocity=file[ii].velocity
                 clump[ii+kk].ch_ra=file[jj].ra
                 clump[ii+kk].ch_dec=file[jj].dec
                 clump[ii+kk].ch_g=file[jj].g
                 clump[ii+kk].ch_r=file[jj].r
                 clump[ii+kk].ch_u=file[jj].u
                 clump[ii+kk].ch_i=file[jj].i
                 clump[ii+kk].ch_velocity=file[jj].velocity
                 clump[ii+kk].veldiff=veldiff[jj]
                 clump[ii+kk].raddiff=sphere[jj]
                 kk=kk+1     
             endif
         endif
     endif
 endfor
endfor

;****************remove all empty entries
clumps=where(clump.par_g ne 0.00) 
clump=clump[clumps]
;****************save structure and file
mwrfits,clump,'/global/data/scr/adi/clumps.fits',/create
return,clump
end
