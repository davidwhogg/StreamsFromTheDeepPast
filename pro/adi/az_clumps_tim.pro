function az_clumps_tim,deltav=deltav,radius=radius
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
readcol,'/global/data/scr/adi/timBHB.dat',id, Teff,eTEff,logg,elogg,fe,efe,velocity,dis,edis,ra,dec
nelem=2490
kk=0
;********************define distance to stars using M=0.8, and x,y,z**********
rad=dis*1d-3  ;distance to stars in kpc

x=rad*cos((ra*!DPI)/180)*cos((dec*!DPI)/180)
y=rad*sin((ra*!DPI)/180)*cos((dec*!DPI)/180)
z=rad*sin((dec*!DPI)/180)
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
str1=create_struct('parent',0.,'child',0.,'par_ra',0.D,'par_dec',0.D,'par_velocity',0.D,'ch_ra',0.D,'ch_dec',0.D,'ch_velocity',0.D,'veldiff',0.,'raddiff',0.,'par_dist',0.D,'ch_dist',0.D)
clump=replicate(str1,100000)

;*************** pick star to be parent************
for ii=0L,nelem-1 do begin
 cen[ii]=ii
; rcen[ii]=10^((file[ii].g+4.3)/5.)
 rcen[ii]=rad[ii]
 xcen[ii]=rcen[ii]*cos((ra[ii]*!DPI)/180)*cos((dec[ii]*!DPI)/180)
 ycen[ii]=rcen[ii]*sin((ra[ii]*!DPI)/180)*cos((dec[ii]*!DPI)/180)
 zcen[ii]=rcen[ii]*sin((dec[ii]*!DPI)/180)
 cvec[ii]=[xcen[ii],ycen[ii],zcen[ii]]
 cenvel[ii]=velocity[ii]

;***********************velocity pairs******************************
 for jj=0L,nelem-1 do begin
     if (jj ne ii) then begin
         veldiff[jj]=velocity[jj]-cenvel[ii]
         if (abs(veldiff[jj]) le deltav) then begin
             
;**************** define radius/sphere around parent star to be searched******
             
             sphere[jj]=(x[jj]-xcen[ii])^2+(y[jj]-ycen[ii])^2+(z[jj]-zcen[ii])^2-(radius)^2.
             if sphere[jj] le 0 then begin
                 clump[ii+kk].parent=ii
                 clump[ii+kk].child=jj
                 clump[ii+kk].par_ra=ra[ii]
                 clump[ii+kk].par_dec=dec[ii]
     ;            clump[ii+kk].par_r=file[ii].r
      ;           clump[ii+kk].par_u=file[ii].u
       ;          clump[ii+kk].par_i=file[ii].i
        ;         clump[ii+kk].par_g=file[ii].g
                 clump[ii+kk].par_velocity=velocity[ii]
                 clump[ii+kk].par_dist=rad[ii]
                 clump[ii+kk].ch_ra=ra[jj]
                 clump[ii+kk].ch_dist=rad[jj]
                 clump[ii+kk].ch_dec=dec[jj]
         ;        clump[ii+kk].ch_g=file[jj].g
          ;       clump[ii+kk].ch_r=file[jj].r
           ;      clump[ii+kk].ch_u=file[jj].u
            ;     clump[ii+kk].ch_i=file[jj].i
                 clump[ii+kk].ch_velocity=velocity[jj]
                 clump[ii+kk].veldiff=veldiff[jj]
                 clump[ii+kk].raddiff=sphere[jj]
                 kk=kk+1     
             endif
         endif
     endif
 endfor
endfor

;****************remove all empty entries
clumps=where(clump.par_ra ne 0.00) 
clump=clump[clumps]
;****************save structure and file
mwrfits,clump,'/global/data/scr/adi/clumps.fits',/create
return,clump
end
