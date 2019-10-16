pro vector
;+
;
;AZ
;
;
;-

file=mrdfits('/global/data/scr/adi/spectro.fits',1)
nelem=n_elements(file)
jj=0L
kk=0L




print, nelem
x=cos((file.ra*!DPI)/180)*cos((file.dec*!DPI)/180)
y=sin((file.ra*!DPI)/180)*cos((file.dec*!DPI)/180)
z=sin((file.dec*!DPI)/180)

;;ra and dec of field Belokurov A7

x1=cos((185.0*!DPI)/180)*cos((11.2*!DPI)/180)
y1=sin((185.0*!DPI)/180)*cos((11.2*!DPI)/180)
z1=sin((11.2*!DPI)/180)

;; ra and dec of Belokurov A16
x2=cos((140.0*!DPI)/180)*cos((19.4*!DPI)/180)
y2=sin((140.0*!DPI)/180)*cos((19.4*!DPI)/180)
z2=sin((19.4*!DPI)/180)



v1=[x1,y1,z1]
v2=[x2,y2,z2]
result=crossp(v1,v2)
norm=sqrt((result[0]*result[0])+(result[1]*result[1])+(result[2]*result[2]))
unit=result/norm

rr=atan(unit[1],unit[0])
rr_deg=(rr*180)/!dpi



p_dec=asin(unit[2])
poledec=(p_dec*180)/!dpi



if (unit[0] lt 0) then begin

polera=rr_deg+180
endif

if ((unit[0] ge 0) and (unit[1] lt 0)) then begin
    polera=rr_deg+360
endif


print, polera,poledec 

for ii=0L,nelem-1 do begin
    adist=djs_diff_angle(polera,poledec,file[ii].ra,file[ii].dec,units=degrees)
    
    
    cuts=where ((adist ge 87) and (adist lt 93),ncuts)
    if (ncuts ge 1) then begin
        
        poles=file[cuts]
        st=create_struct('g',0.,'i',0.,'u',0.,'r',0.,'ra',0.D,'dec',0.D,'velocity',0.D)
        
        mags=replicate(st,ncuts)
        mags.ra=poles.ra
        mags.dec=poles.dec
        mags.g=22.5-(2.5*alog10(poles.psfflux[1]))-poles.extinction[1]
        mags.r=22.5-(2.5*alog10(poles.psfflux[2]))-poles.extinction[2]
        mags.u=22.5-(2.5*alog10(poles.psfflux[0]))-poles.extinction[0]
        mags.i=22.5-(2.5*alog10(poles.psfflux[3]))-poles.extinction[3]
        mags.velocity=poles.z*(3d5)
        
        if (not keyword_set(everything)) then begin
            everything=mags
            everything=replicate(mags[0],3000000L)
            nevery=n_elements(everything)
        endif
        kk=jj+ncuts
        if (kk gt nevery)then begin
            everything=[everything, replicate(mags[0],10000000L)]
            nevery=n_elements(everything)
        endif
        everything[jj:kk-1]=mags
        jj=kk
    endif

    
    
    
    

endfor  


mwrfits,everything[0:jj-1],'/global/data/scr/adi/polecuts_spectro_data.fits'







end
