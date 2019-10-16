pro makecuts
;+
;
;
;
;
;-

file=mrdfits('/global/data/scr/adi/newdata.fits',1)
nelem=n_elements(file)
cut=where((file.ra gt 124) and $  ;124
          (file.ra lt 170) and $  ;170
          (file.dec gt 0) and $
          (file.dec lt 65),ncut) ;65
if (ncut ge 1) then begin
    file=file[cut]
    nc=n_elements(file)
    blue=0.2+((0.1-0.2)*((file.gmag-20.4)/2.0))
    red=0.6+((0.4-0.6)*((file.gmag-20.4)/2.0))
    cutst=create_struct('ra',0.D,'dec',0.D,'gmag',0.,'gmi',0.)
    new=replicate(cutst,nc)
    new.ra=file.ra
    new.dec=file.dec
    new.gmag=file.gmag
    new.gmi=file.gmag-file.imag
    magcut=where((new.gmi gt blue)$
                 and (new.gmi lt red) $
                 and (new.gmag lt 22.5) $
                 and (new.gmag gt 20.2) $
                 ,nmag)  
    if nmag ge 1 then begin
        new=new[magcut]

                
    endif 
endif
mwrfits,new,'/global/data/scr/adi/grill_cuts.fits',/create
end
