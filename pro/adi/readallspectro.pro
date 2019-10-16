pro readallspectro
;+
;By :AZ
;
;Read in all spectro files and save only class=stars
;
;June 28,06
;+


nrows=1d5
nchunks=12
ii=0L
kk=0L

for jj=0,nchunks-1 do begin
    files=mrdfits('/global/data/sdss/spectro/spAll.fits',1,rows=(lindgen(nrows)+jj*nrows))
    stars=where (files.class eq 'STAR  ',nstars)
    if (nstars ge 1) then begin
        files=files[stars]
    endif
    nelem=n_elements(files)

    if (jj eq 0) then begin
        
        everything=files


      endif

      kk=ii+nelem
      if (jj ge 1) then begin



          everything=[everything,files]
      endif

      ii=kk
      
      
  endfor
mwrfits2,everything,'/global/data/scr/adi/spectro.fits',/create    
end
