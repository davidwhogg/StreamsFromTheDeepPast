pro makeimage
;+
;
;
;
;-


file=mrdfits('grillmaircut.fits',1);starstreams.fits',1)
nelem=n_elements(file)
;blah=where((file.ra gt 124 and file.ra lt 170) and (file.dec gt 0 and file.dec lt 65) and (file.gmag gt 17 and file.gmag lt 24))
;new=create_struct('ra',0.D,'dec',0.D,'gmag',0.,'imag',0.)
;struct=replicate(new,nelem)
;file=file[blah]
;struct.ra= where (file.ra gt 124. and file.ra lt 170.) 
;struct.dec= where (file.dec gt 0. and file.dec lt 65.)
;populate_image,

set_plot,'ps'
device,filename='plot.ps'
hogg_scatterplot, file.ra,file.dec,xrange=[170,130],xnpix=20,ynpix=20,grid=grid
device,/close

atv,grid
;mwrfits,file,'grillmaircut.fits',/create
end
