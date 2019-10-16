pro smosaic_pal5
ra0= 229.022                    ; deg J2000
dec0= -0.111                    ; deg J2000
ra_size= 3.0/60.0
dec_size= 3.0/60.0
rerun= lindgen(1000)
pixscale= 0.4/3600.0
prefix= 'pal5'
title= 'Palomar 5'
all= 1
smosaic_make_jpg, ra0,dec0,ra_size,dec_size, $
 rerun=rerun,prefix=prefix,pixscale=pixscale,smfwhm=pixscale*2.0*3600.0, $
 minscore=0.0,all=all,nskyframe=5,title=title,npixround=8, $
 quality=90,/dontdelete,/ivarout,/nogrid
return
end
