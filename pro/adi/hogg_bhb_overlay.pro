pro hogg_bhb_overlay
naxis1= 1000
naxis2= 1000
data= hogg_image_overlay('/home/users/az481/streams/pro/adi/brancha_sag.ps',$
                         naxis1,naxis2)
bhbmodel= hogg_image_overlay('/home/users/az481/streams/pro/adi/models_bhb.ps',$
                          naxis1,naxis2)
bmsmodel= hogg_image_overlay('/home/users/az481/streams/pro/adi/models_bs.ps',$
                          naxis1,naxis2)
; yellow
bhbmodel[*,*,0]= 1.0
bhbmodel[*,*,1]= 1.0

; cyan
bmsmodel[*,*,1]= 1.0
bmsmodel[*,*,2]= 1.0

data= data < bhbmodel
data= data < bmsmodel
data= nw_float_to_byte(data)
WRITE_JPEG,'hogg_bhb_overlay.jpg',temporary(data),TRUE=3,QUALITY=95

return
end
