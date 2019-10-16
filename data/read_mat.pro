pro read_mat

d2r=3.14159265358979D/180.D
data=load_mat_data('blantonv4.mat', $
                   ['centres','assignments','vl','vb','ll','bb'], $
                   centres,assignments,vl,vb,ll,bb)

; make points
seed=-10l
nsample=10
llout=d2r*replicate(1.,nsample)#ll
bbout=d2r*replicate(1.,nsample)#bb
vlout=replicate(1.,nsample)#vl
vbout=replicate(1.,nsample)#vb
vrout=(500.+10.*randomn(seed,nsample*n_elements(ll)))*((replicate(1.,nsample*n_elements(ll))))
;vrout=(((2.*dindgen(nsample)/double(nsample))-1.)#replicate(1.,n_elements(ll)))*((replicate(1.,nsample))#vr)
vout=fltarr(3,n_elements(vlout))
vout[0,*]=-sin(llout)*vlout-cos(llout)*sin(bbout)*vbout+ $
  cos(llout)*cos(bbout)*vrout
vout[1,*]= cos(llout)*vlout-sin(llout)*sin(bbout)*vbout+ $
  sin(llout)*cos(bbout)*vrout
vout[2,*]= cos(bbout)*vbout+sin(bbout)*vrout
openw,unit,'vline.pts',/get_lun
writeu,unit,vout
close,unit
free_lun,unit
openw,11,'pick.dat'
out=fltarr(5,n_elements(vlout))
out[0,*]=vlout
out[1,*]=vbout
out[2,*]=vout[0,*]
out[3,*]=vout[1,*]
out[4,*]=vout[2,*]
printf,11,format='(f14.8,x,f14.8,x,f14.8,x,f14.8,x,f14.8)',out
close,11

; make vertices of lines
;llout=d2r*[ll,ll]
;bbout=d2r*[bb,bb]
;vlout=[vl,vl]
;vbout=[vb,vb]
;vrout=[vr,-vr]
;vout=fltarr(3,n_elements(vlout))
;vout[0,*]=-sin(llout)*vlout-cos(llout)*sin(bbout)*vbout+ $
;  cos(llout)*cos(bbout)*vrout
;vout[1,*]= cos(llout)*vlout-sin(llout)*sin(bbout)*vbout+ $
;  sin(llout)*cos(bbout)*vrout
;vout[2,*]= cos(bbout)*vbout+sin(bbout)*vrout
;openw,unit,'vline.vtx',/get_lun
;writeu,unit,vout
;close,unit
;free_lun,unit
;
;; make indices 
;iout=lonarr(3,n_elements(ll))
;iout[0,*]=2
;iout[1,*]=lindgen(n_elements(ll))
;iout[2,*]=n_elements(ll)+lindgen(n_elements(ll))
;openw,unit,'vline.indx',/get_lun
;writeu,unit,iout
;close,unit
;free_lun,unit

end

