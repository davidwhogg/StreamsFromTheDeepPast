;+
; NAME:
;   testhalo
; PURPOSE:
;   Test and validate ex-max codes
;-
pro testhalo, ngauss

; NEED TO CHECK PROJECTION MATRIX

; read SDSS galaxy colors
;  filename= '../data/emcolors.sample7b60.dat'
;  nfield= 7
;  template= {version:1.0, datastart: 3, $
;    delimiter: byte(32), missingvalue: 0.0, commentsymbol: '#', $
;    fieldcount: long(nfield), fieldtypes: lonarr(nfield)+5, $
;    fieldnames: ['id','u','g','r','i','z','M_r'], $
;    fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
;  data= read_ascii(filename,template=template)
;  seed= long(1)
;  ug= data.u-data.g
;  gr= data.g-data.r
;  ri= data.r-data.i
;  iz= data.i-data.z
;  Mr= data.M_r
;  ndata= n_elements(gr)
;  ddata= 4
;  if ddata EQ 2 then xx= transpose([[gr],[ri]])
;  if ddata EQ 3 then xx= transpose([[ug],[gr],[ri]])
;  if ddata EQ 4 then xx= transpose([[ug],[gr],[iz],[ri]])
;  if ddata EQ 5 then xx= transpose([[Mr],[ug],[gr],[iz],[ri]])

; make fake kinematic data
;  seed= 4L
;  ndata= 2000
;  fake_catalog, seed,ndata,100.0,120.0,0.5,10.0,qa='fake.ps', $
;    ra=ra,dec=dec,vr=vr,rr=rr,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz, $
;    dotra=radot,dotdec=decdot
;  dgauss= 3
;  yy= transpose([[vx],[vy],[vz]])

; read Hipparcos data
;  hip= read_hipparcos(rootdir='/home/users/mb144/hipparcos')
;  save, hip, filename='hip.sav'
  restore, 'hip.sav'
  parsec= 3.08568D13              ; km
  degree= double(!PI)/180D        ; radian
  arcsec= degree/3600D            ; radian
  year= 3.15569D7                 ; s
  rr= parsec*1000.0/hip.plx       ; km
; get into l,b
  radec_to_lb, hip.ra,hip.dec,hip.pmra,hip.pmdec,ll,bb,pmll,pmbb,1991.25
  ra= ll*degree                   ; radian
  dec= bb*degree                  ; radian
  radot= pmll*arcsec/1000.0/year  ; radian/s
  decdot= pmbb*arcsec/1000.0/year ; radian/s
  vra= rr*radot                   ; km/s
  vdec= rr*decdot                 ; km/s
  index= where(abs(vra) LT 450.0 AND abs(vdec) LT 450.0,ndata)
  xx= transpose([[vra[index]],[vdec[index]]])
  ra= ra[index]
  dec= dec[index]
  ddata= 2
  pos=fltarr(3,n_elements(ra))
  pos[0,*]=rr*cos(ra)*cos(dec)/parsec
  pos[1,*]=rr*sin(ra)*cos(dec)/parsec
  pos[2,*]=rr*sin(dec)/parsec
  openw,11,'testhalo.pts'
  writeu,11,pos
  close,11
  
help, ddata,ndata

; setup projection matrix
  dgauss= 3
  proj2= dblarr(dgauss,dgauss,ndata)
; this matrix takes radial and tangential velocities to x,y,z
  proj2[0,0,*]=  cos(dec)*cos(ra)
  proj2[1,0,*]= -sin(ra)
  proj2[2,0,*]= -sin(dec)*cos(ra)
  proj2[0,1,*]=  cos(dec)*sin(ra)
  proj2[1,1,*]=  cos(ra)
  proj2[2,1,*]= -sin(dec)*sin(ra)
  proj2[0,2,*]=  sin(dec)
  proj2[1,2,*]=  0
  proj2[2,2,*]=  cos(dec)
  proj= proj2
  for i=0,ndata-1 do proj[*,*,i]= invert(proj2[*,*,i],/double)
  if ddata EQ 2 then proj= proj[*,1:2,*]

; setup first guess at gaussians
  if (not keyword_set(ngauss)) then ngauss= 5
	print,ngauss
  ;restore, 'testhalo.g'+strtrim(string(ngauss),2)+'.sav'
  amp= dblarr(ngauss)+1.0
  mean= dblarr(dgauss,ngauss)
  var= dblarr(dgauss,dgauss,ngauss)
  yy= dblarr(dgauss,ndata)
  for i=0,ndata-1 do yy[*,i]= proj[*,*,i]#xx[*,i]
  for d=0,dgauss-1 do begin
    mmm= djs_mean(yy[d,*])
    vvv= djsig(yy[d,*])
    mean[d,*]= mmm+vvv*(2.0*randomu(seed,ngauss)-1.0)
    var[d,d,*]= (vvv*(1.0*randomu(seed,ngauss)+0.5))^2
  endfor

; project to observed space
;  xx= dblarr(ddata,ndata)
;  for i=0,ndata-1 do begin
;    xx[*,i]= transpose(yy[*,i]#proj[*,*,i])
;  endfor

; run E-M
  qa='halo2.g'+strtrim(string(ngauss),2)+'.ps'
  ;if ddata EQ 2 then qa='halo2.ps'
  ex_max_proj, xx,proj,amp,mean,var,weight=weight,qa=qa,maxiterate=1000

; save results before exiting
  save, index, xx,proj,amp,mean,var,ddata,dgauss,ndata,ngauss,weight, $
    filename='testhalo.g'+strtrim(string(ngauss),2)+'.sav'
end

