;+
; NAME:
;   testem
; PURPOSE:
;   Test and validate ex-max codes
;-
pro testsimple

ndata=1000l
ddata=3l
ngauss=1l
dgauss=3l

seed=99l
xx=dblarr(ddata,ndata)
xx[0,*]=10.+randomn(seed,ndata)
xx[1,*]=10.+5.*randomn(seed,ndata)
xx[2,*]=10.+randomn(seed,ndata)

proj2= dblarr(dgauss,dgauss,ndata)
proj2[0,0,*]=  1.
proj2[1,0,*]=  0.
proj2[2,0,*]=  0.
proj2[0,1,*]=  0.
proj2[1,1,*]=  1.
proj2[2,1,*]=  0.
proj2[0,2,*]=  0.
proj2[1,2,*]=  0
proj2[2,2,*]=  1.
proj=proj2
for i=0,ndata-1 do proj[*,*,i]= invert(proj2[*,*,i],/double)

;  restore, 'testsimple.sav'
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

; run E-M
  qa='testsimple.ps'
  ;if ddata EQ 2 then qa='em2.ps'
  ex_max_proj, xx,proj,amp,mean,var,weight=weight,qa=qa,miniterate=100, $
    maxiterate=200

; save results before exiting
  save, xx,proj,amp,mean,var,ddata,dgauss,ndata,ngauss,weight, $
    filename='testsimple.sav'
end

