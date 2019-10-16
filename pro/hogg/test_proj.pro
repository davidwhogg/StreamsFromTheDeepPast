pro test_proj

restore,'hip-10-200-1.sav'

all=lindgen(n_elements(hip))
firsthalf=lindgen(n_elements(hip)/2)
secondhalf=n_elements(hip)/2+lindgen(n_elements(hip)-n_elements(hip)/2)

for iexp=0, 5 do begin 
    ngauss=2^iexp
    indx=all
    projection=hip[indx].nsm
    ydata=dblarr(2,n_elements(indx))
    ydata[0,*]=hip[indx].vl
    ydata[1,*]=hip[indx].vb
    ycovar=hip[indx].vlvbc
; fix to use arg_present
    xguess=dblarr(3,n_elements(indx),ngauss)
    xcovarguess=dblarr(3,3,n_elements(indx),ngauss)
    projected_gauss_mixtures, ngauss, ydata, ycovar, projection, $
      amp, xmean, xcovar, avgloglikedata=avgloglikedata,quiet=quiet, $
      /debug,/plot3d,rotation=hip[indx].sm,xguess=xguess, $
      xcovarguess=xcovarguess
    save,filename='gauss'+strtrim(string(ngauss),2)+'-10-200-1.sav'
endfor
stop

indx=secondhalf
projection=hip[indx].nsm
ydata=dblarr(2,n_elements(indx))
ydata[0,*]=hip[indx].vl
ydata[1,*]=hip[indx].vb
ycovar=hip[indx].vlvbc
projected_gauss_analyze, ydata, ycovar, projection, $
  amp, xmean, xcovar, avgloglikedata=avgloglikedata,quiet=quiet

covar=[[1.,0.3,0.5], $
       [0.3,2.,0.9], $
       [0.5,0.9,30.]]
mean=dblarr(3)+50.5
nsig=20.
threemesh, mean, covar, nsig, linepos, lineindx, nptsring=nptsring, $
  nrings=nrings
openw,unit,'test.pts',/get_lun
writeu,unit,linepos
free_lun,unit
openw,unit,'test.indx',/get_lun
writeu,unit,lineindx
free_lun,unit

end
