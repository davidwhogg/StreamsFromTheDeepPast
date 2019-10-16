pro gauss_mixtures_test

; setup
ngauss= 3
ndim= 3
ndata= 1000
seed= -1L
amp= double(ndata)*replicate(1d0,ngauss)/double(ngauss)
mean= double(randomn(seed,ndim,ngauss))
var= dblarr(ndim,ndim,ngauss)
for ii=0L,ndim-1 do begin
    var[ii,ii,*]= 1d-2
endfor
splog, 'amp=',amp
for ii=0L,ngauss-1 do splog, 'mean=',mean[*,ii]

; make fake data
data= fake_catalog(seed,amp,mean,var,gauss=gauss)
ndata= n_elements(data)/ndim
weight= dblarr(ndata)+1d0
zindx= where(gauss EQ 0)
weight[zindx]= weight[zindx]/2D0
tweight= total(weight,/double)

; fit
saveseed= seed
gauss_mixtures, ngauss,data,amp1,mean1,var1,$
  weight=2D0*weight,tol=(2d0*tweight/double(ndata))*1d-6,/quiet,seed=seed
splog, 'amp1=',amp1
for ii=0L,ngauss-1 do splog, 'mean1=',mean1[*,ii]

; fit again
seed= saveseed
gauss_mixtures, ngauss,data,amp2,mean2,var2,$
  weight=weight,tol=(tweight/double(ndata))*1d-6,/quiet,seed=seed
splog, 'amp2=',amp2
for ii=0L,ngauss-1 do splog, 'mean2=',mean2[*,ii]

end
