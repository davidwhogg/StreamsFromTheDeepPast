num1=1000
num2=600
var=0.01
xy1=[[3.+0.1*randomn(seed,num1)],[3.+2.5*randomn(seed,num1)],[3.+2.5*randomn(seed,num1)]]
xy2=[[3.+2.5*randomn(seed,num2)],[3.+0.1*randomn(seed,num2)],[3.+2.5*randomn(seed,num2)]]
xy=dblarr(3,num1+num2)
xy[0,*]=[xy1[*,0],xy2[*,0]]
xy[1,*]=[xy1[*,1],xy2[*,1]]
xy[2,*]=[xy1[*,2],xy2[*,2]]
covar=dblarr(3,3,num1+num2)
covar[0,0,*]=var
covar[1,1,*]=var
covar[2,2,*]=var

ngroup=12
kmeans,ngroup,xy,group
plot,xy[0,*],xy[1,*],psym=3
for i=0, ngroup-1 do begin & $
    indx=where(group eq i) & $
    oplot,xy[0,indx],xy[1,indx],psym=3+i & $
endfor

gauss_mixtures, ngroup, xy, gauss_amp, gauss_mean, gauss_covar, $
  loglikedata=loglikedata, peg_dims=peg_dims
