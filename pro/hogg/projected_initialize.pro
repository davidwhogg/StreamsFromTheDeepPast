;+
; NAME:
;   projected_initialize
; PURPOSE:
;   initialize means for projected_kmeans or projected_gauss_mixtures
; INPUTS:
;   num - number of groups
;   y - [2,N] data points 
; OUTPUTS:
; BUGS:
;   weights not used
;   not actually writtem
;   initialization not general
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro projected_initialize, ngroup, ydata, ycovar, projection, amp, $
                          xmean, xcovar, quiet=quiet

ndimy=n_elements(ycovar)/n_elements(ydata)
ndata=n_elements(ydata)/ndimy
ndimx=n_elements(projection)/ndata/ndimy

amp=dblarr(ngroup)+1.D/double(ngroup)

xdata=dblarr(ndimx,ndata)
for ii=0L, ndata-1L do $
  xdata[*,ii]=projection[*,*,ii]#ydata[*,ii]

xcovar=dblarr(ndimx,ndimx,ngroup)
for kk=0L, ndimx-1L do begin
    xcovar[kk,kk,*]=variance(xdata[kk,*])/double(ngroup)
endfor

randata=lonarr(ngroup)
while(n_elements(uniq(randata,sort(randata))) lt ngroup) do $
  randata=long(double(ndata)*randomu(seed,ngroup))
xmean=xdata[*,randata]

;xmean=20.*randomn(seed,ndimx,ngroup)
;xcovar=dblarr(ndimx,ndimx,ngroup)
;for jj=0L, ngroup-1L do $
;  xcovar[*,*,jj]=(10.)^2*identity(ndimx)

end
