;+
;
;-
pro test_blob_mixtures_oned

ndata= 100000
mm= 0.5
SS= dblarr(ndata)+1.0
VVtrue= 0.01
TT= SS+VVtrue

ww= dblarr(ndata)+1
yy= mm+randomn(seed,ndata)*sqrt(TT)

blob_mixtures_oned, ww,yy,SS

end

