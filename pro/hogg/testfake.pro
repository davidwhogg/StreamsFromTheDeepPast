;+
; NAME:
;   testfake
; PURPOSE:
;   Test and validate ex-max code on a simple galaxy model
;-
pro testfake, ngauss

; restore hipparcos data
indir= '/global/data/hipparcos'
suffix= '-10-200-1.sav'
restore, '/usr/local/data/dwh2/hip.sav'
nhip= n_elements(hip)
splog, 'nhip = ',nhip

; set random seed
seed= -1L

; loop over models to make
nmodel= 1
for model= 0L,nmodel-1 do begin
    modelstr= strtrim(string(model),2)

; set thin disk, thick disk, and halo amplitudes, means, variances
    if model EQ 0 then begin
        amp= [1.0,0.3,0.01]
        mean= dblarr(3,3)
        mean[*,2]= [0,-220,0]
        var= dblarr(3,3,3)
        var[*,*,0]= [[20^2,0,0],[0,10^2,0],[0,0,10^2]]
        var[*,*,1]= [[30^2,0,0],[0,20^2,0],[0,0,15^2]]
        var[*,*,2]= [[160^2,0,0],[0,100^2,0],[0,0,120^2]]
    endif

; make fake catalog
    v3d= fake_catalog(seed,amp,mean,var,ndata=nhip,gauss=gauss)

; loop over points projecting to vl,vb
    for ii=0L,nhip-1 do begin
        v2d= hip[ii].nsm ## v3d[*,ii]
        hip[ii].vl= v2d[0]
        hip[ii].vb= v2d[1]

; add gaussian noise
        noise2d= fake_catalog(seed,1.0,[0.0,0.0],hip[ii].vlvbc,/quiet)
        hip[ii].vl= hip[ii].vl+noise2d[0]
        hip[ii].vb= hip[ii].vb+noise2d[1]
    endfor

; make QA plots
    weight= fltarr(nhip)+1.0
    datablock= transpose([[transpose(v3d)],[hip.vl],[hip.vb]])
    label= ['U','V','W','v_l','v_b']
    vrange= [-70,70]
    range= [[vrange],[vrange],[vrange],[vrange],[vrange]]
    hogg_manyd_scatterplot, weight,datablock,'testfake-'+modelstr+'.ps', $
      label=label,range=range

; make savesets
    smalli= (shuffle(seed,hip))[0:999]
    save, hip,v3d,gauss,smalli, $
      filename=indir+'/fake-'+modelstr+suffix

endfor
end
