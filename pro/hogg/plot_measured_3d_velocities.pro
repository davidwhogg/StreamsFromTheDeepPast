;+
; NAME:
;   plot_measured_3d_velocities
; PURPOSE:
;   match up and plot fully measured velocities in Hipparcos data
; BUGS:
;   match algorithm is stoopid
; REVISION HISTORY:
;   2003-02-20  written - Hogg
;-
pro plot_measured_3d_velocities

; read in data
restore, '/global/data/hipparcos/hip-10-200-1.sav'
restore, '/global/data/hipparcos/hip-radv.sav'
nhip= n_elements(hip)

; match
for ii=0L,nhip-1 do begin
    indx= where(radv.hip EQ hip[ii].hip)
    if indx[0] GE 0 then begin
        if keyword_set(hipindx) then hipindx= [hipindx,ii] else $
          hipindx= ii
        if keyword_set(radvindx) then radvindx= [radvindx,indx[0]] else $
          radvindx= ii
    endif
endfor
nmatch= n_elements(hipindx)

; loop over matches, transforming to UVW
v3d= dblarr(3,nmatch)
for ii=0L,nmatch-1 do begin
    jj= hipindx[ii]
    kk= radvindx[ii]
    vtan= hip[jj].nsm # [hip[jj].vl,hip[jj].vb]
    ll1= hip[jj].l*!DPI/1.80D2
    bb1= hip[jj].b*!DPI/1.80D2
    vrad= radv[kk].vr * [cos(bb1)*cos(ll1),cos(bb1)*sin(ll1),sin(bb1)]
    v3d[*,ii]= vrad+vtan
endfor

; make plots
weight= fltarr(nmatch)+1.0
label= ['U','V','W']
vrange= [-70,70]
range= [[vrange],[vrange],[vrange]]
hogg_manyd_scatterplot,weight,v3d,'measured_3d_velocities.ps', $
  label=label,range=range

end
