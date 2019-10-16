;+
; NAME:
;   check_predicted_v3d
; PURPOSE:
;   plot predicted and measured velocities
; BUGS:
;   Not totally sure I am transforming the covariance matrix correctly - Hogg.
;-
pro plot_one_star_v3d, ampp,vp,varp,vo,varo,hipno=hipno
rrange= [-1,1]*40.0+vo
frange= [-0.1,1.1]*1.0/sqrt(2.0*!DPI)
plot, [0],[0],/nodata, $
  xrange=rrange,xstyle=9,xtickinterval=40.0,xticklen=-0.1, $
  yrange=frange,ystyle=5
tweak= 0.8
oplot_1d_gaussian, ampp*tweak*sqrt(total(ampp*varp)/total(ampp)),vp,varp, $
  thick=1,noclip=1
if varo GE 0.0 then begin
    oplot_1d_gaussian, sqrt(varo),vo,varo, $
      noclip=1
endif
if keyword_set(hipno) then begin
    xyouts, !X.CRANGE[0],!Y.CRANGE[0],strtrim(string(hipno),2), $
      charsize=0.5,orientation=90
endif
end

pro test_radial_velocities, ngauss

; restore data
ngstr= strtrim(string(round(ngauss)),2)
restore, '/global/data/hipparcos/gauss'+ngstr+'-10-200-1.sav'
amp= amp
mean= xguess
covar= xcovarguess
restore, '/global/data/hipparcos/hip-radv.sav'

; initialize stuff
nhip= n_elements(hip)
ngauss= n_elements(amp)
vrlb= dblarr(3,ngauss)
vrlbc= dblarr(3,3,ngauss)

; make hip - radv match
radvindx= lonarr(nhip)
for ii=0L,nhip-1 do begin
    radvindx[ii]= (where(radv.hip EQ hip[ii].hip))[0]
endfor
matchindx= radvindx[where(radvindx GE 0)]

; setup plotting
set_plot, "PS"
psfilename= 'gauss'+ngstr+'_radial_velocities.ps'
xsize= 8.0
ysize= 11.0
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
hogg_plot_defaults
!X.MARGIN= [1,0]
!Y.MARGIN= [2.5,0]

; plot HR diagrams
!P.MULTI=[0,2,1]
distance= 1D3/hip.plx  ; pc
absmag= hip.vmag-5D0*alog10(distance/1D1)
magrange= [-2,15]
colorrange= [-0.2,1.7]
magtitle='M!dV!n  (mag)'
colortitle='(B-V)  (mag)'
;plot, absmag,hip.bvcolor, $
;  xtitle=magtitle,ytitle=colortitle, $
;  xrange=magrange,yrange=colorrange,psym=3
;plot, absmag[matchindx],hip[matchindx].bvcolor, $
;  xtitle=magtitle,ytitle=colortitle, $
;  xrange=magrange,yrange=colorrange,psym=3;

; shuffle hip and start loop
seed= -1L
indx= shuffle(seed,hip)
nplot= 600
!P.MULTI= [0,15,25]
!Y.CHARSIZE= 0.001
for ii=0L,nplot-1 do begin
    hip1= hip[indx[ii]]

; do we have a radial velocity?
    jj= radvindx[indx[ii]]
    if jj GE 0 then begin
        vr= radv[jj].vr
        vr_var= (radv[jj].vr_err)^2
    endif else begin
        vr= 0.0
        vr_var= -1.0
    endelse

; project mean and covar to r,l,b space
    for kk=0L,ngauss-1 do begin
        vrlb[*,kk]= transpose(hip1.sm) # mean[*,indx[ii],kk]
        vrlbc[*,*,kk]= transpose(hip1.sm) # covar[*,*,indx[ii],kk] # hip1.sm
    endfor

; make plots
    plot_one_star_v3d,amp,vrlb[0,*],vrlbc[0,0,*],vr,vr_var,hipno=hip1.hip
    plot_one_star_v3d,amp,vrlb[1,*],vrlbc[1,1,*],hip1.vl,hip1.vlvbc[0,0]
    plot_one_star_v3d,amp,vrlb[2,*],vrlbc[2,2,*],hip1.vb,hip1.vlvbc[1,1]
    if (ii MOD 4) NE 3 then plot,[0],[0],/nodata,xstyle=4,ystyle=4

endfor
device,/close
end
