;
;
;
;
;
;
pro hogg_make_poly, xin,yin,xout,yout
sindx= sort(xin)
xx= xin[sindx]
yy= yin[sindx]
dx= (xx-shift(xx,1)) > (shift(xx,-1)-xx)
npt= 2*n_elements(yy)
xx= reform(transpose([[xx-0.5*dx],[xx+0.5*dx]]),npt)
yy= reform(transpose([[yy],[yy]]),npt)
xout= [!X.CRANGE[0],!X.CRANGE[0],xx,!X.CRANGE[1],!X.CRANGE[1]]
yout= [!Y.CRANGE[0],yy[0],yy,yy[npt-1],!Y.CRANGE[0]]
return
end

pro az_kin_comp
mw1=mrdfits('/global/data/scr/adi2/obsmw1.fits',1)
tunit=13.7/.3331
mw1.accrtime3=mw1.accrtime3*tunit
mw1.accrtime1=mw1.accrtime1*tunit

nbins=26
ntotal=dblarr(1,nbins)
in_z=dblarr(1,nbins)
avg_z=dblarr(1,nbins)
out_z=dblarr(1,nbins)
indisk=dblarr(1,nbins)
inhalo=dblarr(1,nbins)
inbulge=dblarr(1,nbins)
inrest=dblarr(1,nbins)
shalo=dblarr(1,nbins)
srest=dblarr(1,nbins)
sdisk=dblarr(1,nbins)

xsize=6.5
ysize=6.5

for ii=1,nbins-1 do begin
out_z[ii]=ii*0.5
in_z[ii]=(ii*0.5)-0.5
avg_z[ii]=(in_z[ii]+out_z[ii])/2.

all=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii])) 
ntotal[ii]=n_elements(all)

;disk
idisk=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
            and (mw1.comp eq 1)) 
indisk[ii]=n_elements(idisk)
plotsym,0,1.,/fill
smdisk=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
             and (mw1.comp eq 1) and (mw1.maxhalo eq 0.))
sdisk[ii]=n_elements(smdisk)

;halo
ihalo=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
            and (mw1.comp eq 2)) 
inhalo[ii]=n_elements(ihalo)

smhalo=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
             and (mw1.comp eq 2) and (mw1.maxhalo eq 0.))
;bulge
ibulge=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
             and (mw1.comp eq 3)) 
inbulge[ii]=n_elements(ibulge)

;rest
irest=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
            and (mw1.comp eq 4)) 
inrest[ii]=n_elements(irest)
smrest=where((abs(mw1.z) gt in_z[ii]) and (abs(mw1.z) lt out_z[ii]) $
             and (mw1.comp eq 4) and (mw1.maxhalo eq 0.))
srest[ii]=n_elements(smrest)

shalo[ii]=n_elements(smhalo)+n_elements(idisk)+n_elements(irest)
endfor

;*****************************PLOT*********************************
set_plot,'ps'
device,file='brown_fraction_try.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/encapsulated,/color,bits=8
device,/times,/isolatin1
!X.THICK= 2
!Y.THICK= 2
plot,[0],[0],/nodata,xrange=[0,10],yrange=[0,1],ytitle='fraction',charsize=1.2,xtitle='z [kpc]'
;halo
lowx=0.25
highx=9.75
hal=shalo/ntotal
lowy=min(shalo/ntotal)
highy=shalo[1]/ntotal[1]
indices=value_locate(avg_z,[lowx,highx])
lowindex=indices[0]
highindex=indices[1]
if avg_z(lowindex) lt lowx then lowindex=lowindex+1
if avg_z(highindex) gt highx then highindex=highindex-1
hogg_make_poly, avg_z[lowindex:highindex],hal[lowindex:highindex],xpoly,ypoly
polyfill,!X.CRANGE[[0,0,1,1]],!Y.CRANGE[[0,1,1,0]],col=0.53*!d.n_colors
polyfill,xpoly,ypoly,col=0.5*!d.n_colors
;linestyle=2,spacing=3.,/line_fill;
;whiteout
rest=(indisk+inrest)/ntotal
lowrest=min(rest)
highrest=rest[1]
hogg_make_poly, avg_z[lowindex:highindex],rest[lowindex:highindex], $
  xpoly,yrestpoly
polyfill,xpoly,yrestpoly,col=0.73*!d.n_colors
;rest
res=(indisk+srest)/ntotal
lowres=min(res)
highres=res[1]
hogg_make_poly, avg_z[lowindex:highindex],res[lowindex:highindex], $
  xpoly,yrpoly
polyfill,xpoly,yrpoly,col=0.7*!d.n_colors
;linestyle=3,spacing=2.,/line_fill;
;whiteout
disk=(indisk)/ntotal
lowdisk=min(disk)
highdisk=disk[1]
hogg_make_poly, avg_z[lowindex:highindex],disk[lowindex:highindex], $
  xpoly,ydiskpoly
polyfill,xpoly,ydiskpoly,col=0.93*!d.n_colors
;disk
di=(sdisk)/ntotal
lowdi=min(di)
highdi=di[1]
hogg_make_poly, avg_z[lowindex:highindex],di[lowindex:highindex], $
  xpoly,ydipoly
polyfill,xpoly,ydipoly,col=0.9*!d.n_colors
;linestyle=1,spacing=0.5,/line_fill;

okay= lindgen(n_elements(xpoly)-4)+2
plots,xpoly[okay],ydiskpoly[okay],linestyle=0,thick=10.
plots,xpoly[okay],ydipoly[okay],linestyle=0,thick=2
plots,xpoly[okay],ypoly[okay],linestyle=0,thick=2
plots,xpoly[okay],yrestpoly[okay],linestyle=0,thick=10
plots,xpoly[okay],yrpoly[okay],linestyle=0,thick=2

ylabel= 0.3
xyouts,1.5,ylabel,'disk',charthick=4.,charsize=2.,align=0.5
xyouts,4.0,ylabel,'rest',charthick=4.,charsize=2.,align=0.5
xyouts,7.0,ylabel,'halo',charthick=4.,charsize=2.,align=0.5

device,/close

end
