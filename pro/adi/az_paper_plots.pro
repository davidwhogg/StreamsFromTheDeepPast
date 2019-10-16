pro az_paper_plots
;
;
;
;
;
;
mw1=mrdfits('/global/data/scr/adi2/gal1.fits',1)
nelem1=n_elements(mw1)
tunit=13.7/.3331
mw1.accrtime3=mw1.accrtime3*tunit
mw1.tform=mw1.tform*tunit

rvir=383.-150.

halo=where((mw1.haloid eq 1) and (abs(mw1.z) gt 7))
mw1=mw1[halo]
nelem=n_elements(mw1)

xsize=7.5
ysize=7.5
;********
nbins=15
ntotal=dblarr(1,nbins)
out_dist=dblarr(1,nbins)
in_dist=dblarr(1,nbins)
avg_dist=dblarr(1,nbins)
nacc1=dblarr(1,nbins)
nacc2=dblarr(1,nbins)
nacc3=dblarr(1,nbins)
nacc=dblarr(1,nbins)
mass1=dblarr(1,nbins)
mass2=dblarr(1,nbins)
mass3=dblarr(1,nbins)
mass3=dblarr(1,nbins)
pmass1=dblarr(1,nbins)
pmass2=dblarr(1,nbins)
time1=dblarr(1,nbins)
time2=dblarr(1,nbins)

set_plot,'ps'

device,file='radial_mass_gal1.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/encapsulated,/color,bits=8
device,/times,/isolatin1
erase & multiplot,[1,2]
plot,[0],[0],/nodata,xrange=[0,150],yrange=[0,1],ytitle='accreted fraction',charsize=1.2,title='gal1',xticks=8,xtickn=['','0.04r','0.09r','0.14r','0.2r','0.25r','0.3r','0.35r','0.4r']
for ii=1,nbins-1 do begin

out_dist[ii]=10*(ii)
in_dist[ii]=((10*ii)-10)
avg_dist[ii]=(out_dist[ii]+in_dist[ii])/2.

alltotal=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
ntotal[ii]=n_elements(alltotal)

;******
acc=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.delay le 0))
nacc[ii]=n_elements(acc)

acc1=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.accrtime3 lt 4.)and (mw1.delay le 0))    
nacc1[ii]=n_elements(acc1)
m1=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.accrtime3 lt 4.)and (mw1.delay le 0) and (mw1.maxhalo ne 0.)) 
if (n_elements(m1) gt 1) then mass1[ii]=alog10(avg(mw1[m1].maxhalo))

;*********
acc2=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.delay le 0) and (mw1.accrtime3 ge 4.0) $
           and (mw1.accrtime3 lt 7.))
nacc2[ii]=n_elements(acc2)
m2=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.accrtime3 lt 7.)and (mw1.delay le 0) and $
         (mw1.accrtime3 ge 4) and (mw1.maxhalo ne 0.)) 
if (n_elements(m2) gt 1) then mass2[ii]=alog10(avg(mw1[m2].maxhalo))

;*******
acc3=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.accrtime3 ge 7.) and (mw1.delay le 0))
nacc3[ii]=n_elements(acc3)

m3=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and $
           (mw1.accrtime3 ge 7.)and (mw1.delay le 0) and (mw1.maxhalo ne 0.)) 
if (n_elements(m3) gt 1) then mass3[ii]=alog10(avg(mw1[m3].maxhalo))


if (mass1[ii] lt 10.3) then begin
plotsym,8,0.6,/fill ;square
plots,avg_dist[ii],(nacc1[ii]/ntotal[ii]),psym=8,color=djs_icolor('cyan')
endif
if ((mass1[ii] ge 10.3) and (mass1[ii] lt 10.8))then begin
plotsym,8,0.9,/fill ;square
plots,avg_dist[ii],(nacc1[ii]/ntotal[ii]),psym=8,color=djs_icolor('cyan')
endif
if ((mass1[ii] ge 10.8))then begin
plotsym,8,1.3,/fill ;square
plots,avg_dist[ii],(nacc1[ii]/ntotal[ii]),psym=8,color=djs_icolor('cyan')
endif

if (mass2[ii] lt (10.3)) then begin
plotsym,4,0.6,/fill;traiangle
plots,avg_dist[ii],(nacc2[ii]/ntotal[ii]),psym=8,color=djs_icolor('magenta')
endif
if ((mass2[ii] ge 10.3) and (mass2[ii] lt 10.8))then begin
plotsym,4,0.9,/fill;traiangle
plots,avg_dist[ii],(nacc2[ii]/ntotal[ii]),psym=8,color=djs_icolor('magenta')
endif
if (mass2[ii] ge (10.8)) then begin
plotsym,4,1.3,/fill;traiangle
plots,avg_dist[ii],(nacc2[ii]/ntotal[ii]),psym=8,color=djs_icolor('magenta')
endif

if (mass3[ii] lt 10.3) then begin
plotsym,3,0.6,/fill ;star
plots,avg_dist[ii],(nacc3[ii]/ntotal[ii]),psym=8,color=djs_icolor('green')
endif
if ((mass3[ii] ge 10.3) and (mass3[ii] lt 10.8))then begin
plotsym,3,0.9,/fill ;star
plots,avg_dist[ii],(nacc3[ii]/ntotal[ii]),psym=8,color=djs_icolor('green')
endif
if (mass3[ii] ge 10.8) then begin
plotsym,3,1.3,/fill ;star
plots,avg_dist[ii],(nacc3[ii]/ntotal[ii]),psym=8,color=djs_icolor('green')
endif
;*******************mass of progenitor
lmass=alog10(mw1.maxhalo)
prmass1=where((mw1.haloid_orig ne 0) and $
              (mw1.haloid_orig ne 1) and $
              (mw1.delay lt 0.) and $
             (mw1.d gt in_dist[ii]) and $
              (lmass le 10.4) $
              and (mw1.d lt out_dist[ii]))

              ;and (mw1.delay lt 0))
;totalmass[ii]=total(mw1[prmass1].maxhalo)
pmass1[ii]=n_elements(prmass1)
time1[ii]=avg(mw1[prmass1].accrtime3)
prmass2=where((mw1.haloid_orig ne 0) and $
              (mw1.haloid_orig ne 1) and $
              (mw1.delay lt 0) and $
              (mw1.d gt in_dist[ii]) and $
              (lmass gt 10.4) $
              and (lmass le 12.5) $
              and (mw1.d lt out_dist[ii]))
pmass2[ii]=n_elements(prmass2)
time2[ii]=avg(mw1[prmass2].accrtime3)

;multiplot
;plots,avg_dist[ii],(pmass1[ii]/ntotal[ii]),psym=4
;xmultiplot,/reset
endfor
plots,avg_dist,(nacc/ntotal),thick=2

per1=strn((total(nacc1)/total(ntotal))*100,length=2)
per2=strn((total(nacc2)/total(ntotal))*100,length=2)
per3=strn((total(nacc3)/total(ntotal))*100,length=2)
totper=strn((total(nacc)/total(ntotal))*100,length=2)

legend,['t < 4 Gyr '+per1+'% of halo','4 < t < 7 Gyr '+per2+'% of halo','7< t Gyr '+per3+' % of halo',totper+' % accreted fraction'],psym=[6,5,2,0],/fill,color=[djs_icolor('cyan'),djs_icolor('magenta'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2],charsize=1.1
multiplot
t1_1=where(time1[*] lt 7)
t1_2=where((time1[*] ge 7) and (time1[*] lt 9))
t1_3=where(time1[*] ge 9)

t2_1=where(time2[*] lt 7)
t2_2=where((time2[*] ge 7) and (time2[*] lt 9))
t2_3=where(time2[*] ge 9)

plotsym,8,0.5,/fill ;square
plot,avg_dist,(pmass1[t1_1]/nacc[t1_1]),psym=8,xrange=[0,150],yrange=[0,1.0],ytitle='total fraction',xtitle='distance (kpc)',color=djs_icolor('indigo')
if (n_elements(t1_2) gt 1) then begin
oplot,avg_dist[t1_2],(pmass1[t1_2]/nacc[t1_2]),psym=8,color=djs_icolor('blue')
endif
if (n_elements(t1_3) gt 1) then begin
oplot,avg_dist[t1_3],(pmass1[t1_3]/nacc[t1_3]),psym=8,color=djs_icolor('cyan')
endif
if (n_elements(t2_1) gt 1) then begin
plotsym,4,0.5,/fill;traiangle
oplot,avg_dist[t2_1],(pmass2[t2_1]/nacc[t2_1]),psym=8,color=djs_icolor('red')
endif
if (n_elements(t2_2) gt 1) then begin
oplot,avg_dist[t2_2],(pmass2[t2_2]/nacc[t2_2]),psym=8,color=djs_icolor('magenta')
endif
if (n_elements(t2_3) gt 1) then begin
oplot,avg_dist[t2_3],(pmass2[t2_3]/nacc[t2_3]),psym=8,color=djs_icolor('pink')
endif
legend,['m<10.4','10.4<m'],psym=[6,5],/fill,color=[djs_icolor('blue'),djs_icolor('red')],symsize=[1.2,1.2],charsize=1.1
multiplot,/reset
device,/close

stop
;*******************fe/h, time,mass
lowearly=where((mw1.accrtime3 lt 5) and (mw1.delay lt 0) and $
               (alog10(mw1.maxhalo) lt 10.2))
lowlate=where((mw1.accrtime3 gt 5) and (mw1.delay lt 0) and $
              (alog10(mw1.maxhalo) lt 10.2))
highearly=where((mw1.accrtime3 lt 5) and (mw1.delay lt 0) and $
               (alog10(mw1.maxhalo) gt 10.2))
highlate=where((mw1.accrtime3 gt 5) and (mw1.delay lt 0) and $
              (alog10(mw1.maxhalo) gt 10.2))
begplot,'font.ps'
erase & multiplot,[2,2]
plot,mw1[lowearly].d,mw1[lowearly].met,psym=4,yrange=[-5,0],xrange=[0,100]
xyouts,'lowearly',5,-4.5
multiplot
plot,mw1[lowlate].d,mw1[lowlate].met,psym=4,yrange=[-5,0],xrange=[0,100]
xyouts,'lowlate',12,-4.5
multiplot
plot,mw1[highearly].d,mw1[highearly].met,psym=4,yrange=[-5,0],xrange=[0,100]
multiplot
plot,mw1[highlate].d,mw1[highlate].met,psym=4,yrange=[-5,0],xrange=[0,100]
multiplot,/reset
endplot
in1=where(mw1.d lt 20)
in2=where((mw1.d gt 20 ) and(mw1.d lt 40))
in3=where((mw1.d gt 40) and (mw1.d lt 100))
begplot,'AMR.ps'
erase & multiplot,[1,3]
plot,mw1[in1].tform,mw1[in1].met,yrange=[-5,0],xrange=[0,13]
multiplot
plot,mw1[in2].tform,mw1[in2].met,yrange=[-5,0],xrange=[0,13]
multiplot
plot,mw1[in3].tform,mw1[in3].met,yrange=[-5,0],xrange=[0,13]
multiplot,/reset
endplot
begplot,'obs_vel.ps'
erase & multiplot,[3,3]
plothist,mw1[in1].vx,xrange=[-300,300]
multiplot
plothist,mw1[in1].vy,xrange=[-300,300]
multiplot
plothist,mw1[in1].vz,xrange=[-300,300]
multiplot
plothist,mw1[in2].vx,xrange=[-300,300]
multiplot
plothist,mw1[in2].vy,xrange=[-300,300]
multiplot
plothist,mw1[in2].vz,xrange=[-300,300]
multiplot
plothist,mw1[in3].vx,xrange=[-300,300]
multiplot
plothist,mw1[in3].vy,xrange=[-300,300]
multiplot
plothist,mw1[in3].vz,xrange=[-300,300]
endplot
stop

end
     
