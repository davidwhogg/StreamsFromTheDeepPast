pro az_halo_origin
;
;
;
;
mw1=mrdfits('/global/data/scr/adi2/mw1.fits',1)
nelem1=n_elements(mw1)
tunit=13.7/.3331
mw1.accrtime3=mw1.accrtime3*tunit
mw1.accrtime1=mw1.accrtime1*tunit


mw1.tform=mw1.tform*tunit


;***************angle up
vec=[mw1.x,mw1.y,mw1.z]
un=[0,0,1]
dot=[0,0,mw1.z]
norm=sqrt((mw1.x*mw1.x)+(mw1.y*mw1.y)+(mw1.z*mw1.z))
an_1=dot/norm
angle=acos(an_1)
tele=where(angle lt 1.05)
mw_ob=mw1[tele]
begplot,'obseverme.ps'
plot,mw_ob.x,mw_ob.z,psym=4,xrange=[-10,10],yrange=[-10,10]
plot,mw_ob.y,mw_ob.z,psym=4,xrange=[-10,10],yrange=[-10,10]
plot,mw_ob.d,mw_ob.z,psym=4,xrange=[-10,10],yrange=[-10,10]
;plot,mw1.d,mw1.z,psym=4,xrange=[-50,50],yrange=[-50,50]
endplot
stop
progen=where((mw1.haloid eq 1) and (abs(mw1.z) gt 4) and (mw1.accrmode eq 4) $
             and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
insitu=where((mw1.haloid eq 1) and (abs(mw1.z) gt 4) and (mw1.accrmode eq 4)$
             and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
diskfeat=where((mw1.haloid eq 1) and (abs(mw1.z) gt 4) and (mw1.accrmode eq 4)$
               and (mw1.comp eq 1))
restfeat=where((mw1.haloid eq 1) and (abs(mw1.z) gt 4) and (mw1.accrmode eq 4)$
               and (mw1.comp eq 4))
  
openw,lun,'particles_prog.dat',/get_lun,width=1
printf,lun,progen+1
free_lun,lun

openw,lun,'particles_insitu.dat',/get_lun,width=1
printf,lun,insitu+1
free_lun,lun

openw,lun,'disk_feature.dat',/get_lun,width=1
printf,lun,diskfeat+1
free_lun,lun

openw,lun,'rest_feature.dat',/get_lun,width=1
printf,lun,restfeat+1
free_lun,lun

begplot,'diskfeature.ps'
plot,mw1[diskfeat].x,mw1[diskfeat].y,psym=4,xrange=[-15,15],yrange=[-15,15], $
  xtitle='x',ytitle='y',title='disk'
plot,mw1[diskfeat].x,mw1[diskfeat].z,psym=4,xrange=[-15,15],yrange=[-15,15], $
xtitle='x',ytitle='z',title='disk'
plot,mw1[diskfeat].y,mw1[diskfeat].z,psym=4,xrange=[-15,15],yrange=[-15,15], $
xtitle='y',ytitle='z',title='disk'
endplot

begplot,'restfeature.ps'
plot,mw1[restfeat].x,mw1[restfeat].y,psym=4,xrange=[-15,15],yrange=[-15,15], $
  xtitle='x',ytitle='y',title='rest'
plot,mw1[restfeat].x,mw1[restfeat].z,psym=4,xrange=[-15,15],yrange=[-15,15], $
xtitle='x',ytitle='z',title='rest'
plot,mw1[restfeat].y,mw1[restfeat].z,psym=4,xrange=[-15,15],yrange=[-15,15], $
xtitle='y',ytitle='z',title='rest'
endplot


bethgas=where((mw1.delay gt 0) and ((mw1.haloid_orig eq 1) or (mw1.haloid_orig eq 0)) and (mw1.comp eq 2))
begplot,'shf_comp.ps'
autohist,mw1[bethgas].tform,title='accreted as stars'
endplot
halo=where((mw1.haloid eq 1) and (abs(mw1.z) gt 7))
mw1=mw1[halo]
nelem=n_elements(mw1)
stars=where((mw1.delay lt 0.) and(mw1.d lt 100))
gas=where((mw1.delay gt 0) and ((mw1.haloid_orig eq 1) or (mw1.haloid_orig eq 0)) and (mw1.d lt 100))
 

begplot,'sfh_hist.ps'
erase & multiplot,[1,2]
autohist,mw1[gas].tform,title='accreted as gas w/o progenitor'
multiplot
autohist,mw1[stars].tform,xtitle='t(Gyr)'
xyouts,5,1000,'accreted as stars'
multiplot,/reset
;plot,mw1[gas].tform,mw1[gas].d,psym=4,yrange=[0,100]
;oplot,mw1[out].tform,mw1[out].d,psym=5
endplot
in1=where(mw1.d lt 10)
inter=where((mw1.d gt 10) and (mw1.d lt 20))
inter2=where((mw1.d gt 20) and (mw1.d lt 40))
out=where((mw1.d gt 40) and (mw1.d lt 100))
begplot,'velo_4.ps'
erase & multiplot,[2,2]
plothist,mw1[in1].rv
oplot,[0,0],!y.crange
multiplot
plothist,mw1[inter].rv
oplot,[0,0],!y.crange
multiplot
plothist,mw1[inter2].rv
oplot,[0,0],!y.crange
multiplot
plothist,mw1[out].rv
oplot,[0,0],!y.crange
multiplot,/default
endplot

xsize=6.
ysize=6.
;********
nbins=20
ntotal=dblarr(1,nbins)
in_dist=dblarr(1,nbins)
out_dist=dblarr(1,nbins)
avg_dist=dblarr(1,nbins)
;these are for delay<0
nm1=dblarr(1,nbins)
nm2=dblarr(1,nbins)
nm3=dblarr(1,nbins)
nm4=dblarr(1,nbins)
nc1=dblarr(1,nbins)
nc2=dblarr(1,nbins)
nc3=dblarr(1,nbins)
nc4=dblarr(1,nbins)
;these are for delay>0
gnm1=dblarr(1,nbins)
gnm2=dblarr(1,nbins)
gnm3=dblarr(1,nbins)
gnm4=dblarr(1,nbins)
gnc1=dblarr(1,nbins)
gnc2=dblarr(1,nbins)
gnc3=dblarr(1,nbins)
gnc4=dblarr(1,nbins)
;these are for progenitors
pnm1=dblarr(1,nbins)
pnm2=dblarr(1,nbins)
pnm3=dblarr(1,nbins)
pnm4=dblarr(1,nbins)
pnc1=dblarr(1,nbins)
pnc2=dblarr(1,nbins)
pnc3=dblarr(1,nbins)
pnc4=dblarr(1,nbins)
;these are for smooth/insitu
snm1=dblarr(1,nbins)
snm2=dblarr(1,nbins)
snm3=dblarr(1,nbins)
snm4=dblarr(1,nbins)
snc1=dblarr(1,nbins)
snc2=dblarr(1,nbins)
snc3=dblarr(1,nbins)
snc4=dblarr(1,nbins)
;features
ft11=dblarr(1,nbins)
ft12=dblarr(1,nbins)
ft13=dblarr(1,nbins)
ft14=dblarr(1,nbins)

ft21=dblarr(1,nbins)
ft22=dblarr(1,nbins)
ft23=dblarr(1,nbins)
ft24=dblarr(1,nbins)

ft31=dblarr(1,nbins)
ft32=dblarr(1,nbins)
ft33=dblarr(1,nbins)
ft34=dblarr(1,nbins)

ft41=dblarr(1,nbins)
ft42=dblarr(1,nbins)
ft43=dblarr(1,nbins)
ft44=dblarr(1,nbins)

pmass1=dblarr(1,nbins)
pmass2=dblarr(1,nbins)
pmass3=dblarr(1,nbins)
pmass4=dblarr(1,nbins)
totalmass=dblarr(1,nbins)
mass=dblarr(1,nbins)
starmass=dblarr(1,nbins)
ind=dblarr(1,nbins)

;acc1=dblarr(1,nbins)
;acc2=dblarr(1,nbins)
;acc3=dblarr(1,nbins)

for ii=1,nbins-1 do begin
out_dist[ii]=10*(ii)
in_dist[ii]=((10*ii)-10)
avg_dist[ii]=(out_dist[ii]+in_dist[ii])/2.

alltotal=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) and (mw1.delay lt 0))
ntotal[ii]=n_elements(alltotal)

;stars accreted after formation 
;what is stars accrmode?
mode1=where((mw1.accrmode eq 1.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nm1[ii]=n_elements(mode1)
mode2=where((mw1.accrmode eq 2.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nm2[ii]=n_elements(mode2)
mode3=where((mw1.accrmode eq 3.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nm3[ii]=n_elements(mode3)
mode4=where((mw1.accrmode eq 4.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nm4[ii]=n_elements(mode4)
;what is stars comp?
comp1=where((mw1.comp eq 1.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nc1[ii]=n_elements(comp1)
comp2=where((mw1.comp eq 2.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nc2[ii]=n_elements(comp2)
comp3=where((mw1.comp eq 3.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
nc3[ii]=n_elements(comp3)
comp4=where((mw1.comp eq 4.) and (mw1.delay lt 0.) and $
            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]))
nc4[ii]=n_elements(comp4)

;stars accreted after before formation
;what is there accmode?
gmode1=where((mw1.accrmode eq 1) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnm1[ii]=n_elements(gmode1)
gmode2=where((mw1.accrmode eq 2) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnm2[ii]=n_elements(gmode2)
gmode3=where((mw1.accrmode eq 3) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnm3[ii]=n_elements(gmode3)
gmode4=where((mw1.accrmode eq 4) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnm4[ii]=n_elements(gmode4)
;what is stars comp?
gcomp1=where((mw1.comp eq 1) and (mw1.delay gt 0.)and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnc1[ii]=n_elements(gcomp1)
gcomp2=where((mw1.comp eq 2) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnc2[ii]=n_elements(gcomp2)
gcomp3=where((mw1.comp eq 3) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnc3[ii]=n_elements(gcomp3)
gcomp4=where((mw1.comp eq 4) and (mw1.delay gt 0.) and $
             (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])) 
gnc4[ii]=n_elements(gcomp4)

;***********************question of origin/progenitor

;accreted from proginator
;what is stars accrmode?
;pmode1=where((mw1.accrmode eq 1.) and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
  ;          and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnm1[ii]=n_elements(pmode1)
;pmode2=where((mw1.accrmode eq 2.)and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
;and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnm2[ii]=n_elements(pmode2)
;pmode3=where((mw1.accrmode eq 3.) and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
  ;          and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnm3[ii]=n_elements(pmode3)
;pmode4=where((mw1.accrmode eq 4.) and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
  ;          and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnm4[ii]=n_elements(pmode4)
;what is stars comp?
;pcomp1=where((mw1.comp eq 1.) and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
  ;          and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnc1[ii]=n_elements(pcomp1)
;pcomp2=where((mw1.comp eq 2.) and $
  ;          (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;           and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnc2[ii]=n_elements(pcomp2)
;pcomp3=where((mw1.comp eq 3.) and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
   ;         and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnc3[ii]=n_elements(pcomp3)
;pcomp4=where((mw1.comp eq 4.) and $
 ;           (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$
  ;          and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
;pnc4[ii]=n_elements(pcomp4)

;insitu and smooth
;what is there accmode?
;smode1=where((mw1.accrmode eq 1) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snm1[ii]=n_elements(smode1)
;smode2=where((mw1.accrmode eq 2) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snm2[ii]=n_elements(smode2)
;smode3=where((mw1.accrmode eq 3) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snm3[ii]=n_elements(smode3)
;smode4=where((mw1.accrmode eq 4) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snm4[ii]=n_elements(smode4)
;what is stars comp?
;scomp1=where((mw1.comp eq 1) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snc1[ii]=n_elements(scomp1)
;scomp2=where((mw1.comp eq 2) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snc2[ii]=n_elements(scomp2)
;scomp3=where((mw1.comp eq 3) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snc3[ii]=n_elements(scomp3)
;scomp4=where((mw1.comp eq 4) and $
 ;            (mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii])$ 
  ;          and ((mw1.haloid_orig eq 0) or (mw1.haloid_orig eq 1)))
;snc4[ii]=n_elements(scomp4)

;****************features at d~ 10kpc
;acc as gas (mode =4)
feat41=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
           and (mw1.accrmode eq 4) and (mw1.comp eq 1))
ft41[ii]=n_elements(feat41)
feat42=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
           and (mw1.accrmode eq 4) and (mw1.comp eq 2))
ft42[ii]=n_elements(feat42)
feat43=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
           and (mw1.accrmode eq 4) and (mw1.comp eq 3))
ft43[ii]=n_elements(feat43)
feat44=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
           and (mw1.accrmode eq 4) and (mw1.comp eq 4))
ft44[ii]=n_elements(feat44)

;acc as stars (mode=2)
;feat21=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 2) and (mw1.comp eq 1))
;ft21[ii]=n_elements(feat21)
;feat22=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
;           and (mw1.accrmode eq 2) and (mw1.comp eq 2))
;ft22[ii]=n_elements(feat22)
;feat23=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
;           and (mw1.accrmode eq 2) and (mw1.comp eq 3))
;ft23[ii]=n_elements(feat23)
;feat24=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
;           and (mw1.accrmode eq 2) and (mw1.comp eq 4))
;ft24[ii]=n_elements(feat24)

;acc during (mode=3)
;feat31=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 3) and (mw1.comp eq 1))
;ft31[ii]=n_elements(feat31)
;feat32=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 3) and (mw1.comp eq 2))
;ft32[ii]=n_elements(feat32)
;feat33=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 3) and (mw1.comp eq 3))
;ft33[ii]=n_elements(feat33)
;feat34=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 3) and (mw1.comp eq 4))
;ft34[ii]=n_elements(feat34)
;acc early (mode=1)
;feat11=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 1) and (mw1.comp eq 1))
;ft11[ii]=n_elements(feat11)
;feat12=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 1) and (mw1.comp eq 2))
;ft12[ii]=n_elements(feat12)
;feat3=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 1) and (mw1.comp eq 3))
;ft13[ii]=n_elements(feat13)
;feat14=where((mw1.d gt in_dist[ii]) and (mw1.d lt out_dist[ii]) $
 ;          and (mw1.accrmode eq 1) and (mw1.comp eq 4))
;ft14[ii]=n_elements(feat14)

;******mass of progenitor
lmass=alog10(mw1.maxhalo)
prmass1=where((mw1.haloid_orig ne 0) and $
              (mw1.haloid_orig ne 1) and $
              (mw1.delay lt 0.) and $
             (mw1.d gt in_dist[ii]) and $
              (lmass le 10.2) $
              and (mw1.d lt out_dist[ii]))
              ;and (mw1.delay lt 0))
;totalmass[ii]=total(mw1[prmass1].maxhalo)
pmass1[ii]=n_elements(prmass1)
prmass2=where((mw1.haloid_orig ne 0) and $
              (mw1.haloid_orig ne 1) and $
              (mw1.delay lt 0) and $
              (mw1.d gt in_dist[ii]) and $
              (lmass gt 10.2) $
              and (lmass le 12.5) $
              and (mw1.d lt out_dist[ii]))
;and (mw1.delay lt 0))
pmass2[ii]=n_elements(prmass2)
prmass3=where((mw1.haloid_orig ne 0) and $
              (mw1.haloid_orig ne 1) and $
              (mw1.delay lt 0) and $
             (mw1.d gt in_dist[ii]) and $
             (lmass gt 10.2)$
              and (lmass le 10.8) $
              and (mw1.d lt out_dist[ii]))
;and (mw1.delay lt 0))
pmass3[ii]=n_elements(prmass3)
prmass4=where((mw1.haloid_orig ne 0) and $
              (mw1.haloid_orig ne 1) and $
              (mw1.delay lt 0) and $
              (mw1.d gt in_dist[ii]) and $
              (lmass gt 10.8) $
              and (mw1.d lt out_dist[ii]))
;and (mw1.delay lt 0))
pmass4[ii]=n_elements(prmass4)


endfor

begplot,'mw1_maxhalo.ps'
autohist,mw1.maxhalo
endplot
set_plot,'ps'
acc1=where((mw1.haloid_orig ne 0) and (mw1.delay lt 0.) and $
           (mw1.accrtime3 lt 4.0))
acc2=where((mw1.haloid_orig ne 0) and (mw1.delay lt 0.) and $
           (mw1.accrtime3 ge 4.)and (mw1.accrtime3 lt 7.))
acc3=where((mw1.haloid_orig ne 0) and (mw1.delay lt 0.) and $
           (mw1.accrtime3 ge 7))
mhalo1=mw1[acc1].maxhalo
mhalo2=mw1[acc2].maxhalo
mhalo3=mw1[acc3].maxhalo

ahalo1_sort=mhalo1(sort(mhalo1))
afrac1=findgen(n_elements(ahalo1_sort))/n_elements(ahalo1_sort)
ahalo2_sort=mhalo2(sort(mhalo2))
afrac2=findgen(n_elements(ahalo2_sort))/n_elements(ahalo2_sort)
ahalo3_sort=mhalo3(sort(mhalo3))
afrac3=findgen(n_elements(ahalo3_sort))/n_elements(ahalo3_sort)

device,file='prog_mass_acc.ps',/color,bits=8
plot,alog10(ahalo1_sort),afrac1,yrange=[0,1],xrange=[8.5,11.2],thick=3,xtitle='Log(donor mass)',ytitle='cumulative fraction'
oplot,alog10(ahalo2_sort),afrac2,color=djs_icolor('blue'),thick=3,linestyle=1
oplot,alog10(ahalo3_sort),afrac3,color=djs_icolor('red'),thick=3,linestyle=2
legend,['t (acc) <4 Gyr','4<t<7','7<t'],color=[djs_icolor('black'),djs_icolor('blue'),djs_icolor('red')],psym=[0,0,0],linestyle=[0,1,2],thick=[3,3,3]
device,/close
;metals MDF
metals1=where((mw1.d gt 0) and (mw1.d le 20) and (mw1.delay lt 0))
;and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
metals2=where((mw1.d gt 20) and (mw1.d le 40) and (mw1.delay lt 0))
;and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
metals3=where((mw1.d gt 40) and (mw1.d le 60) and (mw1.delay lt 0))
;and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
metals4=where((mw1.d gt 60) and (mw1.d le 100)and (mw1.delay lt 0))
;and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))
metals5=where((mw1.d gt 100) and (mw1.d le 150) and (mw1.delay lt 0))
;and (mw1.haloid_orig ne 0) and (mw1.haloid_orig ne 1))

methalo1=mw1[metals1].met
methalo2=mw1[metals2].met
methalo3=mw1[metals3].met
methalo4=mw1[metals4].met
methalo5=mw1[metals5].met

methalo1_sort=methalo1(sort(methalo1))
metfrac1=findgen(n_elements(methalo1_sort))/n_elements(methalo1_sort)

methalo2_sort=methalo2(sort(methalo2))
metfrac2=findgen(n_elements(methalo2_sort))/n_elements(methalo2_sort)

methalo3_sort=methalo3(sort(methalo3))
metfrac3=findgen(n_elements(methalo3_sort))/n_elements(methalo3_sort)

methalo4_sort=methalo4(sort(methalo4))
metfrac4=findgen(n_elements(methalo4_sort))/n_elements(methalo4_sort)

methalo5_sort=methalo5(sort(methalo5))
metfrac5=findgen(n_elements(methalo5_sort))/n_elements(methalo5_sort)

device,file='MDF_stars.ps',/color,bits=8
plot,methalo1_sort,metfrac1,yrange=[0,1],xrange=[-5,0.5],thick=3,xtitle='[Fe/H]',ytitle='cumulative fraction',title='MDF for accreted stars in halo'
oplot,methalo2_sort,metfrac2,linestyle=1,thick=3,color=djs_icolor('blue')
oplot,methalo3_sort,metfrac3,linestyle=2,thick=3,color=djs_icolor('green')
oplot,methalo4_sort,metfrac4,linestyle=3,thick=3,color=djs_icolor('red')
oplot,methalo5_sort,metfrac5,linestyle=4,thick=3,color=djs_icolor('yellow')
legend,['0<d<20','20<d<40','40<d<60','60<d<100','100<d<150'],psym=[0,0,0,0,0],linestyle=[0,1,2,3,4],thick=[3,3,3,3,3],color=[djs_icolor('black'),djs_icolor('blue'),djs_icolor('green'),djs_icolor('red'),djs_icolor('yellow')]
device,/close
;********progenitor mass a function of distance 
m1=where((mw1.d gt 0) and (mw1.d le 10) and (mw1.haloid_orig ne 0) and (mw1.delay lt 0))
m2=where((mw1.d gt 10) and (mw1.d le 20) and (mw1.haloid_orig ne 0) and (mw1.delay lt 0))
m3=where((mw1.d gt 20) and (mw1.d le 40) and (mw1.haloid_orig ne 0) and (mw1.delay lt 0))
m4=where((mw1.d gt 40) and (mw1.d le 100) and (mw1.haloid_orig ne 0) and (mw1.delay lt 0))
mshalo1=mw1[m1].maxhalo
mshalo2=mw1[m2].maxhalo
mshalo3=mw1[m3].maxhalo
mshalo4=mw1[m4].maxhalo
mhalo1_sort=mshalo1(sort(mshalo1))
frac1=findgen(n_elements(mhalo1_sort))/n_elements(mhalo1_sort)

mhalo2_sort=mshalo2(sort(mshalo2))
frac2=findgen(n_elements(mhalo2_sort))/n_elements(mhalo2_sort)

mhalo3_sort=mshalo3(sort(mshalo3))
frac3=findgen(n_elements(mhalo3_sort))/n_elements(mhalo3_sort)

mhalo4_sort=mshalo4(sort(mshalo4))
frac4=findgen(n_elements(mhalo4_sort))/n_elements(mhalo4_sort)

set_plot,'ps'

;device,file='hist_mass.ps',/inches,xsize=xsize,ysize=ysize, $
;  xoffset=0.5,/color,bits=8
;device,/times,/isolatin1
;autohist,alog10(mw1[m1].maxhalo),xrange=[9,11.5],yrange=[0,6d3]
;xyouts,9.3,1000,'0<d<10'
;autohist,alog10(mw1[m2].maxhalo),xrange=[9,11.5],yrange=[0,6d3]
;xyouts,9.3,1000,'10<d<20'
;autohist,alog10(mw1[m3].maxhalo),xrange=[9,11.5],yrange=[0,6d3]
;xyouts,9.3,1000,'20<d<40'
;autohist,alog10(mw1[m4].maxhalo),xrange=[9,11.5],yrange=[0,6d3]
;xyouts,9.3,1000,'40<d<100'
;device,/close
stars=where(mw1.comp eq 2)
nstars=n_elements(stars)

device,file='prog_mass_dist.ps'
plot,alog10(mhalo1_sort),frac1,yrange=[0,1],xrange=[8.5,11.2],thick=3,ytitle='cumulative fraction',xtitle='Log (donor mass)'
oplot,alog10(mhalo2_sort),frac2,color=djs_icolor('blue'),thick=3.,linestyle=1
oplot,alog10(mhalo3_sort),frac3,color=djs_icolor('green'),thick=3.,linestyle=2
oplot,alog10(mhalo4_sort),frac4,color=djs_icolor('red'),thick=3.,linestyle=3
legend,['0<d (kpc) <10','10<d<20','20<d<40','40<d<100'],color=[djs_icolor('black'),djs_icolor('blue'),djs_icolor('green'),djs_icolor('red')],psym=[0,0,0,0],linestyle=[0,1,2,3],thick=[3,3,3,3]
device,/close


set_plot,'ps'
device,file='prog_mass_stars2.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,200],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='stars with progenitors and delay<0'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(pmass1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(pmass2/ntotal),psym=8,color=djs_icolor('red')
;plotsym,3,0.6,/fill ;star
;plots,avg_dist,(pmass3/ntotal),psym=8,color=djs_icolor('green')
;plotsym,0,0.4;circle
;plots,avg_dist,(pmass4/ntotal),psym=8
legend,['m<10.2','10.2<m'],psym=[6,5],/fill,color=[djs_icolor('blue'),djs_icolor('red')],symsize=[1.2,1.2]
device,/close

;make plots!
set_plot,'ps'
device,file='halo_origin_mw1_mode_stars.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='mode of acc delay < 0'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(nm1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(nm2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(nm3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(nm4/ntotal),thick=8
legend,['early','accr as star','during merger','acc as gas'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

device,file='halo_origin_mw1_comp_stars.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='accreted fraction',charsize=1.2,title='comp delay <0'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(nc1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(nc2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(nc3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(nc4/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close
stop
device,file='halo_origin_mw1_comp_gas.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='accreted fraction',charsize=1.2,title='comp delay >0'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(gnc1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(gnc2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(gnc3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(gnc4/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

set_plot,'ps'
device,file='halo_origin_mw1_mode_gas.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='mode delay > 0'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(gnm1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(gnm2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(gnm3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(gnm4/ntotal),thick=8
legend,['early','accr as star','during merger','acc as gas'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close



;**********insitu/progenitor

set_plot,'ps'
device,file='halo_origin_mw1_mode_pro.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='mode with a progenitor'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(pnm1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(pnm2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(pnm3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(pnm4/ntotal),thick=8
legend,['early','accr as star','during merger','acc as gas'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

set_plot,'ps'
device,file='halo_origin_mw1_mode_insitu.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='mode insitu'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(snm1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(snm2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(snm3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(snm4/ntotal),thick=8
legend,['early','accr as star','during merger','acc as gas'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close


device,file='halo_origin_mw1_comp_pro.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='accreted fraction',charsize=1.2,title='comp with prog'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(pnc1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(pnc2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(pnc3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(pnc4/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

device,file='halo_origin_mw1_comp_insitu.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='accreted fraction',charsize=1.2,title='comp insitu'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(snc1/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(snc2/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(snc3/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(snc4/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

device,file='mw1_feature_early.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='all things accreted earlys (mode=1)'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(ft11/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(ft12/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(ft13/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(ft14/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

device,file='mw1_feature_stars.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='all things accreted as stars (mode=2)'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(ft21/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(ft22/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(ft23/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(ft24/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close


device,file='mw1_feature_gas.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='all things accreted as gas (mode=4)'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(ft41/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(ft42/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(ft43/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(ft44/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close

device,file='mw1_feature_merger.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=0.5,/color,bits=8
device,/times,/isolatin1
plot,[0],[0],/nodata,xrange=[0,100],yrange=[0,1],xtitle='distance(kpc)',ytitle='fraction',charsize=1.2,title='all things accreted during merger (mode=3)'
plotsym,8,0.6,/fill ;square
plots,avg_dist,(ft31/ntotal),psym=8,color=djs_icolor('blue')
plotsym,4,0.6,/fill;traiangle
plots,avg_dist,(ft32/ntotal),psym=8,color=djs_icolor('red')
plotsym,3,0.6,/fill ;star
plots,avg_dist,(ft33/ntotal),psym=8,color=djs_icolor('green')
plotsym,0,0.4;circle
plots,avg_dist,(ft34/ntotal),thick=8
legend,['disk','halo','bulge','rest'],psym=[6,5,2,0],/fill,color=[djs_icolor('blue'),djs_icolor('red'),djs_icolor('green'),djs_icolor('black')],symsize=[1.2,1.2,1.2,1.2]
device,/close


stop
end
