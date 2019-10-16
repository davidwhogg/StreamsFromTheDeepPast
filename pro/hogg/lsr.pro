;+
; NAME:
;   lsr
; PURPOSE:
;   Determine the LSR a la Dehnen & Binney (1998) MNRAS 298 387
; OPTIONAL INPUTS:
;   nboot          - number of bootstrap trials; default 10
;   njack          - number of jackknife trials; default 5
;   ngaussfordisk  - number of gaussians to use on the disk; default 1
;   nsubsamp       - number of color subsamples of the stars
;   tol            - convergence tolerance
; KEYWORDS:
;   clobber        - re-make the saveset even if it exists
;   varfreehalo    - don't fix the variance of the halo gaussian
;   meanfreehalo   - don't fix the mean of the halo gaussian
;   bighalo        - use a larger halo velocity dispersion
;   clobber        - re-analyze and overwrite the saveset
;   fakedata       - work on fake data, not real data
;   bestngauss     - use best value of ngauss-for-disk for each data
;                    point; requires many working savesets; brittle
; BUGS:
;   Doesn't do sensible things with the output when ngaussfordisk>1.
;   Relies on all sorts of stuff to run safely with bestngauss!
; REVISION HISTORY:
;   2003-03-11  started - Hogg
;   2003-08-20  modified to use "tied" and jacknife - Hogg & Roweis
;-
pro twod_fit_by_gauss_mixtures, x,y,ycovar,mean,eigenvector,covar
nx= n_elements(x)
ydata= dblarr(2,nx)
ydata[0,*]= x
ydata[1,*]= y
projection= dblarr(2,2,nx)
projection[0,0,*]= 1.0
projection[1,1,*]= 1.0
projected_gauss_mixtures, 1, ydata, ycovar, projection, $
  fitamp, fitxmean, fitxcovar,/quiet
eval= eigenql(0.5*(fitxcovar+transpose(fitxcovar)),eigenvectors=evec)
foo= max(eval,jj)
mean= fitxmean
covar= fitxcovar
eigenvector= evec[*,jj]
return
end

pro lsr, nboot=nboot,njack=njack, $
         varfreehalo=varfreehalo,meanfreehalo=meanfreehalo, $
         bighalo=bighalo,clobber=clobber,fakedata=fakedata, $
         ngaussfordisk=ngaussfordisk,nsubsamp=nsubsamp,tol=tol, $
         bestngauss=bestngauss

; make filenames
if NOT keyword_set(ngaussfordisk) then ngaussfordisk = 1L
prefix= 'lsr_'+strtrim(string(ngaussfordisk),2)
if keyword_set(bighalo) then prefix= prefix+'bighalo'
if keyword_set(fakedata) then prefix= prefix+'fake'
if keyword_set(bestngauss) then prefix= 'lsr_bestngauss'
savefilename= prefix+'.sav'

if (((NOT keyword_set(clobber)) AND $
     (file_test(savefilename))) OR $
    (keyword_set(bestngauss))) then begin
    splog, 'reading '+savefilename
    restore, savefilename
endif else begin

; set defaults
    if NOT keyword_set(nboot) then nboot= 20L
    if NOT keyword_set(njack) then njack= 5L
    if NOT keyword_set(tol) then tol= 1D-7
    ngauss= ngaussfordisk+1
    fixmean= bytarr(ngauss)
    fixvar= bytarr(ngauss)
    if NOT keyword_set(varfreehalo) then $
      fixvar[ngaussfordisk]= 1
    if NOT keyword_set(meanfreehalo) then $
      fixmean[ngaussfordisk]= 1
    seed= -1L

; get and sort sample
    dehnenfilename= '/global/data/hipparcos/hip-dehnen.sav'
    if file_test(dehnenfilename) then begin
        restore, dehnenfilename
    endif else begin
        hip= read_hipparcos(/dehnen,/savefile,/qa)
    endelse
    nhip= n_elements(hip)
    sindx= sort(hip.bvcolor)

; set halo properties
    halomean= [0D0,-2.2D2,0D0]
    halovar= dblarr(3,3)
    halosigma= 1.0D2
    if keyword_set(bighalo) then halosigma= 1.5D2
    for ii=0,2 do halovar[ii,ii]= halosigma^2

; create arrays for outputs
    if(n_elements(nsubsamp) eq 0) then nsubsamp= 20L
    color= dblarr(nsubsamp)
    dcolor= dblarr(nsubsamp)
    colorunit= ' {mag}'
    parameter= dblarr(10,nsubsamp,nboot+1)
    dparameter= dblarr(10,10,nsubsamp)
    relamp= dblarr(nsubsamp)
    drelamp= dblarr(nsubsamp)
    mean= dblarr(3,nsubsamp)
    dmean= dblarr(3,nsubsamp)
    vdisp2= dblarr(nsubsamp)
    dvdisp2= dblarr(nsubsamp)
    vertexdev= dblarr(nsubsamp,nboot+1)
    dvertexdev= dblarr(nsubsamp)
    jacklike= dblarr(nsubsamp,njack)
    hyperamp= dblarr(ngauss,nsubsamp,nboot+njack+1)
    hyperxmean= dblarr(3,ngauss,nsubsamp,nboot+njack+1)
    hyperxcovar= dblarr(3,3,ngauss,nsubsamp,nboot+njack+1)

; create arrays for fake info
    if keyword_set(fakedata) then begin
        fakeparameter= dblarr(10,nsubsamp)
        fakeamp= double([0.990,0.010])
        fakevar= dblarr(3,3,2)
        fakemean= dblarr(3,2)
        fakevdisp2= dblarr(nsubsamp)
        fakelsr= double([-10.00,-5.25,-7.17])
        fakeslope= double([0.0,-0.0086,0.0])
        fakelovar= double([[ 250, 80,-15],[ 80,120,  0],[-15,  0, 30]])
        fakehivar= double([[1340, 96, 11],[ 96,478, 32],[ 11, 32,421]])
        locolor= 0.10
        hicolor= 0.60
    endif

; loop over subsamples and bootstraps
    nstar= round(nhip/nsubsamp)
    for boot=0L,(nboot+njack) do for ii=0L,nsubsamp-1 do begin
        subindx= sindx[ii*nstar:((ii+1L)*nstar < (nhip-1L))]
        nsub= n_elements(subindx)
        splog, 'working on boot',boot,' of subsample',ii
        splog, 'ngaussfordisk is',ngaussfordisk

; check colors
        if boot EQ 0 then begin
            crange= minmax(hip[subindx].bvcolor)
            color[ii]=  (crange[1]+crange[0])/2.0
            dcolor[ii]= (crange[1]-crange[0])/2.0
        endif

; make fake velocities
        if boot EQ 0 and keyword_set(fakedata) then begin
            fakevar[*,*,0]= fakelovar + (fakehivar-fakelovar) $
              *((((color[ii]-locolor)/(hicolor - locolor)) > 0.0) < 1.0)^2
            fakevar[*,*,1]= halovar
            fakevdisp2[ii]= trace(fakevar[*,*,0])
            fakemean[*,0]= fakelsr+fakeslope*fakevdisp2[ii]
            fakemean[*,1]= halomean
            data= fake_catalog(seed,fakeamp,fakemean,fakevar,ndata=nsub)

; transform and add fake measurement errors
            for jj=0L,nsub-1 do begin
                kk= subindx[jj]
                ytrue= transpose(hip[kk].nsm ## data[*,jj])
                ydata= fake_catalog(seed,[1.0],ytrue,hip[kk].vlvbc,/quiet)
                hip[kk].vl= ydata[0]
                hip[kk].vb= ydata[1]
            endfor
        endif

; bootstrap resample?
        if boot EQ 0 then begin
            splog, 'this is the primary fit'
        endif else if boot LE nboot then begin
            splog, 'this is a bootstrap trial'
            bootindx= floor(randomu(seed,nsub)*double(nsub))
            subindx= subindx[bootindx]
        endif else begin
            splog, 'this is a jackknife trial'
            jackindx= shuffle_indx(nsub,seed=seed)
            splitindx= ceil(0.9*nsub)+1
            testsubindx= subindx[jackindx[splitindx:nsub-1]]
            subindx=     subindx[jackindx[0:splitindx-1]]
            nsub= n_elements(subindx)
        endelse
        splog, 'subsample contains',nsub,' stars'

; load first guess
        amp= dblarr(ngauss)
        xmean= dblarr(3,ngauss)
        xcovar= dblarr(3,3,ngauss)
        for gg=0L,ngaussfordisk-1 do begin
            amp[gg]= 0.99/double(ngaussfordisk)
            xmean[*,gg]= 0.0
            for jj=0,2 do xcovar[jj,jj,gg]= 1D2*2.0^double(gg)
        endfor
        amp[ngaussfordisk]= 0.01
        xmean[*,ngaussfordisk]= halomean
        xcovar[*,*,ngaussfordisk]= halovar

; run projected gaussian mixtures
        ydata= transpose([[hip[subindx].vl],[hip[subindx].vb]])
        ycovar= hip[subindx].vlvbc
        projection= hip[subindx].nsm
        projected_tied_mixtures, 2, [ngaussfordisk,1], $
          ydata, ycovar, projection, $
          amp, xmean, xcovar, fixmean=fixmean, fixcovar=fixvar, tol=tol, $
          /quiet
        hyperamp[*,ii,boot]= amp
        hyperxmean[*,*,ii,boot]= xmean
        hyperxcovar[*,*,*,ii,boot]= xcovar

; compute total weighted disk covar
        diskcovar= dblarr(3,3)
        diskamp= 0D0
        for gg=0,ngaussfordisk-1 do begin
            diskamp= diskamp+amp[gg]
            diskcovar= diskcovar+amp[gg]*xcovar[*,*,gg]
        endfor
        diskmean= xmean[*,0]

; compute vertex deviation
        diskcovar= 5d-1*(diskcovar+transpose(diskcovar))
        eval= eigenql(diskcovar,eigenvec=evec)
        foo= max(eval,jj)
        evec1= evec[*,jj]
        vdev1= 1.8d2*atan(evec1[1]/evec1[0])/!DPI

; load up arrays with primary and bootstrap outputs
        if boot LE nboot then begin
            parameter[[0,1,2],ii,boot]= diskmean
            parameter[3,ii,boot]= diskcovar[0,0]
            parameter[4,ii,boot]= diskcovar[1,1]
            parameter[5,ii,boot]= diskcovar[2,2]
            parameter[6,ii,boot]= diskcovar[0,1]
            parameter[7,ii,boot]= diskcovar[0,2]
            parameter[8,ii,boot]= diskcovar[1,2]
            parameter[9,ii,boot]= amp[ngaussfordisk]/total(amp)
            vertexdev[ii,boot]= vdev1
            if boot EQ 0 then begin
                relamp[ii]= amp[ngaussfordisk]/total(amp)
                mean[*,ii]= diskmean
                vdisp2[ii]= trace(diskcovar)
            endif else begin
                p1= reform((parameter[*,ii,0]-parameter[*,ii,boot]),10)
                dparameter[*,*,ii]= dparameter[*,*,ii]+(p1 # p1)
                drelamp[ii]= drelamp[ii]+(relamp[ii]-amp[1]/total(amp))^2
                dmean[*,ii]= dmean[*,ii]+(mean[*,ii]-diskmean)^2
                dvdisp2[ii]= dvdisp2[ii]+(vdisp2[ii]-trace(diskcovar))^2
                dvertexdev[ii]= dvertexdev[ii]+(vertexdev[ii,0]-vdev1)^2
            endelse
        endif

; save fake inputs
        if boot EQ 0 and keyword_set(fakedata) then begin
            fakeparameter[[0,1,2],ii]= fakemean[*,0]
            fakeparameter[3,ii]= fakevar[0,0,0]
            fakeparameter[4,ii]= fakevar[1,1,0]
            fakeparameter[5,ii]= fakevar[2,2,0]
            fakeparameter[6,ii]= fakevar[0,1,0]
            fakeparameter[7,ii]= fakevar[0,2,0]
            fakeparameter[8,ii]= fakevar[1,2,0]
            fakeparameter[9,ii]= fakeamp[1]/total(fakeamp)
        endif

; compute and save jacknife likelihoods
        if boot GT nboot then begin
            splog, 'computing the jackknife likelihood...'
            tydata= transpose([[hip[testsubindx].vl],[hip[testsubindx].vb]])
            tycovar= hip[testsubindx].vlvbc
            tprojection= hip[testsubindx].nsm
            avgloglikedata= 0D0
            projected_tied_mixtures, 2, [ngaussfordisk,1], $
              tydata, tycovar, tprojection, $
              amp, xmean, xcovar, /likeonly, avgloglikedata=avgloglikedata, $
              /quiet
            jacklike[ii,boot-nboot-1]= avgloglikedata
        endif
    endfor

; finish error estimates
    dparameter= dparameter/double(nboot)
    drelamp= sqrt(drelamp/double(nboot))
    dmean= sqrt(dmean/double(nboot))
    dvdisp2= sqrt(dvdisp2/double(nboot))
    dvertexdev= sqrt(dvertexdev/double(nboot))

; make saveset
    save, filename=savefilename
endelse

; stop here if there aren't enough bootstraps to make plots
if nboot LT 3 then return

; setup arrays for making plots and tables
parametertexname= $
  ['\eex\T\,\vvdisk','\eey\T\,\vvdisk','\eez\T\,\vvdisk', $
   '\eex\T\,\VVdisk\,\eex','\eey\T\,\VVdisk\,\eey','\eez\T\,\VVdisk\,\eez', $
   '\eex\T\,\VVdisk\,\eey','\eex\T\,\VVdisk\,\eez','\eey\T\,\VVdisk\,\eez', $
   '\alphahalo']
parametername= $
  ['!8v!dx!n!3','!8v!dy!n!3','!8v!dz!n!3', $
   '!8V!dxx!n!3','!8V!dyy!n!3','!8V!dzz!n!3', $
   '!8V!dxy!n!3','!8V!dxz!n!3','!8V!dyz!n!3', $
   '!7a!3!dhalo!n']
kmstr= '$\mathrm{km\,s^{-1}}$'
kmstr2= '$\mathrm{km^2\,s^{-2}}$'
parametertexunit= [kmstr,kmstr,kmstr, $
                   kmstr2,kmstr2,kmstr2,kmstr2,kmstr2,kmstr2,'']
kmstr= ' {km s!u-1!n}'
kmstr2= kmstr+'!u2!n'
parameterunit= [kmstr,kmstr,kmstr, $
                kmstr2,kmstr2,kmstr2,kmstr2,kmstr2,kmstr2,'']
parameterrange= dblarr(2,10)
vrange= [-14,14]
parameterrange[*,0]= vrange-10
parameterrange[*,1]= vrange-16
parameterrange[*,2]= vrange-5
parameterrange[*,3]= [-0.1,1.1]*1600.0
parameterrange[*,4]= [-0.1,1.1]*700.0
parameterrange[*,5]= [-0.1,1.1]*500.0
rrange= [-0.25,0.25]
parameterrange[*,6]= (rrange+0.15)*sqrt(1600.0*700.0)
parameterrange[*,7]= (rrange+0.0)*sqrt(1600.0*500.0)
parameterrange[*,8]= (rrange+0.0)*sqrt(700.0*500.0)
parameterrange[*,9]= [-0.1,1.1]*0.03
parameterformat= strarr(10)
parameterformat[0:2]= '(F5.1)'
parameterformat[3:8]= '(F5.0)'
parameterformat[9]= '(F6.4)'

; output covariance matrix
ii= 16
openw, wlun,prefix+'_commands.tex',/get_lun
colortexname= '(B-V)'
colorformat= '(F6.3)'
printf, wlun, '\newcommand{\subsamplecolor}{' $
  +string(color[ii]-dcolor[ii],format=colorformat)+'<'+colortexname+'<' $
  +string(color[ii]+dcolor[ii],format=colorformat)+'}'
close, wlun
free_lun, wlun
openw, wlun,prefix+'_parameters_data.tex',/get_lun
dpar1= dblarr(10)
for pp=0L,9L do dpar1[pp]= sqrt(dparameter[pp,pp,ii])
rmatrix= dparameter[*,*,ii]/(dpar1 # dpar1)
for pp=0L,9L do begin
    printf, wlun,parametertexname[pp],parameter[pp,ii,0],dpar1[pp], $
      parametertexunit[pp], $
      format='("$",(A22),"$ & $",'+parameterformat[pp]+',"\pm",'+parameterformat[pp]+',"$ & ",(A24)," &")'
    printf, wlun,rmatrix[pp,*], $
      format='($,10("$",(F5.2),"$",:," & "))'
    printf, wlun, ' \\'
endfor
close, wlun
free_lun, wlun

; plot parameters as a function of color
set_plot, "PS"
psfilename= prefix+'_parameters.ps'
xsize= 3.375  ; actual width of ApJ column!
ysize= 6.0
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
acs= 1.2
hogg_plot_defaults, axis_char_scale=acs
!P.PSYM=  3
!P.MULTI= [0,1,10]
!P.CHARTHICK= 1
!X.OMARGIN= [8,1.5]
!X.THICK= 1
!Y.OMARGIN= 0.4*!X.OMARGIN
!Y.THICK= !X.THICK
colorrange= [-0.2,1.2]
!X.RANGE= colorrange
colorname= '!8(B-V)!3'
!X.TITLE= colorname+colorunit
for pp=0L,9L do begin
    xcharsize= 0.001
    if pp EQ 9 then xcharsize=!X.CHARSIZE
    plot, [0],[0],/nodata, $
      xcharsize=xcharsize, $
      ytitle=parametername[pp]+parameterunit[pp],yrange=parameterrange[*,pp]
    if keyword_set(fakedata) then $
      djs_oploterr, color,fakeparameter[pp,*],linestyle=1, $
      xerr=dcolor
    djs_oploterr, color,parameter[pp,*,0], $
      xerr=dcolor,yerr=sqrt(dparameter[pp,pp,*])
endfor
device,/close

; plot vertex deviation
set_plot, "PS"
psfilename= prefix+'_vertex.ps'
ysize= 2.0
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
hogg_plot_defaults, axis_char_scale=0.75*acs
!P.MULTI= [0,1,1]
!P.CHARTHICK= 1
!X.OMARGIN= [7,1.5]
!X.RANGE= colorrange
!Y.OMARGIN= 0.5*!X.OMARGIN
plot, [0],[0],/nodata, $
  ytitle='vertex deviation  {deg}',yrange=[-10,40]
djs_oploterr, color,vertexdev[*,0], $
  xerr=dcolor,yerr=dvertexdev[*]
device,/close

; make mean-dispersion covariance matrix
Rmatrix= double([[1,0,0,0,0,0,0,0,0,0], $
                 [0,1,0,0,0,0,0,0,0,0], $
                 [0,0,1,0,0,0,0,0,0,0], $
                 [0,0,0,1,1,1,0,0,0,0]])
smallcovar= dblarr(4,4,nsubsamp)
for ii=0L,nsubsamp-1 do begin
    smallcovar[*,*,ii]= Rmatrix##dparameter[*,*,ii]##transpose(Rmatrix)
endfor

; drop blue stars
fitindx= where(color GT 0)
plottingcolor= strarr(nsubsamp)+'grey'
plottingcolor[fitindx]= 'default'

; get LSR and bootstrap it too
S2fit= [0.0,1D4]
dVfit= 0d0
dUfit= 0d0
dWfit= 0d0
for boot=0L,nboot do begin
    nfit= n_elements(fitindx)
    if boot EQ 0 then begin
        tfitindx= fitindx
    endif else begin
        tfitindx= fitindx[floor(double(nfit)*randomu(seed,nfit))]
    endelse
; fit line using gauss mixtures
    twod_fit_by_gauss_mixtures, vdisp2[tfitindx],mean[1,tfitindx], $
      (smallcovar[[3,1],[3,1],*])[*,*,tfitindx], $ ; don't ask
      fitmean,fitvec,fitcovar
    tmpVfit= fitmean[1]+(S2fit-fitmean[0])*fitvec[1]/fitvec[0]
    tmpUfit= total(mean[0,tfitindx]/dmean[0,tfitindx]^2) $
      /total(1D0/dmean[0,tfitindx]^2)
    tmpWfit= total(mean[2,tfitindx]/dmean[2,tfitindx]^2) $
      /total(1D0/dmean[2,tfitindx]^2)
; deal with bootstrapping
    if boot EQ 0 then begin
        Vfit= tmpVfit
        Ufit= replicate(tmpUfit,2)
        Wfit= replicate(tmpWfit,2)
    endif else begin
        dVfit= dVfit+(Vfit[0]-tmpVfit[0])^2
        dUfit= dUfit+(Ufit[0]-tmpUfit)^2
        dWfit= dWfit+(Wfit[0]-tmpWfit)^2
    endelse
endfor
dVfit= sqrt(dVfit/double(nboot))
dUfit= sqrt(dUfit/double(nboot))
dWfit= sqrt(dWfit/double(nboot))

; plot things as a function of dispersion
set_plot, "PS"
psfilename= prefix+'_lsr.ps'
ysize= 5.0
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
hogg_plot_defaults, axis_char_scale=1.25*acs
!P.CHARTHICK= 1
!P.MULTI= [0,1,3]
!X.RANGE= [0.0,2900]
!X.TITLE= 'tr(!8V!3!ddisk!n)  '+kmstr2
!X.OMARGIN= [7,1]
!Y.OMARGIN= [3,0.5]
plot, [0],[0],/nodata, $
  xcharsize=0.001, $
  ytitle=parametername[0]+parameterunit[0],yrange=vrange-10.0
if keyword_set(fakedata) then $
  oplot, S2fit,fakelsr[0]+fakeslope[0]*S2fit, $
  psym=0,linestyle=1,thick=!P.THICK+1
oplot, S2fit,Ufit,psym=0
csz= 1.0
xyouts, S2fit[0],!Y.CRANGE[0]+2.0, $
  '  intercept at '+string(Ufit[0],format='(F5.1)') $
  +' +/- '+string(dUfit,format='(F3.1)'), $
  charsize=csz
bngcsz= 1.0*csz
bngthk= 5.0
bngshft= -0.9 ; hard-coded!
if keyword_set(bestngauss) then $
  xyouts, vdisp2,parameter[0,*,0]+bngshft, $
  strtrim(string(bestngaussfordisk),2), $
  align=0.5,charsize=bngcsz,charthick=bngthk
hogg_oplot_ellipse, vdisp2,parameter[0,*,0],smallcovar[[3,0],[3,0],*], $
  color=plottingcolor
if keyword_set(bestngauss) then $
  xyouts, vdisp2,parameter[0,*,0]+bngshft, $
  strtrim(string(bestngaussfordisk),2), $
  align=0.5,charsize=bngcsz,charthick=bngthk/3.0,color=djs_icolor('white')

plot, [0],[0],/nodata, $
  xcharsize=0.001, $
  ytitle=parametername[1]+parameterunit[1],yrange=vrange-17.0
if keyword_set(fakedata) then $
  oplot, S2fit,fakelsr[1]+fakeslope[1]*S2fit, $
  psym=0,linestyle=1,thick=!P.THICK+1
oplot, S2fit,Vfit,psym=0
;hogg_oplot_ellipse, fitmean[0],fitmean[1],fitcovar,color='red'
xyouts, S2fit[0],!Y.CRANGE[0]+2.0, $
  '  intercept at '+string(Vfit[0],format='(F4.1)') $
  +' +/- '+string(dVfit,format='(F3.1)'), $
  charsize=csz
if keyword_set(bestngauss) then $
  xyouts, vdisp2,parameter[1,*,0]+bngshft, $
  strtrim(string(bestngaussfordisk),2), $
  align=0.5,charsize=bngcsz,charthick=bngthk
hogg_oplot_ellipse, vdisp2,parameter[1,*,0],smallcovar[[3,1],[3,1],*], $
  color=plottingcolor
if keyword_set(bestngauss) then $
  xyouts, vdisp2,parameter[1,*,0]+bngshft, $
  strtrim(string(bestngaussfordisk),2), $
  align=0.5,charsize=bngcsz,charthick=bngthk/3.0,color=djs_icolor('white')

plot, [0],[0],/nodata, $
  ytitle=parametername[2]+parameterunit[2],yrange=vrange-5.0
if keyword_set(fakedata) then $
  oplot, S2fit,fakelsr[2]+fakeslope[2]*S2fit, $
  psym=0,linestyle=1,thick=!P.THICK+1
oplot, S2fit,Wfit,psym=0
xyouts, S2fit[0],!Y.CRANGE[0]+2.0, $
  '  intercept at '+string(Wfit[0],format='(F4.1)') $
  +' +/- '+string(dWfit,format='(F3.1)'), $
  charsize=csz
if keyword_set(bestngauss) then $
  xyouts, vdisp2,parameter[2,*,0]+bngshft, $
  strtrim(string(bestngaussfordisk),2), $
  align=0.5,charsize=bngcsz,charthick=bngthk
hogg_oplot_ellipse, vdisp2,parameter[2,*,0],smallcovar[[3,2],[3,2],*], $
  color=plottingcolor
if keyword_set(bestngauss) then $
  xyouts, vdisp2,parameter[2,*,0]+bngshft, $
  strtrim(string(bestngaussfordisk),2), $
  align=0.5,charsize=bngcsz,charthick=bngthk/3.0,color=djs_icolor('white')

device, /close

end
