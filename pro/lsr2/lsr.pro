;+
;   NAME:
;      lsr
;   PURPOSE:
;      derive the lsr by removing the moving groups
;   CALLING SEQUENCE:
;   INPUT:
;      groupcut         - if the posterior prob for a data point to be
;                         part of a moving group is larger than this
;                         value, consider it part of a moving group
;      relcut           - if set, do a relative cut, considering a
;                         data point part of a moving groups if its
;                         posterior prob for that group is > relcut
;                         times the second
;      nboot            - number of bootstrap trials to use for error
;                         estimation
;      nsubsamp         - number of color subsamples to use
;      tol              - tolerance for proj_gauss_mixtures
;      halomean         - default [0D0,-2.2D2,0D0]
;      bighalo          - use a halo with a large vdisp (150 km/s)
;      seed             - seed for random number generator
;      prefix           - prefix for plotfilenames
;   KEYWORDS:
;      remove_ngc1901   - If set, consider the Gaussian centered o
;                         NGC1901 a moving group
;      dropblue         - drop the bluest stars when fitting for the
;                         lsr
;      S500             - If set, drop all samples with S < 500 from
;                         the fit
;      ex1000           - If set, drop all samples with S < 1000 from
;                         the fit
;      formal           - use Hipparcos formal errors
;      dontplot
;      silent
;   OUTPUT:
;      out              - [U,sigU,V,sigV,W,sigW,frac]
;   REVISION HISTORY:
;      2009-05-25 - Written - Bovy (partially based on Hogg's lsr code)
;-
pro twod_fit_by_gauss_mixtures, x,y,ycovar,mean,eigenvector,covar
nx= n_elements(x)
ydata= dblarr(2,nx)
ydata[0,*]= x
ydata[1,*]= y
projection= dblarr(2,2,nx)
projection[0,0,*]= 1.0
projection[1,1,*]= 1.0
fitamp= 1
meanx= (moment(x))[0]
meany= (moment(y))[0]
varx= (moment(x))[1]
vary= (moment(y))[1]
fitxmean= [meanx,meany]
fitxcovar= [[varx*varx, 0.],[0.,vary*vary]]
projected_gauss_mixtures, 1, ydata, ycovar, projection, $
  fitamp, fitxmean, fitxcovar,/quiet
eval= eigenql(0.5*(fitxcovar+transpose(fitxcovar)),eigenvectors=evec)
foo= max(eval,jj)
mean= fitxmean
covar= fitxcovar
eigenvector= evec[*,jj]
return
end
PRO LSR, remove_ngc1901=remove_ngc1901,groupcut=groupcut,relcut=relcut, $
         nboot=nboot, nsubsamp=nsubsamp,tol=tol, seed=seed, $
         dropblue=dropblue, prefix=prefix, S500=S500, formal=formal, $
         out=out, dontplot=dontplot, silent=silent, ex1000=ex1000

IF ~keyword_set(groupcut) THEN groupcut=0.5
IF ~keyword_set(nboot) THEN nboot= 20L
IF ~keyword_set(nsubsamp) THEN nsubsamp= 20L
IF ~keyword_set(tol) THEN tol=1D-6
IF ~keyword_set(seed) THEN seed=-1L

IF keyword_set(formal) THEN $
  restore, filename='hip2-aumer_formal.sav' ELSE $
  restore, filename='hip2-aumer.sav'
nhip= n_elements(hip.hip)
ydata= transpose([[hip.vl], [hip.vb]])
ycovar= hip.vlvbc
projection= hip.nsm

restore, filename='fitv10_bestV_V4.0_sample6.sav'
ngauss= n_elements(amp)

;;First construct the sample, determine the data points which aren't
;;part of a moving group (based on a given probability)
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
;;Make the sample cut
IF keyword_set(remove_ngc1901) THEN groups=[0,1,2,3,4,6,7] ELSE $
  groups=[1,2,3,4,6,7]
ngroups= n_elements(groups)
sample_indx= lonarr(nhip)
nsample_indx= 0L
FOR ii=0L, nhip-1 DO BEGIN
    ingroup= 0
    FOR jj=0L, ngroups-1 DO BEGIN
        IF keyword_set(relcut) THEN BEGIN
            IF logpost[groups[jj],ii] NE max(logpost[*,ii]) THEN CONTINUE
            sortindx= sort(logpost[*,ii])
            IF exp(logpost[groups[jj],ii]) GE $
              relcut*exp(logpost[sortindx[ngauss-2],ii]) THEN BEGIN
                ingroup= 1
                BREAK
            ENDIF
        ENDIF ELSE IF exp(logpost[groups[jj],ii]) GE groupcut THEN BEGIN
            ingroup= 1
            BREAK
        ENDIF
    ENDFOR
    IF ingroup EQ 0 THEN BEGIN
        sample_indx[nsample_indx]= ii
        nsample_indx+= 1
    ENDIF
ENDFOR
sample_indx= sample_indx[0:nsample_indx-1]
ydata= transpose([[hip[sample_indx].vl], [hip[sample_indx].vb]])
ycovar= hip[sample_indx].vlvbc
projection= hip[sample_indx].nsm
frac_excluded= (nhip-nsample_indx)/double(nhip)
IF ~keyword_set(silent) THEN splog, strtrim(string(nsample_indx),2)+" stars not considered to be in a moving group"

;;Then derive the lsr for this sample
;;Sample
nhip= nsample_indx
sindx= sort(hip[sample_indx].bvcolor)

; set halo properties
IF ~keyword_set(halomean) THEN halomean= [0D0,-2.2D2,0D0]
halovar= dblarr(3,3)
halosigma= 1.0D2
if keyword_set(bighalo) then halosigma= 1.5D2
for ii=0,2 do halovar[ii,ii]= halosigma^2

;;Set up arrays for the output
color= dblarr(nsubsamp)
dcolor= dblarr(nsubsamp)
colorunit= ' [mag]'
parameter= dblarr(10,nsubsamp,nboot+1)
dparameter= dblarr(10,10,nsubsamp)
mean= dblarr(3,nsubsamp)
dmean= dblarr(3,nsubsamp)
vdisp2= dblarr(nsubsamp)
dvdisp2= dblarr(nsubsamp)
vertexdev= dblarr(nsubsamp,nboot+1)
dvertexdev= dblarr(nsubsamp)


; loop over subsamples and bootstraps
nstar= round(nhip/nsubsamp)
for boot=0L,nboot do for ii=0L,nsubsamp-1 do begin
    subindx= sindx[ii*nstar:((ii+1L)*nstar < (nhip-1L))]
    nsub= n_elements(subindx)
    IF ~keyword_set(silent) THEN splog, 'working on boot',boot,' of subsample',ii
    
; colors
    if boot EQ 0 then begin
        crange= minmax(hip[sample_indx[subindx]].bvcolor)
        color[ii]=  (crange[1]+crange[0])/2.0
        dcolor[ii]= (crange[1]-crange[0])/2.0
    endif

; bootstrap resample?
    if boot EQ 0 then begin
        IF ~keyword_set(silent) THEN splog, 'this is the primary fit'
    endif else begin
        IF ~keyword_set(silent) THEN splog, 'this is a bootstrap trial'
        bootindx= floor(randomu(seed,nsub)*double(nsub))
        subindx= subindx[bootindx]
    endelse
    IF ~keyword_set(silent) THEN splog, 'subsample contains',nsub,' stars'

; load first guess
    ngauss= 2
    xamp= dblarr(ngauss)
    xmean= dblarr(3,ngauss)
    xcovar= dblarr(3,3,ngauss)
    xamp[0]= 0.99
    xmean[*,0]= 0.0
    for jj=0,2 do xcovar[jj,jj,0]= 1D2*2.
    xamp[1]= 0.01
    xmean[*,1]= halomean
    xcovar[*,*,1]= halovar

; run projected gaussian mixtures
    thisydata= transpose([[hip[sample_indx[subindx]].vl], [hip[sample_indx[subindx]].vb]])
    thisycovar= hip[sample_indx[subindx]].vlvbc
    thisprojection= hip[sample_indx[subindx]].nsm
    ;;fix halo except for the amp
    fixamp= bytarr(ngauss)
    fixmean= bytarr(ngauss)
    fixcovar= bytarr(ngauss)
    fixcovar[1]= 1
    fixmean[1]= 1
    logfile=''
    splitnmerge=0
    projected_gauss_mixtures_c, ngauss, thisydata, thisycovar, $
      thisprojection, xamp, xmean, xcovar, fixamp=fixamp, $
      fixmean=fixmean, fixcovar=fixcovar, tol=tol,/quiet, $
      logfile=logfile, splitnmerge=splitnmerge, w=0.
    ;if boot EQ 0 THEN BEGIN
    ;    print, xamp[1]
    ;    wait, 1
    ;ENDIF
; compute vertex deviation
    diskcovar= 5d-1*(xcovar[*,*,0]+transpose(xcovar[*,*,0]))
    eval= eigenql(diskcovar,eigenvec=evec)
    foo= max(eval,jj)
    evec1= evec[*,jj]
    vdev1= 1.8d2*atan(evec1[1]/evec1[0])/!DPI
    
; load up arrays with primary and bootstrap outputs
    parameter[[0,1,2],ii,boot]= xmean[*,0]
    parameter[3,ii,boot]= diskcovar[0,0]
    parameter[4,ii,boot]= diskcovar[1,1]
    parameter[5,ii,boot]= diskcovar[2,2]
    parameter[6,ii,boot]= diskcovar[0,1]
    parameter[7,ii,boot]= diskcovar[0,2]
    parameter[8,ii,boot]= diskcovar[1,2]
    parameter[9,ii,boot]= amp[1]/total(amp)
    vertexdev[ii,boot]= vdev1
    if boot EQ 0 then begin
        mean[*,ii]= xmean[*,0]
        vdisp2[ii]= trace(diskcovar)
    endif else begin
        p1= reform((parameter[*,ii,0]-parameter[*,ii,boot]),10)
        dparameter[*,*,ii]= dparameter[*,*,ii]+(p1 # p1)
        dmean[*,ii]= dmean[*,ii]+(mean[*,ii]-xmean[*,0])^2
        dvdisp2[ii]= dvdisp2[ii]+(vdisp2[ii]-trace(diskcovar))^2
        dvertexdev[ii]= dvertexdev[ii]+(vertexdev[ii,0]-vdev1)^2
    endelse
endfor


; finish error estimates
dparameter= dparameter/double(nboot)
dmean= sqrt(dmean/double(nboot))
dvdisp2= sqrt(dvdisp2/double(nboot))
dvertexdev= sqrt(dvertexdev/double(nboot))


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
IF keyword_set(dropblue) THEN fitindx= where(color GT 0.1 AND vdisp2 GT 225.) ELSE $
  IF keyword_set(s500) THEN fitindx= where(color GT 0.1 AND vdisp2 GT 500.) ELSE $
  IF keyword_set(ex1000) THEN fitindx= where(color GT 0.1 AND vdisp2 GT 1000.) ELSE $
  fitindx= lindgen(nsubsamp)
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
IF ~keyword_set(dontplot) THEN BEGIN
    set_plot, "PS"
    IF keyword_set(remove_ngc1901) THEN prefix+= '_ngc1901'
    IF keyword_set(dropblue) THEN prefix += '_dropblue' ELSE IF keyword_set(S500) THEN prefix+= '_s500' ELSE IF keyword_set(ex1000) THEN prefix+= '_ex1000'
    IF keyword_set(relcut) THEN prefix+= '_relcut'+strtrim(string(relcut,format='(F3.1)'),2) ELSE $
      prefix+= '_groupcut'+strtrim(string(groupcut,format='(F3.1)'),2)
    IF keyword_set(formal) THEN prefix+= '_formal'
    psfilename= prefix+'.ps'
    parametername= $
      ['!8v!dx!n!3','!8v!dy!n!3','!8v!dz!n!3', $
       '!8V!dxx!n!3','!8V!dyy!n!3','!8V!dzz!n!3', $
       '!8V!dxy!n!3','!8V!dxz!n!3','!8V!dyz!n!3', $
       '!7a!3!dhalo!n']
    kmstr= ' [km s!u-1!n]'
    kmstr2= ' [km!u2!n s!u-2!n]'
    parameterunit= [kmstr,kmstr,kmstr, $
                    kmstr2,kmstr2,kmstr2,kmstr2,kmstr2,kmstr2,'']
    vrange= [-14,14]
    xsize= 3.375                ; actual width of ApJ column!
    ysize= 5.0
    acs= 1.2
    device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
      xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
    hogg_plot_defaults, axis_char_scale=1.25*acs
    !P.CHARTHICK= 2
    !X.THICK=2
    !Y.THICK=2
    !P.MULTI= [0,1,3]
    !X.RANGE= [0.0,2900]
    !X.TITLE= 'tr(!8V!3!ddisk!n)  '+kmstr2
    !X.OMARGIN= [7,1]
    !Y.OMARGIN= [3,0.5]
    plot, [0],[0],/nodata, $
      xcharsize=0.001, $
      ytitle=parametername[0]+parameterunit[0],yrange=vrange-10.0
    oplot, S2fit,Ufit,psym=0
    csz= 1.0
    csztitle= .9
    IF groupcut GE 1.0 THEN xyouts, S2fit[0]+400.,!Y.CRANGE[1]+1.0, $
      'No moving-group stars excluded', $
      charsize=csztitle
    IF groupcut EQ 0.5 AND keyword_set(dropblue) THEN xyouts, S2fit[0],!Y.CRANGE[1]+1.0, $
      'Moving-group stars excluded at '+textoidl('p_{MG}')+' = 0.5', $
      charsize=csztitle
    bngcsz= 1.0*csz
    xyouts, S2fit[0],!Y.CRANGE[0]+2.0, $
      '  intercept at '+string(Ufit[0],format='(F5.1)') $
      +' +/- '+string(dUfit,format='(F3.1)'), $
      charsize=csz
    bngcsz= 1.0*csz
    bngthk= 5.0
    bngshft= -0.9               ; hard-coded!
    hogg_oplot_ellipse, vdisp2,parameter[0,*,0],smallcovar[[3,0],[3,0],*], $
      color=plottingcolor
    
    plot, [0],[0],/nodata, $
      xcharsize=0.001, $
      ytitle=parametername[1]+parameterunit[1],yrange=vrange-17.0
    oplot, S2fit,Vfit,psym=0
    xyouts, S2fit[0],!Y.CRANGE[0]+2.0, $
      '  intercept at '+string(Vfit[0],format='(F4.1)') $
      +' +/- '+string(dVfit,format='(F3.1)'), $
      charsize=csz
    hogg_oplot_ellipse, vdisp2,parameter[1,*,0],smallcovar[[3,1],[3,1],*], $
      color=plottingcolor
    
    plot, [0],[0],/nodata, $
      ytitle=parametername[2]+parameterunit[2],yrange=vrange-5.0
    oplot, S2fit,Wfit,psym=0
    xyouts, S2fit[0],!Y.CRANGE[0]+2.0, $
      '  intercept at '+string(Wfit[0],format='(F4.1)') $
      +' +/- '+string(dWfit,format='(F3.1)'), $
      charsize=csz
    hogg_oplot_ellipse, vdisp2,parameter[2,*,0],smallcovar[[3,2],[3,2],*], $
      color=plottingcolor
    
    device, /close
ENDIF

out= [Ufit[0],dUfit,Vfit[0],dVfit,Wfit[0],dWfit,frac_excluded]

END
