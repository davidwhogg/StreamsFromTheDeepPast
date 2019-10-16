;+
; NAME:
;   projected_tied_mixtures
; PURPOSE:
;   iterate on projected gaussian mixtures
; INPUTS:
;   ngauss - number of desired gaussians to fit
;   ydata - [ndimy,ndata] observed velocities
;   ycovar - [ndimy,ndimy,ndata] observed velocities' errors
;   projection - [ndimx, ndimy, ndata] non-square matrices
;                implementing the projection
; OPTIONAL INPUTS:
;   fixamp   - [ngauss] list of integers: 0 to update amp, 1 not to
;   fixmean  - [ngauss] list of integers: 0 to update xmean, 1 not to
;   fixcovar - [ngauss] list of integers: 0 to update xcovar, 1 not to
; KEYWORDS:
;   /likeonly - only compute and return the avgloglikedata; do not update
; OUTPUTS:
; BUGS:
;   Code has many undocumented inputs.
;   What is the ordering of means and covars in return / fixamp / etc.??
; REVISION HISTORY:
;   2003-08-18  written - Blanton and Roweis
;-
;function safe_exp, expon
;return,exp(expon > (-(machar(/double)).xmax))
;end
;
pro projected_tied_mixtures_step, ngauss_per_center, ydata, ycovar, $
                                  projection, logamp, xmean, xcovar, $
                                  fixamp, fixcenter, fixcovar, $
                                  avgloglikedata=avgloglikedata, $
                                  quiet=quiet, debug=debug, $
                                  xguess=xguess, xcovarguess=xcovarguess, $
                                  bestposterior=bestposterior, $
                                  xbestguess=xbestguess, likeonly=likeonly

; mincloud legacy
; weight of minimum cloud
; mincloud_weight=1.d-3
; variance of minimum cloud (hard-wired for now)
; mincloud_variance=0.001

ngauss=n_elements(logamp)
ncenter=n_elements(ngauss_per_center)
ndimy=n_elements(ycovar)/n_elements(ydata)
ndimx=n_elements(xcovar)/n_elements(xmean)
ndata=n_elements(ydata)/ndimy
twopiterm=0.5*double(ndimy)*alog(2.*!DPI)

minlogamp=-(machar(/double)).xmax+alog(double(ndata))
logamp=logamp > minlogamp

; make mapping from gaussian to center
new_xmean=dblarr(ndimx,ngauss)
curr_ndata=dblarr(ngauss)
center_curr_xmean=dblarr(ndimx,ncenter)
center_curr_ndata=dblarr(ndimx,ndimx,ncenter)
center=lonarr(ngauss)
jjgauss_start=0L
for jj=0L, ncenter-1L do begin
    jjgauss_end=jjgauss_start+ngauss_per_center[jj]-1L
    center[jjgauss_start:jjgauss_end]=jj
    jjgauss_start=jjgauss_end+1L
endfor

loglikedata=0.D
logposterior=dblarr(ngauss,ndata)
for ii=0L, ndata-1L do begin
;   calculate log likelihood for each gaussian for this data point
    loglike=dblarr(ngauss)
    for jj=0L, ngauss-1L do begin
        tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,jj]# $
                            projection[*,*,ii]+ycovar[*,*,ii])
        delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,jj]
        loglike[jj]=logamp[jj]+0.5*alog(determ(tinv,/double))- $
          0.5*transpose(delta)#tinv#delta-twopiterm
    endfor

;   calculate posterior probability to be in each gaussian
    curr_loglikedata=logsum(loglike,/double)
    loglikedata=loglikedata+curr_loglikedata
    logposterior[*,ii]=loglike-curr_loglikedata
endfor

avgloglikedata=loglikedata/double(ndata)
if(NOT keyword_set(quiet)) then $
  splog,'avgloglikedata= '+strtrim(string(avgloglikedata,format='(f40.9)'),2)

if keyword_set(likeonly) then begin
    if(NOT keyword_set(quiet)) then splog, 'computed likelihood, returning'
    return
endif

if(arg_present(xbestguess)) then begin
    bestposterior=dblarr(ndata)-1.D
    xbestguess=dblarr(ndimx,ndata)
endif

for jj=0L, ngauss-1L do begin

;   accumulate stuff for this gaussian
    curr_xmean=dblarr(ndimx)
    curr_xcovar1=dblarr(ndimx,ndimx)
    curr_xcovar2=dblarr(ndimx,ndimx)
    
    for ii=0L, ndata-1L do begin
        tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,jj]# $
                            projection[*,*,ii]+ycovar[*,*,ii])
        delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,jj]
        posterior=exp(logposterior[jj,ii])
;       curr_xguess corresponds to b_ij in the paper
        curr_xguess=xmean[*,jj]+xcovar[*,*,jj]#projection[*,*,ii]#tinv#delta
        curr_xmean=curr_xmean+posterior*curr_xguess
        curr_xcovar1=curr_xcovar1+posterior*(curr_xguess##curr_xguess)
;       fix this if you fix the below 
        projtinv=(projection[*,*,ii]#tinv#transpose(projection[*,*,ii]))
;       NOTE THE NECESSITY TO MAINTAIN STRICT SYMMETRY
;       IDL sucks and does not guarantee that the above line returns a
;       perfectly symmetric matrix, even if tinv is perfectly
;       symmetric
        projtinv=0.5*(projtinv+transpose(projtinv))
        curr_xcovar2=curr_xcovar2+posterior*projtinv

        if(n_elements(bestposterior) gt 0) then begin
				    if(posterior gt bestposterior[ii]) then begin
						    bestposterior[ii]=posterior
                xbestguess[*,ii]=curr_xguess
            endif
				endif
        if(n_elements(xguess) gt 0) then begin
            xguess[*,ii,jj]=curr_xguess
        endif

        if(n_elements(xcovarguess) gt 0) then begin
;           fix this if you fix the above 
            xcovarguess[*,*,ii,jj]=xcovar[*,*,jj]- $
              xcovar[*,*,jj]#projection[*,*,ii]#tinv# $
              transpose(projection[*,*,ii])#transpose(xcovar[*,*,jj])
        endif
    endfor
;   add the minimum cloud contribution to curr_xmean
; mincloud legacy ---curr_xmean=curr_xmean+xmean[*,jj]*mincloud_weight

;   now update parameters 
    if fixamp[jj] NE 1 then $
      logamp[jj]=logsum(logposterior[jj,*],/double)
; mincloud legacy ---logamp[jj]=logsum([logposterior[jj,*],alog(mincloud_weight)],/double)
    if(logamp[jj] le minlogamp) then begin
        splog,'component '+string(jj)+' has zero data points'
        logamp[jj]=minlogamp
        fixamp[jj]=1
        fixcovar[jj]=1
    endif
    curr_ndata[jj]=exp(logamp[jj])
    logamp[jj]=logamp[jj]-alog(double(ndata))
; mincloud legacy ---     alog(double(ndata)+double(ngauss)*mincloud_weight)
;   update mean the same for everything with the same center
    if(curr_ndata[jj] gt 0.) then $
      new_xmean[*,jj]=curr_xmean/curr_ndata[jj] $
    else $
      new_xmean[*,jj]=!VALUES.D_NAN
    invxcovar=invert(xcovar[*,*,jj])
    center_curr_xmean[*,center[jj]]=center_curr_xmean[*,center[jj]]+ $
      invxcovar#curr_xmean
    center_curr_ndata[*,*,center[jj]]=center_curr_ndata[*,*,center[jj]]+ $
      invxcovar*curr_ndata[jj]
;   we are waiting until all the center_curr_xmeans are created to add the 
;   xmean#xmean terms
    ;if(keyword_set(stopatend)) then stop
    if(fixcovar[jj] NE 1 AND curr_ndata[jj] gt 0.) then begin
        upxcovar=xcovar[*,*,jj]#curr_xcovar2#transpose(xcovar[*,*,jj])
;       see annoyed comments above about idl's suckiness in symmetry
        upxcovar=(0.5D)*(upxcovar+transpose(upxcovar))
        xcovar[*,*,jj]=xcovar[*,*,jj]- $
          upxcovar/curr_ndata[jj]+ $
          curr_xcovar1/curr_ndata[jj]
    endif
    
    nonsym=max(abs(xcovar[*,*,jj]-transpose(xcovar[*,*,jj]))) gt 1.D-16
    if(nonsym) then begin
        splog, 'WARNING: xcovar got non-symmetric'
        splog, 'WARNING: IDL sux'
        splog, string(max(abs(xcovar[*,*,jj]-transpose(xcovar[*,*,jj]))))
        stop
    endif

    if(keyword_set(debug)) then $
      splog,'jj= '+string(jj)+' ; curr_ndata[jj]= '+string(curr_ndata[jj])

    if(keyword_set(stopatend)) then stop

endfor

; now assign means of center to each gaussian
jjgauss_start=0L
for jj=0L, ncenter-1L do begin
    if(fixcenter[jj] NE 1) then begin
        jjgauss_end=jjgauss_start+ngauss_per_center[jj]-1L
        if(trace(center_curr_ndata[*,*,jj]) gt 0.) then begin
            thiscenter=invert(center_curr_ndata[*,*,jj])# $
              center_curr_xmean[*,jj]
            xmean[*,jjgauss_start:jjgauss_end]= $
              thiscenter#replicate(1.D,ngauss_per_center[jj])
        endif
        jjgauss_start=jjgauss_end+1L
    endif
endfor

; Now we add the -xmean#xmean term, now that we have the fully updated
; xmeans for each center.
for jj=0L, ngauss-1L do begin
    if(fixcovar[jj] NE 1) then begin
        if(curr_ndata[jj] gt 0.) then $
          xcovar[*,*,jj]= xcovar[*,*,jj]+xmean[*,jj]##xmean[*,jj]- $
          new_xmean[*,jj]##xmean[*,jj]-xmean[*,jj]##new_xmean[*,jj]
;       resymmetrizing, AGAIN, because IDL sux
        xcovar[*,*,jj]=(0.5D)*(xcovar[*,*,jj]+transpose(xcovar[*,*,jj]))
        nonsym=max(abs(xcovar[*,*,jj]-transpose(xcovar[*,*,jj]))) gt 1.D-16
        if(nonsym) then begin
            splog, 'WARNING: xcovar got non-symmetric'
            splog, 'WARNING: IDL sux'
            splog, string(max(abs(xcovar[*,*,jj]-transpose(xcovar[*,*,jj]))))
            stop
        endif
    endif
endfor

end
;
pro projected_tied_mixtures, ncenter, in_ngauss_per_center, ydata, ycovar, $
                             projection, $
                             amp, xmean, xcovar, $
                             fixamp=in_fixamp,fixmean=in_fixmean, $
                             fixcovar=in_fixcovar, $
                             avgloglikedata=avgloglikedata, $
                             quiet=quiet, tol=tol, maxiter=maxiter, $
                             debug=debug, plot=plot, plot3d=plot3d, $
                             xguess=xguess,xcovarguess=xcovarguess, $
                             stupid_init=stupid_init, likeonly=likeonly, $
                             rotation=rotation, interact=interact

if(n_elements(in_ngauss_per_center) eq 1) then $
  ngauss_per_center=replicate(in_ngauss_per_center, ncenter) $
else $
  ngauss_per_center=in_ngauss_per_center
ngauss=long(total(ngauss_per_center))
if(NOT keyword_set(in_fixamp)) then in_fixamp= bytarr(ngauss)
if(NOT keyword_set(in_fixmean)) then in_fixmean= bytarr(ngauss)
if(NOT keyword_set(in_fixcovar)) then in_fixcovar= bytarr(ngauss)
fixamp=in_fixamp
fixmean=in_fixmean
fixcovar=in_fixcovar
if(NOT keyword_set(tol)) then tol=1.D-6
if(NOT keyword_set(maxiter)) then maxiter=1000000000L
if(NOT keyword_set(covarinitfactor)) then covarinitfactor=2.

; convert fixmean to fix center
; if any gaussian is fixed, all gaussians for that center are fixed
fixcenter=bytarr(ncenter)
jjgauss_start=0L
for jj=0L, ncenter-1L do begin
    jjgauss_end=jjgauss_start+ngauss_per_center[jj]-1L
    totcenter=long(total(fixmean[jjgauss_start:jjgauss_end]))
    if(totcenter ne 0 AND totcenter ne ngauss_per_center[jj]) then $
      splog,'WARNING: some gaussians fixed, some not, for center '+string(jj)
    fixcenter[jj]=totcenter gt 0
    jjgauss_start=jjgauss_end+1L
endfor

; these *appear* to do nothing
;if(arg_present(xguess) and n_elements(xguess) eq 0) then begin
;endif

ndimy=n_elements(ycovar)/n_elements(ydata)
ndata=n_elements(ydata)/ndimy
ndimx=n_elements(projection)/ndata/ndimy

; set initial conditions if necessary
; it fails in general
if(n_elements(amp) ne ngauss or $
   n_elements(xmean) eq 0 or $
   n_elements(xcovar) eq 0) $
  then begin
    if(NOT keyword_set(stupid_init)) then begin
        projected_kmean, ncenter, ydata, ycovar, projection, center_amp, $
          center_xmean, center_xcovar
    endif else begin
        projected_initialize, ncenter, ydata, ycovar, projection, center_amp, $
          center_xmean, center_xcovar
    endelse
    amp=dblarr(ngauss)
    xmean=dblarr(ndimx,ngauss)
    xcovar=dblarr(ndimx,ndimx,ngauss)
    jjgauss_start=0L
    for jj=0L, ncenter-1L do begin
        jjgauss_end=jjgauss_start+ngauss_per_center[jj]-1L
        amp[jjgauss_start:jjgauss_end]=center_amp/double(ngauss_per_center[jj])
        xmean[*,jjgauss_start:jjgauss_end]= $
          center_xmean[*,jj]#replicate(1.,ngauss_per_center[jj])
        eval=eigenql(center_xcovar[*,*,jj],/double,eigenvectors=evec)
        for jjgauss= jjgauss_start, jjgauss_end do begin
            curr_factor= $
              exp(-alog(covarinitfactor)+2.*alog(covarinitfactor)* $
                  (double(jjgauss-jjgauss_start)+0.5)/ $
                  double(ngauss_per_center[jj]))
            curr_eval=eval*curr_factor
            diag=dblarr(ndimx,ndimx)
            diag[lindgen(ndimx),lindgen(ndimx)]=curr_eval
            xcovar[*,*,jjgauss]=evec#diag#transpose(evec)
        endfor
        jjgauss_start=jjgauss_end+1L
    endfor
endif else begin
    jjgauss_start=0L
    for jj=0L, ncenter-1L do begin
        jjgauss_end=jjgauss_start+ngauss_per_center[jj]-1L
        for jjgauss=jjgauss_start, jjgauss_end-1L do begin
            if(abs(max(xmean[*,jjgauss]-xmean[*,jjgauss+1L])) gt 0.) then begin
                message, 'not all means of tied gaussians are equal on input'
            endif
        endfor
        jjgauss_start=jjgauss_end+1L
    endfor
endelse 
  
diff=2.*tol
niter=0
if(keyword_set(plot3d)) then begin
		projected_gauss_3dplot, ydata, ycovar, projection, amp, xmean, $
      xcovar, rotation=rotation, linelength=linelength, $
		  basename='qa3dinit',/mkinfile,interact=interact, $
      xbestguess=xbestguess
endif
logamp=alog(amp)
while(niter lt maxiter and diff gt tol) do begin
    if(keyword_set(plot)) then begin
        projected_gauss_2dplot, ydata, ycovar, projection, exp(logamp), $
          xmean, xcovar,/nogreyscale
;,filename='qa.'+string(niter,format='(i5.5)')+'.ps'
    endif

    projected_tied_mixtures_step, ngauss_per_center, ydata, ycovar, $
      projection, logamp, xmean, xcovar, $
      fixamp, fixcenter, fixcovar, $
      avgloglikedata=avgloglikedata, $
      quiet=quiet, debug=debug, xbestguess=xbestguess, $
      bestposterior=bestposterior, likeonly=likeonly

    if keyword_set(likeonly) then begin
        if(NOT keyword_set(quiet)) then splog, 'computed likelihood, returning'
        return
    endif

    if(keyword_set(plot3d)) then begin
		    projected_gauss_3dplot, ydata, ycovar, projection, exp(logamp), $
          xmean, xcovar, rotation=rotation, linelength=linelength, $
          xbestguess=xbestguess, /mkinfile, interact=interact, $
		      basename='qa3d.'+string(niter,format='(i5.5)')
    endif
      
;    if(keyword_set(debug)) then begin
;        splog,'xmean= '+string(xmean)
;        splog,'xcovar= '+string(xcovar)
;    endif
    if(niter gt 0) then $
      diff=(avgloglikedata-oldavgloglikedata)
    if(diff lt 0.) then begin
        splog,'WARNING: log likelihood decreased.'
        splog,'WARNING: Sam cannot do math or Mike cannot write code.'
        stop
    endif
    oldavgloglikedata=avgloglikedata
    niter=niter+1L
endwhile

if(n_elements(xguess) gt 0) then begin
    projected_gauss_mixtures_step, ngauss_per_center, ydata, ycovar, $
      projection, logamp, xmean, xcovar, $
      fixamp, fixmean, fixcovar, $
      avgloglikedata=avgloglikedata, $
      quiet=quiet, debug=debug, xguess=xguess, xcovarguess=xcovarguess
endif

amp=exp(logamp)

end
