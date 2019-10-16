;+
; NAME:
;   projected_gauss_mixtures
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
;   likeonly - only compute and return the avgloglikedata; do not update
; OUTPUTS:
; BUGS:
;   Code has many undocumented inputs and keywords.
; REVISION HISTORY:
;   2003-02-18  written - Blanton and Roweis
;-
pro projected_gauss_mixtures_step, ydata, ycovar, projection, $
                                   logamp, xmean, xcovar, $
                                   fixamp, fixmean, fixcovar, $
                                   avgloglikedata=avgloglikedata, $
                                   quiet=quiet, debug=debug, $
                                   xguess=xguess, xcovarguess=xcovarguess, $
                                   bestposterior=bestposterior, $
                                   xbestguess=xbestguess, likeonly=likeonly

ngauss=n_elements(logamp)
ndimy=n_elements(ycovar)/n_elements(ydata)
ndimx=n_elements(xcovar)/n_elements(xmean)
ndata=n_elements(ydata)/ndimy
twopiterm=0.5*double(ndimy)*alog(2.*!DPI)

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
        curr_xcovar2=curr_xcovar2+ $
          posterior*projection[*,*,ii]#tinv#transpose(projection[*,*,ii])

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

;   now update parameters 
    if fixamp[jj] NE 1 then logamp[jj]=logsum(logposterior[jj,*])
    curr_ndata=exp(logamp[jj])
    logamp[jj]=logamp[jj]-alog(double(ndata))
    new_xmean=curr_xmean/curr_ndata
    if fixmean[jj] NE 1 then xmean[*,jj]=new_xmean
;   this next line may be inefficient but it works whether or not mean
;   is fixed ...
    if fixcovar[jj] NE 1 then xcovar[*,*,jj]=xcovar[*,*,jj]- $
      xcovar[*,*,jj]#curr_xcovar2#transpose(xcovar[*,*,jj])/curr_ndata+ $
      curr_xcovar1/curr_ndata+xmean[*,jj]##xmean[*,jj]- $
      xmean[*,jj]##new_xmean-new_xmean##xmean[*,jj]
    nonsym=max(abs(xcovar[*,*,jj]-transpose(xcovar[*,*,jj]))) gt 1.D-10
    if(nonsym) then begin
        splog, 'WARNING: xcovar got non-symmetric'
        splog, 'WARNING: IDL sux'
        stop
    endif

    if(keyword_set(debug)) then $
      splog,'jj= '+string(jj)+' ; curr_ndata= '+string(curr_ndata)
endfor

return
end
;
pro projected_gauss_mixtures, ngauss, ydata, ycovar, projection, $
                              amp, xmean, xcovar, $
                              fixamp=fixamp,fixmean=fixmean,fixcovar=fixcovar, $
                              avgloglikedata=avgloglikedata, $
                              quiet=quiet, tol=tol, maxiter=maxiter, $
                              debug=debug, plot=plot, plot3d=plot3d, $
                              xguess=xguess,xcovarguess=xcovarguess, $
                              stupid_init=stupid_init, $
                              rotation=rotation, interact=interact, $
                              likeonly=likeonly

if(NOT keyword_set(fixamp)) then fixamp= bytarr(ngauss)
if(NOT keyword_set(fixmean)) then fixmean= bytarr(ngauss)
if(NOT keyword_set(fixcovar)) then fixcovar= bytarr(ngauss)
if(NOT keyword_set(tol)) then tol=1.D-6
if(NOT keyword_set(maxiter)) then maxiter=1000000000L

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
        projected_kmean, ngauss, ydata, ycovar, projection, amp, $
          xmean, xcovar, quiet=quiet
    endif else begin
        projected_initialize, ngauss, ydata, ycovar, projection, amp, $
          xmean, xcovar, quiet=quiet
    endelse
endif
  
diff=2.*tol
niter=0
if(keyword_set(plot3d)) then begin
		projected_gauss_3dplot, ydata, ycovar, projection, amp, xmean, $
      xcovar, rotation=rotation, linelength=linelength, $
		  basename='qa3dinit',/mkinfile,interact=interact, $
      xbestguess=xbestguess
endif
while(niter lt maxiter and diff gt tol) do begin
    if(keyword_set(plot)) then begin
        projected_gauss_2dplot, ydata, ycovar, projection, amp, xmean, $
          xcovar,/nogreyscale
;,filename='qa.'+string(niter,format='(i5.5)')+'.ps'
    endif

    logamp=alog(amp)
    projected_gauss_mixtures_step, ydata, ycovar, projection, $
      logamp, xmean, xcovar, $
      fixamp, fixmean, fixcovar, $
      avgloglikedata=avgloglikedata, $
      quiet=quiet, debug=debug, xbestguess=xbestguess, $
      bestposterior=bestposterior, likeonly=likeonly
    amp=exp(logamp)

    if keyword_set(likeonly) then begin
        if(NOT keyword_set(quiet)) then splog, 'computed likelihood, returning'
        return
    endif

    if(keyword_set(plot3d)) then begin
		    projected_gauss_3dplot, ydata, ycovar, projection, amp, xmean, $
          xcovar, rotation=rotation, linelength=linelength, $
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
    projected_gauss_mixtures_step, ydata, ycovar, projection, $
      logamp, xmean, xcovar, $
      fixamp, fixmean, fixcovar, $
      avgloglikedata=avgloglikedata, $
      quiet=quiet, debug=debug, xguess=xguess, xcovarguess=xcovarguess
endif

return
end
