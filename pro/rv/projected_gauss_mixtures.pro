;+
; NAME:
;   projected_gauss_mixtures
; PURPOSE:
;   iterate on projected gaussian mixtures
; INPUTS:
;   ngauss     - number of desired gaussians to fit
;   ydata      - [ndimy,ndata] observed velocities
;   ycovar     - [ndimy,ndimy,ndata] observed velocities' errors
;   projection - [ndimx, ndimy, ndata] non-square matrices
;                implementing the projection
;   amp        - [ngauss] list of relative amplitudes
;   xmean      - [ndimx,ngauss] list of initial guesses for the means
;   xcovar     - [ndimx,ndimx,ngauss] list of initial guesses for the
;                 covar
;
; OPTIONAL INPUTS:
;   fixamp     - [ngauss] list of integers: 0 to update amp, 1 not to
;   fixmean    - [ngauss] list of integers: 0 to update xmean, 1 not to
;   fixcovar   - [ngauss] list of integers: 0 to update xcovar, 1 not
;                to
;   tol        - tolerance (convergence iff difference in avgloglike <
;                tol)
;   maxiter    - maximum number of iterations
;   xguess     - ?
;   xcovarguess- ?
;   constrain  - constrain the any variance to be larger than this
;                value

;
; KEYWORDS:
;   likeonly   - only compute and return the avgloglikedata; do not
;                update
;   quiet      - run quietly
;   debug      - print debugging messages (doesn't do anything now?)
;   plot       - plot something in 2D?
;   plot3d     - plot something in 3D?
;   stupid_init- use stupid initial conditions (see
;                projected_initialize)
;   rotation   - ?
;   interact   - interact with plots?
;
; OUTPUTS:
;   avloglikedata - Average log-likelihood of the data after
;                   convergence
;   +updata xmean, xcovar and amp...
;
; BUGS:
;   Code has many undocumented inputs and keywords.
;
; REVISION HISTORY:
;   2003-02-18  written - Blanton and Roweis
;   2008-06-10  patched up + docs - Bovy
;   2008-07-17  added 'constrain' input - Bovy
;-
pro projected_gauss_mixtures_step, ydata, ycovar, projection, $
                                   logamp, xmean, xcovar, $
                                   fixamp, fixmean, fixcovar, $
                                   avgloglikedata=avgloglikedata, $
                                   quiet=quiet, debug=debug, $
                                   xguess=xguess, xcovarguess=xcovarguess, $
                                   bestposterior=bestposterior, $
                                   xbestguess=xbestguess, likeonly=likeonly, $
                                   constrain=constrain
ON_ERROR, 2


ngauss=n_elements(logamp)
ndimy=n_elements(ycovar)/n_elements(ydata)
ndimx=n_elements(xcovar)/n_elements(xmean)
ndata=n_elements(ydata)/ndimy
twopiterm=0.5*double(ndimy)*alog(2.*!DPI)

minlogamp=-(machar(/double)).xmax+alog(double(ndata))
logamp=logamp > minlogamp

loglikedata=0.D
logposterior=dblarr(ngauss,ndata)
for ii=0L, ndata-1L do begin
;   calculate log likelihood for each gaussian for this data point
    loglike=dblarr(ngauss)
    for jj=0L, ngauss-1L do begin
            tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,jj]# $
                                projection[*,*,ii]+ycovar[*,*,ii])
            delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,jj]
            loglike[jj]=logamp[jj]+0.5*alog(determ(tinv,/double,/check) > (machar(/double)).xmin)- $
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

;Don't do aything if nothing needs to be updated
    IF (fixamp[jj] EQ 1 AND fixmean[jj] EQ 1 AND fixcovar[jj] EQ 1 AND (NOT keyword_set(quiet))) THEN BEGIN
        splog, 'Skipping component '+strtrim(string(jj),2)
        CONTINUE
    ENDIF

;   accumulate stuff for this gaussian
    curr_xmean=dblarr(ndimx)
    curr_xcovar1=dblarr(ndimx,ndimx)
    curr_xcovar2=dblarr(ndimx,ndimx)
    
;AVOID UNDERFLOW
        wignore= 1D-50;see below, ignore weights smaller than this*mean(weight)
        posterior= exp(logposterior[jj,*])
        norm_post= total(posterior,/DOUBLE)
        tolerance= (norm_post > 1)*(machar(/DOUBLE)).xmin
        weight= dblarr(ndata)+ (posterior > tolerance)/norm_post
        weight= weight * (weight GE (wignore * mean(weight)))

;        splog, 'ignored weights '+strtrim(string(n_elements(where(weight LE (wignore * mean(weight))))),2)+' in component '+strtrim(string(jj),2)

    for ii=0L, ndata-1L do begin
        tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,jj]# $
                            projection[*,*,ii]+ycovar[*,*,ii])
        delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,jj]
;       curr_xguess corresponds to b_ij in the paper
        curr_xguess=xmean[*,jj]+xcovar[*,*,jj]#projection[*,*,ii]#tinv#delta
        curr_xmean=curr_xmean+weight[ii]*curr_xguess
        curr_xcovar1=curr_xcovar1+weight[ii]*(curr_xguess##curr_xguess)

;       fix this if you fix the below 
;        curr_xcovar2=curr_xcovar2+ $
;          posterior*projection[*,*,ii]#tinv#transpose(projection[*,*,ii])
;       fixed?
        projtinv=(projection[*,*,ii]#tinv#transpose(projection[*,*,ii]))
;       NOTE THE NECESSITY TO MAINTAIN STRICT SYMMETRY
;       IDL sucks and does not guarantee that the above line returns a
;       perfectly symmetric matrix, even if tinv is perfectly
;       symmetric
        projtinv=0.5*(projtinv+transpose(projtinv))
        curr_xcovar2=curr_xcovar2+weight[ii]*projtinv

        if(n_elements(bestposterior) gt 0) then begin
            if(posterior[ii] gt bestposterior[ii]) then begin
                bestposterior[ii]=posterior[ii]
                xbestguess[*,ii]=curr_xguess
            endif
        endif
        if(n_elements(xguess) gt 0) then begin
            xguess[*,ii,jj]=curr_xguess
        endif

        if(n_elements(xcovarguess) gt 0) then begin
;           fix this if you fix the above 
;           fixed?
            xcovarguess[*,*,ii,jj]=xcovar[*,*,jj]- $
              xcovar[*,*,jj]#projtinv#transpose(xcovar[*,*,jj])
        endif
    endfor

;   now update parameters 
    if fixamp[jj] NE 1 then $
      logamp[jj]=logsum(logposterior[jj,*],/double)-alog(double(ndata))
; mincloud legacy ---logamp[jj]=logsum([logposterior[jj,*],alog(mincloud_weight)],/double)
    if(logamp[jj] le minlogamp) then begin
        splog,'component '+strtrim(string(jj),2)+' has zero data points'
        logamp[jj]=minlogamp
        fixamp[jj]=1
        fixmean[jj]=1
        fixcovar[jj]=1
    endif
    new_xmean=curr_xmean
;;Avoid infinities
    IF ~FINITE(total(curr_xmean)) THEN fixmean[jj]= 1
    if fixmean[jj] NE 1 then xmean[*,jj]=new_xmean
;   this next line may be inefficient but it works whether or not mean
;   is fixed ...
    IF ~FINITE(total(curr_xcovar2)) THEN fixcovar[jj]= 1
    IF ~FINITE(total(curr_xcovar1)) THEN fixcovar[jj]= 1
    if fixcovar[jj] NE 1 then begin
        upxcovar=xcovar[*,*,jj]#curr_xcovar2#transpose(xcovar[*,*,jj])
;       see annoyed comments above about idl's suckiness in symmetry
        upxcovar=(0.5D)*(upxcovar+transpose(upxcovar))
        tmpcovar= xcovar[*,*,jj]- $
          upxcovar+$
          curr_xcovar1+$
          xmean[*,jj]##xmean[*,jj]- $
      xmean[*,jj]##new_xmean-new_xmean##xmean[*,jj]
;;Make sure the variances are larger than constrain if set
;;Explicitly symmetrize
        tmpcovar= 0.5D*(tmpcovar+transpose(tmpcovar))
        IF ~keyword_set(constrain) THEN BEGIN
            xcovar[*,*,jj]=tmpcovar
        ENDIF ELSE BEGIN
            ;;Diagonalize
            l=eigenql(tmpcovar,/double)
            IF (total((l > constrain) - l) NE 0D) THEN BEGIN
                l=eigenql(tmpcovar,/double,eigenvectors=eigenvec)
                l= (l>constrain)
                ldiag=dblarr(ndimx,ndimx)
                FOR ii=0L, ndimx-1 DO ldiag[ii,ii]=l[ii]
                tmpcovar= transpose(eigenvec)##ldiag##eigenvec
;;Explicitly symmetrize
                tmpcovar= 0.5D*(tmpcovar+transpose(tmpcovar))
            ENDIF
            xcovar[*,*,jj]=tmpcovar
        ENDELSE    
    endif

    nonsym=max(abs(xcovar[*,*,jj]-transpose(xcovar[*,*,jj]))) gt 1.D-16
    if(nonsym) then begin
        splog, 'WARNING: xcovar got non-symmetric'
        splog, 'WARNING: IDL sux'
        wait, 10
    endif

    if(keyword_set(debug)) then $
      splog,'jj= '+string(jj)+' ; curr_ndata= '+string(curr_ndata)
endfor

return
end
;
pro projected_gauss_mixtures, ngauss, ydata, ycovar, projection, $
                              amp, xmean, xcovar, $
                              fixamp=fixamp,fixmean=fixmean, $
                              fixcovar=fixcovar, $
                              avgloglikedata=avgloglikedata, $
                              quiet=quiet, tol=tol, maxiter=maxiter, $
                              debug=debug, plot=plot, plot3d=plot3d, $
                              xguess=xguess,xcovarguess=xcovarguess, $
                              stupid_init=stupid_init, $
                              rotation=rotation, interact=interact, $
                              likeonly=likeonly, constrain=constrain

ON_ERROR, 2

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
		projected_gauss_3dplot, ydata, ycovar, projection, $
                  amp, xmean, xcovar, rotation=rotation, $
                  linelength=linelength, basename='qa3dinit', $
                  /mkinfile,interact=interact, $
                  xbestguess=xbestguess
endif

logamp=alog(amp > (machar(/double)).xmin)

while(niter lt maxiter and diff gt tol) do begin
    if(keyword_set(plot)) then begin
        projected_gauss_2dplot, ydata, ycovar, projection, amp, xmean, $
          xcovar,/nogreyscale
;,filename='qa.'+string(niter,format='(i5.5)')+'.ps'
    endif


    projected_gauss_mixtures_step, ydata, ycovar, projection, $
      logamp, xmean, xcovar, $
      fixamp, fixmean, fixcovar, $
      avgloglikedata=avgloglikedata, $
      quiet=quiet, debug=debug, xbestguess=xbestguess, $
      bestposterior=bestposterior, likeonly=likeonly, constrain=constrain


    if keyword_set(likeonly) then begin
        if(NOT keyword_set(quiet)) then splog, 'computed likelihood, returning'
        return
    endif

    if(keyword_set(plot3d)) then begin
        projected_gauss_3dplot, ydata, ycovar, projection, $
          amp, xmean, xcovar, rotation=rotation, $
          linelength=linelength, xbestguess=xbestguess, $
          /mkinfile, interact=interact, $
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

amp=exp(logamp)

if(n_elements(xguess) gt 0) then begin
    projected_gauss_mixtures_step, ydata, ycovar, projection, $
      logamp, xmean, xcovar, $
      fixamp, fixmean, fixcovar, $
      avgloglikedata=avgloglikedata, $
      quiet=quiet, debug=debug, xguess=xguess, xcovarguess=xcovarguess
endif

return
end
