;+
;   NAME:
;      fitv
;
;   PURPOSE:
;      Fits the velocity distribution in the solar neighborhood using
;      Gaussians
;
;   CALLING SEQUENCE:
;      fitv [, ngauss=ngauss, halomean=halomean, halovar=halovar,
;      fixhalo=fixhalo, nsubsamp=nsubsamp, colors=colors, tol=tol, 
;      debug=debug, best_ngauss= best_ngauss,
;      use_best_ngauss=use_best_ngauss, njack=njack, jacksub=jacksub,
;      jackconv=jackconv, sample=sample,maxstag=magstag,
;      initial=initial,maxitt=maxitt, , thinmean=thinmean,
;      thincovar=thincovar, thickmean=thickmean,
;      thickcovar=thickcovar, maxstag=maxstag,
;      Vconstraint=Vconstraint, bestV=bestV, Vjacksub=Vjacksub,
;      Vnjack=Vnjack, justcalc=justcalc, log=log, giants=giants,
;      all=all, minamp=minamp
;
;   INPUTS:
;      ngauss   - number of Gaussians to use (default 3, GE 1), if
;                 array, then it must have the number of gaussians for
;                 each individual subsample
;      halomean - if set, set the mean of the halo (default: -220 y)
;      halovar  - if set, set the var of the halo (default: 100 XX,
;                 YY, ZZ)
;      fixhalo  - if set, fix the halo-Gaussian (default: free)
;      nsubsamp - number of color subsamples to use (default: 20)
;      colors   - alternative way to specify the colors, array that
;                 holds the breakpoints for the color
;                 (i.e. [sup(sample1)=min(sample2), sup(sample2),...]
;      tol      - tolerance for convergence (default: 1D-7)
;      njack    - number of jackknife trials in the determination of
;                 bestngauss (default: 5)
;      jacksub  - fraction of the sample to set aside in the jackknife
;                 determination of bestngauss (default .1)
;      jackconv - number of times the loglikelihood has to decrease to
;                 indicate that we're past it's peak
;      sample   - use only this subsample,
;                 i.e. 1,2,3,4,5=giants,6=all MS, 7=all with giants
;      itttol   - tolerance for convergence of /initial (obsolete)
;      maxitt   - maximal number of iterations for intial procedure
;      maxstag  - maximal number of retries before initial gives up
;      thinmean - sets the initial mean of the thin disk mean
;      thincovar- sets the initial covariance of the thin disk
;      thickmean- sets the initial mean of the thick disk
;      thickcovar- sets the initial covar of the thick disk
;      Vconstraint - sets the constraint on the variance
;                    (i.e. variance > Vconstraint, default 0 )
;      Vjacksub - fraction of the sample to set aside in the jackknife
;                 determination of bestV
;      Vnjack   - number of jackknife trials in the determination of
;                 bestV (default: 5), not effective if /vanilla is set
;      minamp   - minimum value of the amplitude that should be
;                 strived for. In number of datapoints, e.g. 5!!!
;
;   KEYWORDS:
;      best_ngauss     - calculate best number of gaussians to fit with
;                        for each sample.
;      use_best_ngauss - use the best_ngauss output, which can also
;                        just be inputted
;      debug    - print debugging messages and insert breakpoints
;      initial  - try improving the distribution by using the output
;                 of one iteration as the initial condition of the
;                 next, modifying one of the low-amplitude ones into
;                 one of the high amplitude ones
;      bestV     - jackknife the sample, to calculate the likelihood
;                  for each Vconstraint (if Vconstraint is given,
;                  jackknifing is only performed for that Vconstraint)
;      justcalc  - Don't plot or write parameters
;      log       - Write output to log-file
;      giants    - Include GIANTS as a sample
;      all       - Include the union of all samples as a sample
;      realjack  - do a real jackknife, use .1 and 10 subsamples
;      vanilla   - do 'vanilla' jackknife
;      oldhip    - Use old reduction of Hipparcos (default: new van
;                  Leeuwen reduction)
;      tyc       - Use old Hipparcos and Tycho proper motions
;
;   OUTPUTS:
;     Plots the velocity distribution function. + savefile +
;     parameters file
;
;   EXAMPLE:
;      fitv, ngauss=10, colors=[0.0D,0.4D,0.6D], /bestV,
;      Vconstraint=4D,/giants,/all,Vjacksub=0.01, /vanilla, /log, sample=6 
;
;   REVISION HISTORY:
;      06-02-2008  - Started Bovy
;      12-02-2008  - implemented different data sets
;-
PRO FITV, ngauss=ngauss, halomean=halomean, halovar=halovar, fixhalo=fixhalo, $
          nsubsamp=nsubsamp, colors=colors, tol=tol, debug=debug, $
          best_ngauss= best_ngauss, use_best_ngauss=use_best_ngauss, $
          njack=njack, jacksub=jacksub, jackconv=jackconv, sample=sample, $
          initial=initial, itttol=itttol, maxitt=maxitt, thinmean=thinmean, $
          thincovar=thincovar, thickmean=thickmean, thickcovar=thickcovar, $
          maxstag=maxstag, Vconstraint=Vconstraint, bestV=bestV, $
          Vjacksub=Vjacksub, Vnjack=Vnjack, justcalc=justcalc, log=log, $
          giants=giants, all=all, minamp=minamp, realjack=realjack, $
          vanilla=vanilla, oldhip=oldhip, tyc=tyc

;Set dimension
d= 3

;;Use best_ngauss.sav if use_best_ngauss is given and ngauss is not
;;set.
bestngaussfilename= 'bestngauss.sav'
IF (keyword_set(use_best_ngauss)) THEN BEGIN
    IF ~file_test(bestngaussfilename) THEN BEGIN
        IF ~keyword_set(ngauss) THEN BEGIN
            splog, 'Either you give me the numbers of gaussians for each '+ $
              'subsample or the file '+bestngaussfilename+' must exist'
            splog, 'returning...'
            RETURN
        ENDIF
    ENDIF ELSE BEGIN
        restore, bestngaussfilename
    ENDELSE
ENDIF

;Set defaults
IF ~keyword_set(halomean) THEN halomean= [0D0,-2.2D2,0D0]
IF ~keyword_set(halovar) THEN BEGIN
    halovar= dblarr(d,d)
    halosigma= 1.0D2
    FOR ii=0L,d-1 DO halovar[ii,ii]= halosigma^2
ENDIF
IF ~keyword_set(thinmean) THEN thinmean= dblarr(d)
IF ~keyword_set(thincovar) THEN BEGIN
    thincovar= dblarr(d,d)
    thinsigma= 20D
    FOR ii=0L,d-1 DO thincovar[ii,ii]= thinsigma^2
ENDIF
IF ~keyword_set(thickmean) THEN thickmean= [0D0,-40D,0D0]
IF ~keyword_set(thickcovar) THEN BEGIN
    thickcovar= dblarr(d,d)
    thicksigma= 75D
    FOR ii=0L,d-1 DO thickcovar[ii,ii]= thicksigma^2
ENDIF
IF ~keyword_set(fixhalo) THEN fixhalo=0
IF ~keyword_set(nsubsamp) THEN nsubsamp=20L
IF keyword_set(colors) THEN nsubsamp= n_elements(colors)+1
IF keyword_set(giants) THEN nsubsamp++
IF keyword_set(all) THEN nsubsamp+= 2
IF ~keyword_set(ngauss) THEN ngauss=3
ngauss= ngauss > 1
IF (n_elements(ngauss) EQ 1 AND nsubsamp GT 1) THEN $
  ngauss= lonarr(nsubsamp) + ngauss
IF ~keyword_set(tol) THEN tol=1D-6
IF ~keyword_set(debug) THEN debug= 0
IF ~keyword_set(jacksub) THEN jacksub= .1D
IF (keyword_set(best_ngauss) AND ~keyword_set(njack)) THEN njack=5
IF ~keyword_set(jackconv) THEN jackconv= 5L
IF ~keyword_set(best_ngauss) THEN njack= 0L
IF ~keyword_set(sample) THEN sample= 0
IF ~keyword_set(itttol) THEN itttol=1D-7
IF ~keyword_set(maxitt) THEN maxitt=10000L
IF ~keyword_set(initial) THEN maxitt=1L
IF ~keyword_set(maxstag) THEN maxstag= 5L
IF (~keyword_set(Vconstraint) AND ~keyword_set(bestV)) THEN Vconstraint= 0D
IF ~keyword_set(Vnjack) THEN Vnjack= 0L
IF ~keyword_set(Vjacksub) THEN Vjacksub= 0.1D
IF keyword_set(vanilla) THEN Vnjack= floor(1D0/Vjacksub)
IF ~keyword_set(minamp) THEN minamp= 10L

ngauss= long(ngauss)

;;All arrays will have a size equal to the max of the ngauss
maxngauss= max(ngauss)

endngauss= lonarr(nsubsamp)+ngauss
IF keyword_set(best_ngauss) THEN BEGIN
    endngauss= 100
    jacklike= dblarr(nsubsamp,endngauss,njack+1)
    endngauss=lonarr(nsubsamp) + endngauss
ENDIF
down= 0;obsolete

seed= -1L

;;Set up save-file
IF ~keyword_set(best_ngauss) THEN prefix='fitv'+ $
  strtrim(string(maxngauss),2) $
  ELSE prefix= 'fitvbestngauss'
IF keyword_set(bestV) THEN BEGIN 
    prefix= prefix+'_bestV'
    IF keyword_set(Vconstraint) THEN prefix= prefix + '_V'+strtrim(string(Vconstraint,format='(f4.1)'),2)
ENDIF ELSE BEGIN
    prefix= prefix+'_V'+strtrim(string(Vconstraint,format='(f4.1)'),2)
ENDELSE
IF keyword_set(sample) THEN prefix= prefix+'_sample'+strtrim(string(sample),2)
savefilename= prefix+'.sav'
;IF keyword_set(best_ngauss) THEN savefilename=bestngaussfilename
IF file_test(savefilename) THEN BEGIN
    splog, 'reading '+savefilename
    restore, savefilename
ENDIF ELSE BEGIN

;;Open log-file
    IF keyword_set(log) THEN OPENW, loglun, prefix+'.log', /GET_LUN

;;Get sample, sort it by color
    IF keyword_set(oldhip) THEN BEGIN
        dehnen_giants_filename= '/global/data/hipparcos/hip-dehnen-giants.sav'
        splog, 'Using the old Hipparcos reduction'
    ENDIF ELSE IF keyword_set(tyc) THEN BEGIN
        dehnen_giants_filename= '~jb2777/streams/pro/rv/tyc-dehnen-giants.sav'
        splog, 'Using the old Hipparcos reduction with Tycho proper motions'
    ENDIF ELSE BEGIN
        dehnen_giants_filename= '/global/data/hipparcos/hipnewredux/hip2-dehnen-giants.sav'
        splog, 'Using the new Hipparcos reduction'
    ENDELSE
    IF (keyword_set(giants) OR keyword_set(all)) THEN BEGIN
        IF file_test(dehnen_giants_filename) THEN BEGIN
            restore, dehnen_giants_filename
            IF keyword_set(tyc) THEN hip= tyc
        ENDIF ELSE BEGIN
            splog, 'Error: Dehnen-Giants-savefile does not exist in ', $
              dehnen_giants_filename
            splog, 'Dehnen-Giants-savefile must exist, returning...'
            RETURN
        ENDELSE
        hip_giants= hip
    ENDIF
    IF keyword_set(oldhip) THEN BEGIN
        dehnenfilename= '/global/data/hipparcos/hip-dehnen.sav'
        splog, 'Using the old Hipparcos reduction'
    ENDIF ELSE IF keyword_set(tyc) THEN BEGIN
        dehnenfilename= '~jb2777/streams/pro/rv/tyc-dehnen.sav'
        splog, 'Using the old Hipparcos reduction with Tycho proper motions'
    ENDIF ELSE BEGIN
        dehnenfilename= '/global/data/hipparcos/hipnewredux/hip2-dehnen.sav'
        splog, 'Using the new Hipparcos reduction'
    ENDELSE
    IF file_test(dehnenfilename) THEN BEGIN
        restore, dehnenfilename
        IF keyword_set(tyc) THEN hip= tyc
    ENDIF ELSE BEGIN
        splog, 'Error: Dehnen-savefile does not exist in ', $
          dehnenfilename
        splog, 'Dehnen-savefile must exist, returning...'
        RETURN
    ENDELSE
    nhip= n_elements(hip)
    sindx= sort(hip.bvcolor)
    IF keyword_set(colors) THEN BEGIN
        cindx= dblarr(nsubsamp+1)
        cindx[0]=0
        cindx[nsubsamp]= nhip-1
        FOR ii=0L, nsubsamp-2-keyword_set(giants) -2 * keyword_set(all) DO cindx[ii+1]= $
          max(where(hip[sindx].bvcolor LT colors[ii]))
    ENDIF

;;Set constraints to use in bestV
    IF keyword_set(bestV) THEN BEGIN
        IF ~keyword_set(Vconstraint) THEN $
          Vconstraints= [1D,3D,5D,7D,10D,12D,15D,20D,25D,40D,0D,0.1D,0.5D] $
        ELSE Vconstraints= Vconstraint
        nVconstraints= n_elements(Vconstraints)
        IF ~keyword_set(best_ngauss) THEN $
          Vlike= dblarr(nVconstraints,Vnjack+1) $
        ELSE $
          Vlike= dblarr(nVconstraints,Vnjack+1,endngauss[0]+1)
    ENDIF ELSE BEGIN
        Vconstraints= Vconstraint
        nVconstraints= 1L
    ENDELSE


;;Create arrays for the output
    fnlindx= nsubsamp*(~keyword_set(sample))+keyword_set(sample)
    color= dblarr(fnlindx)
    mean= dblarr(d,maxngauss,nVconstraints,fnlindx)
    covar= dblarr(d,d,maxngauss, nVconstraints,fnlindx)
    amp= dblarr(maxngauss,nVconstraints,fnlindx)
    IF keyword_set(best_ngauss) THEN BEGIN
        maxngauss= endngauss[0]
        mean= dblarr(d,maxngauss,nVconstraints,maxngauss,fnlindx)
        covar= dblarr(d,d,maxngauss, nVconstraints,maxngauss,fnlindx)
        amp= dblarr(maxngauss,nVconstraints,maxngauss,fnlindx)
    ENDIF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   MAIN LOOP:
;      loops over subsamples and runs the projected_gauss_mixtures for
;      all of them, also loops over jackknife subsamples
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    nstar= round(nhip/nsubsamp)
    FOR ii=0L, nsubsamp-1 DO WHILE ($;down LE jackconv AND $
                                    ngauss[ii] LE endngauss[ii]) DO BEGIN
        IF (keyword_set(giants) AND (ii EQ (nsubsamp-1-2*keyword_set(all)))) THEN BEGIN
            hip_main= hip
            hip= hip_giants
        ENDIF
        IF (keyword_set(all) AND (ii EQ nsubsamp-1)) THEN BEGIN
            hip= [hip_main, hip_giants]
        ENDIF
        IF (keyword_set(all) AND (ii EQ nsubsamp-2)) THEN BEGIN
            hip= hip_main
        ENDIF
        
        IF (sample NE (ii+1) AND sample NE 0) THEN BEGIN
            ngauss[ii]++
            CONTINUE
        ENDIF

            
;;Loop over Vconstraints samples and best_ngauss jackknife samples
        FOR jj=0L,nVconstraints*(Vnjack+1)-1 DO BEGIN
;;Load the sample
        IF keyword_set(colors) THEN BEGIN
            IF (ii LE (nsubsamp-1-keyword_set(giants) - keyword_set(all))) THEN BEGIN
                subindx=sindx[cindx[ii]:cindx[ii+1]]
            ENDIF ELSE BEGIN
                subindx= lindgen(n_elements(hip.hip))
            ENDELSE
        ENDIF ELSE BEGIN
            subindx= sindx[ii*nstar:((ii+1L)*nstar < (nhip - 1L))]
        ENDELSE
        nsub= n_elements(subindx)
           
        IF jj EQ 0 THEN BEGIN
            splog, 'Working on subsample '+strtrim(string(ii+1),2)+ $
              '/'+strtrim(string(nsubsamp),2)
            splog, 'Subsample has '+ strtrim(string(nsub),2) + ' stars'
            splog, 'Using '+strtrim(string(ngauss[ii]),2)+' Gaussians'
            IF keyword_set(log) THEN BEGIN
                PRINTF, loglun, 'Working on subsample '+strtrim(string(ii+1),2)+ $
                  '/'+strtrim(string(nsubsamp),2)
                PRINTF, loglun, 'Subsample has '+ strtrim(string(nsub),2) + ' stars'
                PRINTF, loglun, 'Using '+strtrim(string(ngauss[ii]),2)+' Gaussians'
                FLUSH, loglun
            ENDIF
        ENDIF
;;Colors
        crange=minmax(hip[subindx].bvcolor)
        color[ii*(~keyword_set(sample))]= (crange[1]+crange[0])/2D
        
;;Set fixmean and fixcovar
        fixmean= bytarr(ngauss[ii])+1*((jj MOD (Vnjack+1)) NE 0)
        fixcovar= bytarr(ngauss[ii])+1*((jj MOD (Vnjack+1)) NE 0)
;        fixmean[ngauss[ii]-1]= fixhalo
;        fixcovar[ngauss[ii]-1]= fixhalo
        
;;If we are doing best_ngauss, split the sample
;            IF (jj GT nVconstraints*Vnjack) THEN BEGIN
;                jackindx= shuffle_indx(nsub,seed=seed)
;                splitindx= ceil((1-jacksub)*nsub)+1
;                testsubindx= subindx[jackindx[splitindx:nsub-1]]
;                subindx=     subindx[jackindx[0:splitindx-1]]
;                nsub= n_elements(subindx)
;                splog, 'Working on (best_ngauss) jackknife-subsample '+strtrim(string(jj-nVconstraints*Vnjack),2)+'/'+$
;                  strtrim(string(njack),2)
;                IF keyword_set(log) THEN  PRINTF, loglun, 'Working on (best_ngauss) jackknife-subsample '+strtrim(string(jj-nVconstraints*Vnjack),2)+'/'+$
;                  strtrim(string(njack),2)
;                IF keyword_set(log) THEN FLUSH, loglun
;            ENDIF

        IF ((jj MOD (Vnjack+1)) EQ 0 AND keyword_set(bestV)) THEN BEGIN
            splog, 'Working on bestV, V= '+strtrim(string(Vconstraints[jj/(Vnjack+1)]),2)
            IF keyword_set(log) THEN PRINTF, loglun, 'Working on bestV, V= '+strtrim(string(Vconstraints[jj/(Vnjack+1)]),2)
            IF keyword_set(log) THEN FLUSH, loglun
        ENDIF
        
;;If we are determining the best variance constraint, split
        IF ((jj MOD (Vnjack+1)) EQ 1 AND keyword_set(realjack)) THEN rndmindx= shuffle_indx(nsub,seed=seed)
        IF ((jj MOD (Vnjack+1)) EQ 1 AND keyword_set(vanilla)) THEN rndmindx= shuffle_indx(nsub,seed=-1L)
                      
        
        IF ((jj MOD (Vnjack+1)) NE 0) THEN BEGIN
            IF (~keyword_set(realjack) AND ~keyword_set(vanilla)) THEN BEGIN
                jackindx= shuffle_indx(nsub,seed=seed)
                splitindx= ceil((1-Vjacksub)*nsub)+1
                testsubindx= subindx[jackindx[splitindx:nsub-1]]
                subindx=     subindx[jackindx[0:splitindx-1]]
            ENDIF ELSE BEGIN
                jjmod= (jj MOD (Vnjack+1))
                nrndm= n_elements(rndmindx)
                testsubindx= subindx[rndmindx[((jjmod-1)*nrndm)/Vnjack:(jjmod*nrndm)/Vnjack-1 < (nrndm-1)]]
                IF jjmod EQ 1 THEN BEGIN
                    subindx= subindx[rndmindx[nrndm/Vnjack:nrndm-1]]
                ENDIF ELSE IF jjmod EQ Vnjack THEN BEGIN
                    subindx= subindx[0:((Vnjack-1)*nrndm)/Vnjack-1]
                ENDIF ELSE BEGIN
                    subindx= [subindx[rndmindx[0:((jjmod-1)*nrndm)/Vnjack-1]],subindx[(jjmod*nrndm)/Vnjack:nrndm-1]]
                ENDELSE
            ENDELSE
            nsub= n_elements(subindx)
            splog, 'Working on (bestV) jackknife-subsample '+strtrim(string((jj MOD (Vnjack+1))),2)+'/'+$
              strtrim(string(Vnjack),2)
            IF keyword_set(log) THEN PRINTF, loglun, 'Working on (bestV) jackknife-subsample '+strtrim(string((jj MOD (Vnjack+1))),2)+'/'+$
              strtrim(string(Vnjack),2)
            IF keyword_set(log) THEN FLUSH, loglun
        ENDIF


;;If we are doing bestV, in a bestngauss determination, load the bestV
;;to do the Vconstraint
        IF keyword_set(Vconstraint) THEN BEGIN
            constraint= Vconstraint
        ENDIF ELSE IF (keyword_set(best_ngauss) AND ~keyword_set(bestV)) THEN BEGIN
            constraint= 0D
        ENDIF ELSE IF keyword_set(bestV) THEN BEGIN
            constraint= Vconstraints[jj/(Vnjack+1)]
            IF ((jj MOD (Vnjack+1)) EQ 0 AND ~keyword_set(bestV)) THEN BEGIN
                splog, 'Using Vconstraint = '+strtrim(string(constraint),2)
                IF keyword_set(log) THEN PRINTF, loglun, 'Using Vconstraint = '+strtrim(string(constraint),2)
                IF keyword_set(log) THEN FLUSH, loglun
            ENDIF
        ENDIF ELSE BEGIN
            constraint= 0D
        ENDELSE


;;Handle errors
        retry=3
        count=0
        CATCH, Error_Status
        IF Error_Status NE 0 THEN BEGIN
            IF keyword_set(debug) THEN PRINTF, -2, 'Error Index: ', Error_Status
;                bad_io1:
            IF keyword_set(debug) THEN PRINTF, -2, 'Error message: ', !ERROR_STATE.MSG
            IF keyword_set(debug) THEN PRINTF, -2, 'Trying again...'
            count=count+1
            IF count GT retry THEN BEGIN
                CATCH, /CANCEL
                IF keyword_set(debug) THEN PRINTF, -2, 'Last try...'
            ENDIF
        ENDIF
        
        
        initmean=dblarr(d,ngauss[ii])
        initcovar= dblarr(d,d,ngauss[ii])
        initamp= dblarr(ngauss[ii])
        IF ((jj MOD (Vnjack+1)) EQ 0) THEN BEGIN
;;Initial 'initial conditions': 1 for thin disk, 1 for thick disk, 1
;;for halo, rest, ~0, growing covariances
;                splog, 'Using initial initial conditions'
;;Thin disk
            initmean[*,0]= thinmean
            initcovar[*,*,0]= thincovar
            initamp[0]= 0.899D/double(ngauss[ii]-2)
            
;;Thick disk
            IF ngauss[ii] NE 1 THEN BEGIN
                initmean[*,1]= thickmean
                initcovar[*,*,1]= thickcovar
                initamp[1]= .1D
            ENDIF ELSE BEGIN
                initamp[0] = 1D
            ENDELSE
;;halo
            IF ngauss[ii] GT 2 THEN BEGIN
                initamp[2]= 0.001D
                initmean[*,2]= halomean
                initcovar[*,*,2]= halovar
            ENDIF ELSE IF ngauss[ii] EQ 2 THEN BEGIN
                initamp[0]= .9D
                initamp[1]= .1D
            ENDIF
;;Rest
            FOR gg=3L, ngauss[ii]-1 DO FOR kk=0,d-1 DO initmean[kk,gg]= $
              randomn(seed)
            FOR gg=3L, ngauss[ii]-1 DO FOR kk=0,d-1 DO initcovar[kk,kk,gg]= $
              1D2*2.0^double(gg < 30)*(gg LE 30) + 1D2*2.0^double(60-gg)*(gg GT 30)
            FOR gg=3L, ngauss[ii]-1 DO initamp[gg]=0.899D/double(ngauss[ii]-2)
        ENDIF ELSE BEGIN
;                splog, 'Using initial conditions from full sample result'
            initmean= currentmean
            initcovar= currentcovar
            initamp= currentamp
        ENDELSE

        diffloglike=2*itttol
        itt=0
        oldavgloglike=-((machar(/double)).xmax)
        didntgetbetter= 0
        WHILE (itt++ LT maxitt AND didntgetbetter LT maxstag) DO BEGIN ;diffloglike GE itttol AND 
            xamp= initamp
            xmean= initmean
            xcovar= initcovar
            IF (keyword_set(initial) AND (jj MOD (Vnjack+1)) EQ 0) THEN BEGIN
                splog, 'Iteration '+strtrim(string(itt),2)
                IF keyword_set(log) THEN PRINTF, loglun, 'Iteration '+strtrim(string(itt),2)
                IF keyword_set(log) THEN FLUSH, loglun
            ENDIF
            
; Projected gauss mixtures
;;Get data (without any radial velocities)
;            ydata= transpose([[hip[subindx].vl], [hip[subindx].vb]])
;            ycovar= hip[subindx].vlvbc
;            projection= hip[subindx].nsm
;;Get data
            ydata= transpose([[dblarr(nsub)],[hip[subindx].vl], [hip[subindx].vb]])
            ycovar_tmp= hip[subindx].vlvbc
            projection= hip[subindx].sm
;;Add radial velocities if those are given as well
            ycovar= dblarr(3,3,nsub)
            FOR gg= 0, nsub-1 DO BEGIN
                ycovar[*,*,gg]=[[10000.^2,0.,0.],[0.,ycovar_tmp[0,0,gg],ycovar_tmp[0,1,gg]],[0.,ycovar_tmp[1,0,gg],ycovar_tmp[1,1,gg]]]
            ENDFOR

;           stop
          
;;Set up logfile (only for the main convergence)
            IF ((jj MOD (Vnjack+1)) EQ 0) THEN BEGIN
                logfile=prefix
;                splitnmerge=
;                ngauss[ii]*(ngauss[ii]-1)*(ngauss[ii]-2)/2
                splitnmerge= ngauss[ii]*3
            ENDIF  ELSE BEGIN
                logfile=''
                splitnmerge=0
            ENDELSE
            projected_gauss_mixtures_c, ngauss[ii], ydata, ycovar, projection, $
              xamp, xmean, xcovar, fixmean=fixmean, fixcovar=fixcovar, $
              tol=tol, avgloglikedata=avgloglikedata, w=constraint, logfile=logfile, $
              splitnmerge=splitnmerge,/quiet

            CATCH, /CANCEL
            
;            stop
            
            IF (keyword_set(initial) AND (jj MOD (Vnjack+1) EQ 0)) THEN BEGIN
                diffloglike= avgloglikedata-oldavgloglike
                IF diffloglike LE 1D-7 THEN BEGIN
                    splog, 'loglikelihood did not decrease'
                    IF keyword_set(log) THEN PRINTF, loglun, 'loglikelihood did not decrease'
                    youwantmetobreak= 0L
                    wait, 5
                    splog, 'Too late!'
                    didntgetbetter++
                    diffloglike= 2*itttol
                    ;;Restore old initial conditions
                    xamp= initamp
                    xmean= initmean
                    xcovar= initcovar
                    IF youwantmetobreak THEN break
                ENDIF ELSE BEGIN
                    didntgetbetter=0
                    oldavgloglike= avgloglikedata
                ENDELSE
                
;;load new initial conditions
                quants= weighted_quantile(xamp,quant=[0.20D,0.7D])
                minindx= where(xamp LE quants[0])
                maxindx= where(xamp GT quants[0])
;;Fix good ones
                IF diffloglike GT 0 THEN BEGIN
                    fixamp= bytarr(ngauss[ii])
;                fixmean= bytarr(ngauss[ii])+1*(xamp GT quants[1])
;                fixcovar= bytarr(ngauss[ii])+1*(xamp GT quants[1])
                    fixmean= bytarr(ngauss[ii])+1*(xamp GT 1D0/double(ngauss[ii]))
                    fixcovar= bytarr(ngauss[ii])+1*(xamp GT 1D0/double(ngauss[ii]))
;                    fixamp[minindx]= 0
;                    fixmean[minindx]= 0
;                    fixcovar[minindx]= 0
                ENDIF
                minindx= where(xamp LE double(minamp)/nsub)
                IF minindx[0] EQ -1 THEN BEGIN
                    splog, 'All Gaussians have large enough amplitudes'
                    IF keyword_set(log) THEN PRINTF, loglun, 'All Gaussians have large enough amplitudes'
                    BREAK
                ENDIF
                maxindx= where(xamp GT quants[1])
                IF maxindx[0] EQ -1 THEN maxindx= (where(xamp EQ max(xamp)))[0]
                nmin= n_elements(minindx)
                nmax= n_elements(maxindx)
                minindx= minindx[floor(randomu(seed)*nmin)]
                maxindx= maxindx[floor(randomu(seed)*nmax)]
                initmean= xmean
                initcovar= xcovar
                initamp= xamp
                FOR dd=0,d-1 DO initmean[dd,minindx]= xmean[dd,maxindx]+randomn(seed)
                initcovar[*,*,minindx]= xcovar[*,*,maxindx]
;            FOR dd=0,d-1 DO initcovar[dd,dd,minindx]= xcovar[dd,dd,maxindx]*1D2
            ENDIF ELSE BEGIN
                BREAK
            ENDELSE

        ENDWHILE

        IF ((jj MOD (Vnjack + 1)) EQ 0) THEN BEGIN
            splog, 'Average loglike after convergence (full sample): '+ $
              strtrim(string(avgloglikedata),2)
            IF keyword_set(log) THEN PRINTF, loglun, 'Average loglike after convergence (full sample): '+ $
              strtrim(string(avgloglikedata),2)
            IF keyword_set(log) THEN FLUSH, loglun
        ENDIF
        
        IF (~keyword_set(best_ngauss) AND (jj MOD (Vnjack + 1)) EQ 0) THEN BEGIN
            mean[*,0:ngauss[ii]-1,jj/(Vnjack+1),ii*(~keyword_set(sample))]= xmean
            covar[*,*,0:ngauss[ii]-1,jj/(Vnjack+1),ii*(~keyword_set(sample))]= xcovar
            amp[0:ngauss[ii]-1,jj/(Vnjack+1),ii*(~keyword_set(sample))]= xamp
        ENDIF
        IF (keyword_set(best_ngauss) AND (jj MOD (Vnjack + 1)) EQ 0) THEN BEGIN
            mean[*,0:ngauss[ii]-1,jj/(Vnjack+1),ngauss[ii],ii]= xmean
            covar[*,*,0:ngauss[ii]-1,jj/(Vnjack+1),ngauss[ii],ii]= xcovar
            amp[0:ngauss[ii]-1,jj/(Vnjack+1),ngauss[ii],ii]= xamp
        ENDIF
        IF ((jj MOD (Vnjack + 1)) EQ 0) THEN BEGIN
            currentmean= xmean
            currentcovar= xcovar
            currentamp= xamp
        ENDIF

;;compute and save jackknife likelihoods
        IF ((jj MOD (Vnjack+1) NE 0) AND (keyword_set(best_ngauss) OR keyword_set(bestV))) THEN BEGIN
;                IF keyword_set(debug) THEN $
;                  splog, 'Computing the jackknife likelihood...'
            tydata= transpose([[hip[testsubindx].vl], $
                               [hip[testsubindx].vb]])
            tycovar= hip[testsubindx].vlvbc
            tprojection= hip[testsubindx].nsm
            avgloglikedata= 0D0
            projected_gauss_mixtures_C, ngauss[ii], tydata, tycovar, $
              tprojection, xamp, xmean, xcovar, /likeonly, $
              avgloglikedata=avgloglikedata, /quiet
;                IF (jj GT nVconstraints*Vnjack) THEN BEGIN
;                    jacklike[ii,ngauss[ii],jj-nVconstraints*Vnjack-1]= avgloglikedata
;                ENDIF ELSE 
            IF (jj MOD (Vnjack+1) NE 0) THEN BEGIN
                Vlike[jj/(Vnjack+1),(jj MOD (Vnjack+1))-1,keyword_set(best_ngauss)*ngauss[ii]]= avgloglikedata*double(n_elements(testsubindx))
            ENDIF
        ENDIF

;;Calculate total loglike of Vlike
        IF (keyword_set(bestV) AND ((jj+1) MOD (Vnjack+1)) EQ 0) THEN BEGIN
            Vlike[jj/(Vnjack+1),Vnjack,keyword_set(best_ngauss)*ngauss[ii]]=$
              total(Vlike[jj/(Vnjack+1),0:Vnjack-1,keyword_set(best_ngauss)*ngauss[ii]]);/double(Vnjack)
            splog, 'Total log-likelihood bestV: '+$
              strtrim(string(Vlike[jj/(Vnjack+1),Vnjack,keyword_set(best_ngauss)*ngauss[ii]]),2);+' +- '+$
;              strtrim(string(stddev(Vlike[jj/(Vnjack+1),0:Vnjack-1,keyword_set(best_ngauss)*ngauss[ii]])),2)
            IF keyword_set(log) THEN PRINTF, loglun, 'Average log-likelihood bestV: '+$
              strtrim(string(Vlike[jj/(Vnjack+1),Vnjack,keyword_set(best_ngauss)*ngauss[ii]]),2);+' +- '+$
;              strtrim(string(stddev(Vlike[jj/(Vnjack+1),0:Vnjack-1,keyword_set(best_ngauss)*ngauss[ii]])),2)
            IF keyword_set(log) THEN FLUSH, loglun
        ENDIF
        
        IF keyword_set(best_ngauss) THEN BEGIN
            IF keyword_set(bestV) THEN $
              save, filename=savefilename, ngauss, mean, amp, covar, Vlike $ ;, jacklike
            ELSE save, filename=savefilename, ngauss, mean, amp, covar ;, jacklike
        ENDIF ELSE IF keyword_set(bestV) THEN BEGIN
            save, filename=savefilename, color, mean, covar, amp, ngauss, Vlike
        ENDIF ELSE BEGIN
            save, filename=savefilename, color, mean, covar, amp, ngauss
        ENDELSE
        ENDFOR
;        IF keyword_set(best_ngauss) THEN BEGIN
;            jacklike[ii,ngauss[ii],njack]= $
;              total(jacklike[ii,ngauss[ii],*])/double(njack)
;            IF ngauss[ii] NE 2 THEN IF (jacklike[ii,ngauss[ii],njack] LT $
;                                        jacklike[ii,ngauss[ii]-1,njack]) THEN $
;              down++ $
;            ELSE down= 0
;            splog, 'Average (best_ngauss) log-likelihood = '+$
;              strtrim(string(jacklike[ii,ngauss[ii],njack]),2)
;            IF keyword_set(log) THEN PRINTF, loglun, 'Average (best_ngauss) log-likelihood = '+$
;              strtrim(string(jacklike[ii,ngauss[ii],njack]),2)
;            IF keyword_set(log) THEN FLUSH, loglun
;    ENDIF
        ngauss[ii]++
        
        IF keyword_set(debug) THEN BEGIN
            splog, 'Waiting 0.5 min...'
            WAIT, 30
        ENDIF

    ENDWHILE
    ngauss--

;    print, 'RESULT'
;    FOR ii=0L, nsubsamp-1 DO FOR gg=0L, ngauss-1 DO BEGIN
;        print, 'Results for Gaussian', gg, ' in subsample ', ii
;        print, 'mean = ', mean[*,gg,ii]
;        print, 'covar= ', covar[*,*,gg,ii]
;        print, 'amp= ', amp[gg,ii]
;    ENDFOR

    IF keyword_set(best_ngauss) THEN BEGIN
;        ngauss= lonarr(nsubsamp)
;        FOR ii=0L, nsubsamp DO BEGIN
;            maxlike= max(jacklike[ii,2:*,njack],max_indx)
;            ngauss[ii]= max_indx
;        ENDFOR
        IF keyword_set(bestV) THEN $
          save, filename=savefilename, ngauss, mean, amp, covar, Vlike $;, jacklike
        ELSE save, filename=savefilename, ngauss, mean, amp, covar;, jacklike
    ENDIF ELSE IF keyword_set(bestV) THEN BEGIN
        save, filename=savefilename, color, mean, covar, amp, ngauss, Vlike
    ENDIF ELSE BEGIN
        save, filename=savefilename, color, mean, covar, amp, ngauss
    ENDELSE
ENDELSE

IF keyword_set(log) THEN FREE_LUN, loglun

IF keyword_set(best_ngauss) THEN RETURN
IF keyword_set(justcalc) THEN RETURN

;Make plots of projected distributions and output parameters to LaTeX file
splog, 'Plotting the projected distributions and writing parameters to file...'
;Set-up LaTeX file
vunit= 'km s$^{-1}$'
v2unit= 'km$^2$ s$^{-2}$'
OPENW, wlun, prefix+'_parameters.tex', /GET_LUN
PRINTF, wlun, '\begin{deluxetable}{lccccccccccc}'
PRINTF, wlun, '\tablecaption{Parameters of the component '+$
  'Gaussians\label{tableparam}}'
PRINTF, wlun, '\tablecolumns{12}'
PRINTF, wlun, '\tablewidth{0pt}'
PRINTF, wlun, '\rotate'
PRINTF, wlun, '\tablehead{\colhead{Sample} & \colhead{\ngauss}'+$
  ' & \colhead{\alphagauss} & \colhead{$\eex\T\,\vv$} & '+$
  '\colhead{$\eey\T\,\vv$} & \colhead{$\eez\T\,\vv$} & '+$
  '\colhead{$\eex\T\,\VV\,\eex$} &\colhead{$\eey\T\,\VV\,\eey$} & '+$
  '\colhead{$\eez\T\,\VV\,\eez$} &\colhead{$\eex\T\,\VV\,\eey$} &'+$
  '\colhead{$\eex\T\,\VV\,\eez$} & \colhead{$\eey\T\,\VV\,\eez$} \\'
PRINTF, wlun, '\colhead{} & \colhead{} & \colhead{} & '+$
  '\colhead{('+ vunit+')} & \colhead{('+ vunit+')} & \colhead{('+ $
  vunit+')} & '+ '\colhead{('+ v2unit+')} & \colhead{('+ $
  v2unit+')} & \colhead{('+ v2unit+')} &'+ ' \colhead{('+ $
  v2unit+')} & \colhead{('+ v2unit+')} & \colhead{('+ v2unit+')}}'
PRINTF, wlun, '\startdata'
;loop over subsamples
parameterformat=strarr(10)
parameterformat[0:2]='(G8.4)'
parameterformat[3:8]='(G13.5)'
parameterformat[9]='(G8.4)'
FOR ii=0L, nsubsamp-1 DO BEGIN
    IF (sample NE (ii+1) AND sample NE 0) THEN CONTINUE
    label='B'+strtrim(string(ii+1),2)
    IF (keyword_set(giants) AND ii EQ (nsubsamp-1-2*keyword_set(all))) THEN BEGIN
        label= 'GI'
    ENDIF
    IF (keyword_set(all) and ii EQ nsubsamp-2) THEN BEGIN
        label= 'ALMS'
    ENDIF
    IF (keyword_set(all) and ii EQ nsubsamp-1) THEN BEGIN
        label= 'AL'
    ENDIF
;Plot
    fnlindx= ii*(~keyword_set(sample))
    plot_projected_gaussians_wrap, mean[*,0:ngauss[ii]-1,fnlindx], $
      covar[*,*,0:ngauss[ii]-1,fnlindx], $
      xrange=[-130,120], yrange=[-120,60], zrange=[-70,70], $
;      xrange=[-50,50], yrange=[-50,50], zrange=[-50,50], $
      amp=amp[0:ngauss[ii]-1,fnlindx], $
      basefilename=label,label=label, quiet=~debug, grid=256
    scatter_projected_gaussians_wrap, mean[*,0:ngauss[ii]-1,fnlindx], $
      covar[*,*,0:ngauss[ii]-1,fnlindx], $
      xrange=[-130,120], yrange=[-120,60], zrange=[-70,70], $
;      xrange=[-50,50], yrange=[-50,50], zrange=[-50,50], $
      amp=amp[0:ngauss[ii]-1,fnlindx], $
      basefilename=label,label=label
;Write parameters to file
    PRINTF, wlun, label+' & 1 & '
    PRINTF, wlun, amp[0,fnlindx],mean[0,0,fnlindx],mean[1,0,fnlindx],mean[2,0,fnlindx],$
      covar[0,0,0,fnlindx],covar[1,1,0,fnlindx],covar[2,2,0,fnlindx],covar[0,1,0,fnlindx],$
      covar[0,2,0,fnlindx],covar[1,2,0,fnlindx], format='('+parameterformat[9]+$
      '," & ", '+parameterformat[0]+'," & ",'+parameterformat[1]+'," & ",'+$
      parameterformat[2]+'," & ",'+parameterformat[3]+'," & ",'+$
      parameterformat[4]+'," & ",'+parameterformat[5]+'," & ",'+$
      parameterformat[6]+'," & ",'+parameterformat[7]+'," & ",'+$
      parameterformat[8]+',"\\")'
    FOR gg=1L, ngauss[ii]-1 DO PRINTF, wlun, gg+1, amp[gg,fnlindx],mean[0,gg,fnlindx],$
      mean[1,gg,fnlindx],mean[2,gg,fnlindx],$
      covar[0,0,gg,fnlindx],covar[1,1,gg,fnlindx],covar[2,2,gg,fnlindx],covar[0,1,gg,fnlindx],$
      covar[0,2,gg,fnlindx],covar[1,2,gg,fnlindx], format='( "&",'+'(I)'+ '," & " '+$
      parameterformat[9]+'," & ",'+parameterformat[0]+'," & ",'+$
      parameterformat[1]+'," & ",'+parameterformat[2]+'," & ",'+$
      parameterformat[3]+'," & ",'+parameterformat[4]+'," & ",'+$
      parameterformat[5]+'," & ",'+parameterformat[6]+'," & ",'+$
      parameterformat[7]+'," & ",'+parameterformat[8]+',"\\")'
    IF ii NE nsubsamp-1 THEN PRINTF, wlun, '& & & & & & & & & & & \\'
ENDFOR
PRINTF, wlun, '\enddata'
PRINTF, wlun, '\tablecomments{}'
PRINTF, wlun, '\end{deluxetable}'
FREE_LUN, wlun




END
