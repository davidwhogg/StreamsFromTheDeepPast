;+
;   NAME:
;      tile_highlowpredict
;   PURPOSE:
;      plot high and low entropy predictions of the GCS radial velocities from
;      the Hipparcos data
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      mplot - array for multiplot (e.g. [0,3,2])
;      filename - filename for plot
;      basedir - basedirectory for density estimate (ends in /)
;      savefilename - save the entropies here
;   KEYWORDS:
;      /colour - plot in color (in the density estimate there is a
;                 member color)
;      /QA     - Print QA info
;      /low    - plot bad predictions (default = high)
;   OUTPUT:
;      plots
;   EXAMPLE:
;       tile_highlowpredict, K=10, w=4, sample=6,mplot=[0,3,2],filename='predict_gcs_high.ps', /QA
;   REVISION HISTORY:
;      2009-01-15 - Written Bovy (NYU)
;-
PRO TILE_HIGHLOWPREDICT, K=K, w=w, sample=sample, mplot=mplot, $
                         filename=filename, basedir=basedir, $
                         colour=colour, QA=QA, low=low, $
                         savefilename=savefilename


;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'


IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;   Calculate entropies
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Restore GCS
    genevafilename= '/global/data/rv/gcs/gcs.sav'
    IF file_test(genevafilename) THEN BEGIN
        restore, genevafilename
    ENDIF ELSE BEGIN
        splog, 'Error: GCS-savefile does not exist in ', $
          genevafilename
        splog, 'GCS-savefile must exist, returning...'
        RETURN
    ENDELSE
    ngcs= n_elements(gcs.hip)
    ;ngcs= 10;;Use for test runs
    
    entropies= dblarr(ngcs)
    altentropies= dblarr(ngcs)
    
;;Restore density estimate
    densityfile=basedir+'fitv'
    densityfile+= strtrim(string(K),2)
    densityfile+='_bestV_V'
    densityfile+=strtrim(string(w,format='(f4.1)'),2)
    densityfile+= '_sample'
    densityfile+= strtrim(string(sample),2)
    densityfile+= '.sav'
    IF file_test(densityfile) THEN BEGIN
        IF keyword_set(QA) THEN splog, 'reading '+densityfile
        restore, densityfile
    ENDIF ELSE BEGIN
        splog, 'Density estimate not found in file '+densityfile
        splog, 'Returning...'
        RETURN
    ENDELSE
    
;;Actual calculation of the entropies
    FOR ii=0L, ngcs-1 DO BEGIN
        predictrv,amp,mean,covar,gcs[ii].sm,vr=gcs[ii].vr,$
          vlvb=[[gcs[ii].vl],[gcs[ii].vb]],cvlvb=gcs[ii].vlvbc,$
          cvr=gcs[ii].vrc,loglike=loglike,altloglike=altloglike, $
          entropy=entropy, altentropy=altentropy
        entropies[ii]= entropy
        altentropies[ii]= altentropy
    ENDFOR
    
;;Sort the likelihoods
    sort_indx= SORT(entropies);;In ascending order
    
    IF keyword_set(savefilename) THEN BEGIN
        splog, "Saving in "+savefilename
        save, filename=savefilename, entropies, altentropies,sort_indx, ngcs
    ENDIF
ENDELSE

IF ~keyword_set(mplot) THEN BEGIN
    IF keyword_set(QA) THEN splog, "Just calculated the entropies, not plotting"
    RETURN
ENDIF

ngood= mplot[1] * mplot[2]

;sort_indx2= SORT(altentropies)

;;Plot the high ones or the low ones
IF ~keyword_set(low) THEN $
  tile_predictgcs, K=K, w=w, sample=sample, mplot=mplot, $
  indices=sort_indx[ngcs-ngood:ngcs-1], filename=filename, $
  basedir=basedir, QA=QA, colour=colour, xsize=6.9, ysize=3.5, $
  yrange=[0,0.032], rrange=[-185,185] $
ELSE $
  tile_predictgcs, K=K, w=w, sample=sample, mplot=mplot, $
  indices=sort_indx[0:ngood-1], filename=filename, $
  basedir=basedir, QA=QA, colour=colour,xsize=6.9, ysize=3.5, $
  rrange=[-55,55], yrange=[0,0.445],offset=.0001D

END
