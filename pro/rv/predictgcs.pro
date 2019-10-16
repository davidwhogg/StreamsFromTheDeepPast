;+
;   NAME:
;      predictgcs
;   PURPOSE:
;      predict the radial velocities from the GCS survey
;   CALLING SEQUENCE:
;      predictgcs, K=K, w=w, sample=sample
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      ii    - catalog index to use (FOR NOW)
;   OPTIONAL INPUT:
;      basedir - base directory for density estimate file (ends in /)
;   OUTPUT:
;      total log likelihood
;   REVISION HISTORY:
;      2008-10-26 - Written Bovy
;-
PRO PREDICTGCS, K=K, w=w,sample=sample, ii=ii, basedir=basedir

;;Default
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMS/'

;;Create the string to restore the density distribution
densityfile=basedir+'fitv'
densityfile+= strtrim(string(K),2)
densityfile+='_bestV_V'
densityfile+=strtrim(string(w,format='(f4.1)'),2)
densityfile+= '_sample'
densityfile+= strtrim(string(sample),2)
densityfile+= '.sav'
IF file_test(densityfile) THEN BEGIN
    splog, 'reading '+densityfile
    restore, densityfile
ENDIF ELSE BEGIN
    splog, 'Density estimate not found in file '+densityfile
    splog, 'Returning...'
    RETURN
ENDELSE

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



;;Run predictrv

predictrv, amp, mean, covar, gcs[ii].sm,gcs[ii].vr, $
  vlvb=[[gcs[ii].vl],[gcs[ii].vb]], cvlvb=gcs[ii].vlvbc,$
  cvr=gcs[ii].vrc,loglike=loglike,altloglike=altloglike, $
  /plot,rrange=[-100,100], $
  identify='HIP'+strtrim(string(gcs[ii].hip),2), /putentropylabels

splog, loglike, altloglike



END
