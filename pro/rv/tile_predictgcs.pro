;+
;   NAME:
;      tile_predictgcs
;   PURPOSE:
;      predict the radial velocities from the GCS survey, output a
;      plot which holds some of them
;   CALLING SEQUENCE:
;      tile_predictgcs, K=K, w=w, sample=sample, 
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      mplot - array for multiplot (e.g. [0,4,4])
;      indices - array which holds the indices in the gcs saveset of
;                the vrs to predict and output
;      filename - filename to plot to
;      seed    - seed for random
;      basedir - basedirectory for density estimate (ends in /)
;      yrange  - prob range
;      rrange  - vr range
;      xsize   - x size for k_print (can be left unset)
;      ysize   - ysize for k_print
;      offset  - tiny offset to add to vr such that it is
;                plotted when it is zero
;   KEYWORDS:
;      /random - pick random gcs observations
;      /colour - plot in color (in the density estimate there is a
;                 member color)
;      /QA     - Print QA info
;      hipnotgcs - don't use the gcs sample, instead use hipparcos
;                  stars that aren't in the GCS catalog
;   OUTPUT:
;      plot
;   EXAMPLE:
;       tile_predictgcs, K=10, w=4, sample=6, mplot=[0,3,4],/random, filename='predict_gcs_ex.ps'
;   REVISION HISTORY:
;      2008-10-27 - Written Bovy
;-
PRO TILE_PREDICTGCS, K=K, w=w,sample=sample, mplot=mplot, offset=offset,$
                     indices=indices, filename=filename, random=random, $
                     seed=seed, colour=colour, basedir=basedir, QA=QA, $
                     yrange=yrange, rrange=rrange, xsize=xsize, ysize=ysize,$
                     hipnotgcs=hipnotgcs

;;Default
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(yrange) THEN yrange=[0,0.055]
IF ~keyword_set(rrange) THEN rrange=[-90,90]
IF ~keyword_set(offset) THEN offset= 0D


;;Create the string to restore the density distribution
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

IF keyword_set(hipnotgcs) THEN BEGIN
    ;;Restore Hipparcos
    dehnenfilename= '/global/data/hipparcos/hipnewredux/hip2-dehnen.sav'
    IF file_test(dehnenfilename) THEN BEGIN
        restore, dehnenfilename
    ENDIF ELSE BEGIN
        splog, 'Data file not found in '+dehnenfilename
        splog, 'Data file must exist!'
        splog, 'Returning...'
        RETURN
    ENDELSE
    nhip= n_elements(hip)
        
    ;;Find the hipparcos stars that aren't in the GCS sample
    ingcs= bytarr(nhip)
    FOR ii= 0L, nhip-1 DO BEGIN
        gcsmember= where(gcs.hip EQ hip[ii].hip)
        IF gcsmember[0] NE -1 THEN ingcs[ii]= 1
    ENDFOR
    notingcs= where(ingcs EQ 0)
    n_notingcs= n_elements(notingcs)
    
    ;;Then swith around the structures
    gcs= hip[notingcs]
    ;;And hope this works
ENDIF

;;check inputs
ntile= mplot[1] * mplot[2]
IF ~keyword_set(random) THEN IF ~keyword_set(indices) THEN BEGIN
    splog, 'Either "/random" or "indices" must be set'
    splog, 'Returning...'
    RETURN
ENDIF ELSE BEGIN
    IF (ntile NE n_elements(indices)) THEN BEGIN
        splog, 'Number of indices does not corresponds to the number of plots requested in mplot'
        splog, 'Returning...'
        RETURN
    ENDIF
ENDELSE

IF keyword_set(filename) THEN k_print, filename=filename, $
  xsize=xsize, ysize=ysize


;!p.multi=mplot
!p.multi=[0,3,3]

;;If random is set, create a random set of "indices"
IF keyword_set(random) THEN BEGIN
    IF ~keyword_set(seed) THEN seed=1L
    indices= lonarr(ntile)
    FOR ii=0L, ntile-1 DO BEGIN
        indices[ii]= floor(n_elements(gcs.hip) * randomu(seed))
    ENDFOR
ENDIF

;;Run predictrv

indexindex= 0L

FOR ii=0L, mplot[2]-1 DO BEGIN
    FOR jj=0L, mplot[1]-1 DO BEGIN
        IF keyword_set(hipnotgcs) THEN vr= 0D ELSE $
          vr= gcs[indices[indexindex]].vr+offset
        IF keyword_set(hipnotgcs) THEN cvr= 0D ELSE $
          cvr= gcs[indices[indexindex]].vrc
        IF jj EQ 0 THEN BEGIN
            IF ii EQ mplot[2]-1 THEN BEGIN
                predictrv, amp, mean, covar, gcs[indices[indexindex]].sm,$
                  vr=vr, $
                  vlvb=[[gcs[indices[indexindex]].vl],$
                        [gcs[indices[indexindex]].vb]], $
                  cvlvb=gcs[indices[indexindex]].vlvbc,$
                  cvr=cvr,$
                  loglike=loglike,altloglike=altloglike, $
                  /plot, rrange=rrange,  yrange=yrange, color=colour,$
                  identify='HIP'+$
                  strtrim(string(gcs[indices[indexindex]].hip),2), $
                  /putentropylabels
            ENDIF ELSE BEGIN
                predictrv, amp, mean, covar, gcs[indices[indexindex]].sm,$
                  vr=vr, $
                  vlvb=[[gcs[indices[indexindex]].vl],$
                        [gcs[indices[indexindex]].vb]], $
                  cvlvb=gcs[indices[indexindex]].vlvbc,$
                  cvr=cvr,$
                  loglike=loglike,altloglike=altloglike, $
                  /plot, rrange=rrange, yrange=yrange, /noxlabels, $
                  color=colour, $
                  identify='HIP'+$
                  strtrim(string(gcs[indices[indexindex]].hip),2), $
                  /putentropylabels
              ENDELSE
          ENDIF ELSE BEGIN
            IF ii EQ mplot[2]-1 THEN BEGIN
                predictrv, amp, mean, covar, gcs[indices[indexindex]].sm,$
                  vr=vr, $
                  vlvb=[[gcs[indices[indexindex]].vl],$
                        [gcs[indices[indexindex]].vb]], $
                  cvlvb=gcs[indices[indexindex]].vlvbc,$
                  cvr=cvr,$
                  loglike=loglike,altloglike=altloglike, $
                  /plot, rrange=rrange, /noylabels, yrange=yrange, $
                  color=colour, $
                  identify='HIP'+$
                  strtrim(string(gcs[indices[indexindex]].hip),2), $
                  /putentropylabels
            ENDIF ELSE BEGIN
                predictrv, amp, mean, covar, gcs[indices[indexindex]].sm,$
                  vr=vr, $
                  vlvb=[[gcs[indices[indexindex]].vl],$
                        [gcs[indices[indexindex]].vb]], $
                  cvlvb=gcs[indices[indexindex]].vlvbc,$
                  cvr=cvr,$
                  loglike=loglike,altloglike=altloglike, $
                  /plot, rrange=rrange,  /noxlabels, /noylabels, $
                  yrange=yrange, color=colour, $
                  identify='HIP'+$
                  strtrim(string(gcs[indices[indexindex]].hip),2), $
                  /putentropylabels
            ENDELSE
        ENDELSE
        IF keyword_set(QA) AND ~keyword_set(hipnotgcs) THEN BEGIN
            splog, 'For '+'HIP'+$
              strtrim(string(gcs[indices[indexindex]].hip),2)
            splog, 'v_r = '+ strtrim(string(gcs[indices[indexindex]].vr,$
                                            format='(F6.1)'),2)
            splog, 'log likelihood= '+ $
              strtrim(string(loglike,format='(F7.2)'),2)
        ENDIF
        indexindex += 1
    ENDFOR
ENDFOR

IF keyword_set(filename) THEN k_end_print

END
