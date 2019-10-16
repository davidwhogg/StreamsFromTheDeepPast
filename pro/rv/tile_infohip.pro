;+
;   NAME:
;      tile_infohip
;   PURPOSE:
;      plot predictions for Hipparcos stars that aren't in GCS that
;      are particularly informative (high/low entropy)
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      mplot - array for multiplot (e.g. [0,3,2])
;      filename - filename for plot
;      basedir - basedirectory for density estimate (ends in /)
;      savefilename - save the entropies here
;      ntable  - number of records per table
;   KEYWORDS:
;      /colour - plot in color (in the density estimate there is a
;                 member color)
;      /QA     - Print QA info
;      /low    - plot bad predictions (default = high)
;      /table  - print table with high entropy predictions to filename
;      /alttable - print table with low entropy predictions to filename
;   OUTPUT:
;      plots/tables
;   EXAMPLE:
;       tile_infohip, K=10, w=4,
;       sample=6,mplot=[0,3,2],filename='info_hip_high.ps', /QA, savefilename='hip_entropies.sav'
;   REVISION HISTORY:
;      2009-01-18 - Written Bovy (NYU)
;-
PRO TILE_INFOHIP, K=K, w=w, sample=sample, mplot=mplot, filename=filename, $
                  basedir=basedir, savefilename=savefilename, $
                  colour=colour, QA=QA, low=low, table=table, $
                  alttable=alttable, ntable=ntable

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(ntable) THEN ntable= 25L

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

;;Restore Hipparcos
dehnenfilename= '/global/data/hipparcos/hipnewredux/hip2-dehnen.sav'
IF keyword_set(QA) THEN splog, 'QA: Using the new Hipparcos reduction'
IF file_test(dehnenfilename) THEN BEGIN
    restore, dehnenfilename
ENDIF ELSE BEGIN
    splog, 'Data file not found in '+dehnenfilename
    splog, 'Data file must exist!'
    splog, 'Returning...'
    RETURN
ENDELSE
nhip= n_elements(hip)
IF keyword_set(QA) THEN splog, 'QA: Data sample has '+$
  strtrim(string(nhip),2)+' members'

;;Find the hipparcos stars that aren't in the GCS sample
ingcs= bytarr(nhip)
FOR ii= 0L, nhip-1 DO BEGIN
    gcsmember= where(gcs.hip EQ hip[ii].hip)
    IF gcsmember[0] NE -1 THEN ingcs[ii]= 1
ENDFOR
notingcs= where(ingcs EQ 0)
n_notingcs= n_elements(notingcs)
IF keyword_set(QA) THEN splog, "Found "+strtrim(string(n_notingcs),2)+$
  " not in GCS"

;;Testing
;n_notingcs= 10



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



IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;   Calculate entropies
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;Now calculate the predictions + their entropies
    entropies= dblarr(n_notingcs)
    altentropies= dblarr(n_notingcs)
    
    ;;Actual calculation of the entropies
    FOR ii=0L, n_notingcs-1 DO BEGIN
        predictrv,amp,mean,covar,hip[notingcs[ii]].sm,$
          vlvb=[[hip[notingcs[ii]].vl],[hip[notingcs[ii]].vb]],$
          cvlvb=hip[notingcs[ii]].vlvbc,$
          entropy=entropy, altentropy=altentropy
        entropies[ii]= entropy
        altentropies[ii]= altentropy
    ENDFOR
    
;;Sort the likelihoods
    sort_indx= SORT(entropies);;In ascending order
    
    IF keyword_set(savefilename) THEN BEGIN
        splog, "Saving in "+savefilename
        save, filename=savefilename, entropies, altentropies,$
          sort_indx, n_notingcs
    ENDIF
ENDELSE


;;Now plot or print what's relevant
IF keyword_set(table) OR keyword_set(alttable) THEN BEGIN
    IF keyword_set(table) THEN BEGIN
        indices= lindgen(ntable)+n_notingcs-ntable
        indices= reverse(indices)
;        tablehip=
;        hip[notingcs[sort_indx[n_notingcs-ntable:n_notingcs-1]]]
        tablehip= hip[notingcs[sort_indx[indices]]]
        rrange= [-250,250]
    ENDIF ELSE BEGIN
        tablehip= hip[notingcs[sort_indx[0:ntable-1]]]
        rrange= [-250,250]
    ENDELSE
    IF keyword_set(alttable) THEN cntr= 'two' ELSE cntr= 'one'
    OPENW, wlun, filename, /GET_LUN
    PRINTF, wlun, '\clearpage'
    PRINTF, wlun, '\begin{deluxetable}{lrrrrrrrrrc}'
    IF keyword_set(table) THEN PRINTF, wlun, "\tablecaption{High entropy predictions for Hipparcos stars' radial velocities\label{table:hip_high}}" $
      ELSE PRINTF, wlun, "\tablecaption{Low entropy predictions for Hipparcos stars' radial velocities\label{table:hip_low}}"
    PRINTF, wlun, '\tablecolumns{11}'
    PRINTF, wlun, '\tablewidth{0pt}'
    PRINTF, wlun, '\tabletypesize{\small}'
;    PRINTF, wlun, '\rotate'
    PRINTF, wlun, '\tablehead{\colhead{ID\tablenotemark{\protect{\ref{HIPNUMBER'+cntr+'}}}} & \multicolumn{3}{c}{\ra\ (1991.25)}& \multicolumn{3}{c}{\dec\ (1991.25)}{} &\colhead{$v_r$\tablenotemark{\protect{\ref{VR'+cntr+'}}}}&\colhead{$v_r > (95\%)$\tablenotemark{\protect{\ref{VRLOWER'+cntr+'}}}}& \colhead{$v_r < (95\%)$\tablenotemark{\protect{\ref{VRUPPER'+cntr+'}}}} & \colhead{\phantom{bbb}$H$\tablenotemark{\protect{\ref{ENTROPY'+cntr+'}}}\phantom{bbb}}\\'
    PRINTF, wlun, '\colhead{} & \colhead{(h} & \colhead{m} & \colhead{s)} & \colhead{($\degree$} & \colhead{$\prime$} & \colhead{$\prime\prime$)} & \colhead{(km s$^{-1}$)} & \colhead{(km s$^{-1}$)} & \colhead{(km s$^{-1}$)} & \colhead{}}'
    PRINTF, wlun, '\startdata'
    FOR ii=0L, ntable-1 DO BEGIN
        predictrv,amp,mean,covar,tablehip[ii].sm,$
          vlvb=[[tablehip[ii].vl],[tablehip[ii].vb]],$
          cvlvb=tablehip[ii].vlvbc,rrange=rrange,$
          entropy=entropy, maxvr=maxvr, cul=cul, cll=cll
        ;;Write parameters to file
        ;;Dot-padding, to 8
        nhipnumber= alog10(tablehip[ii].hip)
        dotpadstr= ''
        FOR jj= 0L, 10-nhipnumber-1 DO dotpadstr+='..'
        IF bovy_sign(tablehip[ii].dec) EQ -1 THEN sgn_str= '-' ELSE sgn_str=''
        PRINTF, wlun, 'HIP'+strtrim(string(tablehip[ii].hip),2) +dotpadstr+'&' +$
          strtrim(string(bovy_ra_to_hms(tablehip[ii].ra),format='(I02)'),2)+'&'+$
          strtrim(string(bovy_ra_to_hms(tablehip[ii].ra,/min),format='(I02)'),2)+'&'+$
          strtrim(string(bovy_ra_to_hms(tablehip[ii].ra,/sec),format='(F05.2)'),2)+'&'+$
          sgn_str+$
          strtrim(string(bovy_dec_to_damas(abs(tablehip[ii].dec)),format='(I02)'),2)+'&'+$
          strtrim(string(bovy_dec_to_damas(tablehip[ii].dec,/arcmin),format='(I02)'),2)+'&'+$
          strtrim(string(bovy_dec_to_damas(tablehip[ii].dec,/arcsec),format='(F04.1)'),2)+'&'+$
          strtrim(string(maxvr,format='(I)'),2)+'&'+$
          strtrim(string(cll,format='(I)'),2)+'&'+$
          strtrim(string(cul,format='(I)'),2)+'&'+$
          strtrim(string(entropy,format='(F5.3)'),2)+'\\'
    ENDFOR
    PRINTF, wlun, '\enddata'
    PRINTF, wlun, '\tablecomments{Predictions for the radial velocity are based on the velocity density estimate using '+strtrim(string(K),2)+' Gaussians with regularization parameter $w$ = '+strtrim(string(long(w)),2)+" km$^2$ s$^{-2}$, marginalized for the position on the sky of the star and conditioned on the stars' tangential velocity.}"
    IF keyword_set(alttable) THEN cntr= 'two' ELSE cntr= 'one'
    PRINTF, wlun, '\setcounter{table'+cntr+'}{1}'
    PRINTF, wlun, '\makeatletter'
    PRINTF, wlun, '\let\@currentlabel\oldlabel'
    PRINTF, wlun, '\newcommand{\@currentlabel}{\thetable'+cntr+'}'
    PRINTF, wlun, '\makeatother'
    PRINTF, wlun, '\renewcommand{\thetable'+cntr+'}{\alph{table'+cntr+'}}'
    PRINTF, wlun, '\tablenotetext{\thetable'+cntr+'}{\label{HIPNUMBER'+cntr+'}'
    PRINTF, wlun, '\Hipparcos\ number.\stepcounter{table'+cntr+'}}'
    PRINTF, wlun, '\tablenotetext{\thetable'+cntr+'}{\label{VR'+cntr+'}'
    PRINTF, wlun, 'Maximum likelihood estimate of the radial velocity.\stepcounter{table'+cntr+'}}'
    PRINTF, wlun, '\tablenotetext{\thetable'+cntr+'}{\label{VRLOWER'+cntr+'}'
    PRINTF, wlun, 'Lower limit on the radial velocity (95-percent confidence).\stepcounter{table'+cntr+'}}'
    PRINTF, wlun, '\tablenotetext{\thetable'+cntr+'}{\label{VRUPPER'+cntr+'}'
    PRINTF, wlun, 'Upper limit on the radial velocity (95-percent confidence).\stepcounter{table'+cntr+'}}'
    PRINTF, wlun, '\tablenotetext{\thetable'+cntr+'}{\label{ENTROPY'+cntr+'}'
    PRINTF, wlun, 'Entropy of the probability distribution of the radial velocity.\stepcounter{table'+cntr+'}}'
    PRINTF, wlun, '\end{deluxetable}'
    FREE_LUN, wlun
ENDIF ELSE IF ~keyword_set(low) THEN BEGIN
    IF ~keyword_set(mplot) THEN BEGIN
        IF keyword_set(QA) THEN splog, "Just calculated entropes, not plotting"
        RETURN
    ENDIF
    ngood= mplot[1] * mplot[2] 
    tile_predictgcs, K=K, w=w, sample=sample, mplot=mplot, $ 
      indices=sort_indx[n_notingcs-ngood:n_notingcs-1], filename=filename, $
      basedir=basedir, QA=QA, colour=colour, xsize=6.9, ysize=3.5, $
      yrange=[0,0.032], rrange=[-230,230], /hipnotgcs
ENDIF ELSE BEGIN
    IF ~keyword_set(mplot) THEN BEGIN
        IF keyword_set(QA) THEN splog, "Just calculated entropes, not plotting"
        RETURN
    ENDIF
    ngood= mplot[1] * mplot[2]
    tile_predictgcs, K=K, w=w, sample=sample, mplot=mplot, $
      indices=sort_indx[0:ngood-1], filename=filename, $
      basedir=basedir, QA=QA, colour=colour,xsize=6.9, ysize=3.5, $
      rrange=[-80,80], yrange=[0,0.445], /hipnotgcs    
ENDELSE

END
