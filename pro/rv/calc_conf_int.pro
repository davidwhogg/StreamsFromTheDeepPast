;+
;   NAME:
;      calc_conf_int
;   PURPOSE:
;      calculate a confidence interval for a distribution
;   CALLING SEQUENCE:
;   INPUT:
;      ds         - distribution (vector of values)
;      confidence - e.g. .95
;   KEYWORDS:
;      single     - compute a single confidence interval that is the
;                   narrowst possible (DEFAULT)
;      symm       - calculate a single confidence interval that has
;                   (1-confidence)/2 to the left
;      multimodal - calculate a confidence interval such that each
;                   point in it has a higher probability density than
;                   the points outside the interval. For multimodal
;                   distributions this means that there will generally
;                   be more than a single interval. All upper and
;                   lower limts are returned in the cul and cll.
;   OUTPUT:
;      cul        - upper limit(s) (index in the array)
;      cll        - lower limit(s) (index in the array)
;      confintwith- width of the computed confidence
;                   interval (in unit step size)
;   REVISION HISTORY:
;      2008-09-04 - Written Bovy (NYU)
;      2008-11-15 - Added symm keyword - Bovy
;-
PRO CALC_CONF_INT, ds, confidence, cul=cul, cll=cll, confintwidth=confintwidth, $
                   symm=symm, multimodal=multimodal, single=single
;;Search for highest likelihood interval that holds confidence% of the
;;probability
IF keyword_set(multimodal) THEN BEGIN
    cumindx= reverse(sort(ds))
    grid= n_elements(ds)
    cumds= dblarr(grid)
    cumds[cumindx]= total(ds[cumindx],/cumulative)
    cumds= cumds/max(cumds)
    confindx= where(cumds LE confidence)
    nconf= n_elements(confindx)
    cul= bytarr(nconf)
    cll= bytarr(nconf)
    IF (cumds[0] LE confidence) THEN BEGIN
        print, cumds[0]
        splog, 'ERROR: interval is not wide enough, increase size and run again'
        return
    ENDIF
    cll[0]= 1B
    cul[nconf-1]= 1B
    up= 1B
    FOR ii= 1, nconf-1 DO BEGIN
        IF confindx[ii] NE (confindx[ii-1]+1) THEN BEGIN
            cul[ii-1]= 1B
            cll[ii]= 1B
        ENDIF
    ENDFOR
    cul= confindx[where(cul EQ 1B)]
    cll= confindx[where(CLL EQ 1B)]
    IF arg_present(confintwidth) THEN confintwidth= nconf
ENDIF ELSE IF keyword_set(symm) THEN BEGIN
;;Search for a 'symmetrical' confidence interval (see header)
    cumds= total(ds,/cumulative)
    cumds= cumds/max(cumds)
    conf2= (1D0-confidence)/2D0
    FOR ii=1L, n_elements(ds)-1 DO BEGIN
        IF cumds[ii-1] LT conf2 AND cumds[ii] GE conf2 THEN cul= ii
        IF cumds[ii-1] LT (1-conf2) AND cumds[ii] GE (1-conf2) Then cll= ii
    ENDFOR
    IF arg_present(confintwidth) THEN confintwidth= cul-cll
ENDIF ELSE BEGIN
;;Search for minimal interval that holds confidence% probability
;;For every point, look ahead to where confidence% is reached
    cds= total(ds,/cumulative,/double)/total(ds,/double)
    candidates= where(cds LE 1D0-confidence)
    ncand= n_elements(candidates)
    pairs= lonarr(3,ncand)
    FOR ii=0, ncand-1 DO BEGIN
        candidates2= where(cds GT cds[ii])
        ncand2= n_elements(candidates2)
        tmpcds= cds[candidates2]- cds[candidates[ii]]
        pairs[0,ii]= ii
        pairs[1,ii]= min(where(tmpcds GE confidence))+ii
        pairs[2,ii]= pairs[1,ii]-pairs[0,ii]
    ENDFOR
    mindist= min(pairs[2,*],minii)
    cll= pairs[0,minii]
    cul= pairs[1,minii]
;;Return the interval?
    IF arg_present(confintwidth) THEN confintwidth= cul-cll
ENDELSE
END
