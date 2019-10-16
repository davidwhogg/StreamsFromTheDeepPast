;+
;   NAME:
;      plot_kVloglike
;
;   PURPOSE:
;      visual representation of (k,V) loglikelihood matrix which is
;      the result of the jackknife-cross validation
;
;   INPUTS:
;      filename   - name of savefile that holds k, V, Vlike
;      sample     - number of the sample
;
;   OUTPUTS:
;      plot of the matrix
;
;   REVISION HISTORY:
;      2008-08-01 - Written Bovy
;-
PRO PLOT_KVLOGLIKE, filename=filename, sample=sample

restore, filename

ks= n_elements(k)
Vs= n_elements(V)


krange= [k[0],k[ks-1]]
Vrange= [V[0], V[Vs-1]]

k_print, filename='kVlike_sample'+strtrim(string(sample),2)+'.ps'

!X.RANGE= krange
!Y.RANGE= Vrange
nticks= 5
xinterval= hogg_interval(!X.RANGE,nticks=nticks)
yinterval= hogg_interval(!Y.RANGE,nticks=nticks)


djs_plot, [0], [1], xtickinterval=xinterval, ytickinterval=yinterval

;;Make image
tvrange= minmax(Vlike)

;tvrange= [tvrange[0]-(tvrange[1]-tvrange[0]),tvrange[1]]
;tvrange= [tvrange[0]-5,tvrange[1]]
tvrange= [-3540D,tvrange[1]]
Vlike[where(Vlike LE -3530D)]= -3530D

nfindx= where(~FINITE(Vlike))
IF nfindx[0] NE -1 THEN Vlike[nfindx]= tvrange[0]


;;Resample
tvimage= dblarr(max(krange),ceil(max(Vrange)/.25))
nk= max(krange)
nV= ceil(max(Vrange)/.25)
kitt=0
FOR kk=0L, nk-1 DO BEGIN
    IF kitt NE ks-1 THEN IF (kk GT k[kitt]) THEN kitt++
    Vitt=0
    FOR vv=0L, nV-1 DO BEGIN
        IF Vitt NE Vs-1 THEN IF (VV*.25 GT V[Vitt]) THEN Vitt++
        tvimage[kk,VV]= Vlike[kitt,Vitt]
    ENDFOR
ENDFOR

tvimage= 255-((floor(256*(tvimage-tvrange[0])/ $
                         (tvrange[1]-tvrange[0])) > 0) $
                  < 255)
    tv, tvimage, !X.RANGE[0],!Y.RANGE[0],/data, $
      xsize=(!X.CRANGE[1]-!X.CRANGE[0]), $
      ysize=(!Y.CRANGE[1]-!Y.CRANGE[0])

;;Remake axes, and put axis-labels
    !P.MULTI[0]= !P.MULTI[0]+1
    djs_plot, [0],[1],xtickinterval=xinterval,ytickinterval=yinterval, $
      xtitle='k', ytitle='V'

k_end_print

END
