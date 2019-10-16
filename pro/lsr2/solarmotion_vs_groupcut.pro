;+
;   NAME:
;      solarmotion_vs_groupcut
;   PURPOSE:
;      plot the derived value of the Solar motion and its uncertainty
;      as a function of the groupcut employed
;   CALLING SEQUENCE:
;   INPUT: (most of these go directly to the 'lsr' function)
;      nboot            - number of bootstrap trials to use for error
;                         estimation
;      nsubsamp         - number of color subsamples to use
;      tol              - tolerance for proj_gauss_mixtures
;      halomean         - default [0D0,-2.2D2,0D0]
;      bighalo          - use a halo with a large vdisp (150 km/s)
;      groupcut         - [min,max,step]
;      plotfilename
;   KEYWORDS:
;      remove_ngc1901   - If set, consider the Gaussian centered o
;                         NGC1901 a moving group
;      dropblue         - drop the bluest stars when fitting for the
;                         lsr
;      S500             - If set, drop all samples with S < 500 from
;                         the fit
;      ex1000           - If set, drop all samples with S < 1000 from
;                         the fit
;      formal           - use Hipparcos formal errors
;      dontplot             
;   OUTPUT:
;   REVISION HISTORY:
;      2009-08-04 - Written - Bovy
;-
PRO SOLARMOTION_VS_GROUPCUT, remove_ngc1901=remove_ngc1901,groupcut=groupcut, $
                             nboot=nboot, nsubsamp=nsubsamp,tol=tol, $
                             dropblue=dropblue, prefix=prefix, S500=S500, formal=formal,$
                             plotfilename=plotfilename, savefilename=savefilename, $
                             ex1000=ex1000

restore, filename='hip2-aumer.sav'
nhip= n_elements(hip.hip)

IF keyword_set(savefilename) and file_test(savefilename) THEN BEGIN
    restore, filename=savefilename
ENDIF ELSE BEGIN
    ngroupcuts= (groupcut[1]-groupcut[0])/groupcut[2]+1
    UVWs= dblarr(ngroupcuts,7)
    groupcuts= dindgen(ngroupcuts)*groupcut[2]+groupcut[0]
    
    ;;Calculate all of the Solar motions
    FOR ii=0L, ngroupcuts-1 DO BEGIN
        splog, 'Working on groupcut '+strtrim(string(ii),2)+': '+strtrim(string(groupcuts[ii]),2)
        lsr, remove_ngc1901=remove_ngc1901,groupcut=groupcuts[ii], $
          nboot=nboot, nsubsamp=nsubsamp,tol=tol, seed=-1L, $
          dropblue=dropblue, S500=S500, formal=formal, $
          out=out, /dontplot, /silent, ex1000=ex1000
        UVWs[ii,*]= out
    ENDFOR
    IF keyword_set(savefilename) THEN BEGIN
        save, filename=savefilename, UVWs, groupcuts, groupcut
    ENDIF
ENDELSE


;;plotting symbol
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill
;;Now plot this
k_print, filename=plotfilename, xsize= 3.375, ysize= 5.0
!P.CHARTHICK= 2
!X.CHARSIZE= 1
!Y.CHARSIZE= 1
!X.THICK=2
!Y.THICK=2
ploterror, groupcuts, -UVWs[*,0], UVWs[*,1], psym=8, $
  position=[.1,.7,1.,1.], symsize=0.5, xtickformat='(A1)', $
  ytitle= textoidl('v_{\odot,x} [km s^{-1}]'), xrange=[groupcut[0]-.01,groupcut[1]+.01]
oplotbarx, 0.4, linestyle=1
ploterror, groupcuts, -UVWs[*,2], UVWs[*,3], psym=8, $
  position=[.1,.4,1.,.7], symsize=0.5, /NOERASE, xtickformat='(A1)', $
  ytitle= textoidl('v_{\odot,y} [km s^{-1}]'), xrange=[groupcut[0]-.01,groupcut[1]+.01]
oplotbarx, 0.4, linestyle=1
;;For each point indicate the fraction excluded
FOR ii=0L, n_elements(UVWs[*,6])-2 DO BEGIN
    IF UVWs[ii,6] GT 0.005 THEN BEGIN
        xyouts, groupcuts[ii]+.01,-UVWs[ii,2]+.1, $
          strtrim(string(UVWs[ii,6],format='(F4.2)'),2), $
          color=djs_icolor('gray'), charsize=.8
    ENDIF ELSE BEGIN
        xyouts, groupcuts[ii]+.01,-UVWs[ii,2]+.1, $
          strtrim(string(UVWs[ii,6]*nhip,format='(I)'),2), $
          color=djs_icolor('gray'), charsize=.8
    ENDELSE
ENDFOR
ploterror, groupcuts, -UVWs[*,4], UVWs[*,5], psym=8, $
  position=[.1,.1,1.,.4], symsize=0.5, /NOERASE, $
  xtitle= textoidl('p_{MG}'), ytitle= textoidl('v_{\odot,z} [km s^{-1}]'), $
  xrange=[groupcut[0]-.01,groupcut[1]+.01]
oplotbarx, 0.4, linestyle=1
k_end_print

END
