;+
;
;-
PRO PLOT_LATTICE_CONST, savefile=savefile, seed=seed

IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(savefile) THEN savefile= "lattice.sav"

IF file_test(savefile) THEN Begin
    splog, "restoring savefile "+savefile
    restore, savefile
ENDIF ELSE BEGIN
    nd= 249
    
    n= lindgen(nd)+1L
    k= dblarr(nd)
    kz= dblarr(nd)
    ku= dblarr(nd)
    
    FOR dd= 0L, nd-1 DO BEGIN
        splog, "Working on dim "+strtrim(string(n[dd]),2)
        IF n[dd] GE 20 THEN extramc= 10 ELSE extramc= 1
        k[dd]= lattice_const(n[dd],seed=seed,nsamples=10000000L*extramc)
        kz[dd]= lattice_const(n[dd],/zador)
        ku[dd]= lattice_const(n[dd],/upper)
        
    ENDFOR
    save, filename=savefile, n, k, kz, ku, nd
ENDELSE

djs_plot, n, k, psym=-3, xtitle="Dimension n", ytitle="\kappa_n", $
  xrange=[0,n[nd-1]+1], yrange=[1D0/(2*!DPI*exp(1)),0.084], linestyle=0
djs_oplot, n, kz, linestyle=2, psym=-3
djs_oplot, n, ku, linestyle=3, psym=-3

END
