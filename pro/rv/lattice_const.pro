;+
;   NAME:
;      lattice_const
;   PURPOSE:
;      calculate the (conjectured lower bound on) the optimal lattice
;      quantizing constant in n dimensions
;   CALLING SEQUENCE:
;      kappa= lattice_const(n)
;   INPUT:
;      n - dimension
;      seed  - seed for random number generator
;      nsamples - number of MC samples
;   KEYWORDS:
;      zador - calculate the Zador lower bound 
;      upper - calculate the Zador upper bound (takes precedence)
;   OUTPUT:
;      kappa - (lower bound on) optimal lattice quantizing constant
;   REVISION HISTORY:
;      2009-01-12 - Written Bovy (NYU)
;-
FUNCTION LATTICE_CONST, n, zador=zador, upper=upper, $
                        nsamples=nsamples, seed=seed

IF ~keyword_set(seed) THEN seed= -1L

IF keyword_set(upper) THEN BEGIN
    kappa= 1D0/(n*!DPI)*(GAMMA(n/2D0+1D0))^(2D0/n)*GAMMA(1D0+2D0/n)
ENDIF ELSE IF keyword_set(zador) THEN BEGIN
    kappa= 1D0/((n+2D0)*!DPI)*(GAMMA(n/2D0+1))^(2D0/n)
ENDIF ELSE BEGIN
    h= 1D0
    FOR kk= 2L, n+2 DO h+= 1D0/kk
    angle= 0.5D0*acos(1D0/n)
    nint= n/2 > 1
;    COMMON integration, startn,ns, prev, lower
;    startn= n
;    ns= bytarr(nint)
;    prev= dblarr(nint)
;    lower= dblarr(nint)
;    FOR nn= 0L, nint-1 DO lower[nn]= 0.5D0*acos(1D0/(n-1D0-2D0*nn))
;    ns[0]= 1
    f= SCHAFLI(n,angle,nsamples=nsamples,seed=seed)
    kappa= (n+3D0-2D0*h)/(4D0*n*(n+1D0))*((n+1D0)*(GAMMA(n+1D0))^4*f^2)^(1D0/n)
ENDELSE

RETURN, kappa

END
