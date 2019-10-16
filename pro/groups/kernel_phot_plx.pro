;+
;   NAME:
;      kernel_phot_plx
;   PURPOSE:
;      calculate the probability of a parallax given photometric data
;   INPUT:
;      plx - parallax to evaluate the probability of (can be an array)
;      bvcolor - B-V color of the plx star
;      vmag - apparent magnitude of the plx star
;      sig_plx - uncertainty on the plx
;      data - set of [B-V,Mv] [N,2] training data
;      width - width parameter of the kernel function
;   KEYWORDS:
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      silent - be silent
;   OUTPUT:
;      p(plx|...)
;   HISTORY:
;      2009-11-05 - Written - Bovy (NYU)
;-
FUNCTION plx_mv_v, mv, v
;;Returns the parallax in mas given the absolute and the apparent
;;magnitude
out= dblarr(n_elements(mv))
FOR ii=0L, n_elements(mv)-1 DO out[ii]= 10.^(-(v-mv[ii])/5.+2.)
RETURN, out
END
FUNCTION tricube_kernel, t
;;No checking whether this is in the range t<0
out= dblarr(n_elements(t))
FOR ii=0L, n_elements(t)-1 DO out[ii]= (1.-t[ii]^3.)^3.
RETURN, out
END
FUNCTION epanechnikov_kernel, t
;;No checking whether this is in the range t<0
out= dblarr(n_elements(t))
for ii=0L, n_elements(t)-1 DO out[ii]= .75*(1.-t[ii]^2)
RETURN, out
END
FUNCTION gaussian_kernel, t
out= dblarr(n_elements(t))
for ii=0L, n_elements(t)-1 DO out[ii]= exp(-t[ii]^2./2.)
RETURN, out
END
FUNCTION kernel_phot_plx, plx, bvcolor, vmag, sig_plx,data, $
                          width=width, kernel=kernel, silent=silent
IF ~keyword_set(width) THEN width= 0.1
IF ~keyword_set(kernel) THEN kernel= 'tricube'
IF kernel EQ 'tricube' THEN BEGIN
    kernel_function= 'tricube_kernel'
    relevant_indx= where(abs(data[*,0]-bvcolor)/width LT 1.)
    if relevant_indx[0] EQ -1 THEN BEGIN
        IF ~keyword_set(silent) THEN splog, "No data point in width range"
        RETURN, -(machar(/double)).xmax
    ENDIF
    relevant_data= data[relevant_indx,*]
ENDIF ELSE IF kernel EQ 'epanechnikov' THEN BEGIN
    kernel_function= 'epanechnikov_kernel'
    relevant_indx= where(abs(data[*,0]-bvcolor)/width LT 1.)
    if relevant_indx[0] EQ -1 THEN BEGIN
        IF ~keyword_set(silent) THEN splog, "No data points in width range"
        RETURN, -(machar(/double)).xmax
    ENDIF
    relevant_data= data[relevant_indx,*]
ENDIF ELSE BEGIN
    kernel_function= 'gaussian_kernel'
    relevant_data= data
ENDELSE
;;Calculate the weigths
ndata= n_elements(relevant_data[*,0])
weights= CALL_FUNCTION(kernel_function,abs(relevant_data[*,0]-bvcolor)/width)
weights= weights/total(weights)
;;and the predicted parallaxes
modelplx= plx_mv_v(relevant_data[*,1],vmag)
vars= dblarr(ndata)+sig_plx^2.
RETURN, oned_sum_gaussians(plx,modelplx,vars,weights)
END
