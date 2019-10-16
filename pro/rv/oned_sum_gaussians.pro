;+
;   NAME:
;      oned_sum_gaussians
;
;   PURPOSE:
;      implements a function that is a sum of 1D Gaussians
;
;   CALLING SEQUENCE:
;      Result = oned_sum_gaussians(z, mean, var [, amp])
;
;   INPUTS:
;      z     : (1D) point to calculate the function @
;      mean  : vector of 1D means (dblarr(ngauss))
;      var : vector of 1D variances (dblarr(ngauss)
;      amp   : relative amplitudes (set to 1 if there is only 1
;              Gaussian)
;
;   OUTPUT:
;      Result : function value
;
;   REVISION HISTORY:
;      06-03-08   - Written Bovy
;-
FUNCTION ONED_SUM_GAUSSIANS, z, mean, var, amp
nz= n_elements(z)
Result=dblarr(nz)
ngauss= n_elements(amp)
FOR zz=0L, nz-1 DO BEGIN
    IF ngauss EQ 1 THEN BEGIN
        Result[zz]=1/sqrt(2D*!DPI*var)*exp(-(z[zz]-mean)^2/(2D*var))
    ENDIF ELSE BEGIN
        FOR ii=0L, ngauss-1 DO BEGIN
            ze= (z[zz]-mean[ii])^2/(2D*var[ii])
            Result[zz]+= amp[ii]/(sqrt(2D*!DPI*var[ii]))*exp(-ze)
        ENDFOR
    ENDELSE
ENDFOR

RETURN, Result

END
