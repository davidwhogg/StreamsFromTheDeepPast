;+
;   NAME:
;      twod_sum_gaussians
;
;   PURPOSE:
;      implements a function that is a sum of 2D Gaussians
;
;   CALLING SEQUENCE:
;      Result = twod_sum_gaussians(z, mean, covar [, amp])
;
;   INPUTS:
;      z     : (2D) point to calculate the function @
;      mean  : vector of 2D means (dblarr(2,ngauss))
;      covar : vector of 2D covariances (dblarr(2,2,ngauss)
;      amp   : relative amplitudes (set to 1 if there is only 1
;              Gaussian
;
;   OUTPUT:
;      Result : function value
;
;   REVISION HISTORY:
;      06-03-08   - Written Bovy
;-
FUNCTION TWOD_SUM_GAUSSIANS, z, mean, covar, amp
Result=0D
ngauss= n_elements(mean[0,*])
IF ngauss EQ 1 THEN BEGIN
    amp=1D
    sigmax= sqrt(covar[0,0])
    sigmay= sqrt(covar[1,1])
    corr= covar[0,1]/(sigmax*sigmay)
    ze= (z[0]-mean[0])^2/covar[0,0]+(z[1]-mean[1])^2/covar[1,1]-2D0*corr*(z[0]-mean[0])*(z[1]-mean[1])/(sigmax*sigmay)
    Result= amp/(2D0*!DPI*sigmax*sigmay*sqrt(1D0-corr^2))*exp(-ze/(2D*(1D0-corr^2)))
ENDIF ELSE BEGIN
    FOR ii=0L, ngauss-1 DO BEGIN
        sigmax= sqrt(covar[0,0,ii])
        sigmay= sqrt(covar[1,1,ii])
        corr= covar[0,1,ii]/(sigmax*sigmay)
        ze= (z[0]-mean[0,ii])^2/covar[0,0,ii]+(z[1]-mean[1,ii])^2/covar[1,1,ii]-2D*corr*(z[0]-mean[0,ii])*(z[1]-mean[1,ii])/(sigmax*sigmay)
        Result+= amp[ii]/(2D0*!DPI*sigmax*sigmay*sqrt(1D0-corr^2))*exp(-ze/(2D0*(1D0-corr^2)))
    ENDFOR
ENDELSE

RETURN, Result

END
