;+
;   NAME:
;      threed_sum_gaussians
;
;   PURPOSE:
;      implements a function that is a sum of 3D Gaussians
;
;   CALLING SEQUENCE:
;      Result = threed_sum_gaussians(z, mean, covar [, amp=amp])
;
;   INPUTS:
;      z     : (3D) point to calculate the function @ (can be array=
;      dblarr(d,N))
;      mean  : vector of 3D means (dblarr(3,ngauss))
;      covar : vector of 3D covariances (dblarr(3,3,ngauss)
;      amp   : relative amplitudes (set to 1 if there is only 1
;              Gaussian, rescaled if they don't add up to one, equal
;              weights when not set)
;
;   OUTPUT:
;      Result : function value
;
;   REVISION HISTORY:
;      06-03-08   - Written Bovy
;-
FUNCTION THREED_SUM_GAUSSIANS, z, mean, covar, amp=amp
d= n_elements(covar)/n_elements(mean)
;;Perhaps jump to oned_sum_gaussians if d=1, rewrite
ngauss= n_elements(mean[0,*])
IF ~keyword_set(amp) THEN amp= dblarr(ngauss)+1D/double(ngauss) $
  ELSE IF TOTAL(amp) NE 1D THEN amp= amp/total(amp)
nz= n_elements(z)/d
Result= dblarr(nz)


;;If this is 1d, use oned_sum_gaussians
IF d EQ 1 THEN BEGIN
    Result= oned_sum_gaussians(z, mean, covar, amp)
    RETURN, Result
ENDIF

;;First diagonalize all covariances, rotate the means, and compute
;;some constants that will be useful
xmean=mean
xcovar=covar
xamp= amp

detj= dblarr(ngauss)+1D
constj= dblarr(ngauss)

var= dblarr(d,ngauss)
FOR gg=0L, ngauss-1 DO BEGIN
    ;;IDL (as we know) sucks, and therefore we need to symmetrize the
    ;;covariances before calculating their eigenvalues
    xcovar[*,*,gg]= 1/2D*(xcovar[*,*,gg]+transpose(xcovar[*,*,gg]))
    var[*,gg]= eigenql(xcovar[*,*,gg],/double,eigenvectors=ev)
    xcovar[*,*,gg]= ev
    xmean[*,gg]= ev##xmean[*,gg]
    FOR dd=0L, d-1 DO BEGIN
        constj[gg]+= 1/var[dd,gg]*xmean[dd,gg]^2
        detj[gg]*= var[dd,gg]
    ENDFOR
    detj[gg]= 1D/sqrt(detj[gg])
ENDFOR

twopid= 1D/(2D0*!DPI)^(d/2D)

FOR zz=0L, nz-1 DO BEGIN
    FOR gg=0L, ngauss-1 DO BEGIN
        tempz= transpose(xcovar[*,*,gg])##z[*,zz]
        arg_exp= constj[gg]
        FOR dd=0L, d-1 DO arg_exp+= 1/var[dd,gg]*$
          tempz[dd]*(tempz[dd]-2*xmean[dd,gg])
        Result[zz]+= amp[gg]*twopid*detj[gg]*$
          exp(-1D/2D*arg_exp)
    ENDFOR
ENDFOR

RETURN, Result

END
