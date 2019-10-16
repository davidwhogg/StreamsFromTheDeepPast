;+
;   NAME:
;      assign_members
;   PURPOSE:
;      calculate posterior probabilities for stars to be in a clump
;   CALLING SEQUENCE:
;      assign_clump_members, logpost, ydata, ycovar, projection, mean, covar, amp
;   INPUT:
;      ydata    - the data [ndimy,ndata]
;      ycovar   - error covariances [ndimy,ndimy,ndata]
;      projection - [ndimx, ndimy, ndata] non-square matrices
;                  implementing the projection
;      xmean     - means of the gaussians [ndimx,ngauss]
;      xcovar    - covariances of the gaussians [ndimx,ndimx,ngauss]
;      xamp      - amplitudes of the gaussians [ngauss]
;   KEYWORDS:
;   OUTPUT:
;      logpost     - matrix of ln of posterior probabilities [ngauss,ndata]
;   REVISION HISTORY:
;      2008-12-16 - Written - Jo Bovy (NYU)
;      2010-02-24 - ReWritten to be stand-alone - Bovy
;-
FUNCTION BOVY_LOGSUM, array
amax= max(array)
RETURN, amax+alog(total(exp(array-amax),/double))
END
FUNCTION SPECIAL_INVERT, matrix
RETURN, invert(matrix,/double)
END
PRO ASSIGN_MEMBERS, logpost, ydata, ycovar, projection, xmean, $
                          xcovar, xamp

ngauss= n_elements(xamp)
ndimx= n_elements(xmean)/ngauss
ndimy= n_elements(ycovar)/n_elements(ydata)
ndata= n_elements(ydata)/ndimy

twopiterm=0.5*double(ndimy)*alog(2.*!DPI)

logpost= dblarr(ngauss,ndata)

;;Loop over data and Gaussians to find posterior probabilities
FOR ii= 0L, ndata-1 DO BEGIN
    loglike= dblarr(ngauss)
    FOR kk= 0L, ngauss-1 DO BEGIN
        tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,kk]# $
                            projection[*,*,ii]+ycovar[*,*,ii])
        delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,kk]
        loglike[kk]=alog(xamp[kk])+0.5*alog(bovy_determ(tinv,/double,/check) > (machar(/double)).xmin)- $
          0.5*transpose(delta)#tinv#delta-twopiterm
    ENDFOR
    ;;normalize the probabilities
    curr_loglikedata=bovy_logsum(loglike)
    logpost[*,ii]=loglike-curr_loglikedata
ENDFOR
END
