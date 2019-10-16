;+
;   NAME:
;      density_estimate_vdist
;   PURPOSE:
;      estimate the velocity distribution using a Gaussian Kernel
;      convolution
;   CALLING SEQUENCE:
;      density_estimate_vdist
;   INPUT:
;   OUTPUT:
;   REVISION HISTORY:
;      2008-10-17 - Started Bovy
;-
PRO DENSITY_ESTIMATE_VDIST, hipparcos=hipparcos, geneva=geneva, savefilename=savefilename

maxitt= 100
tol= 1D-1

IF keyword_set(hipparcos) THEN BEGIN
;;Load the Hipparcos data
    dehnen_giants_filename= '/global/data/hipparcos/hip-dehnen-giants.sav'
;IF (keyword_set(giants) OR keyword_set(all)) THEN BEGIN
    IF file_test(dehnen_giants_filename) THEN BEGIN
        restore, dehnen_giants_filename
    ENDIF ELSE BEGIN
        splog, 'Error: Dehnen-Giants-savefile does not exist in ', $
          dehnen_giants_filename
        splog, 'Dehnen-Giants-savefile must exist, returning...'
        RETURN
    ENDELSE
    hip_giants= hip
;ENDIF
    dehnenfilename= '/global/data/hipparcos/hip-dehnen.sav'
    IF file_test(dehnenfilename) THEN BEGIN
        restore, dehnenfilename
    ENDIF ELSE BEGIN
        splog, 'Error: Dehnen-savefile does not exist in ', $
          dehnenfilename
        splog, 'Dehnen-savefile must exist, returning...'
        RETURN
    ENDELSE
    hip_main= hip
    nhip_main= n_elements(hip_main)
;sindx= sort(hip_main.bvcolor)
;IF keyword_set(colors) THEN BEGIN
;    cindx= dblarr(nsubsamp+1)
;    cindx[0]=0
;    cindx[nsubsamp]= nhip-1
;    FOR ii=0L, nsubsamp-2-keyword_set(giants) -keyword_set(all) DO cindx[ii+1]= $
;      max(where(hip[sindx].bvcolor LT colors[ii]))
;ENDIF
    
;;Merge the two
    hip= [hip_main, hip_giants]
    
    ;;Prepare for take-off
    N= n_elements(hip)
    ydata= transpose([[dblarr(N)],[hip.vl], [hip.vb]])
    ycovar_tmp= hip.vlvbc
    projection= hip.sm
;;Add radial velocities if those are given as well
    ycovar= dblarr(3,3,N)
    FOR gg= 0, N-1 DO BEGIN
        ycovar[*,*,gg]=[[10000.^2,0.,0.],[0.,ycovar_tmp[0,0,gg],ycovar_tmp[0,1,gg]],[0.,ycovar_tmp[1,0,gg],ycovar_tmp[1,1,gg]]]
    ENDFOR
ENDIF ELSE IF keyword_set(geneva) THEN BEGIN
    genevafilename= '/global/data/rv/gcs/gcs.sav'
    IF file_test(genevafilename) THEN BEGIN
        restore, genevafilename
    ENDIF ELSE BEGIN
        splog, 'Error: GCS-savefile does not exist in ', $
          genevafilename
        splog, 'GCS-savefile must exist, returning...'
        RETURN
    ENDELSE
    ;;Prepare for take-off
    gcs= gcs[where((abs(gcs.vl) LE 1000D) AND (abs(gcs.vb) LE 1000D) AND (abs(gcs.vr) LE 1000D))]
    N= n_elements(gcs)
    ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
    ycovar= gcs.vlvbvrc
    projection= gcs.sm
ENDIF

;;now rotate each data point to the U,V,W frame (using projection-1)
FOR ii=0L, N-1 DO BEGIN
    ydata[*,ii]= transpose(projection[*,*,ii])##ydata[*,ii]
    ycovar[*,*,ii]= (transpose(projection[*,*,ii])##ycovar[*,*,ii])##projection[*,*,ii]
ENDFOR


;;Compute the density distribution
;;amplitudes are all equal, so these will be be defaults in the rest
;;of my code (maybe change, think about memory usage)
;;means are just the ydata
;;covariances are obtained by adding a data point dependent smoothing

;;Initial estimate
h= 10D
FOR ii=0L, N-1 DO BEGIN
    ycovar[*,*,ii]= h*IDENTITY(3,/double);ycovar[*,*,ii]+h*IDENTITY(3,/double)
ENDFOR

lambda= dblarr(N)+1D
;;Iterate to find optimal lambda_i s
diff= 2*tol
itt=0L
WHILE (diff GT tol and itt LT maxitt) DO BEGIN
    itt++
    splog, "Iteration "+strtrim(string(itt),2)
    ;;Find lambda s (most time consuming), first find density
    ;;estimates at all the data points
    oldlambda= lambda
    lambda= threed_sum_gaussians(ydata,ydata,ycovar)
    wait, 10
;    stop
    g= total(alog(lambda),/double)/double(N)
    g= exp(g)
    lambda= lambda/g
    lambda= 1/sqrt(lambda)
    diff= max(abs(lambda-oldlambda))
    splog, "Change in lambdas = "+$
      strtrim(string(diff),2)
    ;;reestimate the density
    h= 10D
    FOR ii=0L, N-1 DO BEGIN
        ycovar[*,*,ii]= h*lambda[ii]*IDENTITY(3,/double);ycovar[*,*,ii]+h*lambda[ii]*IDENTITY(3,/double)
    ENDFOR
ENDWHILE

;;Save the result
save, filename=savefilename

;;Plot the density estimate

label= "ALdens"
plot_projected_gaussians_wrap, ydata, ycovar, $
  xrange=[-130,120], yrange=[-120,60], zrange=[-70,70], $
;      xrange=[-50,50], yrange=[-50,50], zrange=[-50,50], $
  basefilename=label,label=label, grid=256
scatter_projected_gaussians_wrap, ydata, ycovar, $
  xrange=[-130,120], yrange=[-120,60], zrange=[-70,70], $
;      xrange=[-50,50], yrange=[-50,50], zrange=[-50,50], $
  basefilename=label,label=label


END
