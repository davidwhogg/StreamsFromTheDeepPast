;+
;   NAME:
;      calc_logpost_fixedpij
;   PURPOSE:
;      calculate the log posterior of a setting of alpha, age, and
;      metallicity given the data and a moving group
;   INPUT:
;      group - which group?
;      alpha - alpha parameter
;      age - log age
;      z - metallicity
;      hip - data structure
;      back_loglike - loglikelihood under the background model
;      sspspread - spread to add to the SSP prediction
;      isodir - directory that holds the isochrones
;      logpost - log of group membership probabilities
;   KEYWORDS:
;      gzip - if set, gzip and gunzip the isochrones
;   OUTPUT:
;      the log of the posterior
;      outage - actual age for which the calculation was performed
;      outz - actual z for which the calculation was performed
;   REVISION HISTORY:
;      2009-11-09 - Written - Bovy
;-
FUNCTION CALC_LOGPOST_FIXEDPIJ, group, alpha, age, z, $
                                back_loglike=back_loglike,hip=hip,$
                                sspspread=sspspread, logpost=logpost,$
                                isodir=isodir, outage=outage, outz=outz, $
                                gzip=gzip
IF ~keyword_set(sspspread) THEN sspspread= 0.
IF ~keyword_set(isodir) THEN isodir='isochrones/'
nhip= n_elements(hip.hip)
;;Read the isochrones
isoname= get_isoname(z,metal=zout)
catch, error;;If the following fails, try again
IF error NE 0 THEN BEGIN  
    PRINT, 'Error index: ', Error_status  
    PRINT, 'Error message: ', !ERROR_STATE.MSG  
    PRINT, 'Waiting and trying again...'
    wait, 10
ENDIF 
IF keyword_set(gzip) THEN BEGIN
    cmd= 'gunzip '+isodir+isoname+'.gz'
    spawn, cmd
ENDIF
iso= read_isochrone(isodir+isoname)
IF keyword_set(gzip) THEN BEGIN
    cmd= 'gzip '+isodir+isoname
    spawn, cmd
ENDIF
catch, /cancel
ages= iso.h01[UNIQ(iso.h01, SORT(iso.h01))]
ii=0L
WHILE ii LT n_elements(ages)-1 DO BEGIN
    IF ages[ii] GE age THEN BREAK
    ii+= 1
ENDWHILE
IF ii GE n_elements(ages)-1 THEN BEGIN
    isoage= ages[ii]
ENDIF ELSE IF ii EQ 0 THEN BEGIN
    isoage= ages[0]
ENDIF ELSE IF abs(age-ages[ii]) GT abs(age-ages[ii-1]) THEN BEGIN
    isoage= ages[ii-1]
ENDIF ELSE BEGIN
    isoage= ages[ii]
ENDELSE
thisiso_indx= where(iso.h01 EQ isoage)
nthisiso= n_elements(thisiso_indx)
thisiso= dblarr(2,nthisiso)
thisiso[0,*]= iso.h09[thisiso_indx]
thisiso[1,*]= iso.h10[thisiso_indx]

out= 0.
FOR ii=0L, nhip-1 DO BEGIN
    BV= hip[ii].bvcolor
    V= hip[ii].vmag
    modelplx= phot_plx(BV,V,thisiso,/ms)
    IF modelplx EQ -1. THEN BEGIN
        out+= exp(logpost[group,ii])*(back_loglike[ii]+alog(alpha))
    ENDIF ELSE BEGIN
        thisspread= hip[ii].e_plx^2.+sspspread^2.*hip[ii].plx^2.
        out+= exp(logpost[group,ii])*alog((1.-alpha)/sqrt(2.*!DPI*thisspread)*exp(-0.5*(hip[ii].plx-modelplx)^2./thisspread)+alpha*exp(back_loglike[ii]))
    ENDELSE
ENDFOR

RETURN, out
END
