;+
;   NAME:
;      calc_entropy_sum_gaussians
;   PURPOSE:
;      calculate the entropy of a distribution that is a sum of
;      Gaussians
;   CALLING SEQUENCE:
;      calc_entropy_sum_gaussians
;   INPUT:
;      xamp   - amplitudes [K]
;      xmean  - means [d,K]
;      xcovar - covar [d,d,K]
;   OPTIONAL INPUTS:
;      nboot  - number of bootstrap samples to use
;      nsamplings - number of samplings per bootstrap sample
;      seed   - random number generator seed
;   KEYWORDS:
;   OUTPUT:
;      entropy  - calculated entropy
;   OPTIONAL OUTPUT:
;      sigma_H  - uncertainty in the entropy (rms)
;   REVISION HISTORY:
;      2008-10-27 - Written Bovy
;-
PRO CALC_ENTROPY_SUM_GAUSSIANS, xamp, xmean, xcovar, $
                                nboot=nboot, nsamplings=nsamplings, $
                                entropy=entropy, sigma_H=sigma_H, seed=seed

;;Check inputs
d= n_elements(xcovar)/n_elements(xmean)
K= n_elements(xamp)
IF ~keyword_set(nboot) THEN nboot= 10
IF ~keyword_set(nsamplings) THEN nsamplings= 100^d
IF ~keyword_set(seed) THEN seed= 1L

boot_H= dblarr(nboot)
FOR bb=0L, nboot-1 DO BEGIN
    ;;sample
    sample_gaussians, nsamples=nsamplings, mean=xmean, covar=xcovar, $
      amp=xamp, sample=sample, seed=seed
    logp= threed_sum_gaussians(sample, xmean, xcovar, amp=xamp)
    logp= alog(logp)
    boot_H[bb]= -mean(logp,/double,/NaN)
ENDFOR

entropy= mean(boot_H,/double)
IF arg_present(sigma_H) THEN sigma_H= stddev(boot_H,/double,/NaN)/sqrt(double(nboot))

END
