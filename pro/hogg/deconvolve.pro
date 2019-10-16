;+
; NAME:
;   deconvolve
; PURPOSE:
;   Deconvolution by transforming to Fourier space.
; INPUTS:
;   data        [N1,N2,N3] array of input data
;   green       [N1,N2,N3] description of the Green's function to deconvolve
; OUTPUTS:
;   data        deconvolved data
; BUGS:
;   
; DEPENDENCIES:
;   idlutils
; REVISION HISTORY:
;   2001-Aug-06  written by Blanton (NYU)
;-
pro deconvolve,data,green

data=fft(data,/overwrite)
greenfft=fft(green)
indx=where(greenfft ne complex(0.,0.))
data[indx]=data[indx]/greenfft[indx]/n_elements(greenfft)
data=fft(data,/overwrite,/inverse)
;data=real(data)

end
