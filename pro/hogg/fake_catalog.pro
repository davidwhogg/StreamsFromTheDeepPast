;+
; NAME:
;   fake_catalog
; PURPOSE:
;   make a fake catalog of data, given outputs of ex_max or ex_max_proj
; CALLING SEQUENCE:
;   data= fake_catalog(seed,amp,mean,var)
; INPUTS:
;   seed         random number seed (should be a long)
;   amp          [M] vector of numbers of stars per gaussian
;   mean         [d,M] array of d-dimensional mean vectors
;   var          [d,d,M] array of dxd variance matrices
; OPTIONAL INPUTS:
;   ndata        number (N) of data points to create (default to total(amp))
;   qa           file for QA plots, if desired
; KEYWORDS:
;   quiet        be quiet about it
; OUTPUTS:
;   data         [d,N] array of fake data
; OPTIONAL OUTPUTS:
;   gauss        [N] array of which gaussian each data point came from
; BUGS:
;   Bombs when d=1 (dimen=1 in code) (on eigenql call).
; REVISION HISTORY:
;   2001-Sep-05  written by Hogg (NYU)
;-
function fake_catalog, seed,amp,mean,var,ndata=ndata,qa=qa,gauss=gauss, $
                       quiet=quiet

  if NOT keyword_set(seed) then seed=-1L

; check dimensions
  ngauss= n_elements(amp)
  dimen= n_elements(mean)/ngauss
  if NOT keyword_set(ndata) then ndata= long(round(total(amp)))
  if NOT keyword_set(quiet) then $
    splog, ndata,' data points,',dimen,' dimensions,',ngauss,' gaussians'

; ram inputs and outputs onto correct format
  amp= reform([amp],ngauss)
  mean= reform([mean],dimen,ngauss)
  var= reform([var],dimen,dimen,ngauss)
  gauss= lonarr(ndata)
  data= dblarr(dimen,ndata)

; set random seed
  seed= long(seed)

; assign data points to gaussians
  frac= total(amp,/cumulative)/total(amp)
  rand= randomu(seed,ndata)
  for j=ngauss-1,0,-1 do begin
    index= where(rand LT frac[j],nin)
    if nin GT 0 then gauss[index]= j
  endfor
  rand= 0

; assign mean vectors
  for i=0L,ndata-1 do data[*,i]= mean[*,gauss[i]]

; assign unscaled delta vectors
  dvector= randomn(seed,dimen,ndata)

; transform unscaled delta vectors with variance matrices
  eval= dblarr(dimen,ngauss)
  evec= dblarr(dimen,dimen,ngauss)
  for j=0L,ngauss-1 do begin
      tmpvar=0.5d*(var[*,*,j]+transpose(var[*,*,j]))
      tmpeval= eigenql(tmpvar,/double,eigenvectors=tmpevec)
      eval[*,j]= sqrt(tmpeval)
      evec[*,*,j]= tmpevec
  endfor
  for i=0L,ndata-1 do begin
    dvector[*,i]= eval[*,gauss[i]]*dvector[*,i]
    dvector[*,i]= evec[*,*,gauss[i]]#dvector[*,i]
  endfor

; add up!
  data= data+dvector

; make postscript file, if necessary
  if keyword_set(qa) then begin
    !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
    set_plot, "PS"
    xsize= 7.5 & ysize= 7.5
    device, file=qa,/inches,xsize=xsize,ysize=ysize, $
      xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
    !P.THICK= 2.0
    !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
    !P.CHARSIZE= 1.2
    !P.PSYM= 3
    !P.TITLE= ''
    !X.STYLE= 3
    !X.TITLE= ''
    !X.RANGE= 0
    !X.MARGIN= [12,0]
    !X.OMARGIN= [-4,0]
    !X.CHARSIZE= 1.0
    !Y.STYLE= 3
    !Y.TITLE= ''
    !Y.RANGE= 0
    !Y.MARGIN= [6,0]
    !Y.OMARGIN= [-2,0]
    !Y.CHARSIZE= 1.0
    !P.MULTI= [0,2,2]
    xyouts, 0,0,'!3'
    for d1=0,dimen-1 do for d2=d1+1,dimen-1 do begin
      djs_plot, data[d1,*],data[d2,*],xtitle=string(d1),ytitle=string(d2)
    endfor
    device, /close
  endif

; return and end
  return, data
end
