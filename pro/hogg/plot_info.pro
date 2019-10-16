;+
; NAME:
;   plot_info
; PURPOSE:
;   make info.ps, a plot of information in ex_max solutions
;-
pro plot_info

; read data
  for ng= 3L,40 do begin
    filename= 'increment.'+string(ng,format='(I2.2)')+'gauss.sav'
    filename = findfile(filepath(filename,root_dir='.'),count=nfile)
    if nfile EQ 1 then begin
      splog, filename[0]
      restore, filename[0]
      if NOT keyword_set(entropylist) then begin
        entropylist= entropy
        ngausslist= ngauss
      endif else begin
        entropylist= [entropylist,entropy]
        ngausslist= [ngausslist,ngauss]
      endelse
    endif
  endfor

; compute information in model parameters
  bitsperparameter= round(alog(sqrt(ndata))/alog(2))    ; guess?
  parameterspergauss= (dgauss+2L)                 ; for spherical gaussians
;  parameterspergauss= (dgauss+1L)*(dgauss+2L)/2L  ; for arbitrary gaussians
  parameterinfo= bitsperparameter*parameterspergauss*ngausslist

; setup postscript plot
  !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
  set_plot, "PS"
  xsize= 7.5 & ysize= 7.5
  device, file='info.ps',/inches,xsize=xsize,ysize=ysize, $
    xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits_per_pixel=24
  !P.THICK= 4.0
  !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
  !P.CHARSIZE= 1.4
  !P.SYMSIZE= !P.CHARSIZE
  !P.PSYM= 0
  !P.TITLE= ''
  !X.STYLE= 3
  !X.TITLE= 'number of gaussians  !8M!3'
  !X.RANGE= 0
  !X.MARGIN= [4.7,0]*!P.CHARSIZE   ; correct for four-digit y-axis numbers
;  !X.MARGIN= [7.5,0]*!P.CHARSIZE   ; correct for exponential notation
  !X.OMARGIN= [3,1]*!P.CHARSIZE
  !X.CHARSIZE= !P.CHARSIZE
  !X.TICKLEN= 1.0
  !Y.STYLE= 3
  !Y.RANGE= 0
  !Y.MARGIN= [1.5,0]*!P.CHARSIZE
  !Y.OMARGIN= 0.5*!X.OMARGIN
  !Y.CHARSIZE= !P.CHARSIZE
  !Y.TICKLEN= 1.0
  !P.MULTI= 0
  xyouts, 0,0,'!3'

; pick ngauss for comparison
  !Y.TITLE= 'information  !8S!3!dcomp!n  (bits)'
  entropylist= entropylist-entropylist[0]
  parameterinfo= parameterinfo-parameterinfo[0]

; plot entropies
  djs_plot, ngausslist,entropylist-parameterinfo, $
    thick=2.0*!P.THICK,psym=-6
  djs_oplot, ngausslist,entropylist-2.0*parameterinfo, $
    linestyle=1,thick=2.0*!P.THICK

  device,/close
  return
end
