;+
; NAME:
;   color_mag
; PURPOSE:
;   plot color-magnitude diagram from output of TESTEM
;-
pro color_mag

; read in data
  restore, 'hip.sav'
  M_V= hip.vmag-5D *alog10(1000.0D /hip.plx/10.0D )
  BV= hip.bvcolor
  restore, 'testem.sav'

; set up plot
  !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
  set_plot, "PS"
  xsize= 7.5 & ysize= 10.0
  device, file='color_mag.ps',/inches,xsize=xsize,ysize=ysize, $
    xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
  !P.THICK= 4.0
  !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
  !P.CHARSIZE= 1.2
  !P.PSYM= 0
  !P.TITLE= ''
  !X.STYLE= 3
  !X.TITLE= ''
  !X.RANGE= [-0.2,1.7]
  !X.MARGIN= [6,0]
  !X.OMARGIN= [0,0]
  !X.CHARSIZE= 1.0
  !X.TICKLEN= 1.0
  !Y.STYLE= 3
  !Y.TITLE= ''
  !Y.RANGE= [12,-2]
  !Y.MARGIN= [3,0]
  !Y.OMARGIN= [0,0]
  !Y.CHARSIZE= 1.0
  !Y.TICKLEN= 1.0
  !P.MULTI= [0,2,2]
  if dgauss EQ 4 then !P.MULTI= [0,3,5]
  if dgauss EQ 6 then !P.MULTI= [0,5,7]
  xyouts, 0,0,'!3'

; loop over gaussians, identifying likely members
  for j=0,ngauss-1 do begin
    index= where(weight[*,j] GT 0.5,nin)
    if nin GT 0 then begin

; plot cmd
      djs_plot, BV[index],M_V[index],psym=3

    endif
  endfor
  device, /close
end