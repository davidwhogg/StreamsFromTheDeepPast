;+
; NAME:
;   shuffle
; PURPOSE:
;   rearrange 1-dimensional array element order
; CALLING SEQUENCES:
;   index = shuffle(seed, array)
;   array2 = array[shuffle(seed,array)]
; BUGS:
;   Incredibly slow, stupid implementation.
;-
function shuffle, seed, array

  nel= n_elements(array)
  orig= lindgen(nel)
  index= lonarr(nel)
  for i= 0L,nel-2 do begin
    el= orig[floor(randomu(seed)*double(n_elements(orig)))]
    index[i]= el
    orig= orig[where(orig NE el)]
  endfor
  index[i]= orig
  return, index
end
