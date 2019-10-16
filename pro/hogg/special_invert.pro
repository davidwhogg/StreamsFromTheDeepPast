;+
; NAME:
;   special_invert
; PURPOSE:
;   inverts symmetric matrices (with 2D as a special case)
; INPUTS:
;   matrix - input symmetric matrix
; OUTPUTS:
; BUGS:
;   no checking of symmetry or anything else
; DEPENDENCIES:
;   idlutils
; REVISION HISTORY:
;   2003-Feb-21  written by Blanton (NYU)
;-
function special_invert,matrix

if(n_elements(matrix) eq 4) then begin
    inverse=dblarr(2,2)
    determinant=matrix[0,0]*matrix[1,1]-matrix[0,1]^2
    if(determinant ne 0.D) then begin
        inverse[0,0]=matrix[1,1]
        inverse[0,1]=-matrix[0,1]
        inverse[1,0]=-matrix[0,1]
        inverse[1,1]=matrix[0,0]
        inverse=inverse/determinant
    endif
endif else begin
    inverse=invert(matrix,/double)
endelse

return,inverse

end
