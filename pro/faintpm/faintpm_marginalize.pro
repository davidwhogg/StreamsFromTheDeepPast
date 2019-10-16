;+
; BUGS:
;  - No proper header.
;  - As with IDL total(), the dimen input is 1-indexed.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_marginalize, dcs,dimen
tiny= 1d-200
foo= dcs
for ii=1,dimen-1 do begin
    foo= foo-min(foo)
    foo= -2.0*alog(total(exp(-0.5*foo),1,/double)+tiny)
endfor
for ii=dimen+1,5 do begin
    foo= foo-min(foo)
    foo= -2.0*alog(total(exp(-0.5*foo),2,/double)+tiny)
endfor
foo= foo-min(foo)
return, foo
end
