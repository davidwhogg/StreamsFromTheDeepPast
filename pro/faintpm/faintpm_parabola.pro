;+
; BUGS:
;  - No proper header.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_parabola, xx,yy,xxfine,xmin,ymin,sigmax
aa= [[replicate(1D0,3)],[xx],[xx*xx]]
aainv= invert(aa,/double)
pp= aainv#yy
; pp= pp+aainv#(yy-aa#pp)
yyfine= pp[0]+pp[1]*xxfine+pp[2]*xxfine*xxfine
xmin= -0.5*pp[1]/pp[2]
ymin= pp[0]+pp[1]*xmin+pp[2]*xmin*xmin
sigmax= -1.0
if (pp[2] GT 0) then sigmax= 1.0/sqrt(pp[2])
if (sigmax LE 0.0) then begin
    ymin= min(yy,imin)
    xmin= xx[imin]
endif
return, yyfine
end
