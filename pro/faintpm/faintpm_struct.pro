;+
; BUGS:
;  - No proper header.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_struct,nn
foo= {faintpm, $
      nepoch:     0, $
      meant:    0.0, $; MJD
      deltat:   0.0, $; yr
      sigmat:   0.0, $; yr
      epoch:    0.0, $; MJD
      flux:    -1D0, $; nMgy
      ra:      -1D0, $; deg
      dec:     -1D2, $; deg
      radot:    0D0, $; mas/yr
      decdot:   0D0, $; mas/yr
      fluxerr:   -1D0, $; nMgy
      raerr:     -1D0, $; deg
      decerr:    -1D0, $; deg
      radoterr:  -1D0, $; mas/yr
      decdoterr: -1D0, $; mas/yr
      filename: ''}
return, replicate(foo,nn)
end
