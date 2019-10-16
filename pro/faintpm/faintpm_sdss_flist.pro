;+
; INPUTS:
;  RA,Dec       - position (deg)
;  name         - name for this object (no spaces, please)
; OPTIONAL INPUTS:
;  run          - list of runs to consider (otherwise consider all
;                 photometric runs)
; OUTPUTS:
;  [return]     - list of fields, SORTED BY RUN
; BUGS:
;  - No proper header.
;  - Returns only SDSS z-band images.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_sdss_flist, ra,dec,run=run
dra= 1.0/60.0 ; deg
ddec= 1.0/60.0 ; deg
pixscale= 0.396/3600.0 ; deg
bigast = hogg_make_astr(ra,dec,dra,ddec,pixscale=pixscale,npixround=1)
rerun= [137,161]
flist = smosaic_flist(bigast,run,rerun=rerun, $
                      filter=4,minscore=0.5,/allruns)
if (NOT keyword_set(flist)) then begin
    splog, 'no overlapping fields'
    return, 0
endif
flist= flist[sort(flist.run)]
splog, 'found',n_elements(flist),' fields'
return, flist
end
