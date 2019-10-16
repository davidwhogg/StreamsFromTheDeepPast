;+
; NAME:
;   read_radial_velocities
; PURPOSE:
;   read what VRs we can find for Hipparcos stars
; CALLING SEQUENCE:
;   radv= read_radial_velocities(/savefile)
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
;   savefile    - if set, make a saveset
; OUTPUTS:
;   radv        - structure of radial velocity information
; BUGS:
;   Hard-coding atrocities everywhere.
; REVISION HISTORY:
;   2003-02-18  written - Hogg
;-
function radv_struct
    str= {hipparcos_radial_velocity, $
          hip     : 0L,  $ ; HIP number
          ra      : 0D,  $ ; right ascension (2000.00) (deg)
          dec     : 0D,  $ ; declination (2000.00) (deg)
          vmag    : 0D,  $ ; apparent V magnitude or similar (mag)
          vr      : 0D,  $ ; radial velocity (km/s)
          vr_err  : 0D   $ ; error (km/s)
         }
return, str
end

function hogg_cast_double,string
    if stregex(string,'[0-9]') GE 0 then return,double(string) else return,0.0
end

function read_radial_velocities,savefile=savefile
help, savefile

; read Barbier-Brossat et al 2000
; open file
indir= '/global/data/hipparcos/vr/3_3213'
filename= 'catalog.dat'
openr, rlun,indir+'/'+filename,/get_lun

; loop over lines, reading only HIP stars
while NOT eof(rlun) do begin
    line=''
    readf, rlun,line
    hipno= -1
    hipstr= strmid(line,64,6)
    if hipstr NE '      ' then hipno= long(strmid(line,64,6))
    if hipno GE 0 then begin

; fill in radv structure element
        radv1= radv_struct()
        radv1.hip= hipno
        radv1.ra= 1.5D1*(hogg_cast_double(strmid(line,78,2))+ $
                         hogg_cast_double(strmid(line,81,5))/6D1)
        radv1.dec=       hogg_cast_double(strmid(line,86,3))+ $
                         hogg_cast_double(strmid(line,90,2))/6D1
        radv1.vmag=      hogg_cast_double(strmid(line,94,5))
        radv1.vr=        hogg_cast_double(strmid(line,112,7))
        radv1.vr_err=    hogg_cast_double(strmid(line,121,4))

; append to list
	if not keyword_set(radv) then radv= radv1 else radv= [radv,radv1]
    endif
endwhile
close, rlun
free_lun, rlun

; make saveset
if keyword_set(savefile) then begin
    filename= '/global/data/hipparcos/hip-radv.sav'
    splog, 'saving '+filename
    save, radv,filename=filename
endif
return, radv
end
