;+
; NAME:
;  faintpm_stack
; PURPOSE:
;  Stack a multi-exten FITS file on a proper motion value.
; INPUTS:
;  infile    - multi-exten FITS file made by faintpm_sdss_verify or
;              faintpm_makefake or equivalent.
;  radot     - dRA/dts (mas/yr)
;  decdot    - dDec/dts (mas/yr)
;  tname     - name of the mean time of observation in FITS headers
;  outfile   - FITS file for stacked image
; OPTIONAL INPUTS:
;  tunit     - multiplier to get tname time into years (eg, 1.0/365.25
;              if tname time is in days); default 1.0
;  jpegfile  - if set, make a jpeg
;  scales    - input for jpeg making
;  nonlinearity  - input for jpeg
; BUGS:
;  - Currently just a grab-bag of code cut out of other code!
;  - Assumes the header for a single EXTEN is good for the stack.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_stack, infile,radot,decdot,tname,outfile, $
                   tunit=tunit,jpegfile=jpegfile,scales=scales, $
                   nonlinearity=nonlinearity
nepoch= sxpar(headfits(infile),'NEPOCH')
splog, 'nepoch',nepoch
timelist= fltarr(nepoch)
for tt=0L,nepoch-1L do begin
    thishdr= headfits(infile,exten=tt+1)
    timelist[tt]= sxpar(thishdr,tname)*tunit
    extast, thishdr,thisastr
    if (tt EQ 0) then begin
        thisastr= struct_addtags(thisastr,{naxis: [0,0]})
        astr= replicate(thisastr,nepoch)
    endif
    struct_assign, astr[tt],thisastr
    astr[tt].naxis= [sxpar(thishdr,'NAXIS1'),sxpar(thishdr,'NAXIS2')]
endfor
medt= median(timelist)
closest= (sort(abs(timelist-medt)))[0]
bigastr= astr[closest]
numerator= fltarr(bigastr.naxis)
denominator= numerator+1e-9     ; tiny
for ii=0,nepoch-1 do begin
    thisimage= mrdfits(infile,(2*ii+1))
    thisweight= mrdfits(infile,(2*ii+2))
    thisastr= astr[ii]
    thisastr.crval= thisastr.crval $
      -[radot,decdot]*(timelist[ii]-medt)
    smosaic_remap, thisimage,thisastr,bigastr,outimage,offset, $
      weight=thisweight,outweight=outweight,reflimits=reflim, $
      outlimits=outlim
    numerator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]]= $
      numerator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]] $
      +outweight[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]] $
      *outimage[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]]
    denominator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]]= $
      denominator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]] $
      +outweight[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]]
endfor
outimage= numerator/denominator
mkhdr, outhdr,outimage
putast, outhdr,bigastr
mwrfits, outimage,outfile,outhdr,/create
mwrfits, denominator,outfile
if (keyword_set(jpegname)) then begin
    subimage= outimage-median(outimage)
    nw_rgb_make, subimage,subimage,subimage, $
      name=jpegname,scales=scales,nonlinearity=nonlinearity,quality=95
endif
return
end
