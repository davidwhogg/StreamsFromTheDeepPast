;+
;   NAME:
;      read_hyades
;   PURPOSE:
;      read the hyades membership catalog (Perryman+98) and construct
;      the 'hip' object for it
;   CALLING SEQUENCE:
;      hip=read_hyades(savefilename=savefilename)
;   INPUT:
;      savefilename - filename for savefile
;   KEYWORDS:
;      tenpc        - only select stars within 10 pc from the center
;      single       - only select single stars
;   OUTPUT:
;      hyades-structure
;   REVISION HISTORY:
;      2009-09-14 - Written - Bovy (NYU)
;-
FUNCTION READ_HYADES, savefilename=savefilename, single=single, tenpc=tenpc

datafilename='../../data/HyadesMembership.dat'

;;First read the membership catalog
nfield= 33L
; most data types are long
fieldtypes= lonarr(nfield)+3
; some are double
fieldtypes[[20,21,22,25,30,31]]= 5
; some are string
fieldtypes[[3,6,9,12,15,18,23,24,27,28,29,32]]= 7
; give dumb names
fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
template= {version:1.0, datastart: 0, $
           delimiter: '|', missingvalue: 0.0, commentsymbol: '#', $
           fieldcount: nfield, fieldtypes: fieldtypes, $
           fieldnames: fieldnames, $
           fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
perrymandata= read_ascii(datafilename,num_records=nline,template=template)
splog, 'found ',n_elements(perrymandata.h01),' stars in the Perryman+98 Hyades catalog'
IF keyword_set(tenpc) AND keyword_set(single) THEN BEGIN
    hyadesmembers= where(perrymandata.h32 EQ '1' AND perrymandata.h30 LE 10 $
                         AND (perrymandata.h27 EQ string(0.0) OR perrymandata.h27 EQ '') $
                         AND (perrymandata.h28 EQ string(0.0) OR perrymandata.h28 EQ '')$
                         AND (perrymandata.h29 EQ string(0.0) OR perrymandata.h29 EQ ''))
ENDIF ELSE IF keyword_set(tenpc) THEN BEGIN
    hyadesmembers= where(perrymandata.h32 EQ '1' AND perrymandata.h30 LT 10)
ENDIF ELSE IF keyword_set(single) THEN BEGIN
    hyadesmembers= where(perrymandata.h32 EQ '1' $
                         AND (perrymandata.h27 EQ string(0.0) OR perrymandata.h27 EQ '') $
                         AND (perrymandata.h28 EQ string(0.0) OR perrymandata.h28 EQ '')$
                         AND (perrymandata.h29 EQ string(0.0) OR perrymandata.h29 EQ ''))
ENDIF ELSE BEGIN
    hyadesmembers= where(perrymandata.h32 EQ '1')
ENDELSE
nmembers= n_elements(hyadesmembers)
splog, strtrim(string(nmembers),2)+' Hyades members selected'

;;now read the Tycho data
tycho_savefile='../lsr2/tycho.sav'
IF file_test(tycho_savefile) THEN BEGIN
    restore, filename=tycho_savefile
ENDIF ELSE BEGIN
    tycdata=read_tycho(/read_raw)
    save, filename=tycho_savefile, tycdata
ENDELSE

;;Create the output object
str= {hyades, $
      hip     : 0L,  $          ; HIP number
      ra      : 0D,  $          ; right ascension (1991.25) (deg)
      dec     : 0D,  $          ; declination (1991.25) (deg)
      plx     : 0D,  $          ; parallax (mas)
      pmra    : 0D,  $        ; proper motion in RA direction (mas/yr)
      pmdec   : 0D,  $       ; proper motion in Dec direction (mas/yr)
      e_ra      : 0D,  $ ; standard error in RA (mas)
      e_dec     : 0D,  $ ; standard error in DEC (mas)
      e_plx     : 0D,  $ ; standard error in plx (mas)
      e_pmra    : 0D,  $ ; standard error in proper motion RA (mas/yr)
      e_pmdec   : 0D,  $ ; standard error in proper motion DEC (mas/yr)
      c_decra     : 0D,  $ ; correlation between dec and ra
      c_plxra     : 0D,  $ ; correlation between plx and ra
      c_plxdec    : 0D,  $ ; correlation between plx and dec
      c_pmrara    : 0D,  $ ; correlation between pmra and ra
      c_pmradec   : 0D,  $ ; correlation between pmra and dec
      c_pmraplx   : 0D,  $ ; correlation between pmra and plx
      c_pmdecra   : 0D,  $ ; correlation between pmdec and ra
      c_pmdecdec  : 0D,  $ ; correlation between pmdec and dec
      c_pmdecplx  : 0D,  $ ; correlation between pmdec and plx
      c_pmdecpmra : 0D,  $ ; correlation between pmdec and pmra
      vmag    : 0D,  $          ; Vmag from old reduction
      bvcolor : 0D,  $          ; BVcolor from new reduction
      vt      : 0D,  $          ; VT mag (mag)
      bt      : 0D,  $          ; BT mag (mag)
      e_vt    : 0D,  $          ; error in VT mag (mag)
      e_bt    : 0D,   $          ; error in BT mag (mag)
      l       : 0D,  $          ; Galactic longitude (deg)
      b       : 0D,  $          ; Galactic latitude (deg)
      radius  : 0D,  $    ; heliocentric radius (km) inferred from plx
      vra     : 0D,  $          ; velocity in RA direction (km/s)
      vdec    : 0D,  $          ; velocity in Dec direction (km/s)
      vrvdc   : dblarr(2,2), $ ; measurement covariance on vra,vdec ([km/s]^2)
      vl      : 0D,  $          ; velocity in l direction (km/s)
      vb      : 0D,  $          ; velocity in b direction (km/s)
      vlvbc   : dblarr(2,2), $ ; measurement covariance on vl,vb ([km/s]^2)
      nsm     : dblarr(3,2), $  ; non-square matrix for projection
      sm      : dblarr(3,3)  $  ; square matrix
     }
hyades = replicate(str,nmembers)
hyades.hip= perrymandata.h01[hyadesmembers]
hyades.plx= perrymandata.h20[hyadesmembers]
hyades.e_plx= perrymandata.h21[hyadesmembers]
;;Cross-correlate with the original Hipparcos catalog to get V
;;Read the original Hipparcos catalog
hip_orig=read_hipparcos(rootdir='/global/data/scr/jb2777/hipparcos',minsn=-100000000000000000,maxdist=1000000000);;hack to make sure almost the whole catalog is read
nhip_orig= n_elements(hip_orig.hip)
FOR indx= 0L, nhip_orig-1 DO BEGIN
    ii= (where(hip_orig[indx].hip EQ hyades.hip))[0]
    IF ii EQ -1 THEN BEGIN
        IF keyword_set(QA)  THEN splog, 'warning: no Original Hipparcos counterpart found'
        CONTINUE
    ENDIF
    hyades[ii].vmag= hip_orig[indx].vmag
ENDFOR
indx= where(hyades.vmag EQ 0.)
splog, strtrim(string(n_elements(indx)),2) +' stars without V magnitude'
;;Cross-correlate with the new Hipparcos catalog to get BVcolor and
;;all the rest
;;Read the Hipparcos catalog
hip_orig=read_hipparcos2(rootdir='/global/data/scr/jb2777/hipparcos/hipnewredux',/novmag,minsn=-100000000000000000,maxdist=1000000000);;hack to make sure almost the whole catalog is read
nhip_orig= n_elements(hip_orig.hip)
FOR indx= 0L, nhip_orig-1 DO BEGIN
    ii= (where(hip_orig[indx].hip EQ hyades.hip))[0]
    IF ii EQ -1 THEN BEGIN
        IF keyword_set(QA)  THEN splog, 'warning: no Original Hipparcos counterpart found'
        CONTINUE
    ENDIF
    hyades[ii].bvcolor= hip_orig[indx].bvcolor
    hyades[ii].ra=      hip_orig[indx].ra
    hyades[ii].dec=     hip_orig[indx].dec
    hyades[ii].pmra=    hip_orig[indx].pmra
    hyades[ii].pmdec=   hip_orig[indx].pmdec
    hyades[ii].e_ra=    hip_orig[indx].e_ra
    hyades[ii].e_dec=   hip_orig[indx].e_dec
    hyades[ii].e_pmra=  hip_orig[indx].e_pmra
    hyades[ii].e_pmdec= hip_orig[indx].e_pmdec
    hyades[ii].c_decra= hip_orig[indx].c_decra
    hyades[ii].c_plxra= hip_orig[indx].c_plxra
    hyades[ii].c_plxdec=hip_orig[indx].c_plxdec
    hyades[ii].c_pmrara=hip_orig[indx].c_pmrara
    hyades[ii].c_pmradec=hip_orig[indx].c_pmradec
    hyades[ii].c_pmraplx=hip_orig[indx].c_pmraplx
    hyades[ii].c_pmdecra=hip_orig[indx].c_pmdecra
    hyades[ii].c_pmdecdec=hip_orig[indx].c_pmdecdec
    hyades[ii].c_pmdecplx=hip_orig[indx].c_pmdecplx
    hyades[ii].c_pmdecpmra=hip_orig[indx].c_pmdecpmra
    hyades[ii].l=    hip_orig[indx].l
    hyades[ii].b=    hip_orig[indx].b
    hyades[ii].radius=hip_orig[indx].radius
    hyades[ii].vra=   hip_orig[indx].vra
    hyades[ii].vdec=  hip_orig[indx].vdec
    hyades[ii].vrvdc= hip_orig[indx].vrvdc
    hyades[ii].vl=     hip_orig[indx].vl
    hyades[ii].vb=     hip_orig[indx].vb
    hyades[ii].vlvbc=  hip_orig[indx].vlvbc
    hyades[ii].nsm=    hip_orig[indx].nsm
    hyades[ii].sm=     hip_orig[indx].sm

ENDFOR
indx= where(hyades.bvcolor EQ 0.)
splog, strtrim(string(n_elements(indx)),2) +' stars without B-V color'
;;Get VT and BT from Tycho-2
ntyc= n_elements(tycdata.h25)
FOR indx= 0L, ntyc-1 DO BEGIN
    IF tycdata[indx].h25 EQ 0 THEN CONTINUE
    ii= (where(tycdata[indx].h25 EQ hyades.hip))[0]
    IF ii EQ -1 THEN BEGIN
        CONTINUE
    ENDIF
    hyades[ii].vt= tycdata[indx].h21
    hyades[ii].e_vt= tycdata[indx].h22
    hyades[ii].bt= tycdata[indx].h19
    hyades[ii].e_bt= tycdata[indx].h20
ENDFOR

IF keyword_set(savefilename) THEN BEGIN
    splog, 'saving file '+savefilename
    save, hyades,filename=savefilename
ENDIF ELSE BEGIN
    splog, '*Not* saving'
ENDELSE

RETURN, hyades

END
