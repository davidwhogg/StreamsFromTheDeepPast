;+
; NAME:
;   read_hipparcos2
; PURPOSE:
;   read in, trim, and reformat the hipparcos data
;   BOVY: reads in Dehnen sample, might not read in rest of the sample correctly
; USAGE:
;   hip= read_hipparcos(/dehnen,/savefile)
;   hip= read_hipparcos(/giants,/savefile)
; INPUTS:
; OPTIONAL INPUTS:
;   minsn       minimum parallax S/N (defaults to 5)
;   maxdist     minimum radial distance (defaults to 200pc)
;   textfile    name of text file to write catalog; writes none default
;   binaryfile  name of binary file to write some output; writes none by default
;   nline       if set, read only the first nline lines of the catalog
; KEYWORDS:
;   dehnen     use the Dehnen & Binney (1998) MNRAS 298, 387 selection
;                by reading their HIP numbers from the streams/data
;                directory.  This option over-rides minsn, maxdist,
;                single.
;   giants     Use the giants sample (Dehnen 1998 AJ 115, 2384) 
;                by reading their HIP numbers from the streams/data
;                directory.  This option over-rides minsn, maxdist,
;                single.
;   single     if set, reject multiple stars based on multiplicity
;   flag (DEPRECATED)
;   qa         if set, make qa plots hip*.eps (DEPRECATED)
;   savefile   if set, save a savefile in the rootdir (DEPRECATED)
;   formal     if set, use the formal errors
;   novmag     don't bother getting vmag from the original Hipparcos catalog
; OUTPUTS:
;   hip        structure containing catalog and several derived quantities;
;                for more information, "use the source, Luke".
; COMMENTS:
;   Covariance matrix calculation tested by Hogg on 2002-02-20 by inputting
;     fake data made with the hipparcos covariance matrix, and checking output
;     against analytical calculation.
; BUGS:
;   Does not apply Lutz-Kelker corrections.
;   Inverts some matrices which, being unitary, could just be
;     transposed.
;   Filenames and paths hard-wired.
; DEPENDENCIES:
;   idlutils
; REVISION HISTORY:
;   2001-08-08  written - Hogg (NYU)
;   2001-10-14  added QA plots, flags, non-square matrices - Hogg
;   2002-02-20  added covariance matrices - Hogg
;   2003-03-11  read Dehnen & Binney HIP numbers - Hogg
;   2008-07-11  read Dehnen Giants HIP numbers - Bovy
;   2008-11-07  read in new reduction - Bovy
;   2008-11-26  Added observational errors to the output struct - Bovy
;   2009-09-14  Added novmag keyword - Bovy
;-
function read_hipparcos2, dehnen=dehnen,minsn=minsn,maxdist=maxdist, $
                          textfile=textfile,nline=nline,rootdir=rootdir,$
                          savefile=savefile, $
                          giants=giants, $
                          formal=formal, novmag=novmag


; set defaults
  if NOT keyword_set(minsn) then minsn= 1D1
  snstr= strtrim(string(round(minsn)),2)
  if NOT keyword_set(maxdist) then maxdist= 200D ; pc
  diststr= strtrim(string(round(maxdist)),2)
  if keyword_set(single) then singlestr= '1' else singlestr= '0'
  if NOT keyword_set(rootdir) then rootdir='/global/data/hipparcos/hipnewredux'

; read data
  nfield= 41L
; most data types are double
  fieldtypes= lonarr(nfield)+5
; some are long
  fieldtypes[[0,1,2,3,14,16,18,22]]= 3
; none are string
;  fieldtypes[[0,2,3,4,10,36,39,42,43,48,52,53,54,55,56,59,60,61,62, $
;    68,69,70,72,73,74,76,77]]= 7
; give dumb names
  fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
  template= {version:1.0, datastart: 0, $
             delimiter: ' ', missingvalue: 0.0, commentsymbol: '#', $
             fieldcount: nfield, fieldtypes: fieldtypes, $
             fieldnames: fieldnames, $
             fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
  data= read_ascii(rootdir+'/hip2.dat',num_records=nline,template=template)
  splog, 'found ',n_elements(data.h01),' stars in the Hipparcos catalog'

; either read dehnen HIPs and select them...
  if keyword_set(dehnen) then begin
      readcol, getenv('HOME')+'/streams/data/dehnen_binney_sample1.id', $
        hipsamp,format='L'
      nsample= n_elements(hipsamp)
      index= lonarr(nsample)
      for jj=0L,nsample-1 do begin
          index[jj]= (where(data.h00 EQ hipsamp[jj]))[0]
      endfor
      jj= where(index NE -1,nsample)
      index= index[jj]
      splog, 'found ',nsample,' stars in the Dehnen & Binney (1998) sample'
  endif else if keyword_set(giants) then begin
;; Or read the Dehnen Giant's Hips and select them...
      readcol, getenv('HOME')+'/streams/data/dehnen_binney_giants.id', $
        hipsamp,format='L'
      nsample= n_elements(hipsamp)
      index= lonarr(nsample)
      for jj=0L,nsample-1 do begin
          index[jj]= (where(data.h00 EQ hipsamp[jj]))[0]
      endfor
      jj= where(index NE -1,nsample)
      index= index[jj]
      splog, 'found ',nsample,' stars in the Dehnen (1998) Giants sample'
  endif else begin
; ...or cut on parallax s/n and distance and, if 'single', multiplicity
      sn= data.h11/data.h16
      sncut= sn GT minsn
      distance= 1000.0/data.h11
      distcut= distance LT maxdist
      if keyword_set(single) then singlecut= data.h59 EQ '' else singlecut= 1B
      index= where(sncut AND distcut AND singlecut,nsample)
  endelse
  if nsample LT 1 then begin
      splog, 'ERROR: nothing makes sample cuts'
      return, 0
  endif

; make data structure
  str= {hipparcos2, $
        hip     : 0L,  $ ; HIP number
        ra      : 0D,  $ ; right ascension (1991.25) (deg)
        dec     : 0D,  $ ; declination (1991.25) (deg)
        plx     : 0D,  $ ; parallax (mas)
        pmra    : 0D,  $ ; proper motion in RA direction (mas/yr)
        pmdec   : 0D,  $ ; proper motion in Dec direction (mas/yr)
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
        vmag    : 0D,  $ ; apparent V magnitude (mag) NOT IN NEW REDUCTION, from OLD
        bvcolor : 0D,  $ ; B-V color (mag)
;        flags   : 0L,  $ ; flag integer; 1 for LMS, 2 for UMS, 4 for
;        RG (NOT IN NEW REDUCTION)
        l       : 0D,  $ ; Galactic longitude (deg)
        b       : 0D,  $ ; Galactic latitude (deg)
        radius  : 0D,  $ ; heliocentric radius (km) inferred from plx
        vra     : 0D,  $ ; velocity in RA direction (km/s)
        vdec    : 0D,  $ ; velocity in Dec direction (km/s)
        vrvdc   : dblarr(2,2), $ ; measurement covariance on vra,vdec ([km/s]^2)
        vl      : 0D,  $ ; velocity in l direction (km/s)
        vb      : 0D,  $ ; velocity in b direction (km/s)
        vlvbc   : dblarr(2,2), $ ; measurement covariance on vl,vb ([km/s]^2)
        nsm     : dblarr(3,2), $ ; non-square matrix for projection
        sm      : dblarr(3,3)  $ ; square matrix
       }
  hip = replicate(str,nsample)

; load data structure
  hip.hip=     data.h00[index]
  hip.ra=      data.h04[index]/!DPI*180.
  hip.dec=     data.h05[index]/!DPI*180.
  hip.plx=     data.h06[index]
  hip.pmra=    data.h07[index]
  hip.pmdec=   data.h08[index]
  hip.bvcolor= data.h23[index]

; Find the errors
FOR ii=0L, nsample-1 DO BEGIN
    U= dblarr(5,5)
    U[0,0]= data.h26[index[ii]];1
    U[1,0]= data.h27[index[ii]];2
    U[1,1]= data.h28[index[ii]];3
    U[2,0]= data.h29[index[ii]];4
    U[2,1]= data.h30[index[ii]];5
    U[2,2]= data.h31[index[ii]];6
    U[3,0]= data.h32[index[ii]];7
    U[3,1]= data.h33[index[ii]];8
    U[3,2]= data.h34[index[ii]];9
    U[3,3]= data.h35[index[ii]];10
    U[4,0]= data.h36[index[ii]];11
    U[4,1]= data.h37[index[ii]];12
    U[4,2]= data.h38[index[ii]];13
    U[4,3]= data.h39[index[ii]];14
    U[4,4]= data.h40[index[ii]];15
    ;;C^{-1}= U\T U
    C= transpose(U)##U
    C= invert(C,/double);;This is the covariance matrix
    ;;Now put everything in place
    IF keyword_set(formal) THEN BEGIN
        hip[ii].e_ra        = data.h09[index[ii]]
        hip[ii].e_dec       = data.h10[index[ii]]
        hip[ii].e_plx       = data.h11[index[ii]] 
        hip[ii].e_pmra      = data.h12[index[ii]]
        hip[ii].e_pmdec     = data.h13[index[ii]]
    ENDIF ELSE BEGIN
        hip[ii].e_ra        = sqrt(C[0,0])
        hip[ii].e_dec       = sqrt(C[1,1])
        hip[ii].e_plx       = sqrt(C[2,2])
        hip[ii].e_pmra      = sqrt(C[3,3])
        hip[ii].e_pmdec     = sqrt(C[4,4])
    ENDELSE
    hip[ii].c_decra     = C[0,1]/(hip[ii].e_ra*hip[ii].e_dec)
    hip[ii].c_plxra     = C[0,2]/(hip[ii].e_plx*hip[ii].e_ra)
    hip[ii].c_plxdec    = C[2,1]/(hip[ii].e_plx*hip[ii].e_dec) 
    hip[ii].c_pmrara    = C[3,0]/(hip[ii].e_pmra*hip[ii].e_ra) 
    hip[ii].c_pmradec   = C[3,1]/(hip[ii].e_pmra*hip[ii].e_dec)
    hip[ii].c_pmraplx   = C[3,2]/(hip[ii].e_pmra*hip[ii].e_plx)
    hip[ii].c_pmdecra   = C[4,0]/(hip[ii].e_pmdec*hip[ii].e_ra)
    hip[ii].c_pmdecdec  = C[4,1]/(hip[ii].e_pmdec*hip[ii].e_dec)
    hip[ii].c_pmdecplx  = C[4,2]/(hip[ii].e_pmdec*hip[ii].e_plx)
    hip[ii].c_pmdecpmra = C[4,3]/(hip[ii].e_pmdec*hip[ii].e_pmra)
ENDFOR


; get Galactic quantities
  epoch= 1991.25
  radec_to_lb, hip.ra,hip.dec,hip.pmra,hip.pmdec,ll,bb,pmll,pmbb,epoch
  hip.l= ll
  hip.b= bb
  km_for_unit_plx= 3.08568D16
  radius= km_for_unit_plx/hip.plx  ; km
  hip.radius= radius               ; km
  invsec_for_unit_pm= !DPI*1.0d0/1.80D2/3.6D6/3.15569D7
  hip.vra = hip.pmra *radius*invsec_for_unit_pm  ; km/s
  hip.vdec= hip.pmdec*radius*invsec_for_unit_pm  ; km/s
  hip.vl  = pmll     *radius*invsec_for_unit_pm  ; km/s
  hip.vb  = pmbb     *radius*invsec_for_unit_pm  ; km/s

; construct the velocity covariance matrices in ra,dec space
  pmra = hip.pmra *invsec_for_unit_pm            ; 1/s
  pmdec= hip.pmdec*invsec_for_unit_pm            ; 1/s
  radius_err= -radius*hip.e_plx/hip.plx     ; km
  pmra_err = hip.e_pmra*invsec_for_unit_pm  ; 1/s
  pmdec_err= hip.e_pmdec*invsec_for_unit_pm  ; 1/s
  hip.vrvdc[0,0]= pmra^2 *radius_err^2 $
    +radius^2*pmra_err^2 $
    +2d0*radius*pmra *radius_err*pmra_err *hip.c_pmraplx
  hip.vrvdc[1,1]= pmdec^2*radius_err^2 $
    +radius^2*pmdec_err^2 $
    +2d0*radius*pmdec*radius_err*pmdec_err*hip.c_pmdecplx
  hip.vrvdc[0,1]= pmra*pmdec*radius_err^2 $
    +radius*pmra *radius_err*pmdec_err*hip.c_pmdecplx $
    +radius*pmdec*radius_err*pmra_err *hip.c_pmraplx $
    +radius^2*pmra_err*pmdec_err*hip.c_pmdecpmra
  hip.vrvdc[1,0]= hip.vrvdc[0,1]

; rotate the covariance matrices into the l,b space
  zero= dblarr(nsample)
  unity= zero+1
  radec_to_lb, hip.ra,hip.dec,unity,zero,ll,bb,cosine,sine,epoch
  rotmax= dblarr(2,2,nsample)
  rotmax[0,0,*]= cosine
  rotmax[0,1,*]= sine
  rotmax[1,0,*]= -sine
  rotmax[1,1,*]= cosine
  for i=0L,nsample-1L do begin
      u= rotmax[*,*,i]
      tmp= u##hip[i].vrvdc##transpose(u)
; re-symmetrize covariance matrix explicitly (to remove numerical noise)
      hip[i].vlvbc= 0.5d0*(tmp+transpose(tmp))
  endfor

; make and load non-square matrices for ex_max_proj
  ll= hip.l*!DPI/1.80D2
  bb= hip.b*!DPI/1.80D2
  tmpmat= dblarr(3,3)
  for i=0L,nsample-1L do begin
      tmpmat[0,0]=  cos(bb[i])*cos(ll[i])
      tmpmat[1,0]= -sin(ll[i])
      tmpmat[2,0]= -sin(bb[i])*cos(ll[i])
      tmpmat[0,1]=  cos(bb[i])*sin(ll[i])
      tmpmat[1,1]=  cos(ll[i])
      tmpmat[2,1]= -sin(bb[i])*sin(ll[i])
      tmpmat[0,2]=  sin(bb[i])
      tmpmat[1,2]=  0
      tmpmat[2,2]=  cos(bb[i])
      tmpmat= transpose(tmpmat)
      hip[i].sm= tmpmat
      hip[i].nsm= tmpmat[*,1:2]
  endfor



IF ~keyword_set(novmag) THEN BEGIN
    ;;Get Vmag from old reduction
    hip_old= read_hipparcos(dehnen=dehnen,giants=giants,minsn=minsn,maxdist=maxdist)
    ;;These should be in the right order, but let's cross-correlate just
    ;;to make sure
    FOR ii= 0L, n_elements(hip.hip)-1 DO BEGIN
        old_indx= (where(hip_old.hip EQ hip[ii].hip))[0]
        IF old_indx EQ -1 THEN BEGIN
            CONTINUE
        ENDIF
        hip[ii].vmag= hip_old[old_indx].vmag
    ENDFOR
ENDIF


if keyword_set(savefile) then begin
    savefilename= rootdir+'/hip-'+snstr+'-'+diststr+'-'+singlestr+'.sav'
    if keyword_set(giants) then savefilename= rootdir+'/hip2-dehnen-giants.sav'
    if keyword_set(dehnen) then savefilename= rootdir+'/hip2-dehnen.sav'
    splog, 'saving file '+savefilename
    save, hip,filename=savefilename
endif

return, hip
end
