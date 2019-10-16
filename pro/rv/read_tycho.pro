;+
; NAME:
;   read_tycho
; PURPOSE:
;   read in, trim, and reformat the Tycho data - 
;   outputs ra, dec, pmra, pmdec, hip, + errors if hipparcos is not
;   set
;   if hipparcos is set, parallaxes etc. are taken from hipparcos and
;   a structure similar to the 'hip' from read_hipparcos is returned
;   in which Tycho proper motions are used
; USAGE:
;   tyc= read_tycho(/dehnen,/savefile,/giants,/hipparcos)
; INPUTS:
; OPTIONAL INPUTS:
;   minsn       minimum parallax S/N (defaults to 5)
;   maxdist     minimum radial distance (defaults to 200pc)
;   textfile    name of text file to write catalog; writes none default
;   binaryfile  name of binary file to write some output; writes none by default
;   nline       if set, read only the first nline lines of the catalog
; KEYWORDS:
;   read_raw   just read and return the raw data (i.e., the table as
;              it appears in the data files
;   dehnen     use the Dehnen & Binney (1998) MNRAS 298, 387 selection
;                by reading their HIP numbers from the streams/data
;                directory.  This option over-rides minsn, maxdist,
;                single.
;   giants     Use the giants sample (Dehnen 1998 AJ 115, 2384) 
;                by reading their HIP numbers from the streams/data
;                directory.  This option over-rides minsn, maxdist,
;                single.
;   hipparcos  Correlate with Hipparcos to get distances and the like
;   single     if set, reject multiple stars based on multiplicity flag
; OUTPUTS:
;   tyc        structure containing catalog and several derived quantities;
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
;   2008-11-3  Started - adapted from read_hipparcos - Bovy
;-
function read_tycho, dehnen=dehnen,minsn=minsn,maxdist=maxdist, $
  textfile=textfile,nline=nline,rootdir=rootdir,single=single,qa=qa, $
  binaryfile=binaryfile,savefile=savefile, giants=giants, $
                     hipparcos=hipparcos, read_raw=read_raw

; set defaults
  if NOT keyword_set(minsn) then minsn= 1D1
  snstr= strtrim(string(round(minsn)),2)
  if NOT keyword_set(maxdist) then maxdist= 200D ; pc
  diststr= strtrim(string(round(maxdist)),2)
  if keyword_set(single) then singlestr= '1' else singlestr= '0'
  if NOT keyword_set(rootdir) then rootdir='/global/data/hipparcos'


tycdir= rootdir+'/Tycho-2'

;;loop over data files
basename= tycdir+'/tyc2.dat.'
nfiles=20
nlines= lonarr(nfiles)
FOR ii=0L, nfiles-1 DO BEGIN
    datafilename= basename+string(ii,format='(I2.2)')
    gunzipname= datafilename+'.gz'
    cmd= 'gunzip '+gunzipname
    splog, 'Gunzipping '+gunzipname+'...'
    spawn, cmd

;now read data from this part of the catalog
; read data
  nfield= 32L
; most data types are double
  fieldtypes= lonarr(nfield)+5
; some are long
  fieldtypes[[6,7,12,21]]= 3
; some are string
  fieldtypes[[0,1,22,23,30]]= 7
; give dumb names
  fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
  template= {version:1.0, datastart: 0, $
             delimiter: '|', missingvalue: 0.0, commentsymbol: '#', $
             fieldcount: nfield, fieldtypes: fieldtypes, $
             fieldnames: fieldnames, $
             fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
  splog, 'Reading '+datafilename+'...'
  data_tmp= read_ascii(datafilename,num_records=nline,template=template)
  IF ii EQ 0 THEN BEGIN
      data018= REPLICATE(data_tmp,nfiles-1)
      data018[ii]= data_tmp
  ENDIF ELSE IF ii EQ nfiles-1 THEN BEGIN
      data19= data_tmp
  ENDIF ELSE BEGIN
      data018[ii]= data_tmp
  ENDELSE
  nlines[ii]= n_elements(data_tmp.h01)
  
  splog, 'Gzipping '+datafilename+' again...'
  cmd= 'gzip '+datafilename
  spawn, cmd
ENDFOR


;;Make silly structure
datastr= {tychosilly, $
       h00      : 0L,  $
       h01      : 0L,  $
       h02      : 0L,  $
       h03      : '',  $
       h04      : 0D,  $
       h05      : 0D,  $
       h06      : 0D,  $
       h07      : 0D,  $
       h08      : 0L,  $
       h09      : 0L,  $
       h10      : 0D,  $
       h11      : 0D,  $
       h12      : 0D,  $
       h13      : 0D,  $
       h14      : 0L,  $
       h15      : 0D,  $
       h16      : 0D,  $
       h17      : 0D,  $
       h18      : 0D,  $
       h19      : 0D,  $
       h20      : 0D,  $
       h21      : 0D,  $
       h22      : 0D,  $
       h23      : 0L,  $
       h24      : '',  $
       h25      : 0L,  $
       h26      : '',  $
       h27      : 0D,  $
       h28      : 0D,  $
       h29      : 0D,  $
       h30      : 0D,  $
       h31      : 0D,  $
       h32      : 0D,  $
       h33      : '',  $
       h34      : 0D  $
}
data = replicate(datastr,total(nlines))

;;and copy stuff into it
nlines= total(nlines,/cumulative)
FOR ii=0L, nfiles-2 DO BEGIN
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h00= long(strmid(data018[ii].h00,0,4))
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h01= long(strmid(data018[ii].h00,5,5))
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h02= long(strmid(data018[ii].h00,11,1))
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h03= data018[ii].h01
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h04= data018[ii].h02
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h05= data018[ii].h03
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h06= data018[ii].h04
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h07= data018[ii].h05
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h08= data018[ii].h06
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h09= data018[ii].h07
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h10= data018[ii].h08
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h11= data018[ii].h09
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h12= data018[ii].h10
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h13= data018[ii].h11
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h14= data018[ii].h12
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h15= data018[ii].h13
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h16= data018[ii].h14
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h17= data018[ii].h15
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h18= data018[ii].h16
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h19= data018[ii].h17
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h20= data018[ii].h18
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h21= data018[ii].h19
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h22= data018[ii].h20
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h23= data018[ii].h21
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h24= data018[ii].h22
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h25= long(strmid(data018[ii].h23,0,6))
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h26= strmid(data018[ii].h23,6,3)
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h27= data018[ii].h24
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h28= data018[ii].h25
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h29= data018[ii].h26
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h30= data018[ii].h27
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h31= data018[ii].h28
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h32= data018[ii].h29
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h33= data018[ii].h30
    data[(ii NE 0) * nlines[(ii NE 0)*(ii-1)]:nlines[ii]-1].h34= data018[ii].h31
ENDFOR
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h00= long(strmid(data19.h00,0,4))
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h01= long(strmid(data19.h00,5,5))
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h02= long(strmid(data19.h00,11,1))
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h03= data19.h01
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h04= data19.h02
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h05= data19.h03
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h06= data19.h04
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h07= data19.h05
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h08= data19.h06
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h09= data19.h07
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h10= data19.h08
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h11= data19.h09
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h12= data19.h10
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h13= data19.h11
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h14= data19.h12
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h15= data19.h13
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h16= data19.h14
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h17= data19.h15
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h18= data19.h16
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h19= data19.h17
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h20= data19.h18
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h21= data19.h19
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h22= data19.h20
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h23= data19.h21
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h24= data19.h22
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h25= long(strmid(data19.h23,0,6))
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h26= strmid(data19.h23,6,3)
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h27= data19.h24
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h28= data19.h25
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h29= data19.h26
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h30= data19.h27
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h31= data19.h28
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h32= data19.h29
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h33= data19.h30
data[nlines[nfiles-2]:nlines[nfiles-1]-1].h34= data19.h31

splog, 'found ',n_elements(data.h01),' stars in the Tycho catalog'


IF keyword_set(read_raw) THEN RETURN, data


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Now read Hipparcos, THIS IS BASICALLY THE READ_HIPPARCOS CODE
; read data
  nfield= 78L
; most data types are double
  fieldtypes= lonarr(nfield)+5
; some are long
  fieldtypes[[1,6,29,31,47,57,58,63,71]]= 3
; some are string
  fieldtypes[[0,2,3,4,10,36,39,42,43,48,52,53,54,55,56,59,60,61,62, $
    68,69,70,72,73,74,76,77]]= 7
; give dumb names
  fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
  template= {version:1.0, datastart: 0, $
             delimiter: '|', missingvalue: 0.0, commentsymbol: '#', $
             fieldcount: nfield, fieldtypes: fieldtypes, $
             fieldnames: fieldnames, $
             fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
  hipdata= read_ascii(rootdir+'/hip_main.dat',num_records=nline,template=template)
  splog, 'found ',n_elements(hipdata.h01),' stars in the Hipparcos catalog'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



; either read dehnen HIPs and select them...
  if keyword_set(dehnen) then begin
      readcol, getenv('HOME')+'/streams/data/dehnen_binney_sample1.id', $
        hipsamp,format='L'
      nsample= n_elements(hipsamp)
      index= lonarr(nsample)
      inhip= bytarr(nsample)
      for jj=0L,nsample-1 do begin
          index[jj]= (where(data.h25 EQ hipsamp[jj]))[0]
          IF (index[jj] EQ -1) THEN BEGIN;;If it's not in Tycho, look for it in Hipparcos
              splog, 'Warning: no Tycho counterpart found'
              index[jj]= (where(hipdata.h01 EQ hipsamp[jj]))[0]
              inhip[jj]= 1
          ENDIF
      endfor
      splog, strtrim(string(total(inhip)),2)+' withouth Tycho entry found'
      jj= where(index NE -1,nsample)
      index= index[jj]
      splog, 'found ',nsample,' stars in the Dehnen & Binney (1998) sample'
  endif else if keyword_set(giants) then begin
;; Or read the Dehnen Giant's Hips and select them...
      readcol, getenv('HOME')+'/streams/data/dehnen_binney_giants.id', $
        hipsamp,format='L'
      nsample= n_elements(hipsamp)
      index= lonarr(nsample)
      inhip= bytarr(nsample)
      for jj=0L,nsample-1 do begin
          index[jj]= (where(data.h25 EQ hipsamp[jj]))[0]
          IF (index[jj] EQ -1) THEN BEGIN;;If it's not in Tycho, look for it in Hipparcos
              splog, 'Warning: no Tycho counterpart found'
              index[jj]= (where(hipdata.h01 EQ hipsamp[jj]))[0]
              inhip[jj]= 1
          ENDIF
      endfor
      splog, strtrim(string(long(total(inhip))),2)+' without Tycho entry found'
      jj= where(index NE -1,nsample)
      index= index[jj]
      splog, 'found ',nsample,' stars in the Dehnen (1998) Giants sample'
  endif else begin
;THIS DOES NOT APPLY, NO PARALLAXES IN THE TYCHO CATALOG, DON"T USE
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
  str= {tycho, $
        hip     : 0L,  $ ; HIP number
        ra      : 0D,  $ ; right ascension (2000.0) (deg)
        dec     : 0D,  $ ; declination (2000.0) (deg)
        e_ra    : 0D,  $ ; error in right ascension (2000.0) (deg)
        e_dec   : 0D,  $ ; error in declination (2000.0) (deg)
        pmra    : 0D,  $ ; proper motion in RA direction (mas/yr)
        pmdec   : 0D,  $ ; proper motion in Dec direction (mas/yr)
        e_pmra  : 0D,  $ ; error in proper motion in RA direction (mas/yr)
        e_pmdec : 0D  $ ; error in proper motion in Dec direction (mas/yr)
       }
  tyc = replicate(str,nsample)

; load data structure
  tyc.hip=     data[index].h25
  tyc.ra=      data[index].h04
  tyc.dec=     data[index].h05
  tyc.e_ra=    data[index].h08
  tyc.e_dec=   data[index].h09
  tyc.pmra=    data[index].h06
  tyc.pmdec=   data[index].h07
  tyc.e_pmra=  data[index].h10
  tyc.e_pmdec= data[index].h11

;;Ignores correlations, fix this, deal with these 25 stars that are
;;not in Tycho
hipnotyc= where(inhip EQ 1)
tyc[hipnotyc].hip=     hipdata.h01[index[hipnotyc]]
tyc[hipnotyc].ra=      hipdata.h08[index[hipnotyc]]
tyc[hipnotyc].dec=     hipdata.h09[index[hipnotyc]]
tyc[hipnotyc].e_ra=    hipdata.h14[index[hipnotyc]]
tyc[hipnotyc].e_dec=   hipdata.h15[index[hipnotyc]]
tyc[hipnotyc].pmra=    hipdata.h12[index[hipnotyc]]
tyc[hipnotyc].pmdec=   hipdata.h13[index[hipnotyc]]
tyc[hipnotyc].e_pmra=  hipdata.h17[index[hipnotyc]]
tyc[hipnotyc].e_pmdec= hipdata.h18[index[hipnotyc]]


save, filename='tyc.sav', tyc



IF ~keyword_set(hipparcos) THEN RETURN, tyc

indextyc= index
data= hipdata

; either read dehnen HIPs and select them...
  if keyword_set(dehnen) then begin
      readcol, getenv('HOME')+'/streams/data/dehnen_binney_sample1.id', $
        hipsamp,format='L'
      nsample= n_elements(hipsamp)
      index= lonarr(nsample)
      for jj=0L,nsample-1 do begin
          index[jj]= (where(data.h01 EQ hipsamp[jj]))[0]
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
          index[jj]= (where(data.h01 EQ hipsamp[jj]))[0]
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
  str= {hipparcos, $
        hip     : 0L,  $ ; HIP number
        ra      : 0D,  $ ; right ascension (1991.25) (deg)
        dec     : 0D,  $ ; declination (1991.25) (deg)
        plx     : 0D,  $ ; parallax (mas)
        pmra    : 0D,  $ ; proper motion in RA direction (mas/yr)
        pmdec   : 0D,  $ ; proper motion in Dec direction (mas/yr)
      e_ra      : 0D,  $ ; error on RA (mas)
      e_dec     : 0D,  $ ; error on DEC (mas)
      e_plx     : 0D,  $ ; error on PLX (mas)
      e_pmra    : 0D,  $ ; error on pmRA (mas/yr)
      e_pmdec   : 0D,  $ ; error on pmDEC (mas/ytr)
        vmag    : 0D,  $ ; apparent V magnitude (mag)
        bvcolor : 0D,  $ ; B-V color (mag)
        flags   : 0L,  $ ; flag integer; 1 for LMS, 2 for UMS, 4 for RG
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
  hip.hip=     tyc.hip
  hip.ra=      tyc.ra
  hip.dec=     tyc.dec
  hip.plx=     data.h11[index]
  hip.pmra=    tyc.pmra
  hip.pmdec=   tyc.pmdec
  hip.vmag=    data.h05[index]
  hip.bvcolor= data.h37[index]

  hip.e_ra    = tyc.e_ra
  hip.e_dec   = tyc.e_dec
  hip.e_plx   = data.h16[index]
  hip.e_pmra  = tyc.e_pmra
  hip.e_pmdec = tyc.e_pmdec
  
; get Galactic quantities
  epoch= 2000.0
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
  radius_err= radius*data.h16[index]/hip.plx     ; km
  pmra_err = tyc.e_pmra*invsec_for_unit_pm  ; 1/s
  pmdec_err= tyc.e_pmdec*invsec_for_unit_pm  ; 1/s
  hip.vrvdc[0,0]= pmra^2 *radius_err^2 $
    +radius^2*pmra_err^2; $
;    +2d0*radius*pmra *radius_err*pmra_err *data.h24[index]
  hip.vrvdc[1,1]= pmdec^2*radius_err^2 $
    +radius^2*pmdec_err^2; $
;    +2d0*radius*pmdec*radius_err*pmdec_err*data.h27[index]
  hip.vrvdc[0,1]= pmra*pmdec*radius_err^2; $
;    +radius*pmra *radius_err*pmdec_err*data.h27[index] $
;    +radius*pmdec*radius_err*pmra_err *data.h24[index] $
;    +radius^2*pmra_err*pmdec_err*data.h28[index]
  hip.vrvdc[1,0]= hip.vrvdc[0,1]

; rotate the covariance matrices into the l,b space
  zero= dblarr(nsample)
  unity= zero+1
  radec_to_lb, tyc.ra,tyc.dec,unity,zero,ll,bb,cosine,sine,epoch
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


  tyc= hip


  if keyword_set(savefile) then begin
      savefilename= rootdir+'/hip-'+snstr+'-'+diststr+'-'+singlestr+'.sav'
      if keyword_set(giants) then savefilename= '~/streams/pro/rv/tyc-dehnen-giants.sav'
      if keyword_set(dehnen) then savefilename= '~/streams/pro/rv/tyc-dehnen.sav'
;      if keyword_set(giants) then savefilename= rootdir+'/tyc-dehnen-giants.sav'
;      if keyword_set(dehnen) then savefilename= rootdir+'/tyc-dehnen.sav'
      splog, 'saving file '+savefilename
      save, tyc,filename=savefilename
  endif


  return, tyc
end


