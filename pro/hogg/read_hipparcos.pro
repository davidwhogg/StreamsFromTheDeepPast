;+
; NAME:
;   read_hipparcos
; PURPOSE:
;   read in, trim, and reformat the hipparcos data
; USAGE:
;   hip= read_hipparcos(/single,qafile='hip.eps',textfile='hip.txt',$
;        binaryfile='hip.bin',/savefile)
;   hip= read_hipparcos(/dehnen,/savefile,/qa)
;   hip= read_hipparcos(/giants,/savefile,/qa)
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
;   single     if set, reject multiple stars based on multiplicity flag
;   qa         if set, make qa plots hip*.eps
;   savefile   if set, save a savefile in the rootdir
; OUTPUTS:
;   hip        structure containing catalog and several derived quantities;
;                for more information, "use the source, Luke".
; COMMENTS:
;   Covariance matrix calculation tested by Hogg on 2002-02-20 by inputting
;     fake data made with the hipparcos covariance matrix, and checking output
;     against analytical calculation.
; BUGS:
;   Cannot handle the HIPPARCOS data in .gz form, which means we have
;     to keep it on disk uncompressed!  If you fix this, fix download
;     script (in bin directory) too.
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
;   2008-11-26  added observational errors to the ouput structure - Bovy
;-
function read_hipparcos, dehnen=dehnen,minsn=minsn,maxdist=maxdist, $
  textfile=textfile,nline=nline,rootdir=rootdir,single=single,qa=qa, $
  binaryfile=binaryfile,savefile=savefile, giants=giants

; set defaults
  if NOT keyword_set(minsn) then minsn= 1D1
  snstr= strtrim(string(round(minsn)),2)
  if NOT keyword_set(maxdist) then maxdist= 200D ; pc
  diststr= strtrim(string(round(maxdist)),2)
  if keyword_set(single) then singlestr= '1' else singlestr= '0'
  if NOT keyword_set(rootdir) then rootdir='/global/data/hipparcos'

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
  data= read_ascii(rootdir+'/hip_main.dat',num_records=nline,template=template)
  splog, 'found ',n_elements(data.h01),' stars in the Hipparcos catalog'

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
  splog, strtrim(string(nsample),2)+' stars selected'

; make data structure
  str= {hipparcos, $
        hip       : 0L,  $ ; HIP number
        ra        : 0D,  $ ; right ascension (1991.25) (deg)
        dec       : 0D,  $ ; declination (1991.25) (deg)
        plx       : 0D,  $ ; parallax (mas)
        pmra      : 0D,  $ ; proper motion in RA direction (mas/yr)
        pmdec     : 0D,  $ ; proper motion in Dec direction (mas/yr)
      e_ra        : 0D,  $ ; standard error in ra (mas)
      e_dec       : 0D,  $ ; standard error in dec (mas)
      e_plx       : 0D,  $ ; standard error in plx (mas)
      e_pmra      : 0D,  $ ; standard error in pmra (mas/yr)
      e_pmdec     : 0D,  $ ; standard error in pmdec (mas/yr)
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
        vmag      : 0D,  $ ; apparent V magnitude (mag)
        bvcolor   : 0D,  $ ; B-V color (mag)
        flags     : 0L,  $ ; flag integer; 1 for LMS, 2 for UMS, 4 for RG
        l         : 0D,  $ ; Galactic longitude (deg)
        b         : 0D,  $ ; Galactic latitude (deg)
        radius    : 0D,  $ ; heliocentric radius (km) inferred from plx
        vra       : 0D,  $ ; velocity in RA direction (km/s)
        vdec      : 0D,  $ ; velocity in Dec direction (km/s)
        vrvdc     : dblarr(2,2), $ ; measurement covariance on vra,vdec ([km/s]^2)
        vl        : 0D,  $ ; velocity in l direction (km/s)
        vb        : 0D,  $ ; velocity in b direction (km/s)
        vlvbc     : dblarr(2,2), $ ; measurement covariance on vl,vb ([km/s]^2)
        nsm       : dblarr(3,2), $ ; non-square matrix for projection
        sm        : dblarr(3,3)  $ ; square matrix
       }
  hip = replicate(str,nsample)

; load data structure
  hip.hip=     data.h01[index]
  hip.ra=      data.h08[index]
  hip.dec=     data.h09[index]
  hip.plx=     data.h11[index]
  hip.pmra=    data.h12[index]
  hip.pmdec=   data.h13[index]
  hip.vmag=    data.h05[index]
  hip.bvcolor= data.h37[index]

  hip.e_ra        = data.h14[index]
  hip.e_dec       = data.h15[index]
  hip.e_plx       = data.h16[index]
  hip.e_pmra      = data.h17[index]
  hip.e_pmdec     = data.h18[index]
  hip.c_decra     = data.h19[index]
  hip.c_plxra     = data.h20[index]
  hip.c_plxdec    = data.h21[index]
  hip.c_pmrara    = data.h22[index]
  hip.c_pmradec   = data.h23[index]
  hip.c_pmraplx   = data.h24[index]
  hip.c_pmdecra   = data.h25[index]
  hip.c_pmdecdec  = data.h26[index]
  hip.c_pmdecplx  = data.h27[index]
  hip.c_pmdecpmra = data.h28[index]
  
; make subsamples based on absolute magnitude and color and set flags
  absvmag= hip.vmag+5D0*alog10(hip.plx/100.0)
  absvmag_cut= 4.50D
  bvcolor_cut= 0.72D
  LMS_index= where(absvmag LT absvmag_cut,nLMS)
  if nLMS GT 0 then hip[LMS_index].flags= hip[LMS_index].flags+1
  UMS_index= where(absvmag GE absvmag_cut AND hip.bvcolor LT bvcolor_cut,nUMS)
  if nUMS GT 0 then hip[UMS_index].flags= hip[UMS_index].flags+2
  RG_index= where(absvmag GE absvmag_cut AND hip.bvcolor GE bvcolor_cut,nRG)
  if nRG GT 0 then hip[RG_index].flags= hip[RG_index].flags+4

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
  radius_err= -radius*data.h16[index]/hip.plx     ; km
  pmra_err = data.h17[index]*invsec_for_unit_pm  ; 1/s
  pmdec_err= data.h18[index]*invsec_for_unit_pm  ; 1/s
  hip.vrvdc[0,0]= pmra^2 *radius_err^2 $
    +radius^2*pmra_err^2 $
    +2d0*radius*pmra *radius_err*pmra_err *data.h24[index]
  hip.vrvdc[1,1]= pmdec^2*radius_err^2 $
    +radius^2*pmdec_err^2 $
    +2d0*radius*pmdec*radius_err*pmdec_err*data.h27[index]
  hip.vrvdc[0,1]= pmra*pmdec*radius_err^2 $
    +radius*pmra *radius_err*pmdec_err*data.h27[index] $
    +radius*pmdec*radius_err*pmra_err *data.h24[index] $
    +radius^2*pmra_err*pmdec_err*data.h28[index]
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

; make qa plot (if asked to do so)
  if keyword_set(qa) then begin
      set_plot, "PS"
      !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
      !P.THICK= 4.0
      !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
      !P.CHARSIZE= 1.0
      !P.PSYM= 0
      !P.TITLE= ''
      !P.MULTI= 0
      !X.STYLE= 1
      !X.TITLE= '(B-V)  (mag)'
      !X.MARGIN= [6,1]
      !X.OMARGIN= [0,0]
      !X.RANGE= [-0.3,1.6]
      !Y.STYLE= 1
      !Y.TITLE= 'M!dV!n  (mag)'
      !Y.MARGIN= 0.5*!X.MARGIN
      !Y.OMARGIN= 0.5*!X.OMARGIN
      !Y.RANGE= [12,-4.5]
      xsize= 5.0 & ysize= 7.5
      device, file='hip_cmd.eps',/inches,xsize=xsize,ysize=ysize, $
        xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
      xyouts,0,0,'!3'
      djs_plot, hip.bvcolor,absvmag,psym=3
      device,/close
  endif

; write minimal text file (if asked to do so)
  if keyword_set(textfile) then begin
      openw, unit,textfile,/get_lun
      printf, unit, $
        '# HIP RA(deg) Dec(deg) plx(mas) pmra(mas/yr) pmdec(mas/yr) V(mag) B-V(mag) flags'
      for line=0L,nsample-1L do printf, unit,hip[line],format='(I6,2F16.10,3F10.3,2F8.3,I6)'
      close, unit
      free_lun, unit 
  endif

; make binary data output file (if asked to do so)
  if keyword_set(binaryfile) then begin
      openw, unit,binaryfile,/get_lun
      writeu, unit, hip.vl
      writeu, unit, hip.vb
      writeu, unit, hip.vlvbc
      writeu, unit, hip.nsm
      writeu, unit, hip.radius
      writeu, unit, hip.l
      writeu, unit, hip.b
      close, unit
      free_lun, unit
  endif

  if keyword_set(savefile) then begin
      savefilename= rootdir+'/hip-'+snstr+'-'+diststr+'-'+singlestr+'.sav'
      if keyword_set(giants) then savefilename= rootdir+'/hip-dehnen-giants.sav'
      if keyword_set(dehnen) then savefilename= rootdir+'/hip-dehnen.sav'
      splog, 'saving file '+savefilename
      save, hip,filename=savefilename
  endif

  return, hip
end
