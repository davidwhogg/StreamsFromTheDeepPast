;+
;   NAME:
;      read_gcs
;
;   PURPOSE:
;      read in the Geneva-Copenhagen sample of radial velocities
;
;   USAGE:
;      gcs= read_gcs(/savefile)
;
;   OPTIONAL INPUT:
;      rootdir   directory that holds the catalog
;
;   KEYWORDS:
;      savefile     - if set, write a savefile in the rootdir
;      nobinaries   - exclude objects suspected to be binaries
;      nogiants     - exclude objects suspected to be giants
;      hipnewredux  - use the new reduction of the Hipparcos data for
;                     the kinematical data
;      survey       - only use survey stars, i.e., stars from Coravel
;
;   REVISION HISTORY:
;      2008-07-11 written - Bovy (NYU), based on read_hipparcos
;      2009-01-14 Added nobinaries, nogiants, and hipnewredux
;                 keywords - Bovy
;-
FUNCTION READ_GCS, savefile=savefile, rootdir=rootdir, $
                   nobinaries=nobinaries, nogiants=nogiants, $
                   hipnewredux=hipnewredux, survey=survey

IF ~keyword_set(rootdir) THEN rootdir='/global/data/rv/gcs'

;;Read data
nfield= 51L
; most data types are double
fieldtypes= lonarr(nfield)+5
; some are long
fieldtypes[[0,7,8,15,30,31,35,36,37,38,41,42,43]]= 3
; some are string
fieldtypes[[1,2,3,4,5,6,18,19,26,33,34,50]]= 7
; give dumb names
fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
template= {version:1.0, datastart: 0, $
           delimiter: '|', missingvalue: -1000000000.0, commentsymbol: '#', $
           fieldcount: nfield, fieldtypes: fieldtypes, $
           fieldnames: fieldnames, $
           fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
gcsdata= read_ascii(rootdir+'/table1.dat',num_records=nline,template=template)
splog, 'found '+strtrim(string(n_elements(gcsdata.h01)),2)+$
  ' stars in the GCS catalog'

;;First sort the RELEVANT records by Hipparcos number
sortindx= UNIQ(gcsdata.h00,sort(gcsdata.h00))
nuniq= n_elements(sortindx)
gcsdata2= {h00: lonarr(nuniq), h03:strarr(nuniq), h19: strarr(nuniq), h27:dblarr(nuniq),h28:dblarr(nuniq),h33:strarr(nuniq),h34:strarr(nuniq)}
gcsdata2.h00= gcsdata.h00[sortindx]
gcsdata2.h03= gcsdata.h03[sortindx]
gcsdata2.h19= gcsdata.h19[sortindx]
gcsdata2.h27= gcsdata.h27[sortindx]
gcsdata2.h28= gcsdata.h28[sortindx]
gcsdata2.h33= gcsdata.h33[sortindx]
gcsdata2.h34= gcsdata.h34[sortindx]
;;Then cross-correlate with Hipparcos data
hipexists= where(gcsdata2.h00 NE -1000000000)
ngcship= n_elements(hipexists)
;;Then exclude one without radial velocity or without error
hasvr= where((gcsdata2.h27[hipexists] NE -1000000000.0) AND (abs(gcsdata2.h28[hipexists]) LT 100000.))
hipexists= hipexists[hasvr]
ngcship= n_elements(hipexists)
splog, 'found '+strtrim(string(ngcship),2)+' gcs stars that '+$
  'have an hipparcos entry'
;;Find binaries
IF keyword_set(nobinaries) THEN BEGIN
    nobinary= where((gcsdata2.h03[hipexists] NE '*') AND (gcsdata2.h33[hipexists] NE '*'))
    hipexists= hipexists[nobinary]
    ngcship= n_elements(hipexists)
    splog, 'found '+strtrim(string(ngcship),2)+' gcs stars that '+$
      "have an hipparcos entry and aren't suspected to be binaries"
ENDIF
IF keyword_set(nogiants) THEN BEGIN
;;Find giants
    nogiant= where(gcsdata2.h19[hipexists] NE '*')
    hipexists= hipexists[nogiant]
    ngcship= n_elements(hipexists)
    mesg= 'found '+strtrim(string(ngcship),2)+' gcs stars that '+$
      "have an hipparcos entry"
    IF keyword_set(nobinaries) THEN mesg+= ", aren't suspected to be binaries"
    mesg+=", and aren't suspected to be Giants"
    splog, mesg 
ENDIF
IF keyword_set(survey) THEN BEGIN
;;Find survey stars
    surveystars= where(gcsdata2.h34[hipexists] EQ 'C')
    hipexists= hipexists[surveystars]
    ngcship= n_elements(hipexists)
    mesg= 'found '+strtrim(string(ngcship),2)+' gcs stars that '+$
      "have an hipparcos entry"
    IF keyword_set(nobinaries) THEN mesg+= ", aren't suspected to be binaries"
    IF keyword_set(nogiants) THEN mesg+=", aren't suspected to be Giants"
    mesg+= ", and are survey stars"
    splog, mesg 
ENDIF

;;Read hipparcos
IF keyword_set(hipnewredux) THEN BEGIN
    ;read data
    nfield= 41L
    ;most data types are double
    fieldtypes= lonarr(nfield)+5
    ;some are long
    fieldtypes[[0,1,2,3,14,16,18,22]]= 3
    ;none are string
    ;fieldtypes[[0,2,3,4,10,36,39,42,43,48,52,53,54,55,56,59,60,61,62, $
    ;68,69,70,72,73,74,76,77]]= 7
    ;give dumb names
    fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
    template= {version:1.0, datastart: 0, $
               delimiter: ' ', missingvalue: 0.0, commentsymbol: '#', $
               fieldcount: nfield, fieldtypes: fieldtypes, $
               fieldnames: fieldnames, $
               fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
    hipdata= read_ascii('/global/data/hipparcos/hipnewredux/hip2.dat',$
                        num_records=nline,template=template)
    index= lonarr(ngcship)
    FOR jj= 0L, ngcship-1 DO BEGIN
        index[jj]= (where(hipdata.h00 EQ gcsdata2.h00[hipexists[jj]]))[0]
        IF index[jj] EQ -1 THEN print, gcsdata2.h00[hipexists[jj]]
    ENDFOR
    
    ;make data structure
    str= {gcshipparcos, $
          hip     : 0L,  $      ; HIP number
          ra      : 0D,  $      ; right ascension (1991.25) (deg)
          dec     : 0D,  $      ; declination (1991.25) (deg)
          plx     : 0D,  $      ; parallax (mas)
          pmra    : 0D,  $    ; proper motion in RA direction (mas/yr)
          pmdec   : 0D,  $   ; proper motion in Dec direction (mas/yr)
          e_ra      : 0D,  $    ; standard error in RA (mas)
          e_dec     : 0D,  $    ; standard error in DEC (mas)
          e_plx     : 0D,  $    ; standard error in plx (mas)
          e_pmra    : 0D,  $ ; standard error in proper motion RA (mas/yr)
          e_pmdec   : 0D,  $ ; standard error in proper motion DEC (mas/yr)
          c_decra     : 0D,  $  ; correlation between dec and ra
          c_plxra     : 0D,  $  ; correlation between plx and ra
          c_plxdec    : 0D,  $  ; correlation between plx and dec
          c_pmrara    : 0D,  $  ; correlation between pmra and ra
          c_pmradec   : 0D,  $  ; correlation between pmra and dec
          c_pmraplx   : 0D,  $  ; correlation between pmra and plx
          c_pmdecra   : 0D,  $  ; correlation between pmdec and ra
          c_pmdecdec  : 0D,  $  ; correlation between pmdec and dec
          c_pmdecplx  : 0D,  $  ; correlation between pmdec and plx
          c_pmdecpmra : 0D,  $  ; correlation between pmdec and pmra
          vmag    : 0D,  $      ; apparent V magnitude (mag)
          bvcolor : 0D,  $      ; B-V color (mag)
          flags   : 0L,  $ ; flag integer; 1 for LMS, 2 for UMS, 4 for RG
          l       : 0D,  $      ; Galactic longitude (deg)
          b       : 0D,  $      ; Galactic latitude (deg)
          radius  : 0D,  $ ; heliocentric radius (km) inferred from plx
          vra     : 0D,  $      ; velocity in RA direction (km/s)
          vdec    : 0D,  $      ; velocity in Dec direction (km/s)
          vrvdc   : dblarr(2,2), $ ; measurement covariance on vra,vdec ([km/s]^2)
          vl      : 0D,  $      ; velocity in l direction (km/s)
          vb      : 0D,  $      ; velocity in b direction (km/s)
          vlvbc   : dblarr(2,2), $ ; measurement covariance on vl,vb ([km/s]^2)
          nsm     : dblarr(3,2), $ ; non-square matrix for projection
          sm      : dblarr(3,3),  $ ; square matrix
          vr      : 0D,  $      ;radial velocity (km/s)
          vrvlvbc : dblarr(3,3),  $ ; measurement covariance ([km/s]^2)
          vrc     : 0D   $ ;measurement uncertainty on the radial velocity
         }
    gcs = replicate(str,ngcship)
    
    ;load data structure
    gcs.hip=     hipdata.h00[index]
    gcs.ra=      hipdata.h04[index]/!DPI*180.
    gcs.dec=     hipdata.h05[index]/!DPI*180.
    gcs.plx=     hipdata.h06[index]
    gcs.pmra=    hipdata.h07[index]
    gcs.pmdec=   hipdata.h08[index]
    gcs.bvcolor= hipdata.h23[index]
   
    ; Find the errors
    FOR ii=0L, ngcship-1 DO BEGIN
        U= dblarr(5,5)
        U[0,0]= hipdata.h26[index[ii]] ;1
        U[1,0]= hipdata.h27[index[ii]] ;2
        U[1,1]= hipdata.h28[index[ii]] ;3
        U[2,0]= hipdata.h29[index[ii]] ;4
        U[2,1]= hipdata.h30[index[ii]] ;5
        U[2,2]= hipdata.h31[index[ii]] ;6
        U[3,0]= hipdata.h32[index[ii]] ;7
        U[3,1]= hipdata.h33[index[ii]] ;8
        U[3,2]= hipdata.h34[index[ii]] ;9
        U[3,3]= hipdata.h35[index[ii]] ;10
        U[4,0]= hipdata.h36[index[ii]] ;11
        U[4,1]= hipdata.h37[index[ii]] ;12
        U[4,2]= hipdata.h38[index[ii]] ;13
        U[4,3]= hipdata.h39[index[ii]] ;14
        U[4,4]= hipdata.h40[index[ii]] ;15
        ;;C^{-1}= U\T U
        C= transpose(U)##U
        C= invert(C,/double);;This is the covariance matrix
        ;;Now put everything in place
        gcs[ii].e_ra        = sqrt(C[0,0])
        gcs[ii].e_dec       = sqrt(C[1,1])
        gcs[ii].e_plx       = sqrt(C[2,2])
        gcs[ii].e_pmra      = sqrt(C[3,3])
        gcs[ii].e_pmdec     = sqrt(C[4,4])
        gcs[ii].c_decra     = C[0,1]/(gcs[ii].e_ra*gcs[ii].e_dec)
        gcs[ii].c_plxra     = C[0,2]/(gcs[ii].e_plx*gcs[ii].e_ra)
        gcs[ii].c_plxdec    = C[2,1]/(gcs[ii].e_plx*gcs[ii].e_dec) 
        gcs[ii].c_pmrara    = C[3,0]/(gcs[ii].e_pmra*gcs[ii].e_ra) 
        gcs[ii].c_pmradec   = C[3,1]/(gcs[ii].e_pmra*gcs[ii].e_dec)
        gcs[ii].c_pmraplx   = C[3,2]/(gcs[ii].e_pmra*gcs[ii].e_plx)
        gcs[ii].c_pmdecra   = C[4,0]/(gcs[ii].e_pmdec*gcs[ii].e_ra)
        gcs[ii].c_pmdecdec  = C[4,1]/(gcs[ii].e_pmdec*gcs[ii].e_dec)
        gcs[ii].c_pmdecplx  = C[4,2]/(gcs[ii].e_pmdec*gcs[ii].e_plx)
        gcs[ii].c_pmdecpmra = C[4,3]/(gcs[ii].e_pmdec*gcs[ii].e_pmra)
    ENDFOR

    ;get Galactic quantities
    epoch= 1991.25
    radec_to_lb, gcs.ra,gcs.dec,gcs.pmra,gcs.pmdec,ll,bb,pmll,pmbb,epoch
    gcs.l= ll
    gcs.b= bb
    km_for_unit_plx= 3.08568D16
    radius= km_for_unit_plx/gcs.plx ; km
    gcs.radius= radius          ; km
    invsec_for_unit_pm= !DPI*1.0d0/1.80D2/3.6D6/3.15569D7
    gcs.vra = gcs.pmra *radius*invsec_for_unit_pm ; km/s
    gcs.vdec= gcs.pmdec*radius*invsec_for_unit_pm ; km/s
    gcs.vl  = pmll     *radius*invsec_for_unit_pm ; km/s
    gcs.vb  = pmbb     *radius*invsec_for_unit_pm ; km/s

    ;construct the velocity covariance matrices in ra,dec space
    pmra = gcs.pmra *invsec_for_unit_pm ; 1/s
    pmdec= gcs.pmdec*invsec_for_unit_pm ; 1/s
    radius_err= -radius*gcs.e_plx/gcs.plx ; km
    pmra_err = gcs.e_pmra*invsec_for_unit_pm ; 1/s
    pmdec_err= gcs.e_pmdec*invsec_for_unit_pm ; 1/s
    gcs.vrvdc[0,0]= pmra^2 *radius_err^2 $
      +radius^2*pmra_err^2 $
      +2d0*radius*pmra *radius_err*pmra_err *gcs.c_pmraplx
    gcs.vrvdc[1,1]= pmdec^2*radius_err^2 $
      +radius^2*pmdec_err^2 $
      +2d0*radius*pmdec*radius_err*pmdec_err*gcs.c_pmdecplx
    gcs.vrvdc[0,1]= pmra*pmdec*radius_err^2 $
      +radius*pmra *radius_err*pmdec_err*gcs.c_pmdecplx $
      +radius*pmdec*radius_err*pmra_err *gcs.c_pmraplx $
      +radius^2*pmra_err*pmdec_err*gcs.c_pmdecpmra
    gcs.vrvdc[1,0]= gcs.vrvdc[0,1]
ENDIF ELSE BEGIN
    ;read data
    nfield= 78L
    ;most data types are double
    fieldtypes= lonarr(nfield)+5
    ;some are long
    fieldtypes[[1,6,29,31,47,57,58,63,71]]= 3
    ;some are string
    fieldtypes[[0,2,3,4,10,36,39,42,43,48,52,53,54,55,56,59,60,61,62, $
            68,69,70,72,73,74,76,77]]= 7
    ;give dumb names
    fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
    template= {version:1.0, datastart: 0, $
               delimiter: '|', missingvalue: 0.0, commentsymbol: '#', $
               fieldcount: nfield, fieldtypes: fieldtypes, $
               fieldnames: fieldnames, $
               fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
    hipdata= read_ascii('/global/data/hipparcos/hip_main.dat', $
                        num_records=nline,template=template)
    index= lonarr(ngcship)
    FOR jj= 0L, ngcship-1 DO BEGIN
        index[jj]= (where(hipdata.h01 EQ gcsdata2.h00[hipexists[jj]]))[0]
        IF index[jj] EQ -1 THEN print, gcsdata2.h00[hipexists[jj]]
    ENDFOR
    
    ;make data structure
    str= {gcshipparcos, $
          hip     : 0L,  $      ; HIP number
          ra      : 0D,  $      ; right ascension (1991.25) (deg)
          dec     : 0D,  $      ; declination (1991.25) (deg)
          plx     : 0D,  $      ; parallax (mas)
          pmra    : 0D,  $    ; proper motion in RA direction (mas/yr)
          pmdec   : 0D,  $   ; proper motion in Dec direction (mas/yr)
          vmag    : 0D,  $      ; apparent V magnitude (mag)
          bvcolor : 0D,  $      ; B-V color (mag)
          flags   : 0L,  $ ; flag integer; 1 for LMS, 2 for UMS, 4 for RG
          l       : 0D,  $      ; Galactic longitude (deg)
          b       : 0D,  $      ; Galactic latitude (deg)
          radius  : 0D,  $ ; heliocentric radius (km) inferred from plx
          vra     : 0D,  $      ; velocity in RA direction (km/s)
          vdec    : 0D,  $      ; velocity in Dec direction (km/s)
          vrvdc   : dblarr(2,2), $ ; measurement covariance on vra,vdec ([km/s]^2)
          vl      : 0D,  $      ; velocity in l direction (km/s)
          vb      : 0D,  $      ; velocity in b direction (km/s)
          vlvbc   : dblarr(2,2), $ ; measurement covariance on vl,vb ([km/s]^2)
          nsm     : dblarr(3,2), $ ; non-square matrix for projection
          sm      : dblarr(3,3),  $ ; square matrix
          vr      : 0D,  $      ;radial velocity (km/s)
          vrvlvbc : dblarr(3,3),  $ ; measurement covariance ([km/s]^2)
          vrc     : 0D   $ ;measurement covariance on the radial velocity ([km s^-1]^2)
         }
    gcs = replicate(str,ngcship)
    
    ;load data structure
    gcs.hip=     hipdata.h01[index]
    gcs.ra=      hipdata.h08[index]
    gcs.dec=     hipdata.h09[index]
    gcs.plx=     hipdata.h11[index]
    gcs.pmra=    hipdata.h12[index]
    gcs.pmdec=   hipdata.h13[index]
    gcs.vmag=    hipdata.h05[index]
    gcs.bvcolor= hipdata.h37[index]
    
    ;get Galactic quantities
    epoch= 1991.25
    radec_to_lb, gcs.ra,gcs.dec,gcs.pmra,gcs.pmdec,ll,bb,pmll,pmbb,epoch
    gcs.l= ll
    gcs.b= bb
    km_for_unit_plx= 3.08568D16
    radius= km_for_unit_plx/gcs.plx ; km
    gcs.radius= radius          ; km
    invsec_for_unit_pm= !DPI*1.0d0/1.80D2/3.6D6/3.15569D7
    gcs.vra = gcs.pmra *radius*invsec_for_unit_pm ; km/s
    gcs.vdec= gcs.pmdec*radius*invsec_for_unit_pm ; km/s
    gcs.vl  = pmll     *radius*invsec_for_unit_pm ; km/s
    gcs.vb  = pmbb     *radius*invsec_for_unit_pm ; km/s
    
    ;construct the velocity covariance matrices in ra,dec space
    pmra = gcs.pmra *invsec_for_unit_pm ; 1/s
    pmdec= gcs.pmdec*invsec_for_unit_pm ; 1/s
    radius_err= radius*hipdata.h16[index]/gcs.plx ; km
    pmra_err = hipdata.h17[index]*invsec_for_unit_pm ; 1/s
    pmdec_err= hipdata.h18[index]*invsec_for_unit_pm ; 1/s
    gcs.vrvdc[0,0]= pmra^2 *radius_err^2 $
      +radius^2*pmra_err^2 $
      +2d0*radius*pmra *radius_err*pmra_err *hipdata.h24[index]
    gcs.vrvdc[1,1]= pmdec^2*radius_err^2 $
      +radius^2*pmdec_err^2 $
      +2d0*radius*pmdec*radius_err*pmdec_err*hipdata.h27[index]
    gcs.vrvdc[0,1]= pmra*pmdec*radius_err^2 $
      +radius*pmra *radius_err*pmdec_err*hipdata.h27[index] $
      +radius*pmdec*radius_err*pmra_err *hipdata.h24[index] $
      +radius^2*pmra_err*pmdec_err*hipdata.h28[index]
    gcs.vrvdc[1,0]= gcs.vrvdc[0,1]
ENDELSE

; rotate the covariance matrices into the l,b space
zero= dblarr(ngcship)
unity= zero+1
radec_to_lb, gcs.ra,gcs.dec,unity,zero,ll,bb,cosine,sine,epoch
rotmax= dblarr(2,2,ngcship)
rotmax[0,0,*]= cosine
rotmax[0,1,*]= sine
rotmax[1,0,*]= -sine
rotmax[1,1,*]= cosine
for i=0L,ngcship-1L do begin
    u= rotmax[*,*,i]
    tmp= u##gcs[i].vrvdc##transpose(u)
; re-symmetrize covariance matrix explicitly (to remove numerical noise)
    gcs[i].vlvbc= 0.5d0*(tmp+transpose(tmp))
endfor

; make and load non-square matrices for projection
ll= gcs.l*!DPI/1.80D2
bb= gcs.b*!DPI/1.80D2
tmpmat= dblarr(3,3)
for i=0L,ngcship-1L do begin
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
    gcs[i].sm= tmpmat
    gcs[i].nsm= tmpmat[*,1:2]
endfor

; load radial velocity and error
gcs.vr= gcsdata2.h27[hipexists]
FOR jj=0L, ngcship-1 DO BEGIN
    gcs[jj].vrvlvbc[1:2,1:2]= gcs[jj].vlvbc
    gcs[jj].vrvlvbc[0,0]= (gcsdata2.h28[hipexists[jj]])^2
    gcs[jj].vrc=gcs[jj].vrvlvbc[0,0]
ENDFOR

IF keyword_set(savefile) THEN save, filename=rootdir+'/gcs.sav', gcs

RETURN, gcs

END
