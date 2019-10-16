;+
;   NAME:
;      create_phot_gcs_sample
;   PURPOSE:
;      create a photometric sample using the GCS metallicities
;   INPUT:
;      rootdir - root dir for GCS data
;   KEYWORDS:
;   OUTPUT:
;      structure
;   HISTORY:
;      2009-11-20 - Written - Bovy (NYU)
;-
FUNCTION CREATE_PHOT_GCS_SAMPLE, rootdir=rootdir
IF ~keyword_set(rootdir) THEN rootdir='/global/data/rv/gcs'
;;Read data
nfield= 29L
;;most data types are double (?)
fieldtypes= lonarr(nfield)+5
;;some are long
fieldtypes[[0,11,17,18,19]]= 3
;;some are string
fieldtypes[[1,2,3,4,5,6,28]]= 7
;;give dumb names
fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
template= {version:1.0, datastart: 0, $
           delimiter: '|', missingvalue: -1000000000.0, commentsymbol: '#', $
           fieldcount: nfield, fieldtypes: fieldtypes, $
           fieldnames: fieldnames, $
           fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
gcsdata= read_ascii(rootdir+'/GCSIII.dat',num_records=nline,template=template)
splog, 'found '+strtrim(string(n_elements(gcsdata.h01)),2)+$
  ' stars in the GCS III catalog'

;;First sort the RELEVANT records by Hipparcos number
sortindx= UNIQ(gcsdata.h00,sort(gcsdata.h00))
nuniq= n_elements(sortindx)
gcsdata2= {h00: lonarr(nuniq), h10: dblarr(nuniq), h14: dblarr(nuniq), $
           h15:dblarr(nuniq), h16:dblarr(nuniq)}
gcsdata2.h00= gcsdata.h00[sortindx]
gcsdata2.h10= gcsdata.h10[sortindx]
gcsdata2.h14= gcsdata.h14[sortindx]
gcsdata2.h15= gcsdata.h15[sortindx]
gcsdata2.h16= gcsdata.h16[sortindx]


;;Then cross-correlate with Hipparcos data
hipexists= where(gcsdata2.h00 NE -1000000000)
ngcship= n_elements(hipexists)


;;read Hipparcos data
nfield= 41L
;;most data types are double
fieldtypes= lonarr(nfield)+5
;;some are long
fieldtypes[[0,1,2,3,14,16,18,22]]= 3
;;give dumb names
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
    
;;make data structure
str= {gcshipparcos, $
      hip     : 0L,  $          ; HIP number
      ra      : 0D,  $          ; right ascension (1991.25) (deg)
      dec     : 0D,  $          ; declination (1991.25) (deg)
      plx     : 0D,  $          ; parallax (mas)
      pmra    : 0D,  $        ; proper motion in RA direction (mas/yr)
      pmdec   : 0D,  $       ; proper motion in Dec direction (mas/yr)
      e_ra      : 0D,  $        ; standard error in RA (mas)
      e_dec     : 0D,  $        ; standard error in DEC (mas)
      e_plx     : 0D,  $        ; standard error in plx (mas)
      e_plx_formal     : 0D,  $        ; Formal error in plx (mas)
      e_pmra    : 0D,  $ ; standard error in proper motion RA (mas/yr)
      e_pmdec   : 0D,  $ ; standard error in proper motion DEC (mas/yr)
      c_decra     : 0D,  $      ; correlation between dec and ra
      c_plxra     : 0D,  $      ; correlation between plx and ra
      c_plxdec    : 0D,  $      ; correlation between plx and dec
      c_pmrara    : 0D,  $      ; correlation between pmra and ra
      c_pmradec   : 0D,  $      ; correlation between pmra and dec
      c_pmraplx   : 0D,  $      ; correlation between pmra and plx
      c_pmdecra   : 0D,  $      ; correlation between pmdec and ra
      c_pmdecdec  : 0D,  $      ; correlation between pmdec and dec
      c_pmdecplx  : 0D,  $      ; correlation between pmdec and plx
      c_pmdecpmra : 0D,  $      ; correlation between pmdec and pmra
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
      sm      : dblarr(3,3),  $ ; square matrix
      vr      : 0D,  $          ;radial velocity (km/s)
      vrvlvbc : dblarr(3,3),  $ ; measurement covariance ([km/s]^2)
      vrc     : 0D,   $ ;measurement uncertainty on the radial velocity
      vmag    : 0D,  $          ; apparent V magnitude (mag)
      bvcolor : 0D,  $          ; B-V color (mag)
      feh : 0D, $               ; Metallicity (wrt Solar)
      age : 0D, $               ; Age in Gyr
      lage : 0D, $              ; 1-sigma lower on the age limit in Gyr
      uage : 0D  $              ; 1-sigma upper on the age limit in Gyr
}
gcs = replicate(str,ngcship)

;;load data structure
gcs.hip=     hipdata.h00[index]
gcs.ra=      hipdata.h04[index]/!DPI*180.
gcs.dec=     hipdata.h05[index]/!DPI*180.
gcs.pmra=    hipdata.h07[index]
gcs.pmdec=   hipdata.h08[index]
gcs.plx=     hipdata.h06[index]
gcs.bvcolor= hipdata.h23[index]
gcs.e_plx_formal= hipdata.h11[index]

;; Find the errors
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


;;get Galactic quantities
epoch= 1991.25
radec_to_lb, gcs.ra,gcs.dec,gcs.pmra,gcs.pmdec,ll,bb,pmll,pmbb,epoch
gcs.l= ll
gcs.b= bb
km_for_unit_plx= 3.08568D16
radius= km_for_unit_plx/gcs.plx ; km
gcs.radius= radius              ; km
invsec_for_unit_pm= !DPI*1.0d0/1.80D2/3.6D6/3.15569D7
gcs.vra = gcs.pmra *radius*invsec_for_unit_pm ; km/s
gcs.vdec= gcs.pmdec*radius*invsec_for_unit_pm ; km/s
gcs.vl  = pmll     *radius*invsec_for_unit_pm ; km/s
gcs.vb  = pmbb     *radius*invsec_for_unit_pm ; km/s

;;construct the velocity covariance matrices in ra,dec space
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
;;Read data
nfield= 51L
;;most data types are double
fieldtypes= lonarr(nfield)+5
;;some are long
fieldtypes[[0,7,8,15,30,31,35,36,37,38,41,42,43]]= 3
;;some are string
fieldtypes[[1,2,3,4,5,6,18,19,26,33,34,50]]= 7
;;give dumb names
fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
template= {version:1.0, datastart: 0, $
           delimiter: '|', missingvalue: -1000000000.0, commentsymbol: '#', $
           fieldcount: nfield, fieldtypes: fieldtypes, $
           fieldnames: fieldnames, $
           fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
oldgcsdata= read_ascii(rootdir+'/table1.dat',num_records=nline,template=template)
splog, 'found '+strtrim(string(n_elements(gcsdata.h01)),2)+$
  ' stars in the GCS catalog'
;;First sort the RELEVANT records by Hipparcos number
sortindx= UNIQ(oldgcsdata.h00,sort(oldgcsdata.h00))
nuniq= n_elements(sortindx)
FOR ii=0L, n_elements(gcs.hip)-1 DO BEGIN
    old_indx= (where(oldgcsdata.h00[sortindx] EQ gcs[ii].hip))[0]
    IF old_indx EQ -1 THEN BEGIN
        splog, "Warning, no old GCS counterpart found"
        CONTINUE
    ENDIF
    gcs[ii].vr= oldgcsdata.h27[sortindx[old_indx]]
    gcs[ii].vrvlvbc[1:2,1:2]= gcs[ii].vlvbc
    gcs[ii].vrvlvbc[0,0]= (oldgcsdata.h28[sortindx[old_indx]])^2
    gcs[ii].vrc=gcs[ii].vrvlvbc[0,0]
ENDFOR

gcs.feh = gcsdata2.h10[hipexists]
gcs.age = gcsdata2.h14[hipexists]
gcs.lage = gcsdata2.h15[hipexists]
gcs.uage = gcsdata2.h16[hipexists]


;;Get Vmag from old reduction
hip_old= read_hipparcos(minsn=-100000000000000000,maxdist=1000000000);;hack to make sure almost the whole catalog is read
;;These should be in the right order, but let's cross-correlate just
;;to make sure
FOR ii= 0L, n_elements(gcs.hip)-1 DO BEGIN
    old_indx= (where(hip_old.hip EQ gcs[ii].hip))[0]
    IF old_indx EQ -1 THEN BEGIN
        CONTINUE
    ENDIF
    gcs[ii].vmag= hip_old[old_indx].vmag
ENDFOR


;;Some quality cuts
print, n_elements(where(gcs.plx/gcs.e_plx GE 10.))
print, n_elements(where(gcs.plx/gcs.e_plx_formal GE 10.))
good= where(gcs.plx/gcs.e_plx_formal GE 10. AND gcs.vr NE -1000000000.0 AND $
            sqrt(gcs.vrc) LT 100000. AND gcs.feh NE -1000000000.0)
RETURN, gcs[good]
end
