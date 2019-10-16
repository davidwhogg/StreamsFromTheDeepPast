;+
;   NAME:
;      select_sample
;   PURPOSE:
;      select the sample of stars as in Aumer & Binney
;   CALLING SEQUENCE:
;      select_sample
;   INPUT:
;      tycho_savefile   - the raw Tycho data
;      filename         - filename for savefile
;   KEYWORDS:
;      savefile  - write a savefile with the resulting structure to
;                  the current dir
;      cmd       - plot the color magnitude diagram of the selected
;                  stars
;      formal    - use formal Hipparcos errors
;      giants    - also include non-MS stars
;      QA        - plot vlvb correlation histogram
;   OUTPUT:
;      structure containing the sample
;   REVISION HISTORY:
;      2009-05-25 - Written - Bovy
;-
FUNCTION SELECT_SAMPLE, tycho_savefile=tycho_savefile, savefile=savefile, $
                        cmd=cmd, formal=formal, filename=filename, QA=QA, $
                        giants=giants

;;First read the Tycho data and save it
IF file_test(tycho_savefile) THEN BEGIN
    restore, filename=tycho_savefile
ENDIF ELSE BEGIN
    tycdata=read_tycho(/read_raw)
    save, filename=tycho_savefile, tycdata
ENDELSE


;;Now determine in each of 16x16x10 bins in sin(b), l, and 
;;-0.3 < (B-V)_T < 1.5 which Hipparcos stars are brighter than the
;;second brightest star that is in Tycho but not in Hipparcos

;;Actually, just determine the VT mag of that star
;;First rotate the stars into l and b
sample_tempfile= '../lsr2/sample_tmp.sav'
nbins_sinbb= 16
nbins_ll= 16
nbins_BV= 10
if file_test(sample_tempfile) THEN BEGIN
    restore, filename=sample_tempfile
ENDIF ELSE BEGIN
    ntyc= n_elements(tycdata.h00)
    radec_to_lb, tycdata.h27, tycdata.h28,dblarr(ntyc),dblarr(ntyc),$
      ll,bb,pmll,pmbb, 1991.7;;This is not quite right
    sinbb= sin(bb/180.*!DPI)
    hipnumbers= lonarr(120000L)
    nhipnumbers= 0L
    VTmags= dblarr(nbins_sinbb,nbins_ll,nbins_BV)
;    FOR ii=0L, 15 DO FOR jj=0L, 0 DO FOR kk=0L,0 DO BEGIN
    FOR ii=0L, nbins_sinbb-1 DO FOR jj=0L, nbins_ll-1 DO FOR kk=0L,nbins_BV-1 DO BEGIN
        thesestars= where(sinbb GE (-1+2./nbins_sinbb*ii) AND $
                          sinbb LT (-1+2./nbins_sinbb*(ii+1)) AND $
                          ll GE 360./nbins_ll*jj AND $
                          ll LT 360./nbins_ll*(jj+1) AND $
                          (1.*(tycdata.h19-tycdata.h21)) GE (1.8/nbins_BV*kk-0.3) AND $
                          (1.*(tycdata.h19-tycdata.h21)) LT (1.8/nbins_BV*(kk+1)-0.3))
        ;;Find the second brightest star that is not in Tycho
        ;;In V   = VT -0.090*(BT-VT)
        mags= tycdata[thesestars].h21;-0.09*(tycdata[thesestars].h19 -$
                                      ;      tycdata[thesestars].h21)
        sortindx= sort(mags)
        sortedmags= mags[sortindx]
        sortedhipns= tycdata[thesestars[sortindx]].h25
        notinhipindx= where(sortedhipns EQ 0)
        VTmags[ii,jj,kk]= sortedmags[notinhipindx[1]]
        ;print, ii*nbins_ll*nbins_BV+jj*nbins_BV+kk, VTmags[ii,jj,kk], tycdata[thesestars[sortindx[notinhipindx[0]]]].h21, tycdata[thesestars[sortindx[notinhipindx[1]]]].h21
        continue
        ;stop
        notinhip=0L
        indx_thesestars= 0L
        WHILE notinhip LT 2 DO BEGIN
            IF tycdata[thesestars[sortindx[indx_thesestars]]].h25 EQ 0 AND $
              tycdata[thesestars[sortindx[indx_thesestars]]].h33 NE 'D' AND $
              tycdata[thesestars[sortindx[indx_thesestars]]].h33 NE 'P' THEN BEGIN
                ;print, notinhip, mags[sortindx[indx_thesestars]], 7.3+abs(sinbb[thesestars[sortindx[indx_thesestars]]])*1.1, sinbb[thesestars[sortindx[indx_thesestars]]]
                notinhip+= 1
                indx_thesestars+= 1
                CONTINUE
            ENDIF
            hipnumbers[nhipnumbers]= $
              tycdata[thesestars[sortindx[indx_thesestars]]].h25
            indx_thesestars+= 1
            nhipnumbers+= 1
        ENDWHILE
        print, nhipnumbers, ii*nbins_ll*nbins_BV+jj*nbins_BV+kk
    ENDFOR
    save, filename=sample_tempfile, VTmags
;    save, filename=sample_tempfile, hipnumbers,nhipnumbers
ENDELSE

sample_tempfile= '../lsr2/sample_tmp2.sav'
if file_test(sample_tempfile) THEN BEGIN
    restore, filename=sample_tempfile
ENDIF ELSE BEGIN
    ;;Read the old reduction of the Hipparcos data for the VT magnitudes
    rootdir='/global/data/hipparcos'
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
    oldhipdata= read_ascii(rootdir+'/hip_main.dat',num_records=nline,template=template)
;;And find those hipnumbers of the stars brighter than the VTmag in
;;VTmags
    noldhip= n_elements(oldhipdata.h00)
    radec_to_lb, oldhipdata.h08, oldhipdata.h09,dblarr(noldhip),dblarr(noldhip),$
      ll,bb,pmll,pmbb, 1991.25
    sinbb= sin(bb/180.*!DPI)
    hipnumbers= lonarr(noldhip)
    nhipnumbers= 0L
    FOR ii=0L, nbins_sinbb-1 DO FOR jj=0L, nbins_ll-1 DO FOR kk=0L,nbins_BV-1 DO BEGIN
        thesestars= where(sinbb GE (-1+2./nbins_sinbb*ii) AND $
                          sinbb LT (-1+2./nbins_sinbb*(ii+1)) AND $
                          ll GE 360./nbins_ll*jj AND $
                          ll LT 360./nbins_ll*(jj+1) AND $
;                          (1.*(oldhipdata.h37)/0.85) GE (1.8/nbins_BV*kk-0.3) AND $
;                          (1.*(oldhipdata.h37)/0.85) LT (1.8/nbins_BV*(kk+1)-0.3))
                          (1.*(oldhipdata.h32-oldhipdata.h34)) GE (1.8/nbins_BV*kk-0.3) AND $
                          (1.*(oldhipdata.h32-oldhipdata.h34)) LT (1.8/nbins_BV*(kk+1)-0.3))
        IF thesestars[0] EQ -1 THEN CONTINUE
        indx= where(oldhipdata.h34[thesestars] LE VTmags[ii,jj,kk]);; AND oldhipdata.h59[thesestars] EQ '')
        IF indx[0] EQ -1 THEN CONTINUE
        addn= n_elements(indx)
        hipnumbers[nhipnumbers:nhipnumbers+addn-1]= oldhipdata.h01[thesestars[indx]]
        nhipnumbers+= addn
        print, nhipnumbers
    ENDFOR
    save, filename=sample_tempfile, hipnumbers, nhipnumbers
ENDELSE

print, nhipnumbers


;;Now read the new reduction of the Hipparcos data
rootdir='/global/data/hipparcos/hipnewredux'
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
hipdata= read_ascii(rootdir+'/hip2.dat',num_records=nline,template=template)
nhip= n_elements(hipdata.h01)
splog, 'found ',n_elements(hipdata.h01),' stars in the Hipparcos catalog'

finalhipindx= lonarr(nhip)
nfinalhipindx= 0L
FOR ii=0L, nhipnumbers-1 DO BEGIN
    ;;Check for each one that it is single, has a relative plx error <
    ;;10% and that it lies on the main sequence
    indx_inhip= (where(hipdata.h00 EQ hipnumbers[ii]))[0]
    IF indx_inhip EQ -1 THEN CONTINUE
    IF hipdata.h01[indx_inhip] LT 10 AND hipdata.h06[indx_inhip]/$
;      sqrt(C[2,2]) GE 10 THEN BEGIN
      hipdata.h11[indx_inhip] GE 10 THEN BEGIN
         IF keyword_set(giants) THEN BEGIN
             finalhipindx[nfinalhipindx]= indx_inhip
             nfinalhipindx+= 1
         ENDIF ELSE BEGIN
             ;;Check MS here
             BV= hipdata.h23[indx_inhip]
                                ;Abs Mag
             Mhip= hipdata.h19[indx_inhip]-5.*alog10(100./hipdata.h06[indx_inhip])
             IF BV LE .35 THEN BEGIN
                 IF Mhip GT (7.5*BV-3.75) AND Mhip LT (4.62*BV+2.383) THEN BEGIN
                     finalhipindx[nfinalhipindx]= indx_inhip
                     nfinalhipindx+= 1
                 ENDIF
             ENDIF ELSE IF BV GT .35 AND BV LE 0.5 THEN BEGIN
                 IF Mhip GT (7.5*BV-3.75) AND Mhip LT (8.33*BV+1.0845)  THEN BEGIN
                     finalhipindx[nfinalhipindx]= indx_inhip
                     nfinalhipindx+= 1
                 ENDIF
             ENDIF ELSE IF BV GT 0.5 AND BV LE .65 THEN BEGIN
                 IF Mhip GT (15.33*BV-7.665) AND Mhip LT (8.33*BV+1.0845) THEN BEGIN
                     finalhipindx[nfinalhipindx]= indx_inhip
                     nfinalhipindx+= 1
                 ENDIF
             ENDIF ELSE IF BV GT .65 AND BV LE .8 THEN BEGIN
                 IF Mhip GT (15.33*BV-7.665) AND Mhip LT (3.33*BV+4.3375) THEN BEGIN
                     finalhipindx[nfinalhipindx]= indx_inhip
                     nfinalhipindx+= 1
                 ENDIF
             ENDIF ELSE IF BV GT .8 AND BV LE 1.25 THEN BEGIN
                 IF Mhip GT (4.43*BV+1.055) AND Mhip LT (3.33*BV+4.3375) THEN BEGIN
                     finalhipindx[nfinalhipindx]= indx_inhip
                     nfinalhipindx+= 1
                 ENDIF
             ENDIF ELSE IF BV GT 1.25 THEN BEGIN
                 IF Mhip GT (4.43*BV+1.055) AND Mhip LT (6.5*BV+0.375) THEN BEGIN
                     finalhipindx[nfinalhipindx]= indx_inhip
                     nfinalhipindx+= 1
                 ENDIF
             ENDIF
         ENDELSE
    ENDIF            
ENDFOR
print, nfinalhipindx
index= finalhipindx[0:nfinalhipindx-1]
nsample= n_elements(index)


;;Plot the cmd
IF keyword_set(cmd) THEN BEGIN
    k_print, filename='cmd.ps', xsize=3.375, ysize=4  ; actual width of ApJ column!
    plot_these= where(hipdata.h23[index] NE 0.)
    djs_plot, hipdata.h23[index[plot_these]], hipdata.h19[index[plot_these]]-$
      5.*alog10(100./hipdata.h06[index[plot_these]]), psym=3, yrange=[10,-6], $
      xtitle='B-V [mag]', ytitle='M_{Hip} [mag]',xrange=[-0.3,1.6], $
      xcharsize=1,ycharsize=1, position=[0.1,0.1,0.9,0.9]
    k_end_print
ENDIF

;;Now make a structure containing all of the relevant information for
;;this sample
; make data structure
str= {aumer, $
      hip     : 0L,  $          ; HIP number
      ra      : 0D,  $          ; right ascension (1991.25) (deg)
      dec     : 0D,  $          ; declination (1991.25) (deg)
      plx     : 0D,  $          ; parallax (mas)
      pmra    : 0D,  $        ; proper motion in RA direction (mas/yr)
      pmdec   : 0D,  $       ; proper motion in Dec direction (mas/yr)
      e_ra      : 0D,  $        ; standard error in RA (mas)
      e_dec     : 0D,  $        ; standard error in DEC (mas)
      e_plx     : 0D,  $        ; standard error in plx (mas)
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
      hmag    : 0D,  $          ; apparent Hip magnitude (mag)
      vmag    : 0D, $           ; apparent V magnitude (mag)
      bvcolor : 0D,  $          ; B-V color (mag)
      vt      : 0D,  $          ; VT mag (mag)
      bt      : 0D,  $          ; BT mag (mag)
      e_vt      : 0D,  $        ; error in VT mag (mag)
      e_bt      : 0D,  $        ; error in BT mag (mag)
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
hip = replicate(str,nsample)

; load hipdata structure
hip.hip=     hipdata.h00[index]
hip.ra=      hipdata.h04[index]/!DPI*180.
hip.dec=     hipdata.h05[index]/!DPI*180.
hip.plx=     hipdata.h06[index]
hip.pmra=    hipdata.h07[index]
hip.pmdec=   hipdata.h08[index]
hip.hmag=    hipdata.h19[index]
hip.bvcolor= hipdata.h23[index]

; Find the errors
FOR ii=0L, nsample-1 DO BEGIN
    U= dblarr(5,5)
    U[0,0]= hipdata.h26[index[ii]];1
    U[1,0]= hipdata.h27[index[ii]];2
    U[1,1]= hipdata.h28[index[ii]];3
    U[2,0]= hipdata.h29[index[ii]];4
    U[2,1]= hipdata.h30[index[ii]];5
    U[2,2]= hipdata.h31[index[ii]];6
    U[3,0]= hipdata.h32[index[ii]];7
    U[3,1]= hipdata.h33[index[ii]];8
    U[3,2]= hipdata.h34[index[ii]];9
    U[3,3]= hipdata.h35[index[ii]];10
    U[4,0]= hipdata.h36[index[ii]];11
    U[4,1]= hipdata.h37[index[ii]];12
    U[4,2]= hipdata.h38[index[ii]];13
    U[4,3]= hipdata.h39[index[ii]];14
    U[4,4]= hipdata.h40[index[ii]];15
    ;;C^{-1}= U\T U
    C= transpose(U)##U
    C= invert(C,/double);;This is the covariance matrix
    ;;Now put everything in place
    IF keyword_set(formal) THEN BEGIN
        hip[ii].e_ra        = hipdata.h09[index[ii]]
        hip[ii].e_dec       = hipdata.h10[index[ii]]
        hip[ii].e_plx       = hipdata.h11[index[ii]] 
        hip[ii].e_pmra      = hipdata.h12[index[ii]]
        hip[ii].e_pmdec     = hipdata.h13[index[ii]]
    ENDIF ELSE BEGIN
        hip[ii].e_ra        = sqrt(C[0,0])
        hip[ii].e_dec       = sqrt(C[1,1])
        hip[ii].e_plx       = sqrt(C[2,2])
        hip[ii].e_pmra      = sqrt(C[3,3])
        hip[ii].e_pmdec     = sqrt(C[4,4])
    ENDELSE
    hip[ii].c_decra     = C[0,1]/sqrt(C[0,0]*C[1,1])
    hip[ii].c_plxra     = C[0,2]/sqrt(C[0,0]*C[2,2])
    hip[ii].c_plxdec    = C[2,1]/sqrt(C[1,1]*C[2,2])
    hip[ii].c_pmrara    = C[3,0]/sqrt(C[0,0]*C[3,3])
    hip[ii].c_pmradec   = C[3,1]/sqrt(C[1,1]*C[3,3])
    hip[ii].c_pmraplx   = C[3,2]/sqrt(C[2,2]*C[3,3])
    hip[ii].c_pmdecra   = C[4,0]/sqrt(C[4,4]*C[0,0])
    hip[ii].c_pmdecdec  = C[4,1]/sqrt(C[4,4]*C[1,1])
    hip[ii].c_pmdecplx  = C[4,2]/sqrt(C[4,4]*C[2,2])
    hip[ii].c_pmdecpmra = C[4,3]/sqrt(C[4,4]*C[3,3])
ENDFOR


; get Galactic quantities
epoch= 1991.25
radec_to_lb, hip.ra,hip.dec,hip.pmra,hip.pmdec,ll,bb,pmll,pmbb,epoch
hip.l= ll
hip.b= bb
km_for_unit_plx= 3.08568D16
radius= km_for_unit_plx/hip.plx ; km
hip.radius= radius              ; km
invsec_for_unit_pm= !DPI*1.0d0/1.80D2/3.6D6/3.15569D7
hip.vra = hip.pmra *radius*invsec_for_unit_pm ; km/s
hip.vdec= hip.pmdec*radius*invsec_for_unit_pm ; km/s
hip.vl  = pmll     *radius*invsec_for_unit_pm ; km/s
hip.vb  = pmbb     *radius*invsec_for_unit_pm ; km/s

; construct the velocity covariance matrices in ra,dec space
pmra = hip.pmra *invsec_for_unit_pm ; 1/s
pmdec= hip.pmdec*invsec_for_unit_pm ; 1/s
radius_err= -radius*hip.e_plx/hip.plx ; km
pmra_err = hip.e_pmra*invsec_for_unit_pm ; 1/s
pmdec_err= hip.e_pmdec*invsec_for_unit_pm ; 1/s
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

;;Cross-correlate with Tycho to get VT and BT
ntyc= n_elements(tycdata.h25)
FOR indx= 0L, ntyc-1 DO BEGIN
    IF tycdata[indx].h25 EQ 0 THEN CONTINUE
    ii= (where(tycdata[indx].h25 EQ hip.hip))[0]
    IF ii EQ -1 THEN BEGIN
        IF keyword_set(QA)  THEN splog, 'warning: no Hipparcos counterpart found'
        CONTINUE
    ENDIF
    hip[ii].vt= tycdata[indx].h21
    hip[ii].e_vt= tycdata[indx].h22
    hip[ii].bt= tycdata[indx].h19
    hip[ii].e_bt= tycdata[indx].h20
ENDFOR
indx= where(hip.vt EQ 0.)
splog, strtrim(string(n_elements(indx)),2) +' stars without VT'

;;Cross-correlate with the original Hipparcos catalog to get V
;;Read the original Hipparcos catalog
hip_orig=read_hipparcos(minsn=-100000000000000000,maxdist=1000000000);;hack to make sure almost the whole catalog is read
nhip_orig= n_elements(hip_orig.hip)
FOR indx= 0L, nhip_orig-1 DO BEGIN
    ii= (where(hip_orig[indx].hip EQ hip.hip))[0]
    IF ii EQ -1 THEN BEGIN
        IF keyword_set(QA)  THEN splog, 'warning: no Original Hipparcos counterpart found'
        CONTINUE
    ENDIF
    hip[ii].vmag= hip_orig[indx].vmag
ENDFOR
indx= where(hip.vmag EQ 0.)
splog, strtrim(string(n_elements(indx)),2) +' stars without V magnitude'
    

if keyword_set(savefile) then begin
    if ~keyword_set(filename) THEN savefilename= 'hip2-aumer.sav' ELSE savefilename=filename
    splog, 'saving file '+savefilename
    save, hip,filename=savefilename
endif


IF keyword_set(QA) THEN BEGIN
    k_print, filename='QA_rho_pmrapmdec.ps'
    hogg_plothist, hip.c_pmdecpmra, xtitle=textoidl('\rho_{\mu_{\alpha^*}\mu_\delta}'), ytitle='N'
    k_end_print
END

return, hip

END
