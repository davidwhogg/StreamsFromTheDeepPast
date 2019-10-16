;+
;   NAME:
;      read_isochrone
;   PURPOSE:
;      read an isochrone
;   CALLING SEQUENCE:
;      iso=read_isochrone(filename,/BTVT)
;   INPUT:
;      filename - name of the file that holds the isochrone
;   KEYWORDS:
;      BTVT - Read a Tycho isochrone (default: UBVRIJHK)
;   OUTPUT:
;      everything in the isochrone file (structure); REMEMBER: THE
;      FIRST FIELD IS EMPTY
;   REVISION HISTORY:
;      2009-09-14 - Written - Bovy (NYU)
;-
FUNCTION READ_ISOCHRONE, filename, BTVT=BTVT
IF keyword_set(BTVT) THEN BEGIN
    ;;Read the isochrone
    nfield= 16L
    ;;most data types are double
    fieldtypes= lonarr(nfield)+5
    ;;some are long
    fieldtypes[[12]]= 3
    ;;give dumb names
    fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
    template= {version:1.0, datastart: 0, $
               delimiter: string(9B), missingvalue: 0.0, commentsymbol: '#', $
               fieldcount: nfield, fieldtypes: fieldtypes, $
               fieldnames: fieldnames, $
               fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
    iso= read_ascii(filename,num_records=nline,template=template)
ENDIF ELSE BEGIN
    ;;Read the isochrone
    nfield= 21L
    ;;most data types are double
    fieldtypes= lonarr(nfield)+5
    ;;some are long
    fieldtypes[[18]]= 3
    ;;give dumb names
    fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
    template= {version:1.0, datastart: 0, $
               delimiter: string(9B), missingvalue: 0.0, commentsymbol: '#', $
               fieldcount: nfield, fieldtypes: fieldtypes, $
               fieldnames: fieldnames, $
               fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
    iso= read_ascii(filename,num_records=nline,template=template)
ENDELSE
RETURN, iso
END
