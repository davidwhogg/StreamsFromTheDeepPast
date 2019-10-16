;
; function to load MATLAB data.
;
; The MATLAB data is read from a MATLAB save file into a list of variables passed
; as arguments.
; Currently the maximum number of variables which can be read is 25.
;
; The format is :
;   
;      status = load_MAT_data(file, MAT_names, arg1, arg2...)
;
; The return value is an array of error message strings.
; If there is no error the first element is an empty string.
;   where :
;
;	file	    - string,	    name of file containing MATLAB data.
;	MAT_names   - string array, MATLAB names for variables to return
;	arg1,..argn - list of variables in which the data will be returned.
;		      MATLAB data for MAT_names(0) => arg1
;		      MATLAB data for MAT_names(1) => arg2
;		      MATLAB data for MAT_names(n-1) => argn
;
;		      There must be one arg for every element of MAT_names.
;                     If you want to read more than 25 matrices just add more
;		      args to the parameter list.			
;
;  Author: Nigel Wade, Space Plasma Physics Group, Leicester University.
;  Copyright: Leicester University.
;  (Use as you wish, but please leave this notice.)
;


FUNCTION load_MAT_data, file, MAT_names, $
	    arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, $
	    arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17, $
	    arg18, arg19, arg20, arg21, arg22, arg23, arg24, arg25

IF N_PARAMS() LT 3 OR N_PARAMS()-2 LT N_ELEMENTS(MAT_names) THEN $
    RETURN, ['load_mat_data','Bad args to load_mat_data']
    
ON_IOERROR, bad_read
OPENR, MAT_lun, file, /GET_LUN

; define the structure for the MAT-file header
MAT_header = {MAT_header, type:0L, mrows:0L, ncols:0L, $
		imagf:0L, namlen:0L}

; read each matrix from the MAT-file and see if its name matches
; one which is wanted
WHILE (NOT EOF(MAT_lun)) DO BEGIN
    READU, MAT_lun, MAT_header
    
    MAT_name = STRING(BINDGEN(MAT_header.namlen)+32b)
    READU, MAT_lun, MAT_name
    
    ; data type indicates how the data is stored, double, float, long etc.
    data_type = (MAT_header.type MOD 100)/ 10

    ; text is a flag to indicate text stored as numbers.
    text = (MAT_header.type MOD 2)

    ; if the data is not complex create a structure to read the data with only one real array.
    IF MAT_header.imagf NE 1  THEN BEGIN
	IF data_type EQ 0 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/DOUBLE,/NOZERO)} $
	ELSE IF data_type EQ 1 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/FLOAT,/NOZERO)} $
	ELSE IF data_type EQ 2 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/LONG,/NOZERO)} $
	ELSE IF data_type EQ 3 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/INT,/NOZERO)} $
	ELSE IF data_type EQ 4 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/INT,/NOZERO)} $
	ELSE BEGIN
	    RETURN, ['load_mat_data','Error in MATLAB data type.']
	ENDELSE
    ENDIF ELSE BEGIN
        ; if the data is complex create a structure to read the data with only real and imag arrays.
	IF data_type EQ 0 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/DOUBLE,/NOZERO), $
			imag: make_array(MAT_header.mrows,MAT_header.ncols,/DOUBLE,/NOZERO)} $
	ELSE IF data_type EQ 1 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/FLOAT,/NOZERO), $
			imag: make_array(MAT_header.mrows,MAT_header.ncols,/FLOAT,/NOZERO)} $
	ELSE IF data_type EQ 2 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/LONG,/NOZERO), $
			imag: make_array(MAT_header.mrows,MAT_header.ncols,/LONG,/NOZERO)} $
	ELSE IF data_type EQ 3 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/INT,/NOZERO), $
			imag: make_array(MAT_header.mrows,MAT_header.ncols,/INT,/NOZERO)} $
	ELSE IF data_type EQ 4 THEN $
	    MAT_data = {real: make_array(MAT_header.mrows,MAT_header.ncols,/INT,/NOZERO), $
			imag: make_array(MAT_header.mrows,MAT_header.ncols,/INT,/NOZERO)} $
	ELSE BEGIN
	    RETURN, ['load_mat_data','Error in MATLAB data type.']
	ENDELSE
    ENDELSE
	
    ; read the matrix
    READU, MAT_lun, MAT_data

    ; find if the matrix is required.
    argno = WHERE(MAT_names EQ MAT_name)

    IF argno(0) NE -1 THEN BEGIN

	; determine the argument into which this matrix should be stored.
	return_arg = 'arg' + STRTRIM(STRING(argno(0)+1),2)

	IF MAT_header.imagf EQ 1 THEN BEGIN
	    ; complex data
	    IF EXECUTE(return_arg+'=complex(MAT_data.real,MAT_data.imag)') NE 1 THEN $
		RETURN, ['load_mat_data','Error converting MATLAB data to IDL.']
	ENDIF ELSE IF text EQ 1 THEN BEGIN
	    ; text data
	    IF N_ELEMENTS(MAT_data.real) GT 1 THEN BEGIN
		IF EXECUTE(return_arg+'=STRING(BYTE(TRANSPOSE(MAT_data.real)))') NE 1 THEN $
		    RETURN, ['load_mat_data','Error converting MATLAB data to IDL.']
	    ENDIF ELSE BEGIN
		IF EXECUTE(return_arg+'=STRING(BYTE(MAT_data.real))') NE 1 THEN $
		    RETURN, ['load_mat_data','Error converting MATLAB data to IDL.']
	    ENDELSE
	ENDIF ELSE BEGIN
	    ; real data
	    IF EXECUTE(return_arg+'=MAT_data.real') NE 1 THEN $
		RETURN, ['load_mat_data','Error converting MATLAB data to IDL.']
	ENDELSE
    ENDIF
    
ENDWHILE

FREE_LUN,MAT_lun

RETURN, ''

bad_read:
    RETURN, ['load_mat_data','Error reading MATLAB file :',!ERR_STRING]

END
