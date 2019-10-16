;+
;   NAME:
;      bovy_determ
;   PURPOSE:
;      version of determ that can handle scalar inputs
;   INPUT:
;      matrix - matrix to evaluate the determinant of, can be a number
;   OUTPUT:
;      determinant of the matrix
;   HISTORY:
;      2009-11-17 - Written - Bovy
;-
FUNCTION BOVY_DETERM, matrix, double=double, check=check
IF n_elements(matrix) EQ 1 THEN return, matrix ELSE return, determ(matrix,double=double,check=check)
END
