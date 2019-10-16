;+
; NAME:
;   threemesh
; PURPOSE:
;   build a three-d mesh out of a covariance matrix
; INPUTS:
; OUTPUTS:
; BUGS:
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro threemesh, mean, covar, nsig, linepos, lineindx, nptsring=nptsring, $
               nrings=nrings, offset=offset

if(NOT keyword_set(nptsring)) then nptsring=30
if(NOT keyword_set(nrings)) then nrings=30
if(NOT keyword_set(offset)) then offset=0L

symcovar=0.5*(covar+transpose(covar))
eigenval=eigenql(symcovar,eigenvectors=eigenvec)
rotate=fltarr(3,3)
for i=0L, 2L do $
  rotate[*,i]=eigenvec[*,i]*sqrt(eigenval[i])

linepos=fltarr(3,nptsring*nrings)
tmp_linepos=fltarr(3,nptsring*nrings)
lineindx=lonarr((nptsring+1)*nrings)
for k=0L, nrings-1L do begin
    theta=(float(k)+0.5)/float(nrings)*!DPI
    phi=dindgen(nptsring)/float(nptsring-1L)*2.*!DPI
    tmp_linepos[0,k*nptsring+lindgen(nptsring)]=sin(theta)*cos(phi)
    tmp_linepos[1,k*nptsring+lindgen(nptsring)]=sin(theta)*sin(phi)
    tmp_linepos[2,k*nptsring+lindgen(nptsring)]=cos(theta)
    lineindx[k*(nptsring+1)]=nptsring
    lineindx[k*(nptsring+1)+lindgen(nptsring)+1]= $
      k*nptsring+lindgen(nptsring)+offset
endfor
linepos=nsig*rotate#tmp_linepos
for i=0L, 2L do $
  linepos[i,*]=linepos[i,*]+mean[i]

end
