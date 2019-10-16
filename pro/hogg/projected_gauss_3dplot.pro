;+
; NAME:
;   projected_gauss_3dplot
; PURPOSE:
;   three-d projections of n-d gaussians
; INPUTS:
; OUTPUTS:
; BUGS:
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro projected_gauss_3dplot, ydata, ycovar, projection, amp, xmean, xcovar, $
                            basename=basename,nptsring=nptsring, $
                            nrings=nrings, rotation=rotation, $
                            linelength=linelength, xbestguess=xbestguess, $
                            interact=interact,mkinfile=mkinfile, $
                            nogauss=nogauss

if(NOT keyword_set(linelength)) then linelength=1.
if(keyword_set(interact)) then mkinfile=1

ngauss=n_elements(amp)
ndimy=n_elements(ycovar)/n_elements(ydata)
ndimx=n_elements(xcovar)/n_elements(xmean)
ndata=n_elements(ydata)/ndimy

openw,punit,basename+'.pts',/get_lun
openw,iunit,basename+'.indx',/get_lun
nsig=1.
offset=0L

if(NOT keyword_set(nogauss)) then begin
    for jj=0L, ngauss-1L do begin
        threemesh, xmean[*,jj], xcovar[*,*,jj], nsig, linepos, lineindx, $
          nptsring=nptsring,nrings=nrings,offset=offset
        writeu,punit,linepos
        writeu,iunit,lineindx
        offset=offset+nptsring*nrings
    endfor
endif
    
linepos=fltarr(3,2)
linedir=[1.,0.,0.]
for ii=0L, ndata-1L do begin
    if(n_elements(rotation) gt 0) then begin
        linedir=rotation[*,*,ii]#[1.,0.,0.]
    endif
    if(n_elements(xbestguess) gt 0) then $
      zerodatamid=xbestguess[*,ii] $
    else $
      zerodatamid=projection[*,*,ii]#ydata[*,ii]
    linepos[*,0]=zerodatamid+linedir*linelength
    linepos[*,1]=zerodatamid-linedir*linelength
    lineindx=[2, offset+[0, 1]]
    writeu,punit,linepos
    writeu,iunit,lineindx
    offset=offset+2
endfor

free_lun,punit
free_lun,iunit

if(keyword_set(mkinfile)) then begin
    instr=strarr(8)
    i=0
    instr[i]='linepositionfile '+basename+'.pts' & i=i+1
    instr[i]='lineindexfile '+basename+'.indx' & i=i+1
    instr[i]='linecolor 1. 1. 1. 1.' & i=i+1
    instr[i]='scale 100.' & i=i+1
    instr[i]='windowsize 700' & i=i+1
    instr[i]='axis0 1.01    0.0 0.0 0.0     50. 0. 0.      1. 0. 0.' & i=i+1
    instr[i]='axis1 1.01    0.0 0.0 0.0     0. 50. 0.      1. 1. 1.' & i=i+1
    instr[i]='axis2 1.01    0.0 0.0 0.0     0. 0. 50.      0. 0. 1.' & i=i+1
    openw,unit,basename+'.in',/get_lun
    for i=0, n_elements(instr)-1 do $
      printf,unit,instr[i]
    free_lun,unit

    if(keyword_set(interact)) then begin
        spawn,'points '+basename+'.in &'
    endif
endif

end
