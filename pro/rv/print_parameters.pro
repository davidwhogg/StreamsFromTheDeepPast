;+
;   NAME:
;      print_parameters
;   PURPOSE:
;      print the parameters of the Gaussians to a LaTeX table
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      filename - filename for table
;      basedir - basedirectory for density estimate (ends in /)
;   KEYWORDS:
;      /QA     - Print QA info
;   OUTPUT:
;      table
;   REVISION HISTORY:
;      2009-01-18 - Written Bovy (NYU)
;-
PRO PRINT_PARAMETERS, K=K, w=w, sample=sample, filename=filename, basedir=basedir, $
                      QA=QA

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'

;;Restore density estimate
densityfile=basedir+'fitv'
densityfile+= strtrim(string(K),2)
densityfile+='_bestV_V'
densityfile+=strtrim(string(w,format='(f4.1)'),2)
densityfile+= '_sample'
densityfile+= strtrim(string(sample),2)
densityfile+= '.sav'
IF file_test(densityfile) THEN BEGIN
    IF keyword_set(QA) THEN splog, 'reading '+densityfile
    restore, densityfile
ENDIF ELSE BEGIN
    splog, 'Density estimate not found in file '+densityfile
    splog, 'Returning...'
    RETURN
ENDELSE

;;Sort the Gaussians by amplitude
sort_indx= sort(amp)
amp= amp[reverse(sort_indx)]
FOR ii= 0, 2 DO mean[ii,*]= mean[ii,reverse(sort_indx)]
FOR ii= 0, 2 DO FOR jj= 0, 2 DO covar[ii,jj,*]= covar[ii,jj,reverse(sort_indx)]

;Set-up LaTeX file
vunit= 'km s$^{-1}$'
v2unit= 'km$^2$ s$^{-2}$'
OPENW, wlun, filename, /GET_LUN
PRINTF, wlun, '\clearpage'
PRINTF, wlun, '\begin{deluxetable}{cr@{.}lr@{.}lr@{.}lr@{.}lr@{.}lr@{.}lr@{.}lr@{.}lr@{.}lr@{.}l}'
PRINTF, wlun, '\tablecaption{Parameters of the component '+$
  'Gaussians for the density estimate with K = '+strtrim(string(K),2)+' and w = '+strtrim(string(long(w)),2)+' km$^2$ s$^{-2}$.\label{table:param}}'
PRINTF, wlun, '\tablecolumns{21}'
PRINTF, wlun, '\tablewidth{0pt}'
PRINTF, wlun, '\rotate'
PRINTF, wlun, '\tablehead{\colhead{\ngauss}'+$
  ' & \multicolumn{2}{c}{\alphagauss} & \multicolumn{2}{c}{$\eex\T\,\vv$} & '+$
  '\multicolumn{2}{c}{$\eey\T\,\vv$} & \multicolumn{2}{c}{$\eez\T\,\vv$} & '+$
  '\multicolumn{2}{c}{$\eex\T\,\VV\,\eex$} &\multicolumn{2}{c}{$\eey\T\,\VV\,\eey$} & '+$
  '\multicolumn{2}{c}{$\eez\T\,\VV\,\eez$} &\multicolumn{2}{c}{$\eex\T\,\VV\,\eey$} &'+$
  '\multicolumn{2}{c}{$\eex\T\,\VV\,\eez$} & \multicolumn{2}{c}{$\eey\T\,\VV\,\eez$} \\'
PRINTF, wlun, '\colhead{} & \multicolumn{2}{c}{} & '+$
  '\multicolumn{2}{c}{('+ vunit+')} & \multicolumn{2}{c}{('+ vunit+')} & \multicolumn{2}{c}{('+ $
  vunit+')} & '+ '\multicolumn{2}{c}{('+ v2unit+')} & \multicolumn{2}{c}{('+ $
  v2unit+')} & \multicolumn{2}{c}{('+ v2unit+')} &'+ ' \multicolumn{2}{c}{('+ $
  v2unit+')} & \multicolumn{2}{c}{('+ v2unit+')} & \multicolumn{2}{c}{('+ v2unit+')}}'
PRINTF, wlun, '\startdata'
;loop over subsamples
parameterformat=strarr(20)
parameterformat[0]='(I)'
parameterformat[1]='(I04)'
parameterformat[2]='(I)'
parameterformat[3]='(I02)'
parameterformat[4]=parameterformat[2]
parameterformat[5]=parameterformat[3]
parameterformat[6]=parameterformat[2]
parameterformat[7]=parameterformat[3]
parameterformat[8]=parameterformat[2]
parameterformat[9]=parameterformat[3]
parameterformat[10]=parameterformat[2]
parameterformat[11]=parameterformat[3]
parameterformat[12]=parameterformat[2]
parameterformat[13]=parameterformat[3]
parameterformat[14]=parameterformat[2]
parameterformat[15]=parameterformat[3]
parameterformat[16]=parameterformat[2]
parameterformat[17]=parameterformat[3]
parameterformat[18]=parameterformat[2]
parameterformat[19]=parameterformat[3];;Whatever
;Write parameters to file
FOR gg=0L, K-1 DO PRINTF, wlun, gg+1, $
  floor(amp[gg]),$
  round(1D4*(amp[gg]-floor(amp[gg]))),$
  long(mean[0,gg]),$
  round(1D2*abs((mean[0,gg]-long(mean[0,gg])))),$
  long(mean[1,gg]),$
  round(1D2*abs((mean[1,gg]-long(mean[1,gg])))),$
  long(mean[2,gg]),$
  round(1D2*abs((mean[2,gg]-long(mean[2,gg])))),$
  long(covar[0,0,gg]),$
  round(1D2*abs((covar[0,0,gg]-long(covar[0,0,gg])))),$
  long(covar[1,1,gg]),$
  round(1D2*abs((covar[1,1,gg]-long(covar[1,1,gg])))),$
  long(covar[2,2,gg]),$
  round(1D2*abs((covar[2,2,gg]-long(covar[2,2,gg])))),$
  long(covar[0,1,gg]),$
  round(1D2*abs((covar[0,1,gg]-long(covar[0,1,gg])))),$
  long(covar[0,2,gg]),$
  round(1D2*abs((covar[0,2,gg]-long(covar[0,2,gg])))),$
  long(covar[1,2,gg]),$
  round(1D2*abs((covar[1,2,gg]-long(covar[1,2,gg])))),$
  format='( "",'+'(I)'+ '," & " '+$
  parameterformat[0]+'," & ",'+parameterformat[1]+'," & ",'+$
  parameterformat[2]+'," & ",'+parameterformat[3]+'," & ",'+$
  parameterformat[4]+'," & ",'+parameterformat[5]+'," & ",'+$
  parameterformat[6]+'," & ",'+parameterformat[7]+'," & ",'+$
  parameterformat[8]+'," & ",'+parameterformat[9]+'," & ",'+$
  parameterformat[10]+'," & ",'+parameterformat[11]+'," & ",'+$
  parameterformat[12]+'," & ",'+parameterformat[13]+'," & ",'+$
  parameterformat[14]+'," & ",'+parameterformat[15]+'," & ",'+$
  parameterformat[16]+'," & ",'+parameterformat[17]+'," & ",'+$
  parameterformat[18]+'," & ",'+parameterformat[19]+',"\\")'
PRINTF, wlun, '\enddata'
;PRINTF, wlun, '\tablecomments{}'
PRINTF, wlun, '\end{deluxetable}'
FREE_LUN, wlun


END
