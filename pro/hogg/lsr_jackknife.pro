;+
; BUGS:
;   No comment header.
;   The structure of this code SUX.
;-
pro lsr_jackknife

; make savesets
;for ng=1,5 do begin
;    lsr, ngaussfordisk=ng
;endfor

; read savesets -- this is crazy because it changes variable values!
lsrj_ngaussfordisk= lindgen(5)+1
nng= n_elements(lsrj_ngaussfordisk)
for iiii=0,nng-1 do begin
    filename= 'lsr_'+strtrim(string(lsrj_ngaussfordisk[iiii]),2)+'.sav'
    restore, filename
    if not keyword_set(lsrjack) then $
      lsrjack= dblarr([size(jacklike,/dimensions),nng])
    lsrjack[*,*,iiii]= jacklike
endfor
help, lsrjack
dimens= size(lsrjack,/dimensions)
nsubsamp= dimens[0]
njack= dimens[1]
nng= dimens[2]
meanlsrjack= total(lsrjack,2)/double(njack)
for ii=0,nsubsamp-1 do $
  meanlsrjack[ii,*]= meanlsrjack[ii,*]-meanlsrjack[ii,0]
help, meanlsrjack

; plot likelihoods
set_plot, 'PS'
xsize= 6.0
ysize= 6.0
device, file='lsr_jackknife.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
hogg_plot_defaults
plot, [0],[0],/nodata, $
  xstyle=1, $
  xrange= [0,nsubsamp]+[-0.5,-0.5], $
  xtitle='subsample', $
  ystyle=1, $
  yrange= [-1.1,1.1]*max(meanlsrjack), $
  ytitle='log likelihood ratio vs 1-gaussian disk'
for ii=0,nng-1 do begin
    oplot, dindgen(nsubsamp),meanlsrjack[*,ii],psym=10,color=djs_icolor('grey')
    xyouts, dindgen(nsubsamp),meanlsrjack[*,ii], $
      strtrim(string(lsrj_ngaussfordisk[ii]),2),align=0.5
endfor
device,/close

tmp_bestngaussfordisk= lonarr(nsubsamp)
for ii=0,nsubsamp-1 do $
  tmp_bestngaussfordisk[ii]=  $
  lsrj_ngaussfordisk[(reverse(sort(meanlsrjack[ii,*])))[0]]
for iiii=0,nng-1 do begin
    tmp_ngaussfordisk= lsrj_ngaussfordisk[iiii]
    filename= 'lsr_'+strtrim(string(tmp_ngaussfordisk),2)+'.sav'
    restore, filename
    if (iiii EQ 0) then begin
        tmp_parameter=  parameter
        tmp_vdisp2=     vdisp2
        tmp_dparameter= dparameter
        tmp_dcolor=     dcolor
        tmp_vertexdev=  vertexdev
        tmp_dvertexdev= dvertexdev
        tmp_mean=       mean
        tmp_dmean=      dmean
    endif else begin
        indx= where(tmp_bestngaussfordisk EQ tmp_ngaussfordisk,nindx)
        if (nindx GT 0) then begin
            tmp_parameter[indx]=  parameter[indx]
            tmp_vdisp2[indx]=     vdisp2[indx]
            tmp_dparameter[indx]= dparameter[indx]
            tmp_dcolor[indx]=     dcolor[indx]
            tmp_vertexdev[indx]=  vertexdev[indx]
            tmp_dvertexdev[indx]= dvertexdev[indx]
            tmp_mean[indx]=       mean[indx]
            tmp_dmean[indx]=      dmean[indx]
        endif
    endelse
endfor
parameter=  temporary(tmp_parameter)
vdisp2=     temporary(tmp_vdisp2)
dparameter= temporary(tmp_dparameter)
dcolor=     temporary(tmp_dcolor)
vertexdev=  temporary(tmp_vertexdev)
dvertexdev= temporary(tmp_dvertexdev)
mean=       temporary(tmp_mean)
dmean=      temporary(tmp_dmean)
bestngaussfordisk= temporary(tmp_bestngaussfordisk)
prefix= 'lsr_bestngauss'
savefilename= 'lsr_bestngauss.sav'
save, filename=savefilename

end
