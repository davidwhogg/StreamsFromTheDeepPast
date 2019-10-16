pro test_uber

; initialize postcript file
  axis_char_scale= 1.0
  set_plot, "PS"
  !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
  !P.THICK= 2.0
  !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
  !P.CHARSIZE= 1.0
  !P.PSYM= 0
  !P.MULTI= [0,6,6]
  !X.MARGIN= [5,0]*axis_char_scale
  !X.OMARGIN= [0,0]*axis_char_scale
  !X.STYLE= 3
  !X.RANGE= 0
  !X.CHARSIZE= !P.CHARSIZE*axis_char_scale
  !Y.MARGIN= 0.6*!X.MARGIN
  !Y.OMARGIN= 0.6*!X.OMARGIN
  !Y.STYLE= !X.STYLE
  !Y.RANGE= 0
  !Y.CHARSIZE= !X.CHARSIZE
  xsize= 6.0 & ysize= 6.0
  device, file="test_uber.ps",/inches,xsize=xsize,ysize=ysize, $
    xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
  xyouts, 0,0,'!3',/device

; initialize variables/inputs
  seed= 1L
  n_point= 30L
  weight= dblarr(n_point)+1d0
  minvar= 1d0
  maxvar= 1d2
  ntry= 36L

; begin loop over tries
  for ii=0L,ntry-1L do begin

; make external variance and mean
    var1= minvar*exp(alog(maxvar/minvar)*double(ii)/double(ntry-1))
    extvar= [[1D0,0.99999],[0.99999,1.0]]*var1
    extmean= 1d2*randomn(seed,2)

; make "true" data values
    truepoint= fake_catalog(seed,[1],extmean,extvar,ndata=n_point)
    point= dblarr(2,n_point)

; make measurement covariances
    covar= dblarr(2,2,n_point)
    for jj=0L,n_point-1L do begin
      for kk=0L,2 do begin
        x1= randomn(seed,2)
        covar[*,*,jj]= covar[*,*,jj]+transpose(x1)##x1
      endfor

; make measurement error
      point[*,jj]= fake_catalog(seed,[1],truepoint[*,jj],covar[*,*,jj],ndata=1,/quiet)
    endfor

; plot points
    djs_plot, extmean[0]+[-3,3]*sqrt(var1),extmean[1]+[-3,3]*sqrt(var1), $
      /nodata
    hogg_oplot_covar, point[0,*],point[1,*],covar,color='grey'

; uber
    uber, weight,point,covar,2,mean,vec
    print, vec
    hogg_oplot_covar, extmean[0],extmean[1],extvar, $
      thick=3*!P.THICK,color='black'
    hogg_oplot_covar, mean[0],mean[1],transpose(vec)##vec, $
      thick=3*!P.THICK,color='red'

  endfor
end

pro test_uber2

; initialize postcript file
  axis_char_scale= 1.0
  set_plot, "PS"
  !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
  !P.THICK= 4.0
  !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
  !P.CHARSIZE= 1.0
  !P.PSYM= 0
  !P.MULTI= [0,1,1]
  !X.MARGIN= [6,0]*axis_char_scale
  !X.OMARGIN= [0,0]*axis_char_scale
  !X.STYLE= 3
  !X.RANGE= 0
  !X.CHARSIZE= !P.CHARSIZE*axis_char_scale
  !Y.MARGIN= 0.6*!X.MARGIN
  !Y.OMARGIN= 0.6*!X.OMARGIN+[0,2]
  !Y.STYLE= !X.STYLE
  !Y.RANGE= 0
  !Y.CHARSIZE= !X.CHARSIZE
  xsize= 6.0 & ysize= 6.0
  device, file="test_uber.ps",/inches,xsize=xsize,ysize=ysize, $
    xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
  xyouts, 0,0,'!3',/device

; initialize variables/inputs
  seed= 1L
  n_point= 100L
  weight= dblarr(n_point)+1d0
  minvar= 1d-3
  maxvar= 1d2
  ntry= 100L
  nboot= 20L
  priorfactor= 1d6

; make axes and label figure with inputs
  djs_plot, alog10([minvar,maxvar]),alog10([minvar,maxvar]), $
    xstyle=3,xtitle='input log!d10!n external variance', $
    ystyle=3,ytitle='inferred log!d10!n external variance', $
    title='output vs input for 1-dimensional, 1-gaussian uber', $
    thick=0.5*!P.THICK
  xl= !X.CRANGE[0]+0.04*(!X.CRANGE[1]-!X.CRANGE[0])
  yl= !Y.CRANGE[1]-0.05*(!Y.CRANGE[1]-!Y.CRANGE[0])
  djs_xyouts, xl,yl,strtrim(string(n_point),2)+' data points in each trial'
  yl= !Y.CRANGE[1]-0.09*(!Y.CRANGE[1]-!Y.CRANGE[0])
  djs_xyouts, xl,yl,'measurment variances drawn from a tophat 0<!4r!3!u2!n<1'
  yl= !Y.CRANGE[1]-0.13*(!Y.CRANGE[1]-!Y.CRANGE[0])
  djs_xyouts, xl,yl,'prior on variance '+strtrim(string(priorfactor),2)+' times true value'
  yl= !Y.CRANGE[1]-0.17*(!Y.CRANGE[1]-!Y.CRANGE[0])
  djs_xyouts, xl,yl,strtrim(string(nboot),2)+' bootstrap resamples per trial'

; begin loop over tries
  for ii=0L,ntry do begin

; make external variance and mean
    extvar= minvar*exp(alog(maxvar/minvar)*double(ii)/double(ntry))
    extmean= 1d2*randomn(seed)

; make measurement covariances
    covar= randomu(seed,n_point)

; make measurements
    point= extmean+randomn(seed,n_point)*sqrt(extvar) $
      +randomn(seed,n_point)*sqrt(covar)

; run uber algorithm
    uber, weight,point,covar,1,mean,vec, $
      nboot=nboot,bootvec=bootvec
    var= vec^2
    bootvar= bootvec^2
    splog, extmean,mean,extvar,var,status
    hogg_usersym, 2,thick=!P.THICK,scale=0.5
    djs_oplot, dblarr(nboot)+alog10(extvar),alog10(bootvar),psym=8
    hogg_usersym, 12,scale=0.5,thick=!P.THICK
    djs_oplot, [alog10(extvar)],[alog10(var)],psym=8,color='white'
    hogg_usersym, 12,/fill,scale=0.5
    djs_oplot, [alog10(extvar)],[alog10(var)],psym=8
  endfor

  device,/close
end
