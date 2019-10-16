pro hogg_get_data_pal5
prefix= 'pal5_sdss_data'
ra0= 229.022                    ; deg J2000
dec0= -0.111                    ; deg J2000
rerun= [137,161,157]
nrerun= n_elements(rerun)

; find runs
res= [sdss_findimage(ra0,dec0), $
      sdss_findimage(ra0,dec0+0.05), $
      sdss_findimage(ra0,dec0-0.05)]
identifier= res.run
sindx= sort(identifier)
res= res[sindx]
identifier= identifier[sindx]
uindx= uniq(identifier)
res= res[uindx]
identifier= identifier[uindx]
nid= n_elements(identifier)
mjd= sdss_run2mjd(res.run)
size= 1.0                       ; deg
bigast = hogg_make_astr(ra0,dec0,size,size)

; set up plot
set_plot, 'ps'
xsize= 7.5 & ysize= 7.5
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
hogg_plot_defaults
degrange= [-0.5,0.5]*size
!X.RANGE= ra0+reverse(degrange)/cos(!PI*dec0/180.0)
!X.TITLE= 'J2000 RA (deg)'
!Y.RANGE= dec0+degrange
!Y.TITLE= 'J2000 Dec (deg)
!P.MULTI= [0,3,3]
superthick= 4.0*!P.CHARTHICK

; loop over runs, plotting
for ii=0L,nid-1L do begin
    splog, ii,res[ii].run
    plot, [0],[0],/nodata
    titlestr= 'run '+strtrim(string(res[ii].run),2) $
      +', mjd '+strtrim(string(mjd[ii]),2)
    xyouts, !X.CRANGE[0],!Y.CRANGE[0],titlestr,charthick=superthick
    for camcol= ((res[ii].camcol-1)>1),((res[ii].camcol+1)<6) do begin
        for field= ((res[ii].field-3)>1),((res[ii].field+3)) do begin
            objs= sdss_readobj(res[ii].run,camcol,field,rerun=rerun)
            if n_elements(objs) GT 1 then oplot, objs.ra,objs.dec,psym=3
        endfor
    endfor
    hogg_oplot_covar, [ra0],[dec0],[[[(5.0/60)^2,0],[0,(5.0/60)^2]]], $
      thick=superthick
endfor

device,/close
return
end
