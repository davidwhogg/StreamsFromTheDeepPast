;+
; BUGS:
;  - Deg -> mas issues in plots.
;  - Needs proper header.
;  - Switch to fitting parabola to lowest 3 values, not lowest plus
;    each side.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_test
prefix= 'faintpm_test'
nreal= 700

spawn, '\rm -vf *.fits*'
faintpm_makefake, nreal

foo= {faintpm, $
      epoch:      0D0, $ ; yr
      trueflux:   0D0, $ ; image data numbers
      truera:     0D0, $ ; deg
      truedec:    0D0, $ ; deg
      trueradot:  0D0, $ ; mas/yr
      truedecdot: 0D0, $ ; mas/yr
      nepoch:       0, $ ;
      totalsnr:   0.0, $ ;
      meansnr:    0.0, $ ;
      meanfwhm:   0.0, $ ; deg
      meant:      0.0, $ ; yr
      sigmat:     0.0, $ ; yr
      maxdt:      0.0, $ ; yr
      fluxerr:    0D0, $ ; fractional
      fluxunc:   -1D0, $ ; fractional
      raerr:      0D0, $ ; mas
      raunc:     -1D0, $ ; mas
      decerr:     0D0, $ ; mas
      decunc:    -1D0, $ ; mas
      radoterr:   0D0, $ ; mas/yr
      radotunc:  -1D0, $ ; mas/yr
      decdoterr:  0D0, $ ; mas/yr
      decdotunc: -1D0  $ ; mas/yr
     }
tstr= replicate(foo,nreal)

set_plot,'ps'
hogg_plot_defaults
xsize= 6.0
ysize= xsize
device, filename=prefix+'.ps',xsize=xsize,xoffset=0.5*(8.5-xsize), $
  ysize=ysize,yoffset=0.5*(11.0-ysize), $
  /inches
!P.MULTI= [0,7,7]
!P.CHARSIZE= 0.75
!X.MARGIN= [1,0.5]*!X.OMARGIN
!X.OMARGIN= [0,0]
!Y.MARGIN= [1,0.5]*!Y.OMARGIN
!Y.OMARGIN= [0,0]

for ii=0L,nreal-1L do begin
    iistr= string(ii,format='(I3.3)')
    file= (findfile('./fake_'+iistr+'.fits.gz'))[0]

; gather truth
    hdr= headfits(file)
    nepoch= sxpar(hdr,'NEPOCH')
    tstr[ii].nepoch= nepoch
    tstr[ii].trueflux= sxpar(hdr,'FLUX')
    tstr[ii].epoch= sxpar(hdr,'MEANT')
    raref= double(sxpar(hdr,'RAREF'))
    decref= double(sxpar(hdr,'DECREF'))
    truedramean= sxpar(hdr,'DRAMEAN')
    trueddecmean= sxpar(hdr,'DDECMEAN')
    tstr[ii].truera= raref+truedramean
    tstr[ii].truedec= decref+trueddecmean
    trueradot= sxpar(hdr,'RADOT')
    tstr[ii].trueradot= trueradot
    truedecdot= sxpar(hdr,'DECDOT')
    tstr[ii].truedecdot= truedecdot

; gather data statistics
    fwhmlist= fltarr(nepoch)
    snrlist= fltarr(nepoch)
    tlist= fltarr(nepoch)
    scalelist= fltarr(nepoch)
    tname= 'TT'
    dd= 0L
    for jj=0L,nepoch-1L do begin
        exten= 2*jj+1
        hdrj= headfits(file,exten=exten)
        fwhmlist[jj]= sxpar(hdrj,'FWHM')
        snrlist[jj]= sxpar(hdrj,'SNR')
        tlist[jj]= sxpar(hdrj,tname)
        scalelist[jj]= sqrt(abs(sxpar(hdrj,'CD1_1')*sxpar(hdrj,'CD2_2') $
                                -sxpar(hdrj,'CD2_1')*sxpar(hdrj,'CD1_2')))
    endfor
    totalsnr2= total(snrlist^2)
    totalsnr= sqrt(totalsnr2)
    tstr[ii].totalsnr= totalsnr
    tstr[ii].meansnr= sqrt(mean(snrlist^2))
; in next line multiply by scale because fwhm is given in pixels
    tstr[ii].meanfwhm= total(scalelist*fwhmlist*snrlist^2)/totalsnr2
    scale= total(scalelist*snrlist^2)/totalsnr2
    tstr[ii].meant= total(tlist*snrlist^2)/totalsnr2
    sigmat= sqrt(total((tlist-tstr[ii].meant)^2*snrlist^2)/totalsnr2)
    tstr[ii].sigmat= sigmat
    maxdt= max(tlist)-min(tlist)
    tstr[ii].maxdt= maxdt
    theoryfluxerr= tstr[ii].trueflux/tstr[ii].totalsnr
    theoryposerr= tstr[ii].meanfwhm/tstr[ii].totalsnr
    theorydoterr= tstr[ii].meanfwhm/(tstr[ii].totalsnr*tstr[ii].sigmat) $
      *3.6D6                    ; mas
    realabel= string(tstr[ii].nepoch,format='(I2)')+' epochs spanning ' $
      +string(tstr[ii].maxdt,format='(F4.2)')+' years with mean !8S/N!3 =' $
      +string(tstr[ii].meansnr,format='(F5.2)')+' per epoch (combined ' $
      +strtrim(string(tstr[ii].totalsnr,format='(F6.2)'),2)+')'
    splog, realabel

; make small grid around truth
    nhalf= 3
    ns= 2.5
    grid= ns*(dindgen(2*nhalf+1)-nhalf)/nhalf
    ra=   raref+truedramean+theoryposerr*grid
    dec=  decref+trueddecmean+theoryposerr*grid
    radot=  trueradot+theorydoterr*grid
    decdot= truedecdot+theorydoterr*grid
    flux= tstr[ii].trueflux+theoryfluxerr*grid
    tunit= 1.0
    dcs= faintpm_verify_grid(ra,dec,radot,decdot,flux,tstr[ii].epoch, $
                             tname,file,tunit=tunit)

; label and plotting crap
    vec= [[ra],[dec],[radot],[decdot],[flux]]
    true= [tstr[ii].truera,tstr[ii].truedec, $
           tstr[ii].trueradot,tstr[ii].truedecdot,tstr[ii].trueflux]
    scaling= [1.0/3.6E6,1.0/3.6E6,1.0,1.0,tstr[ii].trueflux]
    theory= [theoryposerr,theoryposerr,theorydoterr,theorydoterr,theoryfluxerr]
    fitvec= fltarr(5)
    fitvecerr= fltarr(5)-1.0
    label= ['RA error [mas]', $
            'Dec error [mas]', $
            'dRA/d!8t!3 error [mas/yr]', $
            'dDec/d!8t!3 error [mas/yr]', $
            'fractional source flux error']
    hogg_usersym, 12,/fill,scale=0.25
    !Y.TITLE= '!7Dv!3!u2!n'
    !Y.RANGE= [-1,11]

; plots
    for jj=0,4 do begin
        thisvec= (vec[*,jj]-true[jj])/scaling[jj]
        thisdcs= faintpm_marginalize(dcs,jj+1)
        !X.TITLE= label[jj]
        plot, thisvec,thisdcs,psym=8
        oplot, [-1,-1]*theory[jj]/scaling[jj],!Y.CRANGE,linestyle=2
        oplot, [1,1]*theory[jj]/scaling[jj],!Y.CRANGE,linestyle=2
        three= (sort(thisdcs))[0:2]
        xxfine= min(thisvec)+(max(thisvec)-min(thisvec)) $
          *findgen(101)/100.0
        yyfine= faintpm_parabola(thisvec[three],thisdcs[three], $
                                 xxfine,xmin,ymin,sigmax)
        if (sigmax GT 0) then begin
            fitvec[jj]= xmin
            fitvecerr[jj]= sigmax
            oplot, [1,1]*xmin,!Y.CRANGE
            oplot, xxfine,yyfine
        endif
        if (jj EQ 0) then begin
            xl= !X.CRANGE[0]
            yl= !Y.CRANGE[1]+0.05*(!Y.CRANGE[1]-!Y.CRANGE[0])
            xyouts, xl,yl,realabel
        endif
    endfor

; save results
    tstr[ii].raerr=     fitvec[0]
    tstr[ii].decerr=    fitvec[1]
    tstr[ii].radoterr=  fitvec[2]
    tstr[ii].decdoterr= fitvec[3]
    tstr[ii].fluxerr=   fitvec[4]
    tstr[ii].raunc=     fitvecerr[0]
    tstr[ii].decunc=    fitvecerr[1]
    tstr[ii].radotunc=  fitvecerr[2]
    tstr[ii].decdotunc= fitvecerr[3]
    tstr[ii].fluxunc=   fitvecerr[4]

; perform marginalization over flux
    foo= dcs
    foo= foo-min(foo)
    foo= -2.0*alog(total(exp(-0.5*foo),5,/double))
    dcs= foo-min(foo)

; set up to plot ra,dec by marginalizing over dots
    foo= total(total(exp(-0.5*dcs),4,/double),3,/double)
    posdcs= -2.0*alog(foo)
    posdcs= posdcs-min(posdcs)
    ragrid= (ra#replicate(1D0,n_elements(dec)))-raref
    decgrid= (replicate(1D0,n_elements(ra))#dec)-decref

; plot radot,decdot
    !X.TITLE= '[RA - reference] (pix)'
    !Y.TITLE= '[Dec - reference] (pix)'
    levels= [1,4]
    c_thick= [4,2]
    theoryvar= (theoryposerr/scale)^2
    !X.RANGE= 0
    !Y.RANGE= 0
    contour, posdcs,ragrid/scale,decgrid/scale, $
      levels=levels,c_thick=c_thick
    hogg_oplot_covar, [truedramean/scale],[trueddecmean/scale], $
      [[theoryvar,0],[0,theoryvar]],linestyle=2

; set up to plot radot,decdot by marginalizing over ra,dec
    foo= total(total(exp(-0.5*dcs),1,/double),1,/double)
    dotdcs= -2.0*alog(foo)
    dotdcs= dotdcs-min(dotdcs)
    radotgrid= (radot#replicate(1D0,n_elements(decdot)))
    decdotgrid= (replicate(1D0,n_elements(radot))#decdot)

; plot radot,decdot
    theoryvar= (theorydoterr/(scale/maxdt))^2
    !X.TITLE= 'dRA/d!8t!3 (pix/<!7r!d!8t!3!n>)'
    !X.RANGE= 0
    !Y.TITLE= 'dDec/d!8t!3 (pix/<!7r!d!8t!3!n>)'
    !Y.RANGE= 0
    contour, dotdcs,radotgrid/(scale/sigmat),decdotgrid/(scale/sigmat), $
      levels=levels,c_thick=c_thick
    oplot, [trueradot/(scale/sigmat)],[truedecdot/(scale/sigmat)], $
      psym=1,thick=8,symsize=2.0*!P.SYMSIZE
    theoryvar= (theorydoterr/(scale/sigmat))^2
    hogg_oplot_covar, [trueradot/(scale/sigmat)],[truedecdot/(scale/sigmat)], $
      [[theoryvar,0],[0,theoryvar]],linestyle=2

endfor
device,/close

; write fits file
mwrfits, tstr,prefix+'.fits',/create

; make plots
faintpm_test_plot
return
end
