
pro starsweep_stream_hogg
;+ 
; check out ~dwh2/hoggpt/pro/hoggpt/ptcat_struct.pro
; write starsweep_stream_adi_struct.pro
;  ra, dec (doubles)
;  g r i (floats)
;  use extinction to compute g r i.
; replace everything= replicate... lines with
;         everything= starsweep_stream_adi_struct(10000000L)
; read up on objc_flags, flags, objc_type, etc.  photo.astro.princeton.edu
;******i think syntax if wrong*********
;Adi
;started: 2006-06-08
; BUGS:
;   - No correct comment header!
;-

for camcol= 1,6 do begin
    ccstr= string(camcol,format='(I1)')
    tfiles=findfile('/global/data/sdss/redux/datasweep/full_06oct05/137/calibObj-*-'+ccstr+'-star.fits*',count=bb)
    splog, bb
    if keyword_set(files) then files=[files,tfiles] else files=tfiles
endfor
bb=n_elements(files)

jj=0L
kk=0L
for ii=0L,bb-1L do begin
    allinfo=mrdfits(files[ii],1, $
                    columns=['RA','DEC','OBJC_FLAGS','OBJC_FLAGS2', $
                             'RESOLVE_STATUS','flags','flags2',$
                             'Extinction', $
                             'PSFFlux'])
;,'PSFFlux_Ivar']) 
                                ;'OBJC_TYPE','OBJC_ROWC'])
    goodstars=where(allinfo.resolve_status AND $
                    sdss_flagval('Resolve_status','survey_primary'),ngoodstar)
    
    if (ngoodstar GE 1) then begin
        allinfo= allinfo[goodstars]
        beststars=where ((Not (allinfo.flags[1] AND $
                               sdss_flagval('object1','satur'))) $
                         And (Not (allinfo.flags[2] AND $
                               sdss_flagval('object1','satur'))) $
                         And (NOT (allinfo.flags[3] AND $
                               sdss_flagval('object1','satur'))),nbeststar)
        
        if (nbeststar GE 1)then begin
            allinfo= allinfo[beststars]
            final=where ((allinfo.flags[1] And $
                          sdss_flagval('object1','binned1')) $
                         and (allinfo.flags[2] and $
                              sdss_flagval('object1','binned1')) $
                         and (allinfo.flags[3] and $
                              sdss_flagval('object1','binned1')),nfinal)
            
            if (nfinal GE 1) then begin
                allinfo=allinfo[final]
                decrafinal=where((allinfo.ra gt 110) and $
                                 (allinfo.ra lt 190) and $
                                 (allinfo.dec gt -20) and $
                                 (allinfo.dec lt 85),ndecra)
                if (ndecra GE 1) then begin
                    allinfo=allinfo[decrafinal]
                
                str1=create_struct('ra',0.D,'dec',0.D,'gmag',0.,'rmag',0.,'imag',0.)
                starsweep_struct=replicate(str1,ndecra)
                starsweep_struct.ra=allinfo.ra
                starsweep_struct.dec=allinfo.dec
                
                                ;struct_assign,starsweep_struct,allinfo,/nozero
                starsweep_struct.gmag=22.5-(2.5* $
                  alog10(allinfo.psfflux[1]))-allinfo.extinction[1]           
                starsweep_struct.rmag=22.5-(2.5* $
                  alog10(allinfo.psfflux[2]))-allinfo.extinction[2]             
                starsweep_struct.imag=22.5-(2.5* $
                  alog10(allinfo.psfflux[3]))-allinfo.extinction[3]             
                
                
                                ;  stop
                if (not keyword_set(everything)) then begin
                    
                    everything=starsweep_struct
                    everything= replicate(starsweep_struct[0],3000000L)
                    neverything= n_elements(everything)
                endif
                kk= jj+ndecra
                if kk GT neverything then begin
                    everything= [everything,replicate(starsweep_struct[0],10000000L)]
;allinfo[0],3000000L)]
                    neverything= n_elements(everything)
                endif
                everything[jj:kk-1]=starsweep_struct ; allinfo
                jj= kk
            endif
        endif
    endif
    endif
endfor
;qhelp
mwrfits2,everything[0:jj-1],"/global/data/scr/adi/newdata.fits",/create



end
