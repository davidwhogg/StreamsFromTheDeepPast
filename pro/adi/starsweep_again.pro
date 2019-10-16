
pro starsweep_again
;+ 
; check out ~dwh2/hoggpt/pro/hoggpt/ptcat_struct.pro
; write starsweep_stream_adi_struct.pro
;  ra, dec (doubles)
;  g r i (floats)
;  use extinction to compute g r i.
; replace everything= replicate... lines with
;         everything= starsweep_stream_adi_struct(10000000L)
; read up on objc_flags, flags, objc_type, etc.  photo.astro.princeton.edu
;Adi
;started: 2006-06-08
;-

files=findfile('/global/data/sdss/redux/datasweep/dr4uber/137/calibObj*star*fits*')
bb=n_elements(files)-1246

jj=0
kk=0
for ii=0,bb-1 do begin
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
                         And (Not(allinfo.flags[2] AND $
                               sdss_flagval('object1','satur'))) $
                         And (NOT (allinfo.flags[3] AND $
                               sdss_flagval('object1','satur'))),nbeststars)
        
        if (nbeststars GE 1)then begin
            allinfo= allinfo[beststars]
            final=where ((allinfo.flags[1] And $
                          sdss_flagval('object1','binned1')) $
                         and (allinfo.flags[2] and $
                              sdss_flagval('object1','binned1')) $
                         and (allinfo.flags[3] and $
                              sdss_flagval('object1','binned1')),nfinal)
            if (nfinal GE 1) then begin
                allinfo=allinfo[final]
                
                
                str1=create_struct('ra',0.D,'dec',0.D,'gmag',0.,'rmag',0.,'imag',0.)
                starsweep_struct=replicate(str1,nfinal)
                starsweep_struct.ra=allinfo.ra
                starsweep_struct.dec=allinfo.dec
                
                                ;struct_assign,starsweep_struct,allinfo,/nozero
                starsweep_struct.gmag=22.5-2.5* $
                  alog10(allinfo.psfflux[1]-allinfo.extinction[1])           
                starsweep_struct.rmag=22.5-2.5* $
                  alog10(allinfo.psfflux[2]-allinfo.extinction[2])             
                starsweep_struct.rmag=22.5-2.5* $
                  alog10(allinfo.psfflux[3]-allinfo.extinction[3])             
                
                
                                ;  stop
                if (not keyword_set(everything)) then begin
                    
                    everything=starsweep_struct
                    everything= replicate(starsweep_struct[0],3000000L)
                    neverything= n_elements(everything)
                endif
                kk= jj+nfinal
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
    
endfor
mwrfits,everything[0:jj-1],"allinfo.fits",/create



end
