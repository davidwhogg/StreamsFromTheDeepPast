;+
; NAME:
;   hogg_smosaic_gc
; BUGS:
;   No header written.
;   HORRIBLE ascii-reading code
;-
pro hogg_read_mwgc, name,ra,dec
openr, lun,'../../data/mwgc.dat',/get_lun
line= ''
repeat begin
    readf, lun,line
endrep until ((stregex(line,'ID') GE 0) AND $
              (stregex(line,'Name') GE 0) AND $
              (stregex(line,'RA') GE 0) AND $
              (stregex(line,'DEC') GE 0))
readf, lun,line
readf, lun,line
readf, lun,line
while (stregex(line,'.') GE 0) do begin
    name1= strtrim(strmid(line,0,23),2)
    if keyword_set(name) then name= [name,name1] else name= name1
    rah= double(strmid(line,23,2))
    ram= double(strmid(line,26,2))
    ras= double(strmid(line,29,4))
    ra1= 1.5D1*(rah+ram/6D1+ras/3.6D3)
    if keyword_set(ra) then ra= [ra,ra1] else ra= ra1
    decsgn= strmid(line,35,1)
    decd= double(strmid(line,36,2))
    decm= double(strmid(line,39,2))
    decs= double(strmid(line,42,2))
    dec1= decd+decm/6D1+decs/3.6D3
    if (decsgn EQ '-') then dec1= -1D0*dec1
    if keyword_set(dec) then dec= [dec,dec1] else dec= dec1
    readf, lun,line
endwhile
close, lun
return
end
