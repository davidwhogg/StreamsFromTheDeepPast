;+
; NAME:
;   find_all_midpoints
; PURPOSE:
;   Locate all midpoints between all skew lines, where the midpoint is within
;     thresh of both lines
;-
pro find_all_midpoints, thresh,ilist,vlist

; read in hip data
restore, '/global/data/hipparcos/hip-10-200-1.sav'
nhip= n_elements(hip)

; shuffle (and subsample?)
;seed= -1L
;hip= hip[shuffle(seed,hip)]
;nhip=1000

; loop over all pairs, finding close ones
maxlist= 50000000L
ilist= lonarr(2,maxlist)
vlist= dblarr(3,maxlist)
nlist= 0L
for ii=0L,nhip-1L do begin
    splog, ii
    for jj=0L,nhip-1L do if ii NE jj then begin
        v= find_midpoint(hip[ii],hip[jj],dv=dv)
	if dv LT thresh and nlist LT maxlist then begin
            ilist[*,nlist]= [ii,jj]
            vlist[*,nlist]= v
            nlist= nlist+1L
        endif
    endif
endfor
ilist= ilist[*,0:nlist-1]
vlist= vlist[*,0:nlist-1]
help, ilist,vlist

save, hip,ilist,vlist,filename='find_all_midpoints.sav'
end


