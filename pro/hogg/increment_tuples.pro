;+
; NAME:
;   increment_tuples
; PURPOSE:
;   take (N-1) tuples and promote them to N tuples
; INPUTS:
;   viv       - vector of (mean^T . covar^(-1} . mean)
;   iv        - vector of (covar^(-1} . mean)
;   icov      - vector of (covar^(-1))
;   ntup_list - list of index lists for (N-1)-tuples
;   nsigma    - threshold
; OUTPUTS:
;   ntup_list - updated list of index lists of N-tuples
; OPTIONAL OUTPUTS:
;   failed    - binary yes/no for failure
; BUGS:
;   Uses array concatenation which is SLOW-ass!
; REVISION HISTORY:
;   2003-02-18  written - Roweis and Hogg
;-
function increment_tuples, viv,iv,icov,ntup_list,nsigma, $
                           failed=failed

; figure out ntup etc
dims= size(ntup_list,/dimensions)
if n_elements(dims) EQ 1 then begin
    nlist= n_elements(ntup_list)
    ntup_list= reform(ntup_list,1,nlist)
    dims= size(ntup_list,/dimensions)
endif
ntup= dims[0]+1
nlist= dims[1]

; make a list of all unique indices in ntup_list
uniq_list= ntup_list[uniq(ntup_list, sort(ntup_list))]
help, ntup_list
help, uniq_list

; loop over ntuples and all possible appendices
new_nlist= 0L & new_ntup_list= -1
for jj=0L,nlist-1 do begin
    test_tuple= [ntup_list[*,jj],0]
    index= where(uniq_list GT ntup_list[ntup-2,jj],ngreater)
    for kk=0L,ngreater-1 do begin
        test_tuple[ntup-1]= uniq_list[index[kk]]

; check this tuple
        accept= check_tuple(viv,iv,icov,test_tuple,nsigma)

; append this tuple
        if accept then begin
            if new_nlist EQ 0 then new_ntup_list= test_tuple else $
              new_ntup_list= [[new_ntup_list],[test_tuple]]
            new_nlist= new_nlist+1
        endif
    endfor
endfor

failed= 0B
if n_elements(new_ntup_list) EQ 1 AND new_ntup_list[0] EQ -1 then begin
    new_ntup_list=ntup_list
    failed= 1B
endif
return, new_ntup_list
end

