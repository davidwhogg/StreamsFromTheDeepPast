pro outweights,savefile,outfile

restore,'hip.sav'
restore,savefile

dummy={emstruct, hipid:hip[0].hip, weight:transpose(weight[0,*])}
outstruct=replicate(dummy,(size(weight))[1])
outstruct.hipid=hip[index].hip
for i=0, ngauss-1 do outstruct[*].weight[i]=weight[*,i]
mwrfits,outstruct,outfile,/create

end
