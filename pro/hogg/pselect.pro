; Select on the basis of a cut on probability
function pselect, weight, pgauss

ndata=(size(weight))[1]
ngauss=(size(weight))[2]
indx=lonarr(ndata,ngauss)-1l
for i=0l, ngauss-1l do begin 
    tmpindx=where(weight[*,i] gt pgauss,count)
    indx[tmpindx,i]=i
end 

return, indx

end
