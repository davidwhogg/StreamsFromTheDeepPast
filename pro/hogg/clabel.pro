; Select on the basis of a cut on probability
function clabel, weight, pgauss

ndata=(size(weight,/dimensions))[0]
ngauss=(size(weight,/dimensions))[1]
help,ndata,ngauss

clab=fltarr(4,ndata)+1.
clab[3,*]=0.

cols=fltarr(4,10)
cols[*,0]=[1.,0.,0.,0.7]
cols[*,1]=[0.,1.,0.,0.7]
cols[*,2]=[0.,0.,1.,0.7]
cols[*,3]=[1.,1.,0.,0.7]
cols[*,4]=[0.,1.,1.,0.7]
cols[*,5]=[1.,0.,1.,0.7]
cols[*,6]=[1.,0.,0.,0.7]
cols[*,7]=[1.,0.,0.,0.7]
cols[*,8]=[1.,0.,0.,0.7]
cols[*,9]=[1.,0.,0.,0.7]

for i=0l, ngauss-1l do begin 
    indx=where(weight[*,i] gt pgauss,count)
    if (count gt 0) then begin
        clab[0,indx]=cols[0,i mod 10]
        clab[1,indx]=cols[1,i mod 10]
        clab[2,indx]=cols[2,i mod 10]
        clab[3,indx]=cols[3,i mod 10]
    endif
end 

return, clab

end
