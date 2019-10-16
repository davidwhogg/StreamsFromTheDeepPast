pro testdecon,origdata,decondata,filename

indata=dblarr(30l,30l,31l)
openr,11,filename
readf,11,indata
close,11
origdata=dblarr(100l,100l,100l)
origdata[0:29,0:29,0:29]=indata[0:29,0:29,0:29]
decondata=origdata

green=dblarr(100l,100l,100l)
for i = 0, 14 do begin 
for j = 0, 14 do begin 
for k = 0, 14 do begin 
green[i,j,k]=1.d
if(i ne 0 or j ne 0 or k ne 0) then green[i,j,k]=1.d/(i^2+j^2+k^2)
end
end
end

deconvolve,decondata,green

end
