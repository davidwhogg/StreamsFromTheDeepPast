pro oplotbarx,xxx,_extra=extra

n = n_elements(xxx)

yrange = !y.crange

if !y.type then yrange = 10.^yrange

for i=0,n-1 do $

oplot,[xxx[i],xxx[i]],yrange,_extra=extra



return

end
