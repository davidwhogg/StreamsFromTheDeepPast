; hogg_scatterplot,file.ra,file.dec,ynpix=40,xnpix=40,grid=grid
function az_fit,grid
foo= size(grid,/dimens)
nx= foo[0]
ny= foo[1]
err= sqrt(grid)
good= where(grid GT 1)
xvec= fltarr(nx)+1.0
for xorder=0,2 do begin
    yvec= fltarr(ny)+1.0
    for yorder=0,2 do begin
        if keyword_set(aa) then aa=[[aa],[(xvec#yvec)[good]]] $
        else aa=[[(xvec#yvec)[good]]]
        yvec= yvec*findgen(ny)
    endfor
    xvec= xvec*findgen(nx)
endfor
help, aa
x= findgen(nx)#replicate(1.0,ny)
y= replicate(1.0,nx)#findgen(ny)
aa= transpose([[one[good]],[x[good]],[y[good]]])
invvar= 1.0/err[good]
hogg_iter_linfit, aa,grid[good],invvar,pars
fit= fltarr(nx,ny)
fit[good]= aa##pars
return, fit
end
