; hogg_scatterplot,file.ra,file.dec,grid=grid
; fit= az_fit(grid)
; atv, grid-fit
function az_fit,grid,xorder=xorder,yorder=yorder
if (not keyword_set(xorder)) then xorder=2
if (not keyword_set(yorder)) then yorder=2
foo= size(grid,/dimens)
nx= foo[0]
ny= foo[1]
err= sqrt(grid)
good= where(grid GT 5)
xvec= fltarr(nx)+1.0
for ii=0,xorder do begin
    ;xvec= fltarr(ny)+1.0
    yvec=fltarr(ny) +1.0 ; adi added 
    for jj=0,yorder do begin
        if keyword_set(aa) then aa=[[aa],[(xvec#yvec)[good]]] $
        else aa=[[(xvec#yvec)[good]]]
        yvec= yvec*findgen(ny)
    endfor
    xvec= xvec*findgen(nx)
endfor
aa= transpose(aa)
invvar= 1.0/err[good]
hogg_iter_linfit, aa,grid[good],invvar,pars
fit= fltarr(nx,ny)
fit[good]= aa##pars
return, fit
end
