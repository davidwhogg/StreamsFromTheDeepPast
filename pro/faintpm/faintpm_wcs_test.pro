pro faintpm_wcs_test
astr= hogg_make_astr(265.0,15.0,0.1,0.1,pixscale=0.1/3600.0,orientation=15.0)
splog, astr.cd
astr.cd= astr.cd+[[0.0,1d-6],[0.0,0.0]]
splog, astr.cd
xx= astr.crpix[0]-1.5d0
yy= astr.crpix[1]-1.5d0
dx= 1d0
dy= 1d0
xy2ad, xx,yy,astr,aa,dd
xy2ad, xx+dx,yy,astr,aax,ddx
xy2ad, xx,yy+dy,astr,aay,ddy
dadx= (aax-aa)/dx
dddx= (ddx-dd)/dx
dady= (aay-aa)/dy
dddy= (ddy-dd)/dy
splog, astr.cd-([[dadx,dady],[dddx,dddy]])
splog, astr.cd-transpose([[dadx,dady],[dddx,dddy]])
return
end
