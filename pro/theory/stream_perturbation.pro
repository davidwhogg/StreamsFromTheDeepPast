;+
; USAGE:
;  echo stream_perturbation | idl
; BUGS:
;  - No proper header.
;  - Very brittle.
;  - No license.
;-
pro stream_perturbation
seed= -1L
nstar= 3001L
width= 1D-4
length= 3D1
star= dblarr(3,nstar)
star[0,*]= length*(-1.0+2.0*dindgen(nstar)/double(nstar))
star[1,*]= width*randomn(seed,nstar)
star[2,*]= width*randomn(seed,nstar)
cname= ['x','y','z']
crange= [[-7.0,7.0],[-1.0,3.0],[-2.0,2.0]]
cxlist= [0,2,0]
cylist= [1,1,2]
position= [[0.05,0.2,0.75,0.4],[0.75,0.2,0.95,0.4],[0.05,0.0,0.75,0.2]]
set_plot, 'ps'
loadct, 0
xsize= 3.0
ysize= xsize*2.0
poff= 0.05*xsize/ysize
grey= 127
device, xsize=xsize,ysize=ysize,/inches,/encap,/color, $
        filename='grid.eps'
for angle=3D0*!DPI/8D0,-1D-4,-!DPI/8D0 do begin
   anglestr= string(angle,format='(F5.3)')
   vx= sin(angle)
   vz= cos(angle)
   jindx= lonarr(nstar)
   dist2= dblarr(nstar)
   mass= [0D0,1D0,0D0]#replicate(1D0,nstar) $
         +[vx,0D0,vz]#(length*(-1.0+2.0*dindgen(nstar)/double(nstar)))
   for ii=0L,nstar-1L do begin
      thisdist2= total((star[*,ii]#replicate(1D0,nstar)-mass)^2,1)
      dist2[ii]= min(thisdist2,jj)
      jindx[ii]= jj
   endfor
   starv= (mass[*,jindx]-star)/(replicate(1D0,3)#dist2)
   for cc=0,2 do begin
      cx= cxlist[cc]
      cy= cylist[cc]
      !P.title= ''
      !P.charsize= 0.5
      !X.title= cname[cx]
      !X.style= 1
      !X.charsize= 0.5
      !Y.title= cname[cy]
      !Y.style= !X.style+16
      !Y.charsize= !X.charsize
      if (cc EQ 0) then begin
         !P.title= 'arctan(vz/vx) = '+anglestr+' rad'
         !X.charsize= 0.001
      endif else begin
         if (cc EQ 1) then begin
            !Y.charsize= 0.001
         endif
      endelse
      !P.position= position[*,cc]*([1,0,1,0]+[0,1,0,1]*(xsize/ysize)) $
                   +[0,1,0,1]*poff
      splog, !P.position
      plot, [0],/isotropic,/nodata,/noerase,color=grey, $
            xrange=crange[*,cx],yrange=crange[*,cy]
      for time=0.0,2.01,0.25 do begin
         pos= star+time*starv
         oplot, pos[cx,*],pos[cy,*],psym=3
      endfor
   endfor
   poff= poff+0.25
endfor
device,/close
xsize= (5.0/6.0)*ysize
device, xsize=xsize,ysize=ysize,/inches,/encap, $
        filename='views.eps'
for xx=0,4 do for yy=0,5 do begin
   angle= 2.0*!DPI*randomu(seed)
   time= 0.25+1.0*randomu(seed)
   vx= sin(angle)
   vz= cos(angle)
   jindx= lonarr(nstar)
   dist2= dblarr(nstar)
   mass= [0D0,1D0,0D0]#replicate(1D0,nstar) $
         +[vx,0D0,vz]#(length*(-1.0+2.0*dindgen(nstar)/double(nstar)))
   for ii=0L,nstar-1L do begin
      thisdist2= total((star[*,ii]#replicate(1D0,nstar)-mass)^2,1)
      dist2[ii]= min(thisdist2,jj)
      jindx[ii]= jj
   endfor
   starv= (mass[*,jindx]-star)/(replicate(1D0,3)#dist2)
   repeat begin
      xhat= randomn(seed,3)
      xhat= xhat/sqrt(total(xhat^2))
      yhat= randomn(seed,3)
      yhat= yhat-xhat*total(yhat*xhat)
      yhat= yhat/sqrt(total(yhat^2))
      vpos= star##transpose([[xhat],[yhat]])
   endrep until (max(vpos) gt 3.0)
   !P.title= ''
   !X.charsize= 0.001
   !Y.charsize= !X.charsize
   !P.position= [xx/5.0,yy/6.0,(xx+1.0)/5.0,(yy+1.0)/6.0]
   vpos= star##transpose([[xhat],[yhat]])
   plot, vpos[0,*],vpos[1,*],psym=3,/isotropic,/noerase,color=grey, $
         xrange=[-3,3],yrange=[-3,3]
   pos= star+time*starv
   vpos= pos##transpose([[xhat],[yhat]])
   oplot, vpos[0,*],vpos[1,*],psym=3
endfor
device,/close
return
end
