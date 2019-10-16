;+
; NAME:
;   hogg_plot_stars
; PURPOSE:
;   make plot for Adi
;-
pro hogg_plot_stars

; set constants
prefix= 'hogg_plot_stars'
xsize= 7.5 & ysize= xsize
adi= mrdfits('/global/data/scr/adi/sim_halo.fits',1)
drange=[0,300]
massrange=[0,1.5e11]
timerange=[0,14]

; setup plot
set_plot, "PS"
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
hogg_plot_defaults

; plot satellite mass vs radius
;plot, adi.d,adi.maxhalo,psym=1,symsize=0.01,thick=2, $
;  xrange=drange,yrange=massrange

; group by satellite mass
maxhalolist= adi[sort(adi.maxhalo)].maxhalo
maxhalolist= reverse(maxhalolist[uniq(maxhalolist)])
nmaxhalo= n_elements(maxhalolist)
help, adi,/struc

for ii=0,19 do begin
    title= 'progenitor halo mass '+strtrim(string(maxhalolist[ii]),2)
    inthis= where(adi.maxhalo EQ maxhalolist[ii])
    plot, adi[inthis].d,adi[inthis].acctime, $
      psym=1,symsize=0.1,thick=2, $
      xrange=drange,yrange=timerange, $
      title=title
endfor

device, /close
end
