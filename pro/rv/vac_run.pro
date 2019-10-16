;+
;   Run the AL sample while I'm on vacation
;
;   Bovy 2008-08-11
;-

PRO VAC_RUN, processor= processor, stop= stop

;;Grid
grid= dblarr(2,2,9)

grid[0,*,0]= [10,0.01D]
grid[0,*,1]= [20,0.01D]
grid[0,*,2]= [30,0.01D]
grid[0,*,3]= [40,0.01D]
grid[0,*,4]= [50,0.01D]
grid[0,*,5]= [60,0.01D]
grid[0,*,6]= [10,0.3D]
grid[0,*,7]= [20,0.3D]
grid[0,*,8]= [30,0.3D]


grid[1,*,0]= [60,1D]
grid[1,*,1]= [50,1D]
grid[1,*,2]= [40,1D]
grid[1,*,3]= [30,1D]
grid[1,*,4]= [20,1D]
grid[1,*,5]= [10,1D]
grid[1,*,6]= [40,0.3D]
grid[1,*,7]= [50,0.3D]
grid[1,*,8]= [60,0.3D]

IF keyword_set(stop) THEN stop

FOR ii= 0, 8 DO BEGIN
    IF (processor EQ 0 AND (ii EQ 0 OR ii EQ 1 OR ii EQ 2 OR ii EQ 3)) THEN CONTINUE
    fitv, colors=[0.0D,0.4D,0.6D], sample=6, /bestV, Vconstraint=grid[processor,1,ii],ngauss=grid[processor,0,ii], /initial, Vjacksub=.01D, /vanilla, /log, /giants, /all
ENDFOR




END
