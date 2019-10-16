PRO hesstest,nstart,fname
icall=1
plim=0.5
nbin=40
j=0

!path=EXPAND_PATH('+/home/users/astro1/kvj/hipparcos/idl')+':'+!path
;!path=EXPAND_PATH('+/home/kvj/idl')+':'+!path
; Tests for differences in Hess diagrams from groups identified in
; Hipparcos

!P.MULTI=[0,4,5]

IF(icall EQ 0)THEN BEGIN
; Read in Hipparcos data
; note: NaN in bv column!
readcol,'for_kvj.txt',ihipp,ra,dec,p,pmra,pmdec,v,bv,SKIPLINE=1
d=1000/p
v=v-5*ALOG10(d/10)

; Read in probabilities of being members of particular groups
; readcol,'weights.halo.g5.txt',nhipp,w0,w1,w2,w3,w4,SKIPLINE=2
readcol,'weights.halo.g14.txt',nhipp,w0,w1,w2,w3,w4,w5,w6,w7,w8,w9, $
         w10,w11,w12,w13,SKIPLINE=2
nw=14
n=N_ELEMENTS(w0)
w=FLTARR(nw,n)
w[0,*]=w0
w[1,*]=w1
w[2,*]=w2
w[3,*]=w3
w[4,*]=w4
w[5,*]=w5
w[6,*]=w6
w[7,*]=w7
w[8,*]=w8
w[9,*]=w9
w[10,*]=w10
w[11,*]=w11
w[12,*]=w12
w[13,*]=w13

; reduce v,bv vectors to points in weights file
n=max(ihipp+1)
point=FLTARR(LONG(n))
point(ihipp-1)=v
v=point(nhipp-1)
point(ihipp-1)=bv
bv=point(nhipp-1)

SAVE,bv,v,w,FILENAME='hesstest.dat'
ENDIF ELSE BEGIN
;RESTORE,FILENAME='hesstest.dat'
RESTORE,FILENAME=fname
ENDELSE

; Make hess diagram of everything

ntot=N_ELEMENTS(bv)
mask=INDGEN(ntot)

IF(j GT 0)THEN BEGIN
j=j-1
p=w[j,*]
mask=WHERE(p GT plim,count)
ntot=count
ENDIF

makehess,bv,v,mask,pback,nbin

;makeplots,pback,pback,1

makebars,nbin

FOR j=0,4 DO BEGIN
k=j+nstart
; Make hess diagrams at chosen probability level 
	p=w[k,*]
	mask=WHERE(p GT plim,count)
	IF(count GT 0)THEN makehess,bv,v,mask,ptot,nbin

; Plot the diagrams
;	makeplots,ptot,pback,j

; Find khisq assuming poisson approximation to binomial distribution
; for each grid point
	khisq,ptot,pback,j,k2,count,nbin,cback

	print,k,count,cback,k2,SQRT(k2)
ENDFOR

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO khisq,ptot,pback,j,k2,count,nbin,cback
COMMON bins,xmin,xmax,ymin,ymax

n=10
khi2=FLTARR(n)
alpha=(FINDGEN(n)+1.0)/FIX(n)
mask=WHERE(ptot GT 0)
ka=FLTARR(nbin,nbin)

FOR i=0,n-1 DO BEGIN
	ka(mask)=(ptot(mask)-alpha(i)*pback(mask))/SQRT(ptot(mask))
	khi2[i]=TOTAL(ka)
ENDFOR
khi2=khi2*count/FIX(nbin*nbin)

mask=WHERE(pback GT 0,cback)
khi=FLTARR(nbin,nbin)
khi(mask)=(ptot(mask)-pback(mask))/SQRT(pback(mask))
k2=TOTAL(khi*khi)
khi=khi*SQRT(count)
k2=k2*count

;print,min(khi),max(khi)

khiplots,ptot,khi,alpha,khi2,j,count,pback,nbin

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO khiplots,ptot,khi,alpha,khi2,j,count,pback,nbin
COMMON bins,xmin,xmax,ymin,ymax

x=[xmin-1.]
y=[ymin-1.]
FOR i=1,3 DO BEGIN
;SPECIFY LOCATION OF PLOT
        y0=(1.5-0.27*(j+1))*!D.X_SIZE/!D.Y_SIZE-0.05
        y1=y0+0.25*!D.X_SIZE/!D.Y_SIZE
        x0=0.3*(i-1)+0.05
        x1=x0+0.25
;SCALE IMAGE TO FIT DISPLAY
	PX = [x0,x1]* !D.X_VSIZE
	PY = [y0,y1]* !D.Y_VSIZE
;Desired size of image in pixels.
	SX = PX[1] - PX[0] + 1            
	SY = PY[1] - PY[0] + 1

	IF(i LT 4)THEN BEGIN
		PLOT,x,y,XRANGE=[xmin,xmax],YRANGE=[ymax,ymin], $
			 POS=[x0,y0,x1,y1],XSTYLE=1,YSTYLE=1, CHARSIZE=2.

	
		h=ptot*count
		hmin=0.
		hmax=10.
		IF(i EQ 1)THEN BEGIN
			h=ptot
			hmin=0
			hmax=40./nbin/nbin
		ENDIF ELSE BEGIN		
			IF(i EQ 2)THEN BEGIN
				h=khi
				hmin=-3.
				hmax=6.
			ENDIF ELSE BEGIN
			h=pback
			hmin=0
			hmax=40./nbin/nbin
;				h=ptot*count
;				hmin=0.
;				hmax=5.
			ENDELSE
		ENDELSE
;		print,hmin,hmax

	        map=CONGRID(h, SX, SY)
        	map=BYTSCL(map, MIN=hmin, MAX=hmax, TOP = !D.TABLE_SIZE)
	        TV, map, PX[0], PY[0]
;	        TVSCL, CONGRID(h, SX, SY), PX[0], PY[0]

        ENDIF ELSE BEGIN
		PLOT,alpha,khi2,XRANGE=[0.,1.1],YRANGE=[0.,9.], $
			 POS=[x0,y0,x1,y1]
	ENDELSE

ENDFOR

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO makebars,nbin

LOADCT, 5

FOR i=1,3 DO BEGIN

IF(i EQ 2)THEN BEGIN
	x0=0.35
	xmin=-3.
	xmax=6.	
	label='khi'
ENDIF ELSE BEGIN
	IF(i EQ 1)THEN BEGIN
		x0=0.05
		xmin=0.
		xmax=40./nbin/nbin
		label='probability'
	ENDIF ELSE BEGIN
		x0=.65
		xmin=0.
		xmax=5.
		label='number'
	ENDELSE
ENDELSE

;SPECIFY LOCATION
y0=0.95
y1=y0+.015
x1=x0+0.25
PX = [x0,x1]* !D.X_VSIZE
PY = [y0,y1]* !D.Y_VSIZE
SX = PX[1] - PX[0] + 1
SY = PY[1] - PY[0] + 1

nbox=10
dx=(xmax-xmin)/FIX(nbox)
xave=xmin+dx*FINDGEN(nbox)
bar=FLTARR(nbox,2)
bar(*,0)=xave
bar(*,1)=xave
cex=xave+dx/2.
cwhy=[0,1]
CONTOUR, bar,cex,cwhy, XSTYLE=1, YSTYLE=4, POS=[x0,y0,x1,y1],/NODATA, $
XRANGE=[xmin,xmax], CHARSIZE=2.
scale=CONGRID(bar, SX, SY)
scale=BYTSCL(scale, MIN=xmin, MAX=xmax, TOP = !D.TABLE_SIZE)
TV, scale, PX[0], PY[0]
CONTOUR, bar,cex, cwhy, /NODATA, YSTYLE=4,/NOERASE, $
POS=[x0,y0,x1,y1], XTITLE=label, XSTYLE=1, $
XRANGE=[xmin,xmax], CHARSIZE=2.

ENDFOR
RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO makehess,bv,v,mask,h,nbin
COMMON bins,xmin,xmax,ymin,ymax

x=bv(mask)
y=v(mask)
n=N_ELEMENTS(x)
eps=.001

xmin=-.5
xmax=2.5
xbin=(xmax-xmin)/FLOAT(nbin)
ymin=-6.
ymax=16.
ybin=(ymax-ymin)/FLOAT(nbin)

;print,xmin,xmax,ymin,ymax

ix=FIX((x-xmin)/xbin)
iy=FIX((y-ymin)/ybin)

h=FLTARR(nbin,nbin)
FOR i=0,n-1 DO BEGIN
	h(ix(i),nbin-iy(i)-1)=h(ix(i),nbin-iy(i)-1)+1.0
ENDFOR

h=h/FLOAT(n)

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO makeplots,ptot,pback,j
COMMON bins,xmin,xmax,ymin,ymax
x=[xmin-1.]
y=[ymin-1.]


pratio=FLTARR(nbin,nbin)
mask=WHERE(pback GT 0 )
pratio(mask)=ptot(mask)/pback(mask)

FOR i=1,3 DO BEGIN
;SPECIFY LOCATION OF PLOT
        y0=(1.5-0.3*(j+1))*!D.X_SIZE/!D.Y_SIZE
        y1=y0+0.25*!D.X_SIZE/!D.Y_SIZE
        x0=0.3*(i-1)+0.05
        x1=x0+0.25

PLOT,x,y,XRANGE=[xmin,xmax],YRANGE=[ymax,ymin],POS=[x0,y0,x1,y1]
;SCALE IMAGE TO FIT DISPLAY
        PX = [x0,x1]* !D.X_VSIZE
        PY = [y0,y1]* !D.Y_VSIZE

;Desired size of image in pixels.
        SX = PX[1] - PX[0] + 1            
        SY = PY[1] - PY[0] + 1

	h=pback
	IF(i EQ 2)THEN h=ptot
	IF(i EQ 3)THEN h=pratio

        TVSCL, CONGRID(h, SX, SY), PX[0], PY[0]

ENDFOR

RETURN
END

