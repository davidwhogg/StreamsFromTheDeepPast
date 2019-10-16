; Program to convert (ra,dec,pmra,pmdec) to
; (l,b,pml,pmb), where pmra and pml include
; cos dec, cos b factors...
;     - angles in degrees
;     - proper motions in mas/year

; =======================================================================
PRO gco,x,year

; ROTATE X,Y,Z (u,v,w) TO BE ALIGNED WITH GALACTIC COORDS
; Numbers from Allen's Astrophysical quantities, p 575
; RA of ascending node of Galactic Plane is 282.86
; angle between Celestial Equator and Galactic Plane is 62.87	
; Galactic longitude of ascending node of Galactic Plane is 32.93

      a3=(192.85948123d0)
      decpol = (27.4d0 - 5.44d-3 * (year-1950.0d0))
      
; rotate about Celestial pole to line NGP up with z-y plane
      t=(a3+6.d0*15.d0)
      rotate,x,t,3
;     ------
; rotate about x-axis to make NGP new z-axis
      t=-(decpol-90.d0)	
      rotate,x,t,1
;     ------
; rotate about NGP to get Galactic Center at l=0
      l3=(122.9319186d0)
      t=-(l3-90.d0)	
      rotate,x,t,3
;     ------
      RETURN
      END


; =======================================================================

PRO rotate,x,t,idim
      
R=DBLARR(3,3)

; ROTATION ANGLE
      t = t*!DTOR

IF(idim EQ 1)THEN BEGIN
R(0,0) = 1.0
R(1,1) = COS(t)
R(2,1) = SIN(t)
R(1,2) = -SIN(t)
R(2,2) = COS(t)
ENDIF
      
IF(idim EQ 2)THEN BEGIN
R(0,0) = COS(t)
R(0,2) = -SIN(t)
R(1,1) = 1.0
R(2,0) = SIN(t)
R(2,2) = COS(t)
ENDIF

IF(idim EQ 3)THEN BEGIN
R(0,0) = COS(t)
R(1,1) = COS(t)
R(2,2) = 1.0
R(1,0) = SIN(t)
R(0,1) = -SIN(t)
ENDIF

        xyzp = R##x
	x=xyzp

RETURN
END

PRO radec_to_lb, ra,dec,pmra,pmdec,l,b,pml,pmb,year

; DEFINE THETA AND PHI IN TERMS OF RA (a) AND DEC (d) AND CONVERT TO RADIANS
      th=(90.0d0-dec)*!DTOR
      ph=ra*!DTOR
	n=N_ELEMENTS(ra)
;
	x=DBLARR(n,3)
      x[*,0] = (SIN(th)*COS(ph))
      x[*,1] = (SIN(th)*SIN(ph))
      x[*,2] = (COS(th))

; Use Blake's program to find (x,y,z) in Galactic coordinates
      gco,x,year
;     ---         
; Correct coordinate system
      pmth=-pmdec
; Rotate proper motions to (pml,pmb)
	vx=DBLARR(n,3)
      vx[*,0]=pmth*COS(th)*COS(ph) - pmra*SIN(ph)
      vx[*,1]=pmth*COS(th)*SIN(ph) + pmra*COS(ph)
      vx[*,2]=-pmth*SIN(th)
; Use Blake's program to find (x,y,z) in Galactic coordinates
      gco,vx,year
;     ---         
      th=ACOS(x[*,2])
      ph=ATAN(x[*,1],x[*,0])
      l=(ph/!DTOR)

	mask=WHERE(l LT 0.0d0,COUNT)
      IF(COUNT GT 0)THEN l(mask)=360.0d0+l(mask)
      b=(90.d0-th/!DTOR)

      pmb=-vx[*,2]/SIN(th)
      pml=-(vx[*,0]*SIN(ph)-vx[*,1]*COS(ph))

      pmb=-pmb

RETURN

END











