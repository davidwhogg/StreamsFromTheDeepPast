C Program to convert (ra,dec,pmra,pmdec) to
C (l,b,pml,pmb), where pmra and pml include
C cos dec, cos b factors...
C     - angles in degrees
C     - proper motions in mas/year

      SUBROUTINE trans(ra,dec,pmra,pmdec,l,b,pml,pmb,year)
      REAL*8 ra,dec,pmra,pmdec,pml,pmb,year,x(3),v(3),radeg,th,ph,
     &       l,b

      radpdeg=4.0d0*DATAN(1.0d0)/180.d0

C DEFINE THETA AND PHI IN TERMS OF RA (a) AND DEC (d) AND CONVERT TO RADIANS
      th=(90.0d0-dec)*radpdeg
      ph=ra*radpdeg

C Rotate postions to (l,b)
      x(1) = (DSIN(th)*DCOS(ph))
      x(2) = (DSIN(th)*DSIN(ph))
      x(3) = (DCOS(th))
C Use Blake's program to find (x,y,z) in Galactic coordinates
      CALL gco(x,year)
C          ---         
C Correct coordinate system
      pmth=-pmdec
C Rotate proper motions to (pml,pmb)
      v(1)=pmth*COS(th)*COS(ph) - pmra*SIN(ph)
      v(2)=pmth*COS(th)*SIN(ph) + pmra*COS(ph)
      v(3)=-pmth*SIN(th)
C Use Blake's program to find (x,y,z) in Galactic coordinates
      CALL gco(v,year)
C          ---         
      th=DACOS(x(3))
      ph=DATAN2(x(2),x(1))
      l=(ph/radpdeg)
      IF(l.LT.0.d0)l=(360.d0+l)
      b=(90.d0-th/radpdeg)
C
      pmb=-v(3)/SIN(th)
      pml=-(v(1)*SIN(ph)-v(2)*COS(ph))
C
      pmb=-pmb

      END

C =======================================================================
      SUBROUTINE gco(x,year)

C USE MX TO ROTATE X,Y,Z (u,v,w) TO BE ALIGNED WITH GALACTIC COORDS
C Numbers from Allen's Astrophysical quantities, p 575
C RA of ascending node of Galactic Plane is 282.86
C angle between Celestial Equator and Galactic Plane is 62.87	
C Galactic longitude of ascending node of Galactic Plane is 32.93
      REAL*8 a3,decpol,year,tx,tz,l3,x(3)

      a3=(192.85948123d0)
      decpol = (27.4d0 - 5.44d-3 * (year-1950.0d0))
      

C rotate about Celestial pole to line NGP up with z-y plane
      tz=(a3+6.d0*15.d0)
      CALL rotate(x,tz,3)
C          ------
C rotate about x-axis to make NGP new z-axis
      tx=-(decpol-90.d0)	
      CALL rotate(x,tx,1)
C          ------
C rotate about NGP to get Galactic Center at l=0
      l3=(122.9319186d0)
      tz=-(l3-90.d0)	
      CALL rotate(x,tz,3)
C          ------

      RETURN
      END


C =======================================================================

      SUBROUTINE rotate(x,tx,idim)
      INTEGER i,j,idim
      REAL*8 x(3),xarray(3,3),tx,xr(3),radpdeg
      
      radpdeg=4.0d0*DATAN(1.0d0)/180.0d0
C ROTATION ANGLE
      tx = tx*radpdeg

C MAKE ROTATION MATRIX
      DO 10 i=1,3
         xr(i)=x(i)
         DO 20 j=1,3
            xarray(i,j)=0.d0
 20      CONTINUE
 10   CONTINUE

      IF(idim.EQ.3)THEN
         xarray(1,1) = DCOS(tx)
         xarray(1,2) = DSIN(tx)
         xarray(2,1) = -DSIN(tx)
         xarray(2,2) = DCOS(tx)
         xarray(3,3) = 1.0d0
      ELSE IF(idim.EQ.2)THEN
         xarray(1,1) = DCOS(tx)
         xarray(1,3) = DSIN(tx)
         xarray(2,2) = 1.0d0
         xarray(3,1) = -DSIN(tx)
         xarray(3,3) = DCOS(tx)
      ELSE
         xarray(1,1) = 1.0d0
         xarray(2,2) = DCOS(tx)
         xarray(2,3) = DSIN(tx)
         xarray(3,2) = -DSIN(tx)
         xarray(3,3) = DCOS(tx)
      ENDIF

      DO 30 i=1,3
         x(i)=0.0d0
         DO 40 j=1,3
            x(i)=x(i)+xarray(i,j)*xr(j)
 40      CONTINUE
 30   CONTINUE
C
      RETURN
      END











