C Program to bin particles in velocity bins....
      INTEGER ix,iy,iz,ngrid
      PARAMETER(ngrid=50)
      INTEGER nran(ngrid,ngrid,ngrid),ntot(ngrid,ngrid)
      INTEGER i,nhip,nstars
      REAL*8 ra,dec,p,pmra,pmdec,d,
     &       radpdeg,kmppc,maspdeg,secpyr,conv,
     &       l,b,pml,pmb,year,vl,vb
      PARAMETER(kmppc=3.086d13,maspdeg=60.0d0*60.0d0*1.0d3)
      PARAMETER(secpyr=60.0d0*60.0d0*24.0d0*365.0d0)
      
      year=1991.25d0
      radpdeg=4.0d0*DATAN(1.0d0)/180.0d0
      conv=kmppc*radpdeg/maspdeg/secpyr
      dlim=100.0d0
C
      DO 5 ix=1,ngrid
         DO 6 iy=1,ngrid
            ntot(ix,iy)=0
            DO 7 iz=1,ngrid
               nran(ix,iy,iz)=0
 7          CONTINUE
 6       CONTINUE
 5    CONTINUE

      ix=0
      nstars=20847

      OPEN(UNIT=10,FILE='for_kvj.txt',STATUS='old')
      OPEN(UNIT=11,FILE='kvj.test',STATUS='UNKNOWN')
      READ(10,*)
C
      DO 10 i=1,nstars
         IF(MOD(i,1000).EQ.0) WRITE(6,*) i
C Read in phase-space distribution
         READ(10,*)nhip,ra,dec,p,pmra,pmdec
         d=1000.d0/p
         IF(d.LT.dlim)THEN
C
C Translate to Galactic coordinates
            CALL trans(ra,dec,pmra,pmdec,l,b,pml,pmb,year)
C                -----
C         WRITE(6,100)ra,dec,pmra,pmdec
c         WRITE(6,100)SQRT(pmra*pmra+pmdec*pmdec)

C Define test quantities
            l=l*radpdeg
            b=b*radpdeg
            vl=d*pml*conv
            vb=d*pmb*conv
C
            ix=ix+1
C            WRITE(11,100)l,b,vl,vb,ra,dec
C Loop through velocity grid and add star to cell if its vlos
C falls within a grid point
            CALL vcell(l,b,vl,vb,nran)
C                -----
         ENDIF
 10   CONTINUE
      CLOSE(10)
      CLOSE(11)

      WRITE(6,*)ix
C Write our x-y plane as function of z
      OPEN(UNIT=10,FILE='vbox.dat',STATUS='UNKNOWN')
      WRITE(UNIT=10,nran)
C      DO 20 iz=1,ngrid
C         DO 30 ix=1,ngrid
C            DO 40 iy=1,ngrid
C               WRITE(10,99)nran(ix,iy,iz)
C               ntot(ix,iy)=ntot(ix,iy)+nran(ix,iy,iz)
C 40         CONTINUE
C 30      CONTINUE
C 20   CONTINUE

C      DO 50 ix=1,ngrid
C         DO 60 iy=1,ngrid
C            WRITE(10,99)ntot(ix,iy)
C 60      CONTINUE
C 50   CONTINUE

      CLOSE(10)
 99   FORMAT(2I8)
 100  FORMAT(6(1pe10.2))
      END

C ****************************************************************
C
      SUBROUTINE vcell(ra,dec,vra,vdec,ncell)
C
C ****************************************************************
      INTEGER ix,iy,iz,ngrid
      PARAMETER(ngrid=50)
      INTEGER ncell(ngrid,ngrid,ngrid)
      REAL*8 ra,dec,vra,vdec,vvra,vvdec,
     &       x,y,z,vx,vy,vz,dv,vesc,r,
     &       xra,yra,zra,xdec,ydec,zdec
      PARAMETER(vesc=400.0d0)


C Directional vector
      x=DCOS(dec)*DCOS(ra)
      y=DCOS(dec)*DSIN(ra)
      z=DSIN(dec)
C in ra direction
      r=SQRT(x*x+y*y)
      xra=-y/r
      yra=x/r
      zra=0.0d0
C in dec direction
      xdec=y*zra-z*yra
      ydec=z*xra-x*zra
      zdec=x*yra-y*xra
C
C      WRITE(6,99)SQRT(x*x+y*y+z*z)
C      WRITE(6,99)SQRT(xra*xra+yra*yra+zra*zra)
C      WRITE(6,99)SQRT(xdec*xdec+ydec*ydec+zdec*zdec)
C      WRITE(6,99)(x*xra+y*yra+z*zra)
C      WRITE(6,99)(x*xdec+y*ydec+z*zdec)
C      WRITE(6,99)(xdec*xra+ydec*yra+zdec*zra)
C Cell size
      dv=2.0d0*vesc/DBLE(REAL(ngrid))
C
      DO 10 ix=1,ngrid
         vx=-vesc+dv*(DBLE(REAL(ix))-0.5d0)
         DO 20 iy=1,ngrid
            vy=-vesc+dv*(DBLE(REAL(iy))-0.5d0)
            DO 30 iz=1,ngrid
               vz=-vesc+dv*(DBLE(REAL(iz))-0.5d0)
C
               vvra=vx*xra+vy*yra+vz*zra
               vvdec=vx*xdec+vy*ydec+vz*zdec
C
               IF((ABS(vvra-vra).LT.dv).AND.(ABS(vvdec-vdec).LT.dv))
     &              ncell(ix,iy,iz)=ncell(ix,iy,iz)+1
C
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
 99   FORMAT(5(1pe10.2))
      RETURN
      END

