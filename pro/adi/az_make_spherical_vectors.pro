function az_make_spherical_vectors, nstars=nstars
;+
;Adi Zolotov
;April 2008
;-

seed1=7
seed2=18
seed3=-2
nstars=nstars
vect1=randomn(seed1,nstars)
vect2=randomn(seed2,nstars)
vect3=randomn(seed3,nstars)
unitvec=dblarr(nstars,3)
unitvec[*,0]=vect1[*]
unitvec[*,1]=vect2[*]
unitvec[*,2]=vect3[*]
norm=dblarr(1,nstars)
norm[*]=sqrt((unitvec[*,0]*unitvec[*,0])+(unitvec[*,1]*unitvec[*,1])+$
  (unitvec[*,2]*unitvec[*,2]))
unit=dblarr(nstars,3)
unit[*,0]=unitvec[*,0]/norm[*]
unit[*,1]=unitvec[*,1]/norm[*]
unit[*,2]=unitvec[*,2]/norm[*]


return,unit

end
