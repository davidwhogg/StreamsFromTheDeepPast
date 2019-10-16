-- ----------------------------------------------------------------------------
-- NAME:
--  hogg_proper_motion_gc.sql
-- PURPOSE:
--  select stars near a GC with proper motions from the CAS
-- NOTES:
--  * The "fGet" function apparently is much faster than explicit RA,Dec cuts.
-- REVISION HISTORY:
--  2006-07-08  started - Hogg (NYU)
-- ----------------------------------------------------------------------------
select (f.psfmag_u-f.extinction_u) as u,
       (f.psfmag_g-f.extinction_g) as g,
       (f.psfmag_r-f.extinction_r) as r,
       (f.psfmag_i-f.extinction_i) as i,
       (f.psfmag_z-f.extinction_z) as z,
       f.ra as RA,
       f.dec as Dec,
       p.pmRA as pmRA,
       p.pmDec as pmDec,
       p.pmRAErr as pmRAErr,
       p.pmDecErr as pmDecErr
from ProperMotions p,
     fGetNearbyObjAllEq(250.423,36.460,120.0) r, -- NGC 6205 M 13
     PhotoObjAll f
where p.objid=r.objid
  and p.objid=f.objid
  and p.pmRAErr>0.0
  and p.pmDecErr>0.0
