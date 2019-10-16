-- ----------------------------------------------------------------------------
-- NAME:
--  hogg_proper_motion.sql
-- PURPOSE:
--  select halo stars with proper motions from the CAS
-- NOTES:
--  - The "fGetObjFromRectEq" call apparently is much faster than explicit
--    RA,Dec cuts.
-- REVISION HISTORY:
--  2006-07-05  started - Hogg (NYU)
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
     fGetObjFromRectEq(140.0, 5.0,150.0,29.0) r, -- big area
--     fGetObjFromRectEq(140.0,17.0,150.0,22.0) r, -- area A
--     fGetObjFromRectEq(140.0, 5.0,150.0,10.0) r, -- area B1
--     fGetObjFromRectEq(140.0,24.0,150.0,29.0) r, -- area B2
     PhotoObj f
where p.objid=r.objid
  and p.objid=f.objid
  and p.pmRAErr>0.0
  and p.pmDecErr>0.0
