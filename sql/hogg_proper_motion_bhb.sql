-- ----------------------------------------------------------------------------
-- NAME:
--  hogg_proper_motion_bhb.sql
-- PURPOSE:
--  select a superset of BHB stars with proper motions from the CAS
-- NOTES:
--  - The "fGetObjFromRectEq" call apparently is much faster than explicit
--    RA,Dec cuts.
-- REVISION HISTORY:
--  2006-07-05  works - Hogg (NYU)
-- ----------------------------------------------------------------------------
select * from
( select (f.psfmag_u-f.extinction_u) as u,
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
       fGetObjFromRectEq(130.0,0.0,140.0,10.0) r,
       PhotoObj f
  where p.objid=r.objid
    and p.objid=f.objid
    and p.pmRAErr>0.0
    and p.pmDecErr>0.0
) as pm
where (pm.u-pm.g)>(0.80)
  and (pm.u-pm.g)<(1.35)
  and (pm.g-pm.r)>(-0.40)
  and (pm.g-pm.r)<(0.00)
