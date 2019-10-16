-- ----------------------------------------------------------------------------
-- NAME:
--  hogg_proper_motion_quasars
-- PURPOSE:
--  select quasars with proper motions from the CAS
-- NOTES:
--  - I am not sure about the "p.objid=f.SpecBestObjid"; I think it would be
--    better to match these on (RA,Dec).
--  - The "fGetObjFromRectEq" call apparently is much faster than explicit
--    RA,Dec cuts.
-- REVISION HISTORY:
--  2006-07-04  started - Hogg (NYU)
-- ----------------------------------------------------------------------------
select f.SpecRa as RA,
       f.SpecDec as Dec,
       f.SpecZ as z,
       f.SpecZerr as zErr,
       f.SpecZWarning as zWarn,
       p.pmRA as pmRA,
       p.pmDec as pmDec,
       p.pmRAErr as pmRAErr,
       p.pmDecErr as pmDecErr
from ProperMotions p,
     fGetObjFromRectEq(130.0,0.0,140.0,10.0) r,
     QsoSpec f
where p.objid=r.objid
  and p.objid=f.SpecBestObjid
  and p.pmRAErr>0.0
order by f.SpecZ desc
