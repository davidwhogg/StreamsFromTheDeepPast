;+
;   NAME:
;      group_cmd
;   PURPOSE:
;      plot the cmd of a group
;   INPUT:
;      group - which group?
;      cut   - member prob has to be higher than this to be plotted
;      savefilename - filename that holds the velocity distribution
;   OUTPUT:
;      plot
;   REVISION HISTORY:
;      2009-09-08 - Written Bovy
;-
PRO GROUP_CMD, group, savefilename=savefilename, cut=cut

IF ~keyword_set(savefilename) THEN $
  savefilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
IF ~keyword_set(cut) THEN cut= 0.1
;;Name of the group
basefilename='cmd_'
IF group EQ 0 THEN BEGIN
    name='NGC 1901'
    plotfilename=basefilename+'ngc1901'
ENDIF ELSE IF group EQ 1 THEN BEGIN
    name='Sirius'
    plotfilename=basefilename+'sirius'
ENDIF ELSE IF group EQ 2 THEN BEGIN
    name='Arcturus'
    plotfilename=basefilename+'arcturus'
ENDIF ELSE IF group EQ 3 THEN BEGIN
    name='Pleiades'
    plotfilename=basefilename+'pleiades'
ENDIF ELSE IF group EQ 5 THEN BEGIN
    name='Background'
    plotfilename=basefilename+'back1'
ENDIF ELSE IF group EQ 6 THEN BEGIN
    name='Hyades'
    plotfilename=basefilename+'hyades'
ENDIF ELSE IF group EQ 7 THEN BEGIN
    name='Hercules'
    plotfilename=basefilename+'hercules'
ENDIF ELSE IF group EQ 8 THEN BEGIN
    name='Background'
    plotfilename=basefilename+'back2'
ENDIF ELSE IF group EQ 9 THEN BEGIN
    name='Background'
    plotfilename=basefilename+'back3'
ENDIF
plotfilename+='.ps'

;;Restore data
restore, filename='hip2-aumer-giants.sav'
nhip= n_elements(hip.hip)
ydata= transpose([[hip.vl], [hip.vb]])
ycovar= hip.vlvbc
projection= hip.nsm
;;Restore solution
restore, filename=savefilename

;;First compute all the probabilities
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
;;Sort the stars
IF group EQ 3 THEN BEGIN
    logpost[group,*]= alog(exp(logpost[group,*])+exp(logpost[group+1,*]))
ENDIF
sortindx= sort(logpost[group,*])

;;Now plot the cmd
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill
k_print, filename=plotfilename, xsize=3.375, ysize=4 ; actual width of ApJ column!
ncolors=100
loadct, 0, NColors=ncolors
djs_plot, [0.],[0.], yrange=[10,-6], color=0,$
  xtitle='B-V [mag]', ytitle='M_{V} [mag]',xrange=[-0.3,1.6], $
  xcharsize=1,ycharsize=1, position=[0.1,0.1,0.9,0.9],title=name
FOR ii=0L, nhip-1 DO BEGIN
    IF hip[sortindx[ii]].bvcolor EQ 0. THEN CONTINUE
    IF exp(logpost[group,sortindx[ii]]) GT cut AND hip[sortindx[ii]].bt NE 0 THEN BEGIN
        djs_oplot, [hip[sortindx[ii]].bvcolor], [hip[sortindx[ii]].vmag-$
                                                 5.*alog10(100./hip[sortindx[ii]].plx)], color=long(ncolors*(1.-exp(logpost[group,sortindx[ii]]))),$
          psym=8,symsize=0.25
    ENDIF
ENDFOR
;;overplot the mscuts
plot_mscuts
IF group EQ 3 THEN BEGIN
    datafilename='pleiades_isochrone.dat'
ENDIF ELSE IF group EQ 6 THEN BEGIN
    datafilename='hyades_isochrone.dat'
ENDIF ELSE IF group EQ 1 THEN BEGIN
    datafilename='sirius_isochrone.dat'
ENDIF ELSE IF group EQ 0 THEN BEGIN
    datafilename='ngc1901_isochrone.dat'
ENDIF
IF group EQ 3 or group EQ 6 OR group EQ 1 OR group EQ 0 THEN BEGIN
    iso= read_isochrone(datafilename)
    djs_oplot, iso.h09-iso.h10, iso.h10,psym=-3
ENDIF
;;Also plot the other sirius isochrone
IF group EQ 1 OR group EQ 3 OR group EQ 0 OR group EQ 6 THEN BEGIN
    IF group EQ 1 THEN BEGIN
        datafilename='sirius_isochrone2.dat'
    ENDIF ELSE IF group EQ 3 THEN BEGIN
        datafilename='pleiades_isochrone2.dat'
    ENDIF ELSE IF group EQ 0 THEN BEGIN
        datafilename='ngc1901_isochrone2.dat'
    ENDIF ELSE IF group EQ 6 THEN BEGIN
        datafilename='hyades_isochrone2.dat'
    ENDIF
    iso= read_isochrone(datafilename)
    djs_oplot, iso.h09-iso.h10, iso.h10,psym=-3
ENDIF
k_end_print
END
