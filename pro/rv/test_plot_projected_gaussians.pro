;+
;   NAME:
;      test_plot_projected_gaussians
;
;   PURPOSE:
;      test the routine 'plot_projected_gaussians' by plotting a sum
;      of gaussians with a simple structure
;
;   REVISION HISTORY:
;      2008-06-* - Written (Bovy)
;-
PRO TEST_PLOT_PROJECTED_GAUSSIANS

mean=dblarr(2,3)
mean[*,0]=[0D,0D]
mean[*,1]=[2D,3D]
mean[*,2]=[-2D,-3D]

covar=dblarr(2,2,3)
covar[*,*,0]=[[2,1],[1,3]]
covar[*,*,1]=[[2,.5],[.5,2]]
covar[*,*,2]=[[2,0],[0,4]]

amp=[.5D,.25D,0D]

projection=[[1D,0D],[0D,1D]]

plot_projected_gaussians, mean, covar, projection, xrange=[-5,5], yrange=[-5,5], amp=amp

END
