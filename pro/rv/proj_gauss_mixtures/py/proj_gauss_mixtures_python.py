#Wrapper routine for the proj_gauss_mixtures c code
import projGaussMixtures
import scipy as sc

def proj_gauss_mixtures_python(ydata,ycovar,projection,amp,xmean,
                               xcovar,fixamp='',fixmean='',fixcovar='',
                               tol=0.000001,maxiter=10000000000,likeonly=False,
                               w=0.,quiet=False,logfile='',splitnmerge=0,maxsnm=False):
    """
    NAME:
       proj_gauss_mixtures_python
    PURPOSE:
       iterate on projected gaussian mixtures using the C program
    INPUT:
       ydata      - [ndimy,ndata] observed values
       ycovar     - [ndimy,ndimy,ndata] observed values' uncertainties
       projection - [ndimx, ndimy, ndata] non-square matrices
                    implementing the projection
       amp        - [ngauss] list of relative amplitudes
       xmean      - [ndimx,ngauss] list of initial guesses for the means
       xcovar     - [ndimx,ndimx,ngauss] list of initial guesses for the
                    covar
       
    OPTIONAL INPUTS:
       fixamp     - [ngauss] list of integers: 0 to update amp, 1 not to
       fixmean    - [ngauss] list of integers: 0 to update xmean, 1 not to
       fixcovar   - [ngauss] list of integers: 0 to update xcovar, 1 not
                    to
       tol        - tolerance (convergence iff difference in avgloglike <
                    tol)
       maxiter    - maximum number of iterations
       w          - regularization parameter for the covariances
       logfile    - basename for several logfiles (_c.log has output from
                    the c-routine; _loglike.log has the log likelihood
                    path of all the accepted routes, i.e. only parts
                    which increase the likelihood are included, during
                    splitnmerge)
       splitnmerge- depth to go down the splitnmerge path
    BOOLS:
       maxsnm     - use the maximum number of split 'n' merge steps,
                    K*(K-1)*(K-2)/2
       likeonly   - only compute and return the avgloglikedata; do not
                    update
       quiet      - don't print messages
    OUTPUT:
       (avgloglikedata,xmean,xcovar,xamp)  - Average log-likelihood of the data after
                   convergence  +updated xmean, xcovar and amp...
    HISTORY:
       2009-11-23 - Written - Bovy (NYU)
    """
    #Dimensions
    ndata= ydata.shape[1]
    ngauss= len(amp)
    dy= ydata.shape[0]
    dx= xmean.shape[0]
    #Check inputs
    temp_fixamp=''
    for ii in range(ngauss):
        if fixamp[ii] == 'False' or fixamp[ii] == 0:
            temp_fixamp+='0'
        else:
            temp_fixamp+='1'
    temp_fixmean=''
    for ii in range(ngauss):
        if fixmean[ii] == 'False' or fixmean[ii] == 0:
            temp_fixmean+='0'
        else:
            temp_fixmean+='1'
    temp_fixcovar=''
    for ii in range(ngauss):
        if fixcovar[ii] == 'False' or fixcovar[ii] == 0:
            temp_fixcovar+='0'
        else:
            temp_fixcovar+='1'
    if maxsnm:
        splitnmerge= long(ngauss*(ngauss-1)*(ngauss-2)/2)
    if likeonly:
        temp_likeonly= '1'
    else:
        temp_likeonly= '0'

    #ydata
    temp_ydata= projGaussMixtures.doubleArray(ndata*dy)
    for ii in range(dy):
        for jj in range(ndata):
            temp_ydata[ii*ndata+jj]= ydata[ii,jj]
    #ycovar
    temp_ycovar= projGaussMixtures.doubleArray(ndata*dy*dy)
    temp_ycovar[0:-1]= sc.ravel(ycovar)
    #projection
    temp_projection= projGaussMixtures.doubleArray(ndata*dy*dx)
    temp_projection[0:-1]= sc.ravel(projection)
    #xmean
    temp_xmean= projGaussMixtures.doubleArray(dx*ngauss)
    temp_xmean[0,:-1]= sc.ravel(xmean)
    #xcovar
    temp_xcovar= projGaussMixtures.doubleArray(dx*dx*ngauss)
    temp_xcovar[0:-1]= sc.ravel(xcovar)
    #amp
    temp_amp= projGaussMixtures.doubleArray(ngauss)
    temp_amp[0:-1]= sc.ravel(amp)
    #avgloglikedata
    avgloglikedata= projGaussMixtures.doubleArray(1)
    avgloglikedata[0]= 0.
    
    print temp_ydata, temp_ycovar, temp_projection, temp_xmean, temp_xcovar, temp_amp

    if quiet:
        pass
    else:
        print "Loading C-implementation to run the projected_gauss_mixtures algorithm..."
    projGaussMixtures.proj_gauss_mixtures(temp_ydata,temp_ycovar,temp_projection,ndata,dy,
                                          temp_amp,temp_xmean,temp_xcovar,dx,ngauss,
                                          temp_fixamp,temp_fixmean,temp_fixcovar,avgloglikedata,
                                          tol,maxiter,temp_temp_likeonly,w,
                                          logfilename,splitnmerge,convlogfilename)
    if quiet:
        pass
    else:
        print "Successfully ran the C-code"
    avgloglikedata= avgloglikedata[0]
    xmean= sc.array(xmean)
    xcovar= sc.array(xcovar)
    amp= sc.array(amp)
    
    return(avgloglikedata,xmean,xcovar,amp)
