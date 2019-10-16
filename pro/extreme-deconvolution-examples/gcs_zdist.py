import sys
import numpy as nu
import idlsave
from extreme_deconvolution import extreme_deconvolution

def gcs_zdist(ngauss=2):
    ####################
    #Load the data
    ####################
    i= idlsave.read('test.sav')
    ngcs= len(i.test.feh)

    ####################
    #Prepare for ED
    ####################
    ydata= nu.reshape(i.test.feh,(ngcs,1))
    ycovar= nu.reshape(i.test.e_feh,(ngcs,1))
    projection= nu.ones((ngcs,1,1))
    ydata= ydata.astype(nu.float64)
    ycovar=ycovar.astype(nu.float64)
    weight= nu.ones(ngcs)
    #weight= nu.random.random(ngcs)
    #weight= ydata
    logweight= False

    ####################
    #Initial conditions
    ####################
    if ngauss == 2:
        xamp= nu.array([0.5,0.5])
        xmean= nu.reshape(nu.array([0.,0.]),(ngauss,1))
        xcovar= nu.reshape(nu.array([1.,2.]),(ngauss,1,1))
    elif ngauss == 3:
        xamp= nu.array([0.33,0.33,0.34])
        xmean= nu.reshape(nu.array([0.,0.,0.]),(ngauss,1))
        xcovar= nu.reshape(nu.array([1.,2.,4.]),(ngauss,1,1))

    ####################
    #Run ED
    ####################
    print extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,weight=weight,logweight=logweight)
    print xamp
    print xmean
    print xcovar

if __name__ == '__main__':
    if len(sys.argv) > 1:
        gcs_zdist(int(sys.argv[1]))
    else:
        gcs_zdist()
