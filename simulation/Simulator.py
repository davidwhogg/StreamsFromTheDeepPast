# To do:
#	Make this code prettier and more human readable by:
#
#	- Lets work on comments so that they aren't trivial, 
#	and so they actually explain the code. We'll try to avoid
#	pointing out the obvious like "this function sets the 
#	POTENTIAL variable", because duh, the function says 
#	"POTENTIAL = potential"...
#
#	- Everything is version controlled, so I say we delete
#	anything that is "broken" or "outdated" or "not used"
#
#	- As one of my general coding philosophies, any function
#	that has over ~20 lines should be broken up into more
#	generalized functions (obviously there are exceptions)


import pickle

#import os
import numpy
import matplotlib
import matplotlib.pyplot as plt
#import time
#import math
import logging


#CODE RELATED GLOBAL CONSTANTS
#DARK = -9999999999
KEPLER = 11111111
LN = 22222222
FLATLN = 33333333

# IDENTITY = [ [1,0,0,0,0,0],
#              [0,1,0,0,0,0],
#              [0,0,1,0,0,0],
#              [0,0,0,1,0,0],
#              [0,0,0,0,1,0],
#              [0,0,0,0,0,1] ]

# AJ: make the units known and explicit everywhere.
# distances: kpc
# velocities: kpc / myr
# times: myr
# accelerations: (kpc / Myr^2)
# masses: Solar Masses REDUNDANT ish

#CONVERSTION FACTORS
SECPERMYR=60*60*24*365.25*(10**6)
KMPERKPC = 3.08568025e16
KGPERGIGASOLARMASS = 1.98892e30

# PHYSICAL CONSTANTS
NEWTONG = 6.67300e-11 * (1000*KMPERKPC)**(-3) * KGPERGIGASOLARMASS * SECPERMYR**2 #m^3 kg^-1 s-2 --> kpc^3 gigaSolarMass^-1 Myr^-2
V2POTENTIAL = 224.0 * 224.0 * KMPERKPC**(-2) * SECPERMYR**2 #((km/s)^2 ---> (kpc/ Myr)^2
#V2POTENTIAL = 0
GALMASS = 10**10 #gigasolarmass
VELGAUSSSTDEV = 100 *SECPERMYR/KMPERKPC
VELGAUSSNORM = numpy.sqrt(V2POTENTIAL) #224.0 * SECPERMYR/KMPERKPC # km / s --> KPC/MYR
POSGAUSSSTDEV = 20.0 # kpc
CLUSTERPOSNORM = .2 #kpc
CLUSTERVELNORM = 1 * SECPERMYR/KMPERKPC#km/s --> KPC/MYR
MASPERRAD = 206264806 #mas/rad
GALACTICCENTER = numpy.array([0,0,0]) # kpc
FLATPARAM = 1.0 #.87 # dimensionless
EARTH = numpy.array([-10,0,0,0,0.205,0]) #(kpc,kpc,kpc,kpc/myr,kpc/myr,kpc/myr)

#PARAMETERS
#PERTURBERMASS = 0#10 # GIGASOLARMASS
#STARMASS = 1e-9 # GIGASOLARMASS irrelevant parameter

POTENTIAL = LN # KEPLER



def setPotential(potential):
    global POTENTIAL
    POTENTIAL = potential
#    print POTENTIAL

##this function returns a random velocity vector
def get_velocity():
    return numpy.random.randn(3)*VELGAUSSNORM #kpc/myr
##this function returns a random position vector
def get_position():
    return numpy.random.randn(3)*POSGAUSSSTDEV #kpc


def get_acceleration(pos,flatparam=1,vcircsquared=V2POTENTIAL):
    if POTENTIAL == LN:
        r = GALACTICCENTER - pos[:3] #kpc
        sqnormr = numpy.dot(r,r) #kpc
        return V2POTENTIAL*r/(sqnormr) ##kpc/myr^2
    elif POTENTIAL == FLATLN:
        r = GALACTICCENTER - pos[:3] #kpc
        r[2] /= flatparam
        sqnormr = numpy.dot(r,r) #kpc
        rtilde = r
        rtilde[2] /= flatparam
        return vcircsquared*rtilde/(sqnormr)
    elif POTENTIAL == KEPLER:
        r = GALACTICCENTER - pos[:3] #kpc
        sqnormr = numpy.dot(r,r) #kpc
        return NEWTONG*GALMASS*r/(sqnormr**(1.5)) ## kpc/myr^2
    else:
        raise Exception("INCORRECT POTENTIAL")
        
class ColdStream(object):

### the initializer takes a list of values and depending on the size of the list determines which of the types of approximations to a stream to make

    def __init__(self,numForward=10,numBackward=10,timestep=0,fiducialPoint=None,flatparam=1,vcircsquared=V2POTENTIAL,thicknessEigs=[0,0]):
        self.numForward = numForward
        self.numBackward = numBackward
        self.dt = timestep
        self.thicknessEigs = thicknessEigs
        

        if fiducialPoint == None:
            fiducialPoint = numpy.zeros(6)
            fiducialPoint[:3] = get_position()[:]
            fiducialPoint[3:] = get_velocity()[:]

        self.fiducialPoint = Star(fiducialPoint,thicknessEigs = self.thicknessEigs,covarianceMatrix=None)
    

        self.stars = [self.fiducialPoint]
        self.stars_pos = [self.fiducialPoint.pos]
        self.fiducialPointIndex = 0


        self.potential = POTENTIAL
        self.flatparam = flatparam
        self.vcircsquared = vcircsquared
#        print POTENTIAL

    def forwardEuler(self,dt=None, numForward=None,numBackward=None):
        if dt == None:
            dt = self.dt
        if numForward == None:
            numForward = self.numForward
        if numBackward == None:
            numBackward = self.numBackward

        pos = self.stars_pos[0].copy()

        for i in range(numBackward):
            a = get_acceleration(pos)
            pos[3:] -= a*dt ## kpc/myr^2 * myr = kpc/myr
            pos[:3] -= pos[3:]*dt ## kpc/myr * myr = kpc
            self.insertStar(0,pos.copy())
            self.fiducialPointIndex += 1
        pos = self.stars_pos[-1].copy()

        for i in range(numForward):
            a = get_acceleration(pos)
            pos[3:] += a*dt ## kpc/myr^2 * myr = kpc/myr
            pos[:3] += pos[3:]*dt ## kpc/myr * myr = kpc
            self.insertStar(len(self.stars),pos.copy())

    def integrate(self,dt=None,numForward=None,numBackward=None):
#    def leapFrog(self,dt=None, numForward=None,numBackward=None):
        if dt == None:
            dt = self.dt
        if numForward == None:
            numForward = self.numForward
        if numBackward == None:
            numBackward = self.numBackward


        pos = self.stars_pos[0].copy()
        a_i = get_acceleration(pos,self.flatparam,self.vcircsquared)
        for i in range(numBackward):
            pos[:3] += -pos[3:]*dt + a_i*(dt**2)/(2.0) ## kpc/myr * myr = kpc
            a_iplus1 = get_acceleration(pos,self.flatparam,self.vcircsquared)
            pos[3:] -= (a_i + a_iplus1)*dt/(2.0) ## kpc/myr^2 * myr = kpc/myr
            self.insertStar(0,pos.copy())
            self.fiducialPointIndex += 1
            a_i = a_iplus1


        pos = self.stars_pos[-1].copy()
        a_i = get_acceleration(pos,self.flatparam,self.vcircsquared)
        for i in range(numForward):
            pos[:3] += pos[3:]*dt + a_i*(dt**2)/(2.0) ## kpc/myr^2 * myr = kpc/myr
            a_iplus1 = get_acceleration(pos,self.flatparam,self.vcircsquared)
            pos[3:] += (a_i+a_iplus1)*dt/(2.0) ## kpc/myr * myr = kpc
            self.insertStar(len(self.stars),pos.copy())
            a_i = a_iplus1
    
    def insertStar(self,index,pos):
        self.stars.insert(index,Star(pos,self.thicknessEigs))
        self.stars_pos.insert(index,pos)
## this function genenerates the plot data for the streams 
    def plot(self):
        plot = [[],[],[],[],[],[]] 
        for i in self.stars:
            plot[0].append(i.pos[0])
            plot[1].append(i.pos[1])
            plot[2].append(i.pos[2])
            plot[3].append(i.pos[3])
            plot[4].append(i.pos[4])
            plot[5].append(i.pos[5])

        return plot

    def copy(self):
        return pickle.loads(pickle.dumps(self))

    def get_params(self):
        params = [self.numForward+self.numBackward+1]
        params.extend(self.fiducialPoint.pos)
        params.extend([self.flatparam, self.vcircsquared,self.thicknessEigs[0]])
        return params


    def potentialToString(self,abr=False):
        if abr:
            if self.potential == LN:
                str = 'LN'
            elif self.potential == FLATLN:
                str = 'FLN'
            elif self.potential == KEPLER:
                str = 'KEP'
            else:
                str = 'No'
                print self.potential
        else:
            if self.potential == LN:
                str = 'Spherical Log'
            elif self.potential == FLATLN:
                str = 'Flattened Log'
            elif self.potential == KEPLER:
                str = 'Kepler'
            else:
                str = 'No'
                print self.potential
        return str

class DataColdStream(ColdStream):
    
    def __init__(self,numForward=10,numBackward=10,timestep=0,fiducialPoint=None,numRemaining=-1,covarianceOfFiducialPoint=1,hasNoise=False,flatparam=1, vcircsquared=V2POTENTIAL,thicknessEigs=[.000001,0],ln_rad_pos_err=.1):
        ColdStream.__init__(self,numForward,numBackward,timestep,fiducialPoint,flatparam=flatparam,vcircsquared=vcircsquared,thicknessEigs=thicknessEigs)
        
        self.hasNoise = hasNoise
        self.fiducialPoint=Star(self.fiducialPoint.pos,thicknessEigs=thicknessEigs,covarianceMatrix=1)#,hasNoise=False)
        self.integrate()
        if numRemaining == -1:
            numRemaining=numForward+numBackward+1

        for i in range(numForward+numBackward+1-numRemaining):
            self.stars.pop(numpy.random.randint(0,len(self.stars)))

        if self.hasNoise:
            for i in self.stars:
                i.addNoise(ln_rad_pos_err)
       
    def insertStar(self,index,pos):
        self.stars.insert(index,Star(pos.copy(),self.thicknessEigs,1))#,False))#self.hasNoise))
#        print 'det', numpy.linalg.det(self.stars[index].covarianceMatrix)

class ModelColdStream(ColdStream):
    def __init__(self,numForward=10,numBackward=10,timestep=0,fiducialPoint=None,flatparam=1,vcircsquared=V2POTENTIAL,thicknessEigs=[0,0]):
        ColdStream.__init__(self,numForward,numBackward,timestep,fiducialPoint,flatparam=flatparam,vcircsquared=vcircsquared,thicknessEigs=thicknessEigs)
        self.integrate()

    def changeSampleRate(self,sampleRate=1.0):
        self.dt /= sampleRate*1.0
        self.numForward = int(self.numForward*sampleRate)
        self.numBackward = int(self.numBackward*sampleRate)
        self.stars = [self.fiducialPoint]
        self.fiducialPointIndex = 0

        self.integrate()
        
    def addStars(self, numFront=0,numBack=0):
        self.integrate(numForward=numFront,numBackward=numBack)
        self.numForward += numFront
        self.numBackward += numBack

    def copy(self):
        stream = ModelColdStream(numForward=self.numForward,numBackward=self.numBackward,timestep=self.dt,fiducialPoint=self.fiducialPoint.pos,flatparam=self.flatparam,vcircsquared=self.vcircsquared,thicknessEigs=self.thicknessEigs)
        return stream
        
class Star(object):
    def __init__(self,pos,thicknessEigs=[0,0],covarianceMatrix=None):#,hasNoise=False):
        self.pos = pos
        if covarianceMatrix == 1:
            self.covarianceMatrix = numpy.eye(6)
            self.changeCovarianceMatrix()

#             self.covarianceMatrix = numpy.array([ [1,0,0,0,0,0],
#                                                   [0,1,0,0,0,0],
#                                                   [0,0,1,0,0,0],
#                                                   [0,0,0,.01**2,0,0],
#                                                   [0,0,0,0,.01**2,0],
#                                                  [0,0,0,0,0,.01**2] ])
#            print numpy.dot(self.covarianceMatrix,numpy.array([0,0,0,1,0,0]))
#            self.covarianceMatrix = numpy.array(IDENTITY)
        
        else:
            self.covarianceMatrix = covarianceMatrix
        

        self.thicknessTensor = numpy.diag([thicknessEigs[0],thicknessEigs[0],thicknessEigs[0],thicknessEigs[1],thicknessEigs[1],thicknessEigs[1]])
        

#        if hasNoise:
#            self.addNoise()
    

    def addNoise(self,ln_rad_pos_err=.1):
        old_pos = self.pos
        self.pos= numpy.random.multivariate_normal(self.pos,self.thicknessTensor)
        
        self.changeCovarianceMatrix(ln_rad_pos_err=ln_rad_pos_err)
##### INLINE FUNCTIONAL TEST: undo comment out after use

        old_pos = self.pos
        error = numpy.random.multivariate_normal(numpy.zeros(6),self.covarianceMatrix)
        self.pos = self.pos + error

        
        
    def __str__(self):
        s =  "velocity: "+ str(self.pos[3:]) + "\nPosition: " + str(self.pos[:3]) +'\n'
        return s

    def changeCovarianceMatrix(self, observerPos=EARTH, ang_err=5*10**(-7),prop_mot_err=2.0*(10**-3)*(1.0/3600.0)*(numpy.pi/180.0)*10**6, rad_vel_err= 10*SECPERMYR/KMPERKPC, ln_rad_pos_err=.1):
        self.ln_rad_pos_err = ln_rad_pos_err
        relpos = self.pos-observerPos
        r_hat = relpos[:3]
        distance = numpy.linalg.norm(r_hat)
        r_hat /= distance
        theta_hat = numpy.cross(r_hat,[1,0,0])
        epsilon = numpy.linalg.norm(theta_hat)
        if epsilon  <= 10**(-8):
            raise Exception("data on x axis")
        theta_hat /= epsilon
        phi_hat = numpy.cross(r_hat,theta_hat)
        epsilon = numpy.linalg.norm(phi_hat)
        if numpy.abs(epsilon-1.) > 10**(-15):
            raise Exception("cross products foobar")
        zero = [0,0,0]
        # here are the six unit vectors
        e1 = numpy.append(r_hat,zero)
        e2 = numpy.append(theta_hat,zero)
        e3 = numpy.append(phi_hat,zero)
        e4 = numpy.append(zero,r_hat)
        e5 = numpy.append(zero,theta_hat)
        e6 = numpy.append(zero,phi_hat)
#        print e1,e2,e3,e4,e5,e6
        # compute vperp because of distance error propagation into transverse velocity (note distance*e1+vperp vectors below)
        vperp = numpy.dot(relpos,e5)*e5 + numpy.dot(relpos,e6)*e6
#        print vperp
#        print  numpy.outer(distance*e1,distance*e1) + numpy.outer(vperp,vperp)
#        raise IOError
        theta = numpy.arctan2(self.pos[1],self.pos[0])
        vtheta = numpy.arctan2(self.pos[4],self.pos[3])

        


        self.covarianceMatrix = ((ln_rad_pos_err)**2)*(numpy.outer(distance*e1,distance*e1) + numpy.outer(vperp,vperp)) + ((distance*ang_err)**2)*(numpy.outer(e2,e2) + numpy.outer(e3,e3)) + (rad_vel_err**2)*numpy.outer(e4,e4) + ((distance*prop_mot_err)**2)*(numpy.outer(e5,e5) + numpy.outer(e6,e6) )

#        print '\n\n\n###### Covariance Martrix\n',self.covarianceMatrix
#        print 'd*rad_vel_err', (distance*rad_vel_err)
#        print 'd*ang_err', (distance*ang_err)
#        print 'dprop_mot_err', distance*prop_mot_err
#        print '(D*prop_mot_err)^2',(distance*prop_mot_err)**2
#        print 'dD/D * D*vperp', (ln_rad_pos_err*distance*numpy.linalg.norm(vperp))
#        print 'rad_vel_err',rad_vel_err
#        self.covarianceMatrix = numpy.eye(6)
#        self.covarianceMatrix = ((ln_rad_pos_err)**2)*(numpy.outer(distance*e1,distance*e1)) + ((distance*ang_err)**2)*(numpy.outer(e2,e2) + numpy.outer(e3,e3)) + (rad_vel_err**2)*numpy.outer(e4,e4) + ((distance*prop_mot_err)**2)*(numpy.outer(e5,e5) + numpy.outer(e6,e6) )
        
class Simulation:

    #numForward=10,numBackward=10,timestep=0,fiducialPoint=-1,numRemaining=-1,integrationScheme='forwardEuler',covarianceOfFiducialPoint=None):

    def __init__(self,reseed=True,seedVal=20,numStreams=1,timestep=1,potential=LN, fiducialPoints=None,numForward=10,numBackward=10,integrationScheme='forwardEuler',numRemaining=-1,isData=False,hasNoise=False, flatparam=1, vcircsquared=V2POTENTIAL,thicknessEigs=[0,0]):
        if reseed:
            numpy.random.seed(seedVal)
        self.numStreams = numStreams
        self.dt = timestep
        setPotential(potential)
        self.streams= []
        if fiducialPoints == None:
            if isData:
                for i in range(self.numStreams):
                    self.streams.append(DataColdStream(numForward,numBackward,self.dt,fiducialPoints,numRemaining=numRemaining,hasNoise=hasNoise,flatparam=flatparam,vcircsquared=vcircsquared, thicknessEigs=thicknessEigs))
            else:
                for i in range(self.numStreams):
                    self.streams.append(ModelColdStream(numForward,numBackward,self.dt,fiducialPoints,flatparam=flatparam,vcircsquared=vcircsquared,thicknessEigs=thicknessEigs))
        else:
            if isData:
                for i in range(self.numStreams):
                    self.streams.append(DataColdStream(numForward,numBackward,self.dt,fiducialPoints[i],numRemaining=numRemaining,hasNoise=hasNoise,flatparam=flatparam,vcircsquared=vcircsquared,thicknessEigs=thicknessEigs))
            else:
                 for i in range(self.numStreams):
                    self.streams.append(ModelColdStream(numForward,numBackward,self.dt,fiducialPoints[i],flatparam=flatparam, vcircsquared = vcircsquared,thicknessEigs=thicknessEigs))

        
    def plotH(self,streamnum):
        T = []
        U = []
        H = []
        for i in self.streams[streamnum].stars:
            T.append(.5*numpy.dot(i.pos[3:],i.pos[3:]))
            r =GALACTICCENTER-i.pos[:3]
            if POTENTIAL == LN:
                U.append(V2POTENTIAL*numpy.log(numpy.dot(r,r))/2.0)
            elif POTENTIAL == FLATLN:
                r[2] /= FLATPARAM
                sqnormr = numpy.dot(r,r) #kpc
                U.append(V2POTENTIAL*numpy.log(numpy.dot(r,r))/(2.0))
            else:
                U.append(NEWTONG*GALMASS/numpy.linalg.norm(r))
            H.append(T[-1] + U[-1])

        print "deltaH: %.8f" % (H[0]-H[-1])



#        print "t", T[self.streams[streamnum].fiducialPointIndex]
#        print "u", U[self.streams[streamnum].fiducialPointIndex]
#        print 'h', H[self.streams[streamnum].fiducialPointIndex]
#        print 'H', T[self.streams[streamnum].fiducialPointIndex] + U[self.streams[0].fiducialPointIndex]
#        print 'id', self.streams[streamnum].fiducialPointIndex
                     
        plt.plot(range(len(self.streams[streamnum].stars)),T, label="t")
        plt.plot(range(len(self.streams[streamnum].stars)),U, label='u')
        plt.plot(range(len(self.streams[streamnum].stars)),H, label='h')
#        plt.plot([self.streams[streamnum].fiducialPointIndex],[T[self.streams[streamnum].fiducialPointIndex]],'bo', label="fidpt T")
#        plt.plot([self.streams[streamnum].fiducialPointIndex],[U[self.streams[streamnum].fiducialPointIndex]], label="fidpt U")
#        plt.plot([self.streams[streamnum].fiducialPointIndex],[H[self.streams[streamnum].fiducialPointIndex]], label="fidpt H")
        plt.legend()
        plt.savefig("energies")

        plt.clf()

    def plotAllStreams(self,toPlotAgainst=[],againstPotential = 'Spherical Log',label='',plotfid=False,fid=[]):
        fig = plt.figure(figsize=(15,15))
        fig.subplots_adjust(wspace=.2,hspace=.2)
        print len(self.streams)
        stream_potential_name = self.streams[0].potentialToString()
#        if toPlotAgainst==[]:
#            fig.suptitle("True:" + stream_potential_name, fontsize = 30)
#        else:
#            fig.suptitle("True: " + againstPotential + " Model: " + stream_potential_name, fontsize = 30)
    
        streamsPlot = []
        for i in self.streams:
            streamsPlot.append(i.plot())
        streamsPlot = numpy.array(streamsPlot)
        ax = [[],[],[],[],[],[]]
        for i in range(6):
            for j in range(6):
                ax[i].append( fig.add_subplot(6,6,j+1+6*i))
                if len(toPlotAgainst) > 0:
                    plt.plot(toPlotAgainst[j],toPlotAgainst[i], 'r.')#, alpha=1.0)
                    for stream in streamsPlot:
                        plt.plot(stream[j],stream[i], 'k.')#,alpha=.1)
                else:
                    for stream in streamsPlot:
                        plt.plot(stream[j],stream[i],'r.')#,alpha=1)



                xmin = numpy.amin(streamsPlot[:,j])
                xmax = numpy.amax(streamsPlot[:,j])
                ymin = numpy.amin(streamsPlot[:,i])
                ymax = numpy.amax(streamsPlot[:,i])

                xmin -= numpy.abs(xmin*.01)
                xmax += numpy.abs(xmax*.01)
                ymin -= numpy.abs(ymin*.01)
                ymax += numpy.abs(ymax*.01)

                if plotfid:

                    plt.plot([fid[j]],[fid[i]],c='g',marker='x',ls='None')
                    
 
#                 if i < 3:
#                     ymin=-50
#                     ymax=50
#                 else:
#                     ymin=-1
#                     ymax=1
#                 if j <3:
#                     xmin=-50
#                     xmax=50
#                 else:
#                     xmin=-1
#                     xmax=1
                plt.axis([xmin,xmax,ymin,ymax])
                format_plot(ax,i,j,6)
#                 if j == 0:
#                     unitsy = ''
# #                     if i <3:
# #                         unitsy = '$\text{(kpc)}$'
# #                     else:
# #                         unitsy = '$\text{(kpc/Myr)}'
#                     plt.ylabel(('$x_%i$ ' % (i+1))+unitsy,fontsize = 20)
#                     plt.setp(ax.get_yticklabels(),fontsize = 12)

#                 else:
#                     plt.setp(ax.get_yticklabels(),visible=False) 

#                 if i == 5:
#                     unitsx=''
# #                     if j < 3:
# #                         unitsx = '$\text{(kpc)}$'
# #                     else:
# #                        unitsx = '$\text{(kpc/Myr)}$'
                        
#                     plt.xlabel(('$x_%i$ ' % (j+1))+unitsx,fontsize=20)
#                     plt.setp(ax.get_xticklabels(),fontsize=12)
#                 else:
#                     plt.setp(ax.get_xticklabels(),visible=False)


                

        plt.savefig(label+'StreamsPlot')
        plt.close(fig)
        plt.clf()
        


        
    def plotStream(self,streamnum,label=0):
        stream1Plot = self.streams[streamnum].plot()
        for j in range(3):
            plt.subplot(3,3,j+1)
            plt.plot(stream1Plot[j],stream1Plot[(j+1)%3], 'o', label='data')
#            plt.plot([stream1Plot[j][self.streams[streamnum].fiducialPointIndex]],[stream1Plot[(j+1)%3][self.streams[streamnum].fiducialPointIndex]],label='fiducialPoint')
            plt.xlabel('R_%i' % j)
            plt.ylabel('R_%i' % ((j+1)%3))
    
        for j in range(3):
            plt.subplot(3,3,j+4)
            plt.plot(stream1Plot[j],stream1Plot[j+3],'o',label='data')
            plt.xlabel('R_%i' % j)
            plt.ylabel('V_%i' % j)

        for j in range(3):
            plt.subplot(3,3,j+7)
            plt.plot(stream1Plot[j+3],stream1Plot[((j+1)%3)+3], 'o', label='data')
            plt.xlabel('V_%i' % j)
            plt.ylabel('V_%i' % ((j+1)%3))
        
        plt.savefig('PlotStream-'+str(streamnum)+'-'+str(label))

        plt.clf()
        


def plotMultipleStreams(noisy_streams,clean_streams,label=''):
    fig = plt.figure(figsize=(22,22))
    fig.subplots_adjust(wspace=.1,hspace=.1)   
    fig.suptitle('yo')
    streamsPlot = []
    for i in range(len(noisy_streams)):
#    def __init__(self,numForward=10,numBackward=10,timestep=0,fiducialPoint=None,flatparam=1,vcircsquared=V2POTENTIAL,thicknessEigs=[0,0]):
        clean_stream = clean_streams[i]
        redo = ModelColdStream(numForward=clean_stream.numForward, numBackward=clean_stream.numBackward,timestep = clean_stream.dt, fiducialPoint=clean_stream.fiducialPoint.pos,flatparam=clean_stream.flatparam,vcircsquared=clean_stream.vcircsquared,thicknessEigs=clean_stream.thicknessEigs)
        streamsPlot.append([noisy_streams[i].plot(),redo.plot()])
   #     print noisy_streams[i].plot(), redo.plot()
  #  streamsPlot = numpy.array(streamsPlot)
    
    
    for i in range(6):
        for j in range(6):
            ax = fig.add_subplot(6,6,j+1+6*i)
            for k in range(len(noisy_streams)):
                plt.plot(streamsPlot[k][0][j],streamsPlot[k][0][i],'k.',alpha=(.25+ .25/len(noisy_streams) * (k+1)), ls = 'None')
                plt.plot(streamsPlot[k][1][j],streamsPlot[k][1][i],'k-',alpha=(.25 + .25/len(noisy_streams)* (k+1)))
                plt.plot([clean_streams[k].fiducialPoint.pos[j]],[clean_streams[k].fiducialPoint.pos[i]],'gx',ls='None')
                
                
            if i < 3:
                ymin=-50
                ymax=50
            else:
                ymin=-1
                ymax=1
            if j <3:
                xmin=-50
                xmax=50
            else:
                xmin=-1
                xmax=1
            plt.axis([xmin,xmax,ymin,ymax])
            if j == 0:
                unitsy = ''
                plt.ylabel(('$x_%i$ ' % (i+1))+unitsy,fontsize = 20)
                plt.setp(ax.get_yticklabels(),fontsize=12)

            else:
                plt.setp(ax.get_yticklabels(),visible=False) 

            if i == 5:
                unitsx=''
                plt.xlabel(('$x_%i$ ' % (j+1))+unitsx,fontsize=20)
                plt.setp(ax.get_xticklabels(),fontsize=12)
            else:
                plt.setp(ax.get_xticklabels(),visible=False)


                

    plt.savefig(label+'MStreamsPlot')
    plt.close(fig)
    plt.clf()

    fig = plt.figure(figsize=(22,22))
    fig.subplots_adjust(wspace=.1,hspace=.1)   
    fig.suptitle('yo')
        
    for i in range(6):
        for j in range(6):
            ax = fig.add_subplot(6,6,j+1+6*i)
            for k in range(len(noisy_streams)):
                plt.plot([clean_streams[k].fiducialPoint.pos[j]],[clean_streams[k].fiducialPoint.pos[i]],'go',ms = 10.0/len(noisy_streams) * (k+1),ls='None')
                
                
            if i < 3:
                ymin=-50
                ymax=50
            else:
                ymin=-1
                ymax=1
            if j <3:
                xmin=-50
                xmax=50
            else:
                xmin=-1
                xmax=1
            plt.axis([xmin,xmax,ymin,ymax])
            if j == 0:
                unitsy = ''
                plt.ylabel(('$x_%i$ ' % (i+1))+unitsy,fontsize = 20)
                plt.setp(ax.get_yticklabels(),fontsize=12)

            else:
                plt.setp(ax.get_yticklabels(),visible=False) 

            if i == 5:
                unitsx=''
                plt.xlabel(('$x_%i$ ' % (j+1))+unitsx,fontsize=20)
                plt.setp(ax.get_xticklabels(),fontsize=12)
            else:
                plt.setp(ax.get_xticklabels(),visible=False)


                

    plt.savefig(label+'FidMStreamsPlot')
    plt.close(fig)
    plt.clf()



def copy(stream):
    return pickle.loads(pickle.dumps(stream))
# This test demonstrates that energy is conserved by taking the simple case of a circular orbit in a keplerian potential

def EnergyConservationTest(isData,hasNoise,numRemaining=-1):
    S = Simulation(seedVal=30,numForward=200,numBackward=200, timestep = 1, potential=FLATLN, integrationScheme='forwardEuler',  fiducialPoints = [numpy.array([NEWTONG*GALMASS*(10**2),0,0,0,.1,0])])# fiducialPoints = [numpy.array([10,0,0,0,numpy.sqrt(V2POTENTIAL),0])],isData=isData,hasNoise=hasNoise, numRemaining=numRemaining)
    S.plotStream(0)
    S.plotH(0)
#S = Simulation(seedVal = 2000, numForward=10000, numBackward=10000, timestep=.1,potential=KEPLER)

def addStarTest():
    S = Simulation(seedVal=30,numForward=10000,numBackward=10000, timestep = 1, potential=FLATLN,integrationScheme='leapFrog')#, fiducialPoints = [numpy.array([10,0,0,0,numpy.sqrt(V2POTENTIAL),0])])#, fiducialPoints = [numpy.array([NEWTONG*GALMASS*(10**2),0,0,0,.1,0])])
    S.plotStream(0,1)
    print "0 pos", S.streams[0].stars[0].pos
    S.streams[0].addStars(0,100)
    S.plotStream(0,2)
    print "0 pos", S.streams[0].stars[0].pos

def sampleRateTest():
    S = Simulation(seedVal=20,numForward=100,numBackward=100, timestep = 1, potential=FLATLN,integrationScheme='leapFrog')#, fiducialPoints = [numpy.array([10,0,0,0,numpy.sqrt(V2POTENTIAL),0])])#, fiducialPoints = [numpy.array([NEWTONG*GALMASS*(10**2),0,0,0,.1,0])])
    S.plotStream(0,3)
    S.streams[0].changeSampleRate(.1)
    S.plotStream(0,4)

def format_plot(ax,i,j,dimension=6,label=0):
    units = [['kpc','kpc','kpc','kpc/myr','kpc/myr','kpc/myr'],['kpc','arcsec','arcsec','km/s','mas/yr','mas/yr']]
    if dimension==6:
        if j == 0:
#            plt.ticklabel_format(style='sci',scilimits=(-3,3),axis='y')

            plt.ylabel(('$x_%i$ ' % (i+1)) + units[label][i] ,fontsize=40,weight='bold')
#            plt.tick_params(axis='y',pad=15)
            plt.setp(ax[i][j].get_yticklabels(),fontsize=20,weight='bold')
        else:
            plt.setp(ax[i][j].get_yticklabels(),visible=False) 
        if i == 5:
#            plt.ticklabel_format(style='sci',scilimits=(-2,2),axis='x')
#            plt.tick_params(axis='x',pad=15)
            plt.xlabel(('$x_%i$ ' % (j+1))+units[label][j],fontsize=40,weight='bold')
            plt.setp(ax[i][j].get_xticklabels(),fontsize=20,weight='bold')
        else:
            plt.setp(ax[i][j].get_xticklabels(),visible=False)
        ax[i][j].locator_params(nbins=4)#,integer=True)
#        ax[i][j].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(3))
#        ax[i][j].yaxis.set_major_locator(matplotlib.ticker.LinearLocator(3))#locator_params(tight=True,nbins=4)            
        plt.ticklabel_format(useOffset=False)

    elif dimension==4:
        print 'not yet'
        pass
    
def EllipseTest(num_points,s,label):
    #s = Star([-10,0,20,0,.205,0])
    observer = numpy.array([-10,0,0,0,.205,0])
    s.changeCovarianceMatrix(observer,ang_err = 5*10**(-7))
    print "Covariance Matrix: \n",s.covarianceMatrix
    sample_points  = numpy.random.multivariate_normal(s.pos,s.covarianceMatrix,(num_points))





        
    fig = plt.figure(figsize=(15,15))
    fig.subplots_adjust(wspace=.2,hspace=.2)
    
#    fig.suptitle("Ellipse of covariance tensor for usual noise", fontsize=20)
    ax = [[],[],[],[],[],[]]
    
    print sample_points[:,0]
    for i in range(6):
        for j in range(6):
#             if i == 0:
#                 ax[i].append( fig.add_subplot(6,6,j+1+i*6) )
#             elif i != 0 and j ==0:
#                 ax[i].append(fig.add_subplot(6,6,j+1+i*6,sharex=ax[0][0]))
#             else:
#                 ax[i].append(fig.add_subplot(6,6,j+1+i*6,sharex=ax[0][j],sharey=ax[i][0]))

            ax[i].append(fig.add_subplot(6,6,j+1+i*6))
            plt.plot(sample_points[:,j],sample_points[:,i],'k.')#, alpha=.01)
#            plt.plot([observer[j]],[observer[i]],'rx')
            plt.axhline(observer[i])
            plt.axvline(observer[j])

            if j < 3:
                xmin = s.pos[j] - 10
                xmax = s.pos[j] + 10
            else:
                xmin = s.pos[j] - 1
                xmax = s.pos[j] + 1
            if i < 3:
                ymin = s.pos[i] - 10
                ymax = s.pos[i] + 10
            else:
                ymin = s.pos[i] - 1
                ymax = s.pos[i] + 1
            
#             xmin = numpy.amin(sample_points[:,j])
#             xmax = numpy.amax(sample_points[:,j])
#             ymin = numpy.amin(sample_points[:,i])
#             ymax = numpy.amax(sample_points[:,i])

#            xmin -= numpy.abs(xmin*.01)
#            xmax += numpy.abs(xmax*.01)
#            ymin -= numpy.abs(ymin*.01)
#            ymax += numpy.abs(ymax*.01)

                
#            std = numpy.std(sample_points[:,j])
#            xmin = numpy.amin(sample_points[:,j])
#            xmin -= 3*std
#            xmax = numpy.amax(sample_points[:,j])
#            xmax += 3*std
        
#             std = numpy.std(sample_points[:,i])
#             ymin = numpy.amin(sample_points[:,i])
#             ymin -= 3*std
#             ymax = numpy.amax(sample_points[:,i])
#            ymax += 3*std
            plt.axis(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
            format_plot(ax,i,j,6)

    print len(sample_points)
    print len(sample_points.tolist())
    for i in range(6):
        print "sigma %i: %.8f" % (i,numpy.std(sample_points[:,i]))
    plt.savefig("EllipseTest"+label)
    plt.close()
    Tensor_test(s,observer,num_points,label)

def Tensor_test(star,observerPos,num_points,label):
    ang_err=5*10**(-7) #rad
    prop_mot_err=2.0*(10**-3)*(1.0/3600.0)*(numpy.pi/180.0)*10**6 #rad/myr
    rad_vel_err= 10*SECPERMYR/KMPERKPC #kpc/myr
    ln_rad_pos_err=.1
    pos = star.pos
    relpos = star.pos-observerPos
    relpos_sphr = cartesian_to_spherical(relpos)
    print relpos_sphr.shape
    r = relpos_sphr[0]
    noise_tensor = numpy.diag([r*ln_rad_pos_err,ang_err,ang_err,rad_vel_err,prop_mot_err,prop_mot_err])
    print 'nt: ', noise_tensor
    noise_tensor = numpy.dot(noise_tensor,noise_tensor)
    print 'nt: ', noise_tensor

    points = numpy.random.multivariate_normal(relpos_sphr,noise_tensor,num_points)
    diff = []
    units_conversion_sphr = numpy.array([1,10**(-3)*MASPERRAD,10**(-3)*MASPERRAD,KMPERKPC/SECPERMYR,MASPERRAD*10**(-6),MASPERRAD*10**(-6)])
    for point in points:
#        diff.append(cartesian_to_spherical(spherical_to_cartesian(point)-relpos)*units_conversion_sphr)
        diff.append(point*units_conversion_sphr)

    diff = numpy.array(diff)
    fig = plt.figure(figsize=(15,15))
    fig.subplots_adjust(wspace=.2,hspace=.2)

    relpos_sphr *= units_conversion_sphr

    fig.suptitle("("+str(r)+", "+str(relpos_sphr[1])+", "+str(relpos_sphr[2])+", " + str(relpos_sphr[3]) +", "+ str(relpos_sphr[4])  + ", "+str(relpos_sphr[5])+")")

#    diff[:,4] /= diff[:,0]
#    diff[:,5] /= diff[:,0]
    ax = [[],[],[],[],[],[]]
    
    for i in range(6):
        for j in range(6):
            ax[i].append(fig.add_subplot(6,6,j+1+i*6))
            plt.plot(diff[:,j],diff[:,i],'k.', alpha=.01)

            format_plot(ax,i,j,6,label=1)

 
    plt.savefig("tensorEllipseTest"+label)
    plt.close()
    raise Exception("doneo")


def spherical_to_cartesian(pos):
    r = pos[0]
    theta = pos[1]
    phi = pos[2]
    x = r*numpy.cos(phi)*numpy.sin(theta)
    y = r*numpy.sin(phi)*numpy.sin(theta)
    z = r*numpy.cos(theta)

    y_hat = numpy.array([numpy.sin(theta)*numpy.sin(phi),numpy.cos(theta)*numpy.sin(phi),numpy.cos(phi)])
    
    z_hat = numpy.array([numpy.cos(theta),-numpy.sin(theta),0])
    x_hat = numpy.cross(y_hat,z_hat)

    vel = pos[3:].copy()
    vel[1:] *= r
    v_x = numpy.dot(vel,x_hat)
    v_y = numpy.dot(vel,y_hat)
    v_z = numpy.dot(vel,z_hat)
    answer = numpy.array([x,y,z,v_x,v_y,v_z])
    return answer

def cartesian_to_spherical(pos):
    r_hat = pos[:3].copy()

    distance = numpy.linalg.norm(r_hat)
    r_hat /= distance
 
    #positions in spherical coordinates

    r = distance
    theta = numpy.arccos(pos[2]/r)
    phi = numpy.arctan2(pos[1],pos[0])

   #angle_hat
    phi_hat = numpy.array([-numpy.sin(phi),numpy.cos(phi),0]) 
    theta_hat = numpy.cross(phi_hat,r_hat)

    vel = pos[3:].copy()
    v_r = numpy.dot(vel,r_hat)
    v_theta = numpy.dot(vel,theta_hat)
    v_phi = numpy.dot(vel,phi_hat)
    answer= numpy.array([r,theta, phi, v_r, v_theta/r,v_phi/r])
    

    return answer




def test_polar():
    x = numpy.random.randn(6)
    v1 = cartesian_to_spherical(x)
    v2 = spherical_to_cartesian(v1)

    print 'x',x
    print 'v1',v1
    print 'v2',v2
    for i in range(10**4):
        x = numpy.random.randn(6)
        v1 = cartesian_to_spherical(x)
        v2 = spherical_to_cartesian(v1)
        if (numpy.abs(x-v2) >= (10**(-10))).any():
            print x
            print v2
            raise Exception("I...Hate...Everything...")
        

def Diagnostic():
    test_polar()
    EnergyConservationTest(False,False,-1)
    addStarTest()
    sampleRateTest()
    EllipseTest(10000)
if __name__=="__main__":


    Diagnostic()
####mencoder "mf://*.png" -mf fps=10 -o 3str18.mp4 -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=800 
