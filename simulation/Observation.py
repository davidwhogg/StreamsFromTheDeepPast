import Simulator
from Simulator import numpy as np
from Simulator import os
from Simulator import math
### Simulation(numstreams,numstars,numruns,timestep,seedval,jump,rewoundsteps,simtype)

## the goal of this class is to be as an observer on a planet observing the galaxy
class Observer:
## create a simulation instance, take a random position in the galaxy 
## create a cartesion coordinate system based on this 
    def __init__(self):
        os.mkdir('data')
#        os.chdir('data')
        self.S = Simulator.Simulation(1, 100,10,10,20,10,100*10,3)
        os.chdir('data')
        self.pos = np.random.randn(3)*Simulator.POSGAUSSSTDEV #kpc
        self.xaxis = Simulator.GALACTICCENTER - self.pos # kpc
        self.xaxis /= np.linalg.norm(self.xaxis) # dimensionless
        self.yaxis = np.cross(self.xaxis,np.random.randn(3)) #kpc
        self.yaxis /= np.linalg.norm(self.yaxis) #dimensionless
        self.zaxis = np.cross(self.xaxis,self.yaxis) # dimensionless
        print self.pos #kpc
        
## this takes the coordinates of a star and turns it in to ra and dec, and returns [ra,dec,radial distance]
    def data(self,star):
        diff = star.pos-self.pos #kpc
        norm = np.linalg.norm(diff) #kpc
        l = np.pi/2 - np.arccos(np.dot(diff,self.zaxis)/(np.linalg.norm(diff))) # rad       
        b = np.arctan2(np.dot(diff,self.yaxis),np.dot(diff,self.xaxis)) # rad
        return [l,b,np.linalg.norm(diff)]

## takes the simulator through a step and creates data for an hammeraitoff projection
    def run(self,jump):
        for run in range(self.S.rewoundsteps*2):
            self.S.update(run)
            if run % jump == 0:
                self.formattedplot(run)
## takes the simulator through a step and
    def runPatch(self,runno,starid):
        for run in range(runno+1):
            self.S.update(run)
            if run == runno:
                self.patchPlot(runno,starid)
                break

## creates cartesian plot data for stars with respect to its coordinate system
    def plot(self,run):
        zerosinname = int(math.floor(math.log10(self.S.rewoundsteps*2)))
        filename = ''

        for i in range(zerosinname+1):
            if run < 10**i:
                for j in range(zerosinname+1-i):
                    filename += '0'
                filename += str(run)
                break
        
        if run >= 10**zerosinname:
            filename += str(run)
            

        file = open(filename,'w')
        file.write('ID\tl\tb\tDIST\n')
 
        for stream in self.S.streams:
            for star in stream.stars:
                file.write(str(star.id) + '\t'+str(star.pos[0]) +'\t'+str(star.pos[1])+'\t'+str(star.pos[2])+'\n')
        ret = self.data(self.S.perturber)
        file.write(str(self.S.perturber.id)+ '\t'+str(ret[0])+'\t'+str(ret[1])+'\t'+str(ret[2]))
        file.close()

## takes the data and prints out the data for a hammer-aitoff plot for run run
## output: [number] a correclty formatted data file
    def formattedplot(self,run):
        zerosinname = int(math.floor(math.log10(self.S.rewoundsteps*2)))
        filename = ''

        for i in range(zerosinname+1):
            if run < 10**i:
                for j in range(zerosinname+1-i):
                    filename += '0'
                filename += str(run)
                break
        
        if run >= 10**zerosinname:
            filename += str(run)
            

        file = open(filename,'w')
        file.write('ID\tl\tb\tDIST\n')
 
        for stream in self.S.streams:
            for star in stream.stars:
                ret = self.data(star)
                file.write(str(star.id) + '\t'+str(ret[0]) +'\t'+str(ret[1])+'\t'+str(ret[2])+'\n')
        ret = self.data(self.S.perturber)
        file.write(str(self.S.perturber.id)+ '\t'+str(ret[0])+'\t'+str(ret[1])+'\t'+str(ret[2]))
        file.close()

## creates the data for hammer-aitoff plots of the patch of sky around a certain star
    def patchPlot(self,run, starid):
        zerosinname = int(math.floor(math.log10(self.S.rewoundsteps*2)))
        filename = 'patchof:'+str(starid)+':'
        
        for i in range(zerosinname+1):
            if run < 10**i:
                for j in range(zerosinname+1-i):
                    filename += '0'
                filename += str(run)
                break
        
        if run >= 10**zerosinname:
            filename += str(run)
            

        file = open(filename,'w')
#        file.write('ID\tX\tY\n')
        print len(self.S.streams)
        star = self.S.streams[starid/100].stars[starid%100]
        radial = star.pos - self.pos #kpc
        radial[0] = np.dot(self.xaxis,star.pos - self.pos) #kpc
        radial[1] = np.dot(self.yaxis,star.pos - self.pos) #kpc
        radial[2] = np.dot(self.zaxis,star.pos - self.pos) #kpc
        radial /= np.linalg.norm(radial) #dimensionless

        xtp = np.cross(radial,self.zaxis) #dimensionless
        xtp /= np.linalg.norm(xtp) # dimensionless
        ytp = np.cross(radial,xtp) #dimensionless
        
        for stream in self.S.streams:
            for s in stream.stars:
                pos = s.pos-self.pos
                pos[0] = np.dot(self.xaxis,s.pos-self.pos) #kpc
                pos[1] = np.dot(self.yaxis,s.pos-self.pos) #kpc
                pos[2] = np.dot(self.zaxis,s.pos-self.pos) #kpc
                spos = pos/(np.dot(pos,radial)) #dimensionless
                if np.dot(pos,radial) <= 0:
                    pass
                else:
                    diff = spos - radial #dimensionless
                    file.write(str(s.id)+'\t'+str(np.dot(diff,xtp))+'\t'+str(np.dot(diff,ytp))+'\n')
        file.close()
        
## this function demonstrates the growth of a perturbation as a function of the mass of the perturber.
def runGrowth(runno,starid,startval,endval,step):
    earth = Observer()
    Simulator.PERTURBERMASS = startval
    while Simulator.PERTURBERMASS <= endval:
        dir = str(Simulator.PERTURBERMASS)
        os.mkdir(dir)
        os.chdir(dir)
        earth.runPatch(runno,starid)
        Simulator.PERTURBERMASS += step
        os.rename('patchof:'+str(starid)+':'+str(runno),'../'+str(Simulator.PERTURBERMASS)+'patchof:'+str(starid)+':'+str(runno))
        os.chdir('../')
        earth.S = Simulator.Simulation(10, 100,10,1,20,0,100*10,3)
    

earth = Observer()
#earth.runGrowth(1950,127,1e6,1e12,1e1)
earth.run(50)
#runGrowth(1050,127,1e6,1e12,10)

