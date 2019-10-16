# comment
import numpy
import matplotlib
import matplotlib.pyplot as plt
import time
import math
DARK = -9999999999
GM = 10000
Q = 1000
DIRGAUSSNORM = 20.0
FIDGAUSSNORM = 100.0

class ColdStream:

# Pick a velocity from a 3-d gaussian
# Pick a 3 space fiducial point from a 3-d gaussian
# Generate number of stars specified
    def __init__(self,size,id):
        
        self.dirn = numpy.random.randn(3)*DIRGAUSSNORM
        self.vel = numpy.linalg.norm(self.dirn)
        self.dirn = self.dirn/self.vel
        self.fid = numpy.random.randn(3)*FIDGAUSSNORM

        self.id = id
        self.stars = []


        for i in range(size):
            self.stars.append(Star(self.dirn,self.vel, self.fid, 1, i,id))

    def __init__(self,numstars,id,galacticCenter,numruns,rewoundsteps):
        self.numstars = numstars
        self.numruns = numruns
        self.rewoundsteps = rewoundsteps
        self.id = id
        self.pos = numpy.random.randn(3)*FIDGAUSSNORM
        

        self.galacticCenter = galacticCenter

        self.dirn = numpy.random.randn(3)*DIRGAUSSNORM
        self.vel = numpy.linalg.norm(self.dirn)
        self.dirn = self.dirn/self.vel

        self.times = numpy.random.random_integers(-1*rewoundsteps+1,0,numstars)
        self.unusedstars = []
        self.stars = []
        for i in range(numstars):
            s = Star(self.dirn,self.vel,self.pos,1,i,self.id,1)
            self.unusedstars.append(s)
        

    def __str__(self):
        if Q==0:
            return "I am a coldstream"
        else:
            return "I am a cluster"

    def setDist(self,dist):
        self.dist = dist
    
    def reId(self,newid,numstars):
        self.id = newid
#        print self.id
        if len(self.stars) != 0:
            for i in self.stars:
                i.id = newid*numstars+i.id
                i.ownerid = newid

        else:
            for i in self.unusedstars:
                i.id = newid*numstars+i.id
                i.ownerid = newid

    def update(self,dt,perturber):
        alargest = 0
        ### updates then returns tuple (largest acceleration kick, star id,distance)  
        alargest = -1
        largestid = 0
        distance=0
        for i in self.stars:
            r = perturber.pos - i.pos
            normr = numpy.linalg.norm(r)
            i.distances.append(normr)
            a = GM*r/(normr**3)
            i.update(dt,a)
            acc = numpy.linalg.norm(a)
            if acc > alargest:
                alargest = acc
                largestid = i.id
                distance = normr
        return (alargest,largestid,distance)

    def update(self,run,dt,perturber):
        toDelete = []
        for i in range(len(self.times)):

            if self.times[i] == run:
                self.stars.append(self.unusedstars.pop())
                toDelete.append(i)

        
        self.times = numpy.delete(self.times,toDelete)

        for i in self.stars:
            r = self.galacticCenter - i.pos
            normr = numpy.linalg.norm(r)
            aQ = Q*r/(normr**2*i.mass)
            r2 = perturber.pos - i.pos
            normr2 = numpy.linalg.norm(r2)
            a = GM*r2/(normr2**3) 
            if run < self.numruns - self.rewoundsteps:
                a = 0

            i.updateInPotential(dt,aQ,a,normr2)

    def rewind(self,run,dt):
        toDelete = []
        for i in range(len(self.times)):
            if self.times[i] == -run:
                self.stars.append(self.unusedstars.pop())
                toDelete.append(i)
#                if self.id == 0:
 #                   print "new star:"+ str(run)
        self.times = numpy.delete(self.times,toDelete)
        for i in self.stars:
            r = self.galacticCenter -i.pos
            normr = numpy.linalg.norm(r)
            aQ = Q*r/normr**2
            i.rewindInPotential(dt,aQ,0)

    def largePerturbed(self):
        plot = [[],[],[],[],[],[]]
        for i in self.stars:
            if numpy.linalg.norm(i.pos-i.image.pos)>25:
                plot[0].append(i.pos[0])
                plot[1].append(i.pos[1])
                plot[2].append(i.pos[2])
                plot[3].append(i.image.pos[0])
                plot[4].append(i.image.pos[1])
                plot[5].append(i.image.pos[2])

        return plot

    def plot(self):
        plot = [[],[],[],[],[],[],[]] # first three are for regular stars and second three are for image stars and the last is the id of the stars that need arrows generated
        for i in self.stars:
            plot[0].append(i.pos[0])
            plot[1].append(i.pos[1])
            plot[2].append(i.pos[2])
            plot[3].append(i.image.pos[0])
            plot[4].append(i.image.pos[1])
            plot[5].append(i.image.pos[2])
            if i.showArrow ==1:
                plot[6].append(i.id)
        return plot



        
class Star:
# Give it mass, vel,id, owner id and dirn from stream
# give it pos from fid*a multiplier taken from 0 mean 200 s.d. gaussian
# create an imageStar to go along with it
    def __init__(self, dirn, vel,fid, mass,id,ownerid,potentialExists):
        self.ownerid = ownerid
        self.id = id
        self.mass = mass
        self.vel = vel*dirn
        self.streamvel = self.vel
        self.dirn = dirn
        self.showArrow = 0
        self.closestApproach = 2**15
        if(id==DARK):
            self.poswrtfid=fid
        else:
            self.poswrtfid = numpy.random.normal(0,2000)
        if potentialExists:
            self.pos = fid
        else:
            self.pos = fid + dirn * self.poswrtfid
        
        self.image =ImageStar(dirn,vel,fid,mass,id,ownerid,self.poswrtfid,potentialExists)
        self.distances = []

    def __str__(self):
        s =  "velocity: "+ str(self.vel) + "\nDirection: " + str(self.dirn) + "\nPosition: " + str(self.pos) +'\n'
        return s

    def update(self, dt, a,dist):
        self.vel = self.vel + a*dt
        self.pos = self.pos + self.vel*dt
        self.dirn = self.vel/numpy.linalg.norm(self.vel)
        self.image.update(dt)
        if numpy.linalg.norm(self.pos - self.image.pos) > 10 and self.showArrow != 1:
            self.showArrow = 1
        if self.closestApproach > dist:
            self.closestApproach = dist

    def updateInPotential(self,dt,aQ,a,dist):
        self.vel = self.vel+(aQ+a)*dt
        self.pos = self.pos+self.vel*dt
        self.dirn = self.vel/numpy.linalg.norm(self.vel)
        self.image.updateInPotential(dt,aQ)
        self.distances.append(dist)
        #print numpy.linalg.norm(self.pos-self.image.pos)
        if numpy.linalg.norm(self.pos - self.image.pos) > 10:
            self.showArrow = 1
        if self.closestApproach > dist:
            self.closestApproach = dist
    
    def rewindInPotential(self,dt,aQ,a):
        self.vel = self.vel + (aQ+a)*dt
        self.pos = self.pos-self.vel*dt
        self.dirn = self.vel/numpy.linalg.norm(self.vel)
        self.image.rewindInPotential(dt,aQ)
        if numpy.linalg.norm(self.pos-self.image.pos) >10 and self.showArrow != 1:
            self.showArrow = 1

    def minDistance(self):
        d = self.distances[0]
        for i in self.distances:
            if d > i:
                d = i
        return d



class ImageStar(Star):
##image of a star which will show up if there is a significant velocity kick given to it
    def __init__(self,dirn,vel,fid,mass,id,ownerid,poswrtfid,potentialExists):
        self.ownerid = ownerid
        self.id = id
        self.mass = mass
        self.vel = vel*dirn
        self.dirn = dirn
        self.poswrtfid = poswrtfid
        if potentialExists:
            self.pos = fid
        else:
            self.pos = fid + dirn*poswrtfid
        
    def update(self,dt):
        self.pos = self.pos+self.vel*dt
    def updateInPotential(self,dt, aQ):
        self.vel = self.vel + aQ*dt
        self.pos = self.pos+ self.vel*dt
    
    def rewindInPotential(self,dt,aQ):
        self.vel = self.vel + aQ*dt
        self.pos = self.pos - self.vel*dt
        


class Simulation:
    def __init__(self, numstreams, numstars, numruns, timestep,seedval,jump,rewoundsteps,starsOn,arrowsOn):
        self.starsOn = starsOn
        self.arrowsOn = arrowsOn
        self.jump = jump
        self.numstreams = numstreams
        self.numstars = numstars
        self.numruns = numruns
        self.timestep = timestep
        numpy.random.seed(seedval)
        self.perturber = Star(numpy.array([1,0,0]),100,numpy.array([-250,0,0]),1,DARK, DARK,Q)        
        self.streams = []
        self.rewoundsteps = rewoundsteps


        if Q == 0:
            self.initNoPotSim()
        else:
            self.initPotSim()



    def initPotSim(self):
        self.galacticCenter = [-500,-500,-500]
        for i in range(self.numstreams):
###            self.streams.append(Cluster(self.numstars,i,self.galacticCenter, self.numruns,self.rewoundsteps))
            self.streams.append(ColdStream(self.numstars,i,self.galacticCenter, self.numruns,self.rewoundsteps))
        self.streams = self.quickSortVel(self.streams)
        for i in range(len(self.streams)):
            self.streams[i].reId(i,self.numstars)
            distance = numpy.abs(numpy.dot(numpy.cross(self.perturber.dirn,self.streams[i].dirn), (self.perturber.pos-self.streams[i].unusedstars[0].pos)))
            self.streams[i].dist = distance
            #print self.streams[i]

        #print self.streams[0].times

        for run in range(self.rewoundsteps):
            self.rewindInPotential(run)

        for i in self.streams:
            l=[]
            for j in i.stars:
                l.append(j.id)
            #print len(l)
            for j in i.unusedstars:
                print 'SOMETHING IS WRONG'
                raise Exception("UNUSEDSTAR")
#        raise Exception('spam')

    def quickSortVel(self,list):
        if len(list) <= 1:
            return list
        pivot = list.pop(0)
        less = []
        greater = []
        for i in list:
            if i.vel <= pivot.vel:
                less.append(i)
            else:
                greater.append(i)
        ret = []
        ret.extend(self.quickSortVel(less))
        ret.append(pivot)
        ret.extend(self.quickSortVel(greater))
        return ret

    def initNoPotSim(self):
        ###create streams
        closest = 2**15
        
        for i in range(self.numstreams):
            self.streams.append(ColdStream(self.numstars,i))
            distance = numpy.abs(numpy.dot(numpy.cross(self.perturber.dirn,self.streams[i].dirn), (self.perturber.pos-self.streams[i].stars[0].pos)))
            self.streams[i].setDist(distance)
        self.streams = self.quickSortDist(self.streams)
        for i in range(len(self.streams)):
            self.streams[i].reId(i,self.numstars)
            for j in self.streams[i].stars:
                j.closestApproach = numpy.linalg.norm(self.perturber.pos - j.pos)


            
###sorting streams by distances
    def quickSortDist(self, list):
        if len(list) <= 1:
            return list
        pivot = list.pop(0)
        closer = []
        farther = []
        for i in list:
            if i.dist <= pivot.dist:
                closer.append(i)
            else:
                farther.append(i)
        ret = []
        ret.extend(self.quickSortDist(closer))
        ret.append(pivot)
        ret.extend(self.quickSortDist(farther))
        return ret

    def updateNoPotential(self):

        self.perturber.update(self.timestep,0,0)
        alargest = (0,0,0)
        for i in self.streams:
           a = i.update(self.timestep,self.perturber)
           if a[0] > alargest[0]:
               alargest = a
        return alargest

    def updatePotential(self,run):
        r = self.galacticCenter - self.perturber.pos
        normr = numpy.linalg.norm(r)
        a = Q*r/normr**2

        self.perturber.update(self.timestep,a,0)

        for i in self.streams:
            i.update(run,self.timestep,self.perturber)

    def rewindInPotential(self,run):
        #print run
        if run > math.floor(self.rewoundsteps/2):
            r = self.galacticCenter - self.perturber.pos
            normr = numpy.linalg.norm(r)
            a = Q*r/normr**2
            self.perturber.rewindInPotential(self.timestep,a,0)
        for i in self.streams:
            i.rewind(run,self.timestep)
    
    def distances(self):
        f = open("distances",'w')
        for i in self.streams:
            s = "" + str(i.id) + ":"+ str(i.dist) + "\n"
            f.write(s)
        f.close()


    def FunctionalTest1(self):
        for run in range(self.numruns):
            self.perturber.update(self.timestep,0)
            for i in self.streams:
                i.update(self.timestep, self.perturber)

        toPlot =[[],[]]

        for stream in self.streams:
            b = stream.dist
            v = numpy.linalg.norm(self.perturber.vel - stream.vel)
            for star in stream.stars:
                d = star.minDistance()
                deltaV = numpy.linalg.norm(star.vel-star.image.vel)
                toPlot[0].append(d/b)
                toPlot[1].append(deltaV*b*v/(2*GM))

        plt.scatter(toPlot[0],toPlot[1])
        plt.axis([0,10,0,1.5])
        plt.xlabel("d/b")
        plt.ylabel("deltaV/(2*GM/bv)")
        plt.savefig("FunctionalTest.png")
        plt.clf()

### IMPORTANT THIS TEST IS DEPENDENT ON THE TIME STEP SIZE: E.G. IF TIMESTEP = 10^-1, ROUNDING SHOULD BE TO 1

    def FunctionalTest2(self):
        f = open("functionaltest2.out",'w')
        tmins=[]
        for stream in self.streams:
            for star in stream.stars:
                dx = self.perturber.pos - star.pos
                dv = self.perturber.vel - star.vel
                t = -1*numpy.round(numpy.dot(dx,dv)/numpy.dot(dv,dv),1)
                tmins.append(t)

        for run in range(self.numruns):
            self.perturber.update(self.timestep,0)
            for stream in self.streams:
                for star in stream.stars:
                    if run*self.timestep == tmins[star.id]:
                        deltav = 2*GM*(star.pos-self.perturber.pos)/(numpy.linalg.norm(self.perturber.vel-star.vel)*(numpy.linalg.norm(self.perturber.pos-star.vel)**2))
                        star.update(1,deltav)
                    else:
                        star.update(self.timestep,0)


        for stream in self.streams:
            for star in stream.stars:
                dv = numpy.linalg.norm(star.vel - star.image.vel)
                dx = numpy.linalg.norm(star.pos - star.image.pos)
                s = str(star.id)+':'+str(dv)+','+str(dx)+':'+str(tmins[star.id])+'\n'
                f.write(s)
        f.close()


        

    def runNoPotentialNoPlot(self):
        f = open("largestkicks",'w')
        alargest = (0,0,0)
        largestkicks = []
        for run in range(self.numruns): 
            a = self.updateNoPotential()
            f.write(str(a)+'\n')
            if a[0] > alargest[0]:
                alargest = a
                largestkicks.append((run*self.timestep,a))
        f.write(str(largestkicks))
        f.close()

    def runNoPotentialNoPlot2(self):
        f=open("runnoplot2.out",'w')
        for run in range(self.numruns):
            self.updateNoPotential()
        for stream in self.streams:
            for star in stream.stars:
                dv = numpy.linalg.norm(star.vel - star.image.vel)
                dx = numpy.linalg.norm(star.pos - star.image.pos)
                s = str(star.id)+':'+str(dv) + ','+str(dx) + '\n'
                f.write(s)
        f.close()

    def runNoPotential(self,numclosest):
        f = open("largestkicks",'w')
        alargest = 0
        largestkicks = []
        for run in range(self.numruns):
            if run%self.jump == 0:
                for j in range(numclosest):
                    self.plotStream(j,run)
                self.plot(run,1,0)

            a = self.updateNoPotential()
            if a > alargest:
                alargest = a
                largestkicks.append(run)
        f.write(str(largestkicks))
        f.close()

    def runPotential(self):
        for run in range(self.rewoundsteps*2):
            self.updatePotential(run)
            #print run
            if run % self.jump == 0:
                if self.starsOn or self.arrowsOn:
                    self.plot(run,self.starsOn,self.arrowsOn)
#                    self.plotMostPerturbed(run,self.starsOn,self.arrowsOn)
#                self.plotStream(1,run)
        self.plotClosestApproaches()


    def plotMostPerturbed(self,run,plotStars,plotArrows):
        toPlot = [[],[],[],[],[],[],[]]

        for i in self.streams:
            ret = i.largePerturbed()
            for j in range(6):
                toPlot[j].extend(ret[j])

        if len(toPlot[0]) == 0:
            #print "yes"
            return

        #### Plot XY: all
        if plotStars ==1:
            plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
#            plt.scatter(toPlot[3],toPlot[4],c='b', alpha= 0.25, edgecolors = 'None')

        if plotArrows ==1:
            t = time.time()
            for i in toPlot[6]:
                plt.arrow(toPlot[3][i],toPlot[4][i],toPlot[0][i]-toPlot[3][i],toPlot[1][i]-toPlot[4][i],ls='solid')
            print "time to make arrows: " + str(t-time.time())
        plt.scatter([self.perturber.pos[0]],[self.perturber.pos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.axis([-500,500,-500,500])
        


        zerosinname = int(math.floor(math.log10(self.numruns)))
        filename = 'a'
        if run < 0:
            filename +='.v.'
            run = self.rewoundsteps+run
        else:
            filename += '.w.'
        for i in range(zerosinname+1):
            if run < 10**i:
                for j in range(zerosinname+1-i):
                    filename += '0'
                filename += str(run) +'.png'
                break
        
        if run >= 10**zerosinname:
            filename += str(run) + '.png'


        plt.xlabel("x")
        plt.ylabel('y')
        plt.savefig(filename)
        plt.clf()

        #### PLOT Z X: ALL
        if plotStars==1:
            plt.scatter(toPlot[2],toPlot[0], c='b',alpha=0.5,edgecolors='None')
#            plt.scatter(toPlot[5],toPlot[3],c='b',alpha=0.25, edgecolors='None')
        if plotArrows==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[5][i],toPlot[3][i],toPlot[2][i]-toPlot[5][i],toPlot[0][i]-toPlot[3][i],ls='solid')

        plt.scatter([self.perturber.pos[2]],[self.perturber.pos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)           
            
        
                    
        plt.axis([-500,500,-500,500])

        filename = 'b' + filename[1:]


        plt.xlabel("z")
        plt.ylabel('x')
        plt.savefig(filename)
        plt.clf()

        #### PLOT Y Z: ALL
        if plotStars==1:
            plt.scatter(toPlot[1],toPlot[2], c='b',alpha=0.5, edgecolors='None')
#            plt.scatter(toPlot[4],toPlot[5],c='b',alpha=0.25, edgecolors='None')
        if plotArrows==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[4][i],toPlot[5][i],toPlot[1][i]-toPlot[4][i],toPlot[2][i]-toPlot[5][i],ls='solid')
        plt.scatter([self.perturber.pos[1]],[self.perturber.pos[2]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.axis([-500,500,-500,500])

        filename = 'c' + filename[1:]


        plt.xlabel("y")
        plt.ylabel('z')
        plt.savefig(filename)
        plt.clf()

    
    def plot(self, run,plotStars , plotArrows ):
        toPlot = [[],[],[],[],[],[],[]]


        for i in self.streams:
            ret = i.plot()
            for j in range(7):
                toPlot[j].extend(ret[j])
        #print toPlot[6]

        if len(toPlot[0]) == 0:
            #print "yes"
            return

        #### Plot XY: all
        if plotStars ==1:
            plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
#            plt.scatter(toPlot[3],toPlot[4],c='b', alpha= 0.25, edgecolors = 'None')

        if plotArrows ==1:
            t = time.time()
            for i in toPlot[6]:
                plt.arrow(toPlot[3][i],toPlot[4][i],toPlot[0][i]-toPlot[3][i],toPlot[1][i]-toPlot[4][i],ls='solid')
            print "time to make arrows: " + str(t-time.time())
        plt.scatter([self.perturber.pos[0]],[self.perturber.pos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.axis([-500,500,-500,500])
        


        zerosinname = int(math.floor(math.log10(self.numruns)))
        filename = 'a'
        if run < 0:
            filename +='.v.'
            run = self.rewoundsteps+run
        else:
            filename += '.w.'
        for i in range(zerosinname+1):
            if run < 10**i:
                for j in range(zerosinname+1-i):
                    filename += '0'
                filename += str(run) +'.png'
                break
        
        if run >= 10**zerosinname:
            filename += str(run) + '.png'


        plt.xlabel("x")
        plt.ylabel('y')
        plt.savefig(filename)
        plt.clf()

        #### PLOT Z X: ALL
        if plotStars==1:
            plt.scatter(toPlot[2],toPlot[0], c='b',alpha=0.5,edgecolors='None')
#            plt.scatter(toPlot[5],toPlot[3],c='b',alpha=0.25, edgecolors='None')
        if plotArrows==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[5][i],toPlot[3][i],toPlot[2][i]-toPlot[5][i],toPlot[0][i]-toPlot[3][i],ls='solid')

        plt.scatter([self.perturber.pos[2]],[self.perturber.pos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)           
            
        
                    
        plt.axis([-500,500,-500,500])

        filename = 'b' + filename[1:]


        plt.xlabel("z")
        plt.ylabel('x')
        plt.savefig(filename)
        plt.clf()

        #### PLOT Y Z: ALL
        if plotStars==1:
            plt.scatter(toPlot[1],toPlot[2], c='b',alpha=0.5, edgecolors='None')
#            plt.scatter(toPlot[4],toPlot[5],c='b',alpha=0.25, edgecolors='None')
        if plotArrows==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[4][i],toPlot[5][i],toPlot[1][i]-toPlot[4][i],toPlot[2][i]-toPlot[5][i],ls='solid')
        plt.scatter([self.perturber.pos[1]],[self.perturber.pos[2]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.axis([-500,500,-500,500])

        filename = 'c' + filename[1:]


        plt.xlabel("y")
        plt.ylabel('z')
        plt.savefig(filename)
        plt.clf()

    def plotStream(self,streamnum,run):
        toPlot = [[],[],[],[],[],[],[]]
        stream = self.streams[streamnum]
        ret = stream.plot()
        if len(ret[0]) == 0:
            return
       
        for i in range(7):
            toPlot[i].extend(ret[i])
        
        #### PLOT XY
        plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([self.perturber.pos[0]],[self.perturber.pos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter(toPlot[3],toPlot[4],c='b',alpha=0.25,edgecolors  = 'None')
        for i in toPlot[6]:
            plt.arrow(toPlot[3][i],toPlot[4][i],toPlot[0][i]-toPlot[3][i],toPlot[1][i]-toPlot[4][i],ls='solid')

        plt.axis([-500,500,-500,500])
        if run < 10:
            filename = 'Stream:'+str(stream.id)+":0"+str(run)+':0.png'
        else:
            filename = 'Stream:'+str(stream.id)+":"+str(run) + ':0.png'
        plt.xlabel("x")
        plt.ylabel('y')
        plt.savefig(filename)
        plt.clf()
        
        #### PLOT YZ
        plt.scatter(toPlot[1],toPlot[2],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([self.perturber.pos[1]],[self.perturber.pos[2]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter(toPlot[4],toPlot[5],c='b',alpha=0.25,edgecolors  = 'None')
        for i in range(len(toPlot[0])):
            plt.arrow(toPlot[4][i],toPlot[5][i],toPlot[1][i]-toPlot[4][i],toPlot[2][i]-toPlot[5][i],ls='solid')

        
        plt.axis([-500,500,-500,500])
        if run < 10:
            filename = 'Stream:'+str(stream.id)+":0"+str(run)+':1.png'
        else:
            filename = 'Stream:'+str(stream.id)+":"+str(run) + ':1.png'
        plt.xlabel("y")
        plt.ylabel('z')
        plt.savefig(filename)
        plt.clf()
        
        #### PLOT ZX
        plt.scatter(toPlot[2],toPlot[0],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([self.perturber.pos[2]],[self.perturber.pos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter(toPlot[5],toPlot[3],c='b',alpha=0.25,edgecolors  = 'None')
        for i in range(len(toPlot[0])):
            plt.arrow(toPlot[5][i],toPlot[3][i],toPlot[2][i]-toPlot[5][i],toPlot[0][i]-toPlot[3][i],ls='solid')
        
        plt.axis([-500,500,-500,500])
        if run < 10:
            filename = 'Stream:'+str(stream.id)+":0"+str(run)+':2.png'
        else:
            filename = 'Stream:'+str(stream.id)+":"+str(run) + ':2.png'
        plt.xlabel("Z")
        plt.ylabel('x')
        plt.savefig(filename)
        plt.clf()

    def plotClosestApproaches(self):
        toPlotX = []
        toPlotY = []

        for stream in self.streams:
           b = stream.dist
           v = numpy.linalg.norm(self.perturber.vel - stream.vel)
           for star in stream.stars:
                d = star.closestApproach
                deltaV = numpy.linalg.norm(star.vel-star.image.vel)
                toPlotX.append(star.closestApproach)
                toPlotY.append(numpy.linalg.norm(star.vel - star.image.vel))

        plt.loglog(toPlotX,toPlotY,".")

        plt.xlabel("Closest Approach")
        plt.ylabel("Delta V")
        plt.savefig("closestapproachVdv.png")
        plt.clf()


# def comparator():
#     f = open('functionaltest2.out','r')
#     farray = f.readlines()
#     g = open('runnoplot2.out','r')
#     garray = g.readlines()
#     differs = 0
#     for i in range(len(farray)):
#         if farray[i] != garray[i]:
#             differs+=1
            
#     print differs
        
# def functestforpot():
#     S = Simulation(50,1,100*2,.1,20,5,100*2,1,1)
#     avg=0
#     for stream in S.streams:
#         avg +=numpy.linalg.norm(stream.pos - S.galacticCenter)
#     avg /= 50
#     print avg

#functestforpot()

#w = [1,0,0]
#r = [0,1,0]
#print numpy.cross(r,w)
### Simulation(numstreams,numstars,numruns,timestep,seedval,jump,rewoundsteps,starsOn,arrowsOn)
#S = Simulation(10000,20,500,100,.1,20)



#S = Simulation(10,100,100*1,.1,20,5,100*1,1,0)

#S = Simulation(50,1000,100*2,.1,20,5,100*2,1,1)





#S.run(1,1)
#S.FunctionalTest2()
#S.FunctionalTest1()
#S.run(1,1)
#S.FunctionalTest2()
#S.FunctionalTest()
#S.run(1,50)
#S.runNoPlot()
#S.runNoPlot2()
#S.distances()

#S.runPotential()




#comparator()

#c = Cluster(100,0,[-500,-500,-500],100*20,100*10)

#    def __init__(self,numstars,id,galacticCenter,numruns,rewoundsteps):

    




###plot(20,20,500,100,.1)

####mencoder "mf://*.png" -mf fps=10 -o 3str18.mp4 -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=800 
