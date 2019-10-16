import os
import matplotlib
import matplotlib.pyplot as plt
from numpy import abs

##This reads the specifically formatted data files that are output by Simulator.plot and turns it in to a list format as it was stored by Simulation
def readData(filename):
    print filename
    file = open(filename,'r')
    toPlot = [[],[],[],[],[],[],[]]
    lines = file.readlines()
    line = lines[1]
    perturberpos = line[1:-2].rsplit(',')
    for i in range(len(perturberpos)):
        perturberpos[i]=float(perturberpos[i])
    i = 2
#    print len(lines)
    while i < len(lines)-1:
        line = lines[i].rsplit('\t')

        for j in range(len(toPlot)-1):
            
            toPlot[j].append(float(line[j]))
        i+=1

    line = lines[i]
    line = line[1:-1].rsplit(',')
    if line[0] == '':
        line = []
    for i in range(len(line)):
        toPlot[6].append(int(i))
    return (perturberpos,toPlot)

# this plots all of the stream position data for the simulator from a cartesian view looking at the galaxy
# Arguments: plotStars if 1 plot stars if not dont, plotArrows, plot arrows between stars and image  if 1 , if not don't, ObsLoc the location in space of the observer
#output: [letter][run].png where number is a b or c, denoting xy yz zx plot respectively
def plotMainData(plotStars=1,plotArrows=0,ObsLoc=[0,0,0]):
    STARTDIR=os.getcwd()
    if not('OUT' in os.listdir(STARTDIR)):
        os.mkdir('OUT')
    PLOTSDIR=STARTDIR+'/OUT'
    if not('MAINDATA' in os.listdir(STARTDIR)):
        print 'Error: No Data Folder In this Directory'
        return 1

    os.chdir('MAINDATA')
    MAINDATADIR = os.getcwd()

    for fname in os.listdir(MAINDATADIR):
        print fname
        ret = readData(fname)
        perturberpos = ret[0]
        
        toPlot =ret[1]
        
        os.chdir(PLOTSDIR)

        if len(toPlot[0]) == 0:
            return

        #### Plot XY: all
        if plotStars ==1:
            #print toPlot[1]
            plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')



        if plotArrows ==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[3][i],toPlot[4][i],toPlot[0][i]-toPlot[3][i],toPlot[1][i]-toPlot[4][i],ls='solid')

        plt.scatter([ObsLoc[0]],[ObsLoc[1]],c='g',alpha=0.9,edgecolors='None',pickradius=1000.0)

        plt.scatter([perturberpos[0]],[perturberpos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)


        plt.axis([-500,500,-500,500])



        filename = 'a'+fname + '.png'

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

        plt.scatter([perturberpos[2]],[perturberpos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)           
        plt.scatter([ObsLoc[2]],[ObsLoc[0]],c='g',alpha=0.9,edgecolors='None',pickradius=1000.0)
        
                    
#        plt.axis([-200,200,-200,200]) #kpc
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
        plt.scatter([perturberpos[1]],[perturberpos[2]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter([ObsLoc[1]],[ObsLoc[2]],c='g',alpha=0.9,edgecolors='None',pickradius=1000.0)
#        plt.axis([-200,200,-200,200]) # kpc
        plt.axis([-500,500,-500,500])
        filename = 'c' + filename[1:]


        plt.xlabel("y")
        plt.ylabel('z')
        plt.savefig(filename)
        plt.clf()
        
        os.chdir(MAINDATADIR)
   


    os.chdir(STARTDIR)

#plots the total energies
def plotH():
    CWD = os.getcwd()
    if not('energies' in os.listdir(CWD)):
        return 'Error: no energy'
    f = open('energies','r')
    lines = f.readlines()
    toPlot = [[],[],[],[],[]]
    for line in lines:
        l = line[1:-2].rsplit(',')
        toPlot[0].append(float(l[0]))
        toPlot[1].append(float(l[1]))
        toPlot[2].append(float(l[2]))
        toPlot[3].append(float(l[3]))
        toPlot[4].append(float(l[4]))

#    plt.scatter(toPlot[0],toPlot[1],c='b')
#    plt.scatter(toPlot[0],toPlot[2],c='g')
    plt.scatter(toPlot[0],toPlot[3],c='r')
#    plt.scatter(toPlot[0],toPlot[4],c='y')
    plt.savefig("energies.png")
    plt.clf()

#plots only the most perturbed stars
#arguments: plotStars plot stars if 1 else dont. plotArrows draw arrows between star and image if 1 else dont.
#output: [letter]MostPerturbed[run].png where number is a b or c, denoting xy yz zx plot respectively
def plotMostPerturbed(plotStars,plotArrows):
    CWD = os.getcwd()
    if not('MOSTPERTURBED' in os.listdir(CWD)):
        return 'ERROR: NO MOSTPERTURBED DIRECTORY'
    os.chdir('MOSTPERTURBED')
    for fname in os.listdir(os.getcwd()):
        ret = readData(fname)
        toPlot = ret[1]
        perturberpos = ret[0]
        
        
        
        if len(toPlot[0]) == 0:
            #print "yes"
            return

        #### Plot XY: all
        if plotStars ==1:
            plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
#            plt.scatter(toPlot[3],toPlot[4],c='b', alpha= 0.25, edgecolors = 'None')

        if plotArrows ==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[3][i],toPlot[4][i],toPlot[0][i]-toPlot[3][i],toPlot[1][i]-toPlot[4][i],ls='solid')
        plt.scatter([perturberpos[0]],[perturberpos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.axis([-200,200,-200,200]) # kpc
        


        zerosinname = int(math.floor(math.log10(numruns)))
        filename = 'a'+fname+'.png'

        plt.xlabel("x")
        plt.ylabel('y')
        plt.savefig(filename)
        plt.clf()

        #### PLOT Z X: ALL
        if plotStars==1:
            plt.scatter(toPlot[2],toPlot[0], c='b',alpha=0.5,edgecolors='None')

        if plotArrows==1:
            for i in toPlot[6]:
                plt.arrow(toPlot[5][i],toPlot[3][i],toPlot[2][i]-toPlot[5][i],toPlot[0][i]-toPlot[3][i],ls='solid')

        plt.scatter([perturberpos[2]],[perturberpos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)           
            
        
                    
        plt.axis([-200,200,-200,200]) #kpc

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
        plt.scatter([perturberpos[1]],[perturberpos[2]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.axis([-200,200,-200,200])#kpc

        filename = 'c' + filename[1:]


        plt.xlabel("y")
        plt.ylabel('z')
        plt.savefig(filename)
        plt.clf()




#plot a single stream
#arguments streamnum: stream to plot, plotStars: as previous, plotArrows: as previous
#output: Stream:[run]:[number].png where number is 0 1 or 2, denoting xy yz zx plot respectively
def plotStream(streamnum=0,plotStar=1,plotArrow=0):
    CWD = os.getcwd()
    if not('STREAMDATA' in os.listdir(CWD)):
        return 'ERROR: NO STREAMDATA DIRECTORY'
    os.chdir('STREAMDATA')
    for fname in os.listdir(os.getcwd()):
        ret =readData(fname)
        toPlot=ret[1]
        perturberpos = ret[0]

        if len(toPlot[0]) == 0:
            return
        
        #### PLOT XY
        if plotStars ==1:
            plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([perturberpos[0]],[perturberpos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter(toPlot[3],toPlot[4],c='b',alpha=0.25,edgecolors  = 'None')
        if plotArrows == 1:
            for i in toPlot[6]:
                plt.arrow(toPlot[3][i],toPlot[4][i],toPlot[0][i]-toPlot[3][i],toPlot[1][i]-toPlot[4][i],ls='solid')

        plt.axis([-200,200,-200,200]) #kpc
        filename = fname+':0.png'

        plt.xlabel("x")
        plt.ylabel('y')
        plt.savefig(filename)
        plt.clf()
        
        #### PLOT YZ
        plt.scatter(toPlot[1],toPlot[2],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([perturberpos[1]],[perturberpos[2]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter(toPlot[4],toPlot[5],c='b',alpha=0.25,edgecolors  = 'None')
        for i in range(len(toPlot[0])):
            plt.arrow(toPlot[4][i],toPlot[5][i],toPlot[1][i]-toPlot[4][i],toPlot[2][i]-toPlot[5][i],ls='solid')

        
        plt.axis([-200,200,-200,200]) #kpc
        filename = fname+':1.png'

        plt.xlabel("y")
        plt.ylabel('z')
        plt.savefig(filename)
        plt.clf()
        
        #### PLOT ZX
        plt.scatter(toPlot[2],toPlot[0],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([perturberpos[2]],[perturberpos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
        plt.scatter(toPlot[5],toPlot[3],c='b',alpha=0.25,edgecolors  = 'None')
        for i in range(len(toPlot[0])):
            plt.arrow(toPlot[5][i],toPlot[3][i],toPlot[2][i]-toPlot[5][i],toPlot[0][i]-toPlot[3][i],ls='solid')
        
        plt.axis([-200,200,-200,200])
        filename = fname+':2.png'

        plt.xlabel("Z")
        plt.ylabel('x')
        plt.savefig(filename)
        plt.clf()

#plots the closest approahces as a function of change in velocity
# output: closestapproaches.png
def plotClosestApproaches():
    file = open('closestapproaches','r')
    toPlot = [[],[]]
    lines = file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].rsplit('\t')

        for j in range(len(toPlot)):
            toPlot[j].append(float(line[j]))
        i+=1

    print toPlot[0]
    print '\n\n\n'
    print toPlot[1]
    plt.loglog(toPlot[0],toPlot[1],".")

    plt.xlabel("Closest Approach")
    plt.ylabel("Delta V")
    plt.savefig("closestapproachVdv.png")
    plt.clf()

def plotVelocites():
    for file in os.listdir(os.getcwd()):
        data = []
        f = open(file,'r')
        line = f.readlines()[0].rsplit(',')
        line=line[1:-1]
        for i in line:
            data.append(float(i))
        values = set(data)
        count = [[],[]]
        for i in values:
            val = int(i*100)/100.0
            
            count[0].append(val)
            num = 0
            for j in data:
                if abs(val - j)<.05:
                    num +=1
            count[1].append(num)
        
        plt.scatter(count[0],count[1])
        filename = file + '.png'
        plt.savefig(filename)
        plt.clf()

#plotMainData(1,0,[-118.5187002,-134.21555716,15.13030807])

#plotMainData(1,1,[-118.5187002,-134.21555716,15.13030807])

#plotClosestApproaches()
#plotVelocites()
#plotH()
