import os
import matplotlib
from matplotlib import pyplot as plt
import math

## reads data files put out by Observation.py
def readData(filename):
    file = open(filename,'r')
    toPlot = [[],[]]
    perturberpos = []
    lines = file.readlines()
    line = lines[1].rsplit('\t')
    for i in range(len(lines))[1:]:
        line = lines[i].rsplit('\t')
        for j in range(len(toPlot)):
            toPlot[j].append(float(line[j+1]))
    perturberpos.append(toPlot[0].pop())
    perturberpos.append(toPlot[1].pop())
    return (perturberpos,toPlot)

## plots all of the data using a hammer aitoff projection
def plotMainDataAitoff():
    STARTDIR=os.getcwd()
    os.mkdir('OUT')
    PLOTSDIR=STARTDIR+'/OUT'
    if not('data' in os.listdir(STARTDIR)):
        print 'Error: No Data Folder In this Directory'
        return 1

    os.chdir('data')
    MAINDATADIR = os.getcwd()

    for fname in os.listdir(MAINDATADIR):
        print fname
        ret = readData(fname)
        toPlot = ret[1]
        perturberpos=ret[0]
        os.chdir(PLOTSDIR)

        if len(toPlot[0]) == 0:
            return

        #### Plot XY: all
        fig = plt.figure()
        ax=fig.add_subplot(111,projection='hammer')
        print min(toPlot[1]), max(toPlot[1])
        print min(toPlot[0]), max(toPlot[0])
        ax.scatter(toPlot[1],toPlot[0],c='b',alpha=0.5,edgecolors='None')
        ax.scatter([perturberpos[1]],[perturberpos[0]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)           
        #plt.axis([-math.pi/2,math.pi/2,-math.pi,math.pi])
        
        filename = 'a'+fname + '.png'

#        ax.xlabel("l")
#        ax.ylabel('b')
        fig.savefig(filename)
        fig.clf()

        
        os.chdir(MAINDATADIR)
    print 'happend'


    os.chdir(STARTDIR)

## plots the data using a cartesian projection
def plotMainData():
    STARTDIR=os.getcwd()
    os.mkdir('OUT')
    PLOTSDIR=STARTDIR+'/OUT'
    if not('data' in os.listdir(STARTDIR)):
        print 'Error: No Data Folder In this Directory'
        return 1

    os.chdir('data')
    MAINDATADIR = os.getcwd()

    for fname in os.listdir(MAINDATADIR):
        print fname
        ret = readData(fname)
        toPlot = ret[1]
        perturberpos=ret[0]
        os.chdir(PLOTSDIR)

        if len(toPlot[0]) == 0:
            return

        #### Plot XY: all

        
        plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
        plt.scatter([perturberpos[0]],[perturberpos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)           
        plt.axis([-math.pi/2,math.pi/2,-math.pi,math.pi])
        
        filename = 'm'+fname + '.png'

        plt.xlabel("l")
        plt.ylabel('b')
        plt.savefig(filename)
        plt.clf()

        
        os.chdir(MAINDATADIR)
    print 'happend'


    os.chdir(STARTDIR)


# plots the data around a patch of star using a euclidean projection by assuming that the patch is sufficiently small
def plotPatch(run,star):
    filename = 'patchof:'+str(star)+':'+str(run)
    ret = readData(filename)
    toPlot = ret[1]
    toPlot[0].append(ret[0][0])
    toPlot[1].append(ret[0][1])
    plt.scatter(toPlot[0],toPlot[1],c='b',alpha=0.5,edgecolors='None')
    plt.axis([-1,1,-1,1])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(filename+'.png')
    
plotMainDataAitoff()
#plotPatch('1950',127)
