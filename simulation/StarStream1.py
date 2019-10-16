# comment
import numpy
import matplotlib
import matplotlib.pyplot as plt


DARK = -9999999999


class ColdStream:

# Pick a velocity magnitude (pick from zero-mean unit variance gaussian and multiply by 100)
# Pick a 3 space fiducial point (pick individually from 0-m 1-v gaussian  and multiply by 100)
# Pick a 3 space direction (pick individually from 0-mean 1-variance gaussian and normalize) 
# Generate number of stars specified
    def __init__(self,size,id):
        self.vel = numpy.random.randn()*50

        self.fid = numpy.array([numpy.random.randn(), numpy.random.randn(),numpy.random.randn()]) *100

        self.dirn = numpy.random.randn(3)
        self.dirn = self.dirn/numpy.linalg.norm(self.dirn)
        self.id = id
        self.stars = []


        for i in range(size):
            self.stars.append(Star(self.dirn,self.vel, self.fid, 1, i,id))

class Star:
# Give it mass, vel,id, owner id and dirn from stream
# give it pos from fid*a multiplier taken from 0 mean 200 s.d. gaussian
# create an imageStar to go along with it
    def __init__(self, dirn, vel,fid, mass,id,ownerid):
        self.ownerid = ownerid
        self.id = id
        self.mass = mass
        self.vel = vel*dirn
        self.streamvel = self.vel
#        self.mom = vel*mass
        self.dirn = dirn
        if(id==DARK):
            self.poswrtfid=fid
        else:
            self.poswrtfid = numpy.random.normal(0,2000)
        self.pos = fid + dirn * self.poswrtfid

        self.image =imageStar(dirn,vel,fid,mass,id,ownerid,self.poswrtfid)

    def __str__(self):
        s =  "velocity: "+ str(self.vel) + "\nDirection: " + str(self.dirn) + "\nPosition: " + str(self.pos) +'\n'
        return s
    def update(self, dt, a):
        self.vel = self.vel + a*dt
        self.pos = self.pos + self.vel*dt
        self.dirn = self.vel/numpy.linalg.norm(self.vel)
        self.image.update(dt)
        
        if numpy.linalg.norm(self.vel-self.streamvel) > 10 and self.image.isShowing != 1:
            print "aye"
            print str(self.ownerid) + "." + str(self.id)
            self.image.isShowing = 1

class imageStar(Star):
##image of a star which will show up if there is a significant velocity kick given to it
    def __init__(self,dirn,vel,fid,mass,id,ownerid,poswrtfid):
        self.isShowing = 0
        self.ownerid = ownerid
        self.id = id
        self.mass = mass
        self.vel = vel*dirn
        self.dirn = dirn
        self.poswrtfid = poswrtfid
        self.pos = fid + dirn*poswrtfid
    def update(self,dt):
        self.pos = self.pos+self.vel*dt
        

def plot(seedval,numstreams, numstars,numruns,timestep):
    GM = 10000
    numpy.random.seed(seedval)

    perturber = Star(numpy.array([1,0,0]),100,numpy.array([-250,0,0]),1,DARK,DARK)
#    darkstream = ColdStream(1,1)
    streams = []
    closestd = 10000000000
    closest = -1
    for i in range(numstreams):
        streams.append(ColdStream(numstars,i))
        distance = numpy.abs(numpy.dot(numpy.cross(perturber.dirn,streams[i].dirn), (perturber.pos-streams[i].stars[0].pos))) 
        if distance < closestd:
            closest = streams[i].id
    if closest == -1:
        print "ERROR: There is no closest\n" 
        return 


    for i in range(numruns):    
        toPlotX = []
        toPlotY = []
        toPlotZ = []
        itoPlotX = []
        itoPlotY = []
        itoPlotZ = []
        ctoPlotX = []
        ctoPlotY = []
        ctoPlotZ = []
    
        for j in streams:
            
            for k in j.stars:
                if j.id == closest:
                    ctoPlotX.append(k.pos[0])
                    ctoPlotY.append(k.pos[1])
                    ctoPlotZ.append(k.pos[2])
                    if k.image.isShowing == 1:
                        ctoPlotX.append(k.image.pos[0])
                        ctoPlotY.append(k.image.pos[1])
                        ctoPlotZ.append(k.image.pos[2])
                toPlotX.append(k.pos[0])
                toPlotY.append(k.pos[1])
                toPlotZ.append(k.pos[2])
                if k.image.isShowing == 1:
                    itoPlotX.append(k.image.pos[0])
                    itoPlotY.append(k.image.pos[1])
                    itoPlotZ.append(k.image.pos[2])
                r = perturber.pos - k.pos
#                r = darkstream.stars[0].pos - k.pos
                a = GM*r/(numpy.linalg.norm(r)**3)
         #       print numpy.linalg.norm(a)
          #      print '\n'
          
                k.update(timestep,a)

        perturber.update(timestep,0)


#### PLOT X Y: ALL
#        darkstream.stars[0].update(timestep,0)
        plt.scatter(toPlotX,toPlotY,c='b',alpha=0.5,edgecolors='None')
        if len(itoPlotX) > 0:
            plt.scatter(itoPlotX,itoPlotY,c='g',alpha=0.5, edgecolors='None')
        plt.scatter([perturber.pos[0]],[perturber.pos[1]],s=[100],c='r',alpha=0.9,edgecolors='None', pickradius=1000.0)
#        plt.scatter([darkstream.stars[0].pos[0]],[darkstream.stars[0].pos[1]],s=[100], c='r',alpha=0.9,edgecolors='none', pickradius=1000.0)
        plt.axis([-500,500,-500,500])
        if i < 10:
            filename = 'a0'+str(i)+'.png'
        else:
            filename = 'a'+str(i) + '.png'
        plt.xlabel("x")
        plt.ylabel('y')
        plt.savefig(filename)
        plt.clf()

# #### PLOT Z X: ALL
#         plt.scatter([perturber.pos[2]],[perturber.pos[0]],s=[100],c='r',alpha=0.9,"""edgecolors='none',""" pickradius=1000.0)           
# #        plt.scatter([darkstream.stars[0].pos[2]],[darkstream.stars[0].pos[0]],s=[100], c='r',alpha=0.9,edgecolors='none', pickradius=1000.0)
#         plt.scatter(itoPlotZ,itoPlotX,c='g',alpha=0.5""",edgecolors='none'""")
#         plt.scatter(toPlotZ,toPlotX, c='b',alpha=0.5""",edgecolors='none'""")
                    
#         plt.axis([-500,500,-500,500])
#         if i < 10:
#             filename = 'b0'+str(i)+'.png'
#         else:
#             filename = 'b'+str(i) + '.png'
#         plt.xlabel("z")
#         plt.ylabel('x')
#         plt.savefig(filename)
#         plt.clf()

# #### PLOT Y Z: ALL
#         plt.scatter([perturber.pos[1]],[perturber.pos[2]],s=[100],c='r',alpha=0.9,edgecolors='none', pickradius=1000.0)
# #        plt.scatter([darkstream.stars[0].pos[1]],[darkstream.stars[0].pos[2]],s=[100], c='r',marker='o',alpha=0.9,edgecolors='none')
#         plt.scatter(toPlotY,toPlotZ, c='b',alpha=0.5""", edgecolors='none'""")
#         plt.scatter(itoPlotY,itoPlotZ,c='g',alpha=0.5""",edgecolors='none'""")
#         plt.axis([-500,500,-500,500])
#         if i < 10:
#             filename = 'c0'+str(i)+'.png'
#         else:
#             filename = 'c'+str(i) + '.png'
#         plt.xlabel("y")
#         plt.ylabel('z')
#         plt.savefig(filename)
#         plt.clf()

# #### PLOT XY CLOSEST
#         plt.scatter(ctoPlotX,ctoPlotY,c='b',alpha=0.5""",edgecolors='none'""")
#         plt.scatter([perturber.pos[0]],[perturber.pos[1]],s=[100],c='r',alpha=0.9,"""edgecolors='none',""" pickradius=1000.0)
# #        plt.scatter([darkstream.stars[0].pos[0]],[darkstream.stars[0].pos[1]],s=[100], c='r',alpha=0.9,edgecolors='none', pickradius=1000.0)
#         plt.axis([-500,500,-500,500])
#         if i < 10:
#             filename = 'd0'+str(i)+'.png'
#         else:
#             filename = 'd'+str(i) + '.png'
#         plt.xlabel("x")
#         plt.ylabel('y')
#         plt.savefig(filename)
#         plt.clf()

# #### PLOT ZX CLOSEST
#         plt.scatter([perturber.pos[2]],[perturber.pos[0]],s=[100],c='r',alpha=0.9,"""edgecolors='none',""" pickradius=1000.0)           
# #        plt.scatter([darkstream.stars[0].pos[2]],[darkstream.stars[0].pos[0]],s=[100], c='r',alpha=0.9,edgecolors='none', pickradius=1000.0)
#         plt.scatter(ctoPlotZ,ctoPlotX, c='b',alpha=0.5""",edgecolors='none'""")
#         plt.axis([-500,500,-500,500])
#         if i < 10:
#             filename = 'e0'+str(i)+'.png'
#         else:
#             filename = 'e'+str(i) + '.png'
#         plt.xlabel("z")
#         plt.ylabel('x')
#         plt.savefig(filename)
#         plt.clf()

# #### PLOT YZ CLOSEST
#         plt.scatter([perturber.pos[1]],[perturber.pos[2]],s=[100],c='r',alpha=0.9,"""edgecolors='none',""" pickradius=1000.0)
# #        plt.scatter([darkstream.stars[0].pos[1]],[darkstream.stars[0].pos[2]],s=[100], c='r',marker='o',alpha=0.9,edgecolors='none')
#         plt.scatter(ctoPlotY,ctoPlotZ, c='b',alpha=0.5""", edgecolors='none'""")
#         plt.axis([-500,500,-500,500])
#         if i < 10:
#             filename = 'f0'+str(i)+'.png'
#         else:
#             filename = 'f'+str(i) + '.png'
#         plt.xlabel("y")
#         plt.ylabel('z')
#         plt.savefig(filename)
#         plt.clf()

plot(20,20,500,100,.1)

####mencoder "mf://*.png" -mf fps=10 -o 3str18.mp4 -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=800 
