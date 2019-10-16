""" OLD!
# To Do: - APW: Change method of stepping in parameter space to tune the acceptance ratio down to ~40%
#		 - AJ:  Update notion of 'tolerance' to be applicable to MCMC method

# # # Our goal in this experiment is as follows:
# # #   - create a stream in one potential, stream 1
# # #   - create a random stream in other potential, stream 2
# # #   - slowly perturb stream 2 into stream 1


# # # should be  as follows:
# # #     - create first stream
# # #     - create second stream
# # #     - check stream if perturbation is within desired chi squared
# # #     - perturb
# # #     - see if better:
# # #       - if so keep and repeat
# # #       - if not perturb again
# # #     - when perturbation is within desired similarity quit

# # # Thus the program should be structured as follows:
# # #     - get pert val: center a gaussian around point of std epsilon and repick
# # #     - perturbation function:
# # #        - perturb initial conditions from gaussians that are scaled to be small
# # #     - Main:
# # #        - create stream1 in log potential
# # #        - create stream2 in kep potential
# # #        - evaluate similarity function
# # #        - while similarity > tolerance:
# # #           - perturbation 
# # #           - reevaluate similarity
# # #           - if better stream2 = perturbation
"""

import Similarity
import pickle
from Similarity import Simulator
from Simulator import logging
import scipy.stats as stats
import numpy
import matplotlib
import matplotlib.pyplot as plt
import random
#import logging
random.seed(20)
    
## will perturb a given stream of stars by taking each pos and each vel and reevaluating by sampling from a gaussian around the point of width epsilon
def perturbation(stream, step_info):
    pert = []	
    if step_info[2] == 1:
        for i in range(3):
            pert.append(random.gauss(stream[0][i],step_info[0][i]))
        for j in range(3,6):
            pert.append(stream[0][j])
    else:
        for i in range(3):
            pert.append(stream[0][i])
        for j in range(3,6):
            pert.append(random.gauss(stream[0][j],step_info[0][j]))
	
    pert.extend(stream[0][6:])
#    print "pert[7]: " + str(pert[7])
    logging.debug("Density: " + str(step_info[3]))

    return Similarity.runStream(*Similarity.changeDensity(pert[7],pert[8],pert[9],step_info[3],0,0,pert))#,pert[10]))

def log_posterior_probability(streams, prior_info):
	# Likelihood : Similarity.similarity(streams[0], streams[1])
	# Similarity actually returns a Chi Squared, so we take 1/chisquared as our likelihood
    sim, streams[1] = Similarity.similarity(streams[0],streams[1])
    logging.debug('similarity: %.9f'% sim)
    return sim, streams[1] #Similarity.similarity(streams[0], streams[1]) #prior(prior_info) -> this is not correct now that we are working with logs
        
def prior(prior_info):
	return 1.0
def anneal(run):
    new_temp= initial_temperature*annealing_parameter**(run/int(annealing_dampening))
    return new_temp
def mcmc(streams, prior_info, step_info, nsteps,multiple=0):
	# initially, streams = [stream1, stream2, None]
	chain = []
	probs = []
        dens_vs_sims=[[],[]]
        nsteps_x, nsteps_v = 0,0
	n_accept = 0
        n_accept_v = 0
        n_accept_x = 0
        temperature = initial_temperature
	if plot_streams_on: plot_streams(streams[0],streams[1],-1,-1,multiple)
	for i in range(nsteps): 
		logging.info('MCMC step %i' % (i+1))
		# Generate stream3 from stream2
		step_info[2] = numpy.round(random.uniform(0,1))
                if step_info[2] == 1:
                    nsteps_x +=1
                else:
                    nsteps_v +=1
		streams[2] = perturbation(streams[1], step_info)


		new_streams = [streams[0], streams[2]]
                x = len(streams[1][1].stars)
                y = len(new_streams[1][1].stars)
		log_old_prob, streams[1] = log_posterior_probability(streams, prior_info)
		log_new_prob, streams[2] = log_posterior_probability(new_streams, prior_info)

                dens_vs_sims[0].append(i)
                dens_vs_sims[1].append(log_new_prob)
		
		logging.debug('Log Old Prob: %.9f' % log_old_prob)
		logging.debug('Log New Prob: %.9f' % log_new_prob)
		logging.debug('Log Old Prob - Log New Prob: %.6f' % (log_old_prob - log_new_prob))
		# Compute the proposed step in parameter space of the form [R1,R2,R3,V1,V2,V3]
		#	where the 1,2,3 correspond to the components of the position (R) and 
		#	velocity (V) vectors
		delta_p = list(numpy.array(streams[2][0][:6]) - numpy.array(streams[1][0][:6]))
		if numpy.array(delta_p[3:]).all() == 0:
			logging.debug('Proposed change in initial conditions:\n' + \
			'\tRx = %.6f \n \tRy = %.6f \n \tRz = %.6f' % (delta_p[0],delta_p[1],delta_p[2]))
		else:
			logging.debug('Proposed change in initial conditions:\n' + \
			'\tVx = %.6f \n \tVy = %.6f \n \tVz = %.6f' % (delta_p[3],delta_p[4],delta_p[5]))
		
		# Metropolis-Hastings algorithm
                change = random.uniform(0,1.0)  # numpy.log(0) = -inf
		if (log_old_prob-log_new_prob)/temperature > numpy.log(change): 
			chain.append( (new_streams[1][0][:6], log_old_prob-log_new_prob, delta_p, 1) )
			n_accept += 1
                        if step_info[2] == 1:
                            n_accept_x += 1
                        else:
                            n_accept_v += 1
			probs.append(log_new_prob)
			streams[1] = streams[2]
			logging.debug('Step accepted\n')
                else:
			chain.append( (streams[1][0][:6], log_old_prob-log_new_prob, delta_p, 0) )
			probs.append(log_old_prob)
			logging.debug('Step rejected\n')
                
                if plot_streams_on and i % numpy.floor(nsteps/jump_in_plots) == 0 and i != 0: 
                    plot_streams(streams[0],streams[2],i,i,multiple)

                temperature = anneal(i)
                logging.debug('Temperature = %.9f' % (temperature))

        if plot_streams_on:
            plot_streams(streams[0],streams[2],i,i,multiple)
        density_vs_similarity_plot(dens_vs_sims)
	return chain, probs, n_accept,n_accept_x, n_accept_v,nsteps_x,nsteps_v,streams,step_info

def plot_var_vs_t(raw_chain, probs, acc_frac, acc_frac_x,acc_frac_v,nsteps,input,multiple=0):
    chain = []
    for i in range(nsteps):
        chain.append(raw_chain[:][i][0])

    chain = numpy.array(chain)
    f=open('stuff','w')
    f.write(str(chain))
    f.write(str(jump_in_plots))
    f.close()
    for i in range(3):
        plt.clf()
        plt.plot(numpy.arange(0,nsteps,1), chain[:,i], mec='0.2', marker='.', ls='none', ms=0.5)
        scores = []
        for j in range(jump_in_plots):
            x_0=numpy.floor(nsteps/jump_in_plots)*(j)
            x_1=numpy.floor(nsteps/jump_in_plots)*(j+1)
            ninetyfifth = stats.scoreatpercentile(chain[x_0:x_1,i],95)
            yhi = stats.scoreatpercentile(chain[x_0:x_1,i],97.5)
            fifth = stats.scoreatpercentile(chain[x_0:x_1,i],5)
            ylo = stats.scoreatpercentile(chain[x_0:x_1,i],2.5)
            median = numpy.median(chain[x_0:x_1,i])
#            print "fifth = %.4f, 95=%.4f" % (fifth,ninetyfifth)
#            print "x1=%i,x2=%i" % (x_0,x_1)
            xmin = x_0/nsteps
            xmax = x_1/nsteps
            plt.axhline(fifth,xmin,xmax)
            plt.axhline(ninetyfifth,xmin,xmax)
        plt.axhline(input[i],color='g',label='true value')
#        print [0,nsteps,median-3*(median-ylo),median +3*(yhi-median)]
#        plt.axis([0,nsteps,median-3*(median-ylo),median +3*(yhi-median)])
        plt.xlabel("step")
        plt.ylabel("$R_%i$" % (i+1))
        plt.title("$R_%i$(t)" % (i+1))
        plt.savefig('mcmc-scatter_R%i-%i.png' % (i+1,multiple))

    for i in range(3,6):
        plt.clf()
        plt.plot(numpy.arange(0,nsteps,1), chain[:,i], mec='0.2', marker='.', ls='none', ms=0.5)
        scores = []
        for j in range(jump_in_plots):
            x_0=numpy.floor(nsteps/jump_in_plots)*(j)
            x_1=numpy.floor(nsteps/jump_in_plots)*(j+1)
            ninetyfifth = stats.scoreatpercentile(chain[x_0:x_1,i],95)
            fifth = stats.scoreatpercentile(chain[x_0:x_1,i],5)
            yhi = stats.scoreatpercentile(chain[x_0:x_1,i],97.5)
            ylo = stats.scoreatpercentile(chain[x_0:x_1,i],2.5)
            median = numpy.median(chain[x_0:x_1,i])
            plt.axhline(fifth,x_0/nsteps,x_1/nsteps)
            plt.axhline(ninetyfifth,x_0/nsteps,x_1/nsteps)
        plt.axhline(input[i],color='g',label='true value')
#        plt.axis([0,nsteps,median-3*(median-ylo),median +3*(yhi-median)])
        plt.xlabel("t")
        plt.ylabel("$V_%i$" % (i-2))
        plt.title("$V_%i$(t)" % (i-2))
        plt.savefig('mcmc-scatter_V%i-%i.png' % (i-3,multiple))
        plt.clf()
                    

def plot_scatter(input, raw_chain, probs, acc_frac, acc_frac_x, acc_frac_v, nsteps,multiple=0):
	chain = []
	for i in range(nsteps): chain.append(raw_chain[:][i][0])
	chain = numpy.array(chain)
	
	minProbIndex = probs.index(max(probs))
	plot_var_vs_t(raw_chain,probs,acc_frac,acc_frac_x,acc_frac_v,nsteps,input,multiple)

	for i in range(3):
            plt.clf()
            plt.plot(chain[:,i], chain[:,i+3], mec='0.2', marker='.', ls='none', ms=0.5,label="mcmc step")
            plt.plot([input[i]],[input[i+3]],marker = 'x', label="input")
            #plt.plot([chain[minProbIndex,1]], [chain[minProbIndex,0]], 'r.', ls='none', ms=2.0)
            plt.xlabel("$R_%i$" % (i+1))
            plt.ylabel("$V_%i$" % (i+1))
            plt.legend()
            plt.title("$V_%i$($R_%i$)" % (i+1,i+1))
            plt.savefig('mcmc-scatter_R%i(V%i)-%i.png' % (i+1,i+1,multiple))
        for i in range(3):
            for j in range(3):
                # Plot a scatter plot: R-component vs. V-component
		plt.clf()
                plt.plot(chain[:,j+3], chain[:,i+3], mec='0.2', marker='.', ls='none', ms=0.5,label='mcmc step')
                plt.plot([input[j+3]],[input[i+3]],marker='x',label="input")
		#plt.plot([chain[minProbIndex,1]], [chain[minProbIndex,0]], 'r.', ls='none', ms=2.0)
		plt.xlabel("$V_%i$" % (j+1))
		plt.ylabel("$V_%i$" % (i+1))
                plt.legend()
#		acc_frac = float(n_accept)/float(nsteps)
		plt.title("$V_%i$($V_%i$)" % (i+1,j+1))
		plt.savefig('mcmc-scatter_V%ivsV%i-%i.png' % (i+1,j+1,multiple))
                plt.clf()
                

def plot_histogram(raw_chain, probs, n_accept, nsteps):
	chain = []
	for i in range(nsteps): chain.append(raw_chain[:][i][0])
	chain = numpy.array(chain)
	
	for i in range(3):
		mean_r = numpy.mean(chain[:,i])
		mean_v = numpy.mean(chain[:,i+3])
		sigma_r = numpy.std(chain[:,i])
		sigma_v = numpy.std(chain[:,i+3])
	
		# Plot the histogram for R component
		plt.clf()
		plt.hist( chain[:,i], 50 )
		plt.xlabel("$R_%i$" % (i+1))
		plt.savefig('mcmc-hist-R-%i.png' % (i+1))
	
		# Plot the histogram for V component
		plt.clf()
		n,bins,patches = plt.hist( chain[:,i+3], 50 )
		
		#xx_v = numpy.linspace(numpy.min(bins),numpy.max(bins),150,endpoint=True)
		#yy_v = (1.0 / (numpy.sqrt(2.*numpy.pi)*sigma_b)) * numpy.exp(-(xx_v-mean_v)**2 / (2.*sigma_v**2))
		
		plt.xlabel("$V_%i$" % (i+1))
		#plt.xlim(numpy.min(xx_v),numpy.max(xx_v))
		#plt.plot(xx_v,yy_v*numpy.sum(n)*(bins[1]-bins[0]),'r')
		plt.savefig('mcmc-hist-V-%i.png' % (i+1))
	return None


def density_vs_similarity_plot(dens_sims):
    plt.clf()
    plt.plot(dens_sims[0],dens_sims[1],'o',label='similarity')
    plt.xlabel("density of stars in stream")
    plt.ylabel("Similarity")
    plt.legend()
#    if potential == Similarity.LN:
#        model_pot = "LN"
#    else:
#        model_pot = "KEPLER"
    model_pot = "LN"
    plt.title("Similarity(Density of Stars) for " + model_pot +" fit of LN potential, only varying Density")
    plt.savefig("dens_vs_sim.png")
    plt.clf()

def plot_streams(stream1,stream2,fname,title,multiple=0):
    plt.clf()
    stream1Plot = stream1[1].plot()[:3]
    stream2Plot = stream2[1].plot()[:3]
    fname = str(fname) + '-' + str(multiple)
##    XY
#    plt.scatter(stream1Plot[0],stream1Plot[1],c='b', label='data')
    plt.plot(stream1Plot[0],stream1Plot[1], 'o', label='data')
    plt.plot(stream2Plot[0],stream2Plot[1],'.',label='model')
    plt.xlabel('x')
    plt.ylabel('y')
#    plt.axis([-100,100,-100,100])
#    plt.legend(("data","model"))
    plt.legend()
    plt.title("X-Y Model Curve and Data Curve, Link:%i" % title)
    plt.savefig('a'+str(fname) + '.png')
    plt.clf()
##    YZ
    plt.plot(stream1Plot[1],stream1Plot[2],'o',label='data')
    plt.plot(stream2Plot[1],stream2Plot[2],'.',label='model')
    plt.xlabel('y')
    plt.ylabel('z')
#    plt.legend(("data","model"))
    plt.legend()
    plt.title("Y-Z Model Curve and Data Curve, Link:%i" % title)
    plt.savefig('b'+str(fname) + '.png')
    plt.clf()
##    ZX
    plt.plot(stream1Plot[2],stream1Plot[0],'o',label='data')
    plt.plot(stream2Plot[2],stream2Plot[0],'.',label='model')
    plt.xlabel('z')
    plt.ylabel('x')
#    plt.legend(("data","model"))
    plt.legend()
    plt.title("Z-X Model Curve and Data Curve, Link:%i" % title) 
    plt.savefig('c'+str(fname) + '.png')
    plt.clf()


#### This piece of code is becoming very long
def main(tolerance):
####### set flags and store them in internal variables
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-p', '--plot-streams',
					  dest='plot_streams',action='store_true',
					  help='generate plots of the streams', default=False)
	parser.add_option('-s', '--scatter-plot',
					  dest='scatter', action='store_true',
					  help='generate scatter plots', default=False)
	parser.add_option('-g', '--histograms',
					  dest='histogram', action='store_true',
					  help='generate histograms', default=False)
	parser.add_option('-v', '--verbose',
					  action='store_true', dest='verbose',
					  help='be chatty', default=False)
	parser.add_option('-n', '--nsteps',
					  type='int', dest='nsteps',
					  help='how many mcmc steps', default=10**5)
        parser.add_option('-r','--rho-stars',
                          type='float',dest='density_of_stars',
                          help='how many stars ejected per time step for model orbit',default=1.0)
        parser.add_option('-R', '--Rstep',
                                          type='float', dest='Rstep',
                                          help='size of step in position', default=10.0)
        parser.add_option('-V', '--Vstep',
                                          type='float', dest='Vstep',
                                          help='size of step in velocity', default=10.0)
        parser.add_option('-P', '--Potential',
                                          type='string', dest='potential',
                                          help='potential', default='LN')
        parser.add_option('-j', '--jump', 
                          type='int', dest='jump_in_plots',
                          help='plot and do statistics every nsteps/jump_in_plots steps',default=-1)
        parser.add_option('-a','--annealing_parameter',
                          type='float', dest='annealing_parameter',
                          help='annealing parameter specifies exponential decay rate',default = 1)
        parser.add_option('-d','--dampening_annealing',
                          type='float',dest='annealing_dampening',
                          help='slows down cooling',default=1)
        parser.add_option('-t','--temperature',
                          type='float',dest='initial_temperature',
                          help='initial_temperature',default=1)
        parser.add_option('-e', '--experiment',
                          type='string',dest='experiment_type',
                          help='which experiment you wish to run',default='experiment1')
	(options, args) = parser.parse_args()

	global verbose
        global plot_streams_on
        global jump_in_plots  #jump in plot
        global annealing_parameter
        global initial_temperature
        global annealing_dampening
        global scatter_plots_on
        global histogram_on
        global potential
	verbose = options.verbose
	nsteps = options.nsteps
	Rstep = options.Rstep
	Vstep = options.Vstep
        plot_streams_on = options.plot_streams
        jump_in_plots = options.jump_in_plots
        initial_temperature = options.initial_temperature
        annealing_parameter = options.annealing_parameter
        annealing_dampening = options.annealing_dampening
        density_of_stars = options.density_of_stars
        scatter_plots_on = options.scatter
        histogram_on = options.histogram

        if options.potential=='LN':
            potential = Simulator.LN
        elif options.potential=='KEPLER':
            potential = Simulator.KEPLER
        else:
            raise Exception("Not a valid option for POTENTIAL")

### various variables necessary for the code
        
        prior_info = 1.0

        data_num_stars = 5 #total number of stars
        num_runs_sim = [100,200]  # THIS SHOULD BE VARIABLE



### Setting up verbose mode

	if verbose:
		logging.basicConfig(level=logging.DEBUG, format='%(message)s')
		#logging.raiseExceptions = False
	else:
		logging.basicConfig(level=logging.INFO, format='%(message)s')
	

        timestep_for_integrator = 3
### Setting up the data curve and the first theory curve

#	stream1 = Similarity.createStream(num_runs_sim,1,data_num_stars,'random')
        stream1 = Similarity.runStream(*Similarity.createStream(num_runs_sim,timestep_for_integrator,data_num_stars,'random'))


        times = []
        for t in range(-num_runs_sim[0],0,1):
            times.append(t)
            
        bc = []
        bc.extend(stream1[0])
        bc[6]=numpy.array(times)
        bc[9]=num_runs_sim[0]

	stream2 = Similarity.runStream(*Similarity.createStream(num_runs_sim,timestep_for_integrator,bc[9],bc,potential))
#         print "numruns: " +str(stream2[0][7])
#         print "intermodel distance: ", numpy.linalg.norm(stream2[1].stars[0].pos - stream2[1].stars[1].pos)
#         print "intermodel distance 2: ",numpy.linalg.norm(stream2[1].stars[0].pos - stream2[1].stars[-1].pos)
#         print "data distance: ", numpy.linalg.norm(stream1[1].stars[0].pos - stream1[1].stars[-1].pos)
#         print "data rad:"
#         for i in range(5):
#             print numpy.linalg.norm(stream1[1].stars[i].pos)
#         raise IOError


	init_step_info = [(numpy.array([Rstep,Rstep,Rstep,Vstep,Vstep,Vstep])), None, None,density_of_stars] # 1 for x, 0 for v
#	stream2 = perturbation(stream1, init_step_info) #Similarity.createStream(100,1,10,'random')	

        streams = [stream1,stream2,None]
	if options.experiment_type == "experiment1":
            chain, probs, n_accept,n_accept_x, n_accept_v, nsteps_x, nsteps_v, streams,step_info=  mcmc(streams, \
                                                                                                            prior_info, \
                                                                                                            init_step_info, \
                                                                                                            nsteps)
            acc_frac = n_accept*1.0/nsteps
            acc_frac_x = n_accept_x*1.0/nsteps_x
            acc_frac_v = n_accept_v*1.0/nsteps_v
            logging.info('Total Acceptance Percentage: %.4f' % (acc_frac*100.0))
            logging.info('Position Acceptance Percentage: %.4f' % (acc_frac_x*100.0))
            logging.info('Velocity Acceptance Percentage: %.4f' % (acc_frac_v * 100.0))
            if scatter_plots_on: plot_scatter(stream1[0],chain, probs, acc_frac,acc_frac_x,acc_frac_v, nsteps)
            if histogram_on: plot_histogram(chain, probs, n_accept, nsteps)
#            density_vs_similarity(nsteps,density_of_stars,streams[0],streams[1],step_info)
        elif options.experiment_type == "experiment2":
            Experiment2(streams,prior_info, init_step_info, nsteps)

        else:
            print "Don't you think you're clever!"
            return 1




def density_vs_similarity(nsteps,initial_density,stream1,stream2,init_step_info,prior_info):
    step_info = init_step_info
    dens_vs_sims = [[],[]]

    for i in range(len(step_info[0])):
        step_info[0][i] = 0


    step_info[3] = 0.1
    for i in range(nsteps):
        stream3 = perturbation(stream2,step_info)

        dens_vs_sims[0].append(step_info[3])
        dens_vs_sims[1].append(log_posterior_probability([stream1,stream3],prior_info)[0])
        step_info[3]+= .1
        if plot_streams_on and i % numpy.floor(nsteps/jump_in_plots) == 0 and i != 0: 
            plot_streams(stream1,stream3,i,i)

       
    density_vs_similarity_plot(dens_vs_sims)

def fixDensity(stream1,stream2,initial_density,prior_info):
    functionalTest2(stream1,stream2,100,1)
    sim_old = log_posterior_probability([stream1,stream2],prior_info)[0]
    step_info = [numpy.zeros(6),None,None,1*initial_density]

    functionalTest2(stream1,stream2,101,1)
    stream3 = perturbation(stream2,step_info)
    functionalTest2(stream1,stream3,102,1)
    sim_new = log_posterior_probability([stream1,stream3],prior_info)[0]
    functionalTest2(stream1,stream3,103,1)


    while numpy.abs(sim_new - sim_old) > .1:
        step_info[3] *= 2
        stream3 = perturbation(stream2,step_info)
        sim_old = sim_new
        sim_new = log_posterior_probability([stream1,stream3],prior_info)[0]

    return stream3, step_info[3]
    

##Functional test1: create two streams and evaluate the similarity. plot the points generated by the linear interpolation.

def functionalTest1():
    numruns = [100,200]
    timestep = 1
    numstars = 10
    data = Similarity.runStream(*Similarity.createStream(numruns,timestep,numstars))
    times = []
    for t in range(-numruns[0],0,1):
        times.append(t)
        
    bc = []
    bc.extend(data[0])
    bc[6]=numpy.array(times)
    bc[9]=numruns[0]

    model = Similarity.runStream(*Similarity.createStream(numruns,timestep,bc[9],bc,Simulator.KEPLER))

    print Similarity.similarity(data,model,1)[0]

    closest = Similarity.closestPoints(data,model)

    test = len(closest) - len(data[1].stars)
    if test != 0: 
        print "len(closest) - len(data[1].stars): " + test
    ret = []
    toPlot = []
    
    for j in range(len(closest)):
        i = closest[j][0]
        if i == 0:
            i += 1
        if i == len(model[1].stars)-1:
            i -= 1
        middle,left,right,comparison = numpy.zeros(6),numpy.zeros(6),numpy.zeros(6),numpy.zeros(6)
        middle[:3] = model[1].stars[i].pos[:]
        middle[3:] = model[1].stars[i].vel[:]
        left[:3] = model[1].stars[i-1].pos[:]
        left[3:] = model[1].stars[i-1].vel[:]
        right[:3] = model[1].stars[i+1].pos[:]
        right[3:] = model[1].stars[i+1].vel[:]
        comparison[:3] = data[1].stars[j].pos[:]
        comparison[3:] = data[1].stars[j].vel[:]
        v1 = left - middle
        v2 = right - middle
        comparison = comparison - middle

  
        v1closest = Similarity.vectorProjection(comparison,v1)
        v1closestdist = Similarity.chiSquared(v1closest,comparison)
        v2closest = Similarity.vectorProjection(comparison,v2)
        v2closestdist = Similarity.chiSquared(v2closest,comparison)
        dist = v1closestdist
        vclosest = v1closest+middle
        if v2closestdist < v1closestdist:
            dist = v2closestdist
            vclosest = v2closest+middle
        ret.append(dist)
        toPlot.append(vclosest)


    plt.clf()


    stream1Plot = data[1].plot()[:3]
    stream2Plot = model[1].plot()[:3]
    stream3Plot = [[],[],[]]
    for i in range(3):
        for j in range(len(toPlot)):
            stream3Plot[i].append(toPlot[j][i])


    fname = 1
    title = fname
##    XY
#    plt.scatter(stream1Plot[0],stream1Plot[1],c='b', label='data')
    plt.plot(stream1Plot[0],stream1Plot[1], 'o', label='data')
    plt.plot(stream2Plot[0],stream2Plot[1],'.',label='model')
    plt.plot(stream3Plot[0],stream3Plot[1],'x',label='fit')
    plt.xlabel('x')
    plt.ylabel('y')
#    plt.axis([-100,100,-100,100])
#    plt.legend(("data","model"))
    plt.legend()
    plt.title("X-Y Model Curve and Data Curve, Link:%i" % title)
    plt.savefig('a'+str(fname) + '.png')
    plt.clf()
##    YZ
    plt.plot(stream1Plot[1],stream1Plot[2],'o',label='data')
    plt.plot(stream2Plot[1],stream2Plot[2],'.',label='model')
    plt.plot(stream3Plot[1],stream3Plot[2],'x',label='fit')
    plt.xlabel('y')
    plt.ylabel('z')
#    plt.legend(("data","model"))
    plt.legend()
    plt.title("Y-Z Model Curve and Data Curve, Link:%i" % title)
    plt.savefig('b'+str(fname) + '.png')
    plt.clf()
##    ZX
    plt.plot(stream1Plot[2],stream1Plot[0],'o',label='data')
    plt.plot(stream2Plot[2],stream2Plot[0],'.',label='model')
    plt.plot(stream3Plot[2],stream3Plot[0],'x',label='fit')
    plt.xlabel('z')
    plt.ylabel('x')
#    plt.legend(("data","model"))
    plt.legend()
    plt.title("Z-X Model Curve and Data Curve, Link:%i" % title) 
    plt.savefig('c'+str(fname) + '.png')
    plt.clf()
    

## This is another sanity check for our linear interpolation scheme. We evalutate the two closest points on the model stream and draw a line to them
def functionalTest2(data=0,model=0,multiple=0,ForceEvaluation=0):
    if type(data) == int:
        numruns = [100,200]
        timestep = 1
        numstars = 10
        data = Similarity.runStream(*Similarity.createStream(numruns,timestep,numstars))
        times = []
        for t in range(-numruns[0],0,1):
            times.append(t)
            
            bc = []
            bc.extend(data[0])
            bc[6]=numpy.array(times)
            bc[9]=numruns[0]

        model = Similarity.runStream(*Similarity.createStream(numruns,timestep,bc[9],bc,Simulator.KEPLER))

    
    sim, model = Similarity.similarity(data,model,ForceEvaluation)



    closest = Similarity.closestPoints(data,model,ForceEvaluation)
    
    if type(closest) == type((0,1)):
        raise Exception("bullshit")

    test = len(closest) - len(data[1].stars)
    if test != 0: 
        print "len(closest) - len(data[1].stars): " + test
    ret = []



    ClosestList = []
#    print model[1].stars
    for i in range(len(data[1].stars)):
        closestPoints = ([0,Similarity.chiSquared(data[1].stars[i],model[1].stars[0])],[0,Similarity.chiSquared(data[1].stars[i],model[1].stars[0])])
        for j in range(1,len(model[1].stars)):
            dist = Similarity.chiSquared(data[1].stars[i],model[1].stars[j])
            if closestPoints[0][1] > dist:
                closestPoints = ([j,dist],closestPoints[0])
        ClosestList.append([closestPoints[0][0],closestPoints[1][0]])
        

    toPlot = []
    ret = []
    for j in range(len(closest)):
        i = closest[j][0]

        if (i == 0 or i == len(model[1].stars)-1) and not ForceEvaluation:
            raise Exception("even more bullshit")


        middle,left,right,comparison = numpy.zeros(6),numpy.zeros(6),numpy.zeros(6),numpy.zeros(6)
        middle[:3] = model[1].stars[i].pos[:]
        middle[3:] = model[1].stars[i].vel[:]
        left[:3] = model[1].stars[i-1].pos[:]
        left[3:] = model[1].stars[i-1].vel[:]
        right[:3] = model[1].stars[i+1].pos[:]
        right[3:] = model[1].stars[i+1].vel[:]
        comparison[:3] = data[1].stars[j].pos[:]
        comparison[3:] = data[1].stars[j].vel[:]

        v1 = left - middle
        v2 = right - middle
        comparison = comparison - middle

  
        v1closest = Similarity.vectorProjection(comparison,v1)
        v1closestdist = Similarity.chiSquared(v1closest,comparison)
        v2closest = Similarity.vectorProjection(comparison,v2)
        v2closestdist = Similarity.chiSquared(v2closest,comparison)
        dist = v1closestdist
        closestVect = v1closest

        if v2closestdist < v1closestdist:
            dist = v2closestdist
            closestVect = v2closest

        ret.append(dist)
        toPlot.append(closestVect+middle)


    plt.clf()


    stream1Plot = data[1].plot()
    stream2Plot = model[1].plot()
    stream3Plot = [[],[],[],[],[],[]]
    for i in range(6):
        for j in range(len(toPlot)):
            stream3Plot[i].append(toPlot[j][i])


    fname = '1functionalTest2-'+str(multiple)
    title = 1

    plt.clf()
    plt.title('This is a title')
### Pos Vs Pos
    for j in range(3):
        plt.subplot(3,3,j+1)
        plt.plot(stream1Plot[j],stream1Plot[(j+1)%3], 'o', label='data')
        plt.plot(stream2Plot[j],stream2Plot[(j+1)%3],'.',label='model')
        plt.plot(stream3Plot[j],stream3Plot[(j+1)%3],'x',label='fit')
        for i in range(len(ClosestList)):
            plt.arrow(stream1Plot[j][i],stream1Plot[(j+1)%3][i],stream2Plot[j][ClosestList[i][0]]-stream1Plot[j][i],stream2Plot[(j+1)%3][ClosestList[i][0]]-stream1Plot[(j+1)%3][i],ls='solid')
            plt.arrow(stream1Plot[j][i],stream1Plot[(j+1)%3][i],stream2Plot[j][ClosestList[i][1]]-stream1Plot[j][i],stream2Plot[(j+1)%3][ClosestList[i][1]]-stream1Plot[(j+1)%3][i],ls='solid')
        plt.xlabel('R_%i' % j)
        plt.ylabel('R_%i' % ((j+1)%3))
    
### Pos Vs Vel
    for j in range(3):
        plt.subplot(3,3,j+4)
        plt.plot(stream1Plot[j],stream1Plot[j+3],'o',label='data')
        plt.plot(stream2Plot[j],stream2Plot[j+3],'.',label='model')
        plt.plot(stream3Plot[j],stream3Plot[j+3],'x',label='fit')
        for i in range(len(ClosestList)):
            plt.arrow(stream1Plot[j][i],stream1Plot[j+3][i],stream2Plot[j][ClosestList[i][0]]-stream1Plot[j][i],stream2Plot[j+3][ClosestList[i][0]]-stream1Plot[j+3][i],ls='solid')
            plt.arrow(stream1Plot[j][i],stream1Plot[j+3][i],stream2Plot[j][ClosestList[i][1]]-stream1Plot[j][i],stream2Plot[j+3][ClosestList[i][1]]-stream1Plot[j+3][i],ls='solid')


        plt.xlabel('R_%i' % j)
        plt.ylabel('V_%i' % j)

### Vel vs Vel
    for j in range(3):
        plt.subplot(3,3,j+7)
        plt.plot(stream1Plot[j+3],stream1Plot[((j+1)%3)+3], 'o', label='data')
        plt.plot(stream2Plot[j+3],stream2Plot[((j+1)%3)+3],'.',label='model')
        plt.plot(stream3Plot[j+3],stream3Plot[((j+1)%3)+3],'x',label='fit')
        for i in range(len(ClosestList)):
            plt.arrow(stream1Plot[j+3][i],stream1Plot[((j+1)%3)+3][i],stream2Plot[j+3][ClosestList[i][0]]-stream1Plot[j+3][i],stream2Plot[((j+1)%3)+3][ClosestList[i][0]]-stream1Plot[((j+1)%3)+3][i],ls='solid')
            plt.arrow(stream1Plot[j+3][i],stream1Plot[((j+1)%3)+3][i],stream2Plot[j+3][ClosestList[i][1]]-stream1Plot[j+3][i],stream2Plot[((j+1)%3)+3][ClosestList[i][1]]-stream1Plot[((j+1)%3)+3][i],ls='solid')
        plt.xlabel('V_%i' % j)
        plt.ylabel('V_%i' % ((j+1)%3))
#    plt.show()
    plt.savefig('c'+str(fname) + '.png')
    plt.close()
    plt.clf()





###this is the second experiment we will do. In this we will create a data set, add gaussian noise drawn from a Gaussian with noise specified by the matrix


'''

The goal of this experiment is to show that our fitting system will return similar answers given
a data set and the same data set with some gaussian noise provided that the gaussian noise comes
from a distribution with the same covariance matrix as the inverse of the metric used.

The experiment will be run as follows:
  - create data set.
  - for each point in data set create 6-d noise vector.
  - add nosie to corresponding point.
  - fit to original data set
  - fit to new data set
  - compare fits. (?)


'''

def functionalTest3(streams,prior_info,init_step_info,nsteps):
    originalStreams = pickle.loads(pickle.dumps(streams)) # this is my fix to pass by value

    streams[1],init_step_info[3] = fixDensity(streams[0],streams[1],init_step_info[3],prior_info)

    sim1, trash = Similarity.similarity(streams[0],streams[1])
    for star in streams[0][1].stars:
        epsilon = Similarity.createNoise(Similarity.Metric(Similarity.std(0,0)))
        star.pos[:] += epsilon[:3]
        star.vel[:] += epsilon[3:]

    sim2, trash = Similarity.similarity(streams[0],streams[1])

    print "Similarity1 - Similarity2: ", sim1 - sim2
    print "Similairty1: ", sim1
    print "Similarity2: ", sim2
    

def Experiment2(streams,prior_info,init_step_info,nsteps):
    originalStreams = pickle.loads(pickle.dumps(streams)) # this is my fix to pass by value


#    print "length of stream: ", numpy.linalg.norm(streams[0][1].stars[0].pos - streams[0][1].stars[-1].pos)
#    print "dist from center: ",numpy.linalg.norm(streams[0][1].stars[0].pos)
#    print "velocity: ",numpy.linalg.norm(streams[0][1].stars[0].vel)



    sim1, streams[1] = log_posterior_probability([streams[0],streams[1]],prior_info)

#    print "Trash: ",trash
#    print "streams[1]: ",streams[1]

    
    for star in streams[0][1].stars:
        epsilon = Similarity.createNoise(Similarity.Metric(Similarity.std(0,0)))
        star.pos[:] += epsilon[:3]
        star.vel[:] += epsilon[3:]



    sim2, streams[1] = log_posterior_probability([streams[0],streams[1]],prior_info)



    print "Similarity1 - Similarity2: ", sim1 - sim2
    print "Similairty1: ", sim1
    print "Similarity2: ", sim2
    original_density = init_step_info[3]
    streams[1],init_step_info[3] = fixDensity(streams[0],streams[1],init_step_info[3],prior_info)

#    density_vs_similarity(nsteps,init_step_info[3],streams[0],streams[1],init_step_info,prior_info)
    

    chain, probs, n_accept,n_accept_x, n_accept_v, nsteps_x, nsteps_v, streams,step_info=  mcmc(streams, \
                                                                                                    prior_info, \
                                                                                                    init_step_info, \
                                                                                                    nsteps)

                                                                                                            
    acc_frac = n_accept*1.0/nsteps
    acc_frac_x = n_accept_x*1.0/nsteps_x
    acc_frac_v = n_accept_v*1.0/nsteps_v
    logging.info('Total Acceptance Percentage: %.4f' % (acc_frac*100.0))
    logging.info('Position Acceptance Percentage: %.4f' % (acc_frac_x*100.0))
    logging.info('Velocity Acceptance Percentage: %.4f' % (acc_frac_v * 100.0))

#    print originalStreams[0][0]
#    print streams[0][1]
    if scatter_plots_on:
        plot_scatter(originalStreams[0][0],chain, probs, acc_frac,acc_frac_x,acc_frac_v, nsteps)

    
    functionalTest2(streams[0],streams[1],1)


    originalStreams[1], init_step_info[3] = fixDensity(originalStreams[0], originalStreams[1], original_density,prior_info)
    chain, probs, n_accept,n_accept_x, n_accept_v, nsteps_x, nsteps_v, originalStreams,step_info=  mcmc(originalStreams, \
                                                                                                    prior_info, \
                                                                                                    init_step_info, \
                                                                                                    nsteps,1)

    acc_frac = n_accept*1.0/nsteps
    acc_frac_x = n_accept_x*1.0/nsteps_x
    acc_frac_v = n_accept_v*1.0/nsteps_v
    logging.info('Total Acceptance Percentage: %.4f' % (acc_frac*100.0))
    logging.info('Position Acceptance Percentage: %.4f' % (acc_frac_x*100.0))
    logging.info('Velocity Acceptance Percentage: %.4f' % (acc_frac_v * 100.0))

    functionalTest2(originalStreams[0],originalStreams[1],2)
    if scatter_plots_on:
        plot_scatter(originalStreams[0][0],chain, probs, acc_frac,acc_frac_x,acc_frac_v, nsteps,1)
    
    sim , trash = Similarity.similarity(originalStreams[1],streams[1],1)
    logging.info("Similarity of fits: %.4f" % sim)



if __name__ == '__main__':
#    functionalTest2()
#    functionalTest1()
    main(1)

    
