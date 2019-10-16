import MCMC
import pickle
from MCMC import Similarity
import copy
from Similarity import Simulator
from Simulator import logging
import numpy
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt
import random

SEEDVAL = 30
    

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
					  help='how many mcmc steps', default=50)
        parser.add_option('-r','--rate-of-sample',
                          type='float',dest='sample_rate',
                          help='rate of sampling of model',default=1.0)
        parser.add_option('-R', '--Rstep',
                                          type='float', dest='Rstep',
                                          help='size of step in position', default=0.43)
        parser.add_option('-V', '--Vstep',
                                          type='float', dest='Vstep',
                                          help='size of step in velocity', default=0.0035)
        parser.add_option('-P', '--Potential',
                                          type='string', dest='potential',
                                          help='potential', default='FLATLN')
        parser.add_option('-j', '--jump', 
                          type='int', dest='jump_in_plots',
                          help='plot and do statistics every nsteps/jump_in_plots steps',default=-1)
        parser.add_option('-a','--annealing_parameter',
                          type='float', dest='annealing_parameter',
                          help='annealing parameter specifies exponential decay rate',default = 1)
        parser.add_option('-d','--dampening_annealing',
                          type='float',dest='annealing_dampening',
                          help='slows down cooling',default=1)
        parser.add_option('-K','--temperature',
                          type='float',dest='initial_temperature',
                          help='initial_temperature',default=1)
        parser.add_option('-e', '--experiment',
                          type='string',dest='experiment_type',
                          help='which experiment you wish to run',default='experiment1')
        parser.add_option('-b', '--begin',
                          type='int',dest='begin',
                          help='begin hist',default=0)
        parser.add_option('-E', '--end',
                          type='int',dest='end',
                          help='which experiment you wish to run',default=-1)

	parser.add_option('-f','--filename-pickles',
			  type='string',dest='filename_pickles', nargs = 2,
			  help='filename of pickles to be used for plots', default='')
	parser.add_option('-t', '--thickness--tensor-eigenvalues',
			  type='float', dest='thickness_tensor_eigenvalues', nargs = 2,
			  help='Eigenvalues of Thickness Tensor', default = (0.000001,0.0))
	parser.add_option('-T', '--eigenvalue-step-size',
			  type='float',dest='thickness_step',nargs = 2,
			  help='stepsize in T_11 and T_44', default=(0,0) )
	parser.add_option('-Q','--Qstep',
			  type = 'float', dest = 'qstep',
			  help = 'flattening parameter step size', default = 0)
	parser.add_option('-q','--qflattenning-parameter',
			  type = 'float', dest = 'init_flatparam',
			  help = 'initial choice of flattening parameter', default = 1.0)
	parser.add_option('-C','--vcircsquared-step',
			  type = 'float', dest = 'vcircsquared_step',
			  help = 'step size of Vcircsquared', default = 0)
	parser.add_option('-c','--vcircsquared',
			  type = 'float', dest = 'init_vcircsquared',
			  help = 'initial choice of Vcircsquared', default=Simulator.V2POTENTIAL)
	parser.add_option('-M','--step--size',
			  type = 'int', dest='num_step',
			  help = 'stepsize in M', default=1)
	parser.add_option('-S','--seed',
			  type = 'int', dest='seed',
			  help = 'random seed for Rotate', default=1)

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
	global prior_info
	global data_num_stars
	global num_runs_sim
	global timestep_for_integrator
	global init_step_info
	global nsteps
	global init_flatparam
	global init_vcircsquared
	global thickness_tensor_eigenvalues
	verbose = options.verbose
	nsteps = options.nsteps
	Rstep = options.Rstep
	Vstep = options.Vstep
	qstep = options.qstep
	thickness_step = options.thickness_step
	thickness_tensor_eigenvalues = options.thickness_tensor_eigenvalues
        plot_streams_on = options.plot_streams
        jump_in_plots = options.jump_in_plots
        initial_temperature = options.initial_temperature
        annealing_parameter = options.annealing_parameter
        annealing_dampening = options.annealing_dampening
        sample_rate = options.sample_rate
        scatter_plots_on = options.scatter
        histogram_on = options.histogram
	init_flatparam = options.init_flatparam
	init_vcircsquared = options.init_vcircsquared
	V2Cstep = options.vcircsquared_step

        if options.potential=='LN':
		potential = Simulator.LN
        elif options.potential=='KEPLER':
		potential = Simulator.KEPLER
	elif options.potential=='FLATLN':
		potential = Simulator.FLATLN
        else:
            raise Exception("Not a valid option for POTENTIAL")

### various variables necessary for the code
        

        prior_info = 1.0

        data_num_stars = 20  #total number of stars
        num_runs_sim = [200,200]  # THIS SHOULD BE VARIABLE
        timestep_for_integrator = .1

### Setting up verbose mode

	if verbose:
		logging.basicConfig(level=logging.DEBUG, format='%(message)s')
		#logging.raiseExceptions = False
	else:
		logging.basicConfig(level=logging.INFO, format='%(message)s')
	


### Setting up the data curve and the first theory curve

#        def __init__(self,seedVal=20,numStreams=1,timestep=1,potential=LN, fiducialPoints=None,numForward=10,numBackward=10,integrationScheme='forwardEuler',numRemaining=-1,isData=False,hasNoise=False):



	init_step_info = [None,options.num_step,numpy.array([Rstep,Rstep,Rstep,Vstep,Vstep,Vstep]),qstep,V2Cstep,thickness_step,None]

	if options.experiment_type == "mcmc_fit":
            mcmc_fit()
        elif options.experiment_type == "experiment2":
            Experiment2()
        elif options.experiment_type == 'samplerate_vs_similarity':
            samplerate_vs_similarity()
	elif options.experiment_type == 'makeMcmcPlots':
		makeMcmcPlots(*options.filename_pickles)
	elif options.experiment_type == 'makeExamplePlots':
		makeExamplePlots(*options.filename_pickles)
	elif options.experiment_type == 'keep_on_truckin':
		keep_on_truckin(options.filename_pickles[0])
	elif options.experiment_type == 'autocorelation':
		autocorelation(options.filename_pickles[0])
	elif options.experiment_type == 'QPlots':
		QPlots()
	elif options.experiment_type == 'LogRadPosErr':
		LogRadPosErr()
	elif options.experiment_type == 'Rotation':
		Rotation(options.seed)
	elif options.experiment_type == 'makeMultPlots':
		makeMultPlots(args)
        else:
            print "Nah, you good."
            return 1



## This is not so much an experiment as just running an MCMC fit

# MCMC    def __init__(self,streams, prior_info, step_info,initial_temperature,jump_in_plots=1,chain_id=0,annealing_parameter=1,annealing_dampening=1):

def mcmc_fit(hasNoise=False,begin=0,end=-1):
	data = Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1], numRemaining=data_num_stars, isData=True,hasNoise=hasNoise).streams[0]	
	model = Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[data.fiducialPoint.pos],numForward=num_runs_sim[0],numBackward=num_runs_sim[1],potential=potential).streams[0]    
	mcmc = MCMC.MonteCarloMarkovChain(streams,prior_info,init_step_info,initial_temperature,jump_in_plots,annealing_parameter = annealing_parameter, annealing_dampening=annealing_dampening)

	mcmc.run(nsteps,plot_streams_on=plot_streams_on)

	acc_frac = mcmc.n_accept*1.0/nsteps
	acc_frac_x = mcmc.n_accept_x*1.0/mcmc.nsteps_x
	acc_frac_v = mcmc.n_accept_v*1.0/mcmc.nsteps_v
	logging.info('Total Acceptance Percentage: %.4f' % (acc_frac*100.0))
	logging.info('Position Acceptance Percentage: %.4f' % (acc_frac_x*100.0))
	logging.info('Velocity Acceptance Percentage: %.4f' % (acc_frac_v * 100.0))
	if scatter_plots_on: mcmc.plot_scatter()
	if histogram_on: mcmc.plot_histogram(begin=begin,end=end)





def samplerate_vs_similarity(nsteps=5,initial_samplerate=.4,hasNoise=True):
	data = Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=hasNoise,potential=potential).streams[0]	
	model = Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[data.fiducialPoint.pos],numForward=num_runs_sim[0],numBackward=num_runs_sim[1],potential=potential).streams[0]    


	sr_vs_sim=[[],[]]
	oldsim= 10000.0
	for i in range(nsteps):
		sr = initial_samplerate*(2**i)
		print sr
		stream3 = model.copy()
		stream3.changeSampleRate(sr)
		sr_vs_sim[0].append(numpy.log(sr)/numpy.log(2))
		newsim = -.5*Similarity.similarity(data,stream3)
		sr_vs_sim[1].append(oldsim-newsim)
		oldsim = newsim

	plt.clf()
	plt.plot(sr_vs_sim[0],sr_vs_sim[1],'o',label='Badness-of-fit')
#	plt.axis(ymin=-3,ymax=10)
	plt.xlabel("$log_2(\\text{density})$ of stars in stream (stars per Myr)")
	plt.ylabel("$\\frac{d\\text{Badness-of-fit}}{d log_2(\\text{density})}")
	plt.legend()
	if potential == Simulator.LN:
		model_pot = "LN"
	else:
		model_pot = "KEPLER"

#	plt.title("True: " + data.potentialToString() + " Model: " + model.potentialToString())
	plt.savefig("SampleRate_vs_Similarity")
	plt.clf()


    

# def fixDensity(stream1,stream2,initial_density,prior_info):
#     functionalTest2(stream1,stream2,100,1)
#     sim_old = log_posterior_probability([stream1,stream2],prior_info)[0]
#     step_info = [numpy.zeros(6),None,None,1*initial_density]

#     functionalTest2(stream1,stream2,101,1)
#     stream3 = perturbation(stream2,step_info)
#     functionalTest2(stream1,stream3,102,1)
#     sim_new = log_posterior_probability([stream1,stream3],prior_info)[0]
#     functionalTest2(stream1,stream3,103,1)


#     while numpy.abs(sim_new - sim_old) > .1:
#         step_info[3] *= 2
#         stream3 = perturbation(stream2,step_info)
#         sim_old = sim_new
#         sim_new = log_posterior_probability([stream1,stream3],prior_info)[0]

#     return stream3, step_info[3]
    



# '''

# The goal of this experiment is to show that our fitting system will return similar answers given
# a data set and the same data set with some gaussian noise provided that the gaussian noise comes
# from a distribution with the same covariance matrix as the inverse of the metric used.

# The experiment will be run as follows:
#   - create data set.
#   - for each point in data set create 6-d noise vector.
#   - add nosie to corresponding point.
#   - fit to original data set
#   - fit to new data set
#   - compare fits. (?)


# '''
    
def keep_on_truckin(mcmcpickle):
	f = open(mcmcpickle,'r')
	mcmc = pickle.load(f)

	mcmc.run(nsteps,plot_streams_on=plot_streams_on)
	acc_fracs = numpy.zeros(len(mcmc.n_accept))
	for i in range(len(acc_fracs)):
		acc_fracs[i] = mcmc.n_accept[i]*1.0/mcmc.nsteps[i]

	
	logging.info('Total Acceptance Percentage: %.4f' % (acc_fracs[0]*100.0))
	logging.info('acc_fracs: [ M , R, V, Q, C, Sigma_x^2]')
	logging.info('acc_fracs: ' + str(acc_fracs[1:]*100))
	logging.info('steps: ' + str(mcmc.nsteps))
	logging.info('accsteps: ' + str(mcmc.n_accept))
				       


#	logging.info("True Params:" + str(clean_data.get_params()))
	logging.info("Final Step:" + str(mcmc.chain[-1]))

	f = open("mcmcadded.pickle",'w')
	pickle.dump(mcmc,f)
	f.close()

def autocorelation(mcmcpickle):
	f = open(mcmcpickle,'r')
	mcmc = pickle.load(f)
	mcmc.autocorelation(jump_in_plots)



	

def Experiment2():

#	clean_sim =Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,fiducialPoints = [obspos],numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=init_flatparam, vcircsquared=init_vcircsquared,thicknessEigs = thickness_tensor_eigenvalues)
	clean_sim =Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=init_flatparam, vcircsquared=init_vcircsquared,thicknessEigs = thickness_tensor_eigenvalues)
	clean_data = clean_sim.streams[0]

# 	dists = []
	
# 	for star in clean_data.stars:
# 		dists.append(numpy.linalg.norm(star.pos[:3] - Simulator.EARTH[:3])) 


# 	dists = numpy.array(dists)
# 	print 'closest dist:', dists.min()
# 	print 'closest index:', dists.argmin()
# 	print 'check: ', dists[dists.argmin()]
# 	closest_pos = clean_data.stars[dists.argmin()].pos
# 	L = numpy.cross(closest_pos[:3], closest_pos[3:])
# 	Lhat = L/numpy.linalg.norm(L)
# 	print 'ang_mom of data',numpy.dot(Lhat,[1,0,0]), numpy.dot(Lhat, [0,1,0]), numpy.dot(Lhat,[0,0,1])

	model = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[clean_data.fiducialPoint.pos],numForward=num_runs_sim[0],numBackward=num_runs_sim[1], potential=potential, flatparam=init_flatparam,vcircsquared=init_vcircsquared,thicknessEigs=thickness_tensor_eigenvalues).streams[0]    


	noisy_sim = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[clean_data.fiducialPoint.pos], numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True, hasNoise=True, potential=Simulator.FLATLN, flatparam=init_flatparam, vcircsquared=init_vcircsquared, thicknessEigs = thickness_tensor_eigenvalues)

	noisy_data = noisy_sim.streams[0]

	print 'thickness Eigs:',noisy_data.thicknessEigs
	print 'fiducialPos', noisy_data.fiducialPoint.pos
	print 'covtens', noisy_data.fiducialPoint.covarianceMatrix

        relpos = model.fiducialPoint.pos-Simulator.EARTH
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
	change_basis = numpy.array([e1,e2,e3,e4,e5,e6])
#	print numpy.dot(change_basis, numpy.dot(noisy_data.fiducialPoint.covarianceMatrix,change_basis.T))
#	raise IOerror

	print 'marginalized likelyhood perturbed', -.5*Similarity.similarity(noisy_data,model)
	print 'marginalized likelyhood unperturbed', -.5*Similarity.similarity(clean_data,model)	
	f = open('clean_sim.pickle','w')
	pickle.dump(clean_sim,f)
	f.close()
	f = open('noisy_sim.pickle','w')
	pickle.dump(noisy_sim,f)
	f.close()
#	raise IOError



	mcmc = MCMC.MonteCarloMarkovChain(numpy.random.get_state(),noisy_data,model, prior_info=prior_info,step_info=init_step_info,initial_temperature=initial_temperature,jump_in_plots=jump_in_plots,annealing_parameter=annealing_parameter,annealing_dampening=annealing_dampening)

	mcmc.run(nsteps)
#	mcmc.func_test()
#	raise Exception('Your mom')
#	mcmc.fake_run(nsteps,0)

#	mcmc.run(nsteps,plot_streams_on=plot_streams_on)
# 	acc_fracs = numpy.zeros(len(mcmc.n_accept))
# 	for i in range(len(acc_fracs)):
# 		acc_fracs[i] = mcmc.n_accept[i]*1.0/mcmc.nsteps[i]

	
# 	logging.info('Total Acceptance Percentage: %.4f' % (acc_fracs[0]*100.0))
# 	logging.info('acc_fracs: [ M , R, V, Q, C, Sigma_x^2]')
# 	logging.info('acc_fracs: ' + str(acc_fracs[1:]*100))
# 	logging.info('steps: ' + str(mcmc.nsteps))
#	logging.info('accsteps: ' + str(mcmc.n_accept))
	if scatter_plots_on:
		mcmc.plot_scatter()

	if plot_streams_on:
		toPlot = []
		for i in range(0,len(mcmc.chain),len(mcmc.chain)/jump_in_plots):
			toPlot.append(mcmc.chain[i])
		S = Simulator.Simulation(SEEDVAL,numStreams=len(toPlot),timestep=timestep_for_integrator, fiducialPoints = toPlot, numForward=num_runs_sim[0], numBackward=num_runs_sim[1], numRemaining=data_num_stars)
		S.plotAllStreams(noisy_data.plot(),clean_data.potentialToString())
				       
# 	sim1 = mcmc.log_posterior_probability(mcmc.data,mcmc.model)
# 	sim2 = mcmc.log_posterior_probability(clean_data,mcmc.model)
# 	epsilon1 = -2*sim1
# 	for star in mcmc.data.stars:
# 		epsilon1 -= numpy.log(numpy.linalg.det(star.covarianceMatrix+star.thicknessTensor))
# 	epsilon2 = -2*sim2
# 	for star in clean_data.stars:
# 		epsilon2 -= numpy.log(numpy.linalg.det(star.covarianceMatrix+star.thicknessTensor))



	logging.info("True Params:" + str(clean_data.get_params()))
	logging.info("Final Step:" + str(mcmc.chain[-1]))
#	logging.info("Log Posterior Probability of model for noisy data given noisy data: %.8f" % sim1)
#	logging.info("Log Posterior Probability of model for noisy data given the truth: %.8f "% sim2)
#	logging.info("Epsilon for model of noisy data given noisy data: %.8f" % epsilon1)
#	logging.info("Epsilon for model of noisy data given the truth: %.8f" % epsilon2)
#	logging.info("difference: %.8f" % (sim1-sim2))

	f = open("mcmc.pickle",'w')
	mcmc.brine()
	pickle.dump(mcmc,f)
	f.close()

	q_vals = mcmc.chain[:,7,nsteps/2:].flatten()
	ninetyseventh=(stats.scoreatpercentile(q_vals,97.5))
	twopointfifth=(stats.scoreatpercentile(q_vals,2.5))
	print '97.5', ninetyseventh
	print '2.5', twopointfifth


def QPlots():
	clean_sim =Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=init_flatparam, vcircsquared=init_vcircsquared,thicknessEigs = thickness_tensor_eigenvalues)
	clean_data = clean_sim.streams[0]

	model_sim = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[clean_data.fiducialPoint.pos], numForward=num_runs_sim[0],numBackward=num_runs_sim[1],potential=potential,flatparam=init_flatparam,vcircsquared=init_vcircsquared,thicknessEigs=thickness_tensor_eigenvalues)
	model = model_sim.streams[0]


	noisy_data = []
	print 'data num stars',data_num_stars

	for i in xrange(5):
		num_stars = data_num_stars + i*5

		noisy_sim = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[clean_data.fiducialPoint.pos], numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=num_stars, isData=True, hasNoise=True, potential=Simulator.FLATLN, flatparam=init_flatparam, vcircsquared=init_vcircsquared, thicknessEigs = thickness_tensor_eigenvalues)

		noisy_data.append(noisy_sim.streams[0])
	random_state = numpy.random.get_state()
	chains = [MCMC.MonteCarloMarkovChain(random_state,datum,model, prior_info=prior_info,step_info=init_step_info,jump_in_plots=jump_in_plots) for datum in noisy_data]

	ninetyseventh=[]
	twopointfifth = []
	i = 0
	for chain in chains:
		print 'run', i
		chain.run(nsteps)
		print chain.chain[:,7,nsteps/2:].flatten().shape
		q_vals = chain.chain[:,7,nsteps/2:].flatten()
		ninetyseventh.append(stats.scoreatpercentile(q_vals,97.5))
		twopointfifth.append(stats.scoreatpercentile(q_vals,2.5))
		i+=1
	plt.clf()
	print len(numpy.arange(0,25,5))
	print ninetyseventh
	plt.plot(numpy.arange(0,25,5)+data_num_stars,ninetyseventh,marker='x',ls='None',c='b')
	plt.plot(numpy.arange(0,25,5)+data_num_stars,twopointfifth,marker='x',ls='None',c='k')
	plt.axhline(init_flatparam)
	plt.xlabel('Number of Stars')
	plt.ylabel('q')
	plt.savefig('Q')
	plt.clf()

def LogRadPosErr():
	clean_sim =Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=init_flatparam, vcircsquared=init_vcircsquared,thicknessEigs = thickness_tensor_eigenvalues)
	clean_data = clean_sim.streams[0]

	model_sim = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[clean_data.fiducialPoint.pos], numForward=num_runs_sim[0],numBackward=num_runs_sim[1],potential=potential,flatparam=init_flatparam,vcircsquared=init_vcircsquared,thicknessEigs=thickness_tensor_eigenvalues)
	model = model_sim.streams[0]

#    def __init__(self,numForward=10,numBackward=10,timestep=0,fiducialPoint=None,numRemaining=-1,covarianceOfFiducialPoint=1,hasNoise=False,flatparam=1, vcircsquared=V2POTENTIAL,thicknessEigs=[.000001,0],ln_rad_pos_err=.1):

	noisy_data = []
	print 'data num stars',data_num_stars
	nums = [1,10,15,100]
	for i in xrange(4):
		stream = Simulator.DataColdStream(numForward=num_runs_sim[0],numBackward=num_runs_sim[1],timestep=timestep_for_integrator,fiducialPoint=clean_data.fiducialPoint.pos, hasNoise=True,flatparam=init_flatparam,vcircsquared=init_vcircsquared, thicknessEigs=thickness_tensor_eigenvalues,ln_rad_pos_err=.01*nums[i])
	#def __init__(self,numForward=10,numBackward=10,timestep=0,fiducialPoint=None,numRemaining=-1,covarianceOfFiducialPoint=1,hasNoise=False,flatparam=1, vcircsquared=V2POTENTIAL,thicknessEigs=[.000001,0],ln_rad_pos_err=.1):
	
#		noisy_sim = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[clean_data.fiducialPoint.pos], numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True, hasNoise=True, potential=Simulator.FLATLN, flatparam=init_flatparam, vcircsquared=init_vcircsquared, thicknessEigs = thickness_tensor_eigenvalues,ln_rad_pos_err=.01*nums[i])
#		stream = noisy_sim.streams[0]
#		noisy_sim.plotAllStreams(label='old'+str(i))
#		for star in stream.stars:
#			star.addNoise(.01*nums[i])
		noisy_data.append(stream)
#		noisy_sim.plotAllStreams(label='new'+str(i))
	random_state = numpy.random.get_state()
	chains = [MCMC.MonteCarloMarkovChain(random_state,datum,model, prior_info=prior_info,step_info=init_step_info,jump_in_plots=jump_in_plots) for datum in noisy_data]

	ninetyseventh=[]
	twopointfifth = []
	i = 0
	for chain in chains:
		print 'run', i
		chain.run(nsteps)
		print chain.chain[:,7,nsteps/2:].flatten().shape
		q_vals = chain.chain[:,7,nsteps/2:].flatten()
		ninetyseventh.append(stats.scoreatpercentile(q_vals,97.5))
		twopointfifth.append(stats.scoreatpercentile(q_vals,2.5))
		i+=1
	plt.clf()
	plt.plot(nums,ninetyseventh,marker='x',ls='None',c='b')
	plt.plot(nums,twopointfifth,marker='x',ls='None',c='k')
	plt.axhline(init_flatparam)
	plt.xlabel('Log radial velocity err')
	plt.ylabel('q')
	plt.savefig('log_pos_err')
	plt.clf()
	for chain in chains:
		mcmc1.plot_four_hists(7)
		mcmc1.plot_nonstream_vars(jump_in_plots=jump_in_plots,begin=mcmc1.chain.shape[2]/2)



def makeMcmcPlots(mcmcpickle,truepickle):
	print 'do'
	f1 = open(mcmcpickle,'r')
	mcmc1 = pickle.load(f1)
	mcmc1.debrine()
	f1.close()
	f1 = open(truepickle,'r')
	clean_stream = pickle.load(f1).streams[0]
	f1.close()
	mcmc1.plot_nonstream_vars_ensemble_at_t(-1)
	mcmc1.plot_stream_vars_ensemble_at_t(-1)
#	raise Exception("I poop");
	q_vals = mcmc1.chain[:,7,nsteps/2:].flatten()
	ninetyseventh =(stats.scoreatpercentile(q_vals,97.5))
	twopointfifth =(stats.scoreatpercentile(q_vals,2.5))
	sixteen = (stats.scoreatpercentile(q_vals,16))
	half = (stats.scoreatpercentile(q_vals,50))
	eightyfour = (stats.scoreatpercentile(q_vals,84))

	print '97.5', ninetyseventh
	print '84',eightyfour
	print '50', half
	print '16', sixteen
	print '2.5', twopointfifth

	f = open('data.dat','a')
	f.write(str(ninetyseventh)+'\n')
	f.write(str(eightyfour)+'\n')
	f.write(str(half)+'\n')
	f.write(str(sixteen)+'\n')
	f.write(str(twopointfifth)+'\n')
	f.close()
#	raise TypeError
	
	clean_stream =Simulator.Simulation(SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1], fiducialPoints=[clean_stream.fiducialPoint.pos], isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=clean_stream.flatparam, vcircsquared=clean_stream.vcircsquared,thicknessEigs = clean_stream.thicknessEigs).streams[0]
	print 'begin'
	print mcmc1.chain.shape
	print mcmc1.sampler.acceptance_fraction()
	print mcmc1.sampler.naccepted
	print mcmc1.chain
	the_accepted = numpy.zeros(mcmc1.chain.shape[0])
	x = 0
	for j in range(mcmc1.chain.shape[0]):
		for i in mcmc1.chain[j,3,:]:
			if i != x:
				the_accepted[j] += 1
				x = i
	print 1.*the_accepted/mcmc1.chain.shape[2]
	print 'boom sauce'
#	raise IOError

	stream1 = mcmc1.unpack_model_params(mcmc1.chain[0,:,0])
	
	print 'likelyhood of first stream', mcmc1.log_likelihood(mcmc1.data,stream1)
	print 'likelyhood of last stream', mcmc1.log_likelihood(mcmc1.data,mcmc1.model)


	print 'data velocity norm', numpy.linalg.norm(mcmc1.data.fiducialPoint.pos[3:])
	print 'first model velocity norm', numpy.linalg.norm(stream1.fiducialPoint.pos[3:])
	print 'modelvelocity norm', numpy.linalg.norm(mcmc1.model.fiducialPoint.pos[3:])
	print mcmc1.data.fiducialPoint.covarianceMatrix
	print mcmc1.data.fiducialPoint.thicknessTensor
	x =[]
	y = []
	z = []

	for i in mcmc1.data.stars:
		x.append(i.pos[3])
		y.append(i.pos[4])
		z.append(i.pos[5])
	print numpy.std(x)
	print numpy.std(y)
	print numpy.std(z)

#	mcmc1.plot_many_streams(jump_in_plots,clean_stream)	
#	mcmc1.mcmc_func_test()
#	mcmc1.plot_likelihood()
#	likelihoods = mcmc1.func_test()
	likelihoods = None
	for i in range(10):
#		mcmc1.plot_four_hists(i)
		mcmc1.plot_histogram(i,likelihoods,begin=mcmc1.chain.shape[2]/2)
#		mcmc1.plot_histogram(i,likelihoods,begin=-2)
		#	raise IOError
		#	mcmc1.plot_streams(len(mcmc1.chain),len(mcmc1.chain))
	for i in range(0,20,5):
		mcmc1.plot_stream_vars(jump_in_plots=jump_in_plots,begin=mcmc1.chain.shape[2]/2,label=str(i),walker=i)
		mcmc1.plot_nonstream_vars(jump_in_plots=jump_in_plots,begin=mcmc1.chain.shape[2]/2,label=str(i),walker=i)

	prefix = 't' + mcmc1.data.potentialToString(True) + 'm' + mcmc1.model.potentialToString(True)
	toPlot = []

#	for i in range(0,len(mcmc1.chain),len(mcmc1.chain)/jump_in_plots):
#		toPlot.append(mcmc1.chain[i][0][2])
#	S = Simulator.Simulation(SEEDVAL,numStreams=len(toPlot),timestep=timestep_for_integrator, fiducialPoints = toPlot, numForward=num_runs_sim[0], numBackward=num_runs_sim[1], numRemaining=data_num_stars)
#	S.plotAllStreams(mcmc1.data.plot(),mcmc1.data.potentialToString(),prefix)

def makeExamplePlots(cleanpickle, noisypickle):
	f = open(cleanpickle,'r')
	clean_sim = pickle.load(f)
	f.close()
	f = open(noisypickle,'r')
	noisy_sim = pickle.load(f)
	f.close()

	clean_sim.plotAllStreams(label='clean')
	noisy_sim.plotAllStreams(label='noisy')
	clean_sim.streams[0].integrate(numForward=5000,numBackward=5000)
	clean_sim.plotAllStreams(label='clean2',plotfid=True,fid=clean_sim.streams[0].fiducialPoint.pos)
	i=0
	for star in noisy_sim.streams[0].stars:
		Simulator.EllipseTest(10000,star,str(i))
		i+=1
	
	samplerate_vs_similarity(5)

def makeMultPlots(pickles):
	
	clean_streams = []
	noisy_streams = []
	for i in range(len(pickles)):
		f = open(pickles[i],'r')
		stream = pickle.load(f).streams[0]
		f.close()
		if i % 2==0:
			noisy_streams.append(stream)
		else:
			clean_streams.append(stream)
	Simulator.plotMultipleStreams(noisy_streams,clean_streams)


def random_rotate(vector):
    rotation = numpy.random.randn(3)
    rotation = rotation/numpy.linalg.norm(rotation)
    return rotation*numpy.linalg.norm(vector)


def Rotation(seed_number=1):
	clean_sim =Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=init_flatparam, vcircsquared=init_vcircsquared,thicknessEigs = thickness_tensor_eigenvalues)
	clean_data = clean_sim.streams[0]

	state = numpy.random.get_state()
	numpy.random.seed(seed_number)
	rot_pos = random_rotate(clean_data.fiducialPoint.pos[:3])
	rot_vel = random_rotate(clean_data.fiducialPoint.pos[3:])
	numpy.random.set_state(state)

	print "Old Fiducial Position: ", clean_data.fiducialPoint.pos
	fid_pos = []
	fid_pos.extend(rot_pos)
	fid_pos.extend(rot_vel)
	fid_pos = numpy.array(fid_pos)
	print "New Fiducial Position: ", fid_pos

	model = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[fid_pos],numForward=num_runs_sim[0],numBackward=num_runs_sim[1], potential=potential, flatparam=init_flatparam,vcircsquared=init_vcircsquared,thicknessEigs=thickness_tensor_eigenvalues).streams[0]    

	noisy_sim = Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[fid_pos], numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True, hasNoise=True, potential=Simulator.FLATLN, flatparam=init_flatparam, vcircsquared=init_vcircsquared, thicknessEigs = thickness_tensor_eigenvalues)

	noisy_data = noisy_sim.streams[0]

	clean_sim =Simulator.Simulation(seedVal=SEEDVAL,timestep=timestep_for_integrator,fiducialPoints=[fid_pos],numForward=num_runs_sim[0],numBackward=num_runs_sim[1],numRemaining=data_num_stars, isData=True,hasNoise=False,potential=Simulator.FLATLN,flatparam=init_flatparam, vcircsquared=init_vcircsquared,thicknessEigs = thickness_tensor_eigenvalues)

	f = open('clean_sim.pickle','w')
	pickle.dump(clean_sim,f)
	f.close()
	f = open('noisy_sim.pickle','w')
	pickle.dump(noisy_sim,f)
	f.close()

	mcmc = MCMC.MonteCarloMarkovChain(numpy.random.get_state(),noisy_data,model, prior_info=prior_info,step_info=init_step_info,initial_temperature=initial_temperature,jump_in_plots=jump_in_plots,annealing_parameter=annealing_parameter,annealing_dampening=annealing_dampening)

	mcmc.run(nsteps)


	logging.info("True Params:" + str(clean_data.get_params()))
	logging.info("Final Step:" + str(mcmc.chain[-1]))


	f = open("mcmc.pickle",'w')
	mcmc.brine()
	pickle.dump(mcmc,f)
	f.close()

	q_vals = mcmc.chain[:,7,nsteps/2:].flatten()
	ninetyseventh=(stats.scoreatpercentile(q_vals,97.5))
	twopointfifth=(stats.scoreatpercentile(q_vals,2.5))
	print '97.5', ninetyseventh
	print '2.5', twopointfifth

	

if __name__ == '__main__':
#    functionalTest2()
#    functionalTest1()
    main(1)

    
