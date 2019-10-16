import Similarity
import pickle
from Similarity import Simulator
from Simulator import logging
import scipy.stats as stats
import numpy
import matplotlib

import matplotlib.pyplot as plt
import markovpy

try:
    pass
#    plt.locator_params(tight=True,nbins=4)
except:
    print 'Warning: locator_params does not exist, this code may be in compatible with your current version of matplotlib'

SEEDVAL=30
label_font_size = 24
tick_font_size = 9

class MonteCarloMarkovChain(object):
    def __init__(self,random_state,data, model, prior_info=1.0, step_info=None,initial_temperature=1,jump_in_plots=1,chain_id=0,annealing_parameter=1,annealing_dampening=1,nwalkers=100):
        self.step_info = step_info
        self.prior_info = prior_info
        self.data = data
        self.model = model
        self.dt = self.model.dt
        self.jump_in_plots=jump_in_plots
        self.id = chain_id
        self.chain = []
        self.probs = []
	self.n_accept = numpy.zeros(7)
        self.nsteps = numpy.zeros(7)
        self.temperature = initial_temperature
        self.initial_temperature = initial_temperature
        self.annealing_parameter = annealing_parameter
        self.annealing_dampening = annealing_dampening
        self.nwalkers = nwalkers
        self.sampler = markovpy.EnsembleSampler(self.nwalkers,10,self.log_posterior_probability,a=2., outfile="dantheman.out")
        self.random = numpy.random.mtrand.RandomState()
        self.random.set_state(random_state)

    def brine(self):
        self.sampler = None
    
    def debrine(self):
        self.sampler = markovpy.EnsembleSampler(self.nwalkers,10,self.log_posterior_probability,a=2., outfile="dantheman.out")


    def generate_walkers(self,initial_guess,bounds):
        bounds = numpy.array(bounds)
        the_walkers = []
        for i in xrange(self.nwalkers):
            the_walkers.append((bounds[:,1]-bounds[:,0])*self.random.randn(bounds.shape[0]) + bounds[:,0])
        print the_walkers
        return the_walkers

    def pack_model_params(self,model):
        return model.get_params()

    def unpack_model_params(self,params):
        nf = int(params[0])-1
        nf /=2
        S = Simulator.Simulation(reseed=False,timestep=self.dt,fiducialPoints=[params[1:7]],numForward=nf,numBackward=nf,potential=self.model.potential,flatparam = params[7],vcircsquared=params[8],thicknessEigs=[params[9],0.0]).streams[0]

        return S

        

    def log_likelihood(self,data,model):
	# Likelihood : -.5*Similarity.similarity(streams[0], streams[1])
	# Similarity actually returns a Chi Squared
        sim = Similarity.similarity(data,model)
        logging.debug('log_likelihood: %.9f' % (-0.5*sim))
        return -.5*sim
        
    def log_prior(self,model):
        log_prior = -1*numpy.log(model.numForward+ model.numBackward+1)
        if model.flatparam < 0.5 or model.flatparam > 1.5:
            log_prior = log_prior+(-10.0**20)
        if model.thicknessEigs[0] <= (10**(-12)):
            log_prior += -10.0**20
            print 'this happened'
            print log_prior
	return log_prior

    def dummy_log_posterior_probability(self,model_params):#a_dummy_variable,model_params=0):
#        return -0.5*(a_dummy_variable[0])**2
        return (-0.5 * (model_params[3]-7.15)**2 / (0.02**2))

    def log_posterior_probability(self,model_params):
        model = self.unpack_model_params(model_params)
        logL = self.log_likelihood(self.data,model)
        logpr = self.log_prior(model)
        logpost = logpr + logL

        return logpost

    def anneal(self,run):
        new_temp= self.initial_temperature*self.annealing_parameter**(run/int(self.annealing_dampening))
        self.temperature=new_temp

    def func_test(self):
        params = self.pack_model_params(self.model)
        print params
        to_return = []
        likelihoods = []
        the_range = numpy.arange(numpy.log(params[0])-2, numpy.log(params[0])+4,.05)
        
        for value in the_range:
            new_params = [numpy.exp(value)]
            new_params.extend(params[1:])
            print new_params
            likelihoods.append(self.log_posterior_probability(numpy.array(new_params)))
            print 'gobama'
        plt.clf()
        plt.plot(the_range, likelihoods, '.')
        plt.xlabel('log(M)')
        plt.ylabel('log_posterior')
        plt.axis(xmin=the_range[0],xmax=the_range[-1])
        plt.savefig('logpost_model_param0.png')
        to_return.append([numpy.exp(the_range),likelihoods])
 
        for variable in range(1,10):
            likelihoods = []
            the_range = numpy.arange(params[variable]*.999, params[variable]*1.001,(params[variable]*1.001 - params[variable]*.999)/100)
            print len(the_range)
            for value in the_range:
                new_params = params[:variable]
                new_params.append(value)
                new_params.extend(params[variable+1:])
                likelihoods.append(self.log_posterior_probability(numpy.array(new_params)))

            plt.clf()
            plt.plot(the_range, likelihoods, '.')
            plt.xlabel('%i th model param' % variable)
            plt.ylabel('likelihood')
            plt.savefig('likelihood_model_param%i.png' % variable)

            to_return.append([the_range,likelihoods])
        return to_return

    def fake_run(self,nsteps,index):
        bounds = []
        params = numpy.array(self.pack_model_params(self.model))
        for var in params:
            bounds.append([var,var])
            
        if index != 0:
            bounds[index] = [bounds[index][0]*.99999,bounds[index][0]*1.00001]
        else:
            bounds[index] = [bounds[index][0]-1,bounds[index][0]+1]

        the_walkers = self.generate_walkers(params,bounds)
        old_sampler = self.sampler
        self.sampler = markovpy.EnsembleSampler(self.nwalkers,10,self.log_posterior_probability,a=2., outfile="dantheman.out")
        fixed_indicies = numpy.delete(numpy.arange(10),[index])
        fixed_vals = params[fixed_indicies]
        print fixed_indicies

        self.sampler.fix_parameters(fixed_indicies,fixed_vals)
        final_guess, prob, state = self.sampler.run_mcmc(the_walkers,self.random.get_state(),nsteps)
        self.chain = self.sampler.chain
        self.sampler = old_sampler

        
    
    def run(self,nsteps):
        params = numpy.array(self.pack_model_params(self.model))
        bounds = []

        
        for var in params:
           bounds.append([var*.999999,var*1.000001])


        the_walkers = self.generate_walkers(params, bounds)
        random_state =self.random.get_state()
        
        final_guess, prob, state = self.sampler.run_mcmc(the_walkers,self.random.get_state(),nsteps)


        print 'final_guesses',final_guess
        print 'probs',prob
        self.chain = self.sampler.chain
        print self.sampler.chain.shape
        print self.chain.shape
        

    def autocorelation(self,jump_in_plots=100):
        print len(self.chain), jump_in_plots
        for j in range(6):
            chain = []
            for i in range(0,len(self.chain),jump_in_plots):
                print i
                chain.append(self.chain[i][1+j])
            chain = numpy.array(chain)
            print 'chain',chain
            plot = []
            mu = numpy.mean(chain)

            for delta in range(len(chain)):
                stuff = numpy.zeros(len(chain)-delta)
                stuff[:] = (chain[:len(chain)-delta]-mu)*(chain[delta:]-mu)
                stuff = numpy.sum(stuff)/(1.0 *len(stuff))
                plot.append(stuff)
            plt.clf()
            plt.plot(numpy.arange(len(chain)),plot)
            plt.xlabel("$\delta$")
            plt.ylabel("$x_%i$ autocorrelation" % j)
            plt.savefig('autocorrelation%i.png' %j)
            plt.clf()

### I think this is outdated,but dont delete, it is still useful                
    def plot_var_vs_t(self,begin=0,end=-1,walker=0):
        prefix = 't'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        if end==-1: end = len(self.chain[0,0,:])-1
#        chain = []
#        for i in range(begin,end+1):
#            chain.extend(self.chain[:,:,i])
        chain = self.chain[walker]
        print chain.shape
#        raise IOError
#        chain = numpy.array(chain)
        prefix = 't'+self.data.potentialToString(True)+'m'+self.streams[1].potentialToString(True)
    
        for i in range(3):
            plt.clf()
            plt.plot(numpy.arange(begin,end+1,1), chain[i,:], mec='0.2', marker='.', ls='none', ms=0.5)
            scores = []
            for j in range(self.jump_in_plots):
                x_0=numpy.floor((end-begin)/self.jump_in_plots)*(j)
                x_1=numpy.floor((end-begin)/self.jump_in_plots)*(j+1)
                ninetyfifth = stats.scoreatpercentile(chain[i,x_0:x_1],95)
                yhi = stats.scoreatpercentile(chain[i,x_0:x_1],97.5)
                fifth = stats.scoreatpercentile(chain[i,x_0:x_1],5)
                ylo = stats.scoreatpercentile(chain[i,x_0:x_1],2.5)
                median = numpy.median(chain[i,x_0:x_1])
                #            print "fifth = %.4f, 95=%.4f" % (fifth,ninetyfifth)
                #            print "x1=%i,x2=%i" % (x_0,x_1)
                xmin = x_0/(end-begin)
                xmax = x_1/(end-begin)
                plt.axhline(fifth,xmin,xmax)
                plt.axhline(ninetyfifth,xmin,xmax)
            plt.axhline(self.data.fiducialPoint.pos[i],color='g',label='true value')
            #        print [0,nsteps,median-3*(median-ylo),median +3*(yhi-median)]
            #        plt.axis([0,nsteps,median-3*(median-ylo),median +3*(yhi-median)])
            plt.xlabel("Link Number",fontsize=label_font_size)
            plt.ylabel("$x_%i$ (kpc)" % (i+1),fontsize=label_font_size)
            plt.title("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString())
            plt.savefig(prefix+'mcmc-scatter_R%i-%i.png' % (i+1,self.id))

        for i in range(3,6):
            plt.clf()
            plt.plot(numpy.arange(begin,end+1,1), chain[i,:], mec='0.2', marker='.', ls='none', ms=0.5)
            scores = []
            for j in range(self.jump_in_plots):
                x_0=numpy.floor((end-begin)/self.jump_in_plots)*(j)
                x_1=numpy.floor((end-begin)/self.jump_in_plots)*(j+1)
                ninetyfifth = stats.scoreatpercentile(chain[i,x_0:x_1],95)
                fifth = stats.scoreatpercentile(chain[i,x_0:x_1],5)
                yhi = stats.scoreatpercentile(chain[i,x_0:x_1],97.5)
                ylo = stats.scoreatpercentile(chain[i,x_0:x_1],2.5)
                median = numpy.median(chain[i,x_0:x_1])
                plt.axhline(fifth,x_0/(end-begin),x_1/(end-begin))
                plt.axhline(ninetyfifth,x_0/self.nsteps,x_1/(end-begin))

            plt.axhline(self.data.fiducialPoint.pos[i],color='g',label='true value')
            #        plt.axis([0,nsteps,median-3*(median-ylo),median +3*(yhi-median)])
            plt.xlabel("Link Number")
            plt.ylabel("$x_%i$ (kpc/Myr)" % (i+1))
            plt.title("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString())
            plt.savefig(prefix+'mcmc-scatter_V%i-%i.png' % (i-3,self.id))
            plt.clf()
                    
### I think this is outdated. but dont delete. it is still useful
    def plot_scatter(self,begin=0,end=-1):
        prefix = 't'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        chain = []
        if end==-1: end = len(self.chain)-1
        for i in range(begin,end+1): chain.append(self.chain[i][3])
        chain = numpy.array(chain)
	
#        minProbIndex = self.probs.index(max(self.probs))
        self.plot_var_vs_t(begin,end)

        for i in range(3):
            plt.clf()
            plt.plot(chain[:,i], chain[:,i+3], mec='0.2', marker='.', ls='none', ms=0.5,label="mcmc step")
            plt.plot([self.data.fiducialPoint.pos[i]],[self.data.fiducialPoint.pos[i+3]],marker = 'x', label="input")
            plt.xlabel("$x_%i$ (kpc)" % (i+1))
            plt.ylabel("$x_%i$ (kpc/Myr)" % (i+3))
            plt.legend()
            plt.title("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString())
            plt.savefig(prefix+'mcmc-scatter_R%i(V%i)-%i.png' % (i+1,i+1,self.id))
            
        for i in range(3):
            for j in range(3):
                # Plot a scatter plot: R-component vs. V-component
                plt.clf()
                plt.plot(chain[:,j+3], chain[:,i+3], mec='0.2', marker='.', ls='none', ms=0.5,label='mcmc step')
                plt.plot([self.data.fiducialPoint.pos[j+3]],[self.data.fiducialPoint.pos[i+3]],marker='x',label="input")
		
                plt.xlabel("$x_%i$ (kpc/myr)" % (j+1))
                plt.ylabel("$x_%i$ (kpc/myr)" % (i+1))
                plt.legend()
                #		acc_frac = float(n_accept)/float(nsteps)
                plt.title("True: "+ self.data.potentialToString() + " Model: " + self.model.potentialToString())
                plt.savefig(prefix+'mcmc-scatter_V%ivsV%i-%i.png' % (i+1,j+1,self.id))
                plt.clf()
                
    def plot_stream_vars(self,begin=0,end=-1,label='',jump_in_plots=-1,walker=0):
        if jump_in_plots == -1:
            jump_in_plots = self.jump_in_plots
        fig = plt.figure(figsize=(15,15))
        fig.subplots_adjust(wspace=.2,hspace=.2)
        #fig.suptitle("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString())
        streamsPlot = []
        chain = []
        if end==-1: end = len(self.chain)-1
#        for i in range(begin,end+1,100): chain.append(self.chain[i][1:7].tolist())
        chain = self.chain[walker,1:7]
        chain = numpy.array(chain)


        ax = [[],[],[],[],[],[],[]]
        axbounds = numpy.zeros((6,6))
        axbounds = axbounds.tolist()
        for i in range(6):
            std = numpy.std(chain[i,:])
            ymin = numpy.amin(chain[i,:])
            ymin -= 3*std
            ymax = numpy.amax(chain[i,:])
            ymax += 3*std
            for j in range(i,6):
                std = numpy.std(chain[j,:])
                xmin = numpy.amin(chain[j,:])
                xmin -= 3*std
                xmax = numpy.amax(chain[j,:])
                xmax += 3*std
                axbounds[i][j] = [[xmin,xmax],[ymin,ymax]]
                axbounds[j][i] = [axbounds[i][j][1],axbounds[i][j][0]]

        axbounds = numpy.array(axbounds)


        for i in range(6):
            for j in range(6):
                ax[i].append(fig.add_subplot(6,6,j+1+6*i))
                ax[i][j].locator_params(tight=True,nbins=5)
                if i == j:
                    plt.plot(numpy.arange(chain.shape[-1]), chain[i,:], mec='0.2', marker='.', ls='none', ms=0.5)


                    for k in range(jump_in_plots):
                        x_0=numpy.floor((chain.shape[-1])/jump_in_plots)*(k)
                        x_1=numpy.floor((chain.shape[-1])/jump_in_plots)*(k+1)

                        ninetyfifth = stats.scoreatpercentile(chain[i,x_0:x_1],95)
                        yhi = stats.scoreatpercentile(chain[i,x_0:x_1],97.5)
                        fifth = stats.scoreatpercentile(chain[i,x_0:x_1],5)
                        ylo = stats.scoreatpercentile(chain[i,x_0:x_1],2.5)
                        median = numpy.median(chain[i,x_0:x_1])

                        xmin = x_0/(chain.shape[-1])
                        xmax = x_1/(chain.shape[-1])

                        plt.axhline(fifth,xmin,xmax)
                        plt.axhline(ninetyfifth,xmin,xmax)
                        
                    plt.axhline(self.data.fiducialPoint.pos[i],color='g',label='true value')

                    if i == 5:
                        plt.xlabel("Link Number",fontsize=label_font_size)
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)
                    else:
                        if i == 0:
                            plt.ylabel("$x_1$",fontsize=label_font_size)
                        else:
                            plt.setp(ax[i][j].get_yticklabels(),visible=False)

                        plt.setp(ax[i][j].get_xticklabels(),visible=False)
                        
                    plt.axis(ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])

                else:
                    plt.plot(chain[j,:],chain[i,:], marker='.',ms=0.5, mec='0.2',ls="none")
                    plt.axvline(self.data.fiducialPoint.pos[j],c='g')
                    plt.axhline(self.data.fiducialPoint.pos[i],c='g')
                    unitsx = ''
                    unitsy = ''

                    if j == 0:
                        plt.ylabel(("$x_%i$" %(i+1)) + unitsy,fontsize=label_font_size)
                        plt.setp(ax[i][j].get_yticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)

                    if i == 5:
                        plt.xlabel(("$x_%i$" % (j+1)),fontsize=label_font_size)
                        plt.setp(ax[i][j].get_xticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_xticklabels(), visible=False)
                        


                    plt.axis(xmin=axbounds[i,j][0][0],xmax=axbounds[i,j][0][1],ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])


                
        prefix = 't'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        plt.savefig(prefix+label+'mcmcPlotPos'+'.png', bbox_inches='tight')
        plt.close(fig)
        plt.clf()


    def plot_stream_vars_ensemble_at_t(self,index,label=''):

        fig = plt.figure(figsize=(10,10))
        fig.subplots_adjust(wspace=.1,hspace=.1)
        #fig.suptitle("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString())
        streamsPlot = []
        chain = []

        #        for i in range(begin,end+1,100): chain.append(self.chain[i][1:7].tolist())
        chain = self.chain[:,1:7,index].T
        print chain.shape
        chain = numpy.array(chain)


        ax = [[],[],[],[],[],[],[]]
        axbounds = numpy.zeros((6,6))
        axbounds = axbounds.tolist()
        for i in range(6):
            std = numpy.std(chain[i,:])
            ymin = numpy.amin(chain[i,:])
            ymin -= 3*std
            ymax = numpy.amax(chain[i,:])
            ymax += 3*std
            for j in range(i,6):
                std = numpy.std(chain[j,:])
                xmin = numpy.amin(chain[j,:])
                xmin -= 3*std
                xmax = numpy.amax(chain[j,:])
                xmax += 3*std
                axbounds[i][j] = [[xmin,xmax],[ymin,ymax]]
                axbounds[j][i] = [axbounds[i][j][1],axbounds[i][j][0]]

        axbounds = numpy.array(axbounds)


        for i in range(6):
            for j in range(6):
                ax[i].append(fig.add_subplot(6,6,j+1+6*i))
#                ax[i][j].locator_params(tight=True,nbins=4)
                if i == j:
                    plt.hist(chain[i,:])
                    plt.axvline(self.data.fiducialPoint.pos[i],color='g',label='true value')

                    if i == 5:
                        plt.xlabel("Link Number",fontsize=label_font_size-10)
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)
                        plt.setp(ax[i][j].get_xticklabels(),visible=False)
                        
                    plt.axis(ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])

                else:
                    plt.plot(chain[j,:],chain[i,:], marker='.',ms=0.5, mec='0.2',ls="none")
                    plt.axvline(self.data.fiducialPoint.pos[j],c='g')
                    plt.axhline(self.data.fiducialPoint.pos[i],c='g')
                    unitsx = ''
                    unitsy = ''

                    if j == 0:
                        plt.ylabel(("$x_%i$" %(i+1)) + unitsy,fontsize=label_font_size)
                        plt.setp(ax[i][j].get_yticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)

                    if i == 5:
                        plt.xlabel(("$x_%i$" % (j+1)),fontsize=label_font_size)
                        plt.setp(ax[i][j].get_xticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_xticklabels(), visible=False)
                        


                    plt.axis(xmin=axbounds[i,j][0][0],xmax=axbounds[i,j][0][1],ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])


                
        prefix = 'ensemble'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        plt.savefig(prefix+label+'mcmcPlotPos'+'.png', bbox_inches='tight')
        plt.close(fig)
        plt.clf()


    def plot_likelihood(self,begin=0,end=-1):
        print 'THIS FUNCTIONALITY IS BROKEN:plot_likelihood'
#          if end == -1:
#             end = len(self.chain)-1
#         plt.clf()

#         plt.plot(numpy.arange(begin,end),self.probs[begin:end],'k.',alpha=0.1)

#         plt.savefig('Similarity.png')
#         plt.clf()


    def plot_nonstream_vars(self,begin=0,end=-1,label='',jump_in_plots=-1,walker=0):
        if jump_in_plots == -1:
            jump_in_plots = self.jump_in_plots

        if end==-1: end = len(self.chain)-1

        chain = self.chain[walker,[0,9,8,7]]
        print chain.shape
        if chain.shape[0] != 4:
            print chain.shape
            raise Exception('dimension is wrong')


        chain = numpy.array(chain)

        fig = plt.figure(figsize=(10,10))
        fig.subplots_adjust(wspace=.2,hspace=.2)
        #fig.suptitle("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString(),fontsize=30)

        ax = [[],[],[],[],[]]
        axbounds = numpy.zeros((4,4))
        axbounds = axbounds.tolist()
        for i in range(4):
            std = numpy.std(chain[i,:])
            ymin = numpy.amin(chain[i,:])
            ymin -= 2*std
            ymax = numpy.amax(chain[i,:])
            ymax += 2*std
            for j in range(i,4):
                std = numpy.std(chain[j,:])
                xmin = numpy.amin(chain[j,:])
                xmin -= 2*std
                xmax = numpy.amax(chain[j,:])
                xmax += 2*std
                axbounds[i][j] = [[xmin,xmax],[ymin,ymax]]
                axbounds[j][i] = [axbounds[i][j][1],axbounds[i][j][0]]
        axbounds = numpy.array(axbounds)
        the_truth = numpy.array(self.data.get_params())[numpy.array([0,9,8,7])]
        for i in range(4):
            for j in range(4):
            	axis = fig.add_subplot(4,4,j+1+4*i)
            	# Force scientific notation
            	#axis.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            	#axis.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
                ax[i].append(axis)
                ax[i][j].locator_params(tight=True,nbins=5)
                if i == j:

                    plt.plot(numpy.arange(len(chain[0])), chain[i][:], mec='0.2', marker='.', ls='none', ms=0.5)
                    for k in range(jump_in_plots):
                        x_0=numpy.floor((len(chain[0]))/jump_in_plots)*(k)
                        x_1=numpy.floor((len(chain[0]))/jump_in_plots)*(k+1)
                        print x_0,x_1
                        ninetyfifth = stats.scoreatpercentile(chain[i,x_0:x_1],95)
                        yhi = stats.scoreatpercentile(chain[i,x_0:x_1],97.5)
                        fifth = stats.scoreatpercentile(chain[i,x_0:x_1],5)
                        ylo = stats.scoreatpercentile(chain[i,x_0:x_1],2.5)
                        median = numpy.median(chain[i,x_0:x_1])


                        xmin = x_0/(len(chain[0]))
                        xmax = x_1/(len(chain[0]))
                        print xmin, xmax
                        plt.axhline(fifth,xmin,xmax)
                        plt.axhline(ninetyfifth,xmin,xmax)
                    plt.axhline(the_truth[i],c='g')
                    plt.axis(ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])
#                    plt.axis(xmin=axbounds[i,j][0][0],xmax=axbounds[i,j][0][1],ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])
                    
                    if i == 0:
                        plt.ylabel("$M$", fontsize=label_font_size)
                        plt.setp(ax[i][j].get_yticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)
                    if i == 3:
                        plt.xlabel("Link number", fontsize=label_font_size-10)
                        plt.setp(ax[i][j].get_xticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_xticklabels(),visible=False)

                else:
                    plt.plot(chain[j][:],chain[i][:], marker='.',ms=0.5, mec='0.2',ls="none")
                    
                    if j == 0:
                        if i == 1:
                            plt.ylabel("$\mathbf{\Sigma}_x^2$",fontsize=label_font_size)
                        elif i == 2:
                            plt.ylabel("$V_c$",fontsize=label_font_size)
                        elif i == 3:
                            plt.ylabel("$q$",fontsize=label_font_size)
                        plt.setp(ax[i][j].get_yticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)

                    if i == 3:
                        if j == 0 :
                            plt.xlabel("$M$",fontsize=label_font_size)
                        elif j == 1:
                            plt.xlabel("$\Sigma_x^2$",fontsize=label_font_size)
                        elif j == 2:
                            plt.xlabel("$V_c^2$",fontsize=label_font_size)
                        elif j == 3:
                            plt.xlabel("Link Number",fontsize=label_font_size-10)
                        else:
                            plt.xlabel("$q$",fontsize=label_font_size)
                        plt.setp(ax[i][j].get_xticklabels(),fontsize=tick_font_size)
                    else:
                        plt.setp(ax[i][j].get_xticklabels(), visible=False)


                    plt.axvline(the_truth[j],color='g')
                    plt.axhline(the_truth[i],color='g')
                    plt.axis(xmin=axbounds[i,j][0][0],xmax=axbounds[i,j][0][1],ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])

        prefix = 't'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        plt.savefig(prefix+label+"mcmcPlotPot.png", bbox_inches='tight')
        plt.close()
        plt.clf()

    def plot_nonstream_vars_ensemble_at_t(self,index,label=''):
        chain = self.chain[:,[0,9,8,7],index].T
        
        if chain.shape[0] != 4:
            print chain.shape
            raise Exception('dimension is wrong')


        chain = numpy.array(chain)

        fig = plt.figure(figsize=(25,25))
        fig.subplots_adjust(wspace=.1,hspace=.1)
        #fig.suptitle("True: " + self.data.potentialToString() + " Model: " + self.model.potentialToString(),fontsize=30)

        ax = [[],[],[],[],[]]
        axbounds = numpy.zeros((4,4))
        axbounds = axbounds.tolist()
        for i in range(4):
            std = numpy.std(chain[i,:])
            ymin = numpy.amin(chain[i,:])
            ymin -= 2*std
            ymax = numpy.amax(chain[i,:])
            ymax += 2*std
            for j in range(i,4):
                std = numpy.std(chain[j,:])
                xmin = numpy.amin(chain[j,:])
                xmin -= 2*std
                xmax = numpy.amax(chain[j,:])
                xmax += 2*std
                axbounds[i][j] = [[xmin,xmax],[ymin,ymax]]
                axbounds[j][i] = [axbounds[i][j][1],axbounds[i][j][0]]
        axbounds = numpy.array(axbounds)
        the_truth = numpy.array(self.data.get_params())[numpy.array([0,9,8,7])]
        for i in range(4):
            for j in range(4):
            	axis = fig.add_subplot(4,4,j+1+4*i)
            	# Force scientific notation
            	#axis.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            	#axis.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
                ax[i].append(axis)
#                ax[i][j].locator_params(tight=True,nbins=4)
                if i == j:

                    plt.hist(chain[i][:])
                    plt.axvline(the_truth[i],c='g')
#                    plt.axis(ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])
#                    plt.axis(xmin=axbounds[i,j][0][0],xmax=axbounds[i,j][0][1],ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])
                    
                    if i == 0:
                        plt.setp(ax[i][j].get_yticklabels(),fontsize=12)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)
                    if i == 3:
                        plt.xlabel("Link number", fontsize=label_font_size-4)
                        plt.setp(ax[i][j].get_xticklabels(),fontsize=12)
                    else:
                        plt.setp(ax[i][j].get_xticklabels(),visible=False)

                else:
                    plt.plot(chain[j][:],chain[i][:], marker='.',ms=0.5, mec='0.2',ls="none")
                    
                    if j == 0:
                        if i == 1:
                            plt.ylabel("$\mathbf{\Sigma}_x^2$",fontsize=label_font_size)
                        elif i == 2:
                            plt.ylabel("$V_c$",fontsize=label_font_size)
                        elif i == 3:
                            plt.ylabel("$q$",fontsize=label_font_size)
                        plt.setp(ax[i][j].get_yticklabels(),fontsize=12)
                    else:
                        plt.setp(ax[i][j].get_yticklabels(),visible=False)

                    if i == 3:
                        if j == 0 :
                            plt.xlabel("$M$",fontsize=label_font_size)
                        elif j == 1:
                            plt.xlabel("$\Sigma_x^2$",fontsize=label_font_size)
                        elif j == 2:
                            plt.xlabel("$V_c^2$",fontsize=label_font_size)
                        elif j == 3:
                            plt.xlabel("Link Number",fontsize=label_font_size-4)
                        else:
                            plt.xlabel("$q$",fontsize=label_font_size)
                        plt.setp(ax[i][j].get_xticklabels(),fontsize=12)
                    else:
                        plt.setp(ax[i][j].get_xticklabels(), visible=False)


                    plt.axvline(the_truth[j],color='g')
                    plt.axhline(the_truth[i],color='g')
                    plt.axis(xmin=axbounds[i,j][0][0],xmax=axbounds[i,j][0][1],ymin=axbounds[i,j][1][0],ymax=axbounds[i,j][1][1])

        prefix = 'ensemble'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        plt.savefig(prefix+label+"mcmcPlotPot.png", bbox_inches='tight')
        plt.close()
        plt.clf()



    def plot_histogram(self,index,likelihoods=None,begin=0,end=-1):
        if end==-1: end = self.chain.shape[2]
        if begin ==-2: begin=self.chain.shape[2]-1
        chain = []

        for i in range(begin,end): chain.extend(self.chain[:,index,i])
        print len(chain)

        the_true_value = self.data.get_params()[index]
        
        chain = numpy.array(chain)
        print 'length of chain', len(chain)
        prefix = 't'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        
        mean_r = numpy.mean(chain)
        sigma_r = numpy.std(chain)
        print mean_r,sigma_r
        ninetysevenpointfifth = stats.scoreatpercentile(chain,97.5)
        twopointfifth = stats.scoreatpercentile(chain,2.5)
                      

        # Plot the histogram for R component

        plt.clf()
        fig = plt.figure(figsize=(10,10))
        plt.axvline(the_true_value,color='g')
        plt.axvline(ninetysevenpointfifth,color='b')
        plt.axvline(twopointfifth,color='r')
        
        the_stuff = plt.hist(chain, 1000)#,range=(mean_r - 2*sigma_r,mean_r+2*sigma_r))

        if likelihoods != None:
            likelihood=numpy.exp(numpy.array(likelihoods[index][1]))
            const = the_stuff[0].max()/likelihood.max()
            likelihood = const*likelihood
            plt.plot(likelihoods[index][0],likelihood,'r')

        if index == 0:
            plt.xlabel("M",fontsize=label_font_size)
        if index in range(1,7):
            plt.xlabel("$x_%i$ (kpc)" % (index))
        if index == 7:
            plt.xlabel("q",fontsize=label_font_size)
            plt.axis(xmin=.5,xmax=1.5)
        if index == 8:
            plt.xlabel("$V_c^2$ $(kpc/myr)^2$",fontsize=label_font_size)
        if index == 9:
            plt.xlabel("$\Sigma_x$ $(kpc^2)$",fontsize=label_font_size)
        plt.tick_params(labelsize=tick_font_size+10)
#        plt.axis(xmin=mean_r - 2*sigma_r, xmax=mean_r + 2*sigma_r)
        plt.savefig(prefix+'mcmc-hist-R-%i.png' % (index))
        plt.clf()


    def plot_four_hists(self,index):
        end = self.chain.shape[2]
        half = end/2


        fig = plt.figure(figsize=(25,25))
        
        fig.subplots_adjust(wspace=.1,hspace=.1)
    
        if index == 0:
            fig.suptitle("M")
        if index in range(1,7):
            fig.suptitle("$x_%i$ (kpc)" % (index))
        if index == 7:
            fig.suptitle("q")
        if index == 8:
            fig.suptitle("$V_c^2$ (kpc/myr)^2")
        if index == 9:
            fig.suptitle("$\Sigma_x$ (kpc^2)")
        
        the_true_value = self.data.get_params()[index]

        ## Full histogram 1 walker
        axis1 = fig.add_subplot(2,2,1)

        chain = self.chain[0,index,:]

        plt.axvline(the_true_value,color='g')
        the_stuff = plt.hist(chain, 1000)

        ## last half histogram all walkers every 100th step
        axis = fig.add_subplot(2,2,2,sharex=axis1,sharey=axis1)
        chain = []
        for i in range(half,end):chain.extend(self.chain[:,index,i])
        chain = numpy.array(chain)


        plt.axvline(the_true_value,color='g')
        the_stuff = plt.hist(chain,1000,weights=.01*numpy.ones(len(chain)))

        ## first half histogram 1 walker
        axis = fig.add_subplot(2,2,3,sharex=axis1,sharey=axis1)
        chain = self.chain[0,index,:half]
        plt.axvline(the_true_value,color='g')
        the_stuff = plt.hist(chain, 1000)

        ## second half histogram 1 walker
        axis = fig.add_subplot(2,2,4,sharex=axis1,sharey=axis1)
        chain = self.chain[0,index,half:end]
        plt.axvline(the_true_value,color='g')
        the_stuff = plt.hist(chain, 1000)



        plt.savefig('mcmc-hist-R-%i.png' % (index))
        plt.clf()
        
    def plot_many_streams(self,num_plots,true_stream=None, label=''):
        print self.step_info


        againstPotential = self.data.potentialToString()
        fig = plt.figure(figsize=(22,22))
        fig.subplots_adjust(wspace=.1,hspace=.1)
#        print len(self.streams)
        stream_potential_name = 'FlatLn'#self.streams[0].potentialToString()
#        if toPlotAgainst==[]:
#            fig.suptitle("True:" + stream_potential_name, fontsize = 30)
#        else:
        #fig.suptitle("True: " + againstPotential + " Model: " + stream_potential_name, fontsize = 30)
    
        streamsPlot = []
#        indicies = numpy.random.randint(0,len(self.chain),num_plots)
        indicies = [0]
        original_model = self.unpack_model_params(self.chain[0,:,0])

#        print 'modelfidpt-truefidpt', 10*(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
        print 'modelfidpt-truefidpt', (original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)

#        print 'truefidpt', true_stream.fiducialPoint.pos
#        raise IOError

        print 'modelfidpt-truefidpt', (original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
        observerpos = Simulator.EARTH.copy()#numpy.zeros(3)
#        observerpos[:2] = original_model.fiducialPoint.pos[:2]
        print 'obp',observerpos
        print 'omfdp',original_model.fiducialPoint.pos
        r_hat = original_model.fiducialPoint.pos[:3].copy() - observerpos[:3]#- Simulator.EARTH[:3].copy()
        print 'modelfidpt-truefidpt', (original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
#        r_hat /= numpy.linalg.norm(r_hat)
        theta_hat = original_model.fiducialPoint.pos[3:].copy()
        print 'othetahat',theta_hat
#        theta_hat /= numpy.linalg.norm(theta_hat)
#        print 'modelfidpt-truefidpt', 10*(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
#        r_hat = original_model.fiducialPoint.pos[:3] - Simulator.EARTH[:3]
#        print 'modelfidpt-truefidpt', 10*(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
#        r_hat = r_hat - numpy.dot(theta_hat,r_hat)*theta_hat
        print 'r_hat',r_hat
        r_hat /= numpy.linalg.norm(r_hat)

        theta_hat = theta_hat - numpy.dot(theta_hat,r_hat)*r_hat
#        print 'modelfidpt-truefidpt', 10*(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
        theta_hat /= numpy.linalg.norm(theta_hat)
        print 'theta_hat', theta_hat
#        phi_hat = numpy.cross(theta_hat,r_hat)
        phi_hat = numpy.cross(r_hat,theta_hat)
        if numpy.abs(1-numpy.linalg.norm(phi_hat)) >= 10**(-15):
            raise Exception('cross products foobar')

        zeros = numpy.zeros(3)


        e1 = numpy.append(r_hat,zeros)
        e2 = numpy.append(theta_hat,zeros)
        e3 = numpy.append(phi_hat,zeros)
        e4 = numpy.append(zeros,r_hat)
        e5 = numpy.append(zeros,theta_hat)
        e6 = numpy.append(zeros,phi_hat)
#        print e6
        change_basis = numpy.array([e1,e2,e3,e4,e5,e6])
        print numpy.dot(change_basis,numpy.dot(self.data.fiducialPoint.covarianceMatrix,change_basis.T))
        #raise IOError
        for i in self.data.stars:
            print numpy.diag(numpy.dot(change_basis,numpy.dot(i.covarianceMatrix,change_basis.T)))
            
#        print 'modelfidpt-truefidpt', 10*(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)

        #change_basis = numpy.eye(6)
#        print 'change basis', change_basis
#        indicies = [0]
#        print change_basi
#        print 'index',indicies

        for index in indicies:
            model = self.unpack_model_params(self.chain[index,:,0])
            print model.stars[0].pos
#            model_plot = numpy.array(model.plot())
#            print 'modelplotsize',model_plot.shape

            changed_model_plot = numpy.dot(change_basis,(numpy.array(model.plot()).T-observerpos).T)
            print 'model.stars[0].pos',model.stars[0].pos

#            print 'len(changed_model_plot[0])', changed_model_plot.shape
            streamsPlot.append(changed_model_plot)

        changed_model_plot = numpy.array(changed_model_plot)

#        print 'real_fid_pt',changed_model_plot[:,self.data.fiducialPointIndex]
#        print 'modelfidpt-truefidpt', 10*(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)
        print 'original_modelfidptchaanged', numpy.dot(change_basis,original_model.fiducialPoint.pos-observerpos)
#        print 'cmp',changed_model_plot[0,:]
        print 'cmp.shape', changed_model_plot.shape
        print 'modelfidpt-truefidpt',(original_model.fiducialPoint.pos-true_stream.fiducialPoint.pos)

#            print 'modelfiducialpoint', model.fiducialPoint.pos
#        print numpy.dot(change_basis,numpy.array(model_plot)[:,0])
#        print 'x_5',changed_model_plot[4,:]
#        for i in range(0, len(self.chain),len(self.chain)/num_plots):
#            streamsPlot.append(self.unpack_model_params(self.chain[i][0]).plot())
        print 'len(streamsPlot)',len(streamsPlot)

        print 'len(chain)',len(self.chain)
        original_model_plot = numpy.dot(change_basis,(numpy.array(original_model.plot()).T - observerpos).T)
#        print len(self.chain)/num_plots
        o_fid_pos = numpy.dot(change_basis,original_model.fiducialPoint.pos)
        t_fid_pos = numpy.dot(change_basis,true_stream.fiducialPoint.pos)
        print 'ofid-tfid', o_fid_pos - t_fid_pos
        print 'A(o -t)', numpy.dot(change_basis,original_model.fiducialPoint.pos - true_stream.fiducialPoint.pos) 
        truth = []
        if true_stream != None:
#            print true_stream.stars[0].pos
            truth = numpy.dot(change_basis,(numpy.array(true_stream.plot()).T-observerpos).T)
            print 'truestreamfidpt',numpy.dot(change_basis,true_stream.fiducialPoint.pos-observerpos)
#        streamsPlot = numpy.array(streamsPlot)
#        raise IOError
        toPlotAgainst = numpy.dot(change_basis,(numpy.array(self.data.plot()).T - observerpos).T)
        fullPlot = [[],[],[],[],[],[]]
        for stream in streamsPlot:
            for i in range(6):
                fullPlot[i].extend(stream[i])
        for i in range(6):
            fullPlot[i].extend(toPlotAgainst[i])
            fullPlot[i].extend(truth[i])
            fullPlot[i].extend(original_model_plot[i])
        fullPlot= numpy.array(fullPlot)
        axbounds = numpy.zeros((6,6)).tolist()
        centers = []
        halfwidths = []
        for i in range(6):
                xmin = numpy.amin(fullPlot[i,:])
                xmax = numpy.amax(fullPlot[i,:])
                center = .5*(xmin+xmax)
                halfwidth = .5*(xmax-xmin)
                centers.append(center)
                halfwidths.append(halfwidth)
        halfwidthpos = numpy.array(halfwidths)[:3].max()
#        print 'hwp',halfwidthpos
        halfwidthvel = numpy.array(halfwidths)[3:].max()
        halfwidths = [halfwidthpos,halfwidthpos,halfwidthpos,halfwidthvel,halfwidthvel,halfwidthvel]
        for i in range(6):
            for j in range(6):
                ax = fig.add_subplot(6,6,j+1+6*i)
#                if len(toPlotAgainst) > 0:

                plt.plot(toPlotAgainst[j],toPlotAgainst[i], 'g.', alpha=.40)
                plt.plot(original_model_plot[j],original_model_plot[i],'bx',alpha=.1)
#                plt.plot(truth[j],truth[i],'r.',alpha=.1)

                plt.plot
                if len(truth)>0:
#                    pass
                    plt.plot(truth[j],truth[i],'r.',alpha=.1)

                for stream in streamsPlot:
#                    pass
#                    if i == 4:
 #                       print len(stream[j])
#                        print stream[j]
                    plt.plot(stream[j],stream[i], 'k.',alpha=.1)

#                 else:
#                     print 'Error'
#                     print toPlotAgainst
#                     raise IOError
#                     for stream in streamsPlot:
#                         plt.plot(stream[j],stream[i],'r.',alpha=1)


                
#                 xmin = numpy.amin(fullPlot[j,:])
#                 xmax = numpy.amax(fullPlot[j,:])
#                 ymin = numpy.amin(fullPlot[i,:])
#                 ymax = numpy.amax(fullPlot[i,:])

#                xmin -= numpy.abs(xmin*.01)
#                xmax += numpy.abs(xmax*.01)
#                ymin -= numpy.abs(ymin*.01)
#                ymax += numpy.abs(ymax*.01)

                plt.axis([centers[j]-halfwidths[j],centers[j]+halfwidths[j],centers[i]-halfwidths[i],centers[i]+halfwidths[i]])
                
                if j == 0:
                    unitsy = ''
#                     if i <3:
#                         unitsy = '$\text{(kpc)}$'
#                     else:
#                         unitsy = '$\text{(kpc/Myr)}'
                    plt.ylabel(('$x_%i$ ' % (i+1))+unitsy,fontsize = label_font_size)
                    plt.setp(ax.get_yticklabels(),fontsize=12)
                else:
                    plt.setp(ax.get_yticklabels(),visible=False) 

                if i == 5:
                    unitsx=''
#                     if j < 3:
#                         unitsx = '$\text{(kpc)}$'
#                     else:
#                        unitsx = '$\text{(kpc/Myr)}$'
                        
                    plt.xlabel(('$x_%i$ ' % (j+1))+unitsx,fontsize=label_font_size)
                    plt.setp(ax.get_xticklabels(),fontsize=12)
                else:
                    plt.setp(ax.get_xticklabels(),visible=False)
                
                if i == 5 and j == 5:
                    plt.setp(ax.get_yticklabels(),visible=False)


                
#        plt.show()
        plt.savefig(label+'StreamsPlot'+'.png')
        plt.close(fig)
        plt.clf()


    def plot_streams(self,fname,title):
        plt.clf()
        stream1Plot = self.data.plot()[:3]
        stream2Plot = self.model.plot()[:3]
        fname = str(fname) + '-' + str(self.id)
        prefix = 't'+self.data.potentialToString(True)+'m'+self.model.potentialToString(True)
        plt.plot(stream1Plot[0],stream1Plot[1], 'o', label='data')
        plt.plot(stream2Plot[0],stream2Plot[1],'.',label='model')
        plt.plot([stream2Plot[0][self.model.fiducialPointIndex]],[stream2Plot[1][self.model.fiducialPointIndex]],'x')
        plt.xlabel('$x_1$ (kpc)')
        plt.ylabel('$x_2$ (kpc)')
        #    plt.axis([-100,100,-100,100])
        #    plt.legend(("data","model"))
        plt.legend()
        plt.title("True: "+ self.data.potentialToString() + " Model: " + self.model.potentialToString())
        plt.savefig(prefix+'XY-'+str(title)+'-'+str(fname) + '.png')
        plt.clf()


        plt.plot(stream1Plot[1],stream1Plot[2],'o',label='data')
        plt.plot(stream2Plot[1],stream2Plot[2],'.',label='model')
        plt.xlabel('$x_2$ (kpc)')
        plt.ylabel('$x_3$ (kpc)')
        #    plt.legend(("data","model"))
        plt.legend()

        plt.savefig(prefix+'YZ-'+str(title)+'-'+str(fname) + '.png')
        plt.clf()
        

        plt.plot(stream1Plot[2],stream1Plot[0],'o',label='data')
        plt.plot(stream2Plot[2],stream2Plot[0],'.',label='model')
        plt.xlabel('$x_3$ (kpc)')
        plt.ylabel('$x_1$ (kpc)')
        #    plt.legend(("data","model"))
        plt.legend()
        plt.title("True: "+ self.data.potentialToString() + " Model: " + self.model.potentialToString())

        plt.savefig(prefix+'ZX-'+str(title)+'-'+str(fname) + '.png')
        plt.clf()

    def mcmc_func_test(self):
        fig = plt.figure(figsize=(15,15))
        fig.subplots_adjust(wspace=.1,hspace=.1)
        fig.add_subplot(111)

        
        streamsPlot = [[],[]]
#        indicies = numpy.random.randint(0,len(self.chain),10)
#        for index in indicies:
        for link in self.chain:
            streamsPlot[0].append(numpy.sqrt(link[3]))
            streamsPlot[1].append(link[1][3])
#            model = self.unpack_model_params(self.chain[index][0])
#            vels = model.plot()[3]
#            vels = [model.fiducialPoint.pos[3]]
#            streamsPlot.append([numpy.sqrt(model.vcircsquared), vels])

#            streamsPlot.append([model.vcircsquared,model.fiducialPoint.pos[3]])
            #print 'modelfiducialpoint', model.fiducialPoint.pos



        toPlotAgainst = self.data.fiducialPoint.pos[3]
        print len(streamsPlot)
#        print len(toPlotAgainst)
#        print streamsPlot[0][1]
#        print streamsPlot[0][0]
        print toPlotAgainst
        
        plt.plot([numpy.sqrt(self.data.vcircsquared)],[toPlotAgainst],'rx',alpha=1.0)
#        for stream in streamsPlot:
        plt.plot(streamsPlot[0],streamsPlot[1],'.',alpha=.5)
        plt.axhline(toPlotAgainst,color='g')
        plt.axvline(numpy.sqrt(self.data.vcircsquared),color='g')
        plt.xlabel('$V_c$', fontsize=label_font_size)
        plt.ylabel('$x_3$',fontsize=label_font_size)
        plt.savefig('vc_v1.png')

        plt.close(fig)


#        def __init__(self,seedVal=20,numStreams=1,timestep=1,potential=LN, fiducialPoints=None,numForward=10,numBackward=10,integrationScheme='forwardEuler',numRemaining=-1,isData=False,hasNoise=False):


#this is broken. fix.
def mcmcTest():
    data = Simulator.Simulation(SEEDVAL,numForward=100,numBackward=100,numRemaining=10,isData=True,hasNoise=False,).streams[0]
    model = Simulator.Simulation(SEEDVAL,numForward=100,fiducialPoints=[data.fiducialPoint.pos],numBackward=100,isData=False).streams[0]

    streams = [data,model]
    prior_info=1.0
    step_info= [numpy.array([.45,.45,.45,.0035,.0035,.0035]),None,None,1,0]

    initial_temperature = 1

    mcmc = MonteCarloMarkovChain(streams, prior_info, step_info, initial_temperature,jump_in_plots=10)

    mcmc.run(5000, plot_streams_on=True)
    mcmc.plot_scatter()

    


def DiagnosticTest():
    logging.basicConfig(filename='MCMC_LOG',level=logging.DEBUG, format='%(message)s')

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    format = logging.Formatter('%(message)s')
    console.setFormatter(format)
    logging.getLogger('').addHandler(console)
    mcmcTest()
    



if __name__ == '__main__':
    DiagnosticTest()


