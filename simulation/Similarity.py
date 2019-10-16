import Simulator
from Simulator import numpy
from Simulator import logging
from Simulator import matplotlib
from Simulator import plt

DBL_MIN = 2.2250738585072014e-308
DBL_MAX = 1.7976931348623157e+308
INDICIES = [[],[]]
DISTSMIN = []
VALSMIN = []
MINDISTS = []
MINVALS = []
def Metric(star,thicknessEigs):
    if star.covarianceMatrix == None:
        raise Exception("MetricError: This star is a model point")
    modelThicknessTensor = numpy.diag([thicknessEigs[0],thicknessEigs[0],thicknessEigs[0],thicknessEigs[1],thicknessEigs[1],thicknessEigs[1]])
    return numpy.linalg.inv(star.covarianceMatrix + modelThicknessTensor)

def chiSquaredVector(u,v,metricTensor):
    x = v - u
    return numpy.sum(x*numpy.dot(x,metricTensor),axis=1)

def chiSquaredScalar(u,v,metricTensor):
    x = v - u
    return numpy.dot(x,numpy.dot(metricTensor,x))

# assumes all model points have the same associated dt
def marginalizedChiSquared(data_point, model, metricTensor):
    values = chiSquaredVector(data_point.pos,model.stars_pos,metricTensor)
    underflow_const = numpy.array([numpy.log(DBL_MIN) + .5*values.min(),numpy.log(DBL_MAX) +.5*values.max() - numpy.log(len(values))]).min() 
    ## note to self: the negative sign in front of the log of the norm of the tensor is due to the fact that the tensor is the inverse of the covariance tensor
    ret = -2*(numpy.log(numpy.sum(numpy.exp(-.5*values+underflow_const))) - underflow_const- numpy.log(len(values))) - numpy.log(numpy.linalg.det(metricTensor))
    if numpy.isnan(ret):
        if numpy.linalg.det(metricTensor) < 0:
            ret = -numpy.inf
        else:
            print 'ret',ret
            print "log_sum with underflow_const: ", numpy.log(numpy.sum(numpy.exp(-.5*values + underflow_const)))
            print "underflow_const: ", underflow_const
            print 'model.dt: ', model.dt
            print 'det:', numpy.linalg.det(metricTensor)
            print 'ln(det): ', numpy.log(numpy.linalg.det(metricTensor))
            print 'metric tensor', metricTensor
            print 'star covariance', data_point.covarianceMatrix
            print 'star position', data_point.pos
            print 'model thickness', model.thicknessEigs
            print 'model params', model.get_params()
            raise Exception('MarginalizeChiSquared:underflow')

    return ret

def similarity(data, model):
#    if model.thicknessEigs[0] < 10**(-12):
#            raise Exception('Similarity: thickness too small %f' % model.thicknessEigs[0])

    sum = 0
    i = 0
    plt.clf()
    log = numpy.log
    det = numpy.linalg.det
    for star in data.stars:
        tensor = Metric(star,model.thicknessEigs)
        
        mchi2 = marginalizedChiSquared(star, model,tensor)
        i+=1
        if numpy.isnan(mchi2):
            print star
            print model
            print tensor
            print to_add_0
            raise Exception('MarginalizedChiSquared: more underflow. this should never happen.')
        sum+=mchi2
    return sum


def similarityTest(hasNoise):
    sim1 = Simulator.Simulation(numStreams=1,timestep=1,numForward=100,numBackward=100,numRemaining=10,isData=True,hasNoise=hasNoise)
    stream1 = sim1.streams[0]
    sim2 = Simulator.Simulation(numStreams=1,timestep=1,numForward=100,numBackward=100)
    stream2 = sim2.streams[0]
    sim3 = Simulator.Simulation(numStreams=1,timestep=1,numForward=100,numBackward=100)
    stream3 = sim3.streams[0]

    sim1.plotStream(0,0)
    sim2.plotStream(0,1)

    passed = False
    similarityCase1 = similarity(stream1,stream2,forceEvaluate=True)
    print "Similarity, forcing evaluation: ", similarityCase1
    similarityCase2 = similarity(stream1,stream3)
    print "Similarity, allowing length variation: ", similarityCase2

    if hasNoise:
        if similarityCase1 <= 60 and similarityCase2 <= 60:
            passed = True
    else:
        if similarityCase1 == 0 and similarityCase2 == 0:
            passed = True
    print 'Passed: '+str(passed) + '\n\n'

def Diagnostic():
    print "Running Diagnostic Test....\n\n"
    print "Testing Similairity:\n"
    print "Fitting truth to clean data:\n"
    similarityTest(False)
    print "Fitting truth to noisy data:\n"
    similarityTest(True)


if __name__=="__main__":
    Diagnostic()
