import numpy
from scipy.stats import poisson
from scipy.stats import chi2
import matplotlib.pyplot as plt

bin_boundaries = [0 , 0.2, 0.4, 0.6, 0.8, 1]
norm = 26

# define a log function to handle 0 - will be 0 by the product
def mlog(o):
  if o==0: return 0
  else: return numpy.log(o)

# define the goodness of fit test likelihood ratio 
def neg2logLambda(observations,lambdas): 
  return -2*sum([o*mlog(l) - l - o*mlog(o) + o for o,l in zip(observations,lambdas)])

# calculate the Poisson parameter for a given bin
# the bin is defined by the range x0->x1
a=2.0  
def mean_exp(x0,x1):
  N = a*(1-numpy.exp(-a))**-1
  return (N/a)*(numpy.exp(-a*x0)-numpy.exp(-a*x1))

lambdas = [norm*mean_exp(x,x+0.2)  for x in bin_boundaries[0:-1]]

# scale factor for total number of events
for sf,col in zip([1,3,50],['red','green','blue']):

  # generate the distributuon of the test stat under our hypothesis 
  test_stat_toys = []
  for i in range(10000): 
    toy_observations = [numpy.random.poisson(sf*l) for l in lambdas]
    test_stat_t = neg2logLambda(toy_observations,[sf*l for l in lambdas])
    test_stat_toys.append(test_stat_t)

  plt.hist(test_stat_toys, bins=20,range=(0,16),color=col, density=True,histtype='step')

# Compare to a chi2 function with 5 degrees of freedom
xaxis = numpy.arange(0.01,16,0.01)
chi2function = chi2(5)
plt.plot(xaxis,chi2function.pdf(xaxis),color='black')

plt.xlabel("$-2\ln\Lambda$")
plt.ylabel("$f(-2\ln\Lambda)$")
plt.show()
