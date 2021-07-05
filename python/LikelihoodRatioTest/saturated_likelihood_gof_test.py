import numpy
from scipy.stats import poisson
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

# observed counts in each bin.
# norm is what we'll use to obtain the Poisson means
bin_boundaries = [0 , 0.2, 0.4, 0.6, 0.8, 1]
observations   =  [ 9,  6,  4,  5,  2 ]
norm = sum(observations)

# define a log function to handle 0 - will be 0 by the product
def mlog(o):
  if o==0: return 0
  else: return numpy.log(o)

# define the goodness of fit test likelihood ratio 
def neg2logLambda(observations,lambdas): 
  return -2*sum([o*mlog(l) - l - o*mlog(o) + o \
    for o,l in zip(observations,lambdas)])

# calculate the Poisson parameter for a given bin
# the bin is defined by the range x0->x1
a=2.0  
def mean_exp(x0,x1):
  N = a*(1-numpy.exp(-a))**-1
  return (N/a)*(numpy.exp(-a*x0)-numpy.exp(-a*x1))

lambdas = [norm*mean_exp(x,x+0.2)  for x in bin_boundaries[0:-1]]

test_stat_obs = neg2logLambda(observations,lambdas)

# generate the distributuon of the test stat under our hypothesis 
test_stat_toys = []
pval = 0 
for i in range(50000): 
  toy_observations = [numpy.random.poisson(l) for l in lambdas]
  test_stat_t = neg2logLambda(toy_observations,lambdas)
  test_stat_toys.append(test_stat_t)
  if test_stat_t > test_stat_obs : pval += 1

pval/=len(test_stat_toys)
print("-2ln Lambda = ",test_stat_obs, ", p-value = ", pval)

# plot the data, and the exponential hypothesis
xaxis = numpy.arange(0,1,0.01)
fig, (ax1,ax2) = plt.subplots(1,2)

ax1.plot([b+0.1 for b in bin_boundaries[0:-1]],observations,\
    marker="o",color="black", linestyle='none')
ax1.plot([b+0.1 for b in bin_boundaries[0:-1]],\
    lambdas, marker="_",color="red",linestyle='none',markersize=20)

ax2.hist([t for t in test_stat_toys if t > 1.3], \
    bins=numpy.arange(0,16,0.1),range=(0,16),color='cyan')
ax2.hist(test_stat_toys, bins=numpy.arange(0,16,0.1), \
    range=(0,16),color='black',histtype='step')
ax2.plot([test_stat_obs,test_stat_obs],[0,20],color='red')

ax1.set_xlabel("$x$"); ax1.set_ylabel("Observed/Expected events per bin")
ax2.set_xlabel("$-2\ln\Lambda$"); ax2.set_ylabel("# entries")
plt.show()
