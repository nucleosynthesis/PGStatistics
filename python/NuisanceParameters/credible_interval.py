import numpy 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
plt.rcParams.update({'font.size': 14})

from model import *
import marginalised_likelihood

# find and plot the 68.3\% interval using the `shortest' algorithm
# prior on mu P(mu)
def prior_mu_new(mu):
  if (mu > 50 or mu < 0) : return 0
  else: return 1./50

marginalised_likelihood.prior_mu = prior_mu_new

# could also do the full integral but its slower
def integral_range(mua,mub,n): 
  return integrate.quad(pmu,mua,mub,args=(n))[0]

xaxis = numpy.linspace(-1,20,5000)
normalisation = marginalised_likelihood.norm(n) 
posterior_dist = [ marginalised_likelihood.integral(n,mu)/normalisation for mu in xaxis ]

intervals=[]
for i in range(len(xaxis)):
 x = xaxis[i]
 if x > 8: break
 inte=0
 for j in range(i,len(xaxis)-1):
   y = xaxis[j+1]
   yl = xaxis[j]
   inte += posterior_dist[j]*(y-yl)
   #inte = integral_range(x,y,n)/norm(n)
   if inte >= 0.683: 
    intervals.append([y-x,[x,y],[i,j],inte])
    break
for interval in intervals: print(interval)
# find the smallest one
intervals.sort(); interval = intervals[0]
print("68.3%% interval (%.2f,%.2f)"%(interval[1][0],interval[1][1]))

# plot as a shaded region
i=interval[2][0]; j=interval[2][1]
plt.plot(xaxis,posterior_dist,color='blue')
plt.fill_between(xaxis[i:j],posterior_dist[i:j])
plt.xlabel("$\\mu$")
plt.ylabel("$P(\mu)$")
plt.show()
