import numpy
from scipy import stats
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

# x = cos(theta)
# pure photon exchange ~ 1+x^2
class photon_exchange(stats.rv_continuous):
  def _pdf(self,x) :
    if x < -1 : return 0 
    if x > 1  : return 0 
    N = 8./3
    return (1+x*x)/N
  def _get_support(self):
    return -1,1 

# Z* exchange ~ 1+x^2
class Z_exchange(stats.rv_continuous):
  def _pdf(self,x) :
    if x < -1 : return 0 
    if x > 1  : return 0 
    a = 0.4
    b = 0.6
    N = a*8./3
    return (a*(1+x*x)+b*x)/N
  def _get_support(self):
    return -1,1 

alte = Z_exchange()
null = photon_exchange()

# use scipy cosine p.d.f
random_alte = alte.rvs(size=5000)
random_null = null.rvs(size=5000)

# distribution for lambda under H0, and H1
lambda_H0 = [alte.pdf(x)/null.pdf(x) for x in random_null]
lambda_H1 = [alte.pdf(x)/null.pdf(x) for x in random_alte]

# plot the p.d.fs
fig, (ax1,ax2) = plt.subplots(1,2)
xaxis = numpy.arange(-1,1,0.01)
ax1.plot(xaxis,[alte.pdf(x) for x in xaxis],color='red')
ax1.plot(xaxis,[null.pdf(x) for x in xaxis],color='blue')

# and the distribution of Lambda
ax2.hist(lambda_H0, bins=50,range=(0,2),color='blue', density=True,histtype='step')
ax2.hist(lambda_H1, bins=50,range=(0,2),color='red',  density=True,histtype='step')

# label the plots
ax1.set_xlabel("$x$")
ax1.set_ylabel("f(x)")
ax2.set_xlabel("$\Lambda=f(x|H_{1})/f(x|H_{0})$")
ax2.set_ylabel("$f(\Lambda)$")
plt.show()
