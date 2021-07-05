from compare_coverage import * 
from model import * 

import numpy
from scipy.optimize import minimize
from scipy.stats import chi2
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})


def hist_zeta_mu(mu): 
  eta_profiled = profiled_eta(mu,n,0) 
  ntoys = 10000
  toy_n   = numpy.random.poisson(lamb(mu,eta_profiled),size=ntoys)
  toy_eta = numpy.random.normal(eta_profiled,1,size=ntoys)
  zetamu_dist = [zetamu(np,eta_p,mu) for np,eta_p in zip(toy_n,toy_eta)]
  #zetamu_dist.sort() 
  return zetamu_dist

 
fig, (ax1,ax2) = plt.subplots(1,2)

hist_zeta_mu 

chi2vals = chi2.pdf(numpy.arange(0.05,12.05,0.1),df=1)

ax1.hist(hist_zeta_mu(0.1), \
    bins=numpy.arange(0,12,0.1),range=(0,12),\
    color='black',histtype='step',density=True)
ax1.plot(numpy.arange(0.05,12.05,0.1),chi2vals,color='red',label='$\chi^{2}(1)$')
ax2.hist(hist_zeta_mu(9), \
    bins=numpy.arange(0,12,0.1),range=(0,12),\
    color='black',histtype='step',density=True)
ax2.plot(numpy.arange(0.05,12.05,0.1),chi2vals,color='red',label='$\chi^{2}(1)$')

ax1.set_xlabel("$\zeta_{0.1}$")
ax1.set_ylabel("$f(\zeta_{0.1}|H(0.1))$")
ax2.set_xlabel("$\zeta_{9}$")
ax2.set_ylabel("$f(\zeta_{9}|H(9))$")
ax1.set_yscale('log')
ax2.set_yscale('log')

ax1.legend()
plt.show()
