import numpy as np
import matplotlib.pyplot as plt
# scipy has a lot of common pdfs for us to use
from scipy.stats import binom  

# parameters of the binomial probability distribution
N, eps=5000,0.7

def calc_pval():
  # for a random observation K, we calculate the p-val
  # this time, just use in the handy scipy cdf 
  K = np.random.binomial(N,eps)
  pval = binom.cdf(K,N,eps) 
  return pval

nTests=100000
pvals = [calc_pval() for test in range(nTests)]
y,binEdges = np.histogram(pvals,bins=10)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
menStd     = np.sqrt(y)
plt.bar(bincenters, y, width=0.025, yerr=menStd)
plt.xlabel("$p$-value")
plt.ylabel("Number of trials per bin")
plt.show()

