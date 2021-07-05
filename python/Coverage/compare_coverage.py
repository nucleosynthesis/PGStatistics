import numpy
from model import * 
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

# define the unconstrained and constrained functions 
def q(mu,eta,np,eta_p): 
  la = lamb(mu,eta)
  if np==0: return (eta-eta_p)*(eta-eta_p) + 2*la
  elif la<=0: return 99999
  else: return (eta-eta_p)*(eta-eta_p) + 2*la -2*np*numpy.log(la)
 
def q_unconstrained(x, args):
  mu, eta = x[0], x[1]
  np, eta_p = args[0], args[1]
  return q(mu,eta,np,eta_p)

def q_constrained(x, args):
  eta=x[0]
  mu, np, eta_p = args[0], args[1], args[2]
  return q(mu,eta,np,eta_p)

# define the minimisation routines 
def profiled_eta(mu, np, eta_p):
  init_params = [-3.0]
  bounds = [(-5,5)]
  res = minimize(q_constrained,init_params,args=[mu,np,eta_p],bounds=bounds)
  return res.x[0]

def global_min(np,eta_p): 
  init_params = [0.1,-3.]
  bounds = [(-1,50),(-5,5)]
  mle = minimize(q_unconstrained,init_params,args=[np,eta_p],bounds=bounds)
  return mle.fun,mle.x[0]
   
# calculate test statistic 
def zetamu(np,eta_p,mu):
  q_value        = q(mu,profiled_eta(mu,np,eta_p),np,eta_p)
  q_min,mu_min   = global_min(np,eta_p)
  if mu_min < 0     : return q_value-q(0,profiled_eta(0,np,eta_p),np,eta_p)
  else              : return q_value-q_min

def get_t683_neyman(mu):
  eta_profiled = profiled_eta(mu,n,0) 
  ntoys = 100000
  toy_n   = numpy.random.poisson(lamb(mu,eta_profiled),size=ntoys)
  toy_eta = numpy.random.normal(eta_profiled,1,size=ntoys)
  zetamu_dist = [zetamu(np,eta_p,mu) for np,eta_p in zip(toy_n,toy_eta)]
  zetamu_dist.sort() ; zeta_68 = zetamu_dist[int(ntoys*0.683)]
  return zeta_68

if __name__=='__main__':
    mu_range   = numpy.arange(0,12,0.2)
    n_coverage_toys = 10000

    cov_Neyman, cov_MINOS = [], []
    for mu_test in mu_range: 
       # Generate from the profiled value of the nuisance parameter 
       # we could also test the coverage under other values!
       eta_profiled = profiled_eta(mu_test,n,0)
       toy_n   = numpy.random.poisson(lamb(mu_test,eta_profiled),size=n_coverage_toys)
       toy_eta = numpy.random.normal(eta_profiled,1,size=n_coverage_toys)
       
       t683_Neyman = get_t683_neyman(mu_test) 
       t683_MINOS  = 1 
       zeta_obs_vals = [ zetamu(np,eta_p,mu_test) for np,eta_p in zip(toy_n,toy_eta) ]
       coverage_Neyman = float(len([x for x in zeta_obs_vals if x <= t683_Neyman ]))/n_coverage_toys
       coverage_MINOS  = float(len([x for x in zeta_obs_vals if x <= t683_MINOS  ]))/n_coverage_toys
       cov_Neyman.append(coverage_Neyman)
       cov_MINOS.append(coverage_MINOS)

    plt.plot(mu_range,cov_Neyman,color='red',marker="o",label="Neyman")
    plt.plot(mu_range,cov_MINOS,color='blue',marker="o",label="MINOS")
    plt.plot(mu_range,[0.683 for m in mu_range],color='black',linestyle="--")
    plt.xlabel("True value - $\\mu$")
    plt.ylabel("Coverage")
    plt.legend()
    plt.show()
