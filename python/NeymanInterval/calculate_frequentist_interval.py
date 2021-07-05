import numpy
from model import * 
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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

# return a set of heights for each true value of mu 
# for each value, we would mark on where 68.3 % of the distribution lives 
def histo_zetamu(mu,zeta_obs):
  # find the best (profiled) nuisance parameter values for the data (n,0)
  eta_profiled = profiled_eta(mu,n,0) 
  ntoys = 5000
  toy_n   = numpy.random.poisson(lamb(mu,eta_profiled),size=ntoys)
  toy_eta = numpy.random.normal(eta_profiled,1,size=ntoys)
  zetamu_dist = [zetamu(np,eta_p,mu) for np,eta_p in zip(toy_n,toy_eta)]
  zetamu_dist.sort() ; zeta_68 = zetamu_dist[int(ntoys*0.68)]
  return zetamu_dist, zeta_68

mu_range   = numpy.arange(0,10,0.1)
zeta_range = numpy.arange(0,10,0.1)
zeta_obs_vals = []
zeta_68_vals  = []
mu_interval = []
densities = []
for mu_test in mu_range: 
  zeta_obs = zetamu(n,0,mu_test)
  zetamu_toys,zeta_68 = histo_zetamu(mu_test,zeta_obs)
  density_vals = plt.hist(zetamu_toys,density=True,bins=zeta_range)
  densities.append(density_vals[0])
  zeta_68_vals.append(zeta_68)
  zeta_obs_vals.append(zeta_obs)
  if zeta_obs < zeta_68 : mu_interval.append(mu_test)

# Print out the results
mu_l, mu_u = min(mu_interval),max(mu_interval)
mu_hat = global_min(n,0)[1]
print("interval -> (%.2f,%.2f)"%(mu_l,mu_u),", mu = %.2f + %.2f -%.2f "%(mu_hat,mu_u-mu_hat,mu_hat-mu_l))

plt.clf()
X,Y = numpy.meshgrid(zeta_range,mu_range) 
c = plt.pcolor(X,Y,densities, \
    norm=LogNorm(vmin=0.001, vmax=2.0))
plt.colorbar(c)
plt.plot(zeta_obs_vals,mu_range,color='black',linewidth=3)
plt.plot(zeta_68_vals,mu_range,color='red',linewidth=3)
plt.xlabel("$\zeta_{\mu}$")
plt.ylabel("$\\mu$")
plt.title("$f(\zeta_{\mu}|H(\mu))$")
plt.show()
