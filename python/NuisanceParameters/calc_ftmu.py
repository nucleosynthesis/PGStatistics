import numpy 
from model import * 
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

def q(n,eta_p,mu,eta): 
  la = lamb(mu,eta)
  return (eta-eta_p)*(eta-eta_p) + 2*la -2*n*numpy.log(la)

# dq/deta
def dq(n,eta_p,mu,eta): 
  dl = lamb(mu,eta)*numpy.log(1+k)
  return 2*(eta-eta_p) + 2*dl - 2*n*numpy.log(1+k)

# d^2q/deta^2
def d2q(n,mu,eta):
  log_k = numpy.log(1+k)
  la  = lamb(mu,eta)
  return 2+2*la*log_k*log_k

def profiled_eta(n_p,eta_p,mu):
  # numerical minimumisation via newton method
  tol = 0.01
  init_eta = -3
  eps = 100
  while eps > tol: 
    qp = dq(n,eta_p,mu,init_eta)
    eps = abs(qp)
    init_eta = init_eta - qp/d2q(n,mu,init_eta)
  return init_eta

def global_min(n_p,eta_p): 
  vals = []
  for mu_v in numpy.arange(-1,20,0.1):
   #if lamb(mu_v,eta_p) < 0: continue  # avoid hairy situations
   eta_profiled = profiled_eta(n_p,eta_p,mu_v)
   vals.append([q(n_p,eta_p,mu_v,eta_profiled),mu_v])
  vals.sort()
  return vals[0][0],vals[0][1]
   
# calculate test statistic 
def tmu(n_p,eta_p,mu):
  q_value        = q(n_p,eta_p,mu,profiled_eta(n_p,eta_p,mu))
  q_min,mu_min   = global_min(n_p,eta_p)
  #print("toy",mu,n_p,eta_p,q_value, q_min,mu_min, q_value-q_min)   
  if mu_min < 0     : return q_value-q(n_p,eta_p,0,profiled_eta(n_p,eta_p,0))
  elif mu_min <= mu : return q_value-q_min
  else              : return 0 

def histo_tmu(mu,t_obs):
  # find the best (profiled) nuisance parameter values for the data (n,0)
  eta_profiled = profiled_eta(n,0,mu) 
  ntoys = 5000
  toy_n   = numpy.random.poisson(lamb(mu,eta_profiled),size=ntoys)
  toy_eta = numpy.random.normal(eta_profiled,1,size=ntoys)
  tmu_dist = [tmu(n_p,eta_p,mu) for n_p,eta_p in zip(toy_n,toy_eta)]
  return tmu_dist, float(len([tm for tm in tmu_dist if tm > t_obs]))/len(tmu_dist)

# plot only if running this code 
if __name__=='__main__':
    mu_test=10
    t_obs = tmu(n,0,mu_test)
    tmu_toys,pval = histo_tmu(mu_test,t_obs)
    plt.hist(tmu_toys,density=True,bins=24,color='black',histtype='step')
    plt.plot([t_obs,t_obs],[0,0.5],color='red')
    plt.yscale("log")
    plt.xlabel("$t_{%.1f}$"%mu_test)
    plt.ylabel("$f(t_{%.1f}|H(%.1f)$"%(mu_test,mu_test))
    print("p_mu = ",pval)
    plt.show()
