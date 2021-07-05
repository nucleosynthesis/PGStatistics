import sys
import numpy 
import matplotlib.pyplot as plt
from scipy.optimize import minimize
plt.rcParams.update({'font.size': 14})

from gen import * 

nbins = len(data[0])
def lamb(s1,s2,delta_b,b,k):
  return s1+s2+b*(1+k)**delta_b

def log_poisson(data,i,mu1,mu2,delta_b):
  s1 = signal1_counts[i]
  s2 = signal2_counts[i]
  b  = background_counts[i]
  k  = bkg_uncert[i]
  l = lamb(mu1*s1,mu2*s2,delta_b,b,k)
  n = data[i]
  if n==0: return -l
  else: return n*numpy.log(l) - l 

def q(data,mu1,mu2,delta_b):
  constr = (data[1]-delta_b)**2
  return constr -2*sum([log_poisson(data[0]\
    ,i,mu1,mu2,delta_b) for i in range(nbins)])

def q_unconstrained(x, args):
  mu1, mu2, delta_b = x[0], x[1], x[2]
  return q(args,mu1,mu2,delta_b)

# this function constraints mu1 (scan mu1)
def q_constrained_mu1(x, args):
  delta_b = x[1]
  mu2   = x[0]
  data, mu1 = args[0], args[1]
  return q(data,mu1,mu2,delta_b)

def q_constrained_mu2(x, args):
  delta_b = x[1]
  mu1   = x[0]
  data, mu2 = args[0], args[1]
  return q(data,mu1,mu2,delta_b)

def q_constrained(x, args):
  delta_b = x[0]
  data, mu1, mu2 = args[0], args[1], args[2]
  return q(data,mu1,mu2,delta_b)

def global_min(data): 
  init_params = [0.1,0.1,-3.]
  bounds = [(-10,10),(-10,10),(-5,5)]
  mle = minimize(q_unconstrained,init_params,args=data,bounds=bounds)
  return mle.fun,mle.x[0],mle.x[1]

def profiled(data, mu1):
  init_params = [1.0,-3.0]
  bounds =  [(-4,10),(-5,5)]
  res = minimize(q_constrained_mu1,init_params,args=[data,mu1],bounds=bounds)
  return res.x

def delta_qmu1(data,mu1,q_min): 
  profiled_nuisances = profiled(data,mu1)
  q_value        = q(data,mu1,profiled_nuisances[0],profiled_nuisances[1])
  return q_value-q_min 

def profiled_nuis(data, mu1, mu2):
  init_params = [-3.0]
  bounds =  [(-5,5)]
  res = minimize(q_constrained,init_params,args=[data,mu1,mu2],bounds=bounds)
  return res.x[0]

def delta_qmu1_fix(data,mu1,mu2,q_min): 
  profiled_nuisances = profiled_nuis(data,mu1,mu2)
  q_value        = q(data,mu1,mu2,profiled_nuisances)
  return q_value-q_min 


q_min,mu1_min,mu2_min   = global_min(data)

def return_crossing(x1,y1,x2,y2,K): 
  m = (y2-y1)/(x2-x1)
  c = y1-m*x1
  return (K-c)/m,K

def findIntervals(x,y,conts=[1,4]): 
  xx0,yy0 = x[0],y[0]
  crossing_x = []
  for xx,yy in zip(x[1:],y[1:]): 
    for K in conts: 
      if (yy < K and yy0 > K) or (yy > K and yy0 < K): 
        crossing_x.append(return_crossing(xx0,yy0,xx,yy,K))
    xx0=xx
    yy0=yy
  return crossing_x
  
# plotting
mu1_axis = numpy.linspace(-1,5,50)
z = [ delta_qmu1(data,m1,q_min) for m1 in mu1_axis] 
zf = [ delta_qmu1_fix(data,m1,mu2_min,q_min) for m1 in mu1_axis]
plt.plot(mu1_axis,z)
plt.plot(mu1_axis,zf,linestyle="--")
plt.xlabel("$\\mu_{1}$")
plt.ylabel("$\Delta q(\mu_{1})$")
intervals = findIntervals(mu1_axis,z,[1,4])
[ plt.plot([x[0],x[0]],[0,x[1]],color='red') for x in intervals ]
plt.ylim((0,14))
[ plt.plot([mu1_axis[0],mu1_axis[-1]],[k,k],color='black',linestyle='--') for k in [1,4] ]
#plt.plot([mu1_min],[mu2_min],color="white",marker="P")
plt.show()
