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
  #data = args[0]
  return q(args,mu1,mu2,delta_b)

def q_constrained(x, args):
  delta_b = x[0]
  data, mu1, mu2 = args[0], args[1], args[2]
  return q(data,mu1,mu2,delta_b)

def profiled_nuis(data, mu1, mu2):
  init_params = [-3.0]
  bounds =  [(-5,5)]
  res = minimize(q_constrained,init_params,args=[data,mu1,mu2],bounds=bounds)
  return res.x[0]

def global_min(data): 
  init_params = [0.1,0.1,-3.]
  bounds = [(-10,10),(-10,10),(-5,5)]
  mle = minimize(q_unconstrained,init_params,args=data,bounds=bounds)
  return mle.fun,mle.x[0],mle.x[1]

def delta_q(data,mu1,mu2,q_min): 
  q_value        = q(data,mu1,mu2,profiled_nuis(data,mu1,mu2))
  return q_value-q_min 


q_min,mu1_min,mu2_min   = global_min(data)

# plotting
mu1_axis = numpy.linspace(-1,5,50)
mu2_axis = numpy.linspace(-2,5,50)
z = [ [delta_q(data,m1,m2,q_min) for m1 in mu1_axis] for m2 in mu2_axis] 

X,Y = numpy.meshgrid(mu1_axis,mu2_axis) 
c = plt.pcolor(X,Y,z)
plt.colorbar(c)
conts = plt.contour(X,Y,z,[2.3,5.99],colors="white")
plt.clabel(conts, fontsize=9, inline=1)
plt.xlabel("$\\mu_{1}$")
plt.ylabel("$\\mu_{2}$")
plt.title("$\Delta q(\mu_{1},\mu_{2})$")
plt.plot([mu1_min],[mu2_min],color="white",marker="P")



plt.show()
