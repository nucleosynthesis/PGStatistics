import numpy 
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

from model import * 

# q = -2lnL
def q(n,mu,eta): 
  la = lamb(mu,eta)
  return eta*eta + 2*la -2*n*numpy.log(la)

# dq/deta
def dq(n,mu,eta): 
  dl = lamb(mu,eta)*numpy.log(1+k)
  return 2*eta + 2*dl - 2*n*numpy.log(1+k)

# d^2q/deta^2
def d2q(n,mu,eta):
  log_k = numpy.log(1+k)
  la  = lamb(mu,eta)
  return 2+2*la*log_k*log_k

def profiled_eta(n,mu):
  # numerical minimumisation via newton method
  tol = 0.01
  init_eta =-0.5 if lamb(mu,0)-n > 0 else 0.5 
  init_eta = -10
  eps = 100
  while eps > tol: 
    qp = dq(n,mu,init_eta)
    eps = abs(qp)
    init_eta = init_eta - qp/d2q(n,mu,init_eta)
  return init_eta

# plotting
fig, (ax1,ax2) = plt.subplots(1,2)

# 1. plot q(mu,eta)
xaxis = numpy.linspace(-0.9,15,50)
yaxis = numpy.linspace(-3,3,50)
z = [ [q(n,mu,eta) for mu in xaxis] for eta in yaxis ]
X,Y = numpy.meshgrid(xaxis,yaxis) 
c = ax1.pcolor(X,Y,z)
fig.colorbar(c,ax=ax1)
ax1.set_xlabel("$\\mu$")
ax1.set_ylabel("$\eta$")
ax1.set_title("$q(\\mu,\eta)$")

# 2. plot the profiled value of eta as a function of mu
eta_mu = [ profiled_eta(n,mu) for mu in xaxis ]
ax1.plot(xaxis,eta_mu,color='red')

# 2. plot the profile likelihood
q_mu = [ q(n,mu,eta) for mu,eta in zip(xaxis,eta_mu)]
ax2.plot(xaxis,q_mu)
ax2.set_xlabel("$\\mu$")
ax2.set_ylabel("$q(\\mu)$")

plt.show()
