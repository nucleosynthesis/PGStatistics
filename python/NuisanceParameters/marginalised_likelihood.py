import numpy 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
plt.rcParams.update({'font.size': 14})

from model import *

# Likelihood --> P(n|mu,eta)
def Likelihood(n,mu,eta):
  l = lamb(mu,eta)
  return (l**n)*numpy.exp(-l)/numpy.math.factorial(int(n))

# prior on eta P(eta)
def prior_eta(eta): 
  c = 1./((2*numpy.pi)**0.5)
  return c*numpy.exp(-0.5*eta*eta)

# prior on mu P(mu)
def prior_mu(mu):
  return 1./21 if (mu < 20 and mu > -1) else 0 

# and define the product of them (the numerator in Bayes rule)
def product(eta,mu,n):
  return Likelihood(n,mu,eta)*prior_mu(mu)*prior_eta(eta)

# integrate over eta
def integral(mu,n):
  return integrate.quad(product,-6,6,args=(mu,n),epsabs=1.49e-03)[0]
    
# and then integrate over mu
def norm(n):
  return integrate.dblquad(product,-1,22,lambda x:-6,lambda x:6,args=[n],epsabs=1.49e-03)[0]

#if __name__=='main':
# plotting
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(14,6))

# plot the product of the likelihood and priors
xaxis = numpy.linspace(-2.,22,200)
yaxis = numpy.linspace(-3,3,200)
z = [ [product(eta,mu,n) for mu in xaxis] for eta in yaxis ]
X,Y = numpy.meshgrid(xaxis,yaxis) 
c = ax1.pcolor(X,Y,z)
fig.colorbar(c,ax=ax1)
ax1.set_xlabel("$\\mu$")
ax1.set_ylabel("$\eta$")
ax1.set_title("$p(n|\\mu,\eta)p(\\mu)p(\eta)$")

normalise = norm(n)
# plot the posterior
p_mu = [ integral(mu,n)/normalise for mu in xaxis ]
ax2.plot(xaxis,p_mu)
ax2.set_xlabel("$\\mu$")
ax2.set_ylabel("$p(\\mu$)")

plt.show()

# MCMC integration verison 
def proposal(mu,eta):
    return numpy.random.normal(mu,1),numpy.random.normal(eta,1)

# function to accept new proposal or reject it 
def accept_or_reject(mu_new,eta_new,mu,eta):  
    # check for 0 
    if product(eta,mu,n) > 0: 
      r = product(eta_new,mu_new,n)/product(eta,mu,n)
    else:
      r = 2.
    if r > 1: 
        return mu_new,eta_new,True 
    else:
        alpha = min([r,1])
        rnd = numpy.random.uniform(0,1)
        if rnd < alpha : 
            return mu_new,eta_new,True
        else: return mu,eta,False

init = proposal(1,0)
accepted_mc_MH = []

number_steps = 100000
for i in range(number_steps) :
    
    mu_new,eta_new = proposal(init[0],init[1])
    mu_new,eta_new,accept = accept_or_reject(mu_new,eta_new,init[0],init[1])
    init = [mu_new,eta_new]
    if accept: accepted_mc_MH.append(mu_new)
        
# plot the accepted points, overlay the numerical result
accepted_mc_MH = accepted_mc_MH[1000:-1]

res = plt.hist(accepted_mc_MH,density=True,fill=False,edgecolor='blue',bins=50)

plt.plot(xaxis, p_mu)
plt.xlabel("$\mu$")
plt.ylabel("$p(\mu|n)$")
plt.show()