import numpy 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
plt.rcParams.update({'font.size': 14})

from model import *

# Likelihood P(n|mu,eta)
def Likelihood(n,mu,eta):
  l = lamb(mu,eta)
  return (l**n)*numpy.exp(-l)/numpy.math.factorial(int(n))

# prior on eta P(eta)
def prior_eta(eta): 
  c = 1./(numpy.pi**0.5)
  return c*numpy.exp(-0.5*eta*eta)

# prior on mu P(mu)
def prior_mu(mu):
  return 1./15 if (mu < 14 and mu > -1) else 0 

# and define the product of them
def product(eta,mu,n):
  return Likelihood(n,mu,eta)*prior_mu(mu)*prior_eta(eta)

def integral(n,mu):
  return integrate.quad(product,-4,4,args=(mu,n))[0]

def norm(n):
  return integrate.quad(integral,-0.9,14,args=(n))[0]

#if __name__=='main':
# plotting
fig, (ax1,ax2) = plt.subplots(1,2)

# plot the product of the likelihood and priors
xaxis = numpy.linspace(-0.9,10,50)
yaxis = numpy.linspace(-3,3,50)
z = [ [product(eta,mu,n) for mu in xaxis] for eta in yaxis ]
X,Y = numpy.meshgrid(xaxis,yaxis) 
c = ax1.pcolor(X,Y,z)
fig.colorbar(c,ax=ax1)
ax1.set_xlabel("$\\mu$")
ax1.set_ylabel("$\eta$")
ax1.set_title("$L(n|\\mu,\eta)P(\\mu)P(\eta)$")

# plot the posterior
p_mu = [ integral(n,mu)/norm(n) for mu in xaxis ]
ax2.plot(xaxis,p_mu)
ax2.set_xlabel("$\\mu$")
ax2.set_ylabel("$P(\\mu$)")


plt.show()
