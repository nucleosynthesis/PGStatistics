import numpy 
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

# true values
mu_0     = 4.0
sigma_0  = 2.0

# empty list of data for now
data = [] 

def q(params_list = []):
   mu    = params_list[0]
   sigma = params_list[1]
   N = len(data)
   sum_part = sum( [ ((x-mu)/sigma)**2 for x in data ] )
   return 2*N*numpy.log(sigma*((2*numpy.pi)**0.5)) + sum_part

# inital values for mu, sigma
init_params = [mu_0,sigma_0]
bounds  = [(-10,10),(0.01,10)]
# generate random data
steps = [5,10,20, 50, 100, 200, 500, 1000, 5000, 10000, 20000]
hat_mu = []; hat_sigma = []
for N in steps: 
  data = numpy.random.normal(mu_0,sigma_0,size=N)
  mle = minimize(q,init_params,bounds=bounds)
  hat_mu.append(mle.x[0])
  hat_sigma.append(mle.x[1])

plt.plot(steps,hat_mu,color='blue',marker="o",label="$\hat{\mu}$")
plt.plot(steps,hat_sigma,color='red',marker="o",label="$\hat{\sigma}$")

plt.plot(steps,[mu_0 for s in steps],color='blue',linestyle="--")
plt.plot(steps,[sigma_0 for s in steps],color='red',linestyle="--")

plt.xlabel("N")
plt.ylabel("$\hat{\mu}$ or $\hat{\sigma}$")

plt.xscale('log')
plt.legend()
plt.show()
