import numpy 
from model import *
from calc_ftmu import *
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

p_mu_list = []
mu_test_list = numpy.arange(8,15,0.5) 

# for each value of mu to test, calculate p_mu
for mu_test in mu_test_list : 
  t_obs = tmu(n,0,mu_test)
  tmu_toys,pval = histo_tmu(mu_test,t_obs)
  p_mu_list.append(pval)

# from the graph, we can read of upper limits (eg for 95\% CL)
plt.plot(mu_test_list,p_mu_list,color='black',marker="o")
plt.plot([8,15],[0.05,0.05],color='red')

plt.xlabel("$\mu$")
plt.ylabel("$p_{\mu}$")
plt.show()
