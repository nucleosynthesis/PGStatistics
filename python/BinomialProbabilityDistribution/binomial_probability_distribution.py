import numpy as np
import matplotlib.pyplot as plt
# parameters of the binomial probability distribution
N, p=5000,0.70
# our observed number of sparks
K=3420

nMC=10000
random_k = np.random.binomial(N,p,nMC)
plt.hist(random_k,N,[75,N],density=True,color='green')

# count the number of times we see a value of k 
# less than or equal to K
pval = float(len([x for x in random_k if x<=K]))/nMC
print(pval,N,p)
plt.show()
