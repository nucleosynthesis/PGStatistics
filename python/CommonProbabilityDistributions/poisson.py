import numpy 
from scipy.stats import poisson
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

x = range(0,20) 
# 5 different values of \lambda
lambdas = numpy.arange(0.1,10,2.0)
fig, (ax1,ax2) = plt.subplots(1,2)

[ ax1.plot(x,poisson.pmf(x,i),label="$P(%g)$"%i,marker="o") for i in lambdas ]
[ ax2.plot(x,poisson.cdf(x,i),label="$P(%g)$"%i,marker="o") for i in lambdas ]

ax1.set_xlabel("$k$")
ax1.set_ylabel("$f(k;\lambda)$")
ax2.set_xlabel("$k$")
ax2.set_ylabel("$F(k;\lambda)$")

ax1.legend()
plt.show()

