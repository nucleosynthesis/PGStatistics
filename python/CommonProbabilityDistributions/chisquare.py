import numpy
from scipy.stats import chi2
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

x = numpy.arange(0.1,20,0.05)
ndof = range(1,12,1)

# plotting
fig, (ax1,ax2) = plt.subplots(1,2)

[ ax1.plot(x,chi2.pdf(x,i),label="$\chi^{2}(%d)$"%i) for i in ndof ]
[ ax2.plot(x,chi2.cdf(x,i),label="$\chi^{2}(%d)$"%i) for i in ndof ]

ax1.set_xlabel("$X$")
ax1.set_ylabel("$f(X;n)$")
ax2.set_xlabel("$X$")
ax2.set_ylabel("$F(X;n)$")

ax1.legend()
plt.show()
