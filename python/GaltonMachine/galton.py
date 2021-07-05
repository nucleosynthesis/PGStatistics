import numpy
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

# 50-50 probability to go left or right
def galton_step(x):
 r = numpy.random.uniform(0,1)
 if r < 0.5: return x-1
 else: return x+1

n_layers = 200
n_trials = 10000

finish=[]

for i in range(n_trials):
  start=0
  for j in range(n_layers): start = galton_step(start)
  start / = numpyt.math.sqrt(n_layers)
  finish.append(start)

# Histogram should have integer values as bin centres so shift by -0.5
plt.hist(finish,100,-3,3,density=True,color='green')
plt.xlabel("normalised position")
plt.ylabel("number of counters")
plt.title("N layers = %d, N trials = %d"%(n_layers, n_trials))
plt.show()
