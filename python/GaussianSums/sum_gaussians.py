import numpy 
import matplotlib.pyplot as plt

def gaus(x,mu,sigma): 
  return (1./(sigma*(2*numpy.math.pi)**0.5))*numpy.math.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

muX, sigmaX = 10., 1.
muY, sigmaY = -6., 0.5 

# according to our calculation for the sum of two gaussian variables
mu    = muX + muY 
sigma = ((sigmaX)**2 + (sigmaY)**2)**0.5

toys_X = numpy.random.normal(muX,sigmaX,1000)
toys_Y = numpy.random.normal(muY,sigmaY,1000)
toys_U = toys_X+toys_Y

x = numpy.arange(-10,20,0.1)
y = [gaus(xx,mu,sigma) for xx in x]

plt.hist(toys_X,100,(-10,20),density=True,color='green')
plt.hist(toys_Y,100,(-10,20),density=True,color='red')
plt.hist(toys_U,100,(-10,20),density=True,color='cyan')

#plt.plot(x,[gaus(xx,muX,sigmaX) for xx in x],color='black')
#plt.plot(x,[gaus(xx,muY,sigmaY) for xx in x],color='black')
plt.plot(x,[gaus(xx,mu,sigma) for xx in x],color='black')
print("mu=%.1f, sigma=%.1f"%(mu,sigma))

plt.xlabel("$U$")
plt.ylabel("$\\phi(U)$")
plt.show()

