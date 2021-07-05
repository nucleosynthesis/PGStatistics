sigma_TH = 0.01
A   = 0.5 
eff = 0.9
l   = 100.
k   = 0.1 
B   = 1.0 
n   = 2

# Poisson mean
def lamb(mu,eta):
  return mu*eff*A*l*((1+k)**eta)*sigma_TH + B
