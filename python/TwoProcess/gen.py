import numpy 

signal1_counts    = [0,0,0,0.5,1,2,3,4,5,4.5,4.0,3.5,3.0,2.5,2,1.5,1,0.5,0,0,0,0,0,0,0]
signal2_counts    = [0,0,1,2,3,4,5,6,5,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
background_counts = [60*numpy.exp(-0.1*i) for i in range(len(signal1_counts))]
bkg_uncert        = [0.3,0.3,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.02,0.02,0.02,0.02]

nbins = len(signal1_counts)

def lamb(s1,s2,delta_b,b,k):
  return s1+s2+b*(1+k)**delta_b

def gendata():
    data = [numpy.random.poisson(lamb(signal1_counts[i],signal2_counts[i],0,background_counts[i],bkg_uncert[i])) for i in range(nbins)]

    return [data,0]

# we can generate new data or just stick with the data we have below
data = [[65, 45, 47, 37, 37, 40, 42, 36, 34, 36, 22, 23, 23, 18, 17, 13, 12, 12, 11, 12, 6, 6, 6, 9, 4], 0]
#data = gendata()#[[55, 45, 48, 49, 53, 36, 48, 36, 44, 31, 27, 20, 26, 14, 19, 17, 16, 14, 20, 8, 7, 11, 4, 4, 6], 0]
