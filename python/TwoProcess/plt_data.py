import ROOT
import numpy 
from scipy.optimize import minimize
from gen import * 

#data = gendata()
#data = [[55, 45, 48, 49, 53, 36, 48, 36, 44, 31, 27, 20, 26, 14, 19, 17, 16, 14, 20, 8, 7, 11, 4, 4, 6], 0]
#signal1_counts    = [0,0,0,0,0.5,1,2,3,4,5]
#signal2_counts    = [6,5,4,3,2,1.5,1.0,.5,0,0]
#background_counts = [50,50,50,50,50,50,50,50,50,50,50]
#:wq
#bkg_uncert        = [0.1,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01]

nbins = len(data[0])

def lamb(s1,s2,delta_b,b,k):
  return s1+s2+b*(1+k)**delta_b

def log_poisson(data,i,mu1,mu2,delta_b):
  s1 = signal1_counts[i]
  s2 = signal2_counts[i]
  b  = background_counts[i]
  k  = bkg_uncert[i]
  l = lamb(mu1*s1,mu2*s2,delta_b,b,k)
  n = data[i]
  if n==0: return -l
  else: return n*numpy.log(l) - l 

def q(data,mu1,mu2,delta_b):
  constr = (data[1]-delta_b)**2
  return constr -2*sum([log_poisson(data[0]\
    ,i,mu1,mu2,delta_b) for i in range(nbins)])

def q_unconstrained(x, args):
  mu1, mu2, delta_b = x[0], x[1], x[2]
  #data = args[0]
  return q(args,mu1,mu2,delta_b)

def q_constrained(x, args):
  delta_b = x[0]
  data, mu1, mu2 = args[0], args[1], args[2]
  return q(data,mu1,mu2,delta_b)

def profiled_nuis(data, mu1, mu2):
  init_params = [-3.0]
  bounds =  [(-5,5)]
  res = minimize(q_constrained,init_params,args=[data,mu1,mu2],bounds=bounds)
  return res.x[0]

def global_min(data): 
  init_params = [0.1,0.1,-3.]
  bounds = [(-10,10),(-10,10),(-5,5)]
  mle = minimize(q_unconstrained,init_params,args=data,bounds=bounds)
  return mle.fun,mle.x[0],mle.x[1],mle.x[2]

def delta_q(data,mu1,mu2,q_min): 
  q_value        = q(data,mu1,mu2,profiled_nuis(data,mu1,mu2))
  return q_value-q_min 


q_min,mu1_min,mu2_min,delta_b_min   = global_min(data)

hdata = ROOT.TH1F("data",";Bin;Events",nbins,0,nbins)
hbackground = ROOT.TH1F("data","data",nbins,0,nbins)
hsignal1 = ROOT.TH1F("signal_1","signal_1",nbins,0,nbins)
hsignal2 = ROOT.TH1F("signal_2","signal_2",nbins,0,nbins)

for i in range(nbins):
  hdata.SetBinContent(i+1,data[0][i])
  hdata.SetBinError(i+1,(data[0][i])**0.5)
  hsignal1.SetBinContent(i+1,mu1_min*signal1_counts[i])
  hsignal2.SetBinContent(i+1,mu2_min*signal2_counts[i])
  hbackground.SetBinContent(i+1,background_counts[i]*(1+bkg_uncert[i])**delta_b_min)

hbackground.SetFillColor(ROOT.kGreen+2)
hsignal1.SetFillColor(ROOT.kOrange)
hsignal2.SetFillColor(ROOT.kAzure+10)

hbackground.SetLineColor(1)
hsignal1.SetLineColor(1)
hsignal2.SetLineColor(1)
hdata.SetLineColor(1)
hdata.SetLineWidth(2)
hdata.SetMarkerColor(1)
hdata.SetMarkerStyle(20)

stk = ROOT.THStack()
stk.Add(hbackground)
stk.Add(hsignal1)
stk.Add(hsignal2)
print(mu1_min,mu2_min)
print(data)
ROOT.gStyle.SetOptStat(0)

c=ROOT.TCanvas("c","c",800,600)
c.SetBottomMargin(0.15)
c.SetLeftMargin(0.15)
hdata.GetXaxis().SetTitleSize(0.04)
hdata.GetYaxis().SetTitleSize(0.04)
hdata.Draw("P")
hdata.SetMinimum(0)
hdata.SetMaximum(80)
stk.Draw("histsame")
hdata.Draw("psame")
leg = ROOT.TLegend(0.2,0.74,0.84,0.84)
leg.SetNColumns(2)
leg.SetBorderSize(0)
leg.AddEntry(hdata,"Observed data","PEL")
leg.AddEntry(hbackground,"Fitted background","F")
leg.AddEntry(hsignal1,"#hat{#mu}_{1}#times signal_{1}","F")
leg.AddEntry(hsignal2,"#hat{#mu}_{2}#times signal_{2}","F")
leg.Draw()
c.SetTicky()
c.SetTickx()
c.RedrawAxis()
c.SaveAs("test.pdf")

