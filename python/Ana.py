from Objects import *
from Loader import *

from math import sqrt


dEdxlabel = "<dEdx> [MeV.cm^{2}/g]"
lineWidth = 2
filename = "../bin/deuteron.csv"

cands = LoadCSV(filename)

print(len(cands))


"""
*  Effect in data (pathlenght vs dEdX inclusif then as function of p)
* Comparison data/simu
* Compare with timing expectation ...
* Diminution didx en fonction de la layer [en fonction de p]
* Comparaison des estimateurs barrel/ee - thin/thick - TIB/TOB
"""


def GetIndice(val,bins):
   for i in range(len(bins)):
      if val>bins[i] and val<bins[i+1]:
         return i
   return len(bins)


def Statistics(mylist):
   if len(mylist) == 0: return {}
   mylist.sort()
   size = len(mylist)
   mean = sum(mylist)/N
   stddev = pow(sum([e**2 for e in mylist])-mean**2,0.5)
   median = mylist[N/2]
   return {"mean":mean,"meanError":stddev/sqrt(mysize),"median":median,"stddev":stddev}


layers = range(1,22)
layerLabels = ["TIB1","TIB2","TIB3","TIB4","TOB1","TOB2","TOB3","TOB4","TOB5","TOB6","TID1","TID1","TEC1","TEC2","TEC3","TEC4","TEC5","TEC6","TEC7","TEC8","TEC9"]
pmin = 0.5
pmax = 1.5
pstep = 0.1
npsteps = int((pmax-pmin)/pstep)
pbins = [pmin+i*pstep for i in range(npsteps+1)]

#create dictionnaries
dEdx_P_Layer_list = {}
dEdx_P_Layer_res = {}
for layer in layers:
    for ibin in range(len(pbins)):
       dEdx_P_Layer_list[(layer,ibin)] = []
       dEdx_P_Layer_res[(layer,ibin)] = None
       
#run over data
for cand in cands:
    for c in cand.clusters:
        dEdx_P_Layer_list[(c.layer,GetIndice(c.p,pbins))].append(c.eloss)

#statistics per category
for layer in layers:
    for ibin in range(len(pbins)):
	dEdx_P_Layer_res[(layer,ibin)] = Statistics(dEdx_P_Layer_list[(c.layer,GetIndice(c.p,pbins))])

#plots and layout
dEdx_P_Layer_plots = []
cTimingVsLayer = TCanvas()
legP = TLegend()
color = 1
for ibin in range(len(pbins)):
    h = TH1F("dEdxPerLayer_"+str(ibin),"",21,1,22)
    for layer in layers:
        h.SetBinContent(ibin,dEdx_P_Layer_res[(layer,ibin)]["mean"])
        h.SetBinError(ibin,dEdx_P_Layer_res[(layer,ibin)]["meanError"])
    	h.GetXaxis().SetBinLabel(layerLabels(layer))
    dEdx_P_Layer_plots.append(h)
    h.GetYaxis().SetTitle(dEdxlabel)
    h.SetLineWidth(lineWidth)
    h.SetLineColor(color)
    leg.AddEntry(h,"P = "+str(pbins[ibin])+" GeV",'l')
    if ibin==1: h.Draw()
    else: h.Draw("same")

leg.Draw()


