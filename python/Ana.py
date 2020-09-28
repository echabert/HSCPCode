from Objects import *
from Loader import *

from math import sqrt

from ROOT import TH1F, TH2F, TCanvas, TLegend, TFile, gStyle


gStyle.SetOptStat(0)


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
   return len(bins)-1


def Statistics(mylist):
   if len(mylist) == 0: return {}
   mylist.sort()
   N = len(mylist)
   mean = sum(mylist)/N
   stddev = pow(sum([e**2 for e in mylist])-mean**2,0.5)
   median = mylist[int(N/2)]
   return {"mean":mean,"meanError":stddev/sqrt(N),"median":median,"stddev":stddev}



##########################
# beginning of the main #
##########################

# general parameters

#study as function of the layers
layers = range(1,22)
layerLabels = ["TIB1","TIB2","TIB3","TIB4","TOB1","TOB2","TOB3","TOB4","TOB5","TOB6","TID1","TID2","TEC1","TEC2","TEC3","TEC4","TEC5","TEC6","TEC7","TEC8","TEC9"]

#bins in momenta
pmin = 0.5
pmax = 1.5
pstep = 0.1
npsteps = int((pmax-pmin)/pstep)
pbins = [pmin+i*pstep for i in range(npsteps+1)]

#study estimators
estimators = {}
#tuple with (IsH2,IsFiltered,lowFrac,highFrac)
UseFilter=False
estimators["Mean(0,1)"]=(UseFilter,False,0,1)
estimators["Mean(0.15,1)"]=(UseFilter,False,0.15,1)
estimators["Mean(0,0.85)"]=(UseFilter,False,0,0.85)

# plotting tuning
dEdxlabel = "<dEdx> [MeV.cm^{2}/g]"
lineWidth = 2

#file to be analyzed
#filename = "../bin/deuteron.csv"
filename = "deuteron.csv"

# most important line
# running over the data !!!
cands = LoadCSV(filename)
print("We ran over",len(cands),"candidates")



###############################################
# Study of the timing as function of the layer
###############################################

#create dictionnaries
dEdx_P_Layer_list = {}
dEdx_P_Layer_res = {}
for layer in layers:
    for ibin in range(len(pbins)):
       dEdx_P_Layer_list[(layer,ibin)] = []
       dEdx_P_Layer_res[(layer,ibin)] = None

#plots label
h = TH1F("hLabel","",21,1,22)

#run over data
for cand in cands:
    #element 0: Event
    #element 1: HSCPCand
    for c in cand[1].clusters:
        #if c.layer == 14: print(cand[1].p,GetIndice(cand[1].p,pbins)
        if abs(cand[1].eta)<0.1: dEdx_P_Layer_list[(c.layer,GetIndice(cand[1].p,pbins))].append(c.eloss)
        h.Fill(c.layer)
h.Draw()
print(dEdx_P_Layer_list[(14,4)])


#statistics per category
for layer in layers:
    for ibin in range(len(pbins)):
        dEdx_P_Layer_res[(layer,ibin)] = Statistics(dEdx_P_Layer_list[(layer,ibin)])


print(dEdx_P_Layer_res)

#plots and layout
dEdx_P_Layer_plots = []
cTimingVsLayer = TCanvas()
legP = TLegend()
color = 1
for ibin in range(len(pbins)):
    h = TH1F("dEdxPerLayer_"+str(ibin),"",21,1,22)
    for layer in layers:
        #print(dEdx_P_Layer_res[(layer,ibin)])
        #h.SetBinContent(layer,dEdx_P_Layer_res[(layer,ibin)].get("median",0))
        h.SetBinContent(layer,dEdx_P_Layer_res[(layer,ibin)].get("mean",0))
        #h.SetBinError(layer,dEdx_P_Layer_res[(layer,ibin)].get("meanError",0))
        h.GetXaxis().SetBinLabel(layer,layerLabels[layer-1])
    dEdx_P_Layer_plots.append(h)
    h.GetYaxis().SetTitle(dEdxlabel)
    h.SetLineWidth(lineWidth)
    h.SetLineColor(color)
    legP.AddEntry(h,"P = "+str(pbins[ibin])+" GeV",'l')
    if ibin==1: h.Draw()
    else: h.Draw("same")
    color+=1

legP.Draw()

###############################################
# Study of the estimators
###############################################

#create dictionnaries
Estim_P_list = {}
Estim_P_res = {}
for estim in estimators:
    for ibin in range(len(pbins)):
       Estim_P_list[(estim,ibin)] = []
       Estim_P_res[(estim,ibin)] = None
       
#run over data
for cand in cands:
    #for c in cand[1].clusters:
    for estim in estimators:
        print(estimators[estim][2],estimators[estim][3])
        e = Estimator(cand[1].GetClusterEloss(estimators[estim][1]))
        value = 0
        if estimators[estim][0]: value = e.GetHarmonic(estimators[estim][2],estimators[estim][3])
        else: value = e.GetMean(estimators[estim][2],estimators[estim][3])
        print(value)
        Estim_P_list[(estim,GetIndice(cand[1].p,pbins))].append(value)

#statistics per estimator
for estim in estimators:
    for ibin in range(len(pbins)):
        print(Estim_P_list[(estim,ibin)])
        Estim_P_res[(estim,ibin)] = Statistics([c[0] for c in Estim_P_list[(estim,ibin)]])

#plots and layout
estimMean_P_plots = []
estimReso_P_plots = []
cMeanEstim = TCanvas()
cResoEstim = TCanvas()
legEstim = TLegend()
color = 1
for estim in estimators:
    hMean = TH1F("hMeanEstim_"+estim,"",npsteps+1,pmin,pmax)
    hReso = TH1F("hResoEstim_"+estim,"",npsteps+1,pmin,pmax)
    for ibin in range(len(pbins)):
        print(Estim_P_res[(estim,ibin)]["mean"])
        hMean.SetBinContent(ibin,Estim_P_res[(estim,ibin)]["mean"])
        hMean.SetBinError(ibin,Estim_P_res[(estim,ibin)]["meanError"])
        hMean.SetBinContent(ibin,Estim_P_res[(estim,ibin)]["stddev"])
    estimMean_P_plots.append(hMean)
    estimReso_P_plots.append(hReso)
    hMean.GetXaxis().SetTitle("p [GeV]")
    hReso.GetXaxis().SetTitle("p [GeV]")
    hMean.GetYaxis().SetTitle(dEdxlabel)
    hReso.GetYaxis().SetTitle("std-dev(<dEdx>)/<dEdx>")
    hMean.SetLineWidth(lineWidth)
    hReso.SetLineWidth(lineWidth)
    hMean.SetLineColor(color)
    hReso.SetLineColor(color)
    legEstim.AddEntry(h,estim,'l')
    cMeanEstim.cd()
    if color==1: hMean.Draw()
    else: hMean.Draw("same")
    cResoEstim.cd()
    if color==1: hReso.Draw()
    else: hReso.Draw("same")
    color+=1
    estimMean_P_plots.append(hMean)
    estimReso_P_plots.append(hReso)

cMeanEstim.cd()
cResoEstim.cd()
legEstim.Draw("same")


