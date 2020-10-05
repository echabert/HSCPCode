from Objects import *
from Loader import *

from math import sqrt

from ROOT import TH1F, TH2F, TProfile, TCanvas, TLegend, TFile, gStyle


gStyle.SetOptStat(0)


"""
*  Effect in data (pathlenght vs dEdX inclusif then as function of p)
* Comparison data/simu
* Compare with timing expectation ...
* Diminution didx en fonction de la layer [en fonction de p]
* Comparaison des estimateurs barrel/ee - thin/thick - TIB/TOB
"""


def GetIndice(val,bins):
   #for i in range(len(bins)):
   for i in range(len(bins)-1):
      if val>=bins[i] and val<bins[i+1]:
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


#scale histo to the first bin
def ScaleToFirstBin(h):
    if h.GetBinContent(1)>0:
        h.Scale(1./h.GetBinContent(1),"nosw2")


#check if the track stops in that layer
def IsLastLayer(clusters, layer):
 if layer<=11:
    for c in clusters:
       #print (c.layer,layer)
       if c.layer>layer and c.layer<=11: return False
 if layer>=12:
    for c in clusters:
       if c.layer>=layer and c.layer>=12: return False
 return True
     

def DeltaTobTib(clusters):
  meanTOB = 0
  meanTIB = 0
  elossTIB = [c.eloss for c in clusters if 1<=c.layer<=4]
  elossTOB = [c.eloss for c in clusters if 5<=c.layer<=8]
  if len(elossTIB)>0 : meanTIB = sum(elossTIB)/len(elossTIB)
  if len(elossTOB)>0 : meanTOB = sum(elossTOB)/len(elossTOB)
  if meanTOB!=0 and meanTIB!=0:
     return meanTOB-meanTIB
  else:
     return -99


##########################
# beginning of the main #
##########################

# general parameters

#study as function of the layers
layers = range(1,22)
layerLabels = ["TIB1","TIB2","TIB3","TIB4","TOB1","TOB2","TOB3","TOB4","TOB5","TOB6","TID1","TID2","TEC1","TEC2","TEC3","TEC4","TEC5","TEC6","TEC7","TEC8","TEC9"]


#bins in momenta
pfactor=1   #1000 for stop at 2 TeV
pmin = 0.5*pfactor
pmax = 1.0*pfactor #1.5
pstep = 0.1*pfactor
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
#filename = "stop2000.csv"
#filename = "bkg.csv"

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
hQTOB5_last = TH1F("hQTOB5_last","",50,0,50)
hQTOB5_notlast = TH1F("hQTOB5_notlast","",50,0,50)
hNStripVsDedx_last = TH2F("hNStripVsDedx_last","",10,0,30,10,0,20)
hNStripVsDedx_notlast = TH2F("hNStripVsDedx_nolast","",10,0,30,10,0,20)
#hPathVsDedx = TProfile("hPathVsDedx","",40,8,16,0.3,0.6)
#hPathVsDedx = TProfile("hPathVsDedx","",40,0.3,0.5,0,40)
hPathVsDedx = TH2F("hPathVsDedx","",20,0,40,50,0.3,1.0)
hDeltaTibTobvsP = TH2F("hDeltaTibTobvsP","",10,0.5,1.5,100,-5,5)

#run over data
for cand in cands:
    #element 0: Event
    #element 1: HSCPCand
    passSel = False
    if any([True for c in cand[1].clusters if c.layer == 1 and c.eloss>16]): passSel = True
    isTOB5Last =  IsLastLayer(cand[1].clusters,5)
    for c in cand[1].clusters:
        #if c.layer == 14: print(cand[1].p,GetIndice(cand[1].p,pbins)
	#if c.layer==1 and cand[1].p>0.8 and cand[1].p<1:
	if abs(cand[1].eta)<0.5 and passSel and  cand[1].p<1*pfactor:
	    #hPathVsDedx.Fill(c.eloss,c.pathlength*10)
	    if c.layer>=5:  hPathVsDedx.Fill(c.eloss,c.pathlength*10)
        #if abs(cand[1].eta)<0.5 and passSel : #and not IsLastLayer(cand[1].clusters,c.layer): 
        if abs(cand[1].eta)<0.5 and passSel : #and not IsLastLayer(cand[1].clusters,c.layer): 
	    dEdx_P_Layer_list[(c.layer,GetIndice(cand[1].p,pbins))].append(c.eloss)
	    #if c.layer==5:
	    if c.layer>=1:
	        #if isTOB5Last:
	        if c.eloss> 6:#isTOB5Last:
	            hQTOB5_last.Fill(c.eloss)
	            hNStripVsDedx_last.Fill(c.eloss,c.nstrips)
		else:
		    hQTOB5_notlast.Fill(c.eloss)
	            hNStripVsDedx_notlast.Fill(c.eloss,c.nstrips)
        h.Fill(c.layer)
	hDeltaTibTobvsP.Fill(cand[1].p,DeltaTobTib(cand[1].clusters))
h.Draw()

#print(dEdx_P_Layer_list[(14,4)])


cDelta = TCanvas()
hDeltaTibTobvsP.Draw("COLZ")

cPVsDedx = TCanvas()
hPathVsDedx.Draw("COLZ")

cNstripVsDedx = TCanvas()
cNstripVsDedx.Divide(2)
cNstripVsDedx.cd(1)
hNStripVsDedx_last.GetXaxis().SetTitle("dEdx")
hNStripVsDedx_last.GetYaxis().SetTitle("nstrips")
hNStripVsDedx_last.Draw("COLZ")
cNstripVsDedx.cd(2)
hNStripVsDedx_notlast.GetXaxis().SetTitle("dEdx")
hNStripVsDedx_notlast.GetYaxis().SetTitle("nstrips")
hNStripVsDedx_notlast.Draw("COLZ")


cLast = TCanvas()
if hQTOB5_last.Integral()>0: hQTOB5_last.Scale(1./hQTOB5_last.Integral())
hQTOB5_last.SetLineColor(2)
hQTOB5_last.Draw()
if hQTOB5_notlast.Integral()>0: hQTOB5_notlast.Scale(1./hQTOB5_notlast.Integral())
hQTOB5_notlast.Draw("same")

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
    ScaleToFirstBin(h)
    #h.Sumw2()
    dEdx_P_Layer_plots.append(h)
    h.GetYaxis().SetTitle(dEdxlabel)
    h.SetLineWidth(lineWidth)
    h.SetLineColor(color)
    h.GetXaxis().SetRangeUser(0,22)
    legP.AddEntry(h,"P = "+str(pbins[ibin])+" GeV",'l')
    if ibin==0: h.Draw()
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


