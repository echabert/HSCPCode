from Objects import *
from Loader import *

from math import sqrt

from ROOT import TF1, TH1F, TH2F, TProfile, TCanvas, TLegend, TFile, gStyle


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
   print("In stat",len(mylist))
   if len(mylist) == 0: return {}
   mylist.sort()
   N = len(mylist)
   mean = sum(mylist)/N
   print("mean = ",mean)
   stddev = pow(sum([e**2 for e in mylist])/N-mean**2,0.5)
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



def MassReco(K,C,dedx,p):
   if K==0: return -1
   val = (dedx-C)*p*p/K
   if val<0: return -2
   return sqrt(val)


##########################
# beginning of the main #
##########################

# general parameters

#study as function of the layers
layers = range(1,22)
layerLabels = ["TIB1","TIB2","TIB3","TIB4","TOB1","TOB2","TOB3","TOB4","TOB5","TOB6","TID1","TID2","TEC1","TEC2","TEC3","TEC4","TEC5","TEC6","TEC7","TEC8","TEC9"]

lowmass = False

#bins in momenta
pfactor=1 
if not lowmass:
    pfactor = 1000 #for stop at 2 TeV
pmin = 0.5*pfactor
pmax = 1.4*pfactor #1.5
pstep = 0.05*pfactor
npsteps = int((pmax-pmin)/pstep)
pbins = [pmin+i*pstep for i in range(npsteps+1)]

#study estimators
estimators = {}
#tuple with (IsH2,IsFiltered,lowFrac,highFrac)
UseFilter=False #False
UseHarmonic=False
estimators["Mean(0,1)"]=(UseFilter,UseHarmonic,0,1)
estimators["Mean(0,1)-f"]=(True,UseHarmonic,0,1)
#estimators["Mean(0.15,1)"]=(UseFilter,UseHarmonic,0.15,1)
#estimators["Mean(0,0.9)"]=(UseFilter,UseHarmonic,0,0.9)
#estimators["Mean(0,0.85)"]=(UseFilter,UseHarmonic,0,0.85)
#estimators["Mean(0,0.75)"]=(UseFilter,UseHarmonic,0,0.75)
#estimators["Mean(0,0.65)"]=(UseFilter,UseHarmonic,0,0.65)
estimators["Mean(0.15,0.85)"]=(UseFilter,UseHarmonic,0.15,0.85)
estimators["Mean(0.15,0.85)-f"]=(True,UseHarmonic,0.15,0.85)
estimators["Mean(0.25,0.75)"]=(UseFilter,UseHarmonic,0.25,0.75)
estimators["Mean(0.25,0.75)-f"]=(True,UseHarmonic,0.25,0.75)
#estimators["Mean(0,0.6)"]=(UseFilter,UseHarmonic,0,0.6)
#estimators["Mean(0.15,0.7)"]=(UseFilter,UseHarmonic,0.15,0.7)
#estimators["Mean(0.4,0.7)"]=(UseFilter,UseHarmonic,0.4,0.7)

# plotting tuning
dEdxlabel = "<dEdx> [MeV.cm^{2}/g]"
lineWidth = 2

#other constants
deuteronMass = 1.8756 #GeV


#fit results
fitFile = open("fitRes_deuteron.txt")

#file to be analyzed
#filename = "../bin/deuteron.csv"
#filename = "deuteron_light.csv"
#filename = "stop2000.csv"
filename = ""
if lowmass:
    filename = "deuteron.csv"
else:
    filename = "bkg_light2.csv"

# most important line
# running over the data !!!
cands = LoadCSV(filename)
print("We ran over",len(cands),"candidates")


"""
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

"""

###############################################
# Study of the estimators
###############################################

#create dictionnaries
#Estim_P_list = {}
#Estim_P_res = {}
PlotMass = {}
Estim_fit = {}
fitContent = fitFile.read()
fitFile.close()
for estim in estimators:
    #for ibin in range(len(pbins)):
       #Estim_P_list[(estim,ibin)] = []
      # Estim_P_res[(estim,ibin)] = None
    if lowmass:
        PlotMass[estim] = TH1F("Mass_"+estim,"",100,0.5,3)
    else:
        PlotMass[estim] = TH1F("Mass_"+estim,"",100,0,1000)
    for line in fitContent.split("\n"):
        print(line)
	if len(line.split()) == 0: continue
	if line.split()[0] == estim:
            Estim_fit[estim] = [float(line.split()[1]),float(line.split()[2])]

print(Estim_fit)

#fit function
BBfunc = TF1("BBFunc","[0]+[1]/(x*x)",0,2)

pDistrib = TH1F("pDistrib","",100,pmin,pmax)

#p distribution
for cand in cands:
    pDistrib.Fill(cand[1].p)

#run over data
for cand in cands:
    #for c in cand[1].clusters:
    for estim in estimators:
        #estimation 
 	e = Estimator(cand[1].GetClusterEloss(estimators[estim][0])) 
	dedx = []
	if estimators[estim][1]: dedx = e.GetHarmonic(estimators[estim][2],estimators[estim][3])
        else:  dedx = e.GetMean(estimators[estim][2],estimators[estim][3])
	#reconstruct the mass
	if len(dedx)>0:
	    mass = MassReco(Estim_fit[estim][0],Estim_fit[estim][1],dedx[0],cand[1].p) 
  	    if mass>1000: mass = 999.
	    weight = 1
	    #if pDistrib.GetBinContent(pDistrib.FindBin(cand[1].p)): weight = pDistrib.GetBinContent(pDistrib.FindBin(cand[1].p))
	    PlotMass[estim].Fill(mass,weight)

#statistics per estimator
#for estim in estimators:

#plots and layout
cMassEstim = TCanvas()
legEstim = TLegend()
color = 1
for estim in estimators:
    #PlotMass[estim].Draw()
    """
    hMean = TH1F("hMeanEstim_"+estim,"",npsteps,pmin,pmax)
    hReso = TH1F("hResoEstim_"+estim,"",npsteps,pmin,pmax)
    for ibin in range(len(pbins)):
        #print(Estim_P_res[(estim,ibin)]["mean"])
	print((estim,ibin),Estim_P_res[(estim,ibin)])
	hMean.SetBinContent(ibin,Estim_P_res[(estim,ibin)]["mean"])
        hMean.SetBinError(ibin,Estim_P_res[(estim,ibin)]["meanError"])
        hReso.SetBinContent(ibin,Estim_P_res[(estim,ibin)]["stddev"]/Estim_P_res[(estim,ibin)]["mean"])
    estimMean_P_plots.append(hMean)
    estimReso_P_plots.append(hReso)
    """
    #report
    tail300 = PlotMass[estim].Integral(PlotMass[estim].FindBin(300),PlotMass[estim].GetNbinsX()+1)
    tail500 = PlotMass[estim].Integral(PlotMass[estim].FindBin(500),PlotMass[estim].GetNbinsX()+1)
    integr = PlotMass[estim].Integral(1,PlotMass[estim].GetNbinsX()+1)
    print(estim,"fraction > 300 GeV = ",tail300/integr,"500",tail500/integr)

    PlotMass[estim].GetXaxis().SetTitle("Mass [GeV/c^{2}]")
    PlotMass[estim].SetLineWidth(lineWidth)
    PlotMass[estim].SetLineColor(color)
    legEstim.AddEntry(PlotMass[estim],estim,'l')
    cMassEstim.cd()
    if color==1: PlotMass[estim].Draw()
    else: PlotMass[estim].Draw("same")
    color+=1

cMassEstim.cd()
legEstim.Draw("same")


