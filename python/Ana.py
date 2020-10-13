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
   #print("In stat",len(mylist))
   if len(mylist) == 0: return {}
   mylist.sort()
   N = len(mylist)
   mean = sum(mylist)/N
   #print("mean = ",mean)
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
pmax = 1.4*pfactor #1.5
pstep = 0.1*pfactor
npsteps = int((pmax-pmin)/pstep)
pbins = [pmin+i*pstep for i in range(npsteps+1)]
#bin in eta
etamin = 0
etamax = 1.5
etastep = 0.1
netasteps = int((etamax-etamin)/etastep)
ebins = [etamin+i*etastep for i in range(netasteps+1)]
#nhits
nhitsmin = 5
nhitsmax = 25


#study estimators
estimators = {}
#tuple with (IsH2,IsFiltered,lowFrac,highFrac)
UseFilter=False
UseHarmonic=False
estimators["Mean(0,1)"]=(UseFilter,UseHarmonic,0,1)
estimators["Mean(0.15,1)"]=(UseFilter,UseHarmonic,0.15,1)
#estimators["Mean(0,0.9)"]=(UseFilter,UseHarmonic,0,0.9)
estimators["Mean(0,0.85)"]=(UseFilter,UseHarmonic,0,0.85)
#estimators["Mean(0,0.75)"]=(UseFilter,UseHarmonic,0,0.75)
#estimators["Mean(0,0.65)"]=(UseFilter,UseHarmonic,0,0.65)
estimators["Mean(0.15,0.85)"]=(UseFilter,UseHarmonic,0.15,0.85)
estimators["Mean(0.25,0.75)"]=(UseFilter,UseHarmonic,0.25,0.75)
#estimators["Mean(0,0.6)"]=(UseFilter,UseHarmonic,0,0.6)
#estimators["Mean(0.15,0.7)"]=(UseFilter,UseHarmonic,0.15,0.7)
#estimators["Mean(0.4,0.7)"]=(UseFilter,UseHarmonic,0.4,0.7)

# plotting tuning
dEdxlabel = "<dEdx> [MeV.cm^{2}/g]"
lineWidth = 2

#other constants
deuteronMass = 1.8756 #GeV


#fit results
resFile = open("fitRes.txt",'w')

#file to be analyzed
#filename = "../bin/deuteron.csv"
filename = "deuteron.csv"
#filename = "deuteron_light.csv"
#filename = "stop2000.csv"
#filename = "bkg.csv"

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


hEtaVsP = TH2F("hEtaVsP","",npsteps,pmin,pmax,netasteps,etamin,etamax)
hEta = TH1F("hEta","",netasteps,etamin,etamax)
hP = TH1F("hP","",npsteps,pmin,pmax)
hNhits = TH1F("nHits","",nhitsmax-nhitsmin,nhitsmin,nhitsmax)
for cand in cands:
    hEtaVsP.Fill(cand[1].p, cand[1].eta)
    hEta.Fill(cand[1].eta)
    hP.Fill(cand[1].p)
    hNhits.Fill(len(cand[1].GetClusterEloss()))

cEtaVsP = TCanvas()
hEtaVsP.Draw("COLZ")
print(hEtaVsP.GetCorrelationFactor())
c1Dplots = TCanvas()
c1Dplots.Divide(2,2)
c1Dplots.cd(1)
hP.Draw()
c1Dplots.cd(2)
hEta.Draw()
c1Dplots.cd(3)
hNhits.Draw()


#exit()

###############################################
# Study of the estimators
###############################################



#create dictionnaries
#per momenta bins
Estim_P_list = {}
Estim_P_res = {}
#per momenta, #nhits, eta bins
Estim_all_list = {}
Estim_all_res = {}

for estim in estimators:
    for ibin in range(len(pbins)):
       Estim_P_list[(estim,ibin)] = []
       Estim_P_res[(estim,ibin)] = None
       #add info about #nhits - eta
       for ebin in range(len(ebins)):
          for nhits in range(nhitsmin,nhitsmax+1):
	      Estim_all_list[(estim,ibin,ebin,nhits)] = []
	      Estim_all_res[(estim,ibin,ebin,nhits)] = None

#run over data
for cand in cands:
    #for c in cand[1].clusters:
    for estim in estimators:
        #print("Estimator = ",estimators[estim][2],estimators[estim][3])
        #print("collection size: ",len(cand[1].GetClusterEloss(estimators[estim][1])))
	e = Estimator(cand[1].GetClusterEloss(estimators[estim][0]))
        #print(cand[1].GetClusterEloss(estimators[estim][1]))
	#print(e.GetMean())
	value = 0
	if estimators[estim][1]: value = e.GetHarmonic(estimators[estim][2],estimators[estim][3])
        else: value = e.GetMean(estimators[estim][2],estimators[estim][3])
        #print("value = ",value)
        Estim_P_list[(estim,GetIndice(cand[1].p,pbins))].append(value)
	if 0<=GetIndice(cand[1].eta,ebins)<=len(ebins) and nhitsmin<= len(cand[1].GetClusterEloss(estimators[estim][0])) <nhitsmax:
	    Estim_all_list[(estim,GetIndice(cand[1].p,pbins),GetIndice(cand[1].eta,ebins),len(cand[1].GetClusterEloss(estimators[estim][0])))].append(value)


#print(Estim_all_list)
#exit()

#statistics per estimator
for estim in estimators:
    for ibin in range(len(pbins)):
        #print(Estim_P_list[(estim,ibin)])
        #print(([c[0] for c in Estim_P_list[(estim,ibin)]]))
	Estim_P_res[(estim,ibin)] = Statistics([c[0] for c in Estim_P_list[(estim,ibin)]])
	#print((estim,ibin),Estim_P_res[(estim,ibin)])
        for ebin in range(len(ebins)):
            for nhits in range(nhitsmin, nhitsmax):
	        Estim_all_res[(estim,ibin,ebin,nhits)] = Statistics([c[0] for c in Estim_all_list[(estim,ibin,ebin,nhits)]])
		#print("size = ",len(Estim_all_list[(estim,ibin,ebin,nhits)]))
		#print("THERE ",Estim_all_res[(estim,ibin,ebin,nhits)])

#plots and layout
estimMean_P_plots = []
estimReso_P_plots = []
cMeanEstim = TCanvas()
cResoEstim = TCanvas()
legEstim = TLegend()
color = 1
#fit function
BBfunc = TF1("BBFunc","[0]+[1]/(x*x)",0,2)
for estim in estimators:
    hMean = TH1F("hMeanEstim_"+estim,"",npsteps,pmin,pmax)
    hReso = TH1F("hResoEstim_"+estim,"",npsteps,pmin,pmax)
    for ibin in range(len(pbins)):
        #print(Estim_P_res[(estim,ibin)]["mean"])
	#print((estim,ibin),Estim_P_res[(estim,ibin)])
	hMean.SetBinContent(ibin,Estim_P_res[(estim,ibin)]["mean"])
        hMean.SetBinError(ibin,Estim_P_res[(estim,ibin)]["meanError"])
        hReso.SetBinContent(ibin,Estim_P_res[(estim,ibin)]["stddev"]/Estim_P_res[(estim,ibin)]["mean"])
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
    legEstim.AddEntry(hMean,estim,'l')
    cMeanEstim.cd()
    if color==1: hMean.Draw()
    else: hMean.Draw("same")
    cResoEstim.cd()
    if color==1: hReso.Draw()
    else: hReso.Draw("same")
    color+=1
    estimMean_P_plots.append(hMean)
    estimReso_P_plots.append(hReso)
    #perform BB fit
    resFit = hMean.Fit(BBfunc,"S")
    resFile.write(estim+" "+str(BBfunc.GetParameters()[0])+" "+str(BBfunc.GetParameters()[1]/(deuteronMass*deuteronMass))+" "+str(resFit.Chi2())+"\n")


#check bias due to eta
#assign weight based on  nhist and p bins ...
legBias = TLegend()
plotBias = []
cEtaBias = TCanvas()
color = 1
for estim in estimators:
    hEtaBias = TH1F("hEtaBias_"+estim,"",netasteps,etamin,etamax)
    for ebin in range(len(ebins)):
        #print("eta bin",ebin)
	diff = []
	uncert = []
	totalweight = 0
	for pbin in range(len(pbins)):
            Estim_P_res[(estim,pbin)]  
            for nhits in range(nhitsmin,nhitsmax):
                 if Estim_all_res.get((estim,pbin,ebin,nhits),None) is not None:
		     mean = Estim_all_res[(estim,pbin,ebin,nhits)].get("mean",0)
		     error = Estim_all_res[(estim,pbin,ebin,nhits)].get("meanError",0)
		     weight = len(Estim_all_list[(estim,pbin,ebin,nhits)])
		     if mean>0: 
		         diff.append((mean-Estim_P_res[(estim,pbin)]["mean"])*weight)
		         uncert = [pow(error*weight,2)]
			 totalweight+=weight
	bias = 0
	biasError = 0
	#print(sum(diff),totalweight)
	if totalweight>0:  
	    bias = sum(diff)/totalweight
            biasError = sqrt(sum(uncert)/pow(totalweight,2))
	hEtaBias.SetBinContent(ebin,bias)
	if biasError>0 and biasError<100:	hEtaBias.SetBinError(ebin,biasError)
    hEtaBias.SetLineWidth(lineWidth)
    hEtaBias.SetLineColor(color)
    hEtaBias.GetXaxis().SetTitle("track #eta")
    hEtaBias.GetYaxis().SetTitle("bias of I_{h}")
    legBias.AddEntry(hEtaBias,estim,'l')
    if color==1: hEtaBias.Draw()
    else: hEtaBias.Draw("same")
    color+=1
    plotBias.append(hEtaBias)
legBias.Draw("same")

#check bias due to nhits
cNhitsBias = TCanvas()
color = 1
for estim in estimators:
    nbins = nhitsmax-nhitsmin
    hNhitsBias = TH1F("hNhitsBias_"+estim,"",nbins,nhitsmin,nhitsmax)
    for nhits in range(nhitsmin,nhitsmax):
	diff = []
	uncert = []
	totalweight = 0
	for pbin in range(len(pbins)):
            Estim_P_res[(estim,pbin)]  
            for ebin in range(len(ebins)):
                 if Estim_all_res.get((estim,pbin,ebin,nhits),None) is not None:
		     mean = Estim_all_res[(estim,pbin,ebin,nhits)].get("mean",0)
		     error = Estim_all_res[(estim,pbin,ebin,nhits)].get("meanError",0)
		     weight = len(Estim_all_list[(estim,pbin,ebin,nhits)])
		     if mean>0: 
		         diff.append((mean-Estim_P_res[(estim,pbin)]["mean"])*weight)
		         uncert = [pow(error*weight,2)]
			 totalweight+=weight
	bias = 0
	biasError = 0
	#print(sum(diff),totalweight)
	if totalweight>0:  
	    bias = sum(diff)/totalweight
            biasError = sqrt(sum(uncert)/pow(totalweight,2))
	hNhitsBias.SetBinContent(nhits-nhitsmin+1,bias)
	if biasError>0 and biasError<100:	hNhitsBias.SetBinError(nhits-nhitsmin+1,biasError)
    hNhitsBias.SetLineWidth(lineWidth)
    hNhitsBias.SetLineColor(color)
    hNhitsBias.GetXaxis().SetTitle("#hits")
    hNhitsBias.GetYaxis().SetTitle("bias of I_{h}")
    if color==1: hNhitsBias.Draw()
    else: hNhitsBias.Draw("same")
    color+=1
    plotBias.append(hNhitsBias)
legBias.Draw("same")
resFile.close()
cMeanEstim.cd()
legEstim.Draw("same")
cResoEstim.cd()
legEstim.Draw("same")


