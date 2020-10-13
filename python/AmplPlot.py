from ROOT import TH1F, TCanvas, TText, TPaveText, TLegend, gROOT, gStyle
from time import sleep
gStyle.SetOptStat(0)
gROOT.SetBatch(0)

#filename="../bin/Ampls.csv"
filename="../bin/Ampls_deuteron.csv"
filename="../bin/Ampls_deuteron_all.csv"
filename="../bin/AmplsStop1600.csv"
filename="../bin/bkg.csv"

#(3.61*pow(10,-9)*247)*1000*120/0.03
#about 40 ADC for 100 microns


def AmplPathBin(ampls):
    norm=30.
    #threshold should depends on pitch, calibration
    threshold = 30.
    return [norm*el/threshold if el>threshold else 0 for el in ampls]

def IsANormalTrack(ampls):
    threshold = 0.25
    size = len(ampls)
    total = sum(ampls)
    vmax = max(ampls)
    outputs = {"charge":0,"leading":0,"maxNeigh":0,"ratio1":0,"ratio2":0}
    outputs["charge"] = total
    outputs["leading"] = vmax
    if size<5:
        return -1, outputs
    else: 
        imax = ampls.index(max(ampls))
	vneigh = 0
	if imax>0: 
	    vneigh = ampls[imax-1] 
	if imax<size-1: 
	    if vneigh< ampls[imax+1]: vneigh = ampls[imax+1]
	outputs["maxNeigh"] = vneigh
	outputs["ratio1"] = vneigh/total
	outputs["ratio2"] = (vmax+vneigh)/total
	# target the case of no sharing
	if vneigh<threshold*vmax:
	    if (vmax/total)>(1-threshold): return 1, outputs
	    else: return 0, outputs
	# target the case of sharing
	else:
	    if (vmax+vneigh)/total>(1-threshold): return 2, outputs
	    else: return 0, outputs


def IsMerged(ampls):
    threshold = 0.5
    size = len(ampls)
    if size<=2: return False
    localAmpls = ampls[:] 
    vmax = max(localAmpls)
    imax = localAmpls.index(max(localAmpls))
    #set the max to 0
    localAmpls[imax] = 0
    #search second peak
    vmax2 = max(localAmpls)
    imax2 = localAmpls.index(max(localAmpls))
    # the 2 leading strip are not close to each other
    if abs(imax2-imax)>1:
       # strips between
       midStrips = []
       if imax2>imax:  midStrips = localAmpls[imax+1:imax2]
       if imax2<imax:  midStrips = localAmpls[imax2+1:imax]
       average = sum(midStrips)/len(midStrips)
       if average<threshold*(vmax+vmax2)/2 : return True
       else: return False

"""
def FilterXtalk(ampls, xtalk=0.15):
"""       
    

doStatPlot = True
#mass = 1.875
mass = 1600

nbins = 16
hAmpl = TH1F("hAmpl","",nbins,1,nbins)
hxtAmpl = TH1F("hxtAmpl","",nbins,1,nbins)
hPathAmpl = TH1F("hPathAmpl","",nbins,1,nbins)

#statistical plots
PMRmax = 1
hCleanEff_vs_PoverM = TH1F("hCleanEff_vs_PoverM","",50,0,PMRmax)
hCleanXTEff_vs_PoverM = TH1F("hCleanEffXT_vs_PoverM","",50,0,PMRmax)
hCleanEff_vs_Layer = TH1F("hCleanEff_vs_PoverM","",21,1,21)
hCleanXTEff_vs_Layer = TH1F("hCleanEffXT_vs_PoverM","",21,1,21)
hLayer  = TH1F("hLayer","",21,1,21)
hPoverM = TH1F("hPoverM","",50,0,PMRmax)

canvas = TCanvas()
gifname="Ampls.gif"
with open(filename) as ifile:
   for aline in ifile:
      line = aline.split('|')[0]
      xtline = aline.split('|')[1]
      p = float(line.split(",")[0])
      pt = float(line.split(",")[1])
      eta = float(line.split(",")[2])
      layer= int(line.split(",")[3])
      path = float(line.split(",")[4])
      isClean = bool(int(line.split(",")[5]))
      #print(isClean)
      exitcode = int(line.split(",")[6])
      isCleanXT = bool(int(line.split(",")[7]))
      ampls = [float(el) for el in line.split(",")[8:]]
      xtampls = [float(el) for el in xtline.split(",")]
      pathAmpls = AmplPathBin(xtampls)
      charge = sum(ampls)
      leg = TLegend(0.15,0.65,0.4,0.88)
      
      if doStatPlot:
          if isClean: hCleanEff_vs_PoverM.Fill(p/mass)
          if isCleanXT: hCleanXTEff_vs_PoverM.Fill(p/mass)
          if isClean: hCleanEff_vs_Layer.Fill(layer)
          if isCleanXT: hCleanXTEff_vs_Layer.Fill(layer)
      
      hLayer.Fill(layer)
      hPoverM.Fill(p/mass)

      #if pt>4:
      #if pt>100 and charge>500 and len(ampls)>6and ( isClean or isCleanXT):
      """
      if pt>100 and charge>500 and ( isClean or isCleanXT):
      #if isClean!=isCleanXT:
	 hPathAmpl.Reset()
	 hxtAmpl.Reset()
	 hxtAmpl.SetLineWidth(2)
	 hAmpl.Reset()
	 hAmpl.SetFillColor(0)
	 hAmpl.SetLineColor(2)
	 hAmpl.SetLineWidth(2)
	 #hPathAmpl.SetLineColor(3)
	 hPathAmpl.SetFillColor(1)
	 hPathAmpl.SetFillStyle(3344)
	 #hPathAmpl.SetLineColorAlpha(3,0.9)
	 #hPathAmpl.SetFillStyle(3001)
	 hPathAmpl.SetLineWidth(2)
	 nstrips = len(ampls)
	 Q = sum(ampls)
	 xtnstrips = len(xtampls)
	 xtQ = sum(xtampls)
         for ibin in range(nstrips):
	    hAmpl.SetBinContent(nbins/2-nstrips/2+ibin+1,float(ampls[ibin]))
	 for ibin in range(len(xtampls)):
	    hxtAmpl.SetBinContent(nbins/2-xtnstrips/2+ibin+1,float(xtampls[ibin]))
	 for ibin in range(len(pathAmpls)):
	    hPathAmpl.SetBinContent(nbins/2-xtnstrips/2+ibin+1,float(pathAmpls[ibin]))
	 if nstrips>1: 
	     isNormal, outputs = IsANormalTrack(ampls)
	     print("Cleaning:",isClean,isCleanXT)
	     #if isNormal: hAmpl.SetFillColor(2)
	     #if IsMerged(ampls): hAmpl.SetFillColor(3)
	     hAmpl.GetYaxis().SetRange(0,500)
	     hAmpl.GetYaxis().SetRangeUser(0,500)
	     hAmpl.Draw()
	     hxtAmpl.Draw("same")
	     hPathAmpl.Draw("same")
	     leg.AddEntry(hAmpl,"cluster","l")
	     leg.AddEntry(hxtAmpl,"xtalk inversion","l")
	     leg.AddEntry(hPathAmpl,"charge compatibility","f")
	     leg.Draw("same")
	     message = "R1 = "+str(outputs["ratio1"]) + " R2 = "+str(outputs["ratio2"])
	     text = TPaveText(10,350,15,420)
	     text.AddText("clean   :"+str(isClean))
	     text.AddText("cleanXT:"+str(isCleanXT))
	     #text.SetTextSize(20);
	     text.Draw("same")
	     canvas.Update()
	     #sleep(3)
	     canvas.WaitPrimitive()
	     print("good !")
	     #canvas.Print(gifname+"+")
	 """
#if gROOT.IsBatch(): 
#canvas.Print(gifname+"++")


#Statistical plots

leg = TLegend(0.6,0.15,0.88,0.35)

hCleanEff_vs_PoverM.GetXaxis().SetTitle("p/M")
hCleanEff_vs_PoverM.GetYaxis().SetTitle("Efficiency")
hCleanEff_vs_PoverM.Divide(hPoverM)
hCleanXTEff_vs_PoverM.Divide(hPoverM)
hCleanEff_vs_Layer.Divide(hLayer)
hCleanXTEff_vs_Layer.Divide(hLayer)
leg.AddEntry(hCleanEff_vs_PoverM,"cleaning wo xtInv","l")
leg.AddEntry(hCleanXTEff_vs_PoverM,"cleaning w xtInv","l")

cEffPoverM = TCanvas("cEffPoverM")
hCleanEff_vs_PoverM.SetLineWidth(2)
hCleanEff_vs_PoverM.SetLineColor(1)
hCleanXTEff_vs_PoverM.SetLineWidth(2)
hCleanXTEff_vs_PoverM.SetLineColor(2)
hCleanEff_vs_PoverM.GetYaxis().SetRangeUser(0,1)
hCleanEff_vs_PoverM.Draw()
hCleanXTEff_vs_PoverM.Draw("same")
leg.Draw("same")

cEffLayer = TCanvas("cEffLayer")
hCleanEff_vs_Layer.GetXaxis().SetTitle("layer")
hCleanEff_vs_Layer.GetYaxis().SetTitle("Efficiency")
hCleanEff_vs_Layer.SetLineWidth(2)
hCleanEff_vs_Layer.SetLineColor(1)
hCleanXTEff_vs_Layer.SetLineWidth(2)
hCleanXTEff_vs_Layer.SetLineColor(2)
hCleanEff_vs_Layer.GetYaxis().SetRangeUser(0,1)
hCleanEff_vs_Layer.Draw()
hCleanXTEff_vs_Layer.Draw("same")
leg.Draw("same")


