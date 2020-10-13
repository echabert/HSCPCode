from ROOT import TFile, TF1, TMarker, TGraph

fBB = TF1("fBB","[0]+[1]/(x*x)*[2]",0.3,1);
fBB.SetParameters(5.74295715927,1.05439446652,1.86)
fBB.Draw()
sensors = [0.032]*4+[0.05]*10
p = 0.6
mbudgetFactor = 16
normFactor = 1.4 #para Ih / Eloss
markers = []
gr = TGraph(len(sensors)-1)
for i in range(len(sensors)):
    print(fBB.Eval(p)*sensors[i]*1e-3*mbudgetFactor*normFactor)
    p-=fBB.Eval(p)*sensors[i]*1e-3*mbudgetFactor*normFactor
    #m = TMarker(p,fBB.Eval(p),1)
    print(i,p,fBB.Eval(p))
    gr.SetPoint(i,p,fBB.Eval(p))
    #markers.append(m)
    #m.Draw("same")
    print(p)
gr.Draw("same")

    
