#include <numeric>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iterator>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TText.h"
#include "TLegend.h"
#include "TMath.h"
#include "Rtypes.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TExec.h"

#include "../interface/Builder.h"
#include "../interface/Estimator.h"
#include "../interface/PlotterHisto.h"
#include "../interface/Labellizer.h" 
#include "../interface/Calibration.h"
#include "../interface/CommonFunctions.h"

using namespace std;

const float C = 2.51; //obtenus avec fit du proton dqcorr
const float K = 2.61;

//const float C = 3.20573; //obtenus avec fit eloss susy
//const float K = 1.91189;

//MIP: mvp ~ 2.8 - 
//qref = 3.5 (75% distrib)
//qref = 4.15 (85% distrib)
float FractionLowCharge(const vector<float>& v, float qref = 3.5){
   int n = 0;
   for(unsigned int i=0;i<v.size();i++){
     if(v[i]<qref) n++;
   }
   if (v.size()>0) return n*1./v.size();
   return -1;
}

void Stat(const vector<float>& v,float& mean, float& stddev){
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  mean = sum / v.size();
  //cout<<"sum = "<<v[0]<<endl; 
  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(),
  std::bind2nd(std::minus<double>(), mean));
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  stddev = std::sqrt(sq_sum / v.size());

}

void DumpCSV(ofstream& ofile,  Builder* b1, int track){
 float pt                = b1->GetVectTrack()[track].GetPt();
 float p                 = b1->GetVectTrack()[track].GetP();
 float eta               = b1->GetVectTrack()[track].GetEta();

 ofile<<b1->GetRunNb()<<","<<b1->GetEvtNb()<<",";
 ofile<<b1->GetNPV()<<","<<b1->GetInstLumi()<<",";
 ofile<<p<<","<<pt<<","<<eta<<",";
 for(unsigned int cluster = 0;cluster<b1->GetVectTrack()[track].GetVectClusters().size();cluster++){
   ofile<<",";
   ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel()<<",";
   ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPathLength()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetShape()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].IsClean()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStripBoth()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip()<<",";
    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDetId();
  }
  ofile<<endl;
}


TCanvas* CreateCanvasFor2D(const char*title, const char* xtitle, const char*ytitle,TH2D* hlowpt, TH2D* hhighpt, TH2D* hhighdedx){
  TCanvas* c = new TCanvas(title);
  c->Divide(2,2);
  hlowpt->GetXaxis()->SetTitle(xtitle);
  hhighpt->GetXaxis()->SetTitle(xtitle);
  hhighpt->GetYaxis()->SetTitle(ytitle);
  hlowpt->GetYaxis()->SetTitle(ytitle);
  hhighdedx->GetXaxis()->SetTitle(xtitle);
  hhighdedx->GetYaxis()->SetTitle(ytitle);
  c->cd(1);
  c->cd(1)->SetLogz();
  hlowpt->Draw("COLZ");
  c->cd(2);
  //c->SetLogz();
  hhighpt->Draw("COLZ");
  c->cd(3);
  //c->SetLogz();
  hhighdedx->Draw("COLZ");
  return c;
}

int main(int argc,char** argv){

    TChain chain;
	chain.SetName("stage/ttree");
	for(int i=1;i<argc;i++) chain.Add(argv[i]);

    string s1 = argv[1];
    string s2 = s1.substr(0,s1.find('.'))+"_results.root";
    //string s3 = s1.substr(0,s1.find('.'))+"_FitRes.root";


    Builder* b1 = new Builder(chain);
    b1->SetBranchAdd();
    int nentries = b1->GetEntries();

    b1->SetThresholdPartId(0.); //on veut un minimum de 70% de simhits identiques
    b1->SetThresholdPt(0); //on veut des traces avec un minimum de 60 GeV en pt

    ofstream ofile("deuteron.csv");
    ofstream ofile2("bkg.csv");

    //Create the plots
    TH1F hMean("hMean","",200,0,40);
    TH1F hMeanTrunc("hMeanTrunc","",200,0,40);
    TH1F hMeanH2("hMeanH2","",200,0,40);
    TH2F hMeanVsTrunc_Sig("hMeanVsTrunc_Sig","",200,0,40,200,0,40);
    TH2F hMeanVsTrunc_Bkg("hMeanVsTrunc_Bkg","",200,0,40,200,0,40);
    TH1F hStddev("hStddev","",200,0,10);
    TH2D hStddevVsMean_lowpt("hStddevVsMean_lowpt","",200,0,10,200,0,20);
    TH2D hStddevVsMean_highpt("hStddevVsMean_highpt","",200,0,10,200,0,20);
    TH2D hStddevVsMean_highDeDx("hStddevVsMean_highDeDx","",100,5,10,200,0,20);
    TH2D hfQLowVsMean_lowpt("hfQLowVsMean_lowpt","",200,0,10,20,0,1);
    TH2D hfQLowVsMean_highpt("hfQLowVsMean_highpt","",200,0,10,20,0,1);
    TH2D hfQLowVsMean_highDeDx("hfQLowVsMean_highDeDx","",80,6,10,20,0,1);
    TH2D hfQHighVsMean_lowpt("hfQHighVsMean_lowpt","",200,0,10,20,0,1);
    TH2D hfQHighVsMean_highpt("hfQHighVsMean_highpt","",200,0,10,20,0,1);
    TH2D hfQHighVsMean_highDeDx("hfQHighVsMean_highDeDx","",80,6,10,20,0,1);
    TH1F hP("hP","",200,0,20);
    TH2D hDeDxVsP("hDeDxVsP","",500,0,5,1000,0,50);
    TH2D hDeDxVsP_trunc("hDeDxVsP_trunc","",500,0,5,1000,0,50);
    TH2D hDeDxVsP_h2("hDeDxVsP_h2","",500,0,5,1000,0,50);
    TH1F hQMIP("hQMIP","",1000,0,50);
    TH1F hQHSCP("hHSCP","",1000,0,50);
    TH2D hQHSCPvsPLenght("hQHSCPvsPLenght","",200,0,30,100,0,1);
    TH1F hIasMIP("hIasMIP","",100,0,1);
    TH1F hIasHSCP("hIasHSCP","",100,0,1);
    TH2F hIasVsfQLow("hIasVsfQLow","",10,0,1,10,0,1);
    //TH2F hDeDxVsP("hDeDxVsP","",100,0,5,60,0,30);

    TF1 fBB("fBB","[0]/(x*x)+[1]",0,10);
    fBB.SetParameters(3.5,3.96);

    cout<<"N-entries = "<<nentries<<endl;
    //if(nentries>500000) nentries = 200000;

    float lowpcut = 0.7;
    float highpcut = 200;

    float cut_stddev = 7;
    float cut_fQlow = 0.05;
    float cut_fQhigh = 0.3;
    float cut_mean = 6;
    float cut_ias = 0.2;

    float EffSig_Stddev = 0;
    float EffSig_fQLow = 0;
    float EffSig_fQHigh = 0;
    float EffBkg_Stddev = 0;
    float EffBkg_fQLow = 0;
    float EffBkg_fQHigh = 0;
    
    float Eff_1 = 0;
    float Eff_2 = 0;
    float Eff_3 = 0;
    
    int nsig = 0;
    int nbkg = 0;

    for(int i=0;i<nentries;i++)
    {
        if(i%1000==0) cout<<"Event "<<i<<endl;
        b1->GetEntry(i);
	//cout<<b1->GetNtracks()<<endl;
        for(int track=0;track<b1->GetNtracks();track++)
        {
		//cout<<"there"<<endl;
            vector<int> vect_partID;
            vector<float> vect_eloss;
            vector<float> vect_charge;
            vector<float> vect_dedx;
            vector<float> vect_dqdx;
            vector<float> vect_dqcorrdx;
            float pt                = b1->GetVectTrack()[track].GetPt();
            float p                 = b1->GetVectTrack()[track].GetP();
            float eta               = b1->GetVectTrack()[track].GetEta();
            float phi               = b1->GetVectTrack()[track].GetPhi();
            int NCluster            = b1->GetVectTrack()[track].GetNCluster();
            int NClustSat254        = b1->GetVectTrack()[track].GetNSatCluster(254);
            int NClustSat255        = b1->GetVectTrack()[track].GetNSatCluster(255);
            float RatioNClusterSat254 = (double)NClustSat254/(double)NCluster;
            float RatioNClusterSat255 = (double)NClustSat255/(double)NCluster;
            int id                  = b1->GetVectTrack()[track].GetPartId();
            float PoverM            = GetPoverM(p,id);
	    float Ias = b1->GetVectTrack()[track].GetIas();
            //cout<<b1->GetVectTrack()[track].GetNCluster()<<endl;
	    for(int cluster=0;cluster<b1->GetVectTrack()[track].GetNCluster();cluster++)
            {
                float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
                float eloss         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
                //cout<<"QE = "<<charge<<" "<<eloss<<endl;
		vect_eloss.push_back(eloss);
	    }

	    //if(fabs(eta)<1.4) continue ;
	    Estimator estim(vect_eloss);
	    float meanTrunc = estim.GetTrunc40();
	    float meanH2 = estim.GetHarmonic2();
	    //cout<<"size: "<<vect_eloss.size()<<endl;
	    float mean=0, stddev=0;
	    if(vect_eloss.size()>0){
	    	//cout<<"val = "<<vect_eloss[0]<<endl;
		Stat(vect_eloss,mean,stddev);
	    	float fQlow = FractionLowCharge(vect_eloss);
	    	float fQhigh = 1-FractionLowCharge(vect_eloss,mean*1.25);
		hMean.Fill(mean);
		hMeanTrunc.Fill(meanTrunc);
		hMeanH2.Fill(meanH2);
	    	hStddev.Fill(stddev);
	    	//cout<<mean<<" "<<stddev<<endl;
		float weight = 1;
		if(p<lowpcut) weight = 20;//=37
	    	if(p<lowpcut) {
			hMeanVsTrunc_Sig.Fill(meanTrunc,mean);
			//hDeDxVsP.Fill(mean,p);
			hStddevVsMean_lowpt.Fill(mean,stddev);
	    		hfQLowVsMean_lowpt.Fill(mean,fQlow);
	    		hfQHighVsMean_lowpt.Fill(mean,fQhigh);
			if(mean>cut_mean){
			 if(stddev<cut_stddev) EffSig_Stddev++;
			 if(fQlow<cut_fQlow) EffSig_fQLow++;
			 if(fQhigh<cut_fQhigh) EffSig_fQHigh++;
			 nsig++;
			}
		}
		if(p>highpcut){
			hStddevVsMean_highpt.Fill(mean,stddev);
	    		hfQLowVsMean_highpt.Fill(mean,fQlow);
	    		hfQHighVsMean_highpt.Fill(mean,fQhigh);
			if(mean>cut_mean){
			 if(stddev<cut_stddev) EffBkg_Stddev++;
			 if(fQlow<cut_fQlow) EffBkg_fQLow++;
			 if(fQhigh<cut_fQhigh) EffBkg_fQHigh++;
	    		 nbkg++;
			 if(fQlow<cut_fQlow) Eff_1++;
			 if(Ias>cut_ias) Eff_2++;
			 if(Ias>cut_ias && fQlow<cut_fQlow) Eff_3++;
			}
		  DumpCSV(ofile2,b1,track);
		  /*
		  ofile2<<meanTrunc<<","<<meanH2<<","<<mean<<","<<p<<","<<eta;
		  for(unsigned int cluster = 0;cluster<b1->GetVectTrack()[track].GetVectClusters().size();cluster++){
		    bool cleaned = 
		    ofile2<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPathLength()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetShape()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStripBoth()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip()<<",";
		    ofile2<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDetId();
		  }
		  */
			
		}
		if(mean>5 && (p<lowpcut || p>highpcut)){
			hStddevVsMean_highDeDx.Fill(mean,stddev,weight);
	    		hfQLowVsMean_highDeDx.Fill(mean,fQlow,weight);
	    		hfQHighVsMean_highDeDx.Fill(mean,fQhigh,weight);
	    	}
	//	cout<<mean<<" "<<estim.GetTrunc40()<<" "<<estim.GetHarmonic2()<<endl;
		hDeDxVsP.Fill(p,mean);
		hDeDxVsP_trunc.Fill(p,estim.GetTrunc40());
		hDeDxVsP_h2.Fill(p,estim.GetHarmonic2());
		hP.Fill(p);
		if (p>10 && mean>3 && mean<3.5){
		  hIasMIP.Fill(b1->GetVectTrack()[track].GetIas());
		  for(int i=0;i<vect_eloss.size();i++) 
			hQMIP.Fill(vect_eloss[i]);
		}
		//if (p<0.75 && mean>9){
		if (p<1.5 && estim.GetTrunc40()>fBB.Eval(p)){
		//if (p>200 && mean>6){
		  DumpCSV(ofile,b1,track);
		  /*
		  ofile<<meanTrunc<<","<<meanH2<<","<<mean<<","<<p<<","<<eta;
		  for(unsigned int cluster = 0;cluster<b1->GetVectTrack()[track].GetVectClusters().size();cluster++){
		    ofile<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPathLength()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetShape()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStripBoth()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip()<<",";
		    ofile<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDetId();
		  }
		  ofile<<endl;
		  */
		  hMeanVsTrunc_Bkg.Fill(meanTrunc,mean);
		  hIasVsfQLow.Fill(fQlow,Ias);
		  hIasHSCP.Fill(b1->GetVectTrack()[track].GetIas());
		//if (p<1 && mean>9){
		  for(int i=0;i<vect_eloss.size();i++){ 
			hQHSCP.Fill(vect_eloss[i]);
		        /*
			cout<<b1->GetVectTrack()[track].GetVectClusters()[i].GetEloss()<<"\t";
		        cout<<b1->GetVectTrack()[track].GetVectClusters()[i].GetDedxCharge()<<"\t";
		        cout<<b1->GetVectTrack()[track].GetVectClusters()[i].GetSclusCharge()<<"\t";
		        cout<<b1->GetVectTrack()[track].GetVectClusters()[i].GetPathLength()<<endl;
		  	*/
			hQHSCPvsPLenght.Fill(vect_eloss[i],b1->GetVectTrack()[track].GetVectClusters()[i].GetPathLength()*10);
		  }
		}
	    }
        }
    }
    ofile.close();
    ofile2.close();
    cout<<"#track p<0.75\t"<<hStddevVsMean_lowpt.GetEntries()<<endl;
    cout<<"#track p>10\t"<<hStddevVsMean_highpt.GetEntries()<<endl;

    EffSig_Stddev/=nsig;
    EffSig_fQLow/=nsig;
    EffSig_fQHigh/=nsig;
    EffBkg_Stddev/=nbkg;
    EffBkg_fQLow/=nbkg;
    EffBkg_fQHigh/=nbkg;
    Eff_1/=nbkg;
    Eff_2/=nbkg;
    Eff_3/=nbkg;

    cout<<"##########################"<<endl;
    cout<<"##    REPORT            ##"<<endl;
    cout<<"# Stddev "<<endl;
    cout<<"Eff(sig) = "<<EffSig_Stddev<<endl;
    cout<<"Eff(bkg) = "<<EffBkg_Stddev<<endl;
    cout<<"# fQLow "<<endl;
    cout<<"Eff(sig) = "<<EffSig_fQLow<<endl;
    cout<<"Eff(bkg) = "<<EffBkg_fQLow<<endl;
    cout<<"# fQHigh "<<endl;
    cout<<"Eff(sig) = "<<EffSig_fQHigh<<endl;
    cout<<"Eff(bkg) = "<<EffBkg_fQHigh<<endl;
    cout<<"Eff(Ias) = "<<Eff_1<<endl;
    cout<<"Eff(fQlow) = "<<Eff_2<<endl;
    cout<<"Eff(Ias+fQlow) = "<<Eff_3<<endl;
    cout<<"##########################"<<endl;

    //Write output
    TFile fout("results.root","RECREATE");
    fout.cd();
    hMean.Write();
    hMeanTrunc.Write();
    hMeanH2.Write();
    hMeanVsTrunc_Sig.Write();
    hMeanVsTrunc_Bkg.Write();
    hStddev.Write();
    hStddevVsMean_lowpt.Write();
    hStddevVsMean_highpt.Write();
    hStddevVsMean_highDeDx.Write();
    hfQLowVsMean_lowpt.Write();
    hfQLowVsMean_highpt.Write();
    hfQLowVsMean_highDeDx.Write();
    hfQHighVsMean_lowpt.Write();
    hfQHighVsMean_highpt.Write();
    hfQHighVsMean_highDeDx.Write();
    hP.Write();
    hDeDxVsP.Write();
    hDeDxVsP_trunc.Write();
    hDeDxVsP_h2.Write();
    hQMIP.Write();
    hQHSCP.Write();
    hQHSCPvsPLenght.Write();
    hIasMIP.Write();
    hIasHSCP.Write();
    hIasVsfQLow.Write();   
    //hDeDxVsP.Write();

    TCanvas* cfQLowVsMean = CreateCanvasFor2D("cfQLowVsMean","<dEdx> [MeV.cm#^{2}/g]","fraction cluster dEdX<3.5 MeV/g/cm#^2",&hfQLowVsMean_lowpt,&hfQLowVsMean_highpt,&hfQLowVsMean_highDeDx);
    TCanvas* cfQHighVsMean = CreateCanvasFor2D("cfQHighVsMean","<dEdx> [MeV.cm#^2/g]","fraction cluster dEdX> <dEdX>*1.25 MeV.cm#^2/g",&hfQHighVsMean_lowpt,&hfQHighVsMean_highpt,&hfQHighVsMean_highDeDx);
    TCanvas* cStddevVsMean = CreateCanvasFor2D("cStddevVsMean","<dEdx> [MeV.cm#^2/g]","std-dev(dEdX)",&hStddevVsMean_lowpt,&hStddevVsMean_highpt,&hStddevVsMean_highDeDx);
    cfQLowVsMean->Write(); 
    cfQHighVsMean->Write(); 
    cStddevVsMean->Write(); 
    fout.Write();
    fout.Close();

} 
