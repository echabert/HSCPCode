#include "../interface/Track.h"

#include <vector>

using namespace std;

Track::Track()
{
	pt_			= .0;
	p_			= .0;
	eta_		= .0;
	phi_		= .0;
	nhits_		= .0;
	ndedxhits_	= 0;
	chi2_ = 0;
	missing_ = 0;
	validfraction_ = 0;
	validlast_ = 0;
}

Track::Track(float pt,float p,float eta,float phi,float nhits,int ndedxhits, float chi2, float missing, float validFraction, float validLast,const vector<Cluster> &VectClusters)
{
	pt_		= pt;
	p_		= p;
	eta_	= eta;
	phi_	= phi;
	nhits_	= nhits;
	ndedxhits_	= ndedxhits;
	chi2_ = chi2;
	missing_ = missing;
	validfraction_ = validFraction;
	validlast_ = validLast;
	VectClusters_	= VectClusters;
}

Track::~Track()
{
	pt_			= .0;
	p_			= .0;
	nhits_		= .0;
	ndedxhits_	= 0;
	VectClusters_.clear();
}

float Track::GetPt() const
{
	return pt_;
}

float Track::GetP() const
{
	return p_;
}

float Track::GetEta() const
{
	return eta_;
}

float Track::GetPhi() const
{
	return phi_;
}

int Track::GetNhits() const
{
	return nhits_;
}

int Track::GetNCluster() const
{
	return VectClusters_.size();
}

int Track::GetNSatCluster(int sat) const
{
	if(sat==254)
	{
		int nsat254=0;
		for(int i=0;i<VectClusters_.size();i++)
		{
			if(VectClusters_[i].GetSat254()==true && VectClusters_[i].GetSat255()==false) nsat254++;
		}
		return nsat254;
	}
	if(sat==255)
	{
		int nsat255=0;
		for(int i=0;i<VectClusters_.size();i++)
		{
			if(VectClusters_[i].GetSat255()==true) nsat255++;
		}
		return nsat255;
	}
	else return 0;
}

int Track::GetNSatClusterBoth() const
{
	int nsat=0;
	for(int i=0;i<VectClusters_.size();i++)
	{
		if(VectClusters_[i].GetSat254()==true) nsat++;
	}
	return nsat;
}

const vector<Cluster>& Track::GetVectClusters() const
{
	return VectClusters_;
}

TCanvas& Track::GetProfCluster() const
{
	int size = VectClusters_.size();
	TCanvas* c_track = new TCanvas("c_track","c_track",1200,800*size);
	gStyle->SetOptStat(0);
	TProfile* profDistribStrip = new TProfile();
	c_track->Divide(1,size);
	for(int i=0;i<size;i++)
	{
		string str = to_string(i+1);
		profDistribStrip=&VectClusters_[i].GetDistribStrip();
		c_track->cd(i+1);
		profDistribStrip->GetYaxis()->SetRangeUser(0,300);
		profDistribStrip->SetTitle(str.c_str());
		profDistribStrip->Draw();
	}
	return *c_track;
}

int Track::GetPartId() const
{
	return PartId_;
}

void Track::SetPartId(int PartId)
{
	PartId_ = PartId;
}

float Track::GetIas()const{ return ias_;}

void Track::SetIas(float ias){ias_ = ias;}
