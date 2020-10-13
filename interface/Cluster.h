#include <vector>
#include <iostream>
#include "TProfile.h"

#include "CommonFunctions.h"
#include "ClusterStrip.h"
#include "SimHit.h"

using namespace std;

class Cluster{

	private:
		
		float	dedx_charge_;
		float	sclus_charge_;
		float	pathlength_;
		float eloss_;
		int	nstrips_;
		int	nsimhits_;
		int detid_;
		int subdetid_;
		int firstsclus_;
		bool sat254_;
		bool sat255_;
		bool shape_;
		bool clean_;
		bool cleanXT_;
		uint8_t exitcode_;
		vector<int>	ampls_;
		vector<ClusterStrip>	VectStrips_; //certaintly useless - could have been ampls only
		vector<SimHit>	VectSimHits_;
		int partid_;

	public:	
		
		Cluster();
		Cluster(float dedx_charge,float sclus_charge,float pathlength,float eloss,int nstrips,int nsimhits,int detid, int subdetid,bool sat254,bool sat255,bool shape,int firstsclus,int partId,const vector<ClusterStrip> &VectStrip,const vector<SimHit> &VectSimHit);
		~Cluster();
		bool IsClean() const;
		bool IsCleanXT() const;
		float GetDedxCharge() const;
		float GetSclusCharge() const;
		float GetPathLength() const;
		float GetEloss() const;
		int GetFirstsclus() const;
		int GetLayer() const;
		int GetLayerLabel() const;
		int GetNStripLayer(int layerlabel) const;
		bool GetSat254() const;
		bool GetSat255() const;
		bool GetShape() const;
		int GetNSatStrip(int sat) const;
		int GetNSatStripBoth() const;
		int GetNStrip() const;
		int GetNSimHits() const;
		const vector<ClusterStrip>& GetVectStrips() const;
		vector<int> GetAmpls() const;
		const vector<SimHit>& GetVectSimHits() const;
		TProfile& GetDistribStrip() const;
		bool Edge() const;
		bool Cut() const;
		void SetPartId(int partid);
		int GetPartId() const;
		int GetDetId() const;
		int GetSubDetId() const;
		int GetMaxStrip() const;
		uint8_t GetExitCode() const{return exitcode_;}

		void SetClean(bool clean);
		void SetCleanXT(bool clean);
		void SetExitCode(uint8_t exitcode) {exitcode_ = exitcode;}

};
