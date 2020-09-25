import math as m

class Cluster:

    def __init__(self):
	self.layer = 0
	self.pathlength = 0
	self.eloss = 0
	self.shape = True
	self.edge = False
	self.cut = False
	self.clean = True
	self.nsat = 0
	self.nstrips = 0
	self.detid = 0

    def IsGood(self):
    	return shape and clean and not self.edge and not self.cut 

    def IsSat(self):
    	return self.nsat>0
	
    def IsBarrel(self):
    	return 1<=self.layer<=10

    def IsEndCap(self):
    	return 11<=self.layer<=21
	
    def IsTIB(self):
	return 1<=self.layer<=4

    def IsTOB(self):
	return 5<=self.layer<=10

    def IsTID(self):
	return 11<=self.layer<=12

    def IsTEC(self):
	return 13<=self.layer<=21

    def IsBarrelThinSensor(self):
    	return self.layer in [1,2,5,6]

    def IsBarrelThickSensor(self):
    	return self.layer in [3,4,7,8,9,10]

class Event:
    def __init__(self):
        self.EvtNb = 0
	self.RunNb = 0
	self.Lumi = 0
	self.NPV = 0

class HSCPCand:
    
    def __init__(self,_p,_pt,_eta,_clusters):
        self.clusters = _clusters
	self.p = _p
	self.pt = _pt
	selt.eta = _eta

    def Theta():
       return m.atan(m.exp(-eta))*2

class Estimator:
   
    def Estimator(self,clusters):
	clusters = _clusters
	clusters.sort()

    def Check(lowFrac,highFrac):
        if lowFrac<0 or lowFrac>1: return False
        if highFrac<0 or highFrac>1: return False
	if highFrac<lowFrac: return False
	return True

    def GetMean(self, lowFrac, highFrac):
	 if not Check(lowFrac,highFrac): return 0
	 N = len(cluster)
	 if N<=0: return 0
	 coll = self_.clusters[int(lowFrac/N):int(highFrac/N)]
	 mean = sum(coll)/N
	 stddev = pow(sum([c**2 for c in coll])-mean**2,0.5)
	 return mean,stdev

    def GetHarmonic(self,power):
	 if not Check(lowFrac,highFrac): return 0
	 N = len(cluster)
	 if N<=0: return 0
	 coll = self_.clusters[int(lowFrac/N):int(highFrac/N)]
	 mean = pow(sum([pow(c,-power)for c in coll])/N,-1./pow)
	 stddev = pow(sum([ (c-mean)**2 for c in coll])/N,0.5)
	 return mean,stddev


