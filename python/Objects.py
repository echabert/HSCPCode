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
        self.xtalkVar = -1


    def IsGood(self):
    	#return self.shape and self.clean and not self.edge and not self.cut 
    	#return self.clean and not self.edge and not self.cut 
    	return self.clean 

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
        self.eta = _eta
        self.eloss = [c.eloss for c in self.clusters]
        self.eloss.sort()

    def GetClusterEloss(self,filtered = False):
        if filtered: 
            #print(len(self.clusters),len([c.eloss for c in self.clusters if c.IsGood()])," !!")
	    #return [c.eloss for c in self.clusters if c.IsGood()]
	    return [c.eloss for c in self.clusters if c.xtalkVar<0.2]
        else: 
            return self.eloss

    def Theta():
       return m.atan(m.exp(-eta))*2

class Estimator:
   
    def __init__(self,_clusters):
        self.clusters = _clusters[:]
        self.clusters.sort()

    def Check(self,lowFrac,highFrac):
        if lowFrac<0 or lowFrac>1: return False
        if highFrac<0 or highFrac>1: return False
        if highFrac<lowFrac: return False
        return True

    def GetMean(self, lowFrac=0., highFrac=1.):
        if not self.Check(lowFrac,highFrac): return 0
        N = len(self.clusters)
        if N<=0: return [0,0]
        coll = self.clusters[int(lowFrac*N):int(highFrac*N)]
	N = len(coll)
        if N<=0: return [0,0]
	mean = sum(coll)/N
        stddev = pow(sum([c**2 for c in coll])/N-mean**2,0.5)
        return mean,stddev

    def GetHarmonic(self,lowFrac=0,highFrac=1,power=2):
        if not self.Check(lowFrac,highFrac): 
            return 0
        N = len(self.clusters)
        if N<=0: return [0,0]
        coll = self.clusters[int(lowFrac*N):int(highFrac*N)]
        N = len(coll)
        if N<=0: return [0,0]
        mean = pow(sum([pow(c,-power)for c in coll])/N,-1./power)
        stddev = pow(sum([ (c-mean)**2 for c in coll])/N,0.5)
        return mean,stddev


