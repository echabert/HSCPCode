from Objects import *

def LoadCand(line):
    words = line.split(',')
    # build event
    event = Event()
    event.EvtNb = int(words.pop(0))
    event.RunNb = int(words.pop(0))
    event.NPV = int(words.pop(0))
    event.Lumi = float(words.pop(0))
    # build HSCP cand
    p = float(words.pop(0))
    pt = float(words.pop(0))
    eta = float(words.pop(0))
    clusters = []
    #   collection of clusters
    while len(words)>0:
        cluster = Cluster()
        cluster.layer = int(words.pop(0))
        cluster.pathlength = float(words.pop(0))
        cluster.eloss = float(words.pop(0))
        #important remark: conversion string --> int --> bool. Cannot convert directly string to bool otherwise it is always true !
	cluster.shape = bool(int(words.pop(0)))
        cluster.edge = bool(int(words.pop(0)))
        cluster.cut = bool(int(words.pop(0)))
        cluster.clean =bool(int(words.pop(0)))
        cluster.nsat = int(words.pop(0))
        cluster.nstrips = int(words.pop(0))
        cluster.xtalkVar = float(words.pop(0))
	cluster.detid = int(words.pop(0))
        clusters.append(cluster)
    cand = HSCPCand(p,pt,eta,clusters)
    # return both event and cand
    return event, cand


def LoadCSV(filename):
    cands = []
    with open(filename) as ifile:
        for line in ifile:
            event, cand = LoadCand(line)
            cands.append((event,cand))
    return cands

