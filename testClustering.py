# investigate shower development based on RecHits and SimClusters
import ROOT
import numpy as np
import HGCalHelpers
from HGCalImagingAlgo import *

## basic setup for testing
# 2D clustering
ecut = 0.060 # need to be changed to ecut = 3. in case it's given wrt the noise and not in absolute units
deltac = [2.,2.,2.]
# multi-clustering
multiclusterRadii = [.015,.015,.015] # it's in cartesian coordiantes, per detector, used to be one value for all dets in eta/phi coordinates: 0.015
minClusters = 3
# allowed events/layers for testing/histograming
allowedRangeLayers = [] # layer considered for histograming e.g. [10, 15], empty for none
allowedRangeEvents = list(range(0,8)) # event numbers considered for histograming, e.g. [5,6,7], empty for none
# other
verbosityLevel = 1 # 0 - only basic info (default); 1 - additional info; 2 - detailed info printed

def getRecHitDetIds(rechits):
    recHitsList = []
    for rHit in rechits:
        recHitsList.append(rHit.detid)
    # print "RecHits -"*10
    # print recHitsList
    recHits = np.array(recHitsList)
    return recHits

def getHitList(simClus, recHitDetIds):
    sClusHitsList = []
    for DetId in simClus.hits:
        sClusHitsList.append(DetId)
    # print sClusHitsList
    sClusHits = np.array(sClusHitsList)
    # thanks to http://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array
    recHitIndices = np.nonzero(np.in1d(recHitDetIds, sClusHits))
    # print "indices - "*10
    # print recHitIndices
    # if (len(recHitIndices[0]) != len(sClusHitsList)):
    #     print "Mismatch:", len(recHitIndices[0]), len(sClusHits)
    return recHitIndices

# get list of rechist associated to sim-cluster hits
def getRecHitsSimAssoc(rechits_raw, simcluster):
    # get sim-cluster associations
    nSimClus = 0
    simClusHitAssoc = []
    recHitDetIds = getRecHitDetIds(rechits_raw)
    for simClusIndex, simClus in enumerate(simcluster):
        simClusHitAssoc.append(getHitList(simClus, recHitDetIds))
        nSimClus += 1

    # get list of rechist associated to simhits
    rHitsSimAssoc = [[] for k in range(0,nSimClus)]
    for simClusIndex, simCl in enumerate(simcluster):
        if (verbosityLevel>=1): print "Sim-cluster index: ",simClusIndex, ", pT: ",simCl.pt, ", E: ",simCl.energy, ", phi: ",simCl.phi, ", eta: ",simCl.eta

        # loop over sim clusters and then rechits
        rHitsSimAssocTemp = []
        for hitIndexArray in simClusHitAssoc[simClusIndex]:
            for hitIndex in hitIndexArray:
                thisHit = rechits_raw[hitIndex]
                if(thisHit.energy < ecut): continue
                # independent of sim cluster, after cleaning
                rHitsSimAssocTemp.append(thisHit)
        rHitsSimAssoc[simClusIndex]= rHitsSimAssocTemp

    return rHitsSimAssoc

# histograming of rechist associated to sim-cluster hits
def histRecHitsSimAssoc(rHitsSimAssoc, currentEvent, histDict, tag = "rHitsAssoc_", zoomed = False):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    if (zoomed): # zoomed for testing/convenience around the eta/phi of most energetic hit
        rs = sorted(range(len(rHitsSimAssoc[0])), key=lambda k: rHitsSimAssoc[0][k].energy, reverse=True) # indices sorted by decreasing rho
        c_phi = rHitsSimAssoc[0][rs[0]].phi
        c_eta = rHitsSimAssoc[0][rs[0]].eta
        histDict[tag+"map_lay_phi_eta_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"map_lay_phi_eta_evt{}".format(currentEvent), tag+"map_lay_phi_eta_evt{};z(cm);#phi;#eta".format(currentEvent), 160, 320, 400, 50, c_phi-0.5, c_phi+0.5, 50, c_eta-0.5, c_eta+0.5) # 3D rechists associated to sim-cluster (with ecut cleaning)
    else:
        histDict[tag+"map_lay_phi_eta_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"map_lay_phi_eta_evt{}".format(currentEvent), tag+"map_lay_phi_eta_evt{};z(cm);#phi;#eta".format(currentEvent), 160, 320, 400, 314, -3.14, +3.14, 320, -3.2, 3.2) # 3D rechists associated to sim-cluster (with ecut cleaning)

    for simClusIndex in range(0,len(rHitsSimAssoc)):
        # define sim-cluster-level hists
        histDict[tag+"map_lay_phi_eta_evt{}_sim{}".format(currentEvent, simClusIndex)]  = ROOT.TH3F(tag+"map_lay_phi_eta_evt{}_sim{}".format(currentEvent, simClusIndex), tag+"map_lay_phi_eta_evt{}_sim{};z(cm);#phi;#eta".format(currentEvent, simClusIndex), 160, 320, 400, 314, -3.14, +3.14, 320, -3.2, 3.2)
        # loop over assoc. rec hits
        for thisHit in rHitsSimAssoc[simClusIndex]:
            histDict[tag+"map_lay_phi_eta_evt{}".format(currentEvent)].Fill(abs(thisHit.z), thisHit.phi, thisHit.eta) # for each sim cluster, after cleaning
            if(thisHit.energy < ecut): continue
            histDict[tag+"map_lay_phi_eta_evt{}_sim{}".format(currentEvent, simClusIndex)].Fill(abs(thisHit.z), thisHit.phi, thisHit.eta) # independent of sim cluster, before cleaning

    return histDict

# histograming of rechists
def histRecHits(rHits, currentEvent, histDict, tag = "rHits_", zoomed = False):
    # sanity check
    if (histDict == None): return

    # define hists per layer
    for layer in range(1, 41):
        if (layer in allowedRangeLayers): # testing limitation
            if (zoomed): # zoomed for testing/convenience around the eta/phi of most energetic hit
                rs = sorted(range(len(rHits)), key=lambda k: rHits[k].energy, reverse=True) # indices sorted by decreasing rho
                c_phi = rHits[rs[0]].phi
                c_eta = rHits[rs[0]].eta
                histDict[tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer)]  = ROOT.TH2F(tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer), tag+"eng_eta_phi_evt{}_lay{};#eta;#phi".format(currentEvent, layer), 40, c_eta-0.1, c_eta+0.1, 40, c_phi-0.1, c_phi+0.1) # 2D energy-weighted-map of raw rechists (with ecut cleaning)
            else:
                histDict[tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer)]  = ROOT.TH2F(tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer), tag+"eng_eta_phi_evt{}_lay{};#eta;#phi".format(currentEvent, layer), 320, -3.2, 3.2, 314, -3.14, +3.14) # 2D energy-weighted-map of raw rechists (with ecut cleaning)

    # loop over all raw rechits and fill per layer
    for rHit in rHits:
        if (rHit.layer in allowedRangeLayers): # testing limitation
            if(rHit.energy < ecut): continue
            histDict[tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, rHit.layer)].Fill(rHit.eta, rHit.phi, rHit.energy)

    return histDict

# histograming of clustered rechist with algo
def histHexelsClustered(hexelsClustered, currentEvent, histDict, tag = "clustHex_", zoomed = False):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    if (zoomed): # zoomed for testing/convenience around the eta/phi of most energetic hit
        rs = sorted(range(len(hexelsClustered)), key=lambda k: hexelsClustered[k].weight, reverse=True) # indices sorted by decreasing rho
        c_phi = hexelsClustered[rs[0]].phi
        c_eta = hexelsClustered[rs[0]].eta
        histDict[tag+"eng_phi_eta_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"eng_phi_eta_evt{}".format(currentEvent), tag+"eng_phi_eta_evt{};z(cm);#phi;#eta".format(currentEvent), 160, 320, 400, 80, c_phi-0.8, c_phi-0.8, 80, c_eta-0.8, c_eta-0.8) # 3D rechists clustered with algo (with ecut cleaning)
    else:
        histDict[tag+"eng_phi_eta_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"eng_phi_eta_evt{}".format(currentEvent), tag+"eng_phi_eta_evt{};z(cm);#phi;#eta".format(currentEvent), 160, 320, 400, 314, -3.14, +3.14, 320, -3.2, 3.2) # 3D rechists clustered with algo (with ecut cleaning)

    # loop over all clustered rechist
    for iNode in hexelsClustered:
        histDict[tag+"eng_phi_eta_evt{}".format(currentEvent)].Fill(abs(iNode.z), iNode.phi, iNode.eta, iNode.weight)

    return histDict

# histograming of abs. difference
def histDiff(fValues, histDict, tag = "RelDiff_", axunit = "#DeltaE(GeV)"):
    # sanity check
    if (histDict == None): return
    
    # define event-level hists
    histDict[tag]  = ROOT.TH1F(tag, tag+";"+axunit, 200, -10, 10) # abs differencein energy

    # loop over all values
    for value in fValues:
        histDict[tag].Fill(value)
    
    return histDict


# histograming of relative difference
def histRelativeDiff(fValues, histDict, tag = "RelDiff_"):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    histDict[tag]  = ROOT.TH1F(tag, tag+";#delta(%)", 200, -1, 1) # relative difference

    # loop over all values
    for value in fValues:
        histDict[tag].Fill(value)
    
    return histDict

# print/save histograms
def histPrintSaveAll(histDict, outDir):
    imgType = "pdf"
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    if (verbosityLevel>=2): print "histDict.items(): ", histDict.items()
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        ROOT.gStyle.SetPalette(ROOT.kBird);
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            item.Draw("colz")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    return

def main():
    # init output stuff
    outDir = "testR"
    HGCalHelpers.createOutputDir(outDir)
    histDict = {}

    # get sample/tree (where algo @reco step was run with parms: ecut = 0.060, deltac = 2., kappa = 10, multiclusterRadii = .015, minClusters = 3)
    inFile = ROOT.TFile.Open("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/partGun_predragm_PDGid22_nPart1_E20_cmssw900pre2_20170116_ReReco/NTUP/partGun_PDGid22_x96_E20.0To20.0_NTUP_1.root")
    chain = inFile.Get("ana/hgc")
    sampleEvents = chain.GetEntries()
    print "Opening TTree: ", chain
    print "Number of events: ", sampleEvents

    # prepare for comparions (should run on large num. of events)
    multiClusters_engDiffAbs = []
    multiClusters_etaDiff = []
    multiClusters_phiDiff = []
    multiClusters_nClust2DDiff = []
    # start event loop
    for currentEvent, event in enumerate(chain):
        if (not currentEvent in allowedRangeEvents): continue # testing limitation
        print "\nCurrent event: ", currentEvent
        
        # get list of rechist associated to sim-cluster hits
        rHitsSimAssoc = getRecHitsSimAssoc(event.rechits_raw, event.simcluster)
        # histograming of rechist associated to sim-cluster hits
        histDict = histRecHitsSimAssoc(rHitsSimAssoc, currentEvent, histDict, tag = "rHitsSimAssoc_", zoomed = False)

        # get list of raw rechits with ecut cleaning
        rHitsCleaned = [rechit for rechit in event.rechits_raw if not rechit.energy < ecut]
        # histograming of raw rechist (with ecut cleaning)
        histDict = histRecHits(rHitsCleaned, currentEvent, histDict, tag = "rHitsCleaned_", zoomed = True)

        ### Imaging algo run at RECO step (CMSSW)
        # get list of all clusters 2D produced with algo at RECO step (CMSSW)
        clusters2DList_reco = [cls2D for cls2D in event.cluster2d]
        # get list of all multi-clusters produced with algo at RECO step (CMSSW)
        multiClustersList_reco = [multiCluster for multiCluster in event.multicluster]

        ### Imaging algo run as stand-alone (python)
        HGCalAlgo = HGCalImagingAlgo(ecut = ecut, deltac = deltac, multiclusterRadii = multiclusterRadii, minClusters = minClusters, verbosityLevel = 0)
        # produce 2D clusters with stand-alone algo, out of all raw rechits
        clusters2D_rerun = HGCalAlgo.makeClusters(event.rechits_raw) # hexels per-layer, per 2D cluster
        # produce multi-clusters with stand-alone algo, out of all 2D clusters
        multiClustersList_rerun = HGCalAlgo.makePreClusters(clusters2D_rerun) # flat list of multi-clusters (as basic clusters), without KDTree
        # multiClustersList_rerun = HGCalAlgo.make3DClusters(clusters2D_rerun) # flat list of multi-clusters (as basic clusters), with KDTree

        # get list of hexeles from 2D clusters produced with stand-alone algo
        clusters2DList_rerun = HGCalAlgo.getClusters(clusters2D_rerun, verbosityLevel = 0) # flat list of 2D clusters (as basic clusters)
        hexelsClustered_rerun = [iNode for bClust in clusters2DList_rerun for iNode in bClust.thisCluster if not iNode.isHalo]
        # histograming of clustered hexels
        histDict = histHexelsClustered(hexelsClustered_rerun, currentEvent, histDict, tag = "clustHex_", zoomed = True)

        ### Compare stand-alone clustering and sim-clusters
        rHitsSimAssocDID = [rechit.detid for simClus in rHitsSimAssoc for rechit in simClus] # list of detids for sim-associated rehits (with ecut cleaning)
        rHitsClustdDID = [iNode.detid for iNode in hexelsClustered_rerun] # list of detids for clustered hexels
        # print some info
        if (verbosityLevel>=1):
            print "num of rechits associated with sim-clusters : ", len (rHitsSimAssocDID)
            print "num of rechits clustered with imaging algo. : ", len (rHitsClustdDID)
            print "num of clustered not found in sim-associated:", len(list(set(rHitsClustdDID)-set(rHitsSimAssocDID )))
            print "num of sim-associated not found in clustered:", len(list(set(rHitsSimAssocDID )-set(rHitsClustdDID)))

        ### Compare stand-alone and reco-level clustering
        clusters2DListMultiSelected_rerun = [cls for multiCluster in multiClustersList_rerun for cls in multiCluster.thisCluster]
        # print some info
        if (verbosityLevel>=1):
            ls = sorted(range(len(clusters2DListMultiSelected_rerun)), key=lambda k: clusters2DListMultiSelected_rerun[k].thisCluster[0].layer, reverse=False) # indices sorted by increasing layer number
#            for index in range(len(clusters2DListMultiSelected_rerun)): print "LayerID: ", clusters2DListMultiSelected_rerun[ls[index]].thisCluster[0].layer, "| 2D-cluster index: ", ls[index], ", No. of cells = ", len(clusters2DListMultiSelected_rerun[ls[index]].thisCluster), ", Energy  = ", clusters2DListMultiSelected_rerun[ls[index]].energy, ", Phi = ", clusters2DListMultiSelected_rerun[ls[index]].phi, ", Eta = ", clusters2DListMultiSelected_rerun[ls[index]].eta, ", z = ", clusters2DListMultiSelected_rerun[ls[index]].z
            for index in range(len(multiClustersList_rerun)): print "Multi-cluster (RE-RUN) index: ", index, ", No. of 2D-clusters = ", len(multiClustersList_rerun[index].thisCluster), ", Energy  = ", multiClustersList_rerun[index].energy, ", Phi = ", multiClustersList_rerun[index].phi, ", Eta = ", multiClustersList_rerun[index].eta, ", z = ", multiClustersList_rerun[index].z
            ls = sorted(range(len(clusters2DList_reco)), key=lambda k: clusters2DList_reco[k].layer, reverse=False) # indices sorted by increasing layer number
#            for index in range(len(clusters2DList_reco)): print "LayerID: ", clusters2DList_reco[ls[index]].layer, "| 2D-cluster index: ", ls[index], ", No. of cells = ", clusters2DList_reco[ls[index]].nhitAll, ", Energy  = ", clusters2DList_reco[ls[index]].energy, ", Phi = ", clusters2DList_reco[ls[index]].phi, ", Eta = ", clusters2DList_reco[ls[index]].eta, ", z = ", clusters2DList_reco[ls[index]].z
            for index in range(len(multiClustersList_reco)): print "Multi-cluster (RECO) index: ", index, ", No. of 2D-clusters = ", multiClustersList_reco[index].nclus, ", Energy  = ", multiClustersList_reco[index].energy, ", Phi = ", multiClustersList_reco[index].phi, ", Eta = ", multiClustersList_reco[index].eta, ", z = ", multiClustersList_reco[index].z
            print "num of clusters2D @reco : ", len(clusters2DList_reco)
            print "num of clusters2D re-run: ", len(clusters2DListMultiSelected_rerun)
            print "num of multi-cluster @reco : ", len(multiClustersList_reco)
            print "num of multi-cluster re-run: ", len(multiClustersList_rerun)
        # get some historams
        if(len(multiClustersList_rerun)==len(multiClustersList_reco)): # diff. in energy, eta, phi, z for re-run and @reco multi-clustres
            multiClusters_engDiffAbs.extend([(multiClustersList_rerun[k].energy - multiClustersList_reco[k].energy) for k in range(0,len(multiClustersList_rerun))]) # abs. diff.
            multiClusters_etaDiff.extend([100*(multiClustersList_rerun[k].eta - multiClustersList_reco[k].eta)/multiClustersList_reco[k].eta for k in range(0,len(multiClustersList_rerun))])
            multiClusters_phiDiff.extend([100*(multiClustersList_rerun[k].phi - multiClustersList_reco[k].phi)/multiClustersList_reco[k].phi for k in range(0,len(multiClustersList_rerun))])
        multiClusters_nClust2DDiff.append(100*(len(clusters2DListMultiSelected_rerun) - len(clusters2DList_reco))/(0.+len(clusters2DList_reco))) # relative diff. in number of 2D clusters

    histDict = histDiff(multiClusters_engDiffAbs, histDict, tag = "MultClust_EngAbsDiffRerunReco", axunit = "#Delta E(GeV)")
    histDict = histRelativeDiff(multiClusters_etaDiff, histDict, tag = "MultClust_EtaDiffRerunReco")
    histDict = histRelativeDiff(multiClusters_phiDiff, histDict, tag = "MultClust_PhiDiffRerunReco")
    histDict = histRelativeDiff(multiClusters_nClust2DDiff, histDict, tag = "MultClust_nClust2DDiffRerunReco")

    # print/save histograms
    histPrintSaveAll(histDict, outDir)

if __name__ == '__main__':
    main()
