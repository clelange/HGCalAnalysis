from SampleHelper import SampleManager
import ROOT
import logging
import numpy as np
import HGCalHelpers


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
        # print "Mismatch:", len(recHitIndices[0]), len(sClusHits)
    return recHitIndices


def getHists():
    histDict = {}
    clusters = ["SimClus", "PFClus", "GenPart", "RecHits"]
    for clus in clusters:

        histDict["%s_energy" %clus] = ROOT.TH1F("%s_energy" %clus, "%s_energy;E [GeV]" %clus, 200, 0, 100)
        histDict["%s_pt" %clus] = ROOT.TH1F("%s_pt" %clus, "%s_pt;p_{T} [GeV]" %clus, 100, 0, 20)
        histDict["%s_eta" %clus] = ROOT.TH1F("%s_eta" %clus, "%s_eta;#eta" %clus, 100, -5, 5)
        histDict["%s_phi" %clus] = ROOT.TH1F("%s_phi" %clus, "%s_phi;#phi" %clus, 100, -3.2, 3.2)
        if (clus == "SimClus"):
            histDict["%s_simEnergy" %clus] = ROOT.TH1F("%s_simEnergy" %clus, "%s_simEnergy;simE [GeV]" %clus, 200, 0, 50)
            histDict["%s_layers_energy" %clus] = ROOT.TH2F("%s_layers_energy" %clus, "%s_layers_energy;layers;energy [GeV]" %clus, 30, 0, 30, 200, 0, 50)
            histDict["%s_cells_energy" %clus] = ROOT.TH2F("%s_cells_energy" %clus, "%s_cells_energy;cells;energy [GeV]" %clus, 245, 0, 245, 200, 0, 50)
            histDict["%s_wafers_energy" %clus] = ROOT.TH2F("%s_wafers_energy" %clus, "%s_wafers_energy;wafer;energy [GeV]" %clus, 550, 0, 550, 200, 0, 50)
            histDict["%s_fractions_energy" %clus] = ROOT.TH2F("%s_fractions_energy" %clus, "%s_fractions_energy;fraction;energy [GeV]" %clus, 100, 0, 1, 200, 0, 50)
            histDict["%s_layers_pt" %clus] = ROOT.TH2F("%s_layers_pt" %clus, "%s_layers_pt;layers;p_{T} [GeV]" %clus, 30, 0, 30, 200, 0, 50)
            histDict["%s_cells_pt" %clus] = ROOT.TH2F("%s_cells_pt" %clus, "%s_cells_pt;cells;p_{T} [GeV]" %clus, 245, 0, 245, 200, 0, 50)
            histDict["%s_wafers_pt" %clus] = ROOT.TH2F("%s_wafers_pt" %clus, "%s_wafers_pt;wafer;p_{T} [GeV]" %clus, 550, 0, 550, 200, 0, 50)
            histDict["%s_fractions_pt" %clus] = ROOT.TH2F("%s_fractions_pt" %clus, "%s_fractions_pt;fraction;p_{T} [GeV]" %clus, 100, 0, 1, 200, 0, 50)
            histDict["%s_layers_fractions" %clus] = ROOT.TH2F("%s_layers_fractions" %clus, "%s_layers_fractions;layers;fractions" %clus, 30, 0, 30, 100, 0, 1)
            histDict["%s_cells_fractions" %clus] = ROOT.TH2F("%s_cells_fractions" %clus, "%s_cells_fractions;cells;fractions" %clus, 245, 0, 245, 100, 0, 1)
            histDict["%s_wafers_fractions" %clus] = ROOT.TH2F("%s_wafers_fractions" %clus, "%s_wafers_fractions;wafer;fractions" %clus, 550, 0, 550, 100, 0, 1)

        if (clus == "GenPart"):
            histDict["%s_dvz" %clus] = ROOT.TH1F("%s_dvz" %clus, "%s_dvz;dvz [cm]" %clus, 100, 0, 500)

        if (clus == "RecHits"):
            histDict["%s_layers_energy" %clus] = ROOT.TH2F("%s_layers_energy" %clus, "%s_layers_energy;layers;energy [GeV]" %clus, 40, 0.5, 40.5, 200, 0, 30)
            histDict["%s_layers_pt" %clus] = ROOT.TH2F("%s_layers_pt" %clus, "%s_layers_pt;layers;p_{T} [GeV]" %clus, 40, 0.5, 40.5, 200, 0, 10)

        histDict["%s_dRtoSeed" %clus] = ROOT.TH1F(
            "%s_dRtoSeed" %clus, "%s_dRtoSeed;#Delta R to seed" %clus, 100, -0.2, 0.2)

        if (clus != "GenPart"):
            histDict["multi%s_energy" %clus] = ROOT.TH1F("multi%s_energy" %clus, "multi%s_energy;E [GeV]" %clus, 200, 0, 50)
            histDict["multi%s_pt" %clus] = ROOT.TH1F("multi%s_pt" %clus, "multi%s_pt;p_{T} [GeV]" %clus, 100, 0, 10)
            histDict["multi%s_eta" %clus] = ROOT.TH1F("multi%s_eta" %clus, "multi%s_eta;#eta" %clus, 100, -5, 5)
            histDict["multi%s_phi" %clus] = ROOT.TH1F("multi%s_phi" %clus, "multi%s_phi;#phi" %clus, 100, -3.2, 3.2)

    compClusters = ["SimVsPF", "GenVsPF", "GenVsSim", "SimVsRecHits"]
    for comp in compClusters:
        if comp.find("Gen") < 0:
            histDict["%s_delta_energy" %comp] = ROOT.TH1F("%s_delta_energy" %comp, "%s_delta_energy;#Delta E [GeV]" %comp, 100, -20, 20)
            histDict["%s_delta_pt" %comp] = ROOT.TH1F("%s_delta_pt" %comp, "%s_delta_pt;#Delta p_{T} [GeV]" %comp, 100, -10, 10)
            histDict["%s_deltaover_energy" %comp] = ROOT.TH1F("%s_deltaover_energy" %comp, "%s_delta_overenergy;#Delta E/E" %comp, 100, -1, 1)
            histDict["%s_deltaover_pt" %comp] = ROOT.TH1F("%s_deltaover_pt" %comp, "%s_deltaover_pt;#Delta p_{T}/p_{T}" %comp, 100, -1, 1)
            histDict["%s_delta_eta" %comp] = ROOT.TH1F("%s_delta_eta" %comp, "%s_delta_eta;#Delta #eta" %comp, 100, -1, 1)
            histDict["%s_delta_phi" %comp] = ROOT.TH1F("%s_delta_phi" %comp, "%s_delta_phi;#Delta #phi" %comp, 100, -1, 1)
            histDict["%s_delta_R" %comp] = ROOT.TH1F("%s_delta_R" %comp, "%s_delta_R;#Delta R" %comp, 100, -1, 1)

        histDict["multi%s_delta_energy" %comp] = ROOT.TH1F("multi%s_delta_energy" %comp, "multi%s_delta_energy;#Delta E [GeV]" %comp, 100, -10, 10)
        histDict["multi%s_delta_pt" %comp] = ROOT.TH1F("multi%s_delta_pt" %comp, "multi%s_delta_pt;#Delta p_{T} [GeV]" %comp, 100, -5, 5)
        histDict["multi%s_deltaover_energy" %comp] = ROOT.TH1F("multi%s_deltaover_energy" %comp, "multi%s_deltaover_energy;#Delta E/E" %comp, 100, -5, 5)
        histDict["multi%s_deltaover_pt" %comp] = ROOT.TH1F("multi%s_deltaover_pt" %comp, "multi%s_deltaover_pt;#Delta p_{T}/p_{T}" %comp, 100, -5, 5)
        histDict["multi%s_delta_eta" %comp] = ROOT.TH1F("multi%s_delta_eta" %comp, "multi%s_delta_eta;#Delta #eta" %comp, 100, -1, 1)
        histDict["multi%s_delta_phi" %comp] = ROOT.TH1F("multi%s_delta_phi" %comp, "multi%s_delta_phi;#Delta #phi" %comp, 100, -1, 1)
        histDict["multi%s_delta_R" %comp] = ROOT.TH1F("multi%s_delta_R" %comp, "multi%s_delta_R;#Delta R" %comp, 100, -1, 1)
        histDict["multi%s_eff" %comp] = ROOT.TH1F("multi%s_eff" %comp, "multi%s_eff;eff." %comp, 2, -0.5, 1.5)

        if comp.find("Gen") < 0:
            histDict["multi%s_selected_delta_energy" %comp] = ROOT.TH1F("multi%s_selected_delta_energy" %comp, "multi%s_selected_delta_energy;#Delta E [GeV]" %comp, 100, -10, 10)
            histDict["multi%s_selected_delta_pt" %comp] = ROOT.TH1F("multi%s_selected_delta_pt" %comp, "multi%s_selected_delta_pt;#Delta p_{T} [GeV]" %comp, 100, -5, 5)
            histDict["multi%s_selected_deltaover_energy" %comp] = ROOT.TH1F("multi%s_selected_deltaover_energy" %comp, "multi%s_selected_deltaover_energy;#Delta E/E" %comp, 100, -10, 10)
            histDict["multi%s_selected_deltaover_pt" %comp] = ROOT.TH1F("multi%s_selected_deltaover_pt" %comp, "multi%s_selected_deltaover_pt;#Delta p_{T}/p_{T}" %comp, 100, -5, 5)
            histDict["multi%s_selected_delta_eta" %comp] = ROOT.TH1F("multi%s_selected_delta_eta" %comp, "multi%s_selected_delta_eta;#Delta #eta" %comp, 100, -1, 1)
            histDict["multi%s_selected_delta_phi" %comp] = ROOT.TH1F("multi%s_selected_delta_phi" %comp, "multi%s_selected_delta_phi;#Delta #phi" %comp, 100, -1, 1)
            histDict["multi%s_selected_delta_R" %comp] = ROOT.TH1F("multi%s_selected_delta_R" %comp, "multi%s_selected_delta_R;#Delta R" %comp, 100, -1, 1)
            histDict["multi%s_selected_eff" %comp] = ROOT.TH1F("multi%s_selected_eff" %comp, "multi%s_selected_eff;eff." %comp, 2, -0.5, 1.5)
        if comp == "SimVsRecHits":
            detectors = ["EE", "FH", "both"]
            for detect in detectors:
                histDict["%s_delta_energy_%s" %(comp, detect)] = ROOT.TH1F("%s_delta_energy_%s" %(comp, detect), "%s_delta_energy_%s;#Delta E [GeV]" %(comp, detect), 100, -20, 20)
                histDict["%s_delta_pt_%s" %(comp, detect)] = ROOT.TH1F("%s_delta_pt_%s" %(comp, detect), "%s_delta_pt_%s;#Delta p_{T} [GeV]" %(comp, detect), 100, -10, 10)
                histDict["%s_deltaover_energy_%s" %(comp, detect)] = ROOT.TH1F("%s_deltaover_energy_%s" %(comp, detect), "%s_delta_overenergy_%s;#Delta E/E" %(comp, detect), 100, -1, 1)
                histDict["%s_deltaover_pt_%s" %(comp, detect)] = ROOT.TH1F("%s_deltaover_pt_%s" %(comp, detect), "%s_deltaover_pt_%s;#Delta p_{T}/p_{T}" %(comp, detect), 100, -1, 1)
                histDict["%s_frac_energy_%s" %(comp, detect)] = ROOT.TH1F("%s_frac_energy_%s" %(comp, detect), "%s_frac_energy_%s;E fraction" %(comp, detect), 100, -3, 3)
                histDict["%s_frac_pt_%s" %(comp, detect)] = ROOT.TH1F("%s_frac_pt_%s" %(comp, detect), "%s_frac_pt_%s;p_{T} fraction" %(comp, detect), 100, -3, 3)

    return histDict



def main():

    localTest = False
    nEvents = -1
    dvzCut = 0  # 320
    simClusPtCut = 20
    imgType = "pdf"
    outDir = "singlePions"

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)

    sampleManager = SampleManager()
    sample = sampleManager.getSample("chargedPions_nPart1_Pt20")
    chain = sample.getChain()
    if localTest:
        inFile = ROOT.TFile.Open("/afs/cern.ch/user/c/clange/work/HGCal/ntupliser/CMSSW_8_1_0_pre8/src/RecoNtuples/HGCalAnalysis/test/twogamma_pt5_eta2_nosmear_calib_ntuple.root")
        chain = inFile.Get("ana/hgc")

    HGCalHelpers.createOutputDir(outDir)
    histDict = getHists()
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)

    sampleEvents = chain.GetEntries()
    logger.info("Events: %d" % sampleEvents)

    # start event loop
    for currentEvent, event in enumerate(chain):
        if (currentEvent % 100 == 0):
            logger.info("Event {} of {}".format(currentEvent, sampleEvents))

        # associate RecHits to SimClusters
        simClusHitAssoc = []
        recHitDetIds = getRecHitDetIds(event.rechits_raw)
        for simClusIndex,simClus in enumerate(event.simcluster):
            # print ("SimClus %d - "%simClusIndex)*10
            simClusHitAssoc.append(getHitList(simClus, recHitDetIds))

        # get generator particles applying conversion cut
        vGenParticleTLV = HGCalHelpers.getGenParticles(event, histDict, dvzCut)

        # keep track of SimClusters that are in the same hemisphere as the selected GenParticles
        matchesGen = [False]*len(event.simcluster)

        simClusIndex = 0
        for simCl, pfCl in zip(event.simcluster, event.pfcluster):
            # logging.info("{} {}".format(simCl.pt, pfCl.pt))
            for genPartIndex, genPart in enumerate(vGenParticleTLV):
                # if deltaR2(genPart, simCl) < loose_dRCut:
                if genPart.Eta()*simCl.eta  > 0:
                    matchesGen[simClusIndex] = True
                    histDict["SimVsPF_delta_energy"].Fill(simCl.energy-pfCl.energy)
                    histDict["SimVsPF_delta_pt"].Fill(simCl.pt-pfCl.pt)
                    histDict["SimVsPF_deltaover_energy"].Fill((simCl.energy-pfCl.energy)/simCl.energy)
                    histDict["SimVsPF_deltaover_pt"].Fill((simCl.pt-pfCl.pt)/simCl.pt)
                    histDict["SimVsPF_delta_eta"].Fill(simCl.eta-pfCl.eta)
                    histDict["SimVsPF_delta_phi"].Fill(simCl.phi-pfCl.phi)
                    histDict["SimVsPF_delta_R"].Fill(HGCalHelpers.deltaR(simCl, pfCl))
                    break
            simClusIndex += 1

        # use the SimClusters that "match" the GenParticles and study associated RecHits
        detectors = ["EE", "FH", "both"]
        for simClusIndex, simCl in enumerate(event.simcluster):
            # if (matchesGen[simClusIndex]):
            if simCl.pt >= simClusPtCut:
                recHitVectors = {}
                for detect in detectors:
                    recHitVectors[detect] = ROOT.TLorentzVector()
                for hitIndexArray in simClusHitAssoc[simClusIndex]:
                    for hitIndex in hitIndexArray:
                        thisHit = event.rechits_raw[hitIndex]
                        # print thisHit.energy
                        histDict["RecHits_layers_energy"].Fill(thisHit.layer, thisHit.energy)
                        histDict["RecHits_layers_pt"].Fill(thisHit.layer, thisHit.pt)
                        recHitTLV = ROOT.TLorentzVector()
                        recHitTLV.SetPtEtaPhiE(thisHit.pt, thisHit.eta, thisHit.phi, thisHit.energy)
                        recHitVectors["both"] += recHitTLV
                        if (thisHit.layer < 29):
                            recHitVectors["EE"] += recHitTLV
                        else:
                            recHitVectors["FH"] += recHitTLV
                logger.debug("SimCluster pt, E: {}, {} - RecHitVector pt, E: {}, {}".format(simCl.pt, simCl.energy, recHitVectors["both"].Pt(), recHitVectors["both"].E()))
                # relative pT cut to clean up misreconstructed particles
                if (recHitVectors["both"].Pt() < 0.8*simCl.pt):
                    continue
                histDict["RecHits_energy"].Fill(recHitVectors["both"].E())
                histDict["RecHits_eta"].Fill(recHitVectors["both"].Eta())
                histDict["RecHits_phi"].Fill(recHitVectors["both"].Phi())
                histDict["RecHits_pt"].Fill(recHitVectors["both"].Pt())
                histDict["SimVsRecHits_delta_eta"].Fill(simCl.eta-recHitVectors["both"].Eta())
                histDict["SimVsRecHits_delta_phi"].Fill(simCl.phi-recHitVectors["both"].Phi())
                # histDict["SimVsRecHits_delta_R"].Fill(HGCalHelpers.deltaR(simCl, recHitVector))
                for detect in detectors:
                    histDict["SimVsRecHits_delta_energy_%s" % detect].Fill(simCl.energy-recHitVectors[detect].E())
                    histDict["SimVsRecHits_delta_pt_%s" % detect].Fill(simCl.pt-recHitVectors[detect].Pt())
                    histDict["SimVsRecHits_deltaover_energy_%s" % detect].Fill((simCl.energy-recHitVectors[detect].E())/simCl.energy)
                    histDict["SimVsRecHits_deltaover_pt_%s" % detect].Fill((simCl.pt-recHitVectors[detect].Pt())/simCl.pt)
                    histDict["SimVsRecHits_frac_energy_%s" % detect].Fill(recHitVectors[detect].E()/simCl.energy)
                    histDict["SimVsRecHits_frac_pt_%s" % detect].Fill(recHitVectors[detect].Pt()/simCl.pt)


        if (nEvents > 0 and currentEvent >= nEvents):
            break
    HGCalHelpers.saveHistograms(histDict, canvas, outDir, imgType, logScale=False)


if __name__ == '__main__':
    main()
