#!/usr/bin/env python
import ROOT
from SampleHelper import SampleManager
from NtupleDataFormat import HGCalNtuple
import HGCalHelpers


# Make general clustering validation plots for several samples
rangeFolders = ['minus_z', 'plus_z', 'minus_z_eta', 'plus_z_eta']
thicknessDict = {}
thicknessDict[100] = "100_"
thicknessDict[200] = "200_"
thicknessDict[300] = "300_"
thicknessDict[340282346638528859811704183484516925440] = "other_"
thicknessDict[0] = "other_"

maxEvents = 10


def getHists():
    """function to book all the histograms and return as dictionary."""
    histDict = {}
    histDict["selectedEvents"] = ROOT.TH1F("selectedEvents", "selectedEvents", 1, 0.5, 1.5)
    clusters = ["MultiClus", "LayerClus", "LayerClus_100", "LayerClus_200", "LayerClus_300", "LayerClus_other"]

    maxAxisRanges = {}
    maxAxisRanges["MultiClus"] = 1000
    maxAxisRanges["LayerClus"] = 5000
    maxAxisRanges["LayerClus_100"] = 5000
    maxAxisRanges["LayerClus_200"] = 5000
    maxAxisRanges["LayerClus_300"] = 5000
    maxAxisRanges["LayerClus_other"] = 5000
    nhitCoreMax = 200
    nhitAllMax = 500
    for currentRange in rangeFolders:
        for cluster in clusters:
            if "MultiClus" in cluster:
                histDict['{}_{}_nclus'.format(cluster, currentRange)] = ROOT.TH1F('{}_{}_nclus'.format(cluster, currentRange), '{}_{}_nclus;N 2d clusters'.format(cluster, currentRange), 26, -.5, 25.5)
            histDict['{}_{}_mult'.format(cluster, currentRange)] = ROOT.TH1F('{}_{}_mult'.format(cluster, currentRange), '{}_{}_mult'.format(cluster, currentRange), maxAxisRanges[cluster], -0.5, maxAxisRanges[cluster]-0.5)
            histDict['{}_{}_eta'.format(cluster, currentRange)] = ROOT.TH1F('{}_{}_eta'.format(cluster, currentRange), '{}_{}_eta'.format(cluster, currentRange), 34, 1.4, 3.2)
            histDict['{}_{}_pt'.format(cluster, currentRange)] = ROOT.TH1F('{}_{}_pt'.format(cluster, currentRange), '{}_{}_pt'.format(cluster, currentRange), 100, 0, 5)
            if "LayerClus" in cluster:
                histDict['{}_{}_nhitCore'.format(cluster, currentRange)] = ROOT.TH1F('{}_{}_nhitCore'.format(cluster, currentRange), '{}_{}_nhitCore'.format(cluster, currentRange), 100, 0, nhitCoreMax)
                histDict['{}_{}_nhitAll'.format(cluster, currentRange)] = ROOT.TH1F('{}_{}_nhitAll'.format(cluster, currentRange), '{}_{}_nhitAll'.format(cluster, currentRange), 100, 0, nhitAllMax)
                for layer in range(1, 53):
                    histDict['{0}_{1}_{2:0>2}_mult'.format(cluster, currentRange, layer)] = ROOT.TH1F('{0}_{1}_{2:0>2}_mult'.format(cluster, currentRange, layer), '{0}_{1}_{2:0>2}_mult'.format(cluster, currentRange, layer), maxAxisRanges[cluster]/5, 0-0.5, maxAxisRanges[cluster]/5-0.5)
                    histDict['{0}_{1}_{2:0>2}_eta'.format(cluster, currentRange, layer)] = ROOT.TH1F('{0}_{1}_{2:0>2}_eta'.format(cluster, currentRange, layer), '{0}_{1}_{2:0>2}_eta'.format(cluster, currentRange, layer), 34, 1.4, 3.2)
                    histDict['{0}_{1}_{2:0>2}_pt'.format(cluster, currentRange, layer)] = ROOT.TH1F('{0}_{1}_{2:0>2}_pt'.format(cluster, currentRange, layer), '{0}_{1}_{2:0>2}_pt'.format(cluster, currentRange, layer), 100, 0, 5)
                    histDict['{0}_{1}_{2:0>2}_nhitCore'.format(cluster, currentRange, layer)] = ROOT.TH1F('{0}_{1}_{2:0>2}_nhitCore'.format(cluster, currentRange, layer), '{0}_{1}_{2:0>2}_nhitCore'.format(cluster, currentRange, layer), 100, 0, nhitCoreMax)
                    histDict['{0}_{1}_{2:0>2}_nhitAll'.format(cluster, currentRange, layer)] = ROOT.TH1F('{0}_{1}_{2:0>2}_nhitAll'.format(cluster, currentRange, layer), '{0}_{1}_{2:0>2}_nhitAll'.format(cluster, currentRange, layer), 100, 0, nhitAllMax)

    # event display like plot
    for i in range(1, maxEvents+1):
        histDict["eventDisplay_{}".format(i)] = ROOT.TH3F("eventDisplay_{}".format(i), "eventDisplay;layer;#phi;#eta", 52, 1, 52, 50, -3.1415, 3.1415, 50, -3.2, 3.2)

    return histDict


def main():

    samples2Run = ['FlatRandomPtGunProducer_SinglePion_35GeV_20170523',
                   'FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523',
                   'RelValTTbar_14TeV_CMSSW_9_1_0_pre3-PU25ns_91X_upgrade2023_realistic_v1_D13PU200-v2_GEN-SIM-RECO']

    sampleManager = SampleManager()
    for sampleName in samples2Run:
        sample = sampleManager.getSample(sampleName)
        print "Sample {} has {} files".format(sampleName, len(sample.getFiles()))

        outDir = sampleName
        HGCalHelpers.createOutputDir(outDir)

        rootOnly = False
        imgType = 'png'
        canvas = None
        if not rootOnly:
            canvas = ROOT.TCanvas(outDir, outDir, 600, 600)

        outFile = ROOT.TFile.Open("{}.root".format(sampleName), "RECREATE")

        histDict = getHists()

        currentEvent = 0

        for inFile in sample.getFiles():
            if maxEvents > 0:
                if currentEvent > maxEvents:
                    break

            print inFile
            ntuple = HGCalNtuple(inFile)
            for event in ntuple:
                currentEvent += 1
                if (currentEvent % 10 == 0):
                    print "Event", currentEvent, "of", maxEvents
                if maxEvents > 0:
                    if currentEvent > maxEvents:
                        break

                # # multi clusters
                # multiClusters = event.multiClusters()
                # multiClusterCounter = {}
                # for currentRange in rangeFolders:
                #     multiClusterCounter[currentRange] = 0
                # for multiCluster in multiClusters:
                #     if (multiCluster.z() < 0):
                #         multiClusterCounter['minus_z'] += 1
                #         histDict['{}_{}_eta'.format('MultiClus', currentRange)].Fill(abs(multiCluster.eta()))
                #         histDict['{}_{}_pt'.format('MultiClus', currentRange)].Fill(multiCluster.pt())
                #         histDict['{}_{}_nclus'.format('MultiClus', currentRange)].Fill(len(multiCluster.cluster2d()))
                #         if (1.6 < abs(multiCluster.eta()) < 2.6):
                #             multiClusterCounter['minus_z_eta'] += 1
                #             histDict['{}_{}_eta'.format('MultiClus', currentRange)].Fill(abs(multiCluster.eta()))
                #             histDict['{}_{}_pt'.format('MultiClus', currentRange)].Fill(multiCluster.pt())
                #             histDict['{}_{}_nclus'.format('MultiClus', currentRange)].Fill(len(multiCluster.cluster2d()))
                #         # print multiCluster.pt()
                #     else:
                #         multiClusterCounter['plus_z'] += 1
                #         histDict['{}_{}_eta'.format('MultiClus', currentRange)].Fill(abs(multiCluster.eta()))
                #         histDict['{}_{}_pt'.format('MultiClus', currentRange)].Fill(multiCluster.pt())
                #         histDict['{}_{}_nclus'.format('MultiClus', currentRange)].Fill(len(multiCluster.cluster2d()))
                #         if (1.6 < abs(multiCluster.eta()) < 2.6):
                #             multiClusterCounter['plus_z_eta'] += 1
                #             histDict['{}_{}_eta'.format('MultiClus', currentRange)].Fill(abs(multiCluster.eta()))
                #             histDict['{}_{}_pt'.format('MultiClus', currentRange)].Fill(multiCluster.pt())
                #             histDict['{}_{}_nclus'.format('MultiClus', currentRange)].Fill(len(multiCluster.cluster2d()))

                # layer clusters
                layerClusters = event.layerClusters()
                recHits = event.recHits()
                layerClusterCounter = {}
                for currentRange in rangeFolders:
                    layerClusterCounter[currentRange] = []
                    layerClusterCounter["100_"+currentRange] = []
                    layerClusterCounter["200_"+currentRange] = []
                    layerClusterCounter["300_"+currentRange] = []
                    layerClusterCounter["other_"+currentRange] = []
                    for layer in range(0, 53):
                        # use index zero as overall counter
                        layerClusterCounter[currentRange].append(0)
                        layerClusterCounter["100_"+currentRange].append(0)
                        layerClusterCounter["200_"+currentRange].append(0)
                        layerClusterCounter["300_"+currentRange].append(0)
                        layerClusterCounter["other_"+currentRange].append(0)
                for layerCluster in layerClusters:
                    # if recHits[layerCluster.rechitSeed()].flags() != 0:
                    #     print recHits[layerCluster.rechitSeed()].flags(), layerCluster.pt(), len(layerCluster.rechits())
                    histDict["eventDisplay_{}".format(currentEvent)].Fill(layerCluster.layer(), layerCluster.phi(), layerCluster.eta(), layerCluster.pt())
                    # print layerCluster.rechitSeed()
                    # seedHitThickness = int(recHits[layerCluster.rechitSeed()].thickness())
                    # # print "seedHitThickness", seedHitThickness
                    # histStringsForFilling = []
                    # if (layerCluster.z() < 0):
                    #     histStringsForFilling.append('minus_z')
                    #     histStringsForFilling.append(thicknessDict[seedHitThickness] + 'minus_z')
                    #     if (1.6 < abs(layerCluster.eta()) < 2.6):
                    #         histStringsForFilling.append('minus_z_eta')
                    #         histStringsForFilling.append(thicknessDict[seedHitThickness] + 'minus_z_eta')
                    # else:
                    #     histStringsForFilling.append('plus_z')
                    #     histStringsForFilling.append(thicknessDict[seedHitThickness] + 'plus_z')
                    #     if (1.6 < abs(layerCluster.eta()) < 2.6):
                    #         histStringsForFilling.append('plus_z_eta')
                    #         histStringsForFilling.append(thicknessDict[seedHitThickness] + 'plus_z_eta')
                    # for histString in histStringsForFilling:
                    #     layerClusterCounter[histString][0] += 1
                    #     layerClusterCounter[histString][layerCluster.layer()] += 1
                    #     histDict['{}_{}_eta'.format('LayerClus', histString)].Fill(abs(layerCluster.eta()))
                    #     histDict['{0}_{1}_{2:0>2}_eta'.format('LayerClus', histString, layerCluster.layer())].Fill(abs(layerCluster.eta()))
                    #     histDict['{}_{}_pt'.format('LayerClus', histString)].Fill(layerCluster.pt())
                    #     histDict['{0}_{1}_{2:0>2}_pt'.format('LayerClus', histString, layerCluster.layer())].Fill(layerCluster.pt())
                    #     histDict['{}_{}_nhitCore'.format('LayerClus', histString)].Fill(layerCluster.nhitCore())
                    #     histDict['{0}_{1}_{2:0>2}_nhitCore'.format('LayerClus', histString, layerCluster.layer())].Fill(layerCluster.nhitCore())
                    #     histDict['{}_{}_nhitAll'.format('LayerClus', histString)].Fill(layerCluster.nhitAll())
                    #     histDict['{0}_{1}_{2:0>2}_nhitAll'.format('LayerClus', histString, layerCluster.layer())].Fill(layerCluster.nhitAll())

                # fill counting histograms
                # for currentRange in rangeFolders:
                #     for layer in range(0, 53):
                #         # if (currentEvent == 1):
                #         #     print currentRange, layer, layerClusterCounter[currentRange][layer]
                #         if (layer == 0):
                #             histDict['{}_{}_mult'.format('LayerClus', currentRange)].Fill(layerClusterCounter[currentRange][layer])
                #             histDict['{}_{}_mult'.format('MultiClus', currentRange)].Fill(multiClusterCounter[currentRange])
                #         else:
                #             histDict['{0}_{1}_{2:0>2}_mult'.format('LayerClus', currentRange, layer)].Fill(layerClusterCounter[currentRange][layer])

        outFile.cd()
        # for key, value in histDict.items():
        #     if value.GetEntries() != 0:
        #         value.Scale(1./maxEvents)
        HGCalHelpers.saveHistograms(histDict, canvas, outDir, imgType, logScale=False, rootOnly=rootOnly)
        outFile.Write()
        outFile.Close()


if __name__ == '__main__':
    main()
