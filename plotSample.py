"""plot individual sample by looping over all histograms."""
import ROOT
import HGCalHelpers


def main():
    """plot individual sample by looping over all histograms."""
    samples2Run = ["chargedPions_nPart1_E2_pre15_5k", "chargedPions_nPart1_E5_pre15_5k",
                   "chargedPions_nPart1_E10_pre15_5k", "chargedPions_nPart1_E20_pre15_5k",
                   "chargedPions_nPart1_E40_pre15_5k", "chargedPions_nPart1_E80_pre15_5k",
                   "chargedPions_nPart1_E160_pre15_5k", "chargedPions_nPart1_E320_pre15_5k"
                   ]

    outDir = '/afs/cern.ch/work/c/clange/HGCal/output'
    inDir = '/afs/cern.ch/work/c/clange/HGCal/output'
    imgType = 'pdf'
    ROOT.gROOT.SetBatch(True)
    canvas = ROOT.TCanvas("plot", "plot", 500, 500)
    canvas.SetRightMargin(0.12)
    canvas.SetLeftMargin(0.12)

    for sample in samples2Run:
        selectedEvents = 0
        histDict = {}
        print "Plotting", sample
        sampleOutDir = outDir + '/' + sample
        HGCalHelpers.createOutputDir(sampleOutDir)
        file0 = ROOT.TFile(inDir + "/" + sample + ".root")
        file0.cd()
        selectedEventsHist = file0.Get("selectedEvents")
        selectedEvents = selectedEventsHist.GetBinContent(1)
        dirList = ROOT.gDirectory.GetListOfKeys()
        keysToFind = ["RecHitsClus_layers", "fracEvents_", "layers_N", "radius_frac_energy", "radius_events_eFrac", "RecHitsClusVsRecHits_radius_energy", "RecHitsClusVsRecHits_total_energy", "layers_delta_d2_eWeight", "RecHits_layers"]
        for k1 in dirList:
            hist = k1.ReadObj()
            histName = hist.GetName()
            for findKey in keysToFind:
                if (histName.find(findKey) >= 0):
                    # divide by number of selected SimClusters/events
                    hist.Scale(1. / selectedEvents)
                    # make sure to only divide once
                    break
            histDict[histName] = hist
        HGCalHelpers.saveHistograms(
            histDict, canvas, sampleOutDir, imgType, plotsOnly=True)


if __name__ == '__main__':
    main()
