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
        for k1 in dirList:
            hist = k1.ReadObj()
            histName = hist.GetName()
            if ((histName.find("RecHitsClus_layers") >= 0) or (histName.find("fracEvents_") >= 0) or (histName.find("layers_N") >= 0) or (histName.find("radius_frac_energy") >= 0) or (histName.find("radius_events_eFrac") >= 0)):
                # divide by number of selected SimClusters/events
                hist.Scale(1. / selectedEvents)
            histDict[histName] = hist
        HGCalHelpers.saveHistograms(
            histDict, canvas, sampleOutDir, imgType, plotsOnly=True)


if __name__ == '__main__':
    main()
