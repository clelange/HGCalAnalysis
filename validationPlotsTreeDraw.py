import ROOT
from SampleHelper import SampleManager
import HGCalHelpers

nFiles = 3


plotVariablesMultiClus = ["Length$(multiclus_eta)", "multiclus_eta", "multiclus_pt", "multiclus_z"]
plotVariablesLayerClus = ["Length$(cluster2d_layer)", "cluster2d_eta", "cluster2d_pt", "cluster2d_nhitCore", "cluster2d_nhitAll"]
plotVariablesRecHits = ["Length$(rechit_eta)"]

cutsMultiClus = {"nocuts": "1.", "plus_z": "multiclus_eta > 0", "minus_z": "multiclus_eta < 0", "plus_z_etaUnit": "(1.6 < multiclus_eta) && (multiclus_eta <= 2.6)", "minus_z_etaUnit": "(-2.6 < multiclus_eta) && (multiclus_eta <= -1.6)"}
cutsLayerClus = {"nocuts": "1.", "plus_z": "cluster2d_eta > 0", "minus_z": "cluster2d_eta < 0", "plus_z_etaUnit": "(1.6 < cluster2d_eta) && (cluster2d_eta <= 2.6)", "minus_z_etaUnit": "(-2.6 < cluster2d_eta) && (cluster2d_eta <= -1.6)", "nocuts_thickness100": "rechit_thickness[cluster2d_rechitSeed] == 100", "nocuts_thickness200": "rechit_thickness[cluster2d_rechitSeed] == 200", "nocuts_thickness300": "rechit_thickness[cluster2d_rechitSeed] == 300", "nocuts_thicknessOther": "rechit_thickness[cluster2d_rechitSeed] > 300"}
cutsRecHits = {"nocuts": "1.", "plus_z": "rechit_eta > 0", "minus_z": "rechit_eta < 0", "plus_z_etaUnit": "(1.6 < rechit_eta) && (rechit_eta <= 2.6)", "minus_z_etaUnit": "(-2.6 < rechit_eta) && (rechit_eta <= -1.6)", "nocuts_thickness100": "rechit_thickness == 100", "nocuts_thickness200": "rechit_thickness == 200", "nocuts_thickness300": "rechit_thickness == 300", "nocuts_thicknessOther": "rechit_thickness > 300"}

def main():

    samples2Run = ['FlatRandomPtGunProducer_SinglePion_35GeV_20170523',
                   'FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523',
                   'RelValTTbar_14TeV_CMSSW_9_1_0_pre3-PU25ns_91X_upgrade2023_realistic_v1_D13PU200-v2_GEN-SIM-RECO'
                   ]

    sampleManager = SampleManager()
    for sampleName in samples2Run:
        sample = sampleManager.getSample(sampleName)
        print "Sample {} has {} files".format(sampleName, len(sample.fileList))

        outDir = sampleName
        HGCalHelpers.createOutputDir(outDir)
        for category in (cutsMultiClus.keys()+cutsLayerClus.keys()+cutsRecHits.keys()):
            HGCalHelpers.createOutputDir("{}/{}".format(outDir, category))

        imgType = 'png'
        canvas = ROOT.TCanvas(outDir, outDir, 600, 600)

        chain = sample.getChain(nFiles)
        nEvents = chain.GetEntries()
        print "Considering", nEvents, "events"

        for category, cut in cutsMultiClus.items():
            for variable in plotVariablesMultiClus:
                chain.Draw(variable, "({})*1./{}".format(cut, nEvents))
                outName = variable.replace("$", "_").replace("(", "").replace(")", "").replace("@", "")
                if variable.find("Length") >= 0:
                    outName = variable.split("(", 1)[1].rsplit(")")[0].rsplit("_", 1)[0] + "_N"
                canvas.SaveAs("{}/{}/{}.{}".format(outDir, category, outName, imgType))
                # TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp"); // 1D

        for category, cut in cutsLayerClus.items():
            for variable in plotVariablesLayerClus:
                outName = variable.replace("$", "_").replace("(", "").replace(")", "").replace("@", "")
                if variable.find("Length") >= 0:
                    outName = variable.split("(", 1)[1].rsplit(")")[0].rsplit("_", 1)[0] + "_N"

                chain.Draw(variable, "({})*1./{}".format(cut, nEvents))
                canvas.SaveAs("{}/{}/{}.{}".format(outDir, category, outName, imgType))
                for layer in range(1, 53):
                    if variable.find("Length") >= 0:
                        variableNew = variable.replace("Length$(cluster2d_layer)", "Sum$((cluster2d_layer == {}) && ({}))".format(layer, cut))
                        chain.Draw(variableNew, "1./{}".format(nEvents))
                        canvas.SaveAs("{0}/{1}/{2}_layer{3:0>2}.{4}".format(outDir, category, outName, layer, imgType))
                    else:
                        chain.Draw(variable, "(cluster2d_layer == {})&&({})*1./{}".format(layer, cut, nEvents))
                        canvas.SaveAs("{0}/{1}/{2}_layer{3:0>2}.{4}".format(outDir, category, variable.replace("$", "_").replace("(", "").replace(")", "").replace("@", ""), layer, imgType))

        for category, cut in cutsRecHits.items():
            for variable in plotVariablesRecHits:
                outName = variable.replace("$", "_").replace("(", "").replace(")", "").replace("@", "")
                if variable.find("Length") >= 0:
                    outName = variable.split("(", 1)[1].rsplit(")")[0].rsplit("_", 1)[0] + "_N"
                chain.Draw(variable, "({})*1./{}".format(cut, nEvents))
                canvas.SaveAs("{}/{}/{}.{}".format(outDir, category, outName, imgType))

        # outFile.cd()
        # for key, value in histDict.items():
        #     if value.GetEntries() != 0:
        #         value.Scale(1./maxEvents)
        # HGCalHelpers.saveHistograms(histDict, canvas, outDir, imgType, logScale=False, rootOnly=rootOnly)
        # outFile.Write()
        # outFile.Close()


if __name__ == '__main__':
    main()
