import ROOT
# import logging

def main():
    samples2plot = ["chargedPions_nPart1_Pt5_pre15",
                    "chargedPions_nPart1_Pt10_pre15",
                    "chargedPions_nPart1_Pt20_pre15",
                    "chargedPions_nPart1_Pt35_pre15"]

    colourBlindColours = {}
    colourBlindColours[0] = ROOT.TColor(10000, 0, 0.4470588235, 0.6980392157)
    colourBlindColours[1] = ROOT.TColor(10001, 0.337254902, 0.7058823529, 0.9137254902)
    colourBlindColours[2] = ROOT.TColor(10002, 0.8, 0.4745098039, 0.6549019608)
    colourBlindColours[3] = ROOT.TColor(10003, 0, 0.6196078431, 0.4509803922)
    colourBlindColours[4] = ROOT.TColor(10004, 0.8352941176, 0.368627451, 0)

    sampleColours = {}
    sampleColours["chargedPions_nPart1_Pt5_pre15"] = 10000
    sampleColours["chargedPions_nPart1_Pt10_pre15"] = 10001
    sampleColours["chargedPions_nPart1_Pt20_pre15"] = 10002
    sampleColours["chargedPions_nPart1_Pt35_pre15"] = 10003

    sampleMarkerStyle = {}
    sampleMarkerStyle["chargedPions_nPart1_Pt5_pre15"] = 20
    sampleMarkerStyle["chargedPions_nPart1_Pt10_pre15"] = 21
    sampleMarkerStyle["chargedPions_nPart1_Pt20_pre15"] = 22
    sampleMarkerStyle["chargedPions_nPart1_Pt35_pre15"] = 23

    sampleLabels = {}
    sampleLabels["chargedPions_nPart1_Pt5_pre15"] = "p_{T} = 5 GeV"
    sampleLabels["chargedPions_nPart1_Pt10_pre15"] = "p_{T} = 10 GeV"
    sampleLabels["chargedPions_nPart1_Pt20_pre15"] = "p_{T} = 20 GeV"
    sampleLabels["chargedPions_nPart1_Pt35_pre15"] = "p_{T} = 35 GeV"

    histList2D = ["RecHitsClus_layers_energy_cumulative",
                  "RecHitsClus_layers_energy_relative",
                  "RecHits_layers_energy"]

    histList1D = ["RecHitsClus_layers_energy_1D_cumulative",
                  "RecHitsClus_layers_energy_1D_relative",
                  "RecHits_layers_energy_1D",
                  "RecHits_layers_nHits"]

    histYTitles = {}
    histYTitles["RecHitsClus_layers_energy_cumulative"] = "cumulative relative RecHits/SimCluster energy deposits"
    histYTitles["RecHitsClus_layers_energy_relative"] = "relative RecHits/SimCluster energy deposits per layer"
    histYTitles["RecHits_layers_energy"] = "RecHits energy deposits [a.u.]"
    histYTitles["RecHitsClus_layers_energy_1D_cumulative"] = "cumulative relative RecHits/SimCluster energy deposits"
    histYTitles["RecHitsClus_layers_energy_1D_relative"] = "relative RecHits/SimCluster energy deposits per layer"
    histYTitles["RecHits_layers_energy_1D"] = "RecHits energy deposits [a.u.]"
    histYTitles["RecHits_layers_nHits"] = "RecHits number of hits [a.u.]"

    fileDict = {}
    for sampleName in samples2plot:
        fileDict[sampleName] = ROOT.TFile.Open("%s.root" % sampleName)

    firstPlot = True
    myLegend = ROOT.TLegend(0.16,0.75,0.48,0.95)
    myLegend.SetBorderSize(0)
    myLegend.SetFillStyle(0)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.02)

    c = ROOT.TCanvas("c", "c", 500, 500)
    for hist in histList2D:
        histDict2D = {}
        histStackpfX = ROOT.THStack("%s_stack_pfX" % hist, "")
        histStackpjX = ROOT.THStack("%s_stack_pjX" % hist, "")
        maxYprojX = 0
        for sampleName in samples2plot:
            histDict2D[hist] = fileDict[sampleName].Get(hist)
            profileX = histDict2D[hist].ProfileX("%s_%s_pfX" % (sampleName, hist))
            profileX.SetMarkerColor(sampleColours[sampleName])
            profileX.SetLineColor(sampleColours[sampleName])
            profileX.SetMarkerStyle(sampleMarkerStyle[sampleName])
            projX = histDict2D[hist].ProjectionX("%s_%s_pjX" % (sampleName, hist))
            projX.SetMarkerColor(sampleColours[sampleName])
            projX.SetLineColor(sampleColours[sampleName])
            projX.SetMarkerStyle(sampleMarkerStyle[sampleName])
            if projX.GetMaximum() > maxYprojX:
                maxYprojX = projX.GetMaximum()
            if firstPlot:
                myLegend.AddEntry(profileX, sampleLabels[sampleName], "p")
            for bin in range(0, profileX.GetNbinsX()+1):
                profileX.SetBinError(bin, 0.)
            histStackpfX.Add(profileX)
            histStackpjX.Add(projX)
        histStackpfX.Draw("nostack,hist,p")
        histStackpfX.GetXaxis().SetTitle("layer")
        histStackpfX.GetYaxis().SetTitle(histYTitles[hist])
        histStackpfX.GetYaxis().SetTitleOffset(histStackpfX.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        firstPlot = False
        c.SaveAs("%s_pfX.pdf" % hist)
        histStackpjX.Draw("nostack,hist,p")
        histStackpjX.GetXaxis().SetTitle("layer")
        histStackpjX.SetMaximum(maxYprojX*1.2)
        histStackpjX.GetYaxis().SetTitle(histYTitles[hist])
        histStackpjX.GetYaxis().SetTitleOffset(histStackpjX.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        c.SaveAs("%s_pjX.pdf" % hist)

    for hist in histList1D:
        histDict1D = {}
        histStack = ROOT.THStack("%s_stack" % hist, "")
        maxY = 0
        for sampleName in samples2plot:
            histDict1D[hist] = fileDict[sampleName].Get(hist).Clone("%s%s" %(hist,sampleName))
            histDict1D[hist].SetMarkerColor(sampleColours[sampleName])
            histDict1D[hist].SetLineColor(sampleColours[sampleName])
            histDict1D[hist].SetMarkerStyle(sampleMarkerStyle[sampleName])
            if histDict1D[hist].GetMaximum() > maxY:
                maxY = histDict1D[hist].GetMaximum()
            if firstPlot:
                myLegend.AddEntry(histDict1D[hist], sampleLabels[sampleName], "p")
            for bin in range(0, histDict1D[hist].GetNbinsX()+1):
                histDict1D[hist].SetBinError(bin, 0.)
            histStack.Add(histDict1D[hist])
        histStack.Draw("nostack,hist,p")
        histStack.SetMaximum(maxY*1.2)
        histStack.GetXaxis().SetTitle("layer")
        histStack.GetYaxis().SetTitle(histYTitles[hist])
        histStack.GetYaxis().SetTitleOffset(histStack.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        firstPlot = False
        c.SaveAs("%s.pdf" % hist)



if __name__ == '__main__':
    main()
