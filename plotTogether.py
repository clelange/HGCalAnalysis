import ROOT
# import logging

def main():
    samples2plot = ["chargedPions_nPart1_Pt2_pre15_5k",
                    "chargedPions_nPart1_Pt5_pre15_5k",
                    # "chargedPions_nPart1_Pt10_pre15_5k",
                    # "chargedPions_nPart1_Pt20_pre15_5k",
                    "chargedPions_nPart1_Pt35_pre15_5k",
                    "chargedPions_nPart1_Pt200_pre15_5k"
                    ]

    colourBlindColours = {}
    colourBlindColours[0] = ROOT.TColor(10000, 0, 0.4470588235, 0.6980392157)
    colourBlindColours[1] = ROOT.TColor(10001, 0.337254902, 0.7058823529, 0.9137254902)
    colourBlindColours[2] = ROOT.TColor(10002, 0.8, 0.4745098039, 0.6549019608)
    colourBlindColours[3] = ROOT.TColor(10003, 0, 0.6196078431, 0.4509803922)
    colourBlindColours[4] = ROOT.TColor(10004, 0.8352941176, 0.368627451, 0)

    sampleColours = {}
    sampleColours["chargedPions_nPart1_Pt2_pre15_5k"] = 10000
    sampleColours["chargedPions_nPart1_Pt5_pre15_5k"] = 10001
    sampleColours["chargedPions_nPart1_Pt35_pre15_5k"] = 10002
    sampleColours["chargedPions_nPart1_Pt200_pre15_5k"] = 10003

    sampleMarkerStyle = {}
    sampleMarkerStyle["chargedPions_nPart1_Pt2_pre15_5k"] = 20
    sampleMarkerStyle["chargedPions_nPart1_Pt5_pre15_5k"] = 21
    sampleMarkerStyle["chargedPions_nPart1_Pt35_pre15_5k"] = 22
    sampleMarkerStyle["chargedPions_nPart1_Pt200_pre15_5k"] = 23

    sampleLabels = {}
    sampleLabels["chargedPions_nPart1_Pt2_pre15_5k"] = "p_{T} = 2 GeV"
    sampleLabels["chargedPions_nPart1_Pt5_pre15_5k"] = "p_{T} = 5 GeV"
    sampleLabels["chargedPions_nPart1_Pt35_pre15_5k"] = "p_{T} = 35 GeV"
    sampleLabels["chargedPions_nPart1_Pt200_pre15_5k"] = "p_{T} = 200 GeV"

    histList2D = ["RecHitsClus_layers_energy_cumulative",
                  "RecHitsClus_layers_energy_relative",
                  "RecHits_layers_energy",
                  "SimVsRecHits_frac_energy_EE_FH+BH_fullEta",
                  "RecHits_layers_delta_x",
                  "RecHits_layers_delta_y",
                  "RecHits_layers_delta_d",
                  ]

    histList1D = ["RecHitsClus_layers_energy_1D_cumulative",
                  "RecHitsClus_layers_energy_1D_relative",
                  "RecHitsClus_layers_energy_plain_1D_cumulative",
                  "RecHitsClus_layers_energy_plain_1D_relative",
                  "RecHits_layers_energy_1D",
                  "RecHits_layers_nHits",
                  "SimVsRecHits_frac_energy_EE_fullEta",
                  "SimVsRecHits_frac_energy_BH_fullEta",
                  "SimVsRecHits_frac_energy_BH_fullEta_fracEvents_1D",
                  "SimVsRecHits_frac_energy_BH_fullEta_fracEvents_1D_cumulative",
                  "SimClus_pt",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_all",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_EE",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_BH",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH+BH",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_all",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_EE",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_BH",
                  "RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH+BH",
                  "RecHitsClusVsRecHits_radius_frac_energy_SimCluster_all",
                  "RecHitsClusVsRecHits_radius_frac_energy_SimCluster_EE",
                  "RecHitsClusVsRecHits_radius_frac_energy_SimCluster_FH",
                  "RecHitsClusVsRecHits_radius_frac_energy_SimCluster_BH",
                  "RecHitsClusVsRecHits_radius_frac_energy_SimCluster_FH+BH",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_all",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_EE",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_BH",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH+BH",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_all",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_EE",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_BH",
                  "SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH+BH",
                  "SimVsRecHits_radius_frac_energy_SimCluster_all",
                  "SimVsRecHits_radius_frac_energy_SimCluster_EE",
                  "SimVsRecHits_radius_frac_energy_SimCluster_FH",
                  "SimVsRecHits_radius_frac_energy_SimCluster_BH",
                  "SimVsRecHits_radius_frac_energy_SimCluster_FH+BH",
                  ]

    histYTitles = {}
    histYTitles["RecHitsClus_layers_energy_cumulative"] = "cumulative sum of RecHits energy deposits/SimCluster energy"
    histYTitles["RecHitsClus_layers_energy_relative"] = "relative RecHits/SimCluster energy deposits per layer"
    histYTitles["RecHits_layers_energy"] = "number of RecHits energy deposits [a.u.]"
    histYTitles["SimVsRecHits_frac_energy_EE_FH+BH_fullEta"] = "events"
    histYTitles["RecHitsClus_layers_energy_1D_cumulative"] = "cumulative relative RecHits/SimCluster energy deposits"
    histYTitles["RecHitsClus_layers_energy_1D_relative"] = "relative RecHits/SimCluster energy deposits per layer"
    histYTitles["RecHitsClus_layers_energy_plain_1D_cumulative"] = "cumulative relative RecHits/SimCluster energy deposits"
    histYTitles["RecHitsClus_layers_energy_plain_1D_relative"] = "relative RecHits/SimCluster energy deposits per layer"
    histYTitles["RecHits_layers_energy_1D"] = "RecHits energy deposits [a.u.]"
    histYTitles["RecHits_layers_nHits"] = "RecHits number of hits [a.u.]"
    histYTitles["SimVsRecHits_frac_energy_EE_fullEta"] = "events"
    histYTitles["SimVsRecHits_frac_energy_BH_fullEta"] = "events"
    histYTitles["SimVsRecHits_frac_energy_BH_fullEta_fracEvents_1D"] = "fraction of events"
    histYTitles["SimVsRecHits_frac_energy_BH_fullEta_fracEvents_1D_cumulative"] = "cumulative fraction of events"
    histYTitles["SimClus_pt"] = "entries"
    histYTitles["RecHits_layers_delta_x"] = "#Delta x [mm]"
    histYTitles["RecHits_layers_delta_y"] = "#Delta y [mm]"
    histYTitles["RecHits_layers_delta_d"] = "#Delta d [mm]"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_all"] = "#sum RecHits(radius) / #sum RecHits(total) E for all subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_EE"] = "#sum RecHits(radius) / #sum RecHits(total) E for EE subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH"] = "#sum RecHits(radius) / #sum RecHits(total) E for FH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_BH"] = "#sum RecHits(radius) / #sum RecHits(total) E for BH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH+BH"] = "#sum RecHits(radius) / #sum RecHits(total) E for FH+BH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_all"] = "#sum RecHits(radius) / #sum RecHits(total) E for all subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_EE"] = "#sum RecHits(radius) / #sum RecHits(total) E for EE subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH"] = "#sum RecHits(radius) / #sum RecHits(total) E for FH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_BH"] = "#sum RecHits(radius) / #sum RecHits(total) E for BH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH+BH"] = "#sum RecHits(radius) / #sum RecHits(total) E for FH+BH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_all"] = "#sum RecHits(radius) / #sum RecHits(total) E for all subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_EE"] = "#sum RecHits(radius) / #sum RecHits(total) E for EE subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_FH"] = "#sum RecHits(radius) / #sum RecHits(total) E for FH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_BH"] = "#sum RecHits(radius) / #sum RecHits(total) E for BH subdetector(s)"
    histYTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_FH+BH"] = "#sum RecHits(radius) / #sum RecHits(total) E for FH+BH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_all"] = "#sum RecHits(radius) / SimCluster E for all subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_EE"] = "#sum RecHits(radius) / SimCluster E for EE subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH"] = "#sum RecHits(radius) / SimCluster E for FH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_BH"] = "#sum RecHits(radius) / SimCluster E for BH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH+BH"] = "#sum RecHits(radius) / SimCluster E for FH+BH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_all"] = "#sum RecHits(radius) / SimCluster E for all subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_EE"] = "#sum RecHits(radius) / SimCluster E for EE subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH"] = "#sum RecHits(radius) / SimCluster E for FH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_BH"] = "#sum RecHits(radius) / SimCluster E for BH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH+BH"] = "#sum RecHits(radius) / SimCluster E for FH+BH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_SimCluster_all"] = "#sum RecHits(radius) / SimCluster E for all subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_SimCluster_EE"] = "#sum RecHits(radius) / SimCluster E for EE subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_SimCluster_FH"] = "#sum RecHits(radius) / SimCluster E for FH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_SimCluster_BH"] = "#sum RecHits(radius) / SimCluster E for BH subdetector(s)"
    histYTitles["SimVsRecHits_radius_frac_energy_SimCluster_FH+BH"] = "#sum RecHits(radius) / SimCluster E for FH+BH subdetector(s)"

    for histName in histList2D+histList1D:
        if histName not in histYTitles:
            histYTitles[histName] = "entries"

    histXTitles = {}
    histXTitles["SimVsRecHits_frac_energy_EE_fullEta"] = "summed RecHits energy fraction EE"
    histXTitles["SimVsRecHits_frac_energy_BH_fullEta"] = "summed RecHits energy fraction BH"
    histXTitles["SimVsRecHits_frac_energy_BH_fullEta_fracEvents_1D"] = "summed RecHits energy fraction BH"
    histXTitles["SimVsRecHits_frac_energy_BH_fullEta_fracEvents_1D_cumulative"] = "summed RecHits energy fraction BH"
    histXTitles["SimVsRecHits_frac_energy_EE_FH+BH_fullEta"] = "summed RecHits energy fraction EE vs. FH+BH"
    histXTitles["SimClus_pt"] = "SimCluster p_{T} [GeV]"
    histXTitles["RecHits_layers_delta_x"] = "layer"
    histXTitles["RecHits_layers_delta_y"] = "layer"
    histXTitles["RecHits_layers_delta_d"] = "layer"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_all"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_EE"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_BH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH+BH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_all"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_EE"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_BH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH+BH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_all"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_EE"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_FH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_BH"] = "radius [mm]"
    histXTitles["RecHitsClusVsRecHits_radius_frac_energy_SimCluster_FH+BH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_all"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_EE"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_BH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusFirstLayer_FH+BH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_all"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_EE"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_BH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_RecHitsClusMaxLayer_FH+BH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_SimCluster_all"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_SimCluster_EE"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_SimCluster_FH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_SimCluster_BH"] = "radius [mm]"
    histXTitles["SimVsRecHits_radius_frac_energy_SimCluster_FH+BH"] = "radius [mm]"

    fileDict = {}
    for sampleName in samples2plot:
        fileDict[sampleName] = ROOT.TFile.Open("%s.root" % sampleName)

    firstPlot = True
    myLegend = ROOT.TLegend(0.16, 0.75, 0.48, 0.95)
    myLegend.SetBorderSize(0)
    myLegend.SetFillStyle(0)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.02)

    c = ROOT.TCanvas("c", "c", 500, 500)
    printMeanRMS = False
    for hist in histList2D:
        histDict2D = {}
        histStackpfX = ROOT.THStack("%s_stack_pfX" % hist, "")
        histStackpjX = ROOT.THStack("%s_stack_pjX" % hist, "")
        histStackpfY = ROOT.THStack("%s_stack_pfY" % hist, "")
        histStackpjY = ROOT.THStack("%s_stack_pjY" % hist, "")
        maxYprojX = 0
        maxYprojY = 0
        meanVecPfX = {}
        meanVecPjX = {}
        meanVecPfY = {}
        meanVecPjY = {}
        rmsVecPfX = {}
        rmsVecPjX = {}
        rmsVecPfY = {}
        rmsVecPjY = {}
        for sampleName in samples2plot:
            histDict2D[hist] = fileDict[sampleName].Get(hist)
            profileX = histDict2D[hist].ProfileX("%s_%s_pfX" % (sampleName, hist))
            profileX.SetMarkerColor(sampleColours[sampleName])
            profileX.SetLineColor(sampleColours[sampleName])
            profileX.SetMarkerStyle(sampleMarkerStyle[sampleName])
            meanVecPfX[sampleName] = profileX.GetMean()
            rmsVecPfX[sampleName] = profileX.GetRMS()
            projX = histDict2D[hist].ProjectionX("%s_%s_pjX" % (sampleName, hist))
            projX.SetMarkerColor(sampleColours[sampleName])
            projX.SetLineColor(sampleColours[sampleName])
            projX.SetMarkerStyle(sampleMarkerStyle[sampleName])
            meanVecPjX[sampleName] = projX.GetMean()
            rmsVecPjX[sampleName] = projX.GetRMS()
            # and y
            profileY = histDict2D[hist].ProfileY("%s_%s_pfY" % (sampleName, hist))
            profileY.SetMarkerColor(sampleColours[sampleName])
            profileY.SetLineColor(sampleColours[sampleName])
            profileY.SetMarkerStyle(sampleMarkerStyle[sampleName])
            meanVecPfY[sampleName] = profileY.GetMean()
            rmsVecPfY[sampleName] = profileY.GetRMS()
            projY = histDict2D[hist].ProjectionY("%s_%s_pjY" % (sampleName, hist))
            projY.SetMarkerColor(sampleColours[sampleName])
            projY.SetLineColor(sampleColours[sampleName])
            projY.SetMarkerStyle(sampleMarkerStyle[sampleName])
            meanVecPjY[sampleName] = projY.GetMean()
            rmsVecPjY[sampleName] = projY.GetRMS()
            if projX.GetMaximum() > maxYprojX:
                maxYprojX = projX.GetMaximum()
            projY.Scale(1./projY.Integral())
            if projY.GetMaximum() > maxYprojY:
                maxYprojY = projY.GetMaximum()
            if firstPlot:
                myLegend.AddEntry(profileX, sampleLabels[sampleName], "p")
            for bin in range(0, profileX.GetNbinsX()+1):
                profileX.SetBinError(bin, 0.)
            for bin in range(0, profileY.GetNbinsX()+1):
                profileY.SetBinError(bin, 0.)
            histStackpfX.Add(profileX)
            histStackpjX.Add(projX)
            histStackpfY.Add(profileY)
            histStackpjY.Add(projY)
        # draw profileX
        histStackpfX.Draw("nostack,hist,p")
        if (hist.find("layer") >= 0):
            printMeanRMS = False
            histStackpfX.GetXaxis().SetTitle("layer")
        else:
            printMeanRMS = True
            histStackpfX.GetXaxis().SetTitle(histXTitles[hist])
        histStackpfX.GetYaxis().SetTitle(histYTitles[hist])
        histStackpfX.GetYaxis().SetTitleOffset(histStackpfX.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        if printMeanRMS:
            tBox = ROOT.TLatex()
            tBox.SetTextSize(0.03)
            tBox.SetTextFont(42)
            for index in range(len(samples2plot)):
                tBox.DrawLatexNDC(0.5, 0.9-index*0.04, "mean: %5.2f, RMS: %5.2f" % (meanVecPfX[samples2plot[index]], rmsVecPfX[samples2plot[index]]))
        firstPlot = False
        c.SaveAs("%s_pfX.pdf" % hist)
        # draw projectionX
        meanVec = {}
        rmsVec = {}
        histStackpjX.Draw("nostack,hist,p")
        if (hist.find("layer") >= 0):
            histStackpjX.GetXaxis().SetTitle("layer")
        else:
            histStackpjX.GetXaxis().SetTitle(histXTitles[hist])
        histStackpjX.SetMaximum(maxYprojX*1.2)
        histStackpjX.GetYaxis().SetTitle(histYTitles[hist])
        histStackpjX.GetYaxis().SetTitleOffset(histStackpjX.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        if printMeanRMS:
            tBox = ROOT.TLatex()
            tBox.SetTextSize(0.03)
            tBox.SetTextFont(42)
            for index in range(len(samples2plot)):
                tBox.DrawLatexNDC(0.5, 0.9-index*0.04, "mean: %5.2f, RMS: %5.2f" % (meanVecPjX[samples2plot[index]], rmsVecPjX[samples2plot[index]]))
        c.SaveAs("%s_pjX.pdf" % hist)
        # and y
        histStackpfY.Draw("nostack,hist,p")
        histStackpfY.GetXaxis().SetTitle(histYTitles[hist])
        # histStackpfY.SetMaximum(maxYprojY*1.2)
        histStackpfY.GetYaxis().SetTitle("entries")
        histStackpfY.GetYaxis().SetTitleOffset(histStackpfY.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        if printMeanRMS:
            tBox = ROOT.TLatex()
            tBox.SetTextSize(0.03)
            tBox.SetTextFont(42)
            for index in range(len(samples2plot)):
                tBox.DrawLatexNDC(0.5, 0.9-index*0.04, "mean: %5.2f, RMS: %5.2f" % (meanVecPjX[samples2plot[index]], rmsVecPjX[samples2plot[index]]))
        c.SaveAs("%s_pfY.pdf" % hist)
        histStackpjY.Draw("nostack,hist,p")
        histStackpjY.GetXaxis().SetTitle(histYTitles[hist])
        histStackpjY.SetMaximum(maxYprojY*1.2)
        histStackpjY.GetYaxis().SetTitle("a.u. (unit area)")
        histStackpjY.GetYaxis().SetTitleOffset(histStackpjY.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        if printMeanRMS:
            tBox = ROOT.TLatex()
            tBox.SetTextSize(0.03)
            tBox.SetTextFont(42)
            for index in range(len(samples2plot)):
                tBox.DrawLatexNDC(0.5, 0.9-index*0.04, "mean: %5.2f, RMS: %5.2f" % (meanVecPjX[samples2plot[index]], rmsVecPjX[samples2plot[index]]))
        c.SaveAs("%s_pjY.pdf" % hist)

    for hist in histList1D:
        histDict1D = {}
        histStack = ROOT.THStack("%s_stack" % hist, "")
        maxY = 0
        meanVec = {}
        rmsVec = {}
        for sampleName in samples2plot:
            histDict1D[hist] = fileDict[sampleName].Get(hist).Clone("%s%s" %(hist,sampleName))
            histDict1D[hist].SetMarkerColor(sampleColours[sampleName])
            histDict1D[hist].SetLineColor(sampleColours[sampleName])
            histDict1D[hist].SetMarkerStyle(sampleMarkerStyle[sampleName])
            meanVec[sampleName] = histDict1D[hist].GetMean()
            rmsVec[sampleName] = histDict1D[hist].GetRMS()
            if histDict1D[hist].GetMaximum() > maxY:
                maxY = histDict1D[hist].GetMaximum()
            if firstPlot:
                myLegend.AddEntry(histDict1D[hist], sampleLabels[sampleName], "p")
            for bin in range(0, histDict1D[hist].GetNbinsX()+1):
                histDict1D[hist].SetBinError(bin, 0.)
            histStack.Add(histDict1D[hist])
        histStack.Draw("nostack,hist,p")
        histStack.SetMaximum(maxY*1.2)
        if (hist.find("layer") >= 0):
            printMeanRMS = False
            histStack.GetXaxis().SetTitle("layer")
        else:
            printMeanRMS = True
            histStack.GetXaxis().SetTitle(histXTitles[hist])
        histStack.GetYaxis().SetTitle(histYTitles[hist])
        histStack.GetYaxis().SetTitleOffset(histStack.GetYaxis().GetTitleOffset()*1.95)
        myLegend.Draw()
        if printMeanRMS:
            tBox = ROOT.TLatex()
            tBox.SetTextSize(0.03)
            tBox.SetTextFont(42)
            for index in range(len(samples2plot)):
                tBox.DrawLatexNDC(0.5, 0.9-index*0.04, "mean: %5.2f, RMS: %5.2f" % (meanVec[samples2plot[index]], rmsVec[samples2plot[index]]))
        firstPlot = False
        c.SaveAs("%s.pdf" % hist)



if __name__ == '__main__':
    main()
