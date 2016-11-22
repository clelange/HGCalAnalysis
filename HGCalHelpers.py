"""Analysis helper tools for HGCal ntuples on EOS."""
import ROOT
import os
import math
import logging


class NullHandler(logging.Handler):
    """NullHandler for logging module."""

    def emit(self, record):
        """emit."""
        pass

logging.getLogger(__name__).addHandler(NullHandler())


def createOutputDir(outDir):
    """Create output directory if it does not yet exist."""
    if not os.path.exists(outDir):
        os.makedirs(outDir)


def saveHistograms(histDict, canvas, outDir, imgType, logScale=False, doFit=False):
    """Save all the histograms as ROOT file and image files, optionally fit Gaussian."""
    # also store histograms in ROOT file
    outFileName = "%s.root" % outDir
    outFile = ROOT.TFile(outFileName, "recreate")
    logString = ""
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.02)
    if logScale:
        canvas.SetLogy(True)
        logString = "_log"
    else:
        canvas.SetLogy(False)
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F):
            if item.GetEntries() == 0:
                continue
        # write histogram to file
        if type(item) == ROOT.TH2F or type(item) == ROOT.TH1F:
            item.Sumw2()
        item.Write()
        if type(item) == ROOT.TH2F:
            ROOT.gStyle.SetOptStat(0)
            item.Draw("colz")
            item.GetYaxis().SetTitleOffset(1.5)
        else:
            ROOT.gStyle.SetOptStat("mr")
            if type(item) == ROOT.TH1F:
                item.Draw("hist")
                item.GetYaxis().SetTitleOffset(1.5)
                canvas.SetGrid()
                # move the stats box
                ROOT.gPad.Update()
                ps = item.FindObject("stats")
                ps.SetX1NDC(0.15)
                ps.SetX2NDC(0.35)
                canvas.Modified()
                canvas.Update()
            else:
                item.Draw()
            if (doFit):
                if key.find("delta") >= 0 and key.find("delta_R") < 0 and key.find("deltaover") < 0:
                    ROOT.gStyle.SetOptFit(1)
                    item.Fit("gaus")
        if (item.GetYaxis().GetTitle() == ""):
            item.GetYaxis().SetTitle("a.u.")
        canvas.SaveAs("{}/{}{}.{}".format(outDir, key, logString, imgType))
        if type(item) == ROOT.TH2F:
            ROOT.gStyle.SetOptStat("mr")
            pjX = item.ProjectionX("pjX")
            pjX.Draw()
            canvas.SaveAs("{}/{}{}_projectionX.{}".format(outDir, key, logString, imgType))
            pjX.Delete()
            pjY = item.ProjectionY("pjY")
            pjY.Draw()
            canvas.SaveAs("{}/{}{}_projectionY.{}".format(outDir, key, logString, imgType))
            pjY.Delete()
            pfX = item.ProfileX("pfX")
            pfX.Draw()
            canvas.SaveAs("{}/{}{}_profileX.{}".format(outDir, key, logString, imgType))
            pfX.Delete()
            pfY = item.ProfileY("pfY")
            pfY.Draw()
            canvas.SaveAs("{}/{}{}_profileY.{}".format(outDir, key, logString, imgType))
            pfY.Delete()
    outFile.Write()
    outFile.Close()


def getGenParticles(event, histDict, dvzCut):
    """select GenParticles based on dvzCut and save also in histogram."""
    vGenParticleTLV = []
    for particle in event.particles:
        nonConverted = False
        if abs(particle.dvz) > dvzCut:
            nonConverted = True
            logging.debug("gen particle: {}, {}, {}".format(particle.pt, particle.eta, particle.phi))
        if not nonConverted:
            logging.debug("converted photon: {}".format(particle.dvz))
            continue
        if (abs(particle.eta) < 1.6) or (abs(particle.eta) > 2.8):
            logging.debug("photon outside detector coverage, eta: {}".format(particle.eta))
            continue
        if "GenPart_energy" in histDict:
            histDict["GenPart_energy"].Fill(particle.energy)
            histDict["GenPart_pt"].Fill(particle.pt)
            histDict["GenPart_eta"].Fill(particle.eta)
            histDict["GenPart_phi"].Fill(particle.phi)
            histDict["GenPart_dvz"].Fill(particle.dvz)
        particleTLV = ROOT.TLorentzVector()
        particleTLV.SetPtEtaPhiE(
            particle.pt, particle.eta, particle.phi, particle.energy)
        vGenParticleTLV.append(particleTLV)
    return vGenParticleTLV


def getMultiClusters(clusters, histDict, prefix, dRCut, energyCut, matchesGen):
    """get MultiClusters with option of applying cuts."""
    vMulticlusterTLV = []
    usedCluster = [False] * len(clusters)
    for j, simCl1 in enumerate((clusters)):
        logging.debug(usedCluster)
        logging.debug("j: {} {} {} {}".format(j, simCl1.pt, simCl1.eta, simCl1.phi))
        if (simCl1.energy < energyCut) or (usedCluster[j]) or not matchesGen[j]:
            logging.debug("skip")
            continue
        histDict["%s_energy" % prefix].Fill(simCl1.energy)
        histDict["%s_pt" % prefix].Fill(simCl1.pt)
        histDict["%s_eta" % prefix].Fill(simCl1.eta)
        histDict["%s_phi" % prefix].Fill(simCl1.phi)
        if (prefix == "SimClus"):
            for i, layer in enumerate(simCl1.layers):
                histDict["%s_layers_energy" % prefix].Fill(layer, simCl1.energy*simCl1.fractions[i])
                histDict["%s_cells_energy" % prefix].Fill(simCl1.cells[i], simCl1.energy*simCl1.fractions[i])
                histDict["%s_wafers_energy" % prefix].Fill(simCl1.wafers[i], simCl1.energy*simCl1.fractions[i])
                histDict["%s_fractions_energy" % prefix].Fill(simCl1.fractions[i], simCl1.energy*simCl1.fractions[i])
                histDict["%s_layers_pt" % prefix].Fill(layer, simCl1.pt*simCl1.fractions[i])
                histDict["%s_cells_pt" % prefix].Fill(simCl1.cells[i], simCl1.pt*simCl1.fractions[i])
                histDict["%s_wafers_pt" % prefix].Fill(simCl1.wafers[i], simCl1.pt*simCl1.fractions[i])
                histDict["%s_fractions_pt" % prefix].Fill(simCl1.fractions[i], simCl1.pt*simCl1.fractions[i])
                histDict["%s_layers_fractions" % prefix].Fill(layer, simCl1.fractions[i])
                histDict["%s_cells_fractions" % prefix].Fill(simCl1.cells[i], simCl1.fractions[i])
                histDict["%s_wafers_fractions" % prefix].Fill(simCl1.wafers[i], simCl1.fractions[i])
        usedCluster[j] = True
        multiclusterTLV = ROOT.TLorentzVector()
        multiclusterTLV.SetPtEtaPhiE(
            simCl1.pt, simCl1.eta, simCl1.phi, simCl1.energy)

        # for k, simCl2 in enumerate((clusters)[j + 1:]):
        #     l = j+k+1
        #     logging.debug("l: {} {} {} {}".format(l, simCl2.pt, simCl2.eta, simCl2.phi, clusters[l].pt))
        #     if (simCl2.energy < energyCut):
        #         logging.debug("skip")
        #         continue
        #     if not (usedCluster[l]) and matchesGen[j]:
        #         dR = deltaR(simCl1, simCl2)
        #         logging.debug("unused, DeltaR = {}".format(dR))
        #         if (dR < dRCut):
        #             histDict["%s_energy" % prefix].Fill(simCl1.energy)
        #             histDict["%s_pt" % prefix].Fill(simCl1.pt)
        #             histDict["%s_eta" % prefix].Fill(simCl1.eta)
        #             histDict["%s_phi" % prefix].Fill(simCl1.phi)
        #             if (prefix == "SimClus"):
        #                 histDict["%s_simEnergy" % prefix].Fill(simCl1.simEnergy)
        #                 for i,layer in enumerate(simCl1.layers):
        #                     histDict["%s_layers_energy" % prefix].Fill(layer, simCl1.energy*simCl1.fractions[i])
        #                     histDict["%s_cells_energy" % prefix].Fill(simCl1.cells[i], simCl1.energy*simCl1.fractions[i])
        #                     histDict["%s_wafers_energy" % prefix].Fill(simCl1.wafers[i], simCl1.energy*simCl1.fractions[i])
        #                     histDict["%s_fractions_energy" % prefix].Fill(simCl1.fractions[i], simCl1.energy*simCl1.fractions[i])
        #                     histDict["%s_layers_pt" % prefix].Fill(layer, simCl1.pt*simCl1.fractions[i])
        #                     histDict["%s_cells_pt" % prefix].Fill(simCl1.cells[i], simCl1.pt*simCl1.fractions[i])
        #                     histDict["%s_wafers_pt" % prefix].Fill(simCl1.wafers[i], simCl1.pt*simCl1.fractions[i])
        #                     histDict["%s_fractions_pt" % prefix].Fill(simCl1.fractions[i], simCl1.pt*simCl1.fractions[i])
        #                     histDict["%s_layers_fractions" % prefix].Fill(layer, simCl1.fractions[i])
        #                     histDict["%s_cells_fractions" % prefix].Fill(simCl1.cells[i], simCl1.fractions[i])
        #                     histDict["%s_wafers_fractions" % prefix].Fill(simCl1.wafers[i], simCl1.fractions[i])
        #             histDict["%s_dRtoSeed" % prefix].Fill(dR)
        #             logging.debug("pass cut")
        #             usedCluster[l] = True
        #             tmpTLV = ROOT.TLorentzVector()
        #             tmpTLV.SetPtEtaPhiE(
        #                 simCl2.pt, simCl2.eta, simCl2.phi, simCl2.energy)
        #             multiclusterTLV += tmpTLV
        vMulticlusterTLV.append(multiclusterTLV)
        histDict["multi%s_energy" % prefix].Fill(multiclusterTLV.E())
        histDict["multi%s_pt" % prefix].Fill(multiclusterTLV.Pt())
        histDict["multi%s_eta" % prefix].Fill(multiclusterTLV.Eta())
        histDict["multi%s_phi" % prefix].Fill(multiclusterTLV.Phi())
        logging.debug("multicluster: {} {} {}".format(multiclusterTLV.Pt(), multiclusterTLV.Eta(), multiclusterTLV.Phi()))

    logging.debug("end of loop: {}".format(usedCluster))
    vMulticlusterTLV = sorted(
        vMulticlusterTLV, key=lambda tlv: tlv.Pt(), reverse=True)
    return vMulticlusterTLV


def selectMatchingClusters(refCollection, selCollection, dRcut, histDict, comp):
    """get clusters matching with reference collection within DeltaR."""
    selectedClusters = []
    matchedCluster = [False]*len(refCollection)
    for sel in selCollection:
        for i, ref in enumerate(refCollection):
            if sel.DeltaR(ref) < dRcut:
                selectedClusters.append(sel)
                histDict["multi%s_delta_energy" % comp].Fill(ref.E()-sel.E())
                histDict["multi%s_delta_pt" % comp].Fill(ref.Pt()-sel.Pt())
                histDict["multi%s_deltaover_energy" % comp].Fill((ref.E()-sel.E())/ref.E())
                histDict["multi%s_deltaover_pt" % comp].Fill((ref.Pt()-sel.Pt())/ref.Pt())
                histDict["multi%s_delta_eta" % comp].Fill(ref.Eta()-sel.Eta())
                histDict["multi%s_delta_phi" % comp].Fill(ref.Phi()-sel.Phi())
                histDict["multi%s_delta_R" % comp].Fill(ref.DeltaR(sel))
                matchedCluster[i] = True
    for match in matchedCluster:
        if match:
            histDict["multi%s_eff" % comp].Fill(1)
        else:
            histDict["multi%s_eff" % comp].Fill(0)

    return selectedClusters


def deltaR(p1, p2):
    """calculate DeltaR from ntuple eta-phi values."""
    dphi = p1.phi - p2.phi
    deta = p1.eta - p2.eta
    dR = math.sqrt(dphi * dphi + deta * deta)
    return dR


def deltaR2(tlv1, p2):
    """calculate DeltaR from a TLorentzVector and an ntuple eta-phi value."""
    dphi = tlv1.Phi() - p2.phi
    deta = tlv1.Eta() - p2.eta
    dR = math.sqrt(dphi * dphi + deta * deta)
    return dR


def deltad(p1, p2):
    """calculate Deltad from ntuple x-y values."""
    dx = p1.x - p2.x
    dy = p1.y - p2.y
    dd = math.sqrt(dx * dx + dy * dy)
    return dd


def deltad2(tlv, p2, geometry):
    """calculate Deltad from TLV and ntuple x-y values using geometry."""
    tlvx = geometry.layerEtaPhiToX(p2.layer, tlv.Eta(), tlv.Phi())
    tlvy = geometry.layerEtaPhiToY(p2.layer, tlv.Eta(), tlv.Phi())
    p2x = geometry.layerEtaPhiToX(p2.layer, p2.eta, p2.phi)
    p2y = geometry.layerEtaPhiToY(p2.layer, p2.eta, p2.phi)
    dx = tlvx - p2x
    dy = tlvy - p2y
    dd = math.sqrt(dx * dx + dy * dy)
    return dd


class Geometry(object):
    """Sample class to get ROOT Chain and individual files."""

    def __init__(self):
        """initialise directories."""
        self.z_abs = {}
        self.z_rel = {}
        self.x0_cumulative = {}
        self.x0 = {}
        self.dEdx_cumulative = {}
        self.dEdx = {}
        # drop the b in lambda
        self.lamda_cumulative = {}
        self.lamda = {}

    def addLayer(self, layerGeo):
        """add geometry for a layer following structure in file."""
        layer = int(layerGeo[0])+1
        if layer in self.z_abs:
            logging.error("Layer %i info already added." % layer)
            raise RuntimeError
        self.z_abs[layer] = float(layerGeo[1])
        self.z_rel[layer] = float(layerGeo[2])
        self.x0_cumulative[layer] = float(layerGeo[3])
        self.x0[layer] = float(layerGeo[4])
        self.dEdx_cumulative[layer] = float(layerGeo[5])
        self.dEdx[layer] = float(layerGeo[6])
        self.lamda_cumulative[layer] = float(layerGeo[7])
        self.lamda[layer] = float(layerGeo[8])

    def layerToZ(self, layer, eta):
        """convert layer with eta-information to z value."""
        # if (layer not in self.z_abs):
        #     layer = max(self.z_abs.keys(), key=int)
        z_abs = self.z_abs[layer]
        if (eta < 0):
            z_abs *= -1.
        return z_abs

    def layerEtaPhiToX(self, layer, eta, phi):
        """return absolute X value."""
        z = self.layerToZ(layer, eta)
        t = math.exp(-1. * eta)
        if (t == 1):
            x = 0
        else:
            x = z * 2 * t * math.cos(phi)/(1 - t*t)
        return x

    def layerEtaPhiToY(self, layer, eta, phi):
        """return absolute Y value."""
        z = self.layerToZ(layer, eta)
        t = math.exp(-1. * eta)
        if (t == 1):
            y = 0
        else:
            y = z * 2 * t * math.sin(phi)/(1 - t*t)
        return y


def parseGeometry(geoFilename):
    """use regular expressions to parse geometry file."""
    import re
    fpg = "([0-9]*\.?[0-9]+)"  # floating point group
    regex = ur"S[a-z]\s*(\d+): z=\({}\)\s*({})\s*mm;\s*X0=\(\s*{}\)\s*{}; dEdx=\(\s*{}\)\s*{};\s*l=\(\s*{}\)\s*{}".replace("{}", fpg)
    # test_str = (u"Si 39: z=(4078.1) 46.8 mm; X0=(60.34) 2.81; dEdx=(903.87) 49.54; l=(4.813) 0.261\n"
    #             u"Si  0: z=(3207.5)    0 mm; X0=( 0.88) 0.88; dEdx=( 21.32) 21.32; l=(0.181) 0.181")

    geometry = Geometry()
    with open(geoFilename) as f:
        content = f.readlines()
        for line in content:
            matches = re.finditer(regex, line)
            for matchNum, match in enumerate(matches):
                logging.debug("{}".format(match.groups()))
                geometry.addLayer(match.groups())

    return geometry


def main():
    """main function."""
    print "not implemented."


if __name__ == '__main__':
    main()
