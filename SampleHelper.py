"""Sample helper tools for HGCal ntuples on EOS."""
import commands
import ROOT
import logging
from copy import copy


class NullHandler(logging.Handler):
    """NullHandler for logging module."""

    def emit(self, record):
        """emit."""
        pass

logging.getLogger(__name__).addHandler(NullHandler())


class Sample(object):
    """Sample class to get ROOT Chain and individual files."""

    def __init__(self, name, inDir, treeName="ana/hgc", fileList=[]):
        """basic sample settings."""
        self.name = name
        self.inDir = inDir
        self.treeName = treeName
        # fileList needs copy since it's a reference
        self.fileList = copy(fileList)

    def getFiles(self, numberOfFiles=-1):
        """get input files."""
        fileList = []
        for i, fileName in enumerate(self.fileList):
            fileList.append(fileName)
            if numberOfFiles > 0:
                if (i == numberOfFiles):
                    break
        return fileList

    def addFile(self, fileName):
        """add an individual file to sample."""
        logging.debug("addFile: {} - {}".format(self.name, fileName))
        if fileName not in self.fileList:
            logging.debug("addFile: Adding " + fileName)
            self.fileList.append(fileName)
        else:
            logging.error("file " + fileName +
                          " has already been added before.")

    def getChain(self):
        """create a ROOT Chain from sample."""
        chain = ROOT.TChain(self.treeName)
        for fileName in self.fileList:
            logging.debug("getChain: Adding " + fileName)
            chain.AddFile(fileName)
        return chain


class SampleManager(object):
    """Sample class to get ROOT Chain and individual files."""

    def __init__(self, initialise=True):
        """basic location settings and initialisation.

        mind the slashes in the path names.
        """
        self.eosExec = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'
        self.eosPrefix = "root://eoscms.cern.ch/"
        self.baseDir = "/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/"
        self.pathSuffix = "/NTUP/"
        self.sampleDict = {}
        if initialise:
            self.sampleDict = self.addDefaultSamples()
            tmpList = {}
            for name, sample in self.sampleDict.items():
                tmpList[name] = self.addSampleFiles(sample)
            self.sampleDict = tmpList
            self.printSamples()

    def getFullPath(self, inDir, forRoot=False):
        """return full path."""
        fullPath = "{baseDir}{inDir}{pathSuffix}".format(
            baseDir=self.baseDir, inDir=inDir, pathSuffix=self.pathSuffix)
        if forRoot:
            fullPath = fullPath + self.eosPrefix
        return fullPath

    def addSample(self, name, inDir):
        """add an individual sample."""
        self.sampleDict[name] = Sample(name, inDir)

    def addSampleFiles(self, sample):
        """add files to a Sample."""
        logging.debug(
            "addSampleFiles - sample: {}, inDir: {}".format(sample.name, sample.inDir))
        searchPath = self.getFullPath(sample.inDir)
        eosCmd = "{eosExec} ls {searchPath}".format(
            eosExec=self.eosExec, searchPath=searchPath)
        eosOutput = processCmd(eosCmd)
        for line in eosOutput.split("\n"):
            sample.addFile(self.eosPrefix + searchPath + line)
        return sample

    def addDefaultSamples(self):
        """add default samples."""
        sampleDict = {}
        sampleDict["photons_nPart1_Pt5_pre12"] = Sample(
            "photons_nPart1_Pt5_pre12", "partGun_clange_PDGid22_nPart1_Pt5_20160812")
        sampleDict["chargedPions_nPart1_Pt5_pre12"] = Sample(
            "chargedPions_nPart1_Pt5_pre12", "partGun_clange_PDGid211_nPart1_Pt5_20160901")
        sampleDict["chargedPions_nPart1_Pt10_pre12"] = Sample(
            "chargedPions_nPart1_Pt10_pre12", "partGun_clange_PDGid211_nPart1_Pt10_20160901")
        sampleDict["chargedPions_nPart1_Pt20_pre12"] = Sample(
            "chargedPions_nPart1_Pt20_pre12", "partGun_clange_PDGid211_nPart1_Pt20_20160901")
        sampleDict["chargedPions_nPart1_Pt35_pre12"] = Sample(
            "chargedPions_nPart1_Pt35_pre12", "partGun_clange_PDGid211_nPart1_Pt35_20160901")
        # pre15 samples
        # sampleDict["photons_nPart1_Pt5_pre15"] = Sample(
        #     "photons_nPart1_Pt5_pre15", "partGun_predragm_PDGid22_nPart1_Pt5_pre15_20161024")
        sampleDict["chargedPions_nPart1_Pt5_pre15"] = Sample(
            "chargedPions_nPart1_Pt5_pre15", "partGun_predragm_PDGid211_nPart1_Pt5_pre15_20161024")
        sampleDict["chargedPions_nPart1_Pt10_pre15"] = Sample(
            "chargedPions_nPart1_Pt10_pre15", "partGun_predragm_PDGid211_nPart1_Pt10_pre15_20161024")
        sampleDict["chargedPions_nPart1_Pt20_pre15"] = Sample(
            "chargedPions_nPart1_Pt20_pre15", "partGun_predragm_PDGid211_nPart1_Pt20_pre15_20161024")
        sampleDict["chargedPions_nPart1_Pt35_pre15"] = Sample(
            "chargedPions_nPart1_Pt35_pre15", "partGun_predragm_PDGid211_nPart1_Pt35_pre15_20161024")
        return sampleDict

    def getSample(self, sampleName):
        """return Sample by name."""
        if sampleName not in self.sampleDict:
            logging.error("getSample: Sample " + sampleName + "does not exist")
            self.printSamples()
            return 0
        return self.sampleDict[sampleName]

    def getSamples(self):
        """return Sample name list."""
        sampleNameList = []
        for sampleName in self.sampleDict:
            sampleNameList.append(sampleName)
        return sampleNameList

    def printSamples(self):
        """print available samples."""
        logging.info("Available samples:")
        for sampleName, sample in self.sampleDict.items():
            logging.info(" o " + sampleName +
                         " (" + str(len(sample.fileList)) + " files)")


def processCmd(cmd, quiet=False):
    """processing the external os commands."""
    logging.debug("processCmd: " + cmd)
    status, output = commands.getstatusoutput(cmd)
    if (status != 0 and not quiet):
        logging.error('Error in processing command:\n   [' + cmd + ']')
        logging.error('Output:\n   [' + output + '] \n')
    return output
