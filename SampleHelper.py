"""Sample helper tools for HGCal ntuples on EOS."""
import commands
import ROOT
import logging


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
        self.fileList = fileList

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
        self.sampleList = {}
        if initialise:
            self.sampleList = self.addDefaultSamples()
            print "self.sampleList", self.sampleList
            tmpList = {}
            for name, sample in self.sampleList.items():
                tmpList[name] = self.addSampleFiles(sample)
            self.sampleList = tmpList
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
        self.sampleList[name] = Sample(name, inDir)

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
            break
        return sample

    def addDefaultSamples(self):
        """add default samples."""
        sampleList = {}
        # sampleList["photons_nPart1_Pt5"] = Sample(
        #     "photons_nPart1_Pt5", "partGun_clange_PDGid22_nPart1_Pt5_20160812")
        # sampleList["chargedPions_nPart1_Pt5"] = Sample(
        #     "chargedPions_nPart1_Pt5", "partGun_clange_PDGid211_nPart1_Pt5_20160901")
        # sampleList["chargedPions_nPart1_Pt10"] = Sample(
        #     "chargedPions_nPart1_Pt10", "partGun_clange_PDGid211_nPart1_Pt10_20160901")
        sampleList["chargedPions_nPart1_Pt20"] = Sample(
            "chargedPions_nPart1_Pt20", "partGun_clange_PDGid211_nPart1_Pt20_20160901")
        sampleList["chargedPions_nPart1_Pt35"] = Sample(
            "chargedPions_nPart1_Pt35", "partGun_clange_PDGid211_nPart1_Pt35_20160901")
        return sampleList

    def getSample(self, sampleName):
        """return Sample by name."""
        if sampleName not in self.sampleList:
            logging.error("getSample: Sample " + sampleName + "does not exist")
            self.printSamples()
            return 0
        return self.sampleList[sampleName]

    def getSamples(self):
        """return Sample name list."""
        sampleNameList = []
        for sampleName in self.sampleList:
            sampleNameList.append(sampleName)
        return sampleNameList

    def printSamples(self):
        """print available samples."""
        logging.info("Available samples:")
        for sampleName, sample in self.sampleList.items():
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
