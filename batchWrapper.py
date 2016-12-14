"""batch submission wrapping script"""
from test import *
import os
from SampleHelper import SampleManager
import logging
import commands
import time


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def processCmd(cmd, quiet=True):
    if not quiet:
        print cmd
    status, output = commands.getstatusoutput(cmd)
    if (status != 0 and not quiet):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
    return output


def main():
    """Main function and settings."""

    # save working dir
    outDir = '/afs/cern.ch/work/c/clange/HGCal/output'
    currentDir = os.getcwd()
    CMSSW_BASE = os.getenv('CMSSW_BASE')
    CMSSW_VERSION = os.getenv('CMSSW_VERSION')
    SCRAM_ARCH = os.getenv('SCRAM_ARCH')
    geometryFile = currentDir + "/v33-withBH.txt"

    # nEvents = -1
    defaultFilesPerJob = 4
    queue = '1nh'
    batchScript = currentDir + '/' + 'batchScript.sh'
    # dvzCut = -1000  # 320
    # simClusPtCuts = {}
    # simClusPtCuts["chargedPions_nPart1_Pt2_pre15_5k"] = 1.6
    # simClusPtCuts["chargedPions_nPart1_Pt5_pre15_5k"] = 4
    # simClusPtCuts["chargedPions_nPart1_Pt10_pre15_5k"] = 8
    # simClusPtCuts["chargedPions_nPart1_Pt20_pre15_5k"] = 16
    # simClusPtCuts["chargedPions_nPart1_Pt35_pre15_5k"] = 28
    # simClusPtCuts["chargedPions_nPart1_Pt200_pre15_5k"] = 160
    simClusECuts = {}
    simClusECuts["chargedPions_nPart1_E2_pre15_5k"] = 1.8
    simClusECuts["chargedPions_nPart1_E5_pre15_5k"] = 4.5
    simClusECuts["chargedPions_nPart1_E10_pre15_5k"] = 9
    simClusECuts["chargedPions_nPart1_E20_pre15_5k"] = 18
    simClusECuts["chargedPions_nPart1_E40_pre15_5k"] = 36
    simClusECuts["chargedPions_nPart1_E80_pre15_5k"] = 72
    simClusECuts["chargedPions_nPart1_E160_pre15_5k"] = 144
    simClusECuts["chargedPions_nPart1_E320_pre15_5k"] = 288
    # samples2Run = ["chargedPions_nPart1_Pt2_pre15_5k", "chargedPions_nPart1_Pt200_pre15_5k",
    #                "chargedPions_nPart1_Pt5_pre15_5k", "chargedPions_nPart1_Pt10_pre15_5k",
    #                "chargedPions_nPart1_Pt20_pre15_5k", "chargedPions_nPart1_Pt35_pre15_5k"
    #                ]
    samples2Run = ["chargedPions_nPart1_E2_pre15_5k", "chargedPions_nPart1_E5_pre15_5k",
                   "chargedPions_nPart1_E10_pre15_5k", "chargedPions_nPart1_E20_pre15_5k",
                   "chargedPions_nPart1_E40_pre15_5k", "chargedPions_nPart1_E80_pre15_5k",
                   "chargedPions_nPart1_E160_pre15_5k", "chargedPions_nPart1_E320_pre15_5k"
                   ]

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)

    sampleManager = SampleManager()
    for sampleName in samples2Run:
        # sampleName = "chargedPions_nPart1_Pt20_pre15"
        simClusECut = simClusECuts[sampleName]
        if simClusECut > 100:
            filesPerJob = 1
        elif simClusECut > 50:
            filesPerJob = 2
        else:
            filesPerJob = defaultFilesPerJob
        sample = sampleManager.getSample(sampleName)
        sampleFiles = sample.getFiles()
        for i, chunk in enumerate(chunks(sampleFiles, filesPerJob)):
            logger.info("Submitting %s job %d" % (sampleName, i))
            logger.debug("Job %d with files:" % i)
            logger.debug(chunk)
            # submit job
            if not os.path.exists(outDir + '/' + sampleName):
                os.makedirs(outDir + '/' + sampleName)
            if not os.path.exists(outDir + '/' + sampleName + '/std'):
                os.makedirs(outDir + '/' + sampleName + '/std')
            subSampleName = "{sampleName}_{i}".format(sampleName=sampleName, i=i)
            stdOut = "{outDir}/{sampleName}/std/{subSampleName}.out".format(outDir=outDir, sampleName=sampleName, subSampleName=subSampleName)
            stdErr = "{outDir}/{sampleName}/std/{subSampleName}.err".format(outDir=outDir, sampleName=sampleName, subSampleName=subSampleName)
            fileList = ",".join(chunk)
            cmd = 'bsub -o {stdOut} -e {stdErr} -q {queue} -J {jobName} {batchScript} {p1} {p2} {p3} {p4} {p5} {p6} \"{p7}\" {p8} {p9}'.format(stdOut=stdOut, stdErr=stdErr, queue=queue, jobName=subSampleName, batchScript=batchScript, p1=currentDir, p2=CMSSW_BASE, p3=CMSSW_VERSION, p4=SCRAM_ARCH, p5=geometryFile, p6=subSampleName, p7=fileList, p8=simClusECut, p9=outDir)
            processCmd(cmd)
            time.sleep(1)


if __name__ == '__main__':
    main()
