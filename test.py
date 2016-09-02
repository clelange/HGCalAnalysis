from SampleHelper import SampleManager
import ROOT
import logging


def main():

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)

    sampleManager = SampleManager()
    sample = sampleManager.getSample("chargedPions_nPart1_Pt20")
    chain = sample.getChain()
    nEvents = chain.GetEntries()
    logger.info("Events: %d" % nEvents)
    for i, event in enumerate(chain):
        if (i % 100 == 0):
            logger.info("Event {} of {}".format(i, nEvents))


if __name__ == '__main__':
    main()
