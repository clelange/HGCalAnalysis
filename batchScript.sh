## External vars
currentDir=${1}
CMSSWDIR=${2}
CMSSWVER=${3}
CMSSWARCH=${4}
geometryFile=${5}
sampleName=${6}
fileList=${7}
eCut=${8}
outDir=${9}

##Create Work Area
cd $TMPDIR
export SCRAM_ARCH=${CMSSWARCH}
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 project CMSSW ${CMSSWVER}`
cd ${CMSSWVER}/
rm -rf ./*
cp -r -d ${CMSSWDIR}/* ./
cd src
eval `scramv1 runtime -sh`
edmPluginRefresh -p ../lib/$SCRAM_ARCH

## Execute job and retrieve the outputs
echo "Job running on `hostname` at `date`"

cd RecoNtuples/HGCalAnalysis/test/HGCalAnalysis
pwd
echo python test.py --geometry $geometryFile --sampleName $sampleName --files $fileList --eCut $eCut
python test.py --geometry $geometryFile --sampleName $sampleName --files $fileList --eCut $eCut
ls

# copy to outDir
mkdir -p $outDir/${sampleName%_*}
cp ${sampleName%_*}*.root $outDir/${sampleName%_*}/
