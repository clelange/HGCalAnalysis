sampleDir=$HOME/work/HGCal/output
for sample in `find $sampleDir/*_* -maxdepth 0 -type d`; do
    echo "hadding files in $sample"
    sampleName=`basename $sample`
    hadd -f $sampleDir/$sampleName.root $sample/*.root
done
