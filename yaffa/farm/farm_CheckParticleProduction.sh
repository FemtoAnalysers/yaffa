#!/bin/bash

currentDir=$PWD

firstJob=$1
shift

lastJob=$1
shift

oDir=/scratch5/ge86rim/an/LPi/sim/pythia/Nstar
mkdir -p $oDir || exit 1

cp $YAFFA/yaffa/sim/pythia/CheckParticleProduction.C $oDir || exit 1
cp $YAFFA/yaffa/farm/CheckParticleProduction.sh $oDir || exit 1

cd $oDir
echo ".L CheckParticleProduction.C+" | root -l -b

echo "Writing output to " $oDir
for ((iJob = $firstJob; iJob <= $lastJob; iJob++)); do
    mkdir -p $oDir/$iJob || continue
    cd $oDir/$iJob
    
    echo "job: $iJob sbatch --time=3:00:00 --exclude=brett,dallas,kane,lambert,monk --mem 2000 CheckParticleProduction.sh $iJob"
    sbatch --time=3:00:00 --exclude=brett,dallas,kane,lambert,monk --mem 2000 ../CheckParticleProduction.sh $iJob
done

cd $currentDir
