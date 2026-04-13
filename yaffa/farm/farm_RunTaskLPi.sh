#!/bin/bash

isMC=$1
shift

currentDir=$PWD

if $isMC; then
    files=$(find /tank1/scratch/ge86rim/grid/sim/2021/LHC21b3a/294925/AOD -type f -exec realpath {} \;  | sort -V)
    # oDir=/tank1/scratch/ge86rim/grid/ar/sim/2021/LHC21b3a/294925
    oDir=/tank1/scratch/ge86rim/grid/ar/sim/2021/LHC21b3a/294925_tightpid
else
    files=$(find /tank1/scratch/ge86rim/grid/data/2017/LHC17e/000270824/pass2/AOD -type f -exec realpath {} \;  | sort -V)
    oDir=/tank1/scratch/ge86rim/grid/ar/data/2017/LHC17e/000270824/pass2
fi

mkdir -p $oDir || exit 1

cp $YAFFA/yaffa/runtask/RunTaskLPi.C $oDir || exit 1
cp $YAFFA/yaffa/farm/RunTaskLPi.sh $oDir || exit 1

cd $oDir
echo ".L RunTaskLPi.C++" | root -l -b

echo "Writing output to " $oDir
for file in ${files[@]}; do
    subrun=${file: -15:3}
    mkdir -p $oDir/$subrun || continue
    cd $subrun

    command="sbatch --time=3:00:00 --exclude=brett,dallas,kane,lambert,monk --mem 6000 ../RunTaskLPi.sh $isMC $file"
    
    echo "job: $subrun" $command
    $command
    
    cd ..
done

cd $currentDir
