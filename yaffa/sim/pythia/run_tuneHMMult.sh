#!/bin/bash

# Compile
echo '.L /home/ktas/ge86rim/phsw/yaffa/yaffa/sim/pythia/TuneHMTrigger.cc+' | root -l -b

dir=~/phsw/yaffa/yaffa/secrets/mult/monash
mkdir -p $dir
parallel "root '/home/ktas/ge86rim/phsw/yaffa/yaffa/sim/pythia/TuneHMTrigger.cc+(10000, kMonash, kNonDiffractive, \"'$dir'/MultDistr_pp13TeV_tune-Monash_process-NonDiffractive_{1}.root\", {1})'" ::: {21..500}
find $dir -name 'MultDistr_pp13TeV_tune-Monash_process-NonDiffractive_*.root' -exec hadd -f $dir/MultDistr_pp13TeV_tune-Monash_process-NonDiffractive.root {} \+

dir=~/phsw/yaffa/yaffa/secrets/mult/cr2
mkdir -p $dir
parallel "root '/home/ktas/ge86rim/phsw/yaffa/yaffa/sim/pythia/TuneHMTrigger.cc+(10000, kCRMode2, kNonDiffractive, \"'$dir'/MultDistr_pp13TeV_tune-CRMode2_process-NonDiffractive_{1}.root\", {1})'" ::: {21..500}
find $dir -name 'MultDistr_pp13TeV_tune-CRMode2_process-NonDiffractive_*.root' -exec hadd -f $dir/MultDistr_pp13TeV_tune-CRMode2_process-NonDiffractive.root {} \+
