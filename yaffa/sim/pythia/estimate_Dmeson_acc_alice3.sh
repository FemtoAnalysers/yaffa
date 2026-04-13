#!/bin/bash

echo ".L /home/ktas/ge86rim/phsw/yaffa/yaffa/sim/pythia/EstimateDmesonAccALICE3.C+" | root -l -b

parallel "root -l -b -q '/home/ktas/ge86rim/phsw/yaffa/yaffa/sim/pythia/EstimateDmesonAccALICE3.C+(50000, {2}, \"$YAFFA/yaffa/secrets/alice3/DmesonAccALICE3_{2}_{1}.root\", {1})'" ::: {1..20} ::: 421 413

hadd DmesonAccALICE3_Dzero.root DmesonAccALICE3_421_*.root
hadd DmesonAccALICE3_Dstar.root DmesonAccALICE3_413_*.root

cp