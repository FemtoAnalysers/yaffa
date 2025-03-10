#!/bin/bash

rootcling -f CreateDatasetDict.cxx -I$PYTHIA8/include -I$(root-config --incdir) CreateDataset.h CreateDatasetLinkDef.h

g++ -o MakeDistr MakeDistr.cc CreateDatasetDict.cxx \
    $(root-config --cflags) -I$PYTHIA8/include -I$(root-config --incdir) -lGenVector \
    -I/usr/include/yaml-cpp \
    -L$PYTHIA8/lib -L$(root-config --libdir) \
    $(root-config --libs) -lpythia8 -lEG -lMathCore \
    -lCore -lMathMore -lyaml-cpp \
    -std=c++20 -O0 -Wall -Wextra || exit 1

./MakeDistr "/tank2/scratch/share/sim/pythia8311/HMpp13TeV_CRMode2_NonDiffractive/out/PythiaProd_9187.root" "Distr.root" "cfg_make_distr_LPi.yml"
# ./MakeDistr "" Distr.root /home/ktas/ge86rim/an/LPi/sim/pythia/resonances_mT/templates/Sigma1385plus_indirect/cfg_make_distr_LPi.yml
