#!/bin/bash

rootcling -f CreateDatasetDict.cxx -I$PYTHIA8/include -I$(root-config --incdir) CreateDataset.h CreateDatasetLinkDef.h

g++ -o MakeDistr MakeDistr.cc CreateDatasetDict.cxx \
    -I$PYTHIA8/include -I$(root-config --incdir) \
    -L$PYTHIA8/lib -L$(root-config --libdir) \
    -lpythia8 $(root-config --libs) \
    -std=c++20

./MakeDistr
