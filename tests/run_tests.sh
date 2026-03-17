#!/bin/bash

mkdir -p ../build || exit 1
pushd ../build || exit 1
cmake -DLOG_LEVEL=40 .. || exit 1
make || exit 1
popd

[[ -f ../build/bin/SimulateSource ]] || (echo '[Error] file ../build/bin/SimulateSource does not exist' && exit 1)
time ../build/bin/SimulateSource cfg_test_ceca.yaml || exit 1
../scripts/rootdiff.py source.root source_ref.root || exit 1

time ../build/bin/SimulateSource cfg_test_ceca_3b.yaml || exit 1
../scripts/rootdiff.py source_3b.root source_3b_ref.root || exit 1

