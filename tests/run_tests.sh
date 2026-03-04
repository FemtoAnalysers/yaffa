#!/bin/bash

pushd ../build
cmake -DLOG_LEVEL=40 .. || exit 1
make || exit 1
popd

[[ -f ../build/bin/SimulateSource ]] || (echo '[Error] file ../build/bin/SimulateSource does not exist' && exit 1)
../build/bin/SimulateSource cfg_test_ceca.yaml || exit 1
../scripts/rootdiff.py source.root source_ref.root || exit 1

