#!/bin/bash

pushd ../build
make || exit 1
popd

[[ -f ../build/bin/SimulateSource ]] || (echo '[Error] file ../build/bin/SimulateSource does not exist' && exit 1)
../build/bin/SimulateSource 1 60
../scripts/rootdiff.py source.root source_ref.root

