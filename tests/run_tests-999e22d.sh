#!/bin/bash

mkdir -p ../build || exit 1
pushd ../build || exit 1
cmake -DLOG_LEVEL=40 .. || exit 1
make || exit 1
popd

[[ -f ../build/bin/SimulateSource ]] || (echo '[Error] file ../build/bin/SimulateSource does not exist' && exit 1)

# 1 threads
../build/bin/SimulateSource cfg_test_ceca_test-999e22d_1threads.yaml || exit 1
mv source_test-999e22d_1threads.root source_test-999e22d_1threads_ref.root
../build/bin/SimulateSource cfg_test_ceca_test-999e22d_1threads.yaml || exit 1
../scripts/rootdiff.py source_test-999e22d_1threads.root source_test-999e22d_1threads_ref.root || exit 1

# 32 threads
../build/bin/SimulateSource cfg_test_ceca_test-999e22d_32threads.yaml || exit 1
mv source_test-999e22d_32threads.root source_test-999e22d_32threads_ref.root
../build/bin/SimulateSource cfg_test_ceca_test-999e22d_32threads.yaml || exit 1
../scripts/rootdiff.py source_test-999e22d_32threads.root source_test-999e22d_32threads_ref.root || exit 1
