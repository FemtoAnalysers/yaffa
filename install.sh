#!/bin/bash

. bin/check_env.sh || exit 1

# Set the git hooks
git config --local core.hooksPath .githooks/

# Install the executables
mkdir -p ~/.local/bin

# name:path relative to $YAFFA/scripts
SCRIPTS=(
    "CompareGraphs:CompareGraphs.py"
    "ComputeSource:sim/ceca/ComputeSource.py"
    "FitSource:sim/ceca/FitSource.py"
    "FitCF:FitCF.py"
    "ComputeRawCF:ComputeRawCF.py"
    "ComputeWaveFunction:cats/ComputeWaveFunction.py"
    "KooninPratt:KooninPratt.py"
    "RunQA:RunQA.py"
    "BreakUpMomentum:BreakUpMomentum.py"
    "Smear:Smear.py"
    "rootdiff:rootdiff.py"
)

for entry in "${SCRIPTS[@]}"; do
    name=${entry%%:*}
    path=$YAFFA/scripts/${entry#*:}
    chmod +x "$path"
    ln -sf "$path" ~/.local/bin/"$name"
done

# Install yaffa
pushd $YAFFA
pip3 install -e . || pip3 install --user -e .
popd

if [[ $CATS != '' ]]; then
    mkdir -p build
    pushd build
    cmake -U LOG_LEVEL .. || exit 1
    make || exit 1
    popd
else
    echo -e "\033[33mWARNING: path to CATS is not valid, scripts that depend on it are skipped\033[0m"
fi

echo -e "\033[32mInstallation completed!\033[m"
