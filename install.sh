#!/bin/bash

# Check if environment is properly set
if [[ ! -f .env ]]; then
    echo -e "\033[31mError: .env file not found! Exit!\033[0m"
    exit 1
fi
source .env

# Check that the path to YAFFA is defined
if [[ $YAFFA == '' ]]; then
    echo -e "\033[31mError: path to yaffa not found. Check your .env file. Exit!\033[0m"
    exit 1
fi

if [[ $ALIENVLVL == '' ]]; then
    echo -e "\033[31mError: this script must be run inside the ALICE environment. Exit!\033[0m"
    exit 1
fi

# Set the git hooks
git config --local core.hooksPath .githooks/

# Install the executables
# mkdir -p ~/.local/bin/

# cp yaffa/bin/yaffa-lint.sh ~/.local/bin/yaffa-lint
# chmod +x ~/.local/bin/yaffa-lint

# cp yaffa/bin/yaffa-dr.py ~/.local/bin/yaffa-dr
# chmod +x ~/.local/bin/yaffa-dr

# cp yaffa/bin/yaffa-farm-perf.py ~/.local/bin/yaffa-farm-perf
# chmod +x ~/.local/bin/yaffa-farm-perf

# cp yaffa/bin/yaffa-farm-tui.py ~/.local/bin/yaffa-farm-tui
# chmod +x ~/.local/bin/yaffa-farm-tui

# Install yaffa
pushd $YAFFA
pip3 install -e . || pip3 install --user -e .
popd

if [[ $CATSss != '' ]]; then
    mkdir -p build
    pushd build
    cmake -U LOG_LEVEL .. || exit 1
    make || exit 1
    popd
else
    echo -e "\033[33mWARNING: path to CATS is not valid, scripts that depend on it are skipped\033[0m"
fi

echo -e "\033[32mInstallation completed!\033[m"
