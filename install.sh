#!/bin/bash

. bin/check_env.sh || exit 1

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
