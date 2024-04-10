#!/bin/bash

# Script to install yaffa

if [ $YAFFA == '']; then
    echo -e "\033[31mError: the path to yaffa was not found. Add `export YAFFA=/path/to/yaffa` to your .bashrc. Exit!\033[0m"
    exit 1
fi

currentDir=$PWD
cd $YAFFA

# Set the git hooks
git config --local core.hooksPath .githooks/

# Install yaffa
pip3 install --user -e .

# Install the executables
mkdir -p ~/.local/bin/

cp yaffa/bin/yaffa-lint.sh ~/.local/bin/yaffa-lint
chmod +x ~/.local/bin/yaffa-lint

cp yaffa/bin/yaffa-dr.py ~/.local/bin/yaffa-dr
chmod +x ~/.local/bin/yaffa-dr

cp yaffa/bin/yaffa-farm-perf.py ~/.local/bin/yaffa-farm-perf
chmod +x ~/.local/bin/yaffa-farm-perf

cd $currentDir
