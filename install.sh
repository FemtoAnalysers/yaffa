#!/bin/bash

# Script to install yaffa

# Set the git hooks
git config --local core.hooksPath .githooks/

# Install yaffa
pip3 install -e .

# Install the executables
mkdir -p ~/.local/bin/

cp yaffa/bin/yaffa-lint.sh ~/.local/bin/yaffa-lint
chmod +x ~/.local/bin/yaffa-lint

cp yaffa/bin/yaffa-dr.py ~/.local/bin/yaffa-dr
chmod +x ~/.local/bin/yaffa-dr
