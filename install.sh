#!/bin/bash
# Script to install yaffa

mkdir -p ~/.local/bin/

cp yaffa/bin/yaffa-lint ~/.local/bin/yaffa-lint
chmod +x ~/.local/bin/yaffa-lint

# cp yaffa/utils/dr.py ~/.local/bin/dr
# chmod +x ~/.local/bin/dr

# pip3 install .
