#!/bin/bash

isMC=$1
file=$2

time root -l -b -q '../RunTaskLPi.C+('$isMC', "'$file'")'
