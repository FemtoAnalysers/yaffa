#!/usr/bin/python3

'''
Script to evaluate the performance of some jobs. it makes a plot with the
execution time of the jobs run in the batchfarm.

Usage:
In the bash script that is executed by sbatch, prepend "time" to evaluate the
performance, e.g.:
    [[exe.sh]]
    time root -l -b -q myTask.C

The slurm-XXX.out file will end with the execution time. Then run
`yaffa-farm-perf` from the directory that contains all the output

'''

import os

from ROOT import TH1D, TCanvas, gROOT # pylint: disable=no-name-in-module
gROOT.SetBatch(True)

# pylint: disable=anomalous-backslash-in-string
times = os.popen("find . -wholename './[1-9]*/slurm-*.out' -exec tail -n 3 {} \; | grep real | awk '{print $2}'").read()
times = times.split('\n')
times = [time for time in times if time] # remove empty strings

times = [int(time.split("m")[0]) / 60 + float(time.split("m")[1][:-1]) / 3600 for time in times]

nBins = max(35, int(len(times) ** 0.5))
hTimes = TH1D("hTimes", ";Time (h);Counts", 100, 0, 3.5)
for time in times:
    hTimes.Fill(time)

cTimes = TCanvas('cTimes', '', 1800, 600)
hTimes.SetLineWidth(2)
hTimes.Draw()
cTimes.SaveAs('cTimes.png')
