#!/bin/bash

seed=$1

time root -l -b -q '../CheckParticleProduction.C+(100000, {102212, 102214, 102216, 112214, 122212, 202212, 202216, 212212, 212214}, tunes::kCRMode2, "AnalysisResults.root", '$seed')'
