#!/usr/bin/env python
'''
Download AnalysisResults.root from the grid

Install running install.sh
'''

import os
import subprocess
import argparse
import sys

train2pwg = {
    'CF_PbPb': 'CF',
    'CF_PbPb_ESD': 'CF',
    'CF_PbPb_MC': 'CF',
    'CF_PbPb_MC_AOD': 'CF',
    'CF_PbPb_MC_TR': 'CF',
    'CF_pp': 'CF',
    'CF_pPb': 'CF',
    'CF_pPb_ESD': 'CF',
    'CF_pPb_MC': 'CF',
    'CF_pPb_MC_ESD': 'CF',
    'CF_pp_ESD': 'CF',
    'CF_pp_MC': 'CF',
    'CF_pp_MC_ESD': 'CF',
    'DQ_Dimuons_PbPb': 'DQ',
    'DQ_Dimuons_PbPb_ESD': 'DQ',
    'DQ_Dimuons_pp': 'DQ',
    'DQ_Dimuons_pPb': 'DQ',
    'DQ_Dimuons_pPb_ESD': 'DQ',
    'DQ_Dimuons_pPb_MC': 'DQ',
    'DQ_Dimuons_pp_ESD': 'DQ',
    'DQ_Dimuons_pp_MC': 'DQ',
    'DQ_PbPb_AOD': 'DQ',
    'DQ_PbPb_ESD': 'DQ',
    'DQ_PbPb_MC_AOD': 'DQ',
    'DQ_PbPb_MC_ESD': 'DQ',
    'DQ_pp_AOD': 'DQ',
    'DQ_pPb_AOD': 'DQ',
    'DQ_pPb_ESD': 'DQ',
    'DQ_pPb_MC_AOD': 'DQ',
    'DQ_pPb_MC_ESD': 'DQ',
    'DQ_pp_ESD': 'DQ',
    'DQ_pp_MC_AOD': 'DQ',
    'DQ_pp_MC_ESD': 'DQ',
    'GA_PbPb': 'GA',
    'GA_PbPb_AOD': 'GA',
    'GA_PbPb_MC': 'GA',
    'GA_PbPb_MC_AOD': 'GA',
    'GA_pp': 'GA',
    'GA_pp_AOD': 'GA',
    'GA_pPb': 'GA',
    'GA_pPb_AOD': 'GA',
    'GA_pPb_MC': 'GA',
    'GA_pPb_MC_AOD': 'GA',
    'GA_pp_MC': 'GA',
    'GA_pp_MC_AOD': 'GA',
    'D2H_PbPb': 'HF',
    'D2H_PbPb_MC': 'HF',
    'D2H_pp': 'HF',
    'D2H_pPb': 'HF',
    'D2H_pPb_MC': 'HF',
    'D2H_pp_MC': 'HF',
    'Electrons_PbPb': 'HF',
    'Electrons_PbPb_AOD': 'HF',
    'Electrons_PbPb_MC': 'HF',
    'Electrons_PbPb_MC_AOD': 'HF',
    'Electrons_pp_AOD': 'HF',
    'Electrons_pPb': 'HF',
    'Electrons_pPb_AOD': 'HF',
    'Electrons_pPb_MC': 'HF',
    'Electrons_pp_ESD': 'HF',
    'Electrons_pp_MC_AOD': 'HF',
    'Electrons_pp_MC_ESD': 'HF',
    'HFCJ_PbPb': 'HF',
    'HFCJ_PbPb_MC': 'HF',
    'HFCJ_pp': 'HF',
    'HFCJ_pPb': 'HF',
    'HFCJ_pPb_MC': 'HF',
    'HFCJ_pp_MC': 'HF',
    'HF_TreeCreator': 'HF',
    'Muons_PbPb': 'HF',
    'HM_AOD': 'HM',
    'HM_ESD': 'HM',
    'HM_MC_AOD': 'HM',
    'HM_MC_ESD': 'HM',
    'Jets_EMC_PbPb': 'JE',
    'Jets_EMC_PbPb_MC': 'JE',
    'Jets_EMC_pp': 'JE',
    'Jets_EMC_pPb': 'JE',
    'Jets_EMC_pp_MC': 'JE',
    'Jets_PbPb': 'JE',
    'Jets_PbPb_2011': 'JE',
    'Jets_PbPb_AOD': 'JE',
    'Jets_pp': 'JE',
    'Jets_pp_AOD': 'JE',
    'Jets_pp_MC': 'JE',
    'LF_PbPb': 'LF',
    'LF_PbPb_AOD': 'LF',
    'LF_PbPb_MC': 'LF',
    'LF_PbPb_MC_AOD': 'LF',
    'LF_pp': 'LF',
    'LF_pp_AOD': 'LF',
    'LF_pPb': 'LF',
    'LF_pPb_AOD': 'LF',
    'LF_pPb_MC': 'LF',
    'LF_pPb_MC_AOD': 'LF',
    'LF_pp_MC': 'LF',
    'LF_pp_MC_AOD': 'LF',
    'LF_QA': 'LF',
    'LF_QA_AOD': 'LF',
    'LF_QA_MC': 'LF',
    'LF_QA_MC_AOD': 'LF',
    'MM_MCGen_validation': 'MM',
    'MM_pp_ESD': 'MM',
    'MM_pp_MC_ESD': 'MM',
    'AdHocCalibration': 'PP',
    'AnalysisQA_AOD': 'PP',
    'AnalysisQA_ESD': 'PP',
    'DPG_AOD': 'PP',
    'DPG_ESD': 'PP',
    'DPG_MC_AOD': 'PP',
    'DPG_MC_ESD': 'PP',
    'DPG_pPb': 'PP',
    'PP_EMCAL_Calibration': 'PP',
    'QATrain': 'PP',
    'Upgrade_ITS_MFT': 'PP',
    'UD_PbPb_AOD': 'UD',
    'UD_PbPb_ESD': 'UD',
    'UD_pp_AOD': 'UD',
    'UD_pPb_AOD': 'UD',
    'UD_pp_ESD': 'UD',
    'Devel_1': 'ZZ',
    'Devel_2': 'ZZ',
    'Devel_3': 'ZZ',
    'MCGenDev': 'ZZ',
    'MCGen_PbPb': 'ZZ',
    'MCGen_pp': 'ZZ',
    'MCGen_pPb': 'ZZ',
    'NanoAOD_Analysis': 'ZZ',
    'NanoAOD_Filtering': 'ZZ',
    'PbPb_2015': 'ZZ',
    'PbPb_2015_ESD': 'ZZ',
    'PbPb_2015_MC': 'ZZ',
    'PbPb_2015_MC_ESD': 'ZZ',
    'Run3_Benchmarks': 'ZZ',
    'Run3_Conversio': 'ZZ',
}

parser = argparse.ArgumentParser()
parser.add_argument('train', choices=train2pwg.keys(), help='Name of the train from which file should be downloaded')
parser.add_argument('run', help='Number of the run')
args = parser.parse_args()

train = args.train
run = args.run
pwg = train2pwg[train]

if not os.environ.get('ALIENVLVL'):
    print("This script needs to be run from within the ALICE environment (AliPhysics/O2/O2Physics). Exit!")
    sys.exit(1)

# Search the files
path = f'/alice/cern.ch/user/a/alitrain/PWG{pwg}/{train}/{run}_.*'
try:
    files = subprocess.check_output(f'alien_find -r "{path}/AnalysisResults.root"', shell=True, universal_newlines='\n')
except subprocess.CalledProcessError:
    print('No files found in %s. Exit!', path)
    sys.exit(1)

files = str(files).split('\n')
files.remove('')

children = []
merged = []
for name in files:
    if 'child' in name:
        children.append(name)
    else:
        merged.append(name)

print('Found the following files:\nmerged:')
for name in merged:
    print('   ', name)
print('\nChildren:')
for name in children:
    print('   ', name)

choice = None
while choice not in ['m', 'c']:
    choice = input('Download type: merged (m) or children (c): ')

if choice == 'm':
    if len(merged) > 1:
        sys.exit()
    os.system(f'alien_cp {merged[0]} file:AnalysisResults_{run}.root')

if choice == 'c':
    for name in children:
        start = name.index('child_')+6
        child = name[start:start+name[start:].index('/')]

        os.system(f'alien_cp {name} file:merge_{run}/child_{child}/AnalysisResults.root')
