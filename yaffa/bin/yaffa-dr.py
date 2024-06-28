#!/usr/bin/env python
'''
Download AnalysisResults.root from the grid

Install running install.sh
'''

import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('train', help='Name of the train from which file should be downloaded')
parser.add_argument('run', help='Number of the run')
args = parser.parse_args()

train = args.train
run = args.run

if not os.environ.get('ALIENVLVL'):
    print("This script needs to be run from within the ALICE environment (AliPhysics/O2/O2Physics). Exit!")
    sys.exit(1)

# Check that the train is valid
if 'CF' in train:
    pwg = 'CF'
elif 'D2H' in train:
    pwg = 'HF'
else:
    print('Unknown train. Exit!')
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
