# Test the source functions
# Usage:
#   pytest

import os
import pytest
from dotenv import load_dotenv
from pathlib import Path

env_path = Path(__file__).resolve().parent.parent / ".env"
print(f'Loading env from {env_path}')
if not load_dotenv(dotenv_path=env_path, verbose=True, override=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")

from ROOT import TF1, TF2, gInterpreter
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import SourceGauss, SourceAAA, SourceAAAJC, SourceCountsGauss, SourceCountsAAA, SourceCountsAAAJC

EPSILON = 1.e-12

def test_normalization_SourceAAA():
    fSourceAAA = TF1('fSourceAAA', SourceAAA, 0, 1000, 1)
    fSourceAAA.SetParameter(0, 1)
    fSourceAAA.SetNpx(100000)
    assert abs(fSourceAAA.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceGauss():
    fSourceGauss = TF1('fSourceGauss', SourceGauss, 0, 1000, 1)
    fSourceGauss.SetParameter(0, 1)
    fSourceGauss.SetNpx(100000)
    assert abs(fSourceGauss.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceAAAJC():
    fSourceAAAJC = TF2('fSourceAAAJC', SourceAAAJC, 0, 50, 0, 50, 1)
    fSourceAAAJC.SetParameter(0, 1)
    fSourceAAAJC.SetNpx(10000)
    print("---> ", fSourceAAAJC.Integral(0., 50, 0., 50))
    assert abs(fSourceAAAJC.Integral(0., 50, 0., 50, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAA():
    fSourceCountsAAA = TF1('fSourceCountsAAA', SourceCountsAAA, 0, 1000, 2)
    fSourceCountsAAA.SetParameter(0, 1)
    fSourceCountsAAA.SetParameter(1, 1)
    fSourceCountsAAA.SetNpx(100000)
    assert abs(fSourceCountsAAA.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceCountsGauss():
    fSourceCountsGauss = TF1('fSourceCountsGauss', SourceCountsGauss, 0, 1000, 2)
    fSourceCountsGauss.SetParameter(0, 1)
    fSourceCountsGauss.SetParameter(1, 1)
    fSourceCountsGauss.SetNpx(100000)
    assert abs(fSourceCountsGauss.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceCountsAAAJC():
    fSourceCountsAAAJC = TF2('fSourceCountsAAAJC', SourceCountsAAAJC, 0, 50, 0, 50, 2)
    fSourceCountsAAAJC.SetParameter(0, 1)
    fSourceCountsAAAJC.SetParameter(1, 1)
    fSourceCountsAAAJC.SetNpx(10000)
    print("---> ", fSourceCountsAAAJC.Integral(0., 50, 0., 50))
    assert abs(fSourceCountsAAAJC.Integral(0., 50, 0., 50, EPSILON) - 1) < EPSILON
