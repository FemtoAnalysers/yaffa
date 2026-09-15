# Test the source functions
# Usage:
#   pytest

import os
import pytest
from dotenv import load_dotenv
from pathlib import Path

import numpy as np

env_path = Path(__file__).resolve().parent.parent / ".env"
print(f'Loading env from {env_path}')
if not load_dotenv(dotenv_path=env_path, verbose=True, override=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")

from ROOT import TF1, TF2, gInterpreter
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import (
    SourceGauss,
    SourceAAA,
    SourceAAAJC,
    SourceCountsGauss,
    SourceCountsAAA,
    SourceCountsAAAJC,
    SourceCountsAAAppr,
    SourceCountsAAAHypAngle,
    SourceCountsAAApprHypAngle,
    SourceCountsAAAprrHypAngle,
    SourceCountsAAAGaussResonancesHypAngle,
)

EPSILON = 1.e-12

def test_normalization_SourceAAA():
    fSourceAAA = TF1('fSourceAAA', SourceAAA, 0, 1000, 1)
    fSourceAAA.SetParameter(0, 2.5)
    fSourceAAA.SetNpx(100000)
    assert abs(fSourceAAA.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceGauss():
    fSourceGauss = TF1('fSourceGauss', SourceGauss, 0, 1000, 1)
    fSourceGauss.SetParameter(0, 1.3)
    fSourceGauss.SetNpx(100000)
    assert abs(fSourceGauss.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceAAAJC():
    fSourceAAAJC = TF2('fSourceAAAJC', SourceAAAJC, 0, 50, 0, 50, 1)
    fSourceAAAJC.SetParameter(0, 1.2)
    fSourceAAAJC.SetNpx(10000)
    print("---> ", fSourceAAAJC.Integral(0., 50, 0., 50))
    assert abs(fSourceAAAJC.Integral(0., 50, 0., 50, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAA():
    fSourceCountsAAA = TF1('fSourceCountsAAA', SourceCountsAAA, 0, 1000, 2)
    fSourceCountsAAA.SetParameter(0, 1)
    fSourceCountsAAA.SetParameter(1, 2.5)
    fSourceCountsAAA.SetNpx(100000)
    assert abs(fSourceCountsAAA.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceCountsGauss():
    fSourceCountsGauss = TF1('fSourceCountsGauss', SourceCountsGauss, 0, 1000, 2)
    fSourceCountsGauss.SetParameter(0, 1)
    fSourceCountsGauss.SetParameter(1, 1.3)
    fSourceCountsGauss.SetNpx(100000)
    assert abs(fSourceCountsGauss.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourceCountsAAAJC():
    fSourceCountsAAAJC = TF2('fSourceCountsAAAJC', SourceCountsAAAJC, 0, 50, 0, 50, 2)
    fSourceCountsAAAJC.SetParameter(0, 1)
    fSourceCountsAAAJC.SetParameter(1, 1.2)
    fSourceCountsAAAJC.SetNpx(10000)
    print("---> ", fSourceCountsAAAJC.Integral(0., 50, 0., 50))
    assert abs(fSourceCountsAAAJC.Integral(0., 50, 0., 50, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAAppr():
    fSourceCountsAAAppr = TF2("SourceAAAppr", SourceCountsAAAppr, 0, 30, 0, np.pi / 2,  3)
    fSourceCountsAAAppr.SetParameter(0, 1)
    fSourceCountsAAAppr.SetParameter(1, 1.4)
    fSourceCountsAAAppr.SetParameter(2, 2.3)
    fSourceCountsAAAppr.SetNpx(10000)
    fSourceCountsAAAppr.SetNpy(10000)

    assert abs(fSourceCountsAAAppr.Integral(0, 30, 0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAApprHypAngle():
    fSourceCountsAAApprHypAngle = TF1("SourceAAApprHypAngle", SourceCountsAAApprHypAngle, 0, np.pi / 2, 3)
    fSourceCountsAAApprHypAngle.SetParameter(0, 1)
    fSourceCountsAAApprHypAngle.SetParameter(1, 1.4)
    fSourceCountsAAApprHypAngle.SetParameter(2, 2.3)
    fSourceCountsAAApprHypAngle.SetNpx(100000)

    assert abs(fSourceCountsAAApprHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAAHypAngle():
    fSourceCountsAAAHypAngle = TF1("SourceAAAHypAngle", SourceCountsAAAHypAngle, 0, np.pi / 2, 1)
    fSourceCountsAAAHypAngle.SetParameter(0, 1)
    fSourceCountsAAAHypAngle.SetNpx(100000)

    assert abs(fSourceCountsAAAHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAAprrHypAngle():
    fSourceCountsAAAprrHypAngle = TF1("SourceAAAprrHypAngle", SourceCountsAAAprrHypAngle, 0, np.pi / 2, 3)
    fSourceCountsAAAprrHypAngle.SetParameter(0, 1)
    fSourceCountsAAAprrHypAngle.SetParameter(1, 1.4)
    fSourceCountsAAAprrHypAngle.SetParameter(2, 2.3)
    fSourceCountsAAAprrHypAngle.SetNpx(100000)

    assert abs(fSourceCountsAAAprrHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourceCountsAAAGaussResonancesHypAngle():
    fSourceCountsAAAGaussResonancesHypAngle = TF1(
        "SourceAAAGaussResonancesHypAngle", SourceCountsAAAGaussResonancesHypAngle, 0, np.pi / 2, 4
    )
    fSourceCountsAAAGaussResonancesHypAngle.SetParameter(0, 1)
    fSourceCountsAAAGaussResonancesHypAngle.SetParameter(1, 0.3578)
    fSourceCountsAAAGaussResonancesHypAngle.SetParameter(2, 1.4)
    fSourceCountsAAAGaussResonancesHypAngle.SetParameter(3, 2.3)
    fSourceCountsAAAGaussResonancesHypAngle.SetNpx(100000)

    assert abs(fSourceCountsAAAGaussResonancesHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON
