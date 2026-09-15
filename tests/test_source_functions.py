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
    SourcePdfGauss,
    SourcePdfAAAHypRad,
    SourcePdfAAAJC,
    SourcePdfAAAppr,
    SourcePdfAAAHypAngle,
    SourcePdfAAApprHypAngle,
    SourcePdfAAAprrHypAngle,
    SourcePdfAAAGaussResonancesHypAngle,
)

EPSILON = 1.e-12

def test_normalization_SourcePdfAAAHypRad():
    fSourceAAA = TF1('fSourceAAA', SourcePdfAAAHypRad, 0, 1000, 1)
    fSourceAAA.SetParameter(0, 2.5)
    fSourceAAA.SetNpx(100000)
    assert abs(fSourceAAA.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourcePdfGauss():
    fSourceGauss = TF1('fSourceGauss', SourcePdfGauss, 0, 1000, 1)
    fSourceGauss.SetParameter(0, 1.3)
    fSourceGauss.SetNpx(100000)
    assert abs(fSourceGauss.Integral(0, 1000) - 1) < EPSILON

def test_normalization_SourcePdfAAAJC():
    fSourceAAAJC = TF2('fSourceAAAJC', SourcePdfAAAJC, 0, 50, 0, 50, 1)
    fSourceAAAJC.SetParameter(0, 1.2)
    fSourceAAAJC.SetNpx(10000)
    print("---> ", fSourceAAAJC.Integral(0., 50, 0., 50))
    assert abs(fSourceAAAJC.Integral(0., 50, 0., 50, EPSILON) - 1) < EPSILON

def test_normalization_SourcePdfAAAppr():
    fSourceAAAppr = TF2("SourceAAAppr", SourcePdfAAAppr, 0, 30, 0, np.pi / 2,  2)
    fSourceAAAppr.SetParameter(0, 1.4)
    fSourceAAAppr.SetParameter(1, 2.3)
    fSourceAAAppr.SetNpx(10000)
    fSourceAAAppr.SetNpy(10000)

    assert abs(fSourceAAAppr.Integral(0, 30, 0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourcePdfAAApprHypAngle():
    fSourceAAApprHypAngle = TF1("SourceAAApprHypAngle", SourcePdfAAApprHypAngle, 0, np.pi / 2, 2)
    fSourceAAApprHypAngle.SetParameter(0, 1.4)
    fSourceAAApprHypAngle.SetParameter(1, 2.3)
    fSourceAAApprHypAngle.SetNpx(100000)

    assert abs(fSourceAAApprHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourcePdfAAAHypAngle():
    fSourceAAAHypAngle = TF1("SourceAAAHypAngle", SourcePdfAAAHypAngle, 0, np.pi / 2, 0)
    fSourceAAAHypAngle.SetNpx(100000)

    assert abs(fSourceAAAHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourcePdfAAAprrHypAngle():
    fSourceAAAprrHypAngle = TF1("SourceAAAprrHypAngle", SourcePdfAAAprrHypAngle, 0, np.pi / 2, 2)
    fSourceAAAprrHypAngle.SetParameter(0, 1.4)
    fSourceAAAprrHypAngle.SetParameter(1, 2.3)
    fSourceAAAprrHypAngle.SetNpx(100000)

    assert abs(fSourceAAAprrHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON

def test_normalization_SourcePdfAAAGaussResonancesHypAngle():
    fSourceAAAGaussResonancesHypAngle = TF1(
        "SourceAAAGaussResonancesHypAngle", SourcePdfAAAGaussResonancesHypAngle, 0, np.pi / 2, 3
    )
    fSourceAAAGaussResonancesHypAngle.SetParameter(0, 0.3578)
    fSourceAAAGaussResonancesHypAngle.SetParameter(1, 1.4)
    fSourceAAAGaussResonancesHypAngle.SetParameter(2, 2.3)
    fSourceAAAGaussResonancesHypAngle.SetNpx(100000)

    assert abs(fSourceAAAGaussResonancesHypAngle.Integral(0, np.pi / 2, EPSILON) - 1) < EPSILON
