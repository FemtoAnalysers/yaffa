# Test the source functions
# Usage:
#   pytest

import pytest

from ROOT import TF1, TF2, gInterpreter
gInterpreter.Declare('#include "../src/RootFunctions.hxx"')
from ROOT import SourceGauss, SourceAAA, SourceAAAJC

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
