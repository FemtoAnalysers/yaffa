# Test the source functions
# Usage:
#   pytest

import pytest

from ROOT import TF1, gInterpreter
gInterpreter.Declare('#include "../src/RootFunctions.hxx"')
from ROOT import SourceGauss, SourceAAA

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
