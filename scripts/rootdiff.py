#!/bin/env python3

import ROOT
import math
import sys
import ctypes

def print_error(object, message):
    print(f'\033[31mERROR {object.ClassName()} {object.GetName()} {message}\033[0m')

def print_warning(object, message):
    print(f'\033[33mWARNING {object.ClassName()} {object.GetName()} {message}\033[0m')

def compare_tf1(f1, f2, tol=1e-12):
    if not f1 or not f2:
        print_error(f1, f'Invalid objects: f1={f1} and f2={f2}')
        return False

    # Formula
    if f1.GetExpFormula("p") != f2.GetExpFormula("p"):
        print_error(f1, f'Different formula: {f1.GetExpFormula("p")} vs {f2.GetExpFormula("p")}')
        return False

    # Range
    xmin1, xmax1 = f1.GetXmin(), f1.GetXmax()
    xmin2, xmax2 = f2.GetXmin(), f2.GetXmax()
    if abs(xmin1 - xmin2) > tol or abs(xmax1 - xmax2) > tol:
        print_error(f1, f'Different range: {xmin1}, {xmax1} vs {xmin2}, {xmax2}')
        return False

    # Parameters
    if f1.GetNpar() != f2.GetNpar():
        print_error(f1, f'Different number of parameters: {f1.GetNpar()} vs {f2.GetNpar()}')
        return False

    for i in range(f1.GetNpar()):
        if abs(f1.GetParameter(i) - f2.GetParameter(i)) > tol:
            print_error(f1, f'Different parameter value: {f1.GetParameter(i)} vs {f2.GetParameter(i)}')
            return False
        if f1.GetParName(i) != f2.GetParName(i):
            print_error(f1, f'Different parameter name: {f1.GetParName(i)} vs {f2.GetParName(i)}')
            return False

    return True

def compare_hist(h1, h2, tol=1e-12):
    if h1.GetNbinsX() != h2.GetNbinsX():
        print_error(h1, "Different number of bins X")
        return False
    if h1.GetNbinsY() != h2.GetNbinsY():
        print_error(h1, "Different number of bins Y")
        return False
    if h1.GetNbinsZ() != h2.GetNbinsZ():
        print_error(h1, "Different number of bins Z")
        return False

    for i in range(0, h1.GetNcells()):
        if abs(h1.GetBinContent(i) - h2.GetBinContent(i)) > tol:

            binx = ctypes.c_int()
            biny = ctypes.c_int()
            binz = ctypes.c_int()  # always 0 for TH2

            h2.GetBinXYZ(i, binx, biny, binz)

            ix = binx.value
            iy = biny.value
            
            x = h2.GetXaxis().GetBinCenter(ix)
            y = h2.GetYaxis().GetBinCenter(iy)

            # Recompute errors in case they are zero
            kolmogorov = h1.KolmogorovTest(h2, "")
            if kolmogorov > 0.05:
                print_warning(h1, f"Different bin content, but accodring to Kolmogorov test, histograms are compatible")
                return True
            
            print_error(h1, f"Different bin content, histograms are NOT compatible according to  KolmogorovTest: p={kolmogorov}")
            return False

        if abs(h1.GetBinError(i) - h2.GetBinError(i)) > tol:
            kolmogorov = h1.KolmogorovTest(h2)
            if kolmogorov > 0.05:
                print_warning(h1, f"Different bin error, but accodring to Kolmogorov test, histograms are compatible")
                return True
            
            print_error(h1, f"Different bin error, histograms are NOT compatible according to  KolmogorovTest")
            return False

    return True

def compare_graph(g1, g2, tol=1e-12):
    if g1.GetN() != g2.GetN():
        print_error(g1, 'Different number of data points')
        return False

    for i in range(g1.GetN()):
        if abs(g1.GetX()[i] - g2.GetX()[i]) > tol:
            print_error(g1, 'different X values')
            return False
        if abs(g1.GetY()[i] - g2.GetY()[i]) > tol:
            print_error(g1, 'Different Y values')
            return False
        if g1.ClassName() in ['TGraphErrors', 'TGraphAsymmErrors']:
            if abs(g1.GetErrorX(i) - g2.GetErrorX(i)) > tol:
                print_error(g1, 'Different X errors')
                return False
            if abs(g1.GetErrorY(i) - g2.GetErrorY(i)) > tol:
                print_error(g1, 'Different Y errors')
                return False
        if g1.ClassName() == 'TGraphAsymmErrors':
            if abs(g1.GetErrorXhigh(i) - g2.GetErrorXhigh(i)) > tol:
                print_error(g1, 'Different Xhigh errors')
                return False
            if abs(g1.GetErrorXlow(i) - g2.GetErrorXlow(i)) > tol:
                print_error(g1, 'Different Xlow errors')
                return False
            if abs(g1.GetErrorYhigh(i) - g2.GetErrorYhigh(i)) > tol:
                print_error(g1, 'Different Yhigh errors')
                return False
            if abs(g1.GetErrorYlow(i) - g2.GetErrorYlow(i)) > tol:
                print_error(g1, 'Different Ylow errors')
                return False

    return True

def compare_files(f1name, f2name):
    print(f"Comparing file '{f1name}' vs '{f2name}'")
    f1 = ROOT.TFile.Open(f1name)
    f2 = ROOT.TFile.Open(f2name)

    keys1 = {k.GetName(): k.GetClassName() for k in f1.GetListOfKeys()}
    keys2 = {k.GetName(): k.GetClassName() for k in f2.GetListOfKeys()}

    are_same = keys1 == keys2

    if not are_same:
        only1 = set(keys1) - set(keys2)
        only2 = set(keys2) - set(keys1)
        print_error(f1, f'Different keys: only in file1: {only1}, only in file2: {only2}')


        # Try to detect renames
        for n1 in list(only1):
            o1 = f1.Get(n1)
            for n2 in list(only2):
                o2 = f2.Get(n2)

                same = False

                if o1.InheritsFrom("TH1") and o2.InheritsFrom("TH1"):
                    same = compare_hist(o1, o2)
                elif o1.InheritsFrom("TF1") and o2.InheritsFrom("TF1"):
                    same = compare_tf1(o1, o2)
                elif o1.InheritsFrom("TGraph") and o2.InheritsFrom("TGraph"):
                    same = compare_graph(o1, o2)

                if same:
                    print_warning(o1, f"Object renamed: {n1} -> {n2}")
                    only1.remove(n1)
                    only2.remove(n2)
                    keys1.pop(n1)
                    keys2.pop(n2)
                    break

        if only1:
            print_warning(o1, f'Only in file1: {n1}')
        if only2:
            print_warning(o2, f'Only in file2: {n2}')

    for name in keys1:
        o1 = f1.Get(name)
        o2 = f2.Get(name)
        
        if not o2:
            continue

        if o1.InheritsFrom("TH1"):
            are_same &= compare_hist(o1, o2)
        elif o1.InheritsFrom("TF1"):
            are_same &= compare_tf1(o1, o2)
        elif o1.InheritsFrom("TGraph"):
            are_same &= compare_graph(o1, o2)
        else:
            print_warning(o1, 'Comparison not implemented for this class')
            continue

    if are_same:
        print("Files are equivalent")
        return True

    print("\033[31mFiles are different\033[0m")
    return False


if __name__ == "__main__":
    compare_files(sys.argv[1], sys.argv[2])
