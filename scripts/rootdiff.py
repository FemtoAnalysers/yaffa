#!/bin/env python3

import ROOT
import math
import sys
import ctypes


def compare_tf1(f1, f2, tol=1e-12):
    if not f1 or not f2:
        msg = f'  --> f1={f1} and f2={f2}'
        return False, msg

    # Formula
    if f1.GetExpFormula("p") != f2.GetExpFormula("p"):
        msg = f'  --> Different formula: {f1.GetExpFormula("p")} vs {f2.GetExpFormula("p")}'
        return False, msg

    # Range
    xmin1, xmax1 = f1.GetXmin(), f1.GetXmax()
    xmin2, xmax2 = f2.GetXmin(), f2.GetXmax()
    if abs(xmin1 - xmin2) > tol or abs(xmax1 - xmax2) > tol:
        msg = f'  --> Different range: {xmin1}, {xmax1} vs {xmin2}, {xmax2}'
        return False, msg

    # Parameters
    if f1.GetNpar() != f2.GetNpar():
        msg = f'  --> Different number of parameters: {f1.GetNpar()} vs {f2.GetNpar()}'
        return False, msg

    for i in range(f1.GetNpar()):
        if abs(f1.GetParameter(i) - f2.GetParameter(i)) > tol:
            msg = f'  --> Different parameter value: {f1.GetParameter(i)} vs {f2.GetParameter(i)}'
            return False, msg
        if f1.GetParName(i) != f2.GetParName(i):
            msg = f'  --> Different parameter name: {f1.GetParName(i)} vs {f2.GetParName(i)}'
            return False, msg

    return True, None

def compare_hist(h1, h2, tol=1e-12):
    if h1.GetNbinsX() != h2.GetNbinsX():
        print("  --> Different number of bins X")
        return False, None
    if h1.GetNbinsY() != h2.GetNbinsY():
        print("  --> Different number of bins Y")
        return False, None
    if h1.GetNbinsZ() != h2.GetNbinsZ():
        print("  --> Different number of bins Z")
        return False, None

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

            msg = f"  --> Different bin content for bin {i}: x: {ix}({x}), y: {iy}({y}) bc1 = {h1.GetBinContent(i)}, bc2 = {h2.GetBinContent(i)}"
            return False, msg
        if abs(h1.GetBinError(i) - h2.GetBinError(i)) > tol:
            msg = f"  --> Different bin error for bin {i}"
            return False, msg
    return True, None

def compare_graph(g1, g2, tol=1e-12):
    if g1.GetN() != g2.GetN():
        return False

    for i in range(g1.GetN()):
        if abs(g1.GetX()[i] - g2.GetX()[i]) > tol:
            return False
        if abs(g1.GetY()[i] - g2.GetY()[i]) > tol:
            return False

    return True

def compare_files(f1name, f2name):
    print(f"Comparing file '{f1name}' vs '{f2name}'")
    f1 = ROOT.TFile.Open(f1name)
    f2 = ROOT.TFile.Open(f2name)

    keys1 = {k.GetName(): k.GetClassName() for k in f1.GetListOfKeys()}
    keys2 = {k.GetName(): k.GetClassName() for k in f2.GetListOfKeys()}

    if keys1 != keys2:
        only1 = set(keys1) - set(keys2)
        only2 = set(keys2) - set(keys1)


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
                    print(f"Object renamed: {n1} -> {n2}")
                    only1.remove(n1)
                    only2.remove(n2)
                    keys1.pop(n1)
                    keys2.pop(n2)
                    break

        if only1:
            print("\033[33mOnly in file1:\033[0m", only1)
        if only2:
            print("\033[33mOnly in file2:\033[0m", only2)

    are_same = True
    for name, cls in keys1.items():
        o1 = f1.Get(name)
        o2 = f2.Get(name)
        
        if not o2:
            continue

        if o1.InheritsFrom("TH1"):
            ok, msg = compare_hist(o1, o2)
        elif o1.InheritsFrom("TF1"):
            ok, msg = compare_tf1(o1, o2)
        elif o1.InheritsFrom("TGraph"):
            ok = compare_graph(o1, o2)
        else:
            continue

        if not ok:
            print(f" ! \033[31mObject '{name}' ({o1.ClassName()}) is different!\033[0m")
            print(msg)
            are_same = False

    if are_same:
        print("Files are equivalent")
        return True

    print("Files are different")
    return False


if __name__ == "__main__":
    compare_files(sys.argv[1], sys.argv[2])
