#!/bin/python

import ROOT
import math
import sys
import ctypes


def compare_tf1(f1, f2, tol=1e-12):
    if not f1 or not f2:
        return False

    # Formula
    if f1.GetExpFormula("p") != f2.GetExpFormula("p"):
        return False

    # Range
    xmin1, xmax1 = f1.GetXmin(), f1.GetXmax()
    xmin2, xmax2 = f2.GetXmin(), f2.GetXmax()
    if abs(xmin1 - xmin2) > tol or abs(xmax1 - xmax2) > tol:
        return False

    # Parameters
    if f1.GetNpar() != f2.GetNpar():
        return False

    for i in range(f1.GetNpar()):
        if abs(f1.GetParameter(i) - f2.GetParameter(i)) > tol:
            return False
        if f1.GetParName(i) != f2.GetParName(i):
            return False

    return True

def compare_hist(h1, h2, tol=1e-12):
    if h1.GetNbinsX() != h2.GetNbinsX():
        print("Different number of bins X")
        return False
    if h1.GetNbinsY() != h2.GetNbinsY():
        print("Different number of bins Y")
        return False
    if h1.GetNbinsZ() != h2.GetNbinsZ():
        print("Different number of bins Z")
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

            print(f"[{h1.GetName()}] Different bin content for bin {i}: x: {ix}({x}), y: {iy}({y}) bc1 = {h1.GetBinContent(i)}, bc2 = {h2.GetBinContent(i)}")
            return False
        if abs(h1.GetBinError(i) - h2.GetBinError(i)) > tol:
            print(f"Different bin error for bin {i}")
            return False
    return True

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
    f1 = ROOT.TFile.Open(f1name)
    f2 = ROOT.TFile.Open(f2name)

    keys1 = {k.GetName(): k.GetClassName() for k in f1.GetListOfKeys()}
    keys2 = {k.GetName(): k.GetClassName() for k in f2.GetListOfKeys()}

    if keys1 != keys2:
        only1 = set(keys1) - set(keys2)
        only2 = set(keys2) - set(keys1)

        print("Different object lists!")
        print("Only in file1:", only1)
        print("Only in file2:", only2)

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

        if only1 or only2:
            print("Unmatched objects remain.")
            print("Only in file1:", only1)
            print("Only in file2:", only2)
            return False

    for name, cls in keys1.items():
        o1 = f1.Get(name)
        o2 = f2.Get(name)

        if o1.InheritsFrom("TH1"):
            ok = compare_hist(o1, o2)
        elif o1.InheritsFrom("TF1"):
            ok = compare_tf1(o1, o2)
        elif o1.InheritsFrom("TGraph"):
            ok = compare_graph(o1, o2)
        else:
            print(f"Skipping unsupported object: {name} ({cls})")
            continue

        if not ok:
            print(f"{name} is different")
            sys.exit(1)

    print("Files are equivalent")
    return True

if __name__ == "__main__":
    compare_files(sys.argv[1], sys.argv[2])
