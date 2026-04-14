import ROOT
import sys

def is_compatible(val1, val2, threshold=1.e-15):
    return (val1 - val2) / max(min(val1, val2), threshold) < threshold

def get_keys(directory):
    """ Recursively get all keys in a ROOT directory """
    keys = {}
    for key in directory.GetListOfKeys():
        key_name = key.GetName()
        obj = key.ReadObj()
        if obj.InheritsFrom("TDirectory") or obj.InheritsFrom("TDirectoryFile"):
            keys[key_name] = get_keys(obj)  # Recurse into subdirectory
            obj.Close()  # Close directory after use
        else:
            keys[key_name] = obj.ClassName()
            obj.Delete()  # Free memory
    return keys

def compare_histograms(hist1, hist2, path):
    """ Compare two histograms """
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        print(f"Difference in number of bins for {path}")
        return False
    if not is_compatible(hist1.GetEntries(), hist2.GetEntries()):
        print(f"Different number of entries in {path}")
        return False
    if (hist1.GetXaxis().GetTitle() != hist2.GetXaxis().GetTitle() or
        hist1.GetYaxis().GetTitle() != hist2.GetYaxis().GetTitle()):
        print(f"Axis titles differ for {path}")
        return False
    for bin in range(1, hist1.GetNbinsX() + 1):
        bc1 = hist1.GetBinContent(bin)
        bc2 = hist2.GetBinContent(bin)
        if not is_compatible(bc1, bc2):
            print(f"Bin {bin} differs in {path}, bin contents are: {bc1} and {bc2}")
            return False
        be1 = hist1.GetBinError(bin)
        be2 = hist2.GetBinError(bin)
        if not is_compatible(be1, be2):
            print(f"Bin {bin} error differs in {path}, errors are {be1} and {be2}")
            return False
    return True

def compare_objects(obj1, obj2, path):
    """ Compare two ROOT objects based on their type """
    if obj1.ClassName() != obj2.ClassName():
        print(f"Different types at {path}: {obj1.ClassName()} vs {obj2.ClassName()}")
        return False
    if obj1.InheritsFrom("TH1"):  # Compare histograms
        return compare_histograms(obj1, obj2, path)
    else:
        print(f"Unsupported object type {obj1.ClassName()} at {path}")
        return False

def compare_directories(dir1, dir2, path=""):
    """ Recursively compare directories in two ROOT files """
    keys1 = get_keys(dir1)
    keys2 = get_keys(dir2)
    if set(keys1.keys()) != set(keys2.keys()):
        print(f"Different keys found in {path}")
        print(f"Only in first file: {sorted(set(keys1) - set(keys2))}")
        print(f"Only in second file: {sorted(set(keys2) - set(keys1))}")
        return False
    identical = True
    for key in keys1:
        obj1 = dir1.Get(key)
        obj2 = dir2.Get(key)
        if not obj1 or not obj2:
            print(f"Missing object at {path}/{key}")
            identical = False
            continue
        ROOT.SetOwnership(obj1, False)
        ROOT.SetOwnership(obj2, False)
        sub_path = f"{path}/{key}" if path else key
        if isinstance(keys1[key], dict):  # It's a directory
            if not compare_directories(obj1, obj2, sub_path):
                identical = False
        else:  # It's an object
            if not compare_objects(obj1, obj2, sub_path):
                identical = False
    return identical

def compare_root_files(file1_path, file2_path):
    """ Compare two ROOT files """
    file1 = ROOT.TFile.Open(file1_path)
    file2 = ROOT.TFile.Open(file2_path)
    if not file1 or file1.IsZombie() or not file2 or file2.IsZombie():
        print("Error opening files")
        return
    try:
        if compare_directories(file1, file2):
            print("yes")
        else:
            print("Files differ")
    finally:
        file1.Close()
        file2.Close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_root_files.py file1.root file2.root")
        sys.exit(1)
    compare_root_files(sys.argv[1], sys.argv[2])
