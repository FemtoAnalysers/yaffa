'''
Utility functions for i/o operations on root files
'''

import os

from ROOT import TFile, TDirectoryFile, TList # pylint: disable=import-error

import yaffa

histClasses = [f'TH{n}{t}' for n in [1, 2, 3] for t in ['I', 'F', 'D']]
graphClasses = ['TGraph', 'TGraphErrors', 'TGraphAsymErrors']


def Load(container, path):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TDirecotryFile (TFile)
     - TList

    Parameters
    -----------
    container: TFile of the input file
    path: path of the object inside the root file

    Returns:
    -----------
    obj: target root object
    '''

    # Check that the input file is OK
    path = os.path.normpath(path)
    if container == None:  # pylint: disable=singleton-comparison
        yaffa.logger.critical('The container %s is NULL', container.GetName())

    # Start to extract
    for name in path.split(os.sep):
        yaffa.logger.debug('Trying to load %s:%s. Available keys:', container.GetName(), name)

        for key in GetKeyNames(container):
            yaffa.logger.debug('    %s', key)

        if isinstance(container, TDirectoryFile):
            obj = container.Get(name)
        elif isinstance(container, TList):
            obj = container.FindObject(name)
        else:
            yaffa.logger.critical('The container %s of type %s is not valid', container.GetName(), type(container))

        if obj == None:  # pylint: disable=singleton-comparison
            yaffa.logger.critical('The container %s does not contain an object named %s', container.GetName(), name)
        container = obj

    yaffa.logger.debug('The object %s:%s was succesfully loaded', container.GetName(), path)
    return obj


def LoadCF(pair, **kwargs):
    '''
    Load a correlation function from `secrets`

    Parameters
    ----------
    pair : str
        the name of the partice pair
    author : str, optional
        The (first) author of the publication
    method : str, optional
        The method used to obtain the CF
    version : str, optional
        The version of the calculation

    Returns
    -------
    TGraph
        The loaded CF
    '''

    author = kwargs['author']
    method = kwargs['method']
    version = kwargs['version']
    name = '_'.join([pair, author, method, version])
    inFile = TFile(f'{os.environ.get("YAFFA")}/yaffa/secrets/theory/cf/{pair}.root')
    gCF = inFile.Get(f'g{name}')
    return gCF


def LoadResolutionMatrix(pair, **kwargs):
    '''
    Load a correlation function from `secrets`

    Parameters
    ----------
    pair : str
        the name of the partice pair
    dataset : str, optional
        dataset from which to take the resolution matrix (e.g HMpp13TeV)
    NanoAOD : bool, optional
        use NanoAOD resolution matrix
    version : str, optional
        The version of the calculation

    Returns
    -------
    TGraph
        The loaded CF
    '''

    dataset = kwargs['dataset']
    NanoAOD = 'NanoAOD' if kwargs['NanoAOD'] else ''
    version = kwargs['version']

    name = '_'.join(filter(None, ['ResolutionMatrix', pair, dataset, NanoAOD, version]))
    inFile = TFile(f'{os.environ.get("YAFFA")}/yaffa/secrets/resolution/{name}.root')

    dResMat = {}
    combs = ('02', '13')
    for comb in combs:
        dResMat[f'p{comb}'] = Load(inFile, f'p{comb}/hResolutionMatrixME')
        dResMat[f'p{comb}'].SetDirectory(0)


    return dResMat


def GetHistNames(rdir):
    '''
    Returns the list of histogram names in a TDirectoryFile.
    '''
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetClassName() in histClasses]


def GetObjNames(rdir):
    '''
    Returns the list of object names in a TDirectoryFile.
    '''
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetClassName() != 'TDirectoryFile']


def GetSubdirNames(rdir):
    '''
    Returns the list of subdir names in a TDirectoryFile.
    '''
    if rdir == None: # pylint: disable=singleton-comparison
        return []
    return [key.GetName() for key in list(rdir.GetListOfKeys()) if key.GetClassName() == 'TDirectoryFile']


def GetKeyNames(container):  # pylint: disable=inconsistent-return-statements
    '''
    Returns the list of key names in a TDirectoryFile.
    '''
    if isinstance(container, TDirectoryFile):
        return [key.GetName() for key in list(container.GetListOfKeys())]

    if isinstance(container, TList):
        it = container.MakeIterator()
        names = []
        while True:
            obj = it.Next()
            if obj == None: # pylint: disable=singleton-comparison
                break
            names.append(obj.GetName())

        return names
    yaffa.logger.critical('Unknown container type %s', type(container))


def GetHistsInDir(rdir):
    '''
    Returns the list of histogram objects in a TDirectoryFile.
    '''
    return [rdir.Get(key.GetName()) for key in list(rdir.GetListOfKeys()) if key.GetClassName() in histClasses]


def GetGraphsInDir(rdir):
    '''
    Returns the list of graph objects in a TDirectoryFile.
    '''
    return [rdir.Get(key.GetName()) for key in list(rdir.GetListOfKeys()) if key.GetClassName() in graphClasses]


def GetObjsInDir(rdir):
    '''
    Returns the list of objects in a TDirectoryFile.
    '''
    return [rdir.Get(key.GetName()) for key in list(rdir.GetListOfKeys()) if key.GetClassName() != 'TDirectoryFile']
