'''
A mix of utility functions
'''

def SliceVertically(hist, edges, name=None):
    '''
    Slice a TH2 vertically (ProjectionY) and return the list of slices
    '''

    slices = []
    lowEdges = edges[:-1]
    upEdges = edges[1:]
    if not name:
        name = hist.GetName()

    for lowEdge, upEdge in zip(lowEdges, upEdges):
        firstBin = hist.GetXaxis().FindBin(lowEdge * 1.0001)
        lastBin = hist.GetXaxis().FindBin(upEdge * 0.9999)
        slices.append(hist.ProjectionY(f'{name}_{lowEdge:.0f}_{upEdge:.0f}', firstBin, lastBin))
    return slices
