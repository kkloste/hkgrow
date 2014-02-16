import sys
import igraph

def read_smat(filename, set_undirected = True):
    """ Load a graph in SMAT format.
    
    G = dgg.read_smat(filename,set_undirected=True)
    loads filename as an igraph graph based on the smat format.
    
    
    @param filename the path to the file
    @param set_undirected a flag to check if the graph is undirected and set igraph
      accordingly.  If this is false, the output is always a directed igraph graph.
    
    @return the igraph graph.  
    """
    f = open(filename, 'rU')
    hdr = f.readline()
    parts = hdr.split()
    nverts = int(parts[0])
    ncols = int(parts[1])
    nnz = int(parts[2])
    if (nverts != ncols):
        raise ValueError(('read_smat needs nrows (%i) = ncols (%i) in smat for graph' % (nverts,
         ncols)))
    edgelist = []
    for line in f:
        parts = line.split()
        edgelist.append((int(parts[0]),
         int(parts[1])))

    G = igraph.Graph(nverts, edges=edgelist, directed=True)
    if (G.ecount() != nnz):
        raise ValueError(('graph has %i edges but smat file has %i nnz' % (G.ecount(),
         nnz)))
    if (G.vcount() > nverts):
        raise ValueError(('graph has more vertices (%i) than smat file has rows %i' % (G.vcount(),
         nverts)))
    if set_undirected:
        dir = False
        for medge in G.is_mutual():
            if (medge is False):
                dir = True
                break

        if (dir or G.to_undirected()):
            pass
    if ((G.vcount() < nverts) and G.add_vertices((nverts - G.vcout()))):
        pass
    return G



def write_layout(xy, fname):
    f = open(fname, 'w')
    for i in xy:
        f.write(('%18.16e %18.16e\n' % (i[0],
         i[1])))

    f.close()


if (__name__ == '__main__'):
    graphfile = sys.argv[1]
    imfile = sys.argv[2]
    xyfile = sys.argv[3]
    type = sys.argv[4]
    G = read_smat(graphfile)
    xy = G.layout(type)
    fig = igraph.Plot()
    fig.add(G, layout=xy, vertex_label=None, vertex_size=5)
    fig.save(imfile)
    write_layout(xy, xyfile)
