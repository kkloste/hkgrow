#!/bin/bash

# make sure it's all compiled
make all

# run on a random circle
#bin/hypercluster test/randgeom-circ-1000.smat
bin/hypercluster test/stanford-cs-sym.smat
#bin/hypercluster test/wb-cs.stanford.smat 

#bin/pagerank-containment test/stanford-cs-sym.smat test/stanford-cs-sym.hcluster test/stanford-cs-sym.leavetime
