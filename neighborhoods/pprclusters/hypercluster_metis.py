#!/usr/bin/env python

"""
hypercluster_metis.py
=====================

Construct a hyperclustering of a graph with metis

History
-------
:2010-10-19: Initial version
"""

import sys
import os
import time
import random
import subprocess
from math import ceil, floor
from optparse import OptionParser

from overlap_help import *

""" Return a parser for the command line options. """
def options_setup():
  usage = "usage: %prog [options] graphfilename"
  parser = OptionParser(usage=usage)
  parser.add_option("","--maxvol",dest="maxvol",metavar="int|float",
    help="the maximum ideal volume of any cluster expressed as " + 
    "a number of edges (int) or a fraction of total edges (float) " + 
    "the default is max(0.01*nedges,100)",
    default=None,type="float")
  parser.add_option("","--minvol",dest="minvol",metavar="int|float",
    help="the minimum volume of any cluster produced expressed "+
    "the same way as maxvol, the default is 100", default=100,type="float")
  parser.add_option("-o","--out",dest="outfilename",metavar="FILE",
    help = "the output filename for the hyperclustering, " + 
    "the default value is graphname.hcluster",default=None)
  parser.add_option("-n","--nruns", dest="nruns", metavar="int",
    help="the number of times to run metis, default is 2",
    default=2,type="int")
  parser.add_option("-v","--verbose",action="store_true",dest="verbose")
  parser.add_option("-s","--seed",dest="seed",default=0,type="int")
  return parser
  
def main():
  parser = options_setup()
  (options,args) = parser.parse_args()
  if len(args) != 1:
    parser.error("incorrect number of arguments")
  
  verbose = options.verbose
  
  graphfilename = args[0]
  
  if verbose:
    print "reading %s ...."%(graphfilename)
  graphfile = open(graphfilename,"rt")
  header = graphfile.readline()
  graphfile.close()
  parts = header.split()
  nverts = int(parts[0])
  nedges = int(parts[1])
  if verbose:
    print "nverts=%i; nedges=%i"%(nverts,nedges)
    
  # find the path to metis
  metisexe = os.path.abspath( os.path.join(__file__,'..','..',
                'metis-5.0pre2','build','Linux-x86_64','pmetis5.0pre2'))
                
  if verbose:
    print "metisexe: %s"%(metisexe)
  
  if not os.path.isfile(metisexe):
    print "metis command: %s does not exist or is not a file"%(metisexe)
    sys.exit(1)
    
  if options.seed < 0:
    options.seed = int(time.time())
    
  if options.outfilename is None:
    (base,ext) = os.path.splitext(graphfilename)
    options.outfilename = base + ".hcluster"
    
  if options.maxvol is None:
    options.maxvol = max(ceil(0.01*nedges),100)
  elif options.maxvol < 1:
    options.maxvol = ceil(nedges*options.maxvol)
    
  if options.minvol is None:
    options.minvol = 100
  elif options.minvol < 1:
    options.minvol = ceil(nedges*options.minvol)
  
  maxparts = floor(nedges/float(options.minvol))
  minparts = ceil(nedges/float(options.maxvol))
  maxparts = min(maxparts,minparts*3);
  if maxparts<minparts:
     maxparts=minparts
  
  # print the parameters
  print "     nruns: %i"%(options.nruns)
  print "    maxvol: %i"%(options.maxvol)
  print "    minvol: %i"%(options.minvol)
  print "  maxparts: %i"%(maxparts)
  print "  minparts: %i"%(minparts)
  print "      seed: %i"%(options.seed)
  print "    output: %s"%(options.outfilename)
  
  random.seed(options.seed)
  
  assert(options.nruns > 0)
  
  clusters = []
  nclusters = 0
  
  for run in xrange(options.nruns):
    np = random.randint(minparts,maxparts)
    mseed = random.randint(0,65536)
    if run is 0:
      np = minparts
    elif run is 1:
      np = maxparts
    
    # run metis
    metiscmd = "%s -seed %i %s %i"%(metisexe, mseed, graphfilename, np)
    
    if verbose:
      print "running: %s"%(metiscmd)
    
    rval,out,err = run_command(metiscmd)
    
    if rval != 0:
      print "error running metis with %s"%(metiscmd)
      print "continuing ..."
      continue
    
    partvecname = "%s.part.%i"%(graphfilename,np)
    
    # read in the partition vector
    if verbose:
      print "reading %s"%(partvecname)
      
    partmap = {}
    partfile = open(partvecname,'rt')
    for i,v in enumerate(partfile):
      curcluster = int(v)
      if curcluster in partmap:
        mycluster = partmap[curcluster]
      else:
        mycluster = nclusters
        partmap[curcluster] = nclusters
        nclusters += 1
      clusters.append((i,mycluster))
    partfile.close()
    if verbose:
      print "removing %s"%(partvecname)
    os.remove(partvecname)
    
      
  if verbose:
    print "writing %s"%(options.outfilename)
  # now write the hyper clustering
  hclustfile = open(options.outfilename,'wt')
  hclustfile.write('%i %i %i\n'%(nclusters, nverts, len(clusters)))
  for a in clusters:
    hclustfile.write('%i %i 1\n'%(a[1],a[0]))
  hclustfile.close()
  
  
  
if __name__=='__main__':
  main()
    
  
