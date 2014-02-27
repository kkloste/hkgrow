---
title: "Heat kernel diffusion for local clustering"
layout: project
---

hkgrow: heat-kernel based local clustering
==========================================

### Kyle Kloster
### David F. Gleich

_These are research codes and may not work for you._

Download
--------

* [TO DO](filename) (date)

Synopsis
--------

    compile % compile the mex files
    G = load_graph('dolphins');
    hkgrow_mex(G,seeds,t,eps);
    
Reusable codes
--------------

* `hkgrow_mex.m` C++ MEX code for computing a set of best conductance via
  seeded heat kernel with input graph "A" and input parameters "seeds, eps, t".


Codes from others
-----------------

* `pprgrow.m` from personal page rank clustering

Results from the paper [TO DO: replace with KDD results]
----------------------

To reproduce figure 2 (left), run:

    test_tol_accuracy % generate the data
    plot_tol_accuracy % plot the data
