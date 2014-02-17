---
title: "KDD Paper: Comparing Heatkernel and pagerank in local clustering"
layout: project
---

Heatkernel vs pagerank in Local Clustering
==================================

### Kyle Kloster
### David F. Gleich

_These codes are research prototypes and may not work for you. No promises. But do email if you run into problems._


Overview
--------

Functions for heat kernel / pagerank
------------------------------------

gsqppr_mex - uses the Gauss Southwell approach to compute the page rank vector, using a queue instead of a heap. This is our nexpokit function rewritten to the resolvent times a "seed vector" (i.e. a vector with 1s in the rows corresponding to seed nodes for clustering).

gsqexpmv_mex - the nexpokit function (the queue version) for computing columns of expm(P), adapted to compute expm(P)

hppr_mex - "hp" functions use "heap matvecs", rather than Gauss Southwell. "pr" for page rank.

hpexpmv_mex - "hp" functions use "heap matvecs" (not Gauss Southwell). Computes expm(P)*v, where v is a "seed vector".

Clustering Functions
--------------------

hkgrow.m - this is the "pprgrow.m" function from David's "Neighborhoods are good communities" work, except I replace the page rank computation with a heat kernel computation -- my "hpexpmv" version which uses a heap.


