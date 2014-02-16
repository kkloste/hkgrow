hkgrow
======

Heat-kernel based local clustering

### Kyle Kloster
### David F. Gleich

---------

Main mex functions: hkgrow_sresid_mex.cpp
Given an input set of seed nodes s and sparse symmetric adjacency matrix A,
this will epsilon-approximate exp(t*A)*s (where t and epsilon are decided automatically)
in a degree-weighted infinity norm, then do a sweep over the resulting vector
and return the set of best conductance from that sweep.

hkgrow.m  calls this repeatedly to find a set of best conductance

---------

Experiments:

- waiting for trials on ljournal to finish; these trials will give a better idea how to choose t and eps optimally.
- converting the bigger datasets (twitter, friendster) into symmetrized matrices is troubelsome. The sets webbase and ljournal have been saved in a symmetrized form, so trials will start on those shortly.
- since the new hkgrow_sresid_mex.cpp is done, I'm rerunning past experiments with it (on the smaller datasets, since those were the only ones done)

