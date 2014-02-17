hkgrowresults.mat contains 100 trials of (rand/heavy seed/hood) for each of the 12 small datasets.

hkgrowsmalldatasets.mat (once it's done) will contain the same, but computed using the new hkgrow_sresid_mex.

the script refinehkgrowdata.m simply computes averages of conductances and times over the 100 trials and prepares things to be plotted more easily.