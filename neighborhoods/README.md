---
title: "Neighborhoods are good communities"
layout: project
---

Neighborhoods are good communities
==================================

### David F. Gleich
### C. Seshadhri

_These codes are research prototypes and may not work for you. No promises. But do email if you run into problems._


Download
--------

* [neigh-comm.tar.gz](neigh-comm.tar.gz) ___warning, will unzip into the current directory___
* [neigh-comm-data.tar.gz](neigh-comm-data.tar.gz) ___warning, 2gb file___
{: .nobullets}


Prereqs
-------

* [MatlabBGL](https://github.com/dgleich/matlab-bgl)
* [mcode](https://github.com/dgleich/mcode)
* A working C/C++ compiler
* A working Matlab mex compiler
* _Optional_ Metis

Setup
-----

Start matlab in the directory where you unzipped the neigh-comm.tar.gz file

    $ matlab
    >> setup_paths
    >> compile

This should work on Mac OSX (Lion tested) and Ubuntu linux (10.10 tested) with 
Matlab R2011a.

Once that's all done, the following should work (in Matlab)

    >> ncpneighs(load_graph('email-Enron'))

To run some of the codes, you'll need to compile the pprclusters directory.

In bash:

    $ cd pprclusters
    $ make

Please let me know if you run into any issues.
 
Overview
--------

The package is organized by directory

`/`  
: All of the main matlab codes

`data`
: data files for the experiments, and precomputed data

`test`
: codes to poke around the edges and make sure things work as expected

`experiments`
: experimental codes and figures

`pprclusters`
: standalone C++ code for the ppr clusters 

`web`
: this information and all the figures

Figures
-----------
    
|Experiment|Description|Figure|
|:------------------|:------------------------------------|:------------------|
|`experiments/lesmis/lesmis_figs.m` | The les mis figure sequences | Fig. 1, Fig. 6 |
|`data/print_graph_table` | Write out the data for each graph | Tab. 2 |
|`experiments/neighborhoods/neighborhoods.m` | Make the neighborhood plots from precomputed data | Fig. 2 |
|`experiments/neighborhoods/neighvsncp.m` | Make the neighborhood and pprncp plots from precomputed data | Fig. 3 |
|`experiments/cores/corefigs.m` | Compare the core clusters to the others | Fig. 4 |
|`experiments/grow/localmin_figure.m` | Show the locally minimal communities on the ncp| Fig. 5 |
|`experiments/grow/pprgrowfigs.m` | Compare the pprgrown clusters to the others | Fig. 7 |

Online only figures
-------------------

* [Vertex neighborhoods](figures.html#neigh)
* [Vertex neighborhoods and PageRank NCPs](figures.html#nvsncp)
* [Core neighborhoods](figures.html#core)
* [PageRank grown locally minimal neighborhoods communities](figures.html#pprgrow)

The names of the graphs are:

    Penn94, anony-interactions-oneyearA-cc (fb-A-onyear), arxiv-ubc,
    as-22july06, ca-AstroPh-cc, cond-mat-2005-fix-cc,
    dblp-cc, email-Enron-cc, hollywood-2009-cc,
    itdk0304-cc, networkA-anonymized (fb-A), oregon2_010526,
    p2p-Gnutella25, rand-ff-25000-0.4 (ff-0.4), rand-ff-25000-0.49 (ff-0.49),
    soc-LiveJournal1, web-Google

Any figure can be accessed via 
* [neigh-itdk0304-cc.png](neigh-itdk0304-cc.png)
* [nvsncp-itdk0304-cc.png](nvsncp-itdk0304-cc.png)
* [core-itdk0304-cc.png](core-itdk0304-cc.png)
* [pprgrow-itdk0304-cc.png](pprgrow-itdk0304-cc.png)

The other graphs are:

    Penn94, anony-interactions-oneyearA-cc (fb-A-onyear), arxiv-ubc,
    as-22july06, ca-AstroPh-cc, cond-mat-2005-fix-cc,
    dblp-cc, email-Enron-cc, hollywood-2009-cc,
    itdk0304-cc, networkA-anonymized (fb-A), oregon2_010526,
    p2p-Gnutella25, rand-ff-25000-0.4 (ff-0.4), rand-ff-25000-0.49 (ff-0.49),
    soc-LiveJournal1, web-Google
    
So put the general pattern is `<figtype>-<graph>.png`



