function A=load_external_graph(graphname)
% LOAD__EXTERNAL_GRAPH Loads a graph from an external directory
%
% load_external_graph is a helper function to load a graph outside
% of the repository and normalize its contents. 
%
% Example:
%   A = load_external_graph('cond-m

% David F. Gleich
% Copyright, Purdue University, 2011

% History
% 2011-10-02: Initial coding based on load_gaimc_graph

A=readSMAT(fullfile('~/data/graph-db',[graphname '.smat']));
A=A|A';
A=spones(A); % todo, make this an option?
A = largest_component(A);
A = A - diag(diag(A)); % remove self-loops