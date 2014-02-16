function A=load_graph(graphname)
% LOAD_GRAPH Loads a graph from the data directory
%
% load_graph is a helper function to load a graph provided with the
% regardless of the current working directory.  
%
% Example:
%   A = load_graph('cond-mat-2005-fix-cc');

% David F. Gleich
% Copyright, Purdue University, 2011

% History
% 2011-10-02: Initial coding based on load_gaimc_graph

path=fileparts(mfilename('fullpath'));
A=readSMAT(fullfile(path,'data',[graphname '.smat']));
A=spones(A); % todo, make this an option?
A = A - diag(diag(A));