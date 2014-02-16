function [cond,cut,vol,s,v] = whiskercuts(A)
% WHISKERCUTS Build a set of cuts from one-connected whiskers
%
% [cond,cut,vol,s,f] = whiskercuts(A) returns the evaluation of
% a set of whisker cuts. Each whiskercut is a biconnected component
% connected to the largest biconnected component by a single edge.
% Consequently, for each component, we increase the value of the 
% cut by one, but increase the volume/size by a larger fraction.
% This simplicity lets us build up a set of good conductance
% clusters greedily.
%
% v is an indicator over vertices giving the whisker id of
% each vertex
%
% Example:
%   A = load_graph('ca-AstroPh');
%   [cond,cut,vol,s] = whiskercuts(A);
%   ncpplot(cond,s);

[Abc,p,Af,fcc] = biconncore(A);

% A whisker is any component of Af with size at least two
fccsizes = accumarray(fcc,1);
maxcc = max(fcc);
whiskerind = zeros(maxcc,1);
whiskerind(fccsizes>1)=1; % whiskerind is now 1 for all whiskers
[~,bigcc] = max(fccsizes);
whiskerind(bigcc) = 0;    % a whisker isn't the largest CC
whiskerverts = whiskerind(fcc); 
whiskerverts(whiskerverts > 0) = fcc(whiskerverts>0);
% find the number of edges in each whisker
[i,~] = find(Af);
% whiskerverts maps vertices to whisker component ids
fccedges = accumarray(nonzeros(whiskerverts(i)),1);
% fccedges gives the number of edges in each cc of the whisker graph

Gvol = full(sum(sum(A)));

% form all the whisker cuts
[wi,~,vol] = find(fccedges);
vol = vol+1; % add one for the missing edge
cut = ones(length(wi),1);
s = fccsizes(wi);
cond = cut./min(vol,Gvol-vol);

whiskerind(wi) = 1:length(wi);
v = whiskerverts;
v(v>0) = whiskerind(nonzeros(whiskerverts));

