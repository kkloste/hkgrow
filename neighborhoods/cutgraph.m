function Acut = cutgraph(A,set)
% CUTGRAPH Return the graph of only the edges that were cut
%
% Acut = cutgraph(A,set) only has edges that are cut by the set.
%

Acut = A;
Acut(set,set) = 0;
if numel(set) ~= size(A,1)
    % then this is a list of vertices
    otherset = setdiff(1:size(A,1),set);
else
    otherset = ~set;
end
Acut(otherset,otherset) = 0;


