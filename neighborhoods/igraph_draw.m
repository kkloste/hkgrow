function [xy,myimfile] = igraph_draw(A,layout)



mydir = fileparts(mfilename('fullpath'));
% parse options
%options = struct(varargin{:});
if ~isempty(getenv('PYTHON'))
    python=getenv('PYTHON');
else
    python='python';
end
mycmd = [python ' ' fullfile(mydir,'igraph_draw.py')];
mygraphfile = fullfile(mydir,'igraph_draw-graph.smat');
myimfile = fullfile(mydir,'igraph_draw-graph.png');
myxyfile = fullfile(mydir,'igraph_draw-graph.xy');
if nargin==1
    layout = 'lgl';
end


mycmd = [mycmd ' ' mygraphfile ' ' myimfile ' ' myxyfile ' ' layout];

writeSMAT(mygraphfile,A);
status = system(mycmd);

if status==0
    xy=load(myxyfile);
end