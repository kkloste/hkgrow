% first load a 'results' .mat file

% There are 15 total datsets:
%   12 small
%   3 large: ljournal, twitter, friendster
%

% for each graph, get a column vector, X, of runtimes.
%
%   For SMALL DATA
% we'll use the experiment randseed:
%   set experimenttype = 1;
%   use datax = times(id, trialnum, experimenttype)';

load newsmalldata;


numdata = 3+numel(filename);
percdata = zeros(numdata,3,2);
inputsize = zeros(numdata,1);
experimenttype = 1;
numgraphs = numel(filename);
inputsize(1:numgraphs) = gsize(1:numel(filename),1)+gsize(1:numel(filename),2);


functionid = 1; % for hk, 2 is for ppr
for id=1:numgraphs
    datax = times(id,:,experimenttype)';
    percdata(id,:,functionid) = prctile(datax,[25 50 75],1);
end


load newsmallppr;
functionid = 2; % for hk, 2 is for ppr
for id=1:numgraphs
    datax = times(id,:,experimenttype)';
    percdata(id,:,functionid) = prctile(datax,[25 50 75],1);
end

% now get LARGE DATA
% We have to do this one at a time, since each has
% its own .mat file

% friendstertrials
% twitterptrials
% ljournaltrials

% do HK first
functionid = 1;

load friendstertrials;
id = numgraphs+1;
inputsize(id) = gsize(:,1)+gsize(:,2);
datax = times(:,experimenttype);
percdata(id,:,functionid) = prctile(datax,[25 50 75],1);


load twitterptrials
id = numgraphs+2;
inputsize(id) = gsize(:,1)+gsize(:,2);
datax = times(:,experimenttype);
percdata(id,:,functionid) = prctile(datax,[25 50 75],1);


load ljournaltrials
id = numgraphs+3;
inputsize(id) = gsize(:,1)+gsize(:,2);
datax = times(:,experimenttype);
percdata(id,:,functionid) = prctile(datax,[25 50 75],1);


% now do PPR
functionid = 2;

load pprfriendstertrials;
id = numgraphs+1;
% inputsize(id) = gsize(:,1)+gsize(:,2);
datax = times(:,experimenttype);
percdata(id,:,functionid) = prctile(datax,[25 50 75],1);


load pprtwitterptrials
id = numgraphs+2;
% inputsize(id) = gsize(:,1)+gsize(:,2);
datax = times(:,experimenttype);
percdata(id,:,functionid) = prctile(datax,[25 50 75],1);


load pprljournaltrials
id = numgraphs+3;
% inputsize(id) = gsize(:,1)+gsize(:,2);
datax = times(:,experimenttype);
percdata(id,:,functionid) = prctile(datax,[25 50 75],1);



% HAVE ALL DATA

save(['/scratch2/dgleich/kyle/results/' 'runtimes' '.mat'], 'newindexing','percdata', 'inputsize','-v7.3');
