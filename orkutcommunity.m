% Compare the recall of hk and ppr on ground-truth communities in orkut

% dataset = 'youtube';
% load(dataset,'/scratch/yu163/youtube');
load /scratch/yu163/orkut/orkut;
addpath('/scratch2/dgleich/kyle/kdd/ppr');
n = size(A,1);
C(n,end) = 0;
Ctop(n,end) = 0;

%%

totalcommunities = 100;
bestfmeas = zeros(totalcommunities,2);
bestrecsize = zeros(totalcommunities,2);
condofbestfmeas = zeros(totalcommunities,2);
% find the first community with size > 10
n = size(A,1);
e = ones(n,1);
commsize = e'*C;
comm1 = min(find(commsize>10));

% check every 10th community after the first one
% that has size > 10
testcomms = zeros(totalcommunities,1);
for i=1:totalcommunities
testcomms(i) = comm1 + 10*(i-1);
end

for numcom=1:totalcommunities
comm = testcomms(numcom);
verts = find(C(:,comm));

deg = numel(verts);
recalls = zeros(deg,2); % hk = 1, ppr = 2
precisions = zeros(deg,2);
fmeas = zeros(deg,2);
conds = zeros(deg,2);

for trial = 1:deg
[bset,conds(trial,1),cut,vol,~,~] = hkgrow1(A,verts(trial),'t',5);
recalls(trial,1) = numel(intersect(verts,bset))/numel(verts);
precisions(trial,1) = numel(intersect(verts,bset))/numel(bset);
functionID = 1;
fmeas(trial,functionID) = 2*recalls(trial,functionID)*precisions(trial,functionID)/(recalls(trial,functionID)+precisions(trial,functionID));
if fmeas(trial,functionID) > bestfmeas(numcom,functionID),
bestfmeas(numcom,functionID) = fmeas(trial,functionID);
bestrecsize(numcom,functionID) = numel(bset);
condofbestfmeas(numcom,functionID) = conds(trial,1);
end
[bset,conds(trial,2),cut,vol] = pprgrow(A,verts(trial));
recalls(trial,2) = numel(intersect(verts,bset))/numel(verts);
precisions(trial,2) = numel(intersect(verts,bset))/numel(bset);
functionID = 2;
fmeas(trial,functionID) = 2*recalls(trial,functionID)*precisions(trial,functionID)/(recalls(trial,functionID)+precisions(trial,functionID));
if fmeas(trial,functionID) > bestfmeas(numcom,functionID),
bestfmeas(numcom,functionID) = fmeas(trial,functionID);
bestrecsize(numcom,functionID) = numel(bset);
condofbestfmeas(numcom,functionID) = conds(trial,1);
end
end
fprintf('best hk = %8.4f  setsize=%i cond=%6.4f \t best ppr = %8.4f  setsize =%i cond=%6.4f \n',bestfmeas(numcom,1),bestrecsize(numcom,1), condofbestfmeas(numcom,1), bestfmeas(numcom,2), bestrecsize(numcom,2), condofbestfmeas(numcom,2));
end
fprintf('mean for hk = %f \t mean for ppr = %f \n', sum(bestfmeas(:,1))/totalcommunities, sum(bestfmeas(:,2))/totalcommunities);

save(['/scratch2/dgleich/kyle/kdd/' 'orkutcommunity' '.mat'],'fmeas','conds','recalls','precisions','condofbestfmeas','bestrecsize','-v7.3');