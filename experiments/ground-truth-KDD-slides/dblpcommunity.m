% Compare the recall of hk and ppr on ground-truth communities in amazon

% dataset = 'dblp';
% load(dataset,'/scratch/dlgeich/kyle/dblp/dblp');
load /scratch/dgleich/kyle/dblp/dblp;
addpath('/scratch2/dgleich/kyle/kdd/ppr');
n = size(A,1);
C(n,end) = 0;
Ctop(n,end) = 0;

%%

totalcommunities = 100;
avecommsize = 0;
bestfmeas = zeros(totalcommunities,2);
precofbestfmeas = zeros(totalcommunities,2);
recofbestfmeas = zeros(totalcommunities,2);
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
        avecommsize = avecommsize + deg;
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
                precofbestfmeas(numcom,functionID) = precisions(trial,functionID);
                recofbestfmeas(numcom,functionID) = recalls(trial,functionID);
                condofbestfmeas(numcom,functionID) = conds(trial,functionID);
            end
            [bset,conds(trial,2),cut,vol] = pprgrow(A,verts(trial));
            recalls(trial,2) = numel(intersect(verts,bset))/numel(verts);
            precisions(trial,2) = numel(intersect(verts,bset))/numel(bset);
            functionID = 2;
            fmeas(trial,functionID) = 2*recalls(trial,functionID)*precisions(trial,functionID)/(recalls(trial,functionID)+precisions(trial,functionID));
            if fmeas(trial,functionID) > bestfmeas(numcom,functionID),
                bestfmeas(numcom,functionID) = fmeas(trial,functionID);
                bestrecsize(numcom,functionID) = numel(bset);
                precofbestfmeas(numcom,functionID) = precisions(trial,functionID);
                recofbestfmeas(numcom,functionID) = recalls(trial,functionID);
                condofbestfmeas(numcom,functionID) = conds(trial,functionID);
            end
        end
        fprintf('Trial %i \t best hk = %8.4f  ssize=%i cond=%6.4f \t best ppr = %8.4f  ssize =%i cond=%6.4f \n',numcom,bestfmeas(numcom,1),bestrecsize(numcom,1), condofbestfmeas(numcom,1), bestfmeas(numcom,2), bestrecsize(numcom,2), condofbestfmeas(numcom,2));
    end
fprintf('hk: mean fmeas=%6.4f \t mean setsize=%6.4f \t mean prec=%6.4f \t ppr: mean fmeas=%6.4f \t mean setsize=%6.4f \t mean prec=%6.4f \t mean comm size=%6.4f\n', ...
		sum(bestfmeas(:,1))/totalcommunities, sum(bestrecsize(:,1))/totalcommunities, sum(precofbestfmeas(:,1))/totalcommunities, ...
		sum(bestfmeas(:,2))/totalcommunities, sum(bestrecsize(:,2))/totalcommunities, sum(precofbestfmeas(:,2))/totalcommunities, avecommsize);
avecommsize = avecommsize/totalcommunities;
save(['/scratch2/dgleich/kyle/kdd/' 'dblpcommunity' '.mat'],'bestfmeas','conds','avecommsize','recofbestfmeas','precofbestfmeas', 'bestrecsize', 'condofbestfmeas','bestrecsize','-v7.3');
exit;
