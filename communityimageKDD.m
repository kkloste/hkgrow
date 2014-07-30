% Find ground-truth communities of DBLP roughly 100 in size
% compute fmeas, set size of HK, PPR.

% dataset = 'dblp';
% load(dataset,'/scratch/dlgeich/kyle/dblp/dblp');
load /scratch/dgleich/kyle/dblp/dblp;
addpath('/scratch2/dgleich/kyle/hkgrow/ppr');
addpath('/scratch2/dgleich/kyle/hkgrow');

n = size(A,1);
C(n,end) = 0;

% find communities of size 80 < size < 400
n = size(A,1);
e = ones(n,1);
commsize = e'*C;
comminds = find(commsize>50);
dummy = find(commsize(comminds)<120);
comminds = comminds(dummy);
disp(length(comminds))

% now comminds contains indices of C corresponding to communities
% with size between 80 and 140.

%% Next, run PPR and HK on all communities found this way.
% Then we find an example that is representative and easily visualized.



totalcommunities = length(comminds);
bestfmeas = zeros(totalcommunities,3);   % ratio of fmeas at each starting node
bestrecsize = zeros(totalcommunities,2); % size of best-set returned by each alg
commsizes = zeros(totalcommunities,2);   % size of actual community
bestprecs = zeros(totalcommunities,2);

for numcom=1:totalcommunities
    comm = comminds(numcom);
    verts = find(C(:,comm));
    commsizes(numcom) = length(verts);
 
    deg = numel(verts);
    recalls = zeros(deg,2); % hk = 1, ppr = 2
    precisions = zeros(deg,2);
    fmeas = zeros(deg,2);
    conds = zeros(deg,2);

    for trial = 1:deg
        % get HK
        [bset,conds(trial,1),cut,vol,~,~] = hkgrow1(A,verts(trial),'t',5);
        recalls(trial,1) = numel(intersect(verts,bset))/numel(verts);
        precisions(trial,1) = numel(intersect(verts,bset))/numel(bset);
        functionID = 1; hksize = numel(bset); hkprec = precisions(trial,1);
        fmeas(trial,functionID) = 2*recalls(trial,functionID)*precisions(trial,functionID)/(recalls(trial,functionID)+precisions(trial,functionID));

        % get PPR
        [bset,conds(trial,2),cut,vol] = pprgrowKDD(A,verts(trial));
        recalls(trial,2) = numel(intersect(verts,bset))/numel(verts);
        precisions(trial,2) = numel(intersect(verts,bset))/numel(bset);
        functionID = 2; prsize = numel(bset); prprec = precisions(trial,2);
        fmeas(trial,functionID) = 2*recalls(trial,functionID)*precisions(trial,functionID)/(recalls(trial,functionID)+precisions(trial,functionID));
        
        if fmeas(trial,1) > bestfmeas(numcom,1),
            if recalls(trial,1) > 0.2,
                bestfmeas(numcom,1) = fmeas(trial,1);
                bestfmeas(numcom,2) = fmeas(trial,2);
                bestfmeas(numcom,3) = verts(trial);
                bestprecs(numcom,1) = hkprec;
                bestprecs(numcom,2) = prprec;
                bestrecsize(numcom,1) = hksize;
                bestrecsize(numcom,2) = prsize;
            end
        end

    end
    if bestfmeas(numcom,1) > 0,
        fprintf('CommSize=%i  CommID=%i \t HKfmeas=%8.4f  PRfmeas=%8.4f \t HKsetsize=%i  PRsetsize=%i \t HKprec=%8.4f  PRprec=%8.4f \t seedID=%i \n', ...
             commsizes(numcom),comm,bestfmeas(numcom,1),bestfmeas(numcom,2),bestrecsize(numcom,1),bestrecsize(numcom,2), ...
             bestprecs(numcom,1), bestprecs(numcom,2), bestfmeas(numcom,3));
    end
end


save(['/scratch2/dgleich/kyle/kdd/' 'communityimage' '.mat'], 'bestrecsize', 'bestfmeas','bestrecsize','bestprecs','comminds','-v7.3');