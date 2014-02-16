function dataout = ncpneighs(A,varargin)
% NCPNEIGHS A handy tool for showing NCPs based on neighborhood communities

p = inputParser;
p.addOptional('fiedler',true,@islogical);
p.addOptional('greedy',true,@islogical);
p.addOptional('whisker',true,@islogical);
p.addOptional('pprgrow',true,@islogical);
p.addOptional('pprgrow_mindegree',5,@isnumeric);
p.addOptional('pprgrowcores',true,@islogical);
p.addOptional('corecuts',true,@islogical);
p.parse(varargin{:});


n = size(A,1);
fprintf('Computing neighborhood clusters ... ');
data.d = full(sum(A,2));
d = data.d;
t0 = tic;
[data.cond data.cut data.vol data.cc data.t] = triangleclusters(A);
data.time.neighborhood = toc(t0);
fprintf('%.1f sec\n', data.time.neighborhood);


loglog(d+1,data.cond,'.');
bestcond = accumarray(d+1,data.cond,[],@min);
data.bestcond = bestcond;

holdstate = ishold;
if ~holdstate; hold on; end
bestcondi = find(bestcond);
loglog(bestcondi,bestcond(bestcondi),'k-','LineWidth',2);

xlabel('Cluster size in vertices');
ylabel('Conductance');

if p.Results.fiedler
    fprintf('Computing Fiedler vector ... ');
    t0 = tic;
    [x,forder,cuts,bestset,lam2] = nfiedler(A);
    data.time.fiedler = toc(t0);
    fprintf('%.1f sec\n', data.time.fiedler);
    data.fiedler = x;
    data.fiedler_set = bestset;
    data.fiedler_cuts = cuts;
    bestcond = min(cuts.conductance);
    data.fiedler_cond = bestcond;
    line([1,max(d+1)],[bestcond bestcond],'Color',[0.8,0,0]);
    line([1,max(d+1)],[lam2 lam2],'Color',[0.8,0,0]);
    plot(numel(bestset),bestcond,'ro');
end

if p.Results.greedy
    fprintf('Computing neighborhood + greedy clusters ... ');
    t0 = tic;
    [condg,cutg,volg,sg] = triangleclustersgreedy_mex(A);
    data.time.greedy = toc(t0);
    fprintf('%.1f sec\n', data.time.greedy);
    best = accumarray(sg,condg,[],@min);
    besti = find(best);
    loglog(besti,best(besti),'g-','LineWidth',2);
    data.greedy_cond = condg;
    data.greedy_cut = cutg;
    data.greedy_size = sg;
    data.greedy_vol = volg;
end

if p.Results.whisker
    fprintf('Computing whisker clustering ... ');
    t0 = tic;
    [condw,cutw,volw,sw] = bestwhiskercuts(A);
    data.time.whisker = toc(t0);
    fprintf('%.1f sec\n', data.time.whisker);
    loglog(sw,condw,'m.-');
    data.whisker_cond = condw;
    data.whisker_size = sw;
    data.whisker_cut = cutw;
    data.whisker_vol = volw;
end

if p.Results.corecuts
    fprintf('Computing corecuts clustering ... ');
    t0 = tic;
    ccuts = corecuts(A);
    data.time.corecuts = toc(t0);
    fprintf('%.1f sec\n', data.time.corecuts);
    data.core_cond = ccuts.cond;
    data.core_size = ccuts.size;
    data.core_cut = ccuts.cut;
    data.core_vol = ccuts.vol;
    data.core_set = ccuts.bestset;
    loglog(data.core_size, data.core_cond,'c*-');
    
    if p.Results.pprgrowcores
        nv = numel(ccuts.cond);
        fprintf('Computing pprgrow clustering ... for %i coresets ... \n', nv);
        t0 = tic;
        data.pprgrowcores_cond = zeros(nv,1);
        data.pprgrowcores_cut = zeros(nv,1);
        data.pprgrowcores_vol = zeros(nv,1);
        data.pprgrowcores_size = zeros(nv,1);
        for i=1:nv
            nset = ccuts.cuts.order(1:ccuts.ind(i));
            if numel(nset) > n/2
                nset = setdiff(1:n,nset);
            end
            [curset,data.pprgrowcores_cond(i),data.pprgrowcores_cut(i),data.pprgrowcores_vol(i)] = ...
                pprgrow(A,nset,'nruns',4,'maxexpand',nnz(A)/3);
            data.pprgrowcores_size(i) = min(numel(curset),n-numel(curset));
        end
        data.time.pprgrowcores = toc(t0);
        fprintf('%.1f sec\n', data.time.pprgrowcores);
        best = accumarray(data.pprgrowcores_size,data.pprgrowcores_cond,[],@min);
        besti = find(best);
        loglog(besti,best(besti),'co-');
    end
    
%      if p.Results.pprgrowcore
%         nset = data.core_set;
%         if numel(nset)>n/2,
%             nset = setdiff(1:n,nset);
%         end
%         [coregrowset,coregrowcond,coregrowcut,coregrowvol] = ...
%             pprgrow(A,nset,'nruns',12,'maxexpand',nnz(A)/3);
%         plot(min(numel(coregrowset),n-numel(coregrowset)),coregrowcond,'co',...
%             'MarkerSize',18);
%     end
end

if p.Results.pprgrow
    t0 = tic;
    minverts = neighborhoodmin(A,data.cond);
    minverts = minverts(data.d(minverts)>=p.Results.pprgrow_mindegree);
    data.pprgrow_verts=minverts;
    nv = numel(minverts);
    fprintf('Computing pprgrow clustering ... for %i verts\n', nv);
    data.pprgrow_cond = zeros(nv,1);
    data.pprgrow_cut = zeros(nv,1);
    data.pprgrow_vol = zeros(nv,1);
    data.pprgrow_size = zeros(nv,1);
    for i=1:numel(minverts)
        nset = [minverts(i); find(A(:,minverts(i)))];
        [curset,data.pprgrow_cond(i),data.pprgrow_cut(i),data.pprgrow_vol(i)] = ...
            pprgrow(A,nset,'nruns',12,'maxexpand',nnz(A)/3);
        data.pprgrow_size(i) = min(numel(curset),n-numel(curset));
    end
    data.time.pprgrow = toc(t0);
    fprintf('%.1f sec\n', data.time.pprgrow);
    best = accumarray(data.pprgrow_size,data.pprgrow_cond,[],@min);
    besti = find(best);
    loglog(besti,best(besti),'y.-','LineWidth',2);
end




if ~holdstate; hold off; end
    
if nargout>0
    dataout = data;
end