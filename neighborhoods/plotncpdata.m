function plotncpneighs(data)
% PLOTNCPNEIGHS A handy tool for showing NCPs based on neighborhood communities

loglog(data.d+1,data.cond,'.');

holdstate = ishold;
if ~holdstate; hold on; end
bestcondi = find(data.bestcond);
loglog(bestcondi,data.bestcond(bestcondi),'k-','LineWidth',2);

xlabel('Cluster size in vertices');
ylabel('Conductance');

if isfield(data,'fiedler_set')
    line([1,max(data.d+1)],[data.fiedler_cond data.fiedler_cond],'Color',[0.8,0,0]);
    plot(numel(data.fiedler_set),data.fiedler_cond,'ro');
end

if isfield(data,'greedy_cond')
    sg = data.greedy_size;
    condg = data.greedy_cond;
    best = accumarray(sg,condg,[],@min);
    besti = find(best);
    loglog(besti,best(besti),'g-','LineWidth',2);
end
if ~holdstate; hold off; end
    
