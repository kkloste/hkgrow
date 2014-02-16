function [hpts,hline] = ncpplot(cond,csize)


hpts=loglog(csize,cond,'.');
bestcond = accumarray(csize,cond,[],@min);

holdstate = ishold;
if ~holdstate; hold on; end
bestcondi = find(bestcond);
hline=loglog(bestcondi,bestcond(bestcondi),'k-','LineWidth',2);

xlabel('Cluster size in vertices');
ylabel('Conductance');

if ~holdstate; hold off; end
    
