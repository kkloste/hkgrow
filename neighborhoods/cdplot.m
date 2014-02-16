function cdplot(cond,s)
% CDPLOT Community density plot
%


lcond = log10(cond);
ls = log10(s);
jcond = 10.^(lcond+0.1*randn(size(lcond)));
js = 10.^(ls+0.1*randn(size(lcond)));


opts.alpha = 0.2;
opts.plotcolor = [0 0 0];

switch get(gca,'NextPlot')
    case 'replace'
        cla reset;
    case 'replacechildren'
        cla;
    case 'add'
end

set(gcf,'renderer','openGL')

x = log10(js);
y = log10(jcond);

x = [x(:), NaN*ones(numel(x),1)]'; 
y = [y(:), NaN*ones(numel(y),1)]'; 
h = patch(x(:),y(:),...
            [0,0,0],'EdgeAlpha',opts.alpha,'EdgeColor',opts.plotcolor,...
            'FaceColor',opts.plotcolor);
