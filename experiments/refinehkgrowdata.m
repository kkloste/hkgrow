%% Refining data from hkgrowexperiment
% assumes a .mat file from "experimenthkgrow.m" has be loaded already


% get averages over the 100 trials
avtimes = zeros(numel(filename),4);
avconds = zeros(numel(filename),4);

for fileid=1:numel(filename)
    for etype=1:4
        avtimes(fileid,etype) = sum(times(fileid,:,etype))./length(times(fileid,:,etype));
        avconds(fileid,etype) = sum(conds(fileid,:,etype))./length(conds(fileid,:,etype));
    end
end
