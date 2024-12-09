%bs_check_datasets
%
%  temporary function to find suitable datasets in Neuropixels data
%
% 2024, Alexander Heimel

params = bs_default_params();

d = dir(fullfile(params.datafolder,'*.mat' ))
n_sets = length(d);

for i = 1 %:n_sets
    disp(d(i).name)
    load(fullfile(params.datafolder,d(i).name),'sAP');

    % select V1 clusters
    sAP.sCluster = sAP.sCluster(contains(string({sAP.sCluster.Area}),"Primary visual area"));
    n_clusters = length(sAP.sCluster)
end