function params = bs_default_params(record)
%vbs_default_paramters Contains default parameters for brainstate analysis
%
%  PARAM = bs_default_parameters([RECORD])
%
%  Edit processparams_local for local and temporary edits to these default
%  parameters.
%
% 2024, Alexander Heimel 

if nargin<1 
    record = [];
end

% General
params = struct;

params.reproducible = true; % uses fixed random seed

params.datafolder = 'C:\Users\heimel.HERSENINSTITUUT\OneDrive\Projects\Heimel\Brainstates';
params.binsize = 0.2; % s

params.max_time_since_stim_on = 3; 
params.max_time_since_stim_off = 3; 

params.smooth = false;
params.smooth_window = 2; % samples

params.response_offset = 0.050; % s, delay of visual response

params.only_distinguish_orientations = true; % make two opposing directions the same stimulus type

params.cluster_tsne = false;
params.cluster_method = 'kmeans';

% Load processparams_local. Keep at the end
if exist('processparams_local.m','file')
    params = processparams_local( params );
end

