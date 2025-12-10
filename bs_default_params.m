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

params.empty_label = -1; % label for period without a stimulus

params.reproducible = true; % uses fixed random seed

rootpath = fileparts(which('brainstates.mlx'));
params.datafolder = fullfile(rootpath,'Data');

% Neuropixels dataset with V1 recordings
params.dataset = 'Topo6_20220301_AP';
params.block = 2; %

% params.dataset = 'Exp2019-12-16_MP3_S2L5_AP';
% params.block = 2; %

params.binsize = 0.01; % s
params.square_root_transformation = true; % See Yu et al. J Neurophysiol 2009
params.smooth = true;
params.smooth_function = 'smooth'; % {'smooth' = moving window, 'smoothen' = gaussian
params.smooth_window = 0.1; % s

params.max_time_since_stim_on = 1.5; % s 
params.max_time_since_stim_off = 0.5; % s

params.response_offset = 0.050; % s, delay of visual response

params.only_distinguish_orientations = true; % make two opposing directions the same stimulus type

params.cluster_tsne = false;
params.cluster_method = 'kmeans';

params.hmm_version = 'Mathworks'; % alternative Probml
params.hmm_algorithm = 'BaumWelch'; % BaumWelch is more thorough? 
%params.hmm_algorithm = 'Viterbi'; % Viterbi is quicker?

params.jitter_hint = 0.1;

% Load processparams_local. Keep at the end
if exist('processparams_local.m','file')
    params = processparams_local( params );
end

