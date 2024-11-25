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


%params.datafolder = '\\vs03\VS03-CSF-1\Heimel\Heimel_HMM\Data_collection\Montijn_elife_2021';
params.datafolder = 'C:\Users\heimel.HERSENINSTITUUT\OneDrive\Desktop\Brainstates';
params.binsize = 0.1; % s

params.max_time_since_stim_on = 3; 
params.max_time_since_stim_off = 3; 

params.smooth = true;
params.smooth_window = 4; % samples

% Load processparams_local. Keep at the end
if exist('processparams_local.m','file')
    params = processparams_local( params );
end

