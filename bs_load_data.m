function [spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval] = bs_load_data(params)
%bs_load_data. Load data in bs format
%
%
% spike_times = cell list with spike times for all clusters
% trial_stim_on = vector with stimulus onset times for all trials
% trial_stim_off = vector with stimulus offset times for all trials
% trial_stim_type = vector with stimulus types for all trials
% recording_interval contains 1x2 vector with start and end of stimulus set
%
% 2024, Alexander Heimel


switch params.dataset
    case 'Exp2019-12-16_MP3_S2L5_AP' % block 3 
        load(fullfile(params.datafolder,[params.dataset '.mat']),'sAP');

        % select V1 clusters
        sAP.sCluster = sAP.sCluster(contains(string({sAP.sCluster.Area}),"Primary visual area"));
        n_clusters = length(sAP.sCluster);

        spike_times = cell(n_clusters,1);
        for c = 1:n_clusters
            spike_times{c} = sAP.sCluster(c).SpikeTimes;
        end

        trial_stim_on = sAP.cellStim{params.block}.structEP.vecStimOnTime;
        trial_stim_off = sAP.cellStim{params.block}.structEP.vecStimOffTime;
        trial_stim_type = sAP.cellStim{params.block}.structEP.vecTrialStimTypes;
        n_trials = length(trial_stim_on);

        stimuli = unique(trial_stim_type);
        n_stimuli = length(stimuli);

        recording_interval(1) = trial_stim_on(1) - 1;
        recording_interval(2) = trial_stim_off(end) + 1;

        experiment_type = sAP.cellStim{params.block}.structEP.strFile;


    case 'Topo6_20220301_AP' % block 4
        load(fullfile(params.datafolder,[params.dataset '.mat']),'sAP')

        % select V1 clusters
       sAP.sCluster = sAP.sCluster(contains(string({sAP.sCluster.Area}),"Primary visual area"));

        % select responsive clusters
        sAP.sCluster = sAP.sCluster(arrayfun(@(x) x.ZetaP(2),sAP.sCluster)<0.05);

        n_clusters = length(sAP.sCluster);
        spike_times = cell(n_clusters,1);
        for c = 1:n_clusters
            spike_times{c} = sAP.sCluster(c).SpikeTimes;
        end

        trial_stim_on = sAP.cellBlock{params.block}.vecStimOnTime;
        trial_stim_off = sAP.cellBlock{params.block}.vecStimOffTime;
        trial_stim_type = sAP.cellBlock{params.block}.vecTrialStimTypes;

        % Changing stim type to orientation
        if params.only_distinguish_orientations
            trial_stim_type = mod(trial_stim_type-1,12)+1;
        end

        stimuli = unique(trial_stim_type);
        n_stimuli = length(stimuli);
        n_trials = length(trial_stim_on);

        recording_interval(1) = trial_stim_on(1) - 1;
        recording_interval(2) = trial_stim_off(end) + 1;

        experiment_type = sAP.cellBlock{params.block}.strExpType;

end
disp(['Loaded ' num2str(n_clusters) ' clusters and ' num2str(n_trials) ' trials for ' num2str(n_stimuli) ' different stimuli.'])
disp(['Stimulus type' experiment_type]);
