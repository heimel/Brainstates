function [bin_time,spike_counts,stim_type,time_since_stim_on,time_since_stim_off] = bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params)
%bs_bin_data. Bins stimulus and spike data
%
%     spike_times is cell list with spike times
%
% 2024, Alexander Heimel

n_clusters = length(spike_times);
n_trials = length(trial_stim_on);

n_bins = ceil( (recording_interval(2)-recording_interval(1)) / params.binsize );
bin_time = recording_interval(1) + params.binsize * (1:n_bins) - 0.5* params.binsize;


%% Bin spikes
spike_counts = zeros(n_clusters,n_bins);

for c = 1:n_clusters
    for i = 1:length(spike_times{c})
        bin = ceil( (spike_times{c}(i)-recording_interval(1)) /params.binsize);
        if bin<1
            continue
        end
        if bin>n_bins
            break
        end
        spike_counts(c,bin) = spike_counts(c,bin) + 1;
    end

    % smoothing
    if params.smooth
        spike_counts(c,:) = smooth(spike_counts(c,:),params.smooth_window);
        %spike_counts(c,:) = smoothen(spike_counts(c,:),params.smooth_window);
    end
end

%% Creative relative stim times
time_since_stim_on = params.max_time_since_stim_on + params.binsize * (1:n_bins);

for j = 1:n_trials
    ind = find(bin_time>trial_stim_on(j),1);
    dt = bin_time(ind) - trial_stim_on(j);
    time_since_stim_on(ind:end) = time_since_stim_on(ind:end) - time_since_stim_on(ind) + dt;
end
time_since_stim_on( time_since_stim_on > params.max_time_since_stim_on) = params.max_time_since_stim_on;

time_since_stim_off = params.max_time_since_stim_off + params.binsize * (1:n_bins);
for j = 1:n_trials
    ind = find(bin_time>trial_stim_off(j),1);
    dt = bin_time(ind) - trial_stim_off(j);
    time_since_stim_off(ind:end) = time_since_stim_off(ind:end) - time_since_stim_off(ind) + dt;
end
time_since_stim_off( time_since_stim_off > params.max_time_since_stim_off) = params.max_time_since_stim_off;

%% Create stimulus type vector
stim_type = -1*ones(1,n_bins);

for j = 1:n_trials
    sel = (bin_time>(trial_stim_on(j)+params.response_offset) & bin_time<(trial_stim_off(j)+params.response_offset));
    stim_type(sel) = trial_stim_type(j);
end


