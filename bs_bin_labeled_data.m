function [bin_times,bin_counts,bin_labels,rel_bin_times,time_since_stim_off] = bs_bin_labeled_data(spike_times,label_onsets,label_offset,labels,recording_interval,params)
%bs_bin_labeled_data. Bins spike data with associated labels
%
%    [bin_times,bin_counts,bin_labels,rel_bin_times,time_since_stim_off] = ...
%            bs_bin_labeled_data(spike_times,label_onsets,label_offset,labels,recording_interval,params)
%
%     spike_times is cell list with spike times
%
% 2024, Alexander Heimel

n_clusters = length(spike_times);
n_trials = length(label_onsets);

n_bins = ceil( (recording_interval(2)-recording_interval(1)) / params.binsize );
bin_times = recording_interval(1) + params.binsize * (1:n_bins) - 0.5* params.binsize;

%% Bin spikes
bin_counts = zeros(n_clusters,n_bins);

for c = 1:n_clusters

    % Put spikes in bins
    for i = 1:length(spike_times{c})
        bin = ceil( (spike_times{c}(i)-recording_interval(1)) /params.binsize);
        if bin<1
            continue
        end
        if bin>n_bins
            break
        end
        bin_counts(c,bin) = bin_counts(c,bin) + 1;
    end

    % Transformation
    if params.square_root_transformation
        bin_counts(c,:) = bin_counts(c,:).^.5;
    end

    % Smoothing
    if params.smooth
        switch params.smooth_function
            case 'smooth'
                smooth_window_samples = round(params.smooth_window / params.binsize);
                bin_counts(c,:) = smooth(bin_counts(c,:),smooth_window_samples);
            case 'smoothen'
                smooth_window_samples = round(params.smooth_window / params.binsize);
                bin_counts(c,:) = smoothen(bin_counts(c,:),smooth_window_samples);
            otherwise
                logmsg(['Smoothing function ' params.smooth_function ' is unknown.'])
        end
    end
end

%% Creative relative bin times
rel_bin_times = params.max_time_since_stim_on + params.binsize * (1:n_bins);

for j = 1:n_trials
    ind = find(bin_times>label_onsets(j),1);
    dt = bin_times(ind) - label_onsets(j);
    rel_bin_times(ind:end) = rel_bin_times(ind:end) - rel_bin_times(ind) + dt;
end
rel_bin_times( rel_bin_times > params.max_time_since_stim_on) = params.max_time_since_stim_on;

time_since_stim_off = params.max_time_since_stim_off + params.binsize * (1:n_bins);
for j = 1:n_trials
    ind = find(bin_times>label_offset(j),1);
    dt = bin_times(ind) - label_offset(j);
    time_since_stim_off(ind:end) = time_since_stim_off(ind:end) - time_since_stim_off(ind) + dt;
end
time_since_stim_off( time_since_stim_off > params.max_time_since_stim_off) = params.max_time_since_stim_off;

%% Create bin labels vector
bin_labels = ones(1,n_bins) * params.empty_label;

for j = 1:n_trials
    sel = (bin_times>(label_onsets(j)+params.response_offset) & bin_times<(label_offset(j)+params.response_offset));
    bin_labels(sel) = labels(j);
end


