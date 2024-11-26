
params = bs_default_params();


%% Load data

% Required:
% trial_stim_on = vector with stimulus onset times for all trials
% trial_stim_off = vector with stimulus offset times for all trials
% trial_stim_type = vector with stimulus types for all trials
% n_trials = length(trial_stim_on);
% spike_times = cell list with spike times for all clusters
% n_clusters = length(spike_times);

params.dataset = 'Topo6_20220301_AP';
switch params.dataset
    case 'Exp2019-12-16_MP3_S2L5_AP'
        load(fullfile(params.datafolder,[params.dataset '.mat']))

        % select V1 clusters
        sAP.sCluster = sAP.sCluster(contains(string({sAP.sCluster.Area}),"Primary visual area"));
        n_clusters = length(sAP.sCluster);

        spike_times = cell(n_clusters,1);
        for c = 1:n_clusters
            spike_times{c} = sAP.sCluster(c).SpikeTimes;
        end

        block = 3;
        trial_stim_on = sAP.cellStim{block}.structEP.vecStimOnTime;
        trial_stim_off = sAP.cellStim{block}.structEP.vecStimOffTime;
        trial_stim_type = sAP.cellStim{block}.structEP.vecTrialStimTypes;
        n_trials = length(trial_stim_on);

        stimuli = unique(trial_stim_type);
        n_stimuli = length(stimuli);

        recording_interval(1) = trial_stim_on(1) - 1;
        recording_interval(2) = trial_stim_off(end) + 1;

    case 'Topo6_20220301_AP'
        load(fullfile(params.datafolder,[params.dataset '.mat']))
        block = 4;

        % select V1 clusters
        sAP.sCluster = sAP.sCluster(contains(string({sAP.sCluster.Area}),"Primary visual area"));

        % select responsive clusters
        sAP.sCluster = sAP.sCluster(arrayfun(@(x) x.ZetaP(2),sAP.sCluster)<0.05);

        n_clusters = length(sAP.sCluster);
        spike_times = cell(n_clusters,1);
        for c = 1:n_clusters
            spike_times{c} = sAP.sCluster(c).SpikeTimes;
        end

        trial_stim_on = sAP.cellBlock{block}.vecStimOnTime;
        trial_stim_off = sAP.cellBlock{block}.vecStimOffTime;
        trial_stim_type = sAP.cellBlock{block}.vecTrialStimTypes;

        % Changing stim type to orientation
        if params.only_distinguish_orientations
            trial_stim_type = mod(trial_stim_type-1,12)+1;
        end

        stimuli = unique(trial_stim_type);
        n_stimuli = length(stimuli);
        n_trials = length(trial_stim_on);

        recording_interval(1) = trial_stim_on(1) - 1;
        recording_interval(2) = trial_stim_off(end) + 1;

        recording_interval(2) = 3300;
end

disp(['Loaded ' num2str(n_clusters) ' clusters and ' num2str(n_trials) ' trials for ' num2str(n_stimuli) ' different stimuli.'])

[bin_times,spike_counts,stim_types,times_since_stim_on,time_since_stim_off] = ...
    bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
n_bins = length(bin_times);

%%
ind_bin_subset = round(linspace(1,n_bins,min([5000 n_bins]))); % total runs show change over time in t-sne plot


if 0
    %% Get correlated pairs
    no_correlation = true;
    c1 = 0;

    while c1<n_clusters && no_correlation
        c1 = c1 + 1;
        c2 = c1;
        while c2<n_clusters && no_correlation
            c2 = c2+1;
            x = corrcoef(spike_counts(c1,:),spike_counts(c2,:));
            if x(1,2)>0.5
                disp(['Good correlation for ' num2str(c1) ' and ' num2str(c2)]);
                no_correlation = false;
            end
        end
    end


    %% PCA Example with two neurons
    %c1 = 22; c2 = 27;
    %c1 = 7; c2 = 45;
    %c1 = 38; c2 = 39;
    %c1 = 44; c2 = 45;
    figure
    scatter(spike_counts(c1,ind_bin_subset) + 0.5*(rand(1,length(ind_bin_subset))-0.5),...
        spike_counts(c2,ind_bin_subset) + 0.5*(rand(1,length(ind_bin_subset))-0.5),'filled')
    xlabel('Spike count neuron 1 (jittered)');
    ylabel('Spike count neuron 2 (jittered)');
    hold on
    axis equal

    x = spike_counts([c1 c2],:)';
    [coeff,score,latent] = pca(x);

    mx = mean(x);
    plot(mx(1),mx(2),'or','MarkerFaceColor','r')
    plot(mx(1)+[0 coeff(1,1)],mx(2)+[0 coeff(2,1)],'r-','LineWidth',2)
    plot(mx(1)+[0 coeff(1,2)],mx(2)+[0 coeff(2,2)],'g-','LineWidth',2 )


    x_reconstruct = score*coeff' + repmat(mx,n_bins,1);

    scatter(x_reconstruct(ind_bin_subset,1),x_reconstruct(ind_bin_subset,2),'filled')

end

%% PCA for all clusters
disp('Computing PCA')
[coeff,score,latent] = pca(spike_counts');
mx = mean(spike_counts');


figure
subplot(2,2,1)
bar(latent)
ylabel('Explained (%) ')
xlabel('Component')

subplot(2,2,2)
imagesc(spike_counts)
ylabel('Neuron')
xlabel('Sample')
title('Origina data')

subplot(2,2,3)
hold on
plot(bin_times,score(:,1))
plot(bin_times,score(:,2))
plot(bin_times,score(:,3))
plot(bin_times,times_since_stim_on)
xlabel('Time (s)')
ylabel('Score')


subplot(2,2,4)
%x = score(:,1)*coeff(:,1)' + repmat(mx,n_bins,1);
x = score(:,(1:3))*coeff(:,(1:3))' + repmat(mx,n_bins,1);
% x_norm = x ./ repmat(max(x),n_bins,1);
% imagesc(x_norm')
imagesc(x')
%imagesc(x')
ylabel('Neuron')
xlabel('Sample')
title('First 3 PCs')

%% PC1 vs PC2 with color time since stim on
figure;
scatter(score(ind_bin_subset,1),score(ind_bin_subset,2),30,times_since_stim_on(ind_bin_subset),'filled')
axis square
colorbar
xlabel('PC1')
ylabel('PC2')

%% PC1 vs PC2 with color stim_type
figure;
scatter(score(ind_bin_subset,1),score(ind_bin_subset,2),30,stim_types(ind_bin_subset),'filled')
colormap(parula(max(stim_types(ind_bin_subset))+2))

axis square
colorbar
xlabel('PC1')
ylabel('PC2')


%% t-sne for subset
disp('Computing t-sne')
y = tsne(spike_counts(:,ind_bin_subset)');

% Show t-sne with time as color
figure
h = subplot(2,2,1);
colormap(h,'parula')
scatter(y(:,1),y(:,2),30,bin_times(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colorbar
title('Absolute time')

% Show t-sne with time since stim onset as color
h = subplot(2,2,2);
scatter(y(:,1),y(:,2),30,times_since_stim_on(ind_bin_subset),'filled')
%scatter(y(:,1),y(:,2),30,time_since_stim_off(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,periodic_colormap(128))
colorbar
title('Time since stim on')

% Show t-sne with stim type as color
h = subplot(2,2,3);

scatter(y(:,1),y(:,2),30,stim_types(ind_bin_subset)+2,'filled')
%scatter(y(:,1),y(:,2),30,time_since_stim_off(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,parula(max(stim_types(ind_bin_subset))+2))
colorbar
title('Stim type')


%% Optimal HMM fit to two states
n_emissions = 20;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'states_on_off';
params.binsize = 0.01;
[bin_times,spike_counts,stim_types,times_since_stim_on,~] = bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type,recording_interval,params);
[states,~,trans,emis] = bs_fit_hmm(spike_counts,n_states,n_emissions,stim_types,n_stimuli,params);
%confusion = bs_compute_confusion(states,n_states,stim_types,stimuli,n_stimuli);
bs_results_hmm(trans,emis,bin_times,stim_types,states,times_since_stim_on);
correct = (sum(stim_types==-1 & states==1) + sum(stim_types~=1 & states~=1))/length(stim_types);
disp(['Correctly classified states: ' num2str(round(correct*100)) ' % ' ...
    '(random level: ' num2str( round((1-sum(stim_types==-1)/length(stim_types))*100)) ' %)'])


%% HMM fit to PCA
n_emissions = 100;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'states_on_off';
params.binsize = 0.01;
[bin_times,spike_counts,stim_types,times_since_stim_on,~] = bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
[coeff,score,latent] = pca(spike_counts');
score = score(:,1:10); % take first 10 PCs
[states,~,trans,emis] = bs_fit_hmm(score',n_states,n_emissions,stim_types,n_stimuli,params);
bs_results_hmm(trans,emis,bin_times,stim_types,states,times_since_stim_on);
correct = (sum(stim_types==-1 & states==1) + sum(stim_types~=1 & states~=1))/length(stim_types);
disp(['Correctly classified states: ' num2str(round(correct*100)) ' % ' ...
    '(random level: ' num2str( round((1-sum(stim_types==-1)/length(stim_types))*100)) ' %)'])



%% Overfitting HMM (learning all spike counts)
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'states_on_off';
params.binsize = 0.1;
[bin_times,spike_counts,stim_types,times_since_stim_on,~] = bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
n_emissions = length(bin_times);
[states,~,trans,emis] = bs_fit_hmm(spike_counts,n_states,n_emissions,stim_types,n_stimuli,params);
confusion = bs_compute_confusion(states,n_states,stim_types,stimuli,n_stimuli);
bs_results_hmm(trans,emis,bin_times,stim_types,states,times_since_stim_on,confusion);
correct = (sum(stim_types==-1 & states==1) + sum(stim_types~=1 & states~=1))/length(stim_types);
disp(['Correctly classified states: ' num2str(round(correct*100)) ' % ' ...
    '(random level: ' num2str( round((1-sum(stim_types==-1)/length(stim_types))*100)) ' %)'])

%% Fit HMM with 1 ms bins
n_emissions = n_clusters+1; % (in the limit of binsize->0, a single neuron spiking becomes one emission label)
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'states_on_off';
params.binsize = 0.001;
[bin_times,spike_counts,stim_types,times_since_stim_on,~] = bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
[states,p_states,trans,emis] = bs_fit_hmm(spike_counts,n_states,n_emissions,stim_types,n_stimuli,params);
bs_results_hmm(trans,emis,bin_times,stim_types,states,times_since_stim_on);
correct = (sum(stim_types==-1 & states==1) + sum(stim_types~=1 & states~=1))/length(stim_types);
disp(['Correctly classified states: ' num2str(round(correct*100)) ' % ' ...
    '(random level: ' num2str( round((1-sum(stim_types==-1)/length(stim_types))*100)) ' %)'])


%% Optimal HMM fit to n_stimuli + 1 states
n_emissions = 200;
n_states = n_stimuli + 1;
params = bs_default_params();
params.fit_hmm_prompt = 'states_for_all_stimuli';
params.binsize = 0.1;
[bin_times,spike_counts,stim_types,times_since_stim_on,~] = bs_bin_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);

[states,~,trans,emis] = bs_fit_hmm(spike_counts,n_states,n_emissions,stim_types,n_stimuli,params);
bs_results_hmm(trans,emis,bin_times,stim_types,states,times_since_stim_on);

correct = (sum(stim_types==-1 & states==n_states) + sum(stim_types~=1 & states~=n_states))/length(stim_types);
disp(['Correctly classified as no stimulus: ' num2str(round(correct*100)) ' % ' ...
    '(random level: ' num2str( round((1-sum(stim_types==-1)/length(stim_types))*100)) ' %)'])




%%
%- variational auto-encoders
%- CEBRA (latent embeddings of neural space)


