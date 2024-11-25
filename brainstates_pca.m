
params = bs_default_params();


%% Load data
% d = dir(fullfile(params.datafolder,'*.mat'));
% v1 clusters in 2,9
% stimulus info in 9

load(fullfile(params.datafolder,'Exp2019-12-16_MP3_S2L5_AP.mat'))

% params.datafolder = 'W:\Montijn\MontijnNPX2020\Data_collection\DataPreProcessed';
% load(fullfile(params.datafolder,'Topo6_20220301_AP.mat'))

% select V1 clusters
sAP.sCluster = sAP.sCluster(contains(string({sAP.sCluster.Area}),"Primary visual area"));
n_clusters = length(sAP.sCluster);
disp([num2str(n_clusters) ' clusters in V1'])

%% Parse stimulus info

% sAP.cellStim{1}.structEP.vecTrialStimOnSecs contains stimulus onset times
% sAP.cellStim{1}.structEP.vecTrialStimOffSecs contains stimulus offset times
% sAP.cellStim{1}.structEP.vecTrialStimTypes contains stimulus types per trial
% sAP.cellStim{1}.structEP.sStimObject(t).Orientation contains orientation of stimulus type s


%% Get data into right format
% stim_on = sAP.cellStim{3}.structEP.vecTrialStimOnSecs
% stim_off = sAP.cellStim{3}.structEP.vecTrialStimOffSecs


%% Get first and last spike time
stim_number = 3;


first_spike_time = inf;
last_spike_time = -inf;
for c = 1:n_clusters
    if sAP.sCluster(c).SpikeTimes(1) < first_spike_time
        first_spike_time = sAP.sCluster(c).SpikeTimes(1);
    end
    if sAP.sCluster(c).SpikeTimes(end) > last_spike_time
        last_spike_time = sAP.sCluster(c).SpikeTimes(end);
    end
end
recording_interval = [floor(first_spike_time) ceil(last_spike_time)];

% Take subset 
recording_interval(1) = sAP.cellStim{stim_number}.structEP.vecStimOnTime(1) - 1;
recording_interval(2) = sAP.cellStim{stim_number}.structEP.vecStimOffTime(end) + 1;

% recording_interval(1) = sAP.cellBlock{stim_number}.vecStimOnTime(1) - 1;
% recording_interval(2) = sAP.cellBlock{stim_number}.vecStimOffTime(end) + 1;

%sAP.cellStim = sAP.cellBlock;

%% Bin data
n_bins = ceil( (recording_interval(2)-recording_interval(1)) / params.binsize );
bin_time = recording_interval(1) + params.binsize * (1:n_bins) - 0.5* params.binsize;



spike_counts = zeros(n_clusters,n_bins);
spike_counts_norm = zeros(n_clusters,n_bins);
spike_counts_jittered = zeros(n_clusters,n_bins); % for plotting

for c = 1:n_clusters
    for i = 1:length(sAP.sCluster(c).SpikeTimes)
        bin = ceil( (sAP.sCluster(c).SpikeTimes(i)-recording_interval(1)) /params.binsize);
        if bin<1
            continue
        end
        if bin>n_bins
            break
        end
        spike_counts(c,bin) = spike_counts(c,bin) + 1;
    end

    % smooting
    if params.smooth
        spike_counts(c,:) = smooth(spike_counts(c,:),params.smooth_window);
    end

    spike_counts_norm(c,:) = spike_counts(c,:)/max(spike_counts(c,:));
    spike_counts_jittered(c,:) = spike_counts(c,:) + 0.5*(rand(size(spike_counts(c,:)))-0.5);

end



ind_bin_subset = round(linspace(1,n_bins,5000)); % total runs show change over time in t-sne plot
%ind_bin_subset = round(linspace(find(bin_time>3400,1),find(bin_time>5900,1),5000));  % 2nd drifting gratings with 0,5 ,90 and 95 degrees orientation, cellStim{3}
%ind_bin_subset = round(linspace(find(bin_time>2500,1),find(bin_time>3200,1),5000));  % 1st drifitng gratings with 24 orientations


%% Creative relative stim times
disp('Computing relative stim times')
time_since_stim_on = params.max_time_since_stim_on + params.binsize * (1:n_bins);

for i = 1:length(sAP.cellStim)
    if ~strcmp(sAP.cellStim{i}.structEP.strFile,'RunDriftingGratings')
        continue % only take drifting gratings here
    end
    vecStimOnTime = sAP.cellStim{i}.structEP.vecStimOnTime;
    for j = 1:length(vecStimOnTime)
        ind = find(bin_time>vecStimOnTime(j),1);
        dt = bin_time(ind) - vecStimOnTime(j); 
        time_since_stim_on(ind:end) = time_since_stim_on(ind:end) - time_since_stim_on(ind) + dt;
    end
end
time_since_stim_on( time_since_stim_on > params.max_time_since_stim_on) = params.max_time_since_stim_on;

time_since_stim_off = params.max_time_since_stim_off + params.binsize * (1:n_bins);
for i = 1:length(sAP.cellStim)
    if ~strcmp(sAP.cellStim{i}.structEP.strFile,'RunDriftingGratings')
        continue % only take drifting gratings here
    end

    vecStimOffTime = sAP.cellStim{i}.structEP.vecStimOffTime;
    for j = 1:length(vecStimOffTime)
        ind = find(bin_time>vecStimOffTime(j),1);
        dt = bin_time(ind) - vecStimOffTime(j); 
        time_since_stim_off(ind:end) = time_since_stim_off(ind:end) - time_since_stim_off(ind) + dt;
    end
end
time_since_stim_off( time_since_stim_off > params.max_time_since_stim_off) = params.max_time_since_stim_off;

%% Create stimulus type vector
stim_type = -1 * ones(1,n_bins);
for i = 1:length(sAP.cellStim)
    if ~strcmp(sAP.cellStim{i}.structEP.strFile,'RunDriftingGratings')
        continue % only take drifting gratings here
    end
    vecStimOnTime = sAP.cellStim{i}.structEP.vecStimOnTime;
    vecStimOffTime = sAP.cellStim{i}.structEP.vecStimOffTime;
    for j = 1:length(vecStimOnTime)
        ind = find(bin_time>vecStimOnTime(j) & bin_time<vecStimOffTime(j));
        stim_type(ind) = sAP.cellStim{i}.structEP.vecTrialStimTypes(j);
    end

end



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

% %% Show data
% 
% imagesc(spike_counts_norm)
% ylabel('Neuron')
% xlabel('Sample')

%% PCA Example with two neurons
%c1 = 22; c2 = 27;
%c1 = 7; c2 = 45;
%c1 = 38; c2 = 39;
%c1 = 44; c2 = 45;
figure
scatter(spike_counts_jittered(c1,ind_bin_subset),...
    spike_counts_jittered(c2,ind_bin_subset),'filled')
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
imagesc(spike_counts_norm)
ylabel('Neuron')
xlabel('Sample')
title('Origina data')

subplot(2,2,3)
hold on
plot(bin_time,score(:,1))
plot(bin_time,score(:,2))
plot(bin_time,score(:,3))
plot(bin_time,time_since_stim_on)
xlabel('Time (s)')
ylabel('Score')


subplot(2,2,4)
%x = score(:,1)*coeff(:,1)' + repmat(mx,n_bins,1);
x = score(:,(1:3))*coeff(:,(1:3))' + repmat(mx,n_bins,1);
x_norm = x ./ repmat(max(x),n_bins,1);
imagesc(x_norm')
%imagesc(x')
ylabel('Neuron')
xlabel('Sample')
title('First 3 PCs')

%% PC1 vs PC2 with color time since stim on
figure;
scatter(score(ind_bin_subset,1),score(ind_bin_subset,2),30,time_since_stim_on(ind_bin_subset),'filled')
axis square
colorbar
xlabel('PC1')
ylabel('PC2')

%% PC1 vs PC2 with color stim_type
figure;
scatter(score(ind_bin_subset,1),score(ind_bin_subset,2),30,stim_type(ind_bin_subset),'filled')
colormap(parula(max(stim_type(ind_bin_subset))+2))

axis square
colorbar
xlabel('PC1')
ylabel('PC2')


%% t-sne for subset
disp('Computing t-sne')
y = tsne(spike_counts(:,ind_bin_subset)');

%% Show t-sne results

% Show t-sne with time as color
figure
h = subplot(2,2,1);
colormap(h,'parula')
scatter(y(:,1),y(:,2),30,bin_time(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colorbar
title('Absolute time')

% Show t-sne with time since stim onset as color
h = subplot(2,2,2);
scatter(y(:,1),y(:,2),30,time_since_stim_on(ind_bin_subset),'filled')
%scatter(y(:,1),y(:,2),30,time_since_stim_off(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,periodic_colormap(128))
colorbar
title('Time since stim on')

% Show t-sne with stim type as color
h = subplot(2,2,3);

scatter(y(:,1),y(:,2),30,stim_type(ind_bin_subset)+2,'filled')
%scatter(y(:,1),y(:,2),30,time_since_stim_off(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,parula(max(stim_type(ind_bin_subset))+2))
colorbar
title('Stim type')

%% k-means clustering

%idx = kmeans(spike_counts(:,ind_bin_subset)',10);
n_emissions = 100; % 10
idx = kmeans(spike_counts',n_emissions);
%idx = kmeans(y,5)

% Show t-sne with cluster as color

h = subplot(2,2,4);
scatter(y(:,1),y(:,2),30,idx(ind_bin_subset),'filled')
%scatter(y(:,1),y(:,2),30,time_since_stim_off(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,parula(max(idx)))
colorbar
title('k-means')

% %% Hidden Markov Model Example
% trans = [0.9 0.1 0.0; ...
%          0.0 0.9 0.1; ...
%          0.1 0.0 0.9 ];
% emis = [0.8 0.2 0.0  0;...
%         0.2 0.6 0.2  0;...
%         0.0 0.0 0.8 0.2];
% 
% seq1 = hmmgenerate(5000,trans,emis);
% seq2 = hmmgenerate(5000,trans,emis);
% seqs = {seq1,seq2};
% 
% n_emissions = size(emis,2);
% n_states = 3;
% emis_guess = ones(n_states,n_emissions)/n_emissions;
% trans_guess = [ 0.8 0.1 0.1;...
%                 0.1 0.8 0.1;  ...
%                 0.1 0.1 0.8];
% trans_guess = eye(n_states);
% trans_guess = ones(n_states,n_states)/n_states;
% trans_guess = rand(n_states);
% % [trans_est,emis_est] = hmmestimate(seqs,states)
% [trans_est,emis_est] = hmmtrain(seqs,trans_guess,emis_guess)

%% HMM fit to data
disp('Fitting HMM')
%rng(1)

% no prompting
% n_states = 10;  % 6 is okayish
% emis_guess = ones(n_states,n_emissions)/n_emissions;
% trans_guess = rand(n_states);

% complete prompting with 3 states
n_states = 3; 
emis_guess = zeros(n_states,n_emissions);
for i=1:n_emissions
    s = mode(stim_type(idx==i));
    switch s
        case -1
            emis_guess(1,i) = 1;
        case {2,3}
            emis_guess(2,i) = 1;
        case {4,5}
            emis_guess(3,i) = 1;
    end
end
emis_guess = emis_guess + 0.2 *rand(size(emis_guess));
emis_guess = emis_guess./sum(emis_guess,2);
trans_guess = [0.80 0.10 0.10; ...
               0.20 0.80 0.00; ...
               0.20 0.00 0.80];
trans_guess = trans_guess + 0.2 * rand(size(trans_guess));
trans_guess = trans_guess./sum(trans_guess,2);


% complete prompting with 5 states (two states per stimulus)
n_states = 5; 
emis_guess = zeros(n_states,n_emissions);
for i=1:n_emissions
    s = mode(stim_type(idx==i));
    switch s
        case -1
            emis_guess(1,i) = 1;
        case {2,3}
            emis_guess(2,i) = 1;
            emis_guess(3,i) = 1;
        case {4,5}
            emis_guess(4,i) = 1;
            emis_guess(5,i) = 1;
    end
end
emis_guess = emis_guess + 0.2 *rand(size(emis_guess));
emis_guess = emis_guess./sum(emis_guess,2);
trans_guess = 0.1 * ones(n_states,n_states) + eye(n_states);
trans_guess = trans_guess + 0.2 * rand(size(trans_guess));
trans_guess = trans_guess./sum(trans_guess,2);




%
[trans_est,emis_est] = hmmtrain(idx',trans_guess,emis_guess);

p_states_est = hmmdecode(idx',trans_est,emis_est);
[m,states]=max(p_states_est);
disp('Done fitting HMM')

% Compute confusion matrix
stimuli = unique(stim_type);
n_stimuli = length(stimuli);

confusion = zeros(n_states,n_stimuli);
for state = 1:n_states
    for stim = 1:n_stimuli
        confusion(state,stim) = sum(states==state & stim_type==stimuli(stim));
    end
end
for i=1:n_stimuli
    confusion(:,i) = confusion(:,i)/sum(confusion(:,i));
end

%% HMM figures
figure

subplot(2,2,1)
imagesc('ydata',1:n_states,'xdata',1:n_states,'cdata',trans_est)
axis image
colorbar

subplot(2,2,2)
imagesc('ydata',1:n_states,'xdata',1:n_stimuli,'cdata',confusion)
axis image
colorbar
ylabel('State')
xlabel('Stim type')
title('Confusion matrix')

subplot(4,2,5)
plot(bin_time,stim_type,'.-');
xlim([4000 4030])
subplot(4,2,7)
plot(bin_time,states,'.-');
xlim([4000 4030])

subplot(2,2,4)
hold on
% plot(time_since_stim_on,states+0.5*rand(size(states)),'.')
ind_sel = false(size(stim_type));
ind_sel(round(linspace(1,n_bins,4000))) = true;

ind = (stim_type==-1) & ind_sel;
plot(time_since_stim_on(ind),states(ind)+0.5*rand(size(states(ind))),'.')
ind = (stim_type==1|stim_type==2) & ind_sel;
plot(time_since_stim_on(ind),states(ind)+0.5*rand(size(states(ind))),'.')
ind = (stim_type==3|stim_type==4) & ind_sel;
plot(time_since_stim_on(ind),states(ind)+0.5*rand(size(states(ind))),'.')
xlabel('Time from stim on (s)')
ylabel('State')
hold off

%%
% figure;
% plot(time_since_stim_on,idx+0.5*rand(size(idx)),'.')
% xlabel('Time from stim on (s)')
% ylabel('Idx')

%%
%- variational auto-encoders
%- CEBRA (latent embeddings of neural space)
