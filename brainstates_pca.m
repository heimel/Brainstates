%brainstates_pca
%
% Deprecated. Use brainstates.mlx
%
% 2024, Alexander Heimel

%% To do
%  Apply GPFA (Yu et al. J Neurophys 2009) (run alexander_test with data split into trials)
%  Apply augmentation vs smoothing or binning
%  Compute log-likelihood of observed data over different seeds 
%  Make lecture notes
%  Check out other recordings
%  Variational auto-encoders
%  Test CEBRA datasets
%  CEBRA (latent embeddings of neural space)
%
%  Good reference for HMM fit to neural data: Bagi et al. Curr Biol 2022

%  Python implementations:
%  Extensive implementation of HMM: https://github.com/lindermanlab/ssm
%  Gausssian emission HMM: https://hmmlearn.readthedocs.io/en/latest/api.html#hmmlearn.hmm.GaussianHMM
%
%  Matlab:
%  hmmtrain
%  for gaussian: 
%  https://www.cs.ubc.ca/~murphyk/Software/HMM/hmm.html
%   -> moved to https://github.com/probml/pmtk3/tree/master/toolbox/LatentVariableModels/hmm
%
%  https://nl.mathworks.com/matlabcentral/fileexchange/55866-hidden-markov-model-toolbox-hmm
%  https://nl.mathworks.com/matlabcentral/fileexchange/20712-em-for-hmm-multivariate-gaussian-processes



%% Load data
params = bs_default_params();
addpath(fullfile(fileparts(mfilename("fullpath")),'Mutualinfo'));

 params.dataset = 'Topo6_20220301_AP';
 params.block = 2;
% Best fit no hinst (with hints), two states 22% (45%); all states 19% (34%)

%params.dataset = 'Topo6_20220301_AP';
%params.block = 4;
% Best fit no hints (with hints), two states 12% (41%); all states 15% (26%)

% params.dataset = 'Exp2019-12-16_MP3_S2L5_AP';
% params.block = 3;
% Best fit no hints (with hints), two states XX% (23%); all states XX% (XX%)

[spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval] = bs_load_data(params);


params.smooth = true;
params.smooth_function = 'smoothen';
params.smooth_window = 0.1;
params.binsize = 0.01;

[bin_times,bin_counts,bin_labels,bin_times_rel] = ...
    bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
n_bins = length(bin_times);
n_clusters = length(spike_times);

ind_bin_subset = round(linspace(1,n_bins,min([1000 n_bins]))); % total runs show change over time in t-sne plot

return

%% Get correlated pairs
no_correlation = true; %#ok<UNRCH>
c1 = 0;

while c1<n_clusters && no_correlation
    c1 = c1 + 1;
    c2 = c1;
    while c2<n_clusters && no_correlation
        c2 = c2+1;
        x = corrcoef(bin_counts(c1,ind_bin_subset),bin_counts(c2,ind_bin_subset));
        if x(1,2)>0.5
            disp(['Good correlation for ' num2str(c1) ' and ' num2str(c2)]);
            no_correlation = false;
        end
    end
end


%% PCA Example with two neurons
figure
scatter(bin_counts(c1,ind_bin_subset) + 0.5*(rand(1,length(ind_bin_subset))-0.5),...
    bin_counts(c2,ind_bin_subset) + 0.5*(rand(1,length(ind_bin_subset))-0.5),'filled')
xlabel('Spike count neuron 1 (jittered)');
ylabel('Spike count neuron 2 (jittered)');
hold on
axis equal

x = bin_counts([c1 c2],:)';
[coeff,score,latent] = pca(x);

mx = mean(x);
plot(mx(1),mx(2),'or','MarkerFaceColor','r')
plot(mx(1)+[0 coeff(1,1)],mx(2)+[0 coeff(2,1)],'r-','LineWidth',2)
plot(mx(1)+[0 coeff(1,2)],mx(2)+[0 coeff(2,2)],'g-','LineWidth',2 )


% x_reconstruct = score*coeff' + repmat(mx,n_bins,1);
% scatter(x_reconstruct(ind_bin_subset,1),x_reconstruct(ind_bin_subset,2),'filled')


%% PCA for all clusters
disp('Computing PCA')
[coeff,score,latent] = pca(bin_counts');
mx = mean(bin_counts'); %#ok<UDIM>


figure
subplot(2,2,1)
bar(latent*100)
ylabel('Explained (%) ')
xlabel('Component')
xlim([0.5 10.5])
set(gca,'xtick',1:2:10)

subplot(2,2,2)
imagesc(bin_counts)
ylabel('Neuron')
xlabel('Sample')
title('Origina data')

subplot(2,2,3)
hold on
plot(bin_times,score(:,1))
plot(bin_times,score(:,2))
plot(bin_times,score(:,3))
plot(bin_times,bin_times_rel)
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
scatter(score(ind_bin_subset,1),score(ind_bin_subset,2),30,bin_times_rel(ind_bin_subset),'filled')
axis square
colorbar
xlabel('PC1')
ylabel('PC2')

hold on
[~,ind] = sort(bin_times_rel);
x = movingavr(score(ind,1),2000);
y = movingavr(score(ind,2),2000);
plot(x,y,'k-','linewidth',2)
plot(x(1),y(1),'og','MarkerFaceColor','g')
plot(x(end),y(end),'og','MarkerFaceColor','r')


%% PC1 vs PC2 with color stim_type
figure;
scatter(score(ind_bin_subset,1),score(ind_bin_subset,2),30,bin_labels(ind_bin_subset),'filled')
colormap(parula(max(bin_labels(ind_bin_subset))+2))

axis square
colorbar
xlabel('PC1')
ylabel('PC2')

%% Factor analysis
[~,~,~,~,fa_score] = factoran(bin_counts',10);
figure;
subplot(1,2,1)
scatter(fa_score(ind_bin_subset,1),fa_score(ind_bin_subset,2),30,bin_times_rel(ind_bin_subset),'filled')
axis square
colorbar
xlabel('Score 1')
ylabel('Score 2')
hold on
[~,ind] = sort(bin_times_rel);
x = movingavr(fa_score(ind,1),2000);
y = movingavr(fa_score(ind,2),2000);
plot(x,y,'k-','linewidth',2)
plot(x(1),y(1),'og','MarkerFaceColor','g')
plot(x(end),y(end),'og','MarkerFaceColor','r')

subplot(1,2,2)
scatter(fa_score(ind_bin_subset,1),fa_score(ind_bin_subset,2),30,bin_labels(ind_bin_subset),'filled')
colormap(parula(max(bin_labels(ind_bin_subset))+2))
axis square
colorbar
xlabel('Score 1')
ylabel('Score 2')
hold on
[~,ind] = sort(bin_labels);
x = movingavr(fa_score(ind,1),2000);
y = movingavr(fa_score(ind,2),2000);
plot(x,y,'k-','linewidth',2)
plot(x(1),y(1),'og','MarkerFaceColor','g')
plot(x(end),y(end),'og','MarkerFaceColor','r')



%% t-sne for subset
disp('Computing t-sne')
y = tsne(bin_counts(:,ind_bin_subset)');

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
scatter(y(:,1),y(:,2),30,bin_times_rel(ind_bin_subset),'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,periodic_colormap(128))
colorbar
title('Time since stim on')

% Show t-sne with stim type as color
h = subplot(2,2,3);

scatter(y(:,1),y(:,2),30,bin_labels(ind_bin_subset)+2,'filled')
axis square
xlabel('t-sne 1')
ylabel('t-sne 2')
colormap(h,parula(max(bin_labels(ind_bin_subset))+2))
colorbar
title('Stim type')

%% Overfitting HMM (learning all spike counts) (100%)
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'emission_freqs';
params.binsize = 0.1;
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)), recording_interval,params);
n_emissions = length(bin_times);
[bin_states,~,trans,emis] = bs_fit_hmm(bin_counts,n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);

%% HMM fit to two states (11%) 
n_emissions = 10;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'emission_freqs';
params.hmm_algorithm = 'Viterbi';
params.smooth = false;
params.binsize = 0.1;
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)),recording_interval,params);
[bin_states,~,trans,emis] = bs_fit_hmm(bin_counts,n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);


n_emissions = 100;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'emission_freqs';
params.hmm_algorithm = 'Viterbi';
params.smooth = true;
params.smooth_function = 'smoothen';
params.smooth_window = 0.1;
params.binsize = 0.01;
params.square_root_transformation = true;

[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)),recording_interval,params);
[coeff,score,latent] = pca(bin_counts');
score = score(:,1:30);
[bin_states,~,trans,emis] = bs_fit_hmm(score',n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);


%% HMM fit to two states on 10 FA components (37%)
n_emissions = 100;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'emission_freqs';
params.hmm_algorithm = 'Viterbi';
params.smooth = true;
params.smooth_function = 'smoothen';
params.smooth_window = 0.1;
params.binsize = 0.02;
params.square_root_transformation = true;
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)),recording_interval,params);

%[bin_states,~,trans,emis] = bs_fit_hmm(fa_score',n_states,n_emissions,bin_labels,params);

[~,~,~,~,fa_score] = factoran(bin_counts',10);
[bin_states,~,trans,emis] = bs_fit_hmm(fa_score',n_states,n_emissions,bin_labels,params);

bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);

%% GM-HMM fit to two states (41%)
n_emissions = 10;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'mix_gauss_hints';
params.smooth = true;
params.smooth_function = 'smoothen';
params.smooth_window = 0.05;
params.binsize = 0.01;
params.hmm_version = 'Probml_mixgausstied';
params.square_root_transformation = true;
params.n_mixgauss = 6; 
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)),recording_interval,params);
[~,~,~,~,fa_score] = factoran(bin_counts',10);
[bin_states,~,trans,emis] = bs_fit_hmm(fa_score',n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);


%% HMM fit to two states, no hints (11%)
n_emissions = 10;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'no_hints';
params.hmm_algorithm = 'BaumWelch';
%params.reproducible = false;
params.binsize = 0.1;
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)),recording_interval,params);
[bin_states,~,trans,emis] = bs_fit_hmm(bin_counts,n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);


%% GM-HMM fit to two states, no hints (12%)
n_emissions = 10;
n_states = 2;
params = bs_default_params();
params.fit_hmm_prompt = 'no_mix_gauss_hints';
params.smooth = true;
params.smooth_function = 'smoothen';
params.smooth_window = 0.1;
params.binsize = 0.01;
params.hmm_version = 'Probml_mixgausstied';
params.square_root_transformation = true;
params.n_mixgauss = 12; 
params.reproducible = true;
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,ones(size(trial_stim_type)),recording_interval,params);
[~,~,~,~,fa_score] = factoran(bin_counts',10);
[bin_states,~,trans,emis] = bs_fit_hmm(fa_score',n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);




%% GM-HMM fit to all labels (26%)
n_emissions = 10;
params = bs_default_params();
params.fit_hmm_prompt = 'mix_gauss_hints';
params.smooth = true;
params.smooth_function = 'smoothen';
params.smooth_window = 0.05;
params.binsize = 0.01;
params.hmm_version = 'Probml_mixgausstied';
params.square_root_transformation = true;
params.n_mixgauss = 6; 
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
labels = unique(bin_labels);
n_states = length(labels);
[~,~,~,~,fa_score] = factoran(bin_counts',10);
[bin_states,~,trans,emis] = bs_fit_hmm(fa_score',n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel);




%% HMM fit to all labels, no hints (15%) 
n_emissions = 100; % 200
params = bs_default_params();
params.fit_hmm_prompt = 'no_hints';
params.hmm_algorithm = 'BaumWelch';
params.binsize = 0.1;
[bin_times,bin_counts,bin_labels,bin_times_rel,~] = bs_bin_labeled_data(spike_times,trial_stim_on,trial_stim_off,trial_stim_type, recording_interval,params);
labels = unique(bin_labels);
n_states = length(labels);
[states,~,trans,emis] = bs_fit_hmm(bin_counts,n_states,n_emissions,bin_labels,params);
bs_results_hmm(trans,emis,bin_times,bin_labels,states,bin_times_rel);


