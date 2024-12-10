function bs_results_hmm(trans,emis,bin_times,bin_labels,bin_states,bin_times_rel)

n_states = size(trans,1);
n_emissions = size(emis,2);
labels = unique(bin_labels);
n_labels = length(labels);

confusion_matrix = bs_compute_confusion(bin_states,n_states,bin_labels,labels);

ent = entropy(bin_labels);
mi = mutualinfo(bin_labels,bin_states);
disp(['Entropy in labels = ' num2str(ent,2)]);
disp(['Mutual information in states and labels = ' num2str(mi,2)])
disp(['Relative mutual information = ' num2str(round(mi/ent*100)) , '%'])


% HMM figures
figure

% Transition matrix
subplot(4,2,1) 
imagesc('ydata',1:n_states,'xdata',1:n_states,'cdata',trans)
set(gca,'xtick',1:n_states)
set(gca,'ytick',1:n_states)
axis image
xlabel('To')
ylabel('From')
clim([0 1]);
colorbar

% Emission matrix
subplot(4,2,3)
if ~isstruct(emis)
    imagesc('ydata',1:n_states,'xdata',1:n_emissions,'cdata',emis)
else
    imagesc('ydata',1:n_states,'xdata',1:n_emissions,'cdata',emis.M)
end
    clim([0 1])
ylabel('State')
set(gca,'xtick',[])
set(gca,'ytick',1:n_states)
xlabel('Emission')
axis tight

% Example period
subplot(2,3,3)
xl = mean([bin_times(1) bin_times(end)]) + [-10 10];
bs_plot_labels(bin_times,bin_labels,[0 n_states+1],xl)
bs_plot_states(bin_times,bin_states,xl)
set(gca,'YTick',1:n_states)
ylabel('State');
xlim(xl);
ylim([0 n_states+1])
xlabel('Time (s)')

% Confusion sum_label = 1
subplot(2,3,4)
imagesc('ydata',1:n_states,'xdata',1:n_labels,'cdata',confusion_matrix./sum(confusion_matrix))
clim([0 1]);
axis image
colorbar
ylabel('State')
xlabel('Label')
title('\Sigma_{label} = 1')

% Confusion sum_state = 1
subplot(2,3,5)
imagesc('ydata',1:n_states,'xdata',1:n_labels,'cdata',confusion_matrix./sum(confusion_matrix,2))
clim([0 1]);
axis image
colorbar
ylabel('State')
xlabel('Label')
title('\Sigma_{state} = 1')

% Rel times
subplot(2,3,6)
hold on
xlabel('Rel. time (s)')
ylabel('State')
incidence = zeros(n_states,100-1);
edges = linspace(0,1.5,100);
for s=1:n_states
    incidence(s,:) = histcounts(bin_times_rel(bin_states==s),edges);
end
incidence = incidence./(0.001+sum(incidence));
imagesc('ydata',1:n_states,'xdata',edges(1:end-1),'cdata',incidence); % xdata half sample shifted
axis tight
set(gca,'ytick',1:n_states);
clim([0 1])

