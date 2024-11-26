function bs_results_hmm(trans,emis,bin_times,stim_types,states,time_since_stim_on)

n_bins = length(stim_types);
n_states = size(trans,1);
n_emissions = size(emis,2);
stimuli = unique(stim_types);
n_stimuli = length(stimuli);

[confusion_in_state,confusion_in_label,confusion_matrix] = bs_compute_confusion(states,n_states,stim_types,stimuli);
disp(['Confusion in state = ' num2str(confusion_in_state) '. Lower is better.'])
disp(['Confusion in label = ' num2str(confusion_in_label) '. Lower is better.'])

% HMM figures

figure

subplot(2,2,1)
imagesc('ydata',1:n_states,'xdata',1:n_states,'cdata',trans)
set(gca,'xtick',1:n_states)
set(gca,'ytick',1:n_states)
axis image
title('Transition matrix')
colorbar

subplot(4,2,2)
imagesc('ydata',1:n_states,'xdata',1:n_emissions,'cdata',emis)
axis image
colorbar
ylabel('State')
xlabel('Emission')
title('Emission matrix')

subplot(4,4,7)
imagesc('ydata',1:n_states,'xdata',1:n_stimuli,'cdata',confusion_matrix./sum(confusion_matrix))
axis image
colorbar
ylabel('State')
xlabel('Stim type')
title('Confusion matrix (sum over stim = 1)')

subplot(4,4,8)
imagesc('ydata',1:n_states,'xdata',1:n_stimuli,'cdata',confusion_matrix./sum(confusion_matrix,2))
axis image
colorbar
ylabel('State')
xlabel('Stim type')
title('Confusion matrix (sum over state = 1)')


subplot(4,2,5)
plot(bin_times,stim_types,'.-');
xl = mean([bin_times(1) bin_times(end)]) + [-10 10];
xlim(xl);
subplot(4,2,7)
plot(bin_times,states,'.-');
xlim(xl);

subplot(2,2,4)
hold on
% plot(time_since_stim_on,states+0.5*rand(size(states)),'.')
ind_sel = false(size(stim_types));
ind_sel(round(linspace(1,n_bins,4000))) = true;

for i = 1:n_stimuli
    ind = (stim_types==stimuli(i)) & ind_sel;
    plot(time_since_stim_on(ind),states(ind)+0.5*rand(size(states(ind))),'.')
end
xlabel('Time from stim on (s)')
ylabel('State')

