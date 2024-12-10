function [confusion_matrix,confusion_in_state,confusion_in_label] = bs_compute_confusion(states,n_states,stim_types,stimuli)

% No confusion is if each label completely gets mapped to one state
%
% 2024, Alexander Heimel

n_stimuli = length(stimuli);

confusion_matrix = zeros(n_states,n_stimuli);
for state = 1:n_states
    for stim = 1:n_stimuli
        confusion_matrix(state,stim) = sum(states==state & stim_types==stimuli(stim));
    end
end

confusion_in_state = 1-mean(max(confusion_matrix./sum(confusion_matrix)));  
confusion_in_label = 1-mean(max(confusion_matrix./sum(confusion_matrix,2),[],2));  

