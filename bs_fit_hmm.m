function [states,p_states,trans,emis] = bs_fit_hmm(spike_counts,n_states,n_emissions,stim_type,n_stimuli,params)
disp('Fitting HMM')

if params.reproducible
    rng(1);
end

switch params.cluster_method
    case 'kmeans'
        idx = kmeans(spike_counts',n_emissions); % clustering in original space
    case 'gmmodel'
        options = statset('MaxIter',1000);   %    'RegularizationValue',0.1)
        gmmodel = fitgmdist(spike_counts',n_emissions,'SharedCovariance',true,'Options',options);
        idx = cluster(gmmodel,spike_counts');
end



[trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,idx,stim_type,n_stimuli);
[trans,emis] = hmmtrain(idx',trans_guess,emis_guess,'Algorithm','Viterbi');
%[trans,emis] = hmmtrain(idx',trans_guess,emis_guess);

p_states = hmmdecode(idx',trans,emis);

[~,states] = max(p_states);

