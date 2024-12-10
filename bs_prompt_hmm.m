function [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,bin_clusters,bin_labels,bin_counts)

labels = unique(bin_labels);
n_labels = length(labels);
emis_guess = zeros(n_states,n_emissions);

switch params.fit_hmm_prompt
    case 'emission_freqs'
        hint_matrix = zeros(n_labels,n_emissions);
        for i=1:n_emissions
            for j=1:n_labels
                hint_matrix(j,i) = sum(bin_labels(bin_clusters==i)==labels(j));
            end
        end
        hint_matrix = hint_matrix./sum(hint_matrix,2);
        emis_guess = hint_matrix;
    case 'no_hints'
        emis_guess = ones(n_states,n_emissions)/n_emissions;
    case 'mix_gauss_hints'
        labels = unique(bin_labels);
        n_labels = length(labels);
        if n_states~=n_labels
            disp('Cannot hint emission matrix')
        end

        n_dim = size(bin_counts,1);
        n_mix = n_emissions;
        mu = zeros(n_dim,n_states*n_mix);
        sigma = zeros(n_dim,n_dim,n_states*n_mix);
        mix_weights = zeros(n_states,n_states*n_mix);
        for s = 1:n_states
            data = bin_counts(:,bin_labels==labels(s))';
            mix_gauss = mixGaussFit(data, n_mix);
            ind = (1:n_mix) + (s-1)*n_mix;
            mu(:,ind) = mix_gauss.cpd.mu;
            sigma(:,:,ind) = mix_gauss.cpd.Sigma;
            mix_weights(s,ind) = mix_gauss.mixWeight;
        end
        emis_guess = condMixGaussTiedCpdCreate(mu, sigma, mix_weights);

    case 'no_mix_gauss_hints'
        labels = unique(bin_labels);
        n_labels = length(labels);
        if n_states~=n_labels
            disp('Cannot hint emission matrix')
        end

        n_mix = n_emissions;
        mix_weights = zeros(n_states,n_mix);

        data = bin_counts';
        mix_gauss = mixGaussFit(data, n_mix);
        mu = mix_gauss.cpd.mu;
        sigma = mix_gauss.cpd.Sigma;

        for s=1:n_states
            mix_weights(s,:) = mix_gauss.mixWeight;
        end
        emis_guess = condMixGaussTiedCpdCreate(mu, sigma, mix_weights);
        emis_guess.M = rand(n_states,n_mix);
        emis_guess.M = emis_guess.M ./sum(emis_guess.M ); % just makes everything 1/n_states though

        emis_guess.M = zeros(n_states,n_mix);
        emis_guess.M(1,1:n_mix/2) = ones(1,n_mix/2);
        emis_guess.M(2,n_mix/2+1:end) = ones(1,n_mix/2);

end


% Prompting trans_guess
label_duration = 1.0; %s
trans_guess = ones(n_states,n_states)* 1/(n_states*label_duration/params.binsize) + ...
    eye(n_states)*(1-1/(label_duration/params.binsize) );
trans_guess = trans_guess .* (1 + 0.1* rand(size(trans_guess)));
trans_guess = trans_guess./sum(trans_guess,2);

