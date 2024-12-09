function [bin_states,p_states,trans,emis] = bs_fit_hmm(bin_counts,n_states,n_emissions,bin_labels,params)
disp('Fitting HMM')

if params.reproducible
    rng(1);
end

switch params.cluster_method
    case 'kmeans'
        bin_clusters = kmeans(bin_counts',n_emissions); % clustering in original space
    case 'gmmodel'
        options = statset('MaxIter',1000);   %    'RegularizationValue',0.1)
        gmmodel = fitgmdist(bin_counts',n_emissions,'SharedCovariance',true,'Options',options);
        bin_clusters = cluster(gmmodel,bin_counts');
end

switch params.hmm_version
    case 'Mathworks'
        [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,bin_clusters,bin_labels);
        [trans,emis] = hmmtrain(bin_clusters(:)',trans_guess,emis_guess,'Algorithm',params.hmm_algorithm,'Maxiterations',1000);
        p_states = hmmdecode(bin_clusters',trans,emis);
        [~,bin_states] = max(p_states);
    case 'Probml'
        [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,bin_clusters,bin_labels,bin_counts);

        labels = unique(bin_labels);
        n_labels = length(labels);
        % emis_guess = zeros(n_states,n_emissions);
        if n_states~=n_labels
            disp('Unequal number of labels and states. Cannot hint emission matrix')
        end

        n_dim = size(bin_counts,1);
        mu = zeros(n_dim,n_states);
        sigma = zeros(n_dim,n_dim,n_states);
        for s = 1:n_states
            mu(:,s) = mean(bin_counts(:,bin_labels==labels(s)),2);
            sigma(:,:,s) = cov(bin_counts(:,bin_labels==labels(s))');
        end
        emis_guess = condGaussCpdCreate(mu, sigma);
        [model, loglikHist] = hmmFit(bin_counts, n_states, 'gauss', 'trans0',trans_guess,'emission0',emis_guess);
        bin_states = hmmMap(model, bin_counts);
        trans = model.A;
        emis = model.emission;
        p_states = [];

    
    case 'Probml_mixgausstied'
        [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,bin_clusters,bin_labels,bin_counts);
        model = hmmFit(bin_counts, n_states, 'mixgausstied', 'trans0',trans_guess,'nmix',params.n_mixgauss,'emission0',emis_guess);
        bin_states = hmmMap(model, bin_counts);
        trans = model.A;
        emis = model.emission;
        p_states = [];
    otherwise
        logmsg(['HMM version ' params.hmm_version ' is not implemented.']);
end
disp('Done fitting')
