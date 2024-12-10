function [bin_states,p_states,trans,emis] = bs_fit_hmm(data,n_states,n_emissions,hmm_type,bin_labels,params)

if nargin<4 || isempty(hmm_type)
    hmm_type = 'discrete';
end


disp('Fitting HMM')

if params.reproducible
    rng(1);
end

switch hmm_type
    case 'discrete'
        [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,data,bin_labels);
        [trans,emis] = hmmtrain(data(:)',trans_guess,emis_guess,'Algorithm',params.hmm_algorithm,'Maxiterations',1000);
        p_states = hmmdecode(data',trans,emis);
        [~,bin_states] = max(p_states);
    case 'gm'
        %bin_clusters = kmeans(data',n_emissions);
%        [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,bin_clusters,bin_labels,data);
        [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,[],bin_labels,data);

        model = hmmFit(data, n_states, 'mixgausstied', 'trans0',trans_guess,'nmix',n_emissions,'emission0',emis_guess);
        bin_states = hmmMap(model, data);
        trans = model.A;
        emis = model.emission;
        p_states = [];
    otherwise
        logmsg(['HMM type ' hmm_type ' is not implemented.']);
end
disp('Done fitting')
