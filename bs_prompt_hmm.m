function [trans_guess,emis_guess] = bs_prompt_hmm( params,n_states,n_emissions,idx,stim_type,n_stimuli)

switch params.fit_hmm_prompt
    case 'states_for_all_stimuli'
        emis_guess = zeros(n_states,n_emissions);

        % prompt states (putting stim information in)
        for i=1:n_emissions
            s = mode(stim_type(idx==i));
            switch s
                case -1
                    emis_guess(n_stimuli+1,i) = 1;
                otherwise
                    emis_guess(s,i) = 1;
            end
        end
        emis_guess = emis_guess + 0.2 *rand(size(emis_guess));
        emis_guess = emis_guess./sum(emis_guess,2);
        trans_guess = ones(n_states,n_states)/n_stimuli + eye(n_states);

        trans_guess(:,end)  =trans_guess(:,end)+ 1/(1.0/params.binsize);

        trans_guess = trans_guess + 0.2 * rand(size(trans_guess));
        trans_guess = trans_guess./sum(trans_guess,2);
    case 'states_on_off'
        emis_guess = zeros(n_states,n_emissions);
        for i=1:n_emissions
            s = mode(stim_type(idx==i));
            switch s
                case -1
                    emis_guess(1,i) = 1;
                otherwise
                    emis_guess(2,i) = 1;
            end
        end
        emis_guess = emis_guess +  0.2 *rand(size(emis_guess));
        emis_guess = emis_guess./sum(emis_guess,2);

        trans_guess = [  1-1/(0.5/params.binsize) 1/(0.5/params.binsize);
            1/(1.0/params.binsize) 1-1/(1.0/params.binsize)];
        trans_guess = trans_guess .* (1 + 0.1* rand(size(trans_guess)));
        trans_guess = trans_guess./sum(trans_guess,2);
end

