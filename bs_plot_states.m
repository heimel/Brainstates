function bs_plot_states(bin_times,bin_states,xl)
%bs_plot_states. Plot states
%
% 2024, Alexander Heimel

if nargin<3 || isempty(xl)
    xl = [];
end

states = unique(bin_states);
n_states = length(states);
%clrs = parula(n_states);
clrs = 0.5*ones(n_states,3);
clrs(:,4) = 0.8; % add transparency

if ~isempty(xl)
    ind1 = find(bin_times<xl(1),1,'last');
    ind2 = find(bin_times>xl(2),1,'first');
    bin_times = bin_times(ind1:ind2);
    bin_states = bin_states(ind1:ind2);
end

n_bins = length(bin_times);

hold on

current_state = inf;
xleft = [];
xright = [];
for i=2:(n_bins-1)
    x1 = mean([bin_times(i-1) bin_times(i)]);
    x2 = mean([bin_times(i) bin_times(i+1)]);
    c = find(states==bin_states(i),1);
    if c~=current_state   
        if ~isempty(xleft)
            % current_state
            % area([xleft xright],current_state+[0.5 0.5],current_state-0.5,...
            %     'FaceColor',clrs(current_state,:),'EdgeColor','none');
            rectangle( 'Position',[xleft current_state-0.5 xright-xleft 1],...
                'FaceColor',clrs(current_state,:),'EdgeColor','none');
        end
        xleft = x1;
        xright = x2;
        current_state = c;
    else
        xright = x2;
    end
end
if ~isempty(xleft) 
    rectangle( 'Position',[xleft current_state-0.5 xright-xleft 1],...
        'FaceColor',clrs(current_state,:),'EdgeColor','none');
end
