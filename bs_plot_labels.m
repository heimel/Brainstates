function bs_plot_labels(bin_times,bin_labels,yl,xl)
%bs_plot_labels. Adds labels as vertical bars
%
%  bs_plot_labels(bin_times,bin_labels,yl=[0 1],[xl])
%
% 2024, Alexander Heimel

if nargin<3 || isempty(yl)
    yl = [0 1];
end
if nargin<4 || isempty(xl)
    xl = [];
end

labels = unique(bin_labels);
n_labels = length(labels);
clrs = parula(length(labels));

if ~isempty(xl)
    ind1 = find(bin_times<xl(1),1,'last');
    ind2 = find(bin_times>xl(2),1,'first');
    bin_times = bin_times(ind1:ind2);
    bin_labels = bin_labels(ind1:ind2);
end
n_bins = length(bin_times);

hold on


current_label = inf;
xleft = [];
xright = [];
for i=2:(n_bins-1)
    x1 = mean([bin_times(i-1) bin_times(i)]);
    x2 = mean([bin_times(i) bin_times(i+1)]);
    c = find(labels==bin_labels(i),1);
    if c~=current_label   
        if ~isempty(xleft)
            area([xleft xright],[yl(2) yl(2)],yl(1),...
                'FaceColor',clrs(current_label,:),'EdgeColor','none');
        end
        xleft = x1;
        xright = x2;
        current_label = c;
    else
        xright = x2;
    end
end
if ~isempty(xleft)
    area([xleft xright],[yl(2) yl(2)],yl(1),...
        'FaceColor',clrs(current_label,:),'EdgeColor','none');
end

