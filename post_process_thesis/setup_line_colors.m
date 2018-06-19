N = 2; % how many curves/lines on figure
cust_map = brewermap([],'*Greys');close;
% cust_map = brewermap([],'*PuOr');close;
cmap_stop_pts = length(cust_map);
c_map_factor = 1; % starting from the darkest, what fraction is needed
cmap_new_stop_pts = floor(c_map_factor*cmap_stop_pts);
cust_map_new = cust_map(1:cmap_new_stop_pts,:);
line_colors_indices = round(linspace(1,cmap_new_stop_pts,N));
line_colors = cust_map_new(line_colors_indices,:);
run('custom_colors');
line_colors(2,:) = color_imp_coolgrey;
