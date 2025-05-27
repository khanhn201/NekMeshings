n_top = 8; % Need to be even to match the hub
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% n_leading = 3;

k_inner = 5;
k_outer = 10;
k_cyl = 5;
k_downstream = 6;
mult = 1.3;
first_layer_thickness = 0.05;
delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;
hub_radius = 1.78;
z_shift = 3.5;

slice_spacing = 2;
n_slices = 80;
n_smooth = 100;

R_a = 5; % Bounding rect length in y (y from -R_a to R_a)
R_t = 5; % Bounding rect top
R_b = -5; % Bounding rect bottom (z from R_b to R_t)
R_downstream = -20;
R_downstream = 0;

R_end_caps = [61.6, 61.8, 62.5, 64, 67, 69]; % End cap position (-R_x and R_x)




% Coarse
% n_top = 4; % Need to be even to match the hub
% n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% % n_leading = 3;
%
% k_inner = 2;
% k_outer = 6;
% k_cyl = 5;
% k_downstream = 15;
% mult = 1.3;
% first_layer_thickness = 0.1;
% delta_inner = 0.2;
% diamond_mult = 1.5;
%
% % beta % spacing along airfoil;
% hub_radius = 1.78;
% z_shift = 3.5;
%
% slice_spacing = 2;
% n_slices = 40;
% n_smooth = 100;
%
% R_a = 5; % Bounding rect length in y (y from -R_a to R_a)
% R_t = 5; % Bounding rect top
% R_b = -5; % Bounding rect bottom (z from R_b to R_t)
% R_downstream = -20;
% R_downstream = 0;
% R_end_caps = [61.6, 61.8, 62.5, 64, 67, 69]; % End cap position (-R_x and R_x)

