n_top = 24; % Need to be even to match the hub
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% n_leading = 3;

k_inner = 16;
k_outer = 40;
k_cyl = 35;
k_downstream = 15;
mult = 1.3;
first_layer_thickness = 0.005;
delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;
hub_radius = 1.78;
z_shift = 3.5;

slice_spacing = 2;
n_slices = 120;
n_smooth = 100;

R_a = 4; % Bounding rect length in y (y from -R_a to R_a)
R_t = 4; % Bounding rect top
R_b = -4; % Bounding rect bottom (z from R_b to R_t)
% R_downstream = -80;
R_downstream = 0;

R_end_caps = [61.6, 61.8, 62.5, 63.5, 64.5, 65.5, 67, 68.5, 70, 71.5, 73, 74.5, 76]; % End cap position (-R_x and R_x)




% Coarse
n_top = 24; % Need to be even to match the hub
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom

k_inner = 16;
k_outer = 20;
k_cyl = 10;
k_downstream = 7;
mult = 1.3;
first_layer_thickness = 0.005;
delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;
hub_radius = 1.78;
z_shift = 5;

slice_spacing = 2;
n_slices = 20;
n_smooth = 100;

R_a = 5; % Bounding rect length in y (y from -R_a to R_a)
R_t = 5; % Bounding rect top
R_b = -5; % Bounding rect bottom (z from R_b to R_t)
R_downstream = -80;
R_downstream = 0;
R_end_caps = [61.6, 63, 66, 70, 75]; % End cap position (-R_x and R_x)
