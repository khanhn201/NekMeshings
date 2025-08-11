n_top = 48; % Need to be even to match the hub
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% n_leading = 3;

k_inner = 18;
k_outer = 20;
k_cyl = 30;
k_downstream = 15;

k_bl = 25;
mult = 1.08;
first_layer_thickness = 0.02;

delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;
hub_radius = 1.78;
hub_layers = 8;
z_shift = 3.5;

slice_spacing = 2;
n_slices = 240;
n_smooth = 100;

R_a = 5; % Bounding rect length in y (y from -R_a to R_a)
R_t = 5; % Bounding rect top
R_b = -5; % Bounding rect bottom (z from R_b to R_t)
% R_downstream = -80;
R_downstream = 0;

R_end_caps = [61.6, 61.8, 62.5, 63.5, 64.5, 65.5, 67, 68.5, 70, 71.5, 73, 74.5, 76]; % End cap position (-R_x and R_x)




% Coarse
n_top = 16; % Need to be even to match the hub
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom

k_inner = 12;
k_outer = 12;
k_cyl = 20;
k_downstream = 7;

k_bl = 15;
mult = 1.2;
first_layer_thickness = 0.01;

delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;
hub_radius = 1.8;
hub_layers = 5;
z_shift = 5;

slice_spacing = 2;
n_slices = 80;
n_smooth = 100;

R_a = 5; % Bounding rect length in y (y from -R_a to R_a)
R_t = 5; % Bounding rect top
R_b = -5; % Bounding rect bottom (z from R_b to R_t)
R_downstream = -80;
R_downstream = 0;
R_end_caps = [61.6, 63, 66, 70, 75]; % End cap position (-R_x and R_x)
