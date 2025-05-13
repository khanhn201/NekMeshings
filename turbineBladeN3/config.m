n_top = 8; % Need to be even to match the hub
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% n_leading = 3;

k_inner = 3;
k_outer = 5;
k_cyl = 3;
mult = 1.3;
first_layer_thickness = 0.1;
delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;
hub_radius = 1.78;
z_shift = 5;

slice_spacing = 6;
n_smooth = 100;

R_a = 4; % Bounding rect length in y (y from -R_a to R_a)
R_t = 4; % Bounding rect top
R_b = -4; % Bounding rect bottom (z from R_b to R_t)
R_end_caps = [65.5, 68]; % End cap position (-R_x and R_x)
