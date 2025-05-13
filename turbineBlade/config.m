n_top = 20;
n_top = 8;
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% n_leading = 3;

k_inner = 3;
k_outer = 10;
k_outer = 8;
k_cyl = 10;
mult = 1.3;
first_layer_thickness = 6;
delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;

slice_spacing = 4;
n_smooth = 100;

R_a = 800; % Bounding rect length in y (y from -R_a to R_a)
R_t = 600; % Bounding rect top
R_b = -600; % Bounding rect bottom (z from R_b to R_t)
R_end_caps = [5250, 5500]; % End cap position (-R_x and R_x)
