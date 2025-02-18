n_top = 15;
n_bottom = n_top; % inner meshing won't work without n_top = n_bottom
% n_leading = 3;

k_inner = 6;
k_outer = 4;
mult = 1.3;
delta_inner = 0.2;
diamond_mult = 1.5;

% beta % spacing along airfoil;

slice_spacing = 4;

R_a = 800; % Bounding rect length in y (y from -R_a to R_a)
R_t = 400; % Bounding rect top
R_b = -300; % Bounding rect bottom (z from R_b to R_t)
R_x = 5500; % End cap position (-R_x and R_x)
