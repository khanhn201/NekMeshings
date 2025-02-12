n_segs = 4;
k_inner = 9;
k_outer = 10;
k_upstream = 11;
k_downstream = 11;
k_blade = 6;

mult = 1.7;
delta_inner = 0.4;

Ra = 1250; % Inner square width (from -1250 to 1250)

% Turbine Blade region
R_t = 380; % Bounding rect top
R_b = -280; % Bounding rect bottom (z from R_b to R_t)
R_x = 5400; % End cap position (-R_x and R_x)

% Farfield
R_far_top = 30000
R_far_bottom = -10000
R_far_horizontal = 30000
R_downstream = 30000
R_upstream = 10000


